#include "ctrl.h"

#include <string.h>

_Static_assert(CTRL_MAX_FRAME <= UINT16_MAX, "frame lengths and ctrl_t.have are u16");
_Static_assert(CTRL_MAX_SEQ <= UINT16_MAX, "the word count is u16");

uint16_t ctrl_crc16(const uint8_t *d, size_t n) {
    uint16_t crc = 0xFFFF;
    while (n--) {
        crc ^= (uint16_t)(*d++ << 8);
        for (int i = 0; i < 8; i++)
            crc = (crc & 0x8000) ? (uint16_t)((crc << 1) ^ 0x1021) : (uint16_t)(crc << 1);
    }
    return crc;
}

static uint16_t rd16(const uint8_t *p) { return (uint16_t)(p[0] | p[1] << 8); }

static size_t reply(uint8_t *out, uint8_t cmd, uint8_t status, const uint8_t *payload, uint16_t n) {
    out[0] = CTRL_SOF_REP;
    out[1] = cmd;
    out[2] = status;
    out[3] = (uint8_t)n;
    out[4] = (uint8_t)(n >> 8);
    if (n)
        memcpy(out + 5, payload, n);
    uint16_t crc = ctrl_crc16(out, 5u + n);
    out[5 + n] = (uint8_t)crc;
    out[6 + n] = (uint8_t)(crc >> 8);
    return 7u + n;
}

void ctrl_init(ctrl_t *c, const ctrl_hw_t *hw) {
    memset(c, 0, sizeof *c);
    c->hw = *hw;
    c->ref_ok = c->hw.clock(c->hw.ctx, 0);
}

size_t ctrl_handle(ctrl_t *c, const uint8_t *f, size_t len, uint8_t *out) {
    if (c->running && !c->hw.playing(c->hw.ctx))
        c->running = false; /* a sequence with a repeat count has ended */
    if (len < 6 || f[0] != CTRL_SOF_REQ)
        return reply(out, 0, ST_ERR_CRC, NULL, 0);
    uint16_t n = rd16(f + 2);
    if (len != 6u + n || ctrl_crc16(f, len - 2) != rd16(f + 4 + n))
        return reply(out, 0, ST_ERR_CRC, NULL, 0);
    const uint8_t cmd = f[1], *p = f + 4;

    switch (cmd) {
    case CMD_PING: {
        const uint8_t v = CTRL_FW_VERSION;
        return reply(out, cmd, ST_OK, &v, 1);
    }
    case CMD_SET_WORD:
        if (n != CTRL_WORD_BYTES)
            break;
        memcpy(c->word, p, CTRL_WORD_BYTES);
        c->running = false;
        c->hw.latch(c->hw.ctx, c->word);
        return reply(out, cmd, ST_OK, NULL, 0);
    case CMD_LOAD_SEQ: {
        if (n < CTRL_SEQ_HEAD)
            break;
        const uint16_t slot = rd16(p), rep = rd16(p + 2), count = rd16(p + 4);
        if (n != CTRL_SEQ_HEAD + count * CTRL_WORD_BYTES || count < 1 || count > CTRL_MAX_SEQ ||
            slot < CTRL_MIN_SLOT_US)
            break;
        if (c->running) {
            c->running = false;
            c->hw.latch(c->hw.ctx, c->word);
        }
        memcpy(c->seq, p + CTRL_SEQ_HEAD, (size_t)count * CTRL_WORD_BYTES);
        c->seq_n = count;
        c->slot_us = slot;
        c->repeat = rep;
        return reply(out, cmd, ST_OK, NULL, 0);
    }
    case CMD_RUN:
        if (!c->seq_n)
            break;
        c->running = true;
        c->hw.run(c->hw.ctx, c->seq, c->seq_n, c->slot_us, c->repeat, c->word);
        return reply(out, cmd, ST_OK, NULL, 0);
    case CMD_STOP:
        c->running = false;
        c->hw.latch(c->hw.ctx, c->word);
        return reply(out, cmd, ST_OK, NULL, 0);
    case CMD_STATUS: {
        const uint16_t ma = c->hw.usb_ma(c->hw.ctx);
        const uint8_t s[CTRL_STATUS_BYTES] = {CTRL_FW_VERSION,     c->running,
                                              (uint8_t)c->seq_n,   (uint8_t)(c->seq_n >> 8),
                                              (uint8_t)c->slot_us, (uint8_t)(c->slot_us >> 8),
                                              c->clock_src,        c->ref_ok,
                                              (uint8_t)ma,         (uint8_t)(ma >> 8)};
        return reply(out, cmd, ST_OK, s, sizeof s);
    }
    case CMD_CLOCK_SRC:
        if (n != 1 || p[0] > 1)
            break;
        c->clock_src = p[0];
        c->ref_ok = c->hw.clock(c->hw.ctx, p[0]);
        return reply(out, cmd, ST_OK, NULL, 0);
    default:
        return reply(out, cmd, ST_ERR_CMD, NULL, 0);
    }
    return reply(out, cmd, ST_ERR_ARG, NULL, 0);
}

/* Drop the first buffered byte and anything up to the next start byte. */
static void resync(ctrl_t *c) {
    uint16_t i = 1;
    while (i < c->have && c->buf[i] != CTRL_SOF_REQ)
        i++;
    memmove(c->buf, c->buf + i, c->have - i);
    c->have = (uint16_t)(c->have - i);
}

size_t ctrl_feed(ctrl_t *c, uint8_t b, uint32_t now_ms, uint8_t *out) {
    if (c->have && now_ms - c->last_ms > CTRL_FRAME_TIMEOUT_MS)
        c->have = 0; /* a partial frame went stale */
    c->last_ms = now_ms;
    if (!c->have && b != CTRL_SOF_REQ)
        return 0;
    c->buf[c->have++] = b;
    while (c->have >= 4 && rd16(c->buf + 2) > CTRL_MAX_PAYLOAD)
        resync(c); /* no request is that long: that start byte was noise */
    if (c->have >= 4 && c->have == 6u + rd16(c->buf + 2)) {
        size_t r = ctrl_handle(c, c->buf, c->have, out);
        c->have = 0;
        return r;
    }
    return 0;
}

void ctrl_pack(const uint8_t *w, uint32_t fifo[2]) {
    uint64_t v = 0;
    for (int i = CTRL_WORD_BYTES - 1; i >= 0; i--)
        v = v << 8 | w[i];
    fifo[0] = (uint32_t)(v >> 24);
    fifo[1] = (uint32_t)((v & 0xFFFFFFu) << 8);
}

size_t ctrl_size(void) { return sizeof(ctrl_t); }
