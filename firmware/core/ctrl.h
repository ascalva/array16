/* Beam-controller core: protocol, controller state and stream assembly. No hardware here.
 *
 * The behaviour is specified by array16/fwsim.py (frames) and array16/protocol.py (stream).
 * tests/test_firmware_core.py builds this file for the host and runs it and the model on the
 * same frames. The RP2040 layer (firmware/rp2040/main.c) supplies the hardware calls.
 */
#ifndef ARRAY16_CTRL_H
#define ARRAY16_CTRL_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define CTRL_FW_VERSION 2 /* 2: LOAD_SEQ and STATUS carry the word count as u16 */
#define CTRL_WORD_BYTES 7
#define CTRL_MAX_SEQ 256  /* words after the marker */
#define CTRL_SEQ_HEAD 6   /* LOAD_SEQ: slot_us u16, repeat u16, count u16, then the words */
#define CTRL_MAX_PAYLOAD (CTRL_SEQ_HEAD + CTRL_MAX_SEQ * CTRL_WORD_BYTES) /* a full LOAD_SEQ: 1798 */
#define CTRL_MAX_FRAME (6 + CTRL_MAX_PAYLOAD)
#define CTRL_STATUS_BYTES 10
#define CTRL_MAX_REPLY (7 + CTRL_STATUS_BYTES)
#define CTRL_FRAME_TIMEOUT_MS 50
#define CTRL_MIN_SLOT_US 10

enum { CTRL_SOF_REQ = 0xA5, CTRL_SOF_REP = 0x5A };
enum {
    CMD_PING = 0x01,
    CMD_SET_WORD = 0x10,
    CMD_LOAD_SEQ = 0x11,
    CMD_RUN = 0x12,
    CMD_STOP = 0x13,
    CMD_STATUS = 0x20,
    CMD_CLOCK_SRC = 0x30,
};
enum { ST_OK = 0, ST_ERR_CRC = 1, ST_ERR_CMD = 2, ST_ERR_ARG = 3 };

/* What the core asks of the hardware. Every call comes from the main loop, none from an interrupt.
 * Words are 7 bytes, little endian: bit i of the control word is bit (i % 8) of byte i / 8. */
typedef struct {
    /* Stop any sequence and latch this word at once. */
    void (*latch)(void *ctx, const uint8_t *word);
    /* Latch the blank marker word (all zeros) at once, then words[0..n-1] one slot each, and
     * repeat. After `repeat` periods (0 = until stopped) latch `after` and stop. */
    void (*run)(void *ctx, const uint8_t *words, uint16_t n, uint16_t slot_us, uint16_t repeat,
                const uint8_t *after);
    /* True while the sequence started by run() is still playing. */
    bool (*playing)(void *ctx);
    /* Time the slots from reference `src` (0 = on-board 40 MHz, 1 = HackRF 10 MHz).
     * Returns false if that reference is absent and the controller fell back to its own crystal. */
    bool (*clock)(void *ctx, uint8_t src);
    /* Current the USB-C source advertises on CC, mA: 0, 500, 1500 or 3000. */
    uint16_t (*usb_ma)(void *ctx);
    void *ctx;
} ctrl_hw_t;

typedef struct {
    ctrl_hw_t hw;
    uint8_t word[CTRL_WORD_BYTES]; /* static word: latched whenever no sequence plays */
    uint8_t seq[CTRL_MAX_SEQ * CTRL_WORD_BYTES];
    uint16_t seq_n; /* words loaded, 1 to CTRL_MAX_SEQ; 0 = nothing loaded */
    uint16_t slot_us, repeat;
    bool running;
    uint8_t clock_src;
    bool ref_ok;
    uint8_t buf[CTRL_MAX_FRAME]; /* request being assembled from the stream */
    uint16_t have;               /* at most CTRL_MAX_FRAME */
    uint32_t last_ms;
} ctrl_t;

/* Power-up state: static word all zeros, nothing loaded, slots timed from the on-board reference. */
void ctrl_init(ctrl_t *c, const ctrl_hw_t *hw);

/* Handle one complete request frame, exactly as fwsim.ControllerModel.handle does.
 * Writes the reply to `reply` (CTRL_MAX_REPLY bytes) and returns its length. */
size_t ctrl_handle(ctrl_t *c, const uint8_t *frame, size_t len, uint8_t *reply);

/* Feed one byte from the serial stream. Returns the length of a reply written to `reply`
 * (CTRL_MAX_REPLY bytes) when the byte completes a frame, else 0. */
size_t ctrl_feed(ctrl_t *c, uint8_t byte, uint32_t now_ms, uint8_t *reply);

uint16_t ctrl_crc16(const uint8_t *data, size_t n); /* CRC-16/CCITT-FALSE */

/* A control word as two 32-bit words for the PIO shifter, which sends the most significant bit
 * first: bits 55..24, then bits 23..0 in the top 24 bits. Bit 55 enters the chain first and so
 * ends in the last register. */
void ctrl_pack(const uint8_t *word, uint32_t fifo[2]);

size_t ctrl_size(void); /* sizeof(ctrl_t), for the host tests */

#endif
