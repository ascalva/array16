/* array16 beam controller firmware (RP2040).
 *
 * What the controller does is in ../core/ctrl.c, which follows array16/fwsim.py. This file is
 * the hardware under it: the 74HC595 chain driven by PIO, slot timing locked to the selected
 * reference clock, USB CDC, the HackRF trigger and the USB-C current reading.
 *
 * Slot timing: the reference (40 MHz on GPIN0, or the HackRF's 10 MHz on GPIN1) is divided down
 * to 1 MHz on GPOUT0 (GPIO21, unconnected on the board) and the ticker state machine counts its
 * edges, so slot edges never drift against the radios' sample clocks. Without a reference, or
 * once the selected one has gone, the 1 MHz comes from the 12 MHz crystal and STATUS reports
 * ref_ok = 0.
 */
#include <string.h>

#include "hardware/adc.h"
#include "hardware/clocks.h"
#include "hardware/dma.h"
#include "hardware/irq.h"
#include "hardware/pio.h"
#include "pico/stdlib.h"

#include "ctrl.h"
#include "sequencer.pio.h"

/* Pins: hw/design.py, controller sheet. */
#define PIN_SR_CLK 2
#define PIN_SR_DATA 3
#define PIN_SR_LATCH 4
#define PIN_SR_OE_N 5
#define PIN_TRIG 6
#define PIN_REF40 20   /* GPIN0: on-board 40 MHz reference */
#define PIN_SLOTCLK 21 /* GPOUT0: 1 MHz slot clock, read back by the ticker */
#define PIN_REF10 22   /* GPIN1: HackRF 10 MHz */
#define PIN_CC1 26     /* ADC0 */
#define PIN_CC2 27     /* ADC1 */

#define SEQ_PIO pio0
#define FIRST_TICK_US 12 /* the marker latches this long after TRIG rises in hw_run: more than one shift */

static uint sm_shift, sm_tick, off_shift, off_tick;
static uint dma_ch;

static uint32_t table[2 * (CTRL_MAX_SEQ + 1)]; /* marker, then the sequence, packed for the shifter */
static uint32_t table_words;                    /* one period: the DMA transfer count, 2 (n + 1) <= 514 */
static uint32_t after_packed[2];
static volatile bool forever, finishing;
static volatile uint32_t periods_left;

/* ---- sequence playback ---------------------------------------------------------------------- */

/* One period of the table has gone into the shifter's FIFO, which still holds up to four words
 * (four slots, 40 us at the shortest slot), so there is time to queue the next. */
static void __isr dma_done(void) {
    dma_channel_acknowledge_irq0(dma_ch);
    if (finishing)
        return;
    if (forever || --periods_left)
        dma_channel_transfer_from_buffer_now(dma_ch, table, table_words);
    else {
        finishing = true; /* the static word follows the last period */
        dma_channel_transfer_from_buffer_now(dma_ch, after_packed, 2);
    }
}

/* Stop the ticker and the DMA and put the shifter back at its first instruction, empty. */
static void halt(void) {
    pio_sm_set_enabled(SEQ_PIO, sm_tick, false);
    dma_channel_set_irq0_enabled(dma_ch, false); /* RP2040-E13: an abort can post a stray IRQ */
    dma_channel_abort(dma_ch);
    dma_channel_acknowledge_irq0(dma_ch);
    dma_channel_set_irq0_enabled(dma_ch, true);
    forever = false;
    finishing = true;
    pio_sm_set_enabled(SEQ_PIO, sm_shift, false);
    pio_sm_clear_fifos(SEQ_PIO, sm_shift);
    pio_sm_clear_fifos(SEQ_PIO, sm_tick);
    pio_sm_restart(SEQ_PIO, sm_shift);
    pio_sm_restart(SEQ_PIO, sm_tick);
    pio_sm_exec(SEQ_PIO, sm_shift, pio_encode_set(pio_pins, 0)); /* SR_LATCH low, in case we cut a pulse */
    pio_sm_exec(SEQ_PIO, sm_shift, pio_encode_jmp(off_shift));
    pio_sm_exec(SEQ_PIO, sm_tick, pio_encode_jmp(off_tick));
    SEQ_PIO->irq = 1u << 0; /* drop a pending slot tick */
    pio_sm_set_enabled(SEQ_PIO, sm_shift, true);
}

static void hw_latch(void *ctx, const uint8_t *word) {
    (void)ctx;
    uint32_t w[2];
    halt();
    ctrl_pack(word, w);
    pio_sm_put(SEQ_PIO, sm_shift, w[0]);
    pio_sm_put(SEQ_PIO, sm_shift, w[1]);
    SEQ_PIO->irq_force = 1u << 0; /* latch as soon as it is shifted in */
}

static void hw_run(void *ctx, const uint8_t *words, uint16_t n, uint16_t slot_us, uint16_t repeat,
                   const uint8_t *after) {
    (void)ctx;
    static const uint8_t blank[CTRL_WORD_BYTES];
    halt();
    ctrl_pack(blank, table);
    for (uint i = 0; i < n; i++)
        ctrl_pack(words + i * CTRL_WORD_BYTES, table + 2 + 2 * i);
    table_words = 2u * (n + 1u);
    ctrl_pack(after, after_packed);
    forever = repeat == 0;
    periods_left = repeat;
    finishing = false;
    pio_sm_put(SEQ_PIO, sm_tick, FIRST_TICK_US - 1);
    pio_sm_put(SEQ_PIO, sm_tick, slot_us - 1u);
    dma_channel_transfer_from_buffer_now(dma_ch, table, table_words);
    gpio_put(PIN_TRIG, 1); /* HackRF trigger: rises FIRST_TICK_US before the marker slot */
    pio_sm_set_enabled(SEQ_PIO, sm_tick, true);
    busy_wait_us_32(FIRST_TICK_US);
    gpio_put(PIN_TRIG, 0);
}

/* Done once the static word has been queued, shifted and latched, and the shifter waits for more. */
static bool hw_playing(void *ctx) {
    (void)ctx;
    if (!finishing)
        return true;
    return dma_channel_is_busy(dma_ch) || !pio_sm_is_tx_fifo_empty(SEQ_PIO, sm_shift) ||
           pio_sm_get_pc(SEQ_PIO, sm_shift) != off_shift;
}

/* ---- clocks and USB power --------------------------------------------------------------------- */

/* The reference is measured with the frequency counter FC0. Not with frequency_count_khz(): that
 * waits for DONE, and with no clock on the input the counter may wait for a first edge
 * (FC0_STATUS.WAITING) that never comes. A measurement that has not finished by its deadline
 * means no reference. */
#define FC_INTERVAL 10       /* 0.98 us x 2^10: about 1 ms */
#define FC_DEADLINE_US 2500

static uint8_t ref_src; /* the reference last selected */
static bool fc_busy;    /* a measurement of it is under way */
static absolute_time_t fc_deadline;

static void fc_start(void) {
    fc_hw_t *fc = &clocks_hw->fc0;
    const absolute_time_t t = make_timeout_time_us(FC_DEADLINE_US);
    while ((fc->status & CLOCKS_FC0_STATUS_RUNNING_BITS) && !time_reached(t))
        tight_loop_contents();
    fc->ref_khz = clock_get_hz(clk_ref) / 1000;
    fc->interval = FC_INTERVAL;
    fc->min_khz = 0;
    fc->max_khz = 0xffffffff;
    fc->src = ref_src ? CLOCKS_FC0_SRC_VALUE_CLKSRC_GPIN1 : CLOCKS_FC0_SRC_VALUE_CLKSRC_GPIN0; /* starts it */
    fc_deadline = make_timeout_time_us(FC_DEADLINE_US);
    fc_busy = true;
}

/* False once the measurement is over; *ok then says whether the reference is there, within 1 %. */
static bool fc_measuring(bool *ok) {
    const fc_hw_t *fc = &clocks_hw->fc0;
    const bool done = fc->status & CLOCKS_FC0_STATUS_DONE_BITS;
    if (!done && !time_reached(fc_deadline))
        return true;
    const uint32_t want = ref_src ? 10000 : 40000, got = done ? fc->result >> CLOCKS_FC0_RESULT_KHZ_LSB : 0;
    *ok = got + want / 100 > want && got < want + want / 100;
    fc_busy = false;
    return false;
}

/* The 1 MHz slot clock on GPOUT0: the selected reference divided down, or the 12 MHz crystal. */
static void slot_clock(bool locked) {
    if (locked)
        clock_gpio_init(PIN_SLOTCLK,
                        ref_src ? CLOCKS_CLK_GPOUT0_CTRL_AUXSRC_VALUE_CLKSRC_GPIN1 : CLOCKS_CLK_GPOUT0_CTRL_AUXSRC_VALUE_CLKSRC_GPIN0,
                        ref_src ? 10.0f : 40.0f);
    else
        clock_gpio_init(PIN_SLOTCLK, CLOCKS_CLK_GPOUT0_CTRL_AUXSRC_VALUE_XOSC_CLKSRC, 12.0f);
}

static bool hw_clock(void *ctx, uint8_t src) {
    (void)ctx;
    gpio_set_function(src ? PIN_REF10 : PIN_REF40, GPIO_FUNC_GPCK);
    ref_src = src;
    fc_start();
    bool ok;
    while (fc_measuring(&ok))
        tight_loop_contents();
    slot_clock(ok);
    return ok;
}

/* The selected reference can also go away later (the HackRF unplugged or switched off): GPOUT0
 * would stop and the slots with it, while STATUS still said ref_ok. So it is measured over and
 * over; once it has gone, the slots run from the crystal as if it had been missing at CLOCK_SRC,
 * and STATUS reports ref_ok = 0. CLOCK_SRC locks on to it again. */
static void watch_reference(ctrl_t *c) {
    bool ok = true;
    if (!c->ref_ok || (fc_busy && fc_measuring(&ok)))
        return;
    if (ok)
        fc_start();
    else {
        slot_clock(false);
        c->ref_ok = false;
    }
}

static uint32_t cc_mv(uint input) {
    adc_select_input(input);
    return adc_read() * 3300u / 4095u;
}

/* The voltage across the board's 5.1k CC pull-down tells what the source offers (USB Type-C
 * sink thresholds: 0.2 V attached, 0.66 V 1.5 A, 1.23 V 3 A). Only one CC pin carries it. */
static uint16_t hw_usb_ma(void *ctx) {
    (void)ctx;
    const uint32_t a = cc_mv(0), b = cc_mv(1), mv = a > b ? a : b;
    return mv < 200 ? 0 : mv < 660 ? 500 : mv < 1230 ? 1500 : 3000;
}

/* ---- start-up ----------------------------------------------------------------------------------- */

static void sequencer_init(void) {
    off_shift = (uint)pio_add_program(SEQ_PIO, &shifter_program);
    off_tick = (uint)pio_add_program(SEQ_PIO, &ticker_program);
    sm_shift = (uint)pio_claim_unused_sm(SEQ_PIO, true);
    sm_tick = (uint)pio_claim_unused_sm(SEQ_PIO, true);

    pio_sm_config c = shifter_program_get_default_config(off_shift);
    sm_config_set_out_pins(&c, PIN_SR_DATA, 1);
    sm_config_set_set_pins(&c, PIN_SR_LATCH, 1);
    sm_config_set_sideset_pins(&c, PIN_SR_CLK);
    sm_config_set_out_shift(&c, false, false, 32); /* most significant bit first, explicit pulls */
    sm_config_set_fifo_join(&c, PIO_FIFO_JOIN_TX); /* 8 words deep: four control words */
    const uint32_t pins = 1u << PIN_SR_CLK | 1u << PIN_SR_DATA | 1u << PIN_SR_LATCH;
    pio_sm_set_pins_with_mask(SEQ_PIO, sm_shift, 0, pins);
    pio_sm_set_pindirs_with_mask(SEQ_PIO, sm_shift, pins, pins);
    pio_gpio_init(SEQ_PIO, PIN_SR_CLK);
    pio_gpio_init(SEQ_PIO, PIN_SR_DATA);
    pio_gpio_init(SEQ_PIO, PIN_SR_LATCH);
    pio_sm_init(SEQ_PIO, sm_shift, off_shift, &c);

    pio_sm_config t = ticker_program_get_default_config(off_tick);
    sm_config_set_in_pins(&t, PIN_SLOTCLK);
    pio_sm_init(SEQ_PIO, sm_tick, off_tick, &t);

    dma_ch = (uint)dma_claim_unused_channel(true);
    dma_channel_config d = dma_channel_get_default_config(dma_ch);
    channel_config_set_transfer_data_size(&d, DMA_SIZE_32);
    channel_config_set_read_increment(&d, true);
    channel_config_set_write_increment(&d, false);
    channel_config_set_dreq(&d, pio_get_dreq(SEQ_PIO, sm_shift, true));
    dma_channel_configure(dma_ch, &d, &SEQ_PIO->txf[sm_shift], table, 0, false);
    finishing = true;
    dma_channel_set_irq0_enabled(dma_ch, true);
    irq_set_exclusive_handler(DMA_IRQ_0, dma_done);
    irq_set_enabled(DMA_IRQ_0, true);
    pio_sm_set_enabled(SEQ_PIO, sm_shift, true);
}

int main(void) {
    /* The registers' outputs stay off (SR_OE_N is also pulled high on the board) until the safe
     * all-zeros word is in them. */
    gpio_init(PIN_SR_OE_N);
    gpio_put(PIN_SR_OE_N, 1);
    gpio_set_dir(PIN_SR_OE_N, GPIO_OUT);
    gpio_init(PIN_TRIG);
    gpio_put(PIN_TRIG, 0);
    gpio_set_dir(PIN_TRIG, GPIO_OUT);
    adc_init();
    adc_gpio_init(PIN_CC1);
    adc_gpio_init(PIN_CC2);

    sequencer_init();
    static const uint8_t zeros[CTRL_WORD_BYTES];
    hw_latch(NULL, zeros);
    busy_wait_us_32(50); /* one shift is 7.2 us */
    gpio_put(PIN_SR_OE_N, 0);

    static ctrl_t ctrl;
    const ctrl_hw_t hw = {hw_latch, hw_run, hw_playing, hw_clock, hw_usb_ma, NULL};
    ctrl_init(&ctrl, &hw);
    stdio_init_all(); /* USB CDC only, see CMakeLists.txt; replies go out with putchar_raw */

    uint8_t reply[CTRL_MAX_REPLY];
    for (;;) {
        watch_reference(&ctrl);
        const int ch = getchar_timeout_us(1000);
        if (ch < 0)
            continue;
        const size_t n = ctrl_feed(&ctrl, (uint8_t)ch, to_ms_since_boot(get_absolute_time()), reply);
        for (size_t i = 0; i < n; i++)
            putchar_raw(reply[i]);
        if (n)
            stdio_flush();
    }
}
