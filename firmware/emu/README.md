# The firmware in an emulator

`harness.mjs` boots the compiled firmware in [rp2040js](https://github.com/wokwi/rp2040js)
(Wokwi's RP2040 emulator), talks to it over emulated USB CDC, records SR_CLK, SR_DATA,
SR_LATCH, SR_OE_N and TRIG with emulated timestamps (pin numbers are read from
`../rp2040/main.c`), decodes the 74HC595 chain into latched 56-bit words and prints JSON lines.
`tests/test_firmware_emulated.py` builds scripts from `array16/protocol.py`, runs them and checks
the result against `array16/fwsim.py`. It skips if node, rp2040js or the image is missing.

## Running it

    git clone https://github.com/wokwi/rp2040js ~/opt/src/rp2040js && (cd ~/opt/src/rp2040js && npm install)
    # build the firmware (../README.md); the image is looked for in firmware/build/
    node firmware/emu/harness.mjs --demo             # PING, STATUS, SET_WORD, a short sequence
    node firmware/emu/harness.mjs --script s.json    # or --script - to read stdin
    python -m pytest tests/test_firmware_emulated.py # about 2 s

Options: `--image` (`.uf2`, `.hex` or `.bin`; default `firmware/build/array16_controller.uf2`,
or `$ARRAY16_FW_IMAGE`), `--rp2040js DIR` (or `$RP2040JS`), `--trace` (every pin edge and USB
packet).

A script is `{"scenarios": [...]}`; each scenario boots a fresh chip, waits for USB, then runs
its steps in order:

| Step | Does |
| --- | --- |
| `{"send": "a5...", "packets": [3, 10], "gap_us": 2000}` | queue bytes for the OUT endpoint, split into those USB packets, waiting `gap_us` after each one has arrived (both optional; packets over 64 bytes are split further) |
| `{"wait_us": 500}` | run for that long |
| `{"wait_replies": 4, "timeout_us": 50000}` | run until that many reply frames have come back |
| `{"env": {"ref10_hz": 0}}` | change the board: `ref40_hz`, `ref10_hz` (0 = absent), `cc1_mv`, `cc2_mv`, `fc_dead_clock` |
| `{"dtr": false}` | SET_CONTROL_LINE_STATE |

Events (`ev`): `latch` (word as 14 hex digits, bit 55 first; SR_CLK edges since the last latch;
SR_OE_N), `pin` (SR_OE_N, TRIG), `out` (a packet the device took), `reply` (a whole reply frame,
with the time of its first byte), `gpout` (slot clock source and frequency), `fc` (a frequency
measurement), `env`, `timing` (narrowest pulses and setup times the registers saw), `warnings`
(registers the firmware touched that rp2040js lacks), `end`.

## What the harness adds to rp2040js

The firmware leans on parts rp2040js leaves out or simplifies. The harness supplies them from
outside, without changing rp2040js:

- **PIO timing.** rp2040js runs the PIO in batches on wall-clock timers and ignores delay cycles
  and the clock divider. The harness steps it once per clk_sys cycle in emulated time, with
  delays and the divider. Without that there are no meaningful slot times at all.
- **PIO FIFO join and clear.** rp2040js ignores FJOIN, so the shifter's FIFO would be 4 deep
  instead of 8, and `pio_sm_clear_fifos()` (which toggles FJOIN_RX) would clear nothing: STOP
  during a sequence would then latch a stale sequence word.
- **Frequency counter.** Absent in rp2040js (`frequency_count_khz()` would spin for ever).
  Modelled with the scenario's reference frequencies and the real 0.98 us x 2^interval count time.
  What the counter does with a clock that never starts is not documented (FC0_STATUS has a
  WAITING state); `fc_dead_clock` picks either answer: `done` (0 kHz, DIED) or `wait` (never DONE).
- **GPOUT0 and GPIN0/1.** The registers exist but no clock comes out. Modelled: the selected
  source divided by DIV toggles GPIO21's input, which the ticker state machine counts. A GPIN
  carries its reference only once its pin is set to the clock function. A source change keeps
  the output level and starts a new half period: no glitches.
- **USB host.** rp2040js's CDC host answers every armed OUT endpoint at once, with an empty
  packet if it has nothing. Here an armed endpoint waits for data and the script decides where
  packets split. Each packet takes rp2040js's fixed 10 us.
- **ADC.** CC1/CC2 voltages come from the scenario.

## What it shows and what it does not

Shown: the firmware boots, enumerates and answers every frame as fwsim does; the words reach
the registers bit for bit; the PIO programs, DMA re-arming, interrupts and the slot logic work
together, and slot edges land on the model's boundaries to the nanosecond (the emulated slot
clock is exact and the PIO reacts in a fixed number of cycles); the register pins meet 64 ns
setup and pulse widths at 125 MHz.

Not shown: real CPU timing (rp2040js counts about one cycle per instruction, no flash or bus
waits), so latencies such as "frame received to latch" (75-140 us here) are estimates; DMA is
instantaneous and paced by a simplified DREQ; no GPIO input synchronisers (real edges are seen
2-3 cycles later); glitches when GPOUT0 changes source; what FC0 really does with a dead clock;
erratum RP2040-E13 (rp2040js's DMA abort never raises an interrupt); real USB hosts (packet
timing, DTR, the 1200-baud reboot); anything electrical.
