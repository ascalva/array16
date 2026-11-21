# Beam-controller firmware (RP2040)

Two layers, and a harness that runs them in an emulator:

| Path | What it is |
| --- | --- |
| `core/ctrl.c` | Protocol, controller state, serial stream assembly. No hardware. Its behaviour is fixed by `array16/fwsim.py` |
| `rp2040/sequencer.pio` | Two PIO state machines: `shifter` (74HC595 chain) and `ticker` (slot edges) |
| `rp2040/main.c` | DMA, clocks, USB CDC, HackRF trigger, USB-C current reading |
| `emu/harness.mjs` | Runs the built image in the rp2040js emulator and records the register pins (`emu/README.md`) |

`tests/test_firmware_core.py` builds `core/ctrl.c` with the host's C compiler and runs it
and the model on the same frames: thousands of random and malformed requests, split and
noisy serial streams, and a check that every control bit lands on its register output.
Plain `python -m pytest` runs it. `tests/test_firmware_emulated.py` runs the whole compiled image in
an emulator and checks replies, latched words and slot timing against the model (it skips if node,
rp2040js or the image is missing).

The RP2040 layer compiles cleanly (no warnings with `-Wall -Wextra`; Arm GNU Toolchain 14.3.rel1,
pico-sdk 2.3.1, 2026-11-21): 34.0 KB of code, 9.0 KB of RAM. A full sequence takes 5.7 KB of
that: the 256-word store, its copy packed for the DMA and the buffer a 1804-byte LOAD_SEQ frame
is assembled in. It has not run on an RP2040 yet, only in rp2040js (`emu/`).

## Build

Needs CMake, the Arm GNU Toolchain (`arm-none-eabi`, with newlib), pico-sdk 2.x with its TinyUSB
submodule, and picotool (pico-sdk builds it if it cannot find one; Raspberry Pi's prebuilt
`picotool-2.3.1-mac.zip` from github.com/raspberrypi/pico-sdk-tools saves that step).

    export PICO_SDK_PATH=/path/to/pico-sdk
    export PICO_TOOLCHAIN_PATH=/path/to/arm-gnu-toolchain-14.3.rel1-darwin-arm64-arm-none-eabi
    cmake -S firmware -B firmware/build -DPICO_BOARD=pico -Dpicotool_DIR=/path/to/picotool
    cmake --build firmware/build -j

On this machine the three live in `~/opt` (`~/opt/pico-sdk`, `~/opt/arm-gnu-toolchain-14.3.rel1-darwin-arm64-arm-none-eabi`,
`~/opt/picotool-2.3.1/picotool`).

Hold BOOTSEL, plug the board in and copy `firmware/build/array16_controller.uf2` to the
RPI-RP2 drive.

## How it works

- **Shift registers.** `shifter` clocks each 56-bit word out on SR_DATA (GPIO3) with
  SR_CLK (GPIO2) as side-set, most significant bit first, so bit 55 ends in the last
  register. 16 cycles per bit: SR_CLK at 7.8 MHz, 7.2 us per word at 125 MHz. It pulses
  SR_LATCH (GPIO4, 64 ns) when PIO IRQ flag 0 is set.
- **Slot clock.** The selected reference is divided to 1 MHz on GPOUT0 (GPIO21): 40 MHz / 40
  from GPIN0 (GPIO20), or the HackRF's 10 MHz / 10 from GPIN1 (GPIO22). `ticker` reads that
  pin back and sets IRQ 0 every `slot_us` edges, so slot edges never drift against the
  radios' sample clocks. GPIO21 must stay unconnected on the board. If the reference is
  missing (frequency counter more than 1 % off), the 1 MHz comes from the 12 MHz crystal
  instead and STATUS reports `ref_ok = 0`: sequences still play, unlocked. The selected
  reference is measured over and over (1 ms counts); if it goes away later, the slots stop until
  the count shows it (about 1 ms, at most 4.5 ms), then go on from the crystal, and `ref_ok` stays
  0 until CLOCK_SRC selects it again. A count gives up after 2.5 ms: FC0 may wait for ever for the
  first edge of a clock that is not there (FC0_STATUS.WAITING), which `frequency_count_khz()`
  would hang on.
- **Sequences.** DMA feeds the table (marker word, then the sequence, two FIFO words each:
  at most 514 transfers a period) into the shifter's joined 8-word FIFO, four slots of slack.
  The completion interrupt re-arms it once per period; after the last period of a finite
  sequence it queues the static word, which latches at the next slot edge.
- **RUN.** Whatever plays stops, the marker and the sequence are packed into the table, the DMA
  starts and TRIG (GPIO6) rises. The marker latches 11-12 us after TRIG rises (one shift plus
  margin), as TRIG falls, then one word per slot. In the emulator TRIG rises about 73 us after
  a RUN frame arrives with a few words loaded and 190 us with 256: packing takes 0.46 us a word.
- **SET_WORD, STOP, LOAD_SEQ while running.** Everything stops, and the static word is shifted
  in and latched at once, about 7 us after the firmware handles the frame. It handles a frame
  once it has read it and checked its CRC: in the emulator a STOP latches about 78 us after it
  arrives, a SET_WORD about 141 us, and a full LOAD_SEQ about 0.93 ms after its last packet,
  most of that the CRC over its 1802 bytes. Until then the previous sequence plays on.
- **Boot.** SR_OE_N is held high, zeros are shifted in and latched, then SR_OE_N goes low.
- **USB-C current.** ADC0/ADC1 read CC1/CC2 across the 5.1 kohm pull-downs. The Type-C sink
  thresholds (0.2, 0.66, 1.23 V) give 0, 500, 1500 or 3000 mA in STATUS.
- **Serial stream** (`protocol.py`): hunt for 0xA5; a length field above 1798 (a full LOAD_SEQ,
  6 + 256 x 7) is noise (resume one byte later, no reply); a partial frame is dropped after
  50 ms without a byte; a CRC failure is answered with command 0 and ERR_CRC.
- **USB.** `pico_stdio_usb`, binary safe (no newline translation, replies via `putchar_raw`).
  Opening the port at 1200 baud reboots the RP2040 into BOOTSEL (the SDK's reset hook), so
  the host must use any other rate. The host must also assert DTR (pyserial does): while DTR is
  low stdio_usb neither reads requests nor sends replies. Requests are read a byte at a time,
  about 9 us a byte in the emulator: a full LOAD_SEQ (1804 bytes, 29 USB packets) takes about
  16 ms to arrive, with at most 0.6 ms between packets, and the stream drops a partial frame
  only after 50 ms without a byte.

Limits: 1 to 256 words after the marker (protocol version 2, which PING reports); slots of at
least 10 us (the protocol's minimum; a shift takes 7.2 us). Changing CLOCK_SRC while a
sequence runs can stretch or shorten one slot.

## Bring-up on a bare Pico

The firmware runs unchanged on a Raspberry Pi Pico (same crystal, same flash family). With
nothing on GPIO20 it falls back to the crystal. Put a logic analyser on GPIO2-6 and GPIO21,
or wire up a 74HC595 chain with LEDs, and check:

- GPIO21 at 1.000 MHz; SR_CLK in bursts of 56 edges; one SR_LATCH pulse per word.
- A looping sequence at 50 us slots: latch pulses every 50.000 us, one blank word per period.
- STATUS after a sequence with `repeat = 2`: `running = 0`, and the static word back on the outputs.
- On the board: latch edges locked to CLK40 (trigger the scope on the reference).
