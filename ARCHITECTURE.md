# Architecture

## Layers

    experiments, tests, your scripts
            |
    api.Array            calibrate / point / snapshot / locate / monopulse / calibrate_tx
            |
    control, dsp         which words to play; what to do with the samples
            |
    backend.Backend      controller(frame) -> reply    tune(freq)    capture(n) -> samples
       |         |       four-channel: cal_signal(on) -> code replica, set_tx(weights)
    SimBackend   HardwareBackend (to write): serial port to the RP2040, pyadi-iio or HackRF
       |
    fwsim.ControllerModel + BoardModel + radios + scene

`board.py` sits under everything. The schematic generator imports the same control-word
bit positions, so the firmware model, the simulator and the copper cannot drift apart
silently (`hw/design.py` asserts the mapping).

## Control word

56 bits, 7 x 74HC595. Bit i is output Q(i mod 8) of register i div 8. All zeros is the
power-up state: every LNA off, every element in transmit bypass, nothing transmitting.

| Bits | Name | Meaning when 1 |
| --- | --- | --- |
| 0-15 | RX[e] | element e receives through its LNA |
| 16-31 | P180[e] | long path of the 180-degree stage |
| 32-47 | P90[e] | long path of the 90-degree stage |
| 48 | MODE | quadrants feed the shared combiner (HackRF mode) |
| 49-52 | QRX[q] | quadrant port receives |
| 53 | HRX | HackRF port receives |
| 54 | CAL | CAL input routed to the centre probe |
| 55 | LNA2_EN | second-stage LNAs on |

Element e = 4*row + col, row 0 at the top and col 0 at the left seen from the front.
Quadrants: 0 top-left, 1 top-right, 2 bottom-left, 3 bottom-right.

## Controller protocol

USB CDC, little endian, CRC-16/CCITT-FALSE. Request `A5 cmd len16 payload crc16`,
reply `5A cmd status len16 payload crc16`.

| Cmd | Name | Payload |
| --- | --- | --- |
| 01 | PING | none; reply carries the firmware version |
| 10 | SET_WORD | 7 bytes; latched at once, stops any sequence |
| 11 | LOAD_SEQ | slot_us u16, repeat u16 (0 = until STOP), count u8, then count x 7 bytes |
| 12 | RUN | start the loaded sequence |
| 13 | STOP | |
| 20 | STATUS | reply (9 bytes, `protocol.Status`): version u8, running u8, count u8, slot_us u16, clock_src u8, ref_ok u8, usb_ma u16 |
| 30 | CLOCK_SRC | 0 = on-board 40 MHz, 1 = HackRF 10 MHz |

A running sequence plays one blank marker slot (word 0) and then its words, and repeats.
The marker is a power dip in the samples; the host finds it to learn where the slots are.
A sequence with `repeat = n` ends after n periods exactly as if STOP had come: the static
word (the last SET_WORD) is latched again and STATUS reports `running = 0`.

STATUS also says whether the selected reference clock is present (`ref_ok`; if not, slots are
timed from the controller's own crystal and drift against the radios) and how much current
the USB-C source offers (`usb_ma`; both radio boards need 1500).

On the serial stream the controller hunts for 0xA5, treats a length above 453 bytes (a full
LOAD_SEQ) as noise and resumes one byte later without replying, drops a partial frame after
50 ms of silence, and answers a CRC failure with command 0 and ERR_CRC.

### Firmware (RP2040)

Compiles without warnings and runs in the rp2040js emulator, where 20 scenarios match `fwsim.py` byte for byte and latch for latch (`tests/test_firmware_emulated.py`, harness in `firmware/emu/`); not yet run on a real RP2040. See `firmware/README.md`.

- `firmware/core/ctrl.c` holds everything `fwsim.py` specifies and is tested against it.
- PIO state machine `shifter` shifts 56 bits into the registers (7.2 us) and pulses the latch
  when the `ticker` state machine signals a slot edge.
- The ticker counts a 1 MHz slot clock made by dividing the 40 MHz reference on GPIO20 (or the
  HackRF's 10 MHz on GPIO22) onto GPOUT0, GPIO21, and reading it back. GPIO21 must stay unconnected.
- GPIO6 pulses the HackRF trigger at the start of a sequence, 11-12 us before the marker slot; the marker latches that long after
  RUN (the model latches it at once).
- The selected reference is measured continuously with FC0 (2.5 ms deadline, so a missing clock cannot hang the firmware). If it
  disappears, the slot clock falls back to the 12 MHz crystal within a few milliseconds and STATUS reports ref_ok = 0; sequences keep
  playing. It does not re-lock by itself: the host sends CLOCK_SRC again.
- The host must assert DTR (pyserial does by default): without it the USB stdio neither reads requests nor sends replies.
- At boot: shift zeros, latch, then enable the register outputs (SR_OE_N low).
- ADC0/ADC1 read the USB-C CC pins; STATUS reports the advertised current.

## Snapshots

Patterns invert chosen elements by 180 degrees (Hadamard codes). One complex number per
element comes out of the decode.

The "180-degree" state is a delay line, so it is 180 degrees only at 5.8 GHz: at f it is
180 * f / 5.8 GHz. Off 5.8 GHz the decode uses the actual step (`control.toggle_factor`) and
solves the patterns instead of inverting the Hadamard matrix, at a noise cost (four radios)
of about 2.3 dB at 3.6 GHz, 7 dB at 2 GHz and 13 dB at 1 GHz. Beam steering likewise picks each element's state
by its actual delay (`control.state_lag`). At 5.8 GHz both are exactly the original code.

- HackRF mode: 16 patterns. The signal must be a steady tone; its frequency is estimated
  from the slots and removed.
- Four radios, steady signal: 4 patterns, each quadrant stepping through H4.
- Four radios, any signal: 5 patterns (`BLIND_SCHEDULE`). In every pair of neighbouring
  slots one quadrant holds pattern 0, which ties the unknown signal level of one slot to
  the next. Each slot's four-channel covariance gives the quadrant outputs up to that level.

## Calibration

1. `calibrate()` with a steady source straight ahead and at least 2 m away. Stores the raw
   snapshot as the element reference. In four-radio mode it first measures the channels (step 2)
   and stores that as the channel reference.
2. `align()`, called by `tune()`: the host unmutes the transmit channel cabled to CAL IN, the
   controller routes it to the centre probe, and all four channels receive a known pseudo-noise
   code. Correlating against the replica gives the whole-sample offset between the two radio
   boards and the phase of each channel. Needed after every retune or stream restart.
3. `calibrate_tx(power_meter)`: four-channel transmit only. Each channel is turned through three
   phases against channel 0 while a power meter straight ahead reads the sum. Needed after
   every retune. On the bench the power meter is the HackRF with a small antenna.

## What the simulator models, and what it does not

Modelled: embedded element pattern (cos^1.2, 14.6 dBi for the array), switch, line and
combiner losses, delay-line phase that scales with frequency, LNA and radio noise with
thermal noise topped up after every passive section, part-to-part amplitude and phase
errors, different errors on the LNA and bypass paths, an arbitrary oscillator phase per
radio board per retune (receive and transmit), a whole-sample offset between the boards'
streams, USB latency, switching transients, ADC quantisation, probe coupling, leakage
from the auxiliary patch and point reflectors.

Optional, `BoardModel(antenna=...)`: each element's own embedded pattern and the coupling
between patches, from `antenna.PatchArray` (each patch a transmission-line model, every
radiating edge coupled to every other through space; checked against Balanis and against
energy conservation) or a chosen coupling level (`antenna.CoupledIdeal`). Without it every
element has the same ideal pattern and nothing couples.

Not modelled: oscillator phase noise, coupling through surface waves (it adds E-plane coupling
on FR-4), the finite ground plane, multipath, radio non-linearity, a sub-sample offset between
the radio boards, temperature drift. Part-to-part errors are assumed (0.5 dB, 10 degrees per
element), not measured.

## Redesign requirements (decided 2026-11-15)

The next board revision trades the 5.8 GHz-only design for flexibility:

- Receive 0.4-6 GHz.
- All four jobs across the range: direction finding, receive beamforming, transmit beam
  steering (in bands where transmitting is allowed) and radar experiments.
- Antennas: a built-in wideband panel for the upper bands, plus connectors for external
  antennas placed per band (KrakenSDR-style) below it.
- Receivers: four coherent channels from two AD9361-family boards, and still a mode with one
  HackRF and nothing else.

The architecture for this is being worked out from `experiments/antenna_study.md` and
`experiments/wideband_study.md`; the current schematics stay as they are until it is chosen.

## Open issues

- Transmit alignment needs a far-field power meter after each retune. Tracking the transmit
  phases through coupling between quadrants would remove that, and needs patch coupling in
  the simulator first.
- Channel alignment assumes the probe code is stronger than other signals in the room during
  the alignment capture, and that the cal transmit channel is muted the rest of the time.
- Stream offsets between the radio boards are assumed constant between restarts. A dropped
  buffer breaks alignment silently; `align()` fixes it but nothing detects it yet.
- `dsp.cal_slot` only finds offsets of 0 to 63 samples with board B behind board A, and turns
  anything else into 0 without saying so. Two radio boards whose streams are started separately
  can be thousands of samples apart, either way round. The hardware backend needs a coarse step
  first (each board finds the marker on its own channels; that gives the offset modulo one
  sequence period) before the code correlation refines it.
- Power: the plan was that the host keeps the radio boards unpowered when the USB-C source
  offers less than 1.5 A, but the radio sockets' VBUS is wired straight to +5 V. The controller
  can only report the current (STATUS `usb_ma`). A load switch on the two radio sockets'
  VBUS, enabled from a spare RP2040 pin, would make the plan real.
- The control lines float, not low, while SR_OE_N is high: from power-up until the firmware
  runs, and indefinitely in BOOTSEL mode or with blank flash. "All zeros at power-up" needs
  pull-downs on the 56 register outputs (resistor arrays) or a defined state in the parts themselves.
- Schematic corrections from the datasheet check (2026-11-15), for the next regeneration:
  - LMK1C1104 (TI SNAS791D, table 6-1): pin 5 is Y2, pin 7 Y3, pin 8 Y1 (the script has 5 Y1,
    7 Y2, 8 Y3). Harmless electrically, all outputs carry the same clock, but fix the labels.
    CLKIN needs LVCMOS levels (VIH 0.7 VDD, VIL 0.3 VDD, slew of at least 0.1 V/ns).
  - TCXO: no 0.5 ppm 3.2 x 2.5 mm TCXO found with pad 1 = enable. Pad 1 is NC/VC (Epson
    TG-5006CE) or GND (ECS-TXO-32CSMV, Kyocera KT3225K), so tying it to VDD is wrong and on
    the ECS part shorts VCLK. Those parts give clipped sine (about 0.8 Vpp), which cannot drive
    CLKIN: add AC coupling, bias and a squaring buffer, or use a CMOS-output TCXO (the ones
    found are 2.5 ppm). Route pad 1 through a 0-ohm link to GND.
  - LNA enable: SKY65404-31 V_ENABLE draws 2.3 mA (3.0 mA max). Through the 1k + 100 ohm
    filter it loses about 2.5 V and the LNAs never turn on. Feed V_ENABLE through only the
    datasheet's 100 ohm and keep the 1k / 100 pF filter for the switch CTRL pins (nA).
    LNA2_EN drives five LNAs (12-15 mA) from one 74HC595 output (6 mA rated): buffer it.
  - Floating control lines: neither the BGS12WN6 nor the SKY65404-31 defines a floating state,
    and the 74HC595 leaks up to 5 uA when disabled. Add 10k pull-downs on the 56 outputs.
  - FE1.1s (Rev B12.4): BUSJ low means bus-powered (100 mA per port presumed); pull it high
    (10k to HUB_3V3) for self-powered. Pins 12 and 28 (VD18_O / VD18 in the KiCad symbol) are
    NC on Rev B: remove HUB_1V8 and its capacitors. The crystal load capacitors are on chip:
    fit no external 15 pF, use a crystal with CL 16-20 pF and ESR <= 50 ohm. REXT 2.7k 1 %.
    PWRJ (open-drain port-power enable) and OVCJ could drive a load switch on the radio VBUS.
  - SKY65404-31 footprint has a centre pad: include it, to GND. Supply: 11 mA typ (15 max) at
    3 V plus the enable current, so the 21 LNAs draw about 280 mA (380 mA max) from +3V3_RF.
  - Confirmed as drawn: BGS12WN6 pinout and truth table, SKY65404-31 pinout, LMK1C1104 1G
    (active high), FE1.1s XRSTJ, VBUSM (5 V tolerant), OVCJ.
  - Still open: the LNA input match (0.5 pF shunt, 0.6 nH) and whether the DC blocks are needed.
