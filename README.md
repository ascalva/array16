# array16

Schematics, host software and a simulator for the 4x4 phased-array SDR board at 5.8 GHz
(16 patches in four quadrants, 2-bit phase shifters, four radio channels or a HackRF alone).

Nothing here needs hardware. The simulator stands in for the board, the radios and the
room, and the host library talks to it through the same three calls it will use on the bench.

## Run it

    uv venv .env && uv pip install --python .env/bin/python -e ".[dev]"
    python -m pytest                          # 75 tests: virtual bench, antenna model, optimiser, firmware core, and the whole firmware in an RP2040 emulator; about 5 s
    python experiments/virtual_bringup.py     # statistics over 40 simulated boards
    python experiments/antenna_study.py       # antenna options, coupling, scanning, DF with real patches (~75 s)
    python hw/design.py                       # regenerate the KiCad project in hw/kicad/ (needs KiCad 7+)
    python hw/check.py                        # compare KiCad's netlist with the script's (needs kicad-cli 7+)

The firmware-core tests build `firmware/core/ctrl.c` with the local C compiler and skip
without one. `hw/design.py` finds KiCad's symbol libraries on Linux and macOS, or from
`KICAD_SYMBOL_DIR`. The firmware itself builds with pico-sdk, see `firmware/README.md`.

A first session:

    from array16 import Array
    from array16.sim import SimBackend, Source

    sim = SimBackend("quad")                  # or "hackrf"
    arr = Array(sim)
    sim.sources = [Source(0, 0, -75)]         # test source straight ahead
    arr.calibrate()
    sim.sources = [Source(25, 30, -75, kind="qpsk")]
    print(arr.locate(blind=True))             # -> about (25, 30) degrees

## What is here

| Path | What it is |
| --- | --- |
| `array16/board.py` | Geometry, RF budget and the 56-bit control word. The single description of the hardware |
| `array16/control.py` | Beam steering states and snapshot pattern schedules |
| `array16/protocol.py` | Frames between host and beam controller |
| `array16/fwsim.py` | Reference model of the controller firmware; the firmware's specification |
| `array16/backend.py` | What `Array` needs from a backend (simulator or hardware) |
| `array16/sim.py` | Simulator: board, scene, radios |
| `array16/antenna.py` | Patch design and the 16 patches as one network coupled through space: S-matrix and embedded element patterns for the simulator |
| `array16/dsp.py` | Alignment, snapshot decoding, direction finding |
| `array16/api.py` | `Array`: calibrate, point, snapshot, locate, monopulse, transmit alignment |
| `tests/` | Virtual acceptance tests, one per bench check in the design doc, and the firmware core against `fwsim.py` |
| `firmware/` | RP2040 beam-controller firmware: portable core plus PIO, DMA, clocks and USB |
| `experiments/` | Virtual bring-up statistics (`results.json`); antenna study (`antenna_study.md`, `antenna_results.json`, `antenna/`); wideband study (`wideband_study.md`, `wideband_results.json`, `wideband/`) |
| `hw/design.py` | The circuit, sheet by sheet, as a script |
| `hw/kicad/` | Generated KiCad 7 project (opens in KiCad 7, 8 and 9) |
| `hw/array16_schematics.pdf` | The 26 schematic pages |
| `ARCHITECTURE.md` | How the pieces fit, the protocol, the calibration procedures, known gaps |

## State

- **Redesign decided (2026-11-15)**: 0.4-6 GHz, all four jobs, built-in wideband panel plus
  external antennas, four-radio and HackRF modes (`ARCHITECTURE.md`, "Redesign requirements").
  The evidence is in `experiments/antenna_study.md` and `experiments/wideband_study.md`; the
  architecture is being chosen. The schematics below describe the 5.8 GHz board.
- Schematics (5.8 GHz board): first pass, 614 parts, 491 nets, connectivity by net labels. KiCad's
  netlist matches the script's. The three pinouts marked `verify` have been checked against
  datasheets: the corrections, including two real faults (TCXO pad 1, LNA enable drive), are under
  "Open issues" in `ARCHITECTURE.md`. No footprints for the custom parts and no PCB layout yet.
- Host library and simulator: working, 37 tests pass. The simulator can now use real coupled
  patches (`array16/antenna.py`) and runs at any frequency; snapshot decoding and beam steering use
  the delay lines' actual phase at the tuned frequency.
- Firmware: compiles without warnings (Arm GNU 14.3, pico-sdk 2.3.1) and runs in the rp2040js emulator, where 20
  scenarios match `fwsim.py` exactly (`firmware/emu/`, `tests/test_firmware_emulated.py`, skipped without Node and
  rp2040js). The emulator found and the firmware now fixes two faults: a missing reference could hang it, and a lost
  reference froze sequences while STATUS still reported it present. Not yet run on a real RP2040.
- Not written yet: the hardware backend (bladeRF/AD9361 radios, the HackRF tools for the HackRF).
