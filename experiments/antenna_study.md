# Antenna study: printed patches, or something else?

Scope: the 4 x 4 panel at 5.8 GHz, 26 mm pitch. The redesign decided on 2026-11-15 (0.4-6 GHz,
see `ARCHITECTURE.md`) needs a wideband panel instead; `wideband_study.md` covers that. What
follows still decides any 5.8 GHz panel, and the printed-versus-not-printed trade-offs carry over.

`python experiments/antenna_study.py` (about 50 s) writes `antenna_results.json` and the figures
in `antenna/`.

## Short answer

1. **Patches in the main FR-4 board work only with L3 as their ground.** On JLCPCB's default
   4-layer stack the patch sits 0.21 mm above L2: 14 % efficient, and the 15.7 x 11.9 mm patch
   resonates near 6.0 GHz, off the band. With L2 cleared under the patches (ground on L3, 1.28 mm)
   it is 52 % efficient with 180-220 MHz of bandwidth, and the array reaches 14.6 dBi, which is what
   the simulator always assumed. The catch is the laminate: FR-4's permittivity at 5.8 GHz is
   uncertain by about +-0.2, which moves the patch +-128 MHz against a margin of 15-35 MHz. Order a
   test-coupon panel first, measure, then trim the patch length (about -22 MHz per +0.05 mm).
2. **A separate low-loss antenna board buys about 1.5-2 dB**: 1.5 mm RO4003C/RO4350B (76 %) or PTFE
   (84 %). Their permittivity is held to +-0.05, so a simulated array response stays usable on the
   built board. JLCPCB makes these as 2-layer boards only; join them to the main board with pins.
3. **Suspended patches** (printed on a thin FR-4 sheet about 2 mm above the main board's ground, or
   cut from sheet metal) are the most efficient (about 95 %; a measured hobby build reached 74 %) and
   the most tolerant of laminate scatter. But at 26 mm pitch the roughly 20 mm patches nearly touch:
   E-plane coupling is about -10 dB, the active reflection reaches 0.67 when scanning, and
   direction finding with a boresight-only calibration errs by 5.7 degrees (95th percentile). They
   need a measured array response, and smaller elements would help.
4. **Not worth it here**: catalogue 5 GHz ceramic patches (their -10 dB band ends below 5.8 GHz),
   FPV patches, helices, horns, chip antennas. At this pitch an element's gain inside the array is
   capped at about 5 dBi whatever it is; only efficiency raises the array's gain (limit 17.1 dBi).
5. **For direction finding, calibration matters more than the antenna.** With a boresight-only
   calibration the 95th-percentile error is 1.7-2.6 degrees for every printed option; with the
   array's real response it is about 0.5 degrees for all of them.

## How the array was simulated

1. **The simulator as it was**: identical `cos(theta)^1.2` elements, 14.6 dBi, no coupling.
2. **A coupled patch model** (`array16/antenna.py`): each patch is the transmission-line model
   (a length of wide microstrip between two radiating edges) with dielectric, conductor and
   surface-wave loss and the probe's inductance; every radiating edge couples to every other
   through space, so the 16 patches form one network. It gives the 16 x 16 S-matrix and each
   element's embedded pattern, and `BoardModel(antenna=...)` feeds them to the simulator, so the
   real calibration and direction-finding code runs on top.
3. **Full-wave simulation**: not run yet (see "Next steps").

### What the model was checked against

| Check | Result |
| --- | --- |
| Balanis example 14.1 (closed-form design) | W 1.185 cm (book 1.186), eps_eff 1.972, extension 0.081 cm, L 0.905 cm (0.906) |
| Edge conductance and mutual conductance, Balanis eqs. 14-17, 14-18a | equal to 0.07 % (the book's 120 pi ohm) |
| Far-field power against port power, lossless array | equal to 1 part in 10^6 |
| Efficiency against the Jackson-Alexopoulos cavity formulas | within 1-3 % for every printed option (FR-4 51.9 against 53.5 %) |
| Coupling against a moment method with the grounded-slab Green's function (model review) | FR-4 centre element: E -20.7 against -20.7 dB, H -20.7 against -18.8, diagonal -31.5 against -28.5 |
| Embedded gain of a lossless inner element at broadside against the cell aperture | 4.8 against 5.0 dBi |

An independent review (14 agents, each serious finding re-derived by a second one) found and
fixed four problems: surface-wave loss scaled from one edge instead of the whole patch, a coupling
sweep using the S-matrix phase instead of the receive-coupling phase (they differ by 64 degrees on
FR-4 because of the probe), a "known array response" case that was really an exact-model bound,
and scans reported on one side of an asymmetric array. It also showed where the model does not
hold:

- **Printed options**: the E-plane coupling is about right on FR-4 (surface waves take over from the
  grazing space wave); H-plane coupling is 1.5-2 dB too weak and diagonal coupling 2-3 dB too weak;
  on RO4003C and PTFE the E-plane coupling is overstated by 2-5 dB. Bandwidth is about 20 % narrower
  than the cavity formulas give, so both are listed.
- **Air patches at 26 mm**: the facing edges are 3.7-5.9 mm apart and the model's coupling is 6-7 dB
  too weak. `PatchArray.coupling_reliable` flags this, and the air figures below are the review's
  moment-method results (rooftop MPIE with the exact image Green's function, converged to 0.05 dB).

## Element options

| Option | Bandwidth, -10 dB (cavity VSWR 2) | Covers ISM | Efficiency | Element gain | Array gain | Resonance per +0.1 eps_r |
| --- | --- | --- | --- | --- | --- | --- |
| FR-4 main board, ground L2 (0.21 mm) | 114 (133) MHz | no | 14 % | -2.5 dBi | 9.4 dBi | -66 MHz |
| FR-4 main board, ground L3 (1.28 mm) | 182 (224) MHz | yes | 52 % | 3.1 dBi | 14.6 dBi | -64 MHz |
| FR-4 antenna board 1.55 mm | 204 (263) MHz | yes | 54 % | 3.3 dBi | 14.8 dBi | -63 MHz |
| RO4003C board 0.813 mm | 96 (113) MHz | no | 75 % | 5.0 dBi | 16.1 dBi | -79 MHz |
| RO4003C board 1.524 mm | 166 (214) MHz | yes | 76 % | 5.0 dBi | 16.1 dBi | -77 MHz |
| PTFE (F4B) board 1.5 mm | 176 (221) MHz | yes | 84 % | 5.9 dBi | 16.6 dBi | -103 MHz |
| Sheet metal 3 / 5 mm over ground | 356 / 652 MHz | yes | ~100 % | about 9 dBi alone | about 17 dBi | -51 / -61 MHz per +10 % height |

The 0.21 mm patch can only be matched with the probe 0.25 mm from its edge, which cannot be built;
it would need an edge feed. Probe reactance is about +31 ohm through 1.28 mm FR-4 and is absorbed
by the feed position; through 3-5 mm of air it needs a series capacitor (a capacitive disc or
L-probe feed).

## The array with coupled patches

| | FR-4 (L3) | RO4003C 1.524 | PTFE 1.5 | Air 5 mm (moment method) |
| --- | --- | --- | --- | --- |
| Coupling, centre element, E / H / diagonal (S) | -20.7 / -20.7 / -31.5 dB | -17.6 / -17.1 / -28.2 dB | -17.3 / -16.5 / -27.2 dB | -9.9 / -18.2 / -21.7 dB |
| Embedded gain, centre / edge / corner | 2.5 / 2.3-2.9 / 2.7 dBi | 3.9 / 3.6-4.5 / 4.3 dBi | 4.2 / 4.1-5.0 / 5.0 dBi | |
| Scan loss at 30 / 45 / 60 deg, E-plane (worse side) | 0.7 / 1.9 / 3.4 dB | 0.8 / 2.2 / 4.1 dB | 0.8 / 2.3 / 4.4 dB | 0.1 / 1.6 / 6.6 dB |
| Same, H-plane | 1.1 / 2.6 / 5.6 dB | 0.9 / 2.1 / 4.9 dB | 0.9 / 2.0 / 4.8 dB | |
| Worst active reflection, scans to 45 deg | 0.27 | 0.35 | 0.35 | 0.67 |

(`embedded_patterns.png`, `scan_gain.png`, `calibration_residue.png`)

**2-bit phase states.** Over scans of 0-45 degrees the 2-bit states cost 0.4-0.55 dB of gain on
average (up to about 1 dB). Sidelobes are the real limit: -12.5 to -14 dB at broadside, but
quantisation lobes reach -7.6 dB at 20 degrees and -2.8 to -3.5 dB at 40 degrees for one analog
beam (HackRF mode, or the four-radio receive beam as `Array.beam(word)` forms it). Phasing each
quadrant channel toward the beam (four-radio transmit, monopulse, and now
`Array.beam(word, toward=...)`) halves the extra loss and brings the 40-degree lobe to about -6 dB.
Clean beams beyond about 30 degrees need a third phase bit. Direction finding does not use analog
beams and is unaffected.

## Direction finding: calibration matters more than the antenna

Four-radio mode, a single tone at -75 dBm, 72 directions out to 55 degrees, two boards with the
simulator's part-to-part errors, calibrated on a boresight source as `calibrate()` does
(`df_error.png`). Median / 95th percentile, degrees:

| Antenna | Boresight calibration only | Located with the exact array response (a bound) | Board's eps_r off; located with the simulated nominal response |
| --- | --- | --- | --- |
| Ideal elements (old simulator) | 0.12 / 0.45 | | |
| FR-4 (L3) | 0.37 / 1.66 | 0.12 / 0.47 | eps_r +0.2: 0.31 / 1.15 (boresight only 0.35 / 1.98) |
| RO4003C 1.524 | 0.51 / 2.34 | 0.10 / 0.47 | eps_r +0.05: 0.17 / 0.74 |
| PTFE 1.5 | 0.52 / 2.56 | 0.10 / 0.47 | eps_r +0.05: 0.19 / 0.83 |
| Air 5 mm (moment method) | 2.47 / 5.71 | | |

The floor of about 0.5 degrees is the simulator's part-to-part hardware errors acting through the
decode, not noise. A simulated array response helps most where the laminate is predictable, which
argues for Rogers or PTFE if you want to calibrate from simulation; a measured response (turntable
or several known source positions) works for any of them.

How much coupling a boresight calibration tolerates (`df_vs_coupling.png`; the FR-4 receive-coupling
matrix with all its neighbours scaled to the level and its phase turned through 0-270 degrees,
range over phase):

| Nearest-neighbour coupling | -30 dB | -25 dB | -20 dB | -17 dB | -15 dB | -12 dB |
| --- | --- | --- | --- | --- | --- | --- |
| Boresight calibration, 95th pct | 0.7-1.2 | 1.1-1.8 | 1.9-3.2 | 2.6-4.6 | 3.3-5.8 | 4.3-8.3 deg |
| Exact array response, 95th pct | about 0.5 deg throughout | | | | | |

## Facts that decide the printed options (verified research)

- JLCPCB's default 4-layer stack (JLC04161H-7628): L1-L2 0.21 mm 7628 prepreg, L2-L3 1.065 mm core,
  so a patch on L1 is 0.21 mm above L2 and about 1.29 mm above L3. Its other templates (-2116,
  -3313, -1080) give L1-L3 of 1.34-1.40 mm, slightly better. [1]
- JLCPCB publishes Dk only (7628 prepreg 4.4, core 4.6, stated in different places), with no
  frequency and no loss tangent. With L2 cleared the stack is about 4.56 in series, which puts
  the 11.9 mm patch near 5.70 GHz. FR-4 permittivity also falls with frequency by a
  laminate-dependent 0.05-0.25. [2, 3] One hobby-fab FR-4 board measured 4.80 and loss tangent
  0.033 at 5.1 GHz. [4]
- Solder mask over the patches lowers them by 24-36 MHz; model it or open the mask. On the L4
  feed lines it lowers Z0 to about 48 ohm and raises eps_eff 7 %, so draw any delay-line bits with
  the mask in the model. [5]
- FR-4 microstrip on L4 loses 0.17-0.29 dB/cm at 5.8 GHz; keep via-to-LNA lines short.
- JLCPCB's Rogers and PTFE boards are 2-layer only (RO4350B, ZYF PTFE). PCBWay offers mixed
  Rogers + FR-4 multilayers by quote. RO4350B's design Dk is 3.66, not the 3.48 process value
  JLCPCB lists. [6, 7]
- A hobby-built suspended FR-4/air 5.8 GHz patch (two 1 mm boards 4 mm apart, pin-fed) measured
  74 % radiation efficiency; FR-4 batch-to-batch shifts of up to 100 MHz were "normal". The
  carrier sheet's +-0.1 mm thickness tolerance moves such a patch +-60-87 MHz, more than the
  laminate's permittivity does. [8]
- Catalogue patches (Taoglas WLP.4958.12.4.A.02, Abracon APAKN1304-S5517-T) fit the grid but are
  tuned for 5 GHz Wi-Fi: their -10 dB band ends at about 5.73 and 5.78 GHz. [9]
- Measured H-plane coupling of two FR-4 5.8 GHz patches (1.6 mm, 7.9 mm gap): -14.4 dB, against
  -16.5 dB in HFSS. Array16's gap is 10.3 mm. [10]
- U.FL is rated only to 6 GHz with VSWR up to 1.5 above 3 GHz; pins are the cheapest and shortest
  link between an antenna board and the main board. [11]

## Calibration practice (verified research)

- **Source distance.** A boresight source at 2 m is beyond 2D^2/lambda, yet its wavefront
  curvature leaves 0.6 / 2.9 / 5.3 degrees on the inner, edge and corner elements (2.1 at 5 m). [12]
  `calibrate(range_m=...)` now removes it. In the simulator it does not bias single-source
  direction finding at all (the error is symmetric about the centre, and only odd errors move a
  beam's peak), but it defocuses beams and moves nulls, so correct it anyway.
- **Reflections during calibration end up in the reference**: one at -20 dB perturbs elements by
  up to 0.9 dB / 5.7 degrees, at -10 dB by 3.3 dB / 18 degrees. Calibrate outdoors, with absorber,
  or from two source positions. [13]
- **Angle-dependent (manifold) calibration beats a single coupling matrix** in measured
  comparisons, but needs dense angular sampling; a coupling-matrix fit holds up with fewer
  positions. [14]
- **The 16 x 16 Sylvester Hadamard code leaks** errors in each element's 180-degree state into the
  element on its all-ones column; a 20-state Paley code has no such column. [15]
- **The centre probe** is 18, 41 and 55 mm from inner, edge and corner elements: good for tracking
  channel drift, not for calibrating elements.
- The LNA's input match differs between its on and off states (|Gamma| about 0.32 and 0.89), so
  receive and transmit embedded patterns differ; a full-wave model should embed the real loads.

## Next steps

1. **Full-wave check** with openEMS (FDTD; graded mesh, lumped ports, near-to-far-field). There is
   no working Homebrew formula: build it from source, and against Python 3.14, because h5py has no
   wheels for the default Python 3.15 [16, 17]:

       brew install boost hdf5 cmake
       uv venv --python 3.14 --seed ~/opt/openems-venv && source ~/opt/openems-venv/bin/activate
       git clone --recursive https://github.com/thliebig/openEMS-Project.git
       cd openEMS-Project && ./update_openEMS.sh ~/opt/openEMS --python

   Estimated: minutes for one patch; 4-100 minutes per excited port for the 4 x 4, so 16 runs (8 by
   symmetry) take hours, in 1-4 GB of RAM. Cross-check single patches, pairs and a 2 x 2 with
   EMerge (FEM, `pip install emerge`, Python 3.13 or older). [18]
2. **A coupon panel** before the array: single patches at 2-3 lengths, a 2 x 2, a ring resonator.
3. **Calibration**: add the source-distance correction, use a 20-state code, and plan a manifold
   measurement (several known source positions) rather than relying on boresight alone.

## Sources

1. JLCPCB stack-ups (auto-generated from JLCPCB's data): https://github.com/gsuberland/jlcpcb_autogenerated_stackups
2. JLCPCB capabilities: https://jlcpcb.com/capabilities/pcb-capabilities
3. Djordjevic et al., wideband FR-4 model: https://mtt.etf.bg.ac.rs/Mikrotalasna.Tehnika/Clanci/FR4_EMC_2001.pdf
4. FR-4 measured near 5 GHz: https://arxiv.org/pdf/2310.01171
5. JLCPCB laminate structures and solder mask: https://jlcpcb.com/help/article/multi-layer-pcb-standard-laminated-structures
6. JLCPCB high-frequency boards: https://jlcpcb.com/pcb-fabrication/high-frequency-pcb
7. PCBWay advanced capabilities: https://www.pcbway.com/advanced-pcb-capabilities.html
8. Baert, triple-feed patch: https://www.maartenbaert.be/quadcopters/antennas/triple-feed-patch-antenna
9. Taoglas WLP.4958.12.4.A.02: https://www.taoglas.com/datasheets/WLP.4958.12.4.A.02.pdf
10. FR-4 5.8 GHz patch coupling: https://www.radioeng.cz/fulltexts/2019/19_02_0399_0406.pdf
11. Hirose U.FL: https://docs.rs-online.com/3355/0900766b8105243f.pdf
12. Near-field calibration distance: https://arxiv.org/pdf/2209.03082
13. Array calibration with multipath: https://cris.huji.ac.il/en/publications/array-calibration-in-the-presence-of-multipath-2/
14. EUSIPCO 2019 calibration comparison: https://www.eurasip.org/Proceedings/Eusipco/eusipco2019/Proceedings/papers/1570532588.pdf
15. Hadamard demultiplexing errors: https://researchr.org/publication/Silverstein97-0
16. openEMS requirements: https://docs.openems.de/en/latest/install/requirements.html
17. openEMS CI (macOS ARM, Python 3.14): https://github.com/thliebig/openEMS-Project/blob/master/.github/workflows/ci.yml
18. EMerge: https://pypi.org/project/emerge/
