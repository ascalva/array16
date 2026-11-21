<!-- Proposal for owner approval. Produced 2026-11-15 by a design panel (four independent architectures,
three judges, one synthesis) from the requirements in ARCHITECTURE.md and the evidence in
experiments/antenna_study.md and experiments/wideband_study.md. Nothing in hw/ has been changed for it yet.

STATUS 2026-11-15: DRAFT. A final adversarial check (numbers, parts, internal consistency) was running
when this was written; its corrections are not applied yet. Checked so far: the 4-bit delay set gives a
toggle within 0.028 dB of 180 deg everywhere in 2.3-6 GHz (using state pairs, not single bits; single
bits alone are up to 1.04 dB off), and 352 ps covers the 306 ps that +/-60 deg at any azimuth needs. -->

**Status: draft, pending a final check and the owner's approval.**

# array16 rev B: recommended architecture

**Decision, 2026-11-20: the radio is the HackRF Pro alone.** The owner's decision. Its direct consequences for this proposal, read from the sections below:

- The four-radio mode leaves the plan, and with it the two bladeRF 2.0 micro xA4 ($1,080; sections 3 and 10). Boards and parts stay at about $1,300, and the working, calibrated system at about $1,550.
- Four simultaneous coherent channels become one: the board's MODE star sends the four filtered, polarity-switched channels to the clamped HackRF port, and the HackRF takes CLKIN and trigger from the board (section 5, HackRF-only). The noise figures are section 5's HackRF-mode row, 1.1 dB at 0.433 GHz to 3.7 dB at 5.8 GHz, and the panel draws the HackRF mode's 9.0 W (section 10).
- Direction finding goes by snapshots through the board's states, or by time modulation, and needs steady signals (section 9: HackRF DF needs steady carriers). Simulated on the panel it gives 0.55/1.73°, 0.45/1.13° and 0.29/0.63° at 2.4, 3.5 and 5.8 GHz, against 0.12/0.49° to 0.10/0.46° with four channels. Receive beams through the tile's delays stay (section 9, "HackRF too"); beams and nulls formed across four channels in real time go with them, and one channel makes them in software: any receive weight applies after one capture through 16 switched states, at 1/16 of the static beam's signal-to-noise ratio, down to the state calibration's floor, and for signals narrower than the comb spacing f_p when decoded by comb lines (5.88 kHz on rev A, 3.92 kHz on rev B with a marker), with nulls of strong signals made before the converter by the board's states (`docs/VIRTUAL_ANTENNA.md` sections 6.1 and 6.4).
- Transmit through the panel is the HackRF figure of section 7: +1.5 dBm EIRP at 2.44 GHz and -17.8 dBm at 5.8 GHz at full drive (computed for a HackRF One's output; V1's sessions T measure the Pro's), against +8.1 and +6.2 dBm with the bladeRFs. A periodic sequence of the tile's delay states gives continuous weights on comb line 1, at most 0.98 to 0.99 of the amplitude from 2.3 to 6 GHz, for signals narrower than f_p. The continuous-phase steering of four whips on X0-X3 (+15.9 dBm EIRP at 433 MHz) was four transmit channels and leaves with them; through one channel a periodic sequence of the four channels' polarity states gives each whip a continuous weight on comb line 1, at most 2/π = 0.64 in amplitude with a mirror image on line -1, for signals narrower than f_p (`docs/VIRTUAL_ANTENNA.md` sections 6.2 and 6.3).
- Radar is section 7's HackRF-only form: the LMX2572 sends CW, FSK or ramps through RADAR OUT (+10 dBm EIRP) while the HackRF listens on the panel beam, CW Doppler of 1 m² to 300-700 m in 100 ms. The FDD radar on TX3, and passive radar's one reference with three surveillance channels, took four channels.
- Calibration drops the alignment between radios (section 6: FDD, register writes, 3-5 synthesizer offsets for the inter-board delay, libbladeRF overrun flags). With one LO the HackRF needs no per-retune alignment, the LMX2572 on the shared 10 MHz builds static tables with Paley-20 codes, and transmit is calibrated by REV and the AD8317.
- The parts and code that exist for the four-radio mode are to be re-examined: the bladeRF clock and trigger lines (bladeRF A REFIN, A CLKOUT to B CLKIN, J51-1), the RX and TX paths from each channel to the radios and their USB 3 ports, the preselector's reason in section 3 (the AD9361's 3×LO response), the 4-whip simultaneous mode with its ring per band (section 4), hw_bladerf.py (section 8), and step 0 of the build order (section 12). Section 13's decisions are revisited in that light, 3 (transmit ports) and 9 (radar hardware) first, since their recommendations rest on four transmit channels and on passive radar's four receive channels.
- Rev B will be revised around one coherent channel. The sections below are the proposal of 2026-11-15, drafted for four.

## 1. Recommendation

Build Proposal 3's split as **two 4-layer boards** joined by five SMA jumpers. The **core** holds 16 external SMA ports (PMA2-63LN+ LNAs, Latin-square selector), four channels (preselector that also filters transmit, bypassable amplifier, exact 0/180 switch), a clamped HackRF port, an LMX2572 PN calibration source and the RP2040; alone it is a 16-port 0.4-6 GHz KrakenSDR-class instrument. The **tile** is a 6×6 connected Vivaldi grid at 25 mm (inner 4×4 active, 2.3-6 GHz, unambiguous to 6.0 GHz), each element with an LNA, a 4-bit true-time delay and a switched calibration arm, so conducted calibration never radiates. Two bladeRF 2.0 micro xA4 give four coherent channels; a HackRF One works alone. About $1,300 plus $1,080 of radios, 7.3 W from USB-C, 520 hours.

## 2. Block diagram

```
TILE  36 Vivaldi flares, 25 mm, 40 mm tall, cards in columns (V-pol); inner 16 active; 18 dummies -> 49.9R; 2 = probes (49.9R when idle)
 cell (x16): flare-[SW_CAL]-[TR]-PMA3-83LNW+-[TR]-[TTD 38.5/52.5/117.5/139 ps]-+ 4 cells -> Wilkinson
                 |     '---- TX bypass ----'                                  -> [QAMP|bypass] -> CH0..3
           16-way cal tree <-[tree|probe L|probe R]<- tile CAL               == 5 jumpers + ribbon
CORE  X0..X15: SMA-TVS-(5 V bias tee)-DC block-[SW_CAL]-HPF 330 MHz-(X0-3: TR)-PMA2-63LN+-(TR)-> group SP5T (4 + TERM)
                        '- 16-way cal tree
 channel k (x4): [SEL EXT|TILE]-[B1..B5 = TX filter]-[PS]-+-[AMP|bypass]-[POL 0/180]-[MODE]-+-> RX k (bladeRF)
                                                          '<- TX k (bladeRF)                 '-> HackRF star
 HackRF: star -[driver-PE43711 (boots 31.75 dB)-3 dB-clamp-6 dB-clamp-10 dB | TX thru]- SMA <-> HackRF One
 radar:  TX3 or synth -[915|2440|5800 MHz bank]- RADAR OUT -> cavity BPF -> horn >= 1.5 m aside
 cal:    LMX2572 -[BPSK]-20 dB-PE43711-[inject|monitor]-[TERM|ext tree|tile|CAL-ANT +30 dB]; monitor -> AD8317 or RX0
 clock:  10 MHz TCXO -> LMK1C1106 -> RP2040 | bladeRF A REFIN | HackRF CLKIN | LMX2572; A CLKOUT -> B CLKIN
         RP2040 TRIG -> J51-1 (A, B), HackRF trigger
 ctrl:   RP2040 PIO -> lane 0 (core, 56 b) || lane 1 (tile, 104 b); power: USB-C 5 V/3 A, load switch per amplifier
 radios: A RX1/2 = CH0/1, TX1/2 = TX0/1; B RX1/2 = CH2/3, TX1/2 = TX2/3; each on its own Mac USB 3 port
```

## 3. Element chain and RF path

| Stage | Part | Key spec | $ | Qty | Why |
| --- | --- | --- | --- | --- | --- |
| Panel element | Connected antipodal Vivaldi, 0.8 mm FR-4 cards | 2.3-6 GHz | 0.7 | 36 | low Q, no balun |
| Cal arm, T/R | 3× Infineon BGS12WN6 | IL 0.15-0.62 dB; RF1-RF2 isolation 63→31 dB typ (23 dB min at 6 GHz) | 0.32 | 72 | antenna off during calibration; TX bypass |
| Panel LNA | Mini-Circuits PMA3-83LNW+, 0402CS choke | 21 dB, NF 1.2-1.6 dB, 58 mA | 9.0 | 16 | flat gain before delay loss |
| True-time delay | 38.5/52.5/117.5/139 ps, 8× BGS12WN6 | a state pair within 9° of 180° at every frequency (≤0.03 dB phase-only, about 0.2 dB with the line-loss imbalance; the measured tables decode both); RMS error ≤11.4°; idle-arm second resonances ≥6.6 GHz (≥6.09 GHz with FR-4 spread and pad loading) | 2.6 | 16 | DF toggles, steering |
| Combiner | 2-section Wilkinson, 2 levels | 0.5-0.8 dB excess | 0.2 | 4 | coherent sum, TX split |
| QAMP | PMA3-83LNW+ with bypass | +21 dB | 10.5 | 4 | 5.8 GHz NF 3.4, not 5.0 dB |
| Ext input | edge SMA, onsemi ESD8472MUT5G (LCSC C133346), bias tee, 100-150 pF C0G DC block before SW_CAL, 330 MHz HPF after the arm | ≤0.3 dB | 2.3 | 16 | ESD, FM rejection; the DC block keeps the 5 V antenna bias off the switches and the cal tree |
| Ext LNA | Mini-Circuits PMA2-63LN+ | NF 0.5-1.0 dB, 24→11.5 dB, 44 mA | 6.5 | 16 | lowest NF below 2.3 GHz |
| Group selector | 4× BGS12WN6 | 4 ports + TERM | 1.3 | 4 | rows, columns, diagonals in one slot |
| Preselector = TX filter | B1 LC LPF 0.62 GHz; B2 LC HPF + LFCN-1200+; B3 HFCN-1200+ + LFCN-2750+; B4 HFCN-2700+ + LFCN-6000+; B5 BFCN-5750+ + LFCN-6000+ | each <3:1; 3.3 dB at 2.3 GHz, 5.6/6.2 dB at 5.8 GHz | 24 | 4 | AD9361 3×LO response, harmonics |
| Channel amp | PMA3-83LNW+ with bypass | +21 dB; +9 dBm max input below 0.5 GHz (+16 dBm above) | 10.5 | 4 | robust state, HackRF gain |
| Channel polarity | TCM1-83X+ with **one** BGS12WN6, pads 5/6 tied (footprint also takes TCM1-63AX+) | 180°, ≤13° typ unbalance at 2-4 GHz, calibrated | 7.3 | 4 | Hadamard across channels |
| HackRF port | driver, PE43711B-Z (P/S strapped to serial, DI to GND), 2× BAT15-04W, pads | ≤-10 dBm at SMA | 19 | 1 | -5 dBm rating held in hardware |
| Cal source | TI LMX2572, PE43711B-Z, AD8317ACPZ monitor (the AD8318 has no distributor stock) | 12.5 MHz-6.4 GHz; -96 to -45 dBm per cell | 34 | 1 | coherent PN in both modes |
| Radar bank | LFCN-1000+, LFCN-2750+, BFCN-5750+ + LFCN-6000+ | 915/2440/5800 MHz | 16 | 1 | TX3 or synthesizer to RADAR OUT |
| Radios | Nuand bladeRF 2.0 micro xA4 | AD9361, 47 MHz-6 GHz, 2×2 | 540 | 2 | cheapest 6 GHz AD9361 |

About 300 BGS12WN6 in all, machine-placed.

## 4. Antennas

**Panel (2.3-6 GHz).** Six 150×45 mm cards stand in columns (vertical planes, so the E-field across each flare gap is vertical) on the tile ground: 36 flares, neighbours joined, 40 mm deep (the 2.4-6 GHz thickness bound is 4.2 mm; depth buys low-end match). Pitch 25 mm = λ/2 at 6.0 GHz. The inner 4×4 is active; the 20-flare outer ring ends in 49.9 Ω (two flares are probes, switched to 49.9 Ω when not in use), as the measured TCDA and BAVA arrays needed. Each card's tab sits in a plated slot (or a right-angle joint to a GCPW pad with a via to its cell), with a continuous soldered ground seam (or ground posts every 5 mm, λ/10 at 6 GHz) and a slotted comb frame holding the 25 mm pitch to ±0.1 mm. These joints sit outside the conducted calibration plane, so the transition and the radome go into the openEMS model and the coupon cards. Directivity bound 4πA/λ² is 9.2/12.3/16.7 dBi at 2.44/3.5/5.8 GHz; with efficiency of -3 dB at 2.3 GHz rising to -1.2 dB above 4.5 GHz, realised gain is 6.4/10.7/15.5 dBi, HPBW 62/43/26°. Coupling of -10 to -15 dB and edge VSWR 3-4 below 3 GHz make a measured manifold mandatory.

**External (0.4-2.3 GHz, usable to 6 GHz).** Edge SMAs X0-X15 (p = 4·row + col), groups {0,6,11,13} {1,7,10,12} {2,4,9,15} {3,5,8,14}; X0-X3 also transmit. Layouts were searched and checked in `experiments/evolve_study.md` (azimuth DF, 3° calibration residual, 20 dB, sources up to 15° or 30° above the horizon). With 8 or 16 whips one fixed layout covers all of 0.4-2.3 GHz, so nothing is rebuilt between bands; only the 4-whip simultaneous mode keeps a ring per band, because no 4-whip layout survives a 5.75:1 band:

| Setup | Layout on the 500 mm plate | 95th percentile error, 0.4-1 / 1-2.3 GHz | Gross errors, sources 0-30° up |
| --- | --- | --- | --- |
| 16 whips, 4 slots | two rings of 8, r 225 and 112 mm, inner ring turned 22.5° | 0.6-1.4° / 0.25-0.6° | 0.0 % (0.1 % at 5° and 0 dB) |
| 8 whips, X0-X7, 2 slots | fixed irregular layout, `random8` in `experiments/evolve/results.json` (or rings of 225/112/78 mm rebuilt per band: 1.6° worst, no gross errors in any condition) | 0.9-2.1° / 0.45-0.9° | ≤0.1 % (≤0.2 % at 5° and 0 dB) |
| 4 whips, X0-X3, simultaneous | ring per band: r 133 / 74 / 42 mm for 0.40-0.72 / 0.72-1.3 / 1.3-2.3 GHz (r/λ 0.18-0.32; a 4-whip ring turns ambiguous at 0.354 λ) | | |
| λ/2 4×4, 16 ports (beams) | 655 mm at 915 MHz, 381 mm at 1.575 GHz | | |

A single 225 mm ring of 8 also covers the band for sources within 15° of the horizon (95th percentile ≤1.6°) but fails above 1.6 GHz for sources 20-30° up (up to 30 % gross errors): larger rings trade robustness for accuracy. One drilled 500 mm plate holds all of these and a centre monopole driven from CAL-ANT; for grids the host removes the spherical phase.

## 5. Radios, channels, filters, protection

- **4-radio:** each channel takes its tile quadrant or its external group (one port per slot); mixing is allowed. bladeRF A locks to the 10 MHz on REFIN, A's CLKOUT clocks B, J51-1 starts both; each radio has its own Mac USB 3 port (320 MB/s at 4×20 MSPS), so four-radio mode needs three free USB ports (the controller may sit on a hub).
- **HackRF-only:** MODE sends the four filtered, polarity-switched channels through a star to the clamped port; the HackRF takes CLKIN and trigger from the board.
- **Gain states:** panel QAMP, QAMP+AMP (5-6 GHz), robust (amplifier after the filter), low; external mid or high.

NF referred to the element or jack, dB (bladeRF assumed 5-7 dB, HackRF 16-22 dB with amp off). External: mid gain state at 0.433 and 0.915 GHz, high state above; X0-X3 add one switch before the LNA (+0.15-0.5 dB); 6.0 GHz is worse than 5.8 GHz because the switches' loss steps to 0.62 dB above 5.925 GHz. Panel figures use the LNA's 50 Ω noise figure; with the edge elements' VSWR of 3-4 below 3 GHz the 2.3-3 GHz figures become a centre-to-edge range once openEMS gives the active S11:

| | 0.433 | 0.915 | 1.575 | 2.3 | 2.44 | 3.5 | 5.8 GHz |
| --- | --- | --- | --- | --- | --- | --- | --- |
| External | 1.1 | 1.3 | 1.2 | 1.4 | 1.4 | 2.0 | 5.5 |
| Panel QAMP / robust | | | | 1.9/2.1 | 2.0/2.1 | 2.2/2.6 | 3.4/5.0 |
| HackRF mode | 1.1 | 1.2 | 1.5 | 2.0 | 2.0 | 2.2 | 3.7 |

Panel G/T (antenna 290 K): -20.2/-16.1/-12.5 dB/K at 2.44/3.5/5.8 GHz. Bands under 3:1 put 3×LO in a stopband (12-30 dB). The QAMP compresses at -20 dBm (5.8 GHz) element-referred in-beam; near Wi-Fi use the robust state. The clamps hold the HackRF at ≤-10 dBm (BAT15-04W maximum forward voltage; an 11 dB last pad gives ≤-11 dBm) with the driver saturated and the DSA at 0 dB, whatever software does. Below 0.5 GHz the channel amplifier takes only +9 dBm, so the external gain state defaults to mid there and firmware refuses AMP in B1 while an X port transmits.

## 6. Calibration

**Injection points.** 32 arms behind every flare and jack, fed by two 16-way resistive trees, cover live everything after them, from LNA to radios. Code reaching an antenna is ≤-84 dBm (RFIN-RF1 isolation 53→25 dB; ≤-82 dBm worst case at 6 GHz): nothing measurable radiates (P3's taps gave -38 to -42 dBm EIRP) and the scene stays out of the measurement. Flares, cables and antennas rely on the manifold and in-place radiation (centre monopole, two probe flares) at ≤-60 dBm EIRP, refused in 15.205 bands (48 % of 0.4-6 GHz, 51 % of the panel's 2.3-6 GHz: between 3.6 and 5.15 GHz only 4.4-4.5 GHz may radiate), where stored per-port delay/phase fits apply (+1-2°).

**Level and code.** -85 dBm per cell gives 49 dB correlation SNR (0.15°); the DSA raises it to -59 dBm when an emitter leaks in. BPSK between the LMX2572's complementary outputs, 1023 chips at 2 Mchip/s restarted at the marker, resolves ±255 µs, ending the 0-63-sample limit. Between calibrations the synthesizer is muted; bench-check stray leakage ≥40 dB below the calibration level.

**Triggers.** FDD, manual gain, quadrature tracking frozen by register writes; inject after power-up, stream restart, retune, rate or gain change and any board-state change: 3-5 synthesizer offsets in 5 ms give the inter-board delay (<0.1 ns) and phases. Watchdog: a calibration slot every 100 ms, libbladeRF overrun flags, ΔT >2 °C. HackRF-only needs no per-retune alignment (one LO): the synthesizer, on the shared 10 MHz, builds static tables with Paley-20 codes.

**Tables.** All 16 delay states per element or port, receive and transmit, per band and gain state, measured at first use of a frequency (5-20 ms), cached by frequency, temperature and serial; 10 MHz grid near the idle-arm resonances (3.3 and 3.9 GHz nominal, 3.0-3.9 GHz over FR-4 spread), 50 MHz elsewhere.

**Manifold.** openEMS prior, then an outdoor DIY turntable with the board's synthesizer driving a horn ≥5 m away, two source heights and ground absorber (a -10 dB bounce moves phases ±18°): every 100 MHz, ±60° in 5° steps (a weekend), in a MANIFOLD mode (CAL-ANT under TX_ARM, outside 15.205 bands, ≤-41 dBm EIRP; at the general -60 dBm cap about 100 codes, 50 ms, per point also work). In the restricted slices (3.6-5.15 GHz except 4.4-4.5 GHz, and the smaller ones) the manifold is the openEMS model anchored at the nearest measured bands, and DF there carries the +1-2° of unmeasured bands. The EU profile has no general allowance outside SRD bands, so there it measures only in SRD bands. Target ≈0.5° instead of 3-8° with boresight only.

**Transmit.** Arms on the tree: TX0-3 send orthogonal codes through the tree into RX0 (phases in 5-10 ms) at ≥40 dB below full drive (≥45 dB at 5.8 GHz), so switch leakage radiates ≤-60 dBm EIRP; the AD9361's transmit phase against TX attenuation is measured once so the backed-off result transfers to full drive; HackRF transmit uses REV (the trees driven in reverse into the monitor) and the AD8317; the antenna term comes from the manifold by reciprocity.

## 7. Transmit and radar

- Transmit passes the selected preselector band; RADAR OUT has its own bank, plus a cavity band-pass (≥50 dB at 11.6 GHz; LFCN-6000+ alone gives 31-40 dB) at 15.245 or amateur power.
- EIRP at full drive (bladeRF +6/+5/+3 dBm at ≤1.3/2.44/5.8 GHz): panel +8.1/+6.2 dBm (2.44/5.8); four whips +15.9 (433 MHz), +14.6 (915), +12.9 dBm (1270); HackRF via panel +1.5/-17.8 dBm. Compliance cannot rest on the assumed drive (Nuand quotes up to +8 dBm): hw_bladerf.py and hw_hackrf.py cap TX gain per band and region from a measured EIRP per path. US 15.249 steering needs 9.3 dB back-off at 2.44 GHz and 7.5 dB at 5.8 GHz.

| Band | US licence-free | EU licence-free | Amateur |
| --- | --- | --- | --- |
| 420-450 MHz | none useful | 433.05-434.04: 10 mW ERP at ≤10 % duty; 434.04-434.79: 10 mW, ≤25 kHz channels (EN 300 220; back off 3.8 dB) | 70 cm |
| 863-928 MHz | 902-928: 15.249 -1.25 dBm average; 15.245 18.75 dBm | 863-870: 25 mW ERP, duty cycle | 33 cm, pulse |
| 1-2.3 GHz | none | none | 23 cm (no pulse), 2300-2310 |
| 2.4 GHz | 15.249; 15.245 at 2435-2465 | 10 mW (25 mW radiodetermination) | 2390-2450 |
| 5.8 GHz | 15.249; 15.245 at 5785-5815 | 25 mW, time-out for carriers | 5650-5925, pulse |

- Harmonics: 11.6 GHz is restricted (≈60 dB needed; measure once with a borrowed analyser); LFCN-2750+ and LFCN-1200+ remove 4.88 and 2.745 GHz. Firmware TX_ARM(band, region, timeout) gates every transmit bit; a licence flag adds CW ID. Open: FMCW under 15.245, one-way radar under 97.111; until settled use 15.249 levels (range ×0.32).
- **Radar** (4-radio FDD): TX3 into a horn ≥1.5 m aside (≥60 dB isolation), 20 MHz chirps at 20 MSPS (7.5 m bins), or 30 MHz at 30.72 MSPS with 8-bit samples (5 m); 26 MHz at most in 902-928 MHz. Range for 13 dB SNR in 10 ms at 18.75 dBm EIRP (at full drive that needs a ≥22 dBi horn at 5.8 GHz and a ≥14 dBi Yagi at 915 MHz): noise only / with the AD9361's uncorrelated TX/RX synthesizers (-88 dBc/Hz, 10 dB beam rejection) / with a 20 dB null on the horn:

| Target 1 m² | Noise | Synth floor | Null |
| --- | --- | --- | --- |
| 5.8 GHz, 20 dBi horn (10 m²; 0.01 m²) | 610 m (1.1 km; 190 m) | 340 m (610; 110 m) | 600 m |
| 2.44 GHz, 15 dBi patch | 630 m | 240 m | 570 m |
| 915 MHz, Yagi, 4 whips (50 dB isolation) | 1.2 km | 100 m | 320 m |

- Isolation and nulling matter more than power. Passive radar is licence-free: one reference plus three surveillance channels see 10 m² at 5-10 km on DVB-T. HackRF-only: the LMX2572 sends CW, FSK or ramps through RADAR OUT (+10 dBm EIRP) while the HackRF listens on the panel beam: CW Doppler of 1 m² to 300-700 m in 100 ms.

## 8. Control and software

**Lane 0, core (7× 74HC595):** XPWR[row] 0-3; GSEL[g] 4-15 (0xx TERM, 1mm port m); XCAL 16; XTX[k] 17-20; SEL[k] 21-24; BAND 25-27; PS[k] 28-31; AMP[k] 32-35; POL[k] 36-39; MODE 40; HTX 41; RTX 42; RSYN 43; RBAND 44-45; CDST 46-47; CDIR 48; MON 49-50; SYN 51; DET 52; XBIAS 53.
**Lane 1, tile (13× 74HC595):** RX[e] 0-15 (LNA path and power); CAL[e] 16-31; TTD bit b of element e at 32+4e+b; QAMP[q] 96-99; PRB 100-101.
All zeros is safe; 10k pull-downs on all 160 outputs.

**Firmware** (fwsim.py stays the specification): both lanes shift in parallel ("out pins, 2"): 104 clocks, 13.4 µs, minimum slot 15 µs (FIRST_TICK_US and the protocol's minimum slot become 15 in fwsim, ctrl and protocol); 20-byte words; LOAD_SEQ2 carries a host-defined marker so LNAs never power-cycle, and up to 256 words with a u16 count, as rev A's protocol version 2 does; rev A's SET_WORD and LOAD_SEQ (protocol versions 1 and 2) are answered with ERR_CMD, because rev A's 7-byte words would drive unrelated rev B functions (bit 53 is XBIAS). New: SET_WORD2, LOAD_SEQ2, TRIG, TX_ARM, SYNTH, PN, DSA, DET, TEMP, INFO, EEPROM; STATUS gains tile, lock and temperature fields. Interlocks: transmit bits need TX_ARM; radiating calibration refused in 15.205 bands or above -60 dBm EIRP outside MANIFOLD mode; HTX only with amplifiers bypassed; AMP refused in B1 below 0.5 GHz while an X port transmits; words beyond the USB-C current are refused, counting XBIAS at its 50 mA per-port limit; the host turns bias tees and antenna power off when it opens the board. A bit table with encodings lives in hw/design.py; the HackRF driver is powered when MODE and not HTX; the LMX2572 and both PE43711B-Z share one SPI bus with separate LE pins. USB-C current is read from a sink controller over I2C (the RP2040's four ADC inputs cannot take both ports' CC lines and the detector). CLOCK_SRC 0 becomes the board's 10 MHz ÷10.

**Host and simulator:** board.py splits into core.py and tile.py, asserted by hw/design.py; api.py gains a band/gain/region planner, arm_tx and external-port geometry; protocol.py and backend.py gain the new commands and STATUS fields (frame limit ≥5,146 bytes for LOAD_SEQ2: slot, repeat and count u16, the 20-byte marker and 256 words of 20 bytes); fwsim.py is updated first as the specification; hw/design.py becomes a two-board generator; control.py takes state_lag/toggle_factor from measured tables and adds Paley-20 and Latin-square schedules; dsp.py adds long-code alignment, calibration-slot monitoring, loopback/REV, manifold MUSIC, FMCW and passive ECA; new calib.py, hw_bladerf.py (clocking, J51, register writes, metadata) and hw_hackrf.py. sim.py gains the cascade, delay tables, external geometries, both calibration planes and openEMS patterns.

## 9. Expected performance

| Job | 0.4-1 GHz (external) | 1-2.3 GHz (external) | 2.3-6 GHz (panel) |
| --- | --- | --- | --- |
| DF | fixed layouts, 3° residual, median / 95th: 16 whips 0.2-0.5° / 0.6-1.4°, 8 whips 0.3-0.7° / 0.9-2.1° (5° and 0 dB: 0.8-1.8° / 2.5-5.3°); 4 whips, any modulation: 0.7-1.3° / 2.2-3.7°; outdoors 1-3° | same fixed layouts: 16 whips 0.1-0.2° / 0.25-0.6°, 8 whips 0.15-0.3° / 0.45-0.9°; λ/2 4×4 ≈0.5-0.8° / 0.9-1.5° in 2-D; restricted bands +1-2° | simulated 0.12/0.49, 0.12/0.49, 0.10/0.46° at 2.4/3.5/5.8 GHz (HackRF 0.55/1.73, 0.45/1.13, 0.29/0.63°); with coupling and manifold 0.5-1° (1-2° at 2.4 GHz), +1-2° in restricted slices without a measured manifold (3.6-5.15 GHz except 4.4-4.5) |
| Receive beams | 4 whips: +6 dB, 3 nulls of 20-30 dB; λ/2 line 26° | 4×4 rows: 26° fans in real time | 6.4/10.7/15.5 dBi, 62/43/26°; delay beams (HackRF too) lose ≤0.1 dB, sidelobes -19/-14/-11.5 dB; 3 nulls |
| Transmit steering | X0-X3, continuous phase, ≤+15.9 dBm EIRP | amateur only, ≤+13.4 dBm | 16 elements to ±45°, +8.1 (2.44) / +6.2 dBm (5.8) |
| Radar | 915 MHz 15.245: 100-320 m (1 m²; ×0.32 at 15.249 levels); passive DVB-T 5-10 km | passive LTE; active only at 2300-2310 MHz (amateur) | 5.8 GHz 340-600 m, 2.44 GHz 240-570 m (×0.32 at 15.249 levels) |

**Given up.** 0.4-1 GHz: the 4-whip simultaneous mode still rebuilds its ring per ~2:1 band (8 and 16 whips use fixed layouts); real-time beams and transmit use four elements; HackRF DF needs steady carriers. 1-2.3 GHz: no licence-free transmit or active radar; 55 % restricted, so no radiated calibration. 2.3-6 GHz: front ±60°, one polarisation; 2.3-3 GHz rests on the manifold; panel EIRP ≤+9.4 dBm.

## 10. Cost and power

| Item | $ |
| --- | --- |
| Core: front ends 180, channels 170, other blocks 95, SMAs and passives 40, PCB + PCBA 200 | 685 |
| Tile: cells 233, QAMP 42, cards 25, logic 20, PCB + PCBA 180, frame 25 | 523 |
| Cables | 90 |
| **Boards and parts** | **≈1,300** |
| Needed to calibrate it: LiteVNA 120, turntable 50, 2.3-6 GHz source antenna 20-60, 5-10 m cable 30 | ≈250 |
| **Working, calibrated system** | **≈1,550** (≈1,650 with the ring kit, 125, if no whips are at hand) |
| Optional: cavity filter and radar horn | ≈120 |
| 2× bladeRF 2.0 micro xA4 (HackRF One $340 if not owned) | 1,080 |

| Mode, at the USB-C input | W | A at 5 V |
| --- | --- | --- |
| Panel: QAMP / QAMP+AMP / HackRF | 7.3 / 8.6 / 9.0 | 1.5-1.8 |
| External: 4 whips / 16 ports | 1.7 / 6.0 | 0.3 / 1.2 |
| Everything on (refused on 1.5 A) | 13.3 | 2.7 |

Typical 7.3-9.0 W, 13.3 W at most. 5 V/3 A covers every mode only with a power path rated for it: a TPS2121 power mux (4.5 A) and a ≥3 A 5 V converter, not the 2.5 A LM66200 and the TPS63070 of the first parts list. Never request 12 V over PD. Buy PMA2-63LN+ (LCSC holds 9-20) and TCM1-83X+ (LCSC C879838; DigiKey has none before the end of December 2026) before layout.

## 11. What carries over

The RP2040 firmware core (shift registers, PIO, DMA, framing, marker, triggers) with fwsim.py and its tests; BGS12WN6; switched-line delays; the quadrant-MODE-HackRF topology; the clock tree; hw/design.py; the Array API, Hadamard/blind decoding, PN alignment and the simulator. Every open issue in ARCHITECTURE.md closes. Dropped: patches, SKY65404-31, quarter-wave Wilkinsons, centre probe, U.FL links, USB hub.

## 12. Build order

0. **Before any PCB (40 h):** two bladeRFs, a splitter and a PN source on the Mac: REFIN, CLKOUT→CLKIN, J51, FDD, manual gain, tracking freeze, 4×20 MSPS, overrun detection. Four whips on the RX ports already give DF.
1. **Firmware v2** on a bare Pico with 74HC595 LEDs, against fwsim (30 h).
2. **Simulate (60 h):** EMerge Floquet cell of the connected Vivaldi (active S11 against scan); openEMS 6×6, including the card-to-tile transition and the radome, with the loaded ring (8 runs by symmetry, Python 3.14) for embedded patterns, coupling, efficiency and cross-pol over εr 4.2-4.6 and height 35-45 mm; scikit-rf for delay states, idle-arm resonances and LNA stability into reflective shorts.
3. **Core board (180 h):** 16 ports, live calibration, HackRF mode, transmit, passive radar.
4. **Two coupon cards** (with the real transition) on a LiteVNA against openEMS, then the tile, its 36 hand-soldered card joints, and the manifold (180 h).
6. **Host software** for rev B, across steps 1-5 (about 100 h): fwsim, protocol, backend and api changes, calib.py, hw_bladerf.py, hw_hackrf.py, the two-board generator.
5. **Transmit loopback, 2f0 check, active radar** (40 h).

## 13. Decisions for the owner

1. **External ports:** 16 (recommended) or 8 (-$80, -1.8 W, no λ/2 grids).
2. **QAMP:** fit (recommended; +$42, +1.2 W, 1.6 dB at 5.8 GHz).
3. **Transmit ports:** X0-X3 (recommended) or all 16 (12 more T/R pairs and their control bits).
4. **Tile link:** SMA jumpers (recommended) or blind-mate SMP (+$250-300).
5. **Card:** 2-layer antipodal; balanced 3-layer only if cross-pol fails in openEMS.
6. **Power:** 5 V/3 A with a 3 A-rated power path (recommended). 9 V PD would need a PD sink controller and a 9 V-rated input path.
7. **Licence:** Technician, if radar or transmit power matters.
8. **Region:** US or EU profile.
9. **Radar hardware:** passive first; cavity filter and horn later.

## 14. How the proposals compared

| Proposal | Perf. | Build | Risk | Mean | Taken |
| --- | --- | --- | --- | --- | --- |
| P3 Backplane-and-Tiles | 8 | 6.5 | 8 | **7.5** | core/tile split, preselector as TX filter, channel polarity, EEPROM, TX_ARM, A→B clock |
| P4 digital-first | 6 | 7.5 | 5.5 | 6.3 | PMA2-63LN+, Latin square, BFCN band, LMX2572, calibration slots, synthesizer-noise radar budget |
| P1 split-band hybrid | 6.5 | 4.5 | 7 | 6.0 | switched arms, one-switch polarity, REV/AD8318, TVS + HPF, probes, 16 ports |
| P2 minimal change | 3.5 | 4 | 4 | 3.8 | radiated-calibration caps, passive radar, rev A's frames, 5 V/3 A |

Fixed from P3: radiating taps, +9 dBm at the HackRF, two-switch polarity, 12 V PD, $304 of SMP, the 0.21 dB delay penalty, 1 dB-high low-band gain, thermal-only radar, the sub-noise pilot. Calculations: the design panel's working files (not kept in the repo)