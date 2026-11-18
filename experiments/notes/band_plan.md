<!-- Generated 2026-11-16 from the investigators' reports (workflow wf_6f0fba6f-2a0). UNREVIEWED: the synthesis and the adversarial
check for this study did not run (usage limit). Treat every number as a lead to verify, not as a result. -->

# Cherry-picked bands and lattices: research notes (unreviewed)

## investigate:lattice

Lattice lens: a lattice buys beamforming economy, not direction-finding (DF) performance, and periodic 3-D lattices do not help a 16-element array. Script: experiments/lattice_study.py (full run 363 s); outputs in experiments/lattice/. No existing repository file was modified.

**Principle.** A periodic array's array factor returns to full height at every reciprocal-lattice point G. A grating lobe enters real space when |G| <= (1+sin scan)/lambda. Two directions inside a field of view of half-angle t become indistinguishable when |G| <= 2 sin(t)/lambda. Only the shortest G per element matters.

**2-D panel (16 active elements).**
- **Equilateral triangular lattice:** shortest G is 1.1547/spacing, so spacing can be 15.5 % larger or elements 13.4 % fewer for the same grating-lobe-free cone. This reproduces Sharp 1961.
- **Staggering rev B's cards by half a pitch (25 x 25 mm):** this is an isosceles lattice, and its shortest G is unchanged. Grating-lobe-free frequency on a 60 degree cone stays 6.43 GHz (6.00 GHz at 90 degrees); only scans along the cards improve, to 7.9 GHz. It gives essentially no benefit.
- **True triangular lattice with cards still 25 mm apart:** flare pitch 28.87 mm, alternate cards shifted 14.43 mm. It has the same limits as the square 25 mm panel (6.43 / 6.00 GHz; DF unambiguous to 6.92 GHz within 60 degrees), with 15.5 % more aperture (+0.62 dB) and about 5-8 % smaller DF error (95th percentile 0.64 versus 0.67 degrees at 6 GHz).
- **Square panel of the same area (26.9 mm):** grating-free only to 5.98 / 5.58 GHz, with ambiguity level 0.88 at 6 GHz.
- **Delays:** the triangular 28.9 mm panel needs 363 ps for 60 degrees on its worst azimuth, against 347.5 ps available (the square needs 306 ps). Worst quantisation loss at 5.8 GHz is 0.44 dB versus 0.30 dB.
- **Risk:** Ellgardt and Wikstrom 2009 report that single-polarised tapered slots on a triangular grid support guided modes that cause scan blindness, so the usable enlargement is under 15 %. This is from the abstract; I did not read the full paper.
- **Verdict:** +0.6 dB is not worth a new blindness risk on a hobby panel. Restricting band or scan pays far more: the square pitch allowed is lambda/(1+sin scan), which is 30.0 mm at 5.85 GHz for 45 degree scan and 46.2 mm at 3.8 GHz for 45 degrees. A triangular lattice is 15.5 % larger again.

**16 whips, 0.4-2.3 GHz azimuth DF** (same check model and conditions as the layout search; "clean" means gross errors <= 0.5 % for sources up to 15 and 30 degrees elevation).

| Layout | Clean range | Nominal 95th percentile error | Worst gross errors |
| --- | --- | --- | --- |
| Square 4x4, 113 mm | 0.4-1.30 GHz | 1.45 deg at 0.4, 0.48 deg at 1.2 GHz | 30 % above 1.35 GHz |
| Square 4x4, 85 mm | to 1.75 GHz | 1.9 deg at 0.4 GHz | 24 % |
| Triangular, 120 mm | to 1.40 GHz | 1.37 deg at 0.4 GHz | 35 % |
| Triangular, 75 mm | whole band | 2.19-0.45 deg | 1.5 % (0 dB, 5 deg residual) |
| Two rings of 8, 225/112 mm (rev B) | whole band | 1.42-0.24 deg | 0.0 % |
| Searched, random16 | whole band | 1.26-0.23 deg | 0.3 % |

- Lattices fail at c|G|/2: 1.33 GHz for the 113 mm square and 1.44 GHz for the 120 mm triangular lattice. Below that, the 113 mm lattice only matches the two rings.
- The 75 mm triangular lattice is the one lattice that survives the band (its ambiguity frequency is 2.31 GHz, and 75 mm is also the plate's minimum spacing), but it is about 1.5-1.9 times less accurate than the rings.
- For DF with few elements a lattice never wins. It is useful only for real-time fan beams (the lambda/2 grids already in rev B) or below a chosen top frequency.

**3-D lattices.**
- **Shortest G times nearest-neighbour spacing:** simple cubic 1.000; body-centred cubic, face-centred cubic and hexagonal close-packed 1.2247 (close-packed also has a -6 dB partial lobe at 1.1547).
- **Elements per volume for the same ambiguity-free frequency:** body-centred cubic needs 29.3 % fewer than cubic, face-centred and close-packed 23 % fewer. Body-centred cubic is the 3-D analogue of Sharp's result (Petersen and Middleton 1962); the percentages are my derivation.
- **16-element pieces:** ambiguity level passes 0.9 at spacing/lambda of 0.43-0.45 (cubic), 0.47 (planar square), 0.54 (planar triangular), 0.55-0.56 (body-centred, face-centred, close-packed). The aperture is only 2-2.3 spacings per axis, against 3-3.5 for a planar 4x4.
- **Verdict:** a 3-D lattice aliases exactly as a planar one does and spends elements on depth. The gain from 3-D is non-coplanarity (elevation, full sphere), not periodicity. Whips on a ground plate cannot be stacked, and end-fire Vivaldis would shadow each other.
- **Where periodic lattices pay:** MIMO-radar virtual arrays, where a filled uniform virtual grid is the aim, and large arrays with many elements. The sonar and ultrasound usage was only weakly sourced.

### Findings

- **Sharp's triangular-lattice result** (high): An equilateral triangular arrangement needs 13.4 % fewer elements than a square one when the beam scans within a cone about the normal (5.7-13.4 % for a square pyramid of scan). Equivalent to 15.5 % larger spacing. Reproduced numerically: shortest reciprocal vector 1.1547/spacing, element ratio 0.866. Source: https://doi.org/10.1109/TAP.1961.1144967 (abstract via https://www.bohrium.com/paper-details/a-triangular-arrangement-of-planar-array-elements-that-reduces-the-number-needed/812351807132008449-20293); experiments/lattice/results.json (theory)
- **Staggering rev B's cards by half a pitch** (high): With card spacing 25 mm and flare pitch 25 mm, a half-pitch stagger gives an isosceles lattice whose shortest reciprocal vector is unchanged (1/25 mm). Grating-lobe-free frequency on a cone stays 6.43 GHz at 60 deg and 6.00 GHz at 90 deg; only scans along the cards improve (7.9 GHz at 60 deg). Sidelobe and DF results are the same as the square panel within noise. Source: experiments/lattice/results.json (panel)
- **True triangular panel on 25 mm card spacing** (medium): Flare pitch 28.87 mm with alternate cards shifted 14.43 mm is an equilateral lattice with the same limits as square 25 mm (grating-free 6.43/6.00 GHz at 60/90 deg; DF unambiguous to 6.92 GHz within 60 deg), 15.5 % more aperture (+0.62 dB) and 95th percentile DF error 0.64 vs 0.67 deg at 6 GHz (1.55 vs 1.69 deg at 2.4 GHz). A square panel of the same area (26.9 mm) is grating-free only to 5.98/5.58 GHz and has ambiguity level 0.88 at 6 GHz. Source: experiments/lattice/results.json (panel); isotropic elements, no coupling
- **4-bit delays on a triangular panel** (medium): With rev B's delay set (38.5/52.5/117.5/139 ps, 347.5 ps maximum), the triangular 28.9 mm panel needs 363 ps for 60 deg on its worst azimuth (square 25 mm: 306 ps; staggered 25x25: 332 ps). Worst quantisation loss at 5.8 GHz and 60 deg is 0.44 dB vs 0.30 dB for the square; mean 0.26 vs 0.22 dB. The delay set would need about 15 ps more range, or scan limited to about 56 deg on that diagonal. Source: experiments/lattice/results.json (panel, ttd_5p8GHz)
- **Triangular-grid Vivaldi risk** (medium): Ellgardt and Wikstrom (IEEE TAP 57(9), 2009, 2599-2607) report that single-polarised tapered slots in a triangular grid support guided modes that cause scan blindness, so the unit-cell enlargement is limited to less than the optimal 15 %; their 16-element prototype deviated from infinite-array predictions in the H-plane at large scan. Read from the abstract only. Source: https://doi.org/10.1109/TAP.2009.2027044 (abstract via https://manuscript.isc.ac/Inventory/49/1128376.htm)
- **Pitch allowed by restricting band and scan** (high): Maximum square pitch lambda/(1+sin scan): 6.0 GHz 25.0/26.8/29.3 mm for 90/60/45 deg; 5.85 GHz 25.6/27.5/30.0 mm; 5.0 GHz 30.0/32.1/35.1 mm; 3.8 GHz 39.4/42.3/46.2 mm; 2.4835 GHz 60.4/64.7/70.7 mm. A triangular lattice allows 15.5 % more. Dropping the top band or scan range buys far more aperture per element than the lattice choice. Source: experiments/lattice/results.json (max_square_pitch_mm); formula is standard grating-lobe geometry
- **16-whip lattices versus rings and searched layouts** (high): On the 500 mm plate, lattices give gross DF errors above c|G|/2: square 113 mm clean 0.4-1.30 GHz then up to 30 %; square 85 mm clean to 1.75 GHz (24 %); triangular 120 mm clean to 1.40 GHz (35 %). Two rings of 8 (225/112 mm) and random16 are clean over 0.4-2.3 GHz (0.0 % and 0.3 % worst, including 0 dB / 5 deg residual). Below its limit the 113 mm lattice matches the two rings (95th percentile 1.45/0.63/0.48 deg vs 1.42/0.63/0.47 deg at 0.4/0.9/1.2 GHz); it never beats them. Source: experiments/lattice/results.json (whips); experiments/lattice/whips.png
- **The one lattice that covers 0.4-2.3 GHz** (high): A triangular lattice at 75 mm (the plate's minimum spacing) has its azimuth ambiguity frequency at 2.31 GHz and is clean over the whole band (gross errors 0.2/0.2/1.5 % worst for nominal/elevated/stress), but its 95th percentile error is 2.19-0.45 deg vs 1.42-0.24 deg for the two rings, because its RMS radius is 119 mm vs 178 mm. Source: experiments/lattice/results.json (whips)
- **3-D lattices: reciprocal vectors and density** (medium): Shortest full-height reciprocal vector times nearest-neighbour spacing: simple cubic 1.000, BCC 1.2247, FCC 1.2247, HCP 1.2247 (plus a -6 dB partial lobe at 1.1547). Elements per volume for the same ambiguity-free frequency relative to cubic: BCC 0.707 (29.3 % fewer), FCC and HCP 0.770 (23 % fewer). BCC is the optimal 3-D sampling lattice (Petersen and Middleton 1962); the percentages are computed here, and a 2014 survey notes they have been misreported in the literature. Source: experiments/lattice/results.json (theory); https://doi.org/10.1016/S0019-9958(62)90633-2; https://eprints.cs.univie.ac.at/4183
- **3-D lattices with 16 elements** (high): Ambiguity level (highest false peak over true peak, any two directions on the sphere) passes 0.9 at spacing/lambda of 0.43 (cubic 2x2x4), 0.45 (cubic compact), 0.47 (planar square 4x4), 0.54 (planar triangular), 0.55 (BCC), 0.56 (FCC, HCP). A 16-element 3-D piece spans only 2.0-2.3 spacings per axis (HCP 1.6 in depth) against 3.0-3.5 for a planar 4x4, so it has a smaller aperture in every direction and aliases at essentially the same spacing as the planar triangular lattice. Source: experiments/lattice/results.json (lattices_3d); experiments/lattice/lattices_3d.png
- **How volumetric array patterns are analysed** (medium): A volumetric array's pattern is its 3-D pattern space evaluated on a spherical surface whose position is set by the scan angle (Barott and Steffes 2008); they use it to thin volumetric arrays for low sidelobes, and volumetric arrays can keep beam resolution roughly constant over steering angle. A thesis on random spherical volumetric acoustic arrays argues the non-periodic layout is what removes large grating lobes. Source: https://doi.org/10.1109/LAWP.2008.2005360; https://hull-repository.worktribe.com/output/4219778/design-and-analysis-of-an-acoustic-random-spherical-volumetric-array
- **Where periodic lattices are used** (low): MIMO radar builds a virtual array of N_TX x N_RX elements at the pairwise sums of transmit and receive positions, and layouts are chosen so that this virtual array is a filled uniform grid (TI SWRA554). Ultrasound matrix probes are 2-D lattices whose pitch sets the grating-lobe-free frequency. I found no source showing a periodic 3-D lattice used for a small DF array. Source: https://www.ti.com/lit/swra554; https://arxiv.org/html/2602.15174
- **Array ambiguity and geometry in the DF literature** (medium): Manifold ambiguities arise when distinct directions give linearly dependent steering vectors; uniform arrays have a regular ambiguity structure, and for a fixed aperture and element count the geometry with the highest lower bound on response-vector dissimilarity performs best (not limited to linear arrays). Consistent with the search finding irregular layouts and rings ahead of lattices. Source: https://cris.tau.ac.il/en/publications/array-geometry-for-ambiguity-resolution-in-direction-finding/; https://arxiv.org/pdf/2110.10756.pdf

### Files

- `experiments/lattice_study.py`
- `experiments/lattice/results.json`
- `experiments/lattice/run.log`
- `experiments/lattice/panel.png`
- `experiments/lattice/panel_layouts.png`
- `experiments/lattice/whips.png`
- `experiments/lattice/lattices_3d.png`

### Open questions

- Does a triangular-grid connected Vivaldi panel at 28.9 mm flare pitch show the E-plane scan blindness Ellgardt reports inside 2.3-6 GHz and +/-60 deg? Only an openEMS unit-cell model can say; the script uses isotropic elements with no coupling.
- If the owner restricts the panel's top band or scan range (for example 5.85 GHz and +/-45 deg, which allows 30 mm square or 34.7 mm triangular pitch), is the larger aperture worth redesigning the flare and the delay set (the triangular 28.9 mm panel already needs 363 ps against 347.5 ps available)?
- Is a lattice wanted for the whips at all? It pays only for real-time fan beams below its ambiguity frequency (113 mm square: 1.3 GHz; 75 mm triangular: 2.31 GHz); for DF the two rings match or beat it everywhere.
- The sonar and ultrasound use of periodic volumetric lattices was only weakly sourced here (one arXiv paper on 2-D matrix probes, one thesis on random volumetric arrays); a primary sonar reference was not found.
- The panel DF Monte Carlo searches the whole front hemisphere with sources inside 60 deg and ignores the element pattern and mutual coupling; the 6.5-7 GHz gross-error figures are indicative only.
- Only the abstracts of Sharp 1961, Ellgardt and Wikstrom 2009 and Petersen and Middleton 1962 were reachable; the 29.3 % and 23 % 3-D savings are my own derivation from the reciprocal lattices, not quoted figures.

## investigate:codesign

The band-geometry study was already complete on disk from the earlier run (script 672 lines, full run about 26 min: part 1 232 s, parts 2-3 1310 s). I did not re-run it; I checked results.json for completeness, re-checked key rows with fresh seeds (4000 trials, they agree), and added one azimuth-plus-elevation estimator test in scratch.

**Headline: geometry is not what makes 0.4-6 GHz expensive, and cherry-picking bands buys the whip layout about 10 % at most.** Two fixed geometries already tile the range. All numbers are from the idealised model (isotropic elements, no coupling, 3 deg / 0.5 dB calibration residual, 20 dB SNR, 64 samples), so they are upper bounds.

**1. Geometry -> bands (0.3-8 GHz, 80 frequencies)**

| Geometry | DF p95 <= 1 deg | DF p95 <= 2 deg | No gross errors | Clean beams (sidelobe <= -6 dB, width <= 60 deg) |
|---|---|---|---|---|
| rings16 225/112 mm | 0.56-2.72 GHz | 0.30-2.72 | 0.30-2.72 | 0.61-0.92 |
| random8 | none | 0.44-2.50 | 0.30-2.50 | none |
| half-wave 4x4 whips, 164 mm (915 MHz; does not fit the plate) | 0.40-0.88 | 0.30-0.88 | 0.30-0.88 | 0.42-0.75 |
| half-wave 4x4 whips, 95 mm (1575 MHz) | 0.69-1.52 | 0.34-1.52 | 0.30-1.52 | 0.72-1.28 |
| panel 4x4 at 25 mm | 3.63-6.50 | 1.87-6.50 | 0.58-6.50 | 2.50-6.23 |

- **Tiling:** rings16 plus the panel cover 0.30-6.50 GHz at <= 2 deg with a 1.87-2.72 GHz overlap, so rev B's 2.3 GHz crossover sits inside it. At <= 1 deg there is a gap at 2.72-3.63 GHz and below 0.56 GHz.
- **At the bands (p95 nominal):** rings16 gives 1.32 deg at 433, 0.64 at 915, 0.37 at 1575 and 0.24 at 2440 MHz. The panel gives 1.57 at 2440, 1.07 at 3600 and 0.66 at 5800 MHz. At 2.44 GHz the whip plate is about 6x more accurate than the panel, in azimuth only.
- **Whips are DF layouts, not beam layouts:** above about 1 GHz the sparse layouts have -3 dB sidelobes. Half-wave grids give clean beams over less than an octave and go fully ambiguous just above their design frequency (the 164 mm grid already shows 1.1 % gross errors at 915 MHz).
- **The whips' 2.5-2.7 GHz ceiling is the estimator, not the layout.** Above it, gross errors come from sources 15-30 deg up breaking the horizon assumption (rings16: 13 % at 3.6 GHz, 32 % at 5.8 GHz). With an az+el scan the same layouts show 0.0 % gross errors and p95 about 0.5 deg at 3.6 and 5.8 GHz.

**2. Bands -> geometry (6000-evaluation searches, finalists re-selected on 12 draws, 3000-trial check)**

| Band set, 16 whips | Full-band rings16 | Best cherry-picked | Gain |
|---|---|---|---|
| A {433, 915, 2440}: worst p95 nominal / stress | 1.33 / 3.34 deg | 1.19-1.24 / 2.98-3.14 deg | 6-11 % |
| B {915, 1575, 2440} | 0.63 / 1.59 deg | 0.59-0.61 / 1.48-1.49 deg | 3-7 % |
| C {2440, 5800}: gross errors, sources 0-30 deg up | 34 % | 9.1 % (75 mm rule); 0.2 % (30 mm rule) | large, but only for the az-only estimator |

- **Single-band bound:** at 433 MHz alone the best layout is all 16 whips on the 240 mm rim, 1.00 deg against 1.33 (25 %), but that ring fails at 2.44 GHz (11.5 % gross errors).
- **Set C cost:** the 30 mm-spacing layout that fixes 5.8 GHz has nominal p95 of 1.0 deg against 0.24 at 2.44 GHz.
- **8 whips:** cherry-picking only trades. Set A removes random8's 0.9 % gross errors at 2.44 GHz but worsens 433 MHz stress p95 from 5.15 to 5.55 deg; set B is worse on accuracy (0.99 to 1.28 deg at 915 MHz).
- **Method:** random search beat hill-climbing from rings16 in 4 of 6 cases; the cherry-picked 16-whip layouts still hold p95 <= 1.4 deg across 0.4-2.3 GHz.

**3. Panel pitch**

| Pitch | 6x6 tile | 2.44 GHz: directivity / beamwidth / p95 | 5.8 GHz: directivity / p95 | Usable to |
|---|---|---|---|---|
| 25.0 mm (rev B) | 150 mm | 10.1 dBi / 60 deg / 1.49 deg | 16.7 dBi / 0.70 deg | 6.4 GHz |
| 25.4 mm (half-wave at 5.9 GHz) | 152 mm | 10.2 / 59.5 / 1.48 | 16.9 / 0.69 | 6.3 GHz |
| 27.2 mm (+/-60 deg scan at 5.9 GHz) | 163 mm | 10.7 / 55.9 / 1.39 | 17.5 / 0.63 | 5.9 GHz |
| 60 mm (half-wave at 2.5 GHz) | 360 mm | 16.8 / 26.5 / 0.62 | ambiguous (73 % gross) | 2.7 GHz |

Capping the top at 5.9 GHz changes nothing; 25 mm is already that pitch. Dropping 5-6 GHz is the one geometric lever: 60 mm buys +6.7 dB and 2.4x DF accuracy at 2.44 GHz, for a tile 2.4x wider and nothing above about 2.7 GHz.

The script docstring points to a band_plan/README that does not exist.

### Findings

- **Study status** (high): experiments/band_plan_study.py (new file, 672 lines, imports evolve_layout, wideband_study and array16.evolve read-only) had already run to completion: part 1 took 232 s, parts 2-3 took 1310 s (about 26 min total, OMP_NUM_THREADS=4). results.json contains part1, part2 (sets, single, relaxed; budget 6000) and part3. A fresh-seed re-check (4000 trials) reproduced the key rows to within a few percent. Source: experiments/band_plan/results.json; experiments/band_plan/run_parts23.log; scratch/band_plan/lens_sim/spot.py
- **Geometry -> bands: windows** (high): DF windows (p95 nominal within the limit and gross errors <= 1 % for sources 0-15 and 0-30 deg up): rings16_225_112mm <= 2 deg over 0.30-2.72 GHz and <= 1 deg over 0.56-2.72 GHz; random8 <= 2 deg over 0.44-2.50 GHz, never <= 1 deg; 25 mm 4x4 panel <= 2 deg over 1.87-6.50 GHz, <= 1 deg over 3.63-6.50 GHz, <= 5 deg from 0.78 GHz; half-wave 4x4 whip grids: 164 mm pitch 0.30-0.88 GHz, 95 mm pitch 0.34-1.52 GHz. Source: experiments/band_plan/tables.txt
- **Geometry -> bands: tiling** (high): rings16 plus the 25 mm panel tile 0.30-6.50 GHz at <= 2 deg with a 1.87-2.72 GHz overlap, so rev B's 2.3 GHz crossover is inside the overlap. At <= 1 deg there is a gap at 2.72-3.63 GHz (panel p95 is 1.07 deg at 3.6 GHz) and below 0.56 GHz. Source: experiments/band_plan/tables.txt
- **Geometry -> bands: accuracy at the bands** (high): p95 error, sources 0-15 deg up: rings16 1.32 deg at 433 MHz, 0.64 at 915, 0.37 at 1575, 0.24 at 2440. Panel 8.70 deg at 433, 4.21 at 915, 1.57 at 2440, 1.07 at 3600, 0.66 at 5800. At 2.44 GHz the 500 mm whip plate is about 6x more accurate than the 100 mm panel, in azimuth only. Source: experiments/band_plan/tables.txt
- **Geometry -> bands: beams** (high): Windows with peak sidelobe <= -6 dB and beamwidth <= 60 deg: panel 2.50-6.23 GHz (directivity 9.9 dBi at 2.44, 13.0 at 3.5, 16.7 at 5.8 GHz with a cos element; beamwidth 61, 44, 27 deg); rings16 only 0.61-0.92 GHz (sidelobes about -3.6 dB from 1.5 GHz up); random8 none; half-wave grids 0.42-0.75 GHz (164 mm) and 0.72-1.28 GHz (95 mm). Half-wave grids become fully ambiguous (0 dB grating lobe on the horizon) just above their design frequency; the 164 mm grid already shows 1.1 % gross errors at 915 MHz. Source: experiments/band_plan/tables.txt
- **Whips' upper frequency limit is the estimator** (medium): Above 2.5-2.7 GHz the whip layouts have 0 % gross errors for sources 0-15 deg up (rings16 p95 0.10 deg at 5.8 GHz) but 13 % at 3.6 GHz and 32 % at 5.8 GHz for sources 0-30 deg up, because the azimuth-only estimator assumes the horizon. With an az+el Bartlett scan (1 deg grids, 600 trials) rings16, random8 and random16 all give 0.0 % gross errors and azimuth p95 of 0.47-0.61 deg at 2.44, 3.6 and 5.8 GHz. Idealised: isotropic elements, no coupling, and rev B's whips and calibration are only specified to 2.3 GHz. Source: scratch/band_plan/lens_sim/spot.py (output in this session)
- **Bands -> geometry: 16 whips, sets A and B** (high): Set {433, 915, 2440 MHz}: rings16 worst-band p95 is 1.33 deg nominal / 3.34 stress; the layouts optimised for the set give 1.19-1.24 / 2.98-3.14 (6-11 % better; fresh seed 1.32 vs 1.21). Set {915, 1575, 2440}: 0.63 / 1.59 against 0.59-0.61 / 1.48-1.49 (3-7 %). Gross-error rates are 0.0-0.3 % for both, indistinguishable at 3000 trials. The cherry-picked layouts still hold p95 <= 1.4 deg over 0.4-2.3 GHz. Source: experiments/band_plan/tables.txt
- **Bands -> geometry: single-band bound** (high): For 433 MHz alone the best 16-whip layout is all whips on the 240 mm rim: p95 1.00 deg against 1.33 for rings16 (25 % better), but that ring has 11.5 % gross errors at 2.44 GHz for sources 0-30 deg up. For 2.44 GHz alone the optimised layout (0.27 deg) is no better than rings16 (0.24 deg). So the most cherry-picking can gain in accuracy on this plate is about 25 % in one band. Source: experiments/band_plan/tables.txt
- **Bands -> geometry: set {2440, 5800}** (high): With the azimuth-only estimator rings16 has 34 % gross errors at 5.8 GHz (sources 0-30 deg up). Optimising for the set under the 75 mm spacing rule reaches 9.1 % (fresh seed 9.6 %); relaxing the rule to 30 mm reaches 0.2 % (fresh seed 0.1 %) at the cost of nominal p95 1.0 deg against 0.24 deg at 2.44 GHz. The az+el estimator removes the problem without changing the layout. Source: experiments/band_plan/tables.txt; scratch spot.py
- **Bands -> geometry: 8 whips** (high): Cherry-picking trades rather than gains. Set A: gross errors at 2.44 GHz 0.9 % -> 0.0 % and stress p95 there 1.72 -> 1.05 deg, but 433 MHz stress p95 5.15 -> 5.55 deg. Set B: worst nominal p95 0.99 -> 1.28 deg (worse), gross 0.9 % -> 0.0 %. Random search beat hill-climbing from the full-band layout in 4 of 6 set searches. Source: experiments/band_plan/tables.txt
- **Panel pitch** (high): 25.0 mm (rev B), 25.4 mm (half-wave at 5.9 GHz) and 27.2 mm (grating-lobe-free to +/-60 deg at 5.9 GHz) are equivalent: 16.7 / 16.9 / 17.5 dBi and p95 0.70 / 0.69 / 0.63 deg at 5.8 GHz; 27.2 mm cuts the grating-lobe-free scan at 5.8 GHz to 64 deg. A 60 mm pitch (top band 2.5 GHz) gives 16.8 dBi against 10.1 dBi, beamwidth 26.5 against 60.3 deg and p95 0.62 against 1.49 deg at 2.44 GHz, with a 360 mm 6x6 tile against 150 mm, and is ambiguous above about 2.7 GHz (27.6 % gross errors at 3.5 GHz, 73 % at 5.8 GHz). Source: experiments/band_plan/tables.txt
- **Band edges used in the simulation** (high): FCC 15.247 applies to 902-928 MHz, 2400-2483.5 MHz and 5725-5850 MHz. Source: https://www.law.cornell.edu/cfr/text/47/15.247
- **Band edges used in the simulation** (medium): ITU ISM bands include 433.05-434.79 MHz (Region 1 only), 902-928 MHz (Region 2 only), 2400-2500 MHz and 5725-5875 MHz (worldwide). Source: https://en.wikipedia.org/wiki/ISM_radio_band
- **Band edges used in the simulation** (medium): GPS L1 carrier is 1575.42 MHz; the script's +/-10 MHz window around it is its own choice, not an allocation. Source: https://en.wikipedia.org/wiki/GPS_signals
- **Interpretation** (medium): Restricting or cherry-picking the frequency range does not relax the array geometry much: two fixed geometries already tile 0.3-6.5 GHz, and band-specific whip layouts gain 3-11 %. The geometric lever that restriction does open is panel pitch when 5-6 GHz is dropped (+6.7 dB directivity and 2.4x DF accuracy at 2.44 GHz for a 2.4x wider tile). Source: estimate: synthesis of the tables above; idealised model (isotropic elements, no coupling, no element or preselector bandwidth)

### Files

- `experiments/band_plan_study.py`
- `experiments/band_plan/tables.txt`
- `experiments/band_plan/results.json`
- `experiments/band_plan/fom_vs_frequency.png`
- `experiments/band_plan/cherry_layouts.png`
- `experiments/band_plan/panel_pitch.png`
- `experiments/band_plan/run.log`
- `experiments/band_plan/run_parts23.log`
- `scratch/band_plan/lens_sim/spot.py`

### Open questions

- Should the whip plate's DF estimator scan elevation as well as azimuth? In the idealised model that lifts the fixed layouts' gross-error ceiling from about 2.7 GHz to at least 5.8 GHz, but it is untested with real whip patterns, coupling, or calibration above 2.3 GHz, and was run at only 600 trials on a 1 deg grid.
- Does the owner care about 5-6 GHz enough to keep the 25 mm pitch? Dropping it allows a 60 mm pitch: +6.7 dB and 2.4x DF accuracy at 2.44 GHz for a 360 mm tile; a Vivaldi element for that pitch was not modelled.
- If 2.4 GHz moved to the whip plate (window reaches 2.72 GHz, 6x better azimuth accuracy than the panel), could the Vivaldi's low edge rise from 2.3 GHz? That depends on element and preselector design, not on this simulation.
- Only the band edges for 433, 915, 2440, 5800 MHz and the GPS L1 centre were checked against a source in this session; the other entries in the script's BANDS table (868, 1090, 1296, 1420, 1621, 3600, 5500) are used only in the part 1 at-bands table and are unverified here.
- The script docstring refers to a band_plan/README with the allocation citations, which does not exist; a comment in search_bands refers to a first-attempt log that lives only in the scratchpad.
- Each band-set search was run once per array size; the 3-11 % gains for sets A and B have no replicate runs behind them, so the ranking among search methods is not established.

## investigate:bands

Lens: which bands in 0.3-6 GHz matter, and what narrowing buys compared with rev B. The shared web-search budget was used up before this run started, so the research used direct fetches of known datasheet and regulation URLs, the official ERC 70-03 and LoRaWAN RP002 PDFs, and files earlier agents left in the scratchpad. Several new-part prices could not be checked and are marked as estimates.

CATALOGUE. The transmit bands are all narrow, each 7.3 % wide at most:
- **433 MHz:** 433.05-434.79 MHz in the EU (10 mW e.r.p., 10 % duty). In the US only short control bursts at about -14 dBm EIRP.
- **868 MHz (EU):** 863-870 MHz at 25 mW e.r.p. and 0.1-1 % duty; 869.40-869.65 MHz allows 500 mW.
- **902-928 MHz (US):** 15.249 allows -1.25 dBm EIRP, 15.245 sensors +18.75 dBm, 15.247 1 W.
- **2.4 GHz:** EU 10 mW, 25 mW for radar, 100 mW with listen-before-talk.
- **5.725-5.875 GHz:** EU 25 mW. US limits as for 902-928 MHz.
- **Amateur:** 1.5 kW PEP in the US. Pulse is allowed only at 33, 13 and 5 cm. The 9 cm band is gone from the US table.

REDESIGN ignores 47 CFR 15.247. Digitally modulated transmissions of 500 kHz or more, which includes the board's own PN codes, may use 1 W conducted, or 36 dBm EIRP with the 6 dBi rule, instead of 15.249's -1.25 dBm. Whether a home-built steering experiment can rely on 15.247 is untested.

Things to locate or exploit:
- **Direction-finding targets:** DJI O4 at 2.4/5.8 GHz (under 33 dBm EIRP under FCC rules, 14-20 dBm in the EU), LoRa and ELRS at 433/868/915 MHz, tyre sensors and key fobs at 315 and 433 MHz, Bluetooth direction-finding tags, Remote ID, Wi-Fi.
- **Transmitters at known positions, usable for calibration without radiating:** ADS-B at 978/1090 MHz, radiosondes at 400-406 MHz, Inmarsat (geostationary), GNSS L1/L2/L5.
- **Passive-radar illuminators:** DVB-T at 470-694 MHz, ATSC below 608 MHz, LTE/5G downlinks.

SETS. Each set lists its bands, the share inside 15.205 restricted bands, and which radio chips cover it:
- **S1, 2.4 + 5.8 GHz:** 0 % restricted.
- **S2, 400-450 + 863-928 MHz + 2.4 GHz:** 5 % restricted. The AT86RF215 and AD9363 cover all of it.
- **S3, 433 + 868/915 MHz + 2.4 + 5.8 GHz:** 0 % restricted.
- **S4, S3 plus receive-only bands** (400-406, 470-694, 975-1095, 1164-1300, 1400-1427, 1525-1590 MHz): 29 % restricted.
- **Rev B, 0.4-6 GHz:** 47.5 % restricted.

WHAT NARROWING BUYS:
1. **Radiated calibration becomes legal across the S1 and S3 bands.**
2. **Fixed filters replace the B1-B5 bank.** BFCN-2450+ loses 2.0 dB and BFCN-5750+ 1.85 dB, against rev B's 5.6-6.2 dB at 5.8 GHz. The cavity parts ZVBP-2450-S+ and ZVBP-5800-S+ lose 0.7/0.8 dB and give 50-80 dB of rejection, enough for the 11.6 GHz harmonic.
3. **Phase control:** true-time delay is not needed. Phase-only steering moves a 45° beam by at most 1° within any one ISM band. Digital phase shifters do not get cheaper: HMC649A covers 3-6 GHz at about $169, and no ADI part covers 2.4 GHz. The same four switched lines, re-chosen, do better:

   | Lines | RMS phase error | Worst case |
   |---|---|---|
   | Rev B lines over 2.3-6 GHz | 11.3° | 36° |
   | One set for 2.4 + 5.8 GHz (41/82/106/159 ps) | 8.1° | 18° |
   | Single-band sets | 6.6-7.1° (ideal 4-bit is 6.5°) | |

4. **Noise figure improves only at 5.8 GHz:** from 3.4 dB to 2.2 dB (1.8 dB without the cal arm) with two SKY65404-31 per element. Panel LNA power falls by about 5 W. Sub-GHz and 2.4 GHz are already about 1.1-1.9 dB in rev B.
5. **Aperture:** the 25 mm panel is only 0.8 λ wide at 2.44 GHz (6.4 dBi). A 2.4-only 4x4 at 60.4 mm reaches 16.9 dBi, a 26° beam, and +10 dB G/T (-20.2 to -10.3 dB/K). The sibling simulation gives a 95th-percentile direction error of 0.62° against 1.49°.
6. **Radios: 5.8 GHz is the band that forces $270-per-channel AD9361 radios.** Without it, the AT86RF215 fits: 389.5-510 / 779-1020 / 2400-2483.5 MHz, 13-bit 4 MSPS I/Q, at most 2 MHz bandwidth, noise figure 4.3-4.5 dB, a 26 MHz TCXO input, both bands at once. An AD9363 (325-3800 MHz, 20 MHz) also fits.
7. **The whip layout gains nothing.** The sibling study shows the 16-whip rings already cover 0.30-2.72 GHz, and layouts optimised for chosen bands are barely better.

COST AND EFFORT against rev B ($2,630 all-in, 520 h). All of these are estimates:

| Set | Cost change | Effort change |
|---|---|---|
| S1a: keep the Vivaldi, drop the external ports | -$335 | -70 h |
| S1b: two patch tiles (+10 dB G/T at 2.4 GHz) | -$123 | -100 h |
| S3 | -$228 | -100 h |
| S4 | -$228 | -90 h |
| S2a: core only, AD9363 radios (prices unchecked) | -$1,315 | -200 h |
| S2b: one AT86RF215 per whip, 16 coherent channels | -$1,580 | +100 h (FPGA work) |

Note that rev B's 520 h headline differs from its own step list, which sums to 630 h.

My view: S3 is the natural successor. It keeps all four jobs in every licence-free transmit band and moves to a rev-A-style 5.8 GHz patch tile. S2b is the bold alternative if 5.8 GHz can go. The cost of narrowing is losing general-purpose tuning: cellular uplinks, n78/CBRS passive radar, and 5.1-5.35 GHz.

### Findings

- **EU licence-free power, sub-GHz** (high): ERC/REC 70-03 sets these limits: - 433.05-434.79 MHz: 10 mW e.r.p. at 10 % duty or less, or 1 mW with no restriction; 434.04-434.79 MHz allows 10 mW in channels of 25 kHz or less. - 863-870 MHz: 25 mW e.r.p. at 0.1-1 % duty, or with listen-before-talk plus adaptive frequency agility. - 869.40-869.65 MHz: 500 mW e.r.p. at 10 % duty or less (250 kHz wide: Doppler-only radar). - 915-919.4 MHz: 25 mW e.r.p. at 1 % duty or less, only where implemented nationally. Source: https://docdb.cept.org/download/5013 (ERC/REC 70-03 PDF, Annex 1 g1-g3, h1.x, h3; text checked locally)
- **EU licence-free power, 2.4 and 5.8 GHz** (high): ERC/REC 70-03 allows: - 2400-2483.5 MHz: 10 mW e.i.r.p. (Annex 1 i); 25 mW e.i.r.p. for radiodetermination, i.e. radar (Annex 6 c); 100 mW e.i.r.p. for wideband data with LBT/DAA and at most 10 mW/MHz (Annex 3 b). - 5725-5875 MHz: 25 mW e.i.r.p. (Annex 1 j). Tank level probing radar at 4.5-7 GHz applies only inside closed tanks. Source: https://docdb.cept.org/download/5013
- **US licence-free power (15.249 / 15.245 / 15.231)** (high): 15.249 limits 902-928, 2400-2483.5 and 5725-5875 MHz to 50 mV/m at 3 m, which is -1.25 dBm EIRP average. 15.245 field-disturbance sensors (perimeter protection excluded) get 500 mV/m, or +18.75 dBm EIRP average, in 902-928, 2435-2465 and 5785-5815 MHz. Their harmonics below 17.7 GHz must meet 15.209, so the 11.6 GHz harmonic of a 5.8 GHz radar needs about 60 dB of suppression. 15.231 at 433.92 MHz allows only control bursts, released within 5 s, at about -14.4 dBm EIRP average (10,997 uV/m); at 315 MHz it is about -19.6 dBm. Source: https://www.law.cornell.edu/cfr/text/47/15.245 ; https://www.law.cornell.edu/cfr/text/47/15.249 ; https://www.law.cornell.edu/cfr/text/47/15.231 ; conversions in scratch/band_plan/reg_metrics.json
- **US 15.247 for beam steering (missing from REDESIGN)** (medium): Under 15.247, digitally modulated systems (6 dB bandwidth of at least 500 kHz) in 902-928, 2400-2483.5 and 5725-5850 MHz may use 1 W conducted. - Antenna gain above 6 dBi reduces that power dB for dB, which caps EIRP at 36 dBm. Fixed point-to-point links at 5725-5850 MHz need no reduction. - At 2.4 GHz, transmitters with multiple directional beams count array gain as 10 log N plus the element gain and lose 1 dB per 3 dB above 6 dBi. For 16 elements of 6 dBi that is 26 dBm conducted, 44 dBm EIRP. - Spectral density is capped at 8 dBm per 3 kHz. REDESIGN section 7 assumes 15.249 and backs off 9.3/7.5 dB. Under 15.247, PN/BPSK-coded steering would need no back-off: at the 36 dBm EIRP cap, each of 16 coherent 6 dBi elements may carry +5.9 dBm. Source: https://www.law.cornell.edu/cfr/text/47/15.247 ; per-element powers in scratch/band_plan/band_hw.json
- **US amateur bands for steering and radar** (high): 47 CFR 97.305 authorises pulse emissions on 33 cm (902-928), 13 cm and 5 cm (5650-5925). It does not on 70 cm or 23 cm, and 9 cm is no longer in the table. 97.313 caps power at 1.5 kW PEP, with 50 W PEP on 70 cm in listed areas and on 33 cm within 241 km of White Sands Missile Range. Source: https://www.law.cornell.edu/cfr/text/47/97.305 ; https://www.law.cornell.edu/cfr/text/47/97.313 ; https://www.arrl.org/3-ghz-band
- **EU amateur bands (Germany as example)** (medium): German class A licences (AFuV Anlage 1) allow: - 70 cm, 430-440 MHz: 750 W PEP, primary. - 23 cm: 750 W PEP, secondary; 1247-1263 MHz limited to 3.05 W ERP. - 13 cm 2320-2450, 9 cm 3400-3475 and 6 cm 5650-5850 MHz: 75 W PEP, secondary. ECC/DEC/(25)01, in force since 27 June 2025, adds 23 cm EIRP caps to protect Galileo (e.g. 20 mW EIRP at 1262-1270 MHz for amateur satellite). The UK will not impose the amateur conditions. Source: https://gesetze.legal/bund/afuv/anlage_1 ; https://amsat-uk.org/tag/rnss/ (secondary report of ECC/DEC/(25)01)
- **Restricted-band share per set** (high): 47 CFR 15.205 restricted bands cover 47.5 % of 0.4-6 GHz, where radiated calibration is refused. By set: S1 (2.4+5.8) and S3 (433+868/915+2.4+5.8) 0 %; S2 5.0 % (only the 400-410 MHz radiosonde slice, receive-only); S4 29.1 %. Source: https://www.law.cornell.edu/cfr/text/47/15.205 ; computed in scratch/band_plan/catalogue.json
- **Drone links (DF targets)** (high): DJI O4 on the Mini 4 Pro uses 2.400-2.4835 GHz and 5.725-5.850 GHz, plus 5.170-5.250 GHz where permitted. EIRP is under 33 dBm (FCC) at both bands; under CE rules it is under 20 dBm at 2.4 GHz and under 14 dBm at 5.8 GHz. Analog FPV video uses 5645-5945 MHz; ELRS uses 900 MHz and 2.4 GHz; Remote ID uses Bluetooth 4/5 Long Range and Wi-Fi beacon/NAN. Source: https://www.dji.com/mini-4-pro/specs ; https://oscarliang.com/fpv-channels/ ; https://www.expresslrs.org/hardware/hardware-selection/ ; https://www.faa.gov/sites/faa.gov/files/uas/getting_started/remote_id/industry/Remote_ID_Broadcast_Module.pdf
- **LoRa / ISM device targets** (high): LoRaWAN RP002 channel plans: - EU868: default channels 868.1, 868.3 and 868.5 MHz, under 1 % duty, +16 dBm EIRP default. - US915: uplink 64 x 125 kHz at 902.3-914.9 MHz plus 8 x 500 kHz at 903.0-914.2 MHz; downlink 8 channels at 923.3-927.5 MHz. - EU433: 433.175-434.665 MHz, under 12 dBm EIRP, under 10 % duty. Tyre-pressure sensors are sold at both 315 and 433 MHz. Bluetooth 5.1 adds direction-finding tags (AoA/AoD). Source: https://lora-alliance.org/wp-content/uploads/2021/05/RP002-1.0.3-FINAL-1.pdf ; https://www.autoserviceworld.com/carsmagazine/schrader-adds-new-433-mhz-coverage-to-ez-sensor-line ; https://www.bluetooth.com/learn-about-bluetooth/feature-enhancements/direction-finding/
- **315 MHz is outside rev B** (medium): US key fobs and tyre sensors at 315 MHz sit below rev B's 0.4 GHz edge, below its 330 MHz high-pass on the external inputs, and below the AT86RF215's lowest band (389.5 MHz). Adding them needs a lower high-pass and an AD936x, HackRF or RTL-class receiver. Source: REDESIGN.md section 3 ; https://ww1.microchip.com/downloads/en/DeviceDoc/Atmel-42415-WIRELESS-AT86RF215_Datasheet.pdf ; https://www.autoserviceworld.com/carsmagazine/schrader-adds-new-433-mhz-coverage-to-ez-sensor-line
- **Known-position emitters for calibration** (medium): Several transmitters announce or have known positions: - ADS-B (1090 MHz) and UAT (978 MHz) aircraft broadcast their own GPS positions. - Radiosondes (400.15-406 MHz) transmit their GPS positions. - Cospas-Sarsat beacons at 406.025-406.076 MHz send 5 W, 0.5 s bursts every 50 s. - Inmarsat (1525-1559 MHz) is geostationary. - GNSS L1/E1 1575.42, L2C 1227.60 and L5/E5a 1176.45 MHz arrive from directions known from ephemeris. My inference: these let array16 measure its manifold in place in bands where 15.205 forbids radiated calibration. Source: https://www.faa.gov/air_traffic/technology/equipadsb/installation ; https://github.com/projecthorus/radiosonde_auto_rx ; https://en.wikipedia.org/wiki/Emergency_position-indicating_radiobeacon ; https://archive.gps.gov/systems/gps/modernization/civilsignals/ ; https://gssc.esa.int/navipedia/index.php/Galileo_Signal_Plan
- **Passive-radar illuminators** (medium): The EU keeps 470-694 MHz for DVB-T broadcasting; US TV was repacked below 608 MHz. LTE/5G downlinks such as B20 791-821, B5 869-894, B8 925-960, B3 1805-1880, B7 2620-2690 MHz and n78 3300-3800 MHz also serve. FM is out of range. KrakenRF uses HDTV as its passive-radar illuminator. Source: https://eur-lex.europa.eu/eli/dec/2017/899/oj/eng ; https://www.everycrsreport.com/reports/IF10751.html ; https://en.wikipedia.org/wiki/LTE_frequency_bands ; https://github.com/krakenrf/krakensdr_docs/wiki/08.-Passive-Radar/191d28c9f2a665558d40f99d526fbfe0e13d2361
- **Band sets defined** (high): Four band sets are defined (occupied bandwidth; widest piece as a fraction of its centre frequency): - S1, 2.4 + 5.8 GHz (2400-2483.5 and 5650-5925 MHz): 358 MHz; widest piece 4.8 %. - S2, 400-450 + 863-928 + 2400-2483.5 MHz: 198 MHz. - S3, 420-450 + 863-928 + 2.4 + 5.8 GHz: 454 MHz; widest piece 7.3 %. - S4, S3 plus 400-406, 470-694, 975-1095, 1164-1300, 1400-1427 and 1525-1590 MHz: 1,046 MHz. Rev B covers 5,600 MHz. Source: scratch/band_plan/catalogue.json ; chart band_sets.png
- **Narrowing does not help the whip layout** (high): The sibling geometry study finds rev B's two rings (16 whips at 225/112 mm) good from 0.30 to 2.72 GHz. Their 95th-percentile direction errors are 1.32, 0.66, 0.64 and 0.24° at 433, 868, 915 and 2440 MHz. A layout optimised for 433/915/2440 only is barely better (1.24° against 1.33°). Whips at 5.8 GHz fail for sources 0-30° up (32 % gross errors), so 5.8 GHz needs the panel. Source: experiments/band_plan/tables.txt (sibling simulation)
- **Band-specific panel aperture at 2.4 GHz** (medium): Rev B's 25 mm panel is 100 mm wide, 0.8 wavelength at 2.44 GHz: 9.2 dBi at most, 6.4 dBi realised, 62° beam. A 2.4 GHz-only 4x4 at half-wavelength pitch (60.4 mm, 241 mm aperture) reaches 16.9 dBi (about 16.2 dBi realised, estimate) with a 26° beam. G/T rises from -20.2 to about -10.3 dB/K. The sibling simulation gives a 95th-percentile direction error of 0.62° against 1.49°. A dual-band element on the 25 mm lattice cannot recover this, because the aperture stays 0.8 wavelength. Source: computed in scratch/band_plan/band_hw.json ; REDESIGN.md sections 4-5 ; experiments/band_plan/tables.txt part 3
- **Noise figure and power at 5.8 GHz** (medium): Rev B's panel noise figure at 5.8 GHz is 3.4 dB (5.0 dB in the robust state). A narrowband chain with two SKY65404-31 per element (NF 1.0 dB, 13 dB gain, 10 mA, 4.9-5.9 GHz) gives about 2.2 dB with the conducted cal arm and 1.8 dB without it. The same chain with rev B's PMA3-83LNW+ gives 3.2 dB, so the gain comes from the band LNA, not the filter. G/T at 5.8 GHz improves from -12.5 to about -10.2 to -10.7 dB/K. Panel LNA power falls from about 5.8 W to about 1 W (estimate). Source: Skyworks SKY65404-31 datasheet (local copy in scratchpad ds_sky65404) ; cascade in scratch/band_plan/band_hw.json
- **Noise figure below 2.5 GHz does not improve** (medium): A PSA4-5043+ LNA (NF 0.65/0.75/0.98 dB, gain 22.1/18.4/13.3 dB at 0.5/1/2 GHz, 58 mA at 5 V) ahead of a SAW or LTCC filter and an AT86RF215 (NF 4.3-4.5 dB) gives a system noise figure of 1.05 dB at 433 MHz, 1.26 dB at 915 MHz and 1.5-1.9 dB at 2.44 GHz. Rev B's external path already gives 1.1/1.3/1.4 dB, so narrowing gains almost nothing here. Source: https://www.minicircuits.com/pdfs/PSA4-5043+.pdf ; AT86RF215 datasheet https://ww1.microchip.com/downloads/en/DeviceDoc/Atmel-42415-WIRELESS-AT86RF215_Datasheet.pdf ; computed in band_hw.json
- **Fixed filters instead of the B1-B5 bank** (high): Mini-Circuits datasheets (checked locally) give: - BFCN-2450+ LTCC: 2400-2550 MHz, 2.0 dB (2.9 max); 30 dB below 2.1 GHz, 20 dB at 3.4-12 GHz. - BFCN-5750+: 5650-5850 MHz, 1.85 dB; 25 dB below 4.3 GHz and at 8-9.5 GHz. - BFCN-1575+: 1530-1620 MHz, 2.0 dB. - ZVBP-2450-S+ cavity: 2400-2500 MHz, 0.7 dB; at least 70 dB at DC-2.12 GHz and 2.78-4 GHz. - ZVBP-5800-S+ cavity: 5725-5875 MHz, 0.8 dB; at least 50 dB at DC-5.2 GHz and 6.4-14 GHz, which covers the 11.6 GHz harmonic REDESIGN wants a cavity filter for. Rev B's B5 path costs 5.6-6.2 dB at 5.8 GHz. Prices could not be retrieved. Source: https://www.minicircuits.com/pdfs/BFCN-2450+.pdf ; https://www.minicircuits.com/pdfs/BFCN-5750+.pdf ; https://www.minicircuits.com/pdfs/BFCN-1575+.pdf ; https://www.minicircuits.com/pdfs/ZVBP-2450-S+.pdf ; https://www.minicircuits.com/pdfs/ZVBP-5800-S+.pdf
- **Switched lines re-chosen for narrow bands** (high): Rev B's four delay lines (38.5/52.5/117.5/139 ps) give 11.3° RMS and 36° worst-case phase error over 2.3-6 GHz. In the ISM bands alone they give 9.8°/11.4° RMS and 28-30° worst. - Re-chosen for 2.4+5.8 GHz (41.1/82.1/105.9/159.2 ps): 8.1° RMS and 18° worst at both bands. A state pair within about 7-10° of 180° remains for direction-finding toggles. - Single-band sets approach an ideal 4-bit shifter (6.5° RMS / 11.25° worst): 2.4 GHz 26.1/52.2/104.6/206.4 ps gives 6.6° / 11.7°; 5.8 GHz gives 7.1° / 13.1°. - Three bits give 13-15° RMS. The cost is the same $2.6 per element; only the line lengths change. These are lossless-line bounds. Source: scratch/band_plan/phase_states.json (script scripts/phase_states.py)
- **True-time delay is unnecessary within one band** (high): Steering by phase alone to 45° moves the beam by at most 0.98° at the 2.4 GHz band edges, 0.74° at 5725-5875 MHz, 1.36° at 5650-5925 MHz, 0.81° at 902-928 MHz and 2.1° over 863-928 MHz. Beams are 26-60° wide. Over rev B's 2.3-6 GHz the same steering would move the beam about 25°, which is why rev B needs true-time delay. Source: computed in scratch/band_plan/band_hw.json
- **Digital phase shifters do not pay** (medium): Analog Devices' sub-6 GHz digital phase shifters each cover a narrow band: - HMC936A 1.2-1.4 GHz - HMC647A 2.5-3.1 GHz - HMC648A 2.9-3.9 GHz - HMC649A 3-6 GHz: 6-bit, 8 dB loss, about $169 for one unit - HMC1133 4.8-6.2 GHz None covers 2.4 GHz (pSemi PE44820 is 1.7-2.2 GHz), and the vector modulators are narrowband (HMC631 1.8-2.7 GHz). Re-chosen switched lines at about $2.6 per element remain the practical phase control. Source: https://www.analog.com/en/product-category/phase-shifters.html (wideband_study.md source 12; prices from this session's earlier distributor check, not re-verified)
- **AT86RF215 as a band-specific radio** (high): The Microchip AT86RF215 covers 389.5-510, 779-1020 and 2400-2483.5 MHz: - 13-bit I/Q at up to 4 MHz sampling over LVDS (one TX link, two RX links). - Receive analog band-pass of 160-2000 kHz. - NF 4.3/4.5/4.3 dB. - Up to +14.5 dBm at 900 MHz. - 26 MHz TCXO input, 28 mA receiving, sub-GHz and 2.4 GHz at the same time. The open CaribouLite board pairs it with an iCE40LP1K, a 26 MHz TCXO and an RFFC5072 mixer. The unit price could not be checked; about $6-10 is my estimate. Source: https://ww1.microchip.com/downloads/en/DeviceDoc/Atmel-42415-WIRELESS-AT86RF215_Datasheet.pdf ; https://raw.githubusercontent.com/cariboulabs/cariboulite/main/hardware/README.md ; https://raw.githubusercontent.com/cariboulabs/cariboulite/main/hardware/rev2/assembly/full_bom/cariboulite_bom.xlsx
- **Other cheaper radio ICs** (high): Alternatives to the AD9361, each with limits: - Semtech SX1255 (400-510 MHz) and SX1257 (860-1000 MHz) are I/Q transceivers; Semtech's pages list decimated I/Q over I2S and show no price. - The AD9363 in ADALM-PLUTO covers 325 MHz-3.8 GHz at 0.2-20 MHz bandwidth, 12-bit. It is advertised as one RX/one TX and streams 4 MSPS over USB 2. - RTL-SDR Blog V3 (R860/R820T2, 24-1766 MHz) costs $34.95-39.95; the V4 is end-of-life. - KrakenSDR has 5 coherent R820T2 channels, receive-only, 2.56 MHz bandwidth, for $749. Source: https://www.semtech.com/products/wireless-rf/lora-core/sx1255 ; https://www.semtech.com/products/wireless-rf/lora-core/sx1257 ; https://wiki.analog.com/university/tools/pluto/devs/specs ; https://www.rtl-sdr.com/buy-rtl-sdr-dvb-t-dongles/ ; https://www.krakenrf.com/product-page/krakensdr
- **5.8 GHz is what forces AD9361 radios** (high): Of the radios considered, only the AD9361 and HackRF cover every set; the bladeRF 2.0 micro xA4 costs $540 for 2 channels (47 MHz-6 GHz, 56 MHz). The AT86RF215 and AD9363 cover 100 % of S2 but only 23 % of S1 and 39 % of S3, missing exactly the 5.8 GHz part. QuadRF (MAX2850/MAX2851; 4 RX and 4 TX coherent; 4.9-6.0 GHz; 40 MHz; 1 W per antenna; $499 kit) covers only the 5.8 GHz part. Source: https://www.nuand.com/bladerf-2-0-micro/ ; https://github.com/open-space-sdr/main ; https://scalerf.com/ ; coverage computed in catalogue.json
- **Radar range resolution per band and radio** (high): Range bins (c/2B) for each allocation: - EU 433 MHz (1.74 MHz): 86 m - EU 869.40-869.65 MHz: 600 m - US 15.245 at 902-928 (26 MHz) and 2435-2465 / 5785-5815 MHz (30 MHz): 5.8 m and 5 m - EU 2.4 GHz (83.5 MHz): 1.8 m - EU 5.8 GHz (150 MHz): 1.0 m - US amateur 5650-5925 MHz: 0.5 m The radios set floors: AD9361 at 56 MHz gives 2.7 m, AD9363 at 20 MHz 7.5 m, AT86RF215 at 2 MHz 75 m. Source: computed in scratch/band_plan/band_hw.json
- **Per-element power when steering 16 elements** (high): At each region's EIRP limit, 16 coherent 6 dBi elements may carry per element: - EU 25 mW (5.8 GHz, or 2.4 GHz radar): -16.1 dBm - EU 10 mW at 2.4 GHz: -20.1 dBm - US 15.245: -11.3 dBm - US 15.249: -31.3 dBm - US 15.247: +5.9 dBm Rev B's full drive already gives +8.1 / +6.2 dBm EIRP at 2.44 / 5.8 GHz. Regulation, not hardware, limits licence-free steering. Source: computed in band_hw.json and reg_metrics.json from https://docdb.cept.org/download/5013 and https://www.law.cornell.edu/cfr/text/47/15.247
- **S1, 2.4 + 5.8 GHz: hardware** (medium): S1 keeps rev B's core and four channels but drops the 16 external whip front ends. A two-filter bank (BFCN-2450+ / BFCN-5750+) and cavity filters on RADAR OUT replace B1-B5 and the radar bank. Two variants: - S1a keeps the Vivaldi tile with re-chosen lines. - S1b uses two patch tiles: 5.8 GHz at 25.5 mm (RO4003C/PTFE patches, 76-84 % efficient, 16.1-16.6 dBi in antenna_study) and 2.4 GHz at 60.4 mm. The 5.8 GHz half still needs AD9361 radios. S1 loses everything below 2.4 GHz (LoRa/ELRS/tyre-sensor direction finding, ADS-B/GNSS/sonde calibration, TV passive radar), n78 and 5.1-5.35 GHz. Source: experiments/antenna_study.md ; Mini-Circuits datasheets above ; cost_effort.json
- **S2, sub-GHz + 2.4 GHz: hardware** (medium): S2 needs only whips; one fixed ring layout covers all its bands. Fixed SAW filters for 433 and 863-928 MHz plus BFCN-2450+ replace B1-B5, and no tile or panel phase lines are needed. Two radio routes: - S2a: two 2x2 AD9363 boards. - S2b: one AT86RF215 per whip, i.e. 16 coherent receivers in both bands at once and 16 transmitters, giving full digital beamforming and nulling, behind an ECP5 FPGA and a USB 3 bridge. S2 loses 5.8 GHz and high-gain 2.4 GHz beams, unless a 60 mm 2.4 GHz tile is added. With the AT86RF215 it is also limited to 2 MHz bandwidth, so Wi-Fi and DJI signals are seen as 2 MHz slices; that is still enough for direction finding. Source: AT86RF215 datasheet ; experiments/band_plan/tables.txt ; scratch/band_plan/cost_effort.json
- **S3 and S4: hardware** (medium): S3 (all licence-free transmit bands) keeps rev B's 16 whip ports and swaps B1-B5 for five narrow fixed filters (433, 868 and 915 MHz SAW, BFCN-2450+, BFCN-5750+) at about the same $24 per channel. The Vivaldi becomes a rev-A-style 5.8 GHz patch tile, and 2.4 GHz is served by the whips. S4 keeps a five-position bank re-tuned to groups under 3:1 (400-694, 863-1095, 1164-1590 MHz, 2.4, 5.8 GHz), which keeps the AD9361's 3xLO response in a stopband. Both still need AD9361 radios. S4 can add a KrakenSDR for its receive-only bands below 1.77 GHz, which are 67 % of its spectrum. Source: REDESIGN.md sections 3 and 5 ; catalogue.json ; https://www.minicircuits.com/pdfs/BFCN-1575+.pdf
- **Cost and effort against rev B** (low): All figures are estimates against rev B's $2,630 all-in (boards $1,300 + calibration kit $250 + 2 bladeRFs $1,080) and 520 h: - S1a: about -$335, -70 h - S1b: about -$123, -100 h, with +10 dB G/T at 2.44 GHz and about 5 W less - S2a: about -$1,315 (radio price unverified), -200 h - S2b: about -$1,580, but +100 h of FPGA and radio-board work - S3: about -$228, -100 h, 5 W less - S4: about -$228, -90 h Source: scratch/band_plan/cost_effort.json (itemised; anchors from REDESIGN.md section 10)
- **Rev B effort figure is inconsistent** (high): REDESIGN.md gives 520 hours in its headline (section 1), but its build-order steps (section 12: 40+30+60+180+180+100+40 h) sum to 630 h. The 100 h of host software is described as running 'across steps 1-5'. Source: REDESIGN.md
- **Radiated calibration option in ISM-only sets** (low): My inference: in S1 and S3 every band can radiate at 15.249 levels, so a probe or centre antenna could replace the conducted calibration arm on each element. That saves 16 BGS12WN6 and about 0.46 dB of noise figure at 5.8 GHz. Rev B chose conducted arms partly because 51 % of its panel band is restricted. Source: REDESIGN.md section 6 ; computed NF variants (band_hw.py)
- **Lattice follows the band set** (medium): Half-wavelength pitch at each set's top band is 60.4 mm for 2.4 GHz, 25.3-25.5 mm for 5.8 GHz, 162 mm for 928 MHz and 345 mm for 434 MHz. One 2-D lattice serving 2.4 and 5.8 GHz is set by 5.8 GHz and leaves 2.4 GHz with too little aperture. My view: two sub-lattices, interleaved or stacked in 3-D, answer the 2.4+5.8 aperture problem. This is left to the parallel lattice and 3-D studies. Source: estimate: c/(2 f_top) ; band_hw.json panels

### Files

- `scratch/band_plan/catalogue.json`
- `scratch/band_plan/band_sets.png`
- `scratch/band_plan/phase_states.json`
- `scratch/band_plan/band_hw.json`
- `scratch/band_plan/cost_effort.json`
- `scratch/band_plan/reg_metrics.json`
- `scratch/band_plan/scripts/catalogue2.py`
- `scratch/band_plan/scripts/phase_states.py`
- `scratch/band_plan/scripts/band_hw.py`
- `scratch/band_plan/scripts/cost_effort.py`
- `scratch/band_plan/scripts/reg_metrics.py`
- `scratch/band_plan/scripts/bands.py`
- `scratch/band_plan/dl_at86rf215/at86rf215.txt`
- `scratch/band_plan/dl_mcl`
- `scratch/band_plan/dl_erc7003/rec7003.txt`
- `scratch/band_plan/dl_lora/rp002.txt`
- `experiments/band_plan/tables.txt`

### Open questions

- Prices not verified. Microchip, DigiKey, Mouser, Arrow, LCSC, the Mini-Circuits store and analog.com refused or timed out, and the shared web-search budget was used up before this run. Affected parts: AT86RF215 (estimated $6-10), SX1255/SX1257, PSA4-5043+, SKY65404-31 (about $1.06 from earlier session notes), sub-GHz SAW filters, Mini-Circuits LTCC and cavity filters, HMC649A (about $169), and Pluto-class 2x2 AD9363 boards (sibling survey $180-235, no URL).
- Specific SAW/BAW part numbers for 433.92 MHz and 863-928 MHz, and BAW options at 2.4/5.8 GHz, were not selected or verified.
- Do several AT86RF215s sharing one 26 MHz TCXO need re-calibration after every retune (independent PLLs), and how are their samples aligned when the chip has no sync input? This needs a bench test before betting on S2b.
- Can a home-built device (15.23) do PN/BPSK beam-steering experiments under 15.247 (1 W, 36 dBm EIRP) instead of 15.249? Does FMCW qualify under 15.245 or 15.247? In the EU, does Annex 1 j (5725-5875 MHz, 25 mW) allow radar use?
- Has any FCC proceeding changed 902-928 MHz for Part 15 or amateur use (e.g. NextNav's lower-900 MHz proposal)? Not checked.
- The owner's region and licence decide the sub-GHz bands: EU 433/868 against US 915. A Technician licence opens pulsed 33/13/5 cm at high power.
- Coupling, scan match and manifold of a 60.4 mm 2.4 GHz suspended-patch tile and a 25.5 mm 5.8 GHz patch tile need openEMS runs; the gains quoted here use an assumed -0.7 dB efficiency.
- Is 2 MHz instantaneous bandwidth (AT86RF215, 75 m radar bins) acceptable for the owner's radar ambitions, or is S2a (AD9363, 20 MHz) the better sub-GHz + 2.4 route?
- REDESIGN.md's 520 h headline does not match its 630 h step list; which is intended?
