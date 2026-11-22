<!-- 2026-11-21. Who builds charged-particle tracks, traps, rings and accelerators on boards or chips, what a
generic-track RF accelerator on boards could be used for, and what the fields that build such tracks teach array16.
Four sweeps (prior art, applications, lessons, novelty), each checked by an adversarial verifier (87 items
confirmed or corrected; uncertain ones are marked or left out; refuted ones are left out), then this synthesis. -->

# Prior art: particle tracks on boards and chips

The owner's question, 2026-11-21: "it is a crazy thought: building a generic-track particle accelerator on a pcb,
what other interesting applications could come from that? is this even something people work on?"

Yes, in pieces: printed tracks that carry ions round turns, switches and stacked levels, electrons guided and trapped on
microwave boards at 1 to 1.7 GHz, and RF gaps that accelerate beams on FR-4 all exist, and the ion tracks are sold as an
instrument. What has apparently not been done is the combination array16 is after: a closed RF guide that holds
electrons and speeds them up round a drawn ring or knot, driven from one radio.

## Who builds what, closest first

| Field | What they build | Scale | Loop? | Accelerates? | PCB? |
|---|---|---|---|---|---|
| Ion mobility, SLIM [1–5] | tracks between facing boards: turns, switches, escalators; KiCad layout from a path string; sold since 2021 | 1,094 m in 81 passes; 88 m over 8 levels | yes, in 3-D | no; conveys in 2–4 Torr | yes |
| Electron guides, Hommelhoff group [6, 7] | microwave and static alternating guides, splitters | 1–5 eV at 970 MHz; up to 4 keV between boards 1 mm apart | no | no | yes |
| GHz electron traps [8, 9] | Paul traps in printed resonators | 1.6–1.73 GHz, Q 35–41 | no | no | yes |
| RF quadrupole rings [10, 11] | metal-electrode rings and racetracks | 115 mm ring; Mg⁺ crystalline at 1 eV over 3,000+ turns | yes | yes, by laser (PALLAS) | no |
| Trapped-ion computing [12, 13] | racetracks and rings on chips | 1.25–6.6 mm | yes | no; DC conveyors | no |
| Board linacs [14, 15] | RF gaps on FR-4 wafers; a GHz coplanar resonator | Ar⁺ 7.25 keV per gap (kV drive); electrons 2 → 2.3–3.2 keV in 8 cm | no | yes | yes |
| Laser accelerators [16] | nanostructured chips | 28 → 41 keV in 0.5 mm | no | yes | no |
| Microparticle traps [17, 18] | planar traps with a cross junction; an open 60 Hz trap | 0.4–32 µm grains | no | no; shuttles | yes |

## Open ground: what the searches did not find

One sweep's 57 queries (2026-11-21) found none of these; the search budget ran out before a second pass.

- **A closed guide for electrons**, on one board or several (8 queries on electron storage rings, ring guides and ring
  traps). Ions in gas have closed board tracks [3]; electrons have open guides and linear traps [6–9]. The multi-pass
  resonator in [19] is the lead.
- **An electrode guide along a knot** (8 queries on knotted traps and guides, and stacked, flexible or 3-D-printed
  guides). The nearest are a magnetic knotatron design [20] and a theory of charged particles following knotted nodal
  lines of light [21].
- **RF acceleration round an RF-guided ring at board scale** (8 queries). PALLAS accelerated with laser light [11], and
  board linacs run straight [14, 15]. PALLAS's bunched beams [22] are the lead.
- **The confining RF itself circulating as e^(i(mφ − ωt))** (6 queries). The nearest is a phased RF carpet, whose one
  wave levitates and sweeps ions [23].
- **A software-defined radio driving any guide or trap** (5 queries: HackRF, USRP, Red Pitaya, GNU Radio). Published
  drives are amplifier-fed resonators.
- **An open-source GHz electron guide, or a hobby closed RF guide** (4 queries). The open designs found run at 50–60 Hz
  with spores, such as [18].
- Closed ring traps for charged microparticles were not searched.

## Applications, ranked by what the printed track adds and how near the plan they sit

1. **The ring hears its own beam (N4).** Every electrode is a pickup: 10⁶ electrons at 185 eV on the 932 mm loop
   circulate at 8.66 MHz and give about −78 dBm on a 30 % pickup. A GHz trap heard about 1,260 electrons by their image
   current and cooled them through it in 1.74 s [24]. The HackRF Pro is half duplex [25], so reading pickups while it
   drives the ring takes a receive path of its own.
2. **The beam transmits in the radio's band.** A beam over a printed comb of period d radiates Smith–Purcell waves [26]
   at 0.593·√E[eV]/d[mm] GHz (at 90°): 100 eV over 1 mm gives 5.93 GHz. It takes the beam within 80–200 µm of the comb,
   and a 5.9 GHz bunching gap and drift for coherent output.
3. **The guide is an antenna (V4).** The 277 mm ring carries m = 1 at 1.08 GHz in air, its radiation Q switched from
   about 3 × 10⁸ to 20 by PIN-diode stubs (Hypothesis NFD-H1).
4. **A recirculating electron resonator.** Quantum-electron-microscope designs need each electron to pass the sample
   many times, with a microwave pseudopotential among their couplers [19]; one transverse mode (0.88 µeV spacing at 212
   MHz) takes a matched tip source.
5. **Ions at MHz on the same boards.** In 2–5 mbar, fifty laps of the 277 mm track equal SLIM's 13 m, where resolving
   power grows as √path [2]; near 3 × 10⁻⁷ mbar (estimated) a 100 u ion at 10 eV runs 277 m in 63 ms.
6. **Charged grains on a knot (N7).** Grains closer than 0.65 mm buckle into a zigzag, so about 1,250 fit in single file
   on the 814 mm trefoil; each grain's charge is 125 times a printed detector's single-pass floor [27], so the board
   counts them lap by lap. In vacuum the bend rule lets a 450 V guide 0.5 mm from the axis, at q = 0.5, hold 1.24 kV of
   acceleration round the 44 mm ring: 205 m/s for a 200 nm silica grain with 1,000 e, at 56 kHz.
7. **Rotation sensing and qubits.** Counter-propagating electrons round the 61 cm² ring differ in phase by 211 rad per
   rad/s per lap; the ring is also an Aharonov–Bohm loop, so reading Earth's rate takes the flux through it held to 0.8
   fT [28]. Designs for electron qubits at 4 K project 99.99 % two-qubit fidelity [29].

## What the existing fields teach, as what it takes

- **Every joint, feed, port and gap is a junction.** Shaping cut a Y junction's barrier from over 1 eV to under 2 meV
  [30] and a guide entrance's field from 120 to 6 V/m [6]. An RF null cannot pass straight through a crossing [31];
  rings and knots make none, so each joint is shaped from its solved field before boards are drawn.
- **Vacuum keeps every kick.** Gas damps an ion within a microsecond, yet a mismatched escalator entrance still cost
  SLIM 25 % [32]. In vacuum each board-to-board hand-off changes slowly over several secular periods, or the pickups
  damp it.
- **Stray charge.** At q = 0.3 the well tilts empty at 1.3 kV/m on four rods and about 0.3 kV/m on facing boards, and a
  million lost electrons on 1 mm² of laminate 0.7 mm away make about 2 kV/m. That takes no laminate in sight of the
  null, DC returns, a metal dump and trim electrodes on a measured map; nulling stray fields made a PCB ion trap live
  about a hundred times longer [33].
- **The drive does the heating.** The PCB trap's electrons were lost to resonances at the drive's 6th and 7th
  subharmonics, not to gas [8]; junction ions heated on drive noise and DAC aliases [34]. That takes f_sec clear of Ω/n
  for small n, and nothing in the drive, aliases included, at f_sec, 2f_sec or Ω ± f_sec.
- **Tune.** The 932 mm acceleration crosses 11 integer tunes at q = 0.3 and 25 at q = 0.6, so it crosses them fast or
  with damping; PALLAS ran near a tune of 50 [11].
- **Multipactor.** Copper's lowest two-surface boundary, 27.4 V near 0.70 GHz·mm [35], is where the 1 GHz rods sit: 6.5
  dB of margin at 13 V, 0.8 dB at 25 V. The margin follows well depth and q, not band (derived); facing boards reach it
  at 2.5–3 eV wells, so the tracker takes a secondary-emission model.
- **Resonator and gap.** Built board resonators reach Q 35–41 (one fell from a simulated 124) and about 1 kΩ [8, 9]. 100
  kΩ, 450 V from 1 W, takes Q ≈ 3,000 on a 50 Ω half-wave line, or 300 at a tenth of its capacitance. q goes as 1/gap²,
  so ±10 µm on a 1 mm gap is ±2 % in q.
- **Loading and vacuum.** A closed guide returns each electron to its entrance every lap (165 ns at 8 eV), so electrons
  are made inside, as by 405 nm photoemission from gold [24], or the entrance switches within a lap. The 126 m
  acceleration keeps 73 % unscattered at 10⁻⁶ mbar, against 4.5 % at 10⁻⁵.

## Sources

Every DOI and arXiv id was re-checked for this note; 5 and 15 are vendor and conference reports.

1. Deng 2016, Anal. Chem. 88, 8957: [10.1021/acs.analchem.6b01915](https://doi.org/10.1021/acs.analchem.6b01915)
2. Deng 2017, Anal. Chem. 89, 4628: [10.1021/acs.analchem.7b00185](https://doi.org/10.1021/acs.analchem.7b00185)
3. Huntley 2024, Anal. Chem.: [10.1021/acs.analchem.3c05594](https://doi.org/10.1021/acs.analchem.3c05594)
4. Greer 2023, JASMS 34, 1715: [10.1021/jasms.3c00163](https://doi.org/10.1021/jasms.3c00163);
   [SLIM-Pickins](https://github.com/bhclowers/SLIM-Pickins)
5. MOBILion MOBIE launch, 2021:
   [SelectScience](https://www.selectscience.net/article/mobilion-systems-launches-mobie-the-first-slim-based-high-resolution-ion-mobility-product-to-accelerate-biotherapeutic-drug-development-and-multiomic-biomarker-discovery)
6. Hoffrogge 2011, PRL 106, 193001: [10.1103/PhysRevLett.106.193001](https://doi.org/10.1103/PhysRevLett.106.193001)
7. Zimmermann 2021, Nat. Commun. 12, 390: [10.1038/s41467-020-20592-4](https://doi.org/10.1038/s41467-020-20592-4)
8. Matthiesen 2021, PRX 11, 011019: [10.1103/PhysRevX.11.011019](https://doi.org/10.1103/PhysRevX.11.011019)
9. Luo 2026, preprint: [arXiv:2607.18401](https://arxiv.org/abs/2607.18401)
10. Church 1969, J. Appl. Phys. 40, 3127: [10.1063/1.1658153](https://doi.org/10.1063/1.1658153)
11. Schätz 2001, Nature 412, 717: [10.1038/35089045](https://doi.org/10.1038/35089045); Schramm 2003: [PAC
    TOAA004](https://proceedings.jacow.org/p03/PAPERS/TOAA004.PDF)
12. Moses 2023, PRX 13, 041052: [10.1103/PhysRevX.13.041052](https://doi.org/10.1103/PhysRevX.13.041052)
13. Tabakov 2015, Phys. Rev. Applied 4, 031001:
    [10.1103/PhysRevApplied.4.031001](https://doi.org/10.1103/PhysRevApplied.4.031001)
14. Persaud 2017, RSI 88, 063304: [10.1063/1.4984969](https://doi.org/10.1063/1.4984969); Ji 2021, RSI 92, 103301:
    [10.1063/5.0058175](https://doi.org/10.1063/5.0058175)
15. Hwang, Shi, Lal 2013, EIPBN 5C-3: [PDF](https://eipbn.org/abstracts/2013/papers/5C-3.pdf)
16. Chlouba 2023, Nature 622, 476: [10.1038/s41586-023-06602-7](https://doi.org/10.1038/s41586-023-06602-7)
17. Pearson 2006, PRA 73, 032307: [10.1103/PhysRevA.73.032307](https://doi.org/10.1103/PhysRevA.73.032307)
18. Thomas 2025, Am. J. Phys. 93, 581: [10.1119/5.0243389](https://doi.org/10.1119/5.0243389);
    [qt3-ion-trap](https://github.com/qt3uw/qt3-ion-trap)
19. Kruit 2016, Ultramicroscopy 164, 31:
    [10.1016/j.ultramic.2016.03.004](https://doi.org/10.1016/j.ultramic.2016.03.004)
20. Hudson 2014, Phys. Plasmas 21, 010705: [10.1063/1.4863844](https://doi.org/10.1063/1.4863844)
21. Radożycki 2020, PRA 102, 063101: [10.1103/PhysRevA.102.063101](https://doi.org/10.1103/PhysRevA.102.063101)
22. Schramm 2001, PRL 87, 184801: [10.1103/PhysRevLett.87.184801](https://doi.org/10.1103/PhysRevLett.87.184801)
23. Jones 2022, NIM A 1039, 167000: [10.1016/j.nima.2022.167000](https://doi.org/10.1016/j.nima.2022.167000)
24. Taniguchi 2025, preprint: [arXiv:2502.16578](https://arxiv.org/abs/2502.16578)
25. HackRF Pro documentation: [hackrf.readthedocs.io](https://hackrf.readthedocs.io/en/latest/hackrf_pro.html)
26. Massuda 2018, ACS Photonics 5, 3513: [10.1021/acsphotonics.8b00743](https://doi.org/10.1021/acsphotonics.8b00743)
27. Rozsa 2020, RSI 91, 053302: [10.1063/5.0003020](https://doi.org/10.1063/5.0003020)
28. Campbell 2017, J. Phys. B 50, 064002: [10.1088/1361-6455/aa5a8f](https://doi.org/10.1088/1361-6455/aa5a8f)
29. Yu 2022, PRA 105, 022420: [10.1103/PhysRevA.105.022420](https://doi.org/10.1103/PhysRevA.105.022420)
30. Moehring 2011, NJP 13, 075018: [10.1088/1367-2630/13/7/075018](https://doi.org/10.1088/1367-2630/13/7/075018)
31. Wesenberg 2009, PRA 79, 013416: [10.1103/PhysRevA.79.013416](https://doi.org/10.1103/PhysRevA.79.013416)
32. Ibrahim 2017, Anal. Chem. 89, 1972: [10.1021/acs.analchem.6b04500](https://doi.org/10.1021/acs.analchem.6b04500)
33. Brown 2007, PRA 75, 015401: [10.1103/PhysRevA.75.015401](https://doi.org/10.1103/PhysRevA.75.015401)
34. Wright 2013, NJP 15, 033004: [10.1088/1367-2630/15/3/033004](https://doi.org/10.1088/1367-2630/15/3/033004)
35. ECSS-E-ST-20-01C, multipactor standard, 2020:
    [ecss.nl](https://ecss.nl/standard/ecss-e-st-20-01c-multipactor-design-and-test-15-june-2020/)
