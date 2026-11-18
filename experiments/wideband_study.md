# Wideband study: what it would take to tune the array over a wide range

The owner has since decided on a redesign for 0.4-6 GHz with all four jobs, a built-in panel plus
external antennas, and four-radio or HackRF operation (`ARCHITECTURE.md`, "Redesign
requirements"). This study is the evidence base for it: the simulated limits of the current
architecture, what KrakenSDR actually does, and the parts and antenna approaches that exist
(research verified by independent reviewers, 2026-11-15; sources at the end).

## Short answer

- **KrakenSDR is wideband in its receivers, not its antennas.** Five common-clocked RTL-SDR tuners
  cover 24 MHz-1.77 GHz (DF marketed for 100 MHz-1 GHz); the user rebuilds the antenna ring for each
  band (spacing 0.2-0.5 wavelength, so about 2.5:1 per ring), and a noise source switched onto every
  tuner input calibrates the receivers only: nothing ahead of the SMA jacks is calibrated, hence
  cables matched to 0.5-1 cm. At 5.8 GHz coax moves about 10 degrees per mm, so array16 has to keep
  calibrating at the antenna side. US$749 plus US$249 for antennas. [1-4]
- **The 26 mm 4 x 4 cannot be stretched to 0.4 GHz.** Its pitch caps unambiguous DF at 5.77 GHz;
  its 104 mm aperture is a third of a wavelength at 1 GHz; an element that fits a 26 mm cell is at
  most 3 % wide at 1 GHz (Chu); and a ground-backed array covering 0.4-6 GHz must be at least 39 mm
  thick even in theory (11 mm for 1.2-6 GHz, 4 mm for 2.4-6 GHz). A 4 x 4 is all edge: wideband
  arrays need about 5 x 5 wavelengths before edge effects fade. [5-7]
- **So the flexible design splits the band**: a built-in wideband panel for about 2.3-6 GHz, and
  external antennas, KrakenSDR-style, below that; the board between them must be broadband and
  calibrated at the antenna side.
- **Most of the RF chain has to change** (switches excepted): broadband LNAs (PMA3-83LN+, about 5x
  the power), a true 0/180-degree switch that works at every frequency (a TCM1-63AX+ transformer
  reversed by switches), MMIC or resistive combiners instead of quarter-wave Wilkinsons, and
  preselection. Multi-octave phase shifters below 6 GHz do not exist at hobby prices, so fine beam
  steering has to come from the digital side and from time delays.

`python experiments/wideband_study.py` (about 25 s) writes `wideband_results.json` and the
figures in `wideband/`. The simulator's RF parts are frequency-flat, so away from 5.8 GHz it
shows the array as if the LNAs, combiners and antennas were wideband. What remains is the
architecture: the delay-line phase states, the 104 mm aperture, the 26 mm pitch, and the
signal processing.

## What limits the band, part by part

| Part | Today | Why it limits the band | Simulated |
| --- | --- | --- | --- |
| Patches | about 3 % bandwidth (5.7-5.9 GHz) | a resonant patch is one band | `band_limits.png` |
| Element size | 26 mm cell | Chu limit: anything that fits is at most 3 % wide at 1 GHz, 0.5 % at 0.5 GHz | `band_limits.png` |
| Pitch | 26 mm | grating lobes and ambiguous directions above 5.77 GHz | `band_limits.png` |
| Aperture | 104 mm | accuracy and resolution fall as 1/f; at 1 GHz the array has no main beam to speak of | `accuracy_vs_aperture.png`, `df_vs_frequency.png` |
| LNAs | SKY65404-31, a 4.9-5.9 GHz part | band-limited match and gain | (datasheet) |
| Phase shifters | switched delay lines, 90 / 180 degrees at 5.8 GHz | they are true time delays: 180 degrees only at 5.8 GHz | `decode_noise.png`, `analog_beam.png` |
| Combiners | single-section Wilkinsons at 5.8 GHz | in-phase combining stays within 2 dB to 0.3 GHz, isolation falls to 6-10 dB below 3 GHz | `band_limits.png` |
| Radios | HackRF; AD9361-family boards | already wideband | (datasheets) |

## The delay lines, and the fix already made in the code

The snapshot method switches each element's 180-degree line in Hadamard patterns. The line is
a fixed delay, so at 2.4 GHz it flips the phase by only 74 degrees, and the old decode, which
assumed -1, mixed the elements up. The host code now decodes with the actual step and picks
beam states by their actual delay (`control.toggle_factor`, `control.state_lag`; identical at
5.8 GHz, all earlier results unchanged). Remaining cost: the smaller the step, the more the
decode amplifies noise.

| Frequency | Toggle step | Decode noise penalty (four radios / HackRF) | DF error, four radios, median / 95th pct | Same with the old decode |
| --- | --- | --- | --- | --- |
| 5.8 GHz | 180 deg | 0 / 0 dB | 0.10 / 0.43 deg | same |
| 4.5 GHz | 140 deg | 0.8 / 1.0 dB | 0.13 / 0.47 deg | 0.71 / 2.44 deg |
| 3.5 GHz | 109 deg | 2.5 / 2.9 dB | 0.15 / 0.52 deg | 0.88 / 3.30 deg |
| 2.4 GHz | 74 deg | 5.6 / 6.3 dB | 0.22 / 0.73 deg | 0.96 / 3.91 deg |
| 1.3 GHz | 40 deg | 10.8 / 11.7 dB | 0.31 / 1.09 deg | |
| 0.9 GHz | 28 deg | 14.0 / 15.0 dB | 0.47 / 1.55 deg | |

Single source at -75 dBm, calibrated on a boresight source at each frequency, wideband RF
assumed. HackRF mode, with one receiver and 16 patterns, suffers more: 0.65 deg (95th pct) at
5.8 GHz, 3.1 deg at 2.4 GHz, 9 deg at 0.9 GHz. These are noise-limited numbers; calibration
errors and multipath, which grow with the beamwidth (that is, as 1/f for a fixed aperture),
come on top.

Analog beams (HackRF mode, transmit) with the 2-bit delay lines: broadside and about 30 degrees
(one 90-degree step per element, a true-time-delay beam) are exact at every frequency; other
angles lose up to 1.2 dB and grow quantisation lobes, worst between 4 and 6 GHz.

## Phase control for a 0.4-6 GHz board (`phase_control.png`)

Each scheme gives an element a set of phases; the beam picks the nearest for each element.
External antennas are taken as half a wavelength apart at every band (the KrakenSDR way), and
elements as cos^1.2. Worst sidelobe is over scans of 0-45 degrees.

| Per-element phase control | DF decode penalty | Beam gain lost | Worst sidelobe |
| --- | --- | --- | --- |
| Delay lines 90 / 180 deg at 5.8 GHz (today) | 17.5 dB at 0.4 GHz, 10.5 at 0.9, 2.3 at 2.4 | 7.5 dB at 0.4 GHz to 0.4 dB at 5.8 | -10.6 to -1.9 dB |
| Exact 0/180 switch (wideband balun) | 0 dB at every f | 2.1 dB | 0 dB: real weights make a mirror beam at -theta |
| Exact 0/180 plus one 90-deg-at-5.8-GHz delay | 0 dB | 1.8 dB at 0.4 GHz to 0.4 dB at 5.8 | -0.2 to -3.5 dB |
| 3-bit phase shifter, 45 deg at every f | 0 dB | 0.1 dB | -7.7 dB |
| Vector modulator (continuous) | 0 dB | 0 dB | -10.2 dB |

Direction finding needs only an exact 180-degree flip at every frequency, which a switched
wideband balun gives cheaply. Steered beams, and above all transmit beams, need at least three
bits of phase that hold across the whole band; that part decides the redesign's cost.

## Parts for a broadband chain (research, prices November 2026)

| Function | Today (5.8 GHz) | Broadband option | Cost / catch |
| --- | --- | --- | --- |
| Switches | BGS12WN6 | the same: rated 0.05-9 GHz [8] | $0.32; 0.15-0.46 dB loss |
| LNA | SKY65404-31, 4.9-5.9 GHz, 36 mW | PMA3-83LN+, 0.5-8 GHz, 21 dB flat, NF 1.3-1.9 dB [9] | $15.71, 300 mW each: 16 of them draw about 5 W |
| | | PMA-5451+, 0.05-6 GHz, NF 0.6-2.3 dB [10] | $1.96, but gain falls from 22 to 5 dB: needs a second stage above 3 GHz |
| 0/180 deg for DF | 86 ps delay line (180 deg only at 5.8 GHz) | TCM1-63AX+ 1:1 transformer, 10-6000 MHz, reversed by two SPDTs [11] | about $9; 1.3 dB loss, 8 deg / 0.5 dB unbalance (calibrated out) |
| Fine phase | 2-bit delay lines | no multi-octave digital phase shifter below 6 GHz [12]; vector modulators are narrowband (ADL5390 20 MHz-2.4 GHz) [13]; ADAR4002 true-time delay 0.5-19 GHz [14] | ADAR4002 $112 at 1k and 16-20 dB loss per element |
| Combining | quarter-wave Wilkinsons | MMIC splitters EP2W1+ (0.5-9.5 GHz), EP4RKU+ (DC-18 GHz) [15]; resistive splitters [16] | $29 / $55; resistive costs 12 dB per 4-way but only ~0.2 dB of SNR after 18 dB of gain. A 0.4-6 GHz multi-section Wilkinson needs ~100 mm per arm [17] |
| Preselection | none needed in one band | high-pass near 400 MHz per element; filter bank or a tunable filter (ADMV8809, 0.5-9 GHz) per radio channel; limiter before the HackRF [18, 19] | ADMV8809 $101 at 1k; the HackRF is damaged above -5 dBm |

Neither KrakenSDR nor HackRF has preselection built in: a broadband array with ~40 dB of gain
in front of the radio will need it near strong transmitters. [19]

## Antennas for a wide range (research)

- **Built-in panel**: tightly coupled dipole (TCDA) or PUMA-type PCB elements reach 5:1 in large
  arrays, but a measured 2 x 16 TCDA for 1.2-6 GHz had edge elements at VSWR above 4 at the low end
  until resistive loading was added; a 0.7-6 GHz BAVA DF array is 254 x 254 x 25 mm. For a 4 x 4 at
  26 mm, about 2.3-6 GHz with a 5-10 mm thick element and a loaded dummy ring is realistic. [6, 20, 21]
- **External antennas**: KrakenSDR-style rings of identical vertical whips or dipoles at spacing
  about 0.33 wavelength (ring radius 194 mm at 433 MHz, 92 mm at 915 MHz), calibrated in place by
  radiating a PN code from a central antenna; a 4 x 4 at half-wavelength spacing would be 655 mm at
  915 MHz and 1.4 m at 433 MHz. [3]
- **Interchangeable panels** should be complete band-specific front ends on blind-mate SMP
  connectors (U.FL is rated for about 30 matings), each with its own stored array response. [22]
- Rohde & Schwarz cover 20 MHz-8.5 GHz in one 1.1 m radome by nesting active and passive element
  sets: the professional answer is also several arrays, not one. [23]

## Radios and transmit rules (research)

- **Radios that reach 6 GHz.** A stock Pluto (AD9363) is specified for 325 MHz-3.8 GHz and 20 MHz
  and is guaranteed only as one receive / one transmit channel; its "AD9364 mode" (70 MHz-6 GHz) is
  outside the datasheet, and USB 2.0 sustains about 4 MSPS. The bladeRF 2.0 micro xA4 (AD9361,
  47 MHz-6 GHz, 2 x 2, 56 MHz, USB 3.0) costs US$540: two give four coherent channels for US$1,080.
  Its reference input accepts 5-300 MHz, so the board's 40 MHz oscillator can discipline both, and
  its trigger starts chained boards within one sample. KrakenSDR, LimeSDR and a stock Pluto cannot
  reach 5.8 GHz. [26-28]
- **The AD9361 has no RF phase synchronisation.** The phase between two chips re-randomises on every
  LO change, sample-rate change, some gain changes, in TDD mode and with quadrature tracking; even a
  shared external LO leaves a 0/180-degree ambiguity. So: FDD, fixed gain, quadrature tracking frozen
  after calibration, and the PN calibration after every such event, as array16 already does after a
  retune. Within one chip both receive channels share the LO, so two boards add one unknown phase and
  delay per retune. Radios with shared LOs exist but cost US$2,900-26,000 per channel and still
  need periodic calibration. [29]
- **HackRF-only mode** has one receiver and one LO, so a retune moves all 16 elements together and
  needs no re-alignment: a stored per-frequency element table is enough. The HackRF cannot send the
  calibration code itself (half duplex), and its clock input takes only 10 MHz. [30]
- **Transmitting** is band-specific whatever the receiver does. US licence-free: 902-928,
  2400-2483.5 and 5725-5875 MHz at -1.25 dBm EIRP average (+18.75 dBm peak) under 15.249; motion
  sensors at 5785-5815 MHz at 18.75 dBm EIRP average (about 38.75 dBm peak when pulsed) under 15.245.
  EU: 25 mW EIRP at 5725-5875 MHz, measured during transmission, with a time-out required for
  continuous carriers. These limits include the array's gain (17.1 dBi), so 16 elements may carry only
  about -15 dBm each in the EU. The second harmonic of 5.8 GHz falls in a restricted band and needs
  about 60 dB of filtering. A US amateur licence opens 902-928, 1240-1300, 2300-2310, 2390-2450 and
  5650-5925 MHz (pulse emissions allowed except at 23 cm), the practical route to useful radar or
  beam-steering power. [31-34]

## Aperture: why KrakenSDR moves its antennas

A fixed aperture's noise-limited accuracy scales as wavelength / aperture. The 4 x 4 is 104 mm:
at 5.8 GHz that is 2 wavelengths, at 1 GHz a third of one. KrakenSDR instead asks the user to
place five antennas on a circle sized for the band (spacing about half a wavelength), which keeps
accuracy constant in degrees at the cost of a fixture that grows to metres at VHF. With the same
per-element SNR, a resized 5-element circle beats the fixed 4 x 4 below about 1.7 GHz
(`accuracy_vs_aperture.png`); a circle left at its 1 GHz size is better above 1 GHz but its
directions become ambiguous there.

## Sources

1. KrakenSDR datasheet: https://github.com/krakenrf/krakensdr_docs/blob/main/paper_manual/KrakenSDR%20Datasheet.pdf
2. KrakenSDR frequency range (forum): https://forum.krakenrf.com/t/krakensdr-frequency-operating-range/1468
3. KrakenSDR antenna array setup: https://github.com/krakenrf/krakensdr_docs/wiki/04.-Antenna-Array-Setup
4. KrakenSDR store page: https://www.krakenrf.com/product-page/krakensdr
5. Jonsson, Kolitsidas and Hussain, array thickness bound: https://ar5iv.arxiv.org/html/1307.2925
6. Holter and Steyskal, finite wideband arrays: https://doi.org/10.1109/TAP.2002.1017665
7. Chu limit, as applied in `wideband_study.py` (chu_bandwidth)
8. Infineon BGS12WN6: https://www.infineon.com/dgdl/Infineon-BGS12WN6-DataSheet-v01_00-EN.pdf
9. Mini-Circuits PMA3-83LN+: https://www.minicircuits.com/pdfs/PMA3-83LN+.pdf
10. Mini-Circuits PMA-5451+: https://www.minicircuits.com/pdfs/PMA-5451+.pdf
11. Mini-Circuits TCM1-63AX+: https://www.minicircuits.com/pdfs/TCM1-63AX+.pdf
12. Analog Devices phase shifters: https://www.analog.com/en/product-category/phase-shifters.html
13. Analog Devices ADL5390: https://www.analog.com/en/products/adl5390.html
14. Analog Devices ADAR4002: https://www.analog.com/en/products/adar4002.html
15. Mini-Circuits EP4RKU+: https://www.minicircuits.com/pdfs/EP4RKU+.pdf
16. Resistive splitters: https://www.microwaves101.com/encyclopedias/resistive-power-splitters
17. Multi-section Wilkinsons: https://www.microwaves101.com/encyclopedias/how-many-sections-do-you-need
18. Analog Devices ADMV8809: https://www.analog.com/en/products/admv8809.html
19. KrakenSDR overview (filtering): https://github.com/krakenrf/krakensdr_docs/wiki/01.-Overview
20. 1.2-6 GHz TCDA edge effects: https://arxiv.org/abs/2601.01323v1
21. BAVA DF array patent: https://patents.google.com/patent/US9513361B1/en
22. Samtec SMP: https://suddendocs.samtec.com/productspecs/smp-sm.pdf
23. Rohde & Schwarz ADD597: https://scdn.rohde-schwarz.com/ur/pws/dl_downloads/pdm/cl_brochures_and_datasheets/specifications/3609_5680_22/ADD597_dat-sw_en_3609-5680-22_v0300.pdf
24. Switched UCA DF at 2.4 GHz (Werner et al.): https://arxiv.org/html/2502.08592
25. QuadRF / open-space-sdr: https://github.com/open-space-sdr/main
26. ADALM-Pluto specifications: https://wiki.analog.com/university/tools/pluto/devs/specs
27. bladeRF 2.0 micro: https://www.nuand.com/bladerf-2-0-micro/
28. libbladeRF triggers: https://raw.githubusercontent.com/Nuand/bladeRF/master/host/libraries/libbladeRF/include/libbladeRF.h
29. ADI FMCOMMS5 phase synchronisation: https://analogdevicesinc.github.io/system-level/solutions/reference-designs/fmcomms5/phase-sync/
30. HackRF One: https://hackrf.readthedocs.io/en/latest/hackrf_one.html
31. 47 CFR 15.249: https://www.law.cornell.edu/cfr/text/47/15.249
32. 47 CFR 15.245: https://www.law.cornell.edu/cfr/text/47/15.245
33. ETSI EN 300 440 V2.2.1: https://www.etsi.org/deliver/etsi_en/300400_300499/300440/02.02.01_60/en_300440v020201p.pdf
34. 47 CFR 97.301 (amateur bands): https://www.law.cornell.edu/cfr/text/47/97.301
