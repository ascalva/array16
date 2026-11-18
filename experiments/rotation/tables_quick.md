## 0. Cross-checks

az-only model against evolve_layout.check() (median / 95th / gross)
| case | evolve_study.md | evolve_layout.check now | rotation_study.az_mc, K = 1 |
| --- | --- | --- | --- |
| 4 whips, r 0.25 wl | 0.93 / 2.62 / 0.0 % | 0.95 / 2.9 / 0.0 % | 0.95 / 2.9 / 0.0 % |
| 4 whips, r 0.354 wl | 0.75 / 170.7 / 10.6 % | 0.73 / 163.0 / 7.0 % | 0.73 / 163.0 / 7.0 % |
| 8 whips, r 0.5 wl | 0.33 / 0.93 / 0.0 % | 0.34 / 1.0 / 0.0 % | 0.34 / 1.0 / 0.0 % |
| 8 whips, r 225 mm, 2.3 GHz, 0-30 deg up | (30 % gross worst in band) | 0.13 / 126.7 / 26.7 % | 0.13 / 126.7 / 26.7 % |

sphere model, one orientation: great-circle error median / 95th (deg); REDESIGN section 9 quotes 0.12/0.49, 0.12/0.49, 0.10/0.46 for the panel
| case | one snapshot at SNR + 18 dB | 64 samples, covariance |
| --- | --- | --- |
| panel 2.44 GHz, source within 45 deg az / 30 deg el of boresight | 0.65 / 1.54 | 0.70 / 1.61 |
| panel 3.50 GHz, source within 45 deg az / 30 deg el of boresight | 0.46 / 1.07 | 0.47 / 1.15 |
| panel 5.80 GHz, source within 45 deg az / 30 deg el of boresight | 0.28 / 0.64 | 0.30 / 0.62 |
| plate 1.575 GHz, sources 20-70 deg up | 0.54 / 122.63 | 0.58 / 129.97 |
[section 0: 1 s]

## a. Fragile / ambiguous layouts turned about the vertical axis, K orientations combined

a1. Equal time: 64 samples in total, split over K orientations that fill one symmetry sector (step = 360 / (symmetry x K): 45/K deg for a ring of 8, 90/K deg for 4 whips or a square grid). median / 95th percentile (deg) / gross
| layout | K = 1 | K = 2 | K = 4 | K = 8 |
| --- | --- | --- | --- | --- |
| ring of 8, 225 mm, 2.3 GHz, sources 0-30 deg up | 0.16 / 127.6 / 26.7 % | 0.13 / 120.3 / 9.0 % | 0.11 / 0.4 / 4.0 % | 0.10 / 0.3 / 1.7 % |
| ring of 8, 225 mm, 1.8 GHz, sources 0-30 deg up | 0.16 / 75.5 / 11.0 % | 0.14 / 0.4 / 0.7 % | 0.14 / 0.4 / 0.0 % | 0.14 / 0.4 / 0.0 % |
| 4 whips, r 0.354 wl, sources 0-15 deg up | 0.81 / 171.9 / 8.0 % | 0.67 / 1.8 / 0.0 % | 0.61 / 1.8 / 0.0 % | 0.61 / 1.7 / 0.0 % |
| 4 whips, r 0.45 wl | 0.61 / 103.7 / 6.3 % | 0.50 / 1.4 / 0.0 % | 0.50 / 1.4 / 0.0 % | 0.48 / 1.3 / 0.0 % |
| 4 whips, r 0.60 wl | 0.50 / 113.1 / 11.3 % | 0.38 / 1.1 / 0.0 % | 0.38 / 1.1 / 0.0 % | 0.36 / 1.0 / 0.0 % |
| 4 whips, r 1.0 wl | 0.34 / 105.9 / 18.0 % | 0.22 / 0.7 / 0.0 % | 0.22 / 0.6 / 0.0 % | 0.21 / 0.6 / 0.0 % |
| 4x4 grid, 164 mm pitch at 2.3 GHz (1.26 wl) | 0.07 / 105.9 / 13.7 % | 0.06 / 0.2 / 0.0 % | 0.06 / 0.2 / 0.0 % | 0.05 / 0.2 / 0.0 % |
| ring of 8, 225 mm at 5.8 GHz (spacing 3.3 wl) | 0.06 / 123.6 / 19.3 % | 0.05 / 0.1 / 1.0 % | 0.04 / 0.1 / 0.0 % | 0.04 / 0.1 / 0.0 % |

a2. K times the time: 64 samples at every orientation, same steps
| layout | K = 1 | K = 2 | K = 4 | K = 8 |
| --- | --- | --- | --- | --- |
| ring of 8, 225 mm, 2.3 GHz, sources 0-30 deg up | 0.16 / 127.6 / 26.7 % | 0.12 / 123.0 / 9.0 % | 0.12 / 0.5 / 4.7 % | 0.10 / 0.3 / 0.3 % |
| ring of 8, 225 mm, 1.8 GHz, sources 0-30 deg up | 0.16 / 75.5 / 11.0 % | 0.15 / 0.4 / 0.7 % | 0.14 / 0.4 / 0.0 % | 0.14 / 0.4 / 0.0 % |
| 4 whips, r 0.354 wl, sources 0-15 deg up | 0.81 / 171.9 / 8.0 % | 0.67 / 1.9 / 0.0 % | 0.62 / 1.8 / 0.0 % | 0.58 / 1.7 / 0.0 % |
| 4 whips, r 0.45 wl | 0.61 / 103.7 / 6.3 % | 0.53 / 1.5 / 0.0 % | 0.49 / 1.4 / 0.0 % | 0.45 / 1.4 / 0.0 % |
| 4 whips, r 0.60 wl | 0.50 / 113.1 / 11.3 % | 0.39 / 1.1 / 0.0 % | 0.37 / 1.1 / 0.0 % | 0.35 / 1.0 / 0.0 % |
| 4 whips, r 1.0 wl | 0.34 / 105.9 / 18.0 % | 0.24 / 0.6 / 0.0 % | 0.21 / 0.6 / 0.0 % | 0.21 / 0.6 / 0.0 % |
| 4x4 grid, 164 mm pitch at 2.3 GHz (1.26 wl) | 0.07 / 105.9 / 13.7 % | 0.06 / 0.2 / 0.0 % | 0.06 / 0.2 / 0.0 % | 0.06 / 0.2 / 0.0 % |
| ring of 8, 225 mm at 5.8 GHz (spacing 3.3 wl) | 0.06 / 123.6 / 19.3 % | 0.05 / 0.1 / 0.3 % | 0.04 / 0.1 / 0.0 % | 0.04 / 0.1 / 0.0 % |

a3. Equal time, naive steps of 360 / K deg (a symmetric layout lands on itself)
| layout | K = 1 | K = 2 | K = 4 | K = 8 |
| --- | --- | --- | --- | --- |
| ring of 8, 225 mm, 2.3 GHz, sources 0-30 deg up | 0.16 / 127.6 / 26.7 % | 0.05 / 126.9 / 26.7 % | 0.03 / 126.9 / 26.3 % | 0.03 / 126.9 / 26.3 % |
| ring of 8, 225 mm, 1.8 GHz, sources 0-30 deg up | 0.16 / 75.5 / 11.0 % | 0.05 / 75.7 / 10.7 % | 0.03 / 75.5 / 10.7 % | 0.04 / 75.5 / 10.0 % |
| 4 whips, r 0.354 wl, sources 0-15 deg up | 0.81 / 171.9 / 8.0 % | 0.15 / 175.9 / 7.0 % | 0.15 / 175.5 / 7.0 % | 0.15 / 0.4 / 0.0 % |
| 4 whips, r 0.45 wl | 0.61 / 103.7 / 6.3 % | 0.11 / 0.5 / 4.3 % | 0.11 / 0.4 / 4.0 % | 0.11 / 0.3 / 0.0 % |
| 4 whips, r 0.60 wl | 0.50 / 113.1 / 11.3 % | 0.09 / 74.0 / 7.7 % | 0.09 / 73.2 / 7.7 % | 0.08 / 0.2 / 0.0 % |
| 4 whips, r 1.0 wl | 0.34 / 105.9 / 18.0 % | 0.06 / 105.1 / 14.3 % | 0.06 / 105.2 / 14.3 % | 0.05 / 0.2 / 0.0 % |
| 4x4 grid, 164 mm pitch at 2.3 GHz (1.26 wl) | 0.07 / 105.9 / 13.7 % | 0.02 / 105.8 / 13.7 % | 0.02 / 105.8 / 13.3 % | 0.02 / 0.1 / 0.0 % |
| ring of 8, 225 mm at 5.8 GHz (spacing 3.3 wl) | 0.06 / 123.6 / 19.3 % | 0.02 / 113.1 / 20.3 % | 0.02 / 113.0 / 20.0 % | 0.02 / 123.6 / 20.0 % |

a4. The elevated-source failure of the 225 mm ring is an estimator fault: the same data searched over azimuth and elevation (0-60 deg)
| layout | az-only, K = 1 | az-el search, K = 1 | az-el search, K = 4 |
| --- | --- | --- | --- |
| ring of 8, 225 mm, 2.3 GHz, sources 0-30 deg up | 0.16 / 127.6 / 26.7 % | 0.11 / 0.3 / 0.0 % | 0.10 / 0.3 / 0.0 % |
| ring of 8, 225 mm, 1.8 GHz, sources 0-30 deg up | 0.16 / 75.5 / 11.0 % | 0.14 / 0.4 / 0.0 % | 0.14 / 0.4 / 0.0 % |
[section a: 9 s]

## b. Orientation error and cable-flex phase drift

b. ring of 8, 225 mm, 2.3 GHz, sources 0-30 deg up: 64 samples in total; median / 95th (deg) / gross. The first orientation is the surveyed home (exact, no drift); errors apply to the moves
| condition | K = 1 | K = 2 | K = 4 | K = 8 |
| --- | --- | --- | --- | --- |
| perfect | 0.15 / 127.3 / 25.7 % | 0.11 / 122.4 / 7.3 % | 0.12 / 6.3 / 5.0 % | 0.11 / 0.3 / 0.0 % |
| angle 0.1 deg | 0.15 / 127.3 / 25.7 % | 0.12 / 122.5 / 7.3 % | 0.11 / 0.4 / 4.0 % | 0.11 / 0.3 / 0.3 % |
| angle 0.5 deg | 0.15 / 127.3 / 25.7 % | 0.23 / 122.6 / 7.7 % | 0.17 / 0.6 / 4.0 % | 0.18 / 0.5 / 0.3 % |
| angle 1 deg | 0.15 / 127.3 / 25.7 % | 0.42 / 119.8 / 7.7 % | 0.30 / 1.0 / 4.0 % | 0.28 / 0.8 / 1.0 % |
| flex 1 deg | 0.15 / 127.3 / 25.7 % | 0.12 / 57.0 / 7.0 % | 0.11 / 0.4 / 3.7 % | 0.11 / 0.3 / 1.0 % |
| flex 3 deg | 0.15 / 127.3 / 25.7 % | 0.12 / 57.0 / 7.0 % | 0.11 / 0.4 / 4.0 % | 0.11 / 0.3 / 1.0 % |
| flex 10 deg | 0.15 / 127.3 / 25.7 % | 0.13 / 56.9 / 7.0 % | 0.12 / 0.4 / 4.0 % | 0.12 / 0.3 / 1.0 % |
| angle 0.5 deg + flex 3 deg | 0.15 / 127.3 / 25.7 % | 0.24 / 57.1 / 7.0 % | 0.20 / 0.7 / 4.0 % | 0.15 / 0.5 / 0.0 % |

b. 4 whips, r 0.45 wl: 64 samples in total; median / 95th (deg) / gross. The first orientation is the surveyed home (exact, no drift); errors apply to the moves
| condition | K = 1 | K = 2 | K = 4 | K = 8 |
| --- | --- | --- | --- | --- |
| perfect | 0.63 / 103.3 / 5.3 % | 0.51 / 1.4 / 0.0 % | 0.48 / 1.4 / 0.0 % | 0.49 / 1.3 / 0.0 % |
| angle 0.1 deg | 0.63 / 103.3 / 5.3 % | 0.50 / 1.4 / 0.0 % | 0.48 / 1.3 / 0.0 % | 0.51 / 1.3 / 0.0 % |
| angle 0.5 deg | 0.63 / 103.3 / 5.3 % | 0.55 / 1.6 / 0.0 % | 0.50 / 1.3 / 0.0 % | 0.48 / 1.3 / 0.0 % |
| angle 1 deg | 0.63 / 103.3 / 5.3 % | 0.60 / 1.8 / 0.0 % | 0.55 / 1.5 / 0.0 % | 0.52 / 1.4 / 0.0 % |
| flex 1 deg | 0.63 / 103.3 / 5.3 % | 0.48 / 1.4 / 0.0 % | 0.48 / 1.3 / 0.0 % | 0.51 / 1.3 / 0.0 % |
| flex 3 deg | 0.63 / 103.3 / 5.3 % | 0.50 / 1.5 / 0.0 % | 0.48 / 1.4 / 0.0 % | 0.52 / 1.3 / 0.0 % |
| flex 10 deg | 0.63 / 103.3 / 5.3 % | 0.65 / 2.0 / 0.0 % | 0.68 / 1.9 / 0.0 % | 0.61 / 1.6 / 0.0 % |
| angle 0.5 deg + flex 3 deg | 0.63 / 103.3 / 5.3 % | 0.56 / 1.4 / 0.0 % | 0.54 / 1.6 / 0.0 % | 0.48 / 1.3 / 0.0 % |

b. 4x4 grid, 164 mm pitch at 2.3 GHz (1.26 wl): 64 samples in total; median / 95th (deg) / gross. The first orientation is the surveyed home (exact, no drift); errors apply to the moves
| condition | K = 1 | K = 2 | K = 4 | K = 8 |
| --- | --- | --- | --- | --- |
| perfect | 0.08 / 69.4 / 12.3 % | 0.06 / 0.2 / 0.0 % | 0.06 / 0.2 / 0.0 % | 0.06 / 0.2 / 0.0 % |
| angle 0.1 deg | 0.08 / 69.4 / 12.3 % | 0.07 / 0.2 / 0.0 % | 0.07 / 0.2 / 0.0 % | 0.07 / 0.2 / 0.0 % |
| angle 0.5 deg | 0.08 / 69.4 / 12.3 % | 0.18 / 0.5 / 0.0 % | 0.15 / 0.5 / 0.0 % | 0.13 / 0.3 / 0.0 % |
| angle 1 deg | 0.08 / 69.4 / 12.3 % | 0.35 / 0.9 / 0.0 % | 0.28 / 0.9 / 0.0 % | 0.25 / 0.6 / 0.0 % |
| flex 1 deg | 0.08 / 69.4 / 12.3 % | 0.06 / 0.2 / 0.0 % | 0.06 / 0.2 / 0.0 % | 0.06 / 0.2 / 0.0 % |
| flex 3 deg | 0.08 / 69.4 / 12.3 % | 0.07 / 0.2 / 0.0 % | 0.06 / 0.2 / 0.0 % | 0.07 / 0.2 / 0.0 % |
| flex 10 deg | 0.08 / 69.4 / 12.3 % | 0.08 / 0.2 / 0.0 % | 0.07 / 0.2 / 0.0 % | 0.07 / 0.2 / 0.0 % |
| angle 0.5 deg + flex 3 deg | 0.08 / 69.4 / 12.3 % | 0.18 / 0.5 / 0.0 % | 0.16 / 0.5 / 0.0 % | 0.12 / 0.4 / 0.0 % |

b. two rings of 8 (225 / 112 mm), 1.575 GHz, sources 0-15 deg up: a layout with nothing to fix
| condition | K = 1 | K = 2 | K = 4 | K = 8 |
| --- | --- | --- | --- | --- |
| perfect | 0.12 / 0.3 / 0.0 % | 0.11 / 0.3 / 0.0 % | 0.13 / 0.4 / 0.0 % | 0.12 / 0.4 / 0.0 % |
| angle 0.1 deg | 0.12 / 0.3 / 0.0 % | 0.12 / 0.4 / 0.0 % | 0.13 / 0.4 / 0.0 % | 0.12 / 0.4 / 0.0 % |
| angle 0.5 deg | 0.12 / 0.3 / 0.0 % | 0.18 / 0.7 / 0.0 % | 0.19 / 0.6 / 0.0 % | 0.17 / 0.5 / 0.0 % |
| angle 1 deg | 0.12 / 0.3 / 0.0 % | 0.32 / 1.1 / 0.0 % | 0.32 / 0.9 / 0.0 % | 0.27 / 0.7 / 0.0 % |
| flex 1 deg | 0.12 / 0.3 / 0.0 % | 0.12 / 0.4 / 0.0 % | 0.12 / 0.3 / 0.0 % | 0.12 / 0.4 / 0.0 % |
| flex 3 deg | 0.12 / 0.3 / 0.0 % | 0.12 / 0.4 / 0.0 % | 0.12 / 0.4 / 0.0 % | 0.12 / 0.4 / 0.0 % |
| flex 10 deg | 0.12 / 0.3 / 0.0 % | 0.14 / 0.4 / 0.0 % | 0.13 / 0.4 / 0.0 % | 0.13 / 0.4 / 0.0 % |
| angle 0.5 deg + flex 3 deg | 0.12 / 0.3 / 0.0 % | 0.21 / 0.6 / 0.0 % | 0.17 / 0.6 / 0.0 % | 0.17 / 0.5 / 0.0 % |
[section b: 19 s]

## c. The panel on one axis (azimuth)

c1. Fixed panel: azimuth DF and gain against the source's azimuth off boresight (source 0-15 deg up, V-pol, one capture; gain = REDESIGN boresight gain + 15 log10 cos, a model, with a -20 dB floor behind)
|  | 0 deg | 15 deg | 30 deg | 45 deg | 60 deg | 75 deg | 90 deg | 120 deg | 180 deg |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2.44 GHz: az error median / 95th (deg) / gross | 0.29 / 0.9 / 0.0 % | 0.39 / 1.0 / 0.0 % | 0.36 / 1.2 / 0.0 % | 0.54 / 1.7 / 0.0 % | 0.71 / 2.1 / 0.0 % | 1.51 / 5.6 / 0.0 % | 11.53 / 179.7 / 67.0 % | 95.27 / 159.2 / 100.0 % | 94.01 / 170.4 / 100.0 % |
| 2.44 GHz: gain (dBi) | 6.4 | 6.2 | 5.5 | 4.1 | 1.9 | -2.4 | <= -13.6 | <= -13.6 | <= -13.6 |
| 3.50 GHz: az error median / 95th (deg) / gross | 0.20 / 0.7 / 0.0 % | 0.28 / 0.7 / 0.0 % | 0.24 / 0.8 / 0.0 % | 0.36 / 1.1 / 0.0 % | 0.52 / 1.5 / 0.0 % | 0.96 / 2.8 / 0.0 % | 9.65 / 179.0 / 62.0 % | 113.21 / 174.2 / 100.0 % | 118.47 / 173.9 / 100.0 % |
| 3.50 GHz: gain (dBi) | 10.7 | 10.5 | 9.8 | 8.4 | 6.2 | 1.9 | <= -9.3 | <= -9.3 | <= -9.3 |
| 5.80 GHz: az error median / 95th (deg) / gross | 0.13 / 0.4 / 0.0 % | 0.17 / 0.4 / 0.0 % | 0.17 / 0.4 / 0.0 % | 0.24 / 0.7 / 0.0 % | 0.32 / 0.8 / 0.0 % | 0.65 / 2.1 / 0.0 % | 12.70 / 178.9 / 57.0 % | 118.45 / 164.5 / 100.0 % | 145.32 / 172.3 / 100.0 % |
| 5.80 GHz: gain (dBi) | 15.5 | 15.3 | 14.6 | 13.2 | 11.0 | 6.7 | <= -4.5 | <= -4.5 | <= -4.5 |

c2. Panel stepped to K azimuths (one capture each), sources anywhere in azimuth, 0-15 deg up
| band | K | step | az error median / 95th / gross | azimuths found within 1 deg | worst-azimuth gain |
| --- | --- | --- | --- | --- | --- |
| 2.44 GHz | 1 | 360 deg | 18.08 / 170.5 / 53.0 % | 29 % | none (gaps) |
| 2.44 GHz | 2 | 180 deg | 0.69 / 3.1 / 1.0 % | 67 % | none (gaps) |
| 2.44 GHz | 3 | 120 deg | 0.26 / 0.9 / 0.0 % | 96 % | 1.9 dBi |
| 2.44 GHz | 4 | 90 deg | 0.27 / 0.7 / 0.0 % | 98 % | 4.1 dBi |
| 2.44 GHz | 6 | 60 deg | 0.30 / 0.8 / 0.0 % | 98 % | 5.5 dBi |
| 2.44 GHz | 8 | 45 deg | 0.33 / 0.8 / 0.0 % | 98 % | 5.9 dBi |
| 5.80 GHz | 1 | 360 deg | 37.06 / 174.1 / 54.0 % | 44 % | none (gaps) |
| 5.80 GHz | 2 | 180 deg | 0.26 / 1.5 / 0.5 % | 91 % | none (gaps) |
| 5.80 GHz | 3 | 120 deg | 0.16 / 0.4 / 0.0 % | 100 % | 11.0 dBi |
| 5.80 GHz | 4 | 90 deg | 0.16 / 0.4 / 0.0 % | 100 % | 13.2 dBi |
| 5.80 GHz | 6 | 60 deg | 0.15 / 0.4 / 0.0 % | 100 % | 14.6 dBi |
| 5.80 GHz | 8 | 45 deg | 0.18 / 0.4 / 0.0 % | 100 % | 15.0 dBi |

c3. Time for one 360 deg look (K stops, moving between them at the given rate; acceleration ignored); last column: one continuous turn capturing on the fly
| slew rate | K stops | dwell 0.05 s | 0.5 s | 5 s | continuous |
| --- | --- | --- | --- | --- | --- |
| 6 deg/s | 3 | 40.1 s | 41.5 s | 55.0 s | 60 s |
| 6 deg/s | 4 | 45.2 s | 47.0 s | 65.0 s | 60 s |
| 6 deg/s | 6 | 50.3 s | 53.0 s | 80.0 s | 60 s |
| 6 deg/s | 8 | 52.9 s | 56.5 s | 92.5 s | 60 s |
| 30 deg/s | 3 | 8.2 s | 9.5 s | 23.0 s | 12 s |
| 30 deg/s | 4 | 9.2 s | 11.0 s | 29.0 s | 12 s |
| 30 deg/s | 6 | 10.3 s | 13.0 s | 40.0 s | 12 s |
| 30 deg/s | 8 | 10.9 s | 14.5 s | 50.5 s | 12 s |
| 90 deg/s | 3 | 2.8 s | 4.2 s | 17.7 s | 4 s |
| 90 deg/s | 4 | 3.2 s | 5.0 s | 23.0 s | 4 s |
| 90 deg/s | 6 | 3.6 s | 6.3 s | 33.3 s | 4 s |
| 90 deg/s | 8 | 3.9 s | 7.5 s | 43.5 s | 4 s |
[section c: 4 s]

## d. Two axes: elevation near the horizon

d1. Whip plate (two rings of 8), 1.575 GHz, V-pol source at the given elevation, 4 captures: median elevation error with the sign folded (share of answers below the horizon) / median azimuth error, deg
| plate | el 0 | el 2 | el 5 | el 10 | el 20 | el 30 | el 45 | el 70 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| flat, fixed | 1.83 (34 %) / 0.12 | 1.94 (35 %) / 0.14 | 1.39 (38 %) / 0.10 | 0.71 (54 %) / 0.11 | 0.34 (57 %) / 0.13 | 0.25 (53 %) / 0.16 | 0.22 (52 %) / 0.17 | 0.13 (53 %) / 0.39 |
| flat, 4 azimuths (one axis) | 0.69 (30 %) / 0.01 | 0.31 (46 %) / 0.02 | 0.22 (56 %) / 0.01 | 0.08 (51 %) / 0.02 | 0.04 (50 %) / 0.02 | 0.03 (51 %) / 0.02 | 0.02 (49 %) / 0.03 | 0.04 (50 %) / 0.11 |
| tilted 30 deg, fixed | 14.05 (36 %) / 2.55 | 1.51 (25 %) / 1.22 | 1.88 (20 %) / 1.22 | 3.72 (18 %) / 2.36 | 3.57 (26 %) / 1.34 | 2.98 (22 %) / 0.80 | 1.05 (30 %) / 2.18 | 0.71 (19 %) / 3.77 |
| tilted 30 deg, 4 azimuths | 0.17 (5 %) / 0.02 | 0.22 (0 %) / 0.01 | 0.20 (0 %) / 0.02 | 0.24 (0 %) / 0.03 | 0.16 (0 %) / 0.05 | 0.13 (0 %) / 0.10 | 0.06 (0 %) / 0.14 | 0.13 (0 %) / 0.36 |
| tilted 60 deg, 4 azimuths | 0.11 (0 %) / 0.11 | 0.15 (0 %) / 0.12 | 0.25 (0 %) / 0.09 | 0.33 (0 %) / 0.09 | 0.31 (0 %) / 0.10 | 0.23 (0 %) / 0.19 | 0.26 (0 %) / 0.18 | 0.21 (0 %) / 0.33 |
| on edge (90 deg), 4 azimuths | 27.83 (58 %) / 86.60 | 4.48 (0 %) / 2.64 | 3.26 (0 %) / 2.29 | 1.42 (0 %) / 1.35 | 0.58 (0 %) / 0.55 | 0.29 (0 %) / 0.33 | 0.21 (0 %) / 0.11 | 0.23 (0 %) / 0.21 |
| tilts 0/30/60/90 at 4 azimuths | 0.16 (5 %) / 0.08 | 0.21 (0 %) / 0.09 | 0.20 (0 %) / 0.10 | 0.22 (0 %) / 0.09 | 0.15 (0 %) / 0.09 | 0.17 (0 %) / 0.09 | 0.16 (0 %) / 0.12 | 0.27 (0 %) / 0.42 |

d1. Whip plate (two rings of 8), 0.433 GHz, V-pol source at the given elevation, 4 captures: median elevation error with the sign folded (share of answers below the horizon) / median azimuth error, deg
| plate | el 0 | el 2 | el 5 | el 10 | el 20 | el 30 | el 45 | el 70 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| flat, fixed | 3.26 (31 %) / 0.40 | 1.97 (29 %) / 0.47 | 4.12 (22 %) / 0.43 | 2.57 (36 %) / 0.42 | 1.23 (44 %) / 0.43 | 0.90 (58 %) / 0.57 | 0.78 (50 %) / 0.60 | 0.49 (53 %) / 1.36 |
| flat, 4 azimuths (one axis) | 1.16 (34 %) / 0.05 | 1.23 (34 %) / 0.04 | 0.76 (54 %) / 0.05 | 0.26 (53 %) / 0.04 | 0.11 (46 %) / 0.06 | 0.12 (52 %) / 0.06 | 0.10 (47 %) / 0.08 | 0.11 (54 %) / 0.36 |
| tilted 30 deg, fixed | 4.60 (47 %) / 2.14 | 2.52 (26 %) / 2.13 | 3.25 (20 %) / 2.10 | 3.25 (21 %) / 1.28 | 3.41 (19 %) / 1.42 | 3.81 (20 %) / 1.89 | 3.98 (31 %) / 3.19 | 2.13 (18 %) / 4.78 |
| tilted 30 deg, 4 azimuths | 0.30 (14 %) / 0.05 | 0.45 (0 %) / 0.05 | 0.46 (0 %) / 0.07 | 0.45 (0 %) / 0.11 | 0.39 (0 %) / 0.18 | 0.32 (0 %) / 0.38 | 0.32 (0 %) / 0.50 | 0.50 (0 %) / 1.29 |
| tilted 60 deg, 4 azimuths | 0.16 (3 %) / 1.44 | 0.56 (0 %) / 1.35 | 1.38 (0 %) / 1.07 | 2.42 (0 %) / 0.89 | 2.89 (0 %) / 0.59 | 2.33 (0 %) / 0.70 | 1.01 (0 %) / 0.74 | 1.23 (0 %) / 1.83 |
| on edge (90 deg), 4 azimuths | 39.82 (47 %) / 98.23 | 9.76 (23 %) / 9.38 | 7.67 (0 %) / 8.23 | 5.49 (0 %) / 6.42 | 3.15 (0 %) / 4.24 | 2.38 (0 %) / 2.73 | 1.27 (0 %) / 1.22 | 0.79 (0 %) / 1.07 |
| tilts 0/30/60/90 at 4 azimuths | 0.35 (15 %) / 0.31 | 0.39 (0 %) / 0.33 | 0.37 (0 %) / 0.34 | 0.45 (0 %) / 0.34 | 0.51 (0 %) / 0.34 | 0.70 (0 %) / 0.28 | 1.20 (0 %) / 0.36 | 1.35 (0 %) / 2.86 |

d2. Panel, 2.44 GHz, source within 20 deg of the boresight azimuth, 4 captures: median elevation error (deg) / gross
| panel | el 0 | el 15 | el 30 | el 45 | el 60 | el 75 | el 90 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| boresight on the horizon, fixed | 0.39 / 0 % | 0.38 / 0 % | 0.41 / 0 % | 0.50 / 0 % | 0.70 / 0 % | 1.54 / 0 % | 87.58 / 95 % |
| boresight 45 deg up, fixed | 0.56 / 0 % | 0.42 / 0 % | 0.37 / 0 % | 0.37 / 0 % | 0.39 / 0 % | 0.42 / 0 % | 0.75 / 0 % |
| elevation axis: 0 / 30 / 60 / 90 deg | 0.34 / 0 % | 0.35 / 0 % | 0.29 / 0 % | 0.34 / 0 % | 0.29 / 0 % | 0.35 / 0 % | 0.58 / 0 % |

d2. Panel, 5.80 GHz, source within 20 deg of the boresight azimuth, 4 captures: median elevation error (deg) / gross
| panel | el 0 | el 15 | el 30 | el 45 | el 60 | el 75 | el 90 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| boresight on the horizon, fixed | 0.16 / 0 % | 0.16 / 0 % | 0.17 / 0 % | 0.20 / 0 % | 0.32 / 0 % | 0.79 / 0 % | 83.81 / 99 % |
| boresight 45 deg up, fixed | 0.22 / 0 % | 0.18 / 0 % | 0.16 / 0 % | 0.16 / 0 % | 0.15 / 0 % | 0.17 / 0 % | 0.29 / 0 % |
| elevation axis: 0 / 30 / 60 / 90 deg | 0.17 / 0 % | 0.16 / 0 % | 0.14 / 0 % | 0.15 / 0 % | 0.15 / 0 % | 0.16 / 0 % | 0.26 / 0 % |
[section d: 20 s]

## e. Panel polarisation

e1. Fixed V-pol panel, linearly polarised source turned from vertical (ideal flare: no cross-polar response; one capture, source near boresight)
|  | 0 deg | 30 deg | 45 deg | 60 deg | 75 deg | 85 deg | 90 deg |
| --- | --- | --- | --- | --- | --- | --- | --- |
| polarisation loss (dB) = 20 log10 cos(angle) | -0.0 | 1.2 | 3.0 | 6.0 | 11.7 | 21.2 | total |
| 2.44 GHz: great-circle error median / 95th (deg) / gross | 0.63 / 1.2 / 0.0 % | 0.67 / 1.2 / 0.0 % | 0.64 / 1.3 / 0.0 % | 0.60 / 1.4 / 0.0 % | 0.78 / 1.7 / 0.0 % | 1.43 / 3.1 / 0.0 % | 78.25 / 102.7 / 96.0 % |
| 5.80 GHz: great-circle error median / 95th (deg) / gross | 0.26 / 0.5 / 0.0 % | 0.28 / 0.5 / 0.0 % | 0.28 / 0.6 / 0.0 % | 0.27 / 0.6 / 0.0 % | 0.29 / 0.7 / 0.0 % | 0.64 / 1.2 / 0.0 % | 47.76 / 97.5 / 98.0 % |

e2. Panel rolled about boresight to K angles over 180 deg, 5.8 GHz, source of random linear polarisation near boresight: polarisation-angle error and direction error
| SNR | K | roll step | source power change between captures (1 sigma) | polarisation angle error median / 95th (deg) | direction median / 95th / gross |
| --- | --- | --- | --- | --- | --- |
| 20 dB | 1 | - | 0 dB | 52.1 / 84.5 | 0.30 / 0.7 / 0.0 % |
| 20 dB | 2 | 90 deg | 0 dB | 21.2 / 50.5 | 0.25 / 0.5 / 0.0 % |
| 20 dB | 3 | 60 deg | 0 dB | 0.3 / 0.5 | 0.22 / 0.5 / 0.0 % |
| 20 dB | 4 | 45 deg | 0 dB | 0.3 / 0.6 | 0.19 / 0.4 / 0.0 % |
| 20 dB | 4 | 45 deg | 1 dB | 3.6 / 9.3 | 0.19 / 0.4 / 0.0 % |
| 20 dB | 4 | 45 deg | 3 dB | 11.1 / 26.2 | 0.20 / 0.5 / 0.0 % |
| 20 dB | 8 | 22 deg | 0 dB | 0.3 / 0.5 | 0.21 / 0.5 / 0.0 % |
| 20 dB | 8 | 22 deg | 1 dB | 2.7 / 6.6 | 0.21 / 0.5 / 0.0 % |
| 20 dB | 8 | 22 deg | 3 dB | 8.1 / 20.9 | 0.20 / 0.5 / 0.0 % |
| 0 dB | 1 | - | 0 dB | 52.0 / 83.3 | 0.88 / 5.6 / 3.0 % |
| 0 dB | 2 | 90 deg | 0 dB | 21.7 / 55.0 | 0.48 / 1.1 / 0.0 % |
| 0 dB | 3 | 60 deg | 0 dB | 1.0 / 3.0 | 0.46 / 1.2 / 0.0 % |
| 0 dB | 4 | 45 deg | 0 dB | 0.9 / 2.6 | 0.39 / 0.9 / 0.0 % |
| 0 dB | 8 | 22 deg | 0 dB | 0.7 / 2.0 | 0.31 / 0.7 / 0.0 % |
[section e: 2 s]

## f. Arbitrary orientation paths, sources anywhere on the sphere, random linear polarisation
  plate_1.575 L = 180: search scores (lower is better) evolved 0.83, random 1.17, tumble3 0.83  (40 evaluations each, 9 s)

f1. plate, 1.575 GHz, perfect orientation knowledge, 24 captures spread over the path: great-circle error median / 95th (deg) / gross / median polarisation-angle error (deg); 45 deg = polarisation not recovered
| path \ length L (time at 30 deg/s) | 180 deg (6 s) | 720 deg (24 s) |
| --- | --- | --- |
| fixed | 0.63 / 101.2 / 33 % / 56.2 | 1.18 / 102.2 / 39 % / 59.7 |
| one axis: azimuth | 0.69 / 105.8 / 33 % / 55.3 | 0.15 / 111.9 / 38 % / 56.2 |
| one axis: tilt 0-90 deg | 0.36 / 1.9 / 0 % / 3.7 | 0.44 / 2.2 / 0 % / 3.9 |
| two axes: az-el raster | 0.28 / 1.3 / 0 % / 6.0 | 0.16 / 0.5 / 0 % / 0.7 |
| two axes: spiral | 0.24 / 0.9 / 0 % / 1.5 | 0.14 / 0.4 / 0 % / 0.3 |
| two axes: tumble | 0.20 / 0.6 / 0 % / 1.1 | 0.16 / 0.4 / 0 % / 0.3 |
| three axes: tumble (gyroscopic) | 0.21 / 0.5 / 0 % / 1.1 | 0.11 / 0.3 / 0 % / 0.3 |
| three axes: evolved path | 0.21 / 0.6 / 0 % / 1.1 |  |
| three axes: best random path | 0.24 / 0.6 / 0 % / 1.9 |  |

f2. plate, 1.575 GHz, L = 720 deg: imperfect orientation knowledge (per-axis angle errors unknown to the estimator; flex = per-element phase drift with travel). Same columns as f1
| path | perfect | jitter 0.2 deg | jitter 1 deg | axis zero offset 0.5 deg | flex 3 deg | jitter 1 + offset 0.5 + flex 3 |
| --- | --- | --- | --- | --- | --- | --- |
| fixed | 1.43 / 101.4 / 41 % / 59.6 | 1.25 / 101.5 / 40 % / 57.5 | 1.85 / 102.5 / 41 % / 59.0 | 1.13 / 97.7 / 32 % / 55.9 | 1.43 / 101.4 / 41 % / 59.6 | 1.18 / 100.6 / 34 % / 57.6 |
| one axis: azimuth | 0.19 / 105.9 / 39 % / 55.4 | 0.18 / 101.7 / 36 % / 54.6 | 0.59 / 109.3 / 41 % / 57.7 | 0.70 / 101.1 / 31 % / 56.2 | 2.52 / 110.6 / 43 % / 57.9 | 0.92 / 101.6 / 32 % / 57.3 |
| one axis: tilt 0-90 deg | 0.43 / 1.9 / 0 % / 3.8 | 0.49 / 1.9 / 0 % / 3.9 | 0.61 / 2.2 / 0 % / 4.0 | 0.77 / 2.6 / 0 % / 3.8 | 0.45 / 1.9 / 0 % / 3.8 | 0.95 / 2.4 / 0 % / 4.0 |
| two axes: az-el raster | 0.16 / 0.4 / 0 % / 0.6 | 0.19 / 0.5 / 0 % / 0.6 | 0.54 / 1.1 / 0 % / 0.9 | 0.42 / 1.1 / 0 % / 0.7 | 0.22 / 0.6 / 0 % / 0.7 | 0.75 / 1.3 / 0 % / 1.1 |
| two axes: spiral | 0.12 / 0.4 / 0 % / 0.3 | 0.17 / 0.5 / 0 % / 0.3 | 0.48 / 1.0 / 0 % / 0.5 | 0.39 / 1.0 / 0 % / 0.5 | 0.21 / 0.5 / 0 % / 0.3 | 0.68 / 1.3 / 0 % / 0.6 |
| two axes: tumble | 0.13 / 0.4 / 0 % / 0.4 | 0.15 / 0.5 / 0 % / 0.3 | 0.53 / 1.1 / 0 % / 0.6 | 0.41 / 1.1 / 0 % / 0.5 | 0.22 / 0.6 / 0 % / 0.4 | 0.66 / 1.4 / 0 % / 0.7 |
| three axes: tumble (gyroscopic) | 0.12 / 0.3 / 0 % / 0.4 | 0.14 / 0.3 / 0 % / 0.4 | 0.47 / 1.0 / 0 % / 0.5 | 0.36 / 1.0 / 0 % / 0.5 | 0.17 / 0.4 / 0 % / 0.3 | 0.65 / 1.4 / 0 % / 0.6 |

f3. plate, 1.575 GHz, sources 15-60 deg above the plate, any azimuth (where the fixed array already works). Same columns
| path | L = 180 deg | L = 720 deg |
| --- | --- | --- |
| fixed | 35.89 / 109.5 / 57 % / 66.7 | 37.28 / 113.4 / 56 % / 66.7 |
| one axis: azimuth | 35.24 / 113.0 / 52 % / 68.0 | 36.45 / 112.7 / 57 % / 68.5 |
| one axis: tilt 0-90 deg | 0.34 / 1.4 / 0 % / 3.1 | 0.38 / 1.4 / 0 % / 3.3 |
| two axes: az-el raster | 0.28 / 1.1 / 0 % / 5.1 | 0.15 / 0.4 / 0 % / 0.5 |
| two axes: spiral | 0.25 / 0.7 / 0 % / 1.6 | 0.13 / 0.3 / 0 % / 0.2 |
| two axes: tumble | 0.23 / 0.6 / 0 % / 1.1 | 0.14 / 0.3 / 0 % / 0.3 |
| three axes: tumble (gyroscopic) | 0.20 / 0.5 / 0 % / 1.2 | 0.11 / 0.3 / 0 % / 0.3 |
  panel_5.800 L = 180: search scores (lower is better) evolved 5.32, random 4.68, tumble3 5.41  (40 evaluations each, 10 s)

f1. panel, 5.800 GHz, perfect orientation knowledge, 24 captures spread over the path: great-circle error median / 95th (deg) / gross / median polarisation-angle error (deg); 45 deg = polarisation not recovered
| path \ length L (time at 30 deg/s) | 180 deg (6 s) | 720 deg (24 s) |
| --- | --- | --- |
| fixed | 70.65 / 158.7 / 64 % / 54.2 | 83.92 / 160.8 / 64 % / 53.8 |
| one axis: azimuth | 0.35 / 7.0 / 4 % / 47.4 | 0.31 / 0.9 / 2 % / 47.9 |
| one axis: roll about boresight | 77.34 / 164.2 / 63 % / 25.5 | 79.02 / 155.1 / 63 % / 20.3 |
| two axes: az-el raster | 0.35 / 7.0 / 4 % / 47.4 | 0.37 / 1.2 / 0 % / 3.9 |
| two axes: spiral | 0.47 / 3.1 / 2 % / 16.9 | 0.33 / 1.3 / 0 % / 2.9 |
| two axes: tumble | 0.49 / 120.3 / 11 % / 9.5 | 0.32 / 1.0 / 0 % / 3.2 |
| three axes: tumble (gyroscopic) | 0.44 / 124.7 / 13 % / 4.9 | 0.21 / 1.6 / 0 % / 1.3 |
| three axes: evolved path | 0.46 / 104.0 / 7 % / 5.9 |  |
| three axes: best random path | 0.42 / 5.9 / 5 % / 8.0 |  |

f2. panel, 5.800 GHz, L = 720 deg: imperfect orientation knowledge (per-axis angle errors unknown to the estimator; flex = per-element phase drift with travel). Same columns as f1
| path | perfect | jitter 0.2 deg | jitter 1 deg | axis zero offset 0.5 deg | flex 3 deg | jitter 1 + offset 0.5 + flex 3 |
| --- | --- | --- | --- | --- | --- | --- |
| fixed | 83.00 / 166.2 / 64 % / 56.7 | 82.81 / 166.2 / 64 % / 56.7 | 71.01 / 163.4 / 63 % / 52.8 | 83.02 / 166.2 / 64 % / 56.8 | 83.00 / 166.2 / 64 % / 56.7 | 71.02 / 163.4 / 62 % / 55.1 |
| one axis: azimuth | 0.30 / 1.2 / 1 % / 46.5 | 0.31 / 1.2 / 2 % / 46.5 | 0.58 / 1.6 / 1 % / 45.9 | 0.59 / 1.7 / 2 % / 46.5 | 0.41 / 1.3 / 1 % / 46.5 | 0.75 / 2.3 / 1 % / 46.4 |
| one axis: roll about boresight | 90.43 / 160.2 / 63 % / 19.7 | 90.41 / 160.2 / 63 % / 19.4 | 90.65 / 161.2 / 63 % / 19.9 | 90.36 / 160.2 / 63 % / 19.4 | 90.54 / 156.6 / 63 % / 19.8 | 90.58 / 161.4 / 63 % / 19.8 |
| two axes: az-el raster | 0.30 / 1.1 / 0 % / 3.6 | 0.35 / 1.1 / 0 % / 3.3 | 0.81 / 1.5 / 0 % / 3.7 | 0.58 / 1.5 / 0 % / 3.8 | 0.42 / 1.3 / 0 % / 3.7 | 0.95 / 1.9 / 0 % / 3.7 |
| two axes: spiral | 0.33 / 1.4 / 0 % / 2.8 | 0.39 / 1.4 / 0 % / 3.0 | 0.81 / 1.7 / 0 % / 3.1 | 0.61 / 1.6 / 0 % / 3.0 | 0.38 / 1.3 / 0 % / 2.9 | 0.91 / 2.2 / 0 % / 3.1 |
| two axes: tumble | 0.29 / 1.0 / 0 % / 2.8 | 0.35 / 1.0 / 0 % / 2.8 | 0.76 / 1.7 / 0 % / 4.0 | 0.54 / 1.4 / 0 % / 2.9 | 0.41 / 1.4 / 0 % / 2.8 | 0.95 / 1.9 / 0 % / 3.7 |
| three axes: tumble (gyroscopic) | 0.26 / 1.4 / 0 % / 1.4 | 0.26 / 1.3 / 0 % / 1.4 | 0.69 / 1.7 / 0 % / 1.5 | 0.54 / 1.5 / 0 % / 1.4 | 0.29 / 1.5 / 0 % / 1.4 | 0.85 / 2.0 / 0 % / 1.8 |

f3. panel, 5.800 GHz, sources within 45 deg az / 30 deg el of the home boresight (where the fixed array already works). Same columns
| path | L = 180 deg | L = 720 deg |
| --- | --- | --- |
| fixed | 0.27 / 0.6 / 0 % / 46.1 | 0.34 / 0.6 / 0 % / 46.1 |
| one axis: azimuth | 0.29 / 0.7 / 1 % / 46.1 | 0.33 / 0.6 / 0 % / 46.1 |
| one axis: roll about boresight | 0.19 / 0.5 / 0 % / 0.3 | 0.06 / 0.2 / 0 % / 0.3 |
| two axes: az-el raster | 0.29 / 0.7 / 1 % / 46.1 | 0.36 / 1.1 / 0 % / 7.1 |
| two axes: spiral | 0.47 / 4.8 / 2 % / 16.2 | 0.45 / 2.3 / 0 % / 6.8 |
| two axes: tumble | 0.40 / 2.7 / 0 % / 11.9 | 0.42 / 2.3 / 0 % / 3.9 |
| three axes: tumble (gyroscopic) | 0.40 / 1.9 / 0 % / 9.2 | 0.25 / 0.5 / 0 % / 0.8 |

f1. panel, 2.440 GHz, perfect orientation knowledge, 24 captures spread over the path: great-circle error median / 95th (deg) / gross / median polarisation-angle error (deg); 45 deg = polarisation not recovered
| path \ length L (time at 30 deg/s) | 180 deg (6 s) | 720 deg (24 s) |
| --- | --- | --- |
| fixed | 71.33 / 156.9 / 63 % / 56.6 | 51.69 / 152.4 / 63 % / 56.9 |
| one axis: azimuth | 0.69 / 2.9 / 4 % / 46.5 | 0.58 / 2.3 / 3 % / 46.5 |
| one axis: roll about boresight | 81.21 / 157.8 / 65 % / 29.3 | 74.21 / 161.6 / 64 % / 21.2 |
| two axes: az-el raster | 0.69 / 2.9 / 4 % / 46.5 | 0.92 / 3.9 / 0 % / 3.9 |
| two axes: spiral | 1.47 / 7.4 / 4 % / 17.0 | 0.91 / 3.6 / 0 % / 3.4 |
| two axes: tumble | 1.36 / 110.0 / 12 % / 10.7 | 0.90 / 3.3 / 0 % / 3.8 |
| three axes: tumble (gyroscopic) | 1.06 / 106.9 / 14 % / 5.5 | 0.65 / 4.3 / 1 % / 1.9 |
| three axes: evolved path |  |  |
| three axes: best random path |  |  |

f2. panel, 2.440 GHz, L = 720 deg: imperfect orientation knowledge (per-axis angle errors unknown to the estimator; flex = per-element phase drift with travel). Same columns as f1
| path | perfect | jitter 0.2 deg | jitter 1 deg | axis zero offset 0.5 deg | flex 3 deg | jitter 1 + offset 0.5 + flex 3 |
| --- | --- | --- | --- | --- | --- | --- |
| fixed | 60.70 / 157.2 / 64 % / 57.3 | 60.53 / 157.2 / 64 % / 57.2 | 61.23 / 159.3 / 62 % / 56.2 | 60.51 / 157.2 / 64 % / 56.7 | 60.70 / 157.2 / 64 % / 57.3 | 60.66 / 159.2 / 62 % / 56.3 |
| one axis: azimuth | 0.63 / 2.1 / 2 % / 46.5 | 0.61 / 2.0 / 2 % / 47.6 | 0.84 / 3.2 / 2 % / 46.1 | 0.77 / 3.2 / 3 % / 47.7 | 0.85 / 2.3 / 2 % / 46.4 | 1.16 / 3.5 / 2 % / 47.6 |
| one axis: roll about boresight | 81.49 / 154.9 / 65 % / 21.7 | 81.10 / 154.8 / 67 % / 21.6 | 81.01 / 155.0 / 67 % / 19.9 | 81.13 / 152.7 / 65 % / 23.1 | 81.51 / 155.4 / 66 % / 19.0 | 81.61 / 154.4 / 68 % / 20.8 |
| two axes: az-el raster | 0.69 / 3.9 / 0 % / 3.5 | 0.68 / 3.7 / 0 % / 3.3 | 1.10 / 3.9 / 0 % / 3.7 | 0.92 / 4.1 / 0 % / 3.7 | 1.05 / 4.1 / 0 % / 3.7 | 1.34 / 4.3 / 0 % / 3.5 |
| two axes: spiral | 0.85 / 4.1 / 0 % / 3.1 | 0.83 / 4.0 / 0 % / 3.3 | 1.16 / 3.8 / 0 % / 3.7 | 0.93 / 4.1 / 0 % / 3.1 | 0.97 / 3.8 / 0 % / 3.1 | 1.28 / 3.4 / 0 % / 3.5 |
| two axes: tumble | 0.86 / 2.9 / 0 % / 3.3 | 0.88 / 3.0 / 0 % / 3.5 | 1.04 / 3.4 / 0 % / 4.0 | 0.96 / 3.1 / 0 % / 3.4 | 1.01 / 3.7 / 0 % / 3.4 | 1.32 / 4.0 / 0 % / 4.0 |
| three axes: tumble (gyroscopic) | 0.68 / 4.1 / 1 % / 1.5 | 0.68 / 4.1 / 1 % / 1.5 | 0.88 / 3.9 / 1 % / 1.7 | 0.81 / 4.4 / 1 % / 1.8 | 0.82 / 4.1 / 1 % / 1.7 | 1.15 / 3.9 / 1 % / 1.8 |

f3. panel, 2.440 GHz, sources within 45 deg az / 30 deg el of the home boresight (where the fixed array already works). Same columns
| path | L = 180 deg | L = 720 deg |
| --- | --- | --- |
| fixed | 0.66 / 1.5 / 1 % / 46.1 | 0.83 / 1.4 / 0 % / 46.1 |
| one axis: azimuth | 0.63 / 1.4 / 1 % / 46.1 | 0.71 / 1.2 / 1 % / 45.6 |
| one axis: roll about boresight | 0.52 / 1.5 / 0 % / 0.5 | 0.23 / 1.2 / 0 % / 0.4 |
| two axes: az-el raster | 0.63 / 1.4 / 1 % / 46.1 | 0.75 / 3.6 / 0 % / 7.1 |
| two axes: spiral | 1.16 / 8.3 / 3 % / 16.7 | 1.09 / 6.1 / 3 % / 7.5 |
| two axes: tumble | 1.10 / 6.1 / 1 % / 12.4 | 0.95 / 6.7 / 0 % / 4.3 |
| three axes: tumble (gyroscopic) | 1.15 / 5.8 / 0 % / 9.5 | 0.57 / 1.3 / 0 % / 1.1 |
[section f: 134 s]
[total 189 s]
