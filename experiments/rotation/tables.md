## 0. Cross-checks

az-only model against evolve_layout.check() (median / 95th / gross)
| case | evolve_study.md | evolve_layout.check now | rotation_study.az_mc, K = 1 |
| --- | --- | --- | --- |
| 4 whips, r 0.25 wl | 0.93 / 2.62 / 0.0 % | 0.91 / 2.7 / 0.0 % | 0.91 / 2.7 / 0.0 % |
| 4 whips, r 0.354 wl | 0.75 / 170.7 / 10.6 % | 0.72 / 173.3 / 8.9 % | 0.72 / 173.3 / 8.9 % |
| 8 whips, r 0.5 wl | 0.33 / 0.93 / 0.0 % | 0.33 / 1.0 / 0.0 % | 0.33 / 1.0 / 0.0 % |
| 8 whips, r 225 mm, 2.3 GHz, 0-30 deg up | (30 % gross worst in band) | 0.15 / 128.2 / 26.3 % | 0.15 / 128.2 / 26.3 % |

sphere model, one orientation: great-circle error median / 95th (deg); REDESIGN section 9 quotes 0.12/0.49, 0.12/0.49, 0.10/0.46 for the panel
| case | one snapshot at SNR + 18 dB | 64 samples, covariance |
| --- | --- | --- |
| panel 2.44 GHz, source within 45 deg az / 30 deg el of boresight | 0.70 / 1.46 | 0.69 / 1.45 |
| panel 3.50 GHz, source within 45 deg az / 30 deg el of boresight | 0.49 / 1.01 | 0.47 / 1.03 |
| panel 5.80 GHz, source within 45 deg az / 30 deg el of boresight | 0.30 / 0.60 | 0.28 / 0.62 |
| plate 1.575 GHz, sources 20-70 deg up | 0.56 / 130.54 | 0.79 / 130.93 |
[section 0: 4 s]

## a. Fragile / ambiguous layouts turned about the vertical axis, K orientations combined

a1. Equal time: 64 samples in total, split over K orientations that fill one symmetry sector (step = 360 / (symmetry x K): 45/K deg for a ring of 8, 90/K deg for 4 whips or a square grid). median / 95th percentile (deg) / gross
| layout | K = 1 | K = 2 | K = 4 | K = 8 |
| --- | --- | --- | --- | --- |
| ring of 8, 225 mm, 2.3 GHz, sources 0-30 deg up | 0.15 / 128.3 / 27.3 % | 0.11 / 124.4 / 9.9 % | 0.10 / 0.4 / 4.3 % | 0.11 / 0.3 / 0.6 % |
| ring of 8, 225 mm, 1.8 GHz, sources 0-30 deg up | 0.15 / 75.9 / 10.1 % | 0.13 / 0.4 / 0.3 % | 0.13 / 0.4 / 0.0 % | 0.13 / 0.4 / 0.0 % |
| 4 whips, r 0.354 wl, sources 0-15 deg up | 0.76 / 169.4 / 7.2 % | 0.65 / 1.8 / 0.0 % | 0.62 / 1.8 / 0.0 % | 0.62 / 1.8 / 0.0 % |
| 4 whips, r 0.45 wl | 0.57 / 2.2 / 4.5 % | 0.50 / 1.4 / 0.0 % | 0.49 / 1.4 / 0.0 % | 0.49 / 1.4 / 0.0 % |
| 4 whips, r 0.60 wl | 0.46 / 74.7 / 9.1 % | 0.37 / 1.1 / 0.0 % | 0.36 / 1.1 / 0.0 % | 0.37 / 1.0 / 0.0 % |
| 4 whips, r 1.0 wl | 0.32 / 106.3 / 20.4 % | 0.23 / 0.6 / 0.0 % | 0.21 / 0.6 / 0.0 % | 0.21 / 0.6 / 0.0 % |
| 4x4 grid, 164 mm pitch at 2.3 GHz (1.26 wl) | 0.07 / 125.5 / 14.7 % | 0.05 / 0.2 / 0.0 % | 0.05 / 0.2 / 0.0 % | 0.05 / 0.2 / 0.0 % |
| ring of 8, 225 mm at 5.8 GHz (spacing 3.3 wl) | 0.05 / 124.8 / 21.8 % | 0.04 / 0.1 / 1.3 % | 0.04 / 0.1 / 0.0 % | 0.04 / 0.1 / 0.0 % |

a2. K times the time: 64 samples at every orientation, same steps
| layout | K = 1 | K = 2 | K = 4 | K = 8 |
| --- | --- | --- | --- | --- |
| ring of 8, 225 mm, 2.3 GHz, sources 0-30 deg up | 0.15 / 128.3 / 27.3 % | 0.11 / 123.4 / 9.1 % | 0.10 / 0.4 / 4.4 % | 0.10 / 0.3 / 0.1 % |
| ring of 8, 225 mm, 1.8 GHz, sources 0-30 deg up | 0.15 / 75.9 / 10.1 % | 0.13 / 0.4 / 0.2 % | 0.13 / 0.4 / 0.0 % | 0.13 / 0.4 / 0.0 % |
| 4 whips, r 0.354 wl, sources 0-15 deg up | 0.76 / 169.4 / 7.2 % | 0.61 / 1.8 / 0.0 % | 0.61 / 1.7 / 0.0 % | 0.61 / 1.7 / 0.0 % |
| 4 whips, r 0.45 wl | 0.57 / 2.2 / 4.5 % | 0.49 / 1.4 / 0.0 % | 0.47 / 1.4 / 0.0 % | 0.49 / 1.3 / 0.0 % |
| 4 whips, r 0.60 wl | 0.46 / 74.7 / 9.1 % | 0.36 / 1.1 / 0.0 % | 0.36 / 1.0 / 0.0 % | 0.36 / 1.0 / 0.0 % |
| 4 whips, r 1.0 wl | 0.32 / 106.3 / 20.4 % | 0.22 / 0.6 / 0.0 % | 0.21 / 0.6 / 0.0 % | 0.21 / 0.6 / 0.0 % |
| 4x4 grid, 164 mm pitch at 2.3 GHz (1.26 wl) | 0.07 / 125.5 / 14.7 % | 0.05 / 0.2 / 0.0 % | 0.05 / 0.2 / 0.0 % | 0.05 / 0.2 / 0.0 % |
| ring of 8, 225 mm at 5.8 GHz (spacing 3.3 wl) | 0.05 / 124.8 / 21.8 % | 0.04 / 0.1 / 1.4 % | 0.04 / 0.1 / 0.0 % | 0.04 / 0.1 / 0.0 % |

a3. Equal time, naive steps of 360 / K deg (a symmetric layout lands on itself)
| layout | K = 1 | K = 2 | K = 4 | K = 8 |
| --- | --- | --- | --- | --- |
| ring of 8, 225 mm, 2.3 GHz, sources 0-30 deg up | 0.15 / 128.3 / 27.3 % | 0.04 / 128.7 / 28.0 % | 0.03 / 128.4 / 27.7 % | 0.03 / 128.4 / 27.7 % |
| ring of 8, 225 mm, 1.8 GHz, sources 0-30 deg up | 0.15 / 75.9 / 10.1 % | 0.04 / 75.8 / 10.5 % | 0.03 / 76.0 / 10.7 % | 0.03 / 76.0 / 10.7 % |
| 4 whips, r 0.354 wl, sources 0-15 deg up | 0.76 / 169.4 / 7.2 % | 0.14 / 9.0 / 5.0 % | 0.15 / 169.4 / 5.1 % | 0.13 / 0.4 / 0.0 % |
| 4 whips, r 0.45 wl | 0.57 / 2.2 / 4.5 % | 0.11 / 0.4 / 3.3 % | 0.11 / 0.4 / 3.9 % | 0.10 / 0.3 / 0.0 % |
| 4 whips, r 0.60 wl | 0.46 / 74.7 / 9.1 % | 0.08 / 73.8 / 7.9 % | 0.09 / 73.5 / 7.5 % | 0.08 / 0.2 / 0.0 % |
| 4 whips, r 1.0 wl | 0.32 / 106.3 / 20.4 % | 0.06 / 105.5 / 18.3 % | 0.06 / 105.3 / 17.7 % | 0.05 / 0.1 / 0.0 % |
| 4x4 grid, 164 mm pitch at 2.3 GHz (1.26 wl) | 0.07 / 125.5 / 14.7 % | 0.02 / 109.8 / 14.3 % | 0.02 / 125.4 / 14.3 % | 0.02 / 0.0 / 0.0 % |
| ring of 8, 225 mm at 5.8 GHz (spacing 3.3 wl) | 0.05 / 124.8 / 21.8 % | 0.02 / 124.6 / 22.1 % | 0.02 / 124.8 / 22.0 % | 0.02 / 124.9 / 22.0 % |

a4. The elevated-source failure of the 225 mm ring is an estimator fault: the same data searched over azimuth and elevation (0-60 deg)
| layout | az-only, K = 1 | az-el search, K = 1 | az-el search, K = 4 |
| --- | --- | --- | --- |
| ring of 8, 225 mm, 2.3 GHz, sources 0-30 deg up | 0.15 / 128.3 / 27.3 % | 0.10 / 0.3 / 0.0 % | 0.10 / 0.3 / 0.0 % |
| ring of 8, 225 mm, 1.8 GHz, sources 0-30 deg up | 0.15 / 75.9 / 10.1 % | 0.13 / 0.4 / 0.0 % | 0.12 / 0.4 / 0.0 % |
[section a: 47 s]

## b. Orientation error and cable-flex phase drift

b. ring of 8, 225 mm, 2.3 GHz, sources 0-30 deg up: 64 samples in total; median / 95th (deg) / gross. The first orientation is the surveyed home (exact, no drift); errors apply to the moves
| condition | K = 1 | K = 2 | K = 4 | K = 8 |
| --- | --- | --- | --- | --- |
| perfect | 0.15 / 129.4 / 27.9 % | 0.11 / 124.0 / 9.5 % | 0.11 / 0.5 / 4.7 % | 0.10 / 0.3 / 0.4 % |
| angle 0.1 deg | 0.15 / 129.4 / 27.9 % | 0.12 / 124.0 / 9.6 % | 0.12 / 0.6 / 4.9 % | 0.11 / 0.3 / 0.3 % |
| angle 0.5 deg | 0.15 / 129.4 / 27.9 % | 0.21 / 124.3 / 9.6 % | 0.20 / 0.9 / 4.9 % | 0.16 / 0.5 / 0.3 % |
| angle 1 deg | 0.15 / 129.4 / 27.9 % | 0.39 / 124.1 / 9.7 % | 0.34 / 1.6 / 4.9 % | 0.25 / 0.8 / 0.5 % |
| flex 1 deg | 0.15 / 129.4 / 27.9 % | 0.12 / 123.2 / 9.6 % | 0.11 / 0.5 / 4.5 % | 0.11 / 0.3 / 0.5 % |
| flex 3 deg | 0.15 / 129.4 / 27.9 % | 0.12 / 123.2 / 9.6 % | 0.11 / 0.5 / 4.5 % | 0.11 / 0.3 / 0.5 % |
| flex 10 deg | 0.15 / 129.4 / 27.9 % | 0.13 / 123.2 / 9.6 % | 0.11 / 0.5 / 4.4 % | 0.12 / 0.3 / 0.5 % |
| angle 0.5 deg + flex 3 deg | 0.15 / 129.4 / 27.9 % | 0.21 / 123.5 / 9.5 % | 0.19 / 0.8 / 4.3 % | 0.16 / 0.5 / 0.8 % |

b. 4 whips, r 0.45 wl: 64 samples in total; median / 95th (deg) / gross. The first orientation is the surveyed home (exact, no drift); errors apply to the moves
| condition | K = 1 | K = 2 | K = 4 | K = 8 |
| --- | --- | --- | --- | --- |
| perfect | 0.54 / 1.9 / 3.8 % | 0.49 / 1.4 / 0.0 % | 0.48 / 1.4 / 0.0 % | 0.48 / 1.3 / 0.0 % |
| angle 0.1 deg | 0.54 / 1.9 / 3.8 % | 0.49 / 1.4 / 0.0 % | 0.47 / 1.4 / 0.0 % | 0.48 / 1.4 / 0.0 % |
| angle 0.5 deg | 0.54 / 1.9 / 3.8 % | 0.53 / 1.5 / 0.0 % | 0.51 / 1.4 / 0.0 % | 0.50 / 1.4 / 0.0 % |
| angle 1 deg | 0.54 / 1.9 / 3.8 % | 0.62 / 1.7 / 0.0 % | 0.57 / 1.6 / 0.0 % | 0.54 / 1.5 / 0.0 % |
| flex 1 deg | 0.54 / 1.9 / 3.8 % | 0.49 / 1.4 / 0.0 % | 0.49 / 1.4 / 0.0 % | 0.49 / 1.4 / 0.0 % |
| flex 3 deg | 0.54 / 1.9 / 3.8 % | 0.52 / 1.4 / 0.0 % | 0.50 / 1.4 / 0.0 % | 0.49 / 1.4 / 0.0 % |
| flex 10 deg | 0.54 / 1.9 / 3.8 % | 0.70 / 1.9 / 0.0 % | 0.64 / 1.9 / 0.0 % | 0.60 / 1.7 / 0.0 % |
| angle 0.5 deg + flex 3 deg | 0.54 / 1.9 / 3.8 % | 0.54 / 1.6 / 0.0 % | 0.54 / 1.5 / 0.0 % | 0.51 / 1.4 / 0.0 % |

b. 4x4 grid, 164 mm pitch at 2.3 GHz (1.26 wl): 64 samples in total; median / 95th (deg) / gross. The first orientation is the surveyed home (exact, no drift); errors apply to the moves
| condition | K = 1 | K = 2 | K = 4 | K = 8 |
| --- | --- | --- | --- | --- |
| perfect | 0.07 / 108.7 / 14.9 % | 0.06 / 0.2 / 0.0 % | 0.05 / 0.2 / 0.0 % | 0.06 / 0.2 / 0.0 % |
| angle 0.1 deg | 0.07 / 108.7 / 14.9 % | 0.07 / 0.2 / 0.0 % | 0.06 / 0.2 / 0.0 % | 0.06 / 0.2 / 0.0 % |
| angle 0.5 deg | 0.07 / 108.7 / 14.9 % | 0.17 / 0.5 / 0.0 % | 0.16 / 0.5 / 0.0 % | 0.13 / 0.4 / 0.0 % |
| angle 1 deg | 0.07 / 108.7 / 14.9 % | 0.34 / 1.0 / 0.0 % | 0.29 / 0.9 / 0.0 % | 0.24 / 0.7 / 0.0 % |
| flex 1 deg | 0.07 / 108.7 / 14.9 % | 0.06 / 0.2 / 0.0 % | 0.05 / 0.2 / 0.0 % | 0.06 / 0.2 / 0.0 % |
| flex 3 deg | 0.07 / 108.7 / 14.9 % | 0.06 / 0.2 / 0.0 % | 0.06 / 0.2 / 0.0 % | 0.06 / 0.2 / 0.0 % |
| flex 10 deg | 0.07 / 108.7 / 14.9 % | 0.07 / 0.2 / 0.0 % | 0.07 / 0.2 / 0.0 % | 0.07 / 0.2 / 0.0 % |
| angle 0.5 deg + flex 3 deg | 0.07 / 108.7 / 14.9 % | 0.18 / 0.5 / 0.0 % | 0.16 / 0.4 / 0.0 % | 0.13 / 0.4 / 0.0 % |

b. two rings of 8 (225 / 112 mm), 1.575 GHz, sources 0-15 deg up: a layout with nothing to fix
| condition | K = 1 | K = 2 | K = 4 | K = 8 |
| --- | --- | --- | --- | --- |
| perfect | 0.13 / 0.4 / 0.0 % | 0.13 / 0.4 / 0.0 % | 0.13 / 0.4 / 0.0 % | 0.13 / 0.4 / 0.0 % |
| angle 0.1 deg | 0.13 / 0.4 / 0.0 % | 0.13 / 0.4 / 0.0 % | 0.13 / 0.4 / 0.0 % | 0.13 / 0.4 / 0.0 % |
| angle 0.5 deg | 0.13 / 0.4 / 0.0 % | 0.21 / 0.6 / 0.0 % | 0.19 / 0.5 / 0.0 % | 0.17 / 0.5 / 0.0 % |
| angle 1 deg | 0.13 / 0.4 / 0.0 % | 0.38 / 1.1 / 0.0 % | 0.31 / 0.9 / 0.0 % | 0.26 / 0.8 / 0.0 % |
| flex 1 deg | 0.13 / 0.4 / 0.0 % | 0.12 / 0.4 / 0.0 % | 0.13 / 0.4 / 0.0 % | 0.13 / 0.4 / 0.0 % |
| flex 3 deg | 0.13 / 0.4 / 0.0 % | 0.13 / 0.4 / 0.0 % | 0.13 / 0.4 / 0.0 % | 0.13 / 0.4 / 0.0 % |
| flex 10 deg | 0.13 / 0.4 / 0.0 % | 0.14 / 0.4 / 0.0 % | 0.14 / 0.4 / 0.0 % | 0.13 / 0.4 / 0.0 % |
| angle 0.5 deg + flex 3 deg | 0.13 / 0.4 / 0.0 % | 0.21 / 0.6 / 0.0 % | 0.20 / 0.6 / 0.0 % | 0.17 / 0.5 / 0.0 % |
[section b: 96 s]

## c. The panel on one axis (azimuth)

c1. Fixed panel: azimuth DF and gain against the source's azimuth off boresight (source 0-15 deg up, V-pol, one capture; gain = REDESIGN boresight gain + 15 log10 cos, a model, with a -20 dB floor behind)
|  | 0 deg | 15 deg | 30 deg | 45 deg | 60 deg | 75 deg | 90 deg | 120 deg | 180 deg |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2.44 GHz: az error median / 95th (deg) / gross | 0.38 / 1.0 / 0.0 % | 0.35 / 1.1 / 0.0 % | 0.41 / 1.1 / 0.0 % | 0.54 / 1.5 / 0.0 % | 0.79 / 2.1 / 0.0 % | 1.58 / 4.7 / 0.7 % | 11.59 / 179.3 / 66.7 % | 121.66 / 164.5 / 100.0 % | 96.80 / 171.8 / 100.0 % |
| 2.44 GHz: gain (dBi) | 6.4 | 6.2 | 5.5 | 4.1 | 1.9 | -2.4 | <= -13.6 | <= -13.6 | <= -13.6 |
| 3.50 GHz: az error median / 95th (deg) / gross | 0.26 / 0.7 / 0.0 % | 0.24 / 0.7 / 0.0 % | 0.27 / 0.8 / 0.0 % | 0.36 / 1.1 / 0.0 % | 0.54 / 1.5 / 0.0 % | 1.05 / 3.6 / 0.0 % | 11.43 / 178.8 / 63.0 % | 119.33 / 171.7 / 100.0 % | 127.25 / 173.2 / 100.0 % |
| 3.50 GHz: gain (dBi) | 10.7 | 10.5 | 9.8 | 8.4 | 6.2 | 1.9 | <= -9.3 | <= -9.3 | <= -9.3 |
| 5.80 GHz: az error median / 95th (deg) / gross | 0.16 / 0.4 / 0.0 % | 0.14 / 0.4 / 0.0 % | 0.16 / 0.5 / 0.0 % | 0.23 / 0.6 / 0.0 % | 0.29 / 0.8 / 0.0 % | 0.63 / 2.3 / 0.0 % | 13.56 / 178.7 / 57.0 % | 119.40 / 171.6 / 100.0 % | 144.30 / 175.6 / 100.0 % |
| 5.80 GHz: gain (dBi) | 15.5 | 15.3 | 14.6 | 13.2 | 11.0 | 6.7 | <= -4.5 | <= -4.5 | <= -4.5 |

c2. Panel stepped to K azimuths (one capture each), sources anywhere in azimuth, 0-15 deg up
| band | K | step | az error median / 95th / gross | azimuths found within 1 deg | worst-azimuth gain |
| --- | --- | --- | --- | --- | --- |
| 2.44 GHz | 1 | 360 deg | 7.23 / 165.2 / 49.8 % | 34 % | none (gaps) |
| 2.44 GHz | 2 | 180 deg | 0.61 / 4.2 / 0.8 % | 67 % | none (gaps) |
| 2.44 GHz | 3 | 120 deg | 0.29 / 1.0 / 0.0 % | 95 % | 1.9 dBi |
| 2.44 GHz | 4 | 90 deg | 0.31 / 0.9 / 0.0 % | 98 % | 4.1 dBi |
| 2.44 GHz | 6 | 60 deg | 0.29 / 0.8 / 0.0 % | 98 % | 5.5 dBi |
| 2.44 GHz | 8 | 45 deg | 0.29 / 0.9 / 0.0 % | 98 % | 5.9 dBi |
| 5.80 GHz | 1 | 360 deg | 4.24 / 170.0 / 48.8 % | 44 % | none (gaps) |
| 5.80 GHz | 2 | 180 deg | 0.26 / 2.1 / 0.5 % | 88 % | none (gaps) |
| 5.80 GHz | 3 | 120 deg | 0.17 / 0.5 / 0.0 % | 100 % | 11.0 dBi |
| 5.80 GHz | 4 | 90 deg | 0.17 / 0.5 / 0.0 % | 100 % | 13.2 dBi |
| 5.80 GHz | 6 | 60 deg | 0.16 / 0.4 / 0.0 % | 100 % | 14.6 dBi |
| 5.80 GHz | 8 | 45 deg | 0.16 / 0.5 / 0.0 % | 100 % | 15.0 dBi |

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
[section c: 12 s]

## d. Two axes: elevation near the horizon

d1. Whip plate (two rings of 8), 1.575 GHz, V-pol source at the given elevation, 4 captures: median elevation error with the sign folded (share of answers below the horizon) / median azimuth error, deg
| plate | el 0 | el 2 | el 5 | el 10 | el 20 | el 30 | el 45 | el 70 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| flat, fixed | 2.15 (39 %) / 0.13 | 1.90 (28 %) / 0.12 | 1.36 (45 %) / 0.12 | 0.65 (55 %) / 0.13 | 0.35 (50 %) / 0.14 | 0.24 (53 %) / 0.14 | 0.18 (53 %) / 0.19 | 0.13 (53 %) / 0.35 |
| flat, 4 azimuths (one axis) | 1.12 (34 %) / 0.02 | 0.47 (46 %) / 0.01 | 0.21 (50 %) / 0.02 | 0.08 (50 %) / 0.02 | 0.04 (50 %) / 0.02 | 0.03 (50 %) / 0.02 | 0.02 (51 %) / 0.03 | 0.03 (53 %) / 0.09 |
| tilted 30 deg, fixed | 2.37 (31 %) / 1.03 | 2.50 (29 %) / 1.99 | 1.77 (20 %) / 0.98 | 1.71 (23 %) / 1.06 | 3.14 (25 %) / 1.34 | 3.34 (25 %) / 1.01 | 0.97 (25 %) / 0.73 | 20.26 (20 %) / 5.90 |
| tilted 30 deg, 4 azimuths | 0.25 (8 %) / 0.02 | 0.23 (0 %) / 0.02 | 0.24 (0 %) / 0.02 | 0.20 (0 %) / 0.03 | 0.20 (0 %) / 0.06 | 0.12 (0 %) / 0.09 | 0.06 (0 %) / 0.14 | 0.13 (0 %) / 0.32 |
| tilted 60 deg, 4 azimuths | 0.14 (1 %) / 0.10 | 0.14 (0 %) / 0.11 | 0.24 (0 %) / 0.09 | 0.33 (0 %) / 0.10 | 0.29 (0 %) / 0.14 | 0.24 (0 %) / 0.15 | 0.21 (0 %) / 0.20 | 0.22 (0 %) / 0.31 |
| on edge (90 deg), 4 azimuths | 31.18 (47 %) / 95.87 | 4.20 (0 %) / 2.69 | 3.33 (0 %) / 2.31 | 1.47 (0 %) / 1.36 | 0.46 (0 %) / 0.54 | 0.26 (0 %) / 0.25 | 0.18 (0 %) / 0.12 | 0.20 (0 %) / 0.24 |
| tilts 0/30/60/90 at 4 azimuths | 0.22 (8 %) / 0.09 | 0.21 (0 %) / 0.09 | 0.20 (0 %) / 0.09 | 0.19 (0 %) / 0.09 | 0.18 (0 %) / 0.09 | 0.16 (0 %) / 0.09 | 0.16 (0 %) / 0.13 | 0.27 (0 %) / 0.45 |

d1. Whip plate (two rings of 8), 0.433 GHz, V-pol source at the given elevation, 4 captures: median elevation error with the sign folded (share of answers below the horizon) / median azimuth error, deg
| plate | el 0 | el 2 | el 5 | el 10 | el 20 | el 30 | el 45 | el 70 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| flat, fixed | 3.07 (33 %) / 0.48 | 1.98 (27 %) / 0.46 | 4.63 (31 %) / 0.45 | 2.66 (38 %) / 0.47 | 1.23 (50 %) / 0.48 | 0.91 (52 %) / 0.53 | 0.67 (54 %) / 0.66 | 0.50 (55 %) / 1.23 |
| flat, 4 azimuths (one axis) | 1.60 (36 %) / 0.05 | 1.13 (35 %) / 0.05 | 0.60 (48 %) / 0.05 | 0.25 (51 %) / 0.05 | 0.16 (50 %) / 0.06 | 0.11 (51 %) / 0.07 | 0.09 (52 %) / 0.09 | 0.10 (52 %) / 0.31 |
| tilted 30 deg, fixed | 5.43 (42 %) / 2.68 | 3.91 (34 %) / 3.24 | 3.98 (25 %) / 2.46 | 4.20 (27 %) / 2.28 | 5.03 (23 %) / 2.14 | 3.92 (24 %) / 1.79 | 2.81 (21 %) / 1.82 | 12.35 (20 %) / 6.08 |
| tilted 30 deg, 4 azimuths | 0.48 (23 %) / 0.05 | 0.50 (0 %) / 0.05 | 0.49 (0 %) / 0.07 | 0.40 (0 %) / 0.12 | 0.38 (0 %) / 0.21 | 0.35 (0 %) / 0.31 | 0.34 (0 %) / 0.51 | 0.51 (0 %) / 1.20 |
| tilted 60 deg, 4 azimuths | 0.24 (8 %) / 1.34 | 0.55 (0 %) / 1.38 | 1.30 (0 %) / 1.07 | 2.32 (0 %) / 0.82 | 2.80 (0 %) / 0.60 | 2.18 (0 %) / 0.58 | 0.86 (0 %) / 0.76 | 1.23 (0 %) / 1.82 |
| on edge (90 deg), 4 azimuths | 37.25 (52 %) / 90.67 | 9.18 (20 %) / 8.83 | 8.02 (0 %) / 8.63 | 5.66 (0 %) / 6.58 | 2.85 (0 %) / 4.21 | 1.86 (0 %) / 2.65 | 1.20 (0 %) / 1.14 | 0.77 (0 %) / 1.14 |
| tilts 0/30/60/90 at 4 azimuths | 0.49 (24 %) / 0.33 | 0.37 (0 %) / 0.32 | 0.40 (0 %) / 0.35 | 0.41 (0 %) / 0.35 | 0.61 (0 %) / 0.31 | 0.64 (0 %) / 0.31 | 1.17 (0 %) / 0.50 | 1.38 (0 %) / 2.60 |

d2. Panel, 2.44 GHz, source within 20 deg of the boresight azimuth, 4 captures: median elevation error (deg) / gross
| panel | el 0 | el 15 | el 30 | el 45 | el 60 | el 75 | el 90 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| boresight on the horizon, fixed | 0.32 / 0 % | 0.37 / 0 % | 0.37 / 0 % | 0.48 / 0 % | 0.64 / 0 % | 1.84 / 1 % | 86.09 / 96 % |
| boresight 45 deg up, fixed | 0.48 / 0 % | 0.40 / 0 % | 0.33 / 0 % | 0.37 / 0 % | 0.33 / 0 % | 0.44 / 0 % | 0.76 / 0 % |
| elevation axis: 0 / 30 / 60 / 90 deg | 0.31 / 0 % | 0.31 / 0 % | 0.26 / 0 % | 0.31 / 0 % | 0.34 / 0 % | 0.32 / 0 % | 0.58 / 0 % |

d2. Panel, 5.80 GHz, source within 20 deg of the boresight azimuth, 4 captures: median elevation error (deg) / gross
| panel | el 0 | el 15 | el 30 | el 45 | el 60 | el 75 | el 90 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| boresight on the horizon, fixed | 0.14 / 0 % | 0.15 / 0 % | 0.16 / 0 % | 0.21 / 0 % | 0.30 / 0 % | 0.85 / 0 % | 93.34 / 99 % |
| boresight 45 deg up, fixed | 0.21 / 0 % | 0.17 / 0 % | 0.14 / 0 % | 0.15 / 0 % | 0.14 / 0 % | 0.19 / 0 % | 0.33 / 0 % |
| elevation axis: 0 / 30 / 60 / 90 deg | 0.14 / 0 % | 0.15 / 0 % | 0.13 / 0 % | 0.15 / 0 % | 0.15 / 0 % | 0.17 / 0 % | 0.26 / 0 % |
[section d: 57 s]

## e. Panel polarisation

e1. Fixed V-pol panel, linearly polarised source turned from vertical (ideal flare: no cross-polar response; one capture, source near boresight)
|  | 0 deg | 30 deg | 45 deg | 60 deg | 75 deg | 85 deg | 90 deg |
| --- | --- | --- | --- | --- | --- | --- | --- |
| polarisation loss (dB) = 20 log10 cos(angle) | -0.0 | 1.2 | 3.0 | 6.0 | 11.7 | 21.2 | total |
| 2.44 GHz: great-circle error median / 95th (deg) / gross | 0.64 / 1.3 / 0.0 % | 0.60 / 1.3 / 0.0 % | 0.68 / 1.3 / 0.0 % | 0.67 / 1.4 / 0.0 % | 0.77 / 1.6 / 0.0 % | 1.50 / 2.9 / 0.0 % | 76.78 / 103.8 / 98.0 % |
| 5.80 GHz: great-circle error median / 95th (deg) / gross | 0.27 / 0.5 / 0.0 % | 0.26 / 0.6 / 0.0 % | 0.28 / 0.5 / 0.0 % | 0.29 / 0.6 / 0.0 % | 0.33 / 0.7 / 0.0 % | 0.62 / 1.2 / 0.0 % | 51.36 / 97.4 / 97.7 % |

e2. Panel rolled about boresight to K angles over 180 deg, 5.8 GHz, source of random linear polarisation near boresight: polarisation-angle error and direction error
| SNR | K | roll step | source power change between captures (1 sigma) | polarisation angle error median / 95th (deg) | direction median / 95th / gross |
| --- | --- | --- | --- | --- | --- |
| 20 dB | 1 | - | 0 dB | 46.4 / 85.1 | 0.31 / 0.8 / 0.3 % |
| 20 dB | 2 | 90 deg | 0 dB | 22.1 / 45.2 | 0.24 / 0.5 / 0.0 % |
| 20 dB | 3 | 60 deg | 0 dB | 0.3 / 0.6 | 0.24 / 0.5 / 0.0 % |
| 20 dB | 4 | 45 deg | 0 dB | 0.3 / 0.5 | 0.21 / 0.5 / 0.0 % |
| 20 dB | 4 | 45 deg | 1 dB | 3.3 / 9.5 | 0.22 / 0.5 / 0.0 % |
| 20 dB | 4 | 45 deg | 3 dB | 10.0 / 28.2 | 0.22 / 0.5 / 0.0 % |
| 20 dB | 8 | 22 deg | 0 dB | 0.3 / 0.5 | 0.18 / 0.5 / 0.0 % |
| 20 dB | 8 | 22 deg | 1 dB | 2.7 / 7.1 | 0.19 / 0.5 / 0.0 % |
| 20 dB | 8 | 22 deg | 3 dB | 7.5 / 21.5 | 0.19 / 0.5 / 0.0 % |
| 0 dB | 1 | - | 0 dB | 46.3 / 82.7 | 0.82 / 10.6 / 5.7 % |
| 0 dB | 2 | 90 deg | 0 dB | 22.2 / 45.4 | 0.53 / 1.1 / 0.0 % |
| 0 dB | 3 | 60 deg | 0 dB | 1.1 / 3.2 | 0.49 / 1.0 / 0.0 % |
| 0 dB | 4 | 45 deg | 0 dB | 1.0 / 2.8 | 0.42 / 0.9 / 0.0 % |
| 0 dB | 8 | 22 deg | 0 dB | 0.7 / 1.8 | 0.32 / 0.7 / 0.0 % |
[section e: 7 s]

## f. Arbitrary orientation paths, sources anywhere on the sphere, random linear polarisation
  plate_1.575 L = 180: search scores (lower is better) evolved 0.77, random 0.91, tumble3 0.85  (496 evaluations each, 95 s)
  plate_1.575 L = 720: search scores (lower is better) evolved 0.40, random 0.41, tumble3 0.42  (496 evaluations each, 99 s)

f1. plate, 1.575 GHz, perfect orientation knowledge, 24 captures spread over the path: great-circle error median / 95th (deg) / gross / median polarisation-angle error (deg); 45 deg = polarisation not recovered
| path \ length L (time at 30 deg/s) | 90 deg (3 s) | 180 deg (6 s) | 360 deg (12 s) | 720 deg (24 s) | 1440 deg (48 s) |
| --- | --- | --- | --- | --- | --- |
| fixed | 3.83 / 124.5 / 46 % / 61.6 | 8.70 / 134.3 / 49 % / 62.3 | 4.15 / 122.3 / 47 % / 61.9 | 3.08 / 135.3 / 46 % / 61.2 | 1.40 / 117.4 / 41 % / 61.0 |
| one axis: azimuth | 2.02 / 122.4 / 45 % / 61.0 | 1.91 / 125.0 / 45 % / 60.2 | 3.78 / 126.7 / 46 % / 61.1 | 5.39 / 124.8 / 47 % / 61.3 | 4.99 / 125.3 / 47 % / 61.4 |
| one axis: tilt 0-90 deg | 0.40 / 2.3 / 0 % / 4.2 | 0.40 / 2.2 / 0 % / 3.8 | 0.41 / 2.2 / 0 % / 4.2 | 0.39 / 2.4 / 0 % / 4.0 | 0.40 / 2.2 / 0 % / 4.3 |
| two axes: az-el raster | 0.38 / 2.5 / 1 % / 11.4 | 0.31 / 1.2 / 0 % / 5.1 | 0.21 / 0.6 / 0 % / 1.5 | 0.15 / 0.5 / 0 % / 0.6 | 0.13 / 0.4 / 0 % / 0.3 |
| two axes: spiral | 0.34 / 1.6 / 0 % / 3.1 | 0.26 / 0.9 / 0 % / 1.5 | 0.16 / 0.5 / 0 % / 0.5 | 0.13 / 0.4 / 0 % / 0.3 | 0.13 / 0.4 / 0 % / 0.3 |
| two axes: tumble | 0.35 / 1.6 / 0 % / 7.1 | 0.23 / 0.7 / 0 % / 1.1 | 0.18 / 0.5 / 0 % / 0.6 | 0.13 / 0.4 / 0 % / 0.3 | 0.12 / 0.4 / 0 % / 0.3 |
| three axes: tumble (gyroscopic) | 0.34 / 1.6 / 0 % / 7.3 | 0.24 / 0.6 / 0 % / 1.1 | 0.17 / 0.4 / 0 % / 0.6 | 0.11 / 0.3 / 0 % / 0.3 | 0.08 / 0.2 / 0 % / 0.3 |
| three axes: evolved path |  | 0.23 / 0.6 / 0 % / 1.2 |  | 0.12 / 0.3 / 0 % / 0.3 |  |
| three axes: best random path |  | 0.24 / 0.6 / 0 % / 1.3 |  | 0.12 / 0.3 / 0 % / 0.3 |  |

f2. plate, 1.575 GHz, L = 720 deg: imperfect orientation knowledge (per-axis angle errors unknown to the estimator; flex = per-element phase drift with travel). Same columns as f1
| path | perfect | jitter 0.2 deg | jitter 1 deg | axis zero offset 0.5 deg | flex 3 deg | jitter 1 + offset 0.5 + flex 3 |
| --- | --- | --- | --- | --- | --- | --- |
| fixed | 7.37 / 131.8 / 48 % / 61.0 | 9.16 / 135.3 / 49 % / 61.0 | 5.97 / 136.8 / 47 % / 60.9 | 3.45 / 124.8 / 46 % / 60.4 | 7.37 / 131.8 / 48 % / 61.0 | 2.89 / 134.9 / 45 % / 60.5 |
| one axis: azimuth | 7.66 / 136.5 / 48 % / 61.3 | 4.32 / 136.4 / 47 % / 60.8 | 2.56 / 135.1 / 46 % / 60.1 | 2.71 / 135.5 / 47 % / 60.6 | 7.15 / 135.5 / 48 % / 60.9 | 4.85 / 136.8 / 48 % / 60.7 |
| one axis: tilt 0-90 deg | 0.40 / 2.2 / 0 % / 4.0 | 0.40 / 2.3 / 0 % / 4.1 | 0.67 / 2.4 / 0 % / 4.2 | 0.73 / 2.6 / 0 % / 4.3 | 0.43 / 2.3 / 0 % / 4.1 | 0.91 / 2.5 / 0 % / 4.2 |
| two axes: az-el raster | 0.15 / 0.4 / 0 % / 0.6 | 0.19 / 0.5 / 0 % / 0.6 | 0.52 / 1.1 / 0 % / 0.9 | 0.42 / 1.1 / 0 % / 0.7 | 0.23 / 0.6 / 0 % / 0.7 | 0.70 / 1.6 / 0 % / 1.0 |
| two axes: spiral | 0.13 / 0.4 / 0 % / 0.3 | 0.17 / 0.4 / 0 % / 0.3 | 0.50 / 1.1 / 0 % / 0.5 | 0.40 / 1.0 / 0 % / 0.4 | 0.19 / 0.5 / 0 % / 0.3 | 0.63 / 1.4 / 0 % / 0.6 |
| two axes: tumble | 0.14 / 0.4 / 0 % / 0.3 | 0.18 / 0.4 / 0 % / 0.3 | 0.55 / 1.1 / 0 % / 0.6 | 0.41 / 1.0 / 0 % / 0.5 | 0.21 / 0.5 / 0 % / 0.4 | 0.67 / 1.5 / 0 % / 0.7 |
| three axes: tumble (gyroscopic) | 0.11 / 0.3 / 0 % / 0.3 | 0.15 / 0.4 / 0 % / 0.3 | 0.51 / 1.1 / 0 % / 0.6 | 0.41 / 1.0 / 0 % / 0.5 | 0.17 / 0.4 / 0 % / 0.3 | 0.66 / 1.4 / 0 % / 0.7 |
| three axes: evolved path | 0.11 / 0.3 / 0 % / 0.3 | 0.15 / 0.3 / 0 % / 0.4 | 0.49 / 1.0 / 0 % / 0.5 | 0.40 / 1.0 / 0 % / 0.5 | 0.16 / 0.3 / 0 % / 0.3 | 0.59 / 1.4 / 0 % / 0.6 |

f3. plate, 1.575 GHz, sources 15-60 deg above the plate, any azimuth (where the fixed array already works). Same columns
| path | L = 180 deg | L = 720 deg |
| --- | --- | --- |
| fixed | 34.74 / 111.4 / 52 % / 68.8 | 36.64 / 111.8 / 53 % / 67.8 |
| one axis: azimuth | 32.99 / 109.4 / 51 % / 68.6 | 38.61 / 110.9 / 55 % / 68.5 |
| one axis: tilt 0-90 deg | 0.38 / 2.0 / 0 % / 3.8 | 0.37 / 2.0 / 0 % / 3.8 |
| two axes: az-el raster | 0.32 / 1.1 / 0 % / 5.7 | 0.15 / 0.3 / 0 % / 0.7 |
| two axes: spiral | 0.24 / 0.8 / 0 % / 1.8 | 0.12 / 0.3 / 0 % / 0.2 |
| two axes: tumble | 0.22 / 0.6 / 0 % / 1.1 | 0.13 / 0.3 / 0 % / 0.3 |
| three axes: tumble (gyroscopic) | 0.20 / 0.6 / 0 % / 1.1 | 0.10 / 0.3 / 0 % / 0.3 |
| fixed, estimator told the source is above the horizon | 0.30 / 0.8 / 0 % / 43.8 | 0.29 / 0.7 / 0 % / 43.8 |
| one axis: azimuth, estimator told the source is above the horizon | 0.18 / 0.5 / 0 % / 43.9 | 0.03 / 0.2 / 0 % / 44.0 |
  panel_5.800 L = 180: search scores (lower is better) evolved 5.05, random 5.27, tumble3 6.50  (496 evaluations each, 106 s)
  panel_5.800 L = 720: search scores (lower is better) evolved 0.70, random 0.76, tumble3 1.33  (496 evaluations each, 113 s)

f1. panel, 5.800 GHz, perfect orientation knowledge, 24 captures spread over the path: great-circle error median / 95th (deg) / gross / median polarisation-angle error (deg); 45 deg = polarisation not recovered
| path \ length L (time at 30 deg/s) | 90 deg (3 s) | 180 deg (6 s) | 360 deg (12 s) | 720 deg (24 s) | 1440 deg (48 s) |
| --- | --- | --- | --- | --- | --- |
| fixed | 54.21 / 154.8 / 56 % / 54.1 | 68.23 / 150.8 / 56 % / 55.8 | 63.55 / 153.1 / 56 % / 54.7 | 70.40 / 152.2 / 56 % / 55.0 | 59.78 / 151.6 / 56 % / 53.9 |
| one axis: azimuth | 0.45 / 153.2 / 31 % / 49.9 | 0.35 / 2.5 / 2 % / 44.8 | 0.30 / 0.8 / 1 % / 44.2 | 0.30 / 0.8 / 1 % / 43.1 | 0.30 / 0.7 / 1 % / 44.1 |
| one axis: roll about boresight | 62.77 / 147.8 / 57 % / 20.3 | 68.03 / 149.5 / 56 % / 16.4 | 67.98 / 156.5 / 56 % / 19.0 | 70.48 / 153.8 / 57 % / 16.7 | 62.54 / 156.2 / 56 % / 16.2 |
| two axes: az-el raster | 0.45 / 153.2 / 31 % / 49.9 | 0.35 / 2.5 / 2 % / 44.8 | 0.30 / 0.8 / 1 % / 44.2 | 0.33 / 1.2 / 0 % / 4.4 | 0.31 / 1.3 / 0 % / 1.3 |
| two axes: spiral | 0.74 / 144.3 / 25 % / 24.0 | 0.48 / 4.2 / 3 % / 16.3 | 0.37 / 2.8 / 0 % / 6.1 | 0.33 / 1.4 / 0 % / 2.5 | 0.32 / 1.2 / 0 % / 1.4 |
| two axes: tumble | 0.73 / 148.4 / 29 % / 24.2 | 0.54 / 130.5 / 13 % / 11.1 | 0.37 / 2.6 / 0 % / 5.7 | 0.31 / 1.8 / 0 % / 2.5 | 0.32 / 1.5 / 0 % / 1.6 |
| three axes: tumble (gyroscopic) | 0.75 / 146.3 / 30 % / 18.0 | 0.43 / 129.7 / 14 % / 5.7 | 0.34 / 1.9 / 0 % / 3.9 | 0.25 / 1.4 / 0 % / 1.2 | 0.19 / 0.4 / 0 % / 0.6 |
| three axes: evolved path |  | 0.44 / 105.1 / 7 % / 6.5 |  | 0.25 / 0.6 / 0 % / 0.7 |  |
| three axes: best random path |  | 0.45 / 5.4 / 4 % / 13.3 |  | 0.25 / 0.6 / 0 % / 1.0 |  |

f2. panel, 5.800 GHz, L = 720 deg: imperfect orientation knowledge (per-axis angle errors unknown to the estimator; flex = per-element phase drift with travel). Same columns as f1
| path | perfect | jitter 0.2 deg | jitter 1 deg | axis zero offset 0.5 deg | flex 3 deg | jitter 1 + offset 0.5 + flex 3 |
| --- | --- | --- | --- | --- | --- | --- |
| fixed | 52.73 / 155.4 / 55 % / 53.7 | 53.56 / 154.3 / 55 % / 53.4 | 48.44 / 155.5 / 55 % / 52.4 | 52.59 / 155.4 / 55 % / 53.6 | 52.73 / 155.4 / 55 % / 53.7 | 46.21 / 155.5 / 55 % / 52.2 |
| one axis: azimuth | 0.29 / 0.9 / 1 % / 44.2 | 0.30 / 0.9 / 1 % / 44.2 | 0.57 / 1.5 / 1 % / 44.1 | 0.62 / 1.4 / 0 % / 43.5 | 0.39 / 1.1 / 1 % / 44.3 | 0.83 / 2.0 / 1 % / 43.8 |
| one axis: roll about boresight | 66.79 / 152.0 / 56 % / 16.8 | 67.88 / 151.9 / 56 % / 16.9 | 66.06 / 152.4 / 57 % / 16.6 | 66.79 / 152.4 / 56 % / 16.3 | 61.52 / 151.3 / 56 % / 17.0 | 61.07 / 151.8 / 57 % / 18.2 |
| two axes: az-el raster | 0.34 / 1.2 / 0 % / 4.3 | 0.37 / 1.1 / 0 % / 4.3 | 0.69 / 1.7 / 0 % / 3.9 | 0.61 / 1.5 / 0 % / 4.4 | 0.42 / 1.2 / 0 % / 4.3 | 0.92 / 2.0 / 0 % / 4.2 |
| two axes: spiral | 0.32 / 1.5 / 0 % / 2.5 | 0.35 / 1.4 / 0 % / 2.5 | 0.73 / 1.9 / 0 % / 3.0 | 0.61 / 1.6 / 0 % / 2.6 | 0.40 / 1.5 / 0 % / 2.5 | 0.93 / 2.1 / 0 % / 3.0 |
| two axes: tumble | 0.31 / 1.6 / 0 % / 2.5 | 0.34 / 1.7 / 0 % / 2.6 | 0.71 / 2.1 / 0 % / 2.9 | 0.61 / 1.8 / 0 % / 2.6 | 0.39 / 1.7 / 0 % / 2.5 | 0.90 / 2.2 / 0 % / 2.9 |
| three axes: tumble (gyroscopic) | 0.25 / 1.4 / 0 % / 1.2 | 0.29 / 1.4 / 0 % / 1.2 | 0.69 / 2.1 / 0 % / 1.7 | 0.57 / 1.6 / 0 % / 1.5 | 0.32 / 1.4 / 0 % / 1.3 | 0.88 / 2.4 / 0 % / 2.0 |
| three axes: evolved path | 0.24 / 0.6 / 0 % / 0.8 | 0.27 / 0.6 / 0 % / 0.8 | 0.75 / 1.5 / 0 % / 1.2 | 0.53 / 1.1 / 0 % / 1.0 | 0.28 / 0.7 / 0 % / 0.8 | 0.86 / 1.9 / 0 % / 1.3 |

f3. panel, 5.800 GHz, sources within 45 deg az / 30 deg el of the home boresight (where the fixed array already works). Same columns
| path | L = 180 deg | L = 720 deg |
| --- | --- | --- |
| fixed | 0.28 / 0.6 / 0 % / 44.0 | 0.28 / 0.6 / 0 % / 44.0 |
| one axis: azimuth | 0.29 / 0.6 / 0 % / 44.0 | 0.28 / 0.6 / 0 % / 44.0 |
| one axis: roll about boresight | 0.20 / 0.5 / 0 % / 0.3 | 0.05 / 0.2 / 0 % / 0.3 |
| two axes: az-el raster | 0.29 / 0.6 / 0 % / 44.0 | 0.33 / 1.1 / 0 % / 7.1 |
| two axes: spiral | 0.54 / 4.5 / 2 % / 19.5 | 0.37 / 2.5 / 0 % / 6.4 |
| two axes: tumble | 0.40 / 3.3 / 0 % / 11.7 | 0.37 / 1.9 / 0 % / 4.8 |
| three axes: tumble (gyroscopic) | 0.39 / 2.5 / 0 % / 8.7 | 0.20 / 0.6 / 0 % / 0.7 |

f1. panel, 2.440 GHz, perfect orientation knowledge, 24 captures spread over the path: great-circle error median / 95th (deg) / gross / median polarisation-angle error (deg); 45 deg = polarisation not recovered
| path \ length L (time at 30 deg/s) | 90 deg (3 s) | 180 deg (6 s) | 360 deg (12 s) | 720 deg (24 s) | 1440 deg (48 s) |
| --- | --- | --- | --- | --- | --- |
| fixed | 54.45 / 157.3 / 56 % / 53.8 | 51.81 / 155.7 / 56 % / 55.1 | 47.52 / 145.9 / 57 % / 54.8 | 50.03 / 147.4 / 56 % / 53.6 | 47.33 / 146.9 / 56 % / 57.1 |
| one axis: azimuth | 0.85 / 154.0 / 31 % / 51.8 | 0.62 / 4.2 / 2 % / 44.7 | 0.65 / 2.3 / 1 % / 43.7 | 0.64 / 2.1 / 1 % / 43.1 | 0.64 / 1.6 / 1 % / 43.2 |
| one axis: roll about boresight | 48.12 / 158.4 / 59 % / 20.4 | 59.01 / 155.2 / 61 % / 17.0 | 63.02 / 155.6 / 62 % / 16.1 | 53.37 / 159.9 / 61 % / 17.1 | 62.49 / 158.8 / 60 % / 18.9 |
| two axes: az-el raster | 0.85 / 154.0 / 31 % / 51.8 | 0.62 / 4.2 / 2 % / 44.7 | 0.65 / 2.3 / 1 % / 43.7 | 0.87 / 4.0 / 0 % / 4.5 | 0.77 / 3.4 / 0 % / 1.4 |
| two axes: spiral | 1.94 / 143.3 / 27 % / 24.0 | 1.23 / 8.5 / 5 % / 17.8 | 1.09 / 6.1 / 1 % / 6.9 | 0.88 / 4.3 / 0 % / 2.9 | 0.80 / 3.8 / 0 % / 1.6 |
| two axes: tumble | 1.94 / 148.4 / 31 % / 25.6 | 1.43 / 132.3 / 15 % / 12.4 | 0.96 / 6.5 / 2 % / 6.3 | 0.75 / 4.7 / 1 % / 2.8 | 0.78 / 4.8 / 0 % / 1.8 |
| three axes: tumble (gyroscopic) | 2.26 / 148.9 / 31 % / 17.6 | 1.17 / 120.5 / 14 % / 6.8 | 0.91 / 6.3 / 2 % / 4.4 | 0.65 / 4.4 / 1 % / 1.6 | 0.44 / 0.9 / 0 % / 0.7 |
| three axes: evolved path |  |  |  |  |  |
| three axes: best random path |  |  |  |  |  |

f2. panel, 2.440 GHz, L = 720 deg: imperfect orientation knowledge (per-axis angle errors unknown to the estimator; flex = per-element phase drift with travel). Same columns as f1
| path | perfect | jitter 0.2 deg | jitter 1 deg | axis zero offset 0.5 deg | flex 3 deg | jitter 1 + offset 0.5 + flex 3 |
| --- | --- | --- | --- | --- | --- | --- |
| fixed | 55.65 / 144.6 / 56 % / 54.7 | 55.09 / 145.3 / 56 % / 54.6 | 50.52 / 146.0 / 56 % / 54.4 | 55.72 / 144.7 / 56 % / 54.6 | 55.65 / 144.6 / 56 % / 54.7 | 46.99 / 146.4 / 56 % / 54.0 |
| one axis: azimuth | 0.63 / 1.6 / 1 % / 44.2 | 0.62 / 1.6 / 1 % / 44.1 | 0.79 / 2.4 / 1 % / 43.4 | 0.83 / 2.4 / 1 % / 43.9 | 0.82 / 2.1 / 1 % / 44.4 | 1.16 / 2.9 / 1 % / 44.6 |
| one axis: roll about boresight | 64.91 / 159.5 / 60 % / 17.1 | 64.81 / 159.5 / 60 % / 16.8 | 62.08 / 159.4 / 61 % / 17.8 | 64.75 / 159.5 / 60 % / 17.4 | 67.09 / 158.8 / 60 % / 16.4 | 62.88 / 157.8 / 62 % / 17.4 |
| two axes: az-el raster | 0.87 / 3.8 / 0 % / 4.3 | 0.88 / 3.9 / 0 % / 4.5 | 1.11 / 4.0 / 0 % / 4.2 | 1.05 / 4.1 / 0 % / 4.6 | 1.03 / 3.9 / 0 % / 4.4 | 1.40 / 4.0 / 0 % / 4.5 |
| two axes: spiral | 0.78 / 4.4 / 0 % / 3.0 | 0.78 / 4.5 / 0 % / 3.1 | 1.04 / 4.3 / 0 % / 3.3 | 0.99 / 4.5 / 0 % / 3.2 | 0.97 / 4.5 / 0 % / 3.2 | 1.32 / 4.6 / 0 % / 3.3 |
| two axes: tumble | 0.85 / 4.5 / 1 % / 3.1 | 0.87 / 4.4 / 1 % / 3.1 | 1.07 / 4.6 / 1 % / 3.0 | 0.99 / 4.5 / 1 % / 3.0 | 0.96 / 4.9 / 1 % / 2.9 | 1.33 / 5.0 / 1 % / 3.3 |
| three axes: tumble (gyroscopic) | 0.65 / 4.3 / 2 % / 1.6 | 0.67 / 4.3 / 2 % / 1.6 | 0.96 / 4.3 / 2 % / 2.1 | 0.82 / 4.5 / 2 % / 1.7 | 0.81 / 4.3 / 2 % / 1.6 | 1.12 / 4.5 / 1 % / 2.2 |

f3. panel, 2.440 GHz, sources within 45 deg az / 30 deg el of the home boresight (where the fixed array already works). Same columns
| path | L = 180 deg | L = 720 deg |
| --- | --- | --- |
| fixed | 0.66 / 1.5 / 0 % / 44.0 | 0.67 / 1.5 / 0 % / 43.9 |
| one axis: azimuth | 0.66 / 1.3 / 0 % / 44.0 | 0.58 / 1.3 / 0 % / 44.0 |
| one axis: roll about boresight | 0.55 / 1.5 / 0 % / 0.5 | 0.24 / 1.2 / 0 % / 0.4 |
| two axes: az-el raster | 0.66 / 1.3 / 0 % / 44.0 | 0.74 / 3.8 / 0 % / 7.4 |
| two axes: spiral | 1.38 / 9.4 / 5 % / 20.8 | 0.90 / 6.7 / 1 % / 6.7 |
| two axes: tumble | 1.07 / 7.7 / 2 % / 12.6 | 0.88 / 5.4 / 0 % / 5.4 |
| three axes: tumble (gyroscopic) | 1.08 / 6.6 / 1 % / 9.3 | 0.51 / 1.6 / 0 % / 1.0 |
[section f: 955 s]
[total 1178 s]
