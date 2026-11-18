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
[total 11 s]
