<!-- Written by the `report` step of array16.lab. -->

# The verdicts

132 of 142 shots are on file.

| Integrity |  |
| --- | --- |
| shots | 142 |
| measured | 132 |
| predictions match the code | True |
| records made before registration | 0 |
| records from other predictions | 0 |
| records from uncommitted code | 0 |
| commits | ['6c058dfadd51'] |
| benches | ['local'] |

| Criterion | Passed | Judged |
| --- | --- | --- |
| C1 reach | 54 | 55 |
| C2 held or lost | 92 | 93 |
| C3 closed orbit | 21 | 22 |
| C4 turn per lap | 15 | 22 |
| C5 matched launch | 26 | 34 |
| C6 ring | 10 | 10 |
| C8 scaling | 16 | 16 |
| C10 systematics | 6 | 6 |
| C11 first run, this tracker | 24 | 24 |
| C11 first run, the theory | 23 | 23 |
| C7 a twist without a knot (C1 to C5 on the helix and the figure-eight) | 123 | 131 |
| C9 against the field (C1 to C5 on the shots sent backwards) | 82 | 83 |

## Failed: C1 reach

- trefoil-shallow matched rho=+0.5000mm B=20mT laps=20 steps=20: predicted 6.54 mm, and it was lost

## Failed: C2 held or lost

- trefoil-shallow matched rho=+0.5000mm B=20mT laps=20 steps=20: predicted held (reach 6.54 mm), tracked lost after 5.39 laps

## Failed: C3 closed orbit

- helix axis rho=+2.0000mm B=20mT laps=20 steps=20: predicted 1.31 at -0° mm, fitted 1.64 at -0° mm (0.32 mm apart)

## Failed: C4 turn per lap

- trefoil axis rho=+0.5000mm B=20mT laps=20 steps=20: predicted -163.9°, fitted -152.5°
- helix axis rho=+1.0000mm B=20mT laps=20 steps=20: predicted +123.2°, fitted +134.7°
- helix axis rho=+1.5000mm B=20mT laps=20 steps=20: predicted +123.2°, fitted +137.3°
- helix axis rho=+2.0000mm B=20mT laps=20 steps=20: predicted +123.2°, fitted +139.8°
- helix axis rho=-2.0000mm B=20mT laps=20 steps=20: predicted -123.2°, fitted -115.1°
- trefoil axis rho=+0.5000mm B=40mT laps=20 steps=20: predicted -163.9°, fitted -152.5°
- helix axis rho=+1.0000mm B=40mT laps=20 steps=20: predicted +123.2°, fitted +134.7°

## Failed: C5 matched launch

- trefoil matched rho=+1.0000mm B=20mT laps=20 steps=20: closed orbit 3.56 mm from the axis, crossings stray 1.15 mm from it
- helix matched rho=+1.5000mm B=20mT laps=20 steps=20: closed orbit 0.86 mm from the axis, crossings stray 0.37 mm from it
- helix matched rho=+2.0000mm B=20mT laps=20 steps=20: closed orbit 1.31 mm from the axis, crossings stray 0.59 mm from it
- trefoil matched rho=+1.0000mm B=40mT laps=20 steps=20: closed orbit 3.56 mm from the axis, crossings stray 1.15 mm from it
- trefoil matched rho=+1.0000mm B=20mT laps=20 steps=40: closed orbit 3.56 mm from the axis, crossings stray 1.21 mm from it
- trefoil matched rho=+1.0000mm B=20mT laps=60 steps=20: closed orbit 3.56 mm from the axis, crossings stray 1.17 mm from it
- trefoil-fine-curve matched rho=+1.0000mm B=20mT laps=20 steps=20: closed orbit 3.55 mm from the axis, crossings stray 1.16 mm from it
- trefoil-fine-winding matched rho=+1.0000mm B=20mT laps=20 steps=20: closed orbit 3.56 mm from the axis, crossings stray 1.16 mm from it

## Failed: C7 a twist without a knot (C1 to C5 on the helix and the figure-eight)

- helix axis rho=+1.0000mm B=20mT laps=20 steps=20: predicted +123.2°, fitted +134.7°
- helix axis rho=+1.5000mm B=20mT laps=20 steps=20: predicted +123.2°, fitted +137.3°
- helix axis rho=+2.0000mm B=20mT laps=20 steps=20: predicted 1.31 at -0° mm, fitted 1.64 at -0° mm (0.32 mm apart)
- helix axis rho=+2.0000mm B=20mT laps=20 steps=20: predicted +123.2°, fitted +139.8°
- helix axis rho=-2.0000mm B=20mT laps=20 steps=20: predicted -123.2°, fitted -115.1°
- helix axis rho=+1.0000mm B=40mT laps=20 steps=20: predicted +123.2°, fitted +134.7°
- helix matched rho=+1.5000mm B=20mT laps=20 steps=20: closed orbit 0.86 mm from the axis, crossings stray 0.37 mm from it
- helix matched rho=+2.0000mm B=20mT laps=20 steps=20: closed orbit 1.31 mm from the axis, crossings stray 0.59 mm from it

## Failed: C9 against the field (C1 to C5 on the shots sent backwards)

- helix axis rho=-2.0000mm B=20mT laps=20 steps=20: predicted -123.2°, fitted -115.1°

## Every comparison

| Criterion | Shot | Result | Numbers |
| --- | --- | --- | --- |
| C1 reach | trefoil axis rho=+0.5000mm B=20mT laps=20 steps=20 | pass | predicted 4.30 mm, tracked 4.51 mm |
| C1 reach | trefoil axis rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted 5.64 mm, tracked 5.46 mm |
| C1 reach | trefoil matched rho=+0.5000mm B=20mT laps=20 steps=20 | pass | predicted 2.81 mm, tracked 3.08 mm |
| C1 reach | trefoil matched rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted 2.98 mm, tracked 2.96 mm |
| C1 reach | trefoil matched rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted 5.47 mm, tracked 6.31 mm |
| C1 reach | trefoil matched rho=-1.0000mm B=20mT laps=20 steps=20 | pass | predicted 5.31 mm, tracked 5.31 mm |
| C1 reach | trefoil-shallow axis rho=+0.2500mm B=20mT laps=20 steps=20 | pass | predicted 5.36 mm, tracked 5.80 mm |
| C1 reach | trefoil-shallow matched rho=+0.2500mm B=20mT laps=20 steps=20 | pass | predicted 2.95 mm, tracked 3.41 mm |
| C1 reach | trefoil-shallow matched rho=-0.2500mm B=20mT laps=20 steps=20 | pass | predicted 4.59 mm, tracked 4.59 mm |
| C1 reach | trefoil-shallow matched rho=+0.5000mm B=20mT laps=20 steps=20 | **FAIL** | predicted 6.54 mm, and it was lost |
| C1 reach | helix axis rho=+0.5000mm B=20mT laps=20 steps=20 | pass | predicted 1.66 mm, tracked 1.63 mm |
| C1 reach | helix axis rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted 2.12 mm, tracked 2.12 mm |
| C1 reach | helix matched rho=+0.5000mm B=20mT laps=20 steps=20 | pass | predicted 1.60 mm, tracked 1.63 mm |
| C1 reach | helix matched rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted 1.15 mm, tracked 1.17 mm |
| C1 reach | helix axis rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted 3.25 mm, tracked 3.31 mm |
| C1 reach | helix axis rho=-1.0000mm B=20mT laps=20 steps=20 | pass | predicted 3.58 mm, tracked 3.62 mm |
| C1 reach | helix matched rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted 2.85 mm, tracked 2.92 mm |
| C1 reach | helix matched rho=-1.0000mm B=20mT laps=20 steps=20 | pass | predicted 2.15 mm, tracked 2.16 mm |
| C1 reach | helix axis rho=+1.5000mm B=20mT laps=20 steps=20 | pass | predicted 4.95 mm, tracked 5.07 mm |
| C1 reach | helix axis rho=-1.5000mm B=20mT laps=20 steps=20 | pass | predicted 5.27 mm, tracked 5.35 mm |
| C1 reach | helix matched rho=+1.5000mm B=20mT laps=20 steps=20 | pass | predicted 4.10 mm, tracked 4.25 mm |
| C1 reach | helix matched rho=-1.5000mm B=20mT laps=20 steps=20 | pass | predicted 3.39 mm, tracked 3.47 mm |
| C1 reach | helix axis rho=+2.0000mm B=20mT laps=20 steps=20 | pass | predicted 6.66 mm, tracked 6.83 mm |
| C1 reach | helix matched rho=+2.0000mm B=20mT laps=20 steps=20 | pass | predicted 5.35 mm, tracked 5.59 mm |
| C1 reach | helix matched rho=-2.0000mm B=20mT laps=20 steps=20 | pass | predicted 4.64 mm, tracked 4.82 mm |
| C1 reach | figure-eight axis rho=+0.5000mm B=20mT laps=20 steps=20 | pass | predicted 2.99 mm, tracked 2.91 mm |
| C1 reach | figure-eight axis rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted 2.79 mm, tracked 2.79 mm |
| C1 reach | figure-eight matched rho=+0.5000mm B=20mT laps=20 steps=20 | pass | predicted 1.62 mm, tracked 1.62 mm |
| C1 reach | figure-eight matched rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted 1.40 mm, tracked 1.40 mm |
| C1 reach | figure-eight axis rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted 5.52 mm, tracked 5.38 mm |
| C1 reach | figure-eight axis rho=-1.0000mm B=20mT laps=20 steps=20 | pass | predicted 5.57 mm, tracked 5.61 mm |
| C1 reach | figure-eight matched rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted 2.76 mm, tracked 2.76 mm |
| C1 reach | figure-eight matched rho=-1.0000mm B=20mT laps=20 steps=20 | pass | predicted 2.78 mm, tracked 2.82 mm |
| C1 reach | figure-eight matched rho=+1.5000mm B=20mT laps=20 steps=20 | pass | predicted 4.15 mm, tracked 4.15 mm |
| C1 reach | figure-eight matched rho=-1.5000mm B=20mT laps=20 steps=20 | pass | predicted 4.17 mm, tracked 4.30 mm |
| C1 reach | figure-eight matched rho=+2.0000mm B=20mT laps=20 steps=20 | pass | predicted 5.54 mm, tracked 5.54 mm |
| C1 reach | figure-eight matched rho=-2.0000mm B=20mT laps=20 steps=20 | pass | predicted 5.56 mm, tracked 5.83 mm |
| C1 reach | trefoil axis rho=+0.5000mm B=40mT laps=20 steps=20 | pass | predicted 4.30 mm, tracked 4.51 mm |
| C1 reach | trefoil axis rho=-0.5000mm B=40mT laps=20 steps=20 | pass | predicted 5.64 mm, tracked 5.46 mm |
| C1 reach | trefoil matched rho=+0.5000mm B=40mT laps=20 steps=20 | pass | predicted 2.81 mm, tracked 3.08 mm |
| C1 reach | trefoil matched rho=-0.5000mm B=40mT laps=20 steps=20 | pass | predicted 2.98 mm, tracked 2.96 mm |
| C1 reach | trefoil matched rho=+1.0000mm B=40mT laps=20 steps=20 | pass | predicted 5.47 mm, tracked 6.31 mm |
| C1 reach | trefoil matched rho=-1.0000mm B=40mT laps=20 steps=20 | pass | predicted 5.31 mm, tracked 5.31 mm |
| C1 reach | helix axis rho=+0.5000mm B=40mT laps=20 steps=20 | pass | predicted 1.66 mm, tracked 1.63 mm |
| C1 reach | helix axis rho=-0.5000mm B=40mT laps=20 steps=20 | pass | predicted 2.12 mm, tracked 2.12 mm |
| C1 reach | helix matched rho=+0.5000mm B=40mT laps=20 steps=20 | pass | predicted 1.60 mm, tracked 1.63 mm |
| C1 reach | helix matched rho=-0.5000mm B=40mT laps=20 steps=20 | pass | predicted 1.15 mm, tracked 1.17 mm |
| C1 reach | helix axis rho=+1.0000mm B=40mT laps=20 steps=20 | pass | predicted 3.25 mm, tracked 3.31 mm |
| C1 reach | helix axis rho=-1.0000mm B=40mT laps=20 steps=20 | pass | predicted 3.58 mm, tracked 3.62 mm |
| C1 reach | helix matched rho=+1.0000mm B=40mT laps=20 steps=20 | pass | predicted 2.85 mm, tracked 2.92 mm |
| C1 reach | helix matched rho=-1.0000mm B=40mT laps=20 steps=20 | pass | predicted 2.15 mm, tracked 2.16 mm |
| C1 reach | trefoil matched rho=+1.0000mm B=20mT laps=20 steps=40 | pass | predicted 5.47 mm, tracked 6.37 mm |
| C1 reach | trefoil matched rho=+1.0000mm B=20mT laps=60 steps=20 | pass | predicted 5.47 mm, tracked 6.31 mm |
| C1 reach | trefoil-fine-curve matched rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted 5.46 mm, tracked 6.30 mm |
| C1 reach | trefoil-fine-winding matched rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted 5.47 mm, tracked 6.31 mm |
| C2 held or lost | ring axis rho=+0.5000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 60.58 mm), tracked lost after 3.24 laps |
| C2 held or lost | ring axis rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 65.08 mm), tracked lost after 3.04 laps |
| C2 held or lost | ring axis rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 123.41 mm), tracked lost after 1.59 laps |
| C2 held or lost | ring axis rho=-1.0000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 127.92 mm), tracked lost after 1.54 laps |
| C2 held or lost | ring axis rho=+1.5000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 186.24 mm), tracked lost after 1.07 laps |
| C2 held or lost | ring axis rho=-1.5000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 190.75 mm), tracked lost after 1.05 laps |
| C2 held or lost | ring axis rho=+2.0000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 249.07 mm), tracked lost after 0.80 laps |
| C2 held or lost | ring axis rho=-2.0000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 253.58 mm), tracked lost after 0.79 laps |
| C2 held or lost | ring axis rho=+3.0000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 374.73 mm), tracked lost after 0.53 laps |
| C2 held or lost | ring axis rho=-3.0000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 379.24 mm), tracked lost after 0.53 laps |
| C2 held or lost | trefoil axis rho=+0.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 4.30 mm), tracked held after 19.66 laps |
| C2 held or lost | trefoil axis rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 5.64 mm), tracked held after 20.82 laps |
| C2 held or lost | trefoil matched rho=+0.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 2.81 mm), tracked held after 19.66 laps |
| C2 held or lost | trefoil matched rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 2.98 mm), tracked held after 20.84 laps |
| C2 held or lost | trefoil matched rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 5.47 mm), tracked held after 19.01 laps |
| C2 held or lost | trefoil matched rho=-1.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 5.31 mm), tracked held after 21.37 laps |
| C2 held or lost | trefoil axis rho=+1.5000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 13.75 mm), tracked lost after 0.65 laps |
| C2 held or lost | trefoil axis rho=-1.5000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 14.63 mm), tracked lost after 0.63 laps |
| C2 held or lost | trefoil matched rho=+1.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 8.12 mm), tracked held after 18.39 laps |
| C2 held or lost | trefoil matched rho=-1.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 7.81 mm), tracked held after 21.85 laps |
| C2 held or lost | trefoil axis rho=+2.0000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 18.48 mm), tracked lost after 0.41 laps |
| C2 held or lost | trefoil axis rho=-2.0000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 19.35 mm), tracked lost after 0.41 laps |
| C2 held or lost | trefoil axis rho=+3.0000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 27.94 mm), tracked lost after 0.29 laps |
| C2 held or lost | trefoil axis rho=-3.0000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 28.80 mm), tracked lost after 0.30 laps |
| C2 held or lost | trefoil-shallow axis rho=+0.2500mm B=20mT laps=20 steps=20 | pass | predicted held (reach 5.36 mm), tracked held after 19.18 laps |
| C2 held or lost | trefoil-shallow matched rho=+0.2500mm B=20mT laps=20 steps=20 | pass | predicted held (reach 2.95 mm), tracked held after 19.22 laps |
| C2 held or lost | trefoil-shallow matched rho=-0.2500mm B=20mT laps=20 steps=20 | pass | predicted held (reach 4.59 mm), tracked held after 21.50 laps |
| C2 held or lost | trefoil-shallow axis rho=+0.5000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 12.41 mm), tracked lost after 1.86 laps |
| C2 held or lost | trefoil-shallow axis rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 16.04 mm), tracked lost after 1.60 laps |
| C2 held or lost | trefoil-shallow matched rho=+0.5000mm B=20mT laps=20 steps=20 | **FAIL** | predicted held (reach 6.54 mm), tracked lost after 5.39 laps |
| C2 held or lost | trefoil-shallow matched rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 8.08 mm), tracked held after 22.55 laps |
| C2 held or lost | trefoil-shallow axis rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 26.51 mm), tracked lost after 0.82 laps |
| C2 held or lost | trefoil-shallow axis rho=-1.0000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 29.98 mm), tracked lost after 0.75 laps |
| C2 held or lost | trefoil-shallow axis rho=+1.5000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 40.60 mm), tracked lost after 0.56 laps |
| C2 held or lost | trefoil-shallow axis rho=-1.5000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 43.99 mm), tracked lost after 0.53 laps |
| C2 held or lost | trefoil-shallow axis rho=+2.0000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 54.70 mm), tracked lost after 0.40 laps |
| C2 held or lost | trefoil-shallow axis rho=-2.0000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 58.08 mm), tracked lost after 0.38 laps |
| C2 held or lost | trefoil-shallow axis rho=+3.0000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 82.89 mm), tracked lost after 0.27 laps |
| C2 held or lost | trefoil-shallow axis rho=-3.0000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 86.27 mm), tracked lost after 0.27 laps |
| C2 held or lost | helix axis rho=+0.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 1.66 mm), tracked held after 20.16 laps |
| C2 held or lost | helix axis rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 2.12 mm), tracked held after 20.18 laps |
| C2 held or lost | helix matched rho=+0.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 1.60 mm), tracked held after 20.16 laps |
| C2 held or lost | helix matched rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 1.15 mm), tracked held after 20.19 laps |
| C2 held or lost | helix axis rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 3.25 mm), tracked held after 20.13 laps |
| C2 held or lost | helix axis rho=-1.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 3.58 mm), tracked held after 20.18 laps |
| C2 held or lost | helix matched rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 2.85 mm), tracked held after 20.13 laps |
| C2 held or lost | helix matched rho=-1.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 2.15 mm), tracked held after 20.19 laps |
| C2 held or lost | helix axis rho=+1.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 4.95 mm), tracked held after 20.09 laps |
| C2 held or lost | helix axis rho=-1.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 5.27 mm), tracked held after 20.18 laps |
| C2 held or lost | helix matched rho=+1.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 4.10 mm), tracked held after 20.09 laps |
| C2 held or lost | helix matched rho=-1.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 3.39 mm), tracked held after 20.18 laps |
| C2 held or lost | helix axis rho=+2.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 6.66 mm), tracked held after 20.05 laps |
| C2 held or lost | helix axis rho=-2.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 6.98 mm), tracked held after 20.16 laps |
| C2 held or lost | helix matched rho=+2.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 5.35 mm), tracked held after 20.05 laps |
| C2 held or lost | helix matched rho=-2.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 4.64 mm), tracked held after 20.16 laps |
| C2 held or lost | helix matched rho=+3.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 7.84 mm), tracked held after 19.93 laps |
| C2 held or lost | helix matched rho=-3.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 7.14 mm), tracked held after 20.10 laps |
| C2 held or lost | figure-eight axis rho=+0.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 2.99 mm), tracked held after 20.25 laps |
| C2 held or lost | figure-eight axis rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 2.79 mm), tracked held after 20.03 laps |
| C2 held or lost | figure-eight matched rho=+0.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 1.62 mm), tracked held after 20.25 laps |
| C2 held or lost | figure-eight matched rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 1.40 mm), tracked held after 20.03 laps |
| C2 held or lost | figure-eight axis rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 5.52 mm), tracked held after 20.35 laps |
| C2 held or lost | figure-eight axis rho=-1.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 5.57 mm), tracked held after 19.90 laps |
| C2 held or lost | figure-eight matched rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 2.76 mm), tracked held after 20.36 laps |
| C2 held or lost | figure-eight matched rho=-1.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 2.78 mm), tracked held after 19.91 laps |
| C2 held or lost | figure-eight axis rho=+1.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 8.30 mm), tracked held after 20.45 laps |
| C2 held or lost | figure-eight axis rho=-1.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 8.34 mm), tracked held after 19.77 laps |
| C2 held or lost | figure-eight matched rho=+1.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 4.15 mm), tracked held after 20.45 laps |
| C2 held or lost | figure-eight matched rho=-1.5000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 4.17 mm), tracked held after 19.79 laps |
| C2 held or lost | figure-eight matched rho=+2.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 5.54 mm), tracked held after 20.54 laps |
| C2 held or lost | figure-eight matched rho=-2.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 5.56 mm), tracked held after 19.66 laps |
| C2 held or lost | figure-eight axis rho=+3.0000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 16.61 mm), tracked lost after 0.34 laps |
| C2 held or lost | figure-eight axis rho=-3.0000mm B=20mT laps=20 steps=20 | pass | predicted lost (reach 16.65 mm), tracked lost after 0.35 laps |
| C2 held or lost | figure-eight matched rho=+3.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 8.31 mm), tracked held after 20.71 laps |
| C2 held or lost | figure-eight matched rho=-3.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 8.33 mm), tracked held after 19.38 laps |
| C2 held or lost | trefoil axis rho=+0.5000mm B=40mT laps=20 steps=20 | pass | predicted held (reach 4.30 mm), tracked held after 19.66 laps |
| C2 held or lost | trefoil axis rho=-0.5000mm B=40mT laps=20 steps=20 | pass | predicted held (reach 5.64 mm), tracked held after 20.82 laps |
| C2 held or lost | trefoil matched rho=+0.5000mm B=40mT laps=20 steps=20 | pass | predicted held (reach 2.81 mm), tracked held after 19.66 laps |
| C2 held or lost | trefoil matched rho=-0.5000mm B=40mT laps=20 steps=20 | pass | predicted held (reach 2.98 mm), tracked held after 20.84 laps |
| C2 held or lost | trefoil matched rho=+1.0000mm B=40mT laps=20 steps=20 | pass | predicted held (reach 5.47 mm), tracked held after 19.01 laps |
| C2 held or lost | trefoil matched rho=-1.0000mm B=40mT laps=20 steps=20 | pass | predicted held (reach 5.31 mm), tracked held after 21.37 laps |
| C2 held or lost | helix axis rho=+0.5000mm B=40mT laps=20 steps=20 | pass | predicted held (reach 1.66 mm), tracked held after 20.16 laps |
| C2 held or lost | helix axis rho=-0.5000mm B=40mT laps=20 steps=20 | pass | predicted held (reach 2.12 mm), tracked held after 20.18 laps |
| C2 held or lost | helix matched rho=+0.5000mm B=40mT laps=20 steps=20 | pass | predicted held (reach 1.60 mm), tracked held after 20.16 laps |
| C2 held or lost | helix matched rho=-0.5000mm B=40mT laps=20 steps=20 | pass | predicted held (reach 1.15 mm), tracked held after 20.19 laps |
| C2 held or lost | helix axis rho=+1.0000mm B=40mT laps=20 steps=20 | pass | predicted held (reach 3.25 mm), tracked held after 20.13 laps |
| C2 held or lost | helix axis rho=-1.0000mm B=40mT laps=20 steps=20 | pass | predicted held (reach 3.58 mm), tracked held after 20.18 laps |
| C2 held or lost | helix matched rho=+1.0000mm B=40mT laps=20 steps=20 | pass | predicted held (reach 2.85 mm), tracked held after 20.13 laps |
| C2 held or lost | helix matched rho=-1.0000mm B=40mT laps=20 steps=20 | pass | predicted held (reach 2.15 mm), tracked held after 20.19 laps |
| C2 held or lost | trefoil matched rho=+1.0000mm B=20mT laps=20 steps=40 | pass | predicted held (reach 5.47 mm), tracked held after 19.00 laps |
| C2 held or lost | trefoil matched rho=+1.0000mm B=20mT laps=60 steps=20 | pass | predicted held (reach 5.47 mm), tracked held after 56.99 laps |
| C2 held or lost | trefoil-fine-curve matched rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 5.46 mm), tracked held after 19.01 laps |
| C2 held or lost | trefoil-fine-winding matched rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted held (reach 5.47 mm), tracked held after 19.01 laps |
| C3 closed orbit | trefoil axis rho=+0.5000mm B=20mT laps=20 steps=20 | pass | predicted 1.48 at +0° mm, fitted 1.67 at -0° mm (0.18 mm apart) |
| C4 turn per lap | trefoil axis rho=+0.5000mm B=20mT laps=20 steps=20 | **FAIL** | predicted -163.9°, fitted -152.5° |
| C3 closed orbit | trefoil axis rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted 2.67 at +180° mm, fitted 2.61 at -180° mm (0.06 mm apart) |
| C4 turn per lap | trefoil axis rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted +163.9°, fitted +161.9° |
| C3 closed orbit | trefoil-shallow axis rho=+0.2500mm B=20mT laps=20 steps=20 | pass | predicted 2.41 at +0° mm, fitted 2.71 at -0° mm (0.29 mm apart) |
| C4 turn per lap | trefoil-shallow axis rho=+0.2500mm B=20mT laps=20 steps=20 | pass | predicted -52.3°, fitted -49.2° |
| C3 closed orbit | trefoil-shallow axis rho=-0.2500mm B=20mT laps=20 steps=20 | pass | predicted 4.50 at +180° mm, fitted 4.18 at -180° mm (0.32 mm apart) |
| C4 turn per lap | trefoil-shallow axis rho=-0.2500mm B=20mT laps=20 steps=20 | pass | predicted +52.3°, fitted +53.1° |
| C3 closed orbit | helix axis rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted 0.97 at -180° mm, fitted 0.98 at -180° mm (0.01 mm apart) |
| C4 turn per lap | helix axis rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted -123.2°, fitted -125.8° |
| C3 closed orbit | helix axis rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted 0.40 at -0° mm, fitted 0.49 at -1° mm (0.09 mm apart) |
| C4 turn per lap | helix axis rho=+1.0000mm B=20mT laps=20 steps=20 | **FAIL** | predicted +123.2°, fitted +134.7° |
| C3 closed orbit | helix axis rho=-1.0000mm B=20mT laps=20 steps=20 | pass | predicted 1.42 at -180° mm, fitted 1.40 at -180° mm (0.02 mm apart) |
| C4 turn per lap | helix axis rho=-1.0000mm B=20mT laps=20 steps=20 | pass | predicted -123.2°, fitted -122.5° |
| C3 closed orbit | helix axis rho=+1.5000mm B=20mT laps=20 steps=20 | pass | predicted 0.86 at -0° mm, fitted 1.05 at -1° mm (0.19 mm apart) |
| C4 turn per lap | helix axis rho=+1.5000mm B=20mT laps=20 steps=20 | **FAIL** | predicted +123.2°, fitted +137.3° |
| C3 closed orbit | helix axis rho=-1.5000mm B=20mT laps=20 steps=20 | pass | predicted 1.88 at -180° mm, fitted 1.80 at -180° mm (0.09 mm apart) |
| C4 turn per lap | helix axis rho=-1.5000mm B=20mT laps=20 steps=20 | pass | predicted -123.2°, fitted -118.8° |
| C3 closed orbit | helix axis rho=+2.0000mm B=20mT laps=20 steps=20 | **FAIL** | predicted 1.31 at -0° mm, fitted 1.64 at -0° mm (0.32 mm apart) |
| C4 turn per lap | helix axis rho=+2.0000mm B=20mT laps=20 steps=20 | **FAIL** | predicted +123.2°, fitted +139.8° |
| C3 closed orbit | helix axis rho=-2.0000mm B=20mT laps=20 steps=20 | pass | predicted 2.34 at -180° mm, fitted 2.17 at -180° mm (0.17 mm apart) |
| C4 turn per lap | helix axis rho=-2.0000mm B=20mT laps=20 steps=20 | **FAIL** | predicted -123.2°, fitted -115.1° |
| C3 closed orbit | figure-eight axis rho=+0.5000mm B=20mT laps=20 steps=20 | pass | predicted 1.38 at +90° mm, fitted 1.35 at +90° mm (0.03 mm apart) |
| C4 turn per lap | figure-eight axis rho=+0.5000mm B=20mT laps=20 steps=20 | pass | predicted -142.4°, fitted -140.4° |
| C3 closed orbit | figure-eight axis rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted 1.40 at -90° mm, fitted 1.39 at -90° mm (0.00 mm apart) |
| C4 turn per lap | figure-eight axis rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted +142.4°, fitted +145.2° |
| C3 closed orbit | figure-eight axis rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted 2.76 at +90° mm, fitted 2.69 at +90° mm (0.07 mm apart) |
| C4 turn per lap | figure-eight axis rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted -142.4°, fitted -137.9° |
| C3 closed orbit | figure-eight axis rho=-1.0000mm B=20mT laps=20 steps=20 | pass | predicted 2.78 at -90° mm, fitted 2.80 at -90° mm (0.02 mm apart) |
| C4 turn per lap | figure-eight axis rho=-1.0000mm B=20mT laps=20 steps=20 | pass | predicted +142.4°, fitted +147.7° |
| C3 closed orbit | figure-eight axis rho=+1.5000mm B=20mT laps=20 steps=20 | pass | predicted 4.15 at +90° mm, fitted 4.00 at +90° mm (0.15 mm apart) |
| C4 turn per lap | figure-eight axis rho=+1.5000mm B=20mT laps=20 steps=20 | pass | predicted -142.4°, fitted -135.5° |
| C3 closed orbit | figure-eight axis rho=-1.5000mm B=20mT laps=20 steps=20 | pass | predicted 4.17 at -90° mm, fitted 4.24 at -90° mm (0.07 mm apart) |
| C4 turn per lap | figure-eight axis rho=-1.5000mm B=20mT laps=20 steps=20 | pass | predicted +142.4°, fitted +150.1° |
| C3 closed orbit | trefoil axis rho=+0.5000mm B=40mT laps=20 steps=20 | pass | predicted 1.48 at +0° mm, fitted 1.67 at -0° mm (0.18 mm apart) |
| C4 turn per lap | trefoil axis rho=+0.5000mm B=40mT laps=20 steps=20 | **FAIL** | predicted -163.9°, fitted -152.5° |
| C3 closed orbit | trefoil axis rho=-0.5000mm B=40mT laps=20 steps=20 | pass | predicted 2.67 at +180° mm, fitted 2.61 at -180° mm (0.06 mm apart) |
| C4 turn per lap | trefoil axis rho=-0.5000mm B=40mT laps=20 steps=20 | pass | predicted +163.9°, fitted +161.9° |
| C3 closed orbit | helix axis rho=-0.5000mm B=40mT laps=20 steps=20 | pass | predicted 0.97 at -180° mm, fitted 0.98 at -180° mm (0.01 mm apart) |
| C4 turn per lap | helix axis rho=-0.5000mm B=40mT laps=20 steps=20 | pass | predicted -123.2°, fitted -125.8° |
| C3 closed orbit | helix axis rho=+1.0000mm B=40mT laps=20 steps=20 | pass | predicted 0.40 at -0° mm, fitted 0.49 at -1° mm (0.09 mm apart) |
| C4 turn per lap | helix axis rho=+1.0000mm B=40mT laps=20 steps=20 | **FAIL** | predicted +123.2°, fitted +134.7° |
| C3 closed orbit | helix axis rho=-1.0000mm B=40mT laps=20 steps=20 | pass | predicted 1.42 at -180° mm, fitted 1.40 at -180° mm (0.02 mm apart) |
| C4 turn per lap | helix axis rho=-1.0000mm B=40mT laps=20 steps=20 | pass | predicted -123.2°, fitted -122.5° |
| C5 matched launch | trefoil matched rho=+0.5000mm B=20mT laps=20 steps=20 | pass | closed orbit 1.48 mm from the axis, crossings stray 0.34 mm from it |
| C5 matched launch | trefoil matched rho=-0.5000mm B=20mT laps=20 steps=20 | pass | closed orbit 2.67 mm from the axis, crossings stray 0.04 mm from it |
| C5 matched launch | trefoil matched rho=+1.0000mm B=20mT laps=20 steps=20 | **FAIL** | closed orbit 3.56 mm from the axis, crossings stray 1.15 mm from it |
| C5 matched launch | trefoil matched rho=-1.0000mm B=20mT laps=20 steps=20 | pass | closed orbit 4.74 mm from the axis, crossings stray 0.52 mm from it |
| C5 matched launch | trefoil-shallow matched rho=+0.2500mm B=20mT laps=20 steps=20 | pass | closed orbit 2.41 mm from the axis, crossings stray 0.51 mm from it |
| C5 matched launch | trefoil-shallow matched rho=-0.2500mm B=20mT laps=20 steps=20 | pass | closed orbit 4.50 mm from the axis, crossings stray 0.39 mm from it |
| C5 matched launch | helix matched rho=+0.5000mm B=20mT laps=20 steps=20 | pass | closed orbit 0.06 mm from the axis, crossings stray 0.06 mm from it |
| C5 matched launch | helix matched rho=-0.5000mm B=20mT laps=20 steps=20 | pass | closed orbit 0.97 mm from the axis, crossings stray 0.03 mm from it |
| C5 matched launch | helix matched rho=+1.0000mm B=20mT laps=20 steps=20 | pass | closed orbit 0.40 mm from the axis, crossings stray 0.18 mm from it |
| C5 matched launch | helix matched rho=-1.0000mm B=20mT laps=20 steps=20 | pass | closed orbit 1.42 mm from the axis, crossings stray 0.03 mm from it |
| C5 matched launch | helix matched rho=+1.5000mm B=20mT laps=20 steps=20 | **FAIL** | closed orbit 0.86 mm from the axis, crossings stray 0.37 mm from it |
| C5 matched launch | helix matched rho=-1.5000mm B=20mT laps=20 steps=20 | pass | closed orbit 1.88 mm from the axis, crossings stray 0.12 mm from it |
| C5 matched launch | helix matched rho=+2.0000mm B=20mT laps=20 steps=20 | **FAIL** | closed orbit 1.31 mm from the axis, crossings stray 0.59 mm from it |
| C5 matched launch | helix matched rho=-2.0000mm B=20mT laps=20 steps=20 | pass | closed orbit 2.34 mm from the axis, crossings stray 0.26 mm from it |
| C5 matched launch | figure-eight matched rho=+0.5000mm B=20mT laps=20 steps=20 | pass | closed orbit 1.38 mm from the axis, crossings stray 0.05 mm from it |
| C5 matched launch | figure-eight matched rho=-0.5000mm B=20mT laps=20 steps=20 | pass | closed orbit 1.40 mm from the axis, crossings stray 0.01 mm from it |
| C5 matched launch | figure-eight matched rho=+1.0000mm B=20mT laps=20 steps=20 | pass | closed orbit 2.76 mm from the axis, crossings stray 0.14 mm from it |
| C5 matched launch | figure-eight matched rho=-1.0000mm B=20mT laps=20 steps=20 | pass | closed orbit 2.78 mm from the axis, crossings stray 0.03 mm from it |
| C5 matched launch | figure-eight matched rho=+1.5000mm B=20mT laps=20 steps=20 | pass | closed orbit 4.15 mm from the axis, crossings stray 0.28 mm from it |
| C5 matched launch | figure-eight matched rho=-1.5000mm B=20mT laps=20 steps=20 | pass | closed orbit 4.17 mm from the axis, crossings stray 0.13 mm from it |
| C5 matched launch | figure-eight matched rho=+2.0000mm B=20mT laps=20 steps=20 | pass | closed orbit 5.54 mm from the axis, crossings stray 0.45 mm from it |
| C5 matched launch | figure-eight matched rho=-2.0000mm B=20mT laps=20 steps=20 | pass | closed orbit 5.56 mm from the axis, crossings stray 0.27 mm from it |
| C5 matched launch | trefoil matched rho=+0.5000mm B=40mT laps=20 steps=20 | pass | closed orbit 1.48 mm from the axis, crossings stray 0.34 mm from it |
| C5 matched launch | trefoil matched rho=-0.5000mm B=40mT laps=20 steps=20 | pass | closed orbit 2.67 mm from the axis, crossings stray 0.04 mm from it |
| C5 matched launch | trefoil matched rho=+1.0000mm B=40mT laps=20 steps=20 | **FAIL** | closed orbit 3.56 mm from the axis, crossings stray 1.15 mm from it |
| C5 matched launch | trefoil matched rho=-1.0000mm B=40mT laps=20 steps=20 | pass | closed orbit 4.74 mm from the axis, crossings stray 0.52 mm from it |
| C5 matched launch | helix matched rho=+0.5000mm B=40mT laps=20 steps=20 | pass | closed orbit 0.06 mm from the axis, crossings stray 0.06 mm from it |
| C5 matched launch | helix matched rho=-0.5000mm B=40mT laps=20 steps=20 | pass | closed orbit 0.97 mm from the axis, crossings stray 0.03 mm from it |
| C5 matched launch | helix matched rho=+1.0000mm B=40mT laps=20 steps=20 | pass | closed orbit 0.40 mm from the axis, crossings stray 0.18 mm from it |
| C5 matched launch | helix matched rho=-1.0000mm B=40mT laps=20 steps=20 | pass | closed orbit 1.42 mm from the axis, crossings stray 0.03 mm from it |
| C5 matched launch | trefoil matched rho=+1.0000mm B=20mT laps=20 steps=40 | **FAIL** | closed orbit 3.56 mm from the axis, crossings stray 1.21 mm from it |
| C5 matched launch | trefoil matched rho=+1.0000mm B=20mT laps=60 steps=20 | **FAIL** | closed orbit 3.56 mm from the axis, crossings stray 1.17 mm from it |
| C5 matched launch | trefoil-fine-curve matched rho=+1.0000mm B=20mT laps=20 steps=20 | **FAIL** | closed orbit 3.55 mm from the axis, crossings stray 1.16 mm from it |
| C5 matched launch | trefoil-fine-winding matched rho=+1.0000mm B=20mT laps=20 steps=20 | **FAIL** | closed orbit 3.56 mm from the axis, crossings stray 1.16 mm from it |
| C6 ring | ring axis rho=+0.5000mm B=20mT laps=20 steps=20 | pass | predicted lost at lap 3.27, tracked 3.24 |
| C6 ring | ring axis rho=-0.5000mm B=20mT laps=20 steps=20 | pass | predicted lost at lap 3.04, tracked 3.04 |
| C6 ring | ring axis rho=+1.0000mm B=20mT laps=20 steps=20 | pass | predicted lost at lap 1.60, tracked 1.59 |
| C6 ring | ring axis rho=-1.0000mm B=20mT laps=20 steps=20 | pass | predicted lost at lap 1.55, tracked 1.54 |
| C6 ring | ring axis rho=+1.5000mm B=20mT laps=20 steps=20 | pass | predicted lost at lap 1.06, tracked 1.07 |
| C6 ring | ring axis rho=-1.5000mm B=20mT laps=20 steps=20 | pass | predicted lost at lap 1.04, tracked 1.05 |
| C6 ring | ring axis rho=+2.0000mm B=20mT laps=20 steps=20 | pass | predicted lost at lap 0.80, tracked 0.80 |
| C6 ring | ring axis rho=-2.0000mm B=20mT laps=20 steps=20 | pass | predicted lost at lap 0.78, tracked 0.79 |
| C6 ring | ring axis rho=+3.0000mm B=20mT laps=20 steps=20 | pass | predicted lost at lap 0.53, tracked 0.53 |
| C6 ring | ring axis rho=-3.0000mm B=20mT laps=20 steps=20 | pass | predicted lost at lap 0.52, tracked 0.53 |
| C8 scaling | trefoil axis rho=+0.5000mm B=40mT laps=20 steps=20 | pass | reach 4.508 mm against 4.508 mm; turn -152.55° against -152.55° |
| C8 scaling | trefoil axis rho=-0.5000mm B=40mT laps=20 steps=20 | pass | reach 5.460 mm against 5.460 mm; turn +161.91° against +161.91° |
| C8 scaling | trefoil matched rho=+0.5000mm B=40mT laps=20 steps=20 | pass | reach 3.077 mm against 3.077 mm; turn -152.35° against -152.35° |
| C8 scaling | trefoil matched rho=-0.5000mm B=40mT laps=20 steps=20 | pass | reach 2.957 mm against 2.957 mm |
| C8 scaling | trefoil axis rho=+1.0000mm B=40mT laps=20 steps=20 | pass | reach 9.904 mm against 9.927 mm |
| C8 scaling | trefoil axis rho=-1.0000mm B=40mT laps=20 steps=20 | pass | reach 10.004 mm against 9.954 mm |
| C8 scaling | trefoil matched rho=+1.0000mm B=40mT laps=20 steps=20 | pass | reach 6.309 mm against 6.309 mm; turn -149.71° against -149.71° |
| C8 scaling | trefoil matched rho=-1.0000mm B=40mT laps=20 steps=20 | pass | reach 5.313 mm against 5.313 mm; turn +167.84° against +167.83° |
| C8 scaling | helix axis rho=+0.5000mm B=40mT laps=20 steps=20 | pass | reach 1.631 mm against 1.631 mm |
| C8 scaling | helix axis rho=-0.5000mm B=40mT laps=20 steps=20 | pass | reach 2.121 mm against 2.121 mm; turn -125.77° against -125.77° |
| C8 scaling | helix matched rho=+0.5000mm B=40mT laps=20 steps=20 | pass | reach 1.626 mm against 1.626 mm |
| C8 scaling | helix matched rho=-0.5000mm B=40mT laps=20 steps=20 | pass | reach 1.167 mm against 1.167 mm |
| C8 scaling | helix axis rho=+1.0000mm B=40mT laps=20 steps=20 | pass | reach 3.314 mm against 3.314 mm; turn +134.73° against +134.73° |
| C8 scaling | helix axis rho=-1.0000mm B=40mT laps=20 steps=20 | pass | reach 3.622 mm against 3.622 mm; turn -122.46° against -122.45° |
| C8 scaling | helix matched rho=+1.0000mm B=40mT laps=20 steps=20 | pass | reach 2.924 mm against 2.924 mm; turn +133.12° against +133.51° |
| C8 scaling | helix matched rho=-1.0000mm B=40mT laps=20 steps=20 | pass | reach 2.161 mm against 2.161 mm |
| C10 systematics | trefoil axis rho=+1.0000mm B=20mT laps=20 steps=40 | pass | reach 9.925 mm against 9.927 mm |
| C10 systematics | trefoil-fine-curve axis rho=+1.0000mm B=20mT laps=20 steps=20 | pass | reach 9.900 mm against 9.927 mm |
| C10 systematics | trefoil-fine-winding axis rho=+1.0000mm B=20mT laps=20 steps=20 | pass | reach 9.967 mm against 9.927 mm |
| C10 systematics | trefoil matched rho=+1.0000mm B=20mT laps=20 steps=40 | pass | reach 6.366 mm against 6.309 mm; turn -149.79° against -149.71° |
| C10 systematics | trefoil-fine-curve matched rho=+1.0000mm B=20mT laps=20 steps=20 | pass | reach 6.303 mm against 6.309 mm; turn -149.82° against -149.71° |
| C10 systematics | trefoil-fine-winding matched rho=+1.0000mm B=20mT laps=20 steps=20 | pass | reach 6.312 mm against 6.309 mm; turn -149.68° against -149.71° |
| C11 first run, this tracker | trefoil-single first-run rho=+6.7446mm B=5mT laps=20 steps=20 | pass | first run lost (reach 10.2 mm), now lost (reach 10.2 mm) |
| C11 first run, the theory | trefoil-single first-run rho=+6.7446mm B=5mT laps=20 steps=20 | pass | first run lost, theory lost (reach 44.7 mm) |
| C11 first run, this tracker | trefoil-single first-run rho=+11.6831mm B=5mT laps=20 steps=20 | pass | first run lost (reach 9.9 mm), now lost (reach 9.9 mm) |
| C11 first run, the theory | trefoil-single first-run rho=+11.6831mm B=5mT laps=20 steps=20 | pass | first run lost, theory lost (reach 91.4 mm) |
| C11 first run, this tracker | trefoil-single first-run rho=+21.3377mm B=5mT laps=20 steps=20 | pass | first run lost (reach 11.9 mm), now lost (reach 11.9 mm) |
| C11 first run, the theory | trefoil-single first-run rho=+21.3377mm B=5mT laps=20 steps=20 | pass | first run lost, theory lost (reach 182.7 mm) |
| C11 first run, this tracker | trefoil-single first-run rho=+36.9940mm B=5mT laps=20 steps=20 | pass | first run lost (reach 13.0 mm), now lost (reach 13.0 mm) |
| C11 first run, the theory | trefoil-single first-run rho=+36.9940mm B=5mT laps=20 steps=20 | pass | first run lost, theory lost (reach 330.7 mm) |
| C11 first run, this tracker | trefoil-single first-run rho=+3.3723mm B=10mT laps=20 steps=20 | pass | first run lost (reach 9.9 mm), now lost (reach 9.9 mm) |
| C11 first run, the theory | trefoil-single first-run rho=+3.3723mm B=10mT laps=20 steps=20 | pass | first run lost, theory lost (reach 12.8 mm) |
| C11 first run, this tracker | trefoil-single first-run rho=+5.8416mm B=10mT laps=20 steps=20 | pass | first run lost (reach 10.1 mm), now lost (reach 10.1 mm) |
| C11 first run, the theory | trefoil-single first-run rho=+5.8416mm B=10mT laps=20 steps=20 | pass | first run lost, theory lost (reach 36.1 mm) |
| C11 first run, this tracker | trefoil-single first-run rho=+10.6688mm B=10mT laps=20 steps=20 | pass | first run lost (reach 10.1 mm), now lost (reach 10.1 mm) |
| C11 first run, the theory | trefoil-single first-run rho=+10.6688mm B=10mT laps=20 steps=20 | pass | first run lost, theory lost (reach 81.8 mm) |
| C11 first run, this tracker | trefoil-single first-run rho=+18.4970mm B=10mT laps=20 steps=20 | pass | first run lost (reach 11.5 mm), now lost (reach 11.5 mm) |
| C11 first run, the theory | trefoil-single first-run rho=+18.4970mm B=10mT laps=20 steps=20 | pass | first run lost, theory lost (reach 155.8 mm) |
| C11 first run, this tracker | trefoil-single first-run rho=+1.6861mm B=20mT laps=20 steps=20 | pass | first run held (reach 5.5 mm), now held (reach 5.5 mm) |
| C11 first run, the theory | trefoil-single first-run rho=+1.6861mm B=20mT laps=20 steps=20 | pass | first run held, theory held (reach 5.6 mm) |
| C11 first run, this tracker | trefoil-single first-run rho=+2.9208mm B=20mT laps=20 steps=20 | pass | first run held (reach 8.8 mm), now held (reach 8.8 mm) |
| C11 first run, this tracker | trefoil-single first-run rho=+5.3344mm B=20mT laps=20 steps=20 | pass | first run lost (reach 10.0 mm), now lost (reach 10.0 mm) |
| C11 first run, the theory | trefoil-single first-run rho=+5.3344mm B=20mT laps=20 steps=20 | pass | first run lost, theory lost (reach 31.3 mm) |
| C11 first run, this tracker | trefoil-single first-run rho=+9.2485mm B=20mT laps=20 steps=20 | pass | first run lost (reach 10.2 mm), now lost (reach 10.2 mm) |
| C11 first run, the theory | trefoil-single first-run rho=+9.2485mm B=20mT laps=20 steps=20 | pass | first run lost, theory lost (reach 68.4 mm) |
| C11 first run, this tracker | ring-single first-run rho=+6.7446mm B=5mT laps=20 steps=20 | pass | first run lost (reach 10.0 mm), now lost (reach 10.0 mm) |
| C11 first run, the theory | ring-single first-run rho=+6.7446mm B=5mT laps=20 steps=20 | pass | first run lost, theory lost (reach 631.4 mm) |
| C11 first run, this tracker | ring-single first-run rho=+11.6831mm B=5mT laps=20 steps=20 | pass | first run lost (reach 9.9 mm), now lost (reach 9.9 mm) |
| C11 first run, the theory | ring-single first-run rho=+11.6831mm B=5mT laps=20 steps=20 | pass | first run lost, theory lost (reach 1252.0 mm) |
| C11 first run, this tracker | ring-single first-run rho=+21.3377mm B=5mT laps=20 steps=20 | pass | first run lost (reach 10.3 mm), now lost (reach 10.3 mm) |
| C11 first run, the theory | ring-single first-run rho=+21.3377mm B=5mT laps=20 steps=20 | pass | first run lost, theory lost (reach 2465.2 mm) |
| C11 first run, this tracker | ring-single first-run rho=+36.9940mm B=5mT laps=20 steps=20 | pass | first run lost (reach 13.3 mm), now lost (reach 13.3 mm) |
| C11 first run, the theory | ring-single first-run rho=+36.9940mm B=5mT laps=20 steps=20 | pass | first run lost, theory lost (reach 4432.6 mm) |
| C11 first run, this tracker | ring-single first-run rho=+3.3723mm B=10mT laps=20 steps=20 | pass | first run lost (reach 9.9 mm), now lost (reach 9.9 mm) |
| C11 first run, the theory | ring-single first-run rho=+3.3723mm B=10mT laps=20 steps=20 | pass | first run lost, theory lost (reach 207.7 mm) |
| C11 first run, this tracker | ring-single first-run rho=+5.8416mm B=10mT laps=20 steps=20 | pass | first run lost (reach 9.9 mm), now lost (reach 9.9 mm) |
| C11 first run, the theory | ring-single first-run rho=+5.8416mm B=10mT laps=20 steps=20 | pass | first run lost, theory lost (reach 518.0 mm) |
| C11 first run, this tracker | ring-single first-run rho=+10.6688mm B=10mT laps=20 steps=20 | pass | first run lost (reach 9.9 mm), now lost (reach 9.9 mm) |
| C11 first run, the theory | ring-single first-run rho=+10.6688mm B=10mT laps=20 steps=20 | pass | first run lost, theory lost (reach 1124.6 mm) |
| C11 first run, this tracker | ring-single first-run rho=+18.4970mm B=10mT laps=20 steps=20 | pass | first run lost (reach 10.1 mm), now lost (reach 10.1 mm) |
| C11 first run, the theory | ring-single first-run rho=+18.4970mm B=10mT laps=20 steps=20 | pass | first run lost, theory lost (reach 2108.3 mm) |
| C11 first run, this tracker | ring-single first-run rho=+1.6861mm B=20mT laps=20 steps=20 | pass | first run held (reach 6.5 mm), now held (reach 6.5 mm) |
| C11 first run, the theory | ring-single first-run rho=+1.6861mm B=20mT laps=20 steps=20 | pass | first run held, theory held (reach 4.7 mm) |
| C11 first run, this tracker | ring-single first-run rho=+2.9208mm B=20mT laps=20 steps=20 | pass | first run lost (reach 9.9 mm), now lost (reach 9.9 mm) |
| C11 first run, the theory | ring-single first-run rho=+2.9208mm B=20mT laps=20 steps=20 | pass | first run lost, theory lost (reach 151.0 mm) |
| C11 first run, this tracker | ring-single first-run rho=+5.3344mm B=20mT laps=20 steps=20 | pass | first run lost (reach 9.9 mm), now lost (reach 9.9 mm) |
| C11 first run, the theory | ring-single first-run rho=+5.3344mm B=20mT laps=20 steps=20 | pass | first run lost, theory lost (reach 454.2 mm) |
| C11 first run, this tracker | ring-single first-run rho=+9.2485mm B=20mT laps=20 steps=20 | pass | first run lost (reach 9.9 mm), now lost (reach 9.9 mm) |
| C11 first run, the theory | ring-single first-run rho=+9.2485mm B=20mT laps=20 steps=20 | pass | first run lost, theory lost (reach 946.1 mm) |
