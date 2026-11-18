<!-- Written by `experiments/guide_study.py predict`. Theory only: no electron was tracked to make it. -->

# Experiment 1: the predictions

## The guides

| Guide | Shape | Knotted | Winding | Length (mm) | Tightest bend (mm) | Writhe | Turn per lap | Drift per lap (per ρ) | Field-line wander per lap (mm) | Net turns along the curve | Closed orbit's offset (per ρ) | Matched reach at ρ = 1 mm (mm) | On-axis reach at ρ = 1 mm, 20 laps (mm) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ring | ring | no | return | 932 | 148.4 | +0.000 | +0.0° | 6.28 | 0.11 | +0.00 | none | none | 123.41 |
| trefoil | trefoil | yes | return | 932 | 43.8 | -3.455 | -163.9° | 8.22 | 1.17 | +0.00 | 4.15 | 5.47 | 9.02 |
| trefoil-shallow | shallow trefoil | yes | return | 814 | 51.0 | -3.146 | -52.3° | 12.19 | 0.92 | +0.00 | 13.83 | 13.72 | 26.51 |
| helix | helix | no | return | 897 | 60.1 | -0.658 | +123.2° | 1.61 | 0.90 | +0.00 | 0.91 | 2.85 | 3.25 |
| figure-eight | figure-eight | no | return | 875 | 32.4 | +0.604 | -142.4° | 5.25 | 0.02 | +0.00 | 2.77 | 2.76 | 5.52 |
| trefoil-fine-curve | trefoil | yes | return | 932 | 43.8 | -3.456 | -164.0° | 8.21 | 1.17 | +0.00 | 4.15 | 5.46 | 9.01 |
| trefoil-fine-winding | trefoil | yes | return | 932 | 43.8 | -3.455 | -163.9° | 8.22 | 1.18 | +0.00 | 4.15 | 5.47 | 9.02 |
| ring-single | ring | no | single | 932 | 148.4 | +0.000 | +0.0° | 6.28 | 10.80 | +1.00 | none | none | 90.43 |
| trefoil-single | trefoil | yes | single | 932 | 43.8 | -3.455 | -163.9° | 8.22 | 17.46 | +1.00 | 4.15 | 4.67 | 9.34 |

The wall is 9.9 mm from the curve. ρ is the gyroradius at the guide field; a minus sign
means the electron is sent against the field.

## The main shots (20 mT, 20 laps)

| Shot | Energy (eV) | Predicted | Reach (mm) | Lap of loss | Closed orbit (mm) | Launch to closed orbit (mm) | Turn per lap |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ring axis rho=+0.5000mm B=20mT laps=20 steps=20 | 8.8 | lost | 60.58 | 3.27 |  |  |  |
| ring axis rho=-0.5000mm B=20mT laps=20 steps=20 | 8.8 | lost | 65.08 | 3.04 |  |  |  |
| ring axis rho=+1.0000mm B=20mT laps=20 steps=20 | 35.2 | lost | 123.41 | 1.60 |  |  |  |
| ring axis rho=-1.0000mm B=20mT laps=20 steps=20 | 35.2 | lost | 127.92 | 1.55 |  |  |  |
| ring axis rho=+1.5000mm B=20mT laps=20 steps=20 | 79.1 | lost | 186.24 | 1.06 |  |  |  |
| ring axis rho=-1.5000mm B=20mT laps=20 steps=20 | 79.1 | lost | 190.75 | 1.04 |  |  |  |
| ring axis rho=+2.0000mm B=20mT laps=20 steps=20 | 140.7 | lost | 249.07 | 0.80 |  |  |  |
| ring axis rho=-2.0000mm B=20mT laps=20 steps=20 | 140.7 | lost | 253.58 | 0.78 |  |  |  |
| ring axis rho=+3.0000mm B=20mT laps=20 steps=20 | 316.5 | lost | 374.73 | 0.53 |  |  |  |
| ring axis rho=-3.0000mm B=20mT laps=20 steps=20 | 316.5 | lost | 379.24 | 0.52 |  |  |  |
| trefoil axis rho=+0.5000mm B=20mT laps=20 steps=20 | 8.8 | held | 4.30 |  | 1.48 at +0° | 1.48 | -163.9° |
| trefoil axis rho=-0.5000mm B=20mT laps=20 steps=20 | 8.8 | held | 5.64 |  | 2.67 at +180° | 2.67 | +163.9° |
| trefoil matched rho=+0.5000mm B=20mT laps=20 steps=20 | 8.8 | held | 2.81 |  | 1.48 at +0° | 0.00 | -163.9° |
| trefoil matched rho=-0.5000mm B=20mT laps=20 steps=20 | 8.8 | held | 2.98 |  | 2.67 at +180° | 0.00 | +163.9° |
| trefoil axis rho=+1.0000mm B=20mT laps=20 steps=20 | 35.2 | held | 9.02 |  | 3.56 at +0° | 3.56 | -163.9° |
| trefoil axis rho=-1.0000mm B=20mT laps=20 steps=20 | 35.2 | lost | 10.05 | 1.10 | 4.74 at +180° | 4.74 | +163.9° |
| trefoil matched rho=+1.0000mm B=20mT laps=20 steps=20 | 35.2 | held | 5.47 |  | 3.56 at +0° | 0.00 | -163.9° |
| trefoil matched rho=-1.0000mm B=20mT laps=20 steps=20 | 35.2 | held | 5.31 |  | 4.74 at +180° | 0.00 | +163.9° |
| trefoil axis rho=+1.5000mm B=20mT laps=20 steps=20 | 79.1 | lost | 13.75 | 0.73 | 5.64 at +0° | 5.64 | -163.9° |
| trefoil axis rho=-1.5000mm B=20mT laps=20 steps=20 | 79.1 | lost | 14.63 | 0.54 | 6.82 at +180° | 6.82 | +163.9° |
| trefoil matched rho=+1.5000mm B=20mT laps=20 steps=20 | 79.1 | held | 8.12 |  | 5.64 at +0° | 0.00 | -163.9° |
| trefoil matched rho=-1.5000mm B=20mT laps=20 steps=20 | 79.1 | held | 7.81 |  | 6.82 at +180° | 0.00 | +163.9° |
| trefoil axis rho=+2.0000mm B=20mT laps=20 steps=20 | 140.7 | lost | 18.48 | 0.43 | 7.71 at +0° | 7.71 | -163.9° |
| trefoil axis rho=-2.0000mm B=20mT laps=20 steps=20 | 140.7 | lost | 19.35 | 0.40 | 8.89 at +180° | 8.89 | +163.9° |
| trefoil matched rho=+2.0000mm B=20mT laps=20 steps=20 | 140.7 | lost | 10.77 | 0.14 | 7.71 at +0° | 0.00 | -163.9° |
| trefoil matched rho=-2.0000mm B=20mT laps=20 steps=20 | 140.7 | lost | 10.46 | 0.12 | 8.89 at +180° | 0.00 | +163.9° |
| trefoil axis rho=+3.0000mm B=20mT laps=20 steps=20 | 316.5 | lost | 27.94 | 0.31 | 11.86 at +0° | 11.86 | -163.9° |
| trefoil axis rho=-3.0000mm B=20mT laps=20 steps=20 | 316.5 | lost | 28.80 | 0.28 | 13.04 at +180° | 13.04 | +163.9° |
| trefoil matched rho=+3.0000mm B=20mT laps=20 steps=20 | 316.5 | not launched |  |  | 11.86 at +0° |  | -163.9° |
| trefoil matched rho=-3.0000mm B=20mT laps=20 steps=20 | 316.5 | not launched |  |  | 13.04 at +180° |  | +163.9° |
| trefoil-shallow axis rho=+0.2500mm B=20mT laps=20 steps=20 | 2.2 | held | 5.36 |  | 2.41 at +0° | 2.41 | -52.3° |
| trefoil-shallow axis rho=-0.2500mm B=20mT laps=20 steps=20 | 2.2 | held | 9.09 |  | 4.50 at +180° | 4.50 | +52.3° |
| trefoil-shallow matched rho=+0.2500mm B=20mT laps=20 steps=20 | 2.2 | held | 2.95 |  | 2.41 at +0° | 0.00 | -52.3° |
| trefoil-shallow matched rho=-0.2500mm B=20mT laps=20 steps=20 | 2.2 | held | 4.59 |  | 4.50 at +180° | 0.00 | +52.3° |
| trefoil-shallow axis rho=+0.5000mm B=20mT laps=20 steps=20 | 8.8 | lost | 12.41 | 2.11 | 5.87 at +0° | 5.87 | -52.3° |
| trefoil-shallow axis rho=-0.5000mm B=20mT laps=20 steps=20 | 8.8 | lost | 16.04 | 1.49 | 7.96 at +180° | 7.96 | +52.3° |
| trefoil-shallow matched rho=+0.5000mm B=20mT laps=20 steps=20 | 8.8 | held | 6.54 |  | 5.87 at +0° | 0.00 | -52.3° |
| trefoil-shallow matched rho=-0.5000mm B=20mT laps=20 steps=20 | 8.8 | held | 8.08 |  | 7.96 at +180° | 0.00 | +52.3° |
| trefoil-shallow axis rho=+1.0000mm B=20mT laps=20 steps=20 | 35.2 | lost | 26.51 | 0.81 | 12.79 at +0° | 12.79 | -52.3° |
| trefoil-shallow axis rho=-1.0000mm B=20mT laps=20 steps=20 | 35.2 | lost | 29.98 | 0.73 | 14.88 at +180° | 14.88 | +52.3° |
| trefoil-shallow matched rho=+1.0000mm B=20mT laps=20 steps=20 | 35.2 | not launched |  |  | 12.79 at +0° |  | -52.3° |
| trefoil-shallow matched rho=-1.0000mm B=20mT laps=20 steps=20 | 35.2 | not launched |  |  | 14.88 at +180° |  | +52.3° |
| trefoil-shallow axis rho=+1.5000mm B=20mT laps=20 steps=20 | 79.1 | lost | 40.60 | 0.58 | 19.71 at +0° | 19.71 | -52.3° |
| trefoil-shallow axis rho=-1.5000mm B=20mT laps=20 steps=20 | 79.1 | lost | 43.99 | 0.50 | 21.80 at +180° | 21.80 | +52.3° |
| trefoil-shallow matched rho=+1.5000mm B=20mT laps=20 steps=20 | 79.1 | not launched |  |  | 19.71 at +0° |  | -52.3° |
| trefoil-shallow matched rho=-1.5000mm B=20mT laps=20 steps=20 | 79.1 | not launched |  |  | 21.80 at +180° |  | +52.3° |
| trefoil-shallow axis rho=+2.0000mm B=20mT laps=20 steps=20 | 140.7 | lost | 54.70 | 0.40 | 26.62 at +0° | 26.62 | -52.3° |
| trefoil-shallow axis rho=-2.0000mm B=20mT laps=20 steps=20 | 140.7 | lost | 58.08 | 0.38 | 28.71 at +180° | 28.71 | +52.3° |
| trefoil-shallow matched rho=+2.0000mm B=20mT laps=20 steps=20 | 140.7 | not launched |  |  | 26.62 at +0° |  | -52.3° |
| trefoil-shallow matched rho=-2.0000mm B=20mT laps=20 steps=20 | 140.7 | not launched |  |  | 28.71 at +180° |  | +52.3° |
| trefoil-shallow axis rho=+3.0000mm B=20mT laps=20 steps=20 | 316.5 | lost | 82.89 | 0.28 | 40.46 at +0° | 40.46 | -52.3° |
| trefoil-shallow axis rho=-3.0000mm B=20mT laps=20 steps=20 | 316.5 | lost | 86.27 | 0.27 | 42.55 at +180° | 42.55 | +52.3° |
| trefoil-shallow matched rho=+3.0000mm B=20mT laps=20 steps=20 | 316.5 | not launched |  |  | 40.46 at +0° |  | -52.3° |
| trefoil-shallow matched rho=-3.0000mm B=20mT laps=20 steps=20 | 316.5 | not launched |  |  | 42.55 at +180° |  | +52.3° |
| helix axis rho=+0.5000mm B=20mT laps=20 steps=20 | 8.8 | held | 1.66 |  | 0.06 at -179° | 0.06 | +123.2° |
| helix axis rho=-0.5000mm B=20mT laps=20 steps=20 | 8.8 | held | 2.12 |  | 0.97 at -180° | 0.97 | -123.2° |
| helix matched rho=+0.5000mm B=20mT laps=20 steps=20 | 8.8 | held | 1.60 |  | 0.06 at -179° | 0.00 | +123.2° |
| helix matched rho=-0.5000mm B=20mT laps=20 steps=20 | 8.8 | held | 1.15 |  | 0.97 at -180° | 0.00 | -123.2° |
| helix axis rho=+1.0000mm B=20mT laps=20 steps=20 | 35.2 | held | 3.25 |  | 0.40 at -0° | 0.40 | +123.2° |
| helix axis rho=-1.0000mm B=20mT laps=20 steps=20 | 35.2 | held | 3.58 |  | 1.42 at -180° | 1.42 | -123.2° |
| helix matched rho=+1.0000mm B=20mT laps=20 steps=20 | 35.2 | held | 2.85 |  | 0.40 at -0° | 0.00 | +123.2° |
| helix matched rho=-1.0000mm B=20mT laps=20 steps=20 | 35.2 | held | 2.15 |  | 1.42 at -180° | 0.00 | -123.2° |
| helix axis rho=+1.5000mm B=20mT laps=20 steps=20 | 79.1 | held | 4.95 |  | 0.86 at -0° | 0.86 | +123.2° |
| helix axis rho=-1.5000mm B=20mT laps=20 steps=20 | 79.1 | held | 5.27 |  | 1.88 at -180° | 1.88 | -123.2° |
| helix matched rho=+1.5000mm B=20mT laps=20 steps=20 | 79.1 | held | 4.10 |  | 0.86 at -0° | 0.00 | +123.2° |
| helix matched rho=-1.5000mm B=20mT laps=20 steps=20 | 79.1 | held | 3.39 |  | 1.88 at -180° | 0.00 | -123.2° |
| helix axis rho=+2.0000mm B=20mT laps=20 steps=20 | 140.7 | held | 6.66 |  | 1.31 at -0° | 1.31 | +123.2° |
| helix axis rho=-2.0000mm B=20mT laps=20 steps=20 | 140.7 | held | 6.98 |  | 2.34 at -180° | 2.34 | -123.2° |
| helix matched rho=+2.0000mm B=20mT laps=20 steps=20 | 140.7 | held | 5.35 |  | 1.31 at -0° | 0.00 | +123.2° |
| helix matched rho=-2.0000mm B=20mT laps=20 steps=20 | 140.7 | held | 4.64 |  | 2.34 at -180° | 0.00 | -123.2° |
| helix axis rho=+3.0000mm B=20mT laps=20 steps=20 | 316.5 | lost | 10.07 | 0.83 | 2.23 at -0° | 2.23 | +123.2° |
| helix axis rho=-3.0000mm B=20mT laps=20 steps=20 | 316.5 | lost | 10.39 | 0.81 | 3.25 at -180° | 3.25 | -123.2° |
| helix matched rho=+3.0000mm B=20mT laps=20 steps=20 | 316.5 | held | 7.84 |  | 2.23 at -0° | 0.00 | +123.2° |
| helix matched rho=-3.0000mm B=20mT laps=20 steps=20 | 316.5 | held | 7.14 |  | 3.25 at -180° | 0.00 | -123.2° |
| figure-eight axis rho=+0.5000mm B=20mT laps=20 steps=20 | 8.8 | held | 2.99 |  | 1.38 at +90° | 1.38 | -142.4° |
| figure-eight axis rho=-0.5000mm B=20mT laps=20 steps=20 | 8.8 | held | 2.79 |  | 1.40 at -90° | 1.40 | +142.4° |
| figure-eight matched rho=+0.5000mm B=20mT laps=20 steps=20 | 8.8 | held | 1.62 |  | 1.38 at +90° | 0.00 | -142.4° |
| figure-eight matched rho=-0.5000mm B=20mT laps=20 steps=20 | 8.8 | held | 1.40 |  | 1.40 at -90° | 0.00 | +142.4° |
| figure-eight axis rho=+1.0000mm B=20mT laps=20 steps=20 | 35.2 | held | 5.52 |  | 2.76 at +90° | 2.76 | -142.4° |
| figure-eight axis rho=-1.0000mm B=20mT laps=20 steps=20 | 35.2 | held | 5.57 |  | 2.78 at -90° | 2.78 | +142.4° |
| figure-eight matched rho=+1.0000mm B=20mT laps=20 steps=20 | 35.2 | held | 2.76 |  | 2.76 at +90° | 0.00 | -142.4° |
| figure-eight matched rho=-1.0000mm B=20mT laps=20 steps=20 | 35.2 | held | 2.78 |  | 2.78 at -90° | 0.00 | +142.4° |
| figure-eight axis rho=+1.5000mm B=20mT laps=20 steps=20 | 79.1 | held | 8.30 |  | 4.15 at +90° | 4.15 | -142.4° |
| figure-eight axis rho=-1.5000mm B=20mT laps=20 steps=20 | 79.1 | held | 8.34 |  | 4.17 at -90° | 4.17 | +142.4° |
| figure-eight matched rho=+1.5000mm B=20mT laps=20 steps=20 | 79.1 | held | 4.15 |  | 4.15 at +90° | 0.00 | -142.4° |
| figure-eight matched rho=-1.5000mm B=20mT laps=20 steps=20 | 79.1 | held | 4.17 |  | 4.17 at -90° | 0.00 | +142.4° |
| figure-eight axis rho=+2.0000mm B=20mT laps=20 steps=20 | 140.7 | lost | 11.07 | 0.86 | 5.54 at +90° | 5.54 | -142.4° |
| figure-eight axis rho=-2.0000mm B=20mT laps=20 steps=20 | 140.7 | lost | 11.11 | 0.89 | 5.56 at -90° | 5.56 | +142.4° |
| figure-eight matched rho=+2.0000mm B=20mT laps=20 steps=20 | 140.7 | held | 5.54 |  | 5.54 at +90° | 0.00 | -142.4° |
| figure-eight matched rho=-2.0000mm B=20mT laps=20 steps=20 | 140.7 | held | 5.56 |  | 5.56 at -90° | 0.00 | +142.4° |
| figure-eight axis rho=+3.0000mm B=20mT laps=20 steps=20 | 316.5 | lost | 16.61 | 0.35 | 8.31 at +90° | 8.31 | -142.4° |
| figure-eight axis rho=-3.0000mm B=20mT laps=20 steps=20 | 316.5 | lost | 16.65 | 0.34 | 8.33 at -90° | 8.33 | +142.4° |
| figure-eight matched rho=+3.0000mm B=20mT laps=20 steps=20 | 316.5 | held | 8.31 |  | 8.31 at +90° | 0.00 | -142.4° |
| figure-eight matched rho=-3.0000mm B=20mT laps=20 steps=20 | 316.5 | held | 8.33 |  | 8.33 at -90° | 0.00 | +142.4° |

## The same gyroradii at 40 mT

| Shot | Energy (eV) | Predicted | Reach (mm) | Lap of loss | Closed orbit (mm) | Launch to closed orbit (mm) | Turn per lap |
| --- | --- | --- | --- | --- | --- | --- | --- |
| trefoil axis rho=+0.5000mm B=40mT laps=20 steps=20 | 35.2 | held | 4.30 |  | 1.48 at +0° | 1.48 | -163.9° |
| trefoil axis rho=-0.5000mm B=40mT laps=20 steps=20 | 35.2 | held | 5.64 |  | 2.67 at +180° | 2.67 | +163.9° |
| trefoil matched rho=+0.5000mm B=40mT laps=20 steps=20 | 35.2 | held | 2.81 |  | 1.48 at +0° | 0.00 | -163.9° |
| trefoil matched rho=-0.5000mm B=40mT laps=20 steps=20 | 35.2 | held | 2.98 |  | 2.67 at +180° | 0.00 | +163.9° |
| trefoil axis rho=+1.0000mm B=40mT laps=20 steps=20 | 140.7 | held | 9.02 |  | 3.56 at +0° | 3.56 | -163.9° |
| trefoil axis rho=-1.0000mm B=40mT laps=20 steps=20 | 140.7 | lost | 10.05 | 1.10 | 4.74 at +180° | 4.74 | +163.9° |
| trefoil matched rho=+1.0000mm B=40mT laps=20 steps=20 | 140.7 | held | 5.47 |  | 3.56 at +0° | 0.00 | -163.9° |
| trefoil matched rho=-1.0000mm B=40mT laps=20 steps=20 | 140.7 | held | 5.31 |  | 4.74 at +180° | 0.00 | +163.9° |
| helix axis rho=+0.5000mm B=40mT laps=20 steps=20 | 35.2 | held | 1.66 |  | 0.06 at -179° | 0.06 | +123.2° |
| helix axis rho=-0.5000mm B=40mT laps=20 steps=20 | 35.2 | held | 2.12 |  | 0.97 at -180° | 0.97 | -123.2° |
| helix matched rho=+0.5000mm B=40mT laps=20 steps=20 | 35.2 | held | 1.60 |  | 0.06 at -179° | 0.00 | +123.2° |
| helix matched rho=-0.5000mm B=40mT laps=20 steps=20 | 35.2 | held | 1.15 |  | 0.97 at -180° | 0.00 | -123.2° |
| helix axis rho=+1.0000mm B=40mT laps=20 steps=20 | 140.7 | held | 3.25 |  | 0.40 at -0° | 0.40 | +123.2° |
| helix axis rho=-1.0000mm B=40mT laps=20 steps=20 | 140.7 | held | 3.58 |  | 1.42 at -180° | 1.42 | -123.2° |
| helix matched rho=+1.0000mm B=40mT laps=20 steps=20 | 140.7 | held | 2.85 |  | 0.40 at -0° | 0.00 | +123.2° |
| helix matched rho=-1.0000mm B=40mT laps=20 steps=20 | 140.7 | held | 2.15 |  | 1.42 at -180° | 0.00 | -123.2° |

## The instrument made finer, one part at a time

| Shot | Energy (eV) | Predicted | Reach (mm) | Lap of loss | Closed orbit (mm) | Launch to closed orbit (mm) | Turn per lap |
| --- | --- | --- | --- | --- | --- | --- | --- |
| trefoil axis rho=+1.0000mm B=20mT laps=20 steps=40 | 35.2 | held | 9.02 |  | 3.56 at +0° | 3.56 | -163.9° |
| trefoil axis rho=+1.0000mm B=20mT laps=60 steps=20 | 35.2 | held | 9.03 |  | 3.56 at +0° | 3.56 | -163.9° |
| trefoil-fine-curve axis rho=+1.0000mm B=20mT laps=20 steps=20 | 35.2 | held | 9.01 |  | 3.55 at +0° | 3.55 | -164.0° |
| trefoil-fine-winding axis rho=+1.0000mm B=20mT laps=20 steps=20 | 35.2 | held | 9.02 |  | 3.56 at +0° | 3.56 | -163.9° |
| trefoil matched rho=+1.0000mm B=20mT laps=20 steps=40 | 35.2 | held | 5.47 |  | 3.56 at +0° | 0.00 | -163.9° |
| trefoil matched rho=+1.0000mm B=20mT laps=60 steps=20 | 35.2 | held | 5.47 |  | 3.56 at +0° | 0.00 | -163.9° |
| trefoil-fine-curve matched rho=+1.0000mm B=20mT laps=20 steps=20 | 35.2 | held | 5.46 |  | 3.55 at +0° | 0.00 | -164.0° |
| trefoil-fine-winding matched rho=+1.0000mm B=20mT laps=20 steps=20 | 35.2 | held | 5.47 |  | 3.56 at +0° | 0.00 | -163.9° |

## The first run's cases, on its single-wire winding

| Shot | Energy (eV) | Predicted | Reach (mm) | Lap of loss | Closed orbit (mm) | Launch to closed orbit (mm) | Turn per lap |
| --- | --- | --- | --- | --- | --- | --- | --- |
| trefoil-single first-run rho=+6.7446mm B=5mT laps=20 steps=20 | 100.0 | lost | 44.68 | 0.11 | 19.17 at +0° | 17.17 | -163.9° |
| trefoil-single first-run rho=+11.6831mm B=5mT laps=20 steps=20 | 300.0 | lost | 91.38 | 0.06 | 39.67 at +0° | 37.67 | -163.9° |
| trefoil-single first-run rho=+21.3377mm B=5mT laps=20 steps=20 | 1000.0 | lost | 182.68 | 0.03 | 79.75 at +0° | 77.75 | -163.9° |
| trefoil-single first-run rho=+36.9940mm B=5mT laps=20 steps=20 | 3000.0 | lost | 330.73 | 0.02 | 144.73 at +0° | 142.73 | -163.9° |
| trefoil-single first-run rho=+3.3723mm B=10mT laps=20 steps=20 | 100.0 | lost | 12.78 | 0.78 | 5.18 at +0° | 3.18 | -163.9° |
| trefoil-single first-run rho=+5.8416mm B=10mT laps=20 steps=20 | 300.0 | lost | 36.14 | 0.12 | 15.43 at +0° | 13.43 | -163.9° |
| trefoil-single first-run rho=+10.6688mm B=10mT laps=20 steps=20 | 1000.0 | lost | 81.79 | 0.07 | 35.46 at +0° | 33.46 | -163.9° |
| trefoil-single first-run rho=+18.4970mm B=10mT laps=20 steps=20 | 3000.0 | lost | 155.81 | 0.04 | 67.95 at +0° | 65.95 | -163.9° |
| trefoil-single first-run rho=+1.6861mm B=20mT laps=20 steps=20 | 100.0 | held | 5.64 |  | 1.82 at +179° | 3.82 | -163.9° |
| trefoil-single first-run rho=+2.9208mm B=20mT laps=20 steps=20 | 300.0 | held | 8.52 |  | 3.30 at +0° | 1.30 | -163.9° |
| trefoil-single first-run rho=+5.3344mm B=20mT laps=20 steps=20 | 1000.0 | lost | 31.34 | 0.13 | 13.32 at +0° | 11.32 | -163.9° |
| trefoil-single first-run rho=+9.2485mm B=20mT laps=20 steps=20 | 3000.0 | lost | 68.35 | 0.08 | 29.57 at +0° | 27.57 | -163.9° |
| ring-single first-run rho=+6.7446mm B=5mT laps=20 steps=20 | 100.0 | lost | 631.45 | 0.31 |  |  |  |
| ring-single first-run rho=+11.6831mm B=5mT laps=20 steps=20 | 300.0 | lost | 1252.03 | 0.15 |  |  |  |
| ring-single first-run rho=+21.3377mm B=5mT laps=20 steps=20 | 1000.0 | lost | 2465.23 | 0.08 |  |  |  |
| ring-single first-run rho=+36.9940mm B=5mT laps=20 steps=20 | 3000.0 | lost | 4432.63 | 0.04 |  |  |  |
| ring-single first-run rho=+3.3723mm B=10mT laps=20 steps=20 | 100.0 | lost | 207.69 | 0.94 |  |  |  |
| ring-single first-run rho=+5.8416mm B=10mT laps=20 steps=20 | 300.0 | lost | 517.97 | 0.38 |  |  |  |
| ring-single first-run rho=+10.6688mm B=10mT laps=20 steps=20 | 1000.0 | lost | 1124.57 | 0.17 |  |  |  |
| ring-single first-run rho=+18.4970mm B=10mT laps=20 steps=20 | 3000.0 | lost | 2108.27 | 0.09 |  |  |  |
| ring-single first-run rho=+1.6861mm B=20mT laps=20 steps=20 | 100.0 | held | 4.66 |  |  |  |  |
| ring-single first-run rho=+2.9208mm B=20mT laps=20 steps=20 | 300.0 | lost | 150.95 | 1.28 |  |  |  |
| ring-single first-run rho=+5.3344mm B=20mT laps=20 steps=20 | 1000.0 | lost | 454.25 | 0.43 |  |  |  |
| ring-single first-run rho=+9.2485mm B=20mT laps=20 steps=20 | 3000.0 | lost | 946.09 | 0.20 |  |  |  |
