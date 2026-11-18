<!-- Written by `experiments/guide_posthoc.py`, after the run. Exploratory: not a test. -->

# Experiment 1, read after the run: the turn per lap

Fitted minus predicted turn per lap, for the on-axis shots held in both directions. The shared part is
what the two directions have in common, the reversing part what changes sign with the direction. The
candidate is ρ × (the curve's bending) / 2, from the curve alone, with nothing fitted.

| Guide | ρ (mm) | Along | Against | Shared | Candidate | Shared / candidate | Reversing |
| --- | --- | --- | --- | --- | --- | --- | --- |
| figure-eight | 0.50 | +2.05° | +2.82° | +2.43° | +2.44° | 1.00 | -0.38° |
| figure-eight | 1.00 | +4.47° | +5.27° | +4.87° | +4.88° | 1.00 | -0.40° |
| figure-eight | 1.50 | +6.93° | +7.65° | +7.29° | +7.31° | 1.00 | -0.36° |
| helix | 1.00 | +11.52° | +0.75° | +6.14° | +6.13° | 1.00 | +5.38° |
| helix | 1.50 | +14.09° | +4.38° | +9.23° | +9.19° | 1.00 | +4.86° |
| helix | 2.00 | +16.55° | +8.07° | +12.31° | +12.25° | 1.00 | +4.24° |
| trefoil | 0.50 | +11.33° | -1.97° | +4.68° | +4.65° | 1.01 | +6.65° |
| trefoil-shallow | 0.25 | +3.07° | +0.83° | +1.95° | +1.72° | 1.13 | +1.12° |

| Guide | Length (mm) | Total curvature / 2π | Bending (1/m) | Bending × length / 4π² |
| --- | --- | --- | --- | --- |
| ring | 932 | 1.00 | 42.4 | 1.00 |
| trefoil | 932 | 2.73 | 324.8 | 7.67 |
| trefoil-shallow | 814 | 2.16 | 239.5 | 4.94 |
| helix | 897 | 2.20 | 213.8 | 4.86 |
| figure-eight | 875 | 1.49 | 170.2 | 3.77 |
