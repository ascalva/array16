# Experiment 1: is it the knot, or the twist? The protocol

Written on 2026-11-18, before the tracking, with one exception that section 6 states. This file,
`predictions.json` and `predictions.md` are committed together, and the commit that adds the tracking's output
comes after it: `array16.lab` refuses to measure until they are, and stores with every measurement the commit
it was registered at. The numbers that decide the verdicts are constants at the top of
`experiments/guide_study.py`, in the same commit. Anything done differently from this protocol is listed in the
report.

## 1. The question

The owner's idea is a field whose knotted shape holds a particle's path. The first run
(`experiments/accelerator_study.py`, 2026-11-16) tracked electrons in a solenoid wound on a trefoil knot and in one
wound on a plain ring of the same length, and reported that the knot holds them and the ring loses them. That
comparison was confounded. A single-wire winding carries its current once around its curve, which on the ring
makes a field through the hole that bends one momentum, so the ring held a 100 eV electron at 20 mT by
coincidence.

What property of a guide's axis makes it hold an electron: that it is knotted, or something else?

## 2. The theory, and what it claims

An electron streaming along a curved field line drifts across it. Follow its guiding centre's offset from the
axis as a complex number `z` on axes carried round the curve without twisting. To first order in gyroradius and
offset over bend radius, one lap is the map

    z ↦ exp(iα)·(z + ρG + H)

- `α` is the angle by which the curve turns a frame carried once round it: 2π times the curve's writhe, less
  whole turns (ledger CUR-1, after Călugăreanu, White and Fuller [2]). It is zero for every planar curve.
- `G` is the curvature drift a lap accumulates, per unit gyroradius `ρ = p/(eB)`; `ρ` is negative for an electron
  sent against the field.
- `H` is the wander of the winding's own field line through the axis, from the coil's computed field.

The laws of the map are proved (`docs/LAWS.md`, GDE-1 to GDE-7; `formal/Array16/Guide.lean`). That an electron
obeys the map is the hypothesis under test, GDE-H1. If it does:

1. **A planar ring never holds.** With `α = 0` every lap adds the same step, 2πρ per lap (GDE-6).
2. **Any curve that turns the plane holds what is small enough.** There is one closed orbit, at an offset
   proportional to `ρ`, and every other orbit circles it (GDE-1 to GDE-5). A guide of bore `a` holds the
   electrons whose reach, a multiple of `ρ`, is less than `a`.
3. **The knot does not enter.** Only `α`, `G` and `H` do. An unknotted curve that leaves the plane turns the
   plane too, and the theory expects two such curves here to hold more than the trefoil does: about three times
   the gyroradius on the helix.
4. **Direction matters, energy and field only through ρ.** Sent the other way, an electron circles a closed orbit
   on the other side of the axis, the other way round. Doubling the field and quadrupling the energy changes
   nothing.

Claim 3 is not new physics. Spitzer's figure-eight stellarator of 1958 rests on it: "a rotational transform may be
generated ... by a solenoidal field in a twisted, or figure-eight shaped, tube" [1]. Guides in the shape of a
knot have been proposed for plasma since [3]. What is tested here is whether this code's first-order form of it
predicts this code's tracking, in numbers, for the windings a table-top build would use.

## 3. The apparatus

All numerical. `predictions.md` gives each guide's numbers.

| Guide | Curve | Knotted | Planar |
| --- | --- | --- | --- |
| ring | circle, radius 148.4 mm | no | yes |
| trefoil | (2,3) torus knot on a torus of radii 60 and 28 mm (the first run's) | yes | no |
| trefoil-shallow | the same on radii 60 and 16 mm: a smaller turn per lap | yes | no |
| helix | (1,3) torus curve on radii 110 and 30 mm: once round the hole, three times round the tube | no | no |
| figure-eight | 280 by 140 mm, the strands 80 mm apart where they cross | no | no |

- **Winding.** 180 turns of radius 11 mm, then 180 more of radius 12 mm wound back to the start with the same
  circulation: 360 turns, and no net current along the curve (ledger COI-1, COI-2). Each turn is 36 straight
  filaments. The first run's single-wire winding is kept for one group of shots.
- **Field.** 20 mT on the axis, and 40 mT for one group.
- **Wall.** An electron further than 9.9 mm from the axis is lost, as in the first run.
- **Tracker.** `array16.particles.track_guide`: the relativistic Boris step in the winding's Biot-Savart field,
  20 steps per gyration, 20 laps. It is calibrated against closed forms (ledger INT-1 to INT-3), its field sum
  against the textbook sum (COI-3), and the whole against the drift on a ring over two laps
  (`test_tracked_electron_in_a_ring_coil_drifts_as_the_theory_says`).
- **What is measured.** Whether the electron is still inside after its laps; its largest distance from the axis;
  the lap at which it was lost; and its offset each time it completes a lap, to which a circle is fitted for the
  closed orbit and the turn per lap.

## 4. The shots

142 in all (`guide_study.py::shots`), each one electron. `predictions.md` lists every one with its prediction.

- **Main, 94.** On each of the five guides, gyroradii of 0.5, 1, 1.5, 2 and 3 mm (8.8 to 317 eV at 20 mT),
  and 0.25 mm on the shallow trefoil; launched on the axis and on the predicted closed orbit; sent both ways. The
  ring has no closed orbit to launch on.
- **Scaling, 16.** The trefoil and the helix at 40 mT, gyroradii 0.5 and 1 mm: four times the energy.
- **Systematics, 8.** The trefoil's 1 mm shots repeated with 40 steps per gyration, with 60 laps, with the curve
  drawn from 1200 points, and with 72 filaments per turn.
- **The first run's cases, 24.** Its trefoil and ring, single-wire winding, 5, 10 and 20 mT, 100 to 3000 eV,
  launched 2 mm off the axis as it launched them.

## 5. What counts as agreement

Fixed before the run. "Well inside" means a predicted reach of at most 70 % of the wall, 6.9 mm.

| | Criterion | Passes when |
| --- | --- | --- |
| C1 | Reach | every shot predicted well inside stays inside, and its tracked reach is within 15 % + 0.1 mm of the predicted |
| C2 | Held or lost | every shot whose predicted reach is more than 15 % from the wall, either side, is held or lost as predicted |
| C3 | Closed orbit | for held shots not launched on the closed orbit, the fitted centre is within 15 % of its offset + 0.1 mm of the predicted closed orbit |
| C4 | Turn per lap | for the same shots, the fitted turn is within 8° of the predicted |
| C5 | Matched launch | a shot launched on the predicted closed orbit, well inside, never crosses the start further than 25 % of that orbit's offset + 0.1 mm from it |
| C6 | Ring | every electron on the ring is lost, at the predicted lap within 15 % or 0.3 lap |
| C7 | A twist without a knot | C1 to C5 hold on the helix and the figure-eight |
| C8 | Scaling | each 40 mT shot has the reach of its 20 mT twin within 3 % + 0.02 mm, and its turn within 2° |
| C9 | Against the field | C1 to C5 hold on the shots sent backwards |
| C10 | Systematics | each finer repeat has the reach of the standard shot within 3 % + 0.02 mm, and its turn within 2° |
| C11 | The first run | reported, not a test: see section 6 |

A circle is fitted only to a held orbit with at least 8 crossings spread over more than 0.3 mm.

**Why these sizes.** The theory drops terms of second order in offset over bend radius. Well inside, offsets
reach 7 mm against bend radii of 32 to 60 mm, so those terms are 1 to 5 % with coefficients of order one, and the
winding's ripple grows near the wall. Fifteen percent leaves room for that. It leaves none for a wrong sign or a
wrong turn: launching on the mirror image of the closed orbit, or turning the plane the other way, changes the
predicted reach by a factor near two.

**What a failure would mean.**
- C3 or C4 off on every guide: the map has the wrong form, and GDE-H1 is false as stated.
- C1 to C5 good at small gyroradius and failing as it grows: the first order is right and its range is shorter
  than assumed. The report gives the trend.
- C7 failing while the trefoils pass: the knot matters beyond the turn it gives, and the theory is incomplete.
- C6 failing, an electron held on the ring: the winding or the tracker is at fault, and nothing else here stands.
- C8 or C10 failing: the result depends on the instrument, and nothing else here stands until it does not.

## 6. What was known when this was written

- **The first run's table**, the 24 cases of its single-wire trefoil and ring. Before this protocol was written
  the theory was set against it: it gives held or lost correctly in all 23 cases it judges (one lies within 15 %
  of the wall), with the laps of loss on the ring at 0.94, 1.28, 0.38 and 0.43 against the first run's 0.9, 1.3,
  0.4 and 0.4, and the trefoil's reaches at 5.6 and 8.5 mm against 5.5 and 8.8 mm. This is a retrodiction. The
  theory was not fitted to it, and it has no free parameter to fit, but it was known.
- **One ring**, two laps, gyroradius 0.4 mm: the kernel's own calibration test.
- **One shot of this experiment, by mistake.** A timing test of the faster field sum tracked the trefoil's
  on-axis shot at ρ = 1 mm for two laps and printed its reach: 9.93 mm, so it had reached the wall, where the
  theory puts its reach over 20 laps at 9.02 mm. That shot's predicted reach is within 15 % of the wall, so C1
  does not apply to it and C2 does not judge it; it bears on C3 and C4 only if it is held, which it will not be.
  The protocol's text, the shots and every tolerance were written before that run and were not changed after it.
  The git history cannot show this, since nothing had been committed; it is stated here instead.
- Nothing else. No electron had been tracked in the shallow trefoil, in either unknotted curve, on a closed
  orbit, or against the field.

## 7. What is left out

One electron at a time: no space charge, no collisions, no radiation. Filament windings with no thickness, leads
or supports. No vacuum pipe, so the wall is a distance and not a surface. The launch is exact. A real coil's
winding errors would add to `H`.

## 8. Running it

    uv run python experiments/guide_study.py predict     # the predictions; a few seconds
    uv run python experiments/guide_study.py run         # the tracking; resumable, one guide with --only
    uv run python experiments/guide_study.py report      # the verdicts

## Sources

1. L. Spitzer Jr., "The Stellarator Concept," The Physics of Fluids 1(4), 253-264 (1958),
   doi:10.1063/1.1705883.
2. F. B. Fuller, "The Writhing Number of a Space Curve," Proceedings of the National Academy of Sciences 68(4),
   815-819 (1971), doi:10.1073/pnas.68.4.815.
3. S. R. Hudson, E. Startsev, E. Feibush, "A new class of magnetic confinement device in the shape of a knot,"
   Physics of Plasmas 21, 010705 (2014), doi:10.1063/1.4863844.
