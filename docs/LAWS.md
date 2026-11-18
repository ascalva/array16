# The laws ledger

One row per law the particle code relies on: what it says, which known structure it instantiates, where the code
holds it, how strong the claim is, which test would catch it false, and which Lean declarations state it. The form
is the one the owner's mind-palace uses for its own ledger. The ledger indexes and never restates: proofs live in
`formal/`, enforcement lives in the tests it names.

How to read a row:

- **statement**: the law in one clause, with its arithmetic domain in parentheses.
- **structure**: the known structure the law instantiates.
- **anchor**: `file::symbol`, the code that relies on the law.
- **label**: how strong the claim is, on its evidence.
  - *Proposition*: proved from the stated premises. The proof is machine-checked over exact arithmetic in
    `formal/`, and a test attacks the float implementation.
  - *Claim*: an assertion about the code that the named test supports. It is not proved here.
  - *Hypothesis*: a physical model. An experiment supports it or falsifies it; nothing proves it.
- **falsifier**: the pytest node id, or the experiment, that would catch the law false.
- **formal**: the Lean declarations that state the law, as `file::name` for a file under `formal/Array16/`, or `—`.
  `formal/check.sh` holds every one to a finished proof and the three standard axioms. A Lean proof is about a
  model of the code (`formal/Array16/Guide/Defs.lean`) over exact complex numbers, where the code computes in
  binary64. Read the model against the code to judge whether a row's formal column says what its statement says.

The cross-reference is checked in both directions by `tests/test_laws.py`. Every anchor, falsifier and declaration
a row cites exists. Each cited declaration carries its row's law id in bold at the head of its docstring, and every
declaration so marked is cited by that law's row. A Proposition cites at least one declaration; a Claim or a
Hypothesis cites none. Every test named for a law is one of that law's falsifiers.

## I. The lap map of a guide field

`array16/particles/optics.py`, `GuideOptics`. A magnetic guide field follows a closed curve. A guiding centre's offset from
the curve is a complex number `z` on axes carried along the curve without twisting. One lap adds a step
`ρ·G + H` and then turns the plane by `rot = exp(i·α)`: `lap(z) = rot·(z + step)`. `ρ` is the gyroradius, signed
by the direction of travel. GDE-1 to GDE-7 are laws of that map. GDE-H1 is the claim that an electron obeys it.

| law | statement | structure | anchor | label | falsifier | formal |
|---|---|---|---|---|---|---|
| GDE-1 | when the plane turns, `rot ≠ 1`, the offset `step·rot/(1 − rot)` is a fixed point of the lap map (ℂ: exact; binary64: to 2e-8·(1 + \|z\| + \|step\|) for turns at least 0.05 rad from a whole number) | an affine map of the plane whose linear part is a turn is itself a turn, about one point | `array16/particles/optics.py::GuideOptics.closed_orbit` | Proposition | `tests/test_particles.py::test_gde_1_the_closed_orbit_is_a_fixed_point`; control: `tests/test_particles.py::test_gde_1_control_a_nearby_formula_is_not_a_fixed_point` | `Guide.lean::lap_closedOrbit`; control: `Guide.lean::lap_one_closedOrbit_ne` |
| GDE-2 | a lap moves an offset by `rot − 1` times its offset from the closed orbit, so no other offset closes after one lap (ℂ: exact) | the same turn: only its centre stays | `array16/particles/optics.py::GuideOptics.closed_orbit` | Proposition | `tests/test_particles.py::test_gde_2_only_the_closed_orbit_closes` | `Guide.lean::eq_closedOrbit_of_lap_eq`, `Guide.lean::lap_sub_self`; control: `Guide.lean::lap_one_zero` |
| GDE-3 | after `n` laps the offset from the closed orbit has turned by `rot^n` (ℂ: exact) | iterating a turn | `array16/particles/optics.py::GuideOptics.lap` | Proposition | `tests/test_particles.py::test_gde_3_4_5_laps_turn_the_offset_about_the_closed_orbit` | `Guide.lean::iterate_lap_sub_closedOrbit` |
| GDE-4 | with `\|rot\| = 1` the distance from the closed orbit is the same at the start of every lap: the offsets lie on a circle about it (ℂ: exact) | a turn is an isometry | `array16/particles/optics.py::GuideOptics.lap` | Proposition | `tests/test_particles.py::test_gde_3_4_5_laps_turn_the_offset_about_the_closed_orbit` | `Guide.lean::norm_iterate_lap_sub_closedOrbit`; control: `Guide.lean::norm_lap_two_sub_closedOrbit_ne` |
| GDE-5 | so the reach is bounded whatever the number of laps, by the closed orbit's offset plus the launch's distance from it plus the largest wander within a lap (ℂ: exact) | the triangle inequality on GDE-4 | `array16/particles/optics.py::GuideOptics.excursion` | Proposition | `tests/test_particles.py::test_gde_3_4_5_laps_turn_the_offset_about_the_closed_orbit` | `Guide.lean::norm_iterate_lap_add_le` |
| GDE-6 | without a turn, `rot = 1` (a planar curve), every lap adds the same step: unless the step is zero no orbit closes and every orbit leaves every bounded region (ℂ: exact) | a translation | `array16/particles/optics.py::GuideOptics.lap`, `array16/particles/optics.py::GuideOptics.closes` | Proposition | `tests/test_particles.py::test_gde_6_without_a_turn_every_lap_adds_the_same_step`, `tests/test_particles.py::test_ring_guide_drifts_two_pi_rho_per_lap_and_never_closes` | `Guide.lean::iterate_lap_one`, `Guide.lean::lap_one_ne`, `Guide.lean::exists_lt_norm_iterate_lap_one` |
| GDE-7 | the reach from a fixed launch offset, and the reach of the closed orbit, are convex in `ρ`, so the `ρ` that stay inside a wall form one unbroken range (ℝ → ℝ: exact) | the largest of finitely many lengths of affine functions | `array16/particles/optics.py::GuideOptics.excursion`, `array16/particles/optics.py::GuideOptics.matched_excursion`, `array16/particles/optics.py::window` | Proposition | `tests/test_particles.py::test_gde_7_reach_is_convex_in_rho`; control: `tests/test_particles.py::test_gde_7_control_reach_is_not_concave` | `Guide.lean::reach_convexOn`, `Guide.lean::matchedReach_convexOn`, `Guide.lean::reach_lt_ordConnected`, `Guide.lean::matchedReach_lt_ordConnected` |
| GDE-H1 | an electron's guiding centre in the field of a winding obeys the lap map, with `α` the curve's holonomy, `G` the accumulated curvature drift `−∮ T × dT` and `H` the accumulated slope of the field line through the curve, to first order in gyroradius and offset over bend radius | guiding-centre drift theory; a non-planar magnetic axis turns the transverse plane, as in the figure-eight stellarator | `array16/particles/optics.py::GuideOptics.of`, `array16/particles/coils.py::Coil.optics` | Hypothesis | `tests/test_particles.py::test_tracked_electron_in_a_ring_coil_drifts_as_the_theory_says` (a ring, one lap: drift within 3 %) | — |

## II. Closed curves

`array16/particles/curves.py`, `ClosedCurve`.

| law | statement | structure | anchor | label | falsifier | formal |
|---|---|---|---|---|---|---|
| CUR-1 | a normal carried once round a closed curve without twisting comes back turned by `2π` times the curve's writhe, modulo `2π` (ℝ: a cited theorem; binary64 on an n-point polygon: within 1.5e-3 of a turn at 600 points on the four torus knots tested, and at least 2.5 times closer at 1200) | Călugăreanu, White and Fuller: linking = twist + writhe, for the ribbon a twist-free normal sweeps | `array16/particles/curves.py::ClosedCurve.holonomy`, `array16/particles/curves.py::ClosedCurve.writhe` | Claim | `tests/test_particles.py::test_holonomy_is_the_writhe_modulo_one_turn`, `tests/test_particles.py::test_planar_curve_has_no_holonomy_and_no_writhe` | — |

## III. The integrator

`array16/particles/push.py`, `push`: the relativistic Boris step.

| law | statement | structure | anchor | label | falsifier | formal |
|---|---|---|---|---|---|---|
| INT-1 | in a uniform magnetic field one step turns `u` by exactly the gyration angle of the step, at any step size, and keeps its length (ℝ: exact; binary64: phase within 1e-9 rad after 50 turns of 40 steps, length to 1e-12) | the Boris rotation with the tangent of the half angle: a rotation written as a Cayley transform | `array16/particles/push.py::push` | Claim | `tests/test_particles.py::test_boris_gyration_radius_phase_and_energy` | — |
| INT-2 | in crossed uniform fields the step's `E × B` drift is the exact drift times `(θ/2)/tan(θ/2)`, `θ` the gyration angle per step (binary64: to 1e-6 at 40 steps per turn, where the factor is 0.998) | the price of INT-1: the electric half-kicks straddle an exact rotation | `array16/particles/push.py::push` | Claim | `tests/test_particles.py::test_boris_exb_drift_and_its_known_step_error` | — |
| INT-3 | in an ideal quadrupole the step reproduces Mathieu's stability boundary (bounded at `q = 0.85`, growing at `q = 1.0`; the boundary is 0.908) and the slow frequency `q/√2` at `q = 0.2` to 2 % (binary64, 100 steps per period) | Mathieu's equation | `array16/particles/push.py::push` | Claim | `tests/test_particles.py::test_mathieu_stability_boundary_and_secular_frequency` | — |

## IV. Windings

`array16/particles/coils.py`, `Coil`.

| law | statement | structure | anchor | label | falsifier | formal |
|---|---|---|---|---|---|---|
| COI-1 | a single-layer winding carries its current once along its curve, one net turn through the hole, and a return winding carries none (filaments of 36 segments per turn: 1.00 ± 0.02 and 0.00 ± 0.02 by Ampère's law on a loop round the tube) | Ampère's law | `array16/particles/coils.py::Coil.wound`, `array16/particles/coils.py::Coil.linked_current` | Claim | `tests/test_particles.py::test_single_winding_carries_one_net_turn_and_return_winding_none` | — |
| COI-2 | a return winding on a ring leaves a field through the hole of `ln(1 + gap/a)/(2N)` of the guide field, from its layers' radii `a` and `a + gap` and its `N` turns (filaments: within 5 %) | two opposite ring currents on tubes of different radius | `array16/particles/coils.py::Coil.wound`, `array16/particles/coils.py::Coil.slope` | Claim | `tests/test_particles.py::test_return_winding_leaves_a_small_dipole_from_its_two_radii` | — |
| COI-3 | the fast field sum, `(K·dl) × x − K·(dl × mid)` with `K = \|x − mid\|⁻³`, is the plain Biot-Savart sum over the same filaments (ℝ: an identity, `dl × (x − mid) = dl × x − dl × mid`; binary64: to 1e-11 of the largest field, at points up to 1 mm from the wire, whether taken in one block or several) | the cross product is linear in its second factor | `array16/particles/coils.py::Coil.field`, `array16/particles/coils.py::Coil.field_plain` | Claim | `tests/test_particles.py::test_fast_field_sum_is_the_plain_sum` | — |
