"""Calibration of the particle kernel against results known in closed form, and attacks on the laws of the guide
optics. Each test named test_gde_N attacks the float implementation of the law proved as GDE-N in formal/."""

from __future__ import annotations

import cmath
import math

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from array16 import particles as pt
from array16.particles import ClosedCurve, Coil, GuideOptics, Launch


# ---- kinematics and fields-as-functions -----------------------------------------------------------------------
def test_rigidity_of_a_100_ev_electron_and_its_inverse() -> None:
    assert pt.rigidity(100.0) == pytest.approx(3.3723e-5, rel=1e-4)  # T m
    for w in (1.0, 100.0, 1e4, 1e6):
        assert pt.energy_of(pt.rigidity(w)) == pytest.approx(w, rel=1e-9)


def test_fields_scale_and_add() -> None:
    x = np.zeros((4, 3))
    bz, bx = pt.uniform(np.array([0.0, 0.0, 2.0])), pt.uniform(np.array([1.0, 0.0, 0.0]))
    assert np.allclose(pt.scaled(bz, 0.5)(x), [0.0, 0.0, 1.0])
    assert np.allclose(pt.superpose(bz, bx, bx)(x), [2.0, 0.0, 2.0])


def test_fields_together_add_and_skip_what_is_absent() -> None:
    x = np.zeros((2, 3))
    bz = pt.magnetic(pt.uniform(np.array([0.0, 0.0, 2.0])))
    ex: pt.Fields = lambda x, t: (np.full_like(x, t), None)  # noqa: E731
    e, b = pt.together(bz, ex, bz)(x, 3.0)
    assert e is not None
    assert b is not None
    assert np.allclose(e, 3.0)
    assert np.allclose(b, [0.0, 0.0, 4.0])
    assert pt.together(bz)(x, 0.0)[0] is None  # no electric field: nothing made up for it
    assert pt.together()(x, 0.0) == (None, None)


# ---- the integrator -------------------------------------------------------------------------------------------
def test_boris_gyration_radius_phase_and_energy() -> None:
    v, b = pt.v_of(100.0), 0.01
    gamma = 1 / math.sqrt(1 - (v / pt.C0) ** 2)
    period = 2 * math.pi * gamma / (pt.E / pt.ME * b)
    x, u = np.zeros((1, 3)), np.array([[gamma * v, 0.0, 0.0]])
    qm, ef, bf = np.array([pt.QM_E]), np.zeros((1, 3)), np.array([[0.0, 0.0, b]])
    ys = []
    for _ in range(40 * 50):
        x, u = pt.push(x, u, qm, period / 40, ef, bf)
        ys.append(float(x[0, 1]))
    assert max(ys) - min(ys) == pytest.approx(2 * pt.rigidity(100.0) / b, rel=0.01)
    assert math.atan2(u[0, 1], u[0, 0]) == pytest.approx(0.0, abs=1e-9)  # 50 whole turns: the rotation is exact
    assert float(np.linalg.norm(u)) == pytest.approx(gamma * v, rel=1e-12)


def test_boris_exb_drift_and_its_known_step_error() -> None:
    # Crossed fields: the drift is E/B along E x B. With the exact rotation the scheme's drift is low by the
    # factor (theta/2) / tan(theta/2), theta the gyro-angle per step: 0.2 % at 40 steps per turn.
    e0, b, steps, turns = 100.0, 0.01, 40, 20
    period = 2 * math.pi / (pt.E / pt.ME * b)
    x, u = np.zeros((1, 3)), np.zeros((1, 3))
    qm, ef, bf = np.array([pt.QM_E]), np.array([[e0, 0.0, 0.0]]), np.array([[0.0, 0.0, b]])
    for _ in range(steps * turns):
        x, u = pt.push(x, u, qm, period / steps, ef, bf)
    half = math.pi / steps
    assert x[0, 1] == pytest.approx(-(e0 / b) * turns * period * half / math.tan(half), rel=1e-6)
    assert x[0, 1] == pytest.approx(-(e0 / b) * turns * period, rel=3e-3)


def test_mathieu_stability_boundary_and_secular_frequency() -> None:
    # Ideal quadrupole E = g cos(Omega t) (-x, +y, 0), charge/mass 1: q = 2 g / Omega^2. Stable below q = 0.908;
    # for small q the slow motion has frequency beta Omega / 2 with beta = q / sqrt(2).
    om, qs = 2 * math.pi, np.array([0.2, 0.85, 1.0])
    g = qs * om**2 / 2
    x, u = np.tile([1e-3, 1e-3, 0.0], (3, 1)), np.zeros((3, 3))
    per, nper, hist = 100, 200, []
    for it in range(per * nper):
        ef = np.stack([-g * x[:, 0], g * x[:, 1], 0 * g], 1) * math.cos(om * it / per)
        x, u = pt.push(x, u, np.ones(3), 1 / per, ef, np.zeros((3, 3)))
        if it % per == 0:
            hist.append(x[:, 0].copy())
    h = np.array(hist)
    growth = np.abs(h[-50:]).max(0) / 1e-3
    assert growth[1] < 50
    assert growth[2] > 1e3
    up = [i + h[i, 0] / (h[i, 0] - h[i + 1, 0]) for i in range(len(h) - 1) if h[i, 0] < 0 <= h[i + 1, 0]]
    beta = 2 * (len(up) - 1) / (up[-1] - up[0])
    assert beta == pytest.approx(0.2 / math.sqrt(2), rel=0.02)


# ---- the general loop -----------------------------------------------------------------------------------------
def test_integrate_in_a_uniform_electric_field_gains_qe_per_second() -> None:
    e0, dt, steps = 50.0, 1e-9, 200
    field: pt.Fields = lambda x, t: (np.broadcast_to([e0, 0.0, 0.0], x.shape), None)  # noqa: E731
    flight = pt.integrate(np.zeros((3, 3)), np.zeros((3, 3)), np.full(3, pt.QM_E), dt, steps, field)
    assert flight.steps == steps
    assert bool(flight.alive.all())
    assert np.allclose(flight.u[:, 0], pt.QM_E * e0 * steps * dt, rtol=1e-12)
    assert np.allclose(flight.x[:, 0], 0.5 * pt.QM_E * e0 * (steps * dt) ** 2, rtol=0.01)  # non-relativistic here


def test_integrate_stops_what_an_observer_returns_and_ends_when_none_is_left() -> None:
    seen: list[tuple[int, int]] = []

    def far(step: int, t: float, x: pt.Arr, u: pt.Arr, alive: pt.BoolArr) -> pt.BoolArr:
        seen.append((step, int(alive.sum())))
        stop: pt.BoolArr = x[:, 0] > np.array([2.5e-3, 5.5e-3])
        return stop

    def count(step: int, t: float, x: pt.Arr, u: pt.Arr, alive: pt.BoolArr) -> None:
        assert t == pytest.approx((step + 1) * 1e-3)
        seen.append((-1, int(alive.sum())))  # called after `far`, it sees what `far` stopped

    still: pt.Fields = lambda x, t: (None, None)  # noqa: E731
    x0, u0 = np.zeros((2, 3)), np.tile([1.0, 0.0, 0.0], (2, 1))  # 1 m/s along x, one step a millimetre
    flight = pt.integrate(x0, u0, np.ones(2), 1e-3, 100, still, [far, count])
    assert flight.steps == 6  # the second particle passes 5.5 mm on its sixth step, and the run ends
    assert not flight.alive.any()
    assert np.allclose(flight.x[:, 0], [3e-3, 6e-3])  # each stays where it was stopped
    assert seen[:2] == [(0, 2), (-1, 2)]
    assert (2, 2) in seen
    assert (-1, 1) in seen


def test_integrate_refuses_arrays_of_the_wrong_shape() -> None:
    still: pt.Fields = lambda x, t: (None, None)  # noqa: E731
    with pytest.raises(ValueError, match="shape"):
        pt.integrate(np.zeros(3), np.zeros(3), np.ones(1), 1e-3, 1, still)
    with pytest.raises(ValueError, match="shape"):
        pt.integrate(np.zeros((2, 3)), np.zeros((2, 3)), np.ones(3), 1e-3, 1, still)


# ---- closed curves --------------------------------------------------------------------------------------------
def test_torus_knot_is_closed_and_evenly_spaced() -> None:
    knot = ClosedCurve.torus_knot(0.060, 0.028, 600)
    steps = np.linalg.norm(np.roll(knot.points, -1, 0) - knot.points, axis=1)
    assert len(knot) == 600
    assert float(steps.std() / steps.mean()) < 1e-3
    assert float(steps.sum()) == pytest.approx(knot.length, rel=1e-3)


def test_frame_is_orthonormal_and_closes() -> None:
    knot = ClosedCurve.torus_knot(0.060, 0.028, 600)
    T, fr = knot.tangents, knot.frame(2 * math.pi)
    for a, b in ((T, fr.a), (T, fr.b), (fr.a, fr.b)):
        assert float(np.abs(np.sum(a * b, 1)).max()) < 1e-9
    for a in (T, fr.a, fr.b):
        assert np.allclose(np.linalg.norm(a, axis=1), 1.0, atol=1e-12)
    assert float(np.sum(fr.a[-1] * fr.a[0])) > 0.99  # the last station joins the first


@pytest.mark.parametrize(("p", "q", "r"), [(2, 3, 0.028), (2, 3, 0.016), (2, 5, 0.020), (3, 2, 0.020)])
def test_holonomy_is_the_writhe_modulo_one_turn(p: int, q: int, r: float) -> None:
    # Calugareanu-White-Fuller: a normal carried round a closed curve without twisting comes back turned by
    # 2 pi Wr, modulo 2 pi. Two independent computations: transport along the curve, and the Gauss double integral.
    def gap(n: int) -> float:
        knot = ClosedCurve.torus_knot(0.060, r, n, p, q)
        turns = knot.holonomy / (2 * math.pi) - knot.writhe
        return abs(turns - round(turns))

    assert gap(600) < 1.5e-3
    assert gap(1200) < 0.4 * gap(600)  # second order in the spacing of the points


def test_figure_eight_is_closed_unknotted_and_turns_the_plane_only_when_lifted() -> None:
    flat, lifted = ClosedCurve.figure_eight(0.14, 0.07, 0.0, 600), ClosedCurve.figure_eight(0.14, 0.07, 0.04, 600)
    assert abs(flat.points[:, 2]).max() == 0.0
    assert flat.holonomy == pytest.approx(0.0, abs=1e-9)
    assert abs(lifted.holonomy) > 1.0
    assert abs(lifted.writhe) < 1.0  # under one crossing's worth: a knot would need three
    gap = lifted.holonomy / (2 * math.pi) - lifted.writhe
    assert abs(gap - round(gap)) < 1.5e-3


def test_axes_are_orthonormal_and_come_back_turned_by_the_holonomy() -> None:
    knot = ClosedCurve.torus_knot(0.060, 0.028, 600)
    (e1, e2), T = knot.axes, knot.tangents
    assert float(np.abs(np.sum(e1 * T, 1)).max()) < 1e-9
    assert float(np.abs(np.sum(e1 * e2, 1)).max()) < 1e-9
    # the last point's axes, against the first's: all but one step of the lap's turn
    seen = math.atan2(float(e1[-1] @ e2[0]), float(e1[-1] @ e1[0]))
    assert seen == pytest.approx(knot.holonomy, abs=0.05)


def test_planar_curve_has_no_holonomy_and_no_writhe() -> None:
    ring = ClosedCurve.circle(0.1, 400)
    assert ring.holonomy == pytest.approx(0.0, abs=1e-9)
    assert ring.writhe == pytest.approx(0.0, abs=1e-9)
    assert ring.bend_radius == pytest.approx(0.1, rel=1e-3)


CURVES = {
    "ring": ClosedCurve.circle(0.1484, 600),
    "trefoil": ClosedCurve.torus_knot(0.060, 0.028, 600),
    "shallow trefoil": ClosedCurve.torus_knot(0.060, 0.016, 600),
    "helix": ClosedCurve.torus_knot(0.110, 0.030, 600, 1, 3),
    "figure-eight": ClosedCurve.figure_eight(0.14, 0.07, 0.04, 600),
}


def test_cur_2_bending_is_at_least_the_total_curvature_squared_over_the_length() -> None:
    for name, curve in CURVES.items():
        floor = curve.total_curvature**2 / curve.length
        assert curve.bending >= floor * (1 - 1e-12), name
    ring = CURVES["ring"]  # constant curvature: the floor is reached
    assert ring.bending == pytest.approx(ring.total_curvature**2 / ring.length, rel=1e-9)
    assert ring.bending == pytest.approx(2 * math.pi / 0.1484, rel=1e-3)


def test_cur_3_total_curvature_is_two_pi_on_a_ring_and_past_four_pi_on_a_knot() -> None:
    assert CURVES["ring"].total_curvature == pytest.approx(2 * math.pi, rel=1e-4)
    for name, curve in CURVES.items():
        assert curve.total_curvature >= 2 * math.pi * (1 - 1e-4), name  # Fenchel
    for name in ("trefoil", "shallow trefoil"):  # Fary and Milnor: so its bending is past 16 pi^2 / length
        knot = CURVES[name]
        assert knot.total_curvature > 4 * math.pi
        assert knot.bending > 16 * math.pi**2 / knot.length
    assert CURVES["figure-eight"].total_curvature < 4 * math.pi  # an unknot may stay under the knots' floor


def test_trefoil_has_writhe_of_magnitude_between_three_and_four() -> None:
    assert 3.0 < abs(ClosedCurve.torus_knot(0.060, 0.028, 600).writhe) < 4.0


# ---- coils ----------------------------------------------------------------------------------------------------
def test_field_at_the_centre_of_a_loop() -> None:
    ring = ClosedCurve.circle(0.1, 360)
    nxt = np.roll(ring.points, -1, 0)
    loop = Coil(ring, 0.0, 1, "single", (ring.points + nxt) / 2, nxt - ring.points)
    b = loop.field(np.zeros((1, 3)))
    assert b[0, 2] == pytest.approx(pt.MU0 / (2 * 0.1), rel=1e-3)
    assert float(np.abs(b[0, :2]).max()) < 1e-12


def test_single_winding_carries_one_net_turn_and_return_winding_none() -> None:
    ring = ClosedCurve.circle(0.1, 300)
    radius, turns = 0.011, 90
    single = Coil.wound(ring, radius, turns, "single")
    back = Coil.wound(ring, radius, turns, "return")
    assert abs(single.linked_current(2.5 * radius)) == pytest.approx(1.0, abs=0.02)
    assert abs(back.linked_current(2.5 * radius)) < 0.02
    # both layers of the return winding add to the field along the curve, which points along +T
    assert back.turns == 2 * turns
    assert back.per_amp == pytest.approx(pt.MU0 * back.turns / ring.length, rel=0.01)
    assert float(np.sum(back.axis_field * ring.tangents, 1).min()) > 0
    assert float(np.abs(back.slope).max()) < 1e-3


def test_return_winding_leaves_a_small_dipole_from_its_two_radii() -> None:
    # The two layers carry opposite currents along the curve, at tube radii a and a + gap. On a ring of radius R
    # their fields through the hole differ by (mu0 I / 4 pi R) ln(1 + gap / a): a fraction ln(1 + gap / a) / (2 N)
    # of the guide field, for N turns in all.
    a, gap = 0.015, 1e-3
    coil = Coil.wound(ClosedCurve.circle(0.08, 240), a, 60, "return", per_turn=18, layer_gap=gap)
    assert abs(float(np.mean(coil.slope[:, 2]))) == pytest.approx(math.log(1 + gap / a) / (2 * coil.turns), rel=0.05)


def test_fast_field_sum_is_the_plain_sum() -> None:
    # COI-3. Points all over the bore of a knotted winding, out to 1 mm from the wire, one block and several.
    knot = ClosedCurve.torus_knot(0.060, 0.028, 300)
    coil = Coil.wound(knot, 0.011, 60, "return", per_turn=18)
    rng, fr = np.random.default_rng(7), knot.frame(2 * math.pi)
    at = rng.integers(0, len(knot), 50)
    r, phi = 0.010 * np.sqrt(rng.uniform(0, 1, 50)), rng.uniform(0, 2 * math.pi, 50)
    x = knot.points[at] + (r * np.cos(phi))[:, None] * fr.a[at] + (r * np.sin(phi))[:, None] * fr.b[at]
    plain = coil.field_plain(x)
    scale = float(np.abs(plain).max())
    assert float(np.abs(coil.field(x) - plain).max()) < 1e-11 * scale
    several = np.vstack([coil.field(x[i : i + 7]) for i in range(0, 50, 7)])
    assert float(np.abs(several - plain).max()) < 1e-11 * scale
    assert coil.field(x[:1]).shape == (1, 3)


def test_field_takes_points_a_block_at_a_time(monkeypatch: pytest.MonkeyPatch) -> None:
    from array16.particles import coils

    ring = ClosedCurve.circle(0.1, 120)
    coil = Coil.wound(ring, 0.011, 30, "single", per_turn=12)
    x = ring.points[::5] + 0.003
    whole = coil.field(x)
    monkeypatch.setattr(coils, "_BLOCK", 5 * len(coil.mid))  # five points at a time
    # The same to rounding, and not to the bit: a matrix product may sum in another order for another shape
    # (it does under OpenBLAS on Linux, and does not under Accelerate on a Mac).
    assert float(np.abs(coil.field(x) - whole).max()) < 1e-13 * float(np.abs(whole).max())


def test_a_winding_wider_than_the_bend_is_refused() -> None:
    knot = ClosedCurve.torus_knot(0.060, 0.028, 600)  # tightest bend 43.8 mm
    with pytest.raises(ValueError, match="cannot go round a bend"):
        Coil.wound(knot, 0.044, 10, "single")
    with pytest.raises(ValueError, match="cannot go round a bend"):
        Coil.wound(knot, 0.0432, 10, "return")  # the second layer, 1 mm further out, is the one that crosses
    assert Coil.wound(knot, 0.0432, 10, "single").turns == 10


# ---- first-order guide optics ---------------------------------------------------------------------------------
def test_ring_guide_drifts_two_pi_rho_per_lap_and_never_closes() -> None:
    opt = GuideOptics.of(ClosedCurve.circle(0.1, 600))
    assert opt.alpha == pytest.approx(0.0, abs=1e-9)
    assert not opt.closes
    # axes at the first point: e1 radially outward, e2 = T x e1 = -z. An electron moving along the field drifts
    # towards -z by 2 pi rho per lap.
    assert opt.g[-1] == pytest.approx(2j * math.pi, rel=1e-3)
    rho = np.array([1e-3])
    assert bool(np.isnan(opt.closed_orbit(rho)[0]))
    assert bool(np.isinf(opt.matched_excursion(rho)[0]))
    assert opt.excursion(rho, 0j, laps=5)[0] == pytest.approx(5 * 2 * math.pi * 1e-3, rel=1e-3)


def test_transverse_maps_offsets_onto_the_twist_free_axes() -> None:
    ring = ClosedCurve.circle(0.1, 600)
    v = ring.transverse(np.array([1 + 0j, 1j]))
    assert np.allclose(v[0], [1.0, 0.0, 0.0], atol=1e-9)  # e1: radially outward at the first point
    assert np.allclose(v[1], [0.0, 0.0, -1.0], atol=1e-9)  # e2 = T x e1


def test_knot_guide_has_a_closed_orbit_that_scales_with_rho() -> None:
    opt = GuideOptics.of(ClosedCurve.torus_knot(0.060, 0.028, 600))
    assert opt.closes
    assert abs(opt.transform) > 0.05
    rho = np.array([0.5e-3, 1e-3])
    z = opt.closed_orbit(rho)
    assert z[1] == pytest.approx(2 * z[0], rel=1e-12)
    matched, on_axis = opt.matched_excursion(rho), opt.excursion(rho, 0j, laps=30)
    assert bool(np.all(matched < on_axis))
    assert matched[1] == pytest.approx(2 * matched[0], rel=1e-9)


def test_knot_guide_optics_converge_with_the_number_of_points() -> None:
    coarse, fine, finest = (GuideOptics.of(ClosedCurve.torus_knot(0.060, 0.028, n)) for n in (600, 1200, 2400))
    assert coarse.alpha == pytest.approx(fine.alpha, abs=3e-3)
    assert abs(fine.alpha - finest.alpha) < 0.4 * abs(coarse.alpha - fine.alpha)  # second order
    assert abs(coarse.g[-1]) == pytest.approx(abs(fine.g[-1]), rel=5e-3)
    assert abs(fine.g[-1] - finest.g[-1]) < 0.4 * abs(coarse.g[-1] - fine.g[-1])


def test_single_winding_on_a_ring_bends_exactly_one_momentum() -> None:
    # The one net turn makes a field through the hole, B_z. The lap's step vanishes when the field line's rise
    # cancels the curvature drift: rho = (B_z / B_t) R, which is the bending condition B_z R = p / e.
    radius = 0.1
    ring = ClosedCurve.circle(radius, 300)
    coil = Coil.wound(ring, 0.011, 90, "single")
    opt = coil.optics
    rho_match = -opt.h[-1] / opt.g[-1]
    along = np.sum(coil.axis_field * ring.tangents, 1)
    assert abs(rho_match.imag) < 1e-3 * abs(rho_match.real)
    assert abs(rho_match.real) == pytest.approx(abs(float(np.mean(coil.axis_field[:, 2] / along))) * radius, rel=0.01)
    assert abs(opt.step(np.array([rho_match.real]))[0]) < 1e-3 * abs(opt.h[-1])


# ---- laws of the lap map (formal/Array16/GuideOptics.lean), attacked on the floats ------------------------------
turning = st.floats(0.05, 2 * math.pi - 0.05)
offsets = st.complex_numbers(max_magnitude=5.0, allow_nan=False, allow_infinity=False)
rhos = st.floats(-3.0, 3.0)


def _optics(alpha: float, g_end: complex, h_end: complex) -> GuideOptics:
    """Optics whose drift and wander accumulate along a bent path to the given end values."""
    t = np.linspace(0, 1, 9)
    bend = np.sin(math.pi * t) * (0.3 + 0.2j)
    return GuideOptics(alpha, t * g_end + bend, t * h_end - bend)


def _tol(*values: complex) -> float:
    return 1e-9 * (1 + sum(abs(v) for v in values))


@settings(max_examples=200, deadline=None)
@given(turning, offsets, offsets, rhos)
def test_gde_1_the_closed_orbit_is_a_fixed_point(alpha: float, g_end: complex, h_end: complex, rho: float) -> None:
    opt, r = _optics(alpha, g_end, h_end), np.array([rho])
    z = opt.closed_orbit(r)
    assert abs(opt.lap(r, z)[0] - z[0]) <= 20 * _tol(z[0], opt.step(r)[0])


@settings(max_examples=200, deadline=None)
@given(turning, offsets, offsets, rhos)
def test_gde_1_control_a_nearby_formula_is_not_a_fixed_point(
    alpha: float, g_end: complex, h_end: complex, rho: float
) -> None:
    # The formula without the turn in its numerator, step / (1 - rot), sits exactly one step from the closed orbit,
    # so a lap moves it by |1 - rot| times the step: the test above would catch it.
    opt, r = _optics(alpha, g_end, h_end), np.array([rho])
    step, rot = opt.step(r)[0], cmath.exp(1j * alpha)
    wrong = np.array([step / (1 - rot)])
    tol = 40 * _tol(wrong[0], step)
    assert abs(abs(wrong[0] - opt.closed_orbit(r)[0]) - abs(step)) <= tol
    assert abs(abs(opt.lap(r, wrong)[0] - wrong[0]) - abs(1 - rot) * abs(step)) <= tol


@settings(max_examples=200, deadline=None)
@given(turning, offsets, offsets, rhos, offsets)
def test_gde_2_only_the_closed_orbit_closes(
    alpha: float, g_end: complex, h_end: complex, rho: float, z0: complex
) -> None:
    # What a lap moves an offset by is |1 - exp(i alpha)| times its distance from the closed orbit: zero only there.
    opt, r = _optics(alpha, g_end, h_end), np.array([rho])
    star = opt.closed_orbit(r)[0]
    moved = abs(opt.lap(r, np.array([z0]))[0] - z0)
    assert abs(moved - abs(1 - cmath.exp(1j * alpha)) * abs(z0 - star)) <= 40 * _tol(z0, star)


@settings(max_examples=200, deadline=None)
@given(turning, offsets, offsets, rhos, offsets)
def test_gde_3_4_5_laps_turn_the_offset_about_the_closed_orbit(
    alpha: float, g_end: complex, h_end: complex, rho: float, z0: complex
) -> None:
    opt, r = _optics(alpha, g_end, h_end), np.array([rho])
    star = opt.closed_orbit(r)[0]
    z = np.array([z0])
    for lap in range(1, 40):
        z = opt.lap(r, z)
        # GDE-3: after n laps the offset from the closed orbit has turned by n alpha; GDE-4: its length is kept
        assert abs((z[0] - star) - cmath.exp(1j * alpha * lap) * (z0 - star)) <= 40 * _tol(z0, star)
    # GDE-5: so the excursion is bounded for any number of laps
    within = float(np.abs(r[0] * opt.g + opt.h).max())
    assert opt.excursion(r, z0, laps=40)[0] <= abs(star) + abs(z0 - star) + within + 40 * _tol(z0, star)


@settings(max_examples=200, deadline=None)
@given(offsets, offsets, rhos, offsets)
def test_gde_6_without_a_turn_every_lap_adds_the_same_step(
    g_end: complex, h_end: complex, rho: float, z0: complex
) -> None:
    opt, r = _optics(0.0, g_end, h_end), np.array([rho])
    z = np.array([z0])
    for _ in range(25):
        z = opt.lap(r, z)
    assert abs(z[0] - (z0 + 25 * opt.step(r)[0])) <= 25 * _tol(z0, opt.step(r)[0])


@settings(max_examples=200, deadline=None)
@given(st.floats(0.0, 2 * math.pi), offsets, offsets, rhos, rhos, offsets, st.floats(0.0, 1.0))
def test_gde_7_reach_is_convex_in_rho(
    alpha: float, g_end: complex, h_end: complex, rho_a: float, rho_b: float, z0: complex, t: float
) -> None:
    opt = _optics(alpha, g_end, h_end)
    rho = np.array([rho_a, rho_b, (1 - t) * rho_a + t * rho_b])
    reach = opt.excursion(rho, z0, laps=12)
    assert reach[2] <= (1 - t) * reach[0] + t * reach[1] + 1e-9 * (1 + float(reach.max()))
    if opt.closes and abs(1 - cmath.exp(1j * alpha)) > 0.05:
        matched = opt.matched_excursion(rho)
        assert matched[2] <= (1 - t) * matched[0] + t * matched[1] + 1e-7 * (1 + float(matched.max()))


@settings(max_examples=200, deadline=None)
@given(offsets, rhos, rhos, st.floats(0.05, 0.95))
def test_gde_7_control_reach_is_not_concave(g_end: complex, rho_a: float, rho_b: float, t: float) -> None:
    # The opposite inequality fails: with no wander and an on-axis launch the reach is |rho| times a constant, a
    # V in rho, and the chord across rho = 0 lies above it.
    if abs(g_end) < 0.1 or rho_a * rho_b >= -1e-3:
        return
    opt = _optics(1.0, g_end, 0j)
    rho = np.array([rho_a, rho_b, (1 - t) * rho_a + t * rho_b])
    reach = opt.excursion(rho, 0j, laps=6)
    assert reach[2] < (1 - t) * reach[0] + t * reach[1]


@settings(max_examples=100, deadline=None)
@given(st.floats(0.3, 2 * math.pi - 0.3), offsets, offsets, rhos, offsets)
def test_orbit_circle_recovers_the_closed_orbit_and_the_turn(
    alpha: float, g_end: complex, h_end: complex, rho: float, z0: complex
) -> None:
    # GDE-3 and GDE-4 run backwards: from the offsets at the start of 15 laps, find the centre and the turn.
    opt, r = _optics(alpha, g_end, h_end), np.array([rho])
    star = opt.closed_orbit(r)[0]
    if abs(z0 - star) < 1e-3:
        return
    z, laps = np.array([z0]), [z0]
    for _ in range(15):
        z = opt.lap(r, z)
        laps.append(complex(z[0]))
    fit = pt.orbit_circle(np.array(laps))
    scale = 1 + abs(star) + abs(z0)
    assert abs(fit.centre - star) <= 1e-6 * scale
    assert fit.radius == pytest.approx(abs(z0 - star), abs=1e-6 * scale)
    assert abs(cmath.exp(1j * fit.turn) - cmath.exp(1j * alpha)) <= 1e-6
    assert fit.residual <= 1e-6 * scale


@settings(max_examples=100, deadline=None)
@given(st.floats(0.0, 2 * math.pi), offsets, offsets, rhos, offsets)
def test_unlap_undoes_lap(alpha: float, g_end: complex, h_end: complex, rho: float, z0: complex) -> None:
    opt, r, z = _optics(alpha, g_end, h_end), np.array([rho]), np.array([z0])
    assert abs(opt.unlap(r, opt.lap(r, z))[0] - z0) <= 10 * _tol(z0, opt.step(r)[0])
    assert abs(opt.lap(r, opt.unlap(r, z))[0] - z0) <= 10 * _tol(z0, opt.step(r)[0])


@settings(max_examples=100, deadline=None)
@given(turning, offsets, offsets, st.floats(0.01, 3.0), offsets)
def test_distances_along_the_field_are_the_excursion_and_join_lap_to_lap(
    alpha: float, g_end: complex, h_end: complex, rho: float, z0: complex
) -> None:
    opt = _optics(alpha, g_end, h_end)
    for signed in (rho, -rho):
        far = opt.distances(signed, z0, laps=7)
        assert far.shape == (7, len(opt.g))
        assert far[0, 0] == pytest.approx(abs(z0), abs=1e-9)  # it starts where it was launched
        assert np.allclose(far[:-1, -1], far[1:, 0], atol=1e-8 * (1 + far.max()))  # a lap ends where the next starts
        assert np.allclose(far[:, -1], np.abs(opt.crossings(signed, z0, laps=7)), atol=1e-8 * (1 + far.max()))
    # along the field it is the reach the laws are about
    assert opt.distances(rho, z0, laps=7).max() == pytest.approx(opt.excursion(np.array([rho]), z0, laps=7)[0])


def test_both_directions_share_a_closed_orbit_and_turn_opposite_ways() -> None:
    opt = _optics(2.0, 1.5 + 0.5j, 0.1j)
    for signed in (0.8, -0.8):
        star = complex(opt.closed_orbit(np.array([signed]))[0])
        assert np.allclose(opt.crossings(signed, star, laps=5), star)  # launched on it, it stays on it
        matched = opt.matched_excursion(np.array([signed]))[0]
        assert opt.distances(signed, star, laps=3).max() == pytest.approx(matched)
        fit = pt.orbit_circle(np.concatenate([[0j], opt.crossings(signed, 0j, laps=12)]))
        assert abs(fit.centre - star) < 1e-9
        assert fit.turn == pytest.approx(math.copysign(2.0, signed))


def test_exit_lap_on_a_ring_is_the_wall_over_the_drift_per_lap() -> None:
    opt, rho, wall = GuideOptics.of(ClosedCurve.circle(0.1, 600)), 1e-3, 9.9e-3
    for signed in (rho, -rho):  # the drift does not depend on the direction of travel
        assert opt.exit_lap(signed, 0j, wall, laps=5) == pytest.approx(wall / (2 * math.pi * rho), abs=2e-3)
        assert complex(opt.crossings(signed, 0j, laps=1)[0]) == pytest.approx(2j * math.pi * rho, rel=1e-3)
    assert opt.exit_lap(rho, 0j, wall, laps=1) is None  # not yet, after one lap


def test_window_is_the_range_of_rho_inside_the_wall() -> None:
    rho = np.linspace(-2.0, 2.0, 401)
    assert pt.window(rho, np.abs(rho - 0.5), 1.0) == pytest.approx((-0.49, 1.49))
    assert pt.window(rho, np.abs(rho) + 3.0, 1.0) is None


# ---- tracking, end to end -------------------------------------------------------------------------------------
def test_launch_shares_one_offset_and_defaults_to_moving_along_the_field() -> None:
    launch = Launch.of(np.array([10.0, 20.0, 30.0]), np.array([1e-3, 0.0, 0.0]))
    assert launch.offset.shape == (3, 3)
    assert np.allclose(launch.sense, 1.0)


def test_tracked_electron_in_a_ring_coil_drifts_as_the_theory_says() -> None:
    # Biot-Savart field of a discrete return-wound coil, full-orbit integration, one lap, against GuideOptics.
    radius, bt, rho = 0.08, 0.02, 0.4e-3
    ring = ClosedCurve.circle(radius, 240)
    coil = Coil.wound(ring, 0.015, 60, "return", per_turn=18)
    launch = Launch.of(np.array([pt.energy_of(bt * rho)]), np.zeros(3))
    run = pt.track_guide(
        ring, pt.scaled(coil.field, bt / coil.per_amp), bt, launch, 0.9 * coil.radius, 2.02, record=0, record_every=1
    )
    assert bool(run.held[0])
    assert run.laps[0] == pytest.approx(2.02, abs=0.02)
    # the prediction: curvature drift, 2 pi rho per lap towards -z, less the rise of the winding's residual dipole
    step = coil.optics.step(np.array([rho]))
    predicted = ring.transverse(step * run.laps[0])[0]
    assert predicted[2] == pytest.approx(-2 * math.pi * rho * run.laps[0], rel=0.08)
    assert run.path[-1, 2] == pytest.approx(predicted[2], rel=0.03)
    assert run.reach[0] == pytest.approx(abs(predicted[2]), rel=0.03)
    # and where it crossed the start after each lap: one step, then two, on the curve's axes
    assert len(run.section[0]) == 2
    assert abs(run.section[0][0] - step[0]) < 0.03 * abs(step[0])
    assert abs(run.section[0][1] - 2 * step[0]) < 0.03 * abs(2 * step[0])
