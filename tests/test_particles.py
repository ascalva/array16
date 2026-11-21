"""Calibration of the particle kernel against results known in closed form, and attacks on the laws of the guide
optics. Each test named test_<law>_N is a falsifier of law <LAW>-N in docs/LAWS.md; a test named test_gde_N attacks
the float implementation of the law proved as GDE-N in formal/."""

from __future__ import annotations

import cmath
import heapq
import itertools
import math
from collections.abc import Callable

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


def test_guide_watch_times_laps_and_stops_what_passes_the_wall() -> None:
    # The watch alone, shown two points moving round a circle: one on it, one drifting outward.
    ring = ClosedCurve.circle(0.1, 360)
    watch = pt.GuideWatch(ring, wall=5e-3, count=2)
    alive = np.ones(2, dtype=np.bool_)
    period, steps = 2.0, 1000
    for step in range(2 * steps + 30):
        t = (step + 1) * period / steps
        angle = 2 * math.pi * t / period
        radius = np.array([0.1, 0.1 + 2e-3 * t])  # the second leaves at 1 mm a lap
        x = np.stack([radius * math.cos(angle), radius * math.sin(angle), np.zeros(2)], 1)
        alive &= ~watch(step, t, x, np.zeros((2, 3)), alive)
    assert watch.times[0] == pytest.approx([period, 2 * period], abs=2 * period / 360)
    assert watch.laps[0] == pytest.approx(2.03, abs=0.01)
    # a point on the circle lies off the tangent at the nearest station by up to step^2 / (8 R): the watch's grain
    grain = (ring.length / len(ring)) ** 2 / (8 * 0.1)
    assert abs(watch.crossings[0][0]) <= grain
    assert watch.reach[0] <= grain
    assert bool(alive[0])
    assert not alive[1]  # it passed 5 mm during its second lap
    assert watch.reach[1] == pytest.approx(5e-3, abs=1e-4)
    assert len(watch.times[1]) == 1
    assert watch.crossings[1][0] == pytest.approx(4e-3, abs=1e-4)  # 4 mm outward, along the first axis


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


# ---- radio-frequency confinement ------------------------------------------------------------------------------
def test_pseudo_potential_numbers_for_an_electron() -> None:
    assert pt.ponderomotive_ev(pt.field_for_depth(1.0, 0.4e9), 0.4e9) == pytest.approx(1.0)
    assert pt.field_for_depth(1.0, 0.4e9) == pytest.approx(11.99e3, rel=1e-3)  # V/m
    assert pt.field_for_depth(1.0, 2.45e9) == pytest.approx(73.4e3, rel=1e-3)
    assert pt.field_for_depth(0.025, 5.8e9) == pytest.approx(27.5e3, rel=1e-3)
    # half the frequency, four times the depth; and a proton, 1836 times heavier, gets 1836 times less
    assert pt.ponderomotive_ev(1e3, 1e9) == pytest.approx(4 * pt.ponderomotive_ev(1e3, 2e9))
    assert pt.ponderomotive_ev(1e3, 1e9, pt.QM_E / 1836.15) == pytest.approx(pt.ponderomotive_ev(1e3, 1e9) / 1836.15)


def test_quadrupole_numbers_agree_with_each_other() -> None:
    f, r0 = 1e9, 0.5e-3
    for q in (0.3, 0.6):
        g = pt.gradient_for_q(q, f)
        assert pt.mathieu_q(g, f) == pytest.approx(q)
        assert pt.quadrupole_depth_ev(q, f, r0) == pytest.approx(q * g * r0**2 / 8)  # the textbook q V0 / 8
    assert pt.secular_frequency(0.6, f) == pytest.approx(212.1e6, rel=1e-3)
    assert pt.quadrupole_depth_ev(0.6, f, r0) == pytest.approx(1.263, rel=1e-3)
    assert pt.bend_limit_ev(1.263, 16.6e-3, r0) == pytest.approx(20.97, rel=1e-3)
    # the beam's own field at the limiting density equals the guide's restoring field, at any radius
    n, w_sec = pt.space_charge_limit(0.6, f), 2 * math.pi * pt.secular_frequency(0.6, f)
    assert n * pt.E / (2 * pt.EPS0) == pytest.approx(w_sec**2 / abs(pt.QM_E))
    assert n == pytest.approx(1.12e15, rel=0.01)


@settings(max_examples=200, deadline=None)
@given(st.floats(0.01, 1e4), st.floats(0.05, 0.9), st.floats(1e-4, 1e-2), st.sampled_from([1.0, 1836.15, 1e9]))
def test_rfg_1_volts_set_the_depth_and_frequency_times_distance_follows(
    depth: float, q: float, r0: float, mass: float
) -> None:
    # Whatever the particle and the size: a well `depth` deep at Mathieu's q takes V0 = 8 depth / q, and that
    # drive gives that q where f r0 = sqrt(2 |q/m| V0 / q) / 2 pi.
    qm = pt.QM_E / mass
    volts = pt.volts_for_depth(depth, q)
    f = pt.frequency_times_distance(volts, q, qm) / r0
    assert pt.mathieu_q(volts / r0**2, f, qm) == pytest.approx(q, rel=1e-9)
    assert pt.quadrupole_depth_ev(q, f, r0, qm) == pytest.approx(depth, rel=1e-9)


def test_rfg_1_numbers_for_an_electron() -> None:
    assert pt.volts_for_depth(1.0, 0.3) == pytest.approx(26.7, rel=1e-2)
    assert pt.frequency_times_distance(pt.volts_for_depth(1.0, 0.3), 0.3) == pytest.approx(0.891e6, rel=1e-2)  # GHz mm


def test_rod_1_four_thin_rods_make_a_quadrupole_of_the_stated_gradient() -> None:
    # Across the axis the field changes as g (-a, +b) on the frame's axes, g = 2 lam / (pi eps0 r0^2), with
    # nothing along the axis: by central differences, on rings of several radii.
    lam = 1e-9
    for radius, r0, points in ((0.2, 5e-3, 2000), (0.1, 5e-3, 2000), (0.2, 2e-3, 4000)):
        ring = ClosedCurve.circle(radius, points)
        rods, delta = pt.QuadRods.along(ring, r0), 0.04 * r0
        a, b, T = rods.frame.a[0], rods.frame.b[0], ring.tangents[0]
        e = rods.field_all(ring.points[0] + delta * np.array([a, -a, b, -b]), lam)
        g = rods.gradient(lam)
        assert rods.line_charge(g) == pytest.approx(lam)
        assert -(e[0] - e[1]) @ a / (2 * delta) == pytest.approx(g, rel=2e-3)
        assert +(e[2] - e[3]) @ b / (2 * delta) == pytest.approx(g, rel=2e-3)
        assert float(np.abs(e @ T).max()) < 1e-6 * g * delta
    assert pt.QuadRods.along(ring, r0).volts(lam, 25e-6) == pytest.approx(lam * math.log(80) / (math.pi * pt.EPS0))


def test_rod_2_a_bend_moves_the_null_inward_by_a_quarter_of_r0_squared_over_the_radius() -> None:
    # Bent to radius R the rods leave a field lam / (2 pi eps0 R) on their centre line, toward the centre of
    # the bend, whatever r0. With ROD-1's gradient the null lies r0^2 / (4 R) inside the centre line.
    lam = 1e-9
    for radius, r0, points in ((0.2, 5e-3, 2000), (0.1, 5e-3, 2000), (0.2, 2e-3, 4000), (0.05, 5e-3, 2000)):
        ring = ClosedCurve.circle(radius, points)
        rods = pt.QuadRods.along(ring, r0)
        outward = rods.frame.a[0]  # on a circle the frame's first normal points away from the centre
        on_axis = rods.field_all(ring.points[:1], lam)[0]
        assert on_axis @ outward == pytest.approx(-lam / (2 * math.pi * pt.EPS0 * radius), rel=0.015)
        assert abs(on_axis @ rods.frame.b[0]) < 1e-9 * abs(on_axis @ outward)
        null = (on_axis @ outward) / rods.gradient(lam)
        assert null == pytest.approx(-(r0**2) / (4 * radius), rel=0.015)


def test_rods_whole_field_is_the_sum_over_every_station() -> None:
    knot = ClosedCurve.torus_knot(0.060, 0.016, 301)
    rods = pt.QuadRods.along(knot, 6e-3)
    rng = np.random.default_rng(3)
    at = rng.integers(0, len(knot), 20)
    x = knot.points[at] + 3e-3 * rng.uniform(-1, 1, (20, 1)) * rods.frame.a[at]
    lam = rng.uniform(0.5, 2.0, 20) * 1e-9
    every = np.arange(len(knot)) - len(knot) // 2
    by_station = rods.field(x, at.astype(np.intp), lam, every.astype(np.intp))
    assert float(np.abs(rods.field_all(x, lam) - by_station).max()) < 1e-11 * float(np.abs(by_station).max())


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


# ---- the depth a watt buys (rf.py) ----------------------------------------------------------------------------
masses = st.sampled_from([1.0, 1836.15, 1e9])  # in electron masses: an electron, a proton, a charged grain
frequencies = st.floats(1e7, 1e11)  # Hz
# (depth_share, gradient_share): four rods, two facing boards, one board at its best, a rail pair off it
shares = st.sampled_from([(1.0, 1.0), (0.25, math.pi / 4), (0.009136, 0.3093), (0.0031, 0.2104)])


def test_rfg_2_the_ceiling_is_78_micro_ev_per_watt_for_an_electron() -> None:
    assert pt.free_space_depth_ceiling_ev(1.0) == pytest.approx(78.2e-6, rel=1e-3)
    assert 1.0 / pt.free_space_depth_ceiling_ev(1.0) == pytest.approx(12.78e3, rel=1e-3)  # watts for 1 eV
    assert pt.free_space_depth_ceiling_ev(1.0, pt.QM_E / 1836.15) == pytest.approx(78.2e-6 / 1836.15, rel=1e-3)
    assert pt.dipole_wave_peak_field(1.0, 1e9) == pytest.approx(265.0, rel=1e-3)  # V/m, a watt at 1 GHz


@settings(max_examples=200, deadline=None)
@given(st.floats(1e-6, 1e6), st.floats(1e6, 1e12), masses)
def test_rfg_2_the_dipole_waves_depth_is_the_ceiling_at_every_frequency(power: float, f: float, mass: float) -> None:
    qm = pt.QM_E / mass
    field = pt.dipole_wave_peak_field(power, f)
    assert pt.ponderomotive_ev(field, f, qm) == pytest.approx(pt.free_space_depth_ceiling_ev(power, qm), rel=1e-12)


def _dipole_shapes(p: pt.CArr, x: pt.Arr) -> tuple[pt.Arr, pt.Arr, pt.CArr, pt.CArr, pt.CArr]:
    """At points x (N, 3) round a dipole p (C m) at the origin: the distance r, the direction n, and the shapes
    of the dipole's fields, (n x p) x n = p - n (n . p), 3 n (n . p) - p, and n x p."""
    r = np.linalg.norm(x, axis=1)
    n = x / r[:, None]
    along = (n @ p)[:, None]
    return r, n, p[None, :] - n * along, 3 * n * along - p[None, :], np.cross(n, p)


def _outgoing_power(p: pt.CArr, k: float, radius: float) -> float:
    """The power (W, time-averaged) the outgoing wave of a dipole p (C m) carries through a sphere of the given
    radius: (1/2) Re(E x H*) . n over the sphere, with the dipole's whole fields, near ones included (Jackson
    9.18-9.19). Gauss-Legendre in cos(theta) and the trapezoid rule in phi are exact for their angular terms."""
    cos_t, weights = np.polynomial.legendre.leggauss(12)
    phi = 2 * np.pi * np.arange(24) / 24
    c, ph = np.meshgrid(cos_t, phi, indexing="ij")
    s = np.sqrt(1 - c**2)
    unit = np.stack([s * np.cos(ph), s * np.sin(ph), c], -1).reshape(-1, 3)
    area = np.repeat(weights * 2 * np.pi / len(phi), len(phi)) * radius**2
    r, n, radiation, near, turning = _dipole_shapes(p, radius * unit)
    wave = (np.exp(1j * k * r) / r)[:, None]
    e = (k**2 * radiation + near * (1 / r**2 - 1j * k / r)[:, None]) * wave / (4 * math.pi * pt.EPS0)
    h = pt.C0 * k**2 / (4 * math.pi) * turning * wave * (1 - 1 / (1j * k * r))[:, None]
    flux = 0.5 * np.real(np.sum(np.cross(e, np.conj(h)) * n, axis=1))
    return float(np.sum(area * flux))


def test_rfg_2_the_dipole_wave_built_from_the_dipoles_own_fields_agrees() -> None:
    # Independently of the code's formula. The standing wave regular at the centre is the dipole's field with
    # exp(ikr) replaced by exp(ikr) - exp(-ikr): its converging half carries in what the outgoing half carries out.
    f = 1e9
    k = 2 * math.pi * f / pt.C0
    rng = np.random.default_rng(11)
    dipoles: list[pt.CArr] = [  # C m: along z, turning in the x-y plane, and elliptical in no particular plane
        np.array([0, 0, 1e-12], dtype=np.complex128),
        np.array([1e-12, 1e-12j, 0], dtype=np.complex128),
        np.array([0.3e-12 + 0.1e-12j, -0.2e-12j, 1e-12], dtype=np.complex128),
    ]
    for p in dipoles:
        # near the centre, at kr = 1e-3 in eight directions: i k^3 p / (3 pi eps0)
        n = rng.normal(size=(8, 3))
        r, _, radiation, near, _ = _dipole_shapes(p, 1e-3 / k * n / np.linalg.norm(n, axis=1)[:, None])
        kr = k * r
        e = (
            k**2 * radiation * (2j * np.sin(kr) / r)[:, None]
            + near * (2j * (np.sin(kr) / r**3 - k * np.cos(kr) / r**2))[:, None]
        ) / (4 * math.pi * pt.EPS0)
        centre = 1j * k**3 * p / (3 * math.pi * pt.EPS0)
        assert float(np.abs(e - centre).max()) < 1e-6 * float(np.linalg.norm(centre))
        # the power the outgoing half carries, the same through every sphere: w^4 |p|^2 / (12 pi eps0 c^3)
        power = (2 * math.pi * f) ** 4 * float(np.vdot(p, p).real) / (12 * math.pi * pt.EPS0 * pt.C0**3)
        for kr_sphere in (0.5, 3.0, 300.0):
            assert _outgoing_power(p, k, kr_sphere / k) == pytest.approx(power, rel=1e-12)
        # their ratio is 4 Z0 k^2 / (3 pi), and the code's field for that power is the field at the centre
        squared = np.sum(np.abs(e) ** 2, axis=1) / _outgoing_power(p, k, 3.0 / k)
        assert squared == pytest.approx(np.full(8, 4 * pt.Z0 * k**2 / (3 * math.pi)), rel=1e-6)
        assert pt.dipole_wave_peak_field(power, f) == pytest.approx(float(np.linalg.norm(centre)), rel=1e-6)


def test_rfg_2_control_a_half_wave_spot_falls_short_by_two_pi_over_three() -> None:
    # The estimate that spreads the power over a spot (lambda / 2)^2, E^2 = 2 Z0 P / area, is not the ceiling: the
    # dipole wave's field squared is 2 pi / 3 times it, at every frequency.
    for f in (0.4e9, 2.45e9, 5.8e9):
        spot = 2 * pt.Z0 * 1.0 / (pt.C0 / f / 2) ** 2
        assert pt.dipole_wave_peak_field(1.0, f) ** 2 / spot == pytest.approx(2 * math.pi / 3, rel=1e-12)


@settings(max_examples=200, deadline=None)
@given(st.floats(1.0, 1e7), st.floats(1e-5, 1e-1), st.floats(1e6, 1e10), st.floats(1e-3, 1.0), masses)
def test_rfg_3_electrodes_through_a_resonator_gain_over_the_ceiling(
    resistance: float, distance: float, f: float, share: float, mass: float
) -> None:
    qm = pt.QM_E / mass
    per_watt = pt.electrode_depth_per_watt_ev(resistance, distance, f, share, qm)
    gain = pt.electrode_gain(resistance, distance, f, share)
    assert per_watt / pt.free_space_depth_ceiling_ev(1.0, qm) == pytest.approx(gain, rel=1e-12)
    # the resonator's volts, as RFG-1's drive of four rods at r0 = distance: share of that well, per watt
    power = 0.37
    volts = pt.resonator_volts(power, resistance)
    assert volts**2 / (2 * resistance) == pytest.approx(power, rel=1e-12)
    four_rods = pt.quadrupole_depth_ev(pt.mathieu_q(volts / distance**2, f, qm), f, distance, qm)
    assert per_watt * power == pytest.approx(share * four_rods, rel=1e-12)


def test_rfg_3_numbers_at_1_ghz_half_a_millimetre_and_10_kilohm() -> None:
    f, h, r = 1e9, 0.5e-3, 1e4
    assert pt.electrode_gain(r, h, f) == pytest.approx(1.139e6, rel=1e-3)
    assert pt.electrode_gain(r, h, f, 0.25) == pytest.approx(2.848e5, rel=1e-3)
    assert pt.electrode_depth_per_watt_ev(r, h, f) == pytest.approx(89.10, rel=1e-3)  # eV per watt
    assert pt.electrode_depth_per_watt_ev(r, h, f, 0.25) == pytest.approx(22.28, rel=1e-3)
    assert pt.resonator_volts(3.2e-3, 0.25e6) == pytest.approx(40.0)
    assert pt.resonator_volts(1.0, 1e6) == pytest.approx(1414.2, rel=1e-4)


@settings(max_examples=200, deadline=None)
@given(st.sampled_from(["four rods", "facing", "rails"]), st.floats(1e-4, 1e-2), st.floats(0.1, 20.0), frequencies)
def test_rfg_3_the_depth_at_volts_v_is_share_v_q_over_8_gamma(kind: str, gap: float, ratio: float, f: float) -> None:
    # At volts V a geometry's well is share V q / (8 gamma), q the Mathieu q of its null's own gradient and gamma
    # that gradient over V / h^2; for four rods (share and gamma 1) RFG-1's V0 = 8 D / q.
    volts = 37.0
    if kind == "four rods":
        h, share, gradient = gap, 1.0, volts / gap**2
    else:
        board = pt.FacingBoards(gap) if kind == "facing" else pt.RailPair(gap, ratio * gap)
        h, share, gradient = board.distance, board.depth_share, board.gradient(volts)
    gamma = gradient * h**2 / volts
    q = pt.mathieu_q(gradient, f)
    depth = share * pt.ponderomotive_ev(volts / h, f)
    assert depth == pytest.approx(share * volts * q / (8 * gamma), rel=1e-12)
    if kind == "four rods":
        assert pt.volts_for_depth(depth, q) == pytest.approx(volts, rel=1e-12)


@settings(max_examples=200, deadline=None)
@given(st.floats(1e-3, 1e3), st.floats(0.05, 0.9), st.floats(1e2, 1e7), st.floats(1e-5, 1e-2), shares, masses)
def test_rfg_3_a_depth_at_a_q_takes_volts_power_and_a_frequency(
    depth: float, q: float, resistance: float, distance: float, share_gamma: tuple[float, float], mass: float
) -> None:
    # What a depth takes at a q, and what it takes at the frequency that gives: the same operating point. For four
    # rods it is RFG-1's drive and frequency times distance.
    share, gamma = share_gamma
    qm = pt.QM_E / mass
    at_q = pt.electrode_drive_at_q(depth, q, resistance, distance, share, gamma, qm)
    at_f = pt.electrode_drive(depth, resistance, distance, at_q.f, share, gamma, qm)
    assert at_f == pytest.approx(at_q, rel=1e-9)
    assert pt.resonator_volts(at_q.power, resistance) == pytest.approx(at_q.volts, rel=1e-12)
    per_watt = pt.electrode_depth_per_watt_ev(resistance, distance, at_q.f, share, qm)
    assert per_watt * at_q.power == pytest.approx(depth, rel=1e-9)
    assert pt.mathieu_q(gamma * at_q.volts / distance**2, at_q.f, qm) == pytest.approx(q, rel=1e-12)
    if (share, gamma) == (1.0, 1.0):
        assert at_q.volts == pytest.approx(pt.volts_for_depth(depth, q), rel=1e-12)
        assert at_q.f * distance == pytest.approx(pt.frequency_times_distance(at_q.volts, q, qm), rel=1e-12)


def test_rfg_3_what_1_ev_takes_and_where_the_slope_stops() -> None:
    # 0.5 mm from the null through 10 kOhm. At 1 GHz 1 eV puts four rods at q = 0.53 and two facing boards at
    # 0.84, under Mathieu's 0.908, and one board at its best at 1.73, past it: there it takes a higher frequency.
    # At q = 0.5, where the pseudo-potential holds, 1 eV takes 16 V and 13 mW on four rods at 1.07 GHz, 50 V and
    # 126 mW on facing boards at 1.68 GHz, 542 V and 14.7 W on one board at 3.46 GHz. The depth per watt holds up
    # to q = 0.908: 2.9 eV at 32 mW on four rods, 1.2 eV at 53 mW on facing boards, 0.28 eV at 0.34 W on one board.
    h, r, f = 0.5e-3, 1e4, 1e9
    boards = pt.FacingBoards(2 * h)
    s = (1 + math.sqrt(5) + math.sqrt((1 + math.sqrt(5)) ** 2 - 4)) / 2  # sqrt(1 + 2 rail / gap) at the best
    rails = pt.RailPair(2 * h / s, (s**2 - 1) / s * h)
    assert rails.null_height == pytest.approx(h, rel=1e-12)
    assert rails.gradient_share == pytest.approx(0.3093, abs=1e-4)
    assert boards.gradient_share == pytest.approx(math.pi / 4, rel=1e-15)
    geometries = {"four rods": (1.0, 1.0), "facing": (boards.depth_share, boards.gradient_share)}
    geometries["rails"] = (rails.depth_share, rails.gradient_share)
    at_1ghz = {name: pt.electrode_drive(1.0, r, h, f, *sg) for name, sg in geometries.items()}
    assert [at_1ghz[n].q for n in geometries] == pytest.approx([0.534, 0.839, 1.728], abs=1e-3)
    assert [at_1ghz[n].volts for n in geometries] == pytest.approx([14.98, 29.96, 156.7], rel=1e-3)
    assert at_1ghz["four rods"].power == pytest.approx(11.22e-3, rel=1e-3)
    assert at_1ghz["four rods"].q < 0.908
    assert at_1ghz["facing"].q < 0.908 < at_1ghz["rails"].q
    at_half = {name: pt.electrode_drive_at_q(1.0, 0.5, r, h, *sg) for name, sg in geometries.items()}
    assert [at_half[n].volts for n in geometries] == pytest.approx([16.0, 50.27, 541.7], rel=1e-3)
    assert [at_half[n].power for n in geometries] == pytest.approx([12.8e-3, 126.3e-3, 14.67], rel=1e-3)
    assert [at_half[n].f for n in geometries] == pytest.approx([1.068e9, 1.678e9, 3.456e9], rel=1e-3)
    tops = {"four rods": (2.892, 32.45e-3), "facing": (1.172, 52.61e-3), "rails": (0.2761, 0.3392)}  # eV, W
    for name, (share, gamma) in geometries.items():
        edge = 0.908 / at_1ghz[name].q * at_1ghz[name].volts  # at one frequency q grows as the volts
        assert pt.mathieu_q(gamma * edge / h**2, f) == pytest.approx(0.908, rel=1e-12)
        power = edge**2 / (2 * r)
        assert (pt.electrode_depth_per_watt_ev(r, h, f, share) * power, power) == pytest.approx(tops[name], rel=1e-3)


@pytest.mark.parametrize("ratio", [0.5, 3.68, 8.0, None], ids=["rails 0.5", "rails 3.68", "rails 8", "facing"])
def test_rfg_3_a_boards_share_turns_its_resonators_power_into_its_depth(ratio: float | None) -> None:
    # The depth of a board's well, read from a flood over its field rather than from its own barrier, is the share
    # it states of RFG-3's four-rod well at the resonator's volts, at the frequency that puts its null's Mathieu q
    # at 0.5, where the pseudo-potential holds.
    board: pt.RailPair | pt.FacingBoards
    if ratio is None:
        board = pt.FacingBoards(1e-3)
        h = board.distance
        xs, ys = np.linspace(-6, 6, 241) * h, np.linspace(-1, 1, 81)[1:-1] * h
        start = (120, 39)
    else:
        board = pt.RailPair(1e-3, ratio * 1e-3)
        h = board.distance
        xs, ys = np.linspace(-12, 12, 481) * h, np.linspace(1e-3, 12.0, 240) * h
        start = (240, int(np.argmin(np.abs(ys - h))))
    barrier = _lowest_way_out(np.linalg.norm(board.field(xs[:, None], ys[None, :], 1.0), axis=-1), start)
    for power, resistance in ((1e-3, 1e4), (0.37, 2e3), (2.0, 1e5)):
        volts = pt.resonator_volts(power, resistance)
        f = math.sqrt(2 * abs(pt.QM_E) * board.gradient(volts) / 0.5) / (2 * math.pi)
        assert pt.mathieu_q(board.gradient(volts), f) == pytest.approx(0.5, rel=1e-12)
        per_watt = pt.electrode_depth_per_watt_ev(resistance, h, f, board.depth_share)
        assert pt.ponderomotive_ev(barrier * volts, f) == pytest.approx(per_watt * power, rel=0.02)
        drive = pt.electrode_drive(per_watt * power, resistance, h, f, board.depth_share, board.gradient_share)
        assert (drive.power, drive.volts, drive.q) == pytest.approx((power, volts, 0.5), rel=1e-9)


@settings(max_examples=200, deadline=None)
@given(st.floats(0.01, 0.9), st.floats(1e3, 1e11), masses)
def test_rfg_4_at_the_space_charge_limit_the_plasma_frequency_is_root_two_secular(
    q: float, f: float, mass: float
) -> None:
    qm = pt.QM_E / mass
    plasma = pt.plasma_frequency(pt.space_charge_limit(q, f, qm), qm)
    assert plasma == pytest.approx(math.sqrt(2) * pt.secular_frequency(q, f), rel=1e-12)
    assert plasma == pytest.approx(q * f / 2, rel=1e-12)


def test_rfg_4_numbers_for_an_electron() -> None:
    # the first run's guide: q = 0.6 at 1 GHz holds 1.12e15 electrons per cubic metre, at a plasma frequency of
    # 300 MHz; and 8.98 sqrt(n) Hz in general
    assert pt.plasma_frequency(pt.space_charge_limit(0.6, 1e9)) == pytest.approx(300e6, rel=1e-12)
    assert pt.plasma_frequency(1e15) == pytest.approx(8.98 * math.sqrt(1e15), rel=1e-3)


# ---- board electrodes (electrodes.py) -------------------------------------------------------------------------
gaps = st.floats(1e-5, 1e-2)  # m
rail_ratios = st.floats(0.1, 20.0)  # rail / gap
drives = st.floats(0.1, 1e3)  # V


def _complex_field(board: pt.RailPair | pt.FacingBoards, z: complex, volts: float) -> complex:
    """E_x - i E_y at z = x + iy: the analytic function a field across long electrodes is."""
    ex, ey = board.field(z.real, z.imag, volts)
    return complex(ex, -ey)


def _circle_mean(potential: Callable[[pt.Arr, pt.Arr], pt.Arr], x0: float, y0: float, radius: float) -> float:
    """A potential's mean over a circle, by the trapezoid rule on 64 points. For a harmonic function it is the value
    at the centre, and the rule converges geometrically when the circle stays clear of the boundary."""
    t = 2 * np.pi * np.arange(64) / 64
    return float(np.mean(potential(x0 + radius * np.cos(t), y0 + radius * np.sin(t))))


def _lowest_way_out(magnitude: pt.Arr, start: tuple[int, int]) -> float:
    """Over the paths on a grid from `start` to the grid's edge, the least of the largest value along a path: the
    barrier a particle starting there crosses to leave. A flood in order of value (a minimax Dijkstra)."""
    best = np.full(magnitude.shape, np.inf)
    best[start] = magnitude[start]
    queue = [(float(magnitude[start]), start)]
    last = (magnitude.shape[0] - 1, magnitude.shape[1] - 1)
    while queue:
        level, (i, j) = heapq.heappop(queue)
        if level > best[i, j]:
            continue
        if i in (0, last[0]) or j in (0, last[1]):
            return level
        for k, m in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)):
            higher = max(level, float(magnitude[k, m]))
            if higher < best[k, m]:
                best[k, m] = higher
                heapq.heappush(queue, (higher, (k, m)))
    raise AssertionError("the flood did not reach the edge")


def _golden_max(fn: Callable[[float], float], lo: float, hi: float) -> tuple[float, float]:
    """Where a function with one maximum in [lo, hi] is largest, and its value there, by golden-section search
    run until the bracket stops narrowing."""
    g = (math.sqrt(5) - 1) / 2
    a, b = lo, hi
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = fn(c), fn(d)
    for _ in range(200):
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = fn(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = fn(d)
    middle = (a + b) / 2
    return middle, fn(middle)


@settings(max_examples=200, deadline=None)
@given(gaps, rail_ratios, drives)
def test_ele_1_the_rails_null_stands_at_the_geometric_mean_of_their_edges(
    gap: float, ratio: float, volts: float
) -> None:
    rails = pt.RailPair(gap, ratio * gap)
    h, (c1, c2) = rails.null_height, rails.edges
    assert h == pytest.approx(math.sqrt(c1 * c2), rel=1e-12)
    assert float(np.hypot(*rails.field(0.0, h, volts))) < 1e-9 * volts / h
    # the gradient there is the derivative of the complex field, the same whichever way it is taken
    d, null = 1e-4 * h, 1j * h
    along_x = (_complex_field(rails, null + d, volts) - _complex_field(rails, null - d, volts)) / (2 * d)
    along_y = (_complex_field(rails, null + 1j * d, volts) - _complex_field(rails, null - 1j * d, volts)) / (2j * d)
    assert abs(along_x) == pytest.approx(rails.gradient(volts), rel=1e-6)
    assert abs(along_y - along_x) < 1e-6 * abs(along_x)
    # as a share of four rods' V / h^2: 4 b h / (pi (a + b)^2)
    assert rails.gradient_share == pytest.approx(4 * rails.rail * h / (math.pi * (gap + rails.rail) ** 2), rel=1e-12)


@settings(max_examples=200, deadline=None)
@given(gaps, rail_ratios, drives)
def test_ele_1_control_a_tenth_off_the_null_the_field_is_not_null(gap: float, ratio: float, volts: float) -> None:
    # The test above would catch a misplaced null: a tenth of its height above or below it the field is over a
    # thousandth of V / h, a million times the bound at the null.
    rails = pt.RailPair(gap, ratio * gap)
    h = rails.null_height
    for off in (0.9, 1.1):
        assert float(np.hypot(*rails.field(0.0, off * h, volts))) > 1e-3 * volts / h


@settings(max_examples=200, deadline=None)
@given(rail_ratios, st.floats(-4.0, 4.0), st.floats(0.05, 3.0), drives)
def test_ele_1_the_rails_potential_is_harmonic_holds_their_volts_and_gives_the_field(
    ratio: float, x: float, y: float, volts: float
) -> None:
    rails = pt.RailPair(1e-3, ratio * 1e-3)
    h, (c1, c2) = rails.null_height, rails.edges
    x0, y0 = x * h, y * h

    def phi(xs: pt.Arr | float, ys: pt.Arr | float) -> pt.Arr:
        return rails.potential(xs, ys, volts)

    assert abs(_circle_mean(phi, x0, y0, y0 / 2) - float(phi(x0, y0))) < 1e-12 * volts
    # on the board: the rails' volts on the rails, and ground on the centre strip and beyond the rails
    for at, held in ((0.0, 0.0), ((c1 + c2) / 2, volts), (-(c1 + c2) / 2, volts), (c2 + h, 0.0), (-c2 - h, 0.0)):
        assert float(phi(at, 1e-9 * h)) == pytest.approx(held, abs=1e-6 * volts)
    # the field is minus the potential's gradient, by central differences
    d = 1e-6 * h
    ex, ey = rails.field(x0, y0, volts)
    assert abs(ex + (float(phi(x0 + d, y0)) - float(phi(x0 - d, y0))) / (2 * d)) < 1e-6 * volts / h
    assert abs(ey + (float(phi(x0, y0 + d)) - float(phi(x0, y0 - d))) / (2 * d)) < 1e-6 * volts / h


@settings(max_examples=200, deadline=None)
@given(gaps, st.floats(-3.0, 3.0), st.floats(-0.9, 0.9), st.floats(0.1, 1.0), drives)
def test_ele_2_the_facing_boards_potential_is_harmonic(
    gap: float, x: float, y: float, size: float, volts: float
) -> None:
    boards = pt.FacingBoards(gap)
    h = boards.distance

    def phi(xs: pt.Arr | float, ys: pt.Arr | float) -> pt.Arr:
        return boards.potential(xs, ys, volts)

    # circles reaching at most half way to the nearer board
    x0, y0, radius = x * h, y * h, size * (1 - abs(y)) * h / 2
    assert abs(_circle_mean(phi, x0, y0, radius) - float(phi(x0, y0))) < 1e-12 * volts


@settings(max_examples=200, deadline=None)
@given(st.floats(0.2, 2.0), st.floats(0.2, 0.6))
def test_ele_2_control_a_nearby_potential_is_not_harmonic(x: float, y: float) -> None:
    # Without the tangent, -(V / pi) atan((pi y / 2h) tanh(pi x / 2h)) has the same null and nearly the same values
    # near the centre, and its circle means miss its centre values by more than 1e-6 of V (4e-3 at least here).
    h, volts = 1e-3, 1.0
    k = math.pi / (2 * h)

    def nearby(xs: pt.Arr | float, ys: pt.Arr | float) -> pt.Arr:
        out: pt.Arr = -volts / math.pi * np.arctan(k * np.asarray(ys) * np.tanh(k * np.asarray(xs)))
        return out

    x0, y0, radius = x * h, y * h, (1 - y) * h / 2
    assert abs(_circle_mean(nearby, x0, y0, radius) - float(nearby(x0, y0))) > 1e-6 * volts


@settings(max_examples=200, deadline=None)
@given(st.floats(1e-4, 1e-1), st.floats(0.05, 20.0), drives)
def test_ele_2_the_boards_hold_their_volts(gap: float, x: float, volts: float) -> None:
    boards = pt.FacingBoards(gap)
    h = boards.distance
    near = (1 - 1e-9) * h  # a billionth of h from the boards
    for side in (1.0, -1.0):
        assert float(boards.potential(side * x * h, near, volts)) == pytest.approx(-side * volts / 2, abs=1e-6 * volts)
        assert float(boards.potential(side * x * h, -near, volts)) == pytest.approx(side * volts / 2, abs=1e-6 * volts)


def _fourier_midplane(x: float, h: float) -> float:
    """(1 / pi) times the integral over k from 0 to infinity of sin(k x) / sinh(k h): the facing boards' midplane
    field per volt, from the Fourier integral of the boards' volts. Composite Gauss-Legendre in k h from 0 to 40,
    past which the integrand is below 1e-17 of its start."""
    nodes, weights = np.polynomial.legendre.leggauss(20)
    u = (np.arange(40.0)[:, None] + (nodes[None, :] + 1) / 2).ravel()
    return float(np.sum(np.tile(weights / 2, 40) * np.sin(u * x / h) / np.sinh(u))) / (math.pi * h)


@settings(max_examples=100, deadline=None)
@given(st.floats(1e-4, 1e-1), st.floats(-4.0, 4.0), drives)
def test_ele_2_the_midplane_field_is_v_over_2h_tanh_and_the_fourier_integral(
    gap: float, x: float, volts: float
) -> None:
    boards = pt.FacingBoards(gap)
    h = boards.distance
    ex, ey = boards.field(x * h, 0.0, volts)
    assert abs(ex) < 1e-15 * volts / h  # straight across the gap
    assert ey == pytest.approx(volts / (2 * h) * math.tanh(math.pi * x / 2), rel=1e-12, abs=1e-15 * volts / h)
    assert ey == pytest.approx(volts * _fourier_midplane(x * h, h), rel=1e-8, abs=1e-8 * volts / h)


@settings(max_examples=200, deadline=None)
@given(st.floats(1e-4, 1e-1), st.floats(-3.0, 3.0), st.floats(-0.8, 0.8), drives)
def test_ele_2_the_field_is_minus_the_gradient_of_the_potential(gap: float, x: float, y: float, volts: float) -> None:
    boards = pt.FacingBoards(gap)
    h = boards.distance
    x0, y0, d = x * h, y * h, 1e-6 * h
    ex, ey = boards.field(x0, y0, volts)
    gx = (float(boards.potential(x0 + d, y0, volts)) - float(boards.potential(x0 - d, y0, volts))) / (2 * d)
    gy = (float(boards.potential(x0, y0 + d, volts)) - float(boards.potential(x0, y0 - d, volts))) / (2 * d)
    assert abs(ex + gx) < 1e-7 * volts / h
    assert abs(ey + gy) < 1e-7 * volts / h


def test_ele_2_a_null_of_gradient_pi_v_over_4h2_and_a_quarter_of_four_rods() -> None:
    boards, volts = pt.FacingBoards(2e-3), 10.0
    h = boards.distance
    assert boards.gradient(volts) == pytest.approx(math.pi * volts / (4 * h**2), rel=1e-15)
    assert boards.gradient_share == pytest.approx(math.pi / 4, rel=1e-15)
    d = 1e-5 * h  # near the centre |E| is the gradient times the distance, in every direction
    for angle in np.linspace(0, 2 * np.pi, 9):
        e = boards.field(d * math.cos(angle), d * math.sin(angle), volts)
        assert float(np.hypot(*e)) == pytest.approx(boards.gradient(volts) * d, rel=1e-9)
    # on the midplane the field rises toward V / 2h and stays below it; off the midplane it is stronger, and next to
    # the boards stronger than V / 2h
    barrier = boards.barrier_field(volts)
    assert barrier == pytest.approx(volts / (2 * h), rel=1e-15)
    xs = np.linspace(-5 * h, 5 * h, 2001)
    across = [np.linalg.norm(boards.field(xs, y * h, volts), axis=-1) for y in (0.0, 0.3, 0.6, 0.9, 1 - 1e-9)]
    assert float(across[0].max()) < barrier
    assert float(across[0].max()) == pytest.approx(barrier, rel=1e-6)
    for lower, higher in itertools.pairwise(across):
        assert bool(np.all(higher > lower))
    assert float(across[-1].min()) > barrier
    # the lowest way out, by a flood over the slab from the centre: V / 2h, approached far along the midplane
    grid = np.linalg.norm(
        boards.field(np.linspace(-6, 6, 241)[:, None] * h, np.linspace(-1, 1, 81)[None, 1:-1] * h, volts), axis=-1
    )
    assert _lowest_way_out(grid, (120, 39)) == pytest.approx(barrier, rel=1e-6)
    assert boards.depth_share == pytest.approx(0.25, rel=1e-15)


@settings(max_examples=100, deadline=None)
@given(gaps, rail_ratios, drives)
def test_ele_3_the_escape_is_the_saddle_above_the_null(gap: float, ratio: float, volts: float) -> None:
    rails = pt.RailPair(gap, ratio * gap)
    h, saddle = rails.null_height, rails.escape(volts)

    def on_axis(y: float) -> float:
        return float(np.hypot(*rails.field(0.0, y, volts)))

    # an independent search for the strongest field on the axis above the null
    height, field = _golden_max(on_axis, h, 20 * h)
    assert height == pytest.approx(saddle.height, rel=1e-6)
    assert field == pytest.approx(saddle.field, rel=1e-12)
    # across the axis the field grows: a saddle of |E|, not a peak
    for side in (-0.05 * h, 0.05 * h):
        assert float(np.hypot(*rails.field(side, saddle.height, volts))) > saddle.field * (1 + 1e-5)
    w = (rails.gap + rails.rail) / h
    assert rails.depth_share == pytest.approx(4 / math.pi**2 * (w - 2) / ((w + 2) * w**2), rel=1e-12)
    assert rails.depth_share == pytest.approx((saddle.field * h / volts) ** 2, rel=1e-12)


@pytest.mark.parametrize("ratio", [0.25, 1.0, 3.68, 8.0])
def test_ele_3_the_saddle_is_the_lowest_way_out(ratio: float) -> None:
    # A flood from the null over a grid 12 null heights to either side and up, down to a thousandth of h above the
    # board: no path to the board or to the grid's edge crosses a weaker field than the saddle.
    rails = pt.RailPair(1e-3, ratio * 1e-3)
    h = rails.null_height
    ys = np.linspace(1e-3, 12.0, 240) * h
    grid = np.linalg.norm(rails.field(np.linspace(-12, 12, 481)[:, None] * h, ys[None, :], 1.0), axis=-1)
    start = (240, int(np.argmin(np.abs(ys - h))))
    assert _lowest_way_out(grid, start) == pytest.approx(rails.escape(1.0).field, rel=0.01)


def test_ele_3_one_board_holds_at_most_a_hundredth_of_four_rods() -> None:
    ratios = np.linspace(0.25, 8.0, 3101)  # rail / gap
    boards = [pt.RailPair(1e-3, r * 1e-3) for r in ratios]
    shares = np.array([rails.depth_share for rails in boards])
    best = float(shares.max())
    assert best <= 0.0092
    assert best == pytest.approx((5 * math.sqrt(5) - 11) / (2 * math.pi**2), rel=1e-6)
    assert float(ratios[shares.argmax()]) == pytest.approx(3.68, abs=0.01)  # near rail / gap = 4
    assert pt.RailPair(1e-3, 4e-3).depth_share == pytest.approx(0.0091, abs=5e-5)
    heights = np.array([rails.escape(1.0).height / rails.null_height for rails in boards])
    assert bool(np.all(np.diff(heights) > 0))  # the saddle rises with rail / gap, from 1.74 to 2.32 null heights
    assert (float(heights[0]), float(heights[-1])) == pytest.approx((1.744, 2.316), abs=1e-3)
    assert pt.FacingBoards(2e-3).depth_share >= 27 * best
