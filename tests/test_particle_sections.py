"""A solved cross-section swept along a path (array16/particles/sections.py), held to what it must reproduce: the board
electrodes' closed-form fields on a straight path to the interpolation error (ELE-1, ELE-2) and turned by the framing
angle; its null on the path round R1's rail ring, the flat ring and T915-S's trefoil framed k = 0, -1 and -3, within
0.05h (DR-3); Mathieu's exact secular frequency for a tracked electron to 0.1 % at 120 steps an RF period (DR-0 of
docs/KNOTTED_GUIDE.md section 8.2); the travelling wave's phase velocity and its magnetic field T x E / v; nothing on
what is off the grid; and the same field whatever the search remembers between calls."""

from __future__ import annotations

import functools
import math
from collections.abc import Callable

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from array16 import particles as pt
from array16.particles import ClosedCurve, FacingBoards, RailPair, SectionField, SectionGrid

H915 = 0.915e-3  # T915-S: null to copper (m), the gap 1.83 mm
F915 = 915e6
L915 = 0.318  # T915-S's trefoil, and the flat ring of its length (m)
R1_RAILS = (0.40e-3, 1.46e-3)  # R1: a 0.40 mm grounded strip between rails 1.46 mm wide, the null 0.576 mm up
R1_RADIUS = 0.032  # R1's rail ring, one guided wave round at 915 MHz (m)
OMEGA = 2 * math.pi * F915


# ---- sections and paths ---------------------------------------------------------------------------------------------
def _even(lo: float, hi: float, step: float) -> pt.Arr:
    out: pt.Arr = np.linspace(lo, hi, round((hi - lo) / step) + 1)
    return out


def _graded(centre: float, end: float, fine: float, stretch: float = 1.0) -> pt.Arr:
    """Nodes from centre to end, `fine` apart at the centre and stretching away from it: sinh-spaced."""
    count = max(2, math.ceil(abs(end - centre) * stretch / (math.sinh(stretch) * fine)))
    out: pt.Arr = centre + (end - centre) * np.sinh(stretch * np.linspace(0.0, 1.0, count + 1)) / math.sinh(stretch)
    return out


@functools.cache
def _facing(div: int, across: float = 2.0, toward: float = 0.95) -> SectionGrid:
    """ELE-2's potential per volt round the facing boards' null, h / div apart, `across` h either side of the split
    and `toward` h up and down toward the boards."""
    boards = FacingBoards(2 * H915)
    u = _even(-across * H915, across * H915, H915 / div)
    v = _even(-toward * H915, toward * H915, H915 / div)
    return SectionGrid.sampled(lambda a, b: boards.potential(a, b, 1.0), u, v, (0.0, 0.0))


@functools.cache
def _rails(div: int, graded: bool = False) -> SectionGrid:
    """ELE-1's potential per volt over R1's rails, from 0.2h to 3h above the board and 3h either side, h / div apart
    (at the null, stretching away from it when graded)."""
    rails = RailPair(*R1_RAILS)
    h = rails.null_height
    if graded:
        side = _graded(0.0, 3 * h, h / div)
        u = np.concatenate([-side[:0:-1], side])
        v = np.concatenate([_graded(h, 0.2 * h, h / div)[:0:-1], _graded(h, 3 * h, h / div)])
    else:
        u, v = _even(-3 * h, 3 * h, h / div), _even(0.2 * h, 3 * h, h / div)
    return SectionGrid.sampled(lambda a, b: rails.potential(a, b, 1.0), u, v, (0.0, h))


def _facing_curvature(z: pt.CArr) -> pt.Arr:
    """|f''| for f = E_x - i E_y = -i (1/2h) tanh(pi z / 2h), the facing boards' field per volt: (k^2 / h)
    |sech^2(kz) tanh(kz)|, k = pi / 2h."""
    k = math.pi / (2 * H915)
    out: pt.Arr = k * k / H915 * np.abs(np.tanh(k * z) / np.cosh(k * z) ** 2)
    return out


def _rails_curvature(z: pt.CArr) -> pt.Arr:
    """|f''| for f = (i / pi) [1/(z - c2) - 1/(z - c1) + 1/(z + c1) - 1/(z + c2)], the rails' field per volt."""
    c1, c2 = RailPair(*R1_RAILS).edges
    total = sum(sign / (z - c) ** 3 for c, sign in ((c2, 1), (c1, -1), (-c1, 1), (-c2, -1)))
    out: pt.Arr = 2 / math.pi * np.abs(total)
    return out


def _trefoil(n: int) -> ClosedCurve:
    """T915-S's r/R = 0.45 trefoil at 318 mm."""
    knot = ClosedCurve.torus_knot(1.0, 0.45, n)
    return ClosedCurve(knot.points * (L915 / knot.length), L915)


def _closing(curve: ClosedCurve, k: int, theta0: float = 0.3) -> pt.Arr:
    """The framing angle that closes after a lap with k half turns, the knotted study's rule (k pi - holonomy) / L a
    metre on the twist-free axes, turned a quarter so that u runs across the boards' plane."""
    s = np.arange(len(curve)) * curve.length / len(curve)
    out: pt.Arr = theta0 + math.pi / 2 + (k * math.pi - curve.holonomy) / curve.length * s
    return out


def _e(field: SectionField, x: pt.Arr, t: float = 0.0) -> pt.Arr:
    e, _ = field(x, t)
    assert e is not None
    return e


# ---- a straight sweep against the closed forms ----------------------------------------------------------------------
@pytest.mark.parametrize("case", ["facing boards, even grid", "rails, graded grid"])
def test_a_straight_sweep_is_the_closed_form_field_to_the_interpolation_error(case: str) -> None:
    # ELE-2's and ELE-1's potentials sampled on grids and swept along a straight line: within 0.6h of the null, in
    # every period of the line, the field is the closed form's to what the grid's differences and interpolation
    # allow, (5/12) step^2 max|f''| a component (f = E_x - i E_y), and has nothing along the path. The bound is a few
    # thousandths of the largest field, so a frame turned by 2 mrad, or a null misplaced by 0.002h, would break it.
    if case.startswith("facing"):
        boards = FacingBoards(2 * H915)
        grid, h, curvature = _facing(80), H915, _facing_curvature
        closed: Callable[[pt.Arr, pt.Arr], pt.Arr] = lambda a, b: boards.field(a, b, 1.0)  # noqa: E731
    else:
        rails = RailPair(*R1_RAILS)
        grid, h, curvature = _rails(80, graded=True), rails.null_height, _rails_curvature
        closed = lambda a, b: rails.field(a, b, 1.0)  # noqa: E731
    u0, v0 = grid.null
    field = SectionField.straight(grid, 0.02, pt.InPhase(1.0, 0.0))
    rng = np.random.default_rng(3)
    a, b = rng.uniform(-0.6 * h, 0.6 * h, (2, 4000))
    x = np.stack([a, b, rng.uniform(-0.05, 0.07, 4000)], 1)
    e, magnetic = field(x, 0.0)
    assert e is not None
    assert magnetic is None
    assert not field.off_grid.any()
    assert np.all(e[:, 2] == 0.0)
    # the bound: the widest cell within the stencils' reach of the region, and |f''| over that reach
    margin = 0.6 * h + 2 * max(float(np.diff(grid.u).max()), float(np.diff(grid.v).max()))
    near_u = grid.u[np.abs(grid.u - u0) <= margin]
    near_v = grid.v[np.abs(grid.v - v0) <= margin]
    step = max(float(np.diff(near_u).max()), float(np.diff(near_v).max()))
    reach = 0.6 * h + 2 * step
    zz = np.linspace(u0 - reach, u0 + reach, 301)[:, None] + 1j * np.linspace(v0 - reach, v0 + reach, 301)[None, :]
    bound = 5 / 12 * step**2 * float(curvature(zz).max())
    error = float(np.abs(e[:, :2] - closed(u0 + a, v0 + b)).max())
    assert error <= bound
    assert bound < 4e-3 * float(np.abs(e).max())


@settings(max_examples=40, deadline=None)
@given(
    st.floats(-math.pi, math.pi),
    st.one_of(st.integers(-3, 3).map(float), st.floats(-3.0, 3.0)),
    st.integers(-2, 2),
    st.one_of(st.floats(0.0, 1.0, exclude_max=True), st.floats(399 / 400, 1.0, exclude_max=True)),
    st.floats(-0.55, 0.55),
    st.floats(-0.55, 0.55),
)
def test_the_frame_turns_the_section_by_the_framing_angle(
    theta0: float, turns: float, lap: int, along: float, a: float, b: float
) -> None:
    # A straight sweep twisted by theta(s) = theta0 + turns pi s / L: at z, a fraction `along` of the way through any
    # period of the line, the field is the section's turned by theta0 + turns pi along about the path, also in the
    # last sample's segment (half the draws), where the frame turns on toward the frame it comes back with, and seam
    # reads that turn, turns pi. A whole number of half turns closes the frame up to a turn the facing boards' section
    # (odd across the split and toward the boards) takes into itself.
    length, samples = 0.02, 400
    grid = _facing(80)
    s = np.arange(samples) * length / samples
    field = SectionField.straight(
        grid, length, pt.InPhase(1.0, 0.0), samples=samples, theta=theta0 + turns * math.pi * s / length
    )
    assert abs(math.remainder(field.seam - turns * math.pi, 2 * math.pi)) < 1e-9
    z = (lap + along) * length
    th = theta0 + turns * math.pi * along
    first = np.array([math.cos(th), math.sin(th), 0.0])
    second = np.array([-math.sin(th), math.cos(th), 0.0])
    x = (a * H915 * first + b * H915 * second + [0.0, 0.0, z])[None, :]
    across, toward = FacingBoards(2 * H915).field(a * H915, b * H915, 1.0)
    expected = across * first + toward * second
    # V/m per volt: inside the h/80 grid's bound of 0.85 (the straight sweep's) on fields up to V/2h = 546
    assert float(np.abs(_e(field, x)[0] - expected).max()) < 0.6


# ---- the null on the path -------------------------------------------------------------------------------------------
def _null_offsets(field: SectionField, centre: pt.Arr, first: pt.Arr, second: pt.Arr, h: float) -> pt.Arr:
    """For each centre on the path, the null of the swept field in the plane through it spanned by the frame there,
    by Newton's method on (E . first, E . second): its distance from the centre (m)."""
    ab = np.tile([0.15 * h, -0.1 * h], (len(centre), 1))
    step = 1e-4 * h

    def normal_field(ab: pt.Arr) -> pt.Arr:
        e = _e(field, centre + ab[:, :1] * first + ab[:, 1:] * second)
        out: pt.Arr = np.stack([(e * first).sum(1), (e * second).sum(1)], 1)
        return out

    for _ in range(20):
        here = normal_field(ab)
        jac = np.stack(
            [
                (normal_field(ab + np.array([step, 0.0])) - here) / step,
                (normal_field(ab + np.array([0.0, step])) - here) / step,
            ],
            2,
        )
        move = np.linalg.solve(jac, -here[:, :, None])[:, :, 0]
        ab = ab + move
        if float(np.abs(move).max()) < 1e-12 * h:
            break
    out: pt.Arr = np.hypot(ab[:, 0], ab[:, 1])
    return out


def _ring_stations(radius: float, theta: float, count: int) -> tuple[pt.Arr, pt.Arr, pt.Arr]:
    """Points on a flat ring and its frame turned by theta from the twist-free axes (e1 the outward radius, e2 = T x
    e1 down the ring's axis), at angles that fall on samples and between them."""
    angle = np.linspace(0.0, 2 * math.pi, count, endpoint=False) + 0.37 * 2 * math.pi / count
    zero = np.zeros(count)
    radial = np.stack([np.cos(angle), np.sin(angle), zero], 1)
    down = np.tile([0.0, 0.0, -1.0], (count, 1))
    return (
        radius * radial,
        math.cos(theta) * radial + math.sin(theta) * down,
        math.cos(theta) * down - math.sin(theta) * radial,
    )


@pytest.mark.parametrize(
    "case", ["R1's rail ring", "flat facing-board ring", "trefoil k = 0", "trefoil k = -1", "trefoil k = -3"]
)
def test_swept_round_a_ring_or_a_knot_the_null_stays_on_the_path(case: str) -> None:
    # DR-3's null within 0.05h: wherever the frame closes, the swept field vanishes on the path, between samples and
    # across the seam. R1's rails (one board, the null 0.576 mm up) bend round 32 mm with u inward and v up; the flat
    # ring of 318 mm and T915-S's trefoil carry the facing boards, the trefoil framed by the knotted study's rule with
    # k = 0 (seam 0) and k = -1 and -3 (seam pi, a turn the boards' section takes into itself).
    if case.startswith("R1"):
        grid, h = _rails(40), RailPair(*R1_RAILS).null_height
        field = SectionField.along(grid, ClosedCurve.circle(R1_RADIUS, 2000), pt.InPhase(1.0, 0.0), theta=math.pi)
        centre, first, second = _ring_stations(R1_RADIUS, math.pi, 97)
    elif case.startswith("flat"):
        grid, h, radius = _facing(40), H915, L915 / (2 * math.pi)
        field = SectionField.along(grid, ClosedCurve.circle(radius, 2000), pt.InPhase(1.0, 0.0))
        centre, first, second = _ring_stations(radius, 0.0, 97)
    else:
        k = int(case.split("= ")[1])
        knot = _trefoil(3000)
        grid, h = _facing(40), H915
        field = SectionField.along(grid, knot, pt.InPhase(1.0, 0.0), theta=_closing(knot, k))
        assert abs(math.remainder(field.seam - k * math.pi, 2 * math.pi)) < 1e-9
        pick = np.unique(np.concatenate([np.arange(0, 3000, 31), [2998, 2999, 0, 1]]))
        centre, first, second = knot.points[pick], field.first[pick], field.second[pick]
    if case.endswith("ring"):
        assert abs(field.seam) < 1e-12  # a flat ring's twist-free axes close
    offsets = _null_offsets(field, centre, first, second, h)
    assert float(offsets.max()) < 0.05 * h
    if case.startswith("R1"):
        # the rails' null is where the grid's own differences put it, 3.2e-4 h above the closed form's, and the swept
        # null follows it: misplacing the section's null by 0.1h moves the swept one off the path by as much
        grid_null = grid.find_null()
        assert offsets == pytest.approx(np.full(len(offsets), abs(grid_null[1] - h)), abs=1e-3 * h)
        moved = SectionGrid(grid.u, grid.v, grid.phi, (0.1 * h, h))
        wrong = SectionField.along(moved, ClosedCurve.circle(R1_RADIUS, 2000), pt.InPhase(1.0, 0.0), theta=math.pi)
        assert float(_null_offsets(wrong, centre, first, second, h).min()) > 0.05 * h


def test_the_field_crosses_the_seam_of_a_half_turn_framing_without_a_jump() -> None:
    # T915-S's trefoil framed k = -1 comes back after a lap turned by pi: the segment from the last sample to the first
    # turns as the one before it does, so the field is continuous across the seam as it is between any two samples,
    # changing in proportion to the step taken along the path (a jump would stay as the step shrinks).
    knot = _trefoil(3000)
    field = SectionField.along(_facing(40), knot, pt.InPhase(1.0, 0.0), theta=_closing(knot, -1))
    for at in (0, 1500):  # the seam, and a sample half a lap away
        base = knot.points[at] + 0.3 * H915 * field.first[at] - 0.2 * H915 * field.second[at]
        change = []
        for half in (1e-7, 1e-9):
            e = _e(field, base + np.array([[-half], [half]]) * knot.tangents[at])
            change.append(float(np.linalg.norm(e[1] - e[0]) / np.linalg.norm(e[0])))
        assert change[1] == pytest.approx(change[0] / 100, rel=0.05)
        assert change[0] < 3e-5


def test_the_seam_is_an_ordinary_segment_of_the_path() -> None:
    # Started half a lap round, with the normals the half turn reversed flipped back, the same k = -1 framing puts the
    # old seam's segment in the middle of the path and the half turn at the new start. Through the old seam's segment,
    # and the new one's, both fields are the same to rounding: the segment that crosses a seam turns as any other.
    knot = _trefoil(3000)
    grid, drive = _facing(40), pt.InPhase(1.0, 0.0)
    field = SectionField.along(grid, knot, drive, theta=_closing(knot, -1))
    order = np.roll(np.arange(3000), -1500)
    flip = np.where(order < 1500, -1.0, 1.0)[:, None]  # old samples 0 to 1499 now follow old 2999
    rolled = SectionField(
        grid,
        knot.points[order],
        knot.tangents[order],
        flip * field.first[order],
        flip * field.second[order],
        L915,
        drive,
    )
    assert abs(math.remainder(rolled.seam - math.pi, 2 * math.pi)) < 1e-9
    f = np.linspace(0.05, 0.95, 7)[:, None]
    points = []
    for last, first in ((2999, 0), (1499, 1500)):  # the old seam's segment, and the new one's
        chord = (1 - f) * knot.points[last] + f * knot.points[first]
        for a, b in ((0.3, -0.2), (-0.45, 0.35)):
            points.append(chord + a * H915 * field.first[last] + b * H915 * field.second[last])
    x = np.concatenate(points)
    e = _e(field, x)
    assert np.allclose(e, _e(rolled, x), rtol=0.0, atol=1e-9 * float(np.abs(e).max()))
    assert not field.off_grid.any()
    # the curve's own twist-free axes come back turned by its holonomy, -153 degrees on this trefoil
    assert SectionField.along(grid, knot, drive).seam == pytest.approx(knot.holonomy, abs=1e-9)


def test_a_jump_onto_the_strand_passing_nearest_is_searched_in_full() -> None:
    # Where T915-S's strands pass closest, about 18 mm apart, a particle near one strand that jumps next to the other
    # finds, in the window round its old sample, a nearest sample that is a local minimum of its distance but 18 mm
    # off: past the curve's thickness, so it is searched in full and put on the other strand, as a field that
    # remembers nothing puts it.
    knot = _trefoil(1200)
    field = SectionField.along(_facing(20), knot, pt.InPhase(1.0, 0.0), theta=_closing(knot, 0))
    index = np.arange(1200)
    apart = np.abs(index[:, None] - index[None, :])
    apart = np.minimum(apart, 1200 - apart) * (L915 / 1200)
    gaps = np.linalg.norm(knot.points[:, None, :] - knot.points[None, :, :], axis=2)
    i, j = np.unravel_index(int(np.argmin(np.where(apart > math.pi * knot.bend_radius, gaps, np.inf))), gaps.shape)
    assert 0.015 < float(gaps[i, j]) < 0.021
    before = knot.points[i] + 0.2 * H915 * field.first[i]
    after = knot.points[j] + 0.2 * H915 * field.first[j]
    window = np.linalg.norm(after - knot.points[(i + np.arange(-2, 3)) % 1200], axis=1)
    assert 0 < int(np.argmin(window)) < 4  # what the window round the old sample sees: a minimum inside it
    _e(field, before[None, :])
    e = _e(field, after[None, :])
    assert not field.off_grid.any()
    fresh, _ = field.field_at(field.locate(after[None, :]), 0.0)
    assert np.array_equal(e, fresh)
    assert float(np.linalg.norm(e)) > 0


# ---- DR-0: the secular frequency -------------------------------------------------------------------------------------
def _mathieu_beta(q: float, terms: int = 40) -> float:
    """Mathieu's characteristic exponent at a = 0: the beta of the Floquet solutions e^{i beta tau} P(tau), P of period
    pi, of y'' + 2q cos(2 tau) y = 0. The coefficients c_2j of P e^{i beta tau} satisfy
    (2j + beta)^2 c_2j = q (c_2j-2 + c_2j+2), so c_2/c_0 and c_-2/c_0 are continued fractions, and beta is the root in
    (0, 1) of beta^2 = q (c_2 + c_-2) / c_0 (McLachlan; Abramowitz and Stegun 20.3), found by bisection."""

    def balance(beta: float) -> float:
        up = down = 0.0
        for j in range(terms, 0, -1):
            up = q / ((2 * j + beta) ** 2 - q * up)
            down = q / ((beta - 2 * j) ** 2 - q * down)
        return beta * beta - q * (up + down)

    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if balance(mid) < 0 else (lo, mid)
    return (lo + hi) / 2


def _scheme_beta(q: float, steps: int) -> float:
    """The exponent the integrator itself gives an ideal quadrupole at Mathieu's q and `steps` steps an RF period: the
    trace of its one-period matrix of kicks (the field taken at each step's start) and drifts, as the knotted study's
    scheme_secular reads it."""
    dt, omega = 1 / steps, 2 * math.pi  # time in RF periods
    kick = q * omega**2 / 2 * dt  # |q/m| g dt, with q = 2 |q/m| g / omega^2
    cols = []
    for x, u in ((1.0, 0.0), (0.0, 1.0)):
        for n in range(steps):
            u -= kick * x * math.cos(omega * n * dt)
            x += u * dt
        cols.append((x, u))
    return math.acos((cols[0][0] + cols[1][1]) / 2) / math.pi


def test_the_reference_exponent_is_mathieus() -> None:
    # The study's RK4 monodromy put beta at 0.3737 at q = 0.5 and f_s at 231.12 MHz at q = 0.6 and 1 GHz (run.log),
    # and small q gives q / sqrt 2; at q = 0.3 the first order is 1.8 % low, so DR-0's 0.1 % tells the exact from it
    assert _mathieu_beta(0.5) == pytest.approx(0.37374412, abs=1e-8)
    assert _mathieu_beta(0.6) * 1e9 / 2 == pytest.approx(231.12e6, abs=0.01e6)
    assert _mathieu_beta(0.01) == pytest.approx(0.01 / math.sqrt(2), rel=1e-4)
    assert _mathieu_beta(0.3) / (0.3 / math.sqrt(2)) - 1 > 0.018


def _tracked_beta(steps_per_period: int, periods: int) -> pt.Arr:
    """Mathieu's beta as electrons at q = 0.3 in a straight SectionField of the facing boards (T915-S: h = 0.915 mm,
    915 MHz; the section sampled h/80 fine) show it, 1 eV along the path and 0.004h off the null along each diagonal
    (the principal axes) and along the split and the board normal: read stroboscopically, x(k+1) + x(k-1) =
    2 cos(mu) x(k) for any linear oscillation, cos mu a least-squares ratio, and beta = mu / pi."""
    volts = 4 * H915**2 * pt.gradient_for_q(0.3, F915) / math.pi  # FacingBoards: gradient pi V / 4h^2
    field = SectionField.straight(_facing(80, across=1.0, toward=0.9), 0.05, pt.InPhase(volts, OMEGA))
    root = 1 / math.sqrt(2)
    starts = 4e-3 * H915 * np.array([[root, root], [root, -root], [1.0, 0.0], [0.0, 1.0]])
    v = pt.v_of(1.0)
    u0 = np.tile([0.0, 0.0, v / math.sqrt(1 - (v / pt.C0) ** 2)], (4, 1))
    seen = [starts.copy()]

    def strobe(step: int, t: float, x: pt.Arr, u: pt.Arr, alive: pt.BoolArr) -> None:
        if (step + 1) % steps_per_period == 0:
            seen.append(x[:, :2].copy())

    flight = pt.integrate(
        np.column_stack([starts, np.zeros(4)]),
        u0,
        np.full(4, pt.QM_E),
        1 / (F915 * steps_per_period),
        steps_per_period * periods,
        field,
        [strobe],
    )
    assert bool(flight.alive.all())
    assert not field.off_grid.any()
    path = np.array(seen)
    inner, outer = path[1:-1], path[2:] + path[:-2]
    cos_mu = (inner * outer).sum(axis=(0, 2)) / (2 * (inner * inner).sum(axis=(0, 2)))
    out: pt.Arr = np.arccos(cos_mu) / math.pi
    return out


def test_dr_0_a_tracked_electron_keeps_mathieus_exact_secular_frequency() -> None:
    # DR-0: at 120 steps an RF period the tracked secular frequency beta f / 2 is Mathieu's exact one to 0.1 % for every
    # start. Along the principal axes, where the boards' potential has no quartic term, it is the integrator's own,
    # 1.2e-4 high, to 3e-5: the sweep, the grid and the 1 eV along the path add no more. At 30 steps an RF period the
    # same run misses by the integrator's 0.19 %, outside the criterion.
    exact = _mathieu_beta(0.3)
    fine = _tracked_beta(120, 60)
    assert float(np.abs(fine / exact - 1).max()) < 1e-3
    assert fine[:2] / exact - 1 == pytest.approx([_scheme_beta(0.3, 120) / exact - 1] * 2, abs=3e-5)
    coarse = _tracked_beta(30, 60)
    assert float(np.abs(coarse[:2] / exact - 1).min()) > 1e-3
    assert coarse[:2] == pytest.approx([_scheme_beta(0.3, 30)] * 2, rel=3e-5)


# ---- the drives -----------------------------------------------------------------------------------------------------
def _on_ring(radius: float, s: pt.Arr, a: pt.Arr, b: pt.Arr) -> tuple[pt.Arr, pt.Arr, pt.Arr, pt.Arr]:
    """Points a across (outward) and b down from a flat ring at arc lengths s, with the frame (theta = 0) and the
    tangent there."""
    angle = s / radius
    zero = np.zeros(len(s))
    radial = np.stack([np.cos(angle), np.sin(angle), zero], 1)
    down = np.tile([0.0, 0.0, -1.0], (len(s), 1))
    tangent = np.stack([-np.sin(angle), np.cos(angle), zero], 1)
    return radius * radial + a[:, None] * radial + b[:, None] * down, radial, down, tangent


@pytest.mark.parametrize("m", [1, -1, 2, -3])
def test_the_travelling_wave_runs_at_its_phase_velocity_and_carries_t_cross_e_over_v(m: int) -> None:
    # On the flat facing-board ring of 318 mm at 915 MHz, V0 cos(m 2 pi s / L - omega t + phase): what a point sees at
    # (s, t) a point v tau further along sees at t + tau, with v = omega L / (2 pi m), the 2.91e8 m/s of T915-S's m = 1
    # (c / 1.031); 0.1 % off that speed the pattern misses by more than a thousandth. B = T x E / v: |B| v = |E|,
    # B normal to E and to the path, and E x B along the wave's direction, T for m > 0 and -T for m < 0.
    radius = L915 / (2 * math.pi)
    drive = pt.TravellingWave(50.0, OMEGA, m, phase=0.4)
    field = SectionField.along(_facing(40), ClosedCurve.circle(radius, 4000), drive)
    speed = drive.phase_velocity(L915)
    assert speed == pytest.approx(F915 * L915 / m, rel=1e-12)
    rng = np.random.default_rng(10 + m)
    s = rng.uniform(0.0, L915, 300)
    a, b = rng.uniform(-0.5 * H915, 0.5 * H915, (2, 300))
    x, radial, down, tangent = _on_ring(radius, s, a, b)
    t0 = 0.3 / F915
    e0, b0 = field(x, t0)
    assert e0 is not None
    assert b0 is not None
    local = np.stack([(e0 * radial).sum(1), (e0 * down).sum(1)], 1)
    scale = float(np.abs(e0).max())
    for tau, factor, close in ((0.1, 1.0, True), (0.37, 1.0, True), (2.2, 1.0, True), (2.2, 1.001, False)):
        x1, radial1, down1, _ = _on_ring(radius, np.mod(s + factor * speed * tau / F915, L915), a, b)
        e1 = _e(field, x1, t0 + tau / F915)
        miss = float(np.abs(np.stack([(e1 * radial1).sum(1), (e1 * down1).sum(1)], 1) - local).max()) / scale
        assert (miss < 1e-6) if close else (miss > 1e-3)
    size_e, size_b = np.linalg.norm(e0, axis=1), np.linalg.norm(b0, axis=1)
    assert size_b * abs(speed) == pytest.approx(size_e, rel=1e-12)
    assert float(np.abs((e0 * b0).sum(1)).max()) < 1e-12 * float((size_e * size_b).max())
    assert float(np.abs((b0 * tangent).sum(1)).max()) < 1e-9 * float(size_b.max())
    flow = (np.cross(e0, b0) * tangent).sum(1) / (size_e**2 / speed)
    assert flow == pytest.approx(np.ones(300), rel=1e-9)


def test_a_standing_wave_keeps_its_nodes_and_an_in_phase_drive_is_the_same_all_along() -> None:
    # V0 cos(m angle + spatial) cos(omega t - phase) vanishes at its nodes at every time and at its antinodes swings
    # as -V grad phi, the section's own field at those volts; V0 cos(omega t - phase) gives every arc length the same.
    # Neither carries B.
    radius = L915 / (2 * math.pi)
    ring = ClosedCurve.circle(radius, 4000)
    grid = _facing(40)
    standing = SectionField.along(grid, ring, pt.StandingWave(50.0, OMEGA, 2, phase=0.3, spatial=0.5))
    rng = np.random.default_rng(5)
    a, b = rng.uniform(-0.5 * H915, 0.5 * H915, (2, 4))
    nodes = (np.array([0.5, 1.5, 2.5, 3.5]) * math.pi - 0.5) / 2 * radius  # 2 angle + 0.5 an odd multiple of pi/2
    tops = (np.array([0.0, 1.0, 2.0, 3.0]) * math.pi - 0.5) / 2 * radius  # ... a multiple of pi
    gu, gv, on = grid.gradient_at(a, b)  # the frame at theta = 0: u outward along a, v down along b
    assert bool(on.all())
    for t in np.array([0.0, 0.21, 0.55]) / F915:
        at_nodes, magnetic = standing(_on_ring(radius, nodes % L915, a, b)[0], t)
        assert at_nodes is not None
        assert magnetic is None
        x, radial, down, _ = _on_ring(radius, tops % L915, a, b)
        swing = 50.0 * math.cos(OMEGA * t - 0.3) * np.array([1.0, -1.0, 1.0, -1.0])
        e = _e(standing, x, t)
        assert float(np.abs(at_nodes).max()) < 1e-7 * float(np.abs(e).max())
        local = np.stack([(e * radial).sum(1), (e * down).sum(1)], 1)
        assert local == pytest.approx(-swing[:, None] * np.stack([gu, gv], 1), rel=1e-7, abs=1e-9)
    steady = SectionField.along(grid, ring, pt.InPhase(50.0, OMEGA, phase=0.7))
    s = rng.uniform(0.0, L915, 50)
    x, radial, down, _ = _on_ring(radius, s, np.full(50, a[0]), np.full(50, b[0]))
    even, magnetic = steady(x, 0.4 / F915)
    assert even is not None
    assert magnetic is None
    volts = 50.0 * math.cos(OMEGA * 0.4 / F915 - 0.7)
    local = np.stack([(even * radial).sum(1), (even * down).sum(1)], 1)
    assert local == pytest.approx(np.tile([-volts * gu[0], -volts * gv[0]], (50, 1)), rel=1e-7, abs=1e-9)


def test_a_travelling_wave_needs_an_order_and_a_frequency() -> None:
    with pytest.raises(ValueError, match="order"):
        pt.TravellingWave(1.0, OMEGA, 0)
    with pytest.raises(ValueError, match="omega"):
        pt.TravellingWave(1.0, 0.0, 1)


# ---- off the grid ---------------------------------------------------------------------------------------------------
def test_what_is_off_the_grid_feels_nothing_is_flagged_and_can_be_stopped() -> None:
    # Particles beyond the grid's ends, under its lowest row, in a solid cell, at NaN, or far from the knot get no
    # field, no B, and the flag; integrate stops them with stop_off_grid as an observer, a step after they leave.
    base = _rails(40)
    h = RailPair(*R1_RAILS).null_height
    cu, cv = np.meshgrid((base.u[1:] + base.u[:-1]) / 2, (base.v[1:] + base.v[:-1]) / 2, indexing="ij")
    copper = (np.abs(cu) > 0.5 * h) & (np.abs(cu) < 1.5 * h) & (cv < 0.6 * h)  # a block standing off the board
    grid = SectionGrid(base.u, base.v, base.phi, base.null, solid=copper)
    field = SectionField.straight(grid, 0.02, pt.TravellingWave(10.0, OMEGA, 1))
    section = np.array(
        [
            [0.0, 0.0],  # the null: on
            [3.5 * h, 0.0],  # beyond the grid's side
            [0.0, -0.9 * h],  # under the grid's lowest row
            [h, -0.6 * h],  # in the solid block
            [0.3 * h, 0.4 * h],  # on
            [math.nan, 0.0],
        ]
    )
    x = np.column_stack([section, np.full(len(section), 0.013)])
    e, magnetic = field(x, 0.1 / F915)
    assert e is not None
    assert magnetic is not None
    off = np.array([False, True, True, True, False, True])
    assert np.array_equal(field.off_grid, off)
    assert np.all(e[off] == 0.0)
    assert np.all(magnetic[off] == 0.0)
    assert np.all(np.linalg.norm(e[~off], axis=1) > 0)
    # far from a knot
    knot = _trefoil(1200)
    swept = SectionField.along(_facing(20), knot, pt.InPhase(1.0, OMEGA), theta=_closing(knot, 0))
    out = _e(swept, np.array([knot.points[7] + [0.2 * H915, 0.0, 0.0], [0.0, 0.0, 0.05]]))
    assert np.array_equal(swept.off_grid, [False, True])
    assert np.all(out[1] == 0.0)
    # stopped: one electron resting at the null stays, one launched across the split leaves the grid and is stopped
    v = pt.v_of(5.0)
    u0 = np.array([[0.0, 0.0, 0.0], [v, 0.0, 0.0]])
    still = SectionField.straight(base, 0.02, pt.InPhase(0.0, OMEGA))
    flight = pt.integrate(np.zeros((2, 3)), u0, np.full(2, pt.QM_E), 1e-11, 400, still, [still.stop_off_grid])
    assert np.array_equal(flight.alive, [True, False])
    assert 3 * h < float(flight.x[1, 0]) < 3 * h + 2 * v * 1e-11


# ---- what the search remembers --------------------------------------------------------------------------------------
@functools.cache
def _memory_case(name: str) -> tuple[ClosedCurve, SectionGrid, pt.Arr]:
    curve = _trefoil(1200) if name == "trefoil" else ClosedCurve.circle(1.5e-3, 300)
    return curve, _facing(20), _closing(curve, 0)


@settings(max_examples=40, deadline=None)
@given(st.sampled_from(["trefoil", "a ring tighter than the grid"]), st.integers(0, 2**32 - 1))
def test_what_the_search_remembers_between_calls_changes_nothing(name: str, seed: int) -> None:
    # Ten calls on a dozen particles that drift within the window, run past it, jump along the curve, fly far off and
    # back, stay put or turn to NaN: each call's field and flags are, to the bit, those a search of every sample gives.
    # On the trefoil the window is trusted out to its thickness (8.7 mm, half its closest approach to itself); on a
    # 1.5 mm ring the grid reaches past the thickness, and the particles beyond it are searched in full.
    curve, grid, theta = _memory_case(name)
    field = SectionField.along(grid, curve, pt.InPhase(1.0, 0.0), theta=theta)
    rng = np.random.default_rng(seed)
    count, step = 12, curve.length / len(curve)
    x = curve.points[rng.integers(0, len(curve), count)] + rng.normal(0.0, 0.5 * H915, (count, 3))
    for _ in range(10):
        kind = rng.integers(0, 6, count)
        moved = x.copy()
        moved[kind == 0] += rng.normal(0.0, step, (int((kind == 0).sum()), 3))
        moved[kind == 1] += rng.normal(0.0, 6 * step, (int((kind == 1).sum()), 3))
        jumps = int((kind == 2).sum())
        moved[kind == 2] = curve.points[rng.integers(0, len(curve), jumps)] + rng.normal(0.0, 1.5 * H915, (jumps, 3))
        moved[kind == 3] = rng.normal(0.0, 0.2, (int((kind == 3).sum()), 3))
        moved[(kind == 4) & (rng.uniform(size=count) < 0.3)] = math.nan
        x = moved
        e = _e(field, x)
        where = field.locate(x)
        fresh, _ = field.field_at(where, 0.0)
        assert np.array_equal(e, fresh)
        assert np.array_equal(field.off_grid, ~grid.gradient_at(where.u, where.v)[2])


def test_between_the_strands_of_a_knot_tighter_than_the_grid_the_search_is_in_full() -> None:
    # The trefoil shrunk to 55 mm, half its closest approach to itself 1.5 mm, under the grid's 2.0 mm reach: particles
    # first found near one strand and then carried across, past the midpoint toward the strand passing nearest, are
    # still within the grid's reach of the first, and get the second's field, to the bit, as a search of every sample
    # gives: the window round the first is trusted only inside the curve's thickness
    knot = ClosedCurve.torus_knot(1.0, 0.45, 600)
    curve = ClosedCurve(knot.points * (0.055 / knot.length), 0.055)
    grid = _facing(20)
    field = SectionField.along(grid, curve, pt.InPhase(1.0, 0.0), theta=_closing(curve, 0))
    pts, n = curve.points, len(curve)
    apart = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    apart = np.minimum(apart, n - apart) * curve.length / n
    gap = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=2)
    gap = np.where(apart > math.pi * curve.bend_radius, gap, np.inf)
    i, j = np.unravel_index(int(np.argmin(gap)), gap.shape)
    closest = float(gap[i, j])
    assert closest / 2 < grid.reach < closest
    a, b = pts[i], pts[j]
    near = np.repeat(a[None, :] + 0.1 * (b - a), 8, axis=0)
    across = a[None, :] + np.linspace(0.55, 0.9, 8)[:, None] * (b - a)[None, :]
    _e(field, near)  # remembered on the first strand
    e = _e(field, across)
    where = field.locate(across)
    fresh, _ = field.field_at(where, 0.0)
    assert np.array_equal(e, fresh)
    assert np.abs(e).max() > 0
    steps = np.minimum(np.abs(where.segment - j), n - np.abs(where.segment - j))
    assert int(steps.max()) <= 2  # on the strand passing nearest


# ---- the grid and the frame -----------------------------------------------------------------------------------------
def test_a_section_unfolds_from_its_quarter_and_finds_its_null() -> None:
    # The knotted study solves its facing strips on the quarter (odd in both) and one board's rails on a half (even
    # across): unfolded by their parities they are the whole section, solid cells with them, and find_null puts the
    # rails' null where the grid's differences do, within 4e-4 h of ELE-1's, from a start 0.2h away.
    boards, rails_pair = FacingBoards(2 * H915), RailPair(*R1_RAILS)
    h = rails_pair.null_height
    full = _facing(20)
    corner = np.zeros((40, 19), dtype=bool)
    corner[30:, 15:] = True  # cells near the top board, out along the strip
    u_quarter, v_quarter = _even(0.0, 2 * H915, H915 / 20), _even(0.0, 0.95 * H915, H915 / 20)
    quarter = SectionGrid.mirrored(
        u_quarter,
        v_quarter,
        boards.potential(*np.meshgrid(u_quarter, v_quarter, indexing="ij"), 1.0),
        (0.0, 0.0),
        u_parity=-1,
        v_parity=-1,
        solid=corner,
    )
    assert np.allclose(quarter.u, full.u, rtol=0.0, atol=1e-18)
    assert np.allclose(quarter.v, full.v, rtol=0.0, atol=1e-18)
    assert np.allclose(quarter.phi, full.phi, rtol=0.0, atol=1e-15)
    assert quarter.solid is not None
    assert quarter.solid.shape == (80, 38)
    assert bool(quarter.solid[75, 35])
    assert bool(quarter.solid[4, 2])
    assert bool(quarter.solid[75, 2])
    assert not bool(quarter.solid[40, 19])
    rails = _rails(40)
    u_half = _even(0.0, 3 * h, h / 40)
    half = SectionGrid.mirrored(
        u_half,
        rails.v,
        rails_pair.potential(*np.meshgrid(u_half, rails.v, indexing="ij"), 1.0),
        rails.null,
        u_parity=1,
    )
    assert np.allclose(half.phi, rails.phi, rtol=0.0, atol=1e-15)
    u_null, v_null = rails.find_null((0.15 * h, 1.15 * h))
    assert abs(u_null) < 1e-12 * h
    assert abs(v_null - h) < 4e-4 * h
    gu, gv, on = rails.gradient_at(np.array([u_null]), np.array([v_null]))
    assert bool(on[0])
    assert math.hypot(gu[0], gv[0]) < 1e-9 / h
    assert float(rails.potential_at(np.array([0.0]), np.array([h]))[0]) == pytest.approx(
        float(RailPair(*R1_RAILS).potential(0.0, h, 1.0)), abs=1e-4
    )


def test_a_section_grid_refuses_what_it_cannot_hold() -> None:
    u, v = np.linspace(-1e-3, 1e-3, 5), np.linspace(0.0, 1e-3, 4)
    phi = np.zeros((5, 4))
    with pytest.raises(ValueError, match="phi"):
        SectionGrid(u, v, phi[:, :3], (0.0, 5e-4))
    with pytest.raises(ValueError, match="increase"):
        SectionGrid(u[::-1], v, phi, (0.0, 5e-4))
    with pytest.raises(ValueError, match="null"):
        SectionGrid(u, v, phi, (0.0, 2e-3))
    with pytest.raises(ValueError, match="finite"):
        SectionGrid(u, v, np.full((5, 4), math.nan), (0.0, 5e-4))
    with pytest.raises(ValueError, match="cells"):
        SectionGrid(u, v, phi, (0.0, 5e-4), solid=np.zeros((5, 4), dtype=bool))
    with pytest.raises(ValueError, match="parity"):
        SectionGrid.mirrored(u, v, phi, (0.0, 5e-4), u_parity=-1)  # u does not start at 0


def test_normals_given_per_sample_are_the_framing_they_make() -> None:
    # The first normal given at each sample, with some tangent mixed in and any length, is the framing angle it makes
    # from the twist-free axes; the frame must be unit normals with second = T x first.
    knot = _trefoil(1200)
    theta = _closing(knot, -1)
    e1, e2 = knot.axes
    given_normals = 2.5 * (np.cos(theta)[:, None] * e1 + np.sin(theta)[:, None] * e2) + 0.3 * knot.tangents
    by_angle = SectionField.along(_facing(20), knot, pt.InPhase(1.0, 0.0), theta=theta)
    by_normals = SectionField.along(_facing(20), knot, pt.InPhase(1.0, 0.0), normals=given_normals)
    rng = np.random.default_rng(2)
    x = knot.points[rng.integers(0, 1200, 200)] + rng.normal(0.0, 0.4 * H915, (200, 3))
    assert np.allclose(_e(by_angle, x), _e(by_normals, x), rtol=0.0, atol=1e-9)
    with pytest.raises(ValueError, match="frame"):
        SectionField(
            by_angle.grid, knot.points, knot.tangents, by_angle.second, by_angle.first, L915, pt.InPhase(1.0, 0.0)
        )
