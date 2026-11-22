"""The opt-in bend correction of a swept section (array16/particles/sections.py, BendShift), held to what it must do:
off by default and on a straight path; on a ring of R = 40h, the null moved toward the centre by 2 ln 2/pi^2 h^2/R with
the bend in the boards' plane and by h^2/(pi R) with it toward them, and by their superposition at any angle between;
the field near the null the section's own, moved rigidly, so the gradient there unchanged, as the axisymmetric
reference has it to first order; and on T915-S's trefoil, framed k = 0 and k = -1, the null moved by h^2 (0.1405
kappa_u, 0.318 kappa_v) at every station, across the seam too."""

from __future__ import annotations

import functools
import math

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from array16 import particles as pt
from array16.guide import axisym
from array16.particles import ClosedCurve, FacingBoards, SectionField, SectionGrid
from array16.particles.sections import BEND_IN_PLANE, BEND_TOWARD, BendShift

H915 = 0.915e-3  # T915-S: null to copper (m)
R40 = 40 * H915  # a ring of R = 40h, E2-2's radius (m)
L915 = 0.318  # T915-S's trefoil (m)


@functools.cache
def _facing(div: int = 40) -> SectionGrid:
    """ELE-2's potential per volt round the facing boards' null, h / div apart, 2h either side of the split and 0.95h
    toward the boards; its null is the node at the origin, where the differences of an odd potential vanish."""
    boards = FacingBoards(2 * H915)
    u = np.linspace(-2 * H915, 2 * H915, 4 * div + 1)
    v = np.linspace(-0.95 * H915, 0.95 * H915, round(1.9 * div) + 1)
    return SectionGrid.sampled(lambda a, b: boards.potential(a, b, 1.0), u, v, (0.0, 0.0))


@functools.cache
def _ring(samples: int = 2000) -> ClosedCurve:
    return ClosedCurve.circle(R40, samples)


@functools.cache
def _trefoil(n: int = 3000) -> ClosedCurve:
    knot = ClosedCurve.torus_knot(1.0, 0.45, n)
    return ClosedCurve(knot.points * (L915 / knot.length), L915)


def _closing(curve: ClosedCurve, k: int, theta0: float = 0.3) -> pt.Arr:
    """The framing angle closing after a lap with k half turns (the knotted study's rule), turned a quarter so that u
    runs across the boards' plane."""
    s = np.arange(len(curve)) * curve.length / len(curve)
    out: pt.Arr = theta0 + math.pi / 2 + (k * math.pi - curve.holonomy) / curve.length * s
    return out


def _swept_null(field: SectionField, centre: pt.Arr, first: pt.Arr, second: pt.Arr) -> pt.Arr:
    """For each centre on the path, the swept field's null in the plane spanned by the frame there, by Newton's method
    on (E . first, E . second): its offset (a, b) from the centre along first and second (m), (N, 2)."""
    ab = np.zeros((len(centre), 2))
    step = 1e-5 * H915

    def across(ab: pt.Arr) -> pt.Arr:
        e, _ = field(centre + ab[:, :1] * first + ab[:, 1:] * second, 0.0)
        assert e is not None
        out: pt.Arr = np.stack([(e * first).sum(1), (e * second).sum(1)], 1)
        return out

    for _ in range(30):
        here = across(ab)
        jac = np.stack(
            [(across(ab + np.array([step, 0.0])) - here) / step, (across(ab + np.array([0.0, step])) - here) / step],
            2,
        )
        move = np.linalg.solve(jac, -here[:, :, None])[:, :, 0]
        ab = ab + move
        if float(np.abs(move).max()) < 1e-14 * H915:
            break
    return ab


def _ring_frame(theta: float, count: int = 13) -> tuple[pt.Arr, pt.Arr, pt.Arr, pt.Arr]:
    """Stations on the ring of 40h, on samples and between them: their places, the frame turned by theta from the
    twist-free axes (e1 the outward radius, e2 = T x e1 down the ring's axis), and the outward radius."""
    angle = np.linspace(0.0, 2 * math.pi, count, endpoint=False) + 0.37 * 2 * math.pi / 2000
    zero = np.zeros(count)
    out = np.stack([np.cos(angle), np.sin(angle), zero], 1)
    down = np.tile([0.0, 0.0, -1.0], (count, 1))
    first = math.cos(theta) * out + math.sin(theta) * down
    second = math.cos(theta) * down - math.sin(theta) * out
    return R40 * out, first, second, out


# ---- off by default --------------------------------------------------------------------------------------------------
def test_the_correction_is_off_by_default_and_refuses_what_it_does_not_take() -> None:
    # Without `bend` the swept field is the section's and its null on the path: what stage 1 tracked. A BendShift takes
    # a positive h and finite coefficients, the facing boards' by default; a swept field takes a BendShift or None.
    plain = SectionField.along(_facing(), _ring(), pt.InPhase(1.0, 0.0), theta=math.pi)
    assert plain.bend is None
    centre, first, second, _ = _ring_frame(math.pi)
    du, dv = plain.null_shift(plain.locate(centre))
    assert not du.any()
    assert not dv.any()
    assert float(np.abs(_swept_null(plain, centre, first, second)).max()) < 1e-8 * H915
    shift = BendShift(H915)
    assert (shift.inplane, shift.toward) == (BEND_IN_PLANE, BEND_TOWARD)
    for bad in (0.0, -H915, math.nan, math.inf):
        with pytest.raises(ValueError, match="h > 0"):
            BendShift(bad)
    parts = (plain.grid, plain.points, plain.tangents, plain.first, plain.second, plain.length, plain.drive)
    with pytest.raises(ValueError, match="BendShift or None"):
        SectionField(*parts, bend=0.1)  # type: ignore[arg-type]


def test_a_straight_path_has_no_bend_to_correct() -> None:
    # A straight sweep's curvature is zero, and with a BendShift its field is the uncorrected one to the bit.
    line = SectionField.straight(_facing(), 0.02, pt.InPhase(1.0, 0.0))
    parts = (line.grid, line.points, line.tangents, line.first, line.second, line.length, line.drive, line.axis)
    bent = SectionField(*parts, bend=BendShift(H915))
    assert not bent.curvature_vector.any()
    rng = np.random.default_rng(4)
    x = np.column_stack([rng.uniform(-0.5 * H915, 0.5 * H915, (50, 2)), rng.uniform(0.0, 0.05, 50)])
    e0, _ = line(x, 0.0)
    e1, _ = bent(x, 0.0)
    assert e0 is not None
    assert e1 is not None
    assert np.array_equal(e0, e1)


# ---- a ring of 40h ---------------------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("theta", "way"),
    [(math.pi, "in-plane"), (0.0, "in-plane"), (math.pi / 2, "toward"), (-math.pi / 2, "toward")],
)
def test_on_a_ring_of_forty_h_the_null_moves_toward_the_centre_by_the_first_order_shift(
    theta: float, way: axisym.Way
) -> None:
    # The ring of R = 40h (36.6 mm at h = 0.915 mm) framed with its bend in the boards' plane (theta pi: u inward, v up,
    # stage 1's flat ring; theta 0: u outward) or toward them (theta +-pi/2: u along the axis, v radial): the swept
    # null lies toward the centre by 0.140461 h^2/R = 3.2 um or 0.318310 h^2/R = 7.3 um, the axisymmetric reference's
    # first order, to 1e-5 of that (the ring's curvature by second differences is 1/R to 8e-7, its cubic between
    # samples off the circle by 4e-12 m), and nowhere off the line to the centre; the reference itself, solved at
    # 40h, is within 0.02 % of it.
    field = SectionField.along(_facing(), _ring(), pt.InPhase(1.0, 0.0), theta=theta, bend=BendShift(H915))
    centre, first, second, outward = _ring_frame(theta)
    ab = _swept_null(field, centre, first, second)
    moved = ab[:, :1] * first + ab[:, 1:] * second
    inward = -(moved * outward).sum(1)
    expected = axisym.FIRST_ORDER[way] * H915**2 / R40
    assert expected == pytest.approx((3.2e-6 if way == "in-plane" else 7.3e-6), abs=0.05e-6)
    assert inward == pytest.approx(np.full(len(centre), expected), rel=1e-5)
    assert float(np.abs(moved + inward[:, None] * outward).max()) < 1e-5 * expected
    reference = axisym.bend(40.0, way)
    assert reference.shift * H915 == pytest.approx(expected, rel=2e-4)


@settings(max_examples=30, deadline=None)
@given(st.floats(-math.pi, math.pi))
def test_a_bend_between_the_plane_and_the_boards_moves_the_null_by_the_superposition(theta: float) -> None:
    # Framed at any angle theta, the ring's curvature 1/R toward its centre lies at psi = pi - theta from the boards'
    # plane, (cos psi, sin psi)/R along (u, v): the null moves by h^2/R (0.1405 cos psi, 0.318 sin psi) in (u, v),
    # (0.1405 cos^2 psi + 0.318 sin^2 psi) h^2/R toward the centre and (0.318 - 0.1405) sin psi cos psi h^2/R across
    # that, up the ring's axis.
    field = SectionField.along(_facing(), _ring(), pt.InPhase(1.0, 0.0), theta=theta, bend=BendShift(H915))
    centre, first, second, outward = _ring_frame(theta, count=3)
    ab = _swept_null(field, centre, first, second)
    psi, scale = math.pi - theta, H915**2 / R40
    tolerance = 1e-5 * BEND_TOWARD * scale
    assert ab == pytest.approx(
        np.tile([BEND_IN_PLANE * math.cos(psi) * scale, BEND_TOWARD * math.sin(psi) * scale], (3, 1)), abs=tolerance
    )
    moved = ab[:, :1] * first + ab[:, 1:] * second
    toward_centre = (BEND_IN_PLANE * math.cos(psi) ** 2 + BEND_TOWARD * math.sin(psi) ** 2) * scale
    up_the_axis = (BEND_TOWARD - BEND_IN_PLANE) * math.sin(psi) * math.cos(psi) * scale
    assert -(moved * outward).sum(1) == pytest.approx(np.full(3, toward_centre), abs=tolerance)
    assert moved[:, 2] == pytest.approx(np.full(3, up_the_axis), abs=tolerance)


def test_the_field_near_the_null_is_the_sections_moved_rigidly_and_its_gradient_unchanged() -> None:
    # In each sample's normal plane, the corrected field at x is the uncorrected field at x less the null's
    # displacement, to rounding: the section moved as a whole. So the gradient at the moved null is the section's,
    # pi V/4h^2 to the grid's 1e-3, as the axisymmetric reference has it to first order: there the bend changes it by
    # -5.2e-5 in the boards' plane and +4.3e-5 toward them at 40h, second order in h/R.
    grid = _facing()
    cases: tuple[tuple[float, axisym.Way], ...] = ((math.pi, "in-plane"), (math.pi / 2, "toward"))
    for theta, way in cases:
        plain = SectionField.along(grid, _ring(), pt.InPhase(1.0, 0.0), theta=theta)
        bent = SectionField.along(grid, _ring(), pt.InPhase(1.0, 0.0), theta=theta, bend=BendShift(H915))
        rng = np.random.default_rng(11)
        index = rng.integers(0, 2000, 40)
        a, b = rng.uniform(-0.6 * H915, 0.6 * H915, (2, 40))
        x = bent.points[index] + a[:, None] * bent.first[index] + b[:, None] * bent.second[index]
        du, dv = bent.null_shift(bent.locate(x))
        back = x - du[:, None] * bent.first[index] - dv[:, None] * bent.second[index]
        e_bent, _ = bent(x, 0.0)
        e_plain, _ = plain(back, 0.0)
        assert e_bent is not None
        assert e_plain is not None
        assert np.allclose(e_bent, e_plain, rtol=0.0, atol=1e-12 * float(np.abs(e_plain).max()))
        centre, first, second, _ = _ring_frame(theta, count=1)
        null = centre + _swept_null(bent, centre, first, second) @ np.stack([first[0], second[0]])
        step = 1e-3 * H915
        e_u, _ = bent(null + step * first, 0.0)
        assert e_u is not None
        gradient = float(np.linalg.norm(e_u)) / step
        assert gradient == pytest.approx(math.pi / (4 * H915**2), rel=1e-3)
        assert abs(axisym.bend(40.0, way).gradient - 1) < 6e-5


# ---- T915-S ----------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("k", [0, -1])
def test_on_t915s_the_null_moves_by_the_curvature_in_the_frame_across_the_seam_too(k: int) -> None:
    # T915-S's trefoil framed k = 0 (the frame closes) and k = -1 (it comes back turned by pi): at stations round the
    # lap, between samples and on both sides of the seam, the swept null moves by h^2 (0.1405 kappa_u, 0.318 kappa_v),
    # the curvature vector taken in the frame there, from where it lies without the correction (on the path, 8.5e-8 m
    # from the chord between samples these stations are taken on) to 1e-4 of the largest shift; at the tightest bend
    # (15.3 mm) that is 7.7 to 17.4 um, the bend shared between the boards' plane and the boards as the framing turns.
    knot = _trefoil()
    theta = _closing(knot, k)
    plain = SectionField.along(_facing(), knot, pt.InPhase(1.0, 0.0), theta=theta)
    field = SectionField.along(_facing(), knot, pt.InPhase(1.0, 0.0), theta=theta, bend=BendShift(H915))
    assert abs(math.remainder(field.seam - k * math.pi, 2 * math.pi)) < 1e-9
    kappa = field.curvature_vector
    assert np.linalg.norm(kappa, axis=1) == pytest.approx(knot.curvature, rel=1e-6)
    tight = int(np.argmax(knot.curvature))
    pick = np.unique(np.concatenate([np.arange(0, 3000, 211), [tight, 2998, 2999, 0, 1]]))
    fraction = 0.4  # a station between each picked sample and the next
    nxt = (pick + 1) % 3000
    centre = (1 - fraction) * knot.points[pick] + fraction * knot.points[nxt]
    where = field.locate(centre)
    first, second = where.first, where.second
    moved = _swept_null(field, centre, first, second) - _swept_null(plain, centre, first, second)
    du, dv = field.null_shift(where)
    largest = float(np.hypot(du, dv).max())
    assert moved == pytest.approx(np.stack([du, dv], 1), abs=1e-4 * largest)
    k_u, k_v = (kappa[pick] * first).sum(1), (kappa[pick] * second).sum(1)
    assert du == pytest.approx(BEND_IN_PLANE * H915**2 * k_u, rel=0.02, abs=0.02 * largest)
    assert dv == pytest.approx(BEND_TOWARD * H915**2 * k_v, rel=0.02, abs=0.02 * largest)
    at = int(np.flatnonzero(pick == tight)[0])
    size = math.hypot(du[at], dv[at])
    radius = knot.bend_radius
    assert radius == pytest.approx(15.3e-3, abs=0.05e-3)
    assert BEND_IN_PLANE * H915**2 / radius < size < BEND_TOWARD * H915**2 / radius
