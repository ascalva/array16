"""The geometry of a guide's path (array16/guide/geometry.py), held to the knotted study's numbers and to what stage 2's
paths for H1 must be: T915-S framed k = 0 with its twist of +0.48 deg/mm, its in-plane radius of 18.57 mm, its plates of
329.5 and 306.1 mm and their flat boxes, as experiments/knotted_guide/run.log has them; the framing the twin's to the
bit; KNT-10's plate lengths; the ellipse of 318 mm whose tightest radius is T915-S's in-plane one, the flat ring framed
with its bend wholly toward its boards, and the (1, 3) helix at r/R = 0.1, 0.2 and 0.3, by length and tightest radius;
and the flat plates laid out as they must lie."""

from __future__ import annotations

import math

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from array16 import particles as pt
from array16 import twin
from array16.guide import geometry
from array16.particles import C0, ClosedCurve, SectionField, SectionGrid

H915 = 0.915e-3  # T915-S: null to copper (m)
WIDTH915 = 8 * H915  # its strips' width
L915 = C0 / (1.031 * 915e6)  # the study's T915-S: m = 1 at 915 MHz on index 1.031, 317.79 mm
# the study's full run (n = 4000; experiments/knotted_guide/results.json, strips, k = 0)
STUDY = {
    "twist_deg_per_mm": 0.48151253607767686,
    "inplane_radius_mm": 18.57179826950113,
    "plate_plus_mm": 329.5090481770578,
    "plate_minus_mm": 306.0888940184844,
    "integral_kappa_n": -12.798300754977106,
}
BOXES = {"plus": (21.267443790600844, 314.8974360231517), "minus": (20.711088712607346, 291.85531148672015)}


def _t915s(n: int = 4000) -> ClosedCurve:
    return geometry.torus(geometry.T915S_RATIO, L915, n)


# ---- T915-S as the study framed it ----------------------------------------------------------------------------------
def test_t915s_framed_k_0_is_the_studys() -> None:
    # The r/R = 0.45 trefoil of 317.8 mm at 4000 points, framed k = 0 from the start that keeps the largest in-plane
    # radius: twist +0.48 deg/mm, 18.57 mm in the boards' plane (its tightest bend 15.29 mm), plates 329.5 and 306.1 mm
    # at +-h, laid flat in boxes 21 x 315 and 21 x 292 mm with no net turning: the study's numbers to 1e-9.
    knot = _t915s()
    frame = geometry.framing(knot, 0)
    assert knot.length == pytest.approx(0.3177905, abs=1e-7)
    assert knot.bend_radius == pytest.approx(15.29e-3, abs=0.01e-3)
    assert frame.twist_deg_per_mm == pytest.approx(STUDY["twist_deg_per_mm"], rel=1e-9)
    assert round(frame.twist_deg_per_mm, 2) == 0.48
    assert frame.inplane_radius / 1e-3 == pytest.approx(STUDY["inplane_radius_mm"], rel=1e-9)
    assert round(frame.inplane_radius / 1e-3, 1) == 18.6
    integral = float(frame.kappa_n.sum() * knot.length / len(knot))
    assert integral == pytest.approx(STUDY["integral_kappa_n"], rel=1e-9)
    for sign, name in ((+1, "plus"), (-1, "minus")):
        plate = geometry.plate(knot, frame.normal, sign * H915)
        assert plate.length / 1e-3 == pytest.approx(STUDY[f"plate_{name}_mm"], rel=1e-9)
        assert abs(plate.turning) < 1e-9
        laid = geometry.develop(plate.dl, plate.kg)
        box = geometry.smallest_box(geometry.outline(laid.xy, laid.heading, WIDTH915)[::4])
        short, long = BOXES[name]
        assert (box.short / 1e-3, box.long / 1e-3) == pytest.approx((short, long), rel=1e-9)
    assert round(geometry.plate(knot, frame.normal, H915).length / 1e-3, 1) == 329.5
    assert round(geometry.plate(knot, frame.normal, -H915).length / 1e-3, 1) == 306.1


@pytest.mark.parametrize(
    ("ratio", "p", "n", "k"), [(0.45, 2, 4000, 0), (0.45, 2, 3000, -1), (0.2, 1, 3000, 0), (0.3, 1, 3000, 0)]
)
def test_the_framing_is_the_twins_to_the_bit(ratio: float, p: int, n: int, k: int) -> None:
    # The twin sweeps T915-S and the helix with array16.twin.framing; this one gives the same arrays to the bit, so the
    # twin can take it from here and stage 1 keeps its numbers.
    curve = geometry.torus(ratio, L915, n, p, 3)
    ours, theirs = geometry.framing(curve, k), twin.framing(curve, k)
    assert (ours.k, ours.theta0, ours.twist) == (theirs.k, theirs.theta0, theirs.twist)
    for name in ("normal", "across", "kappa_g", "kappa_n"):
        assert np.array_equal(getattr(ours, name), getattr(theirs, name))
    assert ours.inplane_radius == theirs.inplane_radius


@pytest.mark.parametrize("k", [0, -1, -3])
def test_the_plates_differ_by_2h_times_the_curvature_toward_the_boards(k: int) -> None:
    # KNT-10 (docs/KNOTTED_GUIDE.md section 11; the test takes the law's name, test_knt_10_..., when the law joins
    # docs/LAWS.md, which tests/test_laws.py holds every test so named to): the plates at +-h along the board normal
    # differ in length by -2h oint kappa_n ds, to under a micrometre on T915-S (23.42 mm framed k = 0, 13.77 mm at
    # k = -1, 7 um at k = -3). The h^2 terms cancel in the difference, and what is left is h^3 oint tau_g^2 kappa_n ds,
    # tau_g the normal's turning across the boards' plane: -0.69 um at k = 0, which the plates give to 10 % (the rest is
    # h^5 and the samples' spacing).
    knot = _t915s()
    frame = geometry.framing(knot, k)
    step = knot.length / len(knot)
    integral = float(frame.kappa_n.sum() * step)
    apart = geometry.plate(knot, frame.normal, H915).length - geometry.plate(knot, frame.normal, -H915).length
    assert abs(apart + 2 * H915 * integral) < 1e-6
    if k == 0:
        assert apart == pytest.approx(23.42e-3, abs=0.01e-3)
        turning = (np.roll(frame.normal, -1, 0) - np.roll(frame.normal, 1, 0)) / (2 * step)
        tau_g = -(turning * frame.across).sum(1)
        third = H915**3 * float((tau_g**2 * frame.kappa_n).sum() * step)
        assert third == pytest.approx(-0.69e-6, abs=0.01e-6)
        assert apart + 2 * H915 * integral == pytest.approx(third, rel=0.1)


# ---- stage 2's paths for H1 -----------------------------------------------------------------------------------------
def test_the_ellipse_has_t915s_length_and_its_tightest_in_plane_radius() -> None:
    # H1's planar path with a changing bend: the ellipse of 317.8 mm whose tightest radius is T915-S's in-plane 18.57 mm
    # (semi-axes 64.34 and 34.57 mm, its widest radius 119.7 mm), its 4000 points on the ellipse to 1e-12, equally
    # spaced to 1e-6, its second differences' tightest radius the given one to 1e-5; framed in the boards' plane, its
    # whole bend lies there and none toward the boards, with no twist.
    tightest = geometry.framing(_t915s(), 0).inplane_radius
    curve = geometry.ellipse(L915, tightest, 4000)
    a, b = float(curve.points[:, 0].max()), float(curve.points[:, 1].max())
    assert (a, b) == pytest.approx((64.336e-3, 34.566e-3), abs=1e-6)
    assert b * b / a == pytest.approx(tightest, rel=1e-12)
    assert a * a / b == pytest.approx(119.74e-3, abs=0.01e-3)
    assert curve.length == pytest.approx(L915, rel=1e-14)
    assert float(np.abs((curve.points[:, 0] / a) ** 2 + (curve.points[:, 1] / b) ** 2 - 1).max()) < 1e-12
    steps = np.linalg.norm(np.roll(curve.points, -1, 0) - curve.points, axis=1)
    assert float(steps.max() / steps.min() - 1) < 1e-6
    assert curve.bend_radius == pytest.approx(tightest, rel=1e-5)
    assert int(np.argmax(curve.curvature)) in (0, 1999, 2000, 3999)  # the ends of the major axis
    frame = geometry.flat(curve, "in-plane")
    assert frame.theta0 == pytest.approx(math.pi / 2, abs=1e-12)
    assert frame.inplane_radius == pytest.approx(tightest, rel=1e-5)
    assert float(np.abs(frame.kappa_n).max()) < 1e-9
    assert abs(frame.twist) < 1e-12
    with pytest.raises(ValueError, match="tightest radius"):
        geometry.ellipse(L915, L915 / (2 * math.pi), 400)


@settings(max_examples=25, deadline=None)
@given(st.floats(0.05, 1.0), st.floats(0.08, 0.95))
def test_an_ellipse_is_made_to_any_length_and_tightest_radius(length: float, share: float) -> None:
    # Any perimeter and any tightest radius under the circle's of that perimeter (a share of it from 0.08 up): the curve
    # has both, the semi-axes on its samples (a count divisible by four puts samples at the axes' ends), b^2/a to 1e-12,
    # and its second differences' tightest radius on 2000 points to 1e-3.
    tightest = share * length / (2 * math.pi)
    curve = geometry.ellipse(length, tightest, 2000)
    assert curve.length == pytest.approx(length, rel=1e-13)
    a, b = float(curve.points[0, 0]), float(curve.points[500, 1])
    assert b * b / a == pytest.approx(tightest, rel=1e-12)
    assert curve.bend_radius == pytest.approx(tightest, rel=1e-3)


def test_the_flat_ring_framed_toward_its_boards_bends_wholly_on_their_stiff_side() -> None:
    # H1b's path: the ring of 317.8 mm (R 50.58 mm) framed with its boards two coaxial cylinders, the normal radial and
    # u along the axis: no curvature in the boards' plane, all of it, 1/R, toward them. Swept with normals = across,
    # u runs down the ring's axis and v inward, the frame a design file's theta_deg = 90 gives a circle. Framed in the
    # plane instead, it is stage 1's flat ring (theta_deg = 180 in designs/flat_ring_318.toml: u inward, v up).
    curve = geometry.ring(L915, 4000)
    radius = L915 / (2 * math.pi)
    assert curve.length == pytest.approx(L915, rel=1e-15)
    assert curve.bend_radius == pytest.approx(radius, rel=1e-6)
    toward = geometry.flat(curve, "toward")
    assert toward.theta0 == 0.0
    assert toward.inplane_radius > 1e9
    assert np.abs(toward.kappa_n) == pytest.approx(np.full(4000, 1 / radius), rel=1e-6)
    assert np.abs(toward.across[:, 2]) == pytest.approx(np.ones(4000), abs=1e-12)
    grid = SectionGrid.sampled(lambda u, v: u * v, np.linspace(-1e-3, 1e-3, 5), np.linspace(-1e-3, 1e-3, 5), (0.0, 0.0))
    swept = SectionField.along(grid, curve, pt.InPhase(1.0, 0.0), normals=toward.across)
    assert np.abs(swept.first[:, 2]) == pytest.approx(np.ones(4000), abs=1e-12)
    inward = -curve.points / np.linalg.norm(curve.points, axis=1)[:, None]
    assert (swept.second * inward).sum(1) == pytest.approx(np.ones(4000), abs=1e-12)
    quarter = SectionField.along(grid, curve, pt.InPhase(1.0, 0.0), theta=math.pi / 2)  # a design's theta_deg = 90
    assert np.allclose(swept.first, quarter.first, atol=1e-12)
    assert np.allclose(swept.second, quarter.second, atol=1e-12)
    in_plane = geometry.flat(curve, "in-plane")
    assert in_plane.theta0 == pytest.approx(math.pi / 2, abs=1e-12)
    assert float(np.abs(in_plane.kappa_n).max()) < 1e-9
    stage1 = SectionField.along(grid, curve, pt.InPhase(1.0, 0.0), theta=math.pi)
    ours = SectionField.along(grid, curve, pt.InPhase(1.0, 0.0), normals=in_plane.across)
    assert np.allclose(ours.first, stage1.first, atol=1e-12)
    assert np.allclose(ours.second, stage1.second, atol=1e-12)
    with pytest.raises(ValueError, match="not flat"):
        geometry.flat(_t915s(400), "toward")


@pytest.mark.parametrize(
    ("ratio", "bend_mm", "inplane_mm"), [(0.1, 31.478, 60.085), (0.2, 25.049, 29.144), (0.3, 20.155, 21.062)]
)
def test_the_helix_at_three_ratios_has_t915s_length(ratio: float, bend_mm: float, inplane_mm: float) -> None:
    # H1's paths out of the plane without a knot: the (1, 3) torus curve at r/R = 0.1, 0.2 (stage 1's: tightest bend
    # 25.0 mm, 29.1 mm in the boards' plane framed k = 0, designs/helix_13.toml) and 0.3, each 317.8 mm on 3000 points,
    # the twin's torus path to the bit; the tighter the winding, the tighter its bend.
    curve = geometry.helix(ratio, L915, 3000)
    knot = ClosedCurve.torus_knot(1.0, ratio, 3000, 1, 3)
    assert np.array_equal(curve.points, knot.points * (L915 / knot.length))
    assert curve.length == L915
    assert curve.bend_radius / 1e-3 == pytest.approx(bend_mm, abs=1e-3)
    assert geometry.framing(curve, 0).inplane_radius / 1e-3 == pytest.approx(inplane_mm, abs=1e-3)


# ---- plates laid flat -----------------------------------------------------------------------------------------------
def test_the_flat_rings_plates_lie_as_straight_strips_or_as_rings() -> None:
    # Framed toward its boards, the ring's plates at +-h are circles on cylinders, geodesics there: laid flat they are
    # straight strips 2 pi (R +- h) long. Framed in the boards' plane, each is a planar circle of R: laid flat it closes
    # into a ring, turning 2 pi, in a box 2R + W square.
    curve = geometry.ring(L915, 4000)
    radius = L915 / (2 * math.pi)
    toward = geometry.flat(curve, "toward")
    for sign in (+1, -1):
        plate = geometry.plate(curve, toward.normal, sign * H915)
        assert plate.length == pytest.approx(2 * math.pi * (radius + sign * H915), rel=1e-6)
        assert abs(plate.turning) < 1e-9
        laid = geometry.develop(plate.dl, plate.kg)
        assert float(np.ptp(laid.xy[:, 1])) < 1e-9
    in_plane = geometry.flat(curve, "in-plane")
    plate = geometry.plate(curve, in_plane.normal, H915)
    assert abs(plate.turning) == pytest.approx(2 * math.pi, rel=1e-6)
    laid = geometry.develop(plate.dl, plate.kg)
    box = geometry.smallest_box(geometry.outline(laid.xy, laid.heading, WIDTH915))
    assert (box.short, box.long) == pytest.approx((2 * radius + WIDTH915, 2 * radius + WIDTH915), rel=1e-4)


def test_the_smallest_box_finds_a_turned_rectangle() -> None:
    # A 2 x 5 rectangle's outline turned by 0.3 rad: the box is 2 x 5 at the turn that undoes it, to the search's step
    # of pi/720, which leaves the sides long by at most 5 and 2 times half that step.
    side = np.linspace(0.0, 1.0, 50, endpoint=False)
    corners = np.array([[0.0, 0.0], [5.0, 0.0], [5.0, 2.0], [0.0, 2.0]])
    edge = np.concatenate([a + side[:, None] * (b - a) for a, b in zip(corners, np.roll(corners, -1, 0), strict=True)])
    c, s = math.cos(0.3), math.sin(0.3)
    turned = edge @ np.array([[c, s], [-s, c]])
    box = geometry.smallest_box(turned)
    half = math.pi / 1440
    assert 2.0 <= box.short <= 2.0 + 5 * half + 1e-9
    assert 5.0 <= box.long <= 5.0 + 2 * half + 1e-9
    assert abs(math.remainder(box.angle + 0.3, math.pi / 2)) <= half + 1e-12
