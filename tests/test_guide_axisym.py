"""The axisymmetric bend reference (array16/guide/axisym.py), held to what it must reproduce: its Bessel functions to
their power series in 70-digit arithmetic and to Bessel's and Basset's integrals; the first-order shifts' closed forms,
2 ln 2/pi^2 and 1/pi, as Propositions with their integrals; E2-2's 0.1405 and 0.318 at R = 40h, within 0.1 % of the
first order (DR-0); the shift odd in the curvature and the gradient even; the quadrature converged, and a quadrature too
coarse for the in-plane oscillation, the case that must fail, missing by more than 0.1 %; a finite-volume solve of the
same two bodies of revolution meeting the reference; and the reference reducing to ELE-2's straight boards."""

from __future__ import annotations

import math
from collections.abc import Callable
from decimal import Decimal, localcontext

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from array16.guide import axisym
from array16.guide.sections import laplace2d
from array16.particles import FacingBoards
from array16.particles.base import Arr
from array16.particles.sections import BEND_IN_PLANE, BEND_TOWARD

GAMMA = Decimal("0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917")


# ---- the Bessel functions -------------------------------------------------------------------------------------------
def _series(x: float) -> dict[str, float]:
    """J0, J1, exp(-x) I0, exp(-x) I1, exp(x) K0 and exp(x) K1 at x by their power series in 70-digit decimal
    arithmetic (Abramowitz and Stegun 9.1.10, 9.6.10, 9.6.11): K's series cancel I's growth, which 70 digits keep to
    35 or more up to x = 32."""
    with localcontext() as ctx:
        ctx.prec = 70
        big = Decimal(repr(x))
        quarter = big * big / 4
        j0 = j1 = i0 = i1 = k0_tail = k1_tail = Decimal(0)
        term, harmonic, k = Decimal(1), Decimal(0), 0  # term = (x^2/4)^k / (k!)^2, harmonic = H_k
        while True:
            odd = term * big / 2 / (k + 1)  # (x/2)^(2k + 1) / (k! (k + 1)!)
            sign = 1 if k % 2 == 0 else -1
            j0, i0 = j0 + sign * term, i0 + term
            j1, i1 = j1 + sign * odd, i1 + odd
            k0_tail += term * harmonic
            k1_tail += odd * (2 * harmonic + Decimal(1) / (k + 1) - 2 * GAMMA)  # psi(k + 1) + psi(k + 2)
            k += 1
            harmonic += Decimal(1) / k
            term = term * quarter / (k * k)
            if k > 10 and term < Decimal(10) ** -75:
                break
        log = (big / 2).ln()
        k0 = -(log + GAMMA) * i0 + k0_tail
        k1 = 1 / big + log * i1 - k1_tail / 2
        grow = big.exp()
        return {
            "j0": float(j0),
            "j1": float(j1),
            "i0": float(i0 / grow),
            "i1": float(i1 / grow),
            "k0": float(k0 * grow),
            "k1": float(k1 * grow),
        }


def test_the_bessel_functions_meet_their_power_series() -> None:
    # Below the switch at 30 the trapezoidal rule on Bessel's integrals, and just above it the asymptotic expansions,
    # give J to 1e-15 and the scaled I and K to 2e-15 of their size, against the power series in 70 digits.
    points = [1e-6, 0.01, 0.5, 1.0, 2.404825557695773, 5.0, 12.0, 25.0, 29.999, 30.0, 30.5, 32.0]
    j0, j1 = axisym.bessel_j01(np.array(points))
    i0, i1, k0, k1 = axisym.bessel_ik01(np.array(points))
    for n, x in enumerate(points):
        exact = _series(x)
        assert j0[n] == pytest.approx(exact["j0"], abs=1e-15)
        assert j1[n] == pytest.approx(exact["j1"], abs=1e-15)
        for got, name in ((i0[n], "i0"), (i1[n], "i1"), (k0[n], "k0"), (k1[n], "k1")):
            assert got == pytest.approx(exact[name], rel=2e-15)
    assert abs(float(axisym.bessel_j01(np.array([2.404825557695773]))[0][0])) < 1e-15  # J0's first zero


def _bessel_integral(x: float, n: int) -> float:
    """J_n(x) by the trapezoidal rule on (1/2pi) int cos(n t - x sin t) dt with 2x + 100 points: exact to rounding,
    the aliased J_(2x + 100)(x) being under exp(-0.9x)."""
    count = 2 * int(x) + 100
    t = 2 * math.pi * np.arange(count) / count
    return float(np.cos(n * t - x * np.sin(t)).mean())


def _basset(x: float, n: int) -> float:
    """exp(x) K_n(x) by the trapezoidal rule on int_0^inf exp(-x (cosh t - 1)) cosh(n t) dt, steps of 0.01."""
    t = np.arange(0.0, math.acosh(1 + 80 / x), 0.01)
    f = np.exp(-x * (np.cosh(t) - 1)) * np.cosh(n * t)
    return float(0.01 * (f.sum() - f[0] / 2))


@pytest.mark.parametrize("x", [30.0, 47.3, 100.0, 400.0, 1600.0])
def test_the_asymptotic_expansions_meet_the_integrals(x: float) -> None:
    # From 30 up, where the module sums Hankel's expansions, they meet Bessel's integral for J by the trapezoidal rule
    # on 2x + 100 points, to the rounding of its phases x sin t (x times 1e-16 at most), Basset's integral for K,
    # numpy's i0 for I0, and the Wronskian I0 K1 + I1 K0 = 1/x.
    j0, j1 = axisym.bessel_j01(np.array([x]))
    rounding = 2e-15 * (1 + x / 100)
    assert float(j0[0]) == pytest.approx(_bessel_integral(x, 0), abs=rounding)
    assert float(j1[0]) == pytest.approx(_bessel_integral(x, 1), abs=rounding)
    i0, i1, k0, k1 = (float(f[0]) for f in axisym.bessel_ik01(np.array([x])))
    assert k0 == pytest.approx(_basset(x, 0), rel=1e-13)
    assert k1 == pytest.approx(_basset(x, 1), rel=1e-13)
    if x < 700:
        assert i0 == pytest.approx(float(np.i0(x)) * math.exp(-x), rel=1e-14)
    assert x * (i0 * k1 + i1 * k0) == pytest.approx(1.0, rel=1e-14)


@settings(max_examples=60, deadline=None)
@given(st.floats(1e-8, 1e4))
def test_the_wronskian_holds_everywhere(x: float) -> None:
    # I0 K1 + I1 K0 = 1/x for the scaled functions, whose exponentials cancel in the products, at any argument.
    i0, i1, k0, k1 = (float(f[0]) for f in axisym.bessel_ik01(np.array([x])))
    assert x * (i0 * k1 + i1 * k0) == pytest.approx(1.0, rel=1e-13)


def test_the_bessel_functions_refuse_what_they_do_not_take() -> None:
    with pytest.raises(ValueError, match="x >= 0"):
        axisym.bessel_j01(np.array([-1.0]))
    with pytest.raises(ValueError, match="x > 0"):
        axisym.bessel_ik01(np.array([0.0]))


# ---- the first order ------------------------------------------------------------------------------------------------
def _legendre(f: Callable[[Arr], Arr], top: float, panels: int = 400) -> float:
    """int_0^top f(k) dk by Gauss-Legendre on equal panels, 16 points each."""
    x, w = np.polynomial.legendre.leggauss(16)
    edges = np.linspace(0.0, top, panels + 1)
    k = ((edges[:-1, None] + edges[1:, None]) / 2 + (edges[1] - edges[0]) / 2 * x[None, :]).ravel()
    return float(np.sum(np.tile(w, panels) * f(k)) * (edges[1] - edges[0]) / 2)


def test_the_first_order_shifts_are_two_ln_two_over_pi_squared_and_one_over_pi() -> None:
    # Propositions: the integrals the first order rests on, (1/sinh k - 1/k)/k to -ln 2 (its tail beyond 60, -1/60 to
    # 1e-26, added) and 1/cosh k to pi/2, and so the coefficients 2 ln 2/pi^2 and 1/pi; the particles' BendShift
    # carries the same numbers. The reference's Richardson extrapolation from R = 40h and 80h, its (h/R)^2 term taken
    # out, gives them to 2e-8.
    def mellin(k: Arr) -> Arr:
        out: Arr = (1 / np.sinh(k) - 1 / k) / k
        return out

    assert _legendre(mellin, 60.0) - 1 / 60 == pytest.approx(-math.log(2), abs=1e-12)
    assert _legendre(lambda k: 1 / np.cosh(k), 60.0) == pytest.approx(math.pi / 2, abs=1e-12)
    assert axisym.IN_PLANE == BEND_IN_PLANE == pytest.approx(2 * math.log(2) / math.pi**2, rel=1e-15)
    assert axisym.TOWARD == BEND_TOWARD == pytest.approx(1 / math.pi, rel=1e-15)
    assert abs(axisym.IN_PLANE - 0.140461) < 5e-7
    assert abs(axisym.TOWARD - 0.318310) < 5e-7
    for way in ("in-plane", "toward"):
        near, far = axisym.bend(40.0, way), axisym.bend(80.0, way)
        assert (4 * far.coefficient - near.coefficient) / 3 == pytest.approx(near.first_order, rel=2e-8)


# ---- E2-2 at R = 40h ------------------------------------------------------------------------------------------------
@pytest.mark.parametrize(("way", "stated", "value"), [("in-plane", 0.1405, 0.1404729), ("toward", 0.318, 0.3183644)])
def test_e2_2_at_forty_h_the_shifts_are_the_first_order_to_a_tenth_of_a_percent(
    way: axisym.Way, stated: float, value: float
) -> None:
    # DR-0: at R = 40h the shift times R/h^2 is the first-order coefficient to 0.1 % (to 0.0085 % in the boards' plane
    # and 0.017 % toward them), and E2-2's figures to the three digits they are stated to: 0.1405 and 0.318. 0.318 is
    # 1/pi rounded, and the reference lies 0.115 % above it as written; E2-2's own band is 2 %.
    found = axisym.bend(axisym.E2_2_RADIUS, way)
    assert found.coefficient == pytest.approx(value, abs=1e-7)
    assert abs(found.coefficient / found.first_order - 1) < 1e-3
    assert round(found.coefficient, len(str(stated)) - 2) == stated
    assert abs(found.coefficient / stated - 1) < 0.02
    if way == "toward":
        assert 1.1e-3 < found.coefficient / stated - 1 < 1.2e-3
    # on T915-S's tightest bend, 15.29 mm at h = 0.915 mm: 7.7 and 17.4 um
    h, radius = 0.915e-3, 15.29e-3
    assert found.first_order * h * h / radius == pytest.approx(7.7e-6 if way == "in-plane" else 17.4e-6, abs=0.05e-6)


@pytest.mark.parametrize("way", ["in-plane", "toward"])
def test_the_shift_is_odd_in_the_curvature_and_the_gradient_even(way: axisym.Way) -> None:
    # Proposition: shift x R/h^2 = c + c3 (h/R)^2 and gradient = 1 + g2 (h/R)^2, no term in h/R: from R = 10h to 80h
    # (S - c) R^2 and (G - 1) R^2 hold to 1 %, c3 0.0191 and 0.0872, g2 -0.0838 and +0.0693, and the gradient moves by
    # under 6e-5 at R = 40h (by under 3.1e-4 on T915-S's tightest bend, h/R = 0.0598).
    c3, g2 = {"in-plane": (0.0191, -0.0838), "toward": (0.0872, 0.0693)}[way]
    for radius in (10.0, 20.0, 40.0, 80.0):
        found = axisym.bend(radius, way)
        assert (found.coefficient - found.first_order) * radius**2 == pytest.approx(c3, rel=1e-2)
        assert (found.gradient - 1) * radius**2 == pytest.approx(g2, rel=1e-2)
    assert abs(axisym.bend(40.0, way).gradient - 1) < 6e-5
    assert abs(g2) * (0.915 / 15.29) ** 2 < 3.1e-4


# ---- the quadrature -------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("way", ["in-plane", "toward"])
def test_the_quadrature_is_converged(way: axisym.Way) -> None:
    # k taken to 50/h, 14 points a panel, panels halved, or four a period of the in-plane oscillation: the shift x R/h^2
    # and the gradient move by under 1e-11.
    base = axisym.bend(40.0, way)
    for other in (
        axisym.bend(40.0, way, top=50.0),
        axisym.bend(40.0, way, order=14),
        axisym.bend(40.0, way, panel=0.125),
        axisym.bend(40.0, way, resolve=4.0),
    ):
        assert other.coefficient == pytest.approx(base.coefficient, abs=1e-11)
        assert other.gradient == pytest.approx(base.gradient, abs=1e-11)


def test_a_quadrature_too_coarse_for_the_in_plane_oscillation_must_fail() -> None:
    # The case that must fail: J1(kR) J0(k rho) runs through (R + rho)/2pi periods per unit of k h, and with eight
    # points on panels four periods wide the in-plane reference misses the first order by more than 0.1 % (by -19 %);
    # two panels a period with eight points still meet it.
    coarse = axisym.bend(40.0, "in-plane", resolve=0.25, order=8)
    assert abs(coarse.coefficient / coarse.first_order - 1) > 1e-3
    fine = axisym.bend(40.0, "in-plane", resolve=2.0, order=8)
    assert abs(fine.coefficient / fine.first_order - 1) < 1e-3


def test_the_reference_refuses_what_it_does_not_solve() -> None:
    with pytest.raises(ValueError, match="R >= 2h"):
        axisym.bend(1.5, "toward")
    with pytest.raises(ValueError, match="'in-plane' or 'toward'"):
        axisym.bend(40.0, "across")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="between them"):
        axisym.midline_field(40.0, np.array([1.0]), "toward")


# ---- against a finite-volume solve ----------------------------------------------------------------------------------
def _finite_volume(radius: float, way: axisym.Way, step: float, span: float = 6.0) -> float:
    """The shift x R/h^2 of the same body of revolution solved by finite volumes (h = 1, V = 1): div(rho grad phi) = 0
    on a uniform (rho, z) grid of `step`, rho for permittivity, the boards fixed and the box's other faces left free
    `span` from the split; the null where the central difference of phi across the line through it changes sign."""
    if way == "in-plane":
        r = np.linspace(radius - span, radius + span, round(2 * span / step) + 1)
        z = np.linspace(-1.0, 1.0, round(2 / step) + 1)
    else:
        r = np.linspace(radius - 1, radius + 1, round(2 / step) + 1)
        z = np.linspace(-span, span, round(2 * span / step) + 1)
    eps = np.repeat(((r[1:] + r[:-1]) / 2)[:, None], len(z) - 1, axis=1)
    fixed = np.zeros((len(r), len(z)), dtype=np.bool_)
    value = np.zeros((len(r), len(z)))
    if way == "in-plane":  # the top plate +1/2 inside the split and -1/2 outside it, the bottom one opposite
        top = 0.5 * np.sign(radius - r) * (np.abs(r - radius) > step / 4)
        fixed[:, 0] = fixed[:, -1] = True
        value[:, -1], value[:, 0] = top, -top
    else:  # the inner cylinder +1/2 above z = 0 and -1/2 below it, the outer one opposite
        fixed[0, :] = fixed[-1, :] = True
        value[0, :], value[-1, :] = 0.5 * np.sign(z), -0.5 * np.sign(z)
    phi = laplace2d(r, z, eps, fixed, value)
    mid = len(z) // 2
    across = (phi[:, mid + 1] - phi[:, mid - 1]) / (z[mid + 1] - z[mid - 1])
    turns = np.flatnonzero(np.sign(across[:-1]) != np.sign(across[1:]))
    i = int(turns[np.argmin(np.abs(r[turns] - radius))])
    null = r[i] - across[i] * (r[i + 1] - r[i]) / (across[i + 1] - across[i])
    return float((radius - null) * radius)


@pytest.mark.parametrize("way", ["in-plane", "toward"])
def test_a_finite_volume_solve_of_the_same_body_of_revolution_meets_the_reference(way: axisym.Way) -> None:
    # At R = 10h, where the shift is 0.014h and 0.032h: the finite volumes converge as the step squared, and
    # Richardson's extrapolation from h/40 and h/80 meets the quadrature to 3e-5 of the shift, where the grid of h/80
    # alone is 1e-4 to 4e-4 off it.
    coarse, fine = _finite_volume(10.0, way, 1 / 40), _finite_volume(10.0, way, 1 / 80)
    reference = axisym.bend(10.0, way).coefficient
    assert (fine + (fine - coarse) / 3) == pytest.approx(reference, rel=3e-5)
    assert 5e-5 < abs(fine / reference - 1) < 1e-3


# ---- the straight limit ---------------------------------------------------------------------------------------------
@pytest.mark.parametrize("way", ["in-plane", "toward"])
def test_the_reference_reduces_to_the_straight_boards(way: axisym.Way) -> None:
    # Bent round 40h and 400h, the field across the line through the null is ELE-2's straight boards' plus a part in
    # h/R: (E - E_straight) R is the same profile at both radii to its own next order, under 6 % of its size, and at the
    # circle it is the uniform field the bend adds there, ln 2/(2 pi) V/h in the boards' plane and V/4h toward them,
    # which the shift cancels against the gradient pi V/4h^2: 2 ln 2/pi^2 and 1/pi again.
    offsets = np.linspace(-0.9, 0.9, 7)
    boards = FacingBoards(2.0)
    if way == "in-plane":
        straight = np.array([boards.field(x, 0.0, 1.0)[1] for x in offsets])
    else:
        straight = np.array([boards.field(0.0, y, 1.0)[0] for y in offsets])
    profiles = [(axisym.midline_field(radius, offsets, way) - straight) * radius for radius in (40.0, 400.0)]
    assert np.allclose(profiles[0], profiles[1], atol=1e-2)
    assert float(np.abs(profiles[1] - profiles[0]).max()) < float(np.abs(profiles[1]).max()) * 0.06
    centre = float(profiles[1][3])
    assert centre == pytest.approx(math.log(2) / (2 * math.pi) if way == "in-plane" else 0.25, rel=1e-3)
    assert centre / (math.pi / 4) == pytest.approx(axisym.FIRST_ORDER[way], rel=1e-3)
