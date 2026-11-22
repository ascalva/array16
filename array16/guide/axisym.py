"""The axisymmetric bend reference: where a bend moves the facing boards' null, by quadrature (docs/KNOTTED_GUIDE.md
sections 4.4 and 8.1, E2-2 and DR-0).

Facing boards across a gap 2h, each split into halves at +V/2 and -V/2, opposite on the two boards (ELE-2, gapless and
unbounded), bent round a circle of radius R, are a body of revolution in two ways, and in each Laplace's equation
separates. Lengths are in h and volts in V throughout; rho is the distance from the bend's axis, z runs along it.

  in the boards' plane  annular plates at z = +-h, the top one at +V/2 inside rho = R and -V/2 outside it, the
                        bottom one opposite. Its inner disc is a Hankel transform: phi = -(V/2)(z/h) + V R
                        int_0^inf J1(kR) J0(k rho) sinh(kz)/sinh(kh) dk, and on the midplane, where phi vanishes,
                        E_z = V/2h - V R int J1(kR) J0(k rho) k/sinh(kh) dk.
  toward the boards     coaxial cylinders at rho = R -+ h, the inner one at +V/2 above z = 0 and -V/2 below it, the
                        outer one opposite. A Fourier transform along the axis: phi = (2/pi) int_0^inf sin(kz)/k
                        F(k, rho) dk, F the combination of I0(k rho) and K0(k rho) that is V/2 on the inner cylinder
                        and -V/2 on the outer, and on the plane z = 0, where phi vanishes, E_z = -(2/pi) int F dk.

The null is where E_z vanishes on that line, shifted toward the bend's axis from the circle the boards' split (in their
plane) or their midsurface (toward them) follows, and the gradient is |dE_z/drho| there: near the null the field is the
gradient times the distance from it, as on straight boards, whose gradient is pi V/4h^2. Both integrals decay as
exp(-kh) and are summed by Gauss-Legendre panels to k = 40/h, two panels to each period 2 pi/(R + rho) of the fast
part of the in-plane case's J1(kR) J0(k rho). The Bessel functions are this module's, in numpy: Bessel's integrals by
the trapezoidal rule below an argument of 30 and Hankel's asymptotic expansions above it (bessel_j01, bessel_ik01).

**Proposition (the shift to first order).** The null moves toward the bend's centre by 2 ln 2/pi^2 h^2/R = 0.140461
h^2/R in the boards' plane and by h^2/(pi R) = 0.318310 h^2/R toward them (E2-2's 0.1405 and 0.318; 7.7 and 17.4 um on
T915-S's tightest bend, 15.29 mm at h = 0.915 mm). With r = rho - R and a the section's other coordinate, along the
axis, and epsilon = h/R, Laplace's equation is phi_rr + phi_aa + epsilon phi_r/(1 + epsilon r) = 0 (h = 1), so
phi = phi0 + epsilon phi1, phi0 the straight boards' potential, nabla^2 phi1 = -d phi0/dr and phi1 = 0 on the boards.
Toward the boards r is the board normal: phi1's sine transform along the axis is (V/4)(r sinh(kr)/sinh k -
cosh(kr)/cosh k), its d/dz at the centre -(V/2pi) int dk/cosh k = -V/4 against phi0's gradient -pi V/4, so the null
sits at r = -epsilon/pi. In the boards' plane r runs across the split: phi1's cosine transform across it is
(V/(4k sinh k))(z cosh(kz) - cosh k sinh(kz)/sinh k), its d/dz at the centre (V/2pi) int (1 - k coth k)/(k sinh k) dk
= (V/2pi) int (1/sinh k - 1/k) dk/k = -(V/2pi) ln 2 (the Mellin transform 2(1 - 2^-s) Gamma(s) zeta(s) of 1/sinh at
s = 0), so the null sits at r = -2 ln 2 epsilon/pi^2. Both are linear in the curvature vector, so a bend at an angle psi
to the boards' plane moves the null by h^2/R (0.140461 cos psi, 0.318310 sin psi) across and toward the boards
(particles.sections.BendShift carries this into the swept field).

**Proposition (odd and even orders).** Bending the other way is the mirror image across the circle with the volts
reversed, which keeps the null's place and the gradient's size, so the shift is odd in 1/R and the gradient even:
shift x R/h^2 = c + c3 (h/R)^2 + ..., and the gradient over the straight boards' 1 + g2 (h/R)^2 + ...; the gradient
does not change at first order.

**Claim (the reference, by quadrature).** At R = 40h the shift times R/h^2 is 0.140473 in the boards' plane and
0.318364 toward them, 1.000085 and 1.000171 times the first order; c3 is 0.0191 and 0.0872, g2 -0.0838 and +0.0693, so
the gradient moves by -5.2e-5 and +4.3e-5 at R = 40h. The quadrature is converged to 1e-12 in shift x R/h^2 (panels
halved, the order raised or k taken further change it by no more), Richardson's extrapolation from R = 40h and 80h
gives the first order to 2e-8, and a finite-volume solve of the same two problems (sections.laplace2d with rho for
permittivity, Richardson-extrapolated from h/40 and h/80) meets the reference to 3e-5 at R = 10h
(tests/test_guide_axisym.py). E2-2's 0.318 is 1/pi to three figures: the reference at R = 40h lies 0.017 % above 1/pi
and 0.115 % above 0.318 read as written.

The three-dimensional boundary-element solve (array16/guide/bem.py, stage 2b) is calibrated against `bend` at R = 40h
both ways: its bias there, 2.5 to 4.5 % in the study, is DR-0's record.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np

from array16.particles.base import Arr
from array16.particles.sections import BEND_IN_PLANE, BEND_TOWARD

__all__ = [
    "E2_2_RADIUS",
    "FIRST_ORDER",
    "IN_PLANE",
    "TOWARD",
    "Bend",
    "Way",
    "bend",
    "bessel_ik01",
    "bessel_j01",
    "midline_field",
]

Way = Literal["in-plane", "toward"]
IN_PLANE = BEND_IN_PLANE  # 2 ln 2 / pi^2: the shift x R/h^2 to first order, a bend in the boards' plane
TOWARD = BEND_TOWARD  # 1 / pi: the same toward the boards
FIRST_ORDER: dict[Way, float] = {"in-plane": IN_PLANE, "toward": TOWARD}
E2_2_RADIUS = 40.0  # R/h at which E2-2 states the shifts

SWITCH = 30.0  # the argument from which the Bessel functions take their asymptotic expansions
_TERMS = 30  # the expansions' terms: their remainder at 30 is under 1e-20
_POINTS = 96  # the trapezoidal rule's points on Bessel's integrals below 30: the aliased J_95(30) is under 1e-20
_THETA = 2 * math.pi * np.arange(_POINTS) / _POINTS
_STEP = 0.1  # the trapezoidal rule's step on K's integral in cosh t: its error under exp(-pi^2 / 0.1) x exp(30)
_TAIL = 60.0  # where K's integrand is cut off: exp(-60)
_TOP = 40.0  # the transforms' reach in k h: both decay as exp(-k h)
_ORDER = 10  # Gauss-Legendre points a panel
_PANEL = 0.25  # the widest panel (1/h); the in-plane case takes two a period of its fast oscillation
_GRADED = 25  # panels growing geometrically from 1e-9/h to the first full one: the integrands' k^2 log k at 0


def _hankel(nu: int) -> Arr:
    """Hankel's coefficients a_k(nu) = (4 nu^2 - 1)(4 nu^2 - 9)...(4 nu^2 - (2k - 1)^2) / (k! 8^k), k < _TERMS."""
    out = np.ones(_TERMS)
    for k in range(1, _TERMS):
        out[k] = out[k - 1] * (4 * nu * nu - (2 * k - 1) ** 2) / (8 * k)
    return out


_A = (_hankel(0), _hankel(1))


def _series(coefficients: Arr, inverse: Arr) -> Arr:
    """sum_k coefficients[k] inverse^k, by Horner's rule."""
    out = np.zeros_like(inverse)
    for c in coefficients[::-1]:
        out = out * inverse + c
    return out


def bessel_j01(x: Arr) -> tuple[Arr, Arr]:
    """J0(x) and J1(x) for x >= 0: below 30 by the trapezoidal rule on Bessel's integral (1/2pi) int cos(n t - x sin t)
    dt over a period, 96 points, exact to rounding (its error is the aliased J_{96 - n}(x)); from 30 by Hankel's
    asymptotic expansion, sqrt(2/(pi x)) (P cos chi - Q sin chi) with chi = x - (n/2 + 1/4) pi, to 30 terms."""
    xs = np.asarray(x, dtype=np.float64)
    if np.any(xs < 0) or not np.all(np.isfinite(xs)):
        raise ValueError("bessel_j01 takes finite arguments x >= 0")
    j0, j1 = np.empty_like(xs), np.empty_like(xs)
    low = xs < SWITCH
    if low.any():
        phase = xs[low][:, None] * np.sin(_THETA)[None, :]
        j0[low] = np.cos(phase).mean(axis=1)
        j1[low] = np.cos(_THETA[None, :] - phase).mean(axis=1)
    high = ~low
    if high.any():
        y = xs[high]
        inverse = -1 / (y * y)
        for nu, out in ((0, j0), (1, j1)):
            a = _A[nu]
            p = _series(a[0::2], inverse)
            q = _series(a[1::2], inverse) / y
            chi = y - (nu / 2 + 0.25) * math.pi
            out[high] = np.sqrt(2 / (math.pi * y)) * (p * np.cos(chi) - q * np.sin(chi))
    return j0, j1


def bessel_ik01(x: Arr) -> tuple[Arr, Arr, Arr, Arr]:
    """The modified Bessel functions scaled so that none overflows, for x > 0: exp(-x) I0(x), exp(-x) I1(x),
    exp(x) K0(x) and exp(x) K1(x). Below 30 by the trapezoidal rule on (1/2pi) int exp(x (cos t - 1)) cos(n t) dt over a
    period (96 points) and on int_0^inf exp(-x (cosh t - 1)) cosh(n t) dt (steps of 0.1 to where the exponent reaches
    60), both exact to rounding; from 30 by their asymptotic expansions, (1/sqrt(2 pi x)) sum (-1)^k a_k/x^k and
    sqrt(pi/(2x)) sum a_k/x^k, to 30 terms."""
    xs = np.asarray(x, dtype=np.float64)
    if np.any(xs <= 0) or not np.all(np.isfinite(xs)):
        raise ValueError("bessel_ik01 takes finite arguments x > 0")
    i0, i1, k0, k1 = (np.empty_like(xs) for _ in range(4))
    low = xs < SWITCH
    if low.any():
        y = xs[low]
        bump = np.expm1(y[:, None] * (np.cos(_THETA)[None, :] - 1))  # I1's mean of cos t alone vanishes
        i0[low] = 1 + bump.mean(axis=1)
        i1[low] = (bump * np.cos(_THETA)[None, :]).mean(axis=1)
        t = np.arange(0.0, math.acosh(1 + _TAIL / float(y.min())) + _STEP, _STEP)
        weight = np.full(len(t), _STEP)
        weight[0] = _STEP / 2
        decay = np.exp(-y[:, None] * (np.cosh(t)[None, :] - 1)) * weight[None, :]
        k0[low] = decay.sum(axis=1)
        k1[low] = (decay * np.cosh(t)[None, :]).sum(axis=1)
    high = ~low
    if high.any():
        y = xs[high]
        inverse = 1 / y
        for nu, scaled_i, scaled_k in ((0, i0, k0), (1, i1, k1)):
            a = _A[nu]
            scaled_i[high] = _series(a, -inverse) / np.sqrt(2 * math.pi * y)
            scaled_k[high] = _series(a, inverse) * np.sqrt(math.pi / (2 * y))
    return i0, i1, k0, k1


# ---- the quadrature -------------------------------------------------------------------------------------------------
def _nodes(top: float, panel: float, order: int) -> tuple[Arr, Arr]:
    """Gauss-Legendre nodes and weights in k (1/h) on [0, top]: panels growing geometrically from 1e-9 to `panel`,
    then panels at most `panel` wide."""
    edges = np.concatenate(
        [
            [0.0],
            np.geomspace(1e-9, panel, _GRADED),
            np.linspace(panel, top, max(1, math.ceil((top - panel) / panel)) + 1)[1:],
        ]
    )
    x, w = np.polynomial.legendre.leggauss(order)
    lo, hi = edges[:-1, None], edges[1:, None]
    nodes: Arr = ((lo + hi) / 2 + (hi - lo) / 2 * x[None, :]).ravel()
    weights: Arr = ((hi - lo) / 2 * w[None, :]).ravel()
    return nodes, weights


def _in_plane(k: Arr, radius: float, rho: float, j1_r: Arr) -> tuple[Arr, Arr]:
    """The in-plane integrands per volt at rho: R J1(kR) J0(k rho) k/sinh k, whose integral less 1/2 is phi_z on the
    midplane, and its derivative in rho."""
    j0, j1 = bessel_j01(k * rho)
    kernel = radius * j1_r * k / np.sinh(k)
    return kernel * j0, -kernel * k * j1


def _toward(k: Arr, radius: float, rho: float) -> tuple[Arr, Arr]:
    """The toward-boards integrands per volt at rho: F(k, rho), the combination of I0 and K0 that is 1/2 on the inner
    cylinder (rho = R - 1) and -1/2 on the outer (R + 1), whose integral times 2/pi is phi_z on the plane z = 0, and
    dF/drho. Each product is scaled by exp(-k (b - a)) so that none overflows."""
    a, b = radius - 1, radius + 1
    i0_r, i1_r, k0_r, k1_r = bessel_ik01(k * rho)
    i0_a, _, k0_a, _ = bessel_ik01(k * a)
    i0_b, _, k0_b, _ = bessel_ik01(k * b)
    det = i0_a * k0_b * np.exp(-2 * k * (b - a)) - k0_a * i0_b  # I0(ka) K0(kb) - K0(ka) I0(kb)
    up, down = np.exp(k * (rho - 2 * b + a)), np.exp(k * (a - rho))  # the scaled exponentials of u: zero at b
    left, right = np.exp(k * (rho - b)), np.exp(k * (2 * a - rho - b))  # and of w: zero at a
    value = i0_r * k0_b * up - k0_r * i0_b * down + i0_r * k0_a * left - k0_r * i0_a * right
    slope = k * (i1_r * k0_b * up + k1_r * i0_b * down + i1_r * k0_a * left + k1_r * i0_a * right)
    return 0.5 * value / det, 0.5 * slope / det


@dataclass(frozen=True)
class Bend:
    """The facing boards bent round a circle of radius R (in h), in their plane or toward them, as the reference solves
    them: the null's shift toward the bend's axis (in h), the gradient there over the straight boards' pi V/4h^2, and
    the quadrature's nodes."""

    radius: float
    way: Way
    shift: float
    gradient: float
    nodes: int

    @property
    def coefficient(self) -> float:
        """The shift times R/h^2: E2-2's 0.1405 in the boards' plane and 0.318 toward them, at R = 40h."""
        return self.shift * self.radius

    @property
    def first_order(self) -> float:
        """The coefficient's first-order value, 2 ln 2/pi^2 or 1/pi (Proposition)."""
        return FIRST_ORDER[self.way]


def _check(radius: float, way: str, top: float, order: int, panel: float, resolve: float) -> None:
    if way not in FIRST_ORDER:
        raise ValueError(f"a bend is 'in-plane' or 'toward', not {way!r}")
    if not (math.isfinite(radius) and radius >= 2):
        raise ValueError(f"the reference takes a bend radius R >= 2h, not {radius}")
    if not (top > 0 and order >= 2 and panel > 0 and resolve > 0):
        raise ValueError("the quadrature takes a reach, an order of 2 or more, a panel width and a resolution > 0")


def _grid(way: Way, radius: float, top: float, order: int, panel: float, resolve: float) -> tuple[Arr, Arr]:
    width = panel if way == "toward" else min(panel, math.pi / (resolve * radius))
    return _nodes(top, width, order)


def bend(
    radius: float,
    way: Way,
    *,
    top: float = _TOP,
    order: int = _ORDER,
    panel: float = _PANEL,
    resolve: float = 2.0,
) -> Bend:
    """The facing boards bent round radius R = `radius` h in their plane ('in-plane') or toward them ('toward'): where
    the null is and the gradient there, by Newton's method on phi_z along the line through the null, each step's two
    integrals by Gauss-Legendre panels in k to `top`/h, `order` points a panel, panels at most `panel`/h wide and, in
    the boards' plane, `resolve` a period of J1(kR) J0(k rho)'s oscillation, 2 pi/(R + rho)."""
    _check(radius, way, top, order, panel, resolve)
    k, w = _grid(way, radius, top, order, panel, resolve)
    if way == "in-plane":
        _, j1_r = bessel_j01(k * radius)

        def parts(rho: float) -> tuple[Arr, Arr]:
            return _in_plane(k, radius, rho, j1_r)

        constant, scale = -0.5, 1.0
    else:

        def parts(rho: float) -> tuple[Arr, Arr]:
            return _toward(k, radius, rho)

        constant, scale = 0.0, 2 / math.pi
    rho = radius
    for _ in range(60):
        value, slope = parts(rho)
        step = -(constant + scale * float(w @ value)) / (scale * float(w @ slope))
        rho += step
        if abs(step) < 1e-15 * radius:
            break
    else:
        raise ArithmeticError(f"the null of the {way} bend at R = {radius}h did not settle")
    _, slope = parts(rho)
    gradient = abs(scale * float(w @ slope)) / (math.pi / 4)
    return Bend(radius, way, radius - rho, gradient, len(k))


def midline_field(radius: float, offsets: Arr, way: Way) -> Arr:
    """E_z per volt (in V/h) on the line through the null, at rho = R + offsets (in h): on the midplane in the boards'
    plane, where straight boards give (V/2h) tanh(pi x/2h) (ELE-2's E_y on y = 0), and on the plane z = 0 toward them,
    where they give (V/2h) tan(pi y/2h) (ELE-2's E_x on x = 0, |y| < h). Toward the boards the integrand decays as
    exp(-k d), d the distance to the nearer cylinder, so the transform is taken to k = 40/d there."""
    _check(radius, way, _TOP, _ORDER, _PANEL, 2.0)
    xs = np.atleast_1d(np.asarray(offsets, dtype=np.float64))
    out = np.empty(len(xs))
    if way == "in-plane":
        k, w = _grid(way, radius, _TOP, _ORDER, _PANEL, 2.0)
        _, j1_r = bessel_j01(k * radius)
        for i, x in enumerate(xs):
            value, _ = _in_plane(k, radius, radius + float(x), j1_r)
            out[i] = 0.5 - float(w @ value)
        return out
    if not bool(np.all(np.abs(xs) < 1)):
        raise ValueError("toward the boards the line through the null runs between them: |offset| < h")
    for i, x in enumerate(xs):
        k, w = _grid(way, radius, _TOP / (1 - abs(float(x))), _ORDER, _PANEL, 2.0)
        value, _ = _toward(k, radius, radius + float(x))
        out[i] = -2 / math.pi * float(w @ value)
    return out
