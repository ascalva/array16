"""The figures of docs/CAVITY_KNOTS.md: a cavity's modes, the null curves of its real standing fields, and what a
well round them is worth.

The document answers the owner's question of a spherical containment cavity with knotted fields: which mathematics
describes it (Bessel functions, symmetry, Morse theory, knot theory), where a real TE field vanishes, and a worked
copper cylinder at 2.45 GHz whose field holds three closed null loops. This script computes its figures (section of
the document; figure number):

  models     the owner's sphere with its feeds, pickups and one real field's null diameters; the worked cylinder
             with its feeds, loops, open null lines, folds and escape saddles (2, 8.2; figure 1)
  nullsets   one sphere family's diameters and an n = 2 family's inner zero-field shell; the cylinder's loops on the
             vertical surfaces over the curve C; a random field of the cube's (1, 2, 3) family (3 to 5; 2)
  charts     the sphere's TE and TM roots against l; the cylinder's mode chart and its TE113-TE221 crossing (3; 3)
  folds      a loop coloured by det S, the wall's growth at a fold, the normal gradients along the loop, and a Hopf
             link and a trefoil against height (6; 4)
  sections   log |E|^2 of the worked field on the mid-plane, the plane y = 0 and the meridian of its escape saddle,
             with the escape path (7; 5)
  power      depth against power at 2.45 GHz: the loops, RFG-2's ceiling, facing boards and four rods, with the
             thresholds of section 8.5 (8.1, 8.5; 6)
  tolerance  the tube's floor against the split of the two modes in linewidths (8.3; 7)

    uv run python experiments/cavity_knots_study.py              # every figure, about half a minute
    uv run python experiments/cavity_knots_study.py --quick      # coarser grids and steps: cavity_knots/quick/
    uv run python experiments/cavity_knots_study.py --only charts,power

Model: perfectly conducting walls for the fields, copper's surface resistance (skin depth 1.34 um, as the document
takes it) for the Q's; lengths in units of the cavity's radius a (or the cube's side) inside the computations. A
real TE field is E = grad(Pi) x z_hat in a cylinder or cube and E = grad(psi) x r in a sphere. Bessel functions are
computed here, typed, from their integral (integer order) and from series and recurrence (spherical), and checked
against their known zeros in tests/test_cavity_knots_study.py. Null curves are traced by continuation along
grad Pi_x x grad Pi_y with Newton's method across the curve; wells by a minimax flood over a grid and Newton's method
on grad |E|^2 = 0. Outputs: cavity_knots/{*.png, results.json, run.log}.
"""

from __future__ import annotations

import argparse
import heapq
import json
import math
import sys
import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from array16 import lab
from array16 import particles as pt
from array16.particles import Arr, CArr, IntArr

OUT = Path(__file__).resolve().parent / "cavity_knots"
LOG: list[str] = []
F_DRIVE = 2.45e9
OMEGA = 2 * math.pi * F_DRIVE
K_DRIVE = OMEGA / pt.C0  # 51.348 per metre
SKIN = 1.34e-6  # copper's skin depth at 2.45 GHz, as the document takes it (sigma = 5.8e7 S/m gives 1.335 um)
QM = pt.E / pt.ME

# Ink, surfaces and the categorical order of the reference palette (fixed order, never cycled).
INK, INK2, MUTED, GRID, EDGE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
BLUES = ("#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
SURFACE, COPPER = "#fcfcfb", "#c8874a"


def say(line: str = "") -> None:
    print(line, flush=True)
    LOG.append(line)


# ---- Bessel functions ----------------------------------------------------------------------------------------------
NODES = 128
_ANGLES = 2 * math.pi * np.arange(NODES) / NODES
_SINES = np.sin(_ANGLES)
CHUNK = 8192


def bessel_ladder(orders: Sequence[int], x: Any) -> Arr:
    """J_n(x) for every integer order n in `orders` and real x, shape (len(orders),) + shape of x.

    The trapezoid rule on N equally spaced nodes of the periodic integral J_n(x) = (1/2 pi) int_0^2pi cos(n t -
    x sin t) dt gives exactly the sum of J_{n + kN}(x) over all integers k, so with N = 128 nodes it is J_n(x) to
    rounding for |x| up to 40 and |n| up to 20, where J_{n - 128}(x) is below 1e-30. One pass of cos and sin of
    x sin t serves every order: cos(n t - x sin t) = cos(n t) cos(x sin t) + sin(n t) sin(x sin t)."""
    xs = np.asarray(x, dtype=np.float64)
    flat = xs.reshape(-1)
    n = np.asarray(orders, dtype=np.float64)
    cos_n, sin_n = np.cos(np.outer(n, _ANGLES)), np.sin(np.outer(n, _ANGLES))
    out = np.empty((len(n), flat.size))
    for start in range(0, flat.size, CHUNK):
        arg = np.outer(flat[start : start + CHUNK], _SINES)
        out[:, start : start + CHUNK] = (cos_n @ np.cos(arg).T + sin_n @ np.sin(arg).T) / NODES
    result: Arr = out.reshape((len(n), *xs.shape))
    return result


def bessel_j(order: int, x: Any) -> Arr:
    """J_n(x) for an integer order n."""
    out: Arr = bessel_ladder([order], x)[0]
    return out


def bessel_jp(order: int, x: Any) -> Arr:
    """J_n'(x) = (J_{n-1}(x) - J_{n+1}(x)) / 2."""
    pair = bessel_ladder([order - 1, order + 1], x)
    out: Arr = 0.5 * (pair[0] - pair[1])
    return out


def spherical_j(order: int, x: Any) -> Arr:
    """The spherical Bessel function j_l(x), l >= 0: its power series x^l sum (-x^2/2)^k / (k! (2l + 2k + 1)!!)
    where x <= l + 1, and above that the upward recurrence j_{n+1} = (2n + 1) j_n / x - j_{n-1} from
    j_0 = sin x / x and j_1 = sin x / x^2 - cos x / x, stable while the order stays below x."""
    xs = np.asarray(x, dtype=np.float64)
    out = np.empty_like(xs)
    small = np.abs(xs) <= order + 1
    if small.any():
        z = xs[small]
        term = z**order / float(np.prod(np.arange(1, 2 * order + 2, 2, dtype=np.float64)))
        total = term.copy()
        for k in range(80):
            term = term * (-(z**2) / 2) / ((k + 1) * (2 * order + 2 * k + 3))
            total = total + term
        out[small] = total
    if (~small).any():
        z = xs[~small]
        previous, current = np.sin(z) / z, np.sin(z) / z**2 - np.cos(z) / z
        if order == 0:
            current = previous
        for n in range(1, order):
            previous, current = current, (2 * n + 1) * current / z - previous
        out[~small] = current
    return out


def spherical_jp(order: int, x: Any) -> Arr:
    """j_l'(x) = j_{l-1}(x) - (l + 1) j_l(x) / x, and j_0' = -j_1."""
    xs = np.asarray(x, dtype=np.float64)
    if order == 0:
        return -spherical_j(1, xs)
    out: Arr = spherical_j(order - 1, xs) - (order + 1) * spherical_j(order, xs) / xs
    return out


def roots(fn: Callable[[Arr], Arr], lo: float, hi: float, step: float = 0.01) -> Arr:
    """The roots of fn in (lo, hi): sign changes on a grid of the given step, each refined by 64 bisections."""
    xs = np.arange(lo, hi, step)
    values = fn(xs)
    at = np.flatnonzero(np.sign(values[:-1]) * np.sign(values[1:]) < 0)
    left, right = xs[at].copy(), xs[at + 1].copy()
    f_left = values[at].copy()
    for _ in range(64):
        middle = 0.5 * (left + right)
        f_mid = fn(middle)
        same = np.sign(f_mid) == np.sign(f_left)
        left, f_left = np.where(same, middle, left), np.where(same, f_mid, f_left)
        right = np.where(same, right, middle)
    found: Arr = 0.5 * (left + right)
    return found


def bessel_zeros(order: int, count: int, derivative: bool = False) -> Arr:
    """The first `count` positive zeros of J_m (x_mn, TM modes of a cylinder) or of J_m' (x'_mn, TE modes); the zero
    of J_0' at the origin is left out, as is the convention."""
    top = 4.0 + math.pi * (count + order / 2 + 2)
    fn: Callable[[Arr], Arr] = (lambda z: bessel_jp(order, z)) if derivative else (lambda z: bessel_j(order, z))
    found = roots(fn, 0.5, top)
    return found[:count]


def sphere_zeros(order: int, count: int, tm: bool = False) -> Arr:
    """ka of the first `count` TE families of degree l in a sphere, j_l(ka) = 0, or of the TM families,
    [x j_l(x)]' = x j_{l-1}(x) - l j_l(x) = 0."""
    top = 4.0 + math.pi * (count + order / 2 + 2)

    def tm_condition(z: Arr) -> Arr:
        out: Arr = z * spherical_j(order - 1, z) - order * spherical_j(order, z)
        return out

    fn: Callable[[Arr], Arr] = tm_condition if tm else (lambda z: spherical_j(order, z))
    found = roots(fn, 0.5, top)
    return found[:count]


# ---- harmonic polynomials: the angular parts of a sphere's modes ---------------------------------------------------
@dataclass(frozen=True)
class Harmonics:
    """An orthonormal basis of the real harmonic polynomials of degree l in (x, y, z), 2l + 1 of them: on the unit
    sphere they are the real spherical harmonics of degree l up to a rotation of the basis, with the mean of
    h_i h_j over the sphere equal to delta_ij. exponents (M, 3) lists the monomials, coefficients (2l + 1, M)."""

    degree: int
    exponents: IntArr
    coefficients: Arr

    @staticmethod
    def of(degree: int) -> Harmonics:
        def monomials(n: int) -> list[tuple[int, int, int]]:
            return [(a, b, n - a - b) for a in range(n, -1, -1) for b in range(n - a, -1, -1)]

        top, low = monomials(degree), monomials(degree - 2) if degree >= 2 else []
        where = {m: i for i, m in enumerate(low)}
        laplacian = np.zeros((max(len(low), 1), len(top)))
        for j, (a, b, c) in enumerate(top):
            for axis, power in enumerate((a, b, c)):
                if power >= 2:
                    lowered = [a, b, c]
                    lowered[axis] -= 2
                    laplacian[where[(lowered[0], lowered[1], lowered[2])], j] += power * (power - 1)
        _, _, vt = np.linalg.svd(laplacian)
        null = vt[len(low) :] if degree >= 2 else np.eye(len(top))
        exponents = np.array(top, dtype=np.intp)
        points, weights = sphere_quadrature(2 * degree + 2)
        values = null @ monomial_values(exponents, points)
        gram = (values * weights) @ values.T
        lower = np.linalg.cholesky(gram)
        coefficients = np.linalg.solve(lower, null)
        return Harmonics(degree, exponents, coefficients)

    def gradients(self, points: Arr) -> tuple[Arr, Arr]:
        """Every basis polynomial's value (B, P) and gradient (B, P, 3) at points (P, 3)."""
        e = self.exponents
        value = self.coefficients @ monomial_values(e, points)
        grad = np.zeros((len(self.coefficients), len(np.atleast_2d(points)), 3))
        for i in range(3):
            step = np.zeros(3, dtype=np.intp)
            step[i] = 1
            grad[:, :, i] = (self.coefficients * e[:, i]) @ monomial_values(e - step, points)
        return value, grad

    def values(self, weights: Arr, points: Arr) -> tuple[Arr, Arr, Arr]:
        """The polynomial A = sum_i w_i h_i at points (P, 3): its value (P,), gradient (P, 3) and Hessian (P, 3, 3)."""
        c = weights @ self.coefficients
        e = self.exponents
        value = c @ monomial_values(e, points)
        grad = np.zeros((len(points), 3))
        hess = np.zeros((len(points), 3, 3))
        for i in range(3):
            step = np.zeros(3, dtype=np.intp)
            step[i] = 1
            grad[:, i] = (c * e[:, i]) @ monomial_values(e - step, points)
            for j in range(3):
                step2 = step.copy()
                step2[j] += 1
                factor = e[:, i] * (e[:, j] - (1 if i == j else 0))
                hess[:, i, j] = (c * factor) @ monomial_values(e - step2, points)
        return value, grad, hess


def monomial_values(exponents: IntArr, points: Arr) -> Arr:
    """x^a y^b z^c for each row (a, b, c) of exponents at each point, shape (M, P); negative powers give zero (their
    coefficient is zero wherever they arise)."""
    p = np.atleast_2d(points)
    out = np.ones((len(exponents), len(p)))
    for axis in range(3):
        powers = exponents[:, axis]
        safe = np.maximum(powers, 0)
        out *= np.where(powers[:, None] >= 0, p[None, :, axis] ** safe[:, None], 0.0)
    return out


def sphere_quadrature(order: int) -> tuple[Arr, Arr]:
    """Points on the unit sphere and weights (summing to 1) exact for polynomials up to the given degree:
    Gauss-Legendre in cos theta times the trapezoid rule in phi."""
    nodes, w = np.polynomial.legendre.leggauss(order // 2 + 2)
    phis = 2 * math.pi * np.arange(order + 2) / (order + 2)
    ct, ph = np.meshgrid(nodes, phis, indexing="ij")
    st = np.sqrt(1 - ct**2)
    points = np.c_[(st * np.cos(ph)).ravel(), (st * np.sin(ph)).ravel(), ct.ravel()]
    weights = (w[:, None] * np.ones(len(phis))[None, :]).ravel() / (2 * len(phis))
    return points, weights


def fibonacci_sphere(n: int) -> Arr:
    """n points spread evenly over the unit sphere."""
    i = np.arange(n) + 0.5
    ct = 1 - 2 * i / n
    phi = math.pi * (1 + math.sqrt(5)) * i
    st = np.sqrt(1 - ct**2)
    points: Arr = np.c_[st * np.cos(phi), st * np.sin(phi), ct]
    return points


@dataclass(frozen=True)
class CriticalPoint:
    """A critical point of an angular function A on the unit sphere: where, A there, and its kind."""

    where: Arr
    value: float
    kind: str  # "maximum", "minimum" or "saddle"


def tangent_basis(x: Arr) -> tuple[Arr, Arr]:
    """Two unit vectors spanning the tangent plane of the unit sphere at x."""
    helper = np.array([0.0, 0.0, 1.0]) if abs(x[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(x, helper)
    e1 /= np.linalg.norm(e1)
    return e1, np.cross(x, e1)


def critical_points(basis: Harmonics, weights: Arr, samples: int = 6000) -> list[CriticalPoint]:
    """Every critical point of A = sum w_i h_i on the unit sphere: candidates where |grad_S A| is smallest among
    their neighbours on a Fibonacci grid, each refined by Newton's method in its tangent plane with the surface
    Hessian P Hess(A) P - (x . grad A) P, deduplicated, and classed by that Hessian's signs."""
    grid = fibonacci_sphere(samples)
    value, grad, _ = basis.values(weights, grid)
    surface = grad - (grad * grid).sum(1)[:, None] * grid
    size = np.linalg.norm(surface, axis=1)
    spacing = math.sqrt(4 * math.pi / samples)
    found: list[CriticalPoint] = []
    order = np.argsort(size)
    for index in order[: samples // 8]:
        near = np.linalg.norm(grid - grid[index], axis=1) < 2.5 * spacing
        if size[index] > size[near].min():
            continue
        x = grid[index].copy()
        for _ in range(40):
            _, g, h = basis.values(weights, x[None, :])
            e1, e2 = tangent_basis(x)
            frame = np.c_[e1, e2]
            gs = frame.T @ g[0]
            hs = frame.T @ h[0] @ frame - float(x @ g[0]) * np.eye(2)
            step = -np.linalg.solve(hs, gs)
            x = x + frame @ step
            x /= np.linalg.norm(x)
            if np.linalg.norm(step) < 1e-14:
                break
        if any(np.linalg.norm(c.where - x) < 1e-6 for c in found):
            continue
        a, g, h = basis.values(weights, x[None, :])
        e1, e2 = tangent_basis(x)
        frame = np.c_[e1, e2]
        if np.linalg.norm(frame.T @ g[0]) > 1e-9 * max(1.0, float(np.abs(value).max())):
            continue
        eig = np.linalg.eigvalsh(frame.T @ h[0] @ frame - float(x @ g[0]) * np.eye(2))
        kind = "maximum" if eig.max() < 0 else "minimum" if eig.min() > 0 else "saddle"
        found.append(CriticalPoint(x, float(a[0]), kind))
    return found


# ---- a cylinder's TE modes and their real sums -----------------------------------------------------------------------
Partial = tuple[int, int, int]
FIRST: tuple[Partial, ...] = ((1, 0, 0), (0, 1, 0))
SECOND: tuple[Partial, ...] = ((2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 2, 0), (0, 1, 1))
THIRD: tuple[Partial, ...] = (
    (3, 0, 0),
    (2, 1, 0),
    (2, 0, 1),
    (1, 2, 0),
    (1, 1, 1),
    (1, 0, 2),
    (0, 3, 0),
    (0, 2, 1),
    (0, 1, 2),
)


@dataclass(frozen=True)
class CylinderMode:
    """TE_mnp of a cylinder of radius 1 and height d (lengths in units of the radius): Pi = J_m(x'_mn rho) cos(m phi)
    sin(p pi z / d), or sin(m phi) for the partner, with x'_mn the n-th zero of J_m'."""

    m: int
    n: int
    p: int
    d: float
    x: float
    sine: bool = False

    @staticmethod
    def of(m: int, n: int, p: int, d: float, sine: bool = False) -> CylinderMode:
        return CylinderMode(m, n, p, d, float(bessel_zeros(m, n, derivative=True)[-1]), sine)

    @property
    def beta(self) -> float:
        return self.p * math.pi / self.d

    @property
    def ka(self) -> float:
        return math.hypot(self.x, self.beta)

    @property
    def name(self) -> str:
        return f"TE{self.m}{self.n}{self.p}"


def planar_partials(mode: CylinderMode, xy: Arr, wanted: Sequence[tuple[int, int]]) -> dict[tuple[int, int], Arr]:
    """d^i/dx^i d^j/dy^j of J_m(x' rho) cos(m phi) (or sin) at points xy (P, 2), by the ladder D = d_x + i d_y,
    D [J_m(k rho) e^(i m phi)] = -k J_(m+1) e^(i (m+1) phi), and D-bar the same with k J_(m-1) e^(i (m-1) phi)."""
    top = max(i + j for i, j in wanted)
    rho, phi = np.hypot(xy[:, 0], xy[:, 1]), np.arctan2(xy[:, 1], xy[:, 0])
    shifts = list(range(-top, top + 1))
    ladder = bessel_ladder([mode.m + s for s in shifts], mode.x * rho)
    f = {s: ladder[s + top] * np.exp(1j * (mode.m + s) * phi) for s in shifts}
    k = mode.x
    out: dict[tuple[int, int], Arr] = {}
    for i, j in wanted:
        total: CArr = np.zeros(len(xy), dtype=np.complex128)
        for a in range(i + 1):
            for b in range(j + 1):
                up, down = a + b, (i - a) + (j - b)
                total = (
                    total + math.comb(i, a) * math.comb(j, b) * (-1) ** (j - b) * (-k) ** up * k**down * f[up - down]
                )
        total = total * 0.5 ** (i + j) * (-1j) ** j
        out[(i, j)] = total.imag if mode.sine else total.real
    return out


def axial_partial(mode: CylinderMode, z: Arr, order: int) -> Arr:
    """d^k/dz^k sin(beta z) = beta^k sin(beta z + k pi / 2)."""
    out: Arr = mode.beta**order * np.sin(mode.beta * z + order * math.pi / 2)
    return out


@dataclass(frozen=True)
class CylinderField:
    """A real standing sum of a cylinder's TE modes, Pi = sum_n w_n Pi_n, E = (Pi_y, -Pi_x, 0) (radius 1)."""

    modes: tuple[CylinderMode, ...]
    weights: tuple[float, ...]

    @property
    def d(self) -> float:
        return self.modes[0].d

    def partials(self, points: Arr, wanted: Sequence[Partial]) -> dict[Partial, Arr]:
        p = np.atleast_2d(points)
        planar = sorted({(i, j) for i, j, _ in wanted})
        out: dict[Partial, Arr] = {key: np.zeros(len(p)) for key in wanted}
        for mode, weight in zip(self.modes, self.weights, strict=True):
            flat = planar_partials(mode, p[:, :2], planar)
            for i, j, k in wanted:
                out[(i, j, k)] = out[(i, j, k)] + weight * flat[(i, j)] * axial_partial(mode, p[:, 2], k)
        return out

    def inside(self, points: Arr, margin: float) -> Any:
        p = np.atleast_2d(points)
        return (np.hypot(p[:, 0], p[:, 1]) < 1 - margin) & (p[:, 2] > margin) & (p[:, 2] < self.d - margin)

    def wall_distance(self, points: Arr) -> Arr:
        p = np.atleast_2d(points)
        out: Arr = np.minimum(1 - np.hypot(p[:, 0], p[:, 1]), np.minimum(p[:, 2], self.d - p[:, 2]))
        return out

    def mean_square(self) -> float:
        """The mean of |E|^2 = |grad_perp Pi|^2 over the volume: sum_n w_n^2 x'^2 (1 - m^2/x'^2) J_m(x')^2 / 4, the
        modes being orthogonal (distinct m or p) and each cos^2 m phi and sin^2 p pi z / d averaging 1/2."""
        total = 0.0
        for mode, weight in zip(self.modes, self.weights, strict=True):
            jm = float(bessel_j(mode.m, mode.x))
            azimuth = 1.0 if mode.m == 0 else 0.5
            total += weight**2 * (mode.x**2 - mode.m**2) * jm**2 * azimuth / 2
        return total


def cylinder_q(mode: CylinderMode, radius: float, skin: float = SKIN) -> float:
    """A TE_mnp mode's Q in a copper cylinder: (2 k^2 / delta) int |E|^2 dV / oint |curl E|_t^2 dS, in closed form
    k^2 a^4 d (1 - m^2/x'^2) / (2 x'^2 delta) / [(a d / 2)(1 + (beta a m / x'^2)^2) + (beta a^2 / x')^2
    (1 - m^2/x'^2)], the textbook formula (radius a and height d in metres here)."""
    a, d = radius, mode.d * radius
    beta, x, m = mode.p * math.pi / d, mode.x, mode.m
    k2 = (x / a) ** 2 + beta**2
    share = 1 - (m / x) ** 2
    walls = (a * d / 2) * (1 + (beta * a * m / x**2) ** 2) + (beta * a**2 / x) ** 2 * share
    return k2 * a**4 * d * share / (2 * x**2 * skin) / walls


def mode_energy(mode: CylinderMode) -> float:
    """int |grad_perp Pi|^2 dV for the mode with unit weight, over the cylinder of radius 1 and height d."""
    jm = float(bessel_j(mode.m, mode.x))
    azimuth = 2 * math.pi if mode.m == 0 else math.pi
    return (mode.x**2 - mode.m**2) * jm**2 / 2 * azimuth * mode.d / 2


# ---- a cube's TE_z modes ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class CubeField:
    """A real standing sum of a cube's TE_z modes (side 1): Pi = sum w cos(m pi x) cos(n pi y) sin(p pi z)."""

    indices: tuple[tuple[int, int, int], ...]
    weights: tuple[float, ...]
    d: float = 1.0

    def partials(self, points: Arr, wanted: Sequence[Partial]) -> dict[Partial, Arr]:
        p = np.atleast_2d(points)
        out: dict[Partial, Arr] = {key: np.zeros(len(p)) for key in wanted}
        for (m, n, q), weight in zip(self.indices, self.weights, strict=True):
            km, kn, kq = m * math.pi, n * math.pi, q * math.pi
            for i, j, k in wanted:
                out[(i, j, k)] = out[(i, j, k)] + weight * (
                    km**i
                    * np.cos(km * p[:, 0] + i * math.pi / 2)
                    * kn**j
                    * np.cos(kn * p[:, 1] + j * math.pi / 2)
                    * kq**k
                    * np.sin(kq * p[:, 2] + k * math.pi / 2)
                )
        return out

    def inside(self, points: Arr, margin: float) -> Any:
        p = np.atleast_2d(points)
        return np.all((p > margin) & (p < 1 - margin), axis=1)

    def wall_distance(self, points: Arr) -> Arr:
        p = np.atleast_2d(points)
        out: Arr = np.minimum(p, 1 - p).min(axis=1)
        return out


Field = CylinderField | CubeField


# ---- null curves of a real TE field: Pi_x = Pi_y = 0 ---------------------------------------------------------------
def pair(field: Field, points: Arr) -> Arr:
    """(Pi_x, Pi_y) at points (P, 3): E = (Pi_y, -Pi_x, 0) vanishes where both do."""
    d = field.partials(points, FIRST)
    out: Arr = np.stack([d[(1, 0, 0)], d[(0, 1, 0)]], axis=-1)
    return out


def pair_jacobian(field: Field, points: Arr) -> Arr:
    """The 2 x 3 Jacobian of (Pi_x, Pi_y): rows grad Pi_x and grad Pi_y, shape (P, 2, 3)."""
    d = field.partials(points, SECOND)
    xx, xy, xz, yy, yz = (d[key] for key in SECOND)
    out: Arr = np.stack([np.stack([xx, xy, xz], -1), np.stack([xy, yy, yz], -1)], -2)
    return out


def tangents(field: Field, points: Arr) -> Arr:
    """grad Pi_x x grad Pi_y, the null curve's orienting tangent (unnormalised): its z-part is det S, S the in-slice
    Hessian (CAV-5)."""
    j = pair_jacobian(field, points)
    out: Arr = np.cross(j[:, 0, :], j[:, 1, :])
    return out


def project(field: Field, points: Arr, iterations: int = 12, tolerance: float = 1e-13) -> Arr:
    """Newton's method onto the null set: the least-norm step -J^T (J J^T)^-1 (Pi_x, Pi_y), across the curve."""
    x = np.atleast_2d(points).copy()
    for _ in range(iterations):
        g, j = pair(field, x), pair_jacobian(field, x)
        jjt = j @ np.swapaxes(j, 1, 2)
        det = jjt[:, 0, 0] * jjt[:, 1, 1] - jjt[:, 0, 1] * jjt[:, 1, 0]
        det = np.where(np.abs(det) < 1e-300, 1e-300, det)
        inverse = np.stack(
            [np.stack([jjt[:, 1, 1], -jjt[:, 0, 1]], -1), np.stack([-jjt[:, 1, 0], jjt[:, 0, 0]], -1)], -2
        )
        lam = np.einsum("pij,pj->pi", inverse, g) / det[:, None]
        step = np.einsum("pji,pj->pi", j, lam)
        x = x - step
        if float(np.abs(step).max()) < tolerance:
            break
    return x


@dataclass(frozen=True)
class NullCurve:
    """A traced null curve: its points in order, and whether it closes (a loop) or runs from wall to wall."""

    points: Arr
    closed: bool

    @property
    def length(self) -> float:
        steps = np.diff(np.vstack([self.points, self.points[:1]]) if self.closed else self.points, axis=0)
        return float(np.linalg.norm(steps, axis=1).sum())

    @property
    def arc(self) -> Arr:
        steps = np.linalg.norm(np.diff(self.points, axis=0), axis=1)
        out: Arr = np.r_[0.0, np.cumsum(steps)]
        return out


def trace(field: Field, start: Arr, step: float, limit: int, margin: float, sign: float = 1.0) -> tuple[Arr, str]:
    """Follow the null curve through start, along sign times grad Pi_x x grad Pi_y: a midpoint predictor and Newton's
    method across the curve. Ends 'closed' when it passes back within half a step of start, 'exit' at the walls."""
    x = start.copy()
    points = [x.copy()]
    previous = sign * tangents(field, x[None, :])[0]
    for count in range(limit):
        t1 = tangents(field, x[None, :])[0]
        t1 = t1 / np.linalg.norm(t1) * (1 if t1 @ previous >= 0 else -1)
        t2 = tangents(field, (x + 0.5 * step * t1)[None, :])[0]
        t2 = t2 / np.linalg.norm(t2) * (1 if t2 @ t1 >= 0 else -1)
        new = project(field, (x + step * t2)[None, :], iterations=8)[0]
        if count > 8:
            seg = new - x
            s = float(np.clip((start - x) @ seg / (seg @ seg), 0.0, 1.0))
            if np.linalg.norm(x + s * seg - start) < 0.5 * step:
                return np.array(points), "closed"
        points.append(new.copy())
        if not bool(field.inside(new[None, :], margin)[0]):
            return np.array(points), "exit"
        x, previous = new, t2
    return np.array(points), "limit"


def null_curves(field: Field, starts: Arr, step: float, margin: float, limit: int = 20000) -> list[NullCurve]:
    """Every null curve through a set of starting points: Newton from each start, then a trace from each converged
    point that lies off the curves already traced. Loops are kept oriented by grad Pi_x x grad Pi_y."""
    with np.errstate(over="ignore", invalid="ignore"):
        seeds = project(field, starts, iterations=14)
        scale = float(np.abs(pair(field, starts)).max())
        ok = (
            np.isfinite(seeds).all(axis=1)
            & (np.abs(pair(field, seeds)).max(axis=1) < 1e-10 * scale)
            & field.inside(seeds, 2 * margin)
            & (np.linalg.norm(seeds - starts, axis=1) < 0.3)
        )
    traced: list[NullCurve] = []
    cloud = np.zeros((0, 3))
    for seed in seeds[ok]:
        if len(cloud) and float(np.min(np.linalg.norm(cloud - seed, axis=1))) < 3 * step:
            continue
        with np.errstate(over="ignore", invalid="ignore"):
            forward, how = trace(field, seed, step, limit, margin)
            if how == "closed":
                curve = NullCurve(forward, True)
            else:
                backward, _ = trace(field, seed, step, limit, margin, sign=-1.0)
                curve = NullCurve(np.vstack([backward[::-1], forward[1:]]), False)
        traced.append(curve)
        cloud = np.vstack([cloud, curve.points])
    return traced


def linking_number(a: Arr, b: Arr) -> float:
    """Gauss's linking integral of two closed polygons."""
    da, db = np.roll(a, -1, 0) - a, np.roll(b, -1, 0) - b
    ma, mb = a + da / 2, b + db / 2
    r = ma[:, None, :] - mb[None, :, :]
    cross = np.cross(da[:, None, :], db[None, :, :])
    return float((np.einsum("ijk,ijk->ij", cross, r) / np.linalg.norm(r, axis=2) ** 3).sum() / (4 * math.pi))


def folds(curve: NullCurve, field: Field) -> list[int]:
    """Indices where the curve's height turns: where the z-part of its tangent, det S, changes sign along it."""
    tz = tangents(field, curve.points)[:, 2]
    if curve.closed:
        return [int(i) for i in np.flatnonzero(np.sign(tz) != np.sign(np.roll(tz, -1)))]
    return [int(i) for i in np.flatnonzero(np.sign(tz[:-1]) != np.sign(tz[1:]))]


def fold_point(curve: NullCurve, field: Field, index: int) -> Arr:
    """The fold between points index and index + 1, where det S = 0 on the curve: bisection on the chord's fraction,
    each trial projected onto the curve."""
    a, b = curve.points[index], curve.points[(index + 1) % len(curve.points)]
    lo, hi = 0.0, 1.0
    t_lo = tangents(field, a[None, :])[0, 2]
    x: Arr = a
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        x = project(field, (a + mid * (b - a))[None, :], iterations=10)[0]
        t_mid = tangents(field, x[None, :])[0, 2]
        if np.sign(t_mid) == np.sign(t_lo):
            lo, t_lo = mid, t_mid
        else:
            hi = mid
    return x


# ---- wells: |E|^2, its gradient and Hessian, floods and saddles ----------------------------------------------------
def u_value(field: Field, points: Arr) -> Arr:
    g = pair(field, points)
    out: Arr = (g**2).sum(axis=1)
    return out


def u_gradient_hessian(field: Field, points: Arr) -> tuple[Arr, Arr, Arr]:
    """|E|^2 = Pi_x^2 + Pi_y^2 at points, its gradient 2 (Pi_x grad Pi_x + Pi_y grad Pi_y), and its Hessian."""
    d = field.partials(points, FIRST + SECOND + THIRD)

    def get(i: int, j: int, k: int) -> Arr:
        return d[(i, j, k)]

    px, py = get(1, 0, 0), get(0, 1, 0)
    gx = np.stack([get(2, 0, 0), get(1, 1, 0), get(1, 0, 1)], -1)
    gy = np.stack([get(1, 1, 0), get(0, 2, 0), get(0, 1, 1)], -1)
    hx = np.stack(
        [
            np.stack([get(3, 0, 0), get(2, 1, 0), get(2, 0, 1)], -1),
            np.stack([get(2, 1, 0), get(1, 2, 0), get(1, 1, 1)], -1),
            np.stack([get(2, 0, 1), get(1, 1, 1), get(1, 0, 2)], -1),
        ],
        -2,
    )
    hy = np.stack(
        [
            np.stack([get(2, 1, 0), get(1, 2, 0), get(1, 1, 1)], -1),
            np.stack([get(1, 2, 0), get(0, 3, 0), get(0, 2, 1)], -1),
            np.stack([get(1, 1, 1), get(0, 2, 1), get(0, 1, 2)], -1),
        ],
        -2,
    )
    value = px**2 + py**2
    grad = 2 * (px[:, None] * gx + py[:, None] * gy)
    hess = 2 * (
        gx[:, :, None] * gx[:, None, :]
        + gy[:, :, None] * gy[:, None, :]
        + px[:, None, None] * hx
        + py[:, None, None] * hy
    )
    return value, grad, hess


def critical_newton(field: Field, start: Arr, iterations: int = 60) -> Arr:
    """Newton's method on grad |E|^2 = 0 from start, the step capped at 0.02 of the radius."""
    x = start.copy()
    for _ in range(iterations):
        _, g, h = u_gradient_hessian(field, x[None, :])
        step = -np.linalg.solve(h[0], g[0])
        size = float(np.linalg.norm(step))
        if size > 0.02:
            step *= 0.02 / size
        x = x + step
        if size < 1e-14:
            break
    return x


def minimax_flood(values: Arr, sources: IntArr, target: Any) -> tuple[float, int]:
    """The lowest level at which the sublevel set of a grid's values joins the source voxels to a target voxel, by
    growing the region from the sources through its lowest neighbour first (six-connected); returns the level and the
    voxel at which it was met, the saddle's grid neighbour."""
    shape = values.shape
    flat = values.reshape(-1)
    goal = np.asarray(target).reshape(-1)
    strides = (shape[1] * shape[2], shape[2], 1)
    seen = np.zeros(flat.size, dtype=bool)
    heap = [(float(flat[s]), int(s)) for s in sources if math.isfinite(float(flat[s]))]
    heapq.heapify(heap)
    level, where = -math.inf, -1
    while heap:
        value, index = heapq.heappop(heap)
        if seen[index]:
            continue
        seen[index] = True
        if value > level:
            level, where = value, index
        if goal[index]:
            return level, where
        i, rest = divmod(index, strides[0])
        j, k = divmod(rest, strides[1])
        for axis, (coordinate, size) in enumerate(((i, shape[0]), (j, shape[1]), (k, shape[2]))):
            for move in (-1, 1):
                if 0 <= coordinate + move < size:
                    neighbour = index + move * strides[axis]
                    if not seen[neighbour] and math.isfinite(float(flat[neighbour])):
                        heapq.heappush(heap, (float(flat[neighbour]), neighbour))
    return math.inf, -1


def descent(field: Field, start: Arr, step: float, stop_value: float, limit: int = 4000) -> Arr:
    """Steepest descent of |E|^2 from start, in steps of fixed length, until |E|^2 falls under stop_value (a null) or
    the path meets an end cap or the side wall."""
    x = start.copy()
    path = [x.copy()]
    for _ in range(limit):
        _, g, _ = u_gradient_hessian(field, x[None, :])
        direction = -g[0] / np.linalg.norm(g[0])
        half = x + 0.5 * step * direction
        _, g2, _ = u_gradient_hessian(field, half[None, :])
        x = x - step * g2[0] / np.linalg.norm(g2[0])
        path.append(x.copy())
        if float(u_value(field, x[None, :])[0]) < stop_value or float(field.wall_distance(x[None, :])[0]) < step:
            break
    return np.array(path)


def cylinder_grid(field: CylinderField, xs: Arr, ys: Arr, zs: Arr) -> tuple[Arr, Arr]:
    """(Pi_x, Pi_y) on the grid xs x ys x zs, shape (nx, ny, nz) each, by the modes' separation into a planar part
    and sin(beta z): each planar part is computed once on the plane."""
    gx, gy = np.meshgrid(xs, ys, indexing="ij")
    flat = np.c_[gx.ravel(), gy.ravel()]
    px = np.zeros((len(xs), len(ys), len(zs)))
    py = np.zeros_like(px)
    for mode, weight in zip(field.modes, field.weights, strict=True):
        planar = planar_partials(mode, flat, [(1, 0), (0, 1)])
        axial = axial_partial(mode, zs, 0)
        px += weight * planar[(1, 0)].reshape(gx.shape)[:, :, None] * axial[None, None, :]
        py += weight * planar[(0, 1)].reshape(gx.shape)[:, :, None] * axial[None, None, :]
    return px, py


# ---- the sphere: families, feeds and null diameters ----------------------------------------------------------------
def sphere_q(ka: float, degree: int, radius: float, tm: bool) -> float:
    """A sphere family's copper Q: a / delta for TE, (a / delta)(1 - l(l + 1)/(ka)^2) for TM (Claim CAV-10)."""
    return radius / SKIN * ((1 - degree * (degree + 1) / ka**2) if tm else 1.0)


def wall_point(theta: float, phi: float, turn: float) -> tuple[Arr, Arr]:
    """A point on the unit sphere and a unit tangent there, turned by `turn` from the theta direction."""
    x = np.array([math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)])
    e_theta = np.array([math.cos(theta) * math.cos(phi), math.cos(theta) * math.sin(phi), -math.sin(theta)])
    e_phi = np.array([-math.sin(phi), math.cos(phi), 0.0])
    return x, math.cos(turn) * e_theta + math.sin(turn) * e_phi


def wall_coupling(basis: Harmonics, places: Arr, moments: Arr) -> Arr:
    """C[n, k]: how strongly a small loop at wall point k with moment along the tangent moments[k] couples to the
    family's mode n. A TE mode's tangential H at the wall is along the surface gradient of its angular part,
    grad_S h = grad h - l h x on the unit sphere, so C[n, k] = moments[k] . grad_S h_n(places[k])."""
    value, grad = basis.gradients(places)
    surface = grad - basis.degree * value[:, :, None] * places[None, :, :]
    out: Arr = np.einsum("bpi,pi->bp", surface, moments)
    return out


def inverse_condition(c: Arr) -> float:
    s = np.linalg.svd(c, compute_uv=False)
    return float(s[-1] / s[0]) if c.shape[1] >= c.shape[0] else 0.0


def place_feeds(basis: Harmonics, count: int, pickups: int, seed: int, rounds: int) -> dict[str, Any]:
    """count loop feeds on the sphere's wall placed and turned for a well-conditioned coupling matrix C (modes x
    feeds), by a seeded random search and then local moves that keep any improvement of sigma_min / sigma_max; then
    pickups added one at a time where they most raise the smallest singular value of the modes-by-ports matrix."""
    rng = np.random.default_rng(seed)

    def build(params: Arr) -> tuple[Arr, Arr]:
        pairs = [wall_point(float(t), float(p), float(u)) for t, p, u in params]
        return np.array([q[0] for q in pairs]), np.array([q[1] for q in pairs])

    def draw(n: int) -> Arr:
        out: Arr = np.c_[np.arccos(rng.uniform(-1, 1, n)), rng.uniform(0, 2 * math.pi, n), rng.uniform(0, math.pi, n)]
        return out

    best = draw(count)
    score = inverse_condition(wall_coupling(basis, *build(best)))
    for _ in range(rounds // 4):
        trial = draw(count)
        value = inverse_condition(wall_coupling(basis, *build(trial)))
        if value > score:
            best, score = trial, value
    for r in range(rounds):
        trial = best + rng.normal(0, 0.3 * (1 - r / rounds) + 0.01, best.shape) * (rng.uniform(size=(count, 1)) < 0.4)
        value = inverse_condition(wall_coupling(basis, *build(trial)))
        if value > score:
            best, score = trial, value
    places, moments = build(best)
    coupling = wall_coupling(basis, places, moments)
    picks: list[tuple[Arr, Arr]] = []
    for _ in range(pickups):
        candidates = draw(1500)
        rated = []
        for cand in candidates:
            x, t = wall_point(float(cand[0]), float(cand[1]), float(cand[2]))
            stacked = np.c_[coupling, *(wall_coupling(basis, q[0][None, :], q[1][None, :]) for q in picks)]
            extra = wall_coupling(basis, x[None, :], t[None, :])
            rated.append(float(np.linalg.svd(np.c_[stacked, extra], compute_uv=False)[-1]))
        chosen = candidates[int(np.argmax(rated))]
        picks.append(wall_point(float(chosen[0]), float(chosen[1]), float(chosen[2])))
    singular = np.linalg.svd(coupling, compute_uv=False)
    return {
        "feeds": places,
        "moments": moments,
        "pickups": np.array([q[0] for q in picks]) if picks else np.zeros((0, 3)),
        "pickup_moments": np.array([q[1] for q in picks]) if picks else np.zeros((0, 3)),
        "condition": float(singular[0] / singular[-1]),
    }


def compute_sphere(degree: int, family: int, seed: int, quick: bool, feeds: bool) -> dict[str, Any]:
    """One TE family of a sphere at 2.45 GHz and one real field of it, psi = j_l(kr) A: its radius, its modes, its
    copper Q, the critical points of A (the null diameters, CAV-2), the inner zero-field shells, and loop feeds."""
    ka = float(sphere_zeros(degree, family)[-1])
    radius = ka / K_DRIVE
    shells = sphere_zeros(degree, family)[:-1] / ka
    basis = Harmonics.of(degree)
    weights = np.random.default_rng(seed).normal(size=2 * degree + 1)
    points = critical_points(basis, weights, samples=3000 if quick else 8000)
    kinds = [c.kind for c in points]
    count = {kind: kinds.count(kind) for kind in ("maximum", "saddle", "minimum")}
    out: dict[str, Any] = {
        "degree": degree,
        "family": family,
        "ka": ka,
        "radius_cm": radius * 100,
        "modes": 2 * degree + 1,
        "Q": sphere_q(ka, degree, radius, tm=False),
        "shells": shells,
        "weights": weights,
        "critical": count,
        "morse_sum": count["maximum"] - count["saddle"] + count["minimum"],
        "diameters": len(points) // 2,
        "_points": points,
        "_basis": basis,
    }
    say(
        f"  sphere TE l = {degree}, n = {family}: ka {ka:.4f}, radius {radius * 100:.2f} cm, {2 * degree + 1} modes, "
        f"Q {out['Q']:,.0f}; one real field (seed {seed}): {len(points)} critical points of A ({count['maximum']} "
        f"maxima, {count['saddle']} saddles, {count['minimum']} minima; max - saddles + min = {out['morse_sum']}), so "
        f"{out['diameters']} null diameters"
        + (f"; zero-field shells at r/a = {', '.join(f'{v:.4f}' for v in shells)}" if len(shells) else "")
    )
    if feeds:
        placed = place_feeds(basis, 2 * degree + 1, 3, seed + 1, 400 if quick else 3000)
        out.update({k: v for k, v in placed.items() if k == "condition"})
        out["_feeds"] = placed
        say(
            f"  {2 * degree + 1} loop feeds on its wall, placed and turned by a seeded search: the coupling matrix's "
            f"condition number {placed['condition']:.2f}; three pickups added where they most raise its smallest "
            "singular value"
        )
    return out


# ---- the worked cylinder: TE113 + TE221 --------------------------------------------------------------------------
def worked_cylinder() -> tuple[CylinderField, float]:
    """The worked example of section 7: TE113 and TE221 tuned together by d/a, equal weights; and its radius (m)."""
    x11, x22 = float(bessel_zeros(1, 1, True)[0]), float(bessel_zeros(2, 2, True)[1])
    d = math.pi / math.sqrt((x22**2 - x11**2) / (3**2 - 1**2))
    m1, m2 = CylinderMode.of(1, 1, 3, d), CylinderMode.of(2, 2, 1, d)
    return CylinderField((m1, m2), (1.0, 1.0)), m1.ka / K_DRIVE


def well_depth_ev(q: float, saddle_share: float, power: float, volume: float) -> float:
    """The depth (eV) of a well round a null whose escape saddle has saddle_share of the mean |E|^2, in a cavity of
    volume V (m^3) and quality Q fed with P (W): (e/m) Q s P / (2 eps0 w^3 V) (Proposition CAV-7)."""
    return QM * q * saddle_share * power / (2 * pt.EPS0 * OMEGA**3 * volume)


def to_wall(end: Arr, previous: Arr, d: float) -> float:
    """How far the end of a traced open null line lies from the wall of the cylinder (radius 1, height d) along the
    line's end tangent, the chord from the point before it: the tracer stops within a step of the wall."""
    u = (end - previous) / np.linalg.norm(end - previous)
    a, b, c = u[0] ** 2 + u[1] ** 2, 2 * (end[0] * u[0] + end[1] * u[1]), end[0] ** 2 + end[1] ** 2 - 1
    if c >= 0 or not 0 < end[2] < d:
        return 0.0
    hits = [(-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)] if a > 0 else []
    if u[2] != 0:
        hits.append(((d if u[2] > 0 else 0.0) - end[2]) / u[2])
    return float(min(hits))


def compute_worked(quick: bool) -> dict[str, Any]:
    say("\n== worked: TE113 + TE221 in a copper cylinder at 2.45 GHz (section 7) ==")
    field, radius = worked_cylinder()
    m1, m2 = field.modes
    d = field.d
    height, volume = d * radius, math.pi * radius**2 * d * radius
    q1, q2 = cylinder_q(m1, radius), cylinder_q(m2, radius)
    e1, e2 = mode_energy(m1), mode_energy(m2)
    share = e2 / (e1 + e2)
    q_mix = 1 / ((1 - share) / q1 + share / q2)
    mean = field.mean_square()
    e_rms = math.sqrt(2 * q_mix / (pt.EPS0 * OMEGA * volume))
    say(
        f"  x'11 = {m1.x:.6f}, x'22 = {m2.x:.6f}: d/a = {d:.4f}, ka = {m1.ka:.4f} (TE221 {m2.ka:.4f}); a = "
        f"{radius * 100:.2f} cm, d = {height * 100:.2f} cm, {volume * 1e3:.2f} litres"
    )
    say(
        f"  Q {q1:,.0f} (TE113), {q2:,.0f} (TE221), {q_mix:,.0f} for the mix ({share * 100:.1f} % of its energy in "
        f"TE221); linewidth {F_DRIVE / q_mix / 1e3:.1f} kHz; E_rms at 1 W {e_rms / 1e3:.2f} kV/m"
    )
    # the null curves
    ng, step = (16, 0.008) if quick else (24, 0.004)
    axis = np.linspace(-0.97, 0.97, ng)
    heights = np.linspace(0.03 * d, 0.97 * d, int(ng * d / 2) + 2)
    starts = np.stack(np.meshgrid(axis, axis, heights, indexing="ij"), -1).reshape(-1, 3)
    starts = starts[np.hypot(starts[:, 0], starts[:, 1]) < 0.97]
    curves = null_curves(field, starts, step, 0.004)
    loops = sorted((c for c in curves if c.closed), key=lambda c: (-round(c.length, 3), float(c.points[:, 1].mean())))
    opens = sorted((c for c in curves if not c.closed), key=lambda c: -c.length)
    rows: list[dict[str, Any]] = []
    for i, loop in enumerate(loops):
        rho = np.hypot(loop.points[:, 0], loop.points[:, 1])
        turns = folds(loop, field)
        rows.append(
            {
                "length_cm": loop.length * radius * 100,
                "rho": (float(rho.min()), float(rho.max())),
                "z_over_d": (float(loop.points[:, 2].min() / d), float(loop.points[:, 2].max() / d)),
                "wall_cm": float(field.wall_distance(loop.points).min() * radius * 100),
                "folds": len(turns),
                "in_plane_y0": bool(np.abs(loop.points[:, 1]).max() < 1e-9),
            }
        )
        say(
            f"  loop {i}: {rows[-1]['length_cm']:.2f} cm, rho {rho.min():.3f}-{rho.max():.3f} a, z/d "
            f"{rows[-1]['z_over_d'][0]:.3f}-{rows[-1]['z_over_d'][1]:.3f}, {rows[-1]['wall_cm']:.2f} cm from the "
            f"nearest wall, {len(turns)} folds" + (", in the plane y = 0" if rows[-1]["in_plane_y0"] else "")
        )
    links = {
        f"{i},{j}": linking_number(loops[i].points, loops[j].points)
        for i in range(len(loops))
        for j in range(i + 1, len(loops))
    }
    ends, spans = [], []
    for line in opens:
        kinds = []
        for end in (line.points[0], line.points[-1]):
            kinds.append("cap" if min(end[2], d - end[2]) < 0.02 else "side wall")
        ends.append(" to ".join(kinds))
        spans.append(
            line.length + to_wall(line.points[0], line.points[1], d) + to_wall(line.points[-1], line.points[-2], d)
        )
    say(
        "  open null lines, traced to within a step of the walls and carried on along their end tangents to them: "
        + "; ".join(f"{span * radius * 100:.2f} cm from {how}" for span, how in zip(spans, ends, strict=True))
        + "; linking numbers of the loops "
        + ", ".join(f"({k}) {v:+.1e}" for k, v in links.items())
    )
    # the well: a minimax flood from each loop, then Newton on grad |E|^2 = 0
    n = 61 if quick else 101
    xs = np.linspace(-1, 1, n)
    h = xs[1] - xs[0]
    nz = round(d / h) + 1
    zs = np.linspace(0, d, nz)
    px, py = cylinder_grid(field, xs, xs, zs)
    rho_grid = np.hypot(*np.meshgrid(xs, xs, indexing="ij"))[:, :, None] * np.ones((1, 1, nz))
    grid_u = np.where(rho_grid <= 1, (px**2 + py**2) / mean, np.inf)
    layer = np.arange(nz)[None, None, :] * np.ones((n, n, 1))

    def voxels(points: Arr) -> IntArr:
        i = np.clip(np.round((points[:, 0] + 1) / h).astype(np.intp), 0, n - 1)
        j = np.clip(np.round((points[:, 1] + 1) / h).astype(np.intp), 0, n - 1)
        k = np.clip(np.round(points[:, 2] / (zs[1] - zs[0])).astype(np.intp), 0, nz - 1)
        out: IntArr = np.ravel_multi_index((i, j, k), grid_u.shape).astype(np.intp)
        return out

    walls = ((layer == 0) | (layer == nz - 1) | (rho_grid > 1 - 1.5 * h)).reshape(-1)
    barriers: list[dict[str, Any]] = []
    saddles: list[Arr] = []
    for i, loop in enumerate(loops):
        target = walls.copy()
        for j, other in enumerate(loops + opens):
            if j != i:
                target[voxels(other.points)] = True
        level, where = minimax_flood(grid_u, voxels(loop.points), target)
        a, b, c = np.unravel_index(where, grid_u.shape)
        saddle = critical_newton(field, np.array([xs[a], xs[b], zs[c]]))
        value, grad, hess = u_gradient_hessian(field, saddle[None, :])
        eig = np.linalg.eigvalsh(hess[0])
        barriers.append(
            {"grid": level, "saddle": float(value[0] / mean), "where": saddle, "index": int((eig < 0).sum())}
        )
        saddles.append(saddle)
        say(
            f"  loop {i}: the flood on a {n} x {n} x {nz} grid leaves at {level:.4f} of the mean |E|^2; Newton puts "
            "the "
            f"saddle at {value[0] / mean:.5f}, rho {math.hypot(saddle[0], saddle[1]):.3f} a, phi "
            f"{math.degrees(math.atan2(saddle[1], saddle[0])):.1f} deg, z/d {saddle[2] / d:.3f} (index "
            f"{int((eig < 0).sum())}, |grad| {np.linalg.norm(grad):.0e})"
        )
    s = min(row["saddle"] for row in barriers)
    first = saddles[int(np.argmin([row["saddle"] for row in barriers]))]
    images = []
    for flip_y in (1.0, -1.0):
        for flip_z in (False, True):
            guess = np.array([first[0], flip_y * first[1], d - first[2] if flip_z else first[2]])
            images.append(critical_newton(field, guess))
    image_values = [float(u_value(field, p[None, :])[0] / mean) for p in images]
    depth = well_depth_ev(q_mix, s, 1.0, volume)
    ceiling = pt.free_space_depth_ceiling_ev(1.0)
    top = np.unravel_index(int(np.argmax(np.where(np.isfinite(grid_u), grid_u, -1.0))), grid_u.shape)
    summit = critical_newton(field, np.array([xs[top[0]], xs[top[1]], zs[top[2]]]))
    peak = max(float(grid_u[top]), float(u_value(field, summit[None, :])[0] / mean))
    say(
        "  escape saddles (the lowest and its mirror images in y = 0 and z = d/2): s = "
        + ", ".join(f"{v:.5f}" for v in image_values)
        + f"; depth {depth * 1e6:.1f} micro-eV per watt = {depth / ceiling:.2f} times RFG-2's ceiling; 25 meV takes "
        f"{0.025 / depth:.0f} W, 1 eV {1 / depth / 1e3:.2f} kW; peak |E|^2 {peak:.2f} times the mean, "
        f"{e_rms * math.sqrt(peak) / 1e3:.1f} kV/m at 1 W"
    )
    # the field across each loop: the two normal gradients, Mathieu's q, and B at the null
    gradients, b_ratio = [], []
    for loop in loops:
        jac = pair_jacobian(field, loop.points)
        gradients.append(np.linalg.svd(jac, compute_uv=False) / math.sqrt(mean))
        dd = field.partials(loop.points, SECOND)
        curl = np.stack([dd[(1, 0, 1)], dd[(0, 1, 1)], -(dd[(2, 0, 0)] + dd[(0, 2, 0)])], -1)
        b_ratio.append(np.linalg.norm(curl, axis=1) / (m1.ka * math.sqrt(mean)))
    medians = np.median(gradients[0], axis=0)
    unit = e_rms / radius  # V/m^2 per unit of E_rms / a at 1 W
    q_at = {p: 2 * QM * medians * unit * math.sqrt(p) / OMEGA**2 for p in (1.0, 1e3)}
    # an electron moving at beta c through B in quadrature quivers by e v B / (m w), whose mean energy over the depth
    # at the loop's escape saddle is beta^2 (c|B| / E_rms)^2 / s at most (CAV-18)
    magnetic = max(float(b.max()) ** 2 / row["saddle"] for b, row in zip(b_ratio, barriers, strict=True))
    say(
        f"  across loop 0 the gradients run {gradients[0][:, 0].min():.2f}-{gradients[0][:, 0].max():.2f} and "
        f"{gradients[0][:, 1].min():.2f}-{gradients[0][:, 1].max():.2f} E_rms/a (medians {medians[0]:.3f}, "
        f"{medians[1]:.3f}): Mathieu q {sci(q_at[1.0][0])} and {sci(q_at[1.0][1])} at 1 W, secular "
        f"{q_at[1.0][0] * F_DRIVE / (2 * math.sqrt(2)) / 1e6:.2f} and "
        f"{q_at[1.0][1] * F_DRIVE / (2 * math.sqrt(2)) / 1e6:.2f}"
        f" MHz (q {q_at[1e3][0]:.4f} and {q_at[1e3][1]:.4f} at 1 kW); c|B|/E_rms on the loops "
        + ", ".join(f"{b.min():.2f}-{b.max():.2f}" for b in b_ratio)
        + f"; the magnetic quiver term at most {magnetic:.1f} beta^2 of a loop's depth"
    )
    # the feeds: loops at (+-a, 0, d/2) read H_z, k_c^2 J_m(x') cos(m phi) sin(beta z) on the wall
    feed_z = d / 2
    coupling = {}
    for mode in (m1, m2):
        for sine in (False, True):
            partner = CylinderMode(mode.m, mode.n, mode.p, d, mode.x, sine)
            hz = [
                partner.x**2
                * float(planar_partials(partner, np.array([[sx, 0.0]]), [(0, 0)])[(0, 0)][0])
                * float(axial_partial(partner, np.array([feed_z]), 0)[0])
                for sx in (1.0, -1.0)
            ]
            coupling[f"{partner.name}{' sin' if sine else ' cos'}"] = {
                "sum": hz[0] + hz[1],
                "difference": hz[0] - hz[1],
            }
    say(
        "  two loop feeds at (+-a, 0, d/2), coupling to H_z, driven as sum and difference: "
        + "; ".join(f"{k} {v['sum']:+.3f} / {v['difference']:+.3f}" for k, v in coupling.items())
    )
    return {
        "d_over_a": d,
        "ka": m1.ka,
        "radius_cm": radius * 100,
        "height_cm": height * 100,
        "volume_l": volume * 1e3,
        "Q": {"TE113": q1, "TE221": q2, "mix": q_mix},
        "share_TE221": share,
        "linewidth_kHz": F_DRIVE / q_mix / 1e3,
        "E_rms_1W": e_rms,
        "mean_square": mean,
        "loops": rows,
        "links": links,
        "open_lines": [{"length_cm": span * radius * 100, "ends": how} for span, how in zip(spans, ends, strict=True)],
        "barriers": barriers,
        "saddle_share": s,
        "saddles": images,
        "depth_ueV_per_W": depth * 1e6,
        "ceiling_ratio": depth / ceiling,
        "W_for_25meV": 0.025 / depth,
        "kW_for_1eV": 1 / depth / 1e3,
        "peak_over_mean": peak,
        "gradient_medians": medians,
        "q_1W": q_at[1.0],
        "q_1kW": q_at[1e3],
        "B_ratio": [(float(b.min()), float(b.max())) for b in b_ratio],
        "magnetic_over_beta2": magnetic,
        "feed_coupling": coupling,
        "_field": field,
        "_radius": radius,
        "_loops": loops,
        "_opens": opens,
        "_gradients": gradients,
        "_grid": (xs, zs, grid_u),
    }


# ---- section computations ----------------------------------------------------------------------------------------
def compute_models(quick: bool, worked: dict[str, Any]) -> dict[str, Any]:
    say("\n== models: the owner's sphere, and the worked cylinder (figure 1) ==")
    sphere = compute_sphere(3, 1, 2026, quick, feeds=True)
    return {"sphere": sphere}


def disc_determinant(field: CylinderField, xs: Arr) -> tuple[Arr, Arr]:
    """det(grad Phi_1, grad Phi_2) of the two modes' planar parts on the square grid xs x xs (NaN outside the
    disc), and its scale |grad Phi_1| |grad Phi_2| (CAV-3: every null lies over its zero set C)."""
    gx, gy = np.meshgrid(xs, xs, indexing="ij")
    flat = np.c_[gx.ravel(), gy.ravel()]
    one, two = (planar_partials(mode, flat, [(1, 0), (0, 1)]) for mode in field.modes)
    det = one[(1, 0)] * two[(0, 1)] - one[(0, 1)] * two[(1, 0)]
    scale = np.hypot(one[(1, 0)], one[(0, 1)]) * np.hypot(two[(1, 0)], two[(0, 1)])
    inside = np.hypot(flat[:, 0], flat[:, 1]) <= 1
    out: Arr = np.where(inside, det, np.nan).reshape(gx.shape)
    size: Arr = np.where(inside, scale, np.nan).reshape(gx.shape)
    return out, size


CUBE_INDICES = ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))
CUBE_SEED = 30  # the draw shown; most draws of this family hold only wall-to-wall lines (the census below)


def cube_curves(seed: int, quick: bool) -> tuple[CubeField, list[NullCurve]]:
    weights = np.random.default_rng(seed).normal(size=len(CUBE_INDICES))
    cube = CubeField(CUBE_INDICES, tuple(float(w) for w in weights))
    n, step = (12, 0.008) if quick else (16, 0.005)
    axis = np.linspace(0.03, 0.97, n)
    starts = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), -1).reshape(-1, 3)
    return cube, null_curves(cube, starts, step, 0.003)


def compute_nullsets(quick: bool, worked: dict[str, Any]) -> dict[str, Any]:
    say("\n== nullsets: a sphere family, the two-profile cylinder over C, and a cube field (figure 2) ==")
    sphere = compute_sphere(2, 2, 7, quick, feeds=False)
    field: CylinderField = worked["_field"]
    xs = np.linspace(-1, 1, 161 if quick else 401)
    det, scale = disc_determinant(field, xs)
    on_curves = []
    for curve in worked["_loops"] + worked["_opens"]:
        p = curve.points
        one, two = (planar_partials(mode, p[:, :2], [(1, 0), (0, 1)]) for mode in field.modes)
        value = one[(1, 0)] * two[(0, 1)] - one[(0, 1)] * two[(1, 0)]
        size = np.hypot(one[(1, 0)], one[(0, 1)]) * np.hypot(two[(1, 0)], two[(0, 1)])
        on_curves.append(float(np.max(np.abs(value) / size)))
    say(
        f"  the cylinder's nulls lie over C, det(grad Phi1, grad Phi2) = 0: on the traced curves |det| / (|grad Phi1| "
        f"|grad Phi2|) is at most {max(on_curves):.1e}"
    )
    cube, curves = cube_curves(CUBE_SEED, quick)
    loops = [c for c in curves if c.closed]
    links = [linking_number(a.points, b.points) for i, a in enumerate(loops) for b in loops[i + 1 :]]
    say(
        f"  cube of the (1, 2, 3) family, side sqrt(14) lambda / 2 = {math.sqrt(14) * math.pi / K_DRIVE * 100:.1f} cm, "
        f"six TE_z modes, seed {CUBE_SEED}: {len(loops)} closed loops ("
        + ", ".join(f"{c.length:.3f}" for c in loops)
        + f" sides long, {', '.join(str(len(folds(c, cube))) for c in loops)} folds), {len(curves) - len(loops)} lines "
        f"from wall to wall; linking numbers {', '.join(f'{v:+.1e}' for v in links) or 'none'}"
    )
    census: list[int] = []
    if not quick:
        census = [sum(c.closed for c in cube_curves(seed, True)[1]) for seed in range(12)]
        say(f"  closed loops in the fields of seeds 0 to 11: {census}")
    return {
        "sphere": sphere,
        "on_C": max(on_curves),
        "cube": {
            "seed": CUBE_SEED,
            "side_cm": math.sqrt(14) * math.pi / K_DRIVE * 100,
            "loops": len(loops),
            "lines": len(curves) - len(loops),
            "lengths": [c.length for c in loops],
            "links": links,
            "census": census,
        },
        "_det": (xs, det, scale),
        "_cube": (cube, curves),
    }


def compute_charts(quick: bool) -> dict[str, Any]:
    say("\n== charts: the sphere's families and the cylinder's mode chart (figure 3) ==")
    families = []
    for degree in range(1, 18):
        for tm in (False, True):
            for n, ka in enumerate(sphere_zeros(degree, 8, tm), start=1):
                if ka < 20:
                    families.append((float(ka), "TM" if tm else "TE", degree, n))
    families.sort()
    pairs = sorted(
        (abs(b[0] - a[0]) / ((a[0] + b[0]) / 2), a, b)
        for i, a in enumerate(families)
        for b in families[i + 1 : i + 40]
        if (a[1], a[2]) != (b[1], b[2])
    )
    closest = pairs[0]
    widths = [
        closest[0] * sphere_q(family[0], family[2], family[0] / K_DRIVE, family[1] == "TM") for family in closest[1:]
    ]
    say(
        f"  {len(families)} sphere families below ka = 20; the closest of different type or degree: "
        f"{closest[1][1]}({closest[1][2]},{closest[1][3]}) at {closest[1][0]:.6f} and {closest[2][1]}({closest[2][2]},"
        f"{closest[2][3]}) at {closest[2][0]:.6f}, {closest[0]:.3e} apart: {min(widths):.0f} to {max(widths):.0f} "
        f"linewidths of copper at 2.45 GHz (radius {closest[1][0] / K_DRIVE * 100:.1f} cm)"
    )
    table = {}
    for degree in (1, 2, 3):
        for tm in (False, True):
            ka = float(sphere_zeros(degree, 1, tm)[0])
            radius = ka / K_DRIVE
            table[f"{'TM' if tm else 'TE'} l={degree}"] = {
                "ka": ka,
                "radius_cm": radius * 100,
                "modes": 2 * degree + 1,
                "Q": sphere_q(ka, degree, radius, tm),
            }
    say(
        "  at 2.45 GHz: "
        + "; ".join(f"{k} ka {v['ka']:.4f}, {v['radius_cm']:.2f} cm, Q {v['Q']:,.0f}" for k, v in table.items())
    )
    lines = []
    for m in range(0, 5):
        for n in range(1, 4):
            for tm in (False, True):
                x = float(bessel_zeros(m, n, derivative=not tm)[-1])
                for p in range(0 if tm else 1, 5):
                    if x**2 < 90:
                        lines.append({"kind": "TM" if tm else "TE", "m": m, "n": n, "p": p, "x": x})
    field, _ = worked_cylinder()
    crossing = {"a_over_d_sq": 1 / field.d**2, "ka_sq": field.modes[0].ka ** 2}
    say(
        f"  {len(lines)} cylinder modes drawn; TE113 and TE221 cross at (a/d)^2 = {crossing['a_over_d_sq']:.5f} "
        f"(d/a = {field.d:.4f}), (ka)^2 = {crossing['ka_sq']:.3f}"
    )
    return {
        "families": families,
        "closest": {"relative": closest[0], "a": closest[1], "b": closest[2], "linewidths": widths},
        "table": table,
        "cylinder": lines,
        "crossing": crossing,
    }


def torus_trefoil(samples: int) -> Arr:
    """A trefoil, the (2, 3) torus knot on a torus of radii 3 and 1, turned so that height runs across the torus's
    plane: in that position it has two maxima, its bridge number."""
    t = np.linspace(0, 2 * math.pi, samples, endpoint=False)
    r = 3 + np.cos(3 * t)
    out: Arr = np.c_[r * np.cos(2 * t), np.sin(3 * t), r * np.sin(2 * t)]
    return out


def hopf_link(samples: int) -> tuple[Arr, Arr]:
    """Two round circles that link once, each tilted from the vertical."""
    t = np.linspace(0, 2 * math.pi, samples, endpoint=False)
    one: Arr = np.c_[np.cos(t), np.sin(t) * math.sin(0.5), np.sin(t) * math.cos(0.5)]
    two: Arr = np.c_[1 + np.cos(t), np.sin(t) * math.cos(0.4), -np.sin(t) * math.sin(0.4)]
    return one, two


def height_turns(curve: Arr) -> list[int]:
    z = curve[:, 2]
    up = np.sign(np.roll(z, -1) - z)
    return [int(i) for i in np.flatnonzero(up != np.roll(up, -1))]


def compute_folds(quick: bool, worked: dict[str, Any]) -> dict[str, Any]:
    say("\n== folds: where a null curve turns in height, and the wall there (figure 4) ==")
    field: CylinderField = worked["_field"]
    loop: NullCurve = worked["_loops"][0]
    mean = worked["mean_square"]
    det_s = tangents(field, loop.points)[:, 2]
    places = [fold_point(loop, field, i) for i in folds(loop, field)]
    distances = np.logspace(-3.5, -0.7, 40)
    growth = []
    for place in places:
        t = tangents(field, place[None, :])[0]
        u = t / np.linalg.norm(t)
        v = np.cross(u, [0.0, 0.0, 1.0])
        v /= np.linalg.norm(v)
        w = np.cross(u, v)
        directions = {"across, horizontal": v, "across, the other normal": w, "along the tangent": u}
        curves, local = {}, {}
        for name, e in directions.items():
            curves[name] = u_value(field, place[None, :] + distances[:, None] * e[None, :]) / mean
            e1, e2 = (float(u_value(field, (place + r * e)[None, :])[0]) for r in (0.01, 0.02))
            local[name] = math.log(e2 / e1) / math.log(2)
        small = distances < 3e-3
        near = {
            name: float(np.polyfit(np.log(distances[small]), np.log(values[small]), 1)[0])
            for name, values in curves.items()
        }
        growth.append(
            {"where": place, "z_over_d": float(place[2] / field.d), "curves": curves, "exponents": local, "near": near}
        )
        say(
            f"  fold at z/d {place[2] / field.d:.3f}: |E|^2 grows from it with exponent "
            + ", ".join(f"{local[k]:.2f} {k}" for k in local)
            + " (from 0.01 to 0.02 a); "
            + ", ".join(f"{near[k]:.3f}" for k in near)
            + " below 0.003 a"
        )
    soft = worked["_gradients"][0][:, 1]
    stiff = worked["_gradients"][0][:, 0]
    at_folds = [
        float(np.linalg.svd(pair_jacobian(field, p[None, :]), compute_uv=False)[0, 1] / math.sqrt(mean)) for p in places
    ]
    say(
        f"  along loop 0 the softer normal gradient runs {soft.min():.2f} to {soft.max():.2f} E_rms/a "
        f"({soft.max() / soft.min():.2f} times), {', '.join(f'{v:.2f}' for v in at_folds)} at the folds; the stiffer "
        f"{stiff.min():.2f} to {stiff.max():.2f}"
    )
    trefoil = torus_trefoil(600)
    one, two = hopf_link(400)
    turns = {"trefoil": len(height_turns(trefoil)), "Hopf link": len(height_turns(one)) + len(height_turns(two))}
    hopf = linking_number(one, two)
    say(
        f"  height turns: a trefoil in two-bridge position {turns['trefoil']}, a Hopf link of round circles "
        f"{turns['Hopf link']} (its linking number {hopf:+.3f})"
    )
    return {
        "folds": [{"z_over_d": g["z_over_d"], "exponents": g["exponents"], "near": g["near"]} for g in growth],
        "soft_range": (float(soft.min()), float(soft.max())),
        "stiff_range": (float(stiff.min()), float(stiff.max())),
        "soft_at_folds": at_folds,
        "knot_turns": turns,
        "hopf_linking": hopf,
        "_det_s": det_s / np.abs(det_s).max(),
        "_growth": growth,
        "_distances": distances,
        "_knots": (trefoil, one, two),
    }


def section_values(field: CylinderField, xy: Arr, zs: Arr, mean: float) -> Arr:
    """|E|^2 / mean on the points xy (P, 2) of a line or plane times the heights zs: shape (P, len(zs))."""
    px = np.zeros((len(xy), len(zs)))
    py = np.zeros_like(px)
    for mode, weight in zip(field.modes, field.weights, strict=True):
        planar = planar_partials(mode, xy, [(1, 0), (0, 1)])
        axial = axial_partial(mode, zs, 0)
        px += weight * planar[(1, 0)][:, None] * axial[None, :]
        py += weight * planar[(0, 1)][:, None] * axial[None, :]
    out: Arr = (px**2 + py**2) / mean
    return out


def crossings(curve: NullCurve, axis: int, level: float) -> Arr:
    """Where a curve crosses the plane coordinate[axis] = level, by linear interpolation between its points."""
    p = np.vstack([curve.points, curve.points[:1]]) if curve.closed else curve.points
    c = p[:, axis] - level
    at = np.flatnonzero(np.sign(c[:-1]) * np.sign(c[1:]) < 0)
    s = c[at] / (c[at] - c[at + 1])
    out: Arr = p[at] + s[:, None] * (p[at + 1] - p[at])
    return out


def compute_sections(quick: bool, worked: dict[str, Any]) -> dict[str, Any]:
    say("\n== sections: the worked field in section, and the escape path (figure 5) ==")
    field: CylinderField = worked["_field"]
    mean, d = worked["mean_square"], field.d
    n = 121 if quick else 241
    xs = np.linspace(-1, 1, n)
    zs = np.linspace(0, d, round(n * d / 2))
    gx, gy = np.meshgrid(xs, xs, indexing="ij")
    mid = section_values(field, np.c_[gx.ravel(), gy.ravel()], np.array([d / 2]), mean)[:, 0].reshape(gx.shape)
    mid = np.where(np.hypot(gx, gy) <= 1, mid, np.nan)
    plane = section_values(field, np.c_[xs, np.zeros(n)], zs, mean)
    saddle = worked["saddles"][0]
    azimuth = math.atan2(saddle[1], saddle[0])
    meridian = section_values(field, np.c_[xs * math.cos(azimuth), xs * math.sin(azimuth)], zs, mean)
    value, _, hess = u_gradient_hessian(field, saddle[None, :])
    eig, vec = np.linalg.eigh(hess[0])
    down = vec[:, 0]
    step = 0.006 if quick else 0.003
    paths = [descent(field, saddle + sign * 0.01 * down, step, 1e-7 * mean) for sign in (1.0, -1.0)]
    loop_points = np.vstack([c.points for c in worked["_loops"]])
    ends = []
    for path in paths:
        end = path[-1]
        nearest = float(np.min(np.linalg.norm(loop_points - end, axis=1)))
        if min(end[2], d - end[2]) < 2 * step:
            ends.append("cap")
        elif nearest < 0.02:
            ends.append("loop")
        else:
            ends.append("null line" if float(u_value(field, end[None, :])[0]) < 1e-6 * mean else "wall")
    path_u = [u_value(field, p) / mean for p in paths]
    say(
        f"  the saddle at s = {value[0] / mean:.5f}, azimuth {math.degrees(azimuth):.1f} deg: its Hessian's "
        "eigenvalues "
        f"{', '.join(f'{v / mean:.2f}' for v in eig)} (of mean |E|^2 per a^2); descending either way it reaches "
        f"{' and '.join(ends)}, the path rising nowhere above the saddle ({max(float(v.max()) for v in path_u):.5f})"
    )
    mid_marks = {
        "loops": [crossings(c, 2, d / 2) for c in worked["_loops"]],
        "lines": [crossings(c, 2, d / 2) for c in worked["_opens"]],
    }
    return {
        "azimuth_deg": math.degrees(azimuth),
        "hessian": eig / mean,
        "path_ends": ends,
        "path_peak": max(float(v.max()) for v in path_u),
        "_mid": (xs, mid),
        "_plane": (xs, zs, plane),
        "_meridian": (xs, zs, meridian, azimuth),
        "_paths": paths,
        "_marks": mid_marks,
    }


def compute_power(quick: bool, worked: dict[str, Any]) -> dict[str, Any]:
    say("\n== power: depth against power at 2.45 GHz, and the thresholds (figure 6) ==")
    field: CylinderField = worked["_field"]
    mean, radius = worked["mean_square"], worked["_radius"]
    boards = pt.FacingBoards(1e-3)
    h = boards.distance
    routes = {
        "fields carried in from far away (RFG-2)": (pt.free_space_depth_ceiling_ev(1.0), math.inf),
        "the worked cavity's null loops": (worked["depth_ueV_per_W"] * 1e-6, math.inf),
    }
    for name, share, gamma in (
        ("facing boards, h = 0.5 mm, 10 kOhm (RFG-3)", boards.depth_share, boards.gradient_share),
        ("four rods, h = 0.5 mm, 10 kOhm (RFG-3)", 1.0, 1.0),
    ):
        per_watt = pt.electrode_depth_per_watt_ev(1e4, h, F_DRIVE, share)
        edge = share * (0.908 * OMEGA * h) ** 2 / (16 * QM * gamma**2)
        routes[name] = (per_watt, edge)
    for name, (per_watt, edge) in routes.items():
        say(
            f"  {name}: {sci(per_watt)} eV per watt; 25 meV takes {sci(0.025 / per_watt)} W, 1 eV {sci(1 / per_watt)} W"
            + (f"; the slope holds to {edge:.2f} eV (q = 0.908)" if math.isfinite(edge) else "")
        )
    phi = np.linspace(0, 2 * math.pi, 721)
    z = np.linspace(0, field.d, 241)
    wall = section_values(field, np.c_[np.cos(phi), np.sin(phi)], z, mean)
    wall_peak = float(wall.max())
    up_per_watt = well_depth_ev(worked["Q"]["mix"], wall_peak, 1.0, math.pi * radius**3 * field.d)
    emission = {
        "returning": (30 / (3.17 * up_per_watt), 100 / (3.17 * up_per_watt)),
        "crossing": (30 / (8 * up_per_watt), 100 / (8 * up_per_watt)),
    }
    peak_field_1w = worked["E_rms_1W"] * math.sqrt(worked["peak_over_mean"])
    air = (3e6 / peak_field_1w) ** 2
    torr = ((10e3 / peak_field_1w) ** 2, (30e3 / peak_field_1w) ** 2)
    say(
        f"  the side wall's |E|^2 peaks at {wall_peak:.3f} of the mean: {up_per_watt * 1e3:.1f} meV per watt of "
        f"ponderomotive energy, so copper's 30 to 100 eV crossover at {emission['returning'][0]:.0f} to "
        f"{emission['returning'][1]:.0f} W for returning electrons (3.17 U_p) and {emission['crossing'][0]:.0f} to "
        f"{emission['crossing'][1]:.0f} W for crossing ones (8 U_p); air at 3 MV/m peak at {air / 1e3:.1f} kW; 10 to "
        "30 "
        f"kV/m at the antinodes at {torr[0]:.2f} to {torr[1]:.2f} W"
    )
    return {
        "routes": {k: {"eV_per_W": v[0], "edge_eV": v[1]} for k, v in routes.items()},
        "wall_peak": wall_peak,
        "wall_meV_per_W": up_per_watt * 1e3,
        "emission_W": emission,
        "air_kW": air / 1e3,
        "few_torr_W": torr,
        "board_ratio": routes["facing boards, h = 0.5 mm, 10 kOhm (RFG-3)"][0] / (worked["depth_ueV_per_W"] * 1e-6),
        "rod_ratio": routes["four rods, h = 0.5 mm, 10 kOhm (RFG-3)"][0] / (worked["depth_ueV_per_W"] * 1e-6),
    }


def tube_floor(field: CylinderField, points: Arr, weights: tuple[complex, complex], mean: float) -> Arr:
    """The floor of the tube round a null curve when the two modes take complex amplitudes: at each point, the least
    |E|^2 / mean over the curve's normal plane there, by Gauss-Newton on the four real parts of E."""
    t = tangents(field, points)
    t /= np.linalg.norm(t, axis=1)[:, None]
    helper = np.array([0.3, -0.5, 0.81])
    n1 = np.cross(t, helper)
    n1 /= np.linalg.norm(n1, axis=1)[:, None]
    n2 = np.cross(t, n1)
    offset = np.zeros((len(points), 2))
    singles = [CylinderField((mode,), (1.0,)) for mode in field.modes]

    def complex_e(x: Arr) -> tuple[CArr, CArr]:
        e: CArr = np.zeros((len(x), 2), dtype=np.complex128)
        j: CArr = np.zeros((len(x), 2, 3), dtype=np.complex128)
        for single, weight in zip(singles, weights, strict=True):
            g = pair(single, x)
            jac = pair_jacobian(single, x)
            e = e + weight * np.stack([g[:, 1], -g[:, 0]], -1)
            j = j + weight * np.stack([jac[:, 1], -jac[:, 0]], -2)
        return e, j

    best = np.full(len(points), np.inf)
    for _ in range(16):
        x = points + offset[:, :1] * n1 + offset[:, 1:] * n2
        e, j = complex_e(x)
        best = np.minimum(best, (np.abs(e) ** 2).sum(axis=1) / mean)
        r = np.concatenate([e.real, e.imag], axis=1)
        jn = np.stack([np.einsum("pci,pi->pc", j, n1), np.einsum("pci,pi->pc", j, n2)], -1)
        jr = np.concatenate([jn.real, jn.imag], axis=1)
        lhs = np.einsum("pki,pkj->pij", jr, jr)
        lhs = lhs + 1e-6 * np.trace(lhs, axis1=1, axis2=2)[:, None, None] * np.eye(2)[None, :, :]
        rhs = np.einsum("pki,pk->pi", jr, r)
        step = np.linalg.solve(lhs, rhs[:, :, None])[:, :, 0]
        size = np.linalg.norm(step, axis=1)
        step *= np.minimum(1.0, 0.02 / np.maximum(size, 1e-300))[:, None]
        offset = offset - step
    x = points + offset[:, :1] * n1 + offset[:, 1:] * n2
    e, _ = complex_e(x)
    out: Arr = np.minimum(best, (np.abs(e) ** 2).sum(axis=1) / mean)
    return out


def compute_tolerance(quick: bool, worked: dict[str, Any]) -> dict[str, Any]:
    say("\n== tolerance: the tube's floor against the split of the two modes (figure 7) ==")
    field: CylinderField = worked["_field"]
    q1, q2, q_mix = worked["Q"]["TE113"], worked["Q"]["TE221"], worked["Q"]["mix"]
    loop: NullCurve = worked["_loops"][0]
    points = loop.points[:: 5 if quick else 4]
    s = worked["saddle_share"]
    slope = (3**2 - 1**2) * (math.pi / field.d) ** 2 / field.modes[0].ka ** 2
    radius_m = worked["_radius"]

    def at(split: float) -> dict[str, float]:
        delta = split / q_mix
        d1, d2 = 1 / (1 + 2j * q1 * (-delta / 2)), 1 / (1 + 2j * q2 * (delta / 2))
        common = np.exp(-1j * np.angle(d1 + d2))
        floor = tube_floor(field, points, (complex(d1 * common), complex(d2 * common)), worked["mean_square"])
        return {
            "split": split,
            "phase_deg": math.degrees(abs(np.angle(d1 / d2))),
            "floor_over_s": float(floor.max() / s),
            "radius_nm": delta / slope * radius_m * 1e9,
        }

    splits = np.logspace(math.log10(0.004), 0, 9 if quick else 25)
    curve = [at(float(v)) for v in splits]
    table = [at(v) for v in (0.5, 0.1, 0.05, 0.02)]
    logs = np.log([row["floor_over_s"] for row in curve])
    tenth = float(np.exp(np.interp(math.log(0.1), logs, np.log(splits))))
    say(f"  d ln(f1/f2) = {slope:.4f} d ln(a/d) (CAV-8)")
    for row in table:
        say(
            f"  split {row['split']:.2f} linewidth: relative mode phase {row['phase_deg']:.1f} deg, the tube's floor "
            f"{row['floor_over_s']:.3g} s, {row['radius_nm']:.0f} nm on the radius"
        )
    common = 2 * q1 * q2 * (0.5 / q_mix) / (q2 - q1)
    say(
        f"  a floor of a tenth of the depth takes a split under {tenth:.3f} linewidths, "
        f"{tenth / q_mix / slope * radius_m * 1e9:.0f} nm on the radius; driven where Q1 delta1 = Q2 delta2 a "
        f"half-linewidth split costs {1 + common**2:.2f} times the power"
    )
    return {
        "slope": slope,
        "table": table,
        "tenth_split": tenth,
        "tenth_nm": tenth / q_mix / slope * radius_m * 1e9,
        "shared_phase_power": 1 + common**2,
        "_curve": curve,
    }


# ---- figures ---------------------------------------------------------------------------------------------------------
SUPERSCRIPTS = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def sci(value: float, digits: int = 2) -> str:
    """A number for a title: plain between 0.001 and 9999, else m × 10ⁿ."""
    if value == 0 or not math.isfinite(value):
        return f"{value:g}"
    exponent = math.floor(math.log10(abs(value)))
    if -3 <= exponent <= 3:
        return f"{value:.{max(digits - 1 - exponent, 0)}f}"
    mantissa = f"{value / 10**exponent:.{digits - 1}f}".rstrip("0").rstrip(".")
    power = "10" + str(exponent).translate(SUPERSCRIPTS)
    return power if mantissa == "1" else f"{mantissa} × {power}"


def _pyplot() -> Any:
    """matplotlib's pyplot on the Agg backend, in the house style of experiments/virtual_antenna_study.py."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from cycler import cycler

    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 9.5,
            "axes.labelsize": 9.5,
            "figure.titlesize": 11,
            "legend.fontsize": 8.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "axes.edgecolor": EDGE,
            "axes.labelcolor": INK,
            "axes.titlecolor": INK,
            "text.color": INK,
            "xtick.color": INK2,
            "ytick.color": INK2,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.axisbelow": True,
            "axes.titlelocation": "left",
            "lines.linewidth": 1.8,
            "lines.markersize": 5,
            "legend.frameon": False,
            "axes.prop_cycle": cycler(color=list(SERIES)),
            "savefig.dpi": 140,
            "mathtext.default": "regular",
        }
    )
    return plt


def _ramp(reverse: bool = False) -> Any:
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("ramp", list(BLUES[::-1] if reverse else BLUES))


def _heat() -> Any:
    """A sequential ramp for log |E|^2: dark where the field vanishes, light where it is strong."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("heat", ["#0d366b", "#256abf", "#6da7ec", "#cde2fb", "#fbe3d6", "#eb6834"])


def _height_ramp() -> Any:
    """A sequential ramp for height that stays visible at its low end."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("height", list(BLUES[2:]))


def _diverging() -> Any:
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("split", [SERIES[0], "#cde2fb", "#ecebe6", "#fbe3d6", SERIES[1]])


def _suptitle(fig: Any, text: str, width: int = 150) -> None:
    """A figure's title, made from its data, wrapped so that it never sets the figure's width."""
    import re
    import textwrap

    glued = re.sub(r"(\d) (?=\S)", "\\1\u00a0", text)
    fig.suptitle(textwrap.fill(glued, width), x=0.01, ha="left", fontsize=10.5)


def _save(fig: Any, out: Path, name: str) -> None:
    import matplotlib.pyplot as plt

    fig.savefig(out / f"{name}.png", bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    say(f"  wrote {name}.png")


def _style_3d(ax: Any, box: tuple[float, float, float], labels: tuple[str, str, str] = ("x", "y", "z")) -> None:
    ax.set_box_aspect(box)
    for set_label, label in zip((ax.set_xlabel, ax.set_ylabel, ax.set_zlabel), labels, strict=True):
        set_label(label, fontsize=8, labelpad=-2)
    ax.tick_params(labelsize=7, pad=-1)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor(SURFACE)
        axis.pane.set_edgecolor(GRID)
        axis._axinfo["grid"]["color"] = GRID


def _sphere_wire(ax: Any, radius: float, color: str = EDGE, alpha: float = 0.7) -> None:
    t = np.linspace(0, 2 * math.pi, 120)
    for lat in np.linspace(-60, 60, 5):
        r, z = radius * math.cos(math.radians(lat)), radius * math.sin(math.radians(lat))
        ax.plot(r * np.cos(t), r * np.sin(t), z + 0 * t, color=color, lw=0.6, alpha=alpha)
    for lon in np.linspace(0, 180, 6, endpoint=False):
        c, s = math.cos(math.radians(lon)), math.sin(math.radians(lon))
        ax.plot(radius * np.sin(t) * c, radius * np.sin(t) * s, radius * np.cos(t), color=color, lw=0.6, alpha=alpha)


def _cylinder_wire(ax: Any, radius: float, height: float) -> None:
    t = np.linspace(0, 2 * math.pi, 120)
    for z in (0.0, height):
        ax.plot(radius * np.cos(t), radius * np.sin(t), z + 0 * t, color=EDGE, lw=0.9)
    for angle in np.linspace(0, 2 * math.pi, 12, endpoint=False):
        ax.plot([radius * math.cos(angle)] * 2, [radius * math.sin(angle)] * 2, [0, height], color=GRID, lw=0.6)


def diameters(points: Sequence[CriticalPoint]) -> list[tuple[Arr, str]]:
    """One (direction, kind) per null diameter: a critical point and its antipode share a diameter."""
    out: list[tuple[Arr, str]] = []
    for point in points:
        if any(abs(abs(float(point.where @ seen)) - 1) < 1e-6 for seen, _ in out):
            continue
        out.append((point.where, point.kind))
    return out


LOOP_COLOURS = (SERIES[0], SERIES[2], SERIES[6])


def figure_models(res: dict[str, Any], worked: dict[str, Any], out: Path) -> None:
    """Figure 1: the owner's sphere with feeds, pickups and null diameters; the worked cylinder with its feeds, null
    loops, open null lines, folds and escape saddles."""
    plt = _pyplot()
    fig = plt.figure(figsize=(15.5, 7.4), layout="constrained")
    sphere = res["sphere"]
    radius = sphere["radius_cm"]
    left: Any = fig.add_subplot(1, 2, 1, projection="3d")
    _sphere_wire(left, radius)
    kinds_seen: set[str] = set()
    for direction, kind in diameters(sphere["_points"]):
        colour = MUTED if kind == "saddle" else SERIES[0]
        label = "a diameter over a saddle of A" if kind == "saddle" else "a diameter over a maximum and a minimum of A"
        ends = np.array([-direction, direction]) * radius
        left.plot(
            *ends.T,
            color=colour,
            lw=2.0,
            label=None if kind in kinds_seen or (kind != "saddle" and "x" in kinds_seen) else label,
        )
        kinds_seen.add("x" if kind != "saddle" else kind)
        left.scatter(*ends.T, color=colour, s=10)
    feeds = sphere["_feeds"]
    places, moments = feeds["feeds"] * radius, feeds["moments"]
    left.scatter(*places.T, marker="^", s=70, color=SERIES[1], depthshade=False, label="seven loop feeds")
    left.quiver(*places.T, *(2.6 * moments).T, color=SERIES[1], lw=1.4, arrow_length_ratio=0.3)
    picks, pick_moments = feeds["pickups"] * radius, feeds["pickup_moments"]
    left.scatter(*picks.T, marker="s", s=48, color=SERIES[3], depthshade=False, label="three pickups")
    left.quiver(*picks.T, *(2.6 * pick_moments).T, color=SERIES[3], lw=1.2, arrow_length_ratio=0.3)
    left.scatter([0], [0], [0], color=INK, s=14)
    _style_3d(left, (1, 1, 1), ("x (cm)", "y (cm)", "z (cm)"))
    left.view_init(elev=18, azim=-58)
    left.legend(loc="upper left", fontsize=8, bbox_to_anchor=(0.0, 0.98))
    c = sphere["critical"]
    left.set_title(
        f"the owner's sphere, radius {radius:.2f} cm: TE l = 3 at 2.45 GHz, {sphere['modes']} modes, Q "
        f"{sphere['Q']:,.0f}.\nOne real field's E vanishes on {sphere['diameters']} diameters through the centre, over "
        f"the {c['maximum']} maxima, {c['saddle']} saddles\nand {c['minimum']} minima of A; arrows: the loops' "
        f"moments, the feeds' coupling matrix of condition {sphere['condition']:.2f}",
        fontsize=9,
    )
    right: Any = fig.add_subplot(1, 2, 2, projection="3d")
    a, height = worked["radius_cm"], worked["height_cm"]
    _cylinder_wire(right, a, height)
    field: CylinderField = worked["_field"]
    for line in worked["_opens"]:
        right.plot(*(line.points * a).T, color=INK2, lw=1.3, ls=(0, (4, 2)))
    right.plot([], [], color=INK2, lw=1.3, ls=(0, (4, 2)), label="open null lines, wall to wall and cap to cap")
    for i, loop in enumerate(worked["_loops"]):
        p = np.vstack([loop.points, loop.points[:1]]) * a
        right.plot(*p.T, color=LOOP_COLOURS[i % 3], lw=2.4, label=f"loop {i}: {worked['loops'][i]['length_cm']:.2f} cm")
        for j in folds(loop, field):
            place = fold_point(loop, field, j) * a
            right.scatter(*place, color=INK, s=22, depthshade=False)
    right.scatter([], [], color=INK, s=22, label="folds, where det S = 0")
    saddles = np.array(worked["saddles"]) * a
    right.scatter(
        *saddles.T,
        marker="X",
        s=70,
        color=SERIES[7],
        depthshade=False,
        label=f"escape saddles, s = {worked['saddle_share']:.4f}",
    )
    t = np.linspace(0, 2 * math.pi, 40)
    for sign in (1.0, -1.0):
        centre = np.array([sign * (a - 1.2), 0.0, height / 2])
        right.plot(centre[0] + 1.0 * np.cos(t), 1.0 * np.sin(t), centre[2] + 0 * t, color=COPPER, lw=2.2)
        right.text(sign * (a + 2.5), 0, height / 2 + 1.5, "feed +" if sign > 0 else "feed −", fontsize=8, color=INK2)
    right.plot([], [], color=COPPER, lw=2.2, label="two loop feeds at (±a, 0, d/2)")
    _style_3d(right, (2 * a, 2 * a, height), ("x (cm)", "y (cm)", "z (cm)"))
    right.view_init(elev=20, azim=-48)
    right.legend(loc="upper left", fontsize=8, bbox_to_anchor=(0.0, 1.0))
    right.set_title(
        f"the worked cylinder, a = {a:.2f} cm, d = {height:.2f} cm: TE113 + TE221 at ka = {worked['ka']:.4f}, Q "
        f"{worked['Q']['mix']:,.0f}.\n{len(worked['_loops'])} closed null loops, {len(worked['_opens'])} open null "
        f"lines; {worked['depth_ueV_per_W']:.0f} µeV per watt deep, {worked['ceiling_ratio']:.2f} times RFG-2's "
        "ceiling.\nThe feeds, driven as sum and difference, excite TE221 and TE113 alone, their sin partners not at "
        "all",
        fontsize=9,
    )
    _suptitle(
        fig,
        "Three-dimensional models: the owner's sphere with its transceivers on the wall, and the worked cylinder whose "
        f"real field vanishes on closed loops {min(r['wall_cm'] for r in worked['loops']):.2f} cm or more from every "
        "wall",
    )
    _save(fig, out, "models")


def _shell(ax: Any, radius: float, colour: str) -> None:
    """A sphere with the octant x, y, z > 0 cut away."""
    u = np.linspace(0, 2 * math.pi, 73)
    v = np.linspace(0, math.pi, 37)
    x = radius * np.outer(np.cos(u), np.sin(v))
    y = radius * np.outer(np.sin(u), np.sin(v))
    z = radius * np.outer(np.ones_like(u), np.cos(v))
    cut = (x > 1e-9) & (y > 1e-9) & (z > 1e-9)
    ax.plot_surface(
        np.where(cut, np.nan, x),
        np.where(cut, np.nan, y),
        np.where(cut, np.nan, z),
        color=colour,
        alpha=0.32,
        lw=0,
        shade=True,
    )


def figure_nullsets(res: dict[str, Any], worked: dict[str, Any], out: Path) -> None:
    """Figure 2: one sphere family's null set, the two-profile cylinder's loops over the curve C, a cube field."""
    plt = _pyplot()
    import contourpy
    from mpl_toolkits.mplot3d.art3d import Line3DCollection  # type: ignore[import-untyped]

    fig = plt.figure(figsize=(16.5, 6.4), layout="constrained")
    sphere = res["sphere"]
    radius = sphere["radius_cm"]
    a1: Any = fig.add_subplot(1, 3, 1, projection="3d")
    _sphere_wire(a1, radius, alpha=0.5)
    for shell in sphere["shells"]:
        _shell(a1, float(shell) * radius, BLUES[2])
    for direction, kind in diameters(sphere["_points"]):
        ends = np.array([-direction, direction]) * radius
        a1.plot(*ends.T, color=MUTED if kind == "saddle" else SERIES[0], lw=2.2)
    a1.scatter([0], [0], [0], color=INK, s=12)
    _style_3d(a1, (1, 1, 1), ("x (cm)", "y (cm)", "z (cm)"))
    a1.view_init(elev=22, azim=35)
    a1.set_title(
        f"a sphere family, TE l = {sphere['degree']}, n = {sphere['family']}:\nradius {radius:.2f} cm, "
        f"{sphere['modes']} modes. One real field vanishes\non {sphere['diameters']} diameters (blue over extrema, "
        "grey "
        f"over\nsaddles) and on the zero-field shell at r = {float(sphere['shells'][0]):.4f} a\n(an octant cut away)",
        fontsize=9,
    )
    a2: Any = fig.add_subplot(1, 3, 2, projection="3d")
    xs, det, _ = res["_det"]
    field: CylinderField = worked["_field"]
    d = field.d
    generator = contourpy.contour_generator(xs, xs, np.ma.masked_invalid(det.T))
    paths = [np.asarray(p) for p in generator.lines(0.0) if len(p) > 4]
    zz = np.linspace(0, d, 2)
    for path in paths:
        a2.plot(path[:, 0], path[:, 1], 0 * path[:, 0], color=SERIES[3], lw=2.4)
        a2.plot_surface(
            np.repeat(path[:, :1], 2, axis=1),
            np.repeat(path[:, 1:2], 2, axis=1),
            np.tile(zz, (len(path), 1)),
            color=SERIES[3],
            alpha=0.12,
            lw=0,
            shade=False,
        )
    a2.plot([], [], color=SERIES[3], lw=1.6, label="C: grad Φ1 ∥ grad Φ2, and the vertical surfaces over it")
    t = np.linspace(0, 2 * math.pi, 120)
    for z in (0.0, d):
        a2.plot(np.cos(t), np.sin(t), z + 0 * t, color=EDGE, lw=0.8)
    for line in worked["_opens"]:
        a2.plot(*line.points.T, color=INK2, lw=1.2, ls=(0, (4, 2)))
    for i, loop in enumerate(worked["_loops"]):
        p = np.vstack([loop.points, loop.points[:1]])
        a2.plot(*p.T, color=LOOP_COLOURS[i % 3], lw=2.4)
        a2.plot(p[:, 0], p[:, 1], 0 * p[:, 0], color=LOOP_COLOURS[i % 3], lw=1.0, ls=":")
    a2.plot([], [], color=INK2, lw=1.0, ls=":", label="the loops' shadows on the floor: on C")
    _style_3d(a2, (2, 2, d), ("x/a", "y/a", "z/a"))
    a2.view_init(elev=38, azim=-52)
    a2.legend(loc="upper left", fontsize=7.6)
    a2.set_title(
        "the two-profile cylinder (TE113 + TE221):\nevery null lies over the curve C, |det| / scale ≤ "
        f"{sci(res['on_C'])}\non the traced curves; its closed loops lie on\nuntwisted annuli, so unlinked (CAV-3)",
        fontsize=9,
    )
    a3: Any = fig.add_subplot(1, 3, 3, projection="3d")
    _, curves = res["_cube"]
    side = res["cube"]["side_cm"]
    for corner in ((0, 0), (0, 1), (1, 0), (1, 1)):
        a3.plot([corner[0] * side] * 2, [corner[1] * side] * 2, [0, side], color=EDGE, lw=0.8)
        for z in (0.0, side):
            a3.plot([0, side, side, 0, 0], [0, 0, side, side, 0], [z] * 5, color=EDGE, lw=0.8)
    norm = plt.Normalize(0, side)
    cmap = _height_ramp()
    for curve in curves:
        p = (np.vstack([curve.points, curve.points[:1]]) if curve.closed else curve.points) * side
        segments = np.stack([p[:-1], p[1:]], axis=1)
        lines = Line3DCollection(segments, cmap=cmap, norm=norm, linewidths=3.0 if curve.closed else 1.6)
        lines.set_array(0.5 * (p[:-1, 2] + p[1:, 2]))
        a3.add_collection3d(lines)
    _style_3d(a3, (1, 1, 1), ("x (cm)", "y (cm)", "z (cm)"))
    a3.set_xlim(0, side)
    a3.set_ylim(0, side)
    a3.set_zlim(0, side)
    a3.view_init(elev=20, azim=-55)
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=a3, shrink=0.55, pad=0.08)
    bar.set_label("height (cm)")
    cube_res = res["cube"]
    a3.set_title(
        f"a random real field of the cube's (1, 2, 3)\nfamily: side {side:.1f} cm, six TE_z modes, seed "
        f"{cube_res['seed']};\n{cube_res['loops']} closed loops (thick) and {cube_res['lines']} lines from wall to "
        "wall,\n"
        "coloured by height; linking number " + (", ".join(f"{round(v):d}" for v in cube_res["links"]) or "none"),
        fontsize=9,
    )
    _suptitle(
        fig,
        "Three-dimensional null sets: rays in a sphere's family, unlinked loops over C in two profiles, three "
        "profiles in a cube",
    )
    _save(fig, out, "nullsets")


def figure_charts(res: dict[str, Any], out: Path) -> None:
    """Figure 3: the sphere's TE and TM roots against l, with degeneracy and radius; the cylinder's mode chart."""
    plt = _pyplot()
    fig, (left, right) = plt.subplots(1, 2, figsize=(14.5, 5.8), layout="constrained")
    for ka, kind, degree, _family in res["families"]:
        if degree > 9 or ka > 14:
            continue
        left.scatter(
            degree + (-0.12 if kind == "TE" else 0.12),
            ka,
            s=9 * (2 * degree + 1),
            marker="o" if kind == "TE" else "s",
            color=SERIES[0] if kind == "TE" else SERIES[1],
            alpha=0.85,
            edgecolor=SURFACE,
            lw=0.5,
        )
    left.scatter([], [], marker="o", color=SERIES[0], label="TE, j_l(ka) = 0")
    left.scatter([], [], marker="s", color=SERIES[1], label="TM, [x j_l(x)]\u2032 = 0")
    for name, row in res["table"].items():
        left.annotate(
            f"{row['radius_cm']:.2f} cm",
            (int(name[-1]) + (-0.12 if name.startswith("TE") else 0.12), row["ka"]),
            xytext=(-30 if name.startswith("TE") else 8, 4),
            textcoords="offset points",
            fontsize=7.4,
            color=INK2,
        )
    secondary = left.secondary_yaxis("right", functions=(lambda v: v / K_DRIVE * 100, lambda r: r * K_DRIVE / 100))
    secondary.set_ylabel("radius at 2.45 GHz (cm)")
    left.set_xlabel("degree l (marker area ∝ 2l + 1, the modes sharing the frequency)")
    left.set_ylabel("ka")
    left.set_xticks(range(1, 10))
    left.legend(loc="upper left")
    close = res["closest"]
    left.set_title(
        f"the sphere's families: 2l + 1 modes each, and no two of different type or degree share a frequency\n"
        f"(the closest below ka = 20, {close['a'][1]}({close['a'][2]},{close['a'][3]}) and {close['b'][1]}"
        f"({close['b'][2]},{close['b'][3]}), are {sci(close['relative'], 3)} apart)",
        fontsize=9,
    )
    aspect = np.linspace(0, 1.25, 200)
    for line in res["cylinder"]:
        y = line["x"] ** 2 + (line["p"] * math.pi) ** 2 * aspect
        right.plot(
            aspect, y, color=GRID if line["kind"] == "TM" else EDGE, lw=0.8, ls="--" if line["kind"] == "TM" else "-"
        )
    for (m, n, p), colour in (((1, 1, 3), SERIES[0]), ((2, 2, 1), SERIES[1])):
        x = float(bessel_zeros(m, n, derivative=True)[-1])
        right.plot(
            aspect,
            x**2 + (p * math.pi) ** 2 * aspect,
            color=colour,
            lw=2.2,
            label=f"TE{m}{n}{p}: x\u2032{m}{n} = {x:.4f}",
        )
    cross = res["crossing"]
    right.scatter([cross["a_over_d_sq"]], [cross["ka_sq"]], color=INK, s=36, zorder=5)
    right.annotate(
        f"d/a = {1 / math.sqrt(cross['a_over_d_sq']):.4f}, ka = {math.sqrt(cross['ka_sq']):.4f}",
        (cross["a_over_d_sq"], cross["ka_sq"]),
        xytext=(12, -14),
        textcoords="offset points",
        fontsize=8,
    )
    right.plot([], [], color=EDGE, lw=0.8, label="other TE modes")
    right.plot([], [], color=GRID, lw=0.8, ls="--", label="TM modes")
    right.set_xlim(0, 1.25)
    right.set_ylim(0, 90)
    right.set_xlabel("(a/d)²")
    right.set_ylabel("(ka)² = x² + (pπ)² (a/d)²")
    right.legend(loc="upper left", fontsize=8)
    right.set_title(
        "the cylinder's mode chart: each mode a line of slope (pπ)², so one aspect ratio\nbrings one pair of families "
        "together; TE113 and TE221 at the dot",
        fontsize=9,
    )
    _suptitle(fig, "Mode charts: the sphere's degeneracy is SO(3)'s 2l + 1, the cylinder's coincidences are set by d/a")
    _save(fig, out, "charts")


def figure_folds(res: dict[str, Any], worked: dict[str, Any], out: Path) -> None:
    """Figure 4: a loop coloured by det S, the wall's growth at a fold, the normal gradients, a Hopf link and a
    trefoil against height."""
    plt = _pyplot()
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    fig = plt.figure(figsize=(14.5, 11.0), layout="constrained")
    grid = fig.add_gridspec(2, 2)
    a1: Any = fig.add_subplot(grid[0, 0], projection="3d")
    field: CylinderField = worked["_field"]
    loop: NullCurve = worked["_loops"][0]
    a = worked["radius_cm"]
    p = np.vstack([loop.points, loop.points[:1]]) * a
    shade = np.r_[res["_det_s"], res["_det_s"][:1]]
    segments = np.stack([p[:-1], p[1:]], axis=1)
    norm = plt.Normalize(-1, 1)
    lines = Line3DCollection(segments, cmap=_diverging(), norm=norm, linewidths=3.2)
    lines.set_array(0.5 * (shade[:-1] + shade[1:]))
    a1.add_collection3d(lines)
    for g in res["_growth"]:
        a1.scatter(*(g["where"] * a), color=INK, s=40, depthshade=False)
        a1.text(*(g["where"] * a + np.array([0.8, 0.0, 0.6])), f"fold, z/d = {g['z_over_d']:.3f}", fontsize=8)
    centre, half = 0.5 * (p.min(axis=0) + p.max(axis=0)), 0.5 * float((p.max(axis=0) - p.min(axis=0)).max()) + 1
    a1.set_xlim(centre[0] - half, centre[0] + half)
    a1.set_ylim(centre[1] - half, centre[1] + half)
    a1.set_zlim(centre[2] - half, centre[2] + half)
    _style_3d(a1, (1, 1, 1), ("x (cm)", "y (cm)", "z (cm)"))
    a1.view_init(elev=12, azim=-62)
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=_diverging()), ax=a1, shrink=0.55, pad=0.06)
    bar.set_label("det S / max |det S|: rising where > 0")
    a1.set_title(
        f"loop 0 of the worked example ({worked['loops'][0]['length_cm']:.2f} cm), coloured by det S, the z-part of\n"
        "grad Πx × grad Πy: it rises on extremum strands (det S > 0) and falls on saddle strands, turning at two folds",
        fontsize=9,
    )
    a2 = fig.add_subplot(grid[0, 1])
    distances = res["_distances"]
    first = res["_growth"][0]
    for (name, values), colour in zip(first["curves"].items(), (SERIES[0], SERIES[2], SERIES[1]), strict=True):
        a2.loglog(
            distances,
            values,
            color=colour,
            label=f"{name}: exponent {first['near'][name]:.2f} near the fold, "
            f"{first['exponents'][name]:.2f} at 0.01-0.02 a",
        )
    base = first["curves"]["across, horizontal"][0]
    a2.loglog(distances, base * (distances / distances[0]) ** 2, color=MUTED, ls=":", lw=1.2, label="slope 2")
    tangent0 = first["curves"]["along the tangent"]
    a2.loglog(distances, tangent0[10] * (distances / distances[10]) ** 4, color=MUTED, ls="--", lw=1.2, label="slope 4")
    a2.set_xlabel("distance from the fold (units of a)")
    a2.set_ylabel("|E|² / mean |E|²")
    a2.legend(loc="upper left", fontsize=8)
    nears = [g["near"] for g in res["_growth"]]
    a2.set_title(
        f"the wall at the fold z/d = {first['z_over_d']:.3f}: quadratic across the curve, quartic along its tangent "
        "(CAV-6);\nbelow 0.003 a from the two folds the exponents are "
        + "; ".join(", ".join(f"{e[k]:.2f}" for k in e) for e in nears),
        fontsize=9,
    )
    a3 = fig.add_subplot(grid[1, 0])
    arc = loop.arc * a
    gradients = worked["_gradients"][0]
    a3.plot(arc, gradients[:, 0], color=SERIES[0], label="the stiffer normal direction")
    a3.plot(arc, gradients[:, 1], color=SERIES[1], label="the softer normal direction")
    for i in folds(loop, field):
        a3.axvline(arc[i], color=INK, ls=":", lw=1.0)
    a3.set_ylim(0, None)
    a3.set_xlabel("arc length along loop 0 (cm)")
    a3.set_ylabel("field gradient across the null (E_rms / a)")
    a3.legend(loc="upper right", fontsize=8)
    lo_s, hi_s = res["soft_range"]
    a3.set_title(
        f"the two normal gradients along the loop: the softer runs {lo_s:.2f} to {hi_s:.2f} E_rms/a "
        f"({hi_s / lo_s:.1f} times),\nlowest at the folds (dotted), and never zero: soft, not open",
        fontsize=9,
    )
    a4: Any = fig.add_subplot(grid[1, 1], projection="3d")
    trefoil, one, two = res["_knots"]
    curves_k = [
        trefoil * 0.9 + np.array([-4.5, 0, 0]),
        one * 2.0 + np.array([4.0, 0, 0]),
        two * 2.0 + np.array([4.0, 0, 0]),
    ]
    allz = np.concatenate([c[:, 2] for c in curves_k])
    hnorm = plt.Normalize(float(allz.min()), float(allz.max()))
    for curve in curves_k:
        c = np.vstack([curve, curve[:1]])
        seg = np.stack([c[:-1], c[1:]], axis=1)
        coll = Line3DCollection(seg, cmap=_height_ramp(), norm=hnorm, linewidths=3.0)
        coll.set_array(0.5 * (c[:-1, 2] + c[1:, 2]))
        a4.add_collection3d(coll)
        for i in height_turns(curve):
            a4.scatter(*curve[(i + 1) % len(curve)], color=SERIES[7], s=36, depthshade=False)
    a4.set_xlim(-8.5, 8)
    a4.set_ylim(-4, 4)
    a4.set_zlim(-3.6, 3.6)
    _style_3d(a4, (16.5, 8, 10), ("x", "y", "height"))
    a4.view_init(elev=14, azim=-72)
    turns = res["knot_turns"]
    a4.set_title(
        f"a trefoil (left) and a Hopf link (right) against height, their {turns['trefoil']} and {turns['Hopf link']} "
        "folds in red:\nat least 2b(K) = 4 for either, Schubert's bridge number b = 2 (CAV-5)",
        fontsize=9,
    )
    _suptitle(fig, "Folds and the wall: knotting forces folds, and a fold softens the wall without opening it")
    _save(fig, out, "folds")


def figure_sections(res: dict[str, Any], worked: dict[str, Any], out: Path) -> None:
    """Figure 5: log |E|^2 of the worked field on the mid-plane, the plane y = 0 and the meridian of the escape
    saddle, with the loops' crossings, the open lines, the saddles and the escape path."""
    plt = _pyplot()
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 6.6), layout="constrained", width_ratios=[1.0, 0.82, 0.82])
    field: CylinderField = worked["_field"]
    d = field.d
    cmap, vmin, vmax = _heat(), -3.0, 1.0
    xs, mid = res["_mid"]
    image = axes[0].pcolormesh(
        xs, xs, np.log10(np.maximum(mid.T, 1e-12)), cmap=cmap, vmin=vmin, vmax=vmax, shading="auto"
    )
    for i, marks in enumerate(res["_marks"]["loops"]):
        axes[0].scatter(
            marks[:, 0],
            marks[:, 1],
            color=LOOP_COLOURS[i % 3],
            s=54,
            edgecolor=INK,
            lw=0.9,
            zorder=4,
            label=f"loop {i} crosses",
        )
    lines = (
        np.vstack([m for m in res["_marks"]["lines"] if len(m)])
        if any(len(m) for m in res["_marks"]["lines"])
        else np.zeros((0, 3))
    )
    axes[0].scatter(
        lines[:, 0],
        lines[:, 1],
        color=INK,
        marker="D",
        s=26,
        edgecolor=SURFACE,
        lw=0.6,
        zorder=4,
        label="open lines cross",
    )
    t = np.linspace(0, 2 * math.pi, 200)
    axes[0].plot(np.cos(t), np.sin(t), color=INK2, lw=1.0)
    axes[0].set_aspect("equal")
    axes[0].set_xlabel("x / a")
    axes[0].set_ylabel("y / a")
    axes[0].legend(loc="lower left", fontsize=7.4, ncol=2)
    axes[0].set_title("the mid-plane z = d/2: every null crossing it is a zero of the field", fontsize=9)
    xs2, zs, plane = res["_plane"]
    axes[1].pcolormesh(
        xs2, zs / d, np.log10(np.maximum(plane.T, 1e-12)), cmap=cmap, vmin=vmin, vmax=vmax, shading="auto"
    )
    for i, loop in enumerate(worked["_loops"]):
        if np.abs(loop.points[:, 1]).max() < 1e-9:
            p = np.vstack([loop.points, loop.points[:1]])
            axes[1].plot(p[:, 0], p[:, 2] / d, color=LOOP_COLOURS[i % 3], lw=2.4, label=f"loop {i}, in this plane")
    for line in worked["_opens"]:
        if np.abs(line.points[:, 1]).max() < 1e-9:
            axes[1].plot(
                line.points[:, 0],
                line.points[:, 2] / d,
                color=INK,
                lw=1.6,
                ls=(0, (4, 2)),
                label="the line from cap to cap",
            )
    axes[1].set_xlabel("x / a (y = 0)")
    axes[1].set_ylabel("z / d")
    axes[1].legend(loc="lower left", fontsize=7.4)
    axes[1].set_title(
        "the plane y = 0, the mirror plane of the field:\nthe loop and the line that lie in it", fontsize=9
    )
    xs3, zs3, meridian, azimuth = res["_meridian"]
    axes[2].pcolormesh(
        xs3, zs3 / d, np.log10(np.maximum(meridian.T, 1e-12)), cmap=cmap, vmin=vmin, vmax=vmax, shading="auto"
    )
    normal = np.array([-math.sin(azimuth), math.cos(azimuth), 0.0])
    along = np.array([math.cos(azimuth), math.sin(azimuth), 0.0])
    for i, loop in enumerate(worked["_loops"]):
        p = np.vstack([loop.points, loop.points[:1]])
        side = p @ normal
        at = np.flatnonzero(np.sign(side[:-1]) * np.sign(side[1:]) < 0)
        frac = side[at] / (side[at] - side[at + 1])
        hits = p[at] + frac[:, None] * (p[at + 1] - p[at])
        axes[2].scatter(hits @ along, hits[:, 2] / d, color=LOOP_COLOURS[i % 3], s=54, edgecolor=INK, lw=0.9, zorder=4)
    for path, colour in zip(res["_paths"], (SERIES[3], SERIES[4]), strict=True):
        axes[2].plot(path @ along, path[:, 2] / d, color=colour, lw=2.0)
    axes[2].plot([], [], color=SERIES[3], lw=2.0, label=f"descent to the {res['path_ends'][0]}")
    axes[2].plot([], [], color=SERIES[4], lw=2.0, label=f"descent to the {res['path_ends'][1]}")
    saddles = np.array(worked["saddles"])
    keep = np.abs(saddles @ normal) < 1e-6
    axes[2].scatter(
        saddles[keep] @ along,
        saddles[keep, 2] / d,
        marker="X",
        s=80,
        color=SERIES[7],
        zorder=5,
        label=f"the escape saddles, s = {worked['saddle_share']:.4f}",
    )
    axes[2].set_xlabel(f"distance from the axis along azimuth {math.degrees(azimuth):.1f}° (units of a)")
    axes[2].set_ylabel("z / d")
    axes[2].legend(loc="center left", fontsize=7.4)
    axes[2].set_title(
        "the meridian through the escape saddles: the path from a loop\nover the saddle to the cap, projected on it",
        fontsize=9,
    )
    bar = fig.colorbar(image, ax=axes, shrink=0.8, pad=0.01)
    bar.set_label("log₁₀ |E|² / mean |E|²")
    _suptitle(
        fig,
        f"The worked example in section: the loops' well is {worked['depth_ueV_per_W']:.0f} µeV per watt deep, its "
        f"escape the index-1 saddle at {worked['saddle_share']:.4f} of the mean |E|², ρ = "
        f"{math.hypot(saddles[0, 0], saddles[0, 1]):.3f} a, z/d = {min(saddles[0, 2], d - saddles[0, 2]) / d:.3f} and "
        f"{max(saddles[0, 2], d - saddles[0, 2]) / d:.3f}, between the loops and the caps (where every TE field "
        "vanishes)",
    )
    _save(fig, out, "sections")


def figure_power(res: dict[str, Any], worked: dict[str, Any], out: Path) -> None:
    """Figure 6: depth against power at 2.45 GHz, with the thresholds."""
    plt = _pyplot()
    fig, ax = plt.subplots(figsize=(11.5, 6.6), layout="constrained")
    power = np.logspace(-3, 5, 300)
    styles = (
        (MUTED, "--", 1.4),
        (SERIES[0], "-", 2.6),
        (SERIES[1], "-", 2.0),
        (SERIES[2], "-", 2.0),
    )
    for (name, row), (colour, ls, lw) in zip(res["routes"].items(), styles, strict=True):
        depth = row["eV_per_W"] * power
        held = depth <= row["edge_eV"]
        ax.loglog(
            power[held], depth[held], color=colour, ls=ls, lw=lw, label=f"{name}: {sci(row['eV_per_W'])} eV per watt"
        )
        if (~held).any():
            ax.loglog(power[~held], depth[~held], color=colour, ls=":", lw=1.0, alpha=0.6)
            edge_power = row["edge_eV"] / row["eV_per_W"]
            ax.scatter([edge_power], [row["edge_eV"]], color=colour, s=22, zorder=4)
    for level, text in ((0.025, "25 meV, room temperature"), (1.0, "1 eV")):
        ax.axhline(level, color=INK2, lw=0.9, ls="-.")
        ax.text(1.3e-3, level * 1.15, text, fontsize=8, color=INK2)
    cavity = worked["depth_ueV_per_W"] * 1e-6
    for level in (0.025, 1.0):
        ax.scatter([level / cavity], [level], color=SERIES[0], s=40, zorder=5)
        ax.annotate(
            f"{sci(level / cavity)} W",
            (level / cavity, level),
            xytext=(6, -12),
            textcoords="offset points",
            fontsize=8,
            color=SERIES[0],
        )
    torr = res["few_torr_W"]
    ax.axvspan(*torr, color=SERIES[4], alpha=0.12, lw=0)
    ax.text(math.sqrt(torr[0] * torr[1]), 3e-6, "air near 3 Torr\nbreaks down", ha="center", fontsize=7.6, color=INK2)
    emission = res["emission_W"]
    ax.axvspan(emission["crossing"][0], emission["returning"][1], color=SERIES[3], alpha=0.14, lw=0)
    ax.text(
        math.sqrt(emission["crossing"][0] * emission["returning"][1]),
        3e-6,
        "wall electrons reach\ncopper's secondary-\nemission crossover",
        ha="center",
        fontsize=7.6,
        color=INK2,
    )
    ax.axvline(res["air_kW"] * 1e3, color=SERIES[7], lw=1.2)
    ax.text(res["air_kW"] * 1e3 * 1.15, 3e-6, "air at one\natmosphere\nbreaks down", fontsize=7.6, color=INK2)
    ax.set_xlim(1e-3, 1e5)
    ax.set_ylim(1e-6, 1e2)
    ax.set_xlabel("power into the cavity or the electrodes' resonator at 2.45 GHz (W)")
    ax.set_ylabel("depth of the well (eV)")
    ax.legend(loc="upper left", fontsize=8)
    ax.set_title(
        f"the cavity's loops are {worked['ceiling_ratio']:.2f} times the free-space ceiling; facing boards "
        f"{sci(res['board_ratio'])} and four rods {sci(res['rod_ratio'])} times the loops (dots: q = 0.908, where the "
        "electrodes' slope ends)",
        fontsize=9,
    )
    _suptitle(
        fig,
        f"Depth against power at 2.45 GHz: 25 meV takes {sci(0.025 / cavity)} W and 1 eV {sci(1 / cavity / 1e3)} kW in "
        f"the cavity; the side wall's field gives electrons {res['wall_meV_per_W']:.1f} meV per watt, so secondary "
        f"emission from {emission['crossing'][0]:.0f} W to {emission['returning'][1] / 1e3:.1f} kW (hypothesis)",
    )
    _save(fig, out, "power")


def figure_tolerance(res: dict[str, Any], worked: dict[str, Any], out: Path) -> None:
    """Figure 7: the tube's floor against the split in linewidths."""
    plt = _pyplot()
    fig, ax = plt.subplots(figsize=(10.5, 6.2), layout="constrained")
    curve = res["_curve"]
    splits = np.array([row["split"] for row in curve])
    floors = np.array([row["floor_over_s"] for row in curve])
    ax.loglog(
        splits,
        floors,
        color=SERIES[0],
        lw=2.2,
        label="the tube's floor: the highest along loop 0 of its least |E|² across",
    )
    for row in res["table"]:
        ax.scatter([row["split"]], [row["floor_over_s"]], color=SERIES[1], s=40, zorder=4)
        ax.annotate(
            f"{row['phase_deg']:.1f}°, {row['floor_over_s']:.3g} s",
            (row["split"], row["floor_over_s"]),
            xytext=(-78, 6),
            textcoords="offset points",
            fontsize=8,
            color=INK2,
        )
    ax.axhline(1.0, color=INK, lw=1.0, ls="--")
    ax.text(0.0045, 1.12, "the depth s: the floor reaches the escape saddle and the well is gone", fontsize=8)
    ax.axhline(0.1, color=INK2, lw=0.9, ls=":")
    ax.text(0.0045, 0.112, "a tenth of the depth", fontsize=8, color=INK2)
    ax.axvline(res["tenth_split"], color=INK2, lw=0.9, ls=":")
    ax.set_xlabel("split of TE113 and TE221, in linewidths of the mix (driven midway, real weights)")
    ax.set_ylabel("the tube's floor, in units of the depth s")
    slope, q_mix, radius = res["slope"], worked["Q"]["mix"], worked["_radius"]
    top = ax.secondary_xaxis(
        "top", functions=(lambda v: v / q_mix / slope * radius * 1e9, lambda nm: nm * q_mix * slope / (radius * 1e9))
    )
    top.set_xlabel("the same split as an error on the radius, d ln(f1/f2) = 0.829 d ln(a/d) (nm)")
    ax.legend(loc="upper left", fontsize=8)
    ax.set_title(
        f"the floor grows as the square of the relative mode phase: under a tenth of the depth below a split of "
        f"{res['tenth_split']:.3f} linewidths, {res['tenth_nm']:.0f} nm on the {worked['radius_cm']:.2f} cm radius",
        fontsize=9,
    )
    _suptitle(
        fig,
        "Tolerance: driven between split modes with real weights, the modes take different phases and the quadrature "
        "field lifts the whole tube; the common resonance the drive follows, the split it cannot",
    )
    _save(fig, out, "tolerance")


# ---- the run ---------------------------------------------------------------------------------------------------------
SECTIONS: dict[str, tuple[Callable[..., dict[str, Any]], Callable[..., None], bool]] = {
    "models": (compute_models, figure_models, True),
    "nullsets": (compute_nullsets, figure_nullsets, True),
    "charts": (lambda quick: compute_charts(quick), figure_charts, False),
    "folds": (compute_folds, figure_folds, True),
    "sections": (compute_sections, figure_sections, True),
    "power": (compute_power, figure_power, True),
    "tolerance": (compute_tolerance, figure_tolerance, True),
}
ORDER = tuple(SECTIONS)
KEEP = 64  # results.json keeps arrays and lists of at most this many numbers; the figures hold the rest


def slim(value: Any) -> Any:
    """The results without private entries (keys starting with '_') and large arrays, for results.json."""
    if isinstance(value, np.ndarray):
        return value if value.size <= KEEP else None
    if isinstance(value, dict):
        kept = {str(k): slim(v) for k, v in value.items() if not str(k).startswith("_")}
        return {k: v for k, v in kept.items() if v is not None}
    if isinstance(value, (list, tuple)):
        if len(value) > KEEP:
            return None
        return [slim(v) for v in value]
    if isinstance(value, CriticalPoint):
        return {"where": value.where, "value": value.value, "kind": value.kind}
    return value


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument("--quick", action="store_true", help="coarser grids and steps, into cavity_knots/quick/")
    parser.add_argument("--only", default="", help="comma-separated sections, of: " + ", ".join(ORDER))
    parser.add_argument("--out", type=Path, default=None, help="where the figures, results.json and run.log go")
    args = parser.parse_args(argv)
    chosen = [name for name in args.only.split(",") if name] or list(ORDER)
    unknown = [name for name in chosen if name not in ORDER]
    if unknown:
        print(f"no section {unknown}; there are {list(ORDER)}", file=sys.stderr)
        return 2
    out: Path = args.out if args.out is not None else (OUT / "quick" if args.quick else OUT)
    out.mkdir(parents=True, exist_ok=True)
    LOG.clear()
    began = time.time()
    say(f"cavity knots study: {'quick' if args.quick else 'full'} run of {', '.join(chosen)}")
    results: dict[str, Any] = {}
    if any(SECTIONS[name][2] for name in chosen):
        results["worked"] = compute_worked(args.quick)
    for name in chosen:
        start = time.time()
        compute, _, needs = SECTIONS[name]
        results[name] = compute(args.quick, results["worked"]) if needs else compute(args.quick)
        say(f"  [{name}: {time.time() - start:.1f} s]")
    say("\n== figures ==")
    for name in chosen:
        _, draw, needs = SECTIONS[name]
        if needs:
            draw(results[name], results["worked"], out)
        else:
            draw(results[name], out)
    say(f"\ndone in {time.time() - began:.0f} s")
    if chosen == list(ORDER):
        (out / "results.json").write_text(json.dumps(lab.plain(slim(results)), indent=1) + "\n")
        (out / "run.log").write_text("\n".join(LOG) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
