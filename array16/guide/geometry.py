"""The geometry of a guide's path: frames, plates and their flat outlines, as functions on ClosedCurve, and the paths
stage 2 sets against T915-S (docs/KNOTTED_GUIDE.md sections 3.5, 8.1 and 10; KNT-10).

Taken from experiments/knotted_guide_study.py's typed ports and held to its numbers: T915-S, the r/R = 0.45 trefoil of
317.8 mm, framed k = 0, twists +0.48 deg/mm, keeps 18.57 mm in its boards' plane and has plates 329.5 and 306.1 mm long
at +-h (h = 0.915 mm).

  scaled, torus            a curve at a length; the (p, q) torus knot at r/R and a length (the trefoil is (2, 3))
  curvature_vector         dT/ds by second differences, without its part along the tangent
  Framing, framing         board normals that close with k half turns against the twist-free frame, and the
                           curvature in the boards' plane (kappa_g) and toward them (kappa_n)
  plate, develop, outline  a plate's centre line at an offset along the board normal, its arc-length elements and
                           geodesic curvature; laid flat; the outline of a strip of a width along it
  smallest_box             the turn that fits plane points in their smallest box, and its sides
  ring, ellipse, helix     stage 2's paths for H1 (experiments/twin_rung1/REPORT.md section 6): the flat ring of a
                           length, the ellipse of a length whose tightest radius is given (T915-S's in-plane one), and
                           the (1, 3) torus curve at r/R
  flat                     a planar curve's framing with its bend wholly in the boards' plane ('in-plane') or wholly
                           toward them, on their stiff side ('toward')

A framing's `across` is the first normal SectionField.along takes (normals=framing.across), its v then running along
-normal; the twin bench sweeps T915-S that way (array16.twin.TwinGuide.sweep). framing here is the twin's, to the bit.

Hypothesis H1 (stage 1's report): planar paths keep their whole static limit and paths out of the plane about 0.86 of
it. Its two rivals read the same numbers, and these paths tell them apart at section B and 50 V: on the ellipse, planar
with its bend changing along it, the changing bend (H1a) loses like a knot, about 0.86, where H1 keeps the whole limit;
on the flat ring framed toward its boards, planar with its bend wholly on the stiff side, H1 holds its static limit,
about 4.5 times the flat ring's 32 eV, where the stiff side giving a bend nothing (H1b) holds about the flat ring's
own; and on the (1, 3) helix at r/R = 0.1 and 0.3 beside 0.2, H1 keeps 0.86 at each, where H1a moves with how fast the
bend turns. At the study's 317.8 mm: the ellipse's semi-axes are 64.34 and 34.57 mm, its radius 18.57 mm at the ends of
its major axis and 119.7 mm at the ends of its minor one; the ring's radius is 50.58 mm, and framed toward its boards
it is a circle path at theta_deg = 90 (u down its axis, v inward); the helix's tightest bend is 31.48, 25.05 and
20.15 mm at r/R = 0.1, 0.2 and 0.3, and framed k = 0 its tightest radius in the boards' plane 60.09, 29.14 and
21.06 mm.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal, NamedTuple

import numpy as np

from array16.particles import Arr, ClosedCurve

__all__ = [
    "HELIX_RATIOS",
    "T915S_RATIO",
    "Box",
    "Flat",
    "Framing",
    "Plate",
    "curvature_vector",
    "develop",
    "ellipse",
    "flat",
    "framing",
    "helix",
    "outline",
    "plate",
    "ring",
    "scaled",
    "smallest_box",
    "torus",
]

T915S_RATIO = 0.45  # r/R of T915-S's trefoil
HELIX_RATIOS = (0.1, 0.2, 0.3)  # the (1, 3) helix's r/R for H1: 0.2 is stage 1's


def _periodic_gradient(a: Arr, step: float) -> Arr:
    """da/ds by central differences on a closed curve's samples."""
    out: Arr = (np.roll(a, -1, 0) - np.roll(a, 1, 0)) / (2 * step)
    return out


# ---- curves ---------------------------------------------------------------------------------------------------------
def scaled(curve: ClosedCurve, length: float) -> ClosedCurve:
    """The curve scaled to a length (m)."""
    return ClosedCurve(curve.points * (length / curve.length), length)


def torus(ratio: float, length: float, n: int, p: int = 2, q: int = 3) -> ClosedCurve:
    """The (p, q) torus knot on a torus of radii r/R = ratio, scaled to a length (m), n points: (2, 3) is the trefoil,
    (1, 3) the unknotted helix that winds three times round the torus once round. The twin's torus path, to the bit."""
    if not 0 < ratio < 1:
        raise ValueError(f"a torus knot takes 0 < r/R < 1, not {ratio}")
    return scaled(ClosedCurve.torus_knot(1.0, ratio, n, p, q), length)


def curvature_vector(curve: ClosedCurve) -> Arr:
    """dT/ds at every point (n, 3) (1/m), by second differences, without its part along the tangent."""
    pts, ds = curve.points, curve.length / len(curve)
    k = (np.roll(pts, -1, 0) - 2 * pts + np.roll(pts, 1, 0)) / ds**2
    t = curve.tangents
    out: Arr = k - (k * t).sum(1)[:, None] * t
    return out


# ---- frames ---------------------------------------------------------------------------------------------------------
@dataclass(frozen=True, eq=False)
class Framing:
    """Board normals along a curve: its twist-free axes turned by theta0 + twist s, so that after a lap they come back
    turned by k pi (even k: each plate closes on itself; odd k: both plates are one band of twice the length).

    normal is the board normal n_b (toward the boards), across the direction T x n_b in the boards' plane; kappa_g and
    kappa_n are the curvature in the boards' plane and toward them (1/m)."""

    k: int
    theta0: float
    twist: float
    normal: Arr
    across: Arr
    kappa_g: Arr
    kappa_n: Arr

    @property
    def inplane_radius(self) -> float:
        """The smallest radius of curvature in the boards' plane (m): infinite where the bend is wholly toward them."""
        most = float(np.abs(self.kappa_g).max())
        return 1 / most if most > 0 else math.inf

    @property
    def toward_radius(self) -> float:
        """The smallest radius of curvature toward the boards (m): infinite where the bend is wholly in their plane."""
        most = float(np.abs(self.kappa_n).max())
        return 1 / most if most > 0 else math.inf

    @property
    def twist_deg_per_mm(self) -> float:
        """The twist against the twist-free frame in degrees per millimetre: +0.48 on T915-S framed k = 0."""
        return math.degrees(self.twist) * 1e-3


def framing(curve: ClosedCurve, k: int, theta0: float | None = None, starts: int = 720) -> Framing:
    """The framing that closes with k half turns against the twist-free frame, (k pi - holonomy)/L a metre; without
    theta0, the start among `starts` in [0, pi) that keeps the largest in-plane radius (docs/KNOTTED_GUIDE.md section
    3.5). The knotted study's, the same arithmetic as array16.twin.framing."""
    n, length = len(curve), curve.length
    e1, e2 = curve.axes
    kv = curvature_vector(curve)
    a, b = (kv * e1).sum(1), (kv * e2).sum(1)
    twist = (k * math.pi - curve.holonomy) / length
    s = np.arange(n) * length / n
    if theta0 is None:
        trial = np.linspace(0, math.pi, starts, endpoint=False)
        th = trial[:, None] + twist * s[None, :]
        kg = np.cos(th) * b[None, :] - np.sin(th) * a[None, :]  # kappa . (T x n_b), n_b = cos th e1 + sin th e2
        theta0 = float(trial[int(np.argmin(np.abs(kg).max(1)))])
    th = theta0 + twist * s
    normal = np.cos(th)[:, None] * e1 + np.sin(th)[:, None] * e2
    across = np.cross(curve.tangents, normal)
    return Framing(k, theta0, twist, normal, across, (kv * across).sum(1), (kv * normal).sum(1))


def flat(curve: ClosedCurve, way: Literal["in-plane", "toward"]) -> Framing:
    """A planar curve framed k = 0 with its bend wholly in the boards' plane (the boards parallel to the curve's plane,
    the normal along the plane's) or wholly toward them (the boards two strips of a cylinder through the curve, the
    normal in the plane): the start theta0 in [0, pi) that puts the normal there. Refused for a curve off its plane by
    more than 1e-9 of its size."""
    if way not in ("in-plane", "toward"):
        raise ValueError(f"a flat framing is 'in-plane' or 'toward', not {way!r}")
    centred = curve.points - curve.points.mean(axis=0)
    _, sizes, axes = np.linalg.svd(centred, full_matrices=False)
    if sizes[2] > 1e-9 * sizes[0]:
        raise ValueError(f"the curve leaves its plane by {sizes[2] / sizes[0]:.2g} of its size: it is not flat")
    plane = axes[2]
    e1, e2 = curve.axes
    normal = plane if way == "in-plane" else np.cross(curve.tangents[0], plane)
    theta0 = math.atan2(float(normal @ e2[0]), float(normal @ e1[0])) % math.pi
    return framing(curve, 0, theta0)


# ---- plates ---------------------------------------------------------------------------------------------------------
class Plate(NamedTuple):
    """A plate's centre line: its arc-length elements (m) and its geodesic curvature in the plate (1/m)."""

    dl: Arr
    kg: Arr

    @property
    def length(self) -> float:
        """Its length (m)."""
        return float(self.dl.sum())

    @property
    def turning(self) -> float:
        """Its total geodesic turning (rad): what its outline turns through when laid flat."""
        return float((self.kg * self.dl).sum())


class Flat(NamedTuple):
    """A plate laid flat: its centre line in the plane (m, (n, 2)) and its heading (rad)."""

    xy: Arr
    heading: Arr


class Box(NamedTuple):
    """The turn (rad) that fits plane points in their smallest bounding box, and the box's short and long sides (m)."""

    angle: float
    short: float
    long: float


def plate(curve: ClosedCurve, normal: Arr, offset: float) -> Plate:
    """A plate's centre line, the curve moved by offset along the board normal: its arc-length elements (m) and its
    geodesic curvature in the plate (1/m). KNT-10: the plates at +-h differ in length by -2h oint kappa_n ds +
    h^3 oint tau_g^2 kappa_n ds + O(h^5), tau_g the normal's turning across the boards' plane: the h^2 terms cancel
    (23.42 mm on T915-S framed k = 0, the h^3 term -0.69 um)."""
    ds = curve.length / len(curve)
    g = curve.points + offset * normal
    g1 = _periodic_gradient(g, ds)
    speed = np.linalg.norm(g1, axis=1)
    tg = g1 / speed[:, None]
    dt = _periodic_gradient(tg, ds) / speed[:, None]
    u = np.cross(tg, normal)
    u /= np.linalg.norm(u, axis=1)[:, None]
    return Plate(speed * ds, (dt * u).sum(1))


def develop(dl: Arr, kg: Arr) -> Flat:
    """A plate laid flat: its centre line in the plane (m, (n, 2)) and its heading, from its geodesic curvature."""
    heading = np.concatenate([[0.0], np.cumsum(kg * dl)])[:-1]
    xy = np.cumsum(np.stack([np.cos(heading), np.sin(heading)], 1) * dl[:, None], axis=0)
    return Flat(xy, heading)


def outline(xy: Arr, heading: Arr, width: float) -> Arr:
    """The outline (2n, 2) of a strip of a width along a centre line laid flat."""
    side = np.stack([-np.sin(heading), np.cos(heading)], 1) * width / 2
    out: Arr = np.concatenate([xy + side, (xy - side)[::-1]])
    return out


def smallest_box(points: Arr, steps: int = 720) -> Box:
    """The turn (rad) that puts plane points (m, 2) in their smallest bounding box among `steps` turns in [0, pi), and
    the box's short and long sides."""
    best = (math.inf, 0.0, 0.0, 0.0)
    for angle in np.linspace(0, math.pi, steps, endpoint=False):
        c, s = math.cos(angle), math.sin(angle)
        x, y = points[:, 0] * c - points[:, 1] * s, points[:, 0] * s + points[:, 1] * c
        w, h = float(np.ptp(x)), float(np.ptp(y))
        if w * h < best[0]:
            best = (w * h, float(angle), min(w, h), max(w, h))
    return Box(best[1], best[2], best[3])


# ---- stage 2's paths for H1 -----------------------------------------------------------------------------------------
def ring(length: float, n: int) -> ClosedCurve:
    """The flat ring of a length (m) in the xy plane, n points: the twin's circle path. Framed flat(ring, 'toward'), its
    boards are two coaxial cylinders, the gap radial and the bend wholly on their stiff side; flat(ring, 'in-plane') is
    stage 1's flat ring (designs/flat_ring_318.toml, u inward and v up)."""
    return ClosedCurve.circle(length / (2 * math.pi), n)


def helix(ratio: float, length: float, n: int) -> ClosedCurve:
    """The (1, 3) torus curve at r/R = ratio and a length (m), n points: unknotted and out of the plane, winding three
    times round the torus once round (designs/helix_13.toml at 0.2, framed k = 0)."""
    return torus(ratio, length, n, 1, 3)


def _perimeter(a: float, b: float, points: int = 4096) -> float:
    """The perimeter of the ellipse of semi-axes a and b by the trapezoidal rule on its speed over a period, exact to
    rounding for b/a above 0.1."""
    t = 2 * math.pi * np.arange(points) / points
    return float(2 * math.pi * np.sqrt((a * np.sin(t)) ** 2 + (b * np.cos(t)) ** 2).mean())


def ellipse(length: float, tightest: float, n: int) -> ClosedCurve:
    """The ellipse of perimeter `length` (m) whose least radius of curvature is `tightest` (m), b^2/a at the ends of its
    major axis 2a, centred on the origin in the xy plane with its major axis along x: n points equally spaced in arc
    length, each on the ellipse. Its semi-axes follow from the perimeter, which grows with a at b = sqrt(tightest a)
    from 2 pi tightest (the circle), by bisection; its points from the arc length as a function of the angle t of
    (a cos t, b sin t), the speed's Fourier series integrated term by term on 64 n angles, inverted by interpolation."""
    if not (tightest > 0 and length > 2 * math.pi * tightest):
        raise ValueError(
            f"an ellipse of perimeter {length} has its tightest radius under that of the circle of that perimeter, "
            f"{length / (2 * math.pi)}, and above 0: not {tightest}"
        )
    lo, hi = tightest, length / 4
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if _perimeter(mid, math.sqrt(tightest * mid)) < length else (lo, mid)
        if hi - lo < 1e-15 * hi:
            break
    a = (lo + hi) / 2
    b = math.sqrt(tightest * a)
    count = 64 * n
    t = 2 * math.pi * np.arange(count) / count
    speed = np.sqrt((a * np.sin(t)) ** 2 + (b * np.cos(t)) ** 2)
    c = np.fft.rfft(speed) / count
    harmonics = np.arange(len(c))
    integral = np.zeros_like(c)
    integral[1:] = c[1:] / (1j * harmonics[1:])  # the periodic part of s(t): its Fourier series term by term
    periodic = np.fft.irfft(integral * count, count)
    s = c[0].real * t + periodic - periodic[0]
    total = 2 * math.pi * c[0].real
    angle = np.interp(np.arange(n) * total / n, np.append(s, total), np.append(t, 2 * math.pi))
    points = np.stack([a * np.cos(angle), b * np.sin(angle), np.zeros(n)], 1)
    return ClosedCurve(points, total)
