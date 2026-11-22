"""The figures of docs/KNOTTED_GUIDE.md: a knotted RF guide built from several boards, in theory and simulation.

The document asks what a knot of boards takes: which arrangements of boards carry a knotted path, how fields behave
at their joints, how one channel drives the knot as a travelling-wave ring, and what electrons do in it. This script
computes the figures that carry that theory (section of the document; figure number):

  arrangements  three-dimensional views, at one scale, of the boards that carry a trefoil, and the one-board ring
                (3.2; figure 1)
  knots         the four knots' layouts on two levels; energy held against joint angle and bend radius; length
                against W/D with the bound and the ropelength floors (3.1; 2)
  ramp          a ramp compiled on two slotted boards' four layers, the cut-out that must fail, and the slot guide's
                section with its escape saddle (3.3; 3)
  strips        curvature, torsion and the straight-strip limit; framed strips' curvature in and out of their
                plane; the plates laid flat on a fab panel (3.5; 4)
  joints        an offset joint in three dimensions with its null coloured by axial field; decay by parity and the
                width floor; one segment's four basis fields and the banded operator of a segmented trefoil (4; 5)
  drive         the ring model: |V(s)| on the trefoil for three drives, ripple against the number of feeds, and
                the singular values of the drive operator (5.3; 6)
  laminate      Q against the split and the laminate per band; the power a 1 eV well takes against frequency with
                the bend-limited energy; the multipactor chart with the operating lines (5.4, 5.5; 7)
  held          held and lost against energy on the trefoil and the ring; the facing boards' well and force; the
                energy held against curvature along the knot (6.1; 9)
  stopbands     reach against energy for each joint case with the resonance families' centres, and the fine scan
                round tune 54 on the trefoil against the flat ring (6.2; 10)
  acceleration  energy against lap for the six-gap pi-mode ramp with its band crossings, and the crossing gain
                against the ramp rate (6.4; 11)

Figures 8 and 12 of the document are diagrams drawn in the document itself.

    uv run python experiments/knotted_guide_study.py              # every figure, a few minutes on 12 cores
    uv run python experiments/knotted_guide_study.py --quick      # coarse grids, few draws: knotted_guide/quick/
    uv run python experiments/knotted_guide_study.py --only knots,strips

Models: curves and framings from array16.particles.ClosedCurve; the boards' fields from the gapless facing-board
model (array16.particles.FacingBoards, ELE-2), from Fourier kernels between two planes, and from finite-difference
Laplace solves on grids (vacuum in three dimensions; laminate by finite volumes in two); the drive from a ring of
transmission line with joints and feeds; electrons tracked with array16.particles.integrate in the fields of four
thin rods (QuadRods) or of an ideal quadrupole channel. Outputs: knotted_guide/{*.png, results.json, run.log}.
"""

from __future__ import annotations

import argparse
import heapq
import json
import math
import os
import sys
import time
from collections.abc import Callable, Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from array16 import lab
from array16 import particles as pt
from array16.particles import Arr, BoolArr, CArr, ClosedCurve, FacingBoards, GuideWatch, IntArr, QuadRods

OUT = Path(__file__).resolve().parent / "knotted_guide"
LOG: list[str] = []
MM = 1e-3

# Ink, surfaces and the categorical order of the reference palette (fixed order, never cycled).
INK, INK2, MUTED, GRID, EDGE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
BLUE, ORANGE, GREEN, YELLOW, PINK, DGREEN, PURPLE, RED = SERIES
BLUES = ("#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
SURFACE, WARM, PALE = "#fcfcfb", "#fbe3d6", "#ecebe6"
COPPER, LAMINATE, FLEX = "#c07a3a", "#bfd8a8", "#e8c46a"
SUPERSCRIPTS = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")

# The lead designs (docs/KNOTTED_GUIDE.md sections 3 and 9).
F915 = 915e6
FH = 0.839e6  # f h (Hz m) of a 1 eV well near q = 0.5 on facing boards: 50.3 V
H915 = 0.915e-3  # T915-S: null to copper (gap 1.83 mm)
L915 = 0.318  # T915-S: m = 1 at 915 MHz on 25 um polyimide, n_eff 1.031
RATIO915 = 0.45  # T915-S: r / R of its trefoil
WIDTH915 = 8 * H915  # T915-S: strip width
QUARTER915 = L915 / 4  # a quarter guided wave at m = 1: the two closing joints' and the two feeds' spacing
SCARF = 29.4e-3  # three secular wavelengths at 8 eV (q = 0.5, 915 MHz): the joints' scarf
SPACERS = 24
HA = 0.92e-3  # section A (slotted RO4003C) at 915 MHz
FIRST_RUN = (0.020, 0.006)  # the first run's trefoil, R and r (m): 276.5 mm, r / R = 0.30
PANEL_MM = (234.0, 490.0)  # the flex fab's panel


def say(line: str = "") -> None:
    print(line, flush=True)
    LOG.append(line)


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


def periodic_gradient(a: Arr, step: float) -> Arr:
    """da/ds by central differences on a closed curve's samples."""
    out: Arr = (np.roll(a, -1, 0) - np.roll(a, 1, 0)) / (2 * step)
    return out


# ---- geometry: curves, frames, strips, layouts ---------------------------------------------------------------------
def scaled(curve: ClosedCurve, length: float) -> ClosedCurve:
    """The curve scaled to a length (m)."""
    return ClosedCurve(curve.points * (length / curve.length), length)


def torus(ratio: float, length: float, n: int, p: int = 2, q: int = 3) -> ClosedCurve:
    """The (p, q) torus knot on a torus of radii r / R = ratio, scaled to a length (m). (2, 3) is the trefoil."""
    return scaled(ClosedCurve.torus_knot(1.0, ratio, n, p, q), length)


def first_run_trefoil(n: int = 2001) -> ClosedCurve:
    """The first run's trefoil: torus radii 20 and 6 mm, 276.5 mm long, tightest bend 16.6 mm."""
    return ClosedCurve.torus_knot(FIRST_RUN[0], FIRST_RUN[1], n)


def figure_eight_knot(length: float, n: int) -> ClosedCurve:
    """The figure-eight knot 4_1 in its common parametrisation ((2 + cos 2t) cos 3t, (2 + cos 2t) sin 3t, sin 4t)."""
    t = np.linspace(0, 2 * math.pi, 40 * n, endpoint=False)
    pts = np.stack([(2 + np.cos(2 * t)) * np.cos(3 * t), (2 + np.cos(2 * t)) * np.sin(3 * t), np.sin(4 * t)], 1)
    return scaled(ClosedCurve.from_polygon(pts, n), length)


def curvature_vector(curve: ClosedCurve) -> Arr:
    """dT/ds at every point (n, 3), by second differences, without its part along the tangent."""
    pts, ds = curve.points, curve.length / len(curve)
    k = (np.roll(pts, -1, 0) - 2 * pts + np.roll(pts, 1, 0)) / ds**2
    t = curve.tangents
    out: Arr = k - (k * t).sum(1)[:, None] * t
    return out


def curvature_torsion(curve: ClosedCurve, keep: int = 80) -> tuple[Arr, Arr]:
    """Curvature and torsion (1/m) at every point, from the curve's Fourier series in arc length kept to `keep`
    harmonics a lap (the second and third derivatives of a polygon's samples are too rough for torsion)."""
    n, length = len(curve), curve.length
    c = np.fft.rfft(curve.points, axis=0)
    c[keep + 1 :] = 0
    ik = (2j * math.pi / length * np.arange(c.shape[0]))[:, None]
    d1, d2, d3 = (np.fft.irfft(c * ik**m, n, axis=0) for m in (1, 2, 3))
    cross = np.cross(d1, d2)
    kappa: Arr = np.linalg.norm(cross, axis=1) / np.linalg.norm(d1, axis=1) ** 3
    tau: Arr = (cross * d3).sum(1) / (cross * cross).sum(1)
    return kappa, tau


@dataclass(frozen=True)
class Framing:
    """Board normals along a curve: its twist-free axes turned by theta0 + twist s, so that after a lap they come back
    turned by k pi (even k: each plate closes on itself; odd k: both plates are one band of twice the length).

    normal is the board normal n_b (toward the boards), across the direction T x n_b in the boards' plane;
    kappa_g and kappa_n are the curvature in the boards' plane and toward them (1/m)."""

    k: int
    theta0: float
    twist: float
    normal: Arr
    across: Arr
    kappa_g: Arr
    kappa_n: Arr

    @property
    def inplane_radius(self) -> float:
        """The smallest radius of curvature in the boards' plane (m)."""
        return float(1 / np.abs(self.kappa_g).max())


def framing(curve: ClosedCurve, k: int, theta0: float | None = None, starts: int = 720) -> Framing:
    """The framing that closes with k half turns against the twist-free frame, (k pi - holonomy) / L a metre; without
    theta0, the start that keeps the largest in-plane radius (section 3.5)."""
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


def plate(curve: ClosedCurve, normal: Arr, offset: float) -> tuple[Arr, Arr]:
    """A plate's centre line, the curve moved by offset along the board normal: its arc-length elements (m) and its
    geodesic curvature in the plate (1/m)."""
    ds = curve.length / len(curve)
    g = curve.points + offset * normal
    g1 = periodic_gradient(g, ds)
    speed = np.linalg.norm(g1, axis=1)
    tg = g1 / speed[:, None]
    dt = periodic_gradient(tg, ds) / speed[:, None]
    u = np.cross(tg, normal)
    u /= np.linalg.norm(u, axis=1)[:, None]
    return speed * ds, (dt * u).sum(1)


def develop(dl: Arr, kg: Arr) -> tuple[Arr, Arr]:
    """A plate laid flat: its centre line in the plane (m, 2) and its heading, from its geodesic curvature."""
    heading = np.concatenate([[0.0], np.cumsum(kg * dl)])[:-1]
    xy = np.cumsum(np.stack([np.cos(heading), np.sin(heading)], 1) * dl[:, None], axis=0)
    return xy, heading


def outline(xy: Arr, heading: Arr, width: float) -> Arr:
    """The outline (2m, 2) of a strip of a width along a centre line laid flat."""
    side = np.stack([-np.sin(heading), np.cos(heading)], 1) * width / 2
    out: Arr = np.concatenate([xy + side, (xy - side)[::-1]])
    return out


def smallest_box(points: Arr, steps: int = 720) -> tuple[float, float, float]:
    """The turn (rad) that puts plane points (m, 2) in their smallest bounding box, and the box's short and long
    sides."""
    best = (math.inf, 0.0, 0.0, 0.0)
    for angle in np.linspace(0, math.pi, steps, endpoint=False):
        c, s = math.cos(angle), math.sin(angle)
        x, y = points[:, 0] * c - points[:, 1] * s, points[:, 0] * s + points[:, 1] * c
        w, h = float(np.ptp(x)), float(np.ptp(y))
        if w * h < best[0]:
            best = (w * h, float(angle), min(w, h), max(w, h))
    return best[1], best[2], best[3]


def quintic(u: Arr) -> Arr:
    """A climb from -1 to +1 over u in [-1, 1] with no slope or curvature at its ends: (15u - 10u^3 + 3u^5) / 8."""
    v = np.clip(u, -1, 1)
    out: Arr = (15 * v - 10 * v**3 + 3 * v**5) / 8
    return out


def two_level_trefoil(r0: float, a: float, rise: float, ramp: float, n: int) -> tuple[ClosedCurve, Arr]:
    """S2-C3's layout: rho = r0 + a cos(3 phi / 2) over two turns, on two levels `rise` apart joined by six quintic
    ramps `ramp` long (in the plane) centred on the lobes, where the torus knot's height sin 3t changes sign. Returns
    the curve and its ramps' centres as fractions of a lap (section 3.3)."""
    t = np.linspace(0, 2 * math.pi, 40 * n, endpoint=False)
    rho = r0 + a * np.cos(3 * t)
    speed = np.hypot(-3 * a * np.sin(3 * t), 2 * rho)
    s = np.concatenate([[0.0], np.cumsum(speed[:-1])]) * (t[1] - t[0])
    total = s[-1] + speed[-1] * (t[1] - t[0])
    z = 0.5 * rise * np.sign(np.sin(3 * t))
    centres = []
    for k in range(6):
        sk = s[int(np.argmin(np.abs(t - k * math.pi / 3)))]
        d = (s - sk + total / 2) % total - total / 2
        near = np.abs(d) < ramp / 2
        z[near] = 0.5 * rise * (1 if k % 2 == 0 else -1) * quintic(d[near] / (ramp / 2))
        centres.append(sk / total)
    pts = np.stack([rho * np.cos(2 * t), rho * np.sin(2 * t), z], 1)
    return ClosedCurve.from_polygon(pts, n), np.array(centres)


def six_stick_trefoil(perimeter: float) -> Arr:
    """The D3 six-stick trefoil (Randell's stick number 6): two triangles at heights +-h, the second turned 215
    degrees, scaled to a perimeter (m). Its corners turn 128.8 degrees and its corner planes meet at 46.1 and 20.1
    degrees (section 3.4)."""
    v = []
    for k in range(3):
        a = 2 * math.pi * k / 3
        v.append([math.cos(a), math.sin(a), 0.5])
        v.append([math.cos(a + math.radians(215)), math.sin(a + math.radians(215)), -0.5])
    pts = np.array(v)
    out: Arr = pts * perimeter / float(np.linalg.norm(np.roll(pts, -1, 0) - pts, axis=1).sum())
    return out


def lattice_trefoil(pitch: float) -> Arr:
    """A lattice trefoil of 24 unit steps and 14 sticks (Scharein et al.'s minimum step number), centred, at a pitch
    (m)."""
    corners = np.array(
        [
            (1, 1, -1),
            (1, 1, 0),
            (-2, 1, 0),
            (-2, -1, 0),
            (0, -1, 0),
            (0, 0, 0),
            (0, 0, -1),
            (0, 2, -1),
            (-1, 2, -1),
            (-1, 2, 1),
            (-1, 0, 1),
            (-1, 0, -1),
            (-1, -1, -1),
            (1, -1, -1),
        ],
        dtype=np.float64,
    )
    out: Arr = (corners - corners.mean(0)) * pitch
    return out


def corner_turns(vertices: Arr) -> tuple[Arr, Arr]:
    """At each corner of a closed polygon: the angle (rad) its path turns, and the unit normal of its plane."""
    d = np.roll(vertices, -1, 0) - vertices
    t = d / np.linalg.norm(d, axis=1)[:, None]
    before = np.roll(t, 1, 0)
    turn = np.arccos(np.clip((before * t).sum(1), -1, 1))
    normal = np.cross(before, t)
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    return turn, normal


@dataclass(frozen=True)
class Rounded:
    """A polygon with its corners rounded to arcs in their own planes, sampled densely: points (m, 3), tangents,
    board normals (the corner plane's normal on each corner, carried along the straights with a twist where two
    corners' planes differ), and per point which corner (>= 0) or straight (-1 - index) it is on, and whether the
    board normal turns there."""

    points: Arr
    tangents: Arr
    normals: Arr
    piece: IntArr
    twisting: BoolArr


def rounded(
    vertices: Arr, radius: float, normals: Arr, twist_length: float, step: float = 0.4e-3, swap: bool = False
) -> Rounded:
    """vertices rounded to `radius`, each corner's boards with the given normal; along a straight whose ends' normals
    differ the boards turn about the path over twist_length at its middle (a twisted flex section), or, with swap,
    hand over at its middle in one step (a polarity swap at 90 degrees)."""
    m = len(vertices)
    pts: list[Arr] = []
    tans: list[Arr] = []
    nrm: list[Arr] = []
    piece: list[int] = []
    twisting: list[bool] = []
    ends = []
    for i in range(m):
        a, b, c = vertices[i - 1], vertices[i], vertices[(i + 1) % m]
        t1 = (b - a) / np.linalg.norm(b - a)
        t2 = (c - b) / np.linalg.norm(c - b)
        th = math.acos(float(np.clip(t1 @ t2, -1, 1)))
        d = radius * math.tan(th / 2)
        p1, p2 = b - d * t1, b + d * t2
        inward = (t2 - t1) / np.linalg.norm(t2 - t1)
        centre = b + inward * radius / math.cos(th / 2)
        u1 = p1 - centre
        axis = np.cross(t1, t2)
        axis /= np.linalg.norm(axis)
        steps = max(int(radius * th / step), 4)
        samples = [
            centre + u1 * math.cos(phi) + np.cross(axis, u1) * math.sin(phi) for phi in np.linspace(0, th, steps)
        ]
        ends.append((p1, p2, np.array(samples), t1, t2, axis))
    for i in range(m):
        p1, p2, arc, t1, t2, axis = ends[i]
        for j, p in enumerate(arc):
            frac = j / max(len(arc) - 1, 1)
            tangent = t1 * (1 - frac) + t2 * frac
            pts.append(p)
            tans.append(tangent / np.linalg.norm(tangent))
            nrm.append(normals[i])
            piece.append(i)
            twisting.append(False)
        start, stop = p2, ends[(i + 1) % m][0]
        length = float(np.linalg.norm(stop - start))
        direction = (stop - start) / length
        n0, n1 = normals[i], normals[(i + 1) % m]
        if float(n0 @ n1) < 0:
            n1 = -n1
        angle = math.atan2(float(np.cross(n0, n1) @ direction), float(n0 @ n1))
        count = max(int(length / step), 2)
        for u in np.linspace(0, 1, count, endpoint=False)[1:]:
            along = u * length - length / 2
            if swap:
                frac = 0.0 if along < 0 else 1.0
                turning = abs(along) < step
            else:
                frac = float(np.clip(along / twist_length + 0.5, 0, 1))
                turning = 0 < frac < 1 and abs(angle) > 1e-6
            rot = angle * frac
            normal = n0 * math.cos(rot) + np.cross(direction, n0) * math.sin(rot)
            pts.append(start + u * length * direction)
            tans.append(direction)
            nrm.append(normal)
            piece.append(-1 - i)
            twisting.append(turning)
    return Rounded(np.array(pts), np.array(tans), np.array(nrm), np.array(piece), np.array(twisting))


def plane_runs(normals: Arr) -> list[list[int]]:
    """Consecutive corners of a closed polygon whose planes coincide, grouped into runs, the first run starting where
    the plane changes."""
    m = len(normals)
    same = [abs(float(normals[i] @ normals[(i + 1) % m])) > 1 - 1e-9 for i in range(m)]
    start = next((i + 1) % m for i in range(m) if not same[i])
    runs: list[list[int]] = [[start]]
    for j in range(1, m):
        i = (start + j) % m
        if same[(i - 1) % m]:
            runs[-1].append(i)
        else:
            runs.append([i])
    return runs


@dataclass(frozen=True)
class Crossing:
    """A crossing of a closed polygon's projection on the xy plane: the curve parameters (point index plus fraction)
    of the strand passing over and the one passing under, where it is, and the angle between the strands."""

    over: float
    under: float
    at: tuple[float, float]
    angle: float


def crossings(points: Arr) -> list[Crossing]:
    """Where the projection of a closed polygon on the xy plane crosses itself, the higher strand passing over."""
    p = points[:, :2]
    q = np.roll(p, -1, 0)
    d = q - p
    n = len(p)
    i, j = np.triu_indices(n, 2)
    keep = ~((i == 0) & (j == n - 1))
    i, j = i[keep], j[keep]
    den = d[i, 0] * d[j, 1] - d[i, 1] * d[j, 0]
    with np.errstate(divide="ignore", invalid="ignore"):
        rx, ry = p[j, 0] - p[i, 0], p[j, 1] - p[i, 1]
        ta = (rx * d[j, 1] - ry * d[j, 0]) / den
        tb = (rx * d[i, 1] - ry * d[i, 0]) / den
    hit = (np.abs(den) > 0) & (ta >= 0) & (ta < 1) & (tb >= 0) & (tb < 1)
    out = []
    for a, b, sa, sb in zip(i[hit], j[hit], ta[hit], tb[hit], strict=True):
        za = points[a, 2] + sa * (points[(a + 1) % n, 2] - points[a, 2])
        zb = points[b, 2] + sb * (points[(b + 1) % n, 2] - points[b, 2])
        at = p[a] + sa * d[a]
        cosang = abs(float(d[a] @ d[b])) / float(np.linalg.norm(d[a]) * np.linalg.norm(d[b]))
        angle = math.degrees(math.acos(min(cosang, 1.0)))
        pa, pb = float(a + sa), float(b + sb)
        out.append(Crossing(pa if za > zb else pb, pb if za > zb else pa, (float(at[0]), float(at[1])), angle))
    return out


def ramps_needed(cross: Sequence[Crossing]) -> int:
    """The least number of ramps that put a diagram on two levels: the switches between passing over and passing
    under along the curve (each over-pass on the upper level, each under-pass on the lower)."""
    events = sorted([(c.over, 1) for c in cross] + [(c.under, 0) for c in cross])
    kinds = [kind for _, kind in events]
    return sum(1 for k in range(len(kinds)) if kinds[k] != kinds[k - 1])


def determinant(cross: Sequence[Crossing]) -> int:
    """The knot determinant |Delta(-1)| of a diagram, from its colouring matrix: arcs run between under-passes, and at
    each crossing twice the over-arc less the two under-arcs is zero."""
    c = len(cross)
    if c == 0:
        return 1
    unders = sorted(x.under for x in cross)

    def arc(param: float) -> int:
        """The arc a curve parameter lies on: arc k starts at the k-th under-pass."""
        k = int(np.searchsorted(unders, param, side="right")) - 1
        return k % c

    m = np.zeros((c, c))
    for x in cross:
        k = unders.index(x.under)
        m_in, m_out, m_over = (k - 1) % c, k, arc(x.over)
        row = unders.index(x.under)
        m[row, m_over] += 2
        m[row, m_in] -= 1
        m[row, m_out] -= 1
    return round(abs(float(np.linalg.det(m[1:, 1:])))) if c > 1 else 1


def rotational_order(points: Arr, top: int = 12, tol: float = 2e-3) -> int:
    """The largest k for which turning the curve by 2 pi / k about z maps it onto itself: every turned sample within
    half the polygon's longest side, plus tol of its size, of a sample."""
    size = float(np.ptp(points[:, :2]))
    spacing = float(np.linalg.norm(np.roll(points, -1, 0) - points, axis=1).max())
    for k in range(top, 1, -1):
        c, s = math.cos(2 * math.pi / k), math.sin(2 * math.pi / k)
        turned = points @ np.array([[c, s, 0], [-s, c, 0], [0, 0, 1.0]])
        step = max(len(points) // 400, 1)
        d = np.linalg.norm(turned[::step, None, :] - points[None, :, :], axis=2).min(1)
        if float(d.max()) < 0.75 * spacing + tol * size:
            return k
    return 1


# ---- arrangements (figure 1) --------------------------------------------------------------------------------------
def compute_arrangements(quick: bool) -> dict[str, Any]:
    say("\n== arrangements: boards carrying a trefoil, at one scale ==")
    n = 1200 if quick else 3000
    out: dict[str, Any] = {}
    # S2-C3: two slotted RO4003C boards 0.813 mm thick and 1.0 mm apart; levels at the slots' mid-planes
    board, gap = 0.813e-3, 1.0e-3
    rise = board + gap
    s2c3, centres = two_level_trefoil(33.2e-3, 4e-3, rise, 16e-3, n)
    flat = ClosedCurve.from_polygon(s2c3.points * np.array([1.0, 1.0, 0.0]), len(s2c3))
    cross = crossings(s2c3.points[:: max(n // 1500, 1)])
    out["S2-C3"] = {
        "length_mm": s2c3.length / MM,
        "inplane_bend_mm": flat.bend_radius / MM,
        "bend_mm": s2c3.bend_radius / MM,
        "rise_mm": rise / MM,
        "ramp_mm": 16.0,
        "crossings": len(cross),
        "crossing_angle_deg": float(np.mean([c.angle for c in cross])),
        "total_curvature_pi": s2c3.total_curvature / math.pi,
        "_curve": s2c3,
        "_ramps": centres,
        "_crossings": cross,
    }
    say(
        f"  S2-C3: {s2c3.length / MM:.1f} mm, tightest bend in the plane {flat.bend_radius / MM:.1f} mm "
        f"({s2c3.bend_radius / MM:.1f} mm "
        f"with the 16 mm ramps' climb), {len(cross)} crossings at {out['S2-C3']['crossing_angle_deg']:.1f} deg, six "
        "ramps "
        f"between levels {rise / MM:.3f} mm apart"
    )
    # rounded sticks: the six-stick trefoil whose path, its corners rounded to 8.0 mm, is 291 mm (m = 1 on section A)
    unit = six_stick_trefoil(1.0)
    turn, cnorm = corner_turns(unit)
    cut = float(np.sum(2 * np.tan(turn / 2) - turn))  # the length rounding takes, per metre of radius
    perimeter = 0.291 + 8.0e-3 * cut
    hexagon = unit * perimeter
    sides = np.linalg.norm(np.roll(hexagon, -1, 0) - hexagon, axis=1)
    straight = float(sides.min() - 2 * 8.0e-3 * math.tan(float(turn.max()) / 2))
    largest = (float(sides.min() * 0.291 / perimeter) - 10.4e-3) / (
        2 * math.tan(float(turn.max()) / 2) - float(sides.min() / perimeter) * cut
    )
    plane_angles = [
        math.degrees(math.acos(min(abs(float(cnorm[i] @ cnorm[(i + 1) % 6])), 1.0))) for i in range(len(cnorm))
    ]
    sticks = rounded(hexagon, 8.0e-3, cnorm, 10.4e-3)
    out["sticks"] = {
        "path_mm": float(np.linalg.norm(np.diff(sticks.points, axis=0), axis=1).sum()) / MM,
        "hexagon_mm": perimeter / MM,
        "corner_turn_deg": float(np.degrees(turn).mean()),
        "plane_angles_deg": sorted({round(a, 1) for a in plane_angles}),
        "shortest_straight_mm": straight / MM,
        "largest_corner_mm": largest / MM,
        "_vertices": hexagon,
        "_rounded": sticks,
    }
    say(
        f"  sticks: corners turn {np.degrees(turn).mean():.1f} deg, consecutive corner planes meet at "
        f"{sorted({round(a, 1) for a in plane_angles})} deg; rounded to 8.0 mm the path is "
        f"{out['sticks']['path_mm']:.0f} mm on a "
        f"{perimeter / MM:.0f} mm hexagon, its shortest straight {straight / MM:.1f} mm; at that path the corners "
        "round to "
        f"{largest / MM:.1f} mm before a straight falls under 10.4 mm (a 291 mm hexagon rounded to 8.0 mm keeps "
        f"{(float(sides.min()) * 0.291 / perimeter - 2 * 8.0e-3 * math.tan(float(turn.max()) / 2)) / MM:.1f} mm "
        f"straights and a {(0.291 - 8.0e-3 * cut) / MM:.0f} mm path)"
    )
    # the lattice trefoil, 24 steps of 26.2 mm, corners rounded to 7.9 mm; one facing pair per planar run
    lattice = lattice_trefoil(26.2e-3)
    lturn, lnorm = corner_turns(lattice)
    runs = plane_runs(lnorm)
    run_normal = np.zeros_like(lnorm)
    for run in runs:
        run_normal[run] = lnorm[run[0]]
    lat = rounded(lattice, 7.9e-3, run_normal, 0.0, swap=True)
    steps = float(np.abs(np.roll(lattice, -1, 0) - lattice).sum() / 26.2e-3)
    handovers = [i for i in range(len(lattice)) if abs(float(run_normal[i] @ run_normal[(i + 1) % len(lattice)])) < 0.5]
    out["L26"] = {
        "steps": steps,
        "sticks": len(lattice),
        "runs": len(runs),
        "handovers": len(handovers),
        "corner_turn_deg": float(np.degrees(lturn).mean()),
        "length_mm": float(np.linalg.norm(np.diff(lat.points, axis=0), axis=1).sum()) / MM,
        "_vertices": lattice,
        "_rounded": lat,
        "_runs": runs,
        "_handovers": handovers,
    }
    say(
        f"  L26: {steps:.0f} steps, {len(lattice)} sticks, {len(runs)} planar runs, corners of "
        f"{np.degrees(lturn).mean():.0f} deg"
    )
    # T915-S: the r/R = 0.45 trefoil at 318 mm, strips framed with k = 0
    knot = torus(RATIO915, L915, n)
    fr = framing(knot, 0)
    out["T915-S"] = {
        "length_mm": L915 / MM,
        "bend_mm": knot.bend_radius / MM,
        "twist_deg_per_mm": math.degrees(fr.twist) * MM,
        "inplane_radius_mm": fr.inplane_radius / MM,
        "width_mm": WIDTH915 / MM,
        "_curve": knot,
        "_framing": fr,
    }
    say(
        f"  T915-S: {L915 / MM:.0f} mm, tightest bend {knot.bend_radius / MM:.2f} mm; k = 0 twist "
        f"{math.degrees(fr.twist) * MM:+.3f} deg/mm, in-plane radius >= {fr.inplane_radius / MM:.1f} mm"
    )
    # R1: the channel ring m = 1 round on RO4003C, 210 mm
    out["R1"] = {"length_mm": 210.0, "radius_mm": 210.0 / (2 * math.pi), "posts": 3}
    return out


def _axes3d(fig: Any, where: Any, limit: float, title: str) -> Any:
    """A three-dimensional panel at one scale (mm), the house style's light panes."""
    ax = fig.add_subplot(where, projection="3d", computed_zorder=False)
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_zlim(-limit * 0.5, limit * 0.5)
    ax.set_box_aspect((1, 1, 0.5))
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((0.985, 0.985, 0.98, 1.0))
        axis._axinfo["grid"]["color"] = GRID
        axis._axinfo["grid"]["linewidth"] = 0.5
    ax.tick_params(labelsize=7, pad=-1)
    ax.set_xlabel("mm", fontsize=7.5, labelpad=-6)
    ax.set_title(title, fontsize=8.8, loc="left")
    ax.view_init(elev=26, azim=-58)
    return ax


def _poly3d(polygons: Any, **style: Any) -> Any:
    """matplotlib's Poly3DCollection (its toolkit has no type information)."""
    import importlib

    return importlib.import_module("mpl_toolkits.mplot3d.art3d").Poly3DCollection(polygons, **style)


def _slab(ax: Any, corners: Arr, z0: float, z1: float, colour: str, alpha: float) -> None:
    """A board as a translucent slab: a plane polygon (k, 2) between heights z0 and z1 (mm)."""
    top = [np.column_stack([corners, np.full(len(corners), z1)])]
    bottom = [np.column_stack([corners, np.full(len(corners), z0)])]
    sides = [
        np.array(
            [
                [*corners[i], z0],
                [*corners[(i + 1) % len(corners)], z0],
                [*corners[(i + 1) % len(corners)], z1],
                [*corners[i], z1],
            ]
        )
        for i in range(len(corners))
    ]
    ax.add_collection3d(_poly3d(top + bottom + sides, facecolor=colour, edgecolor=colour, alpha=alpha, lw=0.4))


def _ribbon(ax: Any, centre: Arr, across: Arr, half: float, colour: str, alpha: float, step: int = 3) -> None:
    """A strip of half-width `half` along a centre line (k, 3), spread along `across` (k, 3), as quads."""
    quads = []
    k = len(centre)
    for i in range(0, k - 1, step):
        j = min(i + step, k - 1)
        quads.append(
            [
                centre[i] - half * across[i],
                centre[i] + half * across[i],
                centre[j] + half * across[j],
                centre[j] - half * across[j],
            ]
        )
    ax.add_collection3d(_poly3d(quads, facecolor=colour, edgecolor="none", alpha=alpha))


def _split_ribbon(
    ax: Any, centre: Arr, across: Arr, half: float, gap: float, colours: tuple[str, str], alpha: float
) -> None:
    """A facing-board strip split down its centre line into two halves at opposite volts, each its own colour."""
    _ribbon(ax, centre + across * (half + gap) / 2, across, (half - gap) / 2, colours[0], alpha)
    _ribbon(ax, centre - across * (half + gap) / 2, across, (half - gap) / 2, colours[1], alpha)


def _title3d(ax: Any, text: str) -> None:
    """A three-dimensional panel's title, wrapped to its width."""
    import textwrap

    ax.set_title(textwrap.fill(text, 64), fontsize=8.4, loc="left", pad=-2)


def figure_arrangements(res: dict[str, Any], out: Path) -> None:
    """Figure 1: three-dimensional views at one scale of the arrangements carrying a trefoil, and the one-board ring."""
    plt = _pyplot()
    fig = plt.figure(figsize=(15.5, 10.6))
    grid = fig.add_gridspec(2, 3, hspace=0.08, wspace=0.02)
    limit = 52.0
    # (a) S2-C3
    a = res["S2-C3"]
    curve: ClosedCurve = a["_curve"]
    exaggerate = 8.0
    ax = _axes3d(fig, grid[0, 0], limit, "")
    _title3d(
        ax,
        "S2-C3: two slotted RO4003C boards 0.813 mm thick and 1.0 mm apart, the path on the slots' mid-planes, six 16 "
        "mm "
        f"ramps; {a['length_mm']:.0f} mm (m = 2), tightest bend in the plane {a['inplane_bend_mm']:.1f} mm, "
        f"{a['crossings']} crossings at {a['crossing_angle_deg']:.1f}°; heights ×{exaggerate:g}",
    )
    p = curve.points / MM
    z = p[:, 2] * exaggerate
    square = np.array([[-46, -46], [46, -46], [46, 46], [-46, 46]], dtype=np.float64)
    half_rise = a["rise_mm"] / 2
    for lo, hi in ((-half_rise - 0.4065, -half_rise + 0.4065), (half_rise - 0.4065, half_rise + 0.4065)):
        _slab(ax, square, lo * exaggerate, hi * exaggerate, LAMINATE, 0.14)
    upper = p[:, 2] > 0
    for level, colour in ((upper, BLUE), (~upper, ORANGE)):
        pts = np.where(level[:, None], np.column_stack([p[:, :2], z]), np.nan)
        ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=COPPER, lw=4.0, alpha=0.5, solid_capstyle="butt")
        ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=colour, lw=1.4)
    ramp_idx = (np.asarray(a["_ramps"]) * len(p)).astype(int) % len(p)
    ax.scatter(p[ramp_idx, 0], p[ramp_idx, 1], z[ramp_idx], marker="^", s=40, color=INK, depthshade=False, zorder=9)
    for c in a["_crossings"]:
        ax.scatter([c.at[0] / MM], [c.at[1] / MM], [0.0], marker="x", s=44, color=RED, depthshade=False, zorder=9)
    ax.plot([], [], color=BLUE, label="upper level: the upper board's slot")
    ax.plot([], [], color=ORANGE, label="lower level: the lower board's slot")
    ax.plot([], [], color=COPPER, lw=4, alpha=0.5, label="slot-edge copper strips")
    ax.scatter([], [], marker="^", color=INK, label="ramp, on a lobe (6)")
    ax.scatter([], [], marker="x", color=RED, label="crossing (3)")
    ax.legend(loc="upper left", fontsize=6.9, bbox_to_anchor=(0.0, 0.9))
    # (b) rounded sticks
    b = res["sticks"]
    r: Rounded = b["_rounded"]
    ax = _axes3d(fig, grid[0, 1], limit, "")
    _title3d(
        ax,
        "rounded sticks: six facing pairs on the D3 six-stick trefoil, its corners (turning "
        f"{b['corner_turn_deg']:.1f}°) "
        f"rounded to 8.0 mm for a {b['path_mm']:.0f} mm path (m = 1); corner planes {b['plane_angles_deg'][0]:g}° and "
        f"{b['plane_angles_deg'][1]:g}° apart, so each pair turns in a 10.4 mm twisted flex handover; gap ×4",
    )
    _facing_pairs(ax, r, HA / MM * 4, 2 * HA / MM)
    v = b["_vertices"] / MM
    vv = np.vstack([v, v[:1]])
    ax.plot(vv[:, 0], vv[:, 1], vv[:, 2], color=MUTED, lw=0.7, ls="--")
    corners = np.array([r.points[r.piece == i].mean(0) / MM for i in range(len(v))])
    ax.scatter(corners[:, 0], corners[:, 1], corners[:, 2], marker="o", s=24, color=INK, depthshade=False, zorder=9)
    ax.plot([], [], color=GREEN, lw=4, label="twisted flex handover (6)")
    ax.scatter([], [], marker="o", color=INK, label="rounded corner (6)")
    ax.plot([], [], color=MUTED, lw=0.7, ls="--", label="the sharp six sticks: must fail")
    ax.legend(loc="upper left", fontsize=6.9, bbox_to_anchor=(0.0, 0.9))
    # (c) the lattice trefoil
    c = res["L26"]
    rl: Rounded = c["_rounded"]
    ax = _axes3d(fig, grid[0, 2], limit, "")
    _title3d(
        ax,
        f"L26: the lattice trefoil, {c['steps']:.0f} steps of 26.2 mm, its {c['sticks']} corners of "
        f"{c['corner_turn_deg']:.0f}° rounded to 7.9 mm, a {c['length_mm']:.0f} mm path (m = 2); {c['runs']} planar "
        "runs, "
        f"one facing pair each, handing over by polarity swaps; gap ×4",
    )
    _facing_pairs(ax, rl, HA / MM * 4, 2 * HA / MM)
    cl = np.array([rl.points[rl.piece == i].mean(0) / MM for i in range(c["sticks"])])
    ax.scatter(cl[:, 0], cl[:, 1], cl[:, 2], marker="o", s=14, color=INK, depthshade=False, zorder=9)
    hand = np.array([rl.points[rl.piece == -1 - i].mean(0) / MM for i in c["_handovers"]])
    ax.scatter(hand[:, 0], hand[:, 1], hand[:, 2], marker="D", s=26, color=GREEN, depthshade=False, zorder=9)
    ax.scatter([], [], marker="D", color=GREEN, label=f"polarity-swap handover ({len(hand)})")
    ax.scatter([], [], marker="o", color=INK, label=f"corner in its run's plane ({c['sticks']})")
    ax.legend(loc="upper left", fontsize=6.9, bbox_to_anchor=(0.0, 0.9))
    # (d) T915-S
    d = res["T915-S"]
    knot: ClosedCurve = d["_curve"]
    fr: Framing = d["_framing"]
    pk = knot.points / MM
    gap_draw = 4.0
    ax = _axes3d(fig, grid[1, 0], limit, "")
    _title3d(
        ax,
        f"T915-S: facing flex strips {d['width_mm']:.2f} mm wide on the r/R = 0.45 trefoil, {d['length_mm']:.0f} mm "
        f"(m = 1), framed with k = 0 ({d['twist_deg_per_mm']:+.2f}°/mm), in-plane radius ≥ "
        f"{d['inplane_radius_mm']:.1f} mm; "
        f"two closing joints a quarter wave apart, 24 spacers; gap 1.83 mm ×{gap_draw:g}",
    )
    half = WIDTH915 / MM / 2
    offset = H915 / MM * gap_draw
    _split_ribbon(ax, pk + offset * fr.normal, fr.across, half, 0.27 * gap_draw, (ORANGE, BLUE), 0.5)
    _split_ribbon(ax, pk - offset * fr.normal, fr.across, half, 0.27 * gap_draw, (BLUE, ORANGE), 0.5)
    ax.plot(pk[:, 0], pk[:, 1], pk[:, 2], color=INK, lw=0.9)
    n = len(pk)
    start = int(np.argmin(knot.curvature))
    for k, sign in ((start, +1), ((start + int(QUARTER915 / L915 * n)) % n, -1)):
        span = np.arange(k - int(SCARF / L915 * n / 2), k + int(SCARF / L915 * n / 2)) % n
        line = pk[span] + sign * offset * fr.normal[span]
        ax.plot(line[:, 0], line[:, 1], line[:, 2], color=INK, lw=3.0, zorder=9)
    sp = (np.arange(SPACERS) * n // SPACERS + start) % n
    ax.scatter(pk[sp, 0], pk[sp, 1], pk[sp, 2], s=9, color=INK2, depthshade=False, zorder=8)
    ax.plot([], [], color=INK, lw=3, label="closing joint, scarfed 29.4 mm (top; bottom 79.5 mm on)")
    ax.scatter([], [], s=9, color=INK2, label="spacer station (24)")
    ax.plot([], [], color=ORANGE, lw=4, label="+V/2 half")
    ax.plot([], [], color=BLUE, lw=4, label="−V/2 half")
    ax.legend(loc="upper left", fontsize=6.9, bbox_to_anchor=(0.0, 0.9))
    # (e) T915-R: one ribbon carrying rails
    ax = _axes3d(fig, grid[1, 1], limit, "")
    _title3d(
        ax,
        "T915-R: one flex ribbon carrying rails along the same trefoil, framed the same way: rails 1.46 mm wide either "
        "side of a 0.40 mm grounded strip, the null 0.575 mm above (drawn ×4); one piece, one closing joint",
    )
    _ribbon(ax, pk, fr.across, 3.0, FLEX, 0.45)
    for side in (-1, 1):
        _ribbon(ax, pk + side * fr.across * (0.2 + 0.73), fr.across, 0.73, COPPER, 0.8)
    null = pk + 0.575 * 4 * fr.normal
    ax.plot(null[:, 0], null[:, 1], null[:, 2], color=INK, lw=0.9)
    span = np.arange(start - 40, start + 40) % n
    ax.plot(pk[span, 0], pk[span, 1], pk[span, 2], color=INK, lw=3.0)
    ax.plot([], [], color=COPPER, lw=4, label="rails")
    ax.plot([], [], color=INK, lw=0.9, label="the null above the ribbon")
    ax.plot([], [], color=INK, lw=3, label="closing joint")
    ax.legend(loc="upper left", fontsize=6.9, bbox_to_anchor=(0.0, 0.9))
    # (f) R1: one board, the channel ring, its island on posts to a carrier
    e = res["R1"]
    rad = e["radius_mm"]
    ax = _axes3d(fig, grid[1, 2], limit, "")
    _title3d(
        ax,
        f"R1: one board, the channel ring {e['length_mm']:.0f} mm round (m = 1 on RO4003C, radius {rad:.1f} mm), "
        "copper "
        f"on both faces either side of the channel; its island held by {e['posts']} posts to a carrier board below "
        f"({e['posts']} does not divide 2m); thickness ×8",
    )
    ring_t = np.linspace(0, 2 * math.pi, 241)
    _slab(ax, square, 0.0, 0.813 * 8, LAMINATE, 0.13)
    island = np.column_stack([(rad - 0.5) * np.cos(ring_t[:-1]), (rad - 0.5) * np.sin(ring_t[:-1])])
    _slab(ax, island, 0.0, 0.813 * 8, LAMINATE, 0.2)
    _slab(ax, square * 0.9, -14.0, -12.4, LAMINATE, 0.12)
    for zz in (0.0, 0.813 * 8):
        for rr in (rad - 0.9, rad + 0.9):
            ax.plot(rr * np.cos(ring_t), rr * np.sin(ring_t), np.full_like(ring_t, zz), color=COPPER, lw=1.5)
    ax.plot(rad * np.cos(ring_t), rad * np.sin(ring_t), np.full_like(ring_t, 0.813 * 4), color=INK, lw=1.0)
    for k in range(3):
        for rr, phase in ((rad * 0.55, 0.0), (rad + 9.0, math.pi / 3)):
            ang = 2 * math.pi * k / 3 + phase
            ax.plot([rr * math.cos(ang)] * 2, [rr * math.sin(ang)] * 2, [-12.4, 0.0], color=INK2, lw=2.0)
    feeds = [0.0, math.pi / 2]
    ax.scatter(
        [rad * math.cos(f) for f in feeds],
        [rad * math.sin(f) for f in feeds],
        [0.813 * 4] * 2,
        marker="s",
        s=28,
        color=RED,
    )
    ax.plot([], [], color=COPPER, lw=1.5, label="copper either side of the channel")
    ax.plot([], [], color=INK2, lw=2.0, label="posts to the carrier (3 + 3)")
    ax.scatter([], [], marker="s", color=RED, label="two feeds a quarter wave apart")
    ax.legend(loc="upper left", fontsize=6.9, bbox_to_anchor=(0.0, 0.9))
    _suptitle(
        fig,
        "Boards that carry a trefoil, at one scale (mm): two slotted boards joined by six ramps, rounded sticks, "
        "lattice "
        "boards, facing flex strips framed with k = 0 and one flex ribbon carrying rails, beside the one-board ring",
    )
    _save(fig, out, "arrangements")


def _facing_pairs(ax: Any, r: Rounded, offset: float, half: float) -> None:
    """Facing boards along a rounded path: the two boards at +-offset along the normals, the twisted handovers green."""
    p = r.points / MM
    across = np.cross(r.tangents, r.normals)
    edges = np.flatnonzero(np.diff(r.twisting.astype(np.int8)) != 0) + 1
    pieces = np.split(np.arange(len(p)), edges)
    for idx in pieces:
        if len(idx) < 2:
            continue
        twisted = bool(r.twisting[idx[0]])
        for sign, colour in ((1.0, BLUE), (-1.0, ORANGE)):
            _ribbon(
                ax,
                p[idx] + sign * offset * r.normals[idx],
                across[idx],
                half,
                GREEN if twisted else colour,
                0.5,
                step=2,
            )
    ax.plot(p[:, 0], p[:, 1], p[:, 2], color=INK, lw=0.8)


# ---- knots (figure 2) ---------------------------------------------------------------------------------------------
ROPELENGTH = {"3₁": 16.37, "4₁": 21.04, "5₁": 23.60, "8₁₉": 30.49}


def layout_numbers(curve: ClosedCurve, sample: int = 720) -> dict[str, Any]:
    """A layout's numbers from its projection on the xy plane and its heights: crossings, the determinant, ramps
    (where the height changes sign) and the least the diagram needs, and its rotational symmetry."""
    pts = curve.points[:: max(len(curve) // sample, 1)]
    cross = crossings(pts)
    z = pts[:, 2]
    return {
        "crossings": len(cross),
        "determinant": determinant(cross),
        "ramps": int(np.sum((z > 0) != (np.roll(z, 1) > 0))),
        "least_ramps": ramps_needed(cross),
        "symmetry": rotational_order(pts),
        "_points": pts,
        "_crossings": cross,
    }


def compute_knots(quick: bool) -> dict[str, Any]:
    say("\n== knots: layouts on two levels, energy against joints and bends, length against W/D ==")
    n = 1200 if quick else 2400
    layouts = {
        "3₁ as T(3, 2): torus_knot(2, 3)": torus(0.45, 1.0, n, 2, 3),
        "3₁ as T(2, 3): torus_knot(3, 2)": torus(0.45, 1.0, n, 3, 2),
        "4₁, common form": figure_eight_knot(1.0, n),
        "5₁ = T(5, 2): torus_knot(2, 5)": torus(0.45, 1.0, n, 2, 5),
        "5₁ = T(2, 5): torus_knot(5, 2)": torus(0.45, 1.0, n, 5, 2),
        "8₁₉ = T(4, 3): torus_knot(3, 4)": torus(0.45, 1.0, n, 3, 4),
        "8₁₉ = T(3, 4): torus_knot(4, 3)": torus(0.45, 1.0, n, 4, 3),
    }
    out: dict[str, Any] = {"layouts": {}}
    for name, curve in layouts.items():
        numbers = layout_numbers(curve)
        out["layouts"][name] = numbers
        say(
            f"  {name}: {numbers['crossings']} crossings, determinant {numbers['determinant']}, {numbers['ramps']} "
            f"ramps (least for this diagram {numbers['least_ramps']}), C{numbers['symmetry']}"
        )
    # flat points: the trefoil's curvature on its inner equator vanishes at r/R = 4/13
    ratios = np.linspace(0.20, 0.45, 26 if quick else 101)
    floor = np.array([float(torus(float(rr), 1.0, 1500).curvature.min()) for rr in ratios])
    out["flat"] = {"ratio": ratios, "least_curvature": floor, "predicted": 4 / 13}
    say(
        f"  the trefoil's least curvature (L = 1) is smallest at r/R = {ratios[np.argmin(floor)]:.3f}; 4/13 = "
        f"{4 / 13:.4f}"
    )
    # sudden turns and sticks
    need = 16.6
    sticks = next(m for m in range(6, 200) if 1 / math.sin(math.radians(720 / m)) ** 2 >= need)
    out["sticks_for"] = {"W_over_D": need, "sticks": sticks}
    say(f"  a sudden joint keeps W sin^2 theta < D: holding {need} D takes {sticks} sticks or more")
    # designs: length, tightest bend, h, and the W/D each holds by the bend rule
    designs = design_table()
    out["designs"] = designs
    for name, row in designs.items():
        say(
            f"  {name}: {row['length_mm']:.0f} mm, tightest bend {row['bend_mm']:.1f} mm, h {row['h_mm']:.3f} mm: "
            f"rule R/2h = {row['rule']:.1f}, L / (8 pi h W/D) = "
            f"{row['length_mm'] / (8 * math.pi * row['h_mm'] * row['rule']):.2f}"
        )
    out["per_unit_mm"] = 8 * math.pi * 0.5
    out["facing_per_unit_mm"] = 4 * math.pi * 0.915 / 0.605
    say(
        f"  L > 8 pi h W/D: {out['per_unit_mm']:.1f} mm per unit of W/D at h = 0.5 mm; facing boards' static limit "
        f"{out['facing_per_unit_mm']:.1f} mm per unit at h = 0.915 mm"
    )
    return out


def design_table() -> dict[str, dict[str, float]]:
    """The lead designs' lengths, tightest bends and h (mm), and the W/D the bend rule R / 2h gives each."""
    s2c3 = two_level_trefoil(33.2e-3, 4e-3, 1.813e-3, 16e-3, 3000)[0]
    flat = ClosedCurve.from_polygon(s2c3.points * np.array([1.0, 1.0, 0.0]), len(s2c3))
    s2 = (s2c3.length / MM, flat.bend_radius / MM)  # the levels' bend; the ramps' own limit is section 3.3's
    first = first_run_trefoil()
    t915 = torus(RATIO915, L915, 3000)
    rows = {
        "first run (four rods)": (first.length / MM, first.bend_radius / MM, 0.5),
        "T915-S": (L915 / MM, t915.bend_radius / MM, H915 / MM),
        "S2-C3": (s2[0], s2[1], 0.5),
        "rounded sticks": (291.0, 8.0, HA / MM),
        "L26": (582.0, 7.9, HA / MM),
        "R1 channel ring": (210.0, 210.0 / (2 * math.pi), 0.5),
        "flat ring control": (L915 / MM, L915 / MM / (2 * math.pi), H915 / MM),
    }
    return {
        name: {"length_mm": length, "bend_mm": bend, "h_mm": h, "rule": bend / (2 * h)}
        for name, (length, bend, h) in rows.items()
    }


def figure_knots(res: dict[str, Any], out: Path) -> None:
    """Figure 2: the knots' layouts on two levels; energy held against joint angle and bend radius; length against
    W/D with the bound, the ropelength floors and each design."""
    plt = _pyplot()
    fig = plt.figure(figsize=(15.5, 10.2))
    outer = fig.add_gridspec(2, 1, height_ratios=[0.62, 1.0], hspace=0.16)
    top = outer[0].subgridspec(1, 7, wspace=0.08)
    for k, (name, lay) in enumerate(res["layouts"].items()):
        ax = fig.add_subplot(top[0, k])
        _draw_layout(ax, lay)
        label, _, form = name.partition(": ")
        ax.set_title(
            f"{label}\n{form or 'its common form'}\n{lay['crossings']} crossings, det {lay['determinant']}\n"
            f"{lay['ramps']} ramps, C{lay['symmetry']}",
            fontsize=8.2,
            loc="center",
        )
    bottom = outer[1].subgridspec(1, 3, wspace=0.3)
    # energy against joint angle: the sudden turn
    ax = fig.add_subplot(bottom[0, 0])
    theta = np.linspace(1.0, 179.0, 400)
    s2 = np.sin(np.radians(theta)) ** 2
    sudden = np.where(theta < 90, 1 / s2, 1.0)
    ax.semilogy(theta, sudden, color=INK, label="sudden turn: W·sin²θ < D below 90°, W < D beyond")
    upper = np.where(theta < 90, 1 / (0.4 * s2), 1 / 0.4)
    lower = np.where(theta < 90, 1 / (0.75 * s2), 1 / 0.75)
    ax.fill_between(theta, lower, upper, color=BLUES[1], alpha=0.8, label="a kink handing over 0.4 to 0.75 of W·sin²θ")
    sticks = res["sticks_for"]["sticks"]
    th_s = 720 / sticks
    held_s = 1 / math.sin(math.radians(th_s)) ** 2
    ax.plot([th_s], [held_s], "o", color=ORANGE)
    ax.annotate(
        f"{sticks} sticks: every turn ≤ {th_s:.1f}°,\nholding {res['sticks_for']['W_over_D']:g} D",
        (th_s, held_s),
        xytext=(28, 34),
        textcoords="offset points",
        fontsize=8,
        arrowprops={"arrowstyle": "-", "color": INK2, "lw": 0.7},
    )
    ax.plot([128.8], [1.0], "s", color=RED)
    ax.annotate(
        "six sticks: 128.8°, W < D",
        (128.8, 1.0),
        xytext=(-52, 30),
        textcoords="offset points",
        fontsize=8,
        arrowprops={"arrowstyle": "-", "color": INK2, "lw": 0.7},
    )
    ax.plot([90.0], [1.0], "D", color=GREEN)
    ax.annotate("a lattice corner, 90°", (90.0, 1.0), xytext=(-70, -22), textcoords="offset points", fontsize=8)
    ax.set_xlabel("turn of a sharp joint θ (degrees)")
    ax.set_ylabel("energy held over the well's depth, W/D")
    ax.set_ylim(0.5, 400)
    ax.set_xlim(0, 180)
    ax.legend(fontsize=7.8, loc="upper right")
    ax.set_title("a sharp joint (KNT-4)")
    # energy against bend radius
    ax = fig.add_subplot(bottom[0, 1])
    rh = np.geomspace(4, 200, 200)
    ax.loglog(rh, rh / 2, color=INK, label="bend rule R/2h")
    ax.loglog(rh, 0.605 * rh, color=ORANGE, label="facing boards in their plane: 0.605·R/h")
    ax.loglog(rh, 0.350 * rh, color=BLUE, label="four rods: 0.350·R/r₀")
    ax.fill_between(
        rh,
        0.38 * rh / 2,
        0.77 * rh / 2,
        color=BLUES[1],
        alpha=0.8,
        label="tracked on four rods: 0.38 to 0.77 of the rule",
    )
    for k, (name, row) in enumerate(res["designs"].items()):
        x = row["bend_mm"] / row["h_mm"]
        ax.plot([x], [x / 2], "o", color=SERIES[(k + 2) % len(SERIES)], ms=6)
        ax.annotate(name, (x, x / 2), xytext=(5, -11 if k % 2 else 4), textcoords="offset points", fontsize=7.4)
    ax.set_xlabel("tightest bend over the distance to the copper, R/h")
    ax.set_ylabel("W/D")
    ax.legend(fontsize=7.6, loc="upper left")
    ax.set_title("a smooth bend (static limits, JNT-2c)")
    # length against W/D
    ax = fig.add_subplot(bottom[0, 2])
    wd = np.geomspace(1, 60, 200)
    ax.loglog(wd, 8 * math.pi * wd, color=INK, label="a knot: L > 4πR = 8πh·W/D (by the rule)")
    ax.loglog(wd, 4 * math.pi * wd, color=INK, ls="--", lw=1.0, label="a ring: L = 2πR = 4πh·W/D")
    ax.loglog(wd, 4 * math.pi * wd / 0.605, color=ORANGE, label="a knot of facing boards: 4πh·W/(0.605·D)")
    ax.fill_between(
        wd,
        8 * math.pi * wd / 0.77,
        8 * math.pi * wd / 0.38,
        color=BLUES[1],
        alpha=0.8,
        label="a knot at the tracked fraction",
    )
    for knot, rope in ROPELENGTH.items():
        ax.axhline(2 * rope, color=MUTED, lw=0.8, ls=":")
        ax.text(
            58, 2 * rope * 1.04, f"{knot} tight in a tube 2h wide: {2 * rope:.0f} h", fontsize=7, color=INK2, ha="right"
        )
    for k, (name, row) in enumerate(res["designs"].items()):
        x, y = row["rule"], row["length_mm"] / row["h_mm"]
        ax.plot([x], [y], "o", color=SERIES[(k + 2) % len(SERIES)], ms=6)
        ax.annotate(name, (x, y), xytext=(5, 4 if k % 2 else -11), textcoords="offset points", fontsize=7.4)
    ax.set_xlabel("W/D by the bend rule, R/2h")
    ax.set_ylabel("length over h, L/h")
    ax.set_ylim(20, 3000)
    ax.legend(fontsize=7.4, loc="upper left")
    ax.set_title(f"length (KNT-3): {res['per_unit_mm']:.1f} mm per unit of W/D at h = 0.5 mm")
    _suptitle(
        fig,
        "What knot theory fixes: on two levels a knot of bridge number b takes 2b ramps or more (each layout here from "
        "ClosedCurve with its crossings, determinant, ramps and symmetry); a sharp joint holds W < D/sin²θ, so "
        f"{sticks} sticks hold the first run's 16.6 D; and a knot's length exceeds 8πh·W/D",
    )
    _save(fig, out, "knots")


def _draw_layout(ax: Any, lay: dict[str, Any]) -> None:
    """A layout seen from above: the upper level dark, the lower light, ramps as triangles, the lower strand broken
    at each crossing."""
    p = lay["_points"]
    n = len(p)
    under = np.zeros(n, dtype=bool)
    for c in lay["_crossings"]:
        k = int(c.under)
        under[[(k + j) % n for j in range(-4, 5)]] = True
    upper = p[:, 2] > 0
    for level, colour, width in ((upper, INK, 1.6), (~upper, MUTED, 1.6)):
        y = np.where(level & ~under, p[:, 1], np.nan)
        ax.plot(np.append(p[:, 0], p[0, 0]), np.append(y, y[0]), color=colour, lw=width, solid_capstyle="round")
    ramps = np.flatnonzero((p[:, 2] > 0) != (np.roll(p[:, 2], 1) > 0))
    ax.plot(p[ramps, 0], p[ramps, 1], "^", color=ORANGE, ms=4.5)
    ax.set_aspect("equal")
    ax.axis("off")


# ---- strips (figure 4) --------------------------------------------------------------------------------------------
def compute_strips(quick: bool) -> dict[str, Any]:
    say("\n== strips: curvature, torsion, framings and plates ==")
    n = 1500 if quick else 4000
    out: dict[str, Any] = {"curves": {}}
    for name, curve in (
        ("r/R = 0.30, 276.5 mm (the first run)", first_run_trefoil(n)),
        ("r/R = 0.45, 318 mm (T915-S)", torus(RATIO915, L915, n)),
    ):
        kappa, tau = curvature_torsion(curve)
        ratio = tau / kappa
        fold = 1 / np.abs(periodic_gradient(ratio, curve.length / len(curve)))
        out["curves"][name] = {
            "length_mm": curve.length / MM,
            "bend_mm": curve.bend_radius / MM,
            "total_curvature_pi": curve.total_curvature / math.pi,
            "holonomy_deg": math.degrees(curve.holonomy),
            "writhe": curve.writhe,
            "fold_min_mm": float(fold.min()) / MM,
            "_s": np.arange(len(curve)) / len(curve),
            "_kappa": kappa,
            "_tau": tau,
            "_fold": fold,
        }
        say(
            f"  {name}: tightest bend {curve.bend_radius / MM:.2f} mm, total curvature "
            f"{curve.total_curvature / math.pi:.2f} pi, "
            f"holonomy {math.degrees(curve.holonomy):.1f} deg, writhe {curve.writhe:.3f}; a straight strip folds "
            f"{float(fold.min()) / MM:.2g} mm from the curve"
        )
    knot = torus(RATIO915, L915, n)
    out["framings"] = {}
    for k in (0, -1, -3):
        fr = framing(knot, k)
        kappa = np.linalg.norm(curvature_vector(knot), axis=1)
        tops = np.flatnonzero((kappa >= np.roll(kappa, 1)) & (kappa > np.roll(kappa, -1)))
        tops = tops[np.argsort(-kappa[tops])][:6]
        row: dict[str, Any] = {
            "twist_deg_per_mm": math.degrees(fr.twist) * MM,
            "inplane_radius_mm": fr.inplane_radius / MM,
            "zeros_of_kappa_n": int(np.sum(np.sign(fr.kappa_n) != np.sign(np.roll(fr.kappa_n, 1)))),
            "toward_boards_at_maxima": float(np.mean(np.abs(fr.kappa_n[tops]) / kappa[tops])),
            "integral_kappa_n": float(fr.kappa_n.sum() * knot.length / len(knot)),
            "_framing": fr,
        }
        for name, offset in (("plus", +H915), ("minus", -H915)):
            dl, kg = plate(knot, fr.normal, offset)
            xy, heading = develop(dl, kg)
            angle, short, long = smallest_box(outline(xy, heading, WIDTH915)[::4])
            row[f"plate_{name}_mm"] = float(dl.sum()) / MM
            row[f"box_{name}_mm"] = (short / MM, long / MM)
            row[f"turning_{name}_deg"] = math.degrees(float((kg * dl).sum()))
            row[f"_flat_{name}"] = (xy, heading, angle)
        b = WIDTH915 / 2
        row["edge_strain"] = (fr.twist * b) ** 2 / 3
        row["centre_strain"] = -((fr.twist * b) ** 2) / 6
        row["wrinkle_scale"] = (43e-6 / WIDTH915) ** 2
        out["framings"][k] = row
        say(
            f"  k = {k:+d}: twist {row['twist_deg_per_mm']:+.3f} deg/mm, in-plane radius >= "
            f"{row['inplane_radius_mm']:.1f} mm, "
            f"kappa_n changes sign {row['zeros_of_kappa_n']} times a lap, {row['toward_boards_at_maxima']:.2f} of the "
            f"curvature maxima toward the boards; plates {row['plate_plus_mm']:.1f} and {row['plate_minus_mm']:.1f} "
            "mm, "
            f"boxes {row['box_plus_mm'][0]:.0f} x {row['box_plus_mm'][1]:.0f} and {row['box_minus_mm'][0]:.0f} x "
            f"{row['box_minus_mm'][1]:.0f} mm; edges strain {row['edge_strain']:+.1e}, centre "
            f"{row['centre_strain']:+.1e}"
        )
    return out


def signed(k: int) -> str:
    """An integer with a true minus sign."""
    return f"−{-k}" if k < 0 else f"{k}"


def minus(value: float, spec: str) -> str:
    """A number formatted to spec with a true minus sign."""
    return format(value, spec).replace("-", "−")


def figure_strips(res: dict[str, Any], out: Path) -> None:
    """Figure 4: along the r/R = 0.30 and 0.45 trefoils, curvature, torsion and the straight-strip limit; the k = -1
    and k = 0 strips' curvature in their plane and toward the boards; the plates laid flat on a 234 x 490 mm panel."""
    plt = _pyplot()
    fig = plt.figure(figsize=(15.5, 13.6))
    grid = fig.add_gridspec(3, 3, hspace=0.42, wspace=0.26, height_ratios=[1.0, 1.0, 1.15])
    axes = [fig.add_subplot(grid[0, k]) for k in range(3)]
    for k, (name, row) in enumerate(res["curves"].items()):
        colour = SERIES[k]
        s = row["_s"]
        axes[0].plot(s, row["_kappa"] * MM, color=colour, label=f"{name}: tightest bend {row['bend_mm']:.1f} mm")
        axes[1].plot(s, np.clip(row["_tau"] * MM, -0.2, 0.2), color=colour, lw=1.3, label=name)
        axes[2].semilogy(
            s, np.clip(row["_fold"] / MM, 1e-5, 1e5), color=colour, label=f"{name}: least {sci(row['fold_min_mm'])} mm"
        )
    axes[0].set_ylabel("curvature (1/mm)")
    axes[0].set_title("curvature along the trefoils")
    axes[1].set_ylabel("torsion (1/mm; clipped at ±0.2)")
    axes[1].set_title("torsion: near a flat point (r/R near 4/13) it spikes")
    axes[2].axhline(WIDTH915 / 2 / MM, color=MUTED, lw=0.9, ls=":")
    axes[2].text(0.99, WIDTH915 / 2 / MM * 1.4, "T915-S's half-width, 3.66 mm", fontsize=7.6, color=INK2, ha="right")
    axes[2].set_ylabel("1/|d(τ/κ)/ds| (mm)")
    axes[2].set_title("where a straight (rectifying) strip folds")
    for ax in axes:
        ax.set_xlabel("arc length, fraction of a lap")
        ax.legend(fontsize=7.4, loc="lower right" if ax is axes[2] else "upper right")
    for col, k in enumerate((0, -1)):
        ax = fig.add_subplot(grid[1, col if col == 0 else 1])
        row = res["framings"][k]
        fr: Framing = row["_framing"]
        s = np.arange(len(fr.kappa_g)) / len(fr.kappa_g)
        ax.plot(
            s, fr.kappa_g * MM, color=BLUE, label=f"in the boards' plane: radius ≥ {row['inplane_radius_mm']:.1f} mm"
        )
        ax.plot(
            s, fr.kappa_n * MM, color=ORANGE, label=f"toward the boards: {row['zeros_of_kappa_n']} sign changes a lap"
        )
        ax.axhline(0, color=EDGE, lw=0.8)
        lim = 1.45 * float(max(np.abs(fr.kappa_g).max(), np.abs(fr.kappa_n).max())) * MM
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("arc length, fraction of a lap")
        ax.set_ylabel("curvature (1/mm)")
        ax.legend(fontsize=7.6, loc="upper center", ncol=1)
        ax.set_title(f"T915-S's strips framed with k = {signed(k)}: twist {minus(row['twist_deg_per_mm'], '+.3f')}°/mm")
    ax = fig.add_subplot(grid[1, 2])
    ks = sorted(res["framings"])
    rows = [res["framings"][k] for k in ks]
    xpos = np.arange(len(ks))
    ax.bar(xpos - 0.2, [r["inplane_radius_mm"] for r in rows], width=0.38, color=BLUE, label="in-plane radius (mm)")
    ax.bar(
        xpos + 0.2,
        [abs(r["centre_strain"]) / r["wrinkle_scale"] for r in rows],
        width=0.38,
        color=ORANGE,
        label="the twist's compression at the centre over the wrinkling scale (t/W)²",
    )
    for x, r in zip(xpos, rows, strict=True):
        ax.text(x - 0.2, r["inplane_radius_mm"] + 0.6, f"{r['inplane_radius_mm']:.1f}", ha="center", fontsize=7.6)
        ax.text(
            x + 0.2,
            abs(r["centre_strain"]) / r["wrinkle_scale"] + 0.6,
            f"{abs(r['centre_strain']) / r['wrinkle_scale']:.1f}",
            ha="center",
            fontsize=7.6,
        )
    ax.set_xticks(xpos)
    ax.set_xticklabels(
        [f"k = {signed(k)}\n{minus(r['twist_deg_per_mm'], '+.3f')}°/mm" for k, r in zip(ks, rows, strict=True)]
    )
    ax.set_ylim(0, 52)
    ax.legend(fontsize=7.4, loc="upper right")
    ax.set_title("the closures: what each framing keeps and what its twist strains")
    ax5 = fig.add_subplot(grid[2, :])
    w, h = PANEL_MM
    ax5.add_patch(_rect(0, 0, h, w, PALE))
    placed = [(0, "plus", 8.0, 30.0), (0, "minus", 8.0, 75.0), (-1, "plus", math.nan, 112.0)]
    for k, name, x0, y0 in placed:
        row = res["framings"][k]
        xy, heading, angle = row[f"_flat_{name}"]
        shape = outline(xy, heading, WIDTH915) / MM
        c, s_ = math.cos(angle), math.sin(angle)
        turned = np.column_stack([shape[:, 0] * c - shape[:, 1] * s_, shape[:, 0] * s_ + shape[:, 1] * c])
        if np.ptp(turned[:, 1]) > np.ptp(turned[:, 0]):
            turned = turned[:, ::-1]
        turned -= turned.min(0)
        if not math.isfinite(x0):  # right-aligned in the panel
            x0 = h - 6.0 - float(np.ptp(turned[:, 0]))
        turned += [x0, y0]
        colour = BLUE if k == 0 else ORANGE
        ax5.fill(turned[:, 0], turned[:, 1], color=colour, alpha=0.35, lw=0)
        ax5.plot(turned[:, 0], turned[:, 1], color=colour, lw=0.8)
        box = row[f"box_{name}_mm"]
        label = (
            f"k = {signed(k)}, plate {'+h' if name == 'plus' else '−h'}: {row[f'plate_{name}_mm']:.1f} mm long, "
            f"{box[0]:.0f} × {box[1]:.0f} mm box"
        )
        ax5.text(turned[:, 0].min(), turned[:, 1].max() + 4, label, fontsize=8, color=INK)
    k3 = res["framings"][-3]
    ax5.text(
        8,
        w - 14,
        f"k = −3 keeps C3 with plates of equal length ({k3['plate_plus_mm']:.1f} mm) as one band of twice the length, "
        f"{2 * k3['plate_plus_mm']:.0f} mm: longer than the panel",
        fontsize=8,
        color=INK2,
    )
    ax5.set_xlim(-5, h + 5)
    ax5.set_ylim(-5, w + 5)
    ax5.set_aspect("equal")
    ax5.set_xlabel("mm")
    ax5.set_title(
        f"the plates laid flat on the fab's {w:.0f} × {h:.0f} mm panel: k = 0's two plates, and one lap of k = −1"
    )
    k0 = res["framings"][0]
    _suptitle(
        fig,
        "Strips along the trefoil: framed with k = 0 the strips keep C3, an in-plane radius ≥ "
        f"{k0['inplane_radius_mm']:.1f} mm and a "
        f"twist of {k0['twist_deg_per_mm']:+.2f}°/mm, and their two plates differ in length, {k0['plate_plus_mm']:.1f} "
        "and "
        f"{k0['plate_minus_mm']:.1f} mm, since ∮κ_n ds = {minus(k0['integral_kappa_n'], '.1f')} rad",
    )
    _save(fig, out, "strips")


def _rect(x: float, y: float, w: float, h: float, colour: str) -> Any:
    from matplotlib.patches import Rectangle

    return Rectangle((x, y), w, h, facecolor=colour, edgecolor=EDGE, lw=0.8)


# ---- solvers: Laplace on grids, bounded least squares -------------------------------------------------------------
def _pcg(
    apply: Callable[[Arr], Arr], b: Arr, diag: Arr, free: BoolArr, tol: float, iters: int, axes: tuple[int, ...]
) -> Arr:
    """Jacobi-preconditioned conjugate gradients for A x = b on the free nodes, a batch of problems at once (the axes
    listed are one problem's grid; the leading axes the batch)."""
    x = np.zeros_like(b)
    r = np.where(free, b, 0.0)
    z = r / diag
    p = z.copy()
    rz = (r * z).sum(axis=axes, keepdims=True)
    norm_b = np.sqrt((r * r).sum(axis=axes, keepdims=True)) + 1e-300
    for _ in range(iters):
        ap = np.where(free, apply(p), 0.0)
        alpha = rz / ((p * ap).sum(axis=axes, keepdims=True) + 1e-300)
        x += alpha * p
        r -= alpha * ap
        if bool(np.all(np.sqrt((r * r).sum(axis=axes, keepdims=True)) < tol * norm_b)):
            break
        z = r / diag
        rz_new = (r * z).sum(axis=axes, keepdims=True)
        p = z + (rz_new / (rz + 1e-300)) * p
        rz = rz_new
    return x


def laplace3d(fixed: BoolArr, value: Arr, tol: float = 1e-8, iters: int = 40000) -> Arr:
    """The potential on a uniform grid (the last three axes; leading axes a batch): value where fixed, the discrete
    Laplace equation elsewhere, and zero normal derivative on the box's faces where they are not fixed (the grid
    graph's Laplacian, symmetric, by conjugate gradients)."""
    fixed = np.broadcast_to(fixed, value.shape)
    free = ~fixed
    shape = value.shape

    def lap(u: Arr) -> Arr:
        r = np.zeros_like(u)
        for ax in (-3, -2, -1):
            d = np.diff(u, axis=ax)
            lo = [slice(None)] * u.ndim
            hi = [slice(None)] * u.ndim
            lo[ax], hi[ax] = slice(None, -1), slice(1, None)
            r[tuple(lo)] -= d
            r[tuple(hi)] += d
        return r

    deg = np.zeros(shape[-3:])
    for ax in range(3):
        n = shape[-3 + ax]
        count = np.full(n, 2.0)
        count[0] = count[-1] = 1.0
        deg = deg + count.reshape([-1 if k == ax else 1 for k in range(3)])
    u0 = np.where(fixed, value, 0.0)
    b = -lap(u0)
    x = _pcg(lap, b, np.broadcast_to(deg, shape), free, tol, iters, (-3, -2, -1))
    out: Arr = np.where(fixed, value, x)
    return out


def edge_weights(x: Arr, y: Arr, eps: Arr) -> tuple[Arr, Arr]:
    """Finite-volume weights of a tensor grid's edges: permittivity of the cells beside an edge times their share of
    its control face, over the edge's length. eps (nx - 1, ny - 1) holds each cell's relative permittivity."""
    hx, hy = np.diff(x), np.diff(y)
    ep = np.pad(eps, ((0, 0), (1, 1)))
    hp = np.pad(hy, (1, 1))
    wx = (ep[:, :-1] * hp[None, :-1] / 2 + ep[:, 1:] * hp[None, 1:] / 2) / hx[:, None]
    epx = np.pad(eps, ((1, 1), (0, 0)))
    hpx = np.pad(hx, (1, 1))
    wy = (epx[:-1, :] * hpx[:-1, None] / 2 + epx[1:, :] * hpx[1:, None] / 2) / hy[None, :]
    return wx, wy


def laplace2d(x: Arr, y: Arr, eps: Arr, fixed: BoolArr, value: Arr, tol: float = 1e-10, iters: int = 60000) -> Arr:
    """The potential on a tensor grid's nodes (nx, ny) with cells of relative permittivity eps (nx - 1, ny - 1):
    value where fixed, Gauss's law by finite volumes elsewhere, zero normal flux on faces that are not fixed."""
    wx, wy = edge_weights(x, y, eps)

    def apply(u: Arr) -> Arr:
        r = np.zeros_like(u)
        fx = wx * (u[1:, :] - u[:-1, :])
        r[:-1, :] -= fx
        r[1:, :] += fx
        fy = wy * (u[:, 1:] - u[:, :-1])
        r[:, :-1] -= fy
        r[:, 1:] += fy
        return r

    diag = np.zeros_like(value)
    diag[:-1, :] += wx
    diag[1:, :] += wx
    diag[:, :-1] += wy
    diag[:, 1:] += wy
    free = ~fixed
    u0 = np.where(fixed, value, 0.0)
    xs = _pcg(apply, -apply(u0), np.where(diag > 0, diag, 1.0), free, tol, iters, (-2, -1))
    out: Arr = np.where(fixed, value, xs)
    return out


def energy2d(x: Arr, y: Arr, eps: Arr, phi: Arr) -> float:
    """The integral of eps |grad phi|^2 over the grid, by its edges: twice the energy per length in eps0 units."""
    wx, wy = edge_weights(x, y, eps)
    return float((wx * np.diff(phi, axis=0) ** 2).sum() + (wy * np.diff(phi, axis=1) ** 2).sum())


def nodes(lo: float, hi: float, features: Sequence[float], fine: float, coarse: float, growth: float = 0.3) -> Arr:
    """Grid nodes from lo to hi, `fine` apart at the features (which are nodes), spacing growing by `growth` per unit
    distance from them up to `coarse`."""
    feats = np.array(sorted(f for f in features if lo <= f <= hi), dtype=np.float64)
    xs = [lo]
    while xs[-1] < hi:
        d = float(np.min(np.abs(feats - xs[-1]))) if len(feats) else math.inf
        xs.append(xs[-1] + min(coarse, fine + growth * d))
    arr = np.array(xs, dtype=np.float64)
    arr[-1] = hi
    for f in feats:
        arr[int(np.argmin(np.abs(arr - f)))] = f
    arr = np.unique(arr)
    keep = np.concatenate([[True], np.diff(arr) > fine / 4])
    keep |= np.isin(arr, feats)
    out: Arr = np.unique(arr[keep])
    return out


def bounded_lstsq(a: Arr, b: Arr, bound: float, iters: int = 30000) -> Arr:
    """argmin |a x - b|^2 with every |x_i| <= bound, by accelerated projected gradients (FISTA)."""
    step = 1 / float(np.linalg.norm(a, 2)) ** 2
    x = np.clip(np.linalg.lstsq(a, b, rcond=None)[0], -bound, bound)
    y, t = x.copy(), 1.0
    atb = a.T @ b
    ata = a.T @ a
    for _ in range(iters):
        x_new = np.clip(y - step * (ata @ y - atb), -bound, bound)
        t_new = (1 + math.sqrt(1 + 4 * t * t)) / 2
        y = x_new + ((t - 1) / t_new) * (x_new - x)
        if float(np.abs(x_new - x).max()) < 1e-11:
            x = x_new
            break
        x, t = x_new, t_new
    out: Arr = x
    return out


def zero_crossings(z: Arr, col: Arr) -> list[tuple[float, float]]:
    """Where col(z) changes sign, linearly interpolated, with its slope there."""
    s = np.sign(col)
    idx = np.flatnonzero(s[:-1] * s[1:] < 0)
    return [
        (
            float(z[j] - col[j] * (z[j + 1] - z[j]) / (col[j + 1] - col[j])),
            float((col[j + 1] - col[j]) / (z[j + 1] - z[j])),
        )
        for j in idx
    ]


def bilinear(grid_a: Arr, grid_b: Arr, values: Arr, a: Arr, b: Arr) -> Arr:
    """values (na, nb) on a uniform grid, interpolated at points (a, b); zero outside."""
    da, db = grid_a[1] - grid_a[0], grid_b[1] - grid_b[0]
    fa = (a - grid_a[0]) / da
    fb = (b - grid_b[0]) / db
    inside = (fa >= 0) & (fa <= len(grid_a) - 1) & (fb >= 0) & (fb <= len(grid_b) - 1)
    ia = np.clip(np.floor(fa).astype(np.intp), 0, len(grid_a) - 2)
    ib = np.clip(np.floor(fb).astype(np.intp), 0, len(grid_b) - 2)
    ta, tb = fa - ia, fb - ib
    v = (
        values[ia, ib] * (1 - ta) * (1 - tb)
        + values[ia + 1, ib] * ta * (1 - tb)
        + values[ia, ib + 1] * (1 - ta) * tb
        + values[ia + 1, ib + 1] * ta * tb
    )
    out: Arr = np.where(inside, v, 0.0)
    return out


# ---- ramp (figure 3) ----------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class RampBases:
    """The field on the axial plane of a stack of slotted layers (odd in x across the slot), for a 1 mm segment of
    each layer at one volt and every other copper grounded: g[l] (ns, nz) = dphi/dx on x = 0, whose zeros are the
    null line. Lengths in mm."""

    s: Arr
    z: Arr
    g: Arr
    layers: tuple[float, ...]
    slot: float


def ramp_bases(
    layers: Sequence[float], slot: float, step: float, margin: float = 1.5, width: float = 2.5, reach: float = 7.0
) -> RampBases:
    """Solve the four (or n) layers' unit segments at once in three dimensions: x from the slot's centre plane
    (grounded by symmetry) to `width`, s over +-reach, z between grounded shields `margin` beyond the outer layers."""
    zmax = max(abs(v) for v in layers) + margin
    x = np.arange(0, width + step / 2, step)
    s = np.arange(-reach, reach + step / 2, step)
    z = np.arange(-zmax, zmax + step / 2, step)
    gx, gs, gz = np.meshgrid(x, s, z, indexing="ij")
    fixed = (gx == 0) | (np.abs(gz) >= zmax - step / 4)
    copper = [(gx >= slot / 2 - step / 4) & (np.abs(gz - zl) < step / 4) for zl in layers]
    for c in copper:
        fixed = fixed | c
    segment = (gs >= -0.5 - step / 4) & (gs < 0.5 - step / 4)
    value = np.stack([np.where(c & segment, 1.0, 0.0) for c in copper])
    phi = laplace3d(np.broadcast_to(fixed, value.shape), value)
    g = phi[:, 1, :, :] / step
    return RampBases(s, z, g, tuple(layers), slot)


def compile_ramp(bases: RampBases, length: float, target: float, bound: float, g0: float, h: float) -> dict[str, Any]:
    """Weights in [-bound, bound] on 1 mm segments of every layer that carry the null up a quintic climb from the
    lowest level's slot to the highest's over `length` (mm), holding the gradient near target x G0: bounded least
    squares on the null's position, its gradient and the weights' size (section 3.3)."""
    zc = (bases.layers[1] + bases.layers[0]) / 2
    zc = abs(zc)
    centres = np.arange(-9.0, length + 9.0 + 1e-9, 1.0)
    sg = np.arange(-4.0, length + 4.0 + 1e-9, 0.25)
    path = -zc + zc * (quintic(2 * sg / length - 1) + 1)
    dz = bases.z[1] - bases.z[0]
    cols_g, cols_d = [], []
    for gl in bases.g:
        for c in centres:
            cols_g.append(bilinear(bases.s, bases.z, gl, sg - c, path))
            cols_d.append(
                (bilinear(bases.s, bases.z, gl, sg - c, path + dz) - bilinear(bases.s, bases.z, gl, sg - c, path - dz))
                / (2 * dz)
            )
    ag, ad = np.array(cols_g).T, np.array(cols_d).T
    count = ag.shape[1]
    a = np.vstack([ag / (g0 * h) * 6, ad / g0, 0.03 * np.eye(count)])
    rhs = np.concatenate([np.zeros(len(sg)), -target * np.ones(len(sg)), np.zeros(count)])
    w = bounded_lstsq(a, rhs, bound)
    weights = w.reshape(len(bases.layers), len(centres))
    out: dict[str, Any] = {"centres": centres, "weights": weights, "path_s": sg, "path_z": path}
    for name, wts in (("compiled", weights), ("rounded", np.round(weights))):
        sf = np.arange(-2.0, length + 2.0 + 1e-9, 0.05)
        field = np.zeros((len(sf), len(bases.z)))
        for gl, wl in zip(bases.g, wts, strict=True):
            for wj, cj in zip(wl, centres, strict=True):
                if wj != 0:
                    field += wj * bilinear(
                        bases.s,
                        bases.z,
                        gl,
                        (sf - cj)[:, None] + 0 * bases.z[None, :],
                        np.broadcast_to(bases.z, (len(sf), len(bases.z))),
                    )
        tgt = np.interp(sf, sg, path)
        z0 = np.full(len(sf), np.nan)
        grad = np.full(len(sf), np.nan)
        for i in range(len(sf)):
            cands = [(zz, sl) for zz, sl in zero_crossings(bases.z, field[i]) if sl < 0 and abs(zz - tgt[i]) < 0.25]
            if cands:
                zz, sl = min(cands, key=lambda c: abs(c[0] - tgt[i]))
                z0[i], grad[i] = zz, abs(sl) / g0
        out[name] = {
            "s": sf,
            "null": z0,
            "gradient": grad,
            "followed": float(np.mean(np.isfinite(z0))),
            "path_error_um": float(np.nanmax(np.abs(z0 - tgt))) * 1e3,
            "gradient_min": float(np.nanmin(grad)),
            "gradient_max": float(np.nanmax(grad)),
        }
    return out


def cutout(h: float, t: float, cut: float, step: float) -> dict[str, Any]:
    """The ramp that must fail: a gapless stack of three boards (planes split down x = 0 at z = -zmax and +zmax, the
    middle board's faces at +-t/2), all four faces driven alternately, the middle board cut away for `cut` mm. The
    null of the lower channel jumps to the middle at each end of the cut-out (lengths mm)."""
    zmax = 2 * h + t / 2
    x = np.arange(0, 2.5 + step / 2, step)
    s = np.arange(-5.0, cut + 5.0 + step / 2, step)
    z = np.arange(-zmax, zmax + step / 2, step)
    gx, gs, gz = np.meshgrid(x, s, z, indexing="ij")
    fixed = (gx == 0) | (np.abs(gz) >= zmax - step / 4)
    value = np.where(np.abs(gz + zmax) < step / 4, 1.0, 0.0) + np.where(np.abs(gz - zmax) < step / 4, -1.0, 0.0)
    keep = (gs < 0) | (gs > cut)
    for zl, v in ((-t / 2, -1.0), (t / 2, 1.0)):
        face = (np.abs(gz - zl) < step / 4) & keep & (gx > 0)
        fixed = fixed | face
        value = np.where(face, v, value)
    value = np.where(gx == 0, 0.0, value)
    phi = laplace3d(fixed, value)
    g = phi[1] / step
    g0 = math.pi * 2 / (4 * h * h)
    nulls = []
    for i in range(len(s)):
        inside_board = not (0 <= s[i] <= cut)
        for zz, sl in zero_crossings(z, g[i]):
            if not (inside_board and abs(zz) < t / 2):  # the middle board is solid where it is not cut away
                nulls.append((float(s[i]), zz, abs(sl) / g0))
    mid = [n for n in nulls if cut * 0.3 < n[0] < cut * 0.7]
    outside = [n for n in nulls if n[0] < -3.0]
    return {
        "s": s,
        "z": z,
        "nulls": np.array(nulls),
        "inside_gradient": float(np.mean([n[2] for n in mid])),
        "predicted_inside": (h / (2 * h + t / 2)) ** 2,
        "outside_gradient": float(np.mean([n[2] for n in outside])),
        "levels": (-(h + t / 2), h + t / 2),
        "h": h,
        "t": t,
        "cut": cut,
    }


def slot_section(
    board: float, slot: float, pull: float, er: float, shield: float, other: float | None, fine: float
) -> dict[str, Any]:
    """The slot guide's cross-section (mm): one double-sided board `board` thick with a routed slot `slot` wide, the
    four slot-edge strips' copper set back `pull` from the wall, driven as a quadrupole (halves at +-1, so V = 2),
    laminate of permittivity er; grounded shields `shield` beyond the faces, or another board's copper `other` mm off
    (grounded). Solved on x >= 0 (odd in x). Returns the null's gradient and the escape saddle, in G0 and D0 of
    gapless facing boards with h = slot / 2 at the same volts."""
    h = slot / 2
    zmax = board / 2 + (other if other is not None else shield)
    edge = slot / 2 + pull
    x = nodes(0.0, 6.0, [slot / 2, edge], fine, 0.15, 0.25)
    z = nodes(-zmax, zmax, [-board / 2, board / 2, 0.0], fine, 0.15, 0.25)
    gx, gz = np.meshgrid(x, z, indexing="ij")
    fixed = (gx == 0) | (np.abs(gz) >= zmax - 1e-9)
    value = np.zeros_like(gx)
    for zz, v in ((board / 2, 1.0), (-board / 2, -1.0)):
        cu = (gx >= edge - 1e-9) & (np.abs(gz - zz) < 1e-9)
        fixed |= cu
        value[cu] = v
    xc, zcell = (x[1:] + x[:-1]) / 2, (z[1:] + z[:-1]) / 2
    cx, cz = np.meshgrid(xc, zcell, indexing="ij")
    eps = np.where((cx > slot / 2) & (np.abs(cz) < board / 2), er, 1.0)
    phi = laplace2d(x, z, eps, fixed, value)
    ex, ez = np.gradient(-phi, x, z)
    emag = np.hypot(ex, ez)
    g = phi[1] / x[1]
    j0 = int(np.argmin(np.abs(z)))
    grad = abs((g[j0 + 1] - g[j0 - 1]) / (z[j0 + 1] - z[j0 - 1]))
    g0 = math.pi * 2 / (4 * h * h)
    eb0 = 2 / (2 * h)

    def is_exit(i: int, j: int) -> bool:
        """The slot's wall, copper, open space beyond the board's faces: anywhere but the symmetry plane x = 0."""
        wall = x[i] >= slot / 2 - 1e-9 and abs(z[j]) <= board / 2
        return bool(wall or abs(z[j]) > board / 2 + 0.6 or (fixed[i, j] and x[i] > 0))

    level, exit_at = escape_flood(x, z, emag, (1, j0), is_exit)
    return {
        "x": x,
        "z": z,
        "field": emag / eb0,
        "gradient_G0": grad / g0,
        "depth_D0": (level / eb0) ** 2,
        "saddle": exit_at,
        "board": board,
        "slot": slot,
        "pull": pull,
        "er": er,
    }


def escape_flood(
    x: Arr, z: Arr, emag: Arr, start: tuple[int, int], is_exit: Callable[[int, int], bool]
) -> tuple[float, tuple[float, float]]:
    """The lowest |E| a particle must cross to leave a well: a minimax flood from the null (Dijkstra on the largest
    field met), stopping at the first exit. Returns that field and where the path met it."""
    nx, nz = emag.shape
    seen = np.zeros((nx, nz), dtype=bool)
    heap = [(float(emag[start]), start[0], start[1])]
    level, where = 0.0, (float(x[start[0]]), float(z[start[1]]))
    while heap:
        e, i, j = heapq.heappop(heap)
        if seen[i, j]:
            continue
        seen[i, j] = True
        if e > level:
            level, where = e, (float(x[i]), float(z[j]))
        if is_exit(i, j):
            break
        for a, b in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)):
            if 0 <= a < nx and 0 <= b < nz and not seen[a, b]:
                heapq.heappush(heap, (float(emag[a, b]), a, b))
    return level, where


def compute_ramp(quick: bool) -> dict[str, Any]:
    say("\n== ramp: two slotted boards' four layers, the cut-out that must fail, the slot guide's section ==")
    step = 0.1 if quick else 0.05
    h, v = 0.5, 2.0
    g0 = math.pi * v / (4 * h * h)
    out: dict[str, Any] = {}
    board, gap = 0.8, 1.0
    layers = (-(gap / 2 + board), -gap / 2, gap / 2, gap / 2 + board)
    began = time.time()
    bases = ramp_bases(layers, 1.0, step)
    say(
        f"  bases of four layers (0.8 mm boards 1.0 mm apart, 1.0 mm slot) on a {step} mm grid "
        f"[{time.time() - began:.0f} s]"
    )
    out["stack"] = {"layers": layers, "board": board, "gap": gap, "slot": 1.0}
    for name, target, bound in (("rail", 0.45, 1.0), ("twice the rail", 0.9, 2.0)):
        comp = compile_ramp(bases, 16.0, target, bound, g0, h)
        out[name] = comp
        c = comp["compiled"]
        r = comp["rounded"]
        say(
            f"  16 mm ramp at {name} (|w| <= {bound:g}), target {target} G0: null followed {c['followed'] * 100:.0f} "
            "%, path "
            f"error {c['path_error_um']:.0f} um, gradient {c['gradient_min']:.2f} to {c['gradient_max']:.2f} G0; "
            "rounded "
            f"to -1/0/+1: followed {r['followed'] * 100:.0f} %, gradient {r['gradient_min']:.2f} to "
            f"{r['gradient_max']:.2f} G0"
        )
    began = time.time()
    cut = cutout(0.5, 0.4, 8.0, step)
    out["cutout"] = cut
    say(
        f"  cut-out of the middle board (h 0.5, t 0.4 mm) [{time.time() - began:.0f} s]: gradient "
        f"{cut['outside_gradient']:.2f} G0 "
        f"outside, {cut['inside_gradient']:.3f} G0 inside against (h/(2h + t/2))^2 = {cut['predicted_inside']:.3f}"
    )
    fine = 0.02 if quick else 0.01
    sections = {
        "flush, vacuum, shields 1.5 mm": (0.0, 1.0, 1.5, None),
        "flush, FR-4, shields 1.5 mm": (0.0, 4.4, 1.5, None),
        "set back 0.2 mm, FR-4, shields 1.5 mm": (0.2, 4.4, 1.5, None),
        "set back 0.2 mm, FR-4, the other board 1.0 mm off": (0.2, 4.4, 0.0, 1.0),
    }
    out["sections"] = {}
    for name, (pull, er, shield, other) in sections.items():
        sec = slot_section(0.8, 1.0, pull, er, shield, other, fine)
        out["sections"][name] = sec
        say(
            f"  slot guide, {name}: {sec['depth_D0']:.3f} D0, {sec['gradient_G0']:.3f} G0, saddle at x "
            f"{sec['saddle'][0]:.2f}, z {sec['saddle'][1]:.2f} mm"
        )
    return out


def figure_ramp(res: dict[str, Any], out: Path) -> None:
    """Figure 3: along a ramp of the two-board stack, the four layers coloured by compiled volts with the traced null;
    the failing cut-out; the slot guide's section, flush and set back, with its escape saddle."""
    plt = _pyplot()
    fig = plt.figure(figsize=(15.5, 10.0))
    grid = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.36, height_ratios=[1.0, 1.0])
    rail = res["rail"]
    layers = res["stack"]["layers"]
    ax = fig.add_subplot(grid[0, :2])
    cmap = _diverging()
    centres, weights = rail["centres"], rail["weights"]
    for zl, wl in zip(layers, weights, strict=True):
        for c, w in zip(centres, wl, strict=True):
            ax.add_patch(_rect_signed(c - 0.5, zl - 0.06, 1.0, 0.12, cmap((w + 1) / 2)))
    board = res["stack"]["board"]
    for lo in (layers[0], layers[2]):
        ax.add_patch(_rect(centres[0] - 0.5, lo, centres[-1] - centres[0] + 1, board, LAMINATE))
    c = rail["compiled"]
    ax.plot(rail["path_s"], rail["path_z"], color=MUTED, lw=1.0, ls="--", label="target: a quintic climb over 16 mm")
    ax.plot(c["s"], c["null"], color=INK, lw=1.8, label=f"traced null, rail volts: error ≤ {c['path_error_um']:.0f} µm")
    tw = res["twice the rail"]["compiled"]
    ax.plot(
        tw["s"],
        tw["null"],
        color=GREEN,
        lw=1.2,
        label=f"traced null, twice the rail: error ≤ {tw['path_error_um']:.0f} µm",
    )
    ax.set_xlim(centres[0] - 0.5, centres[-1] + 0.5)
    ax.set_ylim(layers[0] - 0.4, layers[-1] + 0.75)
    ax.set_xlabel("along the path s (mm)")
    ax.set_ylabel("height z (mm)")
    ax.legend(loc="upper left", fontsize=8)
    ax.set_title(
        "the four layers' 1 mm segments coloured by compiled volts (−1 blue to +1 orange, rail volts), boards shaded"
    )
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-1, 1))
    cax = ax.inset_axes((0.7, 0.935, 0.28, 0.03))
    fig.colorbar(sm, cax=cax, orientation="horizontal")
    cax.tick_params(labelsize=7, pad=1)
    ax.text(0.69, 0.95, "compiled volts over the rail ", transform=ax.transAxes, ha="right", va="center", fontsize=7.8)
    ax = fig.add_subplot(grid[0, 2])
    for name, colour in (("rail", INK), ("twice the rail", GREEN)):
        for kind, ls in (("compiled", "-"), ("rounded", ":")):
            d = res[name][kind]
            ax.plot(
                d["s"],
                d["gradient"],
                color=colour,
                ls=ls,
                lw=1.4,
                label=f"{name}, {kind}: least {d['gradient_min']:.2f} G₀",
            )
    ax.set_xlabel("along the path s (mm)")
    ax.set_ylabel("gradient at the null over G₀")
    ax.set_ylim(0, 1.1)
    ax.legend(fontsize=7.6, loc="lower left")
    ax.set_title("the gradient along the ramp: a floor where the null\ncrosses an inner copper plane")
    cut = res["cutout"]
    ax = fig.add_subplot(grid[1, 0])
    nulls = cut["nulls"]
    sc = ax.scatter(nulls[:, 0], nulls[:, 1], c=nulls[:, 2], s=4, cmap=_ramp(), vmin=0, vmax=1.0)
    cax = ax.inset_axes((0.62, 0.6, 0.34, 0.03))
    fig.colorbar(sc, cax=cax, orientation="horizontal")
    cax.tick_params(labelsize=7, pad=1)
    ax.text(0.605, 0.615, "gradient over G₀ ", transform=ax.transAxes, ha="right", va="center", fontsize=7.8)
    h, t = cut["h"], cut["t"]
    zmax = 2 * h + t / 2
    for zz in (-zmax, zmax):
        ax.axhline(zz, color=COPPER, lw=2.2)
    for zz in (-t / 2, t / 2):
        ax.plot([cut["s"][0], 0], [zz, zz], color=COPPER, lw=2.2)
        ax.plot([cut["cut"], cut["s"][-1]], [zz, zz], color=COPPER, lw=2.2)
    ax.set_xlabel("along the path s (mm)")
    ax.set_ylabel("height z (mm)")
    ax.set_title(
        f"must fail: the middle board cut away for {cut['cut']:g} mm. The null jumps half\nthe rise at each end and "
        "sits between at "
        f"{cut['inside_gradient']:.3f} G₀ ((h/(2h + t/2))² = {cut['predicted_inside']:.3f})"
    )
    for k, name in enumerate(("flush, FR-4, shields 1.5 mm", "set back 0.2 mm, FR-4, shields 1.5 mm")):
        sec = res["sections"][name]
        ax = fig.add_subplot(grid[1, 1 + k])
        x, z = sec["x"], sec["z"]
        xs = np.concatenate([-x[::-1], x[1:]])
        field = np.concatenate([sec["field"][::-1], sec["field"][1:]])
        mesh = ax.pcolormesh(xs, z, field.T, cmap=_ramp(), vmin=0, vmax=1.6, shading="auto")
        b, sl, pull = sec["board"], sec["slot"], sec["pull"]
        for sign in (-1, 1):
            ax.add_patch(_rect(sign * sl / 2 if sign > 0 else -2.5, -b / 2, 2.5 - sl / 2, b, LAMINATE))
            for zz in (-b / 2, b / 2):
                ax.plot([sign * (sl / 2 + pull), sign * 2.5], [zz, zz], color=COPPER, lw=2.6)
        ax.plot([0], [0], "+", color=ORANGE, ms=12, mew=2)
        sx, sz = sec["saddle"]
        ax.plot([sx, -sx], [sz, sz], "x", color=INK, ms=9, mew=2)
        ax.set_xlim(-2.0, 2.0)
        ax.set_ylim(-1.4, 1.4)
        ax.set_aspect("equal")
        ax.set_xlabel("across the slot x (mm)")
        ax.set_ylabel("z (mm)")
        ax.set_title(
            f"slot guide, {name}:\n{sec['depth_D0']:.2f} D₀ deep, {sec['gradient_G0']:.2f} G₀; saddle (×) at z = "
            f"{sz:+.2f} mm"
        )
        if k == 1:
            fig.colorbar(mesh, ax=ax, fraction=0.04, pad=0.02, label="|E| over V/2h")
    secs = res["sections"]
    _suptitle(
        fig,
        "Stacked boards: bounded least squares over 1 mm segments of two slotted boards' four layers carries the null "
        "up a "
        f"16 mm climb as a sign pattern, the gradient falling to {rail['compiled']['gradient_min']:.2f} G₀ where it "
        "crosses an "
        f"inner copper plane at rail volts; a cut-out leaves {cut['inside_gradient']:.2f} G₀; the slot guide holds "
        f"{secs['flush, vacuum, shields 1.5 mm']['depth_D0']:.2f} to "
        f"{secs['flush, FR-4, shields 1.5 mm']['depth_D0']:.2f} D₀ flush and "
        f"{secs['set back 0.2 mm, FR-4, shields 1.5 mm']['depth_D0']:.2f} D₀ with its copper set back 0.2 mm",
    )
    _save(fig, out, "ramp")


def _rect_signed(x: float, y: float, w: float, h: float, colour: Any) -> Any:
    from matplotlib.patches import Rectangle

    return Rectangle((x, y), w, h, facecolor=colour, edgecolor="none")


# ---- joints (figure 5) --------------------------------------------------------------------------------------------
def planes_axis(top: Arr, bottom: Arr, step: float, h: float = 1.0) -> dict[str, Arr]:
    """Two parallel planes at y = +-h carrying potentials top(x, s) and bottom(x, s) on a periodic grid (x first, the
    split at the grid's middle column): the field on the axis x = y = 0 along s, by the Fourier kernels of the
    pattern's odd part, k / sinh(kh) for E_y on the midplane, and its even part, 1 / cosh(kh) for the potential
    there (section 4.5)."""
    nx, ns = top.shape
    kx = 2 * math.pi * np.fft.fftfreq(nx, step)
    ks = 2 * math.pi * np.fft.fftfreq(ns, step)
    kxx, kss = np.meshgrid(kx, ks, indexing="ij")
    k = np.hypot(kxx, kss)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        odd_kernel = np.where(k > 0, k / np.sinh(np.minimum(k * h, 700.0)), 1 / h)
    even_kernel = 1 / np.cosh(np.minimum(k * h, 700.0))
    odd = np.fft.fft2((top - bottom) / 2) * odd_kernel
    even = np.fft.fft2((top + bottom) / 2) * even_kernel
    mid = nx // 2
    ey = -np.fft.ifft2(odd).real[mid]
    grad = -np.fft.ifft2(odd * 1j * kxx).real[mid]
    phi = np.fft.ifft2(even).real[mid]
    ex = -np.fft.ifft2(even * 1j * kxx).real[mid]
    es = -np.fft.ifft2(even * 1j * kss).real[mid]
    return {"Ey": ey, "gradient": grad, "phi": phi, "Ex": ex, "Es": es}


def segment_patterns(nx: int, ns: int, step: float, start: float, stop: float) -> dict[str, tuple[Arr, Arr]]:
    """The four drives of one segment from s = start to stop, per volt (patterns top left, top right, bottom left,
    bottom right): Q (+ - - +) confines, C (+ + + +) is common mode, Y (+ + - -) steers along x, X (+ - + -) along y
    (section 4.6). Edges soft over one cell."""
    x = (np.arange(nx) - nx // 2) * step
    s = (np.arange(ns) - ns // 2) * step
    left = 0.5 * np.clip(-x / step, -1, 1) + 0.5  # 1 for x < 0, 0 for x > 0
    span = np.clip((s - start) / step + 0.5, 0, 1) * np.clip((stop - s) / step + 0.5, 0, 1)
    lx, sp = np.meshgrid(left, span, indexing="ij")
    half = 0.5 * sp
    signed = (2 * lx - 1) * half  # +1/2 on the left, -1/2 on the right
    return {
        "Q": (signed, -signed),
        "C": (2 * half, 2 * half),
        "Y": (half, -half),
        "X": (signed, signed),
    }


def offset_joint(cases: Sequence[tuple[float, float]], gap: float, step: float, reach: float = 4.0) -> dict[str, Any]:
    """Facing boards (h = 1, halves at +-1/2, 6 wide) in two segments across a bare gap, the second segment offset by
    dx sideways and dy toward the top board, one case per (dx, dy): solved together in three dimensions inside
    grounded planes at y = +-2.2; the null traced along s with the axial field on it."""
    x = np.arange(-4.0, 4.0 + step / 2, step)
    y = np.arange(-2.2, 2.2 + step / 2, step)
    s = np.arange(-reach, reach + step / 2, step)
    gx, gy, gs = np.meshgrid(x, y, s, indexing="ij")
    fixed_all, value_all = [], []
    for dx, dy in cases:
        fixed = np.abs(gy) >= 2.2 - step / 4
        value = np.zeros_like(gx)
        for lo, hi, ox, oy in ((-np.inf, -gap / 2, 0.0, 0.0), (gap / 2, np.inf, dx, dy)):
            inside = (gs >= lo - step / 4) & (gs <= hi + step / 4) & (np.abs(gx - ox) <= 3.0 + step / 4)
            for plane, sign in ((1.0, 1.0), (-1.0, -1.0)):
                board = inside & (np.abs(gy - (plane + oy)) < step / 4)
                v = sign * 0.5 * np.sign(ox - gx)  # top: +1/2 left of the split, -1/2 right; bottom the opposite
                fixed = fixed | board
                value = np.where(board, v, value)
        fixed_all.append(fixed)
        value_all.append(value)
    phi = laplace3d(np.array(fixed_all), np.array(value_all), tol=1e-9)
    out: dict[str, Any] = {"s": s, "x": x, "y": y, "cases": []}
    for k, (dx, dy) in enumerate(cases):
        p = phi[k]
        ex = -np.gradient(p, step, axis=0)
        ey = -np.gradient(p, step, axis=1)
        es = -np.gradient(p, step, axis=2)
        null = np.full((len(s), 2), np.nan)
        axial = np.full(len(s), np.nan)
        win = (np.abs(x)[:, None] < 0.6) & (np.abs(y)[None, :] < 0.6)
        for i in range(len(s)):
            e2 = ex[:, :, i] ** 2 + ey[:, :, i] ** 2
            a, b = np.unravel_index(int(np.argmin(np.where(win, e2, np.inf))), e2.shape)
            px, py = float(x[a]), float(y[b])
            for _ in range(12):  # Newton on the bilinear interpolant of (E_x, E_y)
                vx, gx = cell_bilinear(x, y, ex[:, :, i], px, py)
                vy, gy = cell_bilinear(x, y, ey[:, :, i], px, py)
                jac = np.array([gx, gy])
                if abs(float(np.linalg.det(jac))) < 1e-300:
                    break
                d = np.linalg.solve(jac, -np.array([vx, vy]))
                px, py = px + float(d[0]), py + float(d[1])
                if abs(d[0]) + abs(d[1]) < 1e-12:
                    break
            null[i] = (px, py)
            axial[i] = cell_bilinear(x, y, es[:, :, i], px, py)[0]
        out["cases"].append(
            {"dx": dx, "dy": dy, "null": null, "axial": axial, "peak_axial": float(np.nanmax(np.abs(axial)))}
        )
    return out


def cell_bilinear(x: Arr, y: Arr, values: Arr, px: float, py: float) -> tuple[float, tuple[float, float]]:
    """values on a grid's nodes, bilinearly interpolated at (px, py) within its cell, with the interpolant's
    gradient."""
    i = int(np.clip(np.searchsorted(x, px) - 1, 0, len(x) - 2))
    j = int(np.clip(np.searchsorted(y, py) - 1, 0, len(y) - 2))
    hx, hy = float(x[i + 1] - x[i]), float(y[j + 1] - y[j])
    u, v = (px - float(x[i])) / hx, (py - float(y[j])) / hy
    f00, f10, f01, f11 = (
        float(values[i, j]),
        float(values[i + 1, j]),
        float(values[i, j + 1]),
        float(values[i + 1, j + 1]),
    )
    val = f00 * (1 - u) * (1 - v) + f10 * u * (1 - v) + f01 * (1 - u) * v + f11 * u * v
    du = ((f10 - f00) * (1 - v) + (f11 - f01) * v) / hx
    dv = ((f01 - f00) * (1 - u) + (f11 - f10) * u) / hy
    return val, (du, dv)


def width_floor(widths: Sequence[float], step: float, box: float = 8.0, reach: float = 14.0) -> dict[str, Any]:
    """A common-mode step on facing boards of finite width W (h = 1): both boards at 1 for s > 0 and 0 for s < 0,
    inside grounded walls `box` beyond the boards' edges and above them. Solved on the quarter x, y >= 0 (the
    drive is even in both); the potential on the axis for s < 0 is the even disturbance's reach."""
    out: dict[str, Any] = {"widths": list(widths), "curves": []}
    for width in widths:
        x = np.arange(0, width / 2 + box + step / 2, step)
        y = np.arange(0, 1 + box + step / 2, step)
        s = np.arange(-reach, 4.0 + step / 2, step)
        gx, gy, gs = np.meshgrid(x, y, s, indexing="ij")
        board = (np.abs(gy - 1) < step / 4) & (gx <= width / 2 + step / 4)
        fixed = board | (gx >= x[-1] - step / 4) | (gy >= y[-1] - step / 4)
        value = np.where(board & (gs > 0), 1.0, np.where(board & (np.abs(gs) < step / 4), 0.5, 0.0))
        phi = laplace3d(fixed, value, tol=1e-11)
        out["curves"].append({"width": width, "s": s, "phi": phi[0, 0, :]})
    return out


def banded_operator(segments: int, length: float, step: float) -> dict[str, Any]:
    """The operator of a segmented guide (h = 1): each of `segments` segments `length` long round a lap carries the
    four drives Q, C, Y, X; G maps them to the gradient, E_y, E_x and E_s on the axis, sampled every h/4 (the
    curvature of a knot enters at order h/R)."""
    lap = segments * length
    ns = round(lap / step)
    nx = 256
    samples = np.arange(0, ns, max(round(0.25 / step), 1))
    cols = []
    names = []
    for j in range(segments):
        start = -lap / 2 + j * length
        pats = segment_patterns(nx, ns, step, start, start + length)
        for drive in ("Q", "C", "Y", "X"):
            top, bottom = pats[drive]
            a = planes_axis(top, bottom, step)
            cols.append(np.concatenate([a[c][samples] for c in ("gradient", "Ey", "Ex", "Es")]))
            names.append(f"{drive}{j}")
    g = np.array(cols).T
    return {"G": g, "names": names, "rows_per_component": len(samples), "segments": segments, "length": length}


def compute_joints(quick: bool) -> dict[str, Any]:
    say("\n== joints: the offset joint in 3-D, decay by parity and the width floor, basis fields and banded G ==")
    out: dict[str, Any] = {}
    step = 1 / 16 if quick else 1 / 32
    nx, ns = 512, 2048 if not quick else 1024
    # the four basis fields of one segment 6h long
    basis = {}
    pats = segment_patterns(nx, ns, step, -3.0, 3.0)
    s = (np.arange(ns) - ns // 2) * step
    for drive, (top, bottom) in pats.items():
        basis[drive] = planes_axis(top, bottom, step)
    out["basis"] = {"s": s, "fields": basis, "length": 6.0}
    centre = ns // 2
    say(
        "  one segment 6h long, per volt at its centre: Q gradient "
        f"{basis['Q']['gradient'][centre] / (math.pi / 4):.4f} of pi/4h^2, "
        f"Y field {abs(basis['Y']['Ey'][centre]):.4f} (1/2h = 0.5), X field {abs(basis['X']['Ex'][centre]):.4f}, "
        f"C axial field at its end {abs(basis['C']['Es'][int(np.argmin(np.abs(s - 3.0)))]):.4f} (1/2h = 0.5)"
    )
    # decay by parity: a Q segment's end (odd) and a C step (even) against exp(-pi s/h) and exp(-pi s/2h)
    long_q = segment_patterns(nx, ns, step, 0.0, ns * step)["Q"]
    qa = planes_axis(*long_q, step)
    c_step = segment_patterns(nx, ns, step, 0.0, ns * step)["C"]
    ca = planes_axis(*c_step, step)
    behind = (s < -1.5) & (s > -5)
    odd_rate = -np.polyfit(-s[behind], np.log(np.abs(qa["gradient"][behind]) / (math.pi / 4)), 1)[0]
    even_rate = -np.polyfit(-s[behind], np.log(np.abs(ca["phi"][behind])), 1)[0]
    closed = 0.5 / np.cosh(math.pi * s / 2)
    out["decay"] = {
        "s": s,
        "odd": np.abs(qa["gradient"]) / (math.pi / 4),
        "even": np.abs(ca["phi"]),
        "even_closed": 2 / math.pi * np.arctan(np.exp(math.pi * s / 2)),
        "odd_rate": float(odd_rate),
        "even_rate": float(even_rate),
        "step_error": float(np.abs(np.abs(ca["Es"]) - closed)[np.abs(s) < 4].max()),
    }
    say(
        f"  decay behind a segment's end: odd (gradient) e-folds at {odd_rate:.3f}/h (pi = {math.pi:.3f}), even "
        "(common-mode "
        f"step's potential) at {even_rate:.3f}/h (pi/2 = {math.pi / 2:.3f}); the step's axial field is (1/2h) sech(pi "
        "s/2h) "
        f"to {out['decay']['step_error']:.1e}"
    )
    began = time.time()
    floor = width_floor((4.0, 6.0), 0.25 if quick else 0.125, box=6.0 if quick else 8.0)
    out["floor"] = floor
    for c in floor["curves"]:
        wide = 2 / math.pi * np.arctan(np.exp(math.pi * c["s"] / 2))
        behind_c = c["s"] < 0
        departed = np.flatnonzero(behind_c & (c["phi"] > 2 * wide))
        knee = int(departed[-1]) if len(departed) else 0
        c["knee_s"] = float(-c["s"][knee])
        c["knee"] = float(c["phi"][knee])
        c["predicted"] = 0.1 * math.exp(-math.pi * c["width"] / 4)
        say(
            f"  boards {c['width']:g}h wide in grounded walls: the common-mode step's potential leaves the wide "
            "boards' "
            f"exp(-pi s/2h) (twice it) {c['knee_s']:.1f}h behind the step at {c['knee']:.1e} of the step; "
            f"0.1 exp(-pi W/4h) = {c['predicted']:.1e}"
        )
    say(f"  [width floor {time.time() - began:.0f} s]")
    began = time.time()
    joint_step = 0.2 if quick else 0.1
    offs = 0.2
    joint = offset_joint([(offs, offs), (offs, 0.0), (0.0, offs), (0.0, 0.0)], 0.2, joint_step)
    out["joint"] = joint
    control = 0.32 * offs * offs
    for c in joint["cases"]:
        say(
            f"  joint offset dx {c['dx']:g}h, dy {c['dy']:g}h across a 0.2h gap: largest axial field on the null "
            f"{c['peak_axial']:.2e} V/h"
        )
    out["control"] = {"predicted": control, "measured": joint["cases"][0]["peak_axial"]}
    say(
        f"  the control reads {joint['cases'][0]['peak_axial']:.4f} against 0.32 V dx dy / h^3 = {control:.4f} "
        f"[{time.time() - began:.0f} s]"
    )
    out["banded"] = banded_operator(12, 4.0, 1 / 8 if quick else 1 / 16)
    g = out["banded"]["G"]
    say(
        f"  the banded operator of 12 segments 4h long: {g.shape[0]} rows, {g.shape[1]} drives; "
        f"{np.mean(np.abs(g) > 1e-3 * np.abs(g).max()):.2f} of its entries above 1e-3 of the largest"
    )
    return out


def figure_joints(res: dict[str, Any], out: Path) -> None:
    """Figure 5: the joint offset both ways in three dimensions with its null coloured by axial field; the decay rates
    and the width floor; the four basis fields of one segment with the banded G of a segmented guide."""
    plt = _pyplot()
    fig = plt.figure(figsize=(15.5, 10.4))
    grid = fig.add_gridspec(2, 3, hspace=0.36, wspace=0.34)
    joint = res["joint"]
    both = joint["cases"][0]
    ax = fig.add_subplot(grid[0, 0], projection="3d", computed_zorder=False)
    for lo, hi, ox, oy in ((-4.0, -0.1, 0.0, 0.0), (0.1, 4.0, both["dx"], both["dy"])):
        for plane, colours in ((1.0, (ORANGE, BLUE)), (-1.0, (BLUE, ORANGE))):
            yy = plane + oy
            for (x0, x1), colour in zip(((-3.0 + ox, ox), (ox, 3.0 + ox)), colours, strict=True):
                quad = [[(x0, lo, yy), (x1, lo, yy), (x1, hi, yy), (x0, hi, yy)]]
                ax.add_collection3d(_poly3d(quad, facecolor=colour, alpha=0.18, edgecolor=colour, lw=0.4))
    s = joint["s"]
    null = both["null"]
    scale = max(float(np.nanmax(np.abs(both["axial"]))), 1e-12)
    sc = ax.scatter(
        null[:, 0], s, null[:, 1], c=both["axial"], cmap=_diverging(), vmin=-scale, vmax=scale, s=10, depthshade=False
    )
    fig.colorbar(
        sc, ax=ax, orientation="horizontal", fraction=0.04, pad=0.06, shrink=0.7, label="axial field on the null (V/h)"
    )
    ax.set_xlabel("x / h", fontsize=8)
    ax.set_ylabel("s / h", fontsize=8)
    ax.set_zlabel("y / h", fontsize=8)
    ax.set_xlim(-3.2, 3.4)
    ax.set_zlim(-1.4, 1.4)
    ax.view_init(elev=18, azim=-62)
    ax.set_title(
        f"the joint offset {both['dx']:g}h sideways and {both['dy']:g}h up across a 0.2h gap:\nthe null jogs both ways "
        "and carries an axial field",
        fontsize=9,
    )
    ax = fig.add_subplot(grid[0, 1])
    for k, c in enumerate(joint["cases"]):
        label = f"dx {c['dx']:g}h, dy {c['dy']:g}h: largest {sci(c['peak_axial'])} V/h"
        ax.plot(s, c["axial"], color=SERIES[k], lw=1.6 if k == 0 else 1.2, label=label)
    ctrl = res["control"]
    ax.axhline(ctrl["predicted"], color=MUTED, ls=":", lw=1.0)
    ax.axhline(-ctrl["predicted"], color=MUTED, ls=":", lw=1.0)
    ax.text(s[-1], ctrl["predicted"] * 0.9, "0.32 V·dx·dy/h³ ", fontsize=8, color=INK2, ha="right", va="top")
    ax.set_xlabel("along the guide s / h")
    ax.set_ylabel("axial field on the null (V/h)")
    ax.legend(fontsize=7.6, loc="lower left")
    ax.set_title(
        "one antisymmetry kept leaves no axial field (SEG-1);\nan offset both ways is the control that must read 0.32 "
        "V dx dy/h³"
    )
    ax = fig.add_subplot(grid[0, 2])
    d = res["decay"]
    sd = d["s"]
    behind = sd < 0
    ax.semilogy(
        -sd[behind],
        d["odd"][behind],
        color=BLUE,
        label=f"odd: a Q segment's gradient behind its end, {d['odd_rate']:.2f}/h",
    )
    ax.semilogy(
        -sd[behind],
        d["even"][behind],
        color=ORANGE,
        label=f"even: a common-mode step's potential, {d['even_rate']:.2f}/h",
    )
    u = np.linspace(0, 12, 100)
    ax.semilogy(u, 0.5 * np.exp(-math.pi * u), color=BLUE, ls=":", lw=1.0, label="exp(−πs/h)")
    ax.semilogy(u, 2 / math.pi * np.exp(-math.pi * u / 2), color=ORANGE, ls=":", lw=1.0, label="exp(−πs/2h)")
    for k, c in enumerate(res["floor"]["curves"]):
        sw = c["s"]
        m = sw < 0
        colour = SERIES[2 + k]
        ax.semilogy(
            -sw[m],
            np.abs(c["phi"][m]),
            color=colour,
            lw=1.3,
            label=f"even, boards {c['width']:g}h wide in grounded walls",
        )
        ax.axhline(c["predicted"], color=colour, ls="--", lw=0.9)
        ax.plot([c["knee_s"]], [c["knee"]], "o", color=colour, ms=5)
        ax.text(
            11.8, c["predicted"] * 1.25, f"0.1·exp(−πW/4h), W = {c['width']:g}h", fontsize=7.4, color=INK2, ha="right"
        )
    ax.set_xlim(0, 12)
    ax.set_ylim(1e-9, 2)
    ax.set_xlabel("distance behind the change, s / h")
    ax.set_ylabel("disturbance on the axis, of the drive")
    ax.legend(fontsize=7.2, loc="lower left")
    ax.set_title("decay by parity (SEG-3) and the width floor")
    ax = fig.add_subplot(grid[1, 0])
    b = res["basis"]
    sb = b["s"]
    view = np.abs(sb) < 7
    fields = b["fields"]
    ax.plot(
        sb[view], fields["Q"]["gradient"][view] / (math.pi / 4), color=BLUE, label="Q (+ − − +): gradient over π/4h²"
    )
    ax.plot(sb[view], fields["C"]["Es"][view] * 2, color=ORANGE, label="C (+ + + +): axial field over 1/2h")
    ax.plot(sb[view], -fields["Y"]["Ey"][view] * 2, color=GREEN, label="Y (+ + − −): E_y over 1/2h")
    ax.plot(sb[view], fields["X"]["Ex"][view] * 2, color=PURPLE, label="X (+ − + −): E_x over 1/2h")
    ax.axvspan(-b["length"] / 2, b["length"] / 2, color=PALE, alpha=0.7, lw=0)
    ax.set_xlabel("along the guide s / h")
    ax.set_ylabel("on the axis, per volt")
    ax.legend(fontsize=7.6, loc="lower right")
    ax.set_title(f"one segment {b['length']:g}h long: its four drives on the axis")
    ax = fig.add_subplot(grid[1, 1:])
    band = res["banded"]
    g = np.abs(band["G"])
    norm = g / g.max(axis=0, keepdims=True)
    img = ax.imshow(
        np.log10(np.maximum(norm, 1e-6)),
        aspect="auto",
        cmap=_ramp(reverse=False),
        vmin=-6,
        vmax=0,
        interpolation="nearest",
    )
    fig.colorbar(img, ax=ax, fraction=0.025, pad=0.01, label="log₁₀ |G|, each drive over its largest")
    rows = band["rows_per_component"]
    for k, name in enumerate(("gradient", "E_y", "E_x", "E_s")):
        ax.axhline(k * rows - 0.5, color=INK2, lw=0.6)
        ax.text(-1.2, k * rows + rows / 2, name, fontsize=8, ha="right", va="center")
    ax.set_xticks(range(0, len(band["names"]), 4))
    ax.set_xticklabels([band["names"][k][1:] for k in range(0, len(band["names"]), 4)])
    ax.set_xlabel("drives: Q, C, Y, X of segments 0 to 11 (four columns each)")
    ax.set_yticks([])
    ax.set_title(
        f"the banded operator G of {band['segments']} segments {band['length']:g}h long round a lap: e = Gw on the "
        "axis, sampled every h/4"
    )
    ax.grid(False)
    _suptitle(
        fig,
        "Fields at joints: a piece that keeps one mirror antisymmetry carries no axial field on its null; an offset "
        "both ways "
        f"puts {sci(res['control']['measured'])} V/h on it (0.32 V·dx·dy/h³ = {res['control']['predicted']:.4f}); odd "
        "disturbances "
        f"fade as exp(−πs/h), even ones as exp(−πs/2h) down to a floor set by the boards' width",
    )
    _save(fig, out, "joints")


# ---- drive (figure 6) ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Ring:
    """A closed guide as a ring of transmission line (section 5.3): length (m), index, Q, impedance z0 (ohm) and
    order m; joints as short sections (centre, length, impedance); feeds at positions, each a Norton current with a
    source admittance, a series capacitor into 50 ohm sized for critical coupling of the set (compensated: its
    susceptance tuned out, a pure conductance). V(s_i) = sum_k G_ik w_k over `samples` positions."""

    length: float
    index: float
    q: float
    z0: float
    m: int
    joints: tuple[tuple[float, float, float], ...]
    feeds: tuple[float, ...]
    samples: int = 240
    compensated: bool = False

    @property
    def shunt(self) -> float:
        """R = Z0 Q / (2 pi m) (RNG-2)."""
        return self.z0 * self.q / (2 * math.pi * self.m)

    @property
    def conductance(self) -> float:
        """Each feed's source conductance for critical coupling of the set."""
        return 1 / (len(self.feeds) * self.shunt)

    @property
    def positions(self) -> Arr:
        out: Arr = (np.arange(self.samples) + 0.5) * self.length / self.samples
        return out

    def resonance(self) -> float:
        """f_m = m c / (n L), pulled down by the feeds' coupling capacitors when they are not compensated."""
        f0 = self.m * pt.C0 / (self.index * self.length)
        if self.compensated:
            return f0
        x = math.sqrt(self.conductance / (50 * (1 - 50 * self.conductance)))
        return f0 * (1 - len(self.feeds) * x * self.z0 / (4 * math.pi * self.m))

    def operator(self, freqs: Arr) -> CArr:
        """G (F, samples, K): volts at the samples per ampere of each feed's Norton current, at each frequency."""
        nf, k = len(freqs), len(self.feeds)
        beta = 2 * math.pi * freqs * self.index / pt.C0
        gamma = beta / (2 * self.q) + 1j * beta
        if self.compensated:
            yf = np.full(nf, self.conductance, dtype=np.complex128)
        else:
            x = math.sqrt(self.conductance / (50 * (1 - 50 * self.conductance)))
            yf = 1 / (50 + 1 / (1j * x * freqs / freqs[nf // 2]))
        events: list[
            tuple[float, int, int, float]
        ] = []  # position, kind (0 joint start, 1 joint end, 2 feed, 3 sample)
        for c, ell, zj in self.joints:
            events.append(((c - ell / 2) % self.length, 0, 0, zj))
            events.append(((c + ell / 2) % self.length, 1, 0, 0.0))
        for j, p in enumerate(self.feeds):
            events.append((p % self.length, 2, j, 0.0))
        for i, p in enumerate(self.positions):
            events.append((float(p), 3, i, 0.0))
        events.sort(key=lambda e: e[0])
        t = np.tile(np.eye(2, dtype=np.complex128), (nf, 1, 1))
        b = np.zeros((nf, 2, k), dtype=np.complex128)
        rec_t = np.zeros((self.samples, nf, 2), dtype=np.complex128)
        rec_b = np.zeros((self.samples, nf, k), dtype=np.complex128)
        pos, z = 0.0, self.z0
        inside = sum(1 for c, ell, _ in self.joints if (c - ell / 2) % self.length > (c + ell / 2) % self.length)
        if inside:
            z = next(zj for c, ell, zj in self.joints if (c - ell / 2) % self.length > (c + ell / 2) % self.length)
        for p, kind, idx, zj in [*events, (self.length, 4, 0, 0.0)]:
            if p > pos:
                gl = gamma * (p - pos)
                ch, sh = np.cosh(gl), np.sinh(gl)
                mat = np.empty((nf, 2, 2), dtype=np.complex128)
                mat[:, 0, 0], mat[:, 0, 1], mat[:, 1, 0], mat[:, 1, 1] = ch, -z * sh, -sh / z, ch
                t, b = mat @ t, mat @ b
                pos = p
            if kind == 0:
                z = zj
            elif kind == 1:
                z = self.z0
            elif kind == 2:
                shunt = np.tile(np.eye(2, dtype=np.complex128), (nf, 1, 1))
                shunt[:, 1, 0] = -yf
                t, b = shunt @ t, shunt @ b
                b[:, 1, idx] += 1.0
            elif kind == 3:
                rec_t[idx], rec_b[idx] = t[:, 0, :], b[:, 0, :]
        x0 = np.linalg.solve(np.eye(2)[None] - t, b)  # the state at s = 0 per feed current
        g: CArr = np.einsum("ifj,fjk->fik", rec_t, x0) + np.transpose(rec_b, (1, 0, 2))
        return g


def ripple(v: CArr) -> float:
    """(max |V| - min |V|) / (max |V| + min |V|): zero for a travelling wave (RNG-1), one for a standing wave."""
    a = np.abs(v)
    return float((a.max() - a.min()) / (a.max() + a.min()))


MUS = (0.0, 1e-3, 1e-2, 3e-2, 1e-1, 3e-1)


def solve_weights(g: CArr, target: CArr, mu: float, iters: int = 8) -> CArr:
    """Regularised least squares to the travelling-wave target, then the target's phase relaxed to the field's (a
    uniform |V| with any phase), repeated."""
    u, s, vh = np.linalg.svd(g, full_matrices=False)
    filt = s / (s * s + mu * float(s[0]) ** 2)
    w: CArr = vh.conj().T @ (filt * (u.conj().T @ target))
    for _ in range(iters):
        w = vh.conj().T @ (filt * (u.conj().T @ np.exp(1j * np.angle(g @ w))))
    return w


def best_drive(
    ring: Ring, fixed: CArr | None = None, span: float = 0.03, count: int = 41, least: float = 0.8
) -> dict[str, Any]:
    """The least ripple over a frequency scan about the resonance and over the regularisations, among drives that
    deliver at least `least` of the available power to the ring (when none does, among those within 0.9 of the most
    efficient), or of the fixed weights alone."""
    freqs = ring.resonance() * (1 + np.linspace(-span, span, count))
    gs = ring.operator(freqs)
    target = np.exp(-2j * math.pi * ring.m * ring.positions / ring.length)
    cands: list[dict[str, Any]] = []
    for i in range(count):
        ws = [fixed] if fixed is not None else [solve_weights(gs[i], target, mu) for mu in MUS]
        for w in ws:
            v = gs[i] @ w
            rms = float(np.sqrt(np.mean(np.abs(v) ** 2)))
            eta = (rms**2 / (2 * ring.shunt)) / (float(np.sum(np.abs(w) ** 2)) / (8 * ring.conductance))
            cands.append({"ripple": ripple(v), "efficiency": eta, "f": float(freqs[i]), "w": w / rms, "v": v / rms})
    ok = [c for c in cands if c["efficiency"] >= least]
    if not ok:
        top = max(c["efficiency"] for c in cands)
        ok = [c for c in cands if c["efficiency"] >= 0.9 * top]
    return min(ok, key=lambda c: float(c["ripple"]))


def feed_layout(k: int, length: float, m: int, joints: int, at_joints: bool) -> tuple[float, ...]:
    """K feeds: two a quarter guided wave apart; K not dividing 2m equally spaced; K dividing 2m equally spaced with
    alternate feeds moved a quarter wave (so they reach -m); or at the joints."""
    lam = length / m
    if at_joints:
        return tuple((j + 0.5) * length / joints for j in range(k))
    if k == 1:
        return (0.03 * length,)
    if k == 2:
        return (0.03 * length, 0.03 * length + lam / 4)
    if (2 * m) % k:
        return tuple(0.03 * length + j * length / k for j in range(k))
    return tuple(0.03 * length + j * length / k + (lam / 4 if j % 2 else 0.0) for j in range(k))


@dataclass(frozen=True)
class RingCase:
    """A ring of the drive study: band, index, Q, impedance, order, joint count and length."""

    label: str
    f: float
    index: float
    q: float
    z0: float
    m: int
    joints: int
    joint_length: float = 3e-3
    joint_z: float = 1.4
    compensated: bool = False

    @property
    def length(self) -> float:
        return self.m * pt.C0 / (self.index * self.f)

    def ring(
        self,
        k: int,
        rng: np.random.Generator | None,
        at_joints: bool = False,
        tol: float = 0.10,
        jitter: float = 0.5e-3,
    ) -> Ring:
        """The ring with K feeds, its joints identical or drawn: impedance +-tol of Z0 and position +-jitter (sd)."""
        joints = tuple(
            (
                (j + 0.5) * self.length / self.joints + (0.0 if rng is None else float(rng.normal(0, jitter))),
                self.joint_length,
                self.z0 * (self.joint_z + (0.0 if rng is None else float(rng.normal(0, tol)))),
            )
            for j in range(self.joints)
        )
        return Ring(
            self.length,
            self.index,
            self.q,
            self.z0,
            self.m,
            joints,
            feed_layout(k, self.length, self.m, self.joints, at_joints),
            compensated=self.compensated,
        )


DRIVE_CASES = (
    RingCase("915 MHz, m = 1", F915, 1.13, 193, 73, 1, 6),
    RingCase("the same, feed susceptance compensated", F915, 1.13, 193, 73, 1, 6, compensated=True),
    RingCase("2450 MHz, m = 2", 2450e6, 1.20, 128, 66, 2, 6),
    RingCase("must fail: m = 3, three joints, 2450 MHz", 2450e6, 1.20, 128, 66, 3, 3),
)


def compute_drive(quick: bool) -> dict[str, Any]:
    say("\n== drive: the knot as a ring resonator, feeds compiled to a travelling wave ==")
    draws = 3 if quick else 20
    count = 41
    rng = np.random.default_rng(16)
    out: dict[str, Any] = {"cases": {}, "draws": draws}
    began = time.time()
    for case in DRIVE_CASES:
        ks: tuple[int, ...] = (1, 2, 3, 4, 6) if case.m == 3 else (1, 2, 3, 4, 6, 8, 12)
        if case.compensated:
            ks = (2, 6)
        rows: dict[str, Any] = {"K": list(ks), "median": [], "p90": [], "efficiency": []}
        for k in ks:
            vals = [best_drive(case.ring(k, rng), count=count) for _ in range(draws)]
            rows["median"].append(float(np.median([v["ripple"] for v in vals])))
            rows["p90"].append(float(np.percentile([v["ripple"] for v in vals], 90)))
            rows["efficiency"].append(float(np.median([v["efficiency"] for v in vals])))
        joints = [best_drive(case.ring(case.joints, rng, at_joints=True), count=count) for _ in range(draws)]
        rows["at_joints"] = float(np.median([v["ripple"] for v in joints]))
        hybrid = [best_drive(case.ring(2, rng), fixed=np.array([1, -1j]), count=count) for _ in range(draws)]
        rows["hybrid"] = float(np.median([v["ripple"] for v in hybrid]))
        rows["length_mm"] = case.length / MM
        rows["shunt_ohm"] = case.ring(1, None).shunt
        out["cases"][case.label] = rows
        say(
            f"  {case.label} ({case.length / MM:.0f} mm, R = Z0 Q/(2 pi m) = {rows['shunt_ohm']:.0f} ohm), medians of "
            f"{draws} draws: "
            + ", ".join(f"K = {k} {r:.3f}" for k, r in zip(ks, rows["median"], strict=True))
            + f"; feeds at the joints {rows['at_joints']:.3f}; untuned hybrid {rows['hybrid']:.3f}"
        )
    say(f"  [{time.time() - began:.0f} s]")
    # the singular values of G for K = 6 at 915 MHz, at the frequency where the largest peaks
    out["singular"] = {}
    for case in DRIVE_CASES[:2]:
        ring = case.ring(6, None)
        freqs = ring.resonance() * (1 + np.linspace(-0.02, 0.02, 201 if quick else 801))
        sv = np.array([np.linalg.svd(g, compute_uv=False) for g in ring.operator(freqs)]) / math.sqrt(ring.samples)
        peak = int(np.argmax(sv[:, 0]))
        out["singular"][case.label] = {
            "values": sv[peak],
            "f_MHz": float(freqs[peak]) / 1e6,
            "predicted": math.sqrt(6) * ring.shunt / 2,
        }
        say(
            f"  K = 6, {case.label}: singular values {np.round(sv[peak], 1).tolist()} V/A at {freqs[peak] / 1e6:.1f} "
            f"MHz; sqrt(K) R/2 = {math.sqrt(6) * ring.shunt / 2:.0f}"
        )
    case = DRIVE_CASES[0]
    # one draw, three drives: one feed, the untuned hybrid, two compiled feeds
    draw = np.random.default_rng(3)
    shown = {}
    joints_draw = case.ring(2, draw)
    for name, k, fixed in (
        ("one feed", 1, None),
        ("untuned hybrid", 2, np.array([1, -1j])),
        ("two compiled feeds", 2, None),
    ):
        ring_k = Ring(
            joints_draw.length,
            joints_draw.index,
            joints_draw.q,
            joints_draw.z0,
            joints_draw.m,
            joints_draw.joints,
            feed_layout(k, joints_draw.length, 1, 6, False),
        )
        res = best_drive(ring_k, fixed=fixed, count=count)
        shown[name] = {
            "V": np.abs(res["v"]),
            "ripple": res["ripple"],
            "efficiency": res["efficiency"],
            "feeds": ring_k.feeds,
        }
        say(f"  one draw, {name}: ripple {res['ripple']:.3f}, efficiency {res['efficiency']:.2f}")
    out["shown"] = shown
    out["shown_length"] = joints_draw.length
    out["shown_joints"] = [j[0] for j in joints_draw.joints]
    return out


def figure_drive(res: dict[str, Any], out: Path) -> None:
    """Figure 6: the trefoil coloured by |V(s)| for one feed, an untuned hybrid and two compiled feeds; ripple against
    K with the m = 3 case; the singular values of G."""
    plt = _pyplot()
    fig = plt.figure(figsize=(15.5, 9.6))
    grid = fig.add_gridspec(2, 3, hspace=0.34, wspace=0.26, height_ratios=[1.0, 1.05])
    length = res["shown_length"]
    knot = torus(RATIO915, length, 240)
    p = knot.points / MM
    from matplotlib.collections import LineCollection

    for k, (name, row) in enumerate(res["shown"].items()):
        ax = fig.add_subplot(grid[0, k])
        v = row["V"]
        seg = np.stack([p[:, :2], np.roll(p[:, :2], -1, 0)], axis=1)
        lc = LineCollection(list(seg), cmap=_ramp(), norm=plt.Normalize(0, 1.6), linewidths=4.5)
        lc.set_array(v)
        ax.add_collection(lc)
        for f in row["feeds"]:
            i = int(f / length * len(p)) % len(p)
            ax.plot(p[i, 0], p[i, 1], "s", color=RED, ms=7)
        for jpos in res["shown_joints"]:
            i = int(jpos / length * len(p)) % len(p)
            ax.plot(p[i, 0], p[i, 1], "o", color=INK, ms=4)
        ax.set_aspect("equal")
        ax.set_xlim(p[:, 0].min() - 4, p[:, 0].max() + 4)
        ax.set_ylim(p[:, 1].min() - 4, p[:, 1].max() + 4)
        ax.set_xlabel("mm")
        ax.set_title(f"{name}: ripple {row['ripple']:.3f}, {row['efficiency']:.2f} of the power in")
        if k == 2:
            fig.colorbar(lc, ax=ax, fraction=0.045, pad=0.02, label="|V(s)| over its rms")
    ax = fig.add_subplot(grid[1, :2])
    slots = [1, 2, 3, 4, 6, 8, 12]
    for k, (name, row) in enumerate(res["cases"].items()):
        colour = SERIES[[0, 2, 1, 7][k]]
        xs = [slots.index(kk) for kk in row["K"]]
        ax.semilogy(
            xs,
            row["median"],
            "o-",
            color=colour,
            label=f"{name} ({row['length_mm']:.0f} mm, R {row['shunt_ohm']:.0f} Ω)",
        )
        ax.fill_between(xs, row["median"], row["p90"], color=colour, alpha=0.12, lw=0)
        ax.plot([len(slots) + 0.2], [row["at_joints"]], "D", color=colour, ms=7)
        ax.plot([len(slots) + 1.4], [row["hybrid"]], "s", color=colour, ms=7, mfc="none", mew=1.6)
    ax.set_xticks([*range(len(slots)), len(slots) + 0.2, len(slots) + 1.4])
    ax.set_xticklabels([*map(str, slots), "feeds at\nthe joints", "untuned\nhybrid"])
    ax.set_xlabel("feeds K (medians of the draws; shaded to the 90th percentile)")
    ax.set_ylabel("ripple of |V| round the ring")
    ax.set_ylim(2e-3, 1.5)
    ax.set_xlim(-0.4, len(slots) + 2.0)
    ax.legend(fontsize=8, loc="lower left")
    ax.set_title(
        f"ripple against the number of feeds, {res['draws']} tolerance draws a point (joints ±10 % in impedance, ±0.5 "
        "mm)"
    )
    ax = fig.add_subplot(grid[1, 2])
    for k, (name, sv) in enumerate(res["singular"].items()):
        colour = SERIES[[0, 2][k]]
        idx = np.arange(1, len(sv["values"]) + 1)
        ax.semilogy(
            idx,
            sv["values"],
            "o-",
            color=colour,
            label=f"{'compensated' if 'compensated' in name else 'capacitor-coupled'} feeds, {sv['f_MHz']:.1f} MHz",
        )
        for i, val in enumerate(sv["values"]):
            ax.annotate(
                f"{val:.1f}",
                (i + 1, val),
                xytext=(4, 4 if k else -11),
                textcoords="offset points",
                fontsize=7.2,
                color=colour,
            )
    pred = next(iter(res["singular"].values()))["predicted"]
    ax.axhline(pred, color=MUTED, ls=":", lw=1.0)
    ax.text(6, pred * 1.2, "√K·R/2", fontsize=8, color=INK2, ha="right")
    ax.set_xlabel("singular value index")
    ax.set_ylabel("V/A (rms over the ring)")
    ax.legend(fontsize=7.6, loc="center right")
    comp = res["singular"]["the same, feed susceptance compensated"]["values"]
    ax.set_title(f"G of six feeds at 915 MHz: two cheap modes, ±m,\n{comp[0] / comp[2]:.0f} times the floor")
    first = res["cases"]["915 MHz, m = 1"]
    _suptitle(
        fig,
        f"Driving the knot: one feed makes a standing wave (ripple {first['median'][0]:.2f}), two compiled feeds a "
        "travelling "
        f"wave ({first['median'][1]:.3f} in the median, an untuned hybrid {first['hybrid']:.2f}); m = 3 on a trefoil "
        "with three "
        "strong joints stays standing whatever the feeds, the case that must fail",
    )
    _save(fig, out, "drive")


# ---- laminate (figure 7) ------------------------------------------------------------------------------------------
COPPER_SIGMA = 5.8e7  # S/m, bare copper
COPPER_T = 18e-6  # m, the flex fab's copper: the cutoff of the strips' edge fields
POLYIMIDE = (3.3, 0.004)
FR4 = (4.4, 0.02)
LAMINATES = {"copper alone": (0.0, 1.0, 0.0), "polyimide 25 µm": (25e-6, *POLYIMIDE), "FR-4 0.8 mm": (0.8e-3, *FR4)}
# ECSS-E-ST-20-01C, the lowest two-surface multipactor threshold for copper: f x d (GHz mm) against peak volts
ECSS_CU = (
    (0.47, 36.6), (0.50, 32.5), (0.56, 29.7), (0.62, 28.1), (0.70, 27.4), (0.82, 28.7), (0.92, 31.9), (1.02, 37.0),
    (1.14, 44.4), (1.28, 54.2), (1.43, 67.0), (1.59, 83.3), (1.68, 91.9), (1.78, 101.4), (1.88, 110.8), (1.99, 120.9),
    (2.10, 130.1), (2.22, 139.7), (2.34, 148.2), (2.48, 156.7), (2.62, 160.2), (2.77, 154.3), (2.92, 149.2),
    (3.09, 150.6), (3.27, 155.2), (3.45, 163.0), (3.65, 171.7), (3.85, 183.4), (4.07, 195.8), (4.30, 211.2),
    (4.55, 227.6), (4.81, 244.6), (5.08, 262.5), (5.37, 277.3), (5.67, 290.7), (5.99, 302.9), (6.33, 314.4),
    (6.69, 329.6), (7.07, 350.3), (7.47, 372.1), (7.90, 394.5), (8.34, 418.0), (8.82, 441.6), (9.32, 466.6),
    (9.85, 493.1), (10.41, 521.0),
)  # fmt: skip


def strip_section(width: float, split: float, laminate: float, er: float, cutoff: float) -> dict[str, float]:
    """Facing strips `width` wide split by `split` at y = +-1 (lengths in h; halves at +-1/2, so V = 1), laminate
    `laminate` thick of permittivity er behind the copper, 2h wider than it; solved on the quarter x, y >= 0 (odd in
    both) inside grounded walls 8h out, by finite volumes on a grid `cutoff` fine at the copper's edges. Returns the
    capacitance per metre over eps0 in air and with the laminate, the laminate's share of the field energy, the
    index and impedance of the quadrupole mode, and the air problem's sum of E_n^2 dl over the copper's faces (per
    V^2/h), which sets the copper's loss (the surface current of the air problem's charge)."""
    xmax, ymax = width / 2 + 8.0, 1.0 + laminate + 8.0
    x = nodes(0.0, xmax, [split / 2, width / 2, width / 2 + 2.0], cutoff, 0.4, 0.25)
    yfeats = [1.0] if laminate == 0 else [1.0, 1.0 + laminate]
    y = nodes(0.0, ymax, yfeats, min(cutoff, laminate / 3) if laminate > 0 else cutoff, 0.4, 0.25)
    gx, gy = np.meshgrid(x, y, indexing="ij")
    copper = (np.abs(gy - 1.0) < 1e-12) & (gx >= split / 2 - 1e-12) & (gx <= width / 2 + 1e-12)
    fixed = copper | (gx == 0) | (gy == 0) | (gx >= xmax - 1e-12) | (gy >= ymax - 1e-12)
    value = np.where(copper, -0.5, 0.0)
    xc, yc = (x[1:] + x[:-1]) / 2, (y[1:] + y[:-1]) / 2
    cx, cy = np.meshgrid(xc, yc, indexing="ij")
    lam = (cy > 1.0) & (cy < 1.0 + laminate) & (cx < width / 2 + 2.0)
    ones = np.ones_like(cx)
    eps = np.where(lam, er, 1.0)
    phi_air = laplace2d(x, y, ones, fixed, value)
    c_air = 4 * energy2d(x, y, ones, phi_air)
    if laminate > 0:
        phi_eps = laplace2d(x, y, eps, fixed, value)
        c_eps = 4 * energy2d(x, y, eps, phi_eps)
        share = 4 * energy2d(x, y, np.where(lam, er, 0.0), phi_eps) / c_eps
    else:
        c_eps, share = c_air, 0.0
    jc = int(np.argmin(np.abs(y - 1.0)))
    on = (x >= split / 2 - 1e-12) & (x <= width / 2 + 1e-12)
    below = (phi_air[on, jc] - phi_air[on, jc - 1]) / (y[jc] - y[jc - 1])
    above = (phi_air[on, jc + 1] - phi_air[on, jc]) / (y[jc + 1] - y[jc])
    xi = x[on]
    dl = np.gradient(xi) if len(xi) > 1 else np.ones(1)
    dl[0] /= 2
    dl[-1] /= 2
    surface = 4 * float(((below**2 + above**2) * dl).sum())
    return {
        "c_air": c_air,
        "c_eps": c_eps,
        "share": share,
        "index": math.sqrt(c_eps / c_air),
        "z0": 1 / (pt.C0 * pt.EPS0 * math.sqrt(c_eps * c_air)),
        "surface": surface,
    }


def section_q(sec: dict[str, float], f: float, tand: float, resistance: float = 1.0) -> dict[str, float]:
    """A section's Q at frequency f with h = f h / f: the copper's, 2 c h / (delta S) for surface resistance
    `resistance` times bare copper's, the laminate's, 1 / (share tan delta), and both."""
    h = FH / f
    delta = math.sqrt(2 / (2 * math.pi * f * pt.MU0 * COPPER_SIGMA))
    qc = 2 * sec["c_air"] * h / (delta * sec["surface"]) / resistance
    qd = 1 / (sec["share"] * tand) if sec["share"] * tand > 0 else math.inf
    return {"Qc": qc, "Qd": qd, "Q": 1 / (1 / qc + 1 / qd)}


@dataclass(frozen=True)
class Yield:
    """Secondary emission yield against impact energy (eV): the modified Vaughan curve with copper's parameters
    (first crossover e1, peak dmax at emax), with a floor d0 below 50 eV for elastic reflection."""

    e1: float = 19.0
    dmax: float = 2.48
    emax: float = 232.0
    d0: float = 0.8

    def _vaughan(self, e: Arr, e0: float) -> Arr:
        v = np.clip((e - e0) / (self.emax - e0), 1e-12, None)
        k = np.where(v < 1, 0.56, 0.25)
        out: Arr = self.dmax * (v * np.exp(1 - v)) ** k
        return out

    @property
    def e0(self) -> float:
        lo, hi = 0.0, self.e1 - 1e-6
        for _ in range(100):
            mid = (lo + hi) / 2
            if float(self._vaughan(np.array([self.e1]), mid)[0]) > 1.0:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2

    def __call__(self, e: Arr) -> Arr:
        out: Arr = np.maximum(self._vaughan(e, self.e0), self.d0 * (e < 50))
        return out


def multipactor_growth(fd: float, volts: float, sey: Yield, n: int, periods: float, seed: int = 1) -> float:
    """Two-surface multipactor between plates d = 1 mm apart at f = fd / d, in one dimension, by Monte Carlo with
    weights: secondaries leave with gamma-distributed energy (mean 3 eV) and a cosine law, each impact multiplies a
    particle's weight by the yield; a flight that meets no plate within three RF periods is lost. log10 of the mean
    weight after `periods` RF periods, floored at -10: above zero the population grows."""
    rng = np.random.default_rng(seed)
    d = 1e-3
    w = 2 * math.pi * fd * 1e6 / d
    qm = pt.E / pt.ME
    amp = qm * volts / (d**2 * w**2)
    vscale = w * d
    phase = rng.uniform(0, 2 * math.pi, n)
    elapsed = np.zeros(n)
    logw = np.zeros(n)
    alive = np.ones(n, dtype=bool)
    tmax = periods * 2 * math.pi
    steps = 0.02 * np.arange(1, int(6 * math.pi / 0.02) + 1)
    while alive.any():
        idx = np.flatnonzero(alive)
        m = len(idx)
        energy = rng.gamma(2.0, 1.5, m)
        v0 = np.sqrt(2 * energy * pt.E / pt.ME) * np.sqrt(rng.uniform(0, 1, m)) / vscale
        u0 = phase[idx]
        y = v0[:, None] * steps[None, :] + amp * (
            steps[None, :] * np.cos(u0)[:, None] - np.sin(u0[:, None] + steps[None, :]) + np.sin(u0)[:, None]
        )
        hit = (y >= 1.0) | (y <= 0.0)
        done = hit.any(axis=1)
        k = np.where(done, hit.argmax(axis=1), 0)
        rows = np.arange(m)
        y1 = y[rows, k]
        y0 = np.where(k > 0, y[rows, np.maximum(k - 1, 0)], 0.0)
        far = y1 >= 1.0
        frac = np.where(np.abs(y1 - y0) > 0, (np.where(far, 1.0, 0.0) - y0) / (y1 - y0), 1.0)
        s_hit = np.where(k > 0, steps[np.maximum(k - 1, 0)], 0.0) + 0.02 * np.clip(frac, 0, 1)
        v_imp = v0 + amp * (np.cos(u0) - np.cos(u0 + s_hit))
        e_imp = 0.5 * pt.ME * (v_imp * vscale) ** 2 / pt.E
        logw[idx] += np.where(done, np.log(np.maximum(sey(e_imp), 1e-300)), -np.inf)
        elapsed[idx] += np.where(done, s_hit, tmax)
        phase[idx] = np.mod(u0 + s_hit + np.where(far, math.pi, 0.0), 2 * math.pi)
        alive[idx] = done & (elapsed[idx] < tmax) & np.isfinite(logw[idx])
    mean = float(np.mean(np.exp(np.clip(logw, -700, 700))))
    return max(math.log10(mean), -10.0) if mean > 0 else -10.0  # -10: every electron lost


def growth_job(job: tuple[float, float, int, float]) -> dict[str, Any]:
    """One multipactor run (f d, volts, particles, periods) on copper's yield, for a process pool."""
    fd, volts, count, periods = job
    return {"growth": multipactor_growth(fd, volts, Yield(), count, periods)}


def compute_laminate(quick: bool) -> dict[str, Any]:
    say("\n== laminate: Q against the split and the laminate, the power a 1 eV well takes, multipactor ==")
    out: dict[str, Any] = {}
    splits = (0.1, 0.3, 0.6, 1.0, 1.5) if quick else (0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.25, 1.5)
    bands = (F915, 2450e6)
    began = time.time()
    curves: dict[str, Any] = {}
    for f in bands:
        h = FH / f
        cutoff = COPPER_T / h
        for name, (thick, er, tand) in LAMINATES.items():
            rows = []
            for g in splits:
                sec = strip_section(8.0, g, thick / h, er, cutoff)
                rows.append({"g": g, **sec, **section_q(sec, f, tand)})
            curves[f"{f / 1e6:.0f} MHz, {name}"] = rows
    out["curves"] = curves
    say(f"  Q against g/h on strips 8h wide [{time.time() - began:.0f} s]:")
    for name, rows in curves.items():
        say("    " + name + ": " + ", ".join(f"g/h {r['g']:g} Q {r['Q']:.0f}" for r in rows))
    b915 = next(r for r in curves["915 MHz, polyimide 25 µm"] if abs(r["g"] - 0.3) < 1e-9)
    out["section_B"] = b915
    say(
        f"  section B at 915 MHz (8h wide, split 0.3h, 25 um polyimide): n_eff {b915['index']:.4f}, polyimide share "
        f"{b915['share']:.3f}, Z0 {b915['z0']:.1f} ohm, Q_c {b915['Qc']:.0f} (bare copper), Q {b915['Q']:.0f}; "
        f"m = 1 takes {pt.C0 / (b915['index'] * F915) / MM:.0f} mm"
    )
    # the power a 1 eV well takes on section B against frequency, per order m, and the energy the bend holds
    fb = pt.FacingBoards(2.0)
    drive = pt.electrode_drive_at_q(1.0, 0.5, 1.0, 1.0, fb.depth_share, fb.gradient_share)
    volts = drive.volts
    freqs = (
        np.array([0.433e9, 0.7e9, F915, 1.5e9, 2.45e9, 4.0e9, 5.8e9])
        if not quick
        else np.array([0.433e9, F915, 2.45e9, 5.8e9])
    )
    power: dict[int, list[float]] = {m: [] for m in range(1, 6)}
    index = []
    for f in freqs:
        h = FH / f
        sec = strip_section(8.0, 0.3, 25e-6 / h, *POLYIMIDE[:1], COPPER_T / h)
        q = section_q(sec, f, POLYIMIDE[1])
        index.append(sec["index"])
        for m in power:
            shunt = sec["z0"] * q["Q"] / (2 * math.pi * m)
            power[m].append(volts**2 / (2 * shunt))
    held = {}
    bend_over_length = torus(RATIO915, 1.0, 3000).bend_radius
    for m in power:
        length = m * pt.C0 / (float(np.mean(index)) * F915)
        ratio = bend_over_length * length / (FH / F915)
        held[m] = (0.61 * 0.605 * ratio, 0.95 * 0.605 * ratio)
    out["power"] = {"f": freqs, "P": power, "volts": volts, "held": held, "index": index}
    say(
        f"  1 eV at q = 0.5 takes {volts:.1f} V; on section B: "
        + ", ".join(
            f"{f / 1e9:g} GHz m = {m}: {power[m][i]:.2f} W"
            for m, f, i in (
                (1, F915, list(freqs).index(F915)),
                (2, 2.45e9, list(freqs).index(2.45e9)),
                (5, 5.8e9, list(freqs).index(5.8e9)),
            )
        )
    )
    say(
        "  the bend holds, by the discounted static limit: "
        + ", ".join(f"m = {m}: {lo:.1f} to {hi:.1f} eV per eV" for m, (lo, hi) in held.items())
    )
    # multipactor: ECSS's copper boundary, the operating lines V = 35.7 q (fd)^2 and the one-dimensional model
    out["ecss"] = np.array(ECSS_CU)
    coef_v = 35.7
    coef_d = 1.42
    fd_1ev = math.sqrt(1.0 / (coef_d * 0.25))
    v_1ev = coef_v * 0.5 * fd_1ev**2
    ecss_at = float(np.interp(fd_1ev, out["ecss"][:, 0], out["ecss"][:, 1]))
    out["mpc"] = {
        "fd_1eV": fd_1ev,
        "V_1eV": v_1ev,
        "ecss_at_1eV": ecss_at,
        "margin_db": 20 * math.log10(ecss_at / v_1ev),
    }
    say(
        f"  MPC-2: 1 eV at q = 0.5 is f d = {fd_1ev:.3f} GHz mm and {v_1ev:.1f} V, {out['mpc']['margin_db']:.1f} dB "
        f"below ECSS's copper boundary ({ecss_at:.1f} V)"
    )
    began = time.time()
    count, periods = (150, 20.0) if quick else (500, 40.0)
    boundary = list(ECSS_CU[2:26:4])
    depths = (0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0) if not quick else (1.0, 3.0)
    points = [(qq, d, math.sqrt(d / (coef_d * qq * qq))) for qq in (0.3, 0.5, 0.7) for d in depths]
    jobs = [(fd, v, count, periods) for fd, v in boundary] + [
        (fd, coef_v * qq * fd * fd, count, periods) for qq, _, fd in points
    ]
    grown = [r["growth"] for r in run_jobs(growth_job, jobs)]
    calib = [(fd, v, g) for (fd, v), g in zip(boundary, grown[: len(boundary)], strict=True)]
    lines: dict[float, list[tuple[float, float, float, float]]] = {0.3: [], 0.5: [], 0.7: []}
    for (qq, d, fd), g in zip(points, grown[len(boundary) :], strict=True):
        lines[qq].append((d, fd, coef_v * qq * fd * fd, g))
    out["growth"] = {"calibration": calib, "lines": lines}
    say(
        "  one-dimensional two-surface model on ECSS's copper boundary (log10 growth over the run): "
        + ", ".join(f"{fd:g}: {g:+.2f}" for fd, _, g in calib)
    )
    for qq, line in lines.items():
        say(f"    q = {qq}: " + ", ".join(f"{d:g} eV ({fd:.2f} GHz mm, {v:.0f} V) {g:+.2f}" for d, fd, v, g in line))
    say(f"  [multipactor {time.time() - began:.0f} s]")
    return out


def figure_laminate(res: dict[str, Any], out: Path) -> None:
    """Figure 7: Q against g/h and laminate per band; power for 1 eV against frequency with the bend-limited energy;
    the ECSS boundary with the operating lines at q = 0.3, 0.5 and 0.7 and the model's growth."""
    plt = _pyplot()
    fig = plt.figure(figsize=(15.5, 6.6))
    grid = fig.add_gridspec(1, 3, wspace=0.3, top=0.82)
    ax = fig.add_subplot(grid[0, 0])
    for k, (name, rows) in enumerate(res["curves"].items()):
        g = [r["g"] for r in rows]
        qs = [r["Q"] for r in rows]
        colour = SERIES[k % 3]
        ls = "-" if name.startswith("915") else "--"
        ax.semilogy(g, qs, ls, color=colour, marker="o", ms=3.5, label=name)
        if name == "915 MHz, polyimide 25 µm":
            enig = [1 / (1 / (r["Qc"] / s) + 1 / r["Qd"]) for r, s in ((r, 1.3) for r in rows)]
            enig_low = [1 / (1 / (r["Qc"] / 5.0) + 1 / r["Qd"]) for r in rows]
            ax.fill_between(
                g, enig_low, enig, color=colour, alpha=0.15, lw=0, label="the same on ENIG: R_s 1.3 to 5 times copper's"
            )
    b = res["section_B"]
    ax.plot([0.3], [b["Q"]], "*", color=INK, ms=11)
    ax.annotate(
        f"T915-S: Q {b['Q']:.0f}, n {b['index']:.3f}, Z₀ {b['z0']:.1f} Ω",
        (0.3, b["Q"]),
        xytext=(10, 12),
        textcoords="offset points",
        fontsize=8,
    )
    ax.set_xlabel("split over h, g/h (strips 8h wide)")
    ax.set_ylabel("Q (bare copper)")
    ax.set_ylim(45, 3000)
    ax.legend(fontsize=7.0, loc="upper center", ncol=2)
    _log_ticks(ax.yaxis, [50, 100, 200, 400, 800])
    ax.set_title("Q of the strips' quadrupole mode by laminate\n(solid 915 MHz, dashed 2450 MHz)")
    ax = fig.add_subplot(grid[0, 1])
    pw = res["power"]
    f = np.asarray(pw["f"]) / 1e9
    for m, colour in zip(range(1, 6), SERIES, strict=False):
        lo, hi = pw["held"][m]
        ax.loglog(
            f, pw["P"][m], "o-", color=colour, ms=3.5, label=f"m = {m}: the bend holds {lo:.0f} to {hi:.0f} eV per eV"
        )
    marks = ((F915, 1), (2.45e9, 2), (5.8e9, 5))
    for fm, m in marks:
        i = int(np.argmin(np.abs(np.asarray(pw["f"]) - fm)))
        ax.plot([f[i]], [pw["P"][m][i]], "*", color=INK, ms=11)
        ax.annotate(
            f"{pw['P'][m][i]:.2f} W", (f[i], pw["P"][m][i]), xytext=(6, -12), textcoords="offset points", fontsize=8
        )
    ax.set_xlabel("frequency (GHz)")
    ax.set_ylabel(f"power for a 1 eV well at q = 0.5 ({pw['volts']:.1f} V)")
    ax.legend(fontsize=7.4, loc="upper left")
    _log_ticks(ax.xaxis, [0.433, 0.915, 2.45, 5.8])
    _log_ticks(ax.yaxis, [0.3, 1, 3, 10])
    ax.set_title("the flex strips (section B, bare copper)\nas a ring of m guided waves")
    ax = fig.add_subplot(grid[0, 2])
    ecss = res["ecss"]
    ax.loglog(ecss[:, 0], ecss[:, 1], color=INK, lw=2.0, label="ECSS-E-ST-20-01C, copper")
    fd = np.geomspace(0.4, 10, 200)
    growth = res["growth"]
    for q, colour in zip((0.3, 0.5, 0.7), (BLUE, ORANGE, GREEN), strict=True):
        ax.loglog(fd, 35.7 * q * fd**2, color=colour, lw=1.3, label=f"q = {q}: V = 35.7·q·(fd)², D = 1.42·q²·(fd)²")
        for depth, x, v, grows in growth["lines"][q]:
            ax.plot([x], [v], "o", color=colour, mfc=colour if grows > 0 else SURFACE, ms=6, mew=1.4)
            if q == 0.5:
                ax.annotate(
                    f"{depth:g} eV", (x, v), xytext=(5, -10), textcoords="offset points", fontsize=7.4, color=INK2
                )
    for x, v, grows in growth["calibration"]:
        ax.plot([x], [v], "s", color=INK, mfc=INK if grows > 0 else SURFACE, ms=5)
    ax.set_xlim(0.4, 10)
    ax.set_ylim(5, 600)
    ax.set_xlabel("f × d (GHz·mm), d = 2h")
    ax.set_ylabel("peak volts across the gap")
    ax.legend(fontsize=7.2, loc="upper left")
    mpc = res["mpc"]
    _log_ticks(ax.xaxis, [0.5, 1, 2, 5, 10])
    _log_ticks(ax.yaxis, [10, 30, 100, 300])
    ax.set_title(
        f"multipactor: at 1 eV and q = 0.5 the line sits {mpc['margin_db']:.1f} dB below\ncopper's boundary; filled, "
        "the 1-D model grows; open, decays"
    )
    _suptitle(
        fig,
        f"Laminate and band: polyimide 25 µm holds {b['share'] * 100:.1f} % of the strips' field at 915 MHz (n_eff "
        f"{b['index']:.3f}, "
        f"Q {b['Q']:.0f} on bare copper), so a 1 eV well takes {pw['P'][1][list(np.asarray(pw['f'])).index(F915)]:.2f} "
        "W at m = 1; the "
        "power grows as m√f while the energy the bend holds grows with m alone",
    )
    _save(fig, out, "laminate")


# ---- tracking: jobs for a process pool ----------------------------------------------------------------------------
WORKERS = 1  # set by main: processes for the independent tracking runs


def run_jobs(func: Callable[[Any], dict[str, Any]], jobs: Sequence[Any]) -> list[dict[str, Any]]:
    """Each job through func, in WORKERS processes when there are several."""
    if WORKERS <= 1 or len(jobs) <= 1:
        return [func(job) for job in jobs]
    with ProcessPoolExecutor(max_workers=min(WORKERS, len(jobs))) as pool:
        return list(pool.map(func, jobs))


def mathieu_beta(q: float, steps: int = 4000) -> float:
    """Mathieu's characteristic exponent beta for a = 0: x'' + 2q cos(2 tau) x = 0 over a period, by RK4, from the
    trace of the monodromy matrix. The secular frequency is beta f / 2 (q / sqrt 2 to first order)."""
    h = math.pi / steps

    def rhs(tau: float, y: Arr) -> Arr:
        out: Arr = np.array([y[1], -2 * q * math.cos(2 * tau) * y[0]])
        return out

    cols = []
    for y0 in (np.array([1.0, 0.0]), np.array([0.0, 1.0])):
        y, tau = y0.copy(), 0.0
        for _ in range(steps):
            k1 = rhs(tau, y)
            k2 = rhs(tau + h / 2, y + h / 2 * k1)
            k3 = rhs(tau + h / 2, y + h / 2 * k2)
            k4 = rhs(tau + h, y + h * k3)
            y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            tau += h
        cols.append(y)
    return math.acos((cols[0][0] + cols[1][1]) / 2) / math.pi


def rods_curve(name: str) -> ClosedCurve:
    """The first run's trefoil, or the flat ring of its length."""
    tre = first_run_trefoil()
    return tre if name == "trefoil" else ClosedCurve.circle(tre.length / (2 * math.pi), len(tre))


def rods_job(job: tuple[str, float, tuple[float, ...], float, float, float, float, float]) -> dict[str, Any]:
    """Electrons through four thin rods 0.5 mm from the first run's trefoil (or the ring of its length) at 1 GHz and
    Mathieu q: started on the bend's static equilibrium (matched) or on the axis, with the rods' line charge dipped
    by `dip` over joints `joint` long (six on one orbit of the knot) and the null stepped `step` aside there, tracked
    `laps` laps at 30 steps an RF period; held, laps and reach per energy."""
    name, q, energies, laps, joint, dip, step, matched = job
    f, r0 = 1e9, 0.5e-3
    curve = rods_curve(name)
    rods = QuadRods.along(curve, r0)
    g0 = pt.gradient_for_q(q, f)
    lam0 = rods.line_charge(g0)
    w = np.array(energies, dtype=np.float64)
    count = len(w)
    v = np.array([pt.v_of(float(e)) for e in w])
    gamma = 1 / np.sqrt(1 - (v / pt.C0) ** 2)
    watch = GuideWatch(curve, 0.8 * r0, count)
    near = np.arange(-60, 61, dtype=np.intp)
    omega = 2 * math.pi * f
    length, ds = curve.length, curve.length / len(curve)
    tan, axis_a = curve.tangents, rods.frame.a
    centres = (np.arange(6) + 0.5) * length / 6

    def field(x: Arr, t: float) -> tuple[Arr, None]:
        st = watch.station
        if joint > 0:
            s = (st * ds + ((x - curve.points[st]) * tan[st]).sum(1)) % length
            dz = (s[:, None] - centres[None, :] + length / 2) % length - length / 2
            prof = np.where(np.abs(dz) < joint / 2, np.cos(math.pi * dz / joint) ** 2, 0.0).sum(1)
            e = rods.field(x, st, lam0 * (1 - dip * prof), near)
            if step:
                e = e + (g0 * step * prof)[:, None] * axis_a[st]
        else:
            e = rods.field(x, st, np.full(count, lam0), near)
        return e * math.cos(omega * t), None

    x0 = np.tile(curve.points[0], (count, 1)) + 1e-5 * axis_a[0]
    if matched:
        depth = pt.quadrupole_depth_ev(q, f, r0)
        kvec = (np.roll(curve.points, -1, 0) - 2 * curve.points + np.roll(curve.points, 1, 0))[0] / ds**2
        x0 = np.tile(curve.points[0], (count, 1)) - (w * r0**2 / depth)[:, None] * kvec[None, :]
    dt = 1 / (f * 30)
    flight = pt.integrate(
        x0,
        (gamma * v)[:, None] * tan[0][None, :],
        np.full(count, pt.QM_E),
        dt,
        int(laps * length / v.min() / dt),
        field,
        [watch],
    )
    return {"energy": w, "held": flight.alive, "laps": watch.laps, "reach_um": watch.reach * 1e6}


# The unrolled first-run guide: an ideal quadrupole channel of the trefoil's length at 1 GHz and q = 0.6, r0 = 0.5 mm
UNROLL_F, UNROLL_R0, UNROLL_Q, UNROLL_L, UNROLL_SPP = 1e9, 0.5e-3, 0.6, 0.277, 30


def scheme_secular(q: float, f: float, spp: int) -> float:
    """The secular frequency (Hz) of the integrator itself in an ideal quadrupole at Mathieu q: the one-period matrix
    of its kick-drift steps (the field taken at each step's start)."""
    dt = 1 / (f * spp)
    omega = 2 * math.pi * f
    g0 = pt.gradient_for_q(q, f)
    k = pt.QM_E * dt
    cols = []
    for x, u in ((1.0, 0.0), (0.0, 1.0)):
        for s in range(spp):
            u = u + k * (-g0 * x) * math.cos(omega * s * dt)
            x = x + u * dt
        cols.append((x, u))
    return math.acos((cols[0][0] + cols[1][1]) / 2) / (2 * math.pi) * f


def joint_profile(centres: Arr, joint: float, length: float, s: Arr) -> tuple[Arr, Arr]:
    """The joints' profile w(s), each a cos^2 bump `joint` long, and its slope."""
    dz = (s[:, None] - centres[None, :] + length / 2) % length - length / 2
    inside = np.abs(dz) < joint / 2
    w = np.where(inside, np.cos(math.pi * dz / joint) ** 2, 0.0).sum(1)
    dw = np.where(inside, -math.pi / joint * np.sin(2 * math.pi * dz / joint), 0.0).sum(1)
    return w, dw


def unrolled_centres(count: int, spacing: str, seed: int) -> Arr:
    """J joint centres along the unrolled guide: equally spaced, or each moved up to 0.3 of a spacing at random."""
    centres = (np.arange(count) + 0.5) * UNROLL_L / count
    if spacing == "irregular":
        centres = np.sort(
            (centres + np.random.default_rng(seed).uniform(-0.3, 0.3, count) * UNROLL_L / count) % UNROLL_L
        )
    out: Arr = centres
    return out


def unrolled_job(job: tuple[int, str, float, float, float, tuple[float, ...], float, int]) -> dict[str, Any]:
    """Electrons along the unrolled first-run guide with J joints, each a gradient dip `dip` and a null step `step`
    over `joint` (cos^2), the axial field the dip and step imply included; reach per energy (section 6.2)."""
    count, spacing, joint, dip, step, energies, laps, seed = job
    centres = unrolled_centres(count, spacing, seed)
    g0 = pt.gradient_for_q(UNROLL_Q, UNROLL_F)
    omega = 2 * math.pi * UNROLL_F
    w_ev = np.array(energies, dtype=np.float64)
    n = len(w_ev)
    v = np.array([pt.v_of(float(e)) for e in w_ev])
    use = count > 0 and spacing != "none"

    def field(x: Arr, t: float) -> tuple[Arr, None]:
        c = math.cos(omega * t)
        if not use:
            return np.stack([-g0 * x[:, 0] * c, g0 * x[:, 1] * c, 0 * x[:, 0]], 1), None
        w, dw = joint_profile(centres, joint, UNROLL_L, x[:, 2] % UNROLL_L)
        g = g0 * (1 - dip * w)
        xs = x[:, 0] - step * w
        ez = g0 * dip * dw * (xs * xs - x[:, 1] ** 2) / 2 + g * xs * step * dw
        return np.stack([-g * xs * c, g * x[:, 1] * c, ez * c], 1), None

    class Reach:
        def __init__(self) -> None:
            self.r2 = np.zeros(n)

        def __call__(self, step: int, t: float, x: Arr, u: Arr, alive: BoolArr) -> BoolArr:
            r2 = x[:, 0] ** 2 + x[:, 1] ** 2
            self.r2 = np.where(alive, np.maximum(self.r2, r2), self.r2)
            out: BoolArr = r2 > (0.8 * UNROLL_R0) ** 2
            return out

    a0 = 10e-6
    x0 = np.zeros((n, 3))
    x0[:, 0], x0[:, 1], x0[:, 2] = a0, a0 / 2, UNROLL_L / (2 * max(count, 1)) * 0.37
    u0 = np.zeros((n, 3))
    u0[:, 2] = v
    watch = Reach()
    dt = 1 / (UNROLL_F * UNROLL_SPP)
    flight = pt.integrate(x0, u0, np.full(n, pt.QM_E), dt, int(laps * UNROLL_L / v.min() / dt), field, [watch])
    return {"energy": w_ev, "reach_um": np.sqrt(watch.r2) * 1e6, "held": flight.alive}


def predicted_bands(
    centres: Arr, joint: float, dip: float, step: float, lo: float, hi: float, fs: float
) -> list[dict[str, Any]]:
    """The resonances first-order theory puts in [lo, hi] eV (section 6.2): with g_k the lap harmonic k of the dip and
    d_k that of the null step, a gradient dip is resonant at k f_rev = 2 f_s and f -+ 2 f_s, a null step at f_s and
    f -+ f_s; each with its strength."""
    n = 1 << 15
    s = np.arange(n) * UNROLL_L / n
    w, _ = joint_profile(centres, joint, UNROLL_L, s)
    cg = np.fft.fft(dip * w) / n
    cd = np.fft.fft(step * w) / n
    families = []
    if dip:
        families += [("2f_s", 2 * fs, cg), ("f − 2f_s", UNROLL_F - 2 * fs, cg), ("f + 2f_s", UNROLL_F + 2 * fs, cg)]
    if step:
        families += [("f_s", fs, cd), ("f − f_s", UNROLL_F - fs, cd), ("f + f_s", UNROLL_F + fs, cd)]
    out = []
    for name, freq, coef in families:
        for k in range(1, n // 2):
            v = freq * UNROLL_L / k
            energy = 0.5 * pt.ME * v * v / pt.E
            if energy < lo:
                break
            if energy <= hi and abs(coef[k]) > 1e-9:
                out.append({"family": name, "k": k, "energy": energy, "strength": float(abs(coef[k]))})
    return out


# ---- held (figure 9) ----------------------------------------------------------------------------------------------
def boards_well(q0: float = 0.5) -> dict[str, Any]:
    """The facing boards' well and force (ELE-2, FacingBoards.field) in their plane, D tanh^2(pi x / 2h), and toward a
    board, D tan^2(pi y / 2h), which has no barrier: the local Mathieu q, q0 sec^2(pi y / 2h), reaches 0.908 first.
    Beside them four rods' well along their weak diagonal, D u^2 / (1 + u^4)^2 (JNT-2). Lengths in h, forces in D/h."""
    fb = FacingBoards(2.0)
    u = np.linspace(0, 3.0, 3001)
    inplane = (fb.field(u, 0 * u, 1.0) ** 2).sum(-1) / fb.barrier_field(1.0) ** 2
    yy = np.linspace(0, 0.6, 1201)
    toward = (fb.field(0 * yy, yy, 1.0) ** 2).sum(-1) / fb.barrier_field(1.0) ** 2
    force_in = np.gradient(inplane, u)
    force_to = np.gradient(toward, yy)
    local_q = q0 / np.cos(math.pi * yy / 2) ** 2
    edge = float(yy[np.argmax(local_q >= 0.908)])
    rods_u = np.linspace(1e-4, 1.2, 200001)
    rods_well = rods_u**2 / (1 + rods_u**4) ** 2
    return {
        "x": u,
        "inplane": inplane,
        "force_in": force_in,
        "y": yy,
        "toward": toward,
        "force_to": force_to,
        "local_q": local_q,
        "edge": edge,
        "force_max_in": float(force_in.max()),
        "force_at_edge": float(np.interp(edge, yy, force_to)),
        "rods_force_max": float(np.gradient(rods_well, rods_u).max()),
    }


def compute_held(quick: bool) -> dict[str, Any]:
    say("\n== held: the energy a bend holds, the facing boards' well, the energy held along the knot ==")
    out: dict[str, Any] = {}
    plan = {
        ("trefoil", 0.3): (2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0),
        ("trefoil", 0.6): (8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0),
        ("ring", 0.3): (4.0, 6.0, 8.0, 10.0, 11.0, 12.0, 14.0),
        ("ring", 0.6): (16.0, 20.0, 24.0, 28.0, 32.0, 36.0, 40.0),
    }
    if quick:
        plan = {key: es[:: max(len(es) // 3, 1)] for key, es in plan.items()}
    laps = 2.0 if quick else 6.0
    began = time.time()
    results = run_jobs(rods_job, [(name, q, es, laps, 0.0, 0.0, 0.0, 1.0) for (name, q), es in plan.items()])
    tracked = {}
    r0 = 0.5e-3
    for ((name, q), _), res in zip(plan.items(), results, strict=True):
        curve = rods_curve(name)
        depth = pt.quadrupole_depth_ev(q, 1e9, r0)
        rule = pt.bend_limit_ev(depth, curve.bend_radius, r0)
        held = [float(e) for e, a in zip(res["energy"], res["held"], strict=True) if a]
        lost = [float(e) for e, a in zip(res["energy"], res["held"], strict=True) if not a]
        row = {
            **res,
            "depth": depth,
            "rule": rule,
            "highest_held": max(held) if held else 0.0,
            "lowest_lost": min(lost) if lost else math.inf,
        }
        tracked[f"{name}, q = {q}"] = row
        say(
            f"  four rods on the {name} ({curve.length / MM:.1f} mm, bend {curve.bend_radius / MM:.1f} mm), q = {q}: "
            f"depth {depth:.2f} eV, "
            f"rule {rule:.1f} eV; held up to {row['highest_held']:g} eV ({row['highest_held'] / rule:.2f} of the rule, "
            f"{row['highest_held'] / (0.70 * rule):.2f} of the static limit), lost from {row['lowest_lost']:g} eV"
        )
    say(f"  [{time.time() - began:.0f} s on {WORKERS} processes]")
    out["tracked"] = tracked
    out["boards"] = boards_well()
    say(
        f"  static force maxima: four rods' weak diagonal {out['boards']['rods_force_max']:.4f} D/r0, facing boards in "
        "their plane "
        f"{out['boards']['force_max_in']:.3f} D/h; toward a board q0 sec^2 reaches 0.908 at y = "
        f"{out['boards']['edge']:.3f} h for "
        f"q0 = 0.5, where the force is {out['boards']['force_at_edge']:.2f} D/h"
    )
    # the energy held against curvature along T915-S, per eV of depth, by the in-plane static limit
    knot = torus(RATIO915, L915, 1500 if quick else 4000)
    kappa = np.linalg.norm(curvature_vector(knot), axis=1)
    along = {"s": np.arange(len(knot)) / len(knot), "whole bend in the plane": 0.605 / (H915 * kappa)}
    for k in (0, -1):
        fr = framing(knot, k)
        along[f"k = {k}"] = 0.605 / (H915 * np.abs(fr.kappa_g))
        along[f"k = {k}, toward"] = 2 * 0.605 / (H915 * np.maximum(np.abs(fr.kappa_n), 1e-9))
    out["along"] = along
    say(
        "  T915-S by the in-plane static limit 0.605 D R/h: "
        + ", ".join(
            f"{name}: least {float(np.min(v)):.1f} eV per eV"
            for name, v in along.items()
            if name != "s" and "toward" not in name
        )
    )
    return out


def figure_held(res: dict[str, Any], out: Path) -> None:
    """Figure 9: held and lost against energy for the trefoil and the ring with the rule and the static and sudden
    lines; the facing-board well and force in their plane and toward a board; the energy held along the knot."""
    plt = _pyplot()
    fig = plt.figure(figsize=(15.5, 9.6))
    grid = fig.add_gridspec(2, 3, hspace=0.36, wspace=0.28)
    ax = fig.add_subplot(grid[0, :2])
    rows = list(res["tracked"].items())
    for k, (name, row) in enumerate(rows):
        e = np.asarray(row["energy"]) / row["rule"]
        held = np.asarray(row["held"], dtype=bool)
        y = np.full(len(e), k)
        ax.plot(e[held], y[held], "o", color=SERIES[k], ms=8, label=None)
        ax.plot(e[~held], y[~held], "o", color=SERIES[k], mfc=SURFACE, ms=8, mew=1.6)
        lost_from = f"lost from {row['lowest_lost']:g} eV" if math.isfinite(row["lowest_lost"]) else "none lost"
        ax.text(
            1.16,
            k,
            f"{name}: rule {row['rule']:.1f} eV, held to {row['highest_held']:g} eV, {lost_from}",
            va="center",
            fontsize=8.2,
        )
    for x, label, colour in (
        (1.0, "bend rule D·R/2r₀", INK),
        (0.70, "static 0.350·D·R/r₀ (0.70)", ORANGE),
        (0.47, "sudden (0.47)", BLUE),
    ):
        ax.axvline(x, color=colour, lw=1.2, ls="--" if x < 1 else "-")
        ax.text(x, len(rows) - 0.4, label, rotation=90, va="top", ha="right", fontsize=8, color=colour)
    ax.set_xlim(0.2, 2.0)
    ax.set_ylim(-0.6, len(rows) - 0.3)
    ax.set_yticks([])
    ax.set_xlabel("energy over the bend rule's estimate (filled: held, open: lost)")
    ax.set_title(
        "four rods at 1 GHz, r₀ = 0.5 mm, on the first run's trefoil and on the ring of its length, started matched to "
        "the bend"
    )
    b = res["boards"]
    ax = fig.add_subplot(grid[0, 2])
    ax.plot(b["x"], b["inplane"], color=BLUE, label="in the plane: D·tanh²(πx/2h)")
    ax.plot(b["x"], b["force_in"], color=BLUE, ls="--", label=f"its force, at most {b['force_max_in']:.3f} D/h")
    ax.plot(b["y"], b["toward"], color=ORANGE, label="toward a board: D·tan²(πy/2h)")
    ax.plot(
        b["y"],
        b["force_to"],
        color=ORANGE,
        ls="--",
        label=f"its force, {b['force_at_edge']:.2f} D/h where q₀sec² = 0.908",
    )
    ax.axvline(b["edge"], color=ORANGE, lw=0.9, ls=":")
    ax.text(b["edge"], 5.6, f" y = {b['edge']:.2f}h", fontsize=8, color=ORANGE)
    ax.set_xlim(0, 2.5)
    ax.set_ylim(0, 6)
    ax.set_xlabel("distance from the null over h")
    ax.set_ylabel("well over D; force over D/h")
    ax.legend(fontsize=7.4, loc="center right")
    ax.set_title("facing boards (FacingBoards.field):\nstiff toward a board")
    a = res["along"]
    ax = fig.add_subplot(grid[1, :])
    for name, colour in (("whole bend in the plane", INK), ("k = 0", BLUE), ("k = -1", ORANGE)):
        v = a[name]
        ax.semilogy(
            a["s"], v, color=colour, lw=1.6, label=f"{name.replace('-', '−')}: least {float(np.min(v)):.1f} eV per eV"
        )
    ax.semilogy(
        a["s"],
        a["k = 0, toward"],
        color=BLUE,
        lw=0.9,
        ls=":",
        label="k = 0, the bend toward the boards at twice the in-plane rule (hypothesis)",
    )
    ax.fill_between(
        a["s"],
        0.61 * a["k = 0"],
        0.95 * a["k = 0"],
        color=BLUES[1],
        alpha=0.7,
        lw=0,
        label="k = 0 discounted as four rods tracked, 0.61 to 0.95",
    )
    ax.set_ylim(1.5, 400)
    ax.set_xlabel("arc length along T915-S, fraction of a lap")
    ax.set_ylabel("energy held per eV of depth, 0.605·R/h")
    ax.legend(fontsize=7.8, loc="lower center", ncol=3)
    ax.set_title(
        "the energy held against curvature along the knot: T915-S's strips framed with k = 0 and k = −1, h = 0.915 mm"
    )
    tre = res["tracked"]["trefoil, q = 0.6"]
    ring = res["tracked"]["ring, q = 0.6"]
    _suptitle(
        fig,
        f"Energy held against curvature: four rods on the first run's trefoil held {tre['highest_held']:g} eV at q = "
        "0.6, "
        f"{tre['highest_held'] / tre['rule']:.2f} of the bend rule, where the ring of the same length held "
        f"{ring['highest_held']:g} eV; "
        f"facing boards' in-plane force peaks at {b['force_max_in']:.3f} D/h, and framed with k = 0 T915-S holds "
        f"{float(np.min(a['k = 0'])):.1f} eV per eV by that limit against "
        f"{float(np.min(a['whole bend in the plane'])):.1f} were its bend in the plane",
    )
    _save(fig, out, "held")


# ---- stopbands (figure 10) ----------------------------------------------------------------------------------------
STOP_CASES = (
    ("no joints (control)", 0, "none", 1e-3, 0.0, 0.0),
    ("six 1 mm joints, 15 % dip", 6, "regular", 1e-3, 0.15, 0.0),
    ("the same, placed irregularly", 6, "irregular", 1e-3, 0.15, 0.0),
    ("three 1 mm joints, 15 % dip", 3, "regular", 1e-3, 0.15, 0.0),
    ("six 4 mm joints, 15 % dip", 6, "regular", 4e-3, 0.15, 0.0),
    ("six 8 mm joints, 15 % dip", 6, "regular", 8e-3, 0.15, 0.0),
)


def compute_stopbands(quick: bool) -> dict[str, Any]:
    say("\n== stopbands: joints as resonances, on the unrolled guide and on the trefoil against the ring ==")
    out: dict[str, Any] = {}
    fs = scheme_secular(UNROLL_Q, UNROLL_F, UNROLL_SPP)
    out["fs_MHz"] = fs / 1e6
    out["fs_first_order_MHz"] = pt.secular_frequency(UNROLL_Q, UNROLL_F) / 1e6
    out["fs_mathieu_MHz"] = mathieu_beta(UNROLL_Q) * UNROLL_F / 2 / 1e6
    say(
        f"  secular frequency at q = 0.6, 1 GHz: the integrator's {out['fs_MHz']:.2f} MHz, Mathieu's "
        f"{out['fs_mathieu_MHz']:.2f}, first order {out['fs_first_order_MHz']:.1f}"
    )
    count, laps = (31, 6.0) if quick else (241, 40.0)
    energies = tuple(float(e) for e in np.linspace(3.9, 4.5, count))
    jobs = []
    seeds = []
    for k, (_, joints, spacing, joint, dip, step) in enumerate(STOP_CASES):
        seed = 3 if spacing == "irregular" else k
        seeds.append(seed)
        jobs.append((joints, spacing, joint, dip, step, energies, laps, seed))
    fine_count, fine_laps = (5, 1.5) if quick else (71, 20.0)
    fine = tuple(float(e) for e in np.linspace(3.93, 4.07, fine_count))
    chunks = max(1, min(len(fine), (WORKERS - 1) // 2))
    parts = [fine[i::chunks] for i in range(chunks)]
    began = time.time()
    rod_jobs = [
        (name, UNROLL_Q, part, fine_laps, 4e-3, 0.0, 1e-6, 0.0) for name in ("trefoil", "ring") for part in parts
    ]
    results = run_jobs(_any_job, [("unrolled", j) for j in jobs] + [("rods", j) for j in rod_jobs])
    say(f"  [{time.time() - began:.0f} s on {WORKERS} processes]")
    unrolled = results[: len(jobs)]
    rods = results[len(jobs) :]
    out["unrolled"] = {}
    for (name, joints, spacing, joint, dip, step), res, seed in zip(STOP_CASES, unrolled, seeds, strict=True):
        centres = unrolled_centres(joints, spacing, seed)
        bands = predicted_bands(centres, joint, dip, step, energies[0], energies[-1], fs) if joints else []
        lost = int(np.sum(~np.asarray(res["held"], dtype=bool)))
        out["unrolled"][name] = {**res, "bands": bands, "lost": lost}
        say(
            f"  {name}: {lost} of {len(energies)} lost in {laps:g} laps, median reach "
            f"{float(np.median(res['reach_um'])):.0f} um"
        )
    out["fine"] = {}
    for k, name in enumerate(("trefoil", "ring")):
        pieces = rods[k * chunks : (k + 1) * chunks]
        e = np.concatenate([p["energy"] for p in pieces])
        order = np.argsort(e)
        reach = np.concatenate([p["reach_um"] for p in pieces])[order]
        held = np.concatenate([p["held"] for p in pieces])[order]
        e = e[order]
        curve = rods_curve(name)
        tune = fs * curve.length / np.array([pt.v_of(float(x)) for x in e])
        base = float(np.median(reach))
        raised = reach > max(1.5 * base, base + 20)
        out["fine"][name] = {"energy": e, "reach_um": reach, "held": held, "tune": tune, "raised": raised, "base": base}
        groups = _runs(np.flatnonzero(raised))
        say(
            f"  fine scan on the {name} (six 4 mm joints, the null stepped 1 um): {len(groups)} raised band(s) at "
            "tunes " + ", ".join(f"{tune[g].min():.2f} to {tune[g].max():.2f}" for g in groups)
        )
    out["laps"] = laps
    out["fine_laps"] = fine_laps
    return out


def _runs(idx: IntArr) -> list[IntArr]:
    """Indices grouped into runs of neighbours (a gap of one allowed)."""
    if len(idx) == 0:
        return []
    return [g for g in np.split(idx, np.flatnonzero(np.diff(idx) > 2) + 1) if len(g)]


def _any_job(job: tuple[str, Any]) -> dict[str, Any]:
    """A tracking job of either kind, so one pool runs them all."""
    kind, args = job
    return unrolled_job(args) if kind == "unrolled" else rods_job(args)


def figure_stopbands(res: dict[str, Any], out: Path) -> None:
    """Figure 10: reach against energy for each joint case with the families' predicted centres, and the fine scan
    round tune 54 on the trefoil against the flat ring."""
    plt = _pyplot()
    fig = plt.figure(figsize=(15.5, 9.8))
    outer = fig.add_gridspec(1, 2, width_ratios=[1.6, 1.0], wspace=0.18)
    grid = outer[0].subgridspec(len(STOP_CASES), 1, hspace=0.12)
    right = outer[1].subgridspec(2, 1, hspace=0.32)
    family_colour = {
        "2f_s": BLUE,
        "f − 2f_s": ORANGE,
        "f + 2f_s": GREEN,
        "f_s": PURPLE,
        "f − f_s": PINK,
        "f + f_s": YELLOW,
    }
    shown: set[str] = set()
    for k, (name, row) in enumerate(res["unrolled"].items()):
        ax = fig.add_subplot(grid[k, 0])
        e = np.asarray(row["energy"])
        reach = np.asarray(row["reach_um"])
        held = np.asarray(row["held"], dtype=bool)
        ax.plot(e, reach, color=INK, lw=0.9)
        ax.plot(e[~held], reach[~held], "o", color=RED, ms=2.6)
        strongest = max((b["strength"] for b in row["bands"]), default=1.0)
        for b in row["bands"]:
            if b["strength"] > 0.08 * strongest:
                colour = family_colour[b["family"]]
                ax.axvline(b["energy"], color=colour, lw=1.0, alpha=0.8, ymax=0.25)
                shown.add(b["family"])
        ax.set_yscale("log")
        ax.set_ylim(5, 600)
        ax.set_xlim(e[0], e[-1])
        ax.text(
            0.005,
            0.86,
            f"{name}: {row['lost']} of {len(e)} lost",
            transform=ax.transAxes,
            fontsize=8.2,
            bbox={"facecolor": SURFACE, "edgecolor": "none", "alpha": 0.85, "pad": 1.5},
            zorder=10,
        )
        if k < len(STOP_CASES) - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("energy (eV)")
        if k == len(STOP_CASES) // 2:
            ax.set_ylabel("reach (µm; the wall at 400)")
    last = fig.axes[len(STOP_CASES) - 1]
    for fam in family_colour:
        if fam in shown:
            last.plot([], [], color=family_colour[fam], label=f"{fam}: predicted centre")
    last.plot([], [], "o", color=RED, ms=3, label="lost")
    last.legend(fontsize=7.4, loc="upper center", ncol=4, bbox_to_anchor=(0.5, -0.55))
    fig.axes[0].set_title(
        f"the unrolled first-run guide, q = 0.6 at 1 GHz: {len(res['unrolled']['no joints (control)']['energy'])} "
        f"energies, {res['laps']:g} laps",
        loc="left",
    )
    for k, name in enumerate(("trefoil", "ring")):
        ax = fig.add_subplot(right[k, 0])
        row = res["fine"][name]
        ax.plot(row["tune"], row["reach_um"], color=INK, lw=1.0, marker="o", ms=2.5)
        raised = np.asarray(row["raised"], dtype=bool)
        ax.plot(np.asarray(row["tune"])[raised], np.asarray(row["reach_um"])[raised], "o", color=ORANGE, ms=3.5)
        ax.axvline(54, color=MUTED, ls=":", lw=1.0)
        ax.invert_xaxis()
        ax.set_ylabel("reach (µm)")
        groups = _runs(np.flatnonzero(raised))
        which = "the first run's trefoil" if name == "trefoil" else "the flat ring of its length"
        ax.set_title(
            f"{which}: six 4 mm joints, the null stepped 1 µm; {len(groups)} band{'s' if len(groups) != 1 else ''}",
            loc="left",
            fontsize=9,
        )
        if k == 1:
            ax.set_xlabel("tune $\\nu$ = f_s L / v (energy rises to the left)")
    _suptitle(
        fig,
        "Joints as stopbands: each joint pattern is resonant where its lap harmonics meet 2f_s, f ∓ 2f_s, f_s or f ∓ "
        "f_s; six "
        f"1 mm dips of 15 % lose {res['unrolled']['six 1 mm joints, 15 % dip']['lost']} electrons of "
        f"{len(res['unrolled']['no joints (control)']['energy'])}, placed irregularly "
        f"{res['unrolled']['the same, placed irregularly']['lost']}, "
        f"three {res['unrolled']['three 1 mm joints, 15 % dip']['lost']}, 4 and 8 mm joints "
        f"{res['unrolled']['six 4 mm joints, 15 % dip']['lost']} and "
        f"{res['unrolled']['six 8 mm joints, 15 % dip']['lost']}; "
        "a pure null step on the trefoil against the ring",
    )
    _save(fig, out, "stopbands")


# ---- acceleration (figure 11) -------------------------------------------------------------------------------------
def pi_mode_ramp(
    length: float,
    va: float,
    phase_deg: float,
    w0: float,
    w1: float,
    count: int,
    seed: int,
    spread: float = 0.02,
    harmonic: int = 3,
    gaps: int = 6,
) -> dict[str, Any]:
    """The joints as gaps (section 6.4): `gaps` equal segments at alternating common-mode volts +-va make 2 va across
    each joint, synchronous at f_a = harmonic f_rev; a kick-drift map of `count` electrons spread in phase over an RF
    period and in energy by `spread`, the synchronous one gaining 2 va sin(phase) a gap, from w0 to w1 eV."""
    rng = np.random.default_rng(seed)
    phis = math.radians(phase_deg)
    w = w0 + spread * rng.standard_normal(count)
    trf = length / pt.v_of(w0) / harmonic
    dt = rng.uniform(-0.5, 0.5, count) * trf
    ws, time_s, laps = w0, 0.0, 0
    history = [w.copy()]
    sync = [ws]
    chirp = [harmonic * pt.v_of(ws) / length]
    while ws < w1:
        for _ in range(gaps):
            vs = pt.v_of(ws)
            fa = harmonic * vs / length
            w = np.maximum(w + 2 * va * np.sin(phis + 2 * math.pi * fa * dt), 1e-3)
            ws += 2 * va * math.sin(phis)
            vs2 = pt.v_of(ws)
            v = np.sqrt(2 * w * pt.E / pt.ME)
            dt = dt + (length / gaps) / v - (length / gaps) / vs2
            time_s += (length / gaps) / vs2
        laps += 1
        history.append(w.copy())
        sync.append(ws)
        chirp.append(harmonic * pt.v_of(ws) / length)
    trf = length / pt.v_of(ws) / harmonic
    captured = (np.abs(dt) < trf / 2) & (np.abs(w - ws) < 1.0)
    return {
        "laps": laps,
        "time_us": time_s * 1e6,
        "captured": float(captured.mean()),
        "history": np.array(history),
        "sync": np.array(sync),
        "chirp_MHz": np.array(chirp) / 1e6,
        "captured_mask": captured,
    }


def crossing_gain(
    eps: float, ramp: float, nu0: float = 54.0, w0: float = 4.01, gaps: int = 6, span: float = 0.05
) -> float:
    """Lambda, the log of the amplitude gain through a gradient band 2 nu = 6n crossed while the energy ramps `ramp`
    eV a lap: a lap map of J cells, each a turn by 2 pi nu / J and a thin lens of strength eps 2 pi nu / J, over +-span
    of the band's energy, which keeps the neighbouring bands at nu0 +- 3 out."""
    lo, hi = w0 * (1 - span), w0 * (1 + span)
    m = np.eye(2)
    for n in range(int((hi - lo) / ramp)):
        w = lo + n * ramp
        nu = nu0 * math.sqrt(w0 / w)
        mu = 2 * math.pi * nu / gaps
        c, s = math.cos(mu), math.sin(mu)
        cell = np.array([[1.0, 0.0], [-eps * 2 * math.pi * nu / gaps, 1.0]]) @ np.array([[c, s], [-s, c]])
        for _ in range(gaps):
            m = cell @ m
    sv = np.linalg.svd(m, compute_uv=False)
    return math.log(float(sv[0]) / math.sqrt(abs(float(np.linalg.det(m)))))


def compute_acceleration(quick: bool) -> dict[str, Any]:
    say("\n== acceleration: the joints as gaps, a chirp, and the bands crossed on the way ==")
    out: dict[str, Any] = {}
    count = 1000 if quick else 4000
    beta = mathieu_beta(0.5)
    fs = beta * F915 / 2
    for name, length, seed in (("T915-S6", L915, 1), ("the first run's guide", first_run_trefoil().length, 2)):
        r = pi_mode_ramp(length, 0.010, 30.0, 1.0, 8.0, count, seed)
        r["length"] = length
        if name == "T915-S6":
            w = r["sync"]
            r["tune"] = fs * length / np.array([pt.v_of(float(x)) for x in w])
        out[name] = r
        say(
            f"  {name} ({length / MM:.0f} mm): 10 mV gaps at 30 deg, 60 meV a lap: 1 -> {r['sync'][-1]:.2f} eV in "
            f"{r['laps']} laps, "
            f"{r['time_us']:.1f} us, chirp 3 f_rev {r['chirp_MHz'][0]:.2f} -> {r['chirp_MHz'][-1]:.2f} MHz, captured "
            f"{r['captured'] * 100:.0f} %"
        )
    t = out["T915-S6"]
    say(
        f"  T915-S6 at q = 0.5 (exact Mathieu f_s {fs / 1e6:.1f} MHz): the tune falls from {t['tune'][0]:.1f} to "
        f"{t['tune'][-1]:.1f}"
    )
    rates = np.geomspace(0.005, 0.1, 10 if quick else 24)  # a band takes several laps to cross
    gains = {}
    for eps in (0.002, 0.005, 0.01):
        gains[eps] = {
            "map": np.array([crossing_gain(eps, float(r)) for r in rates]),
            "formula": np.arccosh(np.exp(math.pi**2 * 54 * eps**2 * 4.01 / (2 * rates))),
        }
    out["gain"] = {"rates": rates, "curves": gains}
    # a 15 % dip over 1 mm at six joints on the first run's guide: eps_k = 2 g_k at k = 108
    centres = unrolled_centres(6, "regular", 0)
    n = 1 << 15
    s = np.arange(n) * UNROLL_L / n
    w, _ = joint_profile(centres, 1e-3, UNROLL_L, s)
    gk = abs(np.fft.fft(0.15 * w)[108]) / n
    eps_design = 2 * gk
    g_design = math.pi**2 * 54 * eps_design**2 * 4.01 / (2 * 0.06)
    out["design"] = {"eps": eps_design, "G": g_design, "Lambda": float(np.arccosh(math.exp(g_design)))}
    say(
        f"  six 1 mm dips of 15 %: eps_k = 2 g_k = {eps_design:.4f} at k = 108, G = pi^2 nu eps^2 W / (2 dW) = "
        f"{g_design:.2f} a crossing at 60 meV a lap"
    )
    return out


def figure_acceleration(res: dict[str, Any], out: Path) -> None:
    """Figure 11: energy against lap for the six-gap pi-mode ramp with its crossings marked, and the crossing gain
    against the ramp rate."""
    plt = _pyplot()
    fig = plt.figure(figsize=(15.5, 6.8))
    grid = fig.add_gridspec(1, 2, width_ratios=[1.5, 1.0], wspace=0.22, top=0.8)
    r = res["T915-S6"]
    ax = fig.add_subplot(grid[0, 0])
    hist = r["history"]
    laps = np.arange(len(hist))
    cap = np.asarray(r["captured_mask"], dtype=bool)
    pick_c = np.flatnonzero(cap)[:30]
    pick_l = np.flatnonzero(~cap)[:12]
    for i in pick_l:
        ax.plot(laps, hist[:, i], color=MUTED, lw=0.6, alpha=0.7)
    for i in pick_c:
        ax.plot(laps, hist[:, i], color=BLUE, lw=0.6, alpha=0.7)
    ax.plot(laps, r["sync"], color=INK, lw=2.0, label="synchronous electron")
    tune = r["tune"]
    for n in range(1, 40):
        for target, colour, ls in ((6 * n, ORANGE, "-"), (3 * n, GREEN, ":")):
            hit = np.flatnonzero((tune[:-1] - target) * (tune[1:] - target) <= 0)
            for h in hit:
                ax.axvline(h, color=colour, lw=0.9, ls=ls, alpha=0.8)
    ax.plot([], [], color=BLUE, lw=0.8, label=f"captured ({r['captured'] * 100:.0f} % of {hist.shape[1]})")
    ax.plot([], [], color=MUTED, lw=0.8, label="left behind")
    ax.plot([], [], color=ORANGE, lw=0.9, label="tune crosses $\\nu$ = 6n: null-step bands of the six cuts")
    ax.plot([], [], color=GREEN, lw=0.9, ls=":", label="$\\nu$ = 3n: their gradient bands, 2$\\nu$ = 6n")
    ax.set_xlabel("lap")
    ax.set_ylabel("kinetic energy (eV)")
    ax.set_ylim(0, 9)
    ax.legend(fontsize=7.8, loc="upper left")
    ax2 = ax.twiny()
    ticks = [0, 20, 40, 60, 80, 100, 116]
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(ticks)
    ax2.set_xticklabels([f"{r['chirp_MHz'][min(t, len(r['chirp_MHz']) - 1)]:.1f}" for t in ticks], fontsize=7.6)
    ax2.set_xlabel("the chirp 3 f_rev (MHz)", fontsize=8.5)
    ax2.grid(False)
    ax.set_title(f"T915-S6: six π-mode gaps at 10 mV, 30°; tune {tune[0]:.0f} → {tune[-1]:.0f}", pad=34)
    ax = fig.add_subplot(grid[0, 1])
    g = res["gain"]
    for colour, (eps, c) in zip((BLUE, ORANGE, GREEN), g["curves"].items(), strict=True):
        ax.loglog(g["rates"] * 1e3, c["map"], "o", color=colour, ms=4)
        ax.loglog(
            g["rates"] * 1e3,
            c["formula"],
            color=colour,
            lw=1.3,
            label=f"ε_k = {eps:g}: map (dots), cosh Λ = e^G (line)",
        )
    d = res["design"]
    ax.plot([60], [d["Lambda"]], "*", color=INK, ms=12)
    ax.annotate(
        f"six 1 mm dips of 15 %: ε_k = {d['eps']:.4f},\nG = {d['G']:.2f} a crossing at 60 meV a lap",
        (60, d["Lambda"]),
        xytext=(-150, 25),
        textcoords="offset points",
        fontsize=8,
        arrowprops={"arrowstyle": "-", "color": INK2, "lw": 0.7},
    )
    ax.set_xlabel("ramp (meV a lap)")
    ax.set_ylabel("Λ: log of the amplitude gain through one band")
    ax.legend(fontsize=7.8, loc="lower left")
    ax.set_title("crossing a gradient band ($\\nu$ = 54) while the energy ramps")
    first = res["the first run's guide"]
    _suptitle(
        fig,
        f"Acceleration by the joints: one chirp of 3 f_rev from {r['chirp_MHz'][0]:.2f} to {r['chirp_MHz'][-1]:.2f} "
        "MHz takes electrons "
        f"from 1 to 8 eV in {r['laps']} laps ({r['time_us']:.1f} µs) on T915-S6 and captures {r['captured'] * 100:.0f} "
        "% "
        f"({first['captured'] * 100:.0f} % on the first run's guide, {first['chirp_MHz'][0]:.2f} to "
        f"{first['chirp_MHz'][-1]:.2f} MHz); "
        "each band crossed is a squeeze cosh Λ = e^G",
    )
    _save(fig, out, "acceleration")


# ---- plotting -----------------------------------------------------------------------------------------------------
def _pyplot() -> Any:
    """matplotlib's pyplot on the Agg backend, set to the house style: one sans face, recessive grid and axes, the
    categorical colours in their fixed order."""
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
    """The sequential blue ramp, light to dark (or dark to light)."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("ramp", list(BLUES[::-1] if reverse else BLUES))


def _log_ticks(axis: Any, ticks: Sequence[float]) -> None:
    """Fixed ticks with plain labels on a logarithmic axis, and none on its minor ticks."""
    from matplotlib.ticker import NullFormatter

    axis.set_ticks(list(ticks))
    axis.set_ticklabels([f"{t:g}" for t in ticks])
    axis.set_minor_formatter(NullFormatter())


def _diverging() -> Any:
    """Blue through the surface to orange, for signed quantities."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("signed", [BLUES[6], BLUES[3], SURFACE, "#f3a77f", ORANGE, "#a8401a"])


def _suptitle(fig: Any, text: str, width: int = 150) -> None:
    """A figure's title, made from its data, wrapped so that it never sets the figure's width."""
    import re
    import textwrap

    glued = re.sub(r"(\d) (?=\S)", "\\1\u00a0", text)  # a number keeps its unit on its line
    fig.suptitle(textwrap.fill(glued, width), x=0.01, ha="left", fontsize=10.5)


def _save(fig: Any, out: Path, name: str) -> None:
    import matplotlib.pyplot as plt

    fig.savefig(out / f"{name}.png", bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    say(f"  wrote {name}.png")


# ---- the run ------------------------------------------------------------------------------------------------------
SECTIONS: dict[str, tuple[Callable[[bool], dict[str, Any]], Callable[[dict[str, Any], Path], None]]] = {
    "arrangements": (compute_arrangements, figure_arrangements),
    "knots": (compute_knots, figure_knots),
    "strips": (compute_strips, figure_strips),
    "ramp": (compute_ramp, figure_ramp),
    "joints": (compute_joints, figure_joints),
    "drive": (compute_drive, figure_drive),
    "laminate": (compute_laminate, figure_laminate),
    "held": (compute_held, figure_held),
    "stopbands": (compute_stopbands, figure_stopbands),
    "acceleration": (compute_acceleration, figure_acceleration),
}
ORDER = tuple(SECTIONS)
KEEP = 64  # results.json keeps arrays and lists of at most this many numbers; the figures hold the rest


def slim(value: Any) -> Any:
    """The results without the large arrays behind the figures and without the objects behind them (keys starting
    with an underscore), for results.json."""
    if isinstance(value, np.ndarray):
        return value if value.size <= KEEP else None
    if isinstance(value, dict):
        kept = {str(k): slim(v) for k, v in value.items() if not str(k).startswith("_")}
        return {k: v for k, v in kept.items() if v is not None}
    if isinstance(value, (list, tuple)):
        if len(value) > KEEP:
            return None
        return [slim(v) for v in value]
    if isinstance(value, (float, int, str, bool, np.generic)) or value is None:
        return value
    return None


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument("--quick", action="store_true", help="coarse grids and few draws, into knotted_guide/quick/")
    parser.add_argument("--only", default="", help="comma-separated sections, of: " + ", ".join(ORDER))
    parser.add_argument("--out", type=Path, default=None, help="where the figures, results.json and run.log go")
    parser.add_argument("--workers", type=int, default=None, help="processes for the tracking runs (quick: 1)")
    args = parser.parse_args(argv)
    global WORKERS
    WORKERS = args.workers if args.workers is not None else (1 if args.quick else max((os.cpu_count() or 2) - 1, 1))
    chosen = [name for name in args.only.split(",") if name] or list(ORDER)
    unknown = [name for name in chosen if name not in ORDER]
    if unknown:
        print(f"no section {unknown}; there are {list(ORDER)}", file=sys.stderr)
        return 2
    out: Path = args.out if args.out is not None else (OUT / "quick" if args.quick else OUT)
    out.mkdir(parents=True, exist_ok=True)
    LOG.clear()
    began = time.time()
    say(f"knotted guide study: {'quick' if args.quick else 'full'} run of {', '.join(chosen)}, {WORKERS} processes")
    results: dict[str, Any] = {}
    for name in chosen:
        start = time.time()
        results[name] = SECTIONS[name][0](args.quick)
        say(f"  [{name}: {time.time() - start:.1f} s]")
    say("\n== figures ==")
    for name in chosen:
        SECTIONS[name][1](results[name], out)
    say(f"\ndone in {time.time() - began:.0f} s")
    if chosen == list(ORDER):
        (out / "results.json").write_text(json.dumps(lab.plain(slim(results)), indent=1) + "\n")
        (out / "run.log").write_text("\n".join(LOG) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
