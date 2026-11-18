"""Charged particles in prescribed fields: the kernel under experiments/accelerator_study.py.

Objects (all immutable):
  ClosedCurve    a closed space curve: tangents, twist-free frames, holonomy, writhe
  QuadRods       four thin rods along a curve, alternating in sign (line charges)
  Coil           a filament winding round a curve, and its Biot-Savart field
  GuideOptics    first-order theory of a guide field that follows a curve: lap map, closed orbit, excursion
  Launch, GuideTrack    what goes into track_guide and what comes out
  OrbitFit       the circle through an orbit's offsets at the start of successive laps

Functions (pure):
  push                         one relativistic Boris step
  track_guide                  electrons round a closed guide in any static magnetic field
  orbit_circle                 the closed orbit and the turn per lap, measured from a tracked orbit
  uniform, scaled, superpose   fields built from fields
  v_of, rigidity, energy_of    kinematics

Models: point charges, thin wires, filament coils, no image charges, no collisions, no space charge.

Conventions: SI units. T is a curve's unit tangent, in the direction of increasing index. A guide field points
along +T.
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Callable
from dataclasses import dataclass
from functools import cached_property, reduce
from typing import Literal, NamedTuple

import numpy as np
from numpy.typing import NDArray

Arr = NDArray[np.float64]
CArr = NDArray[np.complex128]
IntArr = NDArray[np.intp]
BoolArr = NDArray[np.bool_]
Winding = Literal["single", "return"]
BField = Callable[[Arr], Arr]  # positions (N, 3) in metres -> magnetic field (N, 3) in tesla

E = 1.602176634e-19
ME = 9.1093837015e-31
C0 = 299792458.0
EPS0 = 8.8541878128e-12
MU0 = 4e-7 * math.pi
KB = 1.380649e-23
Z0 = MU0 * C0
QM_E = -E / ME


# ---- kinematics ---------------------------------------------------------------------------------------------
def v_of(w_ev: float, m: float = ME) -> float:
    """Speed (m/s) of a particle of mass m (kg) with kinetic energy w_ev (eV)."""
    g = 1 + w_ev * E / (m * C0**2)
    return C0 * math.sqrt(1 - 1 / g**2)


def rigidity(w_ev: float, m: float = ME) -> float:
    """Momentum per elementary charge (T m): the field-times-radius that bends the particle."""
    v = v_of(w_ev, m)
    return m * v / math.sqrt(1 - (v / C0) ** 2) / E


def energy_of(b_rho: float, m: float = ME) -> float:
    """Kinetic energy (eV) of a singly charged particle of mass m with rigidity b_rho (T m). Inverse of rigidity."""
    rest = m * C0**2
    return (math.hypot(E * b_rho * C0, rest) - rest) / E


# ---- fields as functions --------------------------------------------------------------------------------------
def uniform(b: Arr) -> BField:
    """The same field b (3,) everywhere."""
    return lambda x: np.broadcast_to(b, x.shape)


def scaled(field: BField, factor: float) -> BField:
    """A field multiplied by a number: a coil's field per ampere times its current."""
    return lambda x: field(x) * factor


def superpose(*fields: BField) -> BField:
    """The sum of several fields."""
    return lambda x: reduce(lambda a, b: a + b, (f(x) for f in fields))


# ---- integrator -----------------------------------------------------------------------------------------------
def cross(a: Arr, b: Arr) -> Arr:
    """Cross product along the last axis."""
    return np.stack(
        [
            a[..., 1] * b[..., 2] - a[..., 2] * b[..., 1],
            a[..., 2] * b[..., 0] - a[..., 0] * b[..., 2],
            a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0],
        ],
        axis=-1,
    )


def push(
    x: Arr, u: Arr, qm: Arr, dt: float, ef: Arr, bf: Arr, drag: float = 0.0, acc: Arr | None = None
) -> tuple[Arr, Arr]:
    """One relativistic Boris step for N particles: positions x (N, 3), u = gamma v (N, 3), charge/mass qm (N,).

    The magnetic rotation uses the tangent of the exact half angle, so the gyration frequency has no time-step
    error. The price is in crossed fields: the E x B drift comes out low by (theta/2) / tan(theta/2), theta being
    the gyro-angle per step (0.2 % at 40 steps per turn). drag is a Stokes rate (1/s); acc is a uniform
    acceleration (3,) in m/s^2.
    """
    k = (qm * dt / 2)[:, None]
    um = u + k * ef
    if acc is not None:
        um = um + acc * (dt / 2)
    gam = np.sqrt(1 + np.sum(um * um, 1) / C0**2)
    bn = np.sqrt(np.sum(bf * bf, 1))
    tv = bf * (np.tan(qm * dt / 2 * bn / gam) / np.where(bn > 0, bn, 1.0))[:, None]
    up = um + cross(um, tv)
    u = um + cross(up, tv * (2 / (1 + np.sum(tv * tv, 1)))[:, None]) + k * ef
    if acc is not None:
        u = u + acc * (dt / 2)
    if drag:
        u = u * math.exp(-drag * dt)
    gam = np.sqrt(1 + np.sum(u * u, 1) / C0**2)
    return x + u / gam[:, None] * dt, u


# ---- closed curves ----------------------------------------------------------------------------------------------
def _carry(v: Arr, a: Arr, b: Arr) -> Arr:
    """The unit vector v, perpendicular to the unit vector a, after the rotation that takes a to b about their
    common perpendicular: v - (a + b) (v . b) / (1 + a . b)."""
    w: Arr = v - (a + b) * (np.dot(v, b) / (1 + np.dot(a, b)))
    return w / np.linalg.norm(w)


class Frame(NamedTuple):
    """Unit normals (a, b) along a curve, with b = T x a, and the number k of close_mod turns added to close them."""

    a: Arr
    b: Arr
    k: int


@dataclass(frozen=True, eq=False)
class ClosedCurve:
    """A closed space curve: n points equally spaced in arc length (the point after the last is the first)."""

    points: Arr
    length: float

    @classmethod
    def from_polygon(cls, pts: Arr, n: int) -> ClosedCurve:
        """n points equally spaced in arc length along the closed polygon pts."""
        d = np.linalg.norm(np.roll(pts, -1, 0) - pts, axis=1)
        s = np.concatenate([[0], np.cumsum(d)])
        cc = np.vstack([pts, pts[:1]])
        si = np.arange(n) * s[-1] / n
        return cls(np.stack([np.interp(si, s, cc[:, k]) for k in range(3)], 1), float(s[-1]))

    @classmethod
    def torus_knot(cls, R: float, r: float, n: int, p: int = 2, q: int = 3) -> ClosedCurve:
        """A (p, q) torus knot on a torus of radii R and r. (2, 3) is the trefoil."""
        t = np.linspace(0, 2 * np.pi, 40 * n, endpoint=False)
        pts = np.stack(
            [(R + r * np.cos(q * t)) * np.cos(p * t), (R + r * np.cos(q * t)) * np.sin(p * t), r * np.sin(q * t)], 1
        )
        return cls.from_polygon(pts, n)

    @classmethod
    def figure_eight(cls, a: float, b: float, lift: float, n: int) -> ClosedCurve:
        """A figure of eight, 2a long and 2b wide, whose two strands pass `2 lift` apart where they cross: an
        unknotted curve that is not planar. With lift = 0 the strands meet."""
        t = np.linspace(0, 2 * np.pi, 40 * n, endpoint=False)
        return cls.from_polygon(np.stack([a * np.sin(t), b * np.sin(2 * t), lift * np.cos(t)], 1), n)

    @classmethod
    def circle(cls, R: float, n: int) -> ClosedCurve:
        """A circle of radius R in the xy plane, anticlockwise seen from +z."""
        t = np.linspace(0, 2 * np.pi, n, endpoint=False)
        return cls(np.stack([R * np.cos(t), R * np.sin(t), 0 * t], 1), 2 * math.pi * R)

    def __len__(self) -> int:
        return len(self.points)

    @cached_property
    def tangents(self) -> Arr:
        """Unit tangents (n, 3), by central differences."""
        t = np.roll(self.points, -1, 0) - np.roll(self.points, 1, 0)
        t /= np.linalg.norm(t, axis=1)[:, None]
        return t

    @cached_property
    def _carried(self) -> tuple[Arr, float]:
        """A unit normal carried round the curve without twisting (n, 3), and the angle (rad) by which it has
        turned about T when it returns to the start (positive: right-handed about T).

        Each step turns the normal by the rotation that takes one tangent to the next about their common
        perpendicular. The angle then converges as 1/n^2. (Projecting the normal onto the next normal plane, the
        other common rule, turns it by a spurious (step angle)^2 / 4 x sin(2 x its angle to the bend) per step and
        converges only as 1/n: 0.008 of a turn out at 600 points on some torus knots.)"""
        n, T = len(self), self.tangents
        n1 = np.empty_like(self.points)
        v = np.cross(T[0], [0.0, 0.0, 1.0])
        if np.linalg.norm(v) < 1e-6:
            v = np.cross(T[0], [1.0, 0.0, 0.0])
        n1[0] = v / np.linalg.norm(v)
        for i in range(1, n):
            n1[i] = _carry(n1[i - 1], T[i - 1], T[i])
        v = _carry(n1[-1], T[-1], T[0])
        return n1, math.atan2(np.dot(v, np.cross(T[0], n1[0])), np.dot(v, n1[0]))

    @cached_property
    def axes(self) -> tuple[Arr, Arr]:
        """Unit normals (e1, e2 = T x e1) at every point, each (n, 3), carried along the curve without twisting.
        They do not close: after a lap they have turned by the holonomy."""
        e1 = self._carried[0]
        return e1, np.cross(self.tangents, e1)

    @property
    def holonomy(self) -> float:
        """Angle (rad, in (-pi, pi]) by which a normal carried once round the curve without twisting has turned.

        By the Calugareanu-White-Fuller theorem it is 2 pi times the writhe, modulo 2 pi; zero for a planar curve."""
        return self._carried[1]

    @cached_property
    def writhe(self) -> float:
        """The Gauss double integral of the curve with itself, by the midpoint rule on the polygon's segments."""
        d = np.roll(self.points, -1, 0) - self.points
        m = self.points + d / 2
        r = m[:, None, :] - m[None, :, :]
        dist = np.sqrt(np.sum(r * r, axis=2))
        np.fill_diagonal(dist, np.inf)
        trip = np.sum(np.cross(d[:, None, :], d[None, :, :]) * r, axis=2)
        return float(np.sum(trip / dist**3) / (4 * math.pi))

    @cached_property
    def bend_radius(self) -> float:
        """Smallest radius of curvature (m)."""
        c, ds = self.points, self.length / len(self)
        k = np.linalg.norm(np.roll(c, -1, 0) - 2 * c + np.roll(c, 1, 0), axis=1) / ds**2
        return float(1 / k.max())

    def frame(self, close_mod: float) -> Frame:
        """Twist-free normals with a uniform extra twist, so that after one lap they return turned by a multiple
        k of close_mod (pi: a quadrupole's like-signed rods meet; 2 pi: a winding closes)."""
        n, T = len(self), self.tangents
        n1, alpha = self._carried
        k = round(alpha / close_mod)
        tw = -(alpha - k * close_mod) * np.arange(n) / n
        n2 = np.cross(T, n1)
        ct, st = np.cos(tw)[:, None], np.sin(tw)[:, None]
        return Frame(ct * n1 + st * n2, -st * n1 + ct * n2, k)

    def transverse(self, z: CArr) -> Arr:
        """Offsets z = a + ib at the first point, on the twist-free axes (e1, e2 = T x e1), as vectors (K, 3)."""
        e1, e2 = self.axes
        vectors: Arr = np.outer(z.real, e1[0]) + np.outer(z.imag, e2[0])
        return vectors

    def advance(self, x: Arr, st: IntArr, adv: Arr, iters: int = 1) -> tuple[IntArr, Arr]:
        """Move each particle's station index st to a nearer neighbouring point, counting the steps in adv."""
        c, n = self.points, len(self)
        for _ in range(iters):
            d0 = ((x - c[st]) ** 2).sum(1)
            ip, im = (st + 1) % n, (st - 1) % n
            fw = ((x - c[ip]) ** 2).sum(1) < d0
            bw = (~fw) & (((x - c[im]) ** 2).sum(1) < d0)
            st = np.where(fw, ip, np.where(bw, im, st))
            adv = adv + fw - bw
        return st, adv

    def offset(self, x: Arr, st: IntArr) -> tuple[Arr, Arr]:
        """Each particle's displacement from the curve at its station, across the tangent: vectors and lengths."""
        T = self.tangents
        d = x - self.points[st]
        d = d - (d * T[st]).sum(1)[:, None] * T[st]
        return d, np.sqrt((d * d).sum(1))


# ---- electrodes -------------------------------------------------------------------------------------------------
@dataclass(frozen=True, eq=False)
class QuadRods:
    """Four thin rods at r0 round a curve, alternating in sign: segment midpoints (4, n, 3), lengths (4, n), signs."""

    curve: ClosedCurve
    r0: float
    mid: Arr
    dl: Arr
    sign: Arr
    frame: Frame

    @classmethod
    def along(cls, curve: ClosedCurve, r0: float) -> QuadRods:
        fr = curve.frame(math.pi)
        W = np.array(
            [
                curve.points + r0 * (math.cos(j * math.pi / 2) * fr.a + math.sin(j * math.pi / 2) * fr.b)
                for j in range(4)
            ]
        )
        nxt = np.empty_like(W)
        for j in range(4):
            nxt[j, :-1] = W[j, 1:]
            nxt[j, -1] = W[(j + 2 * fr.k) % 4, 0]  # after a lap a rod meets the like-signed rod opposite
        return cls(curve, r0, (W + nxt) / 2, np.linalg.norm(nxt - W, axis=2), np.array([1.0, -1.0, 1.0, -1.0]), fr)

    def field(self, x: Arr, st: IntArr, lam: float | Arr, offs: IntArr) -> Arr:
        """Electric field (N, 3) of the four line charges (sign * lam C/m) from the segments at stations st + offs."""
        idx = (st[:, None] + offs[None, :]) % self.mid.shape[1]
        out = np.zeros_like(x)
        for j in range(4):
            r = x[:, None, :] - self.mid[j][idx]
            d2 = np.sum(r * r, axis=2)
            out += self.sign[j] * np.sum(r * (self.dl[j][idx] / (d2 * np.sqrt(d2)))[:, :, None], axis=1)
        scale: Arr = (np.asarray(lam, dtype=np.float64) / (4 * math.pi * EPS0)).reshape(-1, 1)
        return out * scale


# ---- coils ------------------------------------------------------------------------------------------------------
def _helix(curve: ClosedCurve, radius: float, turns: int, per_turn: int, sense: int) -> tuple[Arr, Arr]:
    """One layer of `turns` turns round the curve: segment midpoints and dl vectors. The wire advances along +T;
    sense = +1 circulates right-handed about T, -1 left-handed."""
    c, fr = curve.points, curve.frame(2 * math.pi)
    n, m = len(curve), turns * per_turn
    f = np.arange(m) / m * n
    i0 = np.floor(f).astype(int)
    w = (f - i0)[:, None]
    i1 = (i0 + 1) % n
    ax, aa, bb = c[i0] * (1 - w) + c[i1] * w, fr.a[i0] * (1 - w) + fr.a[i1] * w, fr.b[i0] * (1 - w) + fr.b[i1] * w
    ph = sense * (2 * math.pi * turns * np.arange(m) / m)
    P = ax + radius * (np.cos(ph)[:, None] * aa + np.sin(ph)[:, None] * bb)
    nxt = np.roll(P, -1, 0)
    return (P + nxt) / 2, nxt - P


@dataclass(frozen=True, eq=False)
class Coil:
    """A filament winding round a closed curve, carrying one ampere, oriented so that its field points along +T."""

    curve: ClosedCurve
    radius: float
    turns: int  # in all layers together
    winding: Winding
    mid: Arr  # segment midpoints (M, 3)
    dl: Arr  # segment vectors (M, 3), along the current

    @classmethod
    def wound(
        cls,
        curve: ClosedCurve,
        radius: float,
        turns: int,
        winding: Winding = "single",
        per_turn: int = 36,
        layer_gap: float = 1e-3,
    ) -> Coil:
        """Wind `turns` turns of the given radius per layer.

        "single": one layer. The wire advances once round the curve, so the winding also carries its current once
        along the curve: one net turn through the hole.
        "return": a second layer, layer_gap outside the first, wound back to the start with the same circulation.
        The two advances cancel: twice the turns, and no net current along the curve.
        """
        mid, dl = _helix(curve, radius, turns, per_turn, 1)
        if winding == "return":
            mid2, dl2 = _helix(curve, radius + layer_gap, turns, per_turn, -1)
            mid, dl = np.vstack([mid, mid2]), np.vstack([dl, -dl2])
        made = cls(curve, radius, turns * (2 if winding == "return" else 1), winding, mid, dl)
        along = float(np.mean(np.sum(made.axis_field * curve.tangents, 1)))
        return made if along >= 0 else cls(curve, radius, made.turns, winding, mid, -dl)

    def field(self, x: Arr) -> Arr:
        """Biot-Savart field (N, 3), tesla per ampere, each segment taken at its midpoint."""
        mid, dl = self.mid, self.dl
        r = x[:, None, :] - mid[None, :, :]
        d2 = np.sum(r * r, axis=2)
        k3 = 1 / (d2 * np.sqrt(d2))
        return (
            MU0
            / (4 * math.pi)
            * np.stack(
                [
                    np.sum((dl[None, :, 1] * r[:, :, 2] - dl[None, :, 2] * r[:, :, 1]) * k3, 1),
                    np.sum((dl[None, :, 2] * r[:, :, 0] - dl[None, :, 0] * r[:, :, 2]) * k3, 1),
                    np.sum((dl[None, :, 0] * r[:, :, 1] - dl[None, :, 1] * r[:, :, 0]) * k3, 1),
                ],
                1,
            )
        )

    @cached_property
    def axis_field(self) -> Arr:
        """The field on the curve itself (n, 3), tesla per ampere."""
        c = self.curve.points
        return np.vstack([self.field(c[i : i + 40]) for i in range(0, len(c), 40)])

    @cached_property
    def per_amp(self) -> float:
        """Mean field strength on the curve, tesla per ampere."""
        return float(np.linalg.norm(self.axis_field, axis=1).mean())

    @cached_property
    def slope(self) -> Arr:
        """The field across the curve divided by the field along it, on the curve (n, 3): the slope at which the
        field line through each point leaves the curve."""
        T = self.curve.tangents
        along = np.sum(self.axis_field * T, 1)[:, None]
        return (self.axis_field - along * T) / along

    @cached_property
    def optics(self) -> GuideOptics:
        """First-order optics of this coil's guide field, including the wander of its field lines."""
        return GuideOptics.of(self.curve, self.slope)

    def linked_current(self, loop_radius: float, points: int = 720) -> float:
        """Ampere's law on a circle of loop_radius round the curve at its first point: the current the winding
        carries along the curve, in amperes per ampere in the wire (its net turns through the hole)."""
        fr, c0 = self.curve.frame(2 * math.pi), self.curve.points[0]
        phi = (np.arange(points) + 0.5) * 2 * math.pi / points
        cos, sin = np.cos(phi)[:, None], np.sin(phi)[:, None]
        pts = c0 + loop_radius * (cos * fr.a[0] + sin * fr.b[0])
        dl = loop_radius * (2 * math.pi / points) * (-sin * fr.a[0] + cos * fr.b[0])
        return float(np.sum(self.field(pts) * dl) / MU0)


# ---- first-order theory of a curved guide field -----------------------------------------------------------------
@dataclass(frozen=True, eq=False)
class GuideOptics:
    """A guide field that follows a closed curve, to first order in gyroradius / bend radius.

    A guiding centre's offset from the curve is written z = a + ib on axes (e1, e2 = T x e1) carried along the
    curve without twisting. Within one lap z(s) = z(0) + rho g(s) + h(s). On returning, the axes have turned by
    alpha, so the lap map is z -> exp(i alpha) (z + rho g[-1] + h[-1]).

    rho = p / (e B) is signed: positive for an electron moving along the field, negative against it. An electron
    moving against the field goes round the other way, so its lap is the inverse of this map at its (negative)
    rho. The inverse has the same closed orbit and turns about it on the same circle, so closed_orbit and
    matched_excursion hold for both directions as they stand; excursion follows the forward map.

    The laws of the lap map are GDE-1 to GDE-7 in docs/LAWS.md, proved in formal/Array16/Guide.lean. That a guiding
    centre obeys it is a hypothesis (GDE-H1), tested by tracking.
    """

    alpha: float  # turn of the transverse plane per lap (rad)
    g: CArr  # (n + 1,) accumulated curvature drift per unit rho
    h: CArr  # (n + 1,) accumulated wander of the field line through the curve (m)

    @classmethod
    def of(cls, curve: ClosedCurve, slope: Arr | None = None) -> GuideOptics:
        """Optics of a field along `curve` that is circularly symmetric about it, so that offsets are carried
        without twisting. The curvature drift of an electron is d(xi)/ds = -rho T x dT/ds. slope (n, 3), if given,
        is the field across the curve over the field along it: Coil.slope."""
        T, (e1, e2) = curve.tangents, curve.axes
        kick = np.cross(T, np.roll(T, -1, 0) - T)  # curvature x binormal x ds
        zero = np.zeros(1, dtype=np.complex128)
        g: CArr = np.concatenate([zero, np.cumsum(-(np.sum(kick * e1, 1) + 1j * np.sum(kick * e2, 1)))])
        h: CArr = np.zeros(len(curve) + 1, dtype=np.complex128)
        if slope is not None:
            ds = np.linalg.norm(np.roll(curve.points, -1, 0) - curve.points, axis=1)
            h = np.concatenate([zero, np.cumsum((np.sum(slope * e1, 1) + 1j * np.sum(slope * e2, 1)) * ds)])
        return cls(curve.holonomy, g, h)

    @property
    def transform(self) -> float:
        """Rotational transform: turns of the transverse plane per lap."""
        return self.alpha / (2 * math.pi)

    @property
    def closes(self) -> bool:
        """Whether the transverse plane turns at all. If it does not, every orbit drifts by the same step each lap."""
        return abs(1 - cmath.exp(1j * self.alpha)) > 1e-9

    def step(self, rho: Arr) -> CArr:
        """What one lap adds to the offset before the plane turns: curvature drift plus field-line wander."""
        added: CArr = rho * self.g[-1] + self.h[-1]
        return added

    def lap(self, rho: Arr, z: CArr) -> CArr:
        """The offsets one lap later: exp(i alpha) (z + step)."""
        after: CArr = cmath.exp(1j * self.alpha) * (z + self.step(rho))
        return after

    def closed_orbit(self, rho: Arr) -> CArr:
        """Offset at the first point of the orbit that closes on itself after one lap, for each rho; nan if the
        plane does not turn. GDE-1: it is a fixed point of the lap map. GDE-2: it is the only one."""
        if not self.closes:
            return np.full(rho.shape, complex(math.nan, math.nan), dtype=np.complex128)
        rot = cmath.exp(1j * self.alpha)
        return self.step(rho) * (rot / (1 - rot))

    def excursion(self, rho: Arr, z0: complex | CArr, laps: int) -> Arr:
        """Largest distance from the curve (m) reached in `laps` laps by guiding centres launched at offset z0.

        GDE-4, GDE-5: if the plane turns, the offsets at the start of each lap lie on a circle about the closed
        orbit, so this stays bounded however many laps. GDE-6: if it does not turn, it grows without limit.
        GDE-7: it is convex in rho."""
        z: CArr = np.broadcast_to(np.asarray(z0, dtype=np.complex128), rho.shape).copy()
        reach: Arr = np.zeros(rho.shape)
        for _ in range(laps):
            along = z[:, None] + rho[:, None] * self.g[None, :] + self.h[None, :]
            reach = np.maximum(reach, np.abs(along).max(axis=1))
            z = self.lap(rho, z)
        return reach

    def matched_excursion(self, rho: Arr) -> Arr:
        """Largest distance from the curve (m) of the closed orbit itself; inf if there is none."""
        if not self.closes:
            return np.full(rho.shape, math.inf)
        return self.excursion(rho, self.closed_orbit(rho), 1)


def window(rho: Arr, reach: Arr, wall: float) -> tuple[float, float] | None:
    """The range of rho whose reach stays inside the wall, or None. GDE-7: reach is convex in rho, so the values
    of rho that stay inside form one unbroken range."""
    inside = rho[reach < wall]
    return (float(inside.min()), float(inside.max())) if inside.size else None


# ---- following electrons round a closed guide ---------------------------------------------------------------------
@dataclass(frozen=True, eq=False)
class Launch:
    """Electrons to send round a guide: kinetic energies (eV), start offsets from the curve's first point (N, 3),
    and senses (+1 along the field, -1 against it)."""

    w_ev: Arr
    offset: Arr
    sense: Arr

    @classmethod
    def of(cls, w_ev: Arr, offset: Arr, sense: Arr | None = None) -> Launch:
        """One electron per energy; a single offset (3,) is shared by all."""
        n = len(w_ev)
        return cls(w_ev, np.broadcast_to(offset, (n, 3)).copy(), np.ones(n) if sense is None else sense)


class GuideTrack(NamedTuple):
    """What track_guide measured, one entry per electron."""

    held: BoolArr  # still inside the wall at the end
    laps: Arr  # laps completed (negative: against T)
    reach: Arr  # largest distance from the curve while held (m)
    path: Arr  # (m, 3) positions of the recorded electron; empty if none was recorded
    section: tuple[CArr, ...]  # per electron, its offset on the curve's axes each time it completes a lap


class _State(NamedTuple):
    x: Arr
    u: Arr
    st: IntArr
    adv: Arr
    alive: BoolArr
    reach: Arr


def track_guide(
    curve: ClosedCurve,
    field: BField,
    bt: float,
    launch: Launch,
    wall: float,
    laps: float,
    *,
    steps_per_turn: int = 20,
    record: int | None = None,
    record_every: int = 3,
) -> GuideTrack:
    """Electrons in a static magnetic field that guides them round a closed curve, by full-orbit integration.

    Each electron starts at the curve's first point plus its offset, moving along the local field (against it if
    its sense is -1), and counts as lost beyond `wall` metres from the curve. The run lasts until the slowest has
    had time for `laps` laps. bt is the guide field strength (T); it sets the time step, 1/steps_per_turn of the
    gyration period.
    """
    N = len(launch.w_ev)
    v = np.array([v_of(float(w)) for w in launch.w_ev])
    x0 = np.tile(curve.points[0], (N, 1)) + launch.offset
    b0 = field(x0)
    u0 = (v / np.sqrt(1 - (v / C0) ** 2))[:, None] * b0 / np.linalg.norm(b0, axis=1)[:, None] * launch.sense[:, None]
    qm, no_e = np.full(N, QM_E), np.zeros((N, 3))
    dt = 2 * math.pi / (E / ME * bt) / steps_per_turn

    def step(s: _State) -> _State:
        xn, un = push(s.x, s.u, qm, dt, no_e, field(s.x))
        x, u = np.where(s.alive[:, None], xn, s.x), np.where(s.alive[:, None], un, s.u)
        st, adv = curve.advance(x, s.st, s.adv, iters=4)
        off = curve.offset(x, st)[1]
        return _State(x, u, st, adv, s.alive & (off < wall), np.where(s.alive, np.maximum(s.reach, off), s.reach))

    s = _State(x0, u0, np.zeros(N, dtype=np.intp), np.zeros(N), np.ones(N, dtype=np.bool_), np.zeros(N))
    path: list[Arr] = []
    (e1, e2), crossings = curve.axes, [list[complex]() for _ in range(N)]
    done = np.zeros(N, dtype=np.intp)
    for it in range(int(laps * curve.length / v.min() / dt)):
        s = step(s)
        now = (np.abs(s.adv) // len(curve)).astype(np.intp)
        for i in np.flatnonzero((now > done) & s.alive):
            d = curve.offset(s.x[i : i + 1], s.st[i : i + 1])[0][0]
            crossings[i].append(complex(float(d @ e1[s.st[i]]), float(d @ e2[s.st[i]])))
        done = np.maximum(done, now)
        if record is not None and it % record_every == 0 and s.alive[record]:
            path.append(s.x[record].copy())
        if not s.alive.any():
            break
    section = tuple(np.array(c, dtype=np.complex128) for c in crossings)
    return GuideTrack(s.alive, s.adv / len(curve), s.reach, np.array(path).reshape(-1, 3), section)


class OrbitFit(NamedTuple):
    """The circle through an orbit's offsets at the start of successive laps."""

    centre: complex  # the closed orbit the offsets turn about (m)
    radius: float  # their distance from it (m)
    turn: float  # the angle between successive offsets, seen from the centre (rad, in (-pi, pi])
    residual: float  # rms distance of the offsets from the circle (m)


def orbit_circle(z: CArr) -> OrbitFit:
    """Fit a circle to the offsets z at the start of successive laps, and measure the turn per lap about its centre.

    GDE-3 and GDE-4 say the offsets of a first-order orbit lie on a circle about the closed orbit and advance by
    the plane's turn each lap, so this measures both from a tracked orbit (GuideTrack.section). It needs at least
    three offsets that are not in line; an orbit launched on the closed orbit gives a point, and no turn."""
    design = np.column_stack([2 * z.real, 2 * z.imag, np.ones(len(z))])
    (a, b, k), *_ = np.linalg.lstsq(design, np.abs(z) ** 2, rcond=None)
    centre = complex(float(a), float(b))
    radius = math.sqrt(max(float(k) + abs(centre) ** 2, 0.0))
    w = z - centre
    turn = cmath.phase(complex(np.sum(w[1:] * np.conj(w[:-1]))))
    return OrbitFit(centre, radius, turn, float(np.sqrt(np.mean((np.abs(w) - radius) ** 2))))
