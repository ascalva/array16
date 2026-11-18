"""Closed space curves: tangents, twist-free frames, holonomy, writhe.

A curve is n points equally spaced in arc length; the point after the last is the first. T is its unit tangent, in
the direction of increasing index."""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import cached_property
from typing import NamedTuple

import numpy as np

from .base import Arr, CArr, IntArr


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
    def curvature(self) -> Arr:
        """The curvature at every point (1/m), by second differences."""
        c, ds = self.points, self.length / len(self)
        k: Arr = np.linalg.norm(np.roll(c, -1, 0) - 2 * c + np.roll(c, 1, 0), axis=1) / ds**2
        return k

    @property
    def bend_radius(self) -> float:
        """Smallest radius of curvature (m)."""
        return float(1 / self.curvature.max())

    @cached_property
    def total_curvature(self) -> float:
        """The curvature summed round the curve (rad). CUR-3: at least 2 pi for any closed curve (Fenchel), and
        more than 4 pi for a knotted one (Fary, Milnor)."""
        return float(self.curvature.sum() * self.length / len(self))

    @cached_property
    def bending(self) -> float:
        """The squared curvature summed round the curve (1/m): its bending energy. CUR-2: at least
        total_curvature^2 / length, so a knot of length L has more than 16 pi^2 / L, four times a circle's."""
        return float(np.sum(self.curvature**2) * self.length / len(self))

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
        """Each particle's displacement from the curve at its station, across the tangent: vectors and lengths.
        It is measured from the tangent line at the station, so it is true to step^2 / (8 x bend radius)."""
        T = self.tangents
        d = x - self.points[st]
        d = d - (d * T[st]).sum(1)[:, None] * T[st]
        return d, np.sqrt((d * d).sum(1))
