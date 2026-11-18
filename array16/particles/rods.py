"""Thin rods along a curve, as line charges."""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import cached_property

import numpy as np

from .base import EPS0, Arr, IntArr
from .curves import ClosedCurve, Frame


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

    def gradient(self, lam: float) -> float:
        """The field gradient on the axis (V/m^2) for rods of line charge +/- lam (C/m): 2 lam / (pi eps0 r0^2).
        Near the axis the field is gradient x (-a, +b) on the frame's axes (ROD-1)."""
        return 2 * lam / (math.pi * EPS0 * self.r0**2)

    def line_charge(self, gradient: float) -> float:
        """The line charge (C/m) that gives a field gradient on the axis. Inverse of gradient."""
        return math.pi * EPS0 * self.r0**2 * gradient / 2

    def volts(self, lam: float, wire_radius: float) -> float:
        """The voltage between neighbouring rods of the given wire radius carrying +/- lam:
        lam ln(r0 / wire_radius) / (pi eps0)."""
        return lam * math.log(self.r0 / wire_radius) / (math.pi * EPS0)

    @cached_property
    def _sources(self) -> tuple[Arr, Arr, Arr, Arr]:
        """Every segment of the four rods laid out for field_all: midpoints by component, and each segment's
        signed length beside that times its midpoint, (4n, 4)."""
        mid = self.mid.reshape(-1, 3)
        w = (self.sign[:, None] * self.dl).reshape(-1, 1)
        mx, my, mz = (np.ascontiguousarray(mid[:, k]) for k in range(3))
        return mx, my, mz, np.ascontiguousarray(np.hstack([w, w * mid]))

    def field_all(self, x: Arr, lam: float | Arr) -> Arr:
        """Electric field (N, 3) of the whole of the four rods. With K[n, m] = |x_n - mid_m|^-3 and w the signed
        segment lengths, the sum of w (x - mid) / |x - mid|^3 is x (K w) - K (w mid): one matrix product. The same
        as `field` over every station, to rounding."""
        x = np.ascontiguousarray(x, dtype=np.float64)
        mx, my, mz, sources = self._sources
        d = x[:, 0:1] - mx
        k = d * d
        for component, m in ((1, my), (2, mz)):
            np.subtract(x[:, component : component + 1], m, out=d)
            d *= d
            k += d
        np.sqrt(k, out=d)
        k *= d
        np.reciprocal(k, out=k)
        s = k @ sources
        out: Arr = x * s[:, 0:1] - s[:, 1:]
        scale: Arr = (np.asarray(lam, dtype=np.float64) / (4 * math.pi * EPS0)).reshape(-1, 1)
        return out * scale

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
