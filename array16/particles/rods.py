"""Thin rods along a curve, as line charges."""

from __future__ import annotations

import math
from dataclasses import dataclass

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
