"""Filament windings round a closed curve, and their Biot-Savart field."""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import cached_property

import numpy as np

from .base import MU0, Arr, Winding
from .curves import ClosedCurve
from .optics import GuideOptics

_BLOCK = 2_000_000  # points x segments evaluated at once: about 50 MB of working arrays


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
        outer = radius + (layer_gap if winding == "return" else 0.0)
        if outer >= curve.bend_radius:
            raise ValueError(
                f"a winding of radius {outer * 1e3:.1f} mm cannot go round a bend of {curve.bend_radius * 1e3:.1f} mm: "
                "its turns would cross on the inside of the bend"
            )
        mid, dl = _helix(curve, radius, turns, per_turn, 1)
        if winding == "return":
            mid2, dl2 = _helix(curve, radius + layer_gap, turns, per_turn, -1)
            mid, dl = np.vstack([mid, mid2]), np.vstack([dl, -dl2])
        made = cls(curve, radius, turns * (2 if winding == "return" else 1), winding, mid, dl)
        along = float(np.mean(np.sum(made.axis_field * curve.tangents, 1)))
        return made if along >= 0 else cls(curve, radius, made.turns, winding, mid, -dl)

    @cached_property
    def _sources(self) -> tuple[Arr, Arr, Arr, Arr]:
        """The segments laid out for the sum in `field`: their midpoints by component, and each segment's dl beside
        its dl x midpoint, (M, 6)."""
        mid, dl = self.mid, self.dl
        mx, my, mz = (np.ascontiguousarray(mid[:, k]) for k in range(3))
        return mx, my, mz, np.ascontiguousarray(np.hstack([dl, np.cross(dl, mid)]))

    def field(self, x: Arr) -> Arr:
        """Biot-Savart field (N, 3), tesla per ampere, each segment taken at its midpoint.

        dl x (x - mid) = dl x x - dl x mid, so with K[n, m] = |x_n - mid_m|^-3 the sum over the segments is
        (K dl) x x - K (dl x mid): one matrix product, and no power function. It is the plain sum of field_plain to
        rounding (COI-3), six to eight times faster. Points are taken a block at a time to bound the memory."""
        x = np.ascontiguousarray(x, dtype=np.float64)
        block = max(1, _BLOCK // len(self.mid))
        if len(x) > block:
            return np.vstack([self.field(x[i : i + block]) for i in range(0, len(x), block)])
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
        out = np.empty_like(x)
        out[:, 0] = s[:, 1] * x[:, 2] - s[:, 2] * x[:, 1]
        out[:, 1] = s[:, 2] * x[:, 0] - s[:, 0] * x[:, 2]
        out[:, 2] = s[:, 0] * x[:, 1] - s[:, 1] * x[:, 0]
        out -= s[:, 3:]
        out *= MU0 / (4 * math.pi)
        return out

    def field_plain(self, x: Arr) -> Arr:
        """The same field by the sum as it is written in a textbook: the reference `field` is tested against."""
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
