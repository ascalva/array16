"""First-order theory of a guide field that follows a closed curve: the lap map, its closed orbit, its reach.

The laws of the lap map are GDE-1 to GDE-7 in docs/LAWS.md, proved in formal/Array16/Guide.lean. That a guiding
centre obeys it is a hypothesis (GDE-H1), tested by tracking."""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from typing import NamedTuple

import numpy as np

from .base import Arr, CArr
from .curves import ClosedCurve


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

    def unlap(self, rho: Arr, z: CArr) -> CArr:
        """The offsets one lap earlier: the inverse of lap. It is also one lap later for an electron going round
        against the field, at its (negative) rho."""
        before: CArr = cmath.exp(-1j * self.alpha) * z - self.step(rho)
        return before

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

    def crossings(self, rho: float, z0: complex, laps: int) -> CArr:
        """Where a guiding centre launched at z0 crosses the start after each of `laps` laps, on the start's axes:
        what GuideTrack.section measures. The sign of rho gives the direction of travel."""
        r, z, out = np.array([rho]), np.array([z0], dtype=np.complex128), []
        for _ in range(laps):
            z = self.lap(r, z) if rho >= 0 else self.unlap(r, z)
            out.append(complex(z[0]))
        return np.array(out, dtype=np.complex128)

    def distances(self, rho: float, z0: complex, laps: int) -> Arr:
        """The guiding centre's distance from the curve (m) along its path, in the order travelled: (laps, n + 1),
        from the start of each lap to its end. The sign of rho gives the direction of travel.

        Along the field a lap that starts at offset z passes station i at z + rho g_i + h_i. Against it, a lap
        that will end at offset w passes the stations in falling order at w + rho g_i + h_i."""
        starts = np.concatenate([[z0], self.crossings(rho, z0, laps)])
        if rho >= 0:
            return np.abs(starts[:-1, None] + rho * self.g[None, :] + self.h[None, :])
        return np.abs(starts[1:, None] + rho * self.g[None, ::-1] + self.h[None, ::-1])

    def exit_lap(self, rho: float, z0: complex, wall: float, laps: int) -> float | None:
        """The lap, with its fraction, at which a guiding centre launched at z0 first comes `wall` metres from the
        curve, or None if it has not in `laps` laps."""
        far = self.distances(rho, z0, laps) >= wall
        if not far.any():
            return None
        lap, station = np.unravel_index(int(np.argmax(far)), far.shape)
        return float(lap) + float(station) / (len(self.g) - 1)

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
