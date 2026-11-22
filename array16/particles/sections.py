"""A solved cross-section swept along a path: the field of a guide whose section stays the same along it, as a Fields
for integrate.

A guide of boards keeps one cross-section along its path, turned with the path's frame. The section's RF potential per
volt, phi(u, v), comes from a two-dimensional solve on a tensor grid, or from a closed form sampled on one
(electrodes.FacingBoards, electrodes.RailPair): u runs across the path in the boards' plane, along the frame's first
normal, and v from the board toward the null, along its second. SectionField sweeps it along a ClosedCurve, or along a
straight line one period of which stands for a lap (a ring unrolled), and drives it with a travelling wave, a standing
wave or the same volts all along. The field is quasi-static: at each arc length it is the section's field at the volts
there, E = -V(s, t) grad phi with grad phi turned into the frame and nothing along the path, as on a line carrying a TEM
wave; a travelling wave carries B = T x E / v, v its phase velocity. On a quasi-TEM line over a laminate that is the
section's E made magnetic, where the line's magnetic pattern is the air line's, (n/c) T x E_air: the two differ in the
vacuum above a board by the laminate's pull on E, a correction to a force that is itself v_e / v of the electric one
(1 % at 9 eV on a line of index 1.55); a design that needs it sweeps the air solve's phi for B.

  SectionGrid    phi per volt on a uniform or graded tensor grid, where its null is, and the cells a particle cannot
                 enter; grad phi at the nodes by second-order differences and bilinearly between them
  SectionField   the grid swept along a path with a frame, the curve's twist-free axes turned by a framing angle
                 theta(s) or normals given per sample, driven by TravellingWave, StandingWave or InPhase (a PathDrive)
  BendShift      opt-in: the facing boards' null moved toward the centre of each bend by its first-order shift,
                 2 ln 2 / pi^2 h^2/R in the boards' plane and h^2/(pi R) toward them (array16/guide/axisym.py)

A particle's place along the path follows from its nearest sample: it lies between that sample's normal plane and a
neighbour's, and its foot is on Hermite's cubic through the two. The nearest sample is found by searching a window round
the one it was nearest at the last call, or every sample where the window does not settle it. The window is trusted
only inside the curve's thickness, where no other strand can be nearer, so the field does not depend on what the
SectionField remembers, only the time it takes. A particle off the grid, or in a solid cell, is pushed by nothing and
flagged (off_grid; stop_off_grid stops it in integrate). SI units.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, field
from functools import cached_property
from typing import NamedTuple, Protocol

import numpy as np

from .base import Arr, BoolArr, IntArr
from .curves import ClosedCurve

_BLOCK = 1 << 21  # entries of a particles-by-samples block in a search of every sample


def _dot(a: Arr, b: Arr) -> Arr:
    """Row-wise dot products of two (N, 3) arrays."""
    out: Arr = np.einsum("ij,ij->i", a, b)
    return out


# ---- drives ---------------------------------------------------------------------------------------------------------
class PathDrive(Protocol):
    """The volts along a swept section's path: V (N,) at angles 2 pi s / L round it (N,) and a time t (s); and the
    slowness 1/v (s/m) of the wave whose magnetic field is T x E / v, or None where the drive carries none."""

    def __call__(self, angle: Arr, t: float) -> Arr: ...

    def slowness(self, length: float) -> float | None: ...


@dataclass(frozen=True)
class TravellingWave:
    """A travelling wave of order m round the path: V0 cos(m angle - omega t + phase), angle = 2 pi s / L, V0 = peak
    (V). It runs along T for m > 0 and against it for m < 0 at the phase velocity omega L / (2 pi m), which is c / n_eff
    when omega is the ring's m-th resonance (T915-S at m = 1: 915 MHz x 0.318 m = 2.91e8 m/s, c / 1.031), and carries
    the magnetic field of a wave on a line, T x E / v."""

    peak: float
    omega: float
    m: int
    phase: float = 0.0

    def __post_init__(self) -> None:
        if self.m == 0 or not self.omega > 0:
            raise ValueError(
                f"a travelling wave needs an order m other than 0 and omega > 0, not m = {self.m}, omega = {self.omega}"
            )

    def __call__(self, angle: Arr, t: float) -> Arr:
        out: Arr = self.peak * np.cos(self.m * angle - self.omega * t + self.phase)
        return out

    def phase_velocity(self, length: float) -> float:
        """omega L / (2 pi m) (m/s) on a path of length L: negative for a wave against T."""
        return self.omega * length / (2 * math.pi * self.m)

    def slowness(self, length: float) -> float | None:
        return 1 / self.phase_velocity(length)


@dataclass(frozen=True)
class StandingWave:
    """A standing wave of order m round the path: V0 cos(m angle + spatial) cos(omega t - phase), with nodes where
    m angle + spatial is an odd multiple of pi / 2. Quasi-static: the magnetic field of its current, a quarter period
    and a quarter wave from E, is left out (B is None); two TravellingWaves of orders m and -m carry it."""

    peak: float
    omega: float
    m: int
    phase: float = 0.0
    spatial: float = 0.0

    def __call__(self, angle: Arr, t: float) -> Arr:
        out: Arr = self.peak * math.cos(self.omega * t - self.phase) * np.cos(self.m * angle + self.spatial)
        return out

    def slowness(self, length: float) -> float | None:
        return None


@dataclass(frozen=True)
class InPhase:
    """The same volts all along the path, a constant amplitude: V0 cos(omega t - phase). omega = 0 is a static drive,
    V0 cos(phase). No magnetic field (B is None)."""

    peak: float
    omega: float
    phase: float = 0.0

    def __call__(self, angle: Arr, t: float) -> Arr:
        return np.full(angle.shape, self.peak * math.cos(self.omega * t - self.phase))

    def slowness(self, length: float) -> float | None:
        return None


# ---- the bend's shift of the null -----------------------------------------------------------------------------------
BEND_IN_PLANE = 2 * math.log(2) / math.pi**2  # 0.140461: the null's shift toward a bend's centre in the boards' plane,
# over h^2/R, to first order in h/R (Proposition, array16/guide/axisym.py; E2-2's 0.1405)
BEND_TOWARD = 1 / math.pi  # 0.318310: the same for a bend toward the boards, out of their plane (E2-2's 0.318)


@dataclass(frozen=True)
class BendShift:
    """The first-order bend correction of a facing-board section swept along a curved path: the null moved toward the
    centre of the bend by h^2 (inplane kappa_u, toward kappa_v) in the section, where kappa_u = kappa . first and
    kappa_v = kappa . second are the path's curvature vector dT/ds (1/m) in the boards' plane and toward them, and h is
    the null's distance to either board (m). The section's pattern moves with its null, copper included: a particle at
    (u, v) feels the section's field at (u - du, v - dv). Near the null, to first order in h/R and in the distance from
    it, that is the bent boards' field: the bend adds a uniform field there, which the shift cancels, and no gradient.
    Farther out the bend's first-order terms of third degree make the well lopsided across the bend, which this leaves
    to the boundary-element solve (array16/guide/bem.py), and the copper is misplaced by the shift itself, 7.7 um in the
    boards' plane and 17 um toward them on a bend of 15.3 mm at h = 0.915 mm, T915-S's tightest.

    The coefficients are the gapless facing boards' (ELE-2), derived and held to the axisymmetric reference by
    quadrature in array16/guide/axisym.py: 2 ln 2 / pi^2 = 0.1405 in the boards' plane and 1/pi = 0.3183 toward them
    (Propositions; E2-2). The first order is linear in the curvature vector, so a bend at an angle psi to the boards'
    plane moves the null by (h^2/R) (inplane cos psi, toward sin psi): (inplane cos^2 psi + toward sin^2 psi) h^2/R
    toward the centre and (toward - inplane) sin psi cos psi h^2/R across that. The gradient at the null is even in the
    curvature, so it does not change at first order: the reference puts its change at -0.084 (h/R)^2 in the boards'
    plane and +0.069 (h/R)^2 toward them, under 3.1e-4 on T915-S's tightest bend, and this correction leaves it out.
    Another section (the slot guide, one board's rails) takes its own coefficients."""

    h: float
    inplane: float = BEND_IN_PLANE
    toward: float = BEND_TOWARD

    def __post_init__(self) -> None:
        if not (math.isfinite(self.h) and self.h > 0 and math.isfinite(self.inplane) and math.isfinite(self.toward)):
            raise ValueError(
                f"a bend correction takes the null's distance to the boards h > 0 (m) and finite coefficients, not "
                f"h = {self.h}, inplane = {self.inplane}, toward = {self.toward}"
            )

    def displacement(self, kappa_u: Arr, kappa_v: Arr) -> tuple[Arr, Arr]:
        """The null's displacement (du, dv) (m) in the section where the path's curvature is kappa_u in the boards'
        plane and kappa_v toward them (1/m, along u and v)."""
        h2 = self.h * self.h
        du: Arr = self.inplane * h2 * np.asarray(kappa_u, dtype=np.float64)
        dv: Arr = self.toward * h2 * np.asarray(kappa_v, dtype=np.float64)
        return du, dv


# ---- the section ----------------------------------------------------------------------------------------------------
def _even(nodes: Arr) -> bool:
    """Whether nodes are evenly spaced, to a billionth of their step."""
    steps = np.diff(nodes)
    return bool(np.ptp(steps) <= 1e-9 * float(steps.mean()))


def _cells(nodes: Arr, p: Arr, even: bool) -> IntArr:
    """The cell [nodes[k], nodes[k + 1]] each finite coordinate p falls in, clipped to the grid's cells: directly on
    an even grid, by bisection on a graded one."""
    k = (
        np.floor((p - nodes[0]) * ((len(nodes) - 1) / (nodes[-1] - nodes[0])))
        if even
        else np.searchsorted(nodes, p, side="right") - 1
    )
    out: IntArr = np.clip(k, 0, len(nodes) - 2).astype(np.intp)
    return out


@dataclass(frozen=True, eq=False)
class SectionGrid:
    """A guide's cross-section: its RF potential per volt phi (Nu, Nv) on the nodes of a tensor grid u (Nu,) by v (Nv,)
    (m, each increasing, evenly spaced or graded; three nodes or more), u across the path in the boards' plane and v
    from the board toward the null, and the null (u0, v0) the path runs through. solid (Nu - 1, Nv - 1), if given,
    marks the cells a particle cannot enter (copper, laminate): one in such a cell has hit something, as one beyond the
    grid has left it.

    The field per volt is -grad phi, taken at the nodes by second-order differences (np.gradient, one-sided at the
    edges) and bilinearly between them, as the knotted study's `bilinear` interpolates on its uniform grids. For a
    harmonic phi it is off by at most (5/12) d^2 max|phi'''| a component on a grid of step d: d^2 / 6 from the
    differences, d^2 / 4 from the interpolation."""

    u: Arr
    v: Arr
    phi: Arr
    null: tuple[float, float]
    solid: BoolArr | None = None

    def __post_init__(self) -> None:
        u, v = np.asarray(self.u, dtype=np.float64), np.asarray(self.v, dtype=np.float64)
        phi = np.asarray(self.phi, dtype=np.float64)
        if u.ndim != 1 or v.ndim != 1 or len(u) < 3 or len(v) < 3 or phi.shape != (len(u), len(v)):
            raise ValueError(
                "a section grid needs nodes u (Nu,) and v (Nv,), three or more each, and phi (Nu, Nv): "
                f"got {u.shape}, {v.shape} and {phi.shape}"
            )
        if not (np.all(np.diff(u) > 0) and np.all(np.diff(v) > 0)):
            raise ValueError("a section grid's nodes must increase along u and along v")
        if not np.all(np.isfinite(phi)):
            raise ValueError("a section grid's potential must be finite at every node")
        u0, v0 = float(self.null[0]), float(self.null[1])
        if not (u[0] <= u0 <= u[-1] and v[0] <= v0 <= v[-1]):
            raise ValueError(f"the null ({u0}, {v0}) lies off the grid [{u[0]}, {u[-1]}] x [{v[0]}, {v[-1]}]")
        object.__setattr__(self, "u", u)
        object.__setattr__(self, "v", v)
        object.__setattr__(self, "phi", phi)
        object.__setattr__(self, "null", (u0, v0))
        if self.solid is not None:
            solid = np.asarray(self.solid, dtype=np.bool_)
            if solid.shape != (len(u) - 1, len(v) - 1):
                raise ValueError(f"solid marks cells, ({len(u) - 1}, {len(v) - 1}), not {solid.shape}")
            object.__setattr__(self, "solid", solid)

    @classmethod
    def sampled(
        cls,
        potential: Callable[[Arr, Arr], Arr],
        u: Arr,
        v: Arr,
        null: tuple[float, float],
        solid: BoolArr | None = None,
    ) -> SectionGrid:
        """A potential per volt known in closed form, sampled on the grid's nodes: potential(U, V) with U and V the
        nodes' coordinates (Nu, Nv). FacingBoards(2h).potential(u, v, 1.0) with the null at (0, 0), or
        RailPair(gap, rail).potential(u, v, 1.0) with the null at (0, null_height) and v > 0."""
        uu, vv = np.meshgrid(np.asarray(u, dtype=np.float64), np.asarray(v, dtype=np.float64), indexing="ij")
        return cls(uu[:, 0], vv[0], np.asarray(potential(uu, vv), dtype=np.float64), null, solid)

    @classmethod
    def mirrored(
        cls,
        u: Arr,
        v: Arr,
        phi: Arr,
        null: tuple[float, float],
        *,
        u_parity: int = 0,
        v_parity: int = 0,
        solid: BoolArr | None = None,
    ) -> SectionGrid:
        """A section solved on a half or a quarter, unfolded by its parities: phi(-u, v) = u_parity phi(u, v) with
        u_parity +1 or -1 (u then starts at 0), and the same in v; 0 leaves an axis as it is. The knotted study solves
        its facing strips on the quarter, odd in both (-1, -1), its slot guide on the half x >= 0, odd across the slot
        (-1, 0); one board's rails are even across (+1, 0)."""
        uu, vv = np.asarray(u, dtype=np.float64), np.asarray(v, dtype=np.float64)
        pp = np.asarray(phi, dtype=np.float64)
        cells = None if solid is None else np.asarray(solid, dtype=np.bool_)
        for axis, parity in ((0, u_parity), (1, v_parity)):
            if parity == 0:
                continue
            nodes = uu if axis == 0 else vv
            if parity not in (-1, 1) or abs(nodes[0]) > 1e-12 * (nodes[-1] - nodes[0]):
                raise ValueError(f"mirroring takes a parity of +1 or -1 and nodes from 0, not {parity} from {nodes[0]}")
            nodes = np.concatenate([-nodes[:0:-1], [0.0], nodes[1:]])
            flip = pp[:0:-1] if axis == 0 else pp[:, :0:-1]
            pp = np.concatenate([parity * flip, pp], axis=axis)
            if cells is not None:
                cells = np.concatenate([cells[::-1] if axis == 0 else cells[:, ::-1], cells], axis=axis)
            if axis == 0:
                uu = nodes
            else:
                vv = nodes
        return cls(uu, vv, pp, null, cells)

    @cached_property
    def _even(self) -> tuple[bool, bool]:
        return _even(self.u), _even(self.v)

    @cached_property
    def _gradient(self) -> Arr:
        """d phi / du and d phi / dv at the nodes (Nu, Nv, 2), per metre."""
        gu, gv = np.gradient(self.phi, self.u, self.v, edge_order=2)
        out: Arr = np.ascontiguousarray(np.stack([gu, gv], axis=-1))
        return out

    @cached_property
    def _corners(self) -> Arr:
        """Each cell's four corners' gradients, ((Nu - 1)(Nv - 1), 8): (i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1),
        so that one gather fetches a cell."""
        g = self._gradient
        out: Arr = np.ascontiguousarray(
            np.concatenate([g[:-1, :-1], g[1:, :-1], g[:-1, 1:], g[1:, 1:]], axis=-1).reshape(-1, 8)
        )
        return out

    @property
    def reach(self) -> float:
        """The farthest corner of the grid from the null (m): no particle farther than this from the path is on it."""
        u0, v0 = self.null
        return max(math.hypot(a - u0, b - v0) for a in (self.u[0], self.u[-1]) for b in (self.v[0], self.v[-1]))

    def _locate(self, pu: Arr, pv: Arr) -> tuple[IntArr, IntArr, Arr, Arr, BoolArr]:
        """Each point's cell (i, j), its fractions across the cell, and whether it is on the grid: inside it, finite,
        and not in a solid cell."""
        u, v = self.u, self.v
        on: BoolArr = np.isfinite(pu) & np.isfinite(pv)
        pu, pv = np.where(on, pu, u[0]), np.where(on, pv, v[0])
        on &= (pu >= u[0]) & (pu <= u[-1]) & (pv >= v[0]) & (pv <= v[-1])
        even_u, even_v = self._even
        i, j = _cells(u, pu, even_u), _cells(v, pv, even_v)
        fu = (pu - u[i]) / (u[i + 1] - u[i])
        fv = (pv - v[j]) / (v[j + 1] - v[j])
        if self.solid is not None:
            on &= ~self.solid[i, j]
        return i, j, fu, fv, on

    def _gradient_at(self, pu: Arr, pv: Arr) -> tuple[Arr, BoolArr]:
        """grad phi (N, 2) at points (u, v), bilinearly between the nodes, zero off the grid; and which are on it."""
        i, j, fu, fv, on = self._locate(pu, pv)
        g = self._corners[i * (len(self.v) - 1) + j]
        a, b = (1 - fu)[:, None], fu[:, None]
        c, d = (1 - fv)[:, None], fv[:, None]
        out: Arr = c * (a * g[:, 0:2] + b * g[:, 2:4]) + d * (a * g[:, 4:6] + b * g[:, 6:8])
        out[~on] = 0.0
        return out, on

    def gradient_at(self, pu: Arr, pv: Arr) -> tuple[Arr, Arr, BoolArr]:
        """d phi / du and d phi / dv (1/m, per volt) at points (u, v), each (N,), bilinearly between the nodes' values
        and zero off the grid; and which points are on it (inside it and not in a solid cell)."""
        g, on = self._gradient_at(np.asarray(pu, dtype=np.float64), np.asarray(pv, dtype=np.float64))
        return g[:, 0], g[:, 1], on

    def potential_at(self, pu: Arr, pv: Arr) -> Arr:
        """phi per volt at points (u, v) (N,), bilinearly between the nodes; zero off the grid."""
        i, j, fu, fv, on = self._locate(np.asarray(pu, dtype=np.float64), np.asarray(pv, dtype=np.float64))
        p, ru, rv = self.phi, 1 - fu, 1 - fv
        out: Arr = ru * (rv * p[i, j] + fv * p[i, j + 1]) + fu * (rv * p[i + 1, j] + fv * p[i + 1, j + 1])
        return np.where(on, out, 0.0)

    def find_null(self, start: tuple[float, float] | None = None, iters: int = 60) -> tuple[float, float]:
        """Where the interpolated gradient vanishes (m), by Newton's method from `start` (the stated null by default)
        on each cell's bilinear field: the null as a swept field sees it. A solve with laminate or copper of some
        thickness moves a section's null off the gapless formula's place; this finds it on the grid."""
        p = np.array(self.null if start is None else start, dtype=np.float64)
        best, best_size = p.copy(), math.inf
        g = self._gradient
        for _ in range(iters):
            i, j, fu, fv, on = self._locate(p[:1], p[1:])
            if not on[0]:
                raise ValueError(f"the search for the null left the grid at ({p[0]}, {p[1]})")
            k, m, a, b = int(i[0]), int(j[0]), float(fu[0]), float(fv[0])
            g00, g10, g01, g11 = g[k, m], g[k + 1, m], g[k, m + 1], g[k + 1, m + 1]
            here = (1 - a) * ((1 - b) * g00 + b * g01) + a * ((1 - b) * g10 + b * g11)
            size = float(np.hypot(here[0], here[1]))
            if size < best_size:
                best, best_size = p.copy(), size
            jac = np.stack(
                [
                    ((1 - b) * (g10 - g00) + b * (g11 - g01)) / (self.u[k + 1] - self.u[k]),
                    ((1 - a) * (g01 - g00) + a * (g11 - g10)) / (self.v[m + 1] - self.v[m]),
                ],
                axis=1,
            )
            step = np.linalg.solve(jac, -here)
            p = p + step
            if float(np.hypot(step[0], step[1])) < 1e-13 * (self.u[-1] - self.u[0]):
                break
        return float(best[0]), float(best[1])


# ---- the swept field ------------------------------------------------------------------------------------------------
class Placement(NamedTuple):
    """Where particles are along a swept section's path: the segment each one's foot lies on (from that sample to the
    next), the foot's arc length s (m, 0 to L), the particle's coordinates (u, v) in the section (m), and the frame
    there, its first and second normals (N, 3)."""

    segment: IntArr
    s: Arr
    u: Arr
    v: Arr
    first: Arr
    second: Arr


class _Path(NamedTuple):
    """The path packed segment by segment, so that one gather fetches a particle's segment, (n, 25): the start sample
    c0 and its tangent t0, the chord to the next sample and that one's tangent t1, the frame's first and second normals
    at the start with their changes to the next sample's (the seam taken out of the last), and the chord's length. And
    the seam."""

    packed: Arr
    seam: float


class _Memory:
    """What a SectionField keeps between calls so that it searches a window rather than every sample: each particle's
    nearest sample; for those found far beyond the grid, where they then were and how far beyond it (m); and which the
    last call found off the grid."""

    def __init__(self) -> None:
        self.nearest: IntArr | None = None
        self.anchor: Arr = np.empty((0, 3))
        self.clearance: Arr = np.empty(0)
        self.off: BoolArr = np.empty(0, dtype=np.bool_)


@dataclass(frozen=True, eq=False)
class SectionField:
    """A SectionGrid swept along a path with a frame and driven by a PathDrive: a Fields for integrate.

    The path is n samples equally spaced in arc length (points and unit tangents, (n, 3)), with the frame's first
    normal (u's direction) and its second, T x first (v's), at each; length is the path's length, a lap. A closed
    curve's samples close on themselves (along); with axis set the path is the straight line through the first sample
    along axis, repeating every `length` (straight). A particle lies between the normal planes of two neighbouring
    samples; its foot is the fraction of the way between them its distances from the planes give, on the cubic through
    the two samples with their tangents (Hermite's), where the frame has turned linearly by that fraction. Planes
    shared by neighbouring segments make a particle's place in the section continuous along the path, and the cubic
    follows the curve to order step^4 / R^3 where the polygon's chord would miss it by step^2 / 8R.

    The frame should come back to itself after a lap, or turned by a symmetry of the section (pi for facing boards,
    odd across and toward them); `seam` is the turn it comes back with, and the last segment turns as the one before
    it, so the field crosses the seam continuously exactly when the section is invariant under that turn.

    The swept null is the section's, on the path wherever the frame closes; the path's curvature and torsion enter
    only through the frame, unless `bend` (a BendShift, off by default) moves the facing boards' null toward each
    bend's centre by its first-order shift, 0.1405 h^2/R in the boards' plane and 0.318 h^2/R toward them, the
    curvature taken at the samples by second differences and in the frame, and linearly between samples as the frame
    turns. A straight path has no bend to correct. A call searches a window of 2 window + 1 samples round each
    particle's last nearest sample, and every sample for any the window does not settle: on the first call, after a
    jump of more than `window` samples, or beyond the curve's thickness (its least bend radius, or half its closest
    approach to itself, whichever is less). With integrate's push, on a 4000-sample trefoil, that is 0.42 us a
    particle-step for 10^4 particles (the push alone 0.08), 0.52 for 10^3 and 0.12 ms a step for one; the first call's
    search of every sample takes 7 us a particle, and a straight path, which needs no search, 0.25 us a
    particle-step."""

    grid: SectionGrid
    points: Arr
    tangents: Arr
    first: Arr
    second: Arr
    length: float
    drive: PathDrive
    axis: Arr | None = None
    window: int = 2
    bend: BendShift | None = None
    _memory: _Memory = field(default_factory=_Memory, init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        n = len(self.points)
        shape = (n, 3)
        if self.points.shape != shape or any(a.shape != shape for a in (self.tangents, self.first, self.second)):
            raise ValueError(
                "a swept section takes points, tangents and both normals as (n, 3) arrays alike, not "
                f"{self.points.shape}, {self.tangents.shape}, {self.first.shape}, {self.second.shape}"
            )
        if n < (2 if self.axis is not None else 3) or not self.length > 0 or self.window < 1:
            raise ValueError(
                f"a swept section needs 3 samples or more (2 on a straight path), a positive length and a window of at "
                f"least 1: got {n}, {self.length} and {self.window}"
            )
        frame_error = max(
            float(np.abs(_dot(self.first, self.first) - 1).max()),
            float(np.abs(_dot(self.first, self.tangents)).max()),
            float(np.abs(self.second - np.cross(self.tangents, self.first)).max()),
        )
        if frame_error > 1e-6:
            raise ValueError(f"the frame must be unit normals with second = T x first; it is off by {frame_error:.2g}")
        if self.bend is not None and not isinstance(self.bend, BendShift):
            raise ValueError(f"a swept section's bend correction is a BendShift or None, not {self.bend!r}")

    @classmethod
    def along(
        cls,
        grid: SectionGrid,
        curve: ClosedCurve,
        drive: PathDrive,
        *,
        theta: float | Arr = 0.0,
        normals: Arr | None = None,
        window: int = 2,
        bend: BendShift | None = None,
    ) -> SectionField:
        """The section swept along a closed curve. The frame is the curve's twist-free axes (e1, e2) turned by the
        framing angle theta, one for all or one per sample: first = cos theta e1 + sin theta e2, second = T x first.
        Or normals (n, 3) give the first normal at each sample, made unit and normal to T. bend, a BendShift, moves the
        facing boards' null toward each bend's centre by its first-order shift (None, the default: the section's null
        stays on the path).

        A framing that closes after a lap turns by (k pi - holonomy) / L a metre, theta = theta0 + that times s. The
        knotted study's framing(curve, k) gives the board normal n_b = cos theta_b e1 + sin theta_b e2 and the
        direction across the boards' plane T x n_b: normals = its `across` (theta = theta_b + pi / 2), and v then runs
        along -n_b. On a flat ring the twist-free e1 is the outward radius and e2 = T x e1 points down the ring's axis,
        so theta = 0 puts u outward and v down, and theta = pi u inward and v up."""
        t = curve.tangents
        if normals is None:
            e1, e2 = curve.axes
            th = np.broadcast_to(np.asarray(theta, dtype=np.float64), (len(curve),))
            first = np.cos(th)[:, None] * e1 + np.sin(th)[:, None] * e2
        else:
            given = np.asarray(normals, dtype=np.float64)
            first = given - _dot(given, t)[:, None] * t
            first /= np.linalg.norm(first, axis=1)[:, None]
        return cls(grid, curve.points, t, first, np.cross(t, first), curve.length, drive, window=window, bend=bend)

    @classmethod
    def straight(
        cls, grid: SectionGrid, length: float, drive: PathDrive, *, samples: int = 64, theta: float | Arr = 0.0
    ) -> SectionField:
        """The section swept along the z axis from the origin, the field repeating every `length` along it: a lap of a
        ring unrolled. u runs along x and v along y, both turned about z by theta, one angle for all or one per sample
        at z = i length / samples (a twisted strip); the frame between samples turns linearly."""
        s = np.arange(samples) * (length / samples)
        th = np.broadcast_to(np.asarray(theta, dtype=np.float64), (samples,))
        zero = np.zeros(samples)
        points = np.stack([zero, zero, s], 1)
        tangents = np.tile([0.0, 0.0, 1.0], (samples, 1))
        first = np.stack([np.cos(th), np.sin(th), zero], 1)
        second = np.stack([-np.sin(th), np.cos(th), zero], 1)
        return cls(grid, points, tangents, first, second, length, drive, axis=np.array([0.0, 0.0, 1.0]))

    # ---- the path, precomputed -----------------------------------------------------------------------------------
    @cached_property
    def _path(self) -> _Path:
        pts, t, n1, n2 = self.points, self.tangents, self.first, self.second
        nxt = np.roll(pts, -1, 0)
        if self.axis is not None:
            nxt[-1] = pts[0] + self.length * self.axis
        chord = nxt - pts
        t_next, n1_next, n2_next = np.roll(t, -1, 0), np.roll(n1, -1, 0), np.roll(n2, -1, 0)
        # each segment's turn about the tangent: the next sample's first normal against this one carried to it without
        # twisting. On the last segment the frame carried round the lap meets the first sample's turned by the seam,
        # beyond the turn of the segment before; that segment ends on the first sample's frame turned back by it.
        carried = n1 - (t + t_next) * (_dot(n1, t_next) / (1 + _dot(t, t_next)))[:, None]
        turn = np.arctan2(-_dot(carried, n2_next), _dot(carried, n1_next))
        seam = math.remainder(float(turn[-2] - turn[-1]), 2 * math.pi)
        c, s = math.cos(seam), math.sin(seam)
        n1_next[-1] = c * n1[0] + s * n2[0]
        n2_next[-1] = c * n2[0] - s * n1[0]
        span = np.linalg.norm(chord, axis=1)[:, None]
        packed = np.hstack([pts, t, chord, t_next, n1, n1_next - n1, n2, n2_next - n2, span])
        return _Path(np.ascontiguousarray(packed), seam)

    @property
    def seam(self) -> float:
        """The turn (rad, -pi to pi) about the tangent with which the frame, carried once round the path with its own
        twist, comes back against the frame it started from: 0 for a frame that closes, pi for one that closes by an
        odd number of half turns, the curve's holonomy for its twist-free axes."""
        return self._path.seam

    @cached_property
    def curvature_vector(self) -> Arr:
        """The closed path's curvature vector dT/ds at each sample (n, 3) (1/m), by second differences of the samples
        (equally spaced in arc length) without its part along the tangent: zero on a straight path."""
        pts, t = self.points, self.tangents
        if self.axis is not None:
            return np.zeros_like(pts)
        step = self.length / len(pts)
        k = (np.roll(pts, -1, 0) - 2 * pts + np.roll(pts, 1, 0)) / step**2
        out: Arr = k - _dot(k, t)[:, None] * t
        return out

    @cached_property
    def _bend_rows(self) -> Arr:
        """Each segment's curvature in the frame, (n, 4): kappa . first and kappa . second at its start, and their
        changes to its end, where the next sample's curvature is taken on the frame the segment ends on (across the
        seam, the first sample's frame turned by the seam)."""
        rows, k = self._path.packed, self.curvature_vector
        first, second = rows[:, 12:15], rows[:, 18:21]
        nxt = np.roll(k, -1, 0)
        start_u, start_v = _dot(k, first), _dot(k, second)
        end_u, end_v = _dot(nxt, first + rows[:, 15:18]), _dot(nxt, second + rows[:, 21:24])
        out: Arr = np.ascontiguousarray(np.stack([start_u, start_v, end_u - start_u, end_v - start_v], 1))
        return out

    def null_shift(self, where: Placement) -> tuple[Arr, Arr]:
        """The bend correction's displacement of the null (du, dv) (m) in the section at particles where a Placement
        puts them: the curvature in the frame at each segment's ends, linearly between them; zero without `bend` or on
        a straight path."""
        if self.bend is None or self.axis is not None:
            zero = np.zeros(len(where.s))
            return zero, zero.copy()
        rows = self._bend_rows[where.segment]
        f = np.clip(where.s * (len(self.points) / self.length) - where.segment, 0.0, 1.0)
        return self.bend.displacement(rows[:, 0] + f * rows[:, 2], rows[:, 1] + f * rows[:, 3])

    @cached_property
    def _cloud(self) -> tuple[Arr, Arr, Arr]:
        """The samples about their centre, their squared lengths, and the centre: for a search of every sample."""
        centre = self.points.mean(axis=0)
        pts = np.ascontiguousarray(self.points - centre)
        return pts, _dot(pts, pts), centre

    @cached_property
    def _padded(self) -> tuple[Arr, IntArr]:
        """The samples with `window` of them repeated beyond either end, and a window's offsets into them: a window
        round sample i is rows i to i + 2 window."""
        w = self.window
        pts = self.points[np.arange(-w, len(self.points) + w) % len(self.points)]
        offsets: IntArr = np.arange(2 * w + 1, dtype=np.intp)
        return np.ascontiguousarray(pts), offsets

    @cached_property
    def _trust(self) -> float:
        """The distance from the path (m) inside which the nearest sample of a window that does not reach its own edge
        is the nearest of all: the curve's thickness, the smaller of its least bend radius and half its closest
        approach to itself between points more than pi bend radii apart along it (on at most 1000 of its samples), less
        a sample's step and the subsample's."""
        pts, n = self.points, len(self.points)
        ds = self.length / n
        bend = float((np.linalg.norm(np.roll(pts, -1, 0) - 2 * pts + np.roll(pts, 1, 0), axis=1) / ds**2).max())
        radius = 1 / bend if bend > 0 else math.inf
        every = max(1, n // 1000)
        sub = pts[::every] - pts.mean(axis=0)
        index = np.arange(len(sub)) * every
        apart = np.abs(index[:, None] - index[None, :])
        apart = np.minimum(apart, n - apart) * ds
        size = _dot(sub, sub)
        d2 = size[:, None] + size[None, :] - 2 * (sub @ sub.T)
        d2 = np.where(apart > math.pi * radius + every * ds, d2, np.inf)
        approach = math.sqrt(max(float(d2.min()), 0.0)) / 2
        return max(0.0, min(radius, approach) - (1 + every) * ds)

    # ---- where the particles are -----------------------------------------------------------------------------------
    def _window_search(self, x: Arr, guess: IntArr) -> tuple[IntArr, Arr, BoolArr, BoolArr]:
        """Each particle's nearest sample in the window round its guess, the squared distance to it (m^2), whether it
        sits on the window's edge (where a nearer one may lie beyond), and whether the particle is ahead of its normal
        plane: nearer the next sample than the one before. With the tangent along the chord from the one before to the
        next, as ClosedCurve's central differences make it, that is the side of the plane through the sample normal to
        it, the samples being equally spaced."""
        pts, offsets = self._padded
        w = self.window
        d = x[:, None, :] - pts[guess[:, None] + offsets[None, :]]
        d2 = np.einsum("ijk,ijk->ij", d, d)
        k = np.argmin(d2, axis=1)
        rows = np.arange(len(x))
        edge: BoolArr = (k == 0) | (k == 2 * w)
        inner = np.clip(k, 1, 2 * w - 1)
        ahead: BoolArr = d2[rows, inner + 1] < d2[rows, inner - 1]
        nearest: IntArr = (guess + k - w) % len(self.points)
        return nearest, d2[rows, k], edge, ahead

    def _nearest(self, x: Arr) -> tuple[IntArr, Arr, BoolArr]:
        """Each (finite) particle's nearest sample of all, the squared distance to it (m^2), and whether it is ahead of
        that sample's normal plane: the least of |p|^2 - 2 x.p over the samples p, a block of particles at a time,
        settled by a window round it."""
        pts, size, centre = self._cloud
        guess = np.empty(len(x), dtype=np.intp)
        block = max(1, _BLOCK // len(pts))
        for start in range(0, len(x), block):
            rows = x[start : start + block] - centre
            guess[start : start + block] = np.argmin(size[None, :] - 2.0 * (rows @ pts.T), axis=1)
        nearest, d2, _, ahead = self._window_search(x, guess)
        return nearest, d2, ahead

    def _remembered(self, x: Arr) -> tuple[IntArr, BoolArr]:
        """Each particle's nearest sample and whether it is ahead of that sample's normal plane, searching a window
        round the last call's nearest where that is trusted, every sample where it is not, and neither for a particle
        still too far beyond the grid to have come back onto it."""
        mem, count = self._memory, len(x)
        finite: BoolArr = np.isfinite(x).all(axis=1)
        ahead: BoolArr = np.zeros(count, dtype=np.bool_)
        if mem.nearest is None or len(mem.nearest) != count:
            mem.nearest = np.zeros(count, dtype=np.intp)
            mem.clearance, mem.anchor = np.full(count, -1.0), np.zeros((count, 3))
            full = finite
        else:
            look = finite.copy()
            far = np.flatnonzero(mem.clearance > 0)
            if len(far):
                moved = np.sqrt(((x[far] - mem.anchor[far]) ** 2).sum(axis=1))
                look[far[moved < mem.clearance[far]]] = False
            rows: IntArr | slice = slice(None) if bool(look.all()) else np.flatnonzero(look)
            idx, d2, edge, side = self._window_search(x[rows], mem.nearest[rows])
            mem.nearest[rows], mem.clearance[rows], ahead[rows] = idx, -1.0, side
            full = np.zeros(count, dtype=np.bool_)
            full[rows] = edge | ~(d2 < self._trust**2)
        if full.any():
            which = np.flatnonzero(full)
            idx, d2, side = self._nearest(x[which])
            beyond = np.sqrt(d2) - (self.grid.reach + self.length / len(self.points))
            mem.nearest[which], mem.clearance[which] = idx, np.where(beyond > 0, beyond, -1.0)
            mem.anchor[which], ahead[which] = x[which], side
        return mem.nearest, ahead

    def _station(self, d: Arr, index: IntArr, f: Arr, rows: Arr) -> Placement:
        """The station of particles displaced d (N, 3) from their feet, a fraction f along segments `index` whose packed
        rows are `rows`: the frame turned linearly by f, and the particles' coordinates on it."""
        ff = f[:, None]
        first = rows[:, 12:15] + ff * rows[:, 15:18]
        first /= np.linalg.norm(first, axis=1)[:, None]
        second = rows[:, 18:21] + ff * rows[:, 21:24]
        second /= np.linalg.norm(second, axis=1)[:, None]
        u0, v0 = self.grid.null
        s = (index + f) * (self.length / len(self.points))
        return Placement(index, s, u0 + _dot(d, first), v0 + _dot(d, second), first, second)

    def _foot(self, x: Arr, nearest: IntArr, ahead: BoolArr) -> Placement:
        """The station of particles from their nearest samples. Each lies between the normal planes of its nearest
        sample and the next (ahead of the nearest's plane) or the one before (behind it), at the fraction
        f = p0 / (p0 - p1) of the way that its distances p0 and p1 from the two planes give; its foot is that fraction
        along Hermite's cubic through the two samples with their tangents."""
        index: IntArr = np.where(ahead, nearest, nearest - 1) % len(self.points)
        rows = self._path.packed[index]
        rel = x - rows[:, 0:3]
        chord, t0, t1 = rows[:, 6:9], rows[:, 3:6], rows[:, 9:12]
        p0 = _dot(rel, t0)
        span = p0 - _dot(rel - chord, t1)
        f = np.clip(np.divide(p0, span, out=np.zeros_like(p0), where=span > 0), 0.0, 1.0)
        ff = f[:, None]
        back = 1 - ff
        # Hermite: c(f) = c0 + f^2 (3 - 2f) chord + l (f (1 - f)^2 t0 - f^2 (1 - f) t1)
        d = rel - ff * ff * (3 - 2 * ff) * chord - rows[:, 24:25] * ff * back * (back * t0 - ff * t1)
        return self._station(d, index, f, rows)

    def _straight(self, x: Arr) -> Placement:
        """The station of particles on a straight path, from their distance along its axis."""
        assert self.axis is not None
        n = len(self.points)
        rel = x - self.points[0]
        along = rel @ self.axis
        pos = np.mod(along, self.length) * (n / self.length)
        pos = np.where(np.isfinite(pos), pos, 0.0)
        index: IntArr = np.minimum(pos.astype(np.intp), n - 1)
        return self._station(rel - along[:, None] * self.axis, index, pos - index, self._path.packed[index])

    def locate(self, x: Arr) -> Placement:
        """Where particles (N, 3) are along the path and in the section, searching every sample: what a call finds,
        without the memory it keeps."""
        x = np.asarray(x, dtype=np.float64)
        if self.axis is not None:
            return self._straight(x)
        nearest = np.zeros(len(x), dtype=np.intp)
        ahead: BoolArr = np.zeros(len(x), dtype=np.bool_)
        finite: BoolArr = np.isfinite(x).all(axis=1)
        if finite.any():
            idx, _, side = self._nearest(x[finite])
            nearest[finite], ahead[finite] = idx, side
        return self._foot(x, nearest, ahead)

    # ---- the fields --------------------------------------------------------------------------------------------------
    def _fields(self, where: Placement, t: float) -> tuple[Arr, Arr | None, BoolArr]:
        if self.bend is not None and self.axis is None:
            du, dv = self.null_shift(where)
            g, on = self.grid._gradient_at(where.u - du, where.v - dv)
        else:
            g, on = self.grid._gradient_at(where.u, where.v)
        volts = np.where(on, self.drive(where.s * (2 * math.pi / self.length), t), 0.0)
        across, toward = -volts * g[:, 0], -volts * g[:, 1]
        e: Arr = across[:, None] * where.first + toward[:, None] * where.second
        slow = self.drive.slowness(self.length)
        b: Arr | None = None
        if slow is not None:  # T x first = second, T x second = -first
            b = slow * (across[:, None] * where.second - toward[:, None] * where.first)
        if not on.all():  # nothing at all, even where the frame is undefined (a particle at NaN)
            off = ~on
            e[off] = 0.0
            if b is not None:
                b[off] = 0.0
        return e, b, on

    def field_at(self, where: Placement, t: float) -> tuple[Arr, Arr | None]:
        """E (V/m) and B (T, or None) at particles where a Placement puts them, at time t (s): -V(s, t) grad phi(u, v)
        in the frame there, and for a travelling wave T x E / v."""
        e, b, _ = self._fields(where, t)
        return e, b

    def __call__(self, x: Arr, t: float) -> tuple[Arr | None, Arr | None]:
        x = np.asarray(x, dtype=np.float64)
        where = self._straight(x) if self.axis is not None else self._foot(x, *self._remembered(x))
        e, b, on = self._fields(where, t)
        self._memory.off = ~on
        return e, b

    @property
    def off_grid(self) -> BoolArr:
        """The particles the last call found off the grid or in a solid cell, which it pushed with nothing."""
        return self._memory.off.copy()

    def stop_off_grid(self, step: int, t: float, x: Arr, u: Arr, alive: BoolArr) -> BoolArr:
        """An Observer for integrate: stops the particles the last evaluation found off the grid, a step after they
        left it."""
        off = self._memory.off
        return off.copy() if off.shape == alive.shape else np.zeros_like(alive)
