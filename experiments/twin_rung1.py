"""Experiment 2, stage 1: the calibration rung of the twin bench (docs/KNOTTED_GUIDE.md section 10).

The first run that tracks electrons through the computed field of a design: each design's cross-section, solved by
array16.guide.sections, swept along its path and driven by its realised ring (array16.twin.TwinGuide.field()), on the
twin bench of benches/twin_r1.toml. Four tests, each a set of shots whose predictions come from the 2-D sections and
formulas alone:

  E2-21      R1's rails' section swept round a ring with an ideal in-phase drive, at 915 MHz: the escape from its well
             at 50 V (q 0.5) as drawn and in the bench's realisation, and at q 0.3 as drawn, the control; and the
             energy it holds round its bend at 50 V, as drawn and on the bench (FR-4's permittivity and Earth's field
             drawn), and at q 0.3 as drawn, the control
  E2-22      the flat facing-board ring of 318 mm and the (1, 3) helix at r/R = 0.2 against T915-S at the same section,
             volts and q: the price of the knot
  slot       the straight slot guide tracked in three dimensions against its 2-D section: secular frequency and depth
  must fail  one board's rails at 915 MHz driven for a 1 eV well (q 1.64), and to q 1.89: no well, nothing held; and
             the same launch at the design's 50 V, the control: every electron held

    uv run python experiments/twin_rung1.py predict                  # theory only: twin_rung1/predictions.{json,md}
    uv run python experiments/twin_rung1.py run --bench twin_r1      # tracking: twin_rung1/runs.jsonl
    uv run python experiments/twin_rung1.py report                   # the verdicts: twin_rung1/results.{json,md}

The theory of a well is the exact motion in the section from the null toward its saddle (section_flights), by RK4 at
480 steps an RF period for the shot's own energies, phases and periods: on R1's rails the saddle stands on the axis of
the section's mirror, where the field has no part across, and the motion stays on it, y'' = |q/m| V(t) dphi/dv; in the
slot guide it is two-dimensional. At q 0.5 on rails it loses electrons far below the pseudopotential's depth D, in
islands of energy: the null's local Mathieu q reaches 0.908 0.20h under the null, where
the pseudopotential is 0.75 D, and an electron swinging back through that zone is pumped out. The theory of the held
energy is the static balance (docs/KNOTTED_GUIDE.md section 6.1, JNT-2c) computed on the solved section rather than
taken from a closed form: under the bend's centrifugal force 2W kappa, and any other uniform force, the electron sits
where the pseudopotential's gradient balances it, and a bend holds W while such an equilibrium exists (the fold of the
tilted well) with the null's local Mathieu q under 0.908. The protocol, with what counts as failure, is
experiments/twin_rung1/PROTOCOL.md. The model is array16.twin's: the section quasi-static and swept along the path, its
null where the 2-D solve puts it, the ring one line mode driven in phase on both rails, one electron at a time, the
bench's residual gas.
"""

from __future__ import annotations

import math
import sys
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, replace
from functools import cache
from pathlib import Path
from typing import Any, Literal

import numpy as np

from array16 import lab
from array16 import particles as pt
from array16.lab import Check
from array16.particles import Arr, BoolArr, SectionGrid
from array16.twin import Design, Driven, RingWave, TwinField, TwinGuide, Watch, curvature_vector

OUT = Path(__file__).resolve().parent / "twin_rung1"
BENCH = "twin_r1"
R1 = "designs/r1_rail_ring.toml"
T915S = "designs/t915s.toml"
FLAT = "designs/flat_ring_318.toml"
HELIX = "designs/helix_13.toml"
SLOT = "experiments/twin_rung1/slot_straight.toml"
MASTER = 20261122  # the master seed: each shot's seed is it plus the shot's place in the list
QM = abs(pt.QM_E)  # an electron's |q/m| (C/kg)

# ---- how it is run ------------------------------------------------------------------------------------------------
STEPS = 120  # steps per RF period (section 10: 120 or more)
LAPS = 20.0  # laps an electron must make to count as held (DR-7)
SCAN = 1.025  # the ratio between neighbouring energies of a scan
PHASES = 2  # RF phases launched at each energy of E2-22's holds
R1_PHASES = 8  # RF phases launched at each energy of R1's holds: 40 electrons a sense in a threshold's window
WELL_PHASES = (math.pi / 2, 3 * math.pi / 2)  # a well's launches: the drive's zero crossings
WELL_PERIODS = 40  # RF periods a well's electron is followed: seven secular periods
BETA_PERIODS = 60  # RF periods the secular frequency is read over, once a period
BETA_OFFSET = 4e-3  # of h: the secular frequency's starts off the null (DR-0's)
FAIL_OFFSET = 0.05  # of h: the must-fail starts, at most this far off the null
FAIL_PERIODS = 200
CONTROL_Q = 0.3  # the controls' drive: the null's Mathieu q at the mean volts
PATCH = (0.8, 1 / 400, 5)  # the well's patch: half-width and node spacing in h, smoothing in nodes
ANGLES = 180  # tilt directions the static strength is computed in, round the circle
MODEL_STEPS = 480  # RK4 steps per RF period of the motion in a section: the exact motion's
WINDOW = 5  # energies in a threshold's window
DRAWS = 13  # the laminate's draws of the tuned band, evenly over +-LAMINATE_SDS, for the bench's well

# ---- what agrees: fixed before the run (PROTOCOL.md, section 5) --------------------------------------------------
DISCOUNT = (0.61, 0.95)  # the tracked fraction of the static limit (section 6.1: four rods, T915-S's design range)
CONTROL_BAND = (0.85, 1.05)  # the tracked fraction of the static limit at q 0.3 (four rods there: 0.96 to 1.03)
RATIO_TOL = 0.15  # E2-22: the tracked ratio within 15 % of the static limits' ratio
WELL_BAND = (0.80, 1.25)  # the slot guide's escape threshold over its predicted depth
WELL_STEPS = 2  # E2-21's wells: the first loss within this many scan steps of the axis model's
BETA_TOL = 2e-3  # the slot guide's tracked secular exponent against Mathieu's
SENSE_STEPS = 2  # E2-21 as drawn: the senses' limits within two scan steps
LAMINATE_SDS = 3.0  # on the bench: the laminate's draw may move the tuned band by this many sds
QCAP = 0.908  # the first stability region of Mathieu's equation at a = 0

Kind = Literal["hold", "well", "beta", "fail", "control"]


def sig(value: float, digits: int = 6) -> float:
    """A number to `digits` significant figures, so that the predictions file does not change with the last bits."""
    return float(f"{value:.{digits}g}") if math.isfinite(value) else value


# ---- theory: Mathieu's exponent ------------------------------------------------------------------------------------
def mathieu_beta(q: float, terms: int = 40) -> float:
    """Mathieu's characteristic exponent at a = 0: the beta of the Floquet solutions e^(i beta tau) P(tau) of
    y'' + 2q cos(2 tau) y = 0, the root in (0, 1) of beta^2 = q (c_2 + c_-2) / c_0 with the two ratios the continued
    fractions of (2j + beta)^2 c_2j = q (c_2j-2 + c_2j+2) (McLachlan; Abramowitz and Stegun 20.3), by bisection. The
    secular frequency is beta f / 2: 0.2161 at q = 0.3 and 0.3737 at q = 0.5, against the first order q / sqrt 2."""

    def balance(beta: float) -> float:
        up = down = 0.0
        for j in range(terms, 0, -1):
            up = q / ((2 * j + beta) ** 2 - q * up)
            down = q / ((beta - 2 * j) ** 2 - q * down)
        return beta * beta - q * (up + down)

    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if balance(mid) < 0 else (lo, mid)
    return (lo + hi) / 2


# ---- theory: what a scan reads -------------------------------------------------------------------------------------
def robust_threshold(
    energies: Sequence[float], held: Sequence[int], lost: Sequence[int], window: int = WINDOW
) -> float:
    """The energy at which a scan starts to lose its electrons, robust to the odd resonant loss: the lowest energy at
    which fewer than half of the electrons judged in the window of `window` energies centred on it (held or lost;
    those the gas stopped or that fell short of their laps count in neither) were held, and its geometric mean with
    the energy below; 0 if that is the lowest energy (below the scan), inf if there is none (above it). A sharp edge
    between two energies gives that edge's geometric mean back."""
    count, half = len(energies), window // 2
    for k in range(count):
        lo, hi = max(0, k - half), min(count, k + half + 1)
        kept = int(sum(held[lo:hi]))
        judged = kept + int(sum(lost[lo:hi]))
        if judged and 2 * kept < judged:
            return 0.0 if k == 0 else math.sqrt(energies[k - 1] * energies[k])
    return math.inf


def lowest_lost(energies: Sequence[float], lost: Sequence[int | bool]) -> float:
    """The lowest energy at which an electron was lost (inf if none was)."""
    return next((float(e) for e, n in zip(energies, lost, strict=True) if n), math.inf)


def edges(energies: Sequence[float], held: Sequence[int], lost: Sequence[int]) -> tuple[float, float]:
    """The highest energy held below the lowest lost, and the lowest lost (eV): 0 if the lowest scanned was lost, inf
    if none was. An energy whose electrons the gas stopped, every one, says nothing."""
    lowest = lowest_lost(energies, lost)
    kept = [e for e, n, gone in zip(energies, held, lost, strict=True) if n and not gone and e < lowest]
    return (max(kept) if kept else 0.0), lowest


def edge_threshold(energies: Sequence[float], held: Sequence[int], lost: Sequence[int]) -> float:
    """The geometric mean of the two edges: 0 below the scan and inf above it."""
    below, above = edges(energies, held, lost)
    if below == 0.0 or not math.isfinite(above):
        return 0.0 if below == 0.0 else math.inf
    return math.sqrt(below * above)


# ---- theory: the static balance on a solved section ---------------------------------------------------------------
def _box(a: Arr, k: int, axis: int) -> Arr:
    """The mean over k nodes centred on each, along an axis, the edges continued."""
    if k <= 1:
        return a
    pad = [(0, 0)] * a.ndim
    pad[axis] = (k // 2 + 1, k // 2)
    total = np.cumsum(np.pad(a, pad, mode="edge"), axis=axis)
    size = total.shape[axis]
    out: Arr = (np.take(total, np.arange(k, size), axis=axis) - np.take(total, np.arange(0, size - k), axis=axis)) / k
    return out


def hessian(grid: SectionGrid) -> tuple[Arr, Arr, Arr]:
    """phi's second derivatives per volt at the grid's nodes (1/m^2): d2/du2, d2/dv2 and d2/du dv, by second-order
    differences of the nodes' second-order gradients."""
    gu, gv = np.gradient(grid.phi, grid.u, grid.v, edge_order=2)
    huu = np.asarray(np.gradient(gu, grid.u, axis=0, edge_order=2), dtype=np.float64)
    hvv = np.asarray(np.gradient(gv, grid.v, axis=1, edge_order=2), dtype=np.float64)
    huv = (np.gradient(gu, grid.v, axis=1, edge_order=2) + np.gradient(gv, grid.u, axis=0, edge_order=2)) / 2
    return huu, hvv, np.asarray(huv, dtype=np.float64)


def stiffness_at(grid: SectionGrid, pu: Arr, pv: Arr) -> Arr:
    """The size of the field's gradient per volt (1/m^2), its largest eigenvalue, at points (u, v): from the nodes'
    second derivatives, bilinear between them, second order everywhere (differencing the bilinear field instead gives
    each cell's mean slope, off by up to half a cell's change); zero off the grid. It sets the local Mathieu q,
    2 |q/m| V stiffness / w^2."""
    huu, hvv, huv = (SectionGrid(grid.u, grid.v, a, grid.null).potential_at(pu, pv) for a in hessian(grid))
    out: Arr = np.hypot((huu - hvv) / 2, huv) + np.abs(huu + hvv) / 2
    return out


@dataclass(frozen=True, eq=False)
class Patch:
    """A section's pseudopotential round its null as the swept field gives it: on a square of nodes (offsets from the
    null, m, the same along u and v), |q/m| |grad phi|^2 / (4 w^2) per volt squared of drive (eV / V^2), smoothed over a
    few nodes either way and NaN off the grid or in its copper and boards; the largest eigenvalue of the field's
    gradient per volt (1/m^2), which sets the local Mathieu q 2 |q/m| V stiff / w^2; and w."""

    nodes: Arr
    well: Arr
    stiff: Arr
    omega: float

    @property
    def step(self) -> float:
        return float(self.nodes[1] - self.nodes[0])

    def sample(self, values: Arr, pu: Arr, pv: Arr) -> Arr:
        """values (on the patch's nodes) bilinearly at offsets (pu, pv) from the null; NaN off the patch."""
        n, first, step = len(self.nodes), float(self.nodes[0]), self.step
        fu, fv = (pu - first) / step, (pv - first) / step
        inside = (fu >= 0) & (fu <= n - 1) & (fv >= 0) & (fv <= n - 1)
        i = np.clip(np.floor(np.where(inside, fu, 0.0)), 0, n - 2).astype(np.intp)
        j = np.clip(np.floor(np.where(inside, fv, 0.0)), 0, n - 2).astype(np.intp)
        a, b = fu - i, fv - j
        out: Arr = (1 - a) * ((1 - b) * values[i, j] + b * values[i, j + 1]) + a * (
            (1 - b) * values[i + 1, j] + b * values[i + 1, j + 1]
        )
        return np.where(inside, out, np.nan)


def well_patch(grid: SectionGrid, f: float, h: float) -> Patch:
    """The patch round a section grid's null at band f, PATCH's size in units of the section's h."""
    half, spacing, smooth = PATCH
    step = spacing * h
    count = round(half / spacing)
    nodes = np.asarray(np.arange(-count, count + 1) * step, dtype=np.float64)
    u0, v0 = grid.null
    uu, vv = np.meshgrid(u0 + nodes, v0 + nodes, indexing="ij")
    du, dv, inside = grid.gradient_at(uu.ravel(), vv.ravel())
    gu, gv, on = du.reshape(uu.shape), dv.reshape(uu.shape), inside.reshape(uu.shape)
    omega = 2 * math.pi * f
    raw = QM * (gu * gu + gv * gv) / (4 * omega * omega)
    stiff = stiffness_at(grid, uu.ravel(), vv.ravel()).reshape(uu.shape)
    weight = on.astype(np.float64)
    total = _box(_box(raw * weight, smooth, 0), smooth, 1)
    cover = _box(_box(weight, smooth, 0), smooth, 1)
    well = np.where(cover > 1 - 1e-9, total / np.maximum(cover, 1e-300), np.nan)
    return Patch(nodes, well, stiff, omega)


@dataclass(frozen=True, eq=False)
class Valley:
    """The equilibria of the well tilted along d = (cos theta, sin theta) in the section (u, v): under the tilt
    U + F d.r, a uniform force F along -d, the electron sits where grad U = -F d. Walking out along -d, the least of U
    across d at each distance is the equilibrium there and F = -dU/dt holds it. u and v (m from the null) are the
    equilibria from the null out, force the tilt each holds (eV/m per V^2; NaN where it is not differenced), stiff the
    field's gradient there (1/m^2 per volt)."""

    theta: float
    u: Arr
    v: Arr
    force: Arr
    stiff: Arr
    omega: float

    def limit(self, volts: float) -> tuple[int, str]:
        """The static limit at `volts`: the index of the fold (the greatest force, where it stops growing), or of the
        greatest force before the local Mathieu q reaches 0.908 (the last equilibrium under the cap, within the
        walk's step of the crossing), or of the greatest force before the walk left the patch; and which of the three
        bound it."""
        q = 2 * QM * abs(volts) * self.stiff / self.omega**2
        best, top = -1, -math.inf
        for i, force in enumerate(self.force):
            if not math.isfinite(force):
                continue
            if q[i] >= QCAP:
                return (best if best >= 0 else i), "q"
            if force > top:
                best, top = i, float(force)
            elif force < 0.98 * top:
                return best, "fold"
        return best, "edge"

    def strength(self, volts: float) -> float:
        """The largest force (eV/m) along -d the well holds at `volts`."""
        return float(self.force[self.limit(volts)[0]]) * volts * volts

    def at(self, force: float, volts: float) -> tuple[float, float]:
        """The equilibrium (u, v) (m from the null) that holds `force` (eV/m) at `volts`; beyond the limit, the
        limit's."""
        end = self.limit(volts)[0]
        known = np.flatnonzero(np.isfinite(self.force[: end + 1]))
        held = np.concatenate([[0.0], self.force[known]])
        u = np.concatenate([[0.0], self.u[known]])
        v = np.concatenate([[0.0], self.v[known]])
        target = force / (volts * volts)
        if target >= held[-1]:
            return float(u[-1]), float(v[-1])
        k = max(1, int(np.searchsorted(held, target)))
        a = (target - held[k - 1]) / (held[k] - held[k - 1])
        return float(u[k - 1] + a * (u[k] - u[k - 1])), float(v[k - 1] + a * (v[k] - v[k - 1]))


def valley(patch: Patch, theta: float, window: int = 5) -> Valley:
    """The valley of `patch` tilted along theta (rad), the force differenced over `window` nodes either way."""
    step = patch.step
    d = np.array([math.cos(theta), math.sin(theta)])
    across = np.array([-d[1], d[0]])
    half = float(patch.nodes[-1])
    span = np.arange(-round(0.5 * half / step), round(0.5 * half / step) + 1) * step
    rows: list[tuple[float, float, float]] = []
    before = 0.0
    for t in -np.arange(0.0, 0.95 * half, step):
        s = before + span
        well = patch.sample(patch.well, t * d[0] + s * across[0], t * d[1] + s * across[1])
        finite = np.isfinite(well)
        if not finite.any():
            break
        k = int(np.argmin(np.where(finite, well, np.inf)))
        if k in (0, len(s) - 1) or not (finite[k - 1] and finite[k + 1]):
            break
        lo, mid, hi = float(well[k - 1]), float(well[k]), float(well[k + 1])
        curve = lo - 2 * mid + hi
        if not curve > 0:
            break
        shift = 0.5 * (lo - hi) / curve
        rows.append((float(t), float(s[k] + shift * step), mid - 0.25 * (lo - hi) * shift))
        before = rows[-1][1]
    table = np.array(rows, dtype=np.float64).reshape(-1, 3)
    t, s, u = table[:, 0], table[:, 1], table[:, 2]
    pu, pv = t * d[0] + s * across[0], t * d[1] + s * across[1]
    force = np.full(len(t), np.nan)
    if len(t) > 2 * window:
        force[window:-window] = -(u[2 * window :] - u[: -2 * window]) / (t[2 * window :] - t[: -2 * window])
    return Valley(theta, pu, pv, force, patch.sample(patch.stiff, pu, pv), patch.omega)


# ---- theory: the motion in a section -------------------------------------------------------------------------------
def section_flights(
    grid: SectionGrid,
    amplitude: float,
    omegas: Arr,
    energies: Sequence[float],
    phases: Sequence[float],
    periods: int,
    aim: tuple[float, float] = (0.0, 1.0),
    steps: int = MODEL_STEPS,
) -> BoolArr:
    """Which electrons leave a well in its section, (len(omegas), len(phases), len(energies)): each starts on the null
    with the speed of its energy along `aim` (a unit vector in (u, v): toward the escape saddle), in the drive
    amplitude cos(w t + phase) at its station, and is followed for `periods` RF periods by RK4, `steps` an RF period,
    on r'' = |q/m| V(t) grad phi (r): the field the swept grid gives (second-order node gradients, bilinear between),
    with nothing along the path. One the grid puts off its nodes or in copper or a board is lost, as the tracker loses
    it. On rails the saddle stands on the axis of the section's mirror, where the field has no part across, so the
    motion stays on the axis and is one-dimensional; in the slot guide the field on the axis points across it. In the
    drive's radians tau = w t, r'' = (|q/m| amplitude / w^2) cos(tau + phase) grad phi, so every band in `omegas`
    (rad/s) shares the steps."""
    u0, v0 = grid.null
    w = np.asarray(omegas, dtype=np.float64)[:, None, None]
    turn = np.asarray(phases, dtype=np.float64)[None, :, None]
    speed = np.array([pt.v_of(e) for e in energies], dtype=np.float64)[None, None, :]
    shape = (w.shape[0], turn.shape[1], speed.shape[2])
    r: Arr = np.stack([np.full(shape, u0), np.full(shape, v0)])
    dr: Arr = np.stack([np.broadcast_to(aim[0] * speed / w, shape), np.broadcast_to(aim[1] * speed / w, shape)])
    scale = np.broadcast_to(QM * amplitude / w**2, shape)
    lost: BoolArr = np.zeros(shape, dtype=np.bool_)
    h = 2 * math.pi / steps

    def pull(at: Arr, tau: float) -> tuple[Arr, BoolArr]:
        du, dv, on = grid.gradient_at(at[0].ravel(), at[1].ravel())
        kick = scale * np.cos(tau + turn)
        return np.stack([kick * du.reshape(shape), kick * dv.reshape(shape)]), on.reshape(shape)

    for n in range(periods * steps):
        tau = n * h
        a1, on = pull(r, tau)
        lost |= ~on
        a2, _ = pull(r + h / 2 * dr, tau + h / 2)
        a3, _ = pull(r + h / 2 * (dr + h / 2 * a1), tau + h / 2)
        a4, _ = pull(r + h * (dr + h / 2 * a2), tau + h)
        r_next = r + h * dr + h * h / 6 * (a1 + a2 + a3)
        dr_next = dr + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
        r, dr = np.where(lost, r, r_next), np.where(lost, dr, dr_next)
    _, on = pull(r, periods * steps * h)
    return lost | ~on


# ---- theory: each design as drawn ---------------------------------------------------------------------------------
def _path(design: str) -> Path:
    path = Path(design)
    return path if path.is_absolute() else lab.ROOT / path


def _wave(field: TwinField) -> RingWave:
    wave = field.swept.drive
    if not isinstance(wave, RingWave):
        raise TypeError("the twin's field is driven by its ring's RingWave")
    return wave


def section_curvature(guide: TwinGuide, field: TwinField) -> tuple[Arr, Arr]:
    """The path's curvature (1/m) at every sample along the swept section's first and second normals, u and v; zero on a
    straight path."""
    swept = field.swept
    if guide.curve is None:
        zero = np.zeros(len(swept.points))
        return zero, zero.copy()
    kv = curvature_vector(guide.curve)
    out: tuple[Arr, Arr] = ((kv * swept.first).sum(1), (kv * swept.second).sum(1))
    return out


def _limits(angles: Arr, strengths: Arr, heading: Arr, kappa: Arr, amplitude: Arr) -> Arr:
    """The static limit (eV) at every path sample: the strength along kappa there, times (|V|/mean)^2, over 2 |kappa|;
    inf where the path is straight."""
    held = np.interp(heading, np.append(angles, 360.0), np.append(strengths, strengths[0]))
    out: Arr = np.where(kappa > 0, held * amplitude**2 / (2 * np.maximum(kappa, 1e-300)), math.inf)
    return out


@dataclass(frozen=True, eq=False)
class Theory:
    """What the theory takes from a design as drawn (realisation 0): the guide and its drive as the design gives it,
    the well's patch, the valleys walked in each tilt direction (deg), the drive's |V| over its mean at every path
    sample, the path's |kappa| (1/m) and kappa's direction in the section (deg) there; and at the mean volts `volts`
    (the design's, or another drive's, `at`): the static strength (eV/m) in each tilt direction, the static limit (eV)
    at every sample, 2W |kappa| against the strength along kappa (inf where the path is straight), and the drive's
    phasor at the start station (V)."""

    design: str
    guide: TwinGuide
    driven: Driven
    patch: Patch
    angles: Arr
    valleys: tuple[Valley, ...]
    amplitude: Arr
    kappa: Arr
    heading: Arr
    volts: float
    strengths: Arr
    limits: Arr
    phasor: complex

    @property
    def static(self) -> float:
        """The path's static limit (eV): the least over its samples."""
        return float(self.limits.min())

    @property
    def where(self) -> int:
        return int(np.argmin(self.limits))

    @property
    def scale(self) -> float:
        """These volts over the design's."""
        return self.volts / self.driven.volts

    @property
    def q(self) -> float:
        """The null's Mathieu q at the mean volts and the tuned band."""
        return self.driven.q * self.scale

    def strength(self, angle: Arr) -> Arr:
        """The static strength (eV/m) at the mean volts along tilt directions (deg), interpolated round the circle."""
        a = np.mod(np.asarray(angle, dtype=np.float64), 360.0)
        out: Arr = np.interp(a, np.append(self.angles, 360.0), np.append(self.strengths, self.strengths[0]))
        return out

    def at(self, volts: float) -> Theory:
        """The same design with its drive's pattern scaled to `volts` mean: the strengths and limits at those volts,
        the local q's cap moving with them."""
        strengths = np.array([v.strength(volts) for v in self.valleys])
        limits = _limits(self.angles, strengths, self.heading, self.kappa, self.amplitude)
        return replace(self, volts=volts, strengths=strengths, limits=limits, phasor=self.phasor * volts / self.volts)


@cache
def theory(design: str, cache: str | None = None) -> Theory:
    """A design as drawn, solved (array16.twin's cache by content) and put through the static balance."""
    guide = TwinGuide(design, 0, cache=cache)
    driven = guide.drive()
    field = guide.field(driven)
    patch = well_patch(guide.grid(), driven.f, guide.solved.distance)
    angles = np.asarray(np.arange(ANGLES) * (360.0 / ANGLES), dtype=np.float64)
    valleys = tuple(valley(patch, math.radians(a)) for a in angles)
    strengths = np.array([v.strength(driven.volts) for v in valleys])
    count = len(field.swept.points)
    wave = _wave(field)
    amplitude: Arr = np.abs(wave.phasor(2 * math.pi * np.arange(count) / count)) / driven.volts
    ku, kv = section_curvature(guide, field)
    kappa = np.hypot(ku, kv)
    heading = np.mod(np.degrees(np.arctan2(kv, ku)), 360.0)
    limits = _limits(angles, strengths, heading, kappa, amplitude)
    phasor = complex(wave.phasor(np.zeros(1))[0])
    return Theory(
        design,
        guide,
        driven,
        patch,
        angles,
        valleys,
        amplitude,
        kappa,
        heading,
        driven.volts,
        strengths,
        limits,
        phasor,
    )


def laminate_sd(guide: TwinGuide) -> float:
    """The relative sd of the ring's tuned band that the design's laminate draw gives: d ln n / d eps times its sd,
    n share / (2 eps) over n."""
    laminate = guide.design.section.laminate
    if laminate is None:
        return 0.0
    return guide.design.tolerance.eps * guide.solved.share / (2 * laminate.er)


@dataclass(frozen=True)
class EarthRange:
    """A flat ring's static limit with Earth's field in every orientation: over the orientations, the least and the
    most of the two senses' geometric mean (eV), of either sense, and the largest |sqrt W- - sqrt W+| (sqrt eV)."""

    mean_low: float
    mean_high: float
    each_low: float
    each_high: float
    split: float


def ring_limit(th: Theory, sense: int, field: Arr) -> float:
    """A flat ring's static limit (eV) for one sense of travel with a uniform field B (T): the most W for which at every
    sample the bend's 2W/R and the force e v x B together stay within the strength along their direction. On the ring
    of ClosedCurve.circle u points inward and v up, so -e v x B adds s e v B_z inward and e s v (B_x cos psi + B_y sin
    psi) upward at angle psi."""
    radius = th.guide.length / (2 * math.pi)
    psi = 2 * math.pi * np.arange(len(th.amplitude)) / len(th.amplitude)
    across = field[0] * np.cos(psi) + field[1] * np.sin(psi)
    amp2 = th.amplitude**2
    lo, hi = 0.0, 4 * th.static
    for _ in range(60):
        w = (lo + hi) / 2
        v = math.sqrt(2 * w * pt.E / pt.ME)
        tilt_u = 2 * w / radius - sense * v * field[2]
        tilt_v = -sense * v * across
        cap = th.strength(np.degrees(np.arctan2(tilt_v, tilt_u))) * amp2
        lo, hi = (w, hi) if bool(np.all(np.hypot(tilt_u, tilt_v) <= cap)) else (lo, w)
    return lo


def earth_range(th: Theory, earth: float, count: int = 400) -> EarthRange:
    """ring_limit over `count` orientations of Earth's field `earth` (T), spread evenly over the sphere."""
    k = np.arange(count) + 0.5
    polar, around = np.arccos(1 - 2 * k / count), math.pi * (1 + math.sqrt(5)) * k
    ways = np.stack([np.sin(polar) * np.cos(around), np.sin(polar) * np.sin(around), np.cos(polar)], 1)
    means, each, split = [], [], []
    for way in ways:
        plus, minus = ring_limit(th, 1, earth * way), ring_limit(th, -1, earth * way)
        means.append(math.sqrt(plus * minus))
        each += [plus, minus]
        split.append(abs(math.sqrt(minus) - math.sqrt(plus)))
    return EarthRange(min(means), max(means), min(each), max(each), max(split))


@cache
def earth(design: str, cache: str | None = None) -> EarthRange:
    """earth_range for a design as drawn, with the Earth's field its tolerances draw."""
    th = theory(design, cache)
    return earth_range(th, th.guide.design.tolerance.earth)


def split_bound(span: EarthRange, ratio: float = SCAN) -> float:
    """E2-21's bound on the bench's |sqrt W- - sqrt W+| (sqrt eV): Earth's most, and SENSE_STEPS scan steps in sqrt W
    at the most either sense holds, each threshold up to a step from its truth (half a step of the scan's own and half
    of its window's)."""
    return span.split + SENSE_STEPS * math.sqrt(span.each_high) * (math.sqrt(ratio) - 1)


# ---- the shots ----------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Shot:
    """One run: its test, the design (a file relative to the repository) by a short label, and its realisation (None:
    the bench's guide, benches/twin_r1.toml's); what is measured; the energies (eV) scanned from `low` to `high` a
    factor `ratio` apart, along the path for a hold, a failure and its control and across it for a well; the RF phases
    launched; how long (laps of the slowest for a hold, RF periods otherwise); the drive (None: the design's; or `depth`
    eV deep, or the null's `q` at the mean volts); and the seed."""

    test: str
    label: str
    design: str
    realisation: int | None
    kind: Kind
    low: float = 0.0
    high: float = 0.0
    ratio: float = SCAN
    phases: int = PHASES
    laps: float = LAPS
    periods: int = 0
    depth: float | None = None
    q: float | None = None
    seed: int = 0

    @property
    def ident(self) -> str:
        where = "on the bench" if self.realisation is None else f"realisation {self.realisation}"
        drive = "" if self.depth is None else f" at a {self.depth:g} eV well"
        drive += "" if self.q is None else f" at q {self.q:g}"
        what = "" if self.kind == "fail" else f" {self.kind}"
        return f"{self.test}{what}, {self.label}, {where}{drive}"

    @property
    def energies(self) -> tuple[float, ...]:
        """The energies scanned (eV), from low up by `ratio` while within high."""
        if not self.low > 0:
            return ()
        count = math.floor(math.log(self.high / self.low) / math.log(self.ratio) + 1e-9) + 1
        return tuple(sig(self.low * self.ratio**k, 8) for k in range(count))

    @property
    def launch(self) -> int:
        """The seed its starts and RF phases are drawn from: the master seed for the must-fail shots and their
        control, which share one launch, and its own seed otherwise."""
        return MASTER if self.kind in ("fail", "control") else self.seed


@cache
def shots() -> tuple[Shot, ...]:
    """Every shot, in the protocol's order (section 4)."""
    drawn = [
        Shot("E2-21", "R1", R1, 0, "well", 0.0045, 0.14, periods=WELL_PERIODS),
        Shot("E2-21", "R1", R1, None, "well", 0.0045, 0.14, periods=WELL_PERIODS),
        Shot("E2-21", "R1", R1, 0, "well", 0.006, 0.05, periods=WELL_PERIODS, q=CONTROL_Q),
        Shot("E2-21", "R1", R1, 0, "hold", 1.0, 13.0, phases=R1_PHASES),
        Shot("E2-21", "R1", R1, None, "hold", 1.0, 13.0, phases=R1_PHASES),
        Shot("E2-21", "R1", R1, 0, "hold", 0.8, 5.0, phases=R1_PHASES, q=CONTROL_Q),
        Shot("E2-22", "T915-S", T915S, 0, "hold", 3.5, 16.0),
        Shot("E2-22", "flat ring 318", FLAT, 0, "hold", 9.5, 44.0),
        Shot("E2-22", "helix (1, 3)", HELIX, 0, "hold", 5.5, 25.5),
        Shot("slot", "straight slot guide", SLOT, 0, "beta", phases=1, periods=BETA_PERIODS),
        Shot("slot", "straight slot guide", SLOT, 0, "well", 0.16, 0.33, periods=WELL_PERIODS),
        Shot("must fail", "R1", R1, None, "fail", 0.5, 4.0, ratio=2.0, phases=4, periods=FAIL_PERIODS, depth=1.0),
        Shot("must fail", "R1", R1, None, "fail", 0.5, 4.0, ratio=2.0, phases=4, periods=FAIL_PERIODS, q=1.89),
        Shot("must fail", "R1", R1, None, "control", 0.5, 4.0, ratio=2.0, phases=4, periods=FAIL_PERIODS),
    ]
    return tuple(replace(shot, seed=MASTER + k) for k, shot in enumerate(drawn))


REFERENCE = {"E2-22": T915S}  # the design the test's ratios are taken against


def shot_volts(guide: TwinGuide, shot: Shot) -> float | None:
    """A shot's drive (V, mean): None for the design's own; else the volts of a well `depth` eV deep at the band, or of
    the null's `q` there, on the solved section."""
    solved = guide.solved
    if shot.depth is not None:
        return math.sqrt(shot.depth / solved.depth_ev(1.0))
    if shot.q is not None:
        return shot.q / solved.mathieu_q(1.0)
    return None


def run_steps(shot: Shot, length: float, dt: float) -> int:
    """The steps a scan runs: a hold until its slowest electron, at its launch speed, has gone its laps 2 % over (on a
    path `length` m round, steps of dt s); anything else its RF periods."""
    if shot.kind == "hold":
        return math.ceil(shot.laps * 1.02 * length / (_speed(min(shot.energies)) * dt))
    return shot.periods * STEPS


# ---- the theory's side --------------------------------------------------------------------------------------------
@dataclass(frozen=True, eq=False)
class Prediction:
    """What theory says a shot measures: the claim in words, the number (eV for a well's first loss or depth and a
    hold's limit, Mathieu's exponent for the secular frequency, the count held for a failure and lost for its control),
    the band the measurement must fall in, its unit, and the numbers it rests on; for E2-22's flat ring and helix, the
    ratio of their static limits to T915-S's and its band."""

    claim: str
    value: float
    low: float
    high: float
    unit: str
    basis: Mapping[str, float]
    ratio: float | None = None
    ratio_low: float | None = None
    ratio_high: float | None = None


def _depth_at(th: Theory, volts: float) -> float:
    """The well's depth (eV) at `volts` and the tuned band: the ponderomotive energy of the escape saddle's field."""
    return pt.ponderomotive_ev(th.guide.solved.escape * volts, th.driven.f)


def _common(th: Theory) -> dict[str, float]:
    solved, driven = th.guide.solved, th.driven
    out = {
        "f_MHz": driven.f / 1e6,
        "volts_mean": th.volts,
        "volts_at_start": abs(th.phasor),
        "q_null": th.q,
        "depth_mean_eV": _depth_at(th, th.volts),
        "depth_weakest_eV": driven.depth_ev * th.scale**2,
        "null_mm": solved.distance * 1e3,
        "n_eff": solved.index,
        "Q": solved.quality,
        "length_mm": th.guide.length * 1e3,
    }
    return {key: sig(value) for key, value in out.items()}


@cache
def section_lost(
    design: str, cache: str | None, volts: float, energies: tuple[float, ...], periods: int, shifts: tuple[float, ...]
) -> BoolArr:
    """section_flights for a design as drawn, driven to `volts` mean, at its start station's amplitude, aimed from the
    solved null at the escape saddle, for a shot's energies, its two zero crossings and `periods`, at the tuned band
    moved by each relative shift."""
    th = theory(design, cache).at(volts)
    solved = th.guide.solved
    du, dv = solved.saddle[0] - solved.null[0], solved.saddle[1] - solved.null[1]
    size = math.hypot(du, dv)
    omegas = 2 * math.pi * th.driven.f * (1 + np.asarray(shifts, dtype=np.float64))
    grid = th.guide.grid()
    return section_flights(grid, abs(th.phasor), omegas, energies, WELL_PHASES, periods, (du / size, dv / size))


@dataclass(frozen=True)
class Escape:
    """What the section model reads off a well's scan: the first loss (eV) at each shift of the band, and at no shift
    the robust threshold and the edges' threshold (eV)."""

    firsts: tuple[float, ...]
    robust: float
    edge: float


def escape(shot: Shot, th: Theory, cache: str | None, shifts: tuple[float, ...] = (0.0,)) -> Escape:
    """The section model on a well shot's own energies and phases, at the theory's volts and each shift of the band."""
    lost = section_lost(shot.design, cache, th.volts, shot.energies, shot.periods, shifts)
    firsts = tuple(lowest_lost(shot.energies, lost[k].any(axis=0).tolist()) for k in range(len(shifts)))
    gone = lost[shifts.index(0.0)].sum(axis=0)
    kept = (len(WELL_PHASES) - gone).tolist()
    return Escape(
        firsts, robust_threshold(shot.energies, kept, gone.tolist()), edge_threshold(shot.energies, kept, gone.tolist())
    )


def qcap_below(th: Theory, reach: float = 1.0) -> tuple[float, float]:
    """How far under the null (in h, along v toward the board) the local Mathieu q at the start station's volts first
    reaches 0.908, and the pseudopotential there over the depth D at the start; inf and 0 if not within `reach` h."""
    grid, h = th.guide.grid(), th.guide.solved.distance
    u0, v0 = grid.null
    ys = v0 - np.linspace(0.0, reach, 4001) * h
    amplitude, omega = abs(th.phasor), 2 * math.pi * th.driven.f
    q = 2 * QM * amplitude * stiffness_at(grid, np.full(len(ys), u0), ys) / omega**2
    past = np.flatnonzero(q >= QCAP)
    if not len(past):
        return math.inf, 0.0
    du, dv, _ = grid.gradient_at(np.full(1, u0), ys[past[:1]])
    there = pt.ponderomotive_ev(amplitude * math.hypot(float(du[0]), float(dv[0])), th.driven.f)
    return float((v0 - ys[past[0]]) / h), there / _depth_at(th, amplitude)


def predict_well(shot: Shot, th: Theory, cache: str | None) -> Prediction:
    depth = _depth_at(th, abs(th.phasor))
    name = th.guide.design.name
    basis = _common(th) | {"saddle_mm": sig(th.guide.solved.saddle[1] * 1e3), "depth_start_eV": sig(depth)}
    if shot.test != "E2-21":  # the slot guide's depth, against the band registered for it
        found = escape(shot, th, cache)
        basis |= {"model_first_D": sig(found.firsts[0] / depth), "model_edge_D": sig(found.edge / depth)}
        claim = (
            f"electrons sent from {name}'s null toward its escape saddle at the drive's zero crossings, with nothing "
            f"along the path, escape at its depth there, {depth:.4g} eV (the section model: {found.edge / depth:.2f} D)"
        )
        return Prediction(claim, sig(depth), sig(WELL_BAND[0] * depth), sig(WELL_BAND[1] * depth), "eV", basis)
    sd = laminate_sd(th.guide) if shot.realisation is None else 0.0
    half = DRAWS // 2 if sd else 0
    shifts = tuple(sig(LAMINATE_SDS * sd * k / max(half, 1), 8) for k in range(-half, half + 1))
    found = escape(shot, th, cache, shifts)
    first = found.firsts[shifts.index(0.0)]
    room = SCAN ** (WELL_STEPS + 0.5)
    low, high = min(found.firsts) / room, max(found.firsts) * room
    below, pseudo = qcap_below(th)
    basis |= {
        "model_first_eV": sig(first),
        "model_first_D": sig(first / depth),
        "model_robust_eV": sig(found.robust),
        "model_robust_D": sig(found.robust / depth),
        "qcap_below_null_h": sig(below, 4),
        "pseudo_at_qcap_D": sig(pseudo, 4),
    }
    claim = (
        f"electrons sent from {name}'s null up its axis toward the escape saddle at the drive's zero crossings are "
        f"first lost at {first:.4g} eV, {first / depth:.2f} of the depth D = {depth:.4g} eV, as the exact motion on "
        f"the solved axis loses them in {shot.periods} RF periods"
    )
    if sd:
        basis |= {
            "model_first_low_eV": sig(min(found.firsts)),
            "model_first_high_eV": sig(max(found.firsts)),
            "laminate_sd": sig(sd),
        }
        claim += f", {min(found.firsts):.4g} to {max(found.firsts):.4g} eV over the laminate's draw of the band"
    return Prediction(claim, sig(first), sig(low), sig(high), "eV", basis)


def _hold_basis(th: Theory) -> dict[str, float]:
    solved, guide = th.guide.solved, th.guide
    h = solved.distance
    depth = _depth_at(th, th.volts)
    out = _common(th) | {
        "static_eV": sig(th.static),
        "static_at": sig(th.where / len(th.limits), 4),
        # the strength in units of the well's depth over h: across the boards' plane (a tilt along u), pushed toward
        # the board below (along -v) and away from it (along +v), as a bend or another force pushes the electron
        "strength_across_D_per_h": sig(float(th.strength(np.zeros(1))[0]) * h / depth),
        "strength_down_D_per_h": sig(float(th.strength(np.full(1, 90.0))[0]) * h / depth),
        "strength_up_D_per_h": sig(float(th.strength(np.full(1, 270.0))[0]) * h / depth),
    }
    curve = guide.curve
    if curve is not None:
        out["bend_mm"] = sig(curve.bend_radius * 1e3)
        out["bend_rule_eV"] = sig(depth * curve.bend_radius / (2 * h))
        if guide.framing is not None:
            out["inplane_radius_mm"] = sig(guide.framing.inplane_radius * 1e3)
    return out


def predict_hold(shot: Shot, th: Theory, cache: str | None) -> Prediction:
    static = th.static
    fraction = CONTROL_BAND if shot.q is not None else DISCOUNT
    low, high = fraction[0] * static, fraction[1] * static
    basis = _hold_basis(th)
    name = th.guide.design.name
    if shot.q is not None:
        claim = (
            f"at q {shot.q:g} the static balance holds: {name} holds electrons for {LAPS:g} laps up to {fraction[0]} "
            f"to {fraction[1]} of its static limit, {static:.4g} eV, the geometric mean of the two senses' limits"
        )
    else:
        claim = (
            f"{name} holds electrons for {LAPS:g} laps up to {fraction[0]} to {fraction[1]} of its static limit, "
            f"{static:.4g} eV, the geometric mean of the two senses' limits"
        )
        if shot.test == "E2-21":
            claim = "Hypothesis: " + claim + ", four rods' discount carried to rails"
    if shot.realisation is None:
        field = th.guide.design.tolerance.earth
        span = earth(shot.design, cache)
        allowance = (1 + LAMINATE_SDS * laminate_sd(th.guide)) ** 2
        low, high = fraction[0] * span.mean_low / allowance, fraction[1] * span.mean_high * allowance
        basis |= {
            "earth_uT": sig(field * 1e6),
            "earth_mean_low_eV": sig(span.mean_low),
            "earth_mean_high_eV": sig(span.mean_high),
            "earth_each_low_eV": sig(span.each_low),
            "earth_each_high_eV": sig(span.each_high),
            "split_bound_sqrt_eV": sig(split_bound(span, shot.ratio)),
            "laminate_sd": sig(laminate_sd(th.guide)),
        }
        claim += f", with Earth's {field * 1e6:g} uT in an orientation the prediction does not know"
    test = REFERENCE.get(shot.test)
    if test is None or test == shot.design:
        return Prediction(claim, sig(static), sig(low), sig(high), "eV", basis)
    reference = theory(test, cache)
    ratio = static / reference.static
    claim += f"; {ratio:.3f} times T915-S's limit, the static limits' ratio"
    return Prediction(
        claim,
        sig(static),
        sig(low),
        sig(high),
        "eV",
        basis,
        sig(ratio),
        sig(ratio * (1 - RATIO_TOL)),
        sig(ratio * (1 + RATIO_TOL)),
    )


def predict_beta(shot: Shot, th: Theory) -> Prediction:
    q = pt.mathieu_q(th.guide.solved.gradient * abs(th.phasor), th.driven.f)
    beta = mathieu_beta(q)
    claim = (
        f"electrons {BETA_OFFSET:g}h off the null, at rest on the path, oscillate at Mathieu's exact exponent for the "
        f"2-D section's q there, {q:.4f}: beta {beta:.6f}"
    )
    basis = _common(th) | {"q_start": sig(q), "first_order_beta": sig(q / math.sqrt(2))}
    return Prediction(claim, sig(beta), sig(beta * (1 - BETA_TOL)), sig(beta * (1 + BETA_TOL)), "beta", basis)


def fail_volts(guide: TwinGuide, shot: Shot) -> float:
    """The must-fail drive (V, mean): the volts of a well `depth` eV deep at the band, or of the null's `q` there."""
    volts = shot_volts(guide, shot)
    if volts is None:
        raise ValueError(f"{shot.ident}: a must-fail shot takes a depth or a q")
    return volts


def predict_fail(shot: Shot, th: Theory) -> Prediction:
    volts = fail_volts(th.guide, shot)
    q = th.guide.solved.mathieu_q(volts)
    claim = (
        f"driven to {volts:.1f} V the null's q is {q:.2f}, past Mathieu's {QCAP}: no well forms, and no electron "
        "launched near the null is held"
    )
    basis = {"volts": sig(volts), "q_null": sig(q), "depth_eV": sig(th.guide.solved.depth_ev(volts))}
    return Prediction(claim, 0.0, 0.0, 0.0, "held", basis)


def predict_control(shot: Shot, th: Theory) -> Prediction:
    top = max(shot.energies)
    claim = (
        f"the must-fail launch at the design's {th.volts:g} V (q {th.q:.2f}): every electron is held for "
        f"{shot.periods} RF periods, {top:g} eV being {top / th.static:.2f} of the static limit, under the "
        f"{DISCOUNT[0]} a hold keeps"
    )
    basis = {
        "volts": sig(th.volts),
        "q_null": sig(th.q),
        "static_eV": sig(th.static),
        "top_over_static": sig(top / th.static),
    }
    return Prediction(claim, 0.0, 0.0, 0.0, "lost", basis)


def predict(shot: Shot, cache: str | None = None) -> Prediction:
    th = theory(shot.design, cache)
    if shot.kind == "fail":
        return predict_fail(shot, th)
    if shot.kind == "control":
        return predict_control(shot, th)
    volts = shot_volts(th.guide, shot)
    th = th if volts is None else th.at(volts)
    if shot.kind == "well":
        return predict_well(shot, th, cache)
    if shot.kind == "hold":
        return predict_hold(shot, th, cache)
    return predict_beta(shot, th)


# ---- the tracking's side ------------------------------------------------------------------------------------------
@dataclass(frozen=True, eq=False)
class Scan:
    """What a scan of energies found: for each sense of travel (+1 with the drive's wave, -1 against it; one row, 0,
    for a well) and each energy (eV, along the path or across it), how many electrons were held, lost to the guide,
    stopped by the gas, and (for a hold) still going but short of their laps, over the RF phases launched (the drive's
    phase at the start station as each batch began, rad), and the fewest laps (or RF periods, for anything but a hold)
    after which one was lost; the run (steps of dt s), the drive (Hz, V mean, the null's q), the guide's identity and
    the seed."""

    kind: str
    guide: Mapping[str, Any]
    seed: int
    energies: tuple[float, ...]
    senses: tuple[int, ...]
    phases: tuple[float, ...]
    held: tuple[tuple[int, ...], ...]
    lost: tuple[tuple[int, ...], ...]
    stopped: tuple[tuple[int, ...], ...]
    short: tuple[tuple[int, ...], ...]
    first_loss: tuple[tuple[float, ...], ...]
    steps: int
    dt: float
    f: float
    volts: float
    q: float

    def edges(self, row: int) -> tuple[float, float]:
        """The highest energy held below the lowest lost, and the lowest lost (eV) in one row (edges())."""
        return edges(self.energies, self.held[row], self.lost[row])

    def lowest_lost(self, row: int) -> float:
        """The first loss: the lowest energy at which an electron was lost in one row (inf if none was)."""
        return lowest_lost(self.energies, self.lost[row])

    def edge_threshold(self, row: int) -> float:
        """The geometric mean of the two edges in one row: 0 below the scan and inf above it. The slot guide's depth
        reads it."""
        return edge_threshold(self.energies, self.held[row], self.lost[row])

    def threshold(self, row: int) -> float:
        """The energy at which electrons start to be lost in one row, by robust_threshold's window: 0 below the scan
        and inf above it."""
        return robust_threshold(self.energies, self.held[row], self.lost[row])

    def held_total(self) -> int:
        return int(sum(sum(row) for row in self.held))

    def lost_total(self) -> int:
        return int(sum(sum(row) for row in self.lost))

    def launched(self) -> int:
        return int(sum(sum(row) for table in (self.held, self.lost, self.stopped, self.short) for row in table))


@dataclass(frozen=True, eq=False)
class Beat:
    """The secular frequency as electrons near the null show it: their starts (u, v, in h), the exponent each shows
    read once an RF period, how many stayed, the run (steps of dt s), the drive (Hz, |V| at the station), the null's q
    there by the device's own section, the guide's identity and the seed."""

    kind: str
    guide: Mapping[str, Any]
    seed: int
    offsets: tuple[tuple[float, float], ...]
    beta: tuple[float, ...]
    held: int
    steps: int
    dt: float
    f: float
    volts: float
    q: float


Measurement = Scan | Beat


@dataclass(frozen=True, eq=False)
class Station:
    """Where a shot starts: the path's first sample, its tangent, the section's first and second normals there (u and
    v), the curvature along them (1/m), and the drive's phasor (V)."""

    point: Arr
    tangent: Arr
    first: Arr
    second: Arr
    kappa: tuple[float, float]
    phasor: complex


def station(guide: TwinGuide, field: TwinField) -> Station:
    swept = field.swept
    ku, kv = section_curvature(guide, field)
    phasor = complex(_wave(field).phasor(np.zeros(1))[0])
    return Station(swept.points[0], swept.tangents[0], swept.first[0], swept.second[0], (ku[0], kv[0]), phasor)


def launch(
    field: TwinField, start: Station, offsets: Arr, velocity: Arr, t0: float, omega: float, dt: float = 0.0
) -> tuple[Arr, Arr]:
    """Electrons at `offsets` (N, 2) (m from the null in the section) with velocities (N, 3) (m/s), in step with the
    drive's micromotion there at time t0: the displacement -(|q/m| / w^2) Re(V e^(i w t0)) grad phi and the velocity
    (|q/m| / w) Im(V e^(i w t)) grad phi added, the velocity half a step dt before t0 (integrate's Boris step takes it
    so, leapfrog), so that each starts on its slow motion. Returns x and gamma v."""
    grid = field.swept.grid
    gu, gv, _ = grid.gradient_at(grid.null[0] + offsets[:, 0], grid.null[1] + offsets[:, 1])
    gradient = gu[:, None] * start.first + gv[:, None] * start.second
    volts = start.phasor * complex(math.cos(omega * t0), math.sin(omega * t0))
    before = start.phasor * complex(math.cos(omega * (t0 - dt / 2)), math.sin(omega * (t0 - dt / 2)))
    x0 = start.point + offsets[:, :1] * start.first + offsets[:, 1:] * start.second
    x0 = x0 - (QM / omega**2) * volts.real * gradient
    v = velocity + (QM / omega) * before.imag * gradient
    gamma = 1 / np.sqrt(1 - (v * v).sum(1) / pt.C0**2)
    u0: Arr = gamma[:, None] * v
    return x0, u0


def _speed(energy: float) -> float:
    return float(pt.v_of(energy))


def _phases(shot: Shot) -> tuple[float, ...]:
    """The RF phases (rad) a shot launches at: for a well the drive's two zero crossings, where an electron's launch
    velocity is its slow velocity; otherwise one in each of `phases` equal parts of the period, drawn from the shot's
    launch seed."""
    if shot.kind == "well":
        return WELL_PHASES
    rng = np.random.default_rng(np.random.SeedSequence([shot.launch, 1]))
    return tuple(float(2 * math.pi * (k + rng.uniform()) / shot.phases) for k in range(shot.phases))


def _shifted(field: TwinField, t0: float) -> pt.Fields:
    """The field from time t0 on: a batch launched at its own RF phase."""

    def fields(x: Arr, t: float) -> tuple[Arr | None, Arr | None]:
        return field(x, t + t0)

    return fields


def _watch(field: TwinField, detector: Any, count: int, start: Arr) -> Watch:
    """The detector's count of held, lost and stopped electrons, without its pickups: no criterion reads their
    currents, which a step-by-step record would hold for every one of a long run's steps."""
    return Watch(field, count, (), -pt.E, float(detector.weight), float(detector.collector), 1, start)


def _guide(shot: Shot, devices: Mapping[str, Any]) -> TwinGuide:
    """The shot's guide: the bench's own for realisation None, else the design opened at its realisation beside it,
    sharing the bench guide's cache of solves."""
    bench = devices["guide"]
    if not isinstance(bench, TwinGuide):
        raise lab.Refusal(f"the bench's guide is a {type(bench).__name__}, and this experiment runs on array16.twin")
    if shot.realisation is None:
        named = Design.load(_path(shot.design)).name
        if bench.design.name != named:
            raise lab.Refusal(f"{shot.ident} runs on {named}, and the bench's guide is {bench.design.name}")
        return bench
    return TwinGuide(shot.design, shot.realisation, cache=bench.cache)


def _counts(table: Sequence[Sequence[Sequence[bool]]]) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(int(sum(cell)) for cell in row) for row in table)


def measure_scan(shot: Shot, devices: Mapping[str, Any]) -> Scan:
    """A hold, a well, a failure or its control: electrons launched at every energy (and sense), at each RF phase,
    each phase one batch of integrate with the bench's gas and the detector's watch."""
    guide = _guide(shot, devices)
    volts = shot_volts(guide, shot)
    driven = guide.drive() if volts is None else guide.drive(volts=volts)
    field = guide.field(driven)
    start = station(guide, field)
    omega, period = 2 * math.pi * driven.f, 1 / driven.f
    dt = period / STEPS
    energies = shot.energies
    senses = (0,) if shot.kind == "well" else (1, -1)
    phases = _phases(shot)
    rng = np.random.default_rng(np.random.SeedSequence([shot.launch, 2]))
    h = guide.solved.distance
    grid = field.swept.grid
    count = len(energies) * len(senses)
    steps = run_steps(shot, guide.length, dt)
    if shot.kind == "hold":
        kappa = math.hypot(*start.kappa)
        bend = valley(well_patch(grid, driven.f, h), math.atan2(start.kappa[1], start.kappa[0]))
        offsets = np.array(
            [bend.at(2 * e * kappa, abs(start.phasor)) for _ in senses for e in energies], dtype=np.float64
        ).reshape(count, 2)
    elif shot.kind == "well":
        aim = np.array(guide.solved.saddle) - np.array(guide.solved.null)
        aim /= np.linalg.norm(aim)
        offsets = np.zeros((count, 2))
    else:
        radius = FAIL_OFFSET * h * np.sqrt(rng.uniform(size=count))
        turn = rng.uniform(0, 2 * math.pi, size=count)
        offsets = np.column_stack([radius * np.cos(turn), radius * np.sin(turn)])
    vacuum, detector = devices["vacuum"], devices["detector"]
    held: list[list[list[bool]]] = [[[] for _ in energies] for _ in senses]
    lost: list[list[list[bool]]] = [[[] for _ in energies] for _ in senses]
    stopped: list[list[list[bool]]] = [[[] for _ in energies] for _ in senses]
    short: list[list[list[bool]]] = [[[] for _ in energies] for _ in senses]
    first: list[list[float]] = [[math.inf for _ in energies] for _ in senses]
    for k, phase in enumerate(phases):
        t0 = float(np.mod((phase - math.atan2(start.phasor.imag, start.phasor.real)) / omega, period))
        velocity = np.zeros((count, 3))
        for i, (s, e) in enumerate((s, e) for s in senses for e in energies):
            if shot.kind == "well":
                velocity[i] = _speed(e) * (aim[0] * start.first + aim[1] * start.second)
            else:
                velocity[i] = s * _speed(e) * start.tangent
        if shot.kind in ("fail", "control"):
            x0 = start.point + offsets[:, :1] * start.first + offsets[:, 1:] * start.second
            gamma = 1 / np.sqrt(1 - (velocity * velocity).sum(1) / pt.C0**2)
            u0: Arr = gamma[:, None] * velocity
        else:
            x0, u0 = launch(field, start, offsets, velocity, t0, omega, dt)
        watch = _watch(field, detector, count, x0)
        gas = vacuum.scatter(seed=shot.seed * 16 + k)
        pt.integrate(x0, u0, np.full(count, pt.QM_E), dt, steps, _shifted(field, t0), [gas, watch])
        seen = watch.reading()
        laps = np.abs(seen.laps)
        if shot.kind == "hold":  # held only once round its laps, whatever came after
            done = laps >= shot.laps
            keep, running = done, (seen.fate == 0) & ~done
            gone, gassed = (seen.fate == 1) & ~done, (seen.fate == 2) & ~done
            when = laps
        else:
            keep, gone, gassed = seen.fate == 0, seen.fate == 1, seen.fate == 2
            running = np.zeros(count, dtype=np.bool_)
            when = seen.time / period
        for i in range(count):
            row, col = divmod(i, len(energies))
            held[row][col].append(bool(keep[i]))
            lost[row][col].append(bool(gone[i]))
            stopped[row][col].append(bool(gassed[i]))
            short[row][col].append(bool(running[i]))
            if gone[i]:
                first[row][col] = min(first[row][col], float(when[i]))
    return Scan(
        kind=shot.kind,
        guide=guide.identity(),
        seed=shot.seed,
        energies=energies,
        senses=senses,
        phases=phases,
        held=_counts(held),
        lost=_counts(lost),
        stopped=_counts(stopped),
        short=_counts(short),
        first_loss=tuple(tuple(row) for row in first),
        steps=steps,
        dt=dt,
        f=driven.f,
        volts=driven.volts,
        q=driven.q,
    )


def beta_of(path: Arr) -> Arr:
    """Mathieu's exponent from positions read once an RF period ((K, N, 2): K reads of N electrons): for a linear
    oscillation x(k + 1) + x(k - 1) = 2 cos(mu) x(k), cos mu a least-squares ratio over both coordinates, beta =
    mu / pi."""
    inner, outer = path[1:-1], path[2:] + path[:-2]
    cos_mu = (inner * outer).sum(axis=(0, 2)) / (2 * (inner * inner).sum(axis=(0, 2)))
    out: Arr = np.arccos(np.clip(cos_mu, -1.0, 1.0)) / math.pi
    return out


def measure_beta(shot: Shot, devices: Mapping[str, Any]) -> Beat:
    """Four electrons BETA_OFFSET h off the null, on the section's diagonals and axes, at rest on the path, read once
    an RF period, with the detector's watch and no gas: at 1e-6 mbar the gas would stop one in 10^5 in 60 RF periods,
    and no draw should decide which of four is read."""
    guide = _guide(shot, devices)
    driven = guide.drive()
    field = guide.field(driven)
    start = station(guide, field)
    h = guide.solved.distance
    root = 1 / math.sqrt(2)
    starts = np.array([[root, root], [root, -root], [1.0, 0.0], [0.0, 1.0]])
    offsets = BETA_OFFSET * h * starts
    x0 = start.point + offsets[:, :1] * start.first + offsets[:, 1:] * start.second
    dt = 1 / (driven.f * STEPS)
    reads = [offsets.copy()]

    def strobe(step: int, t: float, x: Arr, u: Arr, alive: BoolArr) -> None:
        if (step + 1) % STEPS == 0:
            d = x - start.point
            reads.append(np.column_stack([d @ start.first, d @ start.second]))

    watch = _watch(field, devices["detector"], len(x0), x0)
    steps = shot.periods * STEPS
    flight = pt.integrate(x0, np.zeros_like(x0), np.full(len(x0), pt.QM_E), dt, steps, field, [watch, strobe])
    beta = beta_of(np.array(reads))
    volts = abs(start.phasor)
    return Beat(
        kind="beta",
        guide=guide.identity(),
        seed=shot.seed,
        offsets=tuple((float(a), float(b)) for a, b in starts * BETA_OFFSET),
        beta=tuple(float(b) for b in beta),
        held=int(flight.alive.sum()),
        steps=steps,
        dt=dt,
        f=driven.f,
        volts=volts,
        q=pt.mathieu_q(guide.solved.gradient * volts, driven.f),
    )


def measure(shot: Shot, devices: Mapping[str, Any]) -> Measurement:
    return measure_beta(shot, devices) if shot.kind == "beta" else measure_scan(shot, devices)


def load(record: Mapping[str, Any]) -> Measurement:
    """A measurement back from the plain data it was stored as."""

    def rows(name: str) -> tuple[tuple[Any, ...], ...]:
        return tuple(tuple(row) for row in record[name])

    if record["kind"] == "beta":
        return Beat(
            kind="beta",
            guide=record["guide"],
            seed=int(record["seed"]),
            offsets=tuple((float(a), float(b)) for a, b in record["offsets"]),
            beta=tuple(float(b) for b in record["beta"]),
            held=int(record["held"]),
            steps=int(record["steps"]),
            dt=float(record["dt"]),
            f=float(record["f"]),
            volts=float(record["volts"]),
            q=float(record["q"]),
        )
    return Scan(
        kind=str(record["kind"]),
        guide=record["guide"],
        seed=int(record["seed"]),
        energies=tuple(float(e) for e in record["energies"]),
        senses=tuple(int(s) for s in record["senses"]),
        phases=tuple(float(p) for p in record["phases"]),
        held=rows("held"),
        lost=rows("lost"),
        stopped=rows("stopped"),
        short=rows("short"),
        first_loss=tuple(tuple(math.inf if x is None else float(x) for x in row) for row in record["first_loss"]),
        steps=int(record["steps"]),
        dt=float(record["dt"]),
        f=float(record["f"]),
        volts=float(record["volts"]),
        q=float(record["q"]),
    )


# ---- the verdicts: each function is one criterion of PROTOCOL.md, section 5 ---------------------------------------
Row = tuple[Shot, Prediction, Measurement | None]


def _ev(value: float) -> str:
    if value == 0.0:
        return "below the scan"
    return "above the scan" if math.isinf(value) else f"{value:.4g} eV"


def _within(value: float, low: float, high: float) -> bool:
    return math.isfinite(value) and low <= value <= high


def mean_limit(scan: Scan) -> float:
    """The two senses' thresholds' geometric mean (eV): 0 if either is below the scan, inf if either is above it."""
    a, b = (scan.threshold(row) for row in range(len(scan.senses)))
    if a == 0.0 or b == 0.0:
        return 0.0
    return math.inf if math.isinf(a) or math.isinf(b) else math.sqrt(a * b)


def check_wells(rows: Iterable[Row]) -> Iterator[Check]:
    """E2-21's wells: the first loss within the axis model's band; the slot guide's depth: the escape threshold within
    WELL_BAND of the predicted depth."""
    for shot, p, m in rows:
        if shot.kind != "well" or not isinstance(m, Scan):
            continue
        if shot.test == "E2-21":
            found = m.lowest_lost(0)
            name = "E2-21 well" if shot.q is None else "E2-21 well control"
            said = (
                f"first lost at {_ev(found)} against the axis model's {p.value:.4g} eV (band {p.low:.4g} to "
                f"{p.high:.4g}); the robust threshold {_ev(m.threshold(0))}"
            )
        else:
            found = m.edge_threshold(0)
            name = "slot depth"
            said = f"escape at {_ev(found)} against the depth {p.value:.4g} eV (band {p.low:.4g} to {p.high:.4g})"
        yield Check(name, shot.ident, _within(found, p.low, p.high), said)


def check_holds(rows: Iterable[Row]) -> Iterator[Check]:
    """E2-21's and E2-22's held limits: the senses' geometric mean within the band; E2-21's senses as theory says."""
    for shot, p, m in rows:
        if shot.kind != "hold" or not isinstance(m, Scan):
            continue
        found = mean_limit(m)
        plus, minus = m.threshold(0), m.threshold(1)
        first = min(m.lowest_lost(0), m.lowest_lost(1))
        if shot.test == "E2-21":
            name = "E2-21 hold" if shot.q is None else "E2-21 hold control"
        else:
            name = "E2-22 static limit"
        said = (
            f"held to {_ev(found)} (with the wave {_ev(plus)}, against it {_ev(minus)}; the first loss at "
            f"{_ev(first)}) against the static limit {p.value:.4g} eV (band {p.low:.4g} to {p.high:.4g})"
        )
        yield Check(name, shot.ident, _within(found, p.low, p.high), said)
        if shot.test != "E2-21":
            continue
        ok = 0 < plus < math.inf and 0 < minus < math.inf
        if shot.realisation is None:
            bound = p.basis["split_bound_sqrt_eV"]
            split = abs(math.sqrt(minus) - math.sqrt(plus)) if ok else math.inf
            said = f"|sqrt W- - sqrt W+| = {split:.3f} sqrt(eV) against Earth's most and two steps, {bound:.3f}"
            yield Check("E2-21 senses", shot.ident, ok and split <= bound, said)
        else:
            apart = abs(math.log(plus / minus)) if ok else math.inf
            limit = SENSE_STEPS * math.log(shot.ratio)
            said = f"the senses' limits {apart / math.log(shot.ratio):.1f} scan steps apart, {SENSE_STEPS} allowed"
            yield Check("E2-21 senses", shot.ident, ok and apart <= limit + 1e-12, said)


def check_ratios(rows: Sequence[Row]) -> Iterator[Check]:
    """E2-22: the flat ring's and the helix's limits over T915-S's within RATIO_TOL of their static limits' ratio."""
    by_design = {(shot.test, shot.design): m for shot, _, m in rows if shot.kind == "hold"}
    for shot, p, m in rows:
        if shot.kind != "hold" or p.ratio is None or p.ratio_low is None or p.ratio_high is None:
            continue
        base = by_design.get((shot.test, REFERENCE[shot.test]))
        if not isinstance(m, Scan) or not isinstance(base, Scan):
            continue
        ours, theirs = mean_limit(m), mean_limit(base)
        ratio = ours / theirs if 0 < theirs < math.inf and 0 < ours < math.inf else math.nan
        said = (
            f"{_ev(ours)} over T915-S's {_ev(theirs)}: {ratio:.3f} against the static {p.ratio:.3f} (band "
            f"{p.ratio_low:.3f} to {p.ratio_high:.3f})"
        )
        yield Check("E2-22 price of the knot", shot.ident, _within(ratio, p.ratio_low, p.ratio_high), said)


def check_beta(rows: Iterable[Row]) -> Iterator[Check]:
    """The slot guide's secular frequency: every electron's exponent within BETA_TOL of Mathieu's, all of them held."""
    for shot, p, m in rows:
        if shot.kind != "beta" or not isinstance(m, Beat):
            continue
        worst = max(abs(b / p.value - 1) for b in m.beta)
        ok = m.held == len(m.beta) and worst <= BETA_TOL
        said = (
            f"beta {min(m.beta):.6f} to {max(m.beta):.6f} against Mathieu's {p.value:.6f} (at most {worst:.2e} off, "
            f"{BETA_TOL:g} allowed); {m.held} of {len(m.beta)} held"
        )
        yield Check("slot secular frequency", shot.ident, ok, said)


def check_fail(rows: Sequence[Row]) -> Iterator[Check]:
    """The control and the cases that must fail: the control's launch, at the design's volts, loses no electron to the
    guide; each must-fail shot holds none, and passes only beside a control that held."""
    control: bool | None = None
    for shot, _, m in rows:
        if shot.kind != "control" or not isinstance(m, Scan):
            continue
        held, lost = m.held_total(), m.lost_total()
        control = lost == 0 and held > 0
        said = f"q {m.q:.2f} at {m.volts:.1f} V: {held} of {m.launched()} held, {lost} lost to the guide"
        yield Check("must fail control", shot.ident, control, said)
    beside = {None: "no control on file", True: "the control held", False: "the control did not hold"}[control]
    for shot, p, m in rows:
        if shot.kind != "fail" or not isinstance(m, Scan):
            continue
        held, launched = m.held_total(), m.launched()
        losses = [x for row in m.first_loss for x in row if math.isfinite(x)]
        late = f", the last lost after {max(losses):.1f} RF periods" if losses else ""
        said = f"q {m.q:.2f} at {m.volts:.1f} V (theory {p.basis['q_null']:.2f}): {held} of {launched} held{late}"
        said += f"; {beside}"
        yield Check("must fail", shot.ident, held == 0 and launched > 0 and control is True, said)


CHECKS: tuple[Callable[[Sequence[Row]], Iterator[Check]], ...] = (
    check_wells,
    check_holds,
    check_ratios,
    check_beta,
    check_fail,
)


# ---- what the predictions file says of the apparatus -------------------------------------------------------------
def describe(design: str, cache: str | None) -> dict[str, Any]:
    """A design as the theory sees it: its file and solve, its line and well, its ring and drive as drawn, and its
    static balance."""
    th = theory(design, cache)
    guide, solved, driven = th.guide, th.guide.solved, th.driven
    out: dict[str, Any] = {
        "identity": guide.identity(),
        "line": {
            "n_eff": sig(solved.index),
            "z0_ohm": sig(solved.z0),
            "share": sig(solved.share),
            "Qc": sig(solved.qc),
            "Qd": sig(solved.qd),
            "Q": sig(solved.quality),
        },
        "well": {
            "null_mm": sig(solved.distance * 1e3),
            "gradient_per_m2": sig(solved.gradient),
            "saddle_mm": [sig(solved.saddle[0] * 1e3), sig(solved.saddle[1] * 1e3)],
            "escape_per_m": sig(solved.escape),
        },
        "ring": {"length_mm": sig(guide.length * 1e3), "f_MHz": sig(driven.f / 1e6)},
        "drive": {
            "watts": sig(driven.watts),
            "volts": sig(driven.volts),
            "least": sig(driven.least),
            "most": sig(driven.most),
            "backward": sig(driven.purity),
            "depth_ratio": sig(driven.depth_ratio),
            "weakest_eV": sig(driven.depth_ev),
            "q": sig(driven.q),
        },
    }
    if guide.curve is not None:
        out["ring"]["bend_mm"] = sig(guide.curve.bend_radius * 1e3)
        out["static"] = _hold_basis(th)
    if guide.framing is not None:
        out["ring"]["inplane_radius_mm"] = sig(guide.framing.inplane_radius * 1e3)
    return out


def _band(p: Prediction) -> str:
    if p.unit == "held":
        return "none held"
    if p.unit == "lost":
        return "none lost"
    digits = 6 if p.unit == "beta" else 4
    return f"{p.low:.{digits}g} to {p.high:.{digits}g} {p.unit}"


def predictions_markdown(experiment: TwinRung1) -> str:
    """The predictions as tables, for the protocol's reader."""
    lines = [
        "<!-- Written by `experiments/twin_rung1.py predict`. Theory only: no electron was tracked to make it. -->",
        "",
        "# Experiment 2, stage 1: the predictions",
        "",
        "Each prediction is made from the design's solved 2-D section and formulas (array16.guide.sections, the ring",
        "model of array16.guide.ring, the motion on the section's axis and the static balance of",
        "experiments/twin_rung1.py), at the design as drawn. E2-21 is a statement about R1's rails' section swept",
        "round a ring with an ideal in-phase drive: the single-mode ring drives both rails alike, which full wave does",
        "not find R1's ring to do (PROTOCOL.md, section 7).",
        "",
    ]
    rows = []
    for shot in experiment.shots():
        p = experiment.predict(shot)
        ratio = "" if p.ratio is None else f"{p.ratio:.3f} ({p.ratio_low:.3f} to {p.ratio_high:.3f})"
        value = "0" if p.unit in ("held", "lost") else f"{p.value:.6g} {p.unit}"
        rows.append([shot.ident, p.claim, value, _band(p), ratio])
    lines += lab.table(["Shot", "Prediction", "Theory", "Band", "Ratio to T915-S"], rows)
    lines += ["", "## The numbers each rests on", ""]
    for shot in experiment.shots():
        p = experiment.predict(shot)
        said = ", ".join(f"{key} {value:g}" for key, value in p.basis.items())
        lines.append(f"- {shot.ident}: {said}")
    return "\n".join(lines) + "\n"


class TwinRung1(lab.Experiment[Shot, Prediction, Measurement]):
    """The calibration rung as array16.lab runs it, on the twin bench (benches/twin_r1.toml): its guide, detector and
    vacuum. `cache` is the folder of solved sections (None: array16.twin's)."""

    out = OUT
    needs: Mapping[str, Sequence[str]] = {
        "guide": ("tune", "drive", "field", "identity"),
        "detector": ("watch", "identity"),
        "vacuum": ("scatter", "free_path", "identity"),
    }
    registered = (
        OUT / "PROTOCOL.md",
        OUT / "predictions.md",
        OUT / "slot_straight.toml",
        Path(__file__).resolve(),
        lab.ROOT / R1,
        lab.ROOT / T915S,
        lab.ROOT / FLAT,
        lab.ROOT / HELIX,
        lab.BENCHES / f"{BENCH}.toml",
        lab.ROOT / "array16" / "twin",
        lab.ROOT / "array16" / "guide",
        lab.ROOT / "array16" / "particles",
        lab.ROOT / "array16" / "lab.py",
    )

    def __init__(self, cache: str | Path | None = None) -> None:
        self.cache = None if cache is None else str(cache)

    def shots(self) -> Sequence[Shot]:
        return shots()

    def ident(self, shot: Shot) -> str:
        return shot.ident

    def predict(self, shot: Shot) -> Prediction:
        return predict(shot, self.cache)

    def measure(self, batch: Sequence[Shot], devices: Mapping[str, Any]) -> Mapping[str, Measurement]:
        return {shot.ident: measure(shot, devices) for shot in batch}

    def load(self, record: Mapping[str, Any]) -> Measurement:
        return load(record)

    def judge(self, rows: Sequence[Row]) -> Sequence[Check]:
        return [check for criterion in CHECKS for check in criterion(rows)]

    def apparatus(self) -> Mapping[str, Any]:
        designs = {d: describe(d, self.cache) for d in (R1, T915S, FLAT, HELIX, SLOT)}
        criteria = {
            "steps_per_period": STEPS,
            "laps": LAPS,
            "scan_ratio": SCAN,
            "window": WINDOW,
            "model_steps_per_period": MODEL_STEPS,
            "discount": DISCOUNT,
            "control_band": CONTROL_BAND,
            "ratio_tolerance": RATIO_TOL,
            "well_band": WELL_BAND,
            "well_steps": WELL_STEPS,
            "beta_tolerance": BETA_TOL,
            "sense_steps": SENSE_STEPS,
            "laminate_sds": LAMINATE_SDS,
            "laminate_draws": DRAWS,
            "master_seed": MASTER,
        }
        return {"bench": f"benches/{BENCH}.toml", "designs": designs, "criteria": criteria}

    def summary(self, checks: Sequence[Check]) -> Mapping[str, Sequence[Check]]:
        by_ident = {shot.ident: shot for shot in shots()}
        tests = ("E2-21", "E2-22", "slot", "must fail")
        return {
            f"{test}, every check": [c for c in checks if c.ident in by_ident and by_ident[c.ident].test == test]
            for test in tests
        }


def main(argv: Sequence[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    experiment = TwinRung1()
    if args[:1] == ["run"]:  # predictions.md is registered beside predictions.json, so it too must be the code's
        written = OUT / "predictions.md"
        if not written.exists() or written.read_text() != predictions_markdown(experiment):
            print(
                f"refused: {written} is missing or not what the code predicts now: run `predict` first", file=sys.stderr
            )
            return 2
    code = lab.main(experiment, args, about=(__doc__ or "").split("\n", 1)[0])
    if code == 0 and args[:1] == ["predict"]:
        (OUT / "predictions.md").write_text(predictions_markdown(experiment))
    return code


if __name__ == "__main__":
    sys.exit(main())
