"""Experiment 2, stage 2a: what a knot of boards keeps of its static limit, where it loses the rest, and what Earth's
field, the bend's shift of the null and secondary emission do to it, on the twin bench (docs/KNOTTED_GUIDE.md section
10; experiments/twin_rung1/REPORT.md sections 5, 6 and 8).

Stage 1 calibrated the tracking chain and left four hypotheses read from its results. Stage 2a registers the shots that
test them, and stage 2's first predictions, on the bench of benches/twin_r2.toml (T915-S built once), every prediction
made from theory before any electron is tracked:

  H1         leaving the boards' plane costs about 14 % of the static limit, knotted or not, against H1a (the bend's
             changing along the path) and H1b (the boards' stiff side giving a bend nothing): the ellipse of 318 mm, the
             flat ring bent toward its boards, the (1, 3) helix at r/R 0.1, 0.2 and 0.3 and the flat ring, at section B
             and 50 V, and T915-S read from H2's run at 20 laps
  E2-9       facing boards hold their static limit bent in their plane and at least twice that toward them: the flat
             ring and the stiff ring at q 0.30, and at q 0.49 (H1's)
  H2         T915-S's band of late losses and its senses: 100 laps at eight phases on a 1 % scan, at q 0.49 and q 0.30
  H3         every electron a hold out of the plane loses is lost in a lobe or within 0.07 lap past one: the records
  E2-15      Earth's field across T915-S at its worst orientation and at ten times its strength, and the bench's T915-S
             with its own field and with none
  E2-2       the bend's shift of the null against the axisymmetric reference, and T915-S's hold with it on beside off
  H4         R1 on the bench's realisation with Earth's field removed and as drawn: hold and well
  E2-13      multipactor in the boards' own field with Furman and Pivi's secondary emission: T915-S at 50.3 V, the split
             lever at 39.6 V, R1's rails at 50 V, section B overdriven (the control), and the onset on facing strips
  must fail  T915-S's section driven past Mathieu's q 0.908, beside the same launch at 50 V

    uv run python experiments/twin_rung2.py predict                  # theory only: twin_rung2/predictions.{json,md}
    uv run python experiments/twin_rung2.py run --bench twin_r2      # tracking: twin_rung2/runs.jsonl
    uv run python experiments/twin_rung2.py report                   # the verdicts: twin_rung2/results.{json,md}
    uv run python experiments/twin_rung2.py estimate                 # the run's time on this machine, from the shots

The theory of a hold is the exact motion in the section along the path (exact_motion): the section's RF field at the
drive the ring gives each place, the bend's centrifugal force with the curvilinear factor 1/(1 - kappa.r), the frame's
turning (Coriolis, and its centrifugal part), the wave's magnetic force and Earth's, on the same energies, senses, RF
phases, launch and laps as the shot, by the tracker's own leapfrog at 120 steps an RF period; computed after stage 1 on
its five holds, it kept what they kept (PROTOCOL.md section 6). The static balance of stage 1 (twin_rung1.theory),
the bound with every bend counted across the boards' plane (H1b's), and the static balance with a uniform field
(field_limit, E2-15 and H4) give the hypotheses' numbers. The theory of a well is stage 1's exact motion on the axis;
of multipactor, multipactor.py's level-1 chart and its one-dimensional operator. The protocol, with what counts as
failure, is experiments/twin_rung2/PROTOCOL.md. The tracking runs each shot's RF phases in a process pool of the
bench's workers (array16.lab runs a bench with devices in one process; the twin's devices are rebuilt from their
options in each worker), and the exact motion of every hold is kept under .twin_cache/ by the hash of its inputs and
of the code it runs, so that `run` and `report` read it back.
"""

from __future__ import annotations

import hashlib
import inspect
import json
import math
import os
import sys
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, replace
from functools import cache
from pathlib import Path
from typing import Any, Literal

import numpy as np
import twin_rung1 as rung

from array16 import lab
from array16 import particles as pt
from array16.guide import axisym
from array16.guide import emission as em
from array16.guide import multipactor as mp
from array16.guide.sections import section_b
from array16.lab import Check
from array16.particles import Arr, BoolArr, IntArr
from array16.twin import CACHE, Driven, Losses, TwinField, TwinGuide, TwinVacuum, Watch, solver_version
from array16.twin.records import LEFT, STRUCK

OUT = Path(__file__).resolve().parent / "twin_rung2"
BENCH = "twin_r2"
T915S = "designs/t915s.toml"
FLAT = "designs/flat_ring_318.toml"
STIFF = "designs/flat_ring_318_stiff.toml"
ELLIPSE = "designs/ellipse_318.toml"
HELIX_010 = "designs/helix_13_r010.toml"
HELIX = "designs/helix_13.toml"
HELIX_030 = "designs/helix_13_r030.toml"
R1 = "designs/r1_rail_ring.toml"
RING_IN = "experiments/twin_rung2/ring_40h_inplane.toml"
RING_TOWARD = "experiments/twin_rung2/ring_40h_toward.toml"
MASTER = 20261123  # stage 2a's master seed: each shot's seed is it plus the shot's place in the list
QM = abs(pt.QM_E)

# ---- how it is run ------------------------------------------------------------------------------------------------
STEPS = rung.STEPS  # 120 steps an RF period, for the tracker and for the exact motion
LAPS = rung.LAPS  # 20 laps make a hold (DR-7)
LONG = 100.0  # H2's laps
FINE = 1.01  # the ratio between neighbouring energies of a fine scan
COARSE = 1.025  # stage 1's, kept where the hypotheses lie far apart (the stiff rings) and for the wells
PHASES = 8  # RF phases a hold launches: 40 electrons a sense in a threshold's window, as R1's in stage 1
FAIL_PHASES = 4
WELL_PERIODS = rung.WELL_PERIODS
FAIL_PERIODS = rung.FAIL_PERIODS
FAIL_Q = 1.2  # the must-fail drive: the null's Mathieu q at the mean volts, past 0.908
ELLIPSE_START = 0.25  # the ellipse's launch, a quarter lap round: the end of its minor axis, its gentlest bend
SHIFT_STATIONS = (0, 500, 1000, 1500)  # E2-2: the samples of the 2000-sample ring the swept null is found at
MPC_SEEDS = (4, 5, 6, 7)  # E2-13: the emission module read its cases on seeds 0 to 3 while it was built
MPC_SETTINGS: Mapping[str, int] = {  # E2-13: the emission module's own readings' settings
    "electrons": 1000,
    "periods": 30,
    "settle": 10,
    "count": 10,
    "generation_settle": 3,
    "steps": 200,
}
ONSET_DEPTHS = (1.5, 2.0, 2.5, 3.0)  # E2-13: facing strips fitted to these wells (eV) at q 0.5
OVERDRIVE = 160.0  # E2-13's control: section B at these volts, past its gap's parallel-plate threshold
EARTH = 50e-6  # T: Earth's field, as the designs draw it
STRONG = 10.0  # E2-15's control: Earth's field this many times stronger, in the same orientation
WORST = (0.234856, -0.116608, 0.965010)  # E2-15: the orientation of a uniform field that parts T915-S's senses most
# by the static balance (worst_orientation: 1000 orientations over the sphere and a refinement); a unit vector in the
# design's frame, z along the knot's axis

# ---- what agrees: fixed before the run (PROTOCOL.md, section 5) --------------------------------------------------
MOTION = (0.975, 1.035)  # a tracked limit over the exact motion's: stage 1's five holds 0.98 to 1.03, and half a step
PLANE = (0.975, 1.03)  # H1, E2-9 and H4: a planar path's limit over its static limit (stage 1: 0.994 to 1.004)
OUT_OF_PLANE = (0.845, 0.89)  # H1: a path out of the plane over its static limit (stage 1: 0.859 and 0.873)
ACROSS = (0.985, 1.035)  # H1b: a limit over the bound with every bend counted across the plane (1.001 and 1.020)
SENSE_STEPS = 2  # the senses' readings within this many scan steps of each other
BAND_STEPS = 2  # H2: the band's readings within this many scan steps, of each other and of the exact motion's
ONSET = 2  # H2: the band's onset, the lowest energy at which this many of a sense's eight electrons were lost
LOBE_REACH = 0.07  # H3: a lost electron within a lobe or less than this fraction of a lap past one, downstream
SPLIT_BAND = (0.5, 2.0)  # E2-15: half to twice the predicted split of the senses
TOWARD_RATIO = 2.0  # E2-9: the stiff ring holds at least twice what the flat ring holds
SHIFT_TOL = 0.02  # E2-2: the twin's shift within 2 % of the axisymmetric reference at R = 40h
WELL_STEPS = rung.WELL_STEPS  # H4: the well's first loss within two scan steps (and half a step) of the exact motion's
SPREADS = 2.0  # E2-13: a growth's sign judged beyond this many spreads over the seeds

Kind = Literal["hold", "well", "fail", "control", "shift", "multipactor"]
Table3 = tuple[tuple[tuple[int, ...], ...], ...]
Float3 = tuple[tuple[tuple[float, ...], ...], ...]


def sig(value: float, digits: int = 6) -> float:
    return rung.sig(value, digits)


def _path(design: str) -> Path:
    path = Path(design)
    return path if path.is_absolute() else lab.ROOT / path


# ---- theory: a design at a realisation and a drive -----------------------------------------------------------------
@cache
def guide_of(design: str, realisation: int, cache: str | None) -> TwinGuide:
    """A design's guide at a realisation, opened once a process."""
    return TwinGuide(design, realisation, cache=cache)


@cache
def theory(design: str, realisation: int = 0, q: float | None = None, cache: str | None = None) -> rung.Theory:
    """Stage 1's static balance (twin_rung1.Theory) on a design at a realisation, driven as the design says or to the
    null's Mathieu q at the mean volts: the realisation's own band, its drive's ripple and the gap errors its field
    carries, (1 + dh/h)^-2."""
    if realisation == 0 and q is None:
        return rung.theory(design, cache)
    guide = guide_of(design, realisation, cache)
    driven = guide.drive() if q is None else guide.drive(q=q)
    field = guide.field(driven)
    patch = rung.well_patch(guide.grid(), driven.f, guide.solved.distance)
    angles = np.asarray(np.arange(rung.ANGLES) * (360.0 / rung.ANGLES), dtype=np.float64)
    valleys = tuple(rung.valley(patch, math.radians(a)) for a in angles)
    strengths = np.array([v.strength(driven.volts) for v in valleys])
    count = len(field.swept.points)
    wave = rung._wave(field)
    amplitude: Arr = np.abs(wave.phasor(2 * math.pi * np.arange(count) / count)) / driven.volts
    ku, kv = rung.section_curvature(guide, field)
    kappa = np.hypot(ku, kv)
    heading = np.mod(np.degrees(np.arctan2(kv, ku)), 360.0)
    limits = rung._limits(angles, strengths, heading, kappa, amplitude)
    phasor = complex(wave.phasor(np.zeros(1))[0])
    return rung.Theory(
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


@cache
def frame_of(design: str, realisation: int, cache: str | None) -> tuple[Arr, Arr, Arr, Arr, Arr]:
    """The path as the section sees it at every sample: its tangent, first and second normals (n, 3), and the
    curvature along them, kappa_u in the boards' plane and kappa_v toward them (1/m)."""
    guide = guide_of(design, realisation, cache)
    field = guide.field()
    ku, kv = rung.section_curvature(guide, field)
    swept = field.swept
    return swept.tangents, swept.first, swept.second, ku, kv


def across_bound(th: rung.Theory, ku: Arr) -> float:
    """H1b's bound (eV): the static limit were every bend counted across the boards' plane: at every sample the
    strength across them on the bend's own side (along +u where kappa_u >= 0, else -u) against the whole curvature, the
    drive's ripple kept (stage 1's 9.80 eV for T915-S, 16.08 eV for the helix)."""
    held = th.strength(np.where(ku >= 0, 0.0, 180.0)) * th.amplitude**2
    return float(np.min(np.where(th.kappa > 0, held / (2 * np.maximum(th.kappa, 1e-300)), math.inf)))


def field_limit(
    th: rung.Theory, sense: int, field: Sequence[float] | Arr, frame: tuple[Arr, Arr, Arr, Arr, Arr]
) -> float:
    """The static limit (eV) of one sense of travel with a uniform field B (T): the most W for which, at every sample,
    the bend's 2W kappa and the field's force, s v (T x B) in the section (eV/m; s the sense), together stay within the
    strength along their direction (twin_rung1.ring_limit for any path)."""
    tangent, first, second, ku, kv = frame
    tb = np.cross(tangent, np.asarray(field, dtype=np.float64)[None, :])
    bu, bv = (tb * first).sum(1), (tb * second).sum(1)
    amp2 = th.amplitude**2
    lo, hi = 0.0, 4 * th.static
    for _ in range(60):
        w = (lo + hi) / 2
        v = math.sqrt(2 * w * QM)
        tilt_u, tilt_v = 2 * w * ku + sense * v * bu, 2 * w * kv + sense * v * bv
        cap = th.strength(np.degrees(np.arctan2(tilt_v, tilt_u))) * amp2
        lo, hi = (w, hi) if bool(np.all(np.hypot(tilt_u, tilt_v) <= cap)) else (lo, w)
    return lo


def worst_orientation(th: rung.Theory, size: float, frame: tuple[Arr, Arr, Arr, Arr, Arr], count: int = 1000) -> Arr:
    """The unit vector along which a uniform field of `size` (T) parts the two senses' static limits most,
    |ln(W- / W+)|: the best of `count` orientations spread evenly over the sphere, then refined by halving steps."""

    def split(way: Arr) -> float:
        b = size * way / np.linalg.norm(way)
        return abs(math.log(field_limit(th, -1, b, frame) / field_limit(th, 1, b, frame)))

    k = np.arange(count) + 0.5
    polar, around = np.arccos(1 - 2 * k / count), math.pi * (1 + math.sqrt(5)) * k
    ways = np.stack([np.sin(polar) * np.cos(around), np.sin(polar) * np.sin(around), np.cos(polar)], 1)
    scores = [split(w) for w in ways]
    best = ways[int(np.argmax(scores))]
    top, step = max(scores), 0.1
    for _ in range(12):
        moved = False
        for axis in range(3):
            for sign in (1.0, -1.0):
                trial = best.copy()
                trial[axis] += sign * step
                trial /= np.linalg.norm(trial)
                score = split(trial)
                if score > top:
                    best, top, moved = trial, score, True
        if not moved:
            step /= 2
    out: Arr = best / np.linalg.norm(best)
    return out


# ---- theory: the lobes of a path out of the plane (H3) -------------------------------------------------------------
def lobes(ku: Arr) -> tuple[tuple[float, float], ...]:
    """The lobes of a path whose bend turns through the boards' plane (stage 1's report, section 6): the stretches
    between consecutive pairs of maxima of the curvature in the boards' plane, |kappa_u|, above 0.8 of its largest,
    in order from the start (fractions of a lap), a plateau of equal values one maximum; none where that curvature is
    constant, as on a circle, or where the maxima do not come in pairs that repeat round the path's symmetry."""
    kg = np.abs(np.asarray(ku, dtype=np.float64))
    n = len(kg)
    if float(np.ptp(kg)) <= 1e-6 * float(kg.max()):
        return ()
    top = kg >= 0.8 * kg.max()
    peak = top & (kg >= np.roll(kg, 1)) & (kg >= np.roll(kg, -1))
    marks = np.flatnonzero(peak)
    groups: list[list[int]] = []
    for i in marks.tolist():
        if groups and i - groups[-1][-1] <= 2:
            groups[-1].append(i)
        else:
            groups.append([i])
    if len(groups) > 1 and groups[0][0] + n - groups[-1][-1] <= 2:
        groups[0] = [i - n for i in groups.pop()] + groups[0]
    places = [float(np.mean(g)) / n % 1.0 for g in groups]
    if len(places) < 2 or len(places) % 2:
        return ()
    pairs = tuple(sorted(zip(places[0::2], places[1::2], strict=True)))
    widths = [(b - a) % 1.0 for a, b in pairs]
    if max(widths) > 1.5 * min(widths) or max(widths) > 0.5:
        return ()
    return pairs


def past(place: float, sense: int, stretches: Sequence[tuple[float, float]]) -> float:
    """How far a loss at `place` (fraction of a lap) lies past the nearest lobe in its direction of travel (fraction
    of a lap; 0 within one; inf without lobes)."""
    out = math.inf
    for a, b in stretches:
        inside = a <= place <= b if a <= b else (place >= a or place <= b)
        if inside:
            return 0.0
        out = min(out, (place - b) % 1.0 if sense > 0 else (a - place) % 1.0)
    return out


# ---- theory: the exact motion in the section along the path --------------------------------------------------------
@dataclass(frozen=True)
class MotionSpec:
    """One exact motion: the design at a realisation and drive (None: the design's; else the null's q), the uniform
    field it runs in (T: the realisation's own unless given), the launch station (a sample of the path), and the
    shot's energies (eV), senses, RF phases (rad) and laps."""

    design: str
    realisation: int
    q: float | None
    earth: tuple[float, float, float]
    start: int
    energies: tuple[float, ...]
    senses: tuple[int, ...]
    phases: tuple[float, ...]
    laps: float
    steps: int = STEPS


@dataclass(frozen=True, eq=False)
class Table:
    """What became of a scan's electrons, each launched at an energy, a sense and an RF phase: its fate (0 held to the
    end of the run, 1 lost to the guide, 2 stopped by the gas, 3 still going but short of the run's laps), and how far
    it had gone when that was decided (laps for a hold, RF periods otherwise), [phase][sense][energy]; `run` is the
    laps (or periods) the scan was run for. A reading at fewer laps counts an electron held there once it had gone
    them."""

    energies: tuple[float, ...]
    senses: tuple[int, ...]
    phases: tuple[float, ...]
    fate: Table3
    reached: Float3
    run: float

    def _arrays(self) -> tuple[IntArr, Arr]:
        fate = np.array(self.fate, dtype=np.intp).reshape(len(self.phases), len(self.senses), len(self.energies))
        reached = np.array(self.reached, dtype=np.float64).reshape(fate.shape)
        return fate, reached

    def counts(self, at: float | None = None) -> tuple[IntArr, IntArr]:
        """Held and lost by sense and energy, (senses, energies), at `at` laps (the run's by default): held once it had
        gone them, lost when the guide lost it sooner; those the gas stopped sooner and those short of them count in
        neither."""
        when = self.run if at is None else at
        fate, reached = self._arrays()
        done = reached >= when - 1e-9
        held: IntArr = (done | (fate == 0) & (at is None)).sum(axis=0)
        lost: IntArr = ((fate == 1) & ~done).sum(axis=0)
        return held, lost

    def threshold(self, row: int, at: float | None = None) -> float:
        """A sense's robust threshold (eV; twin_rung1.robust_threshold): 0 below the scan, inf above it."""
        held, lost = self.counts(at)
        return rung.robust_threshold(self.energies, held[row].tolist(), lost[row].tolist())

    def limit(self, at: float | None = None) -> float:
        """The senses' robust thresholds' geometric mean (eV): 0 if either is below the scan, inf if either is above."""
        found = [self.threshold(row, at) for row in range(len(self.senses))]
        if min(found) == 0.0:
            return 0.0
        return math.inf if math.isinf(max(found)) else math.exp(sum(math.log(x) for x in found) / len(found))

    def onset(self, row: int, at: float | None = None) -> float:
        """The band's onset in a sense (eV): the lowest energy at which ONSET or more electrons were lost (inf if
        none)."""
        _, lost = self.counts(at)
        return next((e for e, n in zip(self.energies, lost[row].tolist(), strict=True) if n >= ONSET), math.inf)

    def edge(self, row: int, at: float | None = None) -> float:
        """The lowest energy (eV) from which every electron judged, there and at every energy above, was lost (inf if
        none; energies whose electrons the gas stopped every one say nothing)."""
        held, lost = self.counts(at)
        out = math.inf
        for k in range(len(self.energies) - 1, -1, -1):
            if held[row][k]:
                break
            if lost[row][k]:
                out = self.energies[k]
        return out

    def first(self, row: int, at: float | None = None) -> float:
        """The first loss in a sense (eV): the lowest energy at which an electron was lost (inf if none)."""
        _, lost = self.counts(at)
        return rung.lowest_lost(self.energies, lost[row].tolist())

    def launched(self) -> int:
        return len(self.phases) * len(self.senses) * len(self.energies)


@dataclass(frozen=True, eq=False)
class Motion:
    """The exact motion of a shot's electrons: what became of each (a Table), where along the path each lost one was
    lost (fractions of a lap, [phase][sense][energy], nan for the held), and the steps it ran."""

    table: Table
    place: Float3
    steps: int


def _motion_rows(field: TwinField, length: float, earth: Arr) -> tuple[Arr, float | None, float]:
    """The path as the exact motion reads it, sample by sample (n + 1, 9), the first sample repeated past the last:
    kappa_u and kappa_v (1/m), their derivatives along the path (1/m^2), the frame's turning about the tangent against
    the twist-free frame (rad/m), the drive's phasor (V, real and imaginary), and (T x B) along u and v for the uniform
    field B (T); with the wave's slowness (s/m) and omega (rad/s)."""
    swept = field.swept
    n = len(swept.points)
    ds = length / n
    ku = (swept.curvature_vector * swept.first).sum(1)
    kv = (swept.curvature_vector * swept.second).sum(1)
    turn = ((np.roll(swept.first, -1, 0) - np.roll(swept.first, 1, 0)) * swept.second).sum(1) / (2 * ds)
    dku = (np.roll(ku, -1) - np.roll(ku, 1)) / (2 * ds)
    dkv = (np.roll(kv, -1) - np.roll(kv, 1)) / (2 * ds)
    wave = rung._wave(field)
    phasor = wave.phasor(2 * math.pi * np.arange(n) / n)
    tb = np.cross(swept.tangents, earth[None, :])
    eu, ev = (tb * swept.first).sum(1), (tb * swept.second).sum(1)
    rows = np.column_stack([ku, kv, dku, dkv, turn, phasor.real, phasor.imag, eu, ev])
    return np.ascontiguousarray(np.vstack([rows, rows[:1]])), wave.slowness(length), wave.omega


def exact_motion(spec: MotionSpec, cache: str | None = None) -> Motion:
    """The exact motion in the section along the path (PROTOCOL.md section 2.1): every electron of the shot followed in
    two dimensions, (u, v) from the section's null in its frame, with its arc length s and p = s' h^2 (h = 1 - kappa.r,
    the curvilinear metric; p is the canonical momentum along the path over the mass), by the tracker's own leapfrog at
    `steps` an RF period, the Coriolis term turned exactly as Boris turns a magnetic one:

        r'' = |q/m| V(s, t) (1 - v_t / v_ph) grad phi(r) - (v_t^2 / h) kappa - |q/m| v_t (T x B) + 2 Omega J r'
              + Omega^2 r,      p' = -(v_t^2 / h) kappa'.r,      s' = p / h^2,      v_t = p / h,

    V the ring's phasor at s at the phase the electron was launched in, v_ph the wave's phase velocity (its magnetic
    force), kappa the curvature in the frame, B the uniform field, Omega = s' theta' the frame's turning about the
    tangent. Each starts as the tracker launches it (launch_offsets): on its static equilibrium at the station with its
    micromotion matched, its speed along the path that of its energy. It is lost where the section's grid loses it, off
    the grid or in copper or a board; held once it has gone the shot's laps."""
    guide = guide_of(spec.design, spec.realisation, cache)
    driven = guide.drive() if spec.q is None else guide.drive(q=spec.q)
    field = guide.field(driven)
    earth = np.asarray(spec.earth, dtype=np.float64)
    rows, slowness, omega = _motion_rows(field, guide.length, earth)
    slow = 0.0 if slowness is None else slowness
    n = len(rows) - 1
    length = guide.length
    ds = length / n
    grid = field.swept.grid
    u0, v0 = grid.null
    f = omega / (2 * math.pi)
    dt = 1 / (f * spec.steps)
    count = len(spec.phases) * len(spec.senses) * len(spec.energies)
    shape = (len(spec.phases), len(spec.senses), len(spec.energies))
    phase = np.repeat(np.asarray(spec.phases, dtype=np.float64), len(spec.senses) * len(spec.energies))
    sense = np.tile(np.repeat(np.asarray(spec.senses, dtype=np.float64), len(spec.energies)), len(spec.phases))
    energy = np.tile(np.asarray(spec.energies, dtype=np.float64), len(spec.phases) * len(spec.senses))
    start = spec.start % n
    ku0, kv0 = float(rows[start, 0]), float(rows[start, 1])
    p0 = complex(rows[start, 5], rows[start, 6])
    patch = rung.well_patch(guide.grid(), driven.f, guide.solved.distance)
    bend = rung.valley(patch, math.atan2(kv0, ku0))
    kappa = math.hypot(ku0, kv0)
    rest = np.array([bend.at(2 * e * kappa, abs(p0)) for e in spec.energies], dtype=np.float64).reshape(-1, 2)
    r = np.tile(rest, (len(spec.phases) * len(spec.senses), 1))
    t0 = np.mod((phase - math.atan2(p0.imag, p0.real)) / omega, 1 / f)
    gu, gv, _ = grid.gradient_at(u0 + r[:, 0], v0 + r[:, 1])
    g = np.column_stack([gu, gv])
    now = p0 * np.exp(1j * omega * t0)
    before = p0 * np.exp(1j * omega * (t0 - dt / 2))
    r = r - (QM / omega**2) * now.real[:, None] * g
    vr = (QM / omega) * before.imag[:, None] * g
    h = 1 - (ku0 * r[:, 0] + kv0 * r[:, 1])
    p = sense * np.array([pt.v_of(float(e)) for e in energy]) * h
    s = np.full(count, start * ds)
    clock = np.exp(1j * omega * t0)
    tick = complex(math.cos(omega * dt), math.sin(omega * dt))
    steps = math.ceil(spec.laps * 1.02 * length / (float(pt.v_of(min(spec.energies))) * dt))
    idx: IntArr = np.arange(count)
    lost_at = np.full(count, math.nan)
    lost_s = np.full(count, math.nan)
    for k in range(steps):
        x = s / ds
        i = np.floor(x).astype(np.intp)
        fr = (x - i)[:, None]
        i %= n
        row = rows[i] + fr * (rows[i + 1] - rows[i])
        gu, gv, on = grid.gradient_at(u0 + r[:, 0], v0 + r[:, 1])
        if not on.all():
            off = ~on
            lost_at[idx[off]] = np.abs(s[off] - start * ds) / length
            lost_s[idx[off]] = s[off]
            keep = np.flatnonzero(on)
            idx, r, vr, p, s, clock = idx[keep], r[keep], vr[keep], p[keep], s[keep], clock[keep]
            row, gu, gv = row[keep], gu[keep], gv[keep]
            if not len(idx):
                break
        ku, kv = row[:, 0], row[:, 1]
        h = 1 - (ku * r[:, 0] + kv * r[:, 1])
        vt = p / h
        along = vt / h
        volts = row[:, 5] * clock.real - row[:, 6] * clock.imag
        kick = QM * volts * (1 - slow * vt)
        centrifugal = vt * vt / h
        turning = along * row[:, 4]
        au = kick * gu - centrifugal * ku - QM * vt * row[:, 7] + turning * turning * r[:, 0]
        av = kick * gv - centrifugal * kv - QM * vt * row[:, 8] + turning * turning * r[:, 1]
        half_u, half_v = vr[:, 0] + (dt / 2) * au, vr[:, 1] + (dt / 2) * av
        c, sn = np.cos(2 * turning * dt), np.sin(2 * turning * dt)
        vr = np.column_stack([c * half_u + sn * half_v + (dt / 2) * au, c * half_v - sn * half_u + (dt / 2) * av])
        r = r + dt * vr
        p = p - dt * centrifugal * (row[:, 2] * r[:, 0] + row[:, 3] * r[:, 1])
        s = s + dt * along
        clock = clock * tick
        if k & 1023 == 1023:
            clock = clock / np.abs(clock)
    gone = np.full(count, math.nan)
    gone[idx] = np.abs(s - start * ds) / length
    reached = np.where(np.isnan(lost_at), gone, lost_at)
    fate = np.where(reached >= spec.laps - 1e-9, 0, np.where(np.isnan(lost_at), 3, 1))
    place = np.mod(lost_s / length, 1.0)
    table = Table(
        spec.energies,
        spec.senses,
        spec.phases,
        _nest_int(fate.reshape(shape)),
        _nest_float(reached.reshape(shape)),
        spec.laps,
    )
    return Motion(table, _nest_float(place.reshape(shape)), steps)


def _nest_int(a: IntArr) -> Table3:
    return tuple(tuple(tuple(int(x) for x in row) for row in block) for block in a.tolist())


def _nest_float(a: Arr) -> Float3:
    return tuple(tuple(tuple(float(x) for x in row) for row in block) for block in np.asarray(a).tolist())


# ---- the exact motion, kept by the hash of its inputs and code ----------------------------------------------------
_MOTION_CODE = (
    "array16/twin/__init__.py",
    "array16/particles/sections.py",
    "array16/particles/curves.py",
    "array16/guide/ring.py",
    "array16/guide/geometry.py",
    "experiments/twin_rung1.py",
)


@cache
def motion_code() -> str:
    """The hash of the code an exact motion runs: its own functions' source, the solver's version, and the twin's,
    the swept section's, the ring's and stage 1's sources."""
    digest = hashlib.sha256(solver_version().encode())
    for part in (exact_motion, _motion_rows, guide_of, MotionSpec, _nest_int, _nest_float):
        digest.update(inspect.getsource(part).encode())
    for name in _MOTION_CODE:
        digest.update((lab.ROOT / name).read_bytes())
    return digest.hexdigest()


def _motion_file(spec: MotionSpec, cache: str | None) -> Path:
    root = (CACHE if cache is None else Path(cache)) / "rung2"
    key = json.dumps({"format": 1, "code": motion_code(), "spec": lab.plain(spec)}, sort_keys=True)
    return root / f"{hashlib.sha256(key.encode()).hexdigest()}.json"


def _motion_plain(m: Motion) -> dict[str, Any]:
    out: dict[str, Any] = lab.plain(m)
    return out


def _motion_load(data: Mapping[str, Any]) -> Motion:
    return Motion(_table_load(data["table"]), _float3(data["place"]), int(data["steps"]))


def _float3(rows: Any) -> Float3:
    return tuple(tuple(tuple(math.nan if x is None else float(x) for x in r) for r in b) for b in rows)


def _int3(rows: Any) -> Table3:
    return tuple(tuple(tuple(int(x) for x in r) for r in b) for b in rows)


def _table_load(data: Mapping[str, Any]) -> Table:
    return Table(
        tuple(float(e) for e in data["energies"]),
        tuple(int(s) for s in data["senses"]),
        tuple(float(p) for p in data["phases"]),
        _int3(data["fate"]),
        _float3(data["reached"]),
        float(data["run"]),
    )


def _motion_job(args: tuple[MotionSpec, str | None]) -> dict[str, Any]:
    spec, cache = args
    return _motion_plain(exact_motion(spec, cache))


_MOTIONS: dict[tuple[MotionSpec, str | None], Motion] = {}


def motion(spec: MotionSpec, cache: str | None = None) -> Motion:
    """exact_motion, read back from .twin_cache/rung2/ (or the cache folder's rung2/) when this code has computed it
    for these inputs before, else computed and written there."""
    if (spec, cache) in _MOTIONS:
        return _MOTIONS[(spec, cache)]
    path = _motion_file(spec, cache)
    if path.exists():
        try:
            found = _motion_load(json.loads(path.read_text()))
        except (OSError, ValueError, KeyError, TypeError):
            found = None
        if found is not None:
            _MOTIONS[(spec, cache)] = found
            return found
    data = _motion_job((spec, cache))
    _keep(path, data)
    found = _motion_load(data)
    _MOTIONS[(spec, cache)] = found
    return found


def _keep(path: Path, data: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    scratch = path.with_name(f".{path.stem}.{os.getpid()}.tmp")
    scratch.write_text(json.dumps(data))
    os.replace(scratch, path)


def motions(specs: Sequence[MotionSpec], cache: str | None, workers: int) -> None:
    """Compute in a process pool the exact motions not yet kept, longest first, and keep them."""
    todo = [s for s in dict.fromkeys(specs) if not _motion_file(s, cache).exists() and (s, cache) not in _MOTIONS]
    if not todo:
        return
    todo.sort(key=lambda s: -_motion_cost(s))
    if workers <= 1 or len(todo) == 1:
        for spec in todo:
            motion(spec, cache)
        return
    with ProcessPoolExecutor(max_workers=min(workers, len(todo))) as pool:
        for spec, data in zip(todo, pool.map(_motion_job, [(s, cache) for s in todo]), strict=True):
            _keep(_motion_file(spec, cache), data)


def _motion_cost(spec: MotionSpec) -> float:
    """A rough cost, to order the pool's work: the steps times the electrons and the step's own cost."""
    count = len(spec.phases) * len(spec.senses) * len(spec.energies)
    return spec.laps * spec.steps / math.sqrt(min(spec.energies)) * (100 + count)


# ---- the shots ----------------------------------------------------------------------------------------------------
Earth = str | tuple[float, float, float]  # "drawn" (the realisation's own), "none", or a uniform field (T)


@dataclass(frozen=True)
class Shot:
    """One run: its test, the design (a file relative to the repository, or for multipactor a case's name) by a short
    label, and its realisation (None: the bench's guide, benches/twin_r2.toml's); what is measured; the energies (eV)
    scanned from `low` up by `ratio` while within `high`, along the path for a hold and across it for a well; the RF
    phases launched; how long (laps for a hold, RF periods otherwise); the drive (None: the design's; else the null's q
    at the mean volts); Earth's field ("drawn": the realisation's own, none as drawn; "none"; or a field, T); whether
    the bend's shift of the null is in the field; where along the path the electrons are launched (a fraction of a
    lap); the seed its launch is drawn from (None: its own); the tracker's steps an RF period, and its seed."""

    test: str
    label: str
    design: str
    realisation: int | None
    kind: Kind
    low: float = 0.0
    high: float = 0.0
    ratio: float = FINE
    phases: int = PHASES
    laps: float = LAPS
    periods: int = 0
    q: float | None = None
    earth: Earth = "drawn"
    bend: bool = False
    start: float = 0.0
    launch_seed: int | None = None
    steps: int = STEPS
    seed: int = 0

    @property
    def ident(self) -> str:
        where = "on the bench" if self.realisation is None else f"realisation {self.realisation}"
        if self.kind in ("shift", "multipactor"):
            return f"{self.test} {self.kind}, {self.label}"
        what = "" if self.kind == "fail" else f" {self.kind}"
        drive = "" if self.q is None else f" at q {self.q:g}"
        field = ""
        if self.earth == "none":
            field = ", Earth's field removed"
        elif isinstance(self.earth, tuple):
            field = f", {math.hypot(*self.earth) * 1e6:g} uT at the worst orientation"
        bend = ", the bend's shift on" if self.bend else ""
        laps = f", {self.laps:g} laps" if self.kind == "hold" and self.laps != LAPS else ""
        steps = "" if self.steps == STEPS else f", {self.steps} steps an RF period"
        return f"{self.test}{what}, {self.label}, {where}{drive}{field}{bend}{laps}{steps}"

    @property
    def energies(self) -> tuple[float, ...]:
        """The energies scanned (eV), from low up by `ratio` while within high."""
        return rung.Shot("", "", "", 0, "hold", self.low, self.high, self.ratio).energies

    @property
    def launch(self) -> int:
        """The seed its starts, RF phases and gas are drawn from: shared by shots registered as a pair."""
        return self.seed if self.launch_seed is None else self.launch_seed

    @property
    def senses(self) -> tuple[int, ...]:
        return (0,) if self.kind == "well" else (1, -1)

    @property
    def run(self) -> float:
        """What a scan is run for: laps for a hold, RF periods otherwise."""
        return self.laps if self.kind == "hold" else float(self.periods)


def _scan(test: str, label: str, design: str, realisation: int | None, kind: Kind, **more: Any) -> Shot:
    return Shot(test, label, design, realisation, kind, **more)


@cache
def shots() -> tuple[Shot, ...]:
    """Every shot, in the protocol's order (section 4); the energies from each path's static limit (section 2.2)."""
    worst = tuple(EARTH * c for c in WORST)
    strong = tuple(STRONG * EARTH * c for c in WORST)
    drawn = [
        # H1 and E2-9 at q 0.49: section B at 50 V, 20 laps
        _scan("H1", "ellipse 318", ELLIPSE, 0, "hold", low=9.5, high=12.7, start=ELLIPSE_START),
        _scan("H1", "flat ring 318", FLAT, 0, "hold", low=29.0, high=34.5),
        _scan("H1", "flat ring 318, stiff", STIFF, 0, "hold", low=26.0, high=165.0, ratio=COARSE),
        _scan("H1", "helix (1, 3), r/R 0.1", HELIX_010, 0, "hold", low=17.0, high=36.5),
        _scan("H1", "helix (1, 3), r/R 0.2", HELIX, 0, "hold", low=13.5, high=19.5),
        _scan("H1", "helix (1, 3), r/R 0.3", HELIX_030, 0, "hold", low=10.0, high=14.3),
        # E2-9 at q 0.30
        _scan("E2-9", "flat ring 318", FLAT, 0, "hold", low=10.3, high=13.1, q=0.3),
        _scan("E2-9", "flat ring 318, stiff", STIFF, 0, "hold", low=19.0, high=69.0, ratio=COARSE, q=0.3),
        # H2: T915-S for 100 laps
        _scan("H2", "T915-S", T915S, 0, "hold", low=6.4, high=11.9, laps=LONG),
        _scan("H2", "T915-S", T915S, 0, "hold", low=6.4, high=11.9, laps=LONG, steps=2 * STEPS),
        _scan("H2", "T915-S", T915S, 0, "hold", low=2.2, high=4.5, laps=LONG, q=0.3),
        # E2-15: Earth's field across T915-S
        _scan("E2-15", "T915-S", T915S, 0, "hold", low=7.55, high=12.1, earth=worst),
        _scan("E2-15", "T915-S", T915S, 0, "hold", low=5.2, high=12.1, earth=strong),
        _scan("E2-15", "T915-S", T915S, None, "hold", low=7.0, high=11.3),
        _scan("E2-15", "T915-S", T915S, None, "hold", low=7.0, high=11.3, earth="none"),
        # E2-2: the bend's shift of the null
        _scan("E2-2", "section B on a ring of R = 40h, bent in the boards' plane", RING_IN, 0, "shift"),
        _scan("E2-2", "section B on a ring of R = 40h, bent toward the boards", RING_TOWARD, 0, "shift"),
        _scan("E2-2", "T915-S", T915S, 0, "hold", low=8.15, high=12.1, bend=True, launch_seed=MASTER + 100),
        _scan("E2-2", "T915-S", T915S, 0, "hold", low=8.15, high=12.1, launch_seed=MASTER + 100),
        # H4: R1 on the bench's realisation
        _scan("H4", "R1", R1, 1, "hold", low=6.55, high=10.85, earth="none"),
        _scan("H4", "R1", R1, 1, "hold", low=6.1, high=10.85),
        _scan("H4", "R1", R1, 1, "well", low=0.0045, high=0.14, ratio=COARSE, periods=WELL_PERIODS, earth="none"),
        _scan("H4", "R1", R1, 1, "well", low=0.0045, high=0.14, ratio=COARSE, periods=WELL_PERIODS),
        # E2-13: multipactor with secondary emission in the boards' own field
        _scan("E2-13", "T915-S's section B at 50.3 V", "t915s", 0, "multipactor"),
        _scan("E2-13", "section B, split h, refitted for 1 eV at q 0.5, 39.6 V", "lever", 0, "multipactor"),
        _scan("E2-13", "R1's rails at 50 V", "r1", 0, "multipactor"),
        _scan("E2-13", "section B at 160 V, the control", "overdriven", 0, "multipactor"),
        _scan("E2-13", "facing strips at q 0.5 fitted to 1.5 to 3 eV, the onset", "onset", 0, "multipactor"),
        # the case that must fail, beside its control
        _scan("must fail", "T915-S", T915S, 0, "fail", low=0.5, high=4.0, ratio=2.0, phases=FAIL_PHASES,
              periods=FAIL_PERIODS, q=FAIL_Q, launch_seed=MASTER),
        _scan("must fail", "T915-S", T915S, 0, "control", low=0.5, high=4.0, ratio=2.0, phases=FAIL_PHASES,
              periods=FAIL_PERIODS, launch_seed=MASTER),
    ]  # fmt: skip
    return tuple(replace(shot, seed=MASTER + k) for k, shot in enumerate(drawn))


def _phases(shot: Shot) -> tuple[float, ...]:
    """The RF phases (rad) a shot launches at: a well's two zero crossings; otherwise one drawn uniformly in each of
    `phases` equal parts of the period, from the shot's launch seed (twin_rung1._phases)."""
    if shot.kind == "well":
        return rung.WELL_PHASES
    rng = np.random.default_rng(np.random.SeedSequence([shot.launch, 1]))
    return tuple(float(2 * math.pi * (k + rng.uniform()) / shot.phases) for k in range(shot.phases))


def _start_index(shot: Shot, samples: int) -> int:
    return round(shot.start * samples) % samples


def _realisation(shot: Shot) -> int:
    """The realisation a shot's theory takes: the bench's (twin_r2: T915-S realisation 1) for None."""
    return 1 if shot.realisation is None else shot.realisation


def _earth_vector(shot: Shot, guide: TwinGuide) -> tuple[float, float, float]:
    """The uniform field a shot runs in (T)."""
    if isinstance(shot.earth, tuple):
        return (float(shot.earth[0]), float(shot.earth[1]), float(shot.earth[2]))
    if shot.earth == "none":
        return (0.0, 0.0, 0.0)
    b = guide.drawn.earth
    return (float(b[0]), float(b[1]), float(b[2]))


def motion_spec(shot: Shot, cache: str | None) -> MotionSpec:
    """The exact motion that predicts a hold: the shot's own design, realisation, drive, field, station, energies,
    senses, phases and laps."""
    realisation = _realisation(shot)
    guide = guide_of(shot.design, realisation, cache)
    return MotionSpec(
        shot.design,
        realisation,
        shot.q,
        _earth_vector(shot, guide),
        _start_index(shot, guide.design.path.samples),
        shot.energies,
        shot.senses,
        _phases(shot),
        shot.laps,
        shot.steps,
    )


# ---- the theory's side --------------------------------------------------------------------------------------------
@dataclass(frozen=True, eq=False)
class Prediction:
    """What theory says a shot measures: the claim in words, the number and the band the measurement must fall in, its
    unit, and the numbers it rests on (the hypotheses' bands among them, as `<name>_low` and `<name>_high`)."""

    claim: str
    value: float
    low: float
    high: float
    unit: str
    basis: Mapping[str, float]


def _band(value: float, band: tuple[float, float]) -> tuple[float, float]:
    return sig(value * band[0]), sig(value * band[1])


def _planar(design: str) -> bool:
    return guide_of(design, 0, None).design.path.kind in ("circle", "ellipse")


def _motion_basis(m: Motion, static: float, prefix: str = "motion") -> dict[str, float]:
    """The exact motion's readings: its limit, each sense's robust threshold and, for a hold, its band's onset and
    its edge (eV), and each over the static limit."""
    t = m.table
    out = {f"{prefix}_limit_eV": t.limit(), f"{prefix}_fraction": t.limit() / static}
    for row, sense in enumerate(t.senses):
        name = {1: "with", -1: "against", 0: "row"}[sense]
        out[f"{prefix}_{name}_eV"] = t.threshold(row)
        out[f"{prefix}_{name}_onset_eV"] = t.onset(row)
        out[f"{prefix}_{name}_edge_eV"] = t.edge(row)
    return {key: sig(value) for key, value in out.items()}


def _common(th: rung.Theory) -> dict[str, float]:
    solved, driven = th.guide.solved, th.driven
    out = {
        "f_MHz": driven.f / 1e6,
        "volts_mean": th.volts,
        "q_null": th.q,
        "depth_weakest_eV": driven.depth_ev * th.scale**2,
        "null_mm": solved.distance * 1e3,
        "length_mm": th.guide.length * 1e3,
        "static_eV": th.static,
        "static_at": th.where / len(th.limits),
    }
    curve = th.guide.curve
    if curve is not None:
        out["bend_mm"] = curve.bend_radius * 1e3
    if th.guide.framing is not None:
        out["inplane_radius_mm"] = th.guide.framing.inplane_radius * 1e3
    return {key: sig(value, 6 if key != "static_at" else 4) for key, value in out.items()}


def hold_theory(shot: Shot, cache: str | None) -> tuple[rung.Theory, Motion]:
    realisation = _realisation(shot)
    th = theory(shot.design, realisation, shot.q, cache)
    return th, motion(motion_spec(shot, cache), cache)


def predict_hold(shot: Shot, cache: str | None) -> Prediction:
    """A hold: the exact motion's limit and the band MOTION round it (Claim), with the static balance, H1b's bound,
    the field's static balance where there is a field, and the hypotheses' bands in the basis."""
    th, m = hold_theory(shot, cache)
    static = th.static
    t = m.table
    value = t.limit()
    basis = _common(th) | _motion_basis(m, static)
    frame = frame_of(shot.design, _realisation(shot), cache)
    across = across_bound(th, frame[3])
    basis["across_eV"] = sig(across)
    plane = _planar(shot.design)
    hypothesis = PLANE if plane else OUT_OF_PLANE
    basis["H1_low_eV"], basis["H1_high_eV"] = _band(static, hypothesis)
    basis["H1b_low_eV"], basis["H1b_high_eV"] = _band(across, ACROSS)
    basis["planar"] = 1.0 if plane else 0.0
    b = np.asarray(motion_spec(shot, cache).earth, dtype=np.float64)
    if bool(np.any(b)):
        plus, minus = field_limit(th, 1, b, frame), field_limit(th, -1, b, frame)
        basis |= {
            "earth_uT": sig(float(np.linalg.norm(b)) * 1e6),
            "field_with_eV": sig(plus),
            "field_against_eV": sig(minus),
            "field_split": sig(math.log(minus / plus)),
        }
    if not plane and th.guide.curve is not None:
        stretches = lobes(frame[3])
        basis["lobes"] = float(len(stretches))
        if stretches:
            basis["motion_in_lobes"] = sig(_in_lobes(m, stretches))
    curve = th.guide.curve
    if plane and curve is not None:  # E2-9's measure: D R/h, D the weakest well's depth and R the tightest bend
        basis["DRh_eV"] = sig(th.driven.depth_ev * th.scale**2 * curve.bend_radius / th.guide.solved.distance)
    if shot.laps > LAPS:  # its reading at 20 laps, for H1
        twenty = m.table.limit(LAPS)
        basis["motion20_limit_eV"], basis["motion20_fraction"] = sig(twenty), sig(twenty / static)
        for row, sense in enumerate(m.table.senses):
            basis[f"motion20_{'with' if sense > 0 else 'against'}_eV"] = sig(m.table.threshold(row, LAPS))
    name = th.guide.design.name
    laps = shot.laps
    claim = (
        f"{name} holds electrons for {laps:g} laps to {value:.4g} eV ({value / static:.3f} of its static limit, "
        f"{static:.4g} eV), the exact motion in the section on the shot's own launch, within {MOTION[0]} to "
        f"{MOTION[1]} of it"
    )
    if shot.test == "H1":
        claim += (
            f"; H1 has it hold {hypothesis[0]} to {hypothesis[1]} of its static limit, H1b {ACROSS[0]} to {ACROSS[1]} "
            "of "
            f"its bound with every bend counted across the boards' plane, {across:.4g} eV"
        )
    low, high = _band(value, MOTION)
    return Prediction(claim, sig(value), low, high, "eV", basis)


def _in_lobes(m: Motion, stretches: Sequence[tuple[float, float]]) -> float:
    """The share of the exact motion's lost electrons lost within a lobe or LOBE_REACH past one, downstream."""
    t = m.table
    fate = np.array(t.fate)
    place = np.array(m.place)
    count = inside = 0
    for k in range(len(t.phases)):
        for row, sense in enumerate(t.senses):
            for j in range(len(t.energies)):
                if fate[k, row, j] == 1:
                    count += 1
                    inside += past(float(place[k, row, j]), sense, stretches) <= LOBE_REACH
    return inside / count if count else math.nan


def well_theory(shot: Shot, cache: str | None) -> tuple[rung.Theory, float, BoolArr]:
    """The exact motion on the axis (twin_rung1.section_flights) at the realisation's own band and the start station's
    volts: which of the shot's electrons leave the well."""
    th = theory(shot.design, _realisation(shot), shot.q, cache)
    solved = th.guide.solved
    du, dv = solved.saddle[0] - solved.null[0], solved.saddle[1] - solved.null[1]
    size = math.hypot(du, dv)
    omegas = np.array([2 * math.pi * th.driven.f])
    lost = rung.section_flights(
        th.guide.grid(), abs(th.phasor), omegas, shot.energies, rung.WELL_PHASES, shot.periods, (du / size, dv / size)
    )
    return th, abs(th.phasor), lost[0]


def predict_well(shot: Shot, cache: str | None) -> Prediction:
    th, volts, lost = well_theory(shot, cache)
    first = rung.lowest_lost(shot.energies, lost.any(axis=0).tolist())
    depth = pt.ponderomotive_ev(th.guide.solved.escape * volts, th.driven.f)
    room = shot.ratio ** (WELL_STEPS + 0.5)
    basis = _common(th) | {"depth_start_eV": sig(depth), "volts_at_start": sig(volts), "model_first_D": first / depth}
    claim = (
        f"electrons sent from {th.guide.design.name}'s null up its axis toward the saddle at the drive's zero "
        "crossings, "
        f"on the bench's realisation at its own band, {th.driven.f / 1e6:.3f} MHz, are first lost at {first:.4g} eV, "
        f"{first / depth:.2f} of the depth there, as the exact motion on the solved axis loses them in {shot.periods} "
        "RF periods"
    )
    return Prediction(claim, sig(first), sig(first / room), sig(first * room), "eV", basis)


def predict_fail(shot: Shot, cache: str | None) -> Prediction:
    th = theory(shot.design, 0, None, cache)
    if shot.kind == "control":
        top = max(shot.energies)
        claim = (
            f"the must-fail launch at the design's {th.volts:g} V (q {th.q:.2f}): every electron held for "
            f"{shot.periods} RF periods, {top:g} eV being {top / th.static:.2f} of the static limit"
        )
        basis = {"volts": sig(th.volts), "q_null": sig(th.q), "static_eV": sig(th.static)}
        return Prediction(claim, 0.0, 0.0, 0.0, "lost", basis)
    assert shot.q is not None
    volts = shot.q / th.guide.solved.mathieu_q(1.0)
    claim = (
        f"driven to {volts:.1f} V the null's q is {shot.q:g}, past Mathieu's {rung.QCAP}: no well forms, and no "
        "electron launched near the null is held"
    )
    basis = {"volts": sig(volts), "q_null": sig(shot.q), "depth_eV": sig(th.guide.solved.depth_ev(volts))}
    return Prediction(claim, 0.0, 0.0, 0.0, "held", basis)


def ring_radius(shot: Shot) -> float:
    """R/h of an E2-2 ring."""
    guide = guide_of(shot.design, 0, None)
    return guide.length / (2 * math.pi * guide.design.section.geometry.distance)


def _way(shot: Shot) -> axisym.Way:
    """Which way an E2-2 ring bends: in the boards' plane where its frame puts u across the ring's plane's radius
    (theta 0 or 180 degrees), toward them otherwise (theta 90)."""
    theta = guide_of(shot.design, 0, None).design.frame.theta
    return "in-plane" if abs(math.cos(theta)) > 0.5 else "toward"


def predict_shift(shot: Shot) -> Prediction:
    """E2-2: the axisymmetric reference's shift x R/h^2 at the ring's R/h, within SHIFT_TOL."""
    way, radius = _way(shot), ring_radius(shot)
    found = axisym.bend(radius, way)
    claim = (
        f"the twin moves section B's null toward the centre of a ring of R = {radius:.0f}h by {found.coefficient:.6f} "
        f"h^2/R {'in the boards plane' if way == 'in-plane' else 'toward the boards'}, the axisymmetric reference's "
        f"(first order {found.first_order:.6f}; E2-2's {0.1405 if way == 'in-plane' else 0.318})"
    )
    basis = {"radius_h": sig(radius), "first_order": sig(found.first_order), "gradient": sig(found.gradient)}
    low, high = _band(found.coefficient, (1 - SHIFT_TOL, 1 + SHIFT_TOL))
    return Prediction(claim, sig(found.coefficient), low, high, "h^2/R", basis)


# ---- E2-13: the cases, and the level-1 chart's reading of each -------------------------------------------------------
Surface = Literal["bare", "conditioned", "chart"]
SURFACES: Mapping[str, tuple[Surface, ...]] = {
    "t915s": ("bare", "conditioned", "chart"),
    "lever": ("bare", "conditioned", "chart"),
    "r1": ("bare", "conditioned", "chart"),
    "overdriven": ("bare", "chart"),
    "onset": ("bare", "chart"),
}
LAMINATES = ("absorbs", "emits")


def surface_of(name: Surface) -> tuple[em.Surface, float]:
    """A surface and the volts gain it is read through: Furman and Pivi's copper as published and conditioned (gain
    one), and the level-1 chart's own yield (multipactor's GAIN)."""
    if name == "bare":
        return em.FP_COPPER, 1.0
    if name == "conditioned":
        return em.FP_CONDITIONED, 1.0
    return em.CHART, mp.GAIN


@cache
def mpc_cases(name: str) -> tuple[em.Case, ...]:
    """The cases of an E2-13 shot, on the default grid."""
    if name == "onset":
        return tuple(em.facing_case(depth) for depth in ONSET_DEPTHS)
    if name == "overdriven":
        return (em.Case(f"section B at {OVERDRIVE:g} V", section_b(), OVERDRIVE),)
    three = em.cases()
    return ({"t915s": three[0], "lever": three[1], "r1": three[2]}[name],)


def chart_gap(case: em.Case) -> mp.Gap:
    """The gap the level-1 chart reads a case as: facing strips 2h apart at their volts; R1's rails as multipactor's
    bound on their arcs, from the slit to pi times the rail's width (MPC-3)."""
    geometry = case.section.geometry
    if case.name.startswith("R1"):
        return next(g for g in mp.CASES if g.name.startswith("R1"))
    gap = 2 * geometry.distance
    return mp.Gap(case.name, case.section.f, case.volts, gap, gap)


def predict_mpc(shot: Shot) -> Prediction:
    """E2-13 from the level-1 chart: each case's margin (dB) and the one-dimensional operator's growth (Np an RF period)
    and effective yield on the chart's yield at its worst field line; decays where the margin is positive, and R1's
    rails bounded by their parallel-plate growth."""
    basis: dict[str, float] = {}
    lines = []
    for k, case in enumerate(mpc_cases(shot.design)):
        reading = mp.level1(chart_gap(case))
        growth = mp.growth_rate(reading.fd, case.volts)
        effective = float(mp.effective_yield(reading.fd, case.volts)[0])
        tag = "" if len(mpc_cases(shot.design)) == 1 else f"{ONSET_DEPTHS[k]:g}eV_"
        basis |= {
            f"{tag}volts": sig(case.volts),
            f"{tag}fd_GHz_mm": sig(reading.fd),
            f"{tag}margin_dB": sig(reading.margin),
            f"{tag}chart_growth": sig(growth),
            f"{tag}chart_effective_yield": sig(effective),
        }
        lines.append(f"{case.name} at {case.volts:.4g} V: {reading.margin:+.2f} dB, {growth:+.3f} Np a period")
    said = "; ".join(lines)
    if shot.design == "onset":
        margins = [basis[f"{d:g}eV_margin_dB"] for d in ONSET_DEPTHS]
        crossing = next((d for d, m in zip(ONSET_DEPTHS, margins, strict=True) if m <= 0), math.inf)
        where = f"at {crossing:g} eV" if math.isfinite(crossing) else f"above {ONSET_DEPTHS[-1]:g} eV"
        claim = (
            f"the level-1 chart puts the onset on facing strips at q 0.5 {where} ({said}); Hypothesis (MPC-H1, "
            "E2-13): in the strips' own field the population decays at 1.5 eV and grows at 3 eV on copper's yield, in "
            "both laminate models"
        )
        return Prediction(claim, ONSET_DEPTHS[0], ONSET_DEPTHS[0], ONSET_DEPTHS[-1], "eV", basis)
    if shot.design == "overdriven":
        claim = f"the control: {said}; in its own field the population grows on every surface and laminate model"
        return Prediction(claim, basis["chart_growth"], 0.0, math.inf, "Np", basis)
    if shot.design == "r1":
        claim = (
            f"{said} (the parallel-plate bound on the rails' arcs): in the rails' own field the population grows no "
            "faster; Hypothesis (read from the emission module's own readings): it dies, its electrons leaving the "
            "rails' field before they multiply"
        )
        return Prediction(claim, basis["chart_growth"], -math.inf, basis["chart_growth"], "Np", basis)
    claim = f"{said}: positive, so the population decays on every surface and laminate model"
    return Prediction(claim, basis["chart_growth"], -math.inf, 0.0, "Np", basis)


def predict(shot: Shot, cache: str | None = None) -> Prediction:
    if shot.kind == "hold":
        return predict_hold(shot, cache)
    if shot.kind == "well":
        return predict_well(shot, cache)
    if shot.kind in ("fail", "control"):
        return predict_fail(shot, cache)
    if shot.kind == "shift":
        return predict_shift(shot)
    return predict_mpc(shot)


# ---- the tracking's side ------------------------------------------------------------------------------------------
@dataclass(frozen=True, eq=False)
class Scan:
    """What a scan found (a hold, a well, a failure or its control): what became of each electron (a Table), and for
    each RF phase the detector's entries of the electrons lost or stopped (array16.twin.records: how, when, where along
    the path and in the section, with what energy and sense); the run (steps of dt s), the drive (Hz, V mean, the null's
    q), the uniform field it ran in (T), whether the bend's shift was on, the launch station, the guide's identity and
    the seed."""

    kind: str
    guide: Mapping[str, Any]
    seed: int
    table: Table
    losses: tuple[Losses, ...]
    steps: int
    dt: float
    f: float
    volts: float
    q: float
    earth: tuple[float, float, float]
    bend: bool
    start: int


@dataclass(frozen=True, eq=False)
class Shift:
    """E2-2's reading of the twin's bend correction: where the swept field vanishes in the normal planes of a ring of
    R/h (`radius`) at its stations, the null's offsets from the path along the frame's first and second normals (m);
    the mean offset along the bend's axis times R/h^2, and the largest offset across it over that."""

    kind: str
    guide: Mapping[str, Any]
    seed: int
    way: str
    radius: float
    stations: tuple[int, ...]
    offsets: tuple[tuple[float, float], ...]
    coefficient: float
    cross: float


@dataclass(frozen=True)
class MpcReading:
    """One E2-13 reading (array16.guide.emission.read): a case on a surface, the laminate absorbing or emitting as the
    copper, its volts; the tracked growth per RF period over the seeds and its spread (-inf: the population died), the
    rate it fell at while it lived, the effective yield and its spread, a strike's mean energy (eV) and the share of
    the population leaving the region a period."""

    case: str
    surface: str
    laminate: str
    volts: float
    growth: float
    spread: float
    decay: float
    secondaries: float
    secondaries_spread: float
    energy: float
    escaped: float
    died: bool

    @classmethod
    def of(cls, data: Mapping[str, Any]) -> MpcReading:
        """A reading back from plain data, where a number that was not finite is null: a growth of -inf where the
        population died, nan otherwise."""

        def number(key: str, missing: float = math.nan) -> float:
            value = data[key]
            return missing if value is None else float(value)

        died = bool(data["died"])
        return cls(
            str(data["case"]),
            str(data["surface"]),
            str(data["laminate"]),
            number("volts"),
            number("growth", -math.inf if died else math.nan),
            number("spread"),
            number("decay"),
            number("secondaries"),
            number("secondaries_spread"),
            number("energy"),
            number("escaped"),
            died,
        )


@dataclass(frozen=True, eq=False)
class Mpc:
    """E2-13's readings of a shot's cases, on its surfaces and both laminate models, on the seeds."""

    kind: str
    guide: Mapping[str, Any]
    seed: int
    seeds: tuple[int, ...]
    readings: tuple[MpcReading, ...]


Measurement = Scan | Shift | Mpc


@dataclass(frozen=True)
class Job:
    """One batch of integrate: a scan's electrons at one RF phase, with what a worker needs to open the devices again:
    the guide's design, realisation and cache folder, the detector's weight and collector, and the vacuum's pressure
    (mbar), temperature (K) and cross-section (m^2)."""

    shot: Shot
    index: int
    phase: float
    guide: tuple[str, int, str]
    weight: float
    collector: float
    vacuum: tuple[float, float, float]


@dataclass(frozen=True)
class MpcJob:
    """One E2-13 reading: the shot's case by its place, the surface and the laminate model."""

    shot: Shot
    case: int
    surface: Surface
    laminate: str


def station_at(guide: TwinGuide, field: TwinField, index: int) -> rung.Station:
    """Where a shot starts: the path's sample `index`, its tangent and normals, the curvature along them and the
    drive's phasor there (twin_rung1.station at any sample)."""
    swept = field.swept
    n = len(swept.points)
    ku, kv = rung.section_curvature(guide, field)
    phasor = complex(rung._wave(field).phasor(np.array([2 * math.pi * index / n]))[0])
    return rung.Station(
        swept.points[index],
        swept.tangents[index],
        swept.first[index],
        swept.second[index],
        (ku[index], kv[index]),
        phasor,
    )


def launch(
    field: TwinField,
    start: rung.Station,
    offsets: Arr,
    shift: tuple[float, float],
    velocity: Arr,
    t0: float,
    omega: float,
    dt: float,
) -> tuple[Arr, Arr]:
    """twin_rung1.launch with the null moved by `shift` (du, dv, m) at the station, where the bend's shift is on:
    electrons at `offsets` from the null as the field places it, in step with the drive's micromotion there."""
    grid = field.swept.grid
    gu, gv, _ = grid.gradient_at(grid.null[0] + offsets[:, 0], grid.null[1] + offsets[:, 1])
    gradient = gu[:, None] * start.first + gv[:, None] * start.second
    volts = start.phasor * complex(math.cos(omega * t0), math.sin(omega * t0))
    before = start.phasor * complex(math.cos(omega * (t0 - dt / 2)), math.sin(omega * (t0 - dt / 2)))
    place = start.point + (offsets[:, :1] + shift[0]) * start.first + (offsets[:, 1:] + shift[1]) * start.second
    x0 = place - (QM / omega**2) * volts.real * gradient
    v = velocity + (QM / omega) * before.imag * gradient
    gamma = 1 / np.sqrt(1 - (v * v).sum(1) / pt.C0**2)
    u0: Arr = gamma[:, None] * v
    return x0, u0


def shot_field(guide: TwinGuide, shot: Shot) -> tuple[Driven, TwinField]:
    """A shot's drive and field: the design's drive or the null's q, the bend's shift where the shot asks for it,
    and Earth's field as the shot gives it."""
    driven = guide.drive() if shot.q is None else guide.drive(q=shot.q)
    field = guide.field(driven, bend=True if shot.bend else None)
    b = np.asarray(_earth_vector(shot, guide), dtype=np.float64)
    return driven, replace(field, earth=b if bool(np.any(b)) else None)


def track(job: Job) -> dict[str, Any]:
    """One RF phase of a scan, tracked: the electrons launched at every energy and sense (a hold on its static
    equilibrium with its micromotion, a well at the saddle, a failure near the null), followed by integrate with the
    gas and a watch writing the entries of every electron lost, as stage 1 ran them. Plain data back: each electron's
    fate and how far it went, the entries, and the run."""
    shot = job.shot
    guide = guide_of(*job.guide)
    driven, field = shot_field(guide, shot)
    swept = field.swept
    index = _start_index(shot, len(swept.points))
    start = station_at(guide, field, index)
    omega, period = 2 * math.pi * driven.f, 1 / driven.f
    dt = period / shot.steps
    energies, senses = shot.energies, shot.senses
    count = len(energies) * len(senses)
    h = guide.solved.distance
    shift = (0.0, 0.0)
    if shot.bend:
        du, dv = swept.null_shift(swept.locate(start.point[None, :]))
        shift = (float(du[0]), float(dv[0]))
    if shot.kind == "hold":
        kappa = math.hypot(*start.kappa)
        bend = rung.valley(rung.well_patch(swept.grid, driven.f, h), math.atan2(start.kappa[1], start.kappa[0]))
        offsets = np.array(
            [bend.at(2 * e * kappa, abs(start.phasor)) for _ in senses for e in energies], dtype=np.float64
        ).reshape(count, 2)
        steps = math.ceil(shot.laps * 1.02 * guide.length / (float(pt.v_of(min(energies))) * dt))
    else:
        offsets = np.zeros((count, 2))
        steps = shot.periods * shot.steps
    aim = np.array(guide.solved.saddle) - np.array(guide.solved.null)
    aim /= np.linalg.norm(aim)
    if shot.kind in ("fail", "control"):
        rng = np.random.default_rng(np.random.SeedSequence([shot.launch, 2]))
        radius = rung.FAIL_OFFSET * h * np.sqrt(rng.uniform(size=count))
        turn = rng.uniform(0, 2 * math.pi, size=count)
        offsets = np.column_stack([radius * np.cos(turn), radius * np.sin(turn)])
    t0 = float(np.mod((job.phase - math.atan2(start.phasor.imag, start.phasor.real)) / omega, period))
    velocity = np.zeros((count, 3))
    for i, (s, e) in enumerate((s, e) for s in senses for e in energies):
        if shot.kind == "well":
            velocity[i] = float(pt.v_of(e)) * (aim[0] * start.first + aim[1] * start.second)
        else:
            velocity[i] = s * float(pt.v_of(e)) * start.tangent
    if shot.kind in ("fail", "control"):
        x0 = start.point + offsets[:, :1] * start.first + offsets[:, 1:] * start.second
        gamma = 1 / np.sqrt(1 - (velocity * velocity).sum(1) / pt.C0**2)
        u0: Arr = gamma[:, None] * velocity
    else:
        x0, u0 = launch(field, start, offsets, shift, velocity, t0, omega, dt)
    watch = Watch(field, count, (), -pt.E, job.weight, job.collector, 1, x0, losses=True)
    pressure, temperature, sigma = job.vacuum
    gas = TwinVacuum(pressure, temperature, sigma).scatter(seed=shot.launch * 16 + job.index)
    pt.integrate(x0, u0, np.full(count, pt.QM_E), dt, steps, rung._shifted(field, t0), [gas, watch])
    seen = watch.reading()
    if shot.kind == "hold":
        reached = np.abs(seen.laps)
        fate = np.where(reached >= shot.laps, 0, np.where(seen.fate == 0, 3, seen.fate))
    else:
        reached = np.where(seen.fate == 0, float(shot.periods), seen.time / period)
        fate = seen.fate
    b = field.earth if field.earth is not None else np.zeros(3)
    return {
        "fate": fate.reshape(len(senses), len(energies)).tolist(),
        "reached": reached.reshape(len(senses), len(energies)).tolist(),
        "losses": lab.plain(seen.losses),
        "steps": steps,
        "dt": dt,
        "f": driven.f,
        "volts": driven.volts,
        "q": driven.q,
        "earth": [float(x) for x in b],
        "start": index,
        "guide": guide.identity(),
    }


def read_mpc(job: MpcJob) -> dict[str, Any]:
    """One E2-13 reading (emission.read) on the shot's case, surface and laminate model, on MPC_SEEDS, at the
    emission module's own settings (1000 electrons, 30 RF periods of which 10 settle, 10 generations of which 3
    settle, 200 steps an RF period)."""
    case = mpc_cases(job.shot.design)[job.case]
    surface, gain = surface_of(job.surface)
    laminate = None if job.laminate == "absorbs" else surface
    run = MPC_SETTINGS
    found = em.read(
        case,
        surface,
        laminate=laminate,
        gain=gain,
        seeds=MPC_SEEDS,
        electrons=run["electrons"],
        periods=run["periods"],
        settle=run["settle"],
        count=run["count"],
        generation_settle=run["generation_settle"],
        steps=run["steps"],
    )
    runs = found.growth.runs
    energy = float(np.nanmean([r.energy for r in runs]))
    escaped = float(np.nanmean([r.escaped for r in runs]))
    out: dict[str, Any] = lab.plain(
        MpcReading(
            case.name,
            job.surface,
            job.laminate,
            case.volts,
            found.growth.rate,
            found.growth.spread,
            found.decay,
            found.secondaries,
            found.secondaries_spread,
            energy,
            escaped,
            any(r.lived < len(r.factors) or r.factors[-1] <= 0 for r in runs),
        )
    )
    return out


def _work(item: Job | MpcJob) -> dict[str, Any]:
    return track(item) if isinstance(item, Job) else read_mpc(item)


def job_cost(item: Job | MpcJob, cache: str | None = None) -> float:
    """A job's time on one of this machine's performance cores (s), by the costs measured while writing this
    protocol: the tracker's step 230 us and 0.7 us more a particle with the entries kept (section 8), an E2-13 reading
    40 s."""
    if isinstance(item, MpcJob):
        return MPC_SECONDS
    shot = item.shot
    guide = guide_of(item.guide[0], item.guide[1], cache)
    count = len(shot.energies) * len(shot.senses)
    if shot.kind == "hold":
        dt = 1 / (guide.tune().f * shot.steps)
        steps = math.ceil(shot.laps * 1.02 * guide.length / (float(pt.v_of(min(shot.energies))) * dt))
    else:
        steps = shot.periods * shot.steps
    return steps * (STEP_COST[0] + STEP_COST[1] * count)


STEP_COST = (230e-6, 0.7e-6)  # s: the tracker's step on this machine, with the entries kept, and a particle's more
MPC_SECONDS = 40.0  # s: an E2-13 reading on four seeds


def run_items(items: Sequence[Job | MpcJob], workers: int, cache: str | None = None) -> list[dict[str, Any]]:
    """The jobs' results in their order: in a pool of `workers` processes, longest first, or in this process."""
    if workers <= 1 or len(items) <= 1:
        return [_work(item) for item in items]
    order = sorted(range(len(items)), key=lambda k: -job_cost(items[k], cache))
    out: list[dict[str, Any]] = [{} for _ in items]
    with ProcessPoolExecutor(max_workers=min(workers, len(items))) as pool:
        futures = {k: pool.submit(_work, items[k]) for k in order}
        for k, future in futures.items():
            out[k] = future.result()
    return out


def _jobs(shot: Shot, devices: Mapping[str, Any], cache: str | None) -> list[Job]:
    guide = _guide(shot, devices, cache)
    detector, vacuum = devices["detector"], devices["vacuum"]
    where = (str(_path(shot.design)), guide.realisation, str(guide.cache))
    gas = (float(vacuum.pressure()), float(vacuum.temperature), float(vacuum.cross_section))
    return [
        Job(shot, k, phase, where, float(detector.weight), float(detector.collector), gas)
        for k, phase in enumerate(_phases(shot))
    ]


def _guide(shot: Shot, devices: Mapping[str, Any], cache: str | None) -> TwinGuide:
    """The shot's guide: the bench's own for realisation None, else the design opened at its realisation beside it,
    sharing the bench guide's cache of solves."""
    bench = devices["guide"]
    if not isinstance(bench, TwinGuide):
        raise lab.Refusal(f"the bench's guide is a {type(bench).__name__}, and this experiment runs on array16.twin")
    if shot.realisation is None:
        named = guide_of(shot.design, 0, cache).design.name
        if bench.design.name != named:
            raise lab.Refusal(f"{shot.ident} runs on {named}, and the bench's guide is {bench.design.name}")
        return bench
    return guide_of(shot.design, shot.realisation, str(bench.cache))


def _scan_of(shot: Shot, parts: Sequence[Mapping[str, Any]]) -> Scan:
    first = parts[0]
    table = Table(
        shot.energies,
        shot.senses,
        _phases(shot),
        tuple(tuple(tuple(int(x) for x in row) for row in part["fate"]) for part in parts),
        tuple(tuple(tuple(float(x) for x in row) for row in part["reached"]) for part in parts),
        shot.run,
    )
    b = first["earth"]
    return Scan(
        kind=shot.kind,
        guide=first["guide"],
        seed=shot.seed,
        table=table,
        losses=tuple(Losses.load(part["losses"]) for part in parts),
        steps=int(first["steps"]),
        dt=float(first["dt"]),
        f=float(first["f"]),
        volts=float(first["volts"]),
        q=float(first["q"]),
        earth=(float(b[0]), float(b[1]), float(b[2])),
        bend=shot.bend,
        start=int(first["start"]),
    )


def _swept_null(field: Any, i: int) -> tuple[float, float]:
    """Where a swept field vanishes in sample i's normal plane: the null's offsets from the path along the frame's
    first and second normals (m), by Newton's method on E, its derivatives by differences 1e-7 m apart."""
    point, first, second = field.points[i], field.first[i], field.second[i]
    at = np.zeros(2)
    for _ in range(8):
        probes = at + np.array([[0.0, 0.0], [1e-7, 0.0], [0.0, 1e-7]])
        e, _ = field(point + probes[:, :1] * first + probes[:, 1:] * second, 0.0)
        assert e is not None
        here = np.column_stack([e @ first, e @ second])
        jacobian = np.column_stack([here[1] - here[0], here[2] - here[0]]) / 1e-7
        at = at - np.linalg.solve(jacobian, here[0])
    return float(at[0]), float(at[1])


def measure_shift(shot: Shot, cache: str | None) -> Shift:
    """E2-2: section B's solved grid swept round the ring with its [bend] shift on and a static volt, the swept null
    found by Newton's method at SHIFT_STATIONS, its offset along the bend's axis (u in the boards' plane, v toward
    them) times R/h^2."""
    guide = guide_of(shot.design, 0, cache)
    field = guide.sweep(guide.grid(), pt.InPhase(1.0, 0.0))
    radius = ring_radius(shot)
    h = guide.design.section.geometry.distance
    axis = 0 if _way(shot) == "in-plane" else 1
    found = [_swept_null(field, i) for i in SHIFT_STATIONS]
    along = [abs(o[axis]) for o in found]
    across = [abs(o[1 - axis]) for o in found]
    coefficient = float(np.mean(along)) * radius / h
    return Shift(
        "shift",
        guide.identity(),
        shot.seed,
        _way(shot),
        radius,
        SHIFT_STATIONS,
        tuple(found),
        coefficient,
        max(across) / max(min(along), 1e-300),
    )


def measure(batch: Sequence[Shot], devices: Mapping[str, Any], workers: int, cache: str | None) -> dict[str, Any]:
    """A batch's shots: every scan's RF phases and every E2-13 reading as one job each, in the pool together."""
    items: list[Job | MpcJob] = []
    owners: list[str] = []
    for shot in batch:
        if shot.kind in ("hold", "well", "fail", "control"):
            for job in _jobs(shot, devices, cache):
                items.append(job)
                owners.append(shot.ident)
        elif shot.kind == "multipactor":
            for k in range(len(mpc_cases(shot.design))):
                for surface in SURFACES[shot.design]:
                    for laminate in LAMINATES:
                        items.append(MpcJob(shot, k, surface, laminate))
                        owners.append(shot.ident)
    results = run_items(items, workers, cache)
    out: dict[str, Any] = {}
    for shot in batch:
        parts = [r for r, who in zip(results, owners, strict=True) if who == shot.ident]
        if shot.kind == "shift":
            out[shot.ident] = measure_shift(shot, cache)
        elif shot.kind == "multipactor":
            readings = tuple(MpcReading.of(part) for part in parts)
            out[shot.ident] = Mpc("multipactor", {"cases": [c.name for c in mpc_cases(shot.design)]}, shot.seed,
                                  MPC_SEEDS, readings)  # fmt: skip
        else:
            out[shot.ident] = _scan_of(shot, parts)
    return out


def load(record: Mapping[str, Any]) -> Measurement:
    """A measurement back from the plain data it was stored as."""
    kind = record["kind"]
    if kind == "shift":
        return Shift(
            "shift",
            record["guide"],
            int(record["seed"]),
            str(record["way"]),
            float(record["radius"]),
            tuple(int(i) for i in record["stations"]),
            tuple((float(a), float(b)) for a, b in record["offsets"]),
            float(record["coefficient"]),
            float(record["cross"]),
        )
    if kind == "multipactor":
        readings = tuple(MpcReading.of(r) for r in record["readings"])
        return Mpc("multipactor", record["guide"], int(record["seed"]), tuple(record["seeds"]), readings)
    b = record["earth"]
    return Scan(
        kind=str(kind),
        guide=record["guide"],
        seed=int(record["seed"]),
        table=_table_load(record["table"]),
        losses=tuple(Losses.load(r) for r in record["losses"]),
        steps=int(record["steps"]),
        dt=float(record["dt"]),
        f=float(record["f"]),
        volts=float(record["volts"]),
        q=float(record["q"]),
        earth=(float(b[0]), float(b[1]), float(b[2])),
        bend=bool(record["bend"]),
        start=int(record["start"]),
    )


# ---- the verdicts: each function is one criterion of PROTOCOL.md, section 5 ---------------------------------------
Row = tuple[Shot, Prediction, Measurement | None]


def _ev(value: float) -> str:
    if value == 0.0:
        return "below the scan"
    return "above the scan" if math.isinf(value) else f"{value:.4g} eV"


def _within(value: float, low: float, high: float) -> bool:
    return math.isfinite(value) and low <= value <= high


def _steps(a: float, b: float, ratio: float) -> float:
    """How many scan steps apart two energies are (inf if either is off the scan)."""
    if not (0 < a < math.inf and 0 < b < math.inf):
        return math.inf
    return abs(math.log(a / b)) / math.log(ratio)


def _holds(rows: Iterable[Row], test: str | None = None) -> Iterator[tuple[Shot, Prediction, Scan]]:
    for shot, p, m in rows:
        if shot.kind == "hold" and isinstance(m, Scan) and (test is None or shot.test == test):
            yield shot, p, m


def check_motion(rows: Sequence[Row], cache: str | None = None) -> Iterator[Check]:
    """Every hold against its exact motion (Claim): the tracked limit within MOTION of the motion's."""
    for shot, p, m in _holds(rows):
        found = m.table.limit()
        said = (
            f"held to {_ev(found)} (with the wave {_ev(m.table.threshold(0))}, against it {_ev(m.table.threshold(1))}) "
            f"against the exact motion's {p.value:.4g} eV (band {p.low:.4g} to {p.high:.4g}); the static limit "
            f"{p.basis['static_eV']:.4g} eV, {found / p.basis['static_eV']:.3f} of it"
        )
        yield Check(f"{shot.test} exact motion", shot.ident, _within(found, p.low, p.high), said)


def _h1_reading(shot: Shot, p: Prediction, m: Scan) -> tuple[float, Mapping[str, float]]:
    """A hold's limit at 20 laps and the bands H1, H1a and H1b give it."""
    at = LAPS if shot.laps > LAPS else None
    found = m.table.limit(at)
    basis = p.basis
    if at is None:
        a_low, a_high = p.low, p.high
    else:
        a_low, a_high = _band(basis["motion20_limit_eV"], MOTION)
    bands = {
        "H1a": (a_low, a_high),
        "H1": (basis["H1_low_eV"], basis["H1_high_eV"]),
        "H1b": (basis["H1b_low_eV"], basis["H1b_high_eV"]),
    }
    return found, {f"{k}_{side}": v[i] for k, v in bands.items() for i, side in enumerate(("low", "high"))}


def check_h1(rows: Sequence[Row], cache: str | None = None) -> Iterator[Check]:
    """H1 against H1a and H1b: each of H1's holds, and T915-S's at q 0.49 read at 20 laps, in each hypothesis's
    band: H1 (planar paths their static limit, paths out of the plane 0.845 to 0.89 of theirs), H1a (the exact motion
    in the section, the changing bend) and H1b (every bend counted across the boards' plane)."""
    names = {"H1": "H1: the plane", "H1a": "H1a: the changing bend", "H1b": "H1b: the stiff side"}
    for shot, p, m in _holds(rows):
        if not (shot.test == "H1" or (shot.test == "H2" and shot.q is None and shot.steps == STEPS)):
            continue
        found, bands = _h1_reading(shot, p, m)
        static = p.basis["static_eV"]
        for key, name in names.items():
            low, high = bands[f"{key}_low"], bands[f"{key}_high"]
            said = f"{_ev(found)} ({found / static:.3f} of the static limit) against {low:.4g} to {high:.4g} eV"
            yield Check(name, shot.ident, _within(found, low, high), said)


def check_e2_9(rows: Sequence[Row], cache: str | None = None) -> Iterator[Check]:
    """E2-9: facing boards bent in their plane hold their static limit (stage 1's tracked fraction, PLANE), and bent
    toward them at least twice what the flat ring holds at the same q."""
    flats = {shot.q: (p, m) for shot, p, m in _holds(rows) if shot.design == FLAT and shot.test in ("H1", "E2-9")}
    stiffs = {shot.q: (shot, m) for shot, _, m in _holds(rows) if shot.design == STIFF and shot.test in ("H1", "E2-9")}
    for q, (p, m) in flats.items():
        found, static = m.table.limit(), p.basis["static_eV"]
        low, high = _band(static, PLANE)
        said = (
            f"{_ev(found)}, {found / static:.3f} of the static limit {static:.4g} eV (band {low:.4g} to {high:.4g}); "
            f"{found / p.basis['DRh_eV']:.3f} D R/h"
        )
        ident = next(s.ident for s, _, mm in _holds(rows) if mm is m)
        yield Check("E2-9 in the boards' plane", ident, _within(found, low, high), said)
        if q in stiffs:
            stiff, sm = stiffs[q]
            ratio = sm.table.limit() / found if 0 < found < math.inf else math.nan
            said = f"{_ev(sm.table.limit())} over the flat ring's {_ev(found)}: {ratio:.2f}, {TOWARD_RATIO:g} asked"
            yield Check("E2-9 toward the boards", stiff.ident, math.isfinite(ratio) and ratio >= TOWARD_RATIO, said)


def _h2(rows: Sequence[Row]) -> list[tuple[Shot, Prediction, Scan]]:
    return list(_holds(rows, "H2"))


def check_h2(rows: Sequence[Row], cache: str | None = None) -> Iterator[Check]:
    """H2: T915-S's senses within two steps, at 100 laps and at 20; the band in both senses, their edges (all lost
    from) within two steps of each other; between q 0.49 and q 0.30 the hold moving to another fraction of its static
    limit as the exact motion moves it (a resonance), not staying at one fraction (an aperture); and from 120 to 240
    steps an RF period the hold moving as the exact motion moves (the integrator's part in it)."""
    for shot, _p, m in _h2(rows):
        t, ratio = m.table, shot.ratio
        for at in (None, LAPS):
            apart = _steps(t.threshold(0, at), t.threshold(1, at), ratio)
            laps = shot.laps if at is None else at
            said = f"at {laps:g} laps the senses' thresholds {apart:.1f} steps apart, {SENSE_STEPS} allowed"
            yield Check("H2 senses", f"{shot.ident}, at {laps:g} laps", apart <= SENSE_STEPS + 1e-9, said)
        edge = [t.edge(0), t.edge(1)]
        apart = _steps(edge[0], edge[1], ratio)
        said = f"every electron lost from {_ev(edge[0])} with the wave and {_ev(edge[1])} against it: {apart:.1f} steps"
        yield Check("H2 the band in both senses", shot.ident, apart <= BAND_STEPS + 1e-9, said)
    by = {(s.q, s.steps): (s, p, m) for s, p, m in _h2(rows)}
    pairs = ((None, STEPS), (0.3, STEPS), "H2 a resonance, not an aperture"), ((None, STEPS), (None, 2 * STEPS),
             "H2 the integrator's part")  # fmt: skip
    for one, two, name in pairs:
        if one not in by or two not in by:
            continue
        (_, pa, ma), (shot, pb, mb) = by[one], by[two]
        found = [m.table.limit() / p.basis["static_eV"] for p, m in ((pa, ma), (pb, mb))]
        model = [p.basis["motion_fraction"] for p in (pa, pb)]
        moved = _steps(found[1], found[0], shot.ratio) * (1 if found[1] >= found[0] else -1)
        expect = _steps(model[1], model[0], shot.ratio) * (1 if model[1] >= model[0] else -1)
        ok = math.isfinite(moved) and abs(moved - expect) <= BAND_STEPS
        if name.startswith("H2 a resonance"):
            ok = ok and abs(moved) > BAND_STEPS
        said = (
            f"{found[0]:.3f} and {found[1]:.3f} of the static limit: moved {moved:+.1f} steps, the exact motion "
            f"{expect:+.1f} ({model[0]:.3f} to {model[1]:.3f})"
        )
        yield Check(name, shot.ident, ok, said)


def check_h3(rows: Sequence[Row], cache: str | None = None) -> Iterator[Check]:
    """H3: on a path out of the plane with lobes, every electron the hold lost to the guide (struck or left) lost in a
    lobe or less than LOBE_REACH lap past one in its direction of travel."""
    for shot, _p, m in _holds(rows):
        if _planar(shot.design) or math.hypot(*m.earth) > 2 * EARTH:  # a strong field sets its own places
            continue
        stretches = lobes(frame_of(shot.design, _realisation(shot), cache)[3])
        if not stretches:
            continue
        worst, count, inside = 0.0, 0, 0
        for losses in m.losses:
            for k in range(len(losses)):
                if losses.how[k] not in (STRUCK, LEFT):
                    continue
                sense = m.table.senses[int(losses.index[k]) // len(m.table.energies)]
                gone = past(float(losses.s[k]) / losses.length % 1.0, sense, stretches)
                count += 1
                inside += gone <= LOBE_REACH
                worst = max(worst, gone)
        said = f"{inside} of {count} lost within a lobe or {LOBE_REACH} lap past one; the farthest {worst:.3f} lap past"
        yield Check("H3 lost in a lobe", shot.ident, count > 0 and inside == count, said)


def check_e2_15(rows: Sequence[Row], cache: str | None = None) -> Iterator[Check]:
    """E2-15: at 50 uT across T915-S at its worst orientation the senses part by at most the static balance's split
    and two steps; at ten times the field by half to twice the static balance's split; without the field (the bench's
    T915-S with none) within two steps."""
    for shot, p, m in _holds(rows, "E2-15"):
        t = m.table
        plus, minus = t.threshold(0), t.threshold(1)
        split = math.log(minus / plus) if 0 < plus < math.inf and 0 < minus < math.inf else math.nan
        if "field_split" not in p.basis:
            apart = _steps(plus, minus, shot.ratio)
            said = f"no field: the senses {apart:.1f} steps apart, {SENSE_STEPS} allowed"
            yield Check("E2-15 no field, no split", shot.ident, apart <= SENSE_STEPS + 1e-9, said)
            continue
        expect = p.basis["field_split"]
        strong = isinstance(shot.earth, tuple) and math.hypot(*shot.earth) > 2 * EARTH
        if strong:
            ok = math.isfinite(split) and split * expect > 0 and SPLIT_BAND[0] <= split / expect <= SPLIT_BAND[1]
            said = f"ln(W-/W+) {split:+.4f} against the static balance's {expect:+.4f} (half to twice)"
            yield Check("E2-15 the field's split", shot.ident, ok, said)
        else:
            room = abs(expect) + SENSE_STEPS * math.log(shot.ratio)
            ok = math.isfinite(split) and abs(split) <= room + 1e-12
            said = f"ln(W-/W+) {split:+.4f}, the static balance's {expect:+.4f} and two steps allowed: {room:.4f}"
            yield Check("E2-15 Earth's split", shot.ident, ok, said)


def check_e2_2(rows: Sequence[Row], cache: str | None = None) -> Iterator[Check]:
    """E2-2: the twin's shift of the null within SHIFT_TOL of the axisymmetric reference, along the bend's axis
    alone; and T915-S's hold with the shift on within two steps of the same launch with it off, in each sense."""
    for shot, p, m in rows:
        if shot.kind == "shift" and isinstance(m, Shift):
            ok = _within(m.coefficient, p.low, p.high) and m.cross < 1e-3
            said = (
                f"{m.coefficient:.6f} h^2/R against the reference's {p.value:.6f} (band {p.low:.6f} to {p.high:.6f}); "
                f"across it {m.cross:.1e} of that"
            )
            yield Check("E2-2 the bend's shift", shot.ident, ok, said)
    pair = {shot.bend: (shot, m) for shot, _, m in _holds(rows, "E2-2")}
    if True in pair and False in pair:
        (on, a), (_, b) = pair[True], pair[False]
        apart = max(_steps(a.table.threshold(k), b.table.threshold(k), on.ratio) for k in range(2))
        said = (
            f"the senses' thresholds with the shift and without at most {apart:.1f} steps apart, {SENSE_STEPS} allowed"
        )
        yield Check("E2-2 the shift on beside off", on.ident, apart <= SENSE_STEPS + 1e-9, said)


def check_h4(rows: Sequence[Row], cache: str | None = None) -> Iterator[Check]:
    """H4: R1 on the bench's realisation keeps the static balance its band and field give, each sense within PLANE of
    its own; with Earth's field removed the senses together; its well's first loss the exact motion's at its band."""
    for shot, p, m in rows:
        if shot.test != "H4" or not isinstance(m, Scan):
            continue
        t = m.table
        if shot.kind == "well":
            first = t.first(0)
            said = (
                f"first lost at {_ev(first)} against the exact motion's {p.value:.4g} eV (band {p.low:.4g} to "
                f"{p.high:.4g})"
            )
            yield Check("H4 well", shot.ident, _within(first, p.low, p.high), said)
            continue
        expect = (
            (p.basis["field_with_eV"], p.basis["field_against_eV"])
            if "field_with_eV" in p.basis
            else (p.basis["static_eV"], p.basis["static_eV"])
        )
        for row, name in ((0, "with the wave"), (1, "against it")):
            found = t.threshold(row)
            low, high = _band(expect[row], PLANE)
            said = (
                f"{name} {_ev(found)} against the static balance's {expect[row]:.4g} eV (band {low:.4g} to {high:.4g})"
            )
            yield Check("H4 the static balance", f"{shot.ident}, {name}", _within(found, low, high), said)
        if "field_with_eV" not in p.basis:
            apart = _steps(t.threshold(0), t.threshold(1), shot.ratio)
            said = f"the senses {apart:.1f} steps apart, {SENSE_STEPS} allowed"
            yield Check("H4 no field, no split", shot.ident, apart <= SENSE_STEPS + 1e-9, said)


def _decays(r: MpcReading) -> bool:
    return (r.growth == -math.inf or r.growth + SPREADS * r.spread < 0) and r.secondaries < 1


def _grows(r: MpcReading) -> bool:
    return math.isfinite(r.growth) and r.growth - SPREADS * r.spread > 0 and r.secondaries > 1


def _said(r: MpcReading) -> str:
    growth = "dies" if r.growth == -math.inf else f"{r.growth:+.3f} +- {r.spread:.3f} Np a period"
    return f"{r.case}, {r.surface}, laminate {r.laminate}: {growth}, effective yield {r.secondaries:.3f}"


def check_e2_13(rows: Sequence[Row], cache: str | None = None) -> Iterator[Check]:
    """E2-13: T915-S and the split lever decay on every surface and laminate model; R1's rails grow no faster than the
    chart's bound, and die (Hypothesis); section B at 160 V grows (the control); and on facing strips at q 0.5 the
    population decays at 1.5 eV and grows at 3 eV in both laminate models, on each surface."""
    for shot, p, m in rows:
        if shot.kind != "multipactor" or not isinstance(m, Mpc):
            continue
        for r in m.readings:
            ident = f"{shot.ident}: {r.surface}, laminate {r.laminate}"
            if shot.design in ("t915s", "lever"):
                yield Check("E2-13 decays", ident + f", {r.case}", _decays(r), _said(r))
            elif shot.design == "r1":
                bound = r.growth <= p.basis["chart_growth"]
                yield Check("E2-13 R1's rails within the chart's bound", ident, bound, _said(r))
                yield Check("E2-13 R1's rails die (Hypothesis)", ident, _decays(r), _said(r))
            elif shot.design == "overdriven":
                yield Check("E2-13 the control grows", ident, _grows(r), _said(r))
        if shot.design == "onset":
            names = [c.name for c in mpc_cases("onset")]
            for surface in SURFACES["onset"]:
                low = [r for r in m.readings if r.surface == surface and r.case == names[0]]
                high = [r for r in m.readings if r.surface == surface and r.case == names[-1]]
                ok = bool(low) and bool(high) and all(_decays(r) for r in low) and all(_grows(r) for r in high)
                said = "; ".join(_said(r) for r in low + high)
                yield Check("E2-13 the onset between 1.5 and 3 eV", f"{shot.ident}: {surface}", ok, said)


def check_fail(rows: Sequence[Row], cache: str | None = None) -> Iterator[Check]:
    """The control and the case that must fail: the control's launch at the design's volts loses none to the guide
    and holds some; the must-fail shot holds none, beside a control that held; and every electron lost has its
    entry in the records."""
    control: bool | None = None
    for shot, _, m in rows:
        if shot.kind == "control" and isinstance(m, Scan):
            held, lost = m.table.counts()
            control = int(lost.sum()) == 0 and int(held.sum()) > 0
            said = (
                f"q {m.q:.2f} at {m.volts:.1f} V: {int(held.sum())} of {m.table.launched()} held, {int(lost.sum())} "
                "lost"
            )
            yield Check("must fail control", shot.ident, control, said)
    beside = {None: "no control on file", True: "the control held", False: "the control did not hold"}[control]
    for shot, _, m in rows:
        if shot.kind != "fail" or not isinstance(m, Scan):
            continue
        held, lost = m.table.counts()
        entries = sum(int(np.isin(losses.how, (STRUCK, LEFT)).sum()) for losses in m.losses)
        said = (
            f"q {m.q:.2f} at {m.volts:.1f} V: {int(held.sum())} of {m.table.launched()} held, {int(lost.sum())} lost "
            f"with {entries} entries; {beside}"
        )
        ok = int(held.sum()) == 0 and control is True and entries == int(lost.sum()) > 0
        yield Check("must fail", shot.ident, ok, said)


CHECKS: tuple[Callable[[Sequence[Row], str | None], Iterator[Check]], ...] = (
    check_motion,
    check_h1,
    check_e2_9,
    check_h2,
    check_h3,
    check_e2_15,
    check_e2_2,
    check_h4,
    check_e2_13,
    check_fail,
)
TESTS = ("H1", "E2-9", "H2", "E2-15", "E2-2", "H4", "E2-13", "must fail")


# ---- the run's time -----------------------------------------------------------------------------------------------
def batches(todo: Sequence[Shot]) -> list[list[Shot]]:
    """The shots in batches of one test each, in the protocol's order: a batch's jobs share the pool, and its
    records are written when it ends."""
    out: dict[str, list[Shot]] = {}
    for shot in todo:
        out.setdefault(shot.test, []).append(shot)
    return [out[test] for test in TESTS if test in out]


def _makespan(costs: Sequence[float], workers: int, pace: float = 1.0) -> float:
    """The longest-first schedule's span of jobs (s, each its time alone) on `workers` processes running at `pace`."""
    loads = [0.0] * max(workers, 1)
    for cost in sorted(costs, reverse=True):
        k = min(range(len(loads)), key=lambda i: loads[i])
        loads[k] += cost / pace
    return max(loads) if costs else 0.0


PACE = 0.83  # ten trackers at once each step at 0.83 of one alone on this machine (370 against 308 us a step)


def estimate(cache: str | None = None, workers: int | None = None) -> dict[str, Any]:
    """The run's time on this machine (s): each shot's work on one core alone (its jobs' job_cost), and each batch's
    span on the bench's workers, longest first, each worker at PACE when they run together."""
    count = workers or lab.Bench.local().workers
    shots_cost: dict[str, float] = {}
    spans: dict[str, float] = {}
    for batch in batches(shots()):
        costs: list[float] = []
        for shot in batch:
            jobs: list[float] = []
            if shot.kind in ("hold", "well", "fail", "control"):
                guide = guide_of(shot.design, _realisation(shot), cache)
                where = (str(_path(shot.design)), guide.realisation, str(guide.cache))
                jobs = [job_cost(Job(shot, k, ph, where, 1.0, 1e-12, (1e-6, 295.0, 1e-19)), cache) for k, ph in
                        enumerate(_phases(shot))]  # fmt: skip
            elif shot.kind == "multipactor":
                cases = len(ONSET_DEPTHS) if shot.design == "onset" else 1
                jobs = [MPC_SECONDS] * (cases * len(SURFACES[shot.design]) * len(LAMINATES))
            shots_cost[shot.ident] = sig(sum(jobs), 4)
            costs += jobs
        spans[batch[0].test] = sig(_makespan(costs, count, PACE if count > 1 else 1.0), 4)
    return {
        "workers": count,
        "cpu_s": sig(sum(shots_cost.values()), 4),
        "wall_s": sig(sum(spans.values()), 4),
        "batches_s": spans,
        "shots_s": shots_cost,
    }


# ---- what the predictions file says of the apparatus -------------------------------------------------------------
def describe(design: str, realisation: int, cache: str | None) -> dict[str, Any]:
    """A design as the theory sees it at a realisation: its file and solve, its ring and drive, its static limit and
    where, and its bound with every bend counted across the boards' plane."""
    th = theory(design, realisation, None, cache)
    guide, driven = th.guide, th.driven
    ku = frame_of(design, realisation, cache)[3]
    out: dict[str, Any] = {
        "identity": guide.identity(),
        "ring": {"length_mm": sig(guide.length * 1e3), "f_MHz": sig(driven.f / 1e6)},
        "drive": {"volts": sig(driven.volts), "q": sig(driven.q), "weakest_eV": sig(driven.depth_ev)},
        "static": {
            "limit_eV": sig(th.static),
            "at": sig(th.where / len(th.limits), 4),
            "across_eV": sig(across_bound(th, ku)),
        },
        "lobes": [[sig(a, 4), sig(b, 4)] for a, b in lobes(ku)],
    }
    if guide.curve is not None:
        out["ring"]["bend_mm"] = sig(guide.curve.bend_radius * 1e3)
    if guide.framing is not None:
        out["ring"]["inplane_radius_mm"] = sig(guide.framing.inplane_radius * 1e3)
    b = guide.drawn.earth
    if bool(np.any(b)):
        out["earth_uT"] = [sig(float(x) * 1e6) for x in b]
    return out


def _value(p: Prediction) -> str:
    if p.unit in ("held", "lost"):
        return "0"
    return f"{p.value:.6g} {p.unit}"


def _range(p: Prediction) -> str:
    if p.unit == "held":
        return "none held"
    if p.unit == "lost":
        return "none lost"
    if p.unit == "Np":
        if p.low == -math.inf:
            return "decays" if p.high <= 0 else f"at most {p.high:+.3f} Np a period"
        return "grows"
    return f"{p.low:.6g} to {p.high:.6g} {p.unit}"


def predictions_markdown(experiment: TwinRung2) -> str:
    """The predictions as tables, for the protocol's reader."""
    lines = [
        "<!-- Written by `experiments/twin_rung2.py predict`. Theory only: no electron was tracked to make it. -->",
        "",
        "# Experiment 2, stage 2a: the predictions",
        "",
        "Each hold is predicted by the exact motion in the section along its path, on the shot's own energies,",
        "senses, RF phases, launch and laps (`twin_rung2.exact_motion`), with stage 1's static balance, the bound with",
        "every bend counted across the boards' plane (H1b's) and the static balance with a uniform field beside it; a",
        "well by the exact motion on the axis at the realisation's band; the bend's shift by the axisymmetric",
        "reference;",
        "multipactor by the level-1 chart and its one-dimensional operator. PROTOCOL.md section 2 gives each claim.",
        "",
    ]
    rows = []
    for shot in experiment.shots():
        p = experiment.predict(shot)
        extra = ""
        if shot.test == "H1" or (shot.test == "H2" and shot.q is None):
            b = p.basis
            extra = (
                f"H1 {b['H1_low_eV']:.4g} to {b['H1_high_eV']:.4g}; H1b {b['H1b_low_eV']:.4g} to "
                f"{b['H1b_high_eV']:.4g} eV"
            )
        rows.append([shot.ident, p.claim, _value(p), _range(p), extra])
    lines += lab.table(["Shot", "Prediction", "Theory", "Band", "H1 and H1b"], rows)
    lines += ["", "## The numbers each rests on", ""]
    for shot in experiment.shots():
        p = experiment.predict(shot)
        said = ", ".join(f"{key} {value:g}" for key, value in p.basis.items())
        lines.append(f"- {shot.ident}: {said}")
    time = estimate(experiment.cache, experiment.workers)
    lines += [
        "",
        "## The run's time",
        "",
        f"On {time['workers']} workers, by each job's cost on one of this machine's performance cores (`estimate`):",
        f"{time['cpu_s'] / 3600:.2f} hours of work, {time['wall_s'] / 3600:.2f} hours from start to end.",
        "",
        *lab.table(["Batch", "Span (s)"], [[k, f"{v:.0f}"] for k, v in time["batches_s"].items()]),
    ]
    return "\n".join(lines) + "\n"


class TwinRung2(lab.Experiment[Shot, Prediction, Measurement]):
    """Stage 2a as array16.lab runs it, on the twin bench (benches/twin_r2.toml): its guide, detector and vacuum.
    `cache` is the folder of solved sections and kept exact motions (None: array16.twin's); `workers` the processes the
    tracking and the exact motions run in (None: all this machine's cores but two, the bench's default)."""

    out = OUT
    needs: Mapping[str, Sequence[str]] = {
        "guide": ("tune", "drive", "field", "identity"),
        "detector": ("watch", "identity"),
        "vacuum": ("scatter", "free_path", "identity"),
    }
    registered = (
        OUT / "PROTOCOL.md",
        OUT / "predictions.md",
        lab.ROOT / RING_IN,
        lab.ROOT / RING_TOWARD,
        Path(__file__).resolve(),
        lab.ROOT / "experiments" / "twin_rung1.py",
        *(lab.ROOT / d for d in (T915S, FLAT, STIFF, ELLIPSE, HELIX_010, HELIX, HELIX_030, R1)),
        lab.BENCHES / f"{BENCH}.toml",
        lab.ROOT / "array16" / "twin",
        lab.ROOT / "array16" / "guide",
        lab.ROOT / "array16" / "particles",
        lab.ROOT / "array16" / "lab.py",
    )

    def __init__(self, cache: str | Path | None = None, workers: int | None = None) -> None:
        self.cache = None if cache is None else str(cache)
        self.workers = workers if workers is not None else lab.Bench.local().workers
        self._prepared = False

    def shots(self) -> Sequence[Shot]:
        return shots()

    def ident(self, shot: Shot) -> str:
        return shot.ident

    def prepare(self) -> None:
        """Compute in the pool every hold's exact motion not yet kept."""
        if not self._prepared:
            specs = [motion_spec(shot, self.cache) for shot in self.shots() if shot.kind == "hold"]
            motions(specs, self.cache, self.workers)
            self._prepared = True

    def predict(self, shot: Shot) -> Prediction:
        self.prepare()
        return predict(shot, self.cache)

    def batches(self, shots: Sequence[Shot]) -> Sequence[Sequence[Shot]]:
        return batches(shots)

    def measure(self, batch: Sequence[Shot], devices: Mapping[str, Any]) -> Mapping[str, Measurement]:
        found: Mapping[str, Measurement] = measure(batch, devices, self.workers, self.cache)
        return found

    def load(self, record: Mapping[str, Any]) -> Measurement:
        return load(record)

    def judge(self, rows: Sequence[Row]) -> Sequence[Check]:
        return [check for criterion in CHECKS for check in criterion(rows, self.cache)]

    def apparatus(self) -> Mapping[str, Any]:
        designs = {d: describe(d, 0, self.cache) for d in (T915S, FLAT, STIFF, ELLIPSE, HELIX_010, HELIX, HELIX_030)}
        designs["designs/t915s.toml at realisation 1, the bench's guide"] = describe(T915S, 1, self.cache)
        designs["designs/r1_rail_ring.toml at realisation 1"] = describe(R1, 1, self.cache)
        criteria = {
            "steps_per_period": STEPS,
            "laps": LAPS,
            "long_laps": LONG,
            "fine_ratio": FINE,
            "coarse_ratio": COARSE,
            "phases": PHASES,
            "window": rung.WINDOW,
            "motion_band": MOTION,
            "plane_band": PLANE,
            "out_of_plane_band": OUT_OF_PLANE,
            "across_band": ACROSS,
            "sense_steps": SENSE_STEPS,
            "band_steps": BAND_STEPS,
            "onset_count": ONSET,
            "lobe_reach": LOBE_REACH,
            "split_band": SPLIT_BAND,
            "toward_ratio": TOWARD_RATIO,
            "shift_tolerance": SHIFT_TOL,
            "well_steps": WELL_STEPS,
            "spreads": SPREADS,
            "mpc_seeds": MPC_SEEDS,
            "worst_orientation": WORST,
            "master_seed": MASTER,
        }
        return {"bench": f"benches/{BENCH}.toml", "designs": designs, "criteria": criteria}

    def summary(self, checks: Sequence[Check]) -> Mapping[str, Sequence[Check]]:
        by_ident = {shot.ident: shot for shot in shots()}

        def test_of(ident: str) -> str | None:
            shot = by_ident.get(ident) or next((s for i, s in by_ident.items() if ident.startswith(i)), None)
            return None if shot is None else shot.test

        return {f"{test}, every check": [c for c in checks if test_of(c.ident) == test] for test in TESTS}


def _bench_workers(args: Sequence[str]) -> int | None:
    """The workers `run` would use: --workers, else the bench file's (all cores but two where it gives none)."""
    if "--workers" in args:
        return int(args[args.index("--workers") + 1])
    name = args[args.index("--bench") + 1] if "--bench" in args else "local"
    try:
        return lab.Bench.named(name).workers
    except lab.Refusal:
        return None


def main(argv: Sequence[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    experiment = TwinRung2(workers=_bench_workers(args) if args[:1] == ["run"] else None)
    if args[:1] == ["estimate"]:
        print(json.dumps(estimate(experiment.cache, experiment.workers), indent=1))
        return 0
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
