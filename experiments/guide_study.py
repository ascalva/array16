"""Experiment 1: does a guide field hold an electron because its axis is knotted, or because it is twisted?

A winding on a closed curve makes a magnetic field that follows the curve. The first-order theory
(array16.particles.GuideOptics; laws GDE-1 to GDE-7 in docs/LAWS.md) says what that field does to an electron from
three numbers per curve: the turn of the transverse plane per lap, the net curvature drift and the wander of the
field line. The turn is the curve's writhe, less whole turns: zero for any planar ring, and not zero for a knot or
for an unknotted curve that leaves the plane. This script computes the theory's predictions, tracks electrons in
the windings' Biot-Savart fields, and judges the one against the other.

    uv run python experiments/guide_study.py predict     # theory only: guide/predictions.{json,md}
    uv run python experiments/guide_study.py run         # tracking, on every core but two: guide/runs.jsonl
    uv run python experiments/guide_study.py report      # the verdicts: guide/results.{json,md}

The three steps are array16.lab's, shared by every experiment: `run` refuses until the predictions are committed,
keeps each measurement with what made it, and picks up where it stopped.

The protocol, with what counts as failure, is experiments/guide/PROTOCOL.md. It and the predictions were committed
before the first `run`. Models: filament windings, one electron at a time, no collisions, no space charge.
"""

from __future__ import annotations

import cmath
import json
import math
import sys
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from functools import cache, cached_property
from itertools import groupby
from pathlib import Path
from typing import Any, Literal

import numpy as np

from array16 import lab
from array16 import particles as pt
from array16.lab import Check
from array16.particles import ClosedCurve, Coil, GuideOptics, Launch, OrbitFit, Winding

OUT = Path(__file__).resolve().parent / "guide"
TUBE = 0.011  # winding radius (m), as in the first run
WALL = 0.9 * TUBE  # an electron further than this from the curve counts as lost
FIELD = 0.020  # the guide field (T)
LAPS = 20
LADDER = (0.5e-3, 1.0e-3, 1.5e-3, 2.0e-3, 3.0e-3)  # gyroradii at FIELD (m)

LaunchKind = Literal["axis", "matched", "first-run"]
Group = Literal["main", "scaling", "systematics", "first-run"]

# ---- what "agrees" means: fixed before the first run (PROTOCOL.md, section 5) ------------------------------------
INSIDE = 0.70  # a predicted reach below this fraction of the wall is "well inside"
BAND = 0.15  # predicted reaches within this fraction of the wall are reported, and not judged held or lost
REACH_TOL, FLOOR = 0.15, 1e-4  # tracked against predicted reach: 15 % and 0.1 mm
CENTRE_TOL = 0.15  # the fitted centre against the closed orbit: 15 % of its offset, and FLOOR
TURN_TOL = math.radians(8.0)  # the fitted turn per lap against the holonomy
MATCHED_TOL = 0.25  # a matched launch stays within this fraction of the closed orbit's offset, and FLOOR
EXIT_TOL, EXIT_FLOOR = 0.15, 0.3  # the lap of loss on a ring: 15 %, or 0.3 lap
REPEAT_TOL, REPEAT_FLOOR = 0.03, 2e-5  # a repeat with a finer instrument or at twice the field: 3 % and 20 um
REPEAT_TURN = math.radians(2.0)
FIT_POINTS, FIT_RADIUS = 8, 3e-4  # a circle is fitted to at least 8 crossings more than 0.3 mm across


# ---- the apparatus -----------------------------------------------------------------------------------------------
def _trefoil(points: int) -> ClosedCurve:
    return ClosedCurve.torus_knot(0.060, 0.028, points)


SHAPES: dict[str, Callable[[int], ClosedCurve]] = {
    "ring": lambda points: ClosedCurve.circle(_trefoil(points).length / (2 * math.pi), points),
    "trefoil": _trefoil,
    "shallow trefoil": lambda points: ClosedCurve.torus_knot(0.060, 0.016, points),
    "helix": lambda points: ClosedCurve.torus_knot(0.110, 0.030, points, 1, 3),
    "figure-eight": lambda points: ClosedCurve.figure_eight(0.14, 0.07, 0.04, points),
}
KNOTTED = {"ring": False, "trefoil": True, "shallow trefoil": True, "helix": False, "figure-eight": False}


@dataclass(frozen=True)
class Guide:
    """One apparatus: a closed curve and the winding on it."""

    name: str
    shape: str
    winding: Winding = "return"
    points: int = 600
    per_turn: int = 36
    turns: int = 180  # per layer

    @cached_property
    def curve(self) -> ClosedCurve:
        return SHAPES[self.shape](self.points)

    @cached_property
    def coil(self) -> Coil:
        return Coil.wound(self.curve, TUBE, self.turns, self.winding, self.per_turn)

    @property
    def optics(self) -> GuideOptics:
        return self.coil.optics


GUIDES: dict[str, Guide] = {
    g.name: g
    for g in (
        Guide("ring", "ring"),
        Guide("trefoil", "trefoil"),
        Guide("trefoil-shallow", "shallow trefoil"),
        Guide("helix", "helix"),
        Guide("figure-eight", "figure-eight"),
        Guide("trefoil-fine-curve", "trefoil", points=1200),
        Guide("trefoil-fine-winding", "trefoil", per_turn=72),
        Guide("ring-single", "ring", winding="single"),
        Guide("trefoil-single", "trefoil", winding="single"),
    )
}
MAIN = ("ring", "trefoil", "trefoil-shallow", "helix", "figure-eight")


@dataclass(frozen=True)
class Shot:
    """One electron: the guide, its gyroradius at the guide field (negative: sent against the field), where it
    starts, and how it is tracked."""

    guide: str
    rho: float
    launch: LaunchKind
    field: float = FIELD
    laps: int = LAPS
    steps_per_turn: int = 20
    group: Group = "main"

    @property
    def ident(self) -> str:
        return (
            f"{self.guide} {self.launch} rho={self.rho * 1e3:+.4f}mm B={self.field * 1e3:g}mT "
            f"laps={self.laps} steps={self.steps_per_turn}"
        )

    @property
    def w_ev(self) -> float:
        return pt.energy_of(self.field * abs(self.rho))

    @property
    def speed_key(self) -> tuple[str, float, float, int, int]:
        """Shots that share this go round together: one guide, one speed, one time step."""
        return (self.guide, self.field, abs(self.rho), self.laps, self.steps_per_turn)


def _both_ways(guide: str, rho: float, kinds: Iterable[LaunchKind], **how: Any) -> Iterator[Shot]:
    return (Shot(guide, sign * rho, kind, **how) for kind in kinds for sign in (1.0, -1.0))


@cache
def shots() -> tuple[Shot, ...]:
    """Every shot of the experiment, in the order of the protocol's section 4."""
    out: list[Shot] = []
    for name in MAIN:
        kinds: tuple[LaunchKind, ...] = ("axis",) if name == "ring" else ("axis", "matched")
        ladder = (0.25e-3, *LADDER) if name == "trefoil-shallow" else LADDER
        for rho in ladder:
            out += _both_ways(name, rho, kinds)
    for name in ("trefoil", "helix"):  # the same gyroradii at twice the field: four times the energy
        for rho in (0.5e-3, 1.0e-3):
            out += _both_ways(name, rho, ("axis", "matched"), field=2 * FIELD, group="scaling")
    for kind in ("axis", "matched"):  # the same shot with each part of the instrument made finer in turn
        out += [
            Shot("trefoil", 1.0e-3, kind, steps_per_turn=40, group="systematics"),
            Shot("trefoil", 1.0e-3, kind, laps=60, group="systematics"),
            Shot("trefoil-fine-curve", 1.0e-3, kind, group="systematics"),
            Shot("trefoil-fine-winding", 1.0e-3, kind, group="systematics"),
        ]
    for name in ("trefoil-single", "ring-single"):  # the first run's cases, on its single-wire winding
        for field in (0.005, 0.010, 0.020):
            for w_ev in (100.0, 300.0, 1000.0, 3000.0):
                out.append(Shot(name, pt.rigidity(w_ev) / field, "first-run", field=field, group="first-run"))
    return tuple(out)


def start(guide: Guide, shot: Shot) -> complex | None:
    """Where the shot starts, on the axes at the curve's first point. None if there is nowhere to start it: a
    matched launch needs a closed orbit, inside the wall."""
    if shot.launch == "axis":
        return 0j
    if shot.launch == "first-run":
        return complex(2e-3, 0.0)  # 2 mm along the winding frame's first normal, as the first run launched
    star = complex(guide.optics.closed_orbit(np.array([shot.rho]))[0])
    return star if abs(star) < WALL else None


# ---- the theory's side -------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Prediction:
    """What the first-order theory says a shot does."""

    launched: bool
    held: bool  # inside the wall for all its laps
    reach: float | None  # largest distance from the curve over its laps (m); None if not launched
    exit_lap: float | None  # the lap at which it reaches the wall; None if it never does
    centre: complex | None  # the closed orbit; None on a curve that does not turn the plane
    radius: float | None  # the launch's distance from the closed orbit (m)
    turn: float | None  # the turn per lap as the electron sees it (rad)
    step: complex  # what a lap adds before the plane turns (m)


def predict(guide: Guide, shot: Shot) -> Prediction:
    opt, z0 = guide.optics, start(guide, shot)
    step = complex(opt.step(np.array([shot.rho]))[0])
    centre = complex(opt.closed_orbit(np.array([shot.rho]))[0]) if opt.closes else None
    turn = (opt.alpha if shot.rho >= 0 else -opt.alpha) if opt.closes else None
    if z0 is None:
        return Prediction(False, False, None, None, centre, None, turn, step)
    exit_lap = opt.exit_lap(shot.rho, z0, WALL, shot.laps)
    reach = float(opt.distances(shot.rho, z0, shot.laps).max())
    radius = abs(z0 - centre) if centre is not None else None
    return Prediction(True, exit_lap is None, reach, exit_lap, centre, radius, turn, step)


# ---- the tracking's side -----------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Measurement:
    """What tracking in the winding's field found."""

    held: bool
    laps: float  # laps completed, or the lap of loss
    reach: float  # largest distance from the curve while held (m)
    section: tuple[complex, ...]  # its offset at the start, then after each lap, on the start's axes

    def fit(self) -> OrbitFit | None:
        """The circle through its crossings, if there are enough of them and they are not a point."""
        if len(self.section) < FIT_POINTS:
            return None
        z = np.array(self.section, dtype=np.complex128)
        return pt.orbit_circle(z) if float(np.ptp(z.real) + np.ptp(z.imag)) >= FIT_RADIUS else None


def track(guide: Guide, batch: Sequence[Shot]) -> dict[str, Measurement]:
    """Track the shots of one speed together. A shot with nowhere to start is left out."""
    first = batch[0]
    starts = {shot.ident: start(guide, shot) for shot in batch}
    live = [shot for shot in batch if starts[shot.ident] is not None]
    if not live:
        return {}
    z0 = np.array([starts[shot.ident] for shot in live], dtype=np.complex128)
    launch = Launch(
        np.full(len(live), first.w_ev),
        guide.curve.transverse(z0),
        np.array([math.copysign(1.0, shot.rho) for shot in live]),
    )
    field = pt.scaled(guide.coil.field, first.field / guide.coil.per_amp)
    run = pt.track_guide(
        guide.curve, field, first.field, launch, WALL, first.laps + 0.02, steps_per_turn=first.steps_per_turn
    )
    return {
        shot.ident: Measurement(
            bool(run.held[i]),
            float(abs(run.laps[i])),
            float(run.reach[i]),
            (complex(z0[i]), *(complex(z) for z in run.section[i])),
        )
        for i, shot in enumerate(live)
    }


# ---- what the predictions file says of the apparatus -------------------------------------------------------------
def describe(guide: Guide) -> dict[str, Any]:
    """The numbers the theory needs from a guide, and the coil's own."""
    curve, coil, opt = guide.curve, guide.coil, guide.optics
    unit = np.array([1e-3])
    return {
        "shape": guide.shape,
        "knotted": KNOTTED[guide.shape],
        "winding": guide.winding,
        "length_mm": curve.length * 1e3,
        "bend_mm": curve.bend_radius * 1e3,
        "writhe": curve.writhe,
        "turn_deg": math.degrees(opt.alpha),
        "drift_per_rho": abs(opt.g[-1]),
        "wander_mm": abs(opt.h[-1]) * 1e3,
        # on a loop just outside both layers: a wider one would cut through the next strand of a tight knot
        "net_turns_along_curve": coil.linked_current(1.5 * TUBE),
        "mT_per_A": coil.per_amp * 1e3,
        "closed_orbit_per_rho": abs(opt.closed_orbit(unit)[0] - opt.closed_orbit(0 * unit)[0]) / 1e-3,
        "matched_reach_at_1mm_mm": float(opt.matched_excursion(unit)[0]) * 1e3,
        "axis_reach_at_1mm_mm": float(opt.distances(1e-3, 0j, LAPS).max()) * 1e3,
    }


def _mm(value: float | None, digits: int = 2) -> str:
    return "" if value is None or not math.isfinite(value) else f"{value * 1e3:.{digits}f}"


def _polar(value: complex | None) -> str:
    if value is None or not cmath.isfinite(value):
        return ""
    return f"{abs(value) * 1e3:.2f} at {math.degrees(cmath.phase(value)):+.0f}°"


GUIDE_COLUMNS = (
    "Guide",
    "Shape",
    "Knotted",
    "Winding",
    "Length (mm)",
    "Tightest bend (mm)",
    "Writhe",
    "Turn per lap",
    "Drift per lap (per ρ)",
    "Field-line wander per lap (mm)",
    "Net turns along the curve",
    "Closed orbit's offset (per ρ)",
    "Matched reach at ρ = 1 mm (mm)",
    "On-axis reach at ρ = 1 mm, 20 laps (mm)",
)
SHOT_COLUMNS = (
    "Shot",
    "Energy (eV)",
    "Predicted",
    "Reach (mm)",
    "Lap of loss",
    "Closed orbit (mm)",
    "Launch to closed orbit (mm)",
    "Turn per lap",
)


def predictions_markdown() -> str:
    """The predictions as tables, for the protocol's reader."""
    guides = {name: describe(guide) for name, guide in GUIDES.items()}
    lines = [
        "<!-- Written by `experiments/guide_study.py predict`. Theory only: no electron was tracked to make it. -->",
        "",
        "# Experiment 1: the predictions",
        "",
        "## The guides",
        "",
        "| " + " | ".join(GUIDE_COLUMNS) + " |",
        "|" + " --- |" * len(GUIDE_COLUMNS),
    ]
    for name, g in guides.items():
        closes = abs(g["turn_deg"]) > 1e-6
        orbit = [f"{g[key]:.2f}" if closes else "none" for key in ("closed_orbit_per_rho", "matched_reach_at_1mm_mm")]
        cells = [
            name,
            g["shape"],
            "yes" if g["knotted"] else "no",
            g["winding"],
            f"{g['length_mm']:.0f}",
            f"{g['bend_mm']:.1f}",
            f"{g['writhe']:+.3f}",
            f"{g['turn_deg']:+.1f}°",
            f"{g['drift_per_rho']:.2f}",
            f"{g['wander_mm']:.2f}",
            f"{g['net_turns_along_curve']:+.2f}",
            *orbit,
            f"{g['axis_reach_at_1mm_mm']:.2f}",
        ]
        lines.append("| " + " | ".join(cells) + " |")
    lines += [
        "",
        f"The wall is {WALL * 1e3:.1f} mm from the curve. ρ is the gyroradius at the guide field; a minus sign",
    ]
    lines += ["means the electron is sent against the field.", ""]
    for group, title in (
        ("main", "The main shots (20 mT, 20 laps)"),
        ("scaling", "The same gyroradii at 40 mT"),
        ("systematics", "The instrument made finer, one part at a time"),
        ("first-run", "The first run's cases, on its single-wire winding"),
    ):
        lines += [f"## {title}", ""]
        said = []
        for shot in shots():
            if shot.group != group:
                continue
            p = predict(GUIDES[shot.guide], shot)
            verdict = "not launched" if not p.launched else "held" if p.held else "lost"
            lap = "" if p.exit_lap is None else f"{p.exit_lap:.2f}"
            turn = "" if p.turn is None else f"{math.degrees(p.turn):+.1f}°"
            said.append(
                [shot.ident, f"{shot.w_ev:.1f}", verdict, _mm(p.reach), lap, _polar(p.centre), _mm(p.radius), turn]
            )
        lines += [*lab.table(SHOT_COLUMNS, said), ""]
    return "\n".join(lines)


# ---- the verdicts: each function is one criterion of PROTOCOL.md, section 5 ---------------------------------------
Row = tuple[Shot, Prediction, Measurement | None]


def _angle(a: float, b: float) -> float:
    """The smallest angle between two directions (rad)."""
    return abs(cmath.phase(cmath.exp(1j * (a - b))))


def _well_inside(p: Prediction) -> bool:
    return p.launched and p.reach is not None and p.reach <= INSIDE * WALL


def check_reach(rows: Iterable[Row]) -> Iterator[Check]:
    """C1: a shot the theory puts well inside the wall stays inside, and reaches as far as predicted."""
    for shot, p, m in rows:
        if shot.group == "first-run" or not _well_inside(p) or p.reach is None:
            continue
        if m is None or not m.held:
            yield Check("C1 reach", shot.ident, False, f"predicted {_mm(p.reach)} mm, and it was lost")
            continue
        ok = abs(m.reach - p.reach) <= REACH_TOL * p.reach + FLOOR
        yield Check("C1 reach", shot.ident, ok, f"predicted {_mm(p.reach)} mm, tracked {_mm(m.reach)} mm")


def check_held(rows: Iterable[Row]) -> Iterator[Check]:
    """C2: held or lost as predicted, where the predicted reach is not within BAND of the wall."""
    for shot, p, m in rows:
        if shot.group == "first-run" or not p.launched or p.reach is None or m is None:
            continue
        if abs(p.reach / WALL - 1) < BAND:
            continue
        yield Check(
            "C2 held or lost",
            shot.ident,
            m.held == p.held,
            f"predicted {'held' if p.held else 'lost'} (reach {_mm(p.reach)} mm), "
            f"tracked {'held' if m.held else 'lost'} after {m.laps:.2f} laps",
        )


def check_circle(rows: Iterable[Row]) -> Iterator[Check]:
    """C3 and C4: the crossings turn about the predicted closed orbit, by the predicted angle per lap."""
    for shot, p, m in rows:
        fit = m.fit() if m is not None and m.held else None
        if shot.group == "first-run" or fit is None or p.centre is None or p.turn is None or p.radius is None:
            continue
        if p.radius < FIT_RADIUS:
            continue
        off = abs(fit.centre - p.centre)
        yield Check(
            "C3 closed orbit",
            shot.ident,
            off <= CENTRE_TOL * abs(p.centre) + FLOOR,
            f"predicted {_polar(p.centre)} mm, fitted {_polar(fit.centre)} mm ({_mm(off)} mm apart)",
        )
        yield Check(
            "C4 turn per lap",
            shot.ident,
            _angle(fit.turn, p.turn) <= TURN_TOL,
            f"predicted {math.degrees(p.turn):+.1f}°, fitted {math.degrees(fit.turn):+.1f}°",
        )


def check_matched(rows: Iterable[Row]) -> Iterator[Check]:
    """C5: launched on the predicted closed orbit, an electron stays near it."""
    for shot, p, m in rows:
        if shot.launch != "matched" or not _well_inside(p) or p.centre is None or m is None or not m.held:
            continue
        strayed = max(abs(z - p.centre) for z in m.section)
        yield Check(
            "C5 matched launch",
            shot.ident,
            strayed <= MATCHED_TOL * abs(p.centre) + FLOOR,
            f"closed orbit {_mm(abs(p.centre))} mm from the axis, crossings stray {_mm(strayed)} mm from it",
        )


def check_ring(rows: Iterable[Row]) -> Iterator[Check]:
    """C6: the planar ring loses every electron, at the lap the drift predicts."""
    for shot, p, m in rows:
        if shot.guide != "ring" or m is None:
            continue
        if m.held or p.exit_lap is None:
            yield Check("C6 ring", shot.ident, False, f"held for {m.laps:.1f} laps; predicted lost")
            continue
        ok = abs(m.laps - p.exit_lap) <= max(EXIT_TOL * p.exit_lap, EXIT_FLOOR)
        yield Check("C6 ring", shot.ident, ok, f"predicted lost at lap {p.exit_lap:.2f}, tracked {m.laps:.2f}")


def check_repeats(rows: Sequence[Row]) -> Iterator[Check]:
    """C8 and C10: the same gyroradius at twice the field, or with a finer instrument, gives the same orbit."""
    by_shot = {(shot.guide, shot.rho, shot.launch): m for shot, _, m in rows if shot.group == "main"}
    for shot, _, m in rows:
        if shot.group not in ("scaling", "systematics") or shot.laps != LAPS:
            continue
        name = "C8 scaling" if shot.group == "scaling" else "C10 systematics"
        base = by_shot.get(("trefoil" if shot.group == "systematics" else shot.guide, shot.rho, shot.launch))
        if base is None or m is None:
            continue
        ok = m.held == base.held and abs(m.reach - base.reach) <= REPEAT_TOL * base.reach + REPEAT_FLOOR
        said = f"reach {_mm(m.reach, 3)} mm against {_mm(base.reach, 3)} mm"
        fits = m.fit(), base.fit()
        if fits[0] is not None and fits[1] is not None:
            ok = ok and _angle(fits[0].turn, fits[1].turn) <= REPEAT_TURN
            said += f"; turn {math.degrees(fits[0].turn):+.2f}° against {math.degrees(fits[1].turn):+.2f}°"
        yield Check(name, shot.ident, ok, said)


def check_first_run(rows: Iterable[Row]) -> Iterator[Check]:
    """C11, a retrodiction: the theory and this tracker against the first run's table, which was known."""
    table = json.loads((OUT / "first_run.json").read_text())
    for shot, p, m in rows:
        if shot.group != "first-run" or m is None:
            continue
        case = table["cases"][shot.guide][f"{shot.field * 1e3:g} mT"]
        i = min(range(4), key=lambda k: abs(table["energies_eV"][k] - shot.w_ev))
        yield Check(
            "C11 first run, this tracker",
            shot.ident,
            m.held == case["held"][i],
            f"first run {'held' if case['held'][i] else 'lost'} (reach {case['reach_mm'][i]:.1f} mm), "
            f"now {'held' if m.held else 'lost'} (reach {_mm(m.reach, 1)} mm)",
        )
        if p.reach is not None and abs(p.reach / WALL - 1) >= BAND:
            yield Check(
                "C11 first run, the theory",
                shot.ident,
                p.held == case["held"][i],
                f"first run {'held' if case['held'][i] else 'lost'}, theory {'held' if p.held else 'lost'} "
                f"(reach {_mm(p.reach, 1)} mm)",
            )


CHECKS: tuple[Callable[[Sequence[Row]], Iterator[Check]], ...] = (
    check_reach,
    check_held,
    check_circle,
    check_matched,
    check_ring,
    check_repeats,
    check_first_run,
)
SUBSETS: dict[str, Callable[[Shot], bool]] = {
    "C7 a twist without a knot (C1 to C5 on the helix and the figure-eight)": lambda s: (
        s.guide in ("helix", "figure-eight")
    ),
    "C9 against the field (C1 to C5 on the shots sent backwards)": lambda s: s.rho < 0 and s.group == "main",
}


def _load(record: Mapping[str, Any]) -> Measurement:
    section = tuple(complex(re, im) for re, im in record["section"])
    return Measurement(record["held"], record["laps"], record["reach"], section)


class GuideExperiment(lab.Experiment[Shot, Prediction, Measurement]):
    """The experiment as array16.lab runs it. It needs no device: its bench is compute alone."""

    out = OUT
    registered = (
        OUT / "PROTOCOL.md",
        OUT / "predictions.md",
        OUT / "first_run.json",
        Path(__file__).resolve(),
        lab.ROOT / "array16" / "particles",
        lab.ROOT / "array16" / "lab.py",
    )

    def shots(self) -> Sequence[Shot]:
        return shots()

    def ident(self, shot: Shot) -> str:
        return shot.ident

    def predict(self, shot: Shot) -> Prediction:
        return predict(GUIDES[shot.guide], shot)

    def batches(self, todo: Sequence[Shot]) -> Sequence[Sequence[Shot]]:
        """Shots of one speed on one guide go round together, slowest first: they take the longest."""
        ordered = sorted(todo, key=lambda shot: (abs(shot.rho) * shot.field / shot.laps, shot.speed_key))
        return [list(group) for _, group in groupby(ordered, key=lambda shot: shot.speed_key)]

    def measure(self, batch: Sequence[Shot], devices: Mapping[str, Any]) -> Mapping[str, Measurement]:
        return track(GUIDES[batch[0].guide], batch)

    def load(self, record: Mapping[str, Any]) -> Measurement:
        return _load(record)

    def judge(self, rows: Sequence[Row]) -> Sequence[Check]:
        return [check for criterion in CHECKS for check in criterion(rows)]

    def apparatus(self) -> Mapping[str, Any]:
        return {"wall_mm": WALL * 1e3, "guides": {name: describe(guide) for name, guide in GUIDES.items()}}

    def summary(self, checks: Sequence[Check]) -> Mapping[str, Sequence[Check]]:
        by_ident = {shot.ident: shot for shot in shots()}
        first_five = [check for check in checks if check.criterion[:2] in ("C1", "C2", "C3", "C4", "C5")]
        return {title: [c for c in first_five if member(by_ident[c.ident])] for title, member in SUBSETS.items()}


def main(argv: Sequence[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    code = lab.main(GuideExperiment(), args, about=(__doc__ or "").split("\n", 1)[0])
    if code == 0 and args[:1] == ["predict"]:
        (OUT / "predictions.md").write_text(predictions_markdown())
    return code


if __name__ == "__main__":
    sys.exit(main())
