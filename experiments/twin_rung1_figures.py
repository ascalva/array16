"""Experiment 2, stage 1, after the run: the figures of its report (experiments/twin_rung1/REPORT.md), drawn from the
committed records alone.

It reads experiments/twin_rung1/runs.jsonl (the tracking: for each shot, energy and sense, how many electrons were
held, lost to the guide, stopped by the gas or short of their laps, and the fewest laps, or RF periods, after which one
was lost), predictions.json (what was registered: each prediction, its band and the numbers it rests on), results.json
(the verdicts) and the design files (the paths and their frames). It tracks nothing and solves nothing, and it reads a
scan with twin_rung1's own functions (load, Scan.threshold, mean_limit), so that the figures show what the verdicts
read.

  holds.png      each hold: the fraction held against energy in both senses, the static limit and the registered band
  wells.png      each well: the fraction held against energy over the depth, the exact motion's first loss and the
                 band; and the slot guide's secular exponent against Mathieu's
  fractions.png  the fraction of its static limit each path kept, planar paths against those out of the plane, and the
                 price of the knot against its band
  must_fail.png  the cases that must fail: when each energy and sense's first electron was lost, in RF periods,
                 beside the control
  paths.png      R1, the flat ring, T915-S and the helix in three dimensions, the boards drawn as a ribbon shaded by the
                 curvature in their plane, and where along each path the first electron of each lost energy was lost
  along.png      the same places along the path, unrolled, for every hold, under the curvature in and toward the boards

What the records do not hold: any position. A lost electron's place along the path follows from the laps it made
before it was lost (the detector's watch counts them; the record keeps the fewest of each energy's lost electrons),
taken in its sense of travel from the start, modulo one lap. Where in the section it left (which copper, or the
aperture), the places and times of the energy's other lost electrons, and every trajectory are not in these records;
stage 2's will carry each lost electron's place, section position and step (the watch keeps the place and the step).

    uv run python experiments/twin_rung1_figures.py              # the six figures, into experiments/twin_rung1/
"""

from __future__ import annotations

import argparse
import importlib
import json
import math
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import twin_rung1 as rung

from array16.guide.sections import FacingStrips, Rails
from array16.lab import ROOT
from array16.particles import Arr, ClosedCurve
from array16.twin import Design, curvature_vector, framing

SRC = rung.OUT  # the committed records: experiments/twin_rung1/
FIGURES = ("holds.png", "wells.png", "fractions.png", "must_fail.png", "paths.png", "along.png")
MM = 1e-3

# the house style (experiments/knotted_guide_study.py): ink, recessive grid, the categorical colours in fixed order
INK, INK2, MUTED, GRID, EDGE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
SURFACE, PALE = "#fcfcfb", "#ecebe6"
SENSE_COLOUR = {1: SERIES[0], -1: SERIES[1]}  # with the drive's wave, against it
SENSE_WORD = {1: "with the wave", -1: "against it"}


# ---- the records ---------------------------------------------------------------------------------------------------
@dataclass(frozen=True, eq=False)
class Row:
    """One shot as the records hold it: its name, what it was (runs.jsonl's shot), what was measured (read back by
    twin_rung1.load), the prediction registered for it (predictions.json), how long it ran (s) and its provenance."""

    ident: str
    shot: Mapping[str, Any]
    measured: rung.Measurement
    prediction: Mapping[str, Any]
    seconds: float
    provenance: Mapping[str, Any]

    @property
    def kind(self) -> str:
        return str(self.shot["kind"])

    @property
    def scan(self) -> rung.Scan:
        if not isinstance(self.measured, rung.Scan):
            raise TypeError(f"{self.ident} is a {self.kind}, not a scan")
        return self.measured

    @property
    def beat(self) -> rung.Beat:
        if not isinstance(self.measured, rung.Beat):
            raise TypeError(f"{self.ident} is a {self.kind}, not the secular frequency")
        return self.measured

    @property
    def basis(self) -> Mapping[str, float]:
        out: Mapping[str, float] = self.prediction["basis"]
        return out


@dataclass(frozen=True, eq=False)
class Records:
    """The shots on file in the protocol's order, the registered apparatus, and the verdicts."""

    rows: tuple[Row, ...]
    apparatus: Mapping[str, Any]
    checks: tuple[Mapping[str, Any], ...]

    def row(self, ident: str) -> Row:
        for row in self.rows:
            if row.ident == ident:
                return row
        raise KeyError(f"no record of {ident!r}")

    def of(self, *kinds: str) -> tuple[Row, ...]:
        return tuple(row for row in self.rows if row.kind in kinds)

    def verdict(self, criterion: str, ident: str) -> bool:
        """Whether the check `criterion` passed on the shot `ident`, as results.json has it."""
        for check in self.checks:
            if check["criterion"] == criterion and check["ident"] == ident:
                return bool(check["passed"])
        raise KeyError(f"no verdict of {criterion!r} on {ident!r}")


def read(src: Path = SRC) -> Records:
    """The committed records in `src`: runs.jsonl (a later record of a shot replaces an earlier one, as array16.lab
    reads it), predictions.json and results.json."""
    runs = [json.loads(line) for line in (src / "runs.jsonl").read_text().splitlines() if line.strip()]
    latest = {record["ident"]: record for record in runs}
    registered = json.loads((src / "predictions.json").read_text())
    rows = tuple(
        Row(
            ident=entry["ident"],
            shot=latest[entry["ident"]]["shot"],
            measured=rung.load(latest[entry["ident"]]["measurement"]),
            prediction=entry["prediction"],
            seconds=float(latest[entry["ident"]]["seconds"]),
            provenance=latest[entry["ident"]]["provenance"],
        )
        for entry in registered["shots"]
        if entry["ident"] in latest
    )
    checks = tuple(json.loads((src / "results.json").read_text())["checks"])
    return Records(rows, registered["apparatus"], checks)


@cache
def records() -> Records:
    return read()


# ---- what a scan reads ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Sense:
    """One sense of a hold (+1 with the drive's wave, -1 against it): the protocol's robust threshold, the first loss
    (the lowest energy at which an electron was lost), the highest energy held below it, and the lowest energy from
    which every electron judged, there and at every energy above, was lost (eV; inf if there is none)."""

    sense: int
    threshold: float
    first_loss: float
    held_below: float
    all_lost: float


@dataclass(frozen=True)
class Hold:
    """A hold as the records read: the static limit registered and its band (eV), the limit the protocol reads (the
    senses' robust thresholds' geometric mean, twin_rung1.mean_limit), each sense, whether the path is planar (a
    circle) and the verdict."""

    ident: str
    planar: bool
    static: float
    low: float
    high: float
    limit: float
    senses: tuple[Sense, ...]
    passed: bool

    @property
    def fraction(self) -> float:
        """The limit over the static limit."""
        return self.limit / self.static


def all_lost_from(energies: Sequence[float], held: Sequence[int], lost: Sequence[int]) -> float:
    """The lowest energy from which every electron judged was lost, there and at every energy above (inf if none):
    energies where the gas stopped every electron say nothing."""
    out = math.inf
    for k in range(len(energies) - 1, -1, -1):
        if held[k]:
            break
        if lost[k]:
            out = energies[k]
    return out


def hold(rec: Records, row: Row) -> Hold:
    scan, p = row.scan, row.prediction
    senses = []
    for k, sense in enumerate(scan.senses):
        below, first = scan.edges(k)
        senses.append(
            Sense(sense, scan.threshold(k), first, below, all_lost_from(scan.energies, scan.held[k], scan.lost[k]))
        )
    criterion = "E2-22 static limit" if row.shot["test"] == "E2-22" else "E2-21 hold"
    if row.shot["test"] == "E2-21" and row.shot["q"] is not None:
        criterion = "E2-21 hold control"
    planar = Design.load(ROOT / str(row.shot["design"])).path.kind == "circle"
    return Hold(
        ident=row.ident,
        planar=planar,
        static=float(p["value"]),
        low=float(p["low"]),
        high=float(p["high"]),
        limit=rung.mean_limit(scan),
        senses=tuple(senses),
        passed=rec.verdict(criterion, row.ident),
    )


@dataclass(frozen=True)
class Loss:
    """An energy and sense at which electrons were lost: how many of how many judged, and `when` the first of them
    was (the fewest laps for a hold, RF periods otherwise)."""

    sense: int
    energy: float
    when: float
    lost: int
    judged: int

    @property
    def place(self) -> float:
        """Where along the path, as a fraction of a lap from the start, a hold's first-lost electron was lost: its
        laps taken in its sense of travel, modulo one."""
        return (self.sense * self.when) % 1.0


def losses(scan: rung.Scan) -> tuple[Loss, ...]:
    """Every energy and sense of a scan at which an electron was lost."""
    out = []
    for k, sense in enumerate(scan.senses):
        for j, energy in enumerate(scan.energies):
            when = scan.first_loss[k][j]
            if math.isfinite(when):
                judged = scan.held[k][j] + scan.lost[k][j]
                out.append(Loss(sense, energy, when, scan.lost[k][j], judged))
    return tuple(out)


def held_fraction(scan: rung.Scan, k: int) -> Arr:
    """The fraction of the electrons judged (held or lost) that were held, at each energy of a row (nan where the gas
    stopped every one)."""
    held = np.asarray(scan.held[k], dtype=np.float64)
    judged = held + np.asarray(scan.lost[k], dtype=np.float64)
    out: Arr = np.divide(held, judged, out=np.full_like(held, np.nan), where=judged > 0)
    return out


def ratio(rec: Records, row: Row) -> float:
    """A hold's limit over its reference's (T915-S's, for E2-22), as the price of the knot reads it."""
    base = rec.row(next(r.ident for r in rec.of("hold") if r.shot["design"] == rung.REFERENCE[row.shot["test"]]))
    return rung.mean_limit(row.scan) / rung.mean_limit(base.scan)


def title(row: Row) -> str:
    """A shot's name for a panel: the design, as drawn or on the bench, and its drive with the null's q there. A scan
    records the mean volts round the ring; the secular frequency is read at the start station's volts, and says so."""
    where = "as drawn" if row.shot["realisation"] == 0 else "on the bench"
    m = row.measured
    if isinstance(m, rung.Beat):
        return f"{row.shot['label']} {where}, {m.volts:.2f} V at the start station (q {m.q:#.3g})"
    if isinstance(m, rung.Scan):
        return f"{row.shot['label']} {where}, {m.volts:.0f} V (q {m.q:#.3g})"
    return f"{row.shot['label']} {where}"


# ---- the paths -----------------------------------------------------------------------------------------------------
@dataclass(frozen=True, eq=False)
class Route:
    """A design's path as the twin builds it (TwinGuide.curve and .framing), at the registered length: the points (m),
    the direction across the path in the boards' plane and the boards' normal, the curvature in the boards' plane
    and toward them (1/m), the copper's width across the path (m) and the length (m)."""

    points: Arr
    tangents: Arr
    across: Arr
    normal: Arr
    kappa_g: Arr
    kappa_n: Arr
    width: float
    length: float

    @property
    def s(self) -> Arr:
        """Each sample's place along the path, as a fraction of a lap from the start."""
        out: Arr = np.arange(len(self.points)) / len(self.points)
        return out

    def at(self, place: float) -> Arr:
        """The point (m) at a fraction of a lap along the path."""
        out: Arr = self.points[round(place * len(self.points)) % len(self.points)]
        return out


def route(design: str, length: float) -> Route:
    """The path of a design file at `length` m, framed as the twin frames it: a circle in the xy plane, or a torus
    curve scaled to the length; the twist-free axes turned by the frame's theta, or the framing that closes with k
    half turns (array16.twin.framing)."""
    d = Design.load(ROOT / design)
    spec = d.path
    if spec.kind == "circle":
        curve = ClosedCurve.circle(length / (2 * math.pi), spec.samples)
    elif spec.kind == "torus":
        knot = ClosedCurve.torus_knot(1.0, spec.ratio, spec.samples, spec.p, spec.q)
        curve = ClosedCurve(knot.points * (length / knot.length), length)
    else:
        raise ValueError(f"{design}: a {spec.kind} path has no lap to draw")
    tangents = curve.tangents
    if d.frame.k is not None:
        framed = framing(curve, d.frame.k, d.frame.theta0)
        across, normal = framed.across, framed.normal
    else:
        e1, e2 = curve.axes
        across = math.cos(d.frame.theta) * e1 + math.sin(d.frame.theta) * e2
        normal = np.cross(tangents, across)
    bend = curvature_vector(curve)
    geometry = d.section.geometry
    if isinstance(geometry, FacingStrips):
        width = geometry.width
    elif isinstance(geometry, Rails):
        width = geometry.gap + 2 * geometry.rail
    else:
        raise ValueError(f"{design}: no width drawn for a {type(geometry).__name__}")
    return Route(
        points=curve.points,
        tangents=tangents,
        across=across,
        normal=normal,
        kappa_g=np.abs((bend * across).sum(1)),
        kappa_n=np.abs((bend * normal).sum(1)),
        width=width,
        length=curve.length,
    )


def route_of(rec: Records, row: Row) -> Route:
    """A shot's path at the length its prediction rests on (the ring tuned m guided waves round)."""
    return route(str(row.shot["design"]), float(row.basis["length_mm"]) * MM)


def lobes(path: Route) -> tuple[tuple[float, float], ...]:
    """The lobes of a path whose bend turns through the boards' plane: each stretch between the two maxima of the
    curvature in the boards' plane that flank one of its tightest bends (fractions of a lap from the start, in order);
    none where that curvature is constant, as on a circle."""
    kg = path.kappa_g
    if float(np.ptp(kg)) <= 1e-6 * float(kg.max()):
        return ()
    n = len(kg)
    peaks = [i / n for i in range(n) if kg[i] >= kg[i - 1] and kg[i] >= kg[(i + 1) % n] and kg[i] > 0.8 * kg.max()]
    return tuple(zip(peaks[0::2], peaks[1::2], strict=True))


def past(loss: Loss, stretches: Sequence[tuple[float, float]]) -> float:
    """How far a hold's loss lies past the nearest lobe in its direction of travel (fraction of a lap; 0 within one,
    inf without lobes)."""
    out = math.inf
    for a, b in stretches:
        if a <= loss.place <= b:
            return 0.0
        out = min(out, (loss.place - b) % 1.0 if loss.sense == 1 else (a - loss.place) % 1.0)
    return out


# ---- drawing -------------------------------------------------------------------------------------------------------
def _pyplot() -> Any:
    """matplotlib's pyplot on the Agg backend, in the house style."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from cycler import cycler

    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 9.5,
            "axes.labelsize": 9.5,
            "figure.titlesize": 11,
            "legend.fontsize": 8.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "axes.edgecolor": EDGE,
            "axes.labelcolor": INK,
            "axes.titlecolor": INK,
            "text.color": INK,
            "xtick.color": INK2,
            "ytick.color": INK2,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.axisbelow": True,
            "axes.titlelocation": "left",
            "lines.linewidth": 1.6,
            "lines.markersize": 5,
            "legend.frameon": False,
            "axes.prop_cycle": cycler(color=list(SERIES)),
            "savefig.dpi": 140,
            "mathtext.default": "regular",
        }
    )
    return plt


def _save(plt: Any, fig: Any, out: Path, name: str) -> Path:
    path = out / name
    fig.savefig(path)
    plt.close(fig)
    return path


def _sci(value: float) -> str:
    """A small number as the reports write it: 1.4 × 10⁻⁴."""
    mantissa, exponent = f"{value:.1e}".split("e")
    sup = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    return f"{mantissa} × 10{str(int(exponent)).translate(sup)}"


def _verdict(passed: bool, value: float, low: float, high: float) -> str:
    if passed:
        return "pass"
    return "FAIL, above the band" if value > high else "FAIL, below the band"


def figure_holds(rec: Records, out: Path) -> Path:
    """Each hold: the fraction held for 20 laps at each energy in both senses, the band, the static limit, each
    sense's robust threshold and the limit the protocol reads."""
    plt = _pyplot()
    rows = rec.of("hold")
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.8))
    for ax, row in zip(axes.flat, rows, strict=True):
        h, scan = hold(rec, row), row.scan
        energies = np.asarray(scan.energies)
        ax.axvspan(h.low, h.high, color=PALE, lw=0, zorder=0)
        ax.axvline(h.static, color=INK, lw=1.3, zorder=2)
        for k, sense in enumerate(scan.senses):
            colour = SENSE_COLOUR[sense]
            wide = sense == 1  # drawn wider beneath the other, so that where the two senses agree both show
            ax.plot(
                energies,
                held_fraction(scan, k),
                color=colour,
                lw=3.2 if wide else 1.1,
                marker="o",
                ms=5.0 if wide else 2.4,
                zorder=3 if wide else 4,
            )
            reading = h.senses[k].threshold
            if math.isfinite(reading) and reading > 0:
                ax.plot([reading], [0.5], marker="v", ms=9, color=colour, mec=SURFACE, mew=1.0, zorder=5)
        ax.axvline(h.limit, color=INK, lw=1.3, ls=(0, (4, 2.5)), zorder=4)
        ax.text(h.static, 1.045, f" static {h.static:.4g}", fontsize=8, color=INK, ha="left", va="bottom")
        ax.text(h.limit, -0.13, f"held to {h.limit:.4g} ", fontsize=8, color=INK, ha="right", va="top")
        if row.shot["realisation"] is None:
            low, high = row.basis["earth_mean_low_eV"], row.basis["earth_mean_high_eV"]
            ax.plot([low, high], [1.16, 1.16], color=INK2, lw=2.2, solid_capstyle="butt")
            ax.text(high, 1.16, "  Earth's range for the senses' mean", fontsize=7.5, color=INK2, va="center")
        verdict = _verdict(h.passed, h.limit, h.low, h.high)
        ax.set_title(
            f"{title(row)}\nheld to {h.limit:.4g} eV, {h.fraction:.3f} of the static limit\n"
            f"band {h.low:.4g} to {h.high:.4g} eV: {verdict}",
            fontsize=9,
        )
        ax.set_ylim(-0.2, 1.24)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xlim(energies[0] / 1.03, energies[-1] * 1.03)
        ax.set_xlabel("energy along the path (eV)")
        ax.set_ylabel("fraction held for 20 laps")
    handles = [
        plt.Line2D([], [], color=SERIES[0], marker="o", ms=5, lw=3.2, label="held, with the drive's wave"),
        plt.Line2D([], [], color=SERIES[1], marker="o", ms=2.4, lw=1.1, label="held, against it"),
        plt.Line2D([], [], color=MUTED, marker="v", ms=8, lw=0, label="each sense's robust threshold"),
        plt.Line2D([], [], color=INK, lw=1.3, ls=(0, (4, 2.5)), label="the limit read: their geometric mean"),
        plt.Line2D([], [], color=INK, lw=1.3, label="static limit (registered)"),
        plt.Rectangle((0, 0), 1, 1, color=PALE, label="registered band"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.01, 0.94), ncol=6, fontsize=8.5)
    every = [hold(rec, row) for row in rows]
    drawn = [h for h, row in zip(every, rows, strict=True) if row.shot["realisation"] == 0]
    flat = ", ".join(f"{h.fraction:.3f}" for h in drawn if h.planar)
    knotted = ", ".join(f"{h.fraction:.3f}" for h in drawn if not h.planar)
    missed = [h for h in every if not h.passed]
    above = sum(h.limit > h.high for h in missed)
    fig.suptitle(
        "Each hold, electrons launched on the bend's static equilibrium at every energy and RF phase in both senses: "
        f"the planar rings kept their static limits ({flat}),\nthe paths out of the plane less ({knotted}); "
        f"{len(missed)} of {len(every)} holds missed their registered bands, "
        + ("each above it" if above == len(missed) else f"{above} of them above"),
        x=0.01,
        y=0.995,
        ha="left",
        fontsize=10.5,
    )
    fig.subplots_adjust(left=0.05, right=0.985, top=0.81, bottom=0.07, hspace=0.62, wspace=0.2)
    return _save(plt, fig, out, "holds.png")


def figure_wells(rec: Records, out: Path) -> Path:
    """Each well: the fraction held after 40 RF periods at each energy, over the depth at the start station, with
    the exact motion's first loss and the band; and the slot guide's secular exponent against Mathieu's."""
    plt = _pyplot()
    fig = plt.figure(figsize=(15.5, 8.6))
    grid = fig.add_gridspec(2, 6, hspace=0.62, wspace=0.9)
    wells = rec.of("well")
    spots = (grid[0, 0:2], grid[0, 2:4], grid[0, 4:6], grid[1, 0:3])
    for spot, row in zip(spots, wells, strict=True):
        ax = fig.add_subplot(spot)
        scan, p, basis = row.scan, row.prediction, row.basis
        depth = float(basis["depth_start_eV"])
        x = np.asarray(scan.energies) / depth
        ax.axvspan(p["low"] / depth, p["high"] / depth, color=PALE, lw=0, zorder=0)
        ax.axvline(1.0, color=INK, lw=1.3, zorder=2)
        ax.text(1.0, 1.06, " D", fontsize=8.5, va="bottom")
        ax.plot(x, held_fraction(scan, 0), color=SERIES[2], lw=1.1, marker="o", ms=3, zorder=3)
        slot = row.shot["test"] == "slot"
        if slot:
            model, robust = float(basis["model_edge_D"]), math.nan
            found = scan.edge_threshold(0)
            word = "escape, the edges' threshold"
            criterion = "slot depth"
        else:
            model, robust = float(basis["model_first_eV"]) / depth, float(basis["model_robust_eV"]) / depth
            found = scan.lowest_lost(0)
            word = "first lost"
            criterion = "E2-21 well" if row.shot["q"] is None else "E2-21 well control"
        ax.axvline(model, color=INK, lw=1.3, ls=(0, (4, 2.5)), zorder=4)
        if math.isfinite(robust):
            ax.axvline(robust, color=INK2, lw=1.0, ls=":", zorder=4)
        ax.plot([found / depth], [-0.09], marker="^", ms=9, color=INK, zorder=6, clip_on=False)
        verdict = _verdict(rec.verdict(criterion, row.ident), found, float(p["low"]), float(p["high"]))
        model_word = f"{model:.3f} D" if slot else f"{float(basis['model_first_eV']):.4g} eV"
        as_drawn = "" if row.shot["realisation"] == 0 else " as drawn"
        ax.set_title(
            f"{title(row)}, D = {depth:.4g} eV{as_drawn}\n{word} at {found:.4g} eV, {found / depth:.3f} D; "
            f"exact motion {model_word}\nband {float(p['low']):.4g} to {float(p['high']):.4g} eV: {verdict}",
            fontsize=8.8,
        )
        ax.set_ylim(-0.16, 1.16)
        ax.set_yticks([0, 0.5, 1.0])
        ax.set_xlabel("energy across the path over the depth at the start, W / D")
        ax.set_ylabel("fraction held, 40 RF periods")
    handles = [
        plt.Line2D([], [], color=SERIES[2], marker="o", ms=4, label="held (two launches an energy)"),
        plt.Line2D([], [], color=INK, marker="^", ms=8, lw=0, label="tracked reading"),
        plt.Line2D([], [], color=INK, lw=1.3, ls=(0, (4, 2.5)), label="the exact motion's first loss (registered)"),
        plt.Line2D([], [], color=INK2, lw=1.0, ls=":", label="the exact motion's robust threshold"),
        plt.Line2D([], [], color=INK, lw=1.3, label="pseudopotential's depth D at the start"),
        plt.Rectangle((0, 0), 1, 1, color=PALE, label="registered band"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.01, 0.935), ncol=6, fontsize=8.5)

    ax = fig.add_subplot(grid[1, 3:6])
    row = rec.of("beta")[0]
    beat, mathieu = row.beat, float(row.prediction["value"])
    names = []
    for k, (u, v) in enumerate(beat.offsets):
        names.append(f"diagonal u = {'' if u * v > 0 else '−'}v" if u and v else f"along {'u' if u else 'v'}")
        ax.plot([(beat.beta[k] / mathieu - 1) * 1e4], [k], "o", color=SERIES[2], ms=8, zorder=4)
    tol = rung.BETA_TOL * 1e4
    ax.axvspan(-tol, tol, color=PALE, lw=0, zorder=0)
    ax.axvline(0.0, color=INK, lw=1.3)
    first = (float(row.basis["first_order_beta"]) / mathieu - 1) * 1e4
    ax.annotate(
        f"first order q/√2: {first:.0f} × 10⁻⁴".replace("-", "−"),
        xy=(-tol * 1.2, 1.5),
        xytext=(-tol * 0.95, 1.5),
        fontsize=8,
        va="center",
        arrowprops={"arrowstyle": "->", "color": INK2, "lw": 0.9},
    )
    ax.set_xlim(-tol * 1.2, tol * 1.2)
    ax.set_ylim(-0.7, len(names) - 0.3)
    ax.set_yticks(range(len(names)), names)
    ax.set_ylabel(f"start, {rung.BETA_OFFSET:g}h off the null")
    ax.set_xlabel("tracked β over Mathieu's, less one (× 10⁻⁴)")
    worst = max(abs(b / mathieu - 1) for b in beat.beta)
    verdict = "pass" if rec.verdict("slot secular frequency", row.ident) else "FAIL"
    ax.set_title(
        f"{title(row)}\nthe secular exponent over {row.shot['periods']} RF periods, 4 of 4 held: "
        f"β {min(beat.beta):.6f} to {max(beat.beta):.6f}\n"
        f"against Mathieu's {mathieu:.6f} at q {beat.q:.4f}, at most {_sci(worst)} off; "
        f"band ±{_sci(rung.BETA_TOL)} of it: {verdict}",
        fontsize=8.8,
    )
    drawn = [r for r in wells if r.shot["realisation"] == 0 and r.shot["test"] == "E2-21"]
    firsts = " and ".join(
        f"{r.scan.lowest_lost(0) / float(r.basis['depth_start_eV']):.3f} D at q {r.scan.q:.2f}" for r in drawn
    )
    slot = next(r for r in wells if r.shot["test"] == "slot")
    escape = slot.scan.edge_threshold(0) / float(slot.basis["depth_start_eV"])
    fig.suptitle(
        "The wells: electrons sent from the null toward the escape saddle at the drive's zero crossings were first "
        f"lost on R1's rails at {firsts},\nat the exact motion's first losses; the straight slot guide escapes at "
        f"{escape:.3f} D (the exact motion {float(slot.basis['model_edge_D']):.3f} D) and keeps Mathieu's exponent "
        f"to {_sci(worst)}",
        x=0.01,
        y=0.995,
        ha="left",
        fontsize=10.5,
    )
    fig.subplots_adjust(left=0.06, right=0.985, top=0.82, bottom=0.08)
    return _save(plt, fig, out, "wells.png")


def figure_fractions(rec: Records, out: Path) -> Path:
    """The fraction of its static limit each path kept, against its band; and the price of the knot."""
    plt = _pyplot()
    fig = plt.figure(figsize=(15.5, 6.6))
    grid = fig.add_gridspec(1, 2, width_ratios=(1.75, 1.0), wspace=0.3)
    ax = fig.add_subplot(grid[0, 0])
    rows = rec.of("hold")

    def group(row: Row) -> tuple[int, str]:
        if row.shot["realisation"] is None:
            return 2, "laminate and Earth drawn"
        return (0, "planar") if hold(rec, row).planar else (1, "out of the plane")

    ordered = sorted(rows, key=lambda row: group(row)[0])  # stable: the protocol's order within a group
    labels = []
    for place, row in enumerate(ordered):
        h = hold(rec, row)
        y = float(len(ordered) - 1 - place)
        ax.plot([h.low / h.static, h.high / h.static], [y, y], color=PALE, lw=13, solid_capstyle="butt", zorder=1)
        for sense in h.senses:
            dy = 0.17 if sense.sense == 1 else -0.17
            colour = SENSE_COLOUR[sense.sense]
            top = sense.all_lost if math.isfinite(sense.all_lost) else sense.first_loss
            ax.plot([sense.first_loss / h.static, top / h.static], [y + dy, y + dy], color=colour, lw=1.6, zorder=3)
            ax.plot(
                [sense.threshold / h.static],
                [y + dy],
                marker="v" if sense.sense == 1 else "^",
                ms=8,
                color=colour,
                mec=SURFACE,
                zorder=4,
            )
        ax.plot([h.fraction], [y], "o", ms=9, color=INK, zorder=5)
        verdict = "pass" if h.passed else "FAIL, above" if h.limit > h.high else "FAIL, below"
        ax.text(1.115, y, f"{h.fraction:.3f}  {verdict}", va="center", fontsize=8.5)
        labels.append((y, f"{title(row)}\n{group(row)[1]}"))
    ax.axvline(1.0, color=INK, lw=1.3)
    ax.set_yticks([y for y, _ in labels], [text for _, text in labels], fontsize=8)
    ax.set_xlim(0.45, 1.25)
    ax.set_ylim(-0.6, len(labels) - 0.4)
    ax.set_xlabel("energy held over the path's static limit")
    ax.set_title(
        "What each path kept: the limit read (dot), each sense's robust threshold (triangles) and its edge, from the\n"
        "first loss to the energy from which every electron was lost (lines), against the registered band (grey)",
        fontsize=9,
    )
    handles = [
        plt.Line2D([], [], color=INK, marker="o", ms=7, lw=0, label="limit read (the senses' geometric mean)"),
        plt.Line2D([], [], color=SERIES[0], marker="v", ms=7, lw=1.6, label="with the wave: threshold, edge"),
        plt.Line2D([], [], color=SERIES[1], marker="^", ms=7, lw=1.6, label="against it: threshold, edge"),
        plt.Rectangle((0, 0), 1, 1, color=PALE, label="registered band"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=8, ncol=4, bbox_to_anchor=(-0.02, -0.12))

    bx = fig.add_subplot(grid[0, 1])
    priced = [row for row in rec.of("hold") if row.prediction["ratio"] is not None]
    for y, row in enumerate(priced):
        p = row.prediction
        found = ratio(rec, row)
        bx.plot([p["ratio_low"], p["ratio_high"]], [y, y], color=PALE, lw=13, solid_capstyle="butt", zorder=1)
        bx.plot([p["ratio"]], [y], marker="|", ms=18, mew=2.0, color=INK2, zorder=3)
        bx.plot([found], [y], "o", ms=9, color=INK, zorder=5)
        passed = rec.verdict("E2-22 price of the knot", row.ident)
        bx.text(
            found,
            y + 0.22,
            f"{found:.3f}, static {float(p['ratio']):.3f}: {'pass' if passed else 'FAIL'}",
            ha="center",
            fontsize=8.5,
        )
    bx.set_yticks(range(len(priced)), [f"{row.shot['label']} over T915-S" for row in priced])
    bx.set_ylim(-0.6, len(priced) - 0.3)
    bx.set_xlim(1.2, 3.6)
    bx.set_xlabel("held limit over T915-S's")
    bx.set_title(
        "The price of the knot: the tracked ratio (dot) against the static\nlimits' ratio (tick) and its "
        "registered band, 15 % either way (grey)",
        fontsize=9,
    )
    drawn = [(row, hold(rec, row)) for row in rows if row.shot["realisation"] == 0]
    carried = [(row, h) for row, h in drawn if h.planar and row.shot["q"] is None]  # banded on four rods' discount
    controls = [(row, h) for row, h in drawn if h.planar and row.shot["q"] is not None]  # banded on CONTROL_BAND
    knotted = [h.fraction for _, h in drawn if not h.planar]

    def placed(row: Row) -> str:
        if rec.verdict("E2-22 price of the knot", row.ident):
            return "inside its band"
        return "above its band" if ratio(rec, row) > float(row.prediction["ratio_high"]) else "below its band"

    said = ", and ".join(
        f"the {row.shot['label']}'s limit over T915-S's, {ratio(rec, row):.3f}, {placed(row)}" for row in priced
    )
    held = " and ".join(f"{row.shot['label']} at q {row.scan.q:.2f}" for row, _ in carried)
    control = "; ".join(
        f"{row.shot['label']}'s control at q {row.scan.q:.2f} held {h.fraction:.3f}, inside its {rung.CONTROL_BAND[0]} "
        f"to {rung.CONTROL_BAND[1]}"
        for row, h in controls
    )
    fig.suptitle(
        f"{held} held "
        + " and ".join(f"{h.fraction:.3f}" for _, h in carried)
        + f" of their static limits, above the four rods' {rung.DISCOUNT[0]} to {rung.DISCOUNT[1]} their bands "
        f"carried;\n{control}; the paths out of the plane "
        + " and ".join(f"{x:.3f}" for x in knotted)
        + f".\nSo the price of the knot: {said}",
        x=0.01,
        y=0.99,
        ha="left",
        fontsize=10.5,
    )
    fig.subplots_adjust(left=0.165, right=0.98, top=0.79, bottom=0.2)
    return _save(plt, fig, out, "fractions.png")


def figure_must_fail(rec: Records, out: Path) -> Path:
    """The must-fail shots and their control: when the first electron of each energy and sense was lost (RF periods),
    and the control's electrons held through the run."""
    plt = _pyplot()
    rows = rec.of("fail", "control")
    fig, axes = plt.subplots(1, len(rows), figsize=(15.5, 4.9), sharey=True)
    for ax, row in zip(axes, rows, strict=True):
        scan = row.scan
        periods = float(row.shot["periods"])
        energies = np.asarray(scan.energies)
        ax.axvline(periods, color=INK, lw=1.3)
        ax.text(periods, energies[-1] * 1.55, "end of the run ", fontsize=8, ha="right", va="center")
        for k, sense in enumerate(scan.senses):
            colour = SENSE_COLOUR[sense]
            dy = 1.08 if sense == 1 else 1 / 1.08
            for j, energy in enumerate(scan.energies):
                when = scan.first_loss[k][j]
                if math.isfinite(when):
                    ax.plot([when], [energy * dy], "o", ms=7, color=colour, zorder=4)
                    ax.plot([when, periods], [energy * dy, energy * dy], color=colour, lw=0.6, alpha=0.35, zorder=2)
                if scan.held[k][j]:
                    ax.plot([0.3, periods], [energy * dy, energy * dy], color=colour, lw=2.4, zorder=3)
                    ax.text(0.33, energy * dy, f"{scan.held[k][j]} of 4 held", fontsize=7, va="bottom", color=INK2)
        ax.set_xscale("log")
        ax.set_yscale("log", base=2)
        ax.set_xlim(0.3, 400)
        ax.set_xticks([0.5, 1, 2, 5, 10, 20, 50, 100, 200], ["0.5", "1", "2", "5", "10", "20", "50", "100", "200"])
        ax.minorticks_off()
        ax.set_yticks(energies, [f"{e:g} eV" for e in energies])
        ax.set_ylim(energies[0] / 1.6, energies[-1] * 1.9)
        ax.set_xlabel("RF periods after the launch")
        held, launched = scan.held_total(), scan.launched()
        if row.kind == "control":
            verdict = "pass" if rec.verdict("must fail control", row.ident) else "FAIL"
            said = (
                f"{held} of {launched} held, {scan.lost_total()} lost,\nevery one for the {periods:.0f} RF periods: "
                f"{verdict}"
            )
            head = "the control"
        else:
            verdict = "pass" if rec.verdict("must fail", row.ident) else "FAIL"
            losses_ = [x for k in range(len(scan.senses)) for x in scan.first_loss[k] if math.isfinite(x)]
            lost = scan.lost_total()
            every = "every electron" if lost == launched else f"{lost} electrons"
            said = (
                f"{held} of {launched} held; each energy and sense's first loss\nafter {min(losses_):.2f} to "
                f"{max(losses_):.2f} RF periods, {every} lost within the {periods:.0f}: {verdict}"
            )
            head = f"must fail, theory q {float(row.basis['q_null']):.2f}"
        ax.set_title(f"{head}: {title(row)}\n{said}", fontsize=9)
    axes[0].set_ylabel("launch energy along the path (both senses)")
    handles = [
        plt.Line2D([], [], color=SERIES[0], marker="o", ms=6, lw=0, label="first lost, with the wave"),
        plt.Line2D([], [], color=SERIES[1], marker="o", ms=6, lw=0, label="first lost, against it"),
        plt.Line2D([], [], color=MUTED, lw=0.6, label="the other three: lost before the end, when not recorded"),
        plt.Line2D([], [], color=MUTED, lw=2.4, label="held to the end"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.01, 0.875), ncol=4, fontsize=8.5)
    fig.suptitle(
        "The case that must fail: past Mathieu's q of 0.908 no electron is held, the first of each energy and sense "
        "lost within about two RF periods,\nwhile the same 32 starts at the design's 50 V (q 0.51) are all held for "
        "200; four electrons an energy and sense, one at each RF phase",
        x=0.01,
        y=0.99,
        ha="left",
        fontsize=10.5,
    )
    fig.subplots_adjust(left=0.07, right=0.985, top=0.68, bottom=0.13, wspace=0.08)
    return _save(plt, fig, out, "must_fail.png")


def _poly3d(polygons: Any, **style: Any) -> Any:
    """matplotlib's Poly3DCollection (its toolkit has no type information)."""
    return importlib.import_module("mpl_toolkits.mplot3d.art3d").Poly3DCollection(polygons, **style)


def _ribbon(ax: Any, path: Route, colours: Any, step: int = 4) -> None:
    """The boards as a ribbon of the copper's width across the path (mm), each quad shaded by its colour."""
    half = path.width / 2 / MM
    centre, across = path.points / MM, path.across
    n = len(centre)
    quads, faces = [], []
    for i in range(0, n, step):
        j = (i + step) % n
        quads.append(
            [
                centre[i] - half * across[i],
                centre[i] + half * across[i],
                centre[j] + half * across[j],
                centre[j] - half * across[j],
            ]
        )
        faces.append(colours[i])
    ax.add_collection3d(_poly3d(quads, facecolors=faces, edgecolors=faces, linewidths=0.4))


HOLDS_DRAWN = (
    ("designs/r1_rail_ring.toml", 0, None),
    ("designs/flat_ring_318.toml", 0, None),
    ("designs/t915s.toml", 0, None),
    ("designs/helix_13.toml", 0, None),
)


def _drawn_holds(rec: Records) -> list[Row]:
    """The holds drawn in three dimensions: each path's own, as drawn at the design's volts."""
    out: list[Row] = []
    for design, realisation, q in HOLDS_DRAWN:
        out.extend(
            row
            for row in rec.of("hold")
            if row.shot["design"] == design and row.shot["realisation"] == realisation and row.shot["q"] == q
        )
    return out


def figure_paths(rec: Records, out: Path) -> Path:
    """R1, the flat ring, T915-S and the helix in three dimensions: the boards as a ribbon shaded by the curvature in
    their plane, the start and the sense of the wave, where the static limit is weakest, and where the first electron
    of each lost energy was lost."""
    plt = _pyplot()
    from matplotlib import cm, colors

    rows = _drawn_holds(rec)
    routes = [route_of(rec, row) for row in rows]
    top = max(float(r.kappa_g.max()) for r in routes)
    scale = colors.Normalize(vmin=0.0, vmax=top * 1.05)
    shade = colors.LinearSegmentedColormap.from_list("ink", ["#ecebe6", "#898781", "#0b0b0b"])
    fig = plt.figure(figsize=(14.5, 12.4))
    for k, (row, path) in enumerate(zip(rows, routes, strict=True)):
        ax: Any = fig.add_subplot(2, 2, k + 1, projection="3d", computed_zorder=False)
        mm = path.points / MM
        reach = float(np.abs(mm).max()) * 1.08
        _ribbon(ax, path, shade(scale(path.kappa_g)))
        ax.plot(*np.vstack([mm, mm[:1]]).T, color=INK, lw=0.5, zorder=2)
        start, ahead = mm[0], path.tangents[0] * reach * 0.28
        ax.quiver(*start, *ahead, color=INK, lw=1.6, arrow_length_ratio=0.25, zorder=6)
        ax.scatter(*start, marker="s", s=40, color=INK, zorder=6, depthshade=False)
        ax.text(*(start + ahead * 1.2), "start", fontsize=8.5, zorder=7)
        weakest = path.at(float(row.basis["static_at"])) / MM
        ax.scatter(*weakest, marker="*", s=260, color=SERIES[3], edgecolor=INK, lw=0.6, zorder=8, depthshade=False)
        h = hold(rec, row)
        for loss in losses(row.scan):
            where = path.at(loss.place) / MM
            first_lap = loss.when < 1.0
            ax.scatter(
                *where,
                s=34 if first_lap else 60,
                color=SENSE_COLOUR[loss.sense] if first_lap else SURFACE,
                edgecolor=SENSE_COLOUR[loss.sense],
                lw=1.5,
                zorder=9 if first_lap else 10,
                depthshade=False,
            )
        ax.set_xlim(-reach, reach)
        ax.set_ylim(-reach, reach)
        depth = max(float(np.abs(mm[:, 2]).max()) * 1.3, reach * 0.18)
        ax.set_zlim(-depth, depth)
        ax.set_box_aspect((1, 1, depth / reach), zoom=1.3)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.set_pane_color((0.985, 0.985, 0.98, 1.0))
            axis._axinfo["grid"]["color"] = GRID
            axis._axinfo["grid"]["linewidth"] = 0.5
        ax.tick_params(labelsize=7, pad=-1)
        ax.set_xlabel("x (mm)", fontsize=8, labelpad=-4)
        ax.set_ylabel("y (mm)", fontsize=8, labelpad=-4)
        ax.view_init(elev=34 if h.planar else 28, azim=-62)
        lost = losses(row.scan)
        ax.set_title(
            f"{title(row)}: {path.length / MM:.1f} mm round, {'planar' if h.planar else 'out of the plane'}\n"
            f"held to {h.limit:.4g} eV, {h.fraction:.3f} of the static limit; {len(lost)} energies and senses lost",
            fontsize=9.2,
            loc="left",
        )
    bar = fig.colorbar(
        cm.ScalarMappable(norm=scale, cmap=shade), cax=fig.add_axes((0.3, 0.06, 0.4, 0.013)), orientation="horizontal"
    )
    bar.set_label("the boards (a ribbon the copper's width) shaded by the curvature in their plane, |κ_g| (1/m)")
    handles = [
        plt.Line2D(
            [], [], color=SERIES[0], marker="o", ms=6, lw=0, label="first loss, with the wave, on its first lap"
        ),
        plt.Line2D([], [], color=SERIES[1], marker="o", ms=6, lw=0, label="against it, on its first lap"),
        plt.Line2D([], [], color=SURFACE, mec=SERIES[0], mew=1.5, marker="o", ms=8, lw=0, label="after a lap or more"),
        plt.Line2D(
            [], [], color=SERIES[3], mec=INK, marker="*", ms=13, lw=0, label="where the static limit is weakest"
        ),
        plt.Line2D([], [], color=INK, marker="s", ms=6, lw=0, label="the start; the arrow, the drive's wave"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.01, 0.955), ncol=5, fontsize=8.5)
    fig.suptitle(
        "Where electrons were lost along each path: one mark for each energy and sense that lost any, at the place its "
        "first-lost electron left, its laps\ntaken in its sense of travel from the start, modulo one lap. The records "
        "hold no position in the section, nor the other lost electrons' places",
        x=0.01,
        y=0.995,
        ha="left",
        fontsize=10.5,
    )
    fig.subplots_adjust(left=0.0, right=1.0, top=0.88, bottom=0.1, wspace=0.0, hspace=0.1)
    return _save(plt, fig, out, "paths.png")


def figure_along(rec: Records, out: Path) -> Path:
    """For every hold, the place along the path of each lost energy's first loss, against its energy over the static
    limit, under the curvature in the boards' plane and toward them."""
    plt = _pyplot()
    rows = rec.of("hold")
    fig = plt.figure(figsize=(15.5, 9.6))
    outer = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.22)
    for k, row in enumerate(rows):
        inner = outer[k // 3, k % 3].subgridspec(2, 1, height_ratios=(1, 3.2), hspace=0.08)
        track = fig.add_subplot(inner[0])
        ax = fig.add_subplot(inner[1], sharex=track)
        path, h = route_of(rec, row), hold(rec, row)
        track.plot(path.s, path.kappa_g, color=INK, lw=1.1, label="in the boards' plane")
        track.plot(path.s, path.kappa_n, color=MUTED, lw=1.1, label="toward them")
        track.set_ylim(0, max(float(path.kappa_g.max()), float(path.kappa_n.max())) * 1.2)
        track.set_ylabel("κ (1/m)", fontsize=8)
        track.tick_params(labelbottom=False, labelsize=7.5)
        weakest = float(row.basis["static_at"])
        for axis in (track, ax):
            axis.axvline(weakest, color=SERIES[3], lw=1.4, zorder=1)
            for a, b in lobes(path):
                axis.axvspan(a, b, color=PALE, lw=0, zorder=0)
        ax.axhline(1.0, color=INK, lw=1.2)
        ax.axhline(h.fraction, color=INK, lw=1.2, ls=(0, (4, 2.5)))
        for loss in losses(row.scan):
            first_lap = loss.when < 1.0
            colour = SENSE_COLOUR[loss.sense]
            ax.plot(
                [loss.place],
                [loss.energy / h.static],
                "o",
                ms=5 if first_lap else 7,
                color=colour if first_lap else SURFACE,
                mec=colour,
                mew=1.4,
                zorder=4,
            )
            if not first_lap:
                ax.text(
                    loss.place + 0.014,
                    loss.energy / h.static,
                    f"{loss.when:.1f}",
                    fontsize=6.5,
                    color=INK2,
                    va="center",
                )
        energies = np.asarray(row.scan.energies) / h.static
        ax.set_ylim(min(0.7, h.fraction - 0.08), min(float(energies[-1]) * 1.02, 1.45))
        ax.set_xlim(0, 1)
        ax.set_xlabel("place along the path (laps from the start)")
        ax.set_ylabel("energy over the static limit")
        track.set_title(f"{title(row)}: held to {h.fraction:.3f} (dashed)", fontsize=9)
    handles = [
        plt.Line2D([], [], color=SERIES[0], marker="o", ms=5, lw=0, label="with the wave, lost on its first lap"),
        plt.Line2D([], [], color=SERIES[1], marker="o", ms=5, lw=0, label="against it, on its first lap"),
        plt.Line2D([], [], color=SURFACE, mec=MUTED, mew=1.4, marker="o", ms=7, lw=0, label="after laps (labelled)"),
        plt.Line2D([], [], color=SERIES[3], lw=1.4, label="where the static limit is weakest"),
        plt.Rectangle((0, 0), 1, 1, color=PALE, label="lobes: between the two maxima of |κ_g| at a tightest bend"),
        plt.Line2D([], [], color=INK, lw=1.1, label="track: curvature in the boards' plane, |κ_g|"),
        plt.Line2D([], [], color=MUTED, lw=1.1, label="toward them, |κ_n|"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.01, 0.94), ncol=4, fontsize=8.5)
    beyond = []
    for row in rows:
        stretches = lobes(route_of(rec, row))
        if stretches:
            worst = max(past(loss, stretches) for loss in losses(row.scan))
            beyond.append(f"{row.shot['label']} {worst:.3f}")
    rings = []
    for row in rows:
        h = hold(rec, row)
        if row.shot["realisation"] == 0 and h.planar:
            edge = {sense.sense: sense.all_lost for sense in h.senses}
            rings += [loss.when for loss in losses(row.scan) if loss.energy >= edge[loss.sense]]
    fig.suptitle(
        "Where each lost energy's first electron left, unrolled. Out of the plane every one was lost in a lobe or "
        "just past it in its direction of travel (at most "
        + " and ".join(beyond)
        + " lap past);\non the rings as drawn, which have no lobes, every one at or above the edge within "
        f"{min(rings):.2f} to {max(rings):.2f} lap of the start. Read after the run: a hypothesis, not a test",
        x=0.01,
        y=0.995,
        ha="left",
        fontsize=10.5,
    )
    fig.subplots_adjust(left=0.055, right=0.985, top=0.86, bottom=0.07)
    return _save(plt, fig, out, "along.png")


def draw(rec: Records, out: Path) -> tuple[Path, ...]:
    """Every figure, into `out`."""
    out.mkdir(parents=True, exist_ok=True)
    return (
        figure_holds(rec, out),
        figure_wells(rec, out),
        figure_fractions(rec, out),
        figure_must_fail(rec, out),
        figure_paths(rec, out),
        figure_along(rec, out),
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument("--records", type=Path, default=SRC, help="the folder holding the committed records")
    parser.add_argument("--out", type=Path, default=SRC, help="where the figures go")
    args = parser.parse_args(argv)
    rec = read(args.records)
    for path in draw(rec, args.out):
        print(path.name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
