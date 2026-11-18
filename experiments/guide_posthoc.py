"""Experiment 1, read after the run: where the first-order theory missed, and one term that accounts for most of it.

Everything here was written after the tracking of experiments/guide_study.py had been seen. Nothing in it is a
test of anything: it is where the next experiment's hypothesis comes from (GDE-H2 in docs/LAWS.md).

    uv run python experiments/guide_posthoc.py     # guide/posthoc.md and guide/posthoc.png
"""

from __future__ import annotations

import cmath
import math
import sys
from dataclasses import dataclass
from functools import cache

import guide_study as study
import numpy as np

from array16 import lab

Row = tuple[study.Shot, study.Prediction, study.Measurement | None]


@cache
def rows() -> tuple[Row, ...]:
    found = lab.rows(study.GuideExperiment())
    return tuple(found)


def _wrap(angle: float) -> float:
    return cmath.phase(cmath.exp(1j * angle))


@dataclass(frozen=True)
class Pair:
    """One gyroradius on one guide, launched on the axis both ways: how far each fitted turn is from the
    predicted one, and that pair of shifts split into the part the two share and the part that changes sign."""

    guide: str
    rho: float  # m
    along: float  # fitted minus predicted turn, sent along the field (rad)
    against: float  # the same, sent against it
    candidate: float  # rho x bending / 2 (rad)

    @property
    def shared(self) -> float:
        return (self.along + self.against) / 2

    @property
    def reversing(self) -> float:
        return (self.along - self.against) / 2


@cache
def pairs() -> tuple[Pair, ...]:
    """Every gyroradius at which both on-axis shots were held and gave a circle to fit."""
    shift: dict[tuple[str, float], float] = {}
    for shot, p, m in rows():
        fit = m.fit() if m is not None and m.held else None
        if shot.group != "main" or shot.launch != "axis" or fit is None or p.turn is None or p.radius is None:
            continue
        if p.radius >= study.FIT_RADIUS:
            shift[(shot.guide, shot.rho)] = _wrap(fit.turn - p.turn)
    out = []
    for (name, rho), along in sorted(shift.items()):
        if rho > 0 and (name, -rho) in shift:
            bending = study.GUIDES[name].curve.bending
            out.append(Pair(name, rho, along, shift[(name, -rho)], rho * bending / 2))
    return tuple(out)


def markdown() -> str:
    lines = [
        "<!-- Written by `experiments/guide_posthoc.py`, after the run. Exploratory: not a test. -->",
        "",
        "# Experiment 1, read after the run: the turn per lap",
        "",
        "Fitted minus predicted turn per lap, for the on-axis shots held in both directions. The shared part is",
        "what the two directions have in common, the reversing part what changes sign with the direction. The",
        "candidate is ρ × (the curve's bending) / 2, from the curve alone, with nothing fitted.",
        "",
    ]
    said = [
        [
            pair.guide,
            f"{pair.rho * 1e3:.2f}",
            f"{math.degrees(pair.along):+.2f}°",
            f"{math.degrees(pair.against):+.2f}°",
            f"{math.degrees(pair.shared):+.2f}°",
            f"{math.degrees(pair.candidate):+.2f}°",
            f"{pair.shared / pair.candidate:.2f}",
            f"{math.degrees(pair.reversing):+.2f}°",
        ]
        for pair in pairs()
    ]
    columns = ["Guide", "ρ (mm)", "Along", "Against", "Shared", "Candidate", "Shared / candidate", "Reversing"]
    lines += lab.table(columns, said)
    lines += [
        "",
        "| Guide | Length (mm) | Total curvature / 2π | Bending (1/m) | Bending × length / 4π² |",
        "|" + " --- |" * 5,
    ]
    for name in study.MAIN:
        curve = study.GUIDES[name].curve
        lines.append(
            f"| {name} | {curve.length * 1e3:.0f} | {curve.total_curvature / (2 * math.pi):.2f} | {curve.bending:.1f} "
            f"| {curve.bending * curve.length / (4 * math.pi**2):.2f} |"
        )
    return "\n".join(lines) + "\n"


def figure() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colour = {name: f"C{i}" for i, name in enumerate(study.MAIN)}
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    wall = study.WALL * 1e3
    for shot, p, m in rows():
        if shot.group == "main" and m is not None and p.reach is not None and p.reach < 1.6 * study.WALL:
            mark = "o" if m.held else "x"
            ax[0].plot(p.reach * 1e3, m.reach * 1e3, mark, color=colour[shot.guide], ms=5, alpha=0.8)
    line = np.array([0.0, 1.6 * wall])
    ax[0].plot(line, line, "k-", lw=0.8)
    ax[0].fill_between(line, 0.85 * line, 1.15 * line, color="0.85", zorder=0)
    ax[0].axhline(wall, color="k", ls=":", lw=0.8)
    ax[0].axvline(wall, color="k", ls=":", lw=0.8)
    ax[0].set_xlabel("predicted reach, first order (mm)")
    ax[0].set_ylabel("tracked reach (mm)")
    ax[0].set_title("Reach: circles held, crosses lost; band ±15 %; dotted, the wall", fontsize=9)
    for name in study.MAIN:
        ax[0].plot([], [], "o", color=colour[name], label=name)
    ax[0].legend(fontsize=7)

    for name in study.MAIN[1:]:
        mine = [pair for pair in pairs() if pair.guide == name]
        ax[1].plot(
            [pair.rho * 1e3 for pair in mine], [math.degrees(pair.shared) for pair in mine], "o", color=colour[name]
        )
        rho = np.array([0.0, 2.1e-3])
        ax[1].plot(
            rho * 1e3, np.degrees(rho * study.GUIDES[name].curve.bending / 2), "-", color=colour[name], lw=1, label=name
        )
    ax[1].set_xlabel("gyroradius ρ (mm)")
    ax[1].set_ylabel("turn per lap beyond the predicted, shared by both directions (°)")
    ax[1].set_title("Found after the run: points tracked; lines ρ × bending / 2, nothing fitted", fontsize=9)
    ax[1].legend(fontsize=7)

    for shot, p, m in rows():
        if (shot.guide, shot.launch, shot.group) == ("helix", "axis", "main") and abs(shot.rho - 2e-3) < 1e-9:
            assert m is not None
            assert p.centre is not None
            assert p.radius is not None
            z = np.array(m.section) * 1e3
            ax[2].plot(z.real, z.imag, "o-", color="C3", ms=4, lw=0.5, label="tracked, lap by lap")
            turn = np.exp(1j * np.linspace(0, 2 * math.pi, 200))
            first = (p.centre + p.radius * turn) * 1e3
            ax[2].plot(first.real, first.imag, "k--", lw=1, label="first-order circle")
            fit = m.fit()
            assert fit is not None
            fitted = (fit.centre + fit.radius * turn) * 1e3
            ax[2].plot(fitted.real, fitted.imag, "-", color="C3", lw=0.8, label="circle fitted to the crossings")
    ax[2].set_aspect("equal")
    ax[2].set_xlabel("offset along e1 (mm)")
    ax[2].set_ylabel("offset along e2 (mm)")
    ax[2].set_title("Helix, ρ = 2 mm, launched on the axis: where it crosses the start", fontsize=9)
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(study.OUT / "posthoc.png", dpi=130)
    plt.close(fig)


def main() -> int:
    (study.OUT / "posthoc.md").write_text(markdown())
    figure()
    for pair in pairs():
        print(
            f"{pair.guide:16s} rho {pair.rho * 1e3:4.2f} mm: shared {math.degrees(pair.shared):+6.2f} deg, candidate "
            f"{math.degrees(pair.candidate):+6.2f} deg, reversing {math.degrees(pair.reversing):+6.2f} deg"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
