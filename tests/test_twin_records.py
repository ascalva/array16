"""array16.twin.records, the detector's entries for the electrons a run loses (H3 of experiments/twin_rung1/REPORT.md
section 6): on R1's ring with no volts, so that every electron flies straight, an electron aimed at the copper is
recorded on the face it struck, one sent out of the aperture on the grid's edge it crossed and one the gas stopped where
it was, each at its arc length and laps, its section position, step, time, energy and sense, and how; on T915-S's knot
the same in its turning frame; on a straight guide, the whole laps gone; on R1 driven, every entry on the boundary its
step crossed, while everything stage 1 records stays bit for bit what it was; the entries through array16.lab's plain
record and back; and with the option off, a bench's stamp and a watch's reading as stage 1's."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from array16 import particles as pt
from array16.lab import Bench, plain
from array16.particles import Arr, Flight
from array16.twin import Detection, Losses, TwinDetector, TwinField, TwinGuide, Watch
from array16.twin.records import HOW, LEFT, STOPPED, STRUCK

ROOT = Path(__file__).resolve().parents[1]
DESIGNS = ROOT / "designs"
STEPS = 120  # steps an RF period, stage 1's
FR4_FACE = -35e-6  # R1's board under an etched gap: the copper's thickness below its top face (m)
FILM_FACE = -18e-6  # T915-S's polyimide under the split: its copper's thickness behind the face (m)
STRAIGHT = """
name = "straight facing strips"
seed = 1
[band]
f = 915e6
m = 1
[section]
kind = "facing"
gap = 1.83e-3
width = 7.32e-3
split = 0.2745e-3
margin = 1.83e-3
laminate = "polyimide"
thickness = 25e-6
copper = 18e-6
grid = "coarse"
[path]
kind = "straight"
length = 0.1
samples = 64
[drive]
volts = 50.0
"""


@pytest.fixture(scope="module")
def coarse(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """R1 and T915-S on the coarse grid, with one cache of solves for the module."""
    where = tmp_path_factory.mktemp("records")
    for name in ("r1_rail_ring", "t915s"):
        text = (DESIGNS / f"{name}.toml").read_text()
        assert 'grid = "default"' in text
        (where / f"{name}.toml").write_text(text.replace('grid = "default"', 'grid = "coarse"'))
    (where / "straight.toml").write_text(STRAIGHT)
    return where


def _guide(coarse: Path, name: str) -> TwinGuide:
    return TwinGuide(coarse / f"{name}.toml", 0, cache=coarse / "cache")


def _still(guide: TwinGuide) -> TwinField:
    """The design's section swept along its path with no volts on it: every electron flies straight."""
    return TwinField(guide.sweep(guide.grid(), pt.InPhase(0.0, 0.0)), guide.curve)


def _momentum(energies: Arr, headings: Arr) -> Arr:
    """gamma v (N, 3) of electrons of kinetic energies (eV) along unit headings (N, 3)."""
    speed = np.array([pt.v_of(float(w)) for w in energies])
    out: Arr = (speed / np.sqrt(1 - (speed / pt.C0) ** 2))[:, None] * headings
    return out


def _kinetic(u: Arr) -> Arr:
    """Kinetic energies (eV) from gamma v (N, 3)."""
    w = (u * u).sum(axis=1) / pt.C0**2
    out: Arr = w / (np.sqrt(1 + w) + 1) * (pt.ME * pt.C0**2 / pt.E)
    return out


def _speed(u: Arr) -> Arr:
    """Speeds (m/s) from gamma v (N, 3)."""
    size = np.sqrt((u * u).sum(axis=1))
    out: Arr = size / np.sqrt(1 + (size / pt.C0) ** 2)
    return out


def _at(losses: Losses, row: int) -> int:
    """Where electron `row`'s entry is among the losses: once, or the test fails."""
    (at,) = np.flatnonzero(losses.index == row)
    return int(at)


@dataclass(frozen=True, eq=False)
class Still:
    """R1's ring with no volts, its launch and what the watch made of it."""

    guide: TwinGuide
    field: TwinField
    x0: Arr
    u0: Arr
    samples: tuple[int, ...]
    dt: float
    flight: Flight
    seen: Detection


@pytest.fixture(scope="module")
def still(coarse: Path) -> Still:
    """Fifteen electrons on R1's ring with no volts, watched with the entries for 260 steps: four sent down from the
    null onto the strip, a rail, the other rail and the copper beyond the outer gap, one down an etched gap, one up,
    one across, two along the path either way, one down and along it, one a third of a step above the copper, one
    launched in it, two the gas stops at step 40 and one at rest at the null."""
    guide = _guide(coarse, "r1_rail_ring")
    field = _still(guide)
    swept = field.swept
    first, second, along = swept.first, swept.second, swept.tangents
    dt = 1 / (STEPS * guide.tune().f)
    v0 = swept.grid.null[1]
    tilt = math.radians(30.0)
    third = pt.v_of(50.0) * dt / 3
    # (sample, offset across (m), offset up from the null (m), heading, energy (eV))
    plan: list[tuple[int, float, float, Arr, float]] = [
        (0, 0.0, 0.0, -second[0], 50.0),
        (250, 0.8e-3, 0.0, -second[250], 40.0),
        (1000, -1.2e-3, 0.0, -second[1000], 60.0),
        (1600, 2.0e-3, 0.0, -second[1600], 45.0),
        (400, 0.2e-3, 0.0, -second[400], 50.0),
        (500, 0.0, 0.0, second[500], 50.0),
        (1500, 0.0, 0.0, first[1500], 50.0),
        (300, 0.0, 0.0, along[300], 100.0),
        (700, 0.0, 0.0, -along[700], 100.0),
        (1200, 0.0, 0.0, math.cos(tilt) * -second[1200] + math.sin(tilt) * along[1200], 50.0),
        (1800, 0.8e-3, third - v0, -second[1800], 50.0),
        (100, 0.8e-3, -10e-6 - v0, -second[100], 50.0),
        (900, 0.0, 0.0, along[900], 0.01),
        (950, 0.0, 0.0, -along[950], 0.01),
        (1100, 0.0, 0.0, along[1100], 0.0),
    ]
    samples = tuple(i for i, *_ in plan)
    x0 = np.array([swept.points[i] + a * first[i] + b * second[i] for i, a, b, _, _ in plan])
    u0 = _momentum(np.array([w for *_, w in plan]), np.array([heading for _, _, _, heading, _ in plan]))
    count = len(plan)
    gassed = np.zeros(count, dtype=np.bool_)
    gassed[[12, 13]] = True

    def gas(step: int, t: float, x: Arr, u: Arr, alive: pt.BoolArr) -> pt.BoolArr:
        return gassed & (step == 40)

    watch = TwinDetector().watch(field, count, start=x0, losses=True)
    flight = pt.integrate(x0, u0, np.full(count, pt.QM_E), dt, 260, field, [gas, watch])
    return Still(guide, field, x0, u0, samples, dt, flight, watch.reading())


# ---- where the entries put the electrons --------------------------------------------------------------------------
def test_an_electron_aimed_at_the_copper_is_recorded_where_it_struck(still: Still) -> None:
    # R1 as drawn, no volts: an electron sent down from the null strikes the copper's top face (v = 0) square under
    # where it set out, at its own arc length, no laps gone and no sense along the path, the null's height over its
    # speed after it set out (13.5 to 16.5 steps), in the step before the one the watch stopped it at, with the energy
    # it was launched with; down an etched gap, the FR-4 35 um below; down and along the path, 30 degrees from the
    # vertical, where its line meets the face, 0.33 mm on along the path and 1.6 um outward. One launched a third of a
    # step above the copper strikes it in its first step, a third of the way through; one launched in the copper is
    # recorded where it was launched, at step 0 and time 0
    guide, swept, dt, seen = still.guide, still.field.swept, still.dt, still.seen
    losses = seen.losses
    assert losses is not None
    grid, lap = swept.grid, guide.length
    u0, v0 = grid.null
    ds, radius = lap / len(swept.points), lap / (2 * math.pi)
    speed, energy = _speed(still.u0), _kinetic(still.u0)
    assert losses.length == lap
    for row, du, face in ((0, 0.0, 0.0), (1, 0.8e-3, 0.0), (2, -1.2e-3, 0.0), (3, 2.0e-3, 0.0), (4, 0.2e-3, FR4_FACE)):
        k = _at(losses, row)
        fall = (v0 - face) / speed[row]
        assert losses.how[k] == STRUCK
        assert losses.v[k] == pytest.approx(face, abs=1e-11)
        assert losses.u[k] == pytest.approx(u0 + du, abs=1e-11)
        assert losses.s[k] == pytest.approx(still.samples[row] * ds, abs=1e-10)
        assert losses.laps[k] == pytest.approx(0.0, abs=1e-12)
        assert losses.time[k] == pytest.approx(fall, abs=1e-9 * dt)
        assert losses.step[k] == math.floor(fall / dt) == seen.step[row] - 1
        assert losses.energy[k] == pytest.approx(energy[row], rel=1e-12)
        assert losses.sense[k] == 0
        assert np.allclose(
            losses.point[k],
            still.x0[row] + still.u0[row] / np.linalg.norm(still.u0[row]) * speed[row] * fall,
            atol=1e-11,
        )
    assert 13 < min(losses.time[_at(losses, row)] / dt for row in range(5)) < 14
    k = _at(losses, 9)  # down and along: where the line meets the face
    fall = v0 / (speed[9] * math.cos(math.radians(30.0)))
    meet = still.x0[9] + still.u0[9] / np.linalg.norm(still.u0[9]) * speed[9] * fall
    assert losses.how[k] == STRUCK
    assert losses.v[k] == pytest.approx(0.0, abs=1e-11)
    assert losses.time[k] == pytest.approx(fall, abs=1e-9 * dt)
    assert np.allclose(losses.point[k], meet, atol=1e-11)
    on = np.mod(radius * math.atan2(meet[1], meet[0]) - still.samples[9] * ds + lap / 2, lap) - lap / 2
    assert losses.s[k] - still.samples[9] * ds == pytest.approx(on, abs=1e-10)
    assert on == pytest.approx(0.326e-3, rel=0.01)
    assert losses.u[k] == pytest.approx(u0 + radius - math.hypot(meet[0], meet[1]), abs=1e-11)
    assert losses.laps[k] == pytest.approx(on / lap, abs=1e-10)
    assert losses.sense[k] == 1
    k = _at(losses, 10)  # a third of a step above the copper: struck a third of the way through the first step
    assert (losses.how[k], losses.step[k], seen.step[10]) == (STRUCK, 0, 1)
    assert losses.time[k] == pytest.approx(dt / 3, abs=1e-9 * dt)
    assert losses.v[k] == pytest.approx(0.0, abs=1e-11)
    k = _at(losses, 11)  # in the copper from its launch
    assert (losses.how[k], losses.step[k], losses.time[k], seen.step[11]) == (STRUCK, 0, 0.0, 0)
    assert np.array_equal(losses.point[k], still.x0[11])
    assert losses.v[k] == pytest.approx(-10e-6, abs=1e-11)
    assert losses.u[k] == pytest.approx(u0 + 0.8e-3, abs=1e-11)


def test_an_electron_sent_out_of_the_aperture_is_recorded_on_the_edge_it_crossed(still: Still) -> None:
    # up from the null it leaves through the grid's top, 2.20 mm over the null; across, through its inner edge, 2.23 mm
    # in; along the path, off the ring's tangent, through its outer edge where the tangent is 2.23 mm out from the
    # ring, 12.43 mm from where it set out, 11.90 mm on along the path, 0.0565 lap, with the wave's sense (+1) or
    # against it (-1)
    guide, swept, dt, seen = still.guide, still.field.swept, still.dt, still.seen
    losses = seen.losses
    assert losses is not None
    grid, lap = swept.grid, guide.length
    u0, v0 = grid.null
    ds, radius = lap / len(swept.points), lap / (2 * math.pi)
    speed = _speed(still.u0)
    for row, place, flight_length in ((5, (u0, grid.v[-1]), grid.v[-1] - v0), (6, (grid.u[-1], v0), grid.u[-1] - u0)):
        k = _at(losses, row)
        assert losses.how[k] == LEFT
        assert (losses.u[k], losses.v[k]) == pytest.approx(place, abs=1e-11)
        assert losses.s[k] == pytest.approx(still.samples[row] * ds, abs=1e-10)
        assert losses.time[k] == pytest.approx(flight_length / speed[row], abs=1e-9 * dt)
        assert losses.step[k] == seen.step[row] - 1
    out = radius + u0 - grid.u[0]  # u runs inward on R1: the outer edge is this far from the ring's centre
    reach = math.sqrt(out * out - radius * radius)
    on = radius * math.atan2(reach, radius)
    assert (reach, on, on / lap) == pytest.approx((12.43e-3, 11.90e-3, 0.0565), rel=1e-3)
    for row, sense in ((7, 1), (8, -1)):
        k = _at(losses, row)
        assert losses.how[k] == LEFT
        assert (losses.u[k], losses.v[k]) == pytest.approx((grid.u[0], v0), abs=1e-11)
        assert np.mod(losses.s[k] - still.samples[row] * ds - sense * on + lap / 2, lap) - lap / 2 == pytest.approx(
            0.0, abs=1e-10
        )
        assert losses.laps[k] == pytest.approx(sense * on / lap, abs=1e-10)
        assert losses.time[k] == pytest.approx(reach / speed[row], abs=1e-9 * dt)
        assert losses.sense[k] == sense


def test_the_gas_stops_an_electron_where_it_was_and_the_held_have_no_entry(still: Still) -> None:
    # two electrons at 0.01 eV along the path either way, stopped by the gas after step 40, are recorded where they
    # stopped, at the end of that step (41 dt), with their energy and sense; the electron at rest at the null is held
    # and has no entry; one entry an electron, in the order they were lost
    swept, dt, seen, flight = still.field.swept, still.dt, still.seen, still.flight
    losses = seen.losses
    assert losses is not None
    lap = still.guide.length
    for row, sense in ((12, 1), (13, -1)):
        k = _at(losses, row)
        assert (losses.how[k], losses.step[k], seen.step[row], seen.fate[row]) == (STOPPED, 40, 40, 2)
        assert losses.time[k] == pytest.approx(41 * dt, rel=1e-12)
        assert np.array_equal(losses.point[k], flight.x[row])
        assert losses.s[k] == swept.locate(flight.x[row : row + 1]).s[0]
        gone = 41 * dt * _speed(still.u0)[row]
        assert losses.laps[k] == pytest.approx(sense * gone / lap, rel=1e-6)
        assert losses.energy[k] == pytest.approx(_kinetic(still.u0)[row], rel=1e-12)
        assert losses.energy[k] == pytest.approx(0.01, rel=1e-6)
        assert losses.sense[k] == sense
    assert seen.fate[14] == 0
    assert sorted(losses.index.tolist()) == list(range(14))
    assert len(losses) == seen.lost + seen.stopped == 14
    assert set(losses.how.tolist()) == set(HOW)
    assert np.all(np.diff(losses.time) >= 0)
    assert all(not getattr(losses, name).flags.writeable for name in ("index", "how", "time", "point"))


def test_on_a_knot_the_entry_is_where_the_section_met_its_board(coarse: Path) -> None:
    # T915-S, its boards turning with the trefoil (k = 0), no volts: electrons sent from the null square across the
    # section, at four places round the knot, strike the lower board's copper face (v = 0), the upper's (v = 2h,
    # 1.83 mm) or, through the split, the polyimide behind the lower board's copper (v = -18 um), each at its own arc
    # length with no laps gone and no sense along the path, h / speed after it set out (31 steps at 30 eV)
    guide = _guide(coarse, "t915s")
    field = _still(guide)
    swept = field.swept
    grid, lap = swept.grid, guide.length
    u0, v0 = grid.null
    ds = lap / len(swept.points)
    dt = 1 / (STEPS * guide.tune().f)
    plan = [(400, 1.0e-3, -1, 0.0), (1700, -1.5e-3, 1, 2 * v0), (2900, 0.0, -1, FILM_FACE), (3300, 2.5e-3, -1, 0.0)]
    x0 = np.array([swept.points[i] + du * swept.first[i] for i, du, _, _ in plan])
    u_start = _momentum(np.full(len(plan), 30.0), np.array([way * swept.second[i] for i, _, way, _ in plan]))
    watch = TwinDetector(losses=True).watch(field, len(plan), start=x0)
    pt.integrate(x0, u_start, np.full(len(plan), pt.QM_E), dt, 80, field, [watch])
    seen = watch.reading()
    losses = seen.losses
    assert losses is not None
    assert seen.lost == len(losses) == 4
    speed = _speed(u_start)
    for row, (i, du, _, face) in enumerate(plan):
        k = _at(losses, row)
        assert losses.how[k] == STRUCK
        assert (losses.u[k], losses.v[k]) == pytest.approx((u0 + du, face), abs=1e-11)
        assert losses.s[k] == pytest.approx(i * ds, abs=1e-10)
        assert losses.laps[k] == pytest.approx(0.0, abs=1e-12)
        assert losses.time[k] == pytest.approx(abs(v0 - face) / speed[row], abs=1e-9 * dt)
        assert losses.sense[k] == 0
    assert 30 < losses.time[_at(losses, 0)] / dt < 31


def test_on_a_straight_guide_the_entry_counts_the_whole_laps_gone(coarse: Path) -> None:
    # T915-S's section on a straight path a lap of which is 0.1 m (a ring unrolled), no volts: electrons sent from the
    # null 1 mm across, along the path and down toward the lower board at a slope of h in 0.15 m, strike its copper
    # face after 1.5 laps, with the path or against it: halfway round a lap (s = 0.05 m), +1.5 and -1.5 laps gone
    guide = _guide(coarse, "straight")
    assert guide.curve is None
    field = _still(guide)
    swept = field.swept
    grid, lap = swept.grid, guide.length
    u0, v0 = grid.null
    dt = 1 / (STEPS * guide.tune().f)
    slope = math.atan2(v0, 1.5 * lap)
    x0 = np.array([[1e-3, 0.0, 0.0], [1e-3, 0.0, 0.0]])
    headings = np.array([[0.0, -math.sin(slope), math.cos(slope)], [0.0, -math.sin(slope), -math.cos(slope)]])
    u_start = _momentum(np.full(2, 100.0), headings)
    flight_time = math.hypot(1.5 * lap, v0) / float(_speed(u_start)[0])
    watch = TwinDetector().watch(field, 2, start=x0, losses=True)
    pt.integrate(x0, u_start, np.full(2, pt.QM_E), dt, math.ceil(flight_time / dt) + 3, field, [watch])
    losses = watch.reading().losses
    assert losses is not None
    for row, sense in ((0, 1), (1, -1)):
        k = _at(losses, row)
        assert losses.how[k] == STRUCK
        assert (losses.u[k], losses.v[k]) == pytest.approx((u0 + 1e-3, 0.0), abs=1e-11)
        assert losses.s[k] == pytest.approx(0.5 * lap, abs=1e-10)
        assert losses.laps[k] == pytest.approx(1.5 * sense, abs=1e-10)
        assert losses.time[k] == pytest.approx(flight_time, rel=1e-12)  # 2780 steps added up
        assert losses.sense[k] == sense


def test_without_the_launch_a_first_step_loss_is_placed_where_the_field_found_it(coarse: Path) -> None:
    # a watch not given the launch (start) never sees the first step's segment: an electron a third of a step above
    # R1's copper is placed where the field first found it, two thirds of a step into the copper, at the end of step
    # 0; one that crosses later is placed on the face, five and a half steps on, as ever, its laps counted from where
    # the first step left it
    guide = _guide(coarse, "r1_rail_ring")
    field = _still(guide)
    swept = field.swept
    dt = 1 / (STEPS * guide.tune().f)
    v0 = swept.grid.null[1]
    hop = pt.v_of(50.0) * dt
    first, second = swept.first, swept.second
    x0 = np.array(
        [swept.points[i] + 0.8e-3 * first[i] + (b - v0) * second[i] for i, b in ((10, hop / 3), (20, 5.5 * hop))]
    )
    u_start = _momentum(np.full(2, 50.0), -second[[10, 20]])
    watch = TwinDetector().watch(field, 2, losses=True)
    pt.integrate(x0, u_start, np.full(2, pt.QM_E), dt, 20, field, [watch])
    losses = watch.reading().losses
    assert losses is not None
    k = _at(losses, 0)
    assert (losses.how[k], losses.step[k]) == (STRUCK, 0)
    assert losses.time[k] == pytest.approx(dt, rel=1e-12)
    assert losses.v[k] == pytest.approx(-2 * hop / 3, rel=1e-6)
    k = _at(losses, 1)
    assert (losses.how[k], losses.step[k]) == (STRUCK, 5)
    assert losses.time[k] == pytest.approx(5.5 * dt, abs=1e-9 * dt)
    assert losses.v[k] == pytest.approx(0.0, abs=1e-11)
    assert losses.laps[k] == pytest.approx(0.0, abs=1e-12)


def test_a_driven_hold_records_each_loss_on_the_boundary_its_step_crossed(coarse: Path) -> None:
    # R1 as drawn at 50 V for 5 RF periods, stage 1's detector test: 20 electrons at rest within 0.05h of the null are
    # held, 20 sent up at 5 eV leave the aperture's top and 20 sent down at 100 eV strike the copper. Watched with the
    # entries and without, the flights and every reading stage 1 records are the same, bit for bit. Each entry is on the
    # boundary its last step crossed: a nanometre back along the step on the grid, a nanometre on off it, in a solid
    # cell for those that struck and beyond the grid for those that left; in the step before the one the watch stopped
    # it at, within two steps of where it stopped, with the energy it stopped with, its laps within 1e-5 of the watch's
    guide = _guide(coarse, "r1_rail_ring")
    curve = guide.curve
    assert curve is not None
    h = guide.solved.distance
    count, steps = 60, 600
    dt = 1 / (STEPS * guide.tune().f)
    rng = np.random.default_rng(3)
    x0 = curve.points[rng.integers(0, len(curve), count)] + rng.normal(0, 0.03 * h, (count, 3))
    u0 = np.zeros((count, 3))
    u0[20:40, 2], u0[40:, 2] = pt.v_of(5.0), -pt.v_of(100.0)
    runs = []
    for keep in (False, True):
        field = guide.field()
        watch = TwinDetector().watch(field, count, start=x0, losses=keep)
        flight = pt.integrate(x0, u0, np.full(count, pt.QM_E), dt, steps, field, [watch])
        runs.append((field, flight, watch.reading()))
    (_, bare_flight, bare), (field, flight, seen) = runs
    assert bare.losses is None
    for name in ("x", "u", "alive"):
        assert np.array_equal(getattr(flight, name), getattr(bare_flight, name))
    for name in ("held", "lost", "stopped", "fate", "step", "time", "at", "laps", "times", "currents"):
        assert np.array_equal(getattr(seen, name), getattr(bare, name), equal_nan=True)
    assert (seen.dump_charge, seen.dump_volts, seen.loss_current) == (
        bare.dump_charge,
        bare.dump_volts,
        bare.loss_current,
    )
    losses = seen.losses
    assert losses is not None
    assert (seen.held, seen.lost, seen.stopped) == (20, 40, 0)
    assert sorted(losses.index.tolist()) == list(range(20, 60))
    order = np.argsort(losses.index)
    assert losses.how[order].tolist() == [LEFT] * 20 + [STRUCK] * 20
    swept, grid = field.swept, field.swept.grid
    rows = losses.index
    heading = flight.u[rows] / np.linalg.norm(flight.u[rows], axis=1)[:, None]
    for side, expect_on in ((-1e-9, True), (1e-9, False)):
        where = swept.locate(losses.point + side * heading)
        _, _, on = grid.gradient_at(where.u, where.v)
        assert (on == expect_on).all()
        beyond = (where.u < grid.u[0]) | (where.u > grid.u[-1]) | (where.v < grid.v[0]) | (where.v > grid.v[-1])
        if not expect_on:
            assert np.array_equal(beyond, losses.how == LEFT)
    left = losses.how == LEFT
    edge = np.minimum.reduce(
        [
            np.abs(losses.u - grid.u[0]),
            np.abs(losses.u - grid.u[-1]),
            np.abs(losses.v - grid.v[0]),
            np.abs(losses.v - grid.v[-1]),
        ]
    )
    assert (edge[left] < 1e-11).all()
    assert (np.abs(losses.v[~left]) < 1e-11).all()  # every one sent down struck the copper's top face
    assert np.array_equal(losses.step, seen.step[rows] - 1)
    assert ((losses.time >= losses.step * dt) & (losses.time <= (losses.step + 1) * dt)).all()
    assert (np.linalg.norm(losses.point - flight.x[rows], axis=1) <= 2 * _speed(flight.u[rows]) * dt * (1 + 1e-9)).all()
    assert np.allclose(losses.energy, _kinetic(flight.u[rows]), rtol=1e-12, atol=0)
    assert np.abs(losses.laps - seen.laps[rows]).max() < 1e-5


# ---- the record, through array16.lab and back ---------------------------------------------------------------------
_NUMBER = st.one_of(st.floats(-1e3, 1e3, allow_nan=False), st.just(math.nan))


@st.composite
def _losses(draw: st.DrawFn) -> Losses:
    count = draw(st.integers(0, 5))

    def column(values: st.SearchStrategy[Any]) -> list[Any]:
        return draw(st.lists(values, min_size=count, max_size=count))

    return Losses(
        index=np.array(column(st.integers(0, 10**6)), dtype=np.intp),
        how=np.array(column(st.sampled_from(HOW)), dtype=np.str_),
        step=np.array(column(st.integers(0, 10**7)), dtype=np.intp),
        time=np.array(column(_NUMBER)),
        s=np.array(column(_NUMBER)),
        laps=np.array(column(_NUMBER)),
        u=np.array(column(_NUMBER)),
        v=np.array(column(_NUMBER)),
        energy=np.array(column(_NUMBER)),
        sense=np.array(column(st.sampled_from((-1, 0, 1))), dtype=np.intp),
        point=np.array([draw(st.lists(_NUMBER, min_size=3, max_size=3)) for _ in range(count)]).reshape(count, 3),
        length=draw(st.floats(1e-3, 10.0)),
    )


@settings(max_examples=60, deadline=None)
@given(_losses())
def test_any_entries_come_back_from_their_plain_record(losses: Losses) -> None:
    # array16.lab keeps a measurement as plain data, one JSON line of runs.jsonl, and reads it back with the
    # experiment's load: the entries come back as they were, a number that was not finite as nan, each column of its
    # own kind
    back = Losses.load(json.loads(json.dumps(plain(losses))))
    assert plain(back) == plain(losses)
    for name in ("time", "s", "laps", "u", "v", "energy", "point"):
        assert np.array_equal(getattr(back, name), getattr(losses, name), equal_nan=True)
    assert all(getattr(back, name).dtype == np.dtype(np.intp) for name in ("index", "step", "sense"))
    assert back.how.dtype.kind == "U"
    assert back.point.shape == (len(losses), 3)
    assert back.length == losses.length


@dataclass(frozen=True, eq=False)
class Hold:
    """A stand-in for a stage-2 hold's measurement: its energies and each RF phase's entries."""

    energies: tuple[float, ...]
    losses: tuple[Losses, ...]

    @classmethod
    def load(cls, record: Mapping[str, Any]) -> Hold:
        return cls(tuple(float(e) for e in record["energies"]), tuple(Losses.load(r) for r in record["losses"]))


def test_a_holds_entries_come_back_from_its_record(still: Still) -> None:
    # a measurement carrying a run's entries, kept as array16.lab keeps it (plain, a JSON line) and read back as an
    # experiment's load would: every column as it was, the run's 14 entries and a phase that lost none
    losses = still.seen.losses
    assert losses is not None
    none = TwinDetector(losses=True).watch(still.field, 2).reading().losses
    assert none is not None
    assert len(none) == 0
    hold = Hold((1.0, 1.025), (losses, none))
    line = json.dumps({"ident": "a hold", "measurement": plain(hold)})
    back = Hold.load(json.loads(line)["measurement"])
    assert plain(back) == plain(hold)
    assert back.losses[1].point.shape == (0, 3)
    for name in ("index", "how", "step", "time", "s", "laps", "u", "v", "energy", "sense", "point"):
        assert np.array_equal(getattr(back.losses[0], name), getattr(losses, name))


def test_entries_say_what_they_take_when_they_are_refused() -> None:
    one: dict[str, Any] = {
        "index": [0],
        "how": [STRUCK],
        "step": [3],
        "time": [1e-11],
        "s": [0.1],
        "laps": [0.5],
        "u": [0.0],
        "v": [0.0],
        "energy": [1.0],
        "sense": [1],
        "point": [[0.0, 0.0, 0.0]],
        "length": 0.2,
    }
    assert len(Losses.load(one)) == 1
    with pytest.raises(ValueError, match=r"one of \['struck', 'left', 'stopped'\]"):
        Losses.load({**one, "how": ["burnt"]})
    with pytest.raises(ValueError, match="one entry an electron: time"):
        Losses.load({**one, "time": [1e-11, 2e-11]})
    with pytest.raises(ValueError, match=r"\+1, -1 or 0"):
        Losses.load({**one, "sense": [2]})
    with pytest.raises(ValueError, match="positive number of metres"):
        Losses.load({**one, "length": 0.0})


# ---- the option --------------------------------------------------------------------------------------------------
def test_with_the_option_off_the_records_are_stage_1s(coarse: Path) -> None:
    # stage 1's bench stamps its detector as stage 1's records hold it, nothing about entries in it; the detector keeps
    # none by default, nor does the watch twin_rung1 makes; and the reading's other fields are stage 1's (the driven
    # hold above holds them bit for bit with the entries and without)
    record = json.loads((ROOT / "experiments" / "twin_rung1" / "runs.jsonl").read_text().splitlines()[0])
    stamped = record["provenance"]["bench"]["devices"]["detector"]
    bench = Bench.named("twin_r1")
    with bench.session({"detector": ("watch", "identity")}) as devices:
        about = bench.describe(devices)
        detector = devices["detector"]
    assert json.loads(json.dumps(plain(about["devices"]["detector"]))) == stamped
    assert set(stamped["identity"]) == {"pickups", "collector_F", "weight"}
    assert detector.losses is False
    field = _still(_guide(coarse, "r1_rail_ring"))
    start = field.swept.points[:3].copy()
    assert detector.watch(field, 3, start=start).reading().losses is None
    stage1 = Watch(field, 3, (), -pt.E, 1.0, 1e-12, 1, start)  # twin_rung1._watch's own
    assert stage1.reading().losses is None


def test_a_bench_or_a_shot_asks_for_the_entries(tmp_path: Path, coarse: Path) -> None:
    # `losses = true` in a bench file's detector makes every watch keep the entries and stamps the records so; a shot
    # asks for them of its own watch, or declines them, whatever the bench says
    path = tmp_path / "twin_records.toml"
    path.write_text(
        'name = "twin-records"\n\n[device.detector]\ndriver = "array16.twin:TwinDetector"\n'
        "pickups = [0.0, 0.25]\ncollector = 1e-12\nlosses = true\n"
    )
    bench = Bench.load(path)
    with bench.session({"detector": ("watch", "identity")}) as devices:
        about = bench.describe(devices)
        detector = devices["detector"]
    assert about["devices"]["detector"]["identity"]["losses"] is True
    field = _still(_guide(coarse, "r1_rail_ring"))
    swept = field.swept
    start = (swept.points[0] - 0.5 * swept.grid.null[1] * swept.second[0])[None, :]
    down = _momentum(np.array([50.0]), -swept.second[:1])
    dt = 1 / (STEPS * 915e6)
    for asked, kept in ((None, True), (False, False)):
        watch = detector.watch(field, 1, start=start, losses=asked)
        pt.integrate(start, down, np.array([pt.QM_E]), dt, 20, field, [watch])
        losses = watch.reading().losses
        assert (losses is not None) == kept
        if losses is not None:
            assert losses.how.tolist() == [STRUCK]
    shot = TwinDetector().watch(field, 1, start=start, losses=True)
    pt.integrate(start, down, np.array([pt.QM_E]), dt, 20, field, [shot])
    assert shot.reading().losses is not None
    with pytest.raises(ValueError, match="true or false"):
        TwinDetector(losses="yes")  # type: ignore[arg-type]
