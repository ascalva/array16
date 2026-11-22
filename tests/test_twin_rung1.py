"""experiments/twin_rung1.py, experiment 2's calibration rung (docs/KNOTTED_GUIDE.md section 10): its theory held to
closed forms (Mathieu's exponent; the static balance on ELE-2's gapless boards and on RailPair's exact field, ELE-1;
the section model on a harmonic well against Mathieu's stability edge; Earth's field along a flat ring's axis leaving
the two senses' geometric mean where it was, and the bench's bound on the senses recomputed); a scan's thresholds; its
shots against the predictions file, a hold's steps covering its laps; measure() on one tiny shot of each kind, on
coarse copies of the designs, so that the run cannot fail on a bug; and its verdicts on scans made up to pass and to
fail. No test tracks long enough, or near enough a registered threshold, to bear on a verdict: the protocol's section
6 lists what they run."""

from __future__ import annotations

import json
import math
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import twin_rung1 as rung

from array16 import particles as pt
from array16.lab import Refusal, plain
from array16.particles import FacingBoards, RailPair, SectionField, SectionGrid
from array16.twin import TwinDetector, TwinField, TwinGuide, TwinVacuum

ROOT = Path(__file__).resolve().parents[1]
F915 = 915e6


def _coarse(name: str, source: str, where: Path) -> Path:
    """A design file on the coarse grid, beside the test."""
    text = (ROOT / source).read_text()
    assert 'grid = "default"' in text
    out = where / f"{name}.toml"
    out.write_text(text.replace('grid = "default"', 'grid = "coarse"'))
    return out


@pytest.fixture(scope="module")
def coarse(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """R1, T915-S and the straight slot guide on the coarse grid, with one cache for the module."""
    where = tmp_path_factory.mktemp("rung1")
    for name, source in (("r1", rung.R1), ("t915s", rung.T915S), ("slot", rung.SLOT)):
        _coarse(name, source, where)
    return where


# ---- the theory against closed forms ------------------------------------------------------------------------------
def test_mathieus_exponent_is_the_continued_fractions_root() -> None:
    # DR-0's exact beta at q = 0.3, the knotted study's at q = 0.5 (run.log), and the first order q / sqrt 2 for small q
    assert rung.mathieu_beta(0.3) == pytest.approx(0.216059, abs=1e-6)
    assert rung.mathieu_beta(0.5) == pytest.approx(0.3737, abs=1e-4)
    assert rung.mathieu_beta(1e-3) == pytest.approx(1e-3 / math.sqrt(2), rel=1e-6)
    # the first stability region's edge at a = 0, q = 0.90805, where beta reaches 1: the static balance's cap
    assert abs(rung.QCAP - 0.90805) < 1e-3
    assert rung.mathieu_beta(0.90) > 0.9


def _span(lo: float, hi: float, step: float) -> pt.Arr:
    return np.asarray(np.arange(round(lo / step), round(hi / step) + 1) * step, dtype=np.float64)


def test_the_static_balance_meets_the_gapless_boards() -> None:
    # ELE-2's facing boards sampled h/40 apart, as the twin's gapless model is: tilted across their plane the well holds
    # at most 1.209 D/h (D tanh^2(pi x/2h), section 6.1), and pushed toward a board it holds until the local q,
    # q0 sec^2(pi y/2h), reaches 0.908: at y = 0.4677h for q 0.5 at the null, with the force D (pi/h) tan sec^2 there,
    # 5.15 D/h. The walk reads the local q to 0.2 % and stops within its step of the crossing; the equilibria walk out
    # from the null
    h = 0.915e-3
    boards = FacingBoards(2 * h)
    grid = SectionGrid.sampled(
        lambda pu, pv: boards.potential(pu, pv - h, 1.0),
        _span(-4 * h, 4 * h, h / 40),
        _span(0.0, 2 * h, h / 40),
        (0, h),
    )
    patch = rung.well_patch(grid, F915, h)
    volts = 0.5 / pt.mathieu_q(boards.gradient(1.0), F915)
    depth = pt.ponderomotive_ev(boards.barrier_field(volts), F915)
    across = rung.valley(patch, 0.0)
    assert across.strength(volts) * h / depth == pytest.approx(1.209, rel=5e-3)
    assert across.limit(volts)[1] == "fold"
    toward = rung.valley(patch, math.pi / 2)
    end, why = toward.limit(volts)
    assert why == "q"
    crossing = 2 / math.pi * math.acos(math.sqrt(0.5 / 0.908))  # q0 sec^2(pi y / 2h) = 0.908: 0.4677
    assert abs(toward.v[end]) / h == pytest.approx(crossing, abs=0.003)
    angle = math.pi * crossing / 2
    closed = math.pi * math.tan(angle) / math.cos(angle) ** 2  # D (pi/h) tan sec^2 at the crossing, in D/h: 5.15
    assert toward.strength(volts) * h / depth == pytest.approx(closed, rel=0.01)
    y = np.abs(toward.v[: end + 1]) / h
    local = 2 * rung.QM * volts * toward.stiff[: end + 1] / toward.omega**2
    assert np.abs(local / (0.5 / np.cos(math.pi * y / 2) ** 2) - 1).max() < 2e-3
    held = [toward.at(f * toward.strength(volts), volts)[1] for f in (0.2, 0.5, 0.8)]
    assert held[0] > held[1] > held[2] > toward.v[end]  # pushed toward the lower board, the electron moves down


def test_the_static_balance_meets_the_rail_pairs_exact_field() -> None:
    # ELE-1's rails at R1's 0.40 and 1.46 mm sampled h/40 apart, against the valley walked on RailPair's own field with
    # no grid: across the board's plane the well holds 3.546 D/h, where the bend rule takes D/h, within 0.5 %
    pair = RailPair(0.40e-3, 1.46e-3)
    h = pair.null_height
    grid = SectionGrid.sampled(
        lambda pu, pv: pair.potential(pu, pv, 1.0), _span(-4 * h, 4 * h, h / 40), _span(h / 40, 5 * h, h / 40), (0, h)
    )
    depth = pt.ponderomotive_ev(pair.escape(1.0).field, F915)  # per volt squared
    gridded = rung.valley(rung.well_patch(grid, F915, h), 0.0).strength(1.0) * h / depth
    w = 2 * math.pi * F915
    t = -np.arange(0.0, 0.6 * h, h / 1000)
    s = np.linspace(-0.3 * h, 0.3 * h, 601)
    e = pair.field(t[:, None] + 0 * s[None, :], h + s[None, :] + 0 * t[:, None], 1.0)
    well = rung.QM * (e[..., 0] ** 2 + e[..., 1] ** 2) / (4 * w * w)
    force = -np.gradient(well.min(axis=1), t)
    exact = float(force.max()) * h / depth
    assert exact == pytest.approx(3.546, rel=5e-3)
    assert gridded == pytest.approx(exact, rel=5e-3)


def test_the_section_model_keeps_mathieus_stability_edge() -> None:
    # in a harmonic well, phi = G ((v - v0)^2 - (u - u0)^2) / 2, the section model's motion along v is Mathieu's: an
    # electron sent from the null at a hundredth of a meV stays within the grid for 40 RF periods at q 0.8, under the
    # edge 0.908, and leaves it at q 1.0 past it; the grid's field is the linear one exactly (second-order differences
    # of a quadratic, bilinear between)
    nodes = np.linspace(-1e-3, 1e-3, 41)
    gradient = 1e6  # V/m^2 per volt
    grid = SectionGrid.sampled(lambda pu, pv: gradient * (pv * pv - pu * pu) / 2, nodes, nodes, (0.0, 0.0))
    omega = np.array([2 * math.pi * F915])
    for q, gone in ((0.8, False), (1.0, True)):
        volts = q / pt.mathieu_q(gradient, F915)
        lost = rung.section_flights(grid, volts, omega, (1e-5,), rung.WELL_PHASES, 40, steps=240)
        assert lost.shape == (1, 2, 1)
        assert bool(lost.all()) is gone
        assert bool(lost.any()) is gone


def test_a_scans_threshold_is_its_edge_and_shrugs_off_a_lone_resonant_loss() -> None:
    # a sharp edge, eight electrons an energy held below it and lost above: the robust threshold and the edges' both
    # give back its geometric mean; one electron lost far below the edge moves the first loss and the edges' threshold
    # there and leaves the robust one where it was; a scan lost from its first energy reads below the scan, one never
    # lost above it, and energies whose electrons the gas stopped count in neither
    energies = tuple(1.025**k for k in range(40))
    held = [8 if k < 20 else 0 for k in range(40)]
    lost = [8 - n for n in held]
    edge = math.sqrt(energies[19] * energies[20])
    assert rung.robust_threshold(energies, held, lost) == pytest.approx(edge, rel=1e-12)
    assert rung.edge_threshold(energies, held, lost) == pytest.approx(edge, rel=1e-12)
    held[5], lost[5] = 7, 1
    assert rung.lowest_lost(energies, lost) == energies[5]
    assert rung.robust_threshold(energies, held, lost) == pytest.approx(edge, rel=1e-12)
    assert rung.edge_threshold(energies, held, lost) < energies[6]
    assert rung.robust_threshold(energies, [0] * 40, [8] * 40) == 0.0
    assert rung.robust_threshold(energies, [8] * 40, [0] * 40) == math.inf
    gassed = [0, 0, 0, *held[3:]], [0, 0, 0, *lost[3:]]
    assert rung.robust_threshold(energies, *gassed) == pytest.approx(edge, rel=1e-12)


def test_earth_along_a_flat_rings_axis_leaves_the_senses_geometric_mean(coarse: Path) -> None:
    # section 6.1's term: B_z adds s e v B_z to the bend's 2W/R, so the two senses' limits W+ and W- move apart and
    # sqrt(W+ W-) stays at the limit without the field (2W/R +- e B_z sqrt(2W/m) = F, a quadratic in sqrt W)
    th = rung.theory(str(coarse / "r1.toml"), str(coarse / "cache"))
    bare = rung.ring_limit(th, 1, np.zeros(3))
    assert bare == pytest.approx(th.static, rel=1e-5)  # the path's curvature by second differences: 1/R to 1e-6
    up = np.array([0.0, 0.0, 50e-6])
    plus, minus = rung.ring_limit(th, 1, up), rung.ring_limit(th, -1, up)
    assert plus > bare > minus
    assert math.sqrt(plus * minus) == pytest.approx(bare, rel=1e-9)
    radius = th.guide.length / (2 * math.pi)
    assert math.sqrt(plus) - math.sqrt(minus) == pytest.approx(50e-6 * radius * math.sqrt(rung.QM / 2), rel=1e-6)


def test_the_bench_bound_on_the_senses_is_earths_most_and_two_steps(coarse: Path) -> None:
    # E2-21 on the bench: |sqrt W- - sqrt W+| within Earth's most over the orientations and two scan steps in sqrt W
    # at the most either sense holds, each sense's threshold up to a step from its truth (PROTOCOL.md, section 5)
    design, cache = str(coarse / "r1.toml"), str(coarse / "cache")
    th = rung.theory(design, cache)
    shot = rung.Shot("E2-21", "R1", design, None, "hold", 1.0, 13.0, phases=rung.R1_PHASES)
    span = rung.earth(design, cache)
    step = math.sqrt(span.each_high) * (math.sqrt(rung.SCAN) - 1)
    bound = rung.predict_hold(shot, th, cache).basis["split_bound_sqrt_eV"]
    assert bound == pytest.approx(span.split + 2 * step, rel=1e-5)
    assert span.mean_low <= span.mean_high <= th.static * (1 + 1e-9)  # the axis's orientation keeps the mean


# ---- the shots and the predictions file ---------------------------------------------------------------------------
def test_the_shots_are_named_once_and_scan_past_their_bands() -> None:
    # fourteen shots in the protocol's order, seeded from the master seed; each scan reaches below and above its band;
    # the predictions file on disk is these shots' (predict writes it from the code); the must-fail shots and their
    # control share one launch, and the control's launch stays under the discount's foot at a q under the edge
    shots = rung.shots()
    idents = [shot.ident for shot in shots]
    assert len(set(idents)) == len(idents) == 14
    assert [shot.seed for shot in shots] == [rung.MASTER + k for k in range(14)]
    assert {shot.test for shot in shots} == {"E2-21", "E2-22", "slot", "must fail"}
    record = json.loads((rung.OUT / "predictions.json").read_text())
    assert [entry["ident"] for entry in record["shots"]] == idents
    for shot, entry in zip(shots, record["shots"], strict=True):
        assert entry["shot"] == plain(shot)
        p = entry["prediction"]
        if shot.kind in ("hold", "well"):
            energies = shot.energies
            assert energies[0] < p["low"] < p["value"]
            assert p["high"] < energies[-1]
            assert energies[1] / energies[0] == pytest.approx(rung.SCAN, rel=1e-6)
        if shot.kind == "hold" and shot.realisation == 0:  # as drawn, the band is a fraction of the static limit
            band = rung.CONTROL_BAND if shot.q is not None else rung.DISCOUNT
            assert p["low"] == pytest.approx(band[0] * p["value"], rel=1e-5)
            assert p["high"] == pytest.approx(band[1] * p["value"], rel=1e-5)
        if shot.kind == "fail":
            assert p["basis"]["q_null"] > rung.QCAP
        if shot.kind == "control":
            assert p["basis"]["q_null"] < rung.QCAP
            assert p["basis"]["top_over_static"] < rung.DISCOUNT[0]
    launches = [shot for shot in shots if shot.kind in ("fail", "control")]
    assert {shot.launch for shot in launches} == {rung.MASTER}
    assert len({rung._phases(shot) for shot in launches}) == 1
    assert len({shot.energies for shot in launches}) == 1


def test_a_holds_steps_cover_its_laps_at_its_slowest_energy() -> None:
    # DR-7's 20 laps: a hold runs until its slowest electron, at its launch speed, has gone 20.4 laps, and no longer
    record = json.loads((rung.OUT / "predictions.json").read_text())
    holds = [shot for shot in rung.shots() if shot.kind == "hold"]
    assert len(holds) == 6
    for shot in holds:
        ring = record["apparatus"]["designs"][shot.design]["ring"]
        length, dt = ring["length_mm"] * 1e-3, 1 / (ring["f_MHz"] * 1e6 * rung.STEPS)
        steps = rung.run_steps(shot, length, dt)
        slowest = pt.v_of(min(shot.energies))
        assert steps * dt * slowest >= rung.LAPS * 1.02 * length
        assert (steps - 1) * dt * slowest < rung.LAPS * 1.02 * length
    beta = next(shot for shot in rung.shots() if shot.kind == "beta")
    assert rung.run_steps(beta, 0.2, 1e-12) == rung.BETA_PERIODS * rung.STEPS


# ---- the tracking's side ------------------------------------------------------------------------------------------
def test_measure_takes_one_tiny_shot_of_each_kind(coarse: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # each kind of shot on coarse copies of the designs, for a few RF periods or a few thousandths of a lap: the counts
    # add up to what was launched, the seeds and phases are recorded, the bench's own guide is used for realisation
    # None and another opened beside it otherwise, a hold runs until its slowest electron has made its laps, and each
    # measurement comes back from its plain record
    cache = str(coarse / "cache")
    r1, slot, knot = (str(coarse / f"{name}.toml") for name in ("r1", "slot", "t915s"))
    devices: dict[str, Any] = {
        "guide": TwinGuide(r1, 1, cache=cache),
        "detector": TwinDetector(collector=1e-12),
        "vacuum": TwinVacuum(1e-6),
    }
    tiny = [
        rung.Shot("E2-21", "R1", r1, None, "hold", 1.0, 1.03, laps=0.004, phases=1, seed=1),
        rung.Shot("E2-22", "T915-S", knot, 0, "hold", 3.5, 3.6, laps=0.003, phases=1, seed=2),
        rung.Shot("E2-21", "R1", r1, 0, "well", 0.02, 0.0205, periods=2, seed=3),
        rung.Shot("slot", "slot", slot, 0, "beta", phases=1, periods=4, seed=4),
        rung.Shot("must fail", "R1", r1, None, "fail", 0.5, 1.0, ratio=2.0, phases=1, periods=2, depth=1.0, seed=5),
        rung.Shot("must fail", "R1", r1, None, "control", 0.5, 1.0, ratio=2.0, phases=1, periods=2, seed=6),
        rung.Shot("E2-21", "R1", r1, 0, "well", 0.02, 0.0205, periods=2, q=rung.CONTROL_Q, seed=7),
    ]
    for shot in tiny:
        m = rung.measure(shot, devices)
        assert plain(rung.load(plain(m))) == plain(m)
        assert m.seed == shot.seed
        assert m.guide["realisation"] == (1 if shot.realisation is None else 0)
        if isinstance(m, rung.Scan):
            rows = 1 if shot.kind == "well" else 2
            assert len(m.senses) == rows
            assert m.launched() == rows * len(shot.energies) * len(m.phases)
            assert m.steps > 0
            if shot.kind == "hold":
                length = rung._guide(shot, devices).length
                assert m.steps * m.dt * pt.v_of(min(shot.energies)) >= shot.laps * 1.02 * length
                assert not any(sum(row) for row in m.short)
        else:
            assert len(m.beta) == 4
            assert all(0 < b < 1 for b in m.beta)
    well = rung.measure(tiny[2], devices)
    assert isinstance(well, rung.Scan)
    assert well.phases == (math.pi / 2, 3 * math.pi / 2)
    fail, control = rung.measure(tiny[4], devices), rung.measure(tiny[5], devices)
    assert isinstance(fail, rung.Scan)
    assert isinstance(control, rung.Scan)
    assert fail.q == pytest.approx(devices["guide"].solved.mathieu_q(fail.volts) * (F915 / fail.f) ** 2, rel=1e-6)
    assert control.volts == pytest.approx(50.0, rel=1e-12)
    assert fail.phases == control.phases  # one launch for the must-fail shots and their control
    calm = rung.measure(tiny[6], devices)
    assert isinstance(calm, rung.Scan)
    assert calm.q == pytest.approx(rung.CONTROL_Q, rel=2e-3)
    # a hold whose run ends before its electrons have made their laps holds none of them: they are short, not held
    monkeypatch.setattr(rung, "run_steps", lambda shot, length, dt: 20)
    cut = rung.measure(tiny[0], devices)
    assert isinstance(cut, rung.Scan)
    assert cut.held_total() == 0
    assert sum(sum(row) for row in cut.short) + cut.lost_total() + sum(sum(row) for row in cut.stopped) == 4
    other = rung.Shot("E2-21", "T915-S", knot, None, "hold", 1.0, 1.03, laps=0.004, phases=1)
    with pytest.raises(Refusal, match="the bench's guide is R1"):
        rung.measure(other, devices)


def test_a_run_waits_for_predictions_md_to_be_the_codes(monkeypatch: pytest.MonkeyPatch) -> None:
    # predictions.md is registered beside predictions.json: `run` refuses while it is not what the code predicts now
    monkeypatch.setattr(rung, "predictions_markdown", lambda experiment: "not the file on disk\n")
    assert rung.main(["run", "--bench", "twin_r1"]) == 2


def test_a_launch_starts_on_its_slow_motion_in_step_with_the_drive() -> None:
    # in a uniform field E = -V(t) along u, V(t) = Re(V e^(i w t)) per volt per metre, an electron launched by launch()
    # with the micromotion of its place and phase (its velocity half a step early, as Boris's leapfrog takes it) has no
    # slow motion: averaged over its tenth RF period it is where it was put, within 20 nm at 240 steps a period (the
    # step's second order; without the half step 0.4 um, the first), while one put there at rest drifts off at the
    # micromotion's velocity, a third of a millimetre
    nodes = np.linspace(-2e-3, 2e-3, 21)
    grid = SectionGrid(nodes, nodes, np.repeat(nodes[:, None], len(nodes), axis=1), (0.0, 0.0))
    omega = 2 * math.pi * F915
    swept = SectionField.straight(grid, 0.1, pt.InPhase(1000.0, omega, 0.4))
    field = TwinField(swept, None, None)
    phasor = 1000.0 * complex(math.cos(0.4), -math.sin(0.4))  # InPhase: V cos(w t - 0.4)
    start = rung.Station(swept.points[0], swept.tangents[0], swept.first[0], swept.second[0], (0.0, 0.0), phasor)
    t0, steps = 0.3 / F915, 240
    matched = rung.launch(field, start, np.zeros((1, 2)), np.zeros((1, 3)), t0, omega, 1 / (F915 * steps))
    drift = []
    for x0, u0 in (matched, (start.point[None, :], np.zeros((1, 3)))):
        seen: list[pt.Arr] = []

        def keep(step: int, t: float, x: pt.Arr, u: pt.Arr, alive: pt.BoolArr) -> None:
            seen.append(x[0].copy())  # noqa: B023

        pt.integrate(x0, u0, np.array([pt.QM_E]), 1 / (F915 * steps), 10 * steps, rung._shifted(field, t0), [keep])
        drift.append(float(np.linalg.norm(np.mean(seen[-steps:], axis=0) - start.point)))
    assert drift[0] < 2e-8
    assert drift[1] > 1e4 * drift[0]


# ---- the verdicts -------------------------------------------------------------------------------------------------
def _scan(kind: str, energies: tuple[float, ...], lost_from: tuple[float, ...], q: float = 0.5) -> rung.Scan:
    """A made-up scan: in each sense, one electron at each energy, lost from the given energy on."""
    senses = (0,) if kind == "well" else (1, -1)
    lost = tuple(tuple(int(e >= edge) for e in energies) for edge in lost_from)
    held = tuple(tuple(1 - n for n in row) for row in lost)
    zero = tuple(tuple(0 for _ in energies) for _ in senses)
    first = tuple(tuple(math.inf for _ in energies) for _ in senses)
    return rung.Scan(kind, {}, 0, energies, senses, (0.0,), held, lost, zero, zero, first, 1, 1e-12, F915, 50.0, q)


def _prediction(value: float, low: float, high: float, **more: float) -> rung.Prediction:
    ratio = more.pop("ratio", None)
    lo, hi = (None, None) if ratio is None else (ratio * 0.85, ratio * 1.15)
    return rung.Prediction("made up", value, low, high, "eV", dict(more), ratio, lo, hi)


def test_the_verdicts_read_the_scans_by_the_protocols_criteria() -> None:
    shots = list(rung.shots())
    matches: dict[str, Callable[[rung.Shot], bool]] = {
        "well": lambda s: s.kind == "well" and s.test == "E2-21" and s.realisation == 0 and s.q is None,
        "well control": lambda s: s.kind == "well" and s.q is not None,
        "drawn": lambda s: s.kind == "hold" and s.test == "E2-21" and s.realisation == 0 and s.q is None,
        "bench": lambda s: s.kind == "hold" and s.test == "E2-21" and s.realisation is None,
        "hold control": lambda s: s.kind == "hold" and s.q is not None,
        "knot": lambda s: s.kind == "hold" and s.label == "T915-S",
        "flat": lambda s: s.kind == "hold" and s.label == "flat ring 318",
        "beta": lambda s: s.kind == "beta",
        "slot well": lambda s: s.kind == "well" and s.test == "slot",
        "fail": lambda s: s.kind == "fail",
        "control": lambda s: s.kind == "control",
    }
    by_name = {name: next(shot for shot in shots if match(shot)) for name, match in matches.items()}
    grid = tuple(1.0 * rung.SCAN**k for k in range(160))  # 1 to 50 eV
    drawn, bench = by_name["drawn"], by_name["bench"]
    for plus, minus, apart, split in ((7.0, 7.1, True, True), (7.0, 8.0, False, True), (9.5, 6.0, False, False)):
        rows: list[rung.Row] = [
            (drawn, _prediction(9.0, 5.0, 9.0), _scan("hold", grid, (plus, minus))),
            (bench, _prediction(9.0, 5.0, 9.0, split_bound_sqrt_eV=0.6), _scan("hold", grid, (plus, minus))),
            (by_name["hold control"], _prediction(3.3, 2.8, 3.45), _scan("hold", grid, (3.2, 3.2))),
        ]
        checks = {(c.criterion, c.ident): c for c in rung.TwinRung1().judge(rows)}
        assert checks[("E2-21 hold", drawn.ident)].passed  # the geometric means 7.0, 7.5 and 7.5 eV: in the band
        assert checks[("E2-21 hold", bench.ident)].passed
        assert checks[("E2-21 hold control", by_name["hold control"].ident)].passed
        assert checks[("E2-21 senses", drawn.ident)].passed is apart  # within two steps of 2.5 %, or not
        assert checks[("E2-21 senses", bench.ident)].passed is split  # sqrt 9.5 - sqrt 6 = 0.63, past 0.6
    knot, flat = by_name["knot"], by_name["flat"]
    rows = [
        (knot, _prediction(11.6, 7.1, 11.1), _scan("hold", grid, (9.0, 9.0))),
        (flat, _prediction(32.4, 19.8, 30.8, ratio=2.79), _scan("hold", grid, (36.0, 36.0))),
    ]
    checks = {(c.criterion, c.ident): c for c in rung.TwinRung1().judge(rows)}
    assert checks[("E2-22 static limit", knot.ident)].passed
    assert not checks[("E2-22 static limit", flat.ident)].passed  # 36 eV is above 0.95 of the static limit
    assert not checks[("E2-22 price of the knot", flat.ident)].passed  # 36/9 = 4.0, above 2.79 x 1.15
    # the wells: E2-21's first loss against the section model's band, the slot guide's edges against its depth
    sweep = tuple(0.0045 * rung.SCAN**k for k in range(140))
    calm = tuple(0.006 * rung.SCAN**k for k in range(86))
    rows = [
        (by_name["well"], _prediction(0.0118, 0.0111, 0.0125), _scan("well", sweep, (0.012,))),
        (by_name["well control"], _prediction(0.0314, 0.0295, 0.0334), _scan("well", calm, (0.027,))),
        (by_name["slot well"], _prediction(0.229, 0.183, 0.287), _scan("well", grid[:30], (1.2,))),
    ]
    checks = {(c.criterion, c.ident): c for c in rung.TwinRung1().judge(rows)}
    assert checks[("E2-21 well", by_name["well"].ident)].passed
    assert not checks[("E2-21 well control", by_name["well control"].ident)].passed  # lost from 0.027, under 0.0295
    assert not checks[("slot depth", by_name["slot well"].ident)].passed  # the edges at 1.19 eV, above the band
    # the cases that must fail: none held, and only beside a control whose launch, at the design's volts, held
    fail, control = by_name["fail"], by_name["control"]
    energies = (0.5, 1.0, 2.0, 4.0)
    nothing = _scan("fail", energies, (0.0, 0.0), q=1.64)
    for control_scan, ok in (
        (_scan("control", energies, (8.0, 8.0)), True),
        (_scan("control", energies, (3.0, 8.0)), False),
    ):
        rows = [
            (fail, _prediction(0.0, 0.0, 0.0, q_null=1.64), nothing),
            (control, _prediction(0.0, 0.0, 0.0, q_null=0.51), control_scan),
        ]
        checks = {(c.criterion, c.ident): c for c in rung.TwinRung1().judge(rows)}
        assert checks[("must fail control", control.ident)].passed is ok
        assert checks[("must fail", fail.ident)].passed is ok
    rows = [(fail, _prediction(0.0, 0.0, 0.0, q_null=1.64), nothing)]
    checks = {(c.criterion, c.ident): c for c in rung.TwinRung1().judge(rows)}
    assert not checks[("must fail", fail.ident)].passed  # no control on file: the failure says nothing alone
    rows = [
        (fail, _prediction(0.0, 0.0, 0.0, q_null=1.64), _scan("fail", energies, (0.0, 2.0), q=1.64)),
        (control, _prediction(0.0, 0.0, 0.0, q_null=0.51), _scan("control", energies, (8.0, 8.0))),
    ]
    checks = {(c.criterion, c.ident): c for c in rung.TwinRung1().judge(rows)}
    assert not checks[("must fail", fail.ident)].passed  # two electrons held against the wave: the case did not fail
    beta = by_name["beta"]
    exact = rung.mathieu_beta(0.4421)
    for found, held, ok in (
        ((1.0005, 0.9995, 1.0, 1.001), 4, True),
        ((1.0, 1.0, 1.0, 1.003), 4, False),
        ((1.0,) * 4, 3, False),
    ):
        beat = rung.Beat(
            "beta", {}, 0, ((0.0, 0.0),) * 4, tuple(exact * b for b in found), held, 1, 1e-12, F915, 40.0, 0.44
        )
        rows = [(beta, rung.Prediction("made up", exact, exact * 0.998, exact * 1.002, "beta", {}), beat)]
        checks = {(c.criterion, c.ident): c for c in rung.TwinRung1().judge(rows)}
        assert checks[("slot secular frequency", beta.ident)].passed is ok  # within 0.2 %, and every electron held
