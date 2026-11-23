"""experiments/twin_rung2.py, experiment 2's stage 2a (docs/KNOTTED_GUIDE.md section 10): its theory held to stage 1's
and to closed forms (the static balance with a uniform field against twin_rung1.ring_limit on a ring, the bound with
every bend counted across the boards' plane against the static limit, the lobes against stage 1's, the worst
orientation against a sample of the sphere); the exact motion in the section keeping a matched electron round a ring and
losing one past its fold, and losing on a knot what the tracker loses, where it loses it; its shots against the
predictions file, paired shots sharing their launch; measure() on one tiny shot of each kind, on coarse copies of the
designs, in one process and in the pool alike, so that the run cannot fail on a bug; its verdicts on measurements made
up to pass and to fail; and the run's time under eight hours. No test tracks long enough, or near enough a registered
threshold, to bear on a verdict: the protocol's section 6 lists what they run."""

from __future__ import annotations

import json
import math
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import twin_rung1 as rung
import twin_rung2 as r2

from array16 import particles as pt
from array16.guide import axisym
from array16.guide import emission as em
from array16.guide.sections import COARSE, section_b
from array16.lab import Refusal, plain
from array16.twin import Losses, TwinDetector, TwinGuide, TwinVacuum

ROOT = Path(__file__).resolve().parents[1]
NAMES = {
    "t915s": r2.T915S,
    "flat": r2.FLAT,
    "stiff": r2.STIFF,
    "helix": r2.HELIX,
    "r1": r2.R1,
    "ring_in": r2.RING_IN,
}


@pytest.fixture(scope="module")
def coarse(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The designs on the coarse grid, beside the test, with one cache for the module."""
    where = tmp_path_factory.mktemp("rung2")
    for name, source in NAMES.items():
        text = (ROOT / source).read_text()
        assert 'grid = "default"' in text
        (where / f"{name}.toml").write_text(text.replace('grid = "default"', 'grid = "coarse"'))
    return where


def _design(coarse: Path, name: str) -> str:
    return str(coarse / f"{name}.toml")


def _cache(coarse: Path) -> str:
    return str(coarse / "cache")


# ---- the theory against stage 1's and closed forms ----------------------------------------------------------------
def test_a_rings_field_limit_is_stage_1s_ring_limit(coarse: Path) -> None:
    # on a flat ring the static balance with a uniform field is twin_rung1.ring_limit, whatever the field's
    # orientation; with no field it is the static limit, and a field along the ring's axis parts the senses' sqrt W by
    # B R sqrt(e/2m) and keeps their geometric mean (stage 1's Proposition)
    design, cache = _design(coarse, "r1"), _cache(coarse)
    th = r2.theory(design, 0, None, cache)
    frame = r2.frame_of(design, 0, cache)
    assert r2.field_limit(th, 1, np.zeros(3), frame) == pytest.approx(th.static, rel=1e-6)
    rng = np.random.default_rng(1)
    for _ in range(3):
        b = 50e-6 * rng.normal(size=3) / math.sqrt(3)
        for sense in (1, -1):
            assert r2.field_limit(th, sense, b, frame) == pytest.approx(rung.ring_limit(th, sense, b), rel=1e-6)
    up = np.array([0.0, 0.0, 50e-6])
    plus, minus = r2.field_limit(th, 1, up, frame), r2.field_limit(th, -1, up, frame)
    assert math.sqrt(plus * minus) == pytest.approx(th.static, rel=1e-6)


def test_the_across_bound_is_the_static_limit_bent_in_the_plane_and_under_it_out_of_it(coarse: Path) -> None:
    # H1b's bound counts every bend across the boards' plane: on the flat ring, bent in that plane, it is the static
    # limit; on the knot it lies under it, stage 1's 9.80 against 11.64 eV, 0.842 of it; on the stiff ring, bent wholly
    # toward the boards, it is the flat ring's static limit, the same circle and drive
    cache = _cache(coarse)
    flat = r2.theory(_design(coarse, "flat"), 0, None, cache)
    knot = r2.theory(_design(coarse, "t915s"), 0, None, cache)
    stiff = r2.theory(_design(coarse, "stiff"), 0, None, cache)
    ku = {name: r2.frame_of(_design(coarse, name), 0, cache)[3] for name in ("flat", "t915s", "stiff")}
    assert r2.across_bound(flat, ku["flat"]) == pytest.approx(flat.static, rel=1e-6)
    assert r2.across_bound(knot, ku["t915s"]) / knot.static == pytest.approx(0.842, abs=0.01)
    assert r2.across_bound(stiff, ku["stiff"]) == pytest.approx(flat.static, rel=1e-6)
    assert stiff.static > 4 * flat.static


def test_the_lobes_are_stage_1s_and_a_loss_is_placed_downstream(coarse: Path) -> None:
    # stage 1's lobes of T915-S (REPORT.md section 6, H3): between the maxima of the curvature in the boards' plane
    # that pair round each third of the knot; none on a circle; a loss in a lobe is 0 past it, one after it in its
    # sense of travel the distance from the lobe's end, and one before it the distance on to the next lobe's
    cache = _cache(coarse)
    knot = r2.lobes(r2.frame_of(_design(coarse, "t915s"), 0, cache)[3])
    stage1 = ((0.1385, 0.1948), (0.4720, 0.5280), (0.8052, 0.8615))
    assert len(knot) == 3
    for (a, b), (c, d) in zip(knot, stage1, strict=True):
        assert a == pytest.approx(c, abs=1e-3)
        assert b == pytest.approx(d, abs=1e-3)
    assert len(r2.lobes(r2.frame_of(_design(coarse, "helix"), 0, cache)[3])) == 3
    assert r2.lobes(r2.frame_of(_design(coarse, "flat"), 0, cache)[3]) == ()
    assert r2.past(0.15, 1, knot) == 0.0
    assert r2.past(0.20, 1, knot) == pytest.approx(0.20 - knot[0][1], abs=1e-12)
    assert r2.past(0.13, -1, knot) == pytest.approx(knot[0][0] - 0.13, abs=1e-12)
    assert r2.past(0.13, 1, knot) == pytest.approx(0.13 + 1 - knot[2][1], abs=1e-12)
    assert r2.past(0.5, 1, ()) == math.inf


def test_the_worst_orientation_parts_t915s_senses_most() -> None:
    # E2-15's field: the registered orientation parts T915-S's senses at least as much as any of a hundred orientations
    # spread over the sphere, by the static balance with the field, and by about 1.9 % at 50 uT; ten times the field
    # parts them by about ten times as much
    th = r2.theory(r2.T915S, 0, None, None)
    frame = r2.frame_of(r2.T915S, 0, None)
    way = np.array(r2.WORST)

    def split(direction: pt.Arr, size: float = r2.EARTH) -> float:
        b = size * direction / np.linalg.norm(direction)
        return math.log(r2.field_limit(th, -1, b, frame) / r2.field_limit(th, 1, b, frame))

    k = np.arange(100) + 0.5
    polar, around = np.arccos(1 - 2 * k / 100), math.pi * (1 + math.sqrt(5)) * k
    sample = np.stack([np.sin(polar) * np.cos(around), np.sin(polar) * np.sin(around), np.cos(polar)], 1)
    worst = abs(split(way))
    assert worst >= max(abs(split(w)) for w in sample) - 1e-9
    assert worst == pytest.approx(0.019, abs=0.003)
    assert abs(split(way, r2.STRONG * r2.EARTH)) > 5 * worst


def test_the_exact_motion_keeps_a_matched_electron_on_a_ring_and_loses_one_past_its_fold(coarse: Path) -> None:
    # round a ring the bend's force is constant: an electron launched on its static equilibrium with its micromotion
    # matched stays on it, at half the static limit, for a third of a lap; one at 1.3 times the limit has no
    # equilibrium and is lost before a quarter lap, where the ring's drive is weakest
    design, cache = _design(coarse, "flat"), _cache(coarse)
    th = r2.theory(design, 0, None, cache)
    spec = r2.MotionSpec(design, 0, None, (0.0, 0.0, 0.0), 0, (0.5 * th.static, 1.3 * th.static), (1,), (0.4,), 0.3)
    m = r2.exact_motion(spec, cache)
    fate, reached = np.array(m.table.fate), np.array(m.table.reached)
    assert fate[0, 0].tolist() == [0, 1]
    assert reached[0, 0, 0] >= 0.3
    assert reached[0, 0, 1] < 0.25
    assert math.isnan(m.place[0][0][0])
    assert 0 < m.place[0][0][1] < 0.25


def test_the_exact_motion_loses_on_a_knot_what_the_tracker_loses_where_it_loses_it(coarse: Path) -> None:
    # T915-S's first lobe from its start: launched at 0.7 of the static limit both the exact motion and the tracker
    # hold the electron past it, and at 1.3 both lose it there, in the lobe or within 0.07 lap past it, at two phases
    design, cache = _design(coarse, "t915s"), _cache(coarse)
    th = r2.theory(design, 0, None, cache)
    energies = (0.7 * th.static, 1.3 * th.static)
    shot = r2.Shot("H2", "T915-S", design, 0, "hold", energies[0], energies[1] * 1.001, ratio=1.3 / 0.7, phases=2,
                   laps=0.3, seed=5)  # fmt: skip
    assert shot.energies == pytest.approx(energies, rel=1e-6)
    lobe = r2.lobes(r2.frame_of(design, 0, cache)[3])
    m = r2.exact_motion(r2.motion_spec(shot, cache), cache)
    devices = _devices(coarse, "t915s")
    tracked = r2.measure([shot], devices, 1, cache)[shot.ident]
    assert isinstance(tracked, r2.Scan)
    for table in (m.table, tracked.table):
        fate = np.array(table.fate)
        assert (fate[:, :, 0] == 0).all()
        assert (fate[:, :, 1] == 1).all()
    for k in range(2):
        for row, sense in enumerate((1, -1)):
            assert r2.past(m.place[k][row][1], sense, lobe) <= r2.LOBE_REACH
    entries = [r for losses in tracked.losses for r in range(len(losses))]
    assert len(entries) == 4
    for losses in tracked.losses:
        for r in range(len(losses)):
            sense = (1, -1)[int(losses.index[r]) // 2]
            assert r2.past(float(losses.s[r]) / losses.length % 1.0, sense, lobe) <= r2.LOBE_REACH


# ---- the shots and the predictions file ---------------------------------------------------------------------------
def test_the_shots_are_named_once_and_the_pairs_share_their_launch() -> None:
    # thirty shots in the protocol's order, seeded from the master seed; the bend's pair and the must-fail pair each
    # share one launch (starts, phases and gas); the predictions file on disk is these shots' (predict writes it from
    # the code), each hold's scan reaching past its exact motion's limit either way
    shots = r2.shots()
    idents = [shot.ident for shot in shots]
    assert len(set(idents)) == len(idents) == 30
    assert [shot.seed for shot in shots] == [r2.MASTER + k for k in range(30)]
    assert {shot.test for shot in shots} == set(r2.TESTS)
    for name in ("E2-2", "must fail"):
        pair = [shot for shot in shots if shot.test == name and shot.kind in ("hold", "fail", "control")]
        assert len(pair) == 2
        assert len({shot.launch for shot in pair}) == 1
        assert len({r2._phases(shot) for shot in pair}) == 1
    path = r2.OUT / "predictions.json"
    if not path.exists():
        pytest.skip("the predictions file is written by `predict`")
    record = json.loads(path.read_text())
    assert [entry["ident"] for entry in record["shots"]] == idents
    for shot, entry in zip(shots, record["shots"], strict=True):
        assert entry["shot"] == plain(shot)
        p = entry["prediction"]
        if shot.kind == "hold":
            assert shot.energies[0] < p["low"] <= p["value"] <= p["high"] < shot.energies[-1]


def test_the_run_takes_under_eight_hours() -> None:
    # the run's time on this machine's ten workers, by each job's measured cost: under eight hours, H2's batch the
    # longest; and a batch is one test, in the protocol's order
    time = r2.estimate()
    assert time["wall_s"] < 8 * 3600
    assert max(time["batches_s"], key=time["batches_s"].get) == "H2"
    assert [batch[0].test for batch in r2.batches(r2.shots())] == list(r2.TESTS)


def test_a_run_waits_for_predictions_md_to_be_the_codes(monkeypatch: pytest.MonkeyPatch) -> None:
    # predictions.md is registered beside predictions.json: `run` refuses while it is not what the code predicts now
    monkeypatch.setattr(r2, "predictions_markdown", lambda experiment: "not the file on disk\n")
    assert r2.main(["run", "--bench", "twin_r2"]) == 2


# ---- the tracking's side ------------------------------------------------------------------------------------------
def _devices(coarse: Path, bench: str) -> dict[str, Any]:
    return {
        "guide": TwinGuide(_design(coarse, bench), 1, cache=_cache(coarse)),
        "detector": TwinDetector(collector=1e-12, losses=True),
        "vacuum": TwinVacuum(1e-6),
    }


def _tiny_mpc(monkeypatch: pytest.MonkeyPatch) -> None:
    """E2-13's readings cut to a coarse case and a few electrons for a few RF periods."""
    case = em.Case("section B at 50.3 V, coarse", section_b(grid=COARSE), 50.3)
    monkeypatch.setattr(r2, "mpc_cases", lambda name: (case,))
    monkeypatch.setattr(r2, "MPC_SEEDS", (4, 5))
    monkeypatch.setattr(r2, "MPC_SETTINGS", {**r2.MPC_SETTINGS, "electrons": 40, "periods": 4, "settle": 1,
                                             "count": 3, "generation_settle": 1, "steps": 40})  # fmt: skip


def test_measure_takes_one_tiny_shot_of_each_kind(coarse: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # each kind of shot on coarse copies of the designs, for a few RF periods or a few hundredths of a lap: the fates
    # add up to what was launched, the seeds and phases are recorded, each lost electron has its entry, the bench's own
    # guide is used for realisation None and another opened beside it otherwise, a hold with the bend's shift and one
    # with a field run in them, and each measurement comes back from its plain record
    _tiny_mpc(monkeypatch)
    cache = _cache(coarse)
    t915s, r1, ring = (_design(coarse, name) for name in ("t915s", "r1", "ring_in"))
    devices = _devices(coarse, "t915s")
    tiny = [
        r2.Shot("E2-15", "T915-S", t915s, None, "hold", 8.0, 8.2, laps=0.01, phases=1, seed=1),
        r2.Shot("E2-2", "T915-S", t915s, 0, "hold", 8.0, 8.2, laps=0.01, phases=2, bend=True, seed=2),
        r2.Shot("E2-15", "T915-S", t915s, 0, "hold", 8.0, 8.2, laps=0.01, phases=1, earth=(0.0, 0.0, 5e-4), seed=3),
        r2.Shot("H4", "R1", r1, 1, "well", 0.02, 0.0205, periods=2, earth="none", seed=4),
        r2.Shot("must fail", "T915-S", t915s, 0, "fail", 0.5, 1.0, ratio=2.0, phases=1, periods=2, q=1.2, seed=5),
        r2.Shot("must fail", "T915-S", t915s, 0, "control", 0.5, 1.0, ratio=2.0, phases=1, periods=2, seed=6),
        r2.Shot("E2-2", "ring", ring, 0, "shift", seed=7),
        r2.Shot("E2-13", "tiny", "t915s", 0, "multipactor", seed=8),
    ]
    found = r2.measure(tiny, devices, 1, cache)
    for shot in tiny:
        m = found[shot.ident]
        assert plain(r2.load(json.loads(json.dumps(plain(m))))) == plain(m)
        assert m.seed == shot.seed
        if isinstance(m, r2.Scan):
            table = m.table
            assert table.phases == r2._phases(shot)
            assert np.array(table.fate).shape == (len(table.phases), len(shot.senses), len(shot.energies))
            assert m.guide["realisation"] == (1 if shot.realisation is None else shot.realisation)
            lost = int((np.array(table.fate) == 1).sum())
            assert sum(int(np.isin(losses.how, ("struck", "left")).sum()) for losses in m.losses) >= lost
            assert m.bend is shot.bend
    hold = found[tiny[0].ident]
    assert isinstance(hold, r2.Scan)
    assert hold.earth == tuple(float(x) for x in devices["guide"].drawn.earth)
    field = found[tiny[2].ident]
    assert isinstance(field, r2.Scan)
    assert field.earth == (0.0, 0.0, 5e-4)
    shift = found[tiny[6].ident]
    assert isinstance(shift, r2.Shift)
    assert shift.coefficient == pytest.approx(axisym.IN_PLANE, rel=1e-3)
    mpc = found[tiny[7].ident]
    assert isinstance(mpc, r2.Mpc)
    assert len(mpc.readings) == len(r2.SURFACES["t915s"]) * len(r2.LAMINATES)
    other = r2.Shot("E2-15", "R1", r1, None, "hold", 1.0, 1.03, laps=0.004, phases=1)
    with pytest.raises(Refusal, match="the bench's guide is T915-S"):
        r2.measure([other], devices, 1, cache)


def test_the_pool_gives_what_one_process_gives(coarse: Path) -> None:
    # a scan's RF phases tracked in two worker processes, each opening the devices again from their options, give the
    # fates, laps and entries one process gives
    cache = _cache(coarse)
    shot = r2.Shot("H4", "R1", _design(coarse, "r1"), 1, "well", 0.02, 0.0205, periods=2, seed=9)
    devices = _devices(coarse, "t915s")
    one = r2.measure([shot], devices, 1, cache)[shot.ident]
    two = r2.measure([shot], devices, 2, cache)[shot.ident]
    assert plain(one) == plain(two)


def test_a_bent_launch_starts_on_the_moved_null(coarse: Path) -> None:
    # with the bend's shift on the swept null stands toward the bend's centre, and a launch at no offset from it puts
    # the electron where the field vanishes, at the moved null, not on the path
    guide = TwinGuide(_design(coarse, "t915s"), 0, cache=_cache(coarse))
    shot = r2.Shot("E2-2", "T915-S", guide.design.source, 0, "hold", 8.0, 8.2, bend=True)
    _, field = r2.shot_field(guide, shot)
    swept = field.swept
    start = r2.station_at(guide, field, 0)
    du, dv = swept.null_shift(swept.locate(start.point[None, :]))
    assert math.hypot(float(du[0]), float(dv[0])) > 1e-6
    x0, _ = r2.launch(field, start, np.zeros((1, 2)), (float(du[0]), float(dv[0])), np.zeros((1, 3)), 0.0, 1.0, 1e-12)
    e_null, _ = field(x0 + 0 * start.first, 0.0)
    e_path, _ = field(start.point[None, :], 0.0)
    assert e_null is not None
    assert e_path is not None
    assert np.linalg.norm(e_null) < 1e-3 * np.linalg.norm(e_path)


# ---- the verdicts -------------------------------------------------------------------------------------------------
def _table(energies: tuple[float, ...], edges: tuple[float, ...], run: float = 20.0, phases: int = 8) -> r2.Table:
    """A made-up hold: in each sense, `phases` electrons at each energy, lost from the given energy on."""
    fate = tuple(tuple(tuple(int(e >= edge) for e in energies) for edge in edges) for _ in range(phases))
    reached = tuple(tuple(tuple(run if e < edge else 0.5 for e in energies) for edge in edges) for _ in range(phases))
    return r2.Table(energies, (1, -1), tuple(range(phases)), fate, reached, run)


def _hold(table: r2.Table, earth: tuple[float, float, float] = (0.0, 0.0, 0.0)) -> r2.Scan:
    return r2.Scan("hold", {}, 0, table, (), 1, 1e-12, 915e6, 50.0, 0.49, earth, False, 0)


def test_the_verdicts_read_the_measurements_by_the_protocols_criteria() -> None:
    shots = {shot.ident: shot for shot in r2.shots()}
    grid = tuple(1.0 * r2.FINE**k for k in range(400))
    ellipse = next(s for s in shots.values() if s.design == r2.ELLIPSE)
    basis = {"static_eV": 11.9, "H1_low_eV": 11.6, "H1_high_eV": 12.26, "H1b_low_eV": 11.74, "H1b_high_eV": 12.33}
    p = r2.Prediction("made up", 11.4, 11.1, 11.8, "eV", basis)
    rows: list[r2.Row] = [(ellipse, p, _hold(_table(grid, (11.45, 11.45))))]
    checks = {c.criterion: c for c in r2.TwinRung2().judge(rows)}
    assert checks["H1 exact motion"].passed
    assert checks["H1a: the changing bend"].passed
    assert not checks["H1: the plane"].passed  # 11.4 eV is 0.958 of the static limit: under H1's 0.975
    assert not checks["H1b: the stiff side"].passed
    # E2-15: at 50 uT the senses within the static split and two steps; at ten times the field half to twice its split
    weak, strong = (s for s in shots.values() if s.test == "E2-15" and isinstance(s.earth, tuple))
    split = {"static_eV": 11.6, "field_split": -0.02}
    rows = [
        (weak, r2.Prediction("made up", 9.0, 8.8, 9.3, "eV", split), _hold(_table(grid, (9.3, 9.0)))),
        (strong, r2.Prediction("made up", 7.0, 6.8, 7.2, "eV", split | {"field_split": -0.1}),
         _hold(_table(grid, (7.4, 6.6)))),
    ]  # fmt: skip
    pairs = {(c.criterion, c.ident): c for c in r2.TwinRung2().judge(rows)}
    assert pairs[("E2-15 Earth's split", weak.ident)].passed  # ln(9.0/9.3) = -0.033, within 0.02 and two steps
    assert pairs[("E2-15 the field's split", strong.ident)].passed  # -0.114 against -0.1
    # E2-13: decays, the bound and the control
    t915s = next(s for s in shots.values() if s.design == "t915s")
    control = next(s for s in shots.values() if s.design == "overdriven")
    decays = r2.MpcReading("c", "bare", "absorbs", 50.3, -0.17, 0.01, -0.17, 0.80, 0.01, 16.0, 0.05, False)
    grows = replace(decays, growth=0.5, secondaries=1.3)
    rows = [
        (
            t915s,
            r2.Prediction("made up", -0.1, -math.inf, 0.0, "Np", {}),
            r2.Mpc("multipactor", {}, 0, (4,), (decays,)),
        ),
        (
            control,
            r2.Prediction("made up", 1.5, 0.0, math.inf, "Np", {}),
            r2.Mpc("multipactor", {}, 0, (4,), (decays,)),
        ),
    ]
    checks = {c.criterion: c for c in r2.TwinRung2().judge(rows)}
    assert checks["E2-13 decays"].passed
    assert not checks["E2-13 the control grows"].passed
    rows[1] = (rows[1][0], rows[1][1], r2.Mpc("multipactor", {}, 0, (4,), (grows,)))
    assert {c.criterion: c for c in r2.TwinRung2().judge(rows)}["E2-13 the control grows"].passed
    # the case that must fail: none held, its losses each with an entry, beside a control that held
    fail = next(s for s in shots.values() if s.kind == "fail")
    held = next(s for s in shots.values() if s.kind == "control")
    energies = (0.5, 1.0, 2.0, 4.0)
    entries = Losses.load(
        {
            "index": list(range(8)), "how": ["struck"] * 8, "step": [1] * 8, "time": [1e-11] * 8, "s": [0.0] * 8,
            "laps": [0.0] * 8, "u": [0.0] * 8, "v": [0.0] * 8, "energy": [1.0] * 8, "sense": [1] * 8,
            "point": [[0.0, 0.0, 0.0]] * 8, "length": 0.3,
        }
    )  # fmt: skip
    nothing = r2.Scan("fail", {}, 0, _table(energies, (0.0, 0.0), 200.0, 1), (entries,), 1, 1e-12, 915e6, 123.0, 1.2,
                      (0.0, 0.0, 0.0), False, 0)  # fmt: skip
    kept = replace(nothing, kind="control", table=_table(energies, (9.0, 9.0), 200.0, 1), losses=())
    prediction = r2.Prediction("made up", 0.0, 0.0, 0.0, "held", {})
    checks = {c.criterion: c for c in r2.TwinRung2().judge([(fail, prediction, nothing), (held, prediction, kept)])}
    assert checks["must fail control"].passed
    assert checks["must fail"].passed
    checks = {c.criterion: c for c in r2.TwinRung2().judge([(fail, prediction, nothing)])}
    assert not checks["must fail"].passed  # no control on file: the failure says nothing alone
