"""experiments/twin_rung1_figures.py, the figures of experiment 2's stage 1 report, drawn from the committed records:
that it reads a scan as the verdicts did (every hold's limit and the knot's price as results.json states them), builds
each path at its registered length and frame (the in-plane radii and tightest bends the predictions rest on), places a
lost electron by its laps, holds the numbers the report quotes to the records, and draws its six figures from the
records alone, into any folder, from records in any folder."""

from __future__ import annotations

import json
import math
import re
import shutil
from itertools import pairwise
from pathlib import Path

import numpy as np
import pytest
import twin_rung1 as rung
import twin_rung1_figures as figures

HOLDS = {
    "E2-21 hold, R1, realisation 0": (9.344, 1.004, True),
    "E2-21 hold, R1, on the bench": (8.465, 0.910, True),
    "E2-21 hold, R1, realisation 0 at q 0.3": (3.269, 0.997, True),
    "E2-22 hold, T915-S, realisation 0": (9.996, 0.859, False),
    "E2-22 hold, flat ring 318, realisation 0": (32.25, 0.994, True),
    "E2-22 hold, helix (1, 3), realisation 0": (16.10, 0.873, False),
}


def test_the_records_read_as_the_verdicts() -> None:
    # every hold's limit and every price of the knot as results.json words them, from twin_rung1's own readings
    rec = figures.records()
    assert len(rec.rows) == 14
    said = {(c["criterion"], c["ident"]): str(c["said"]) for c in rec.checks}
    for row in rec.of("hold"):
        h = figures.hold(rec, row)
        criterion = next(
            c for c, ident in said if ident == row.ident and c not in ("E2-21 senses", "E2-22 price of the knot")
        )
        assert said[(criterion, row.ident)].startswith(f"held to {h.limit:.4g} eV")
        assert h.passed == rec.verdict(criterion, row.ident)
        limit, fraction, planar = HOLDS[row.ident]
        assert h.limit == pytest.approx(limit, rel=6e-4)
        assert round(h.fraction, 3) == fraction
        assert h.planar is planar
    for row in rec.of("hold"):
        if row.prediction["ratio"] is not None:
            found = figures.ratio(rec, row)
            assert f"{found:.3f} against the static" in said[("E2-22 price of the knot", row.ident)]
    flat = rec.row("E2-22 hold, flat ring 318, realisation 0")
    assert figures.ratio(rec, flat) == pytest.approx(3.226, abs=5e-4)
    assert not rec.verdict("E2-22 price of the knot", flat.ident)


def test_the_edges_of_each_hold() -> None:
    # the planar rings lose every electron one scan step above their first loss; T915-S's senses lose all from 10.90
    # and 11.17 eV, above islands of loss from 9.17 eV
    rec = figures.records()
    drawn = figures.hold(rec, rec.row("E2-21 hold, R1, realisation 0"))
    for sense in drawn.senses:
        assert (sense.first_loss, sense.all_lost) == pytest.approx((9.2289, 9.4596), abs=1e-4)
    knot = figures.hold(rec, rec.row("E2-22 hold, T915-S, realisation 0"))
    with_wave, against = knot.senses
    assert (with_wave.first_loss, with_wave.all_lost) == pytest.approx((9.1685, 11.1709), abs=1e-4)
    assert (against.first_loss, against.all_lost) == pytest.approx((9.3977, 10.8985), abs=1e-4)
    assert (with_wave.threshold, against.threshold) == pytest.approx((9.282, 10.765), abs=1e-3)
    assert figures.all_lost_from([1.0, 2.0, 3.0, 4.0], [2, 0, 1, 0], [0, 2, 1, 2]) == 4.0
    assert figures.all_lost_from([1.0, 2.0], [2, 2], [0, 0]) == math.inf


def test_the_paths_are_the_designs_at_their_registered_lengths() -> None:
    rec = figures.records()
    for row in rec.of("hold"):
        path, basis = figures.route_of(rec, row), row.basis
        assert path.length == pytest.approx(basis["length_mm"] * 1e-3, rel=1e-9)
        assert 1 / float(np.hypot(path.kappa_g, path.kappa_n).max()) == pytest.approx(basis["bend_mm"] * 1e-3, rel=1e-4)
        if "inplane_radius_mm" in basis:  # framed k = 0: the in-plane radius the predictions rest on
            assert 1 / float(path.kappa_g.max()) == pytest.approx(basis["inplane_radius_mm"] * 1e-3, rel=2e-4)
            stretches = figures.lobes(path)
            assert len(stretches) == 3
            assert all(0 < b - a < 1 / 6 for a, b in stretches)
        else:  # a circle in its own plane: the bend wholly in the boards' plane
            assert 1 / float(path.kappa_g.mean()) == pytest.approx(basis["bend_mm"] * 1e-3, rel=1e-4)
            assert float(path.kappa_n.max()) < 1e-6 * float(path.kappa_g.mean())
            assert figures.lobes(path) == ()


def test_a_loss_is_placed_by_its_laps_in_its_sense() -> None:
    assert figures.Loss(1, 9.0, 3.8784, 1, 8).place == pytest.approx(0.8784)
    assert figures.Loss(-1, 9.0, 1.3497, 3, 8).place == pytest.approx(1 - 0.3497)
    rec = figures.records()
    found = {(loss.sense, round(loss.energy, 4)): loss for loss in figures.losses(rec.row(next(iter(HOLDS))).scan)}
    assert found[(1, 9.2289)].place == pytest.approx(0.878, abs=1e-3)
    assert found[(-1, 9.2289)].place == pytest.approx(0.650, abs=1e-3)
    assert (found[(1, 9.2289)].lost, found[(-1, 9.2289)].lost) == (1, 3)


def test_the_places_the_report_quotes() -> None:
    # the records place only each energy and sense's first-lost electron: 40 places on T915-S, which lost 75, and 38 on
    # the helix, which lost 76. Out of the plane every such place lies in a lobe or at most 0.025 lap (T915-S) and
    # 0.067 lap (the helix) past one in its direction of travel, 8 of the helix's 0.050 to 0.067 lap past, at 16.30 to
    # 18.90 eV; on the rings as drawn every first loss at or above the edge came within 0.08 to 0.17 lap, the four
    # below it after 1.3 to 4.0 laps
    rec = figures.records()
    worst, places, beyond = {}, {}, {}
    for row in rec.of("hold"):
        stretches = figures.lobes(figures.route_of(rec, row))
        if stretches:
            lost = figures.losses(row.scan)
            worst[row.shot["label"]] = max(figures.past(loss, stretches) for loss in lost)
            places[row.shot["label"]] = (len(lost), sum(map(sum, row.scan.lost)))
            beyond[row.shot["label"]] = [loss.energy for loss in lost if figures.past(loss, stretches) >= 0.05]
    assert worst["T915-S"] == pytest.approx(0.025, abs=1e-3)
    assert worst["helix (1, 3)"] == pytest.approx(0.067, abs=1e-3)
    assert places == {"T915-S": (40, 75), "helix (1, 3)": (38, 76)}
    assert beyond["T915-S"] == []
    assert len(beyond["helix (1, 3)"]) == 8
    assert (min(beyond["helix (1, 3)"]), max(beyond["helix (1, 3)"])) == pytest.approx((16.30, 18.90), abs=0.01)
    rings, below = [], []
    for row in rec.of("hold"):
        h = figures.hold(rec, row)
        if row.shot["realisation"] == 0 and h.planar:
            edge = {sense.sense: sense.all_lost for sense in h.senses}
            rings += [loss.when for loss in figures.losses(row.scan) if loss.energy >= edge[loss.sense]]
            below += [loss.when for loss in figures.losses(row.scan) if loss.energy < edge[loss.sense]]
    assert (min(rings), max(rings)) == pytest.approx((0.078, 0.172), abs=1e-3)
    assert len(below) == 4
    assert (min(below), max(below)) == pytest.approx((1.3, 4.0), abs=0.05)
    # H4: on the bench, with the wave, each energy's first loss from 7.958 to 9.460 eV came on the first lap, from
    # 0.76 down to 0.49 lap as the energy rose
    bench = rec.row("E2-21 hold, R1, on the bench")
    first_lap = [loss for loss in figures.losses(bench.scan) if loss.sense == 1 and 7.9 < loss.energy < 9.5]
    assert len(first_lap) == 8
    assert (first_lap[0].energy, first_lap[-1].energy) == pytest.approx((7.958, 9.460), abs=1e-3)
    assert [round(loss.when, 2) for loss in (first_lap[0], first_lap[-1])] == [0.76, 0.49]
    assert all(a.when > b.when for a, b in pairwise(first_lap))


def test_the_must_fail_records_hold_each_energys_first_loss() -> None:
    rec = figures.records()
    fails = rec.of("fail")
    assert len(fails) == 2
    for row in fails:
        scan = row.scan
        assert (scan.held_total(), scan.lost_total()) == (0, 32)
        assert all(n == 4 for counts in scan.lost for n in counts)  # four lost an energy and sense, one time each
        firsts = [x for k in range(2) for x in scan.first_loss[k]]
        assert min(firsts) > 1.3
        assert max(firsts) < 2.3
    control = rec.of("control")[0].scan
    assert (control.held_total(), control.lost_total()) == (32, 0)


def test_the_script_draws_its_figures_from_the_records(tmp_path: Path) -> None:
    out = tmp_path / "figures"
    assert figures.main(["--out", str(out)]) == 0
    assert sorted(p.name for p in out.iterdir()) == sorted(figures.FIGURES)
    for name in figures.FIGURES:
        data = (out / name).read_bytes()
        assert data[:8] == b"\x89PNG\r\n\x1a\n"
        assert len(data) > 30_000


def test_it_reads_the_records_where_they_are_and_a_later_record_wins(tmp_path: Path) -> None:
    # records copied elsewhere, the first shot recorded twice: the later record is the one read, as array16.lab reads
    src = tmp_path / "records"
    src.mkdir()
    for name in ("predictions.json", "results.json"):
        shutil.copy(rung.OUT / name, src / name)
    lines = (rung.OUT / "runs.jsonl").read_text().splitlines()
    again = json.loads(lines[0])
    again["seconds"] = 1234.5
    (src / "runs.jsonl").write_text("\n".join([*lines, json.dumps(again)]) + "\n")
    rec = figures.read(src)
    assert len(rec.rows) == 14
    assert rec.row(again["ident"]).seconds == 1234.5
    assert figures.figure_must_fail(rec, tmp_path).name == "must_fail.png"


def test_the_report_shows_every_figure() -> None:
    # REPORT.md links each figure the script draws, and each is committed beside it
    report = (rung.OUT / "REPORT.md").read_text()
    linked = set(re.findall(r"\]\(([a-z_]+\.png)\)", report))
    assert linked == set(figures.FIGURES)
    assert all((rung.OUT / name).exists() for name in figures.FIGURES)
