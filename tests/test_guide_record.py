"""Experiment 1's record, as committed: the verdicts it gives, and what was read from it afterwards."""

from __future__ import annotations

import json
import math

import guide_posthoc
import guide_study

from array16 import lab


def test_gde_h1_the_record_gives_the_verdicts_on_file() -> None:
    # The registered criteria, applied again to the measurements on file, give the counts in results.json:
    # the hypothesis passes whole on the ring, on scaling and on the instrument, and fails on the turn per lap.
    experiment = guide_study.GuideExperiment()
    rows = lab.rows(experiment)
    checks = experiment.judge(rows)
    on_file = json.loads((guide_study.OUT / "results.json").read_text())["summary"]
    tally = {name: lab.Tally.of([c for c in checks if c.criterion == name]) for name in {c.criterion for c in checks}}
    assert {name: (t.passed, t.judged) for name, t in tally.items()} == {
        name: (entry["passed"], entry["judged"]) for name, entry in on_file.items() if not name.startswith(("C7", "C9"))
    }
    assert (tally["C6 ring"].passed, tally["C6 ring"].judged) == (10, 10)
    assert (tally["C8 scaling"].passed, tally["C8 scaling"].judged) == (16, 16)
    assert (tally["C4 turn per lap"].passed, tally["C4 turn per lap"].judged) == (15, 22)


def test_gde_h2_the_shared_shift_of_the_turn_is_half_rho_times_the_bending() -> None:
    # Read from the record after the run, so a description of it and not a test passed: on every pair of
    # on-axis shots held both ways, what the two fitted turns share beyond the prediction is rho x bending / 2.
    pairs = guide_posthoc.pairs()
    assert len(pairs) == 8
    assert {pair.guide for pair in pairs} == {"trefoil", "trefoil-shallow", "helix", "figure-eight"}
    for pair in pairs:
        close = 0.15 if pair.guide == "trefoil-shallow" else 0.02
        assert abs(pair.shared / pair.candidate - 1) < close, pair
        assert pair.shared > math.radians(1.5)  # and it is not small: up to 12 degrees a lap
