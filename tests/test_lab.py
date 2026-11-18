"""The experiment harness on a toy experiment: timing a dropped ball."""

from __future__ import annotations

import json
import math
import subprocess
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from array16 import lab
from array16.lab import Bench, Check, Device, Experiment, Refusal


@dataclass(frozen=True)
class Drop:
    height: float  # m


@dataclass(frozen=True)
class Fall:
    seconds: float


class Drops(Experiment[Drop, Fall, Fall]):
    """Predict the time of a fall from Galileo's law; "measure" it in a world whose gravity can be set."""

    def __init__(self, out: Path, gravity: float = 9.81, theory: float = 9.81, broken: float | None = None) -> None:
        self.out, self.gravity, self.theory, self.broken = out, gravity, theory, broken
        self.registered = (out / "PROTOCOL.md",)

    def shots(self) -> Sequence[Drop]:
        return [Drop(1.0), Drop(2.0), Drop(5.0)]

    def ident(self, shot: Drop) -> str:
        return f"drop from {shot.height:g} m"

    def predict(self, shot: Drop) -> Fall:
        return Fall(math.sqrt(2 * shot.height / self.theory))

    def measure(self, batch: Sequence[Drop], devices: Mapping[str, Any]) -> Mapping[str, Fall]:
        if any(shot.height == self.broken for shot in batch):
            raise OSError("the ball rolled under the bench")
        return {self.ident(shot): Fall(math.sqrt(2 * shot.height / self.gravity)) for shot in batch}

    def load(self, record: Mapping[str, Any]) -> Fall:
        return Fall(record["seconds"])

    def judge(self, rows: Sequence[tuple[Drop, Fall, Fall | None]]) -> Sequence[Check]:
        return [
            Check("C1 time", self.ident(s), abs(m.seconds - p.seconds) < 0.01, f"{p.seconds:.3f} s, {m.seconds:.3f} s")
            for s, p, m in rows
            if m is not None
        ]


class Clock:
    """A stand-in for a device: it times falls, and knows whether it was put away."""

    put_away = 0

    def __init__(self, gravity: float = 9.81) -> None:
        self.gravity = gravity

    def time(self, height: float) -> float:
        return math.sqrt(2 * height / self.gravity)

    def identity(self) -> dict[str, str]:
        return {"serial": "clock-0001"}

    def close(self) -> None:
        Clock.put_away += 1


class TimedDrops(Drops):
    """The same experiment, measured with a device from the bench."""

    needs = {"clock": ("time",)}  # noqa: RUF012

    def measure(self, batch: Sequence[Drop], devices: Mapping[str, Any]) -> Mapping[str, Fall]:
        return {self.ident(shot): Fall(devices["clock"].time(shot.height)) for shot in batch}


def _git(repo: Path, *args: str) -> None:
    subprocess.run(["git", "-C", str(repo), "-c", "commit.gpgsign=false", *args], check=True, capture_output=True)


@pytest.fixture
def repo(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """An empty repository standing in for this one."""
    _git(tmp_path, "init", "-q", "-b", "main")
    _git(tmp_path, "config", "user.name", "a test")
    _git(tmp_path, "config", "user.email", "test@example.invalid")
    (tmp_path / "exp").mkdir()
    (tmp_path / "exp" / "PROTOCOL.md").write_text("what counts as agreement: 10 ms\n")
    (tmp_path / "benches").mkdir()
    (tmp_path / "benches" / "table.toml").write_text(
        'workers = 1\n\n[device.clock]\ndriver = "test_lab:Clock"\ngravity = 9.81\n'
    )
    _git(tmp_path, "add", "exp/PROTOCOL.md", "benches/table.toml")
    _git(tmp_path, "commit", "-q", "-m", "the protocol")
    monkeypatch.setattr(lab, "ROOT", tmp_path)
    monkeypatch.setattr(lab, "BENCHES", tmp_path / "benches")
    return tmp_path


def _register(repo: Path) -> None:
    _git(repo, "add", "exp/predictions.json")
    _git(repo, "commit", "-q", "-m", "the predictions")


# ---- plain data ----------------------------------------------------------------------------------------------------
def test_plain_makes_json_of_dataclasses_complex_numbers_and_arrays() -> None:
    value = {"fall": Fall(1.5), "z": 1 + 2j, "nowhere": complex(math.nan, 0), "far": math.inf, "a": np.arange(3)}
    assert lab.plain(value) == {"fall": {"seconds": 1.5}, "z": [1.0, 2.0], "nowhere": None, "far": None, "a": [0, 1, 2]}
    assert lab.plain((np.float64(2.5), Path("a/b"), None, "s")) == [2.5, "a/b", None, "s"]
    json.dumps(lab.plain(value))


def test_table_is_markdown() -> None:
    assert lab.table(["a", "b"], [["1", "2"]]) == ["| a | b |", "| --- | --- |", "| 1 | 2 |"]


# ---- the bench -----------------------------------------------------------------------------------------------------
def test_a_bench_loads_from_its_file(repo: Path) -> None:
    bench = Bench.named("table")
    assert (bench.name, bench.workers) == ("table", 1)
    assert bench.devices == (Device("clock", "test_lab:Clock", {"gravity": 9.81}),)
    assert Bench.named("table", workers=4).workers == 4
    assert Bench.named("local").devices == ()
    assert Bench.local().workers >= 1


def test_a_bench_file_with_a_mistake_is_refused_by_name(repo: Path) -> None:
    (repo / "benches" / "typo.toml").write_text("wokers = 2\n")
    with pytest.raises(Refusal, match=r"unknown keys.*wokers"):
        Bench.named("typo")
    (repo / "benches" / "bare.toml").write_text("[device.clock]\ngravity = 1.0\n")
    with pytest.raises(Refusal, match="names no driver"):
        Bench.named("bare")
    with pytest.raises(Refusal, match=r"no bench `attic`.*benches: bare, table, typo"):
        Bench.named("attic")


def test_a_session_opens_what_is_needed_checks_it_and_puts_it_away(repo: Path) -> None:
    bench, before = Bench.named("table"), Clock.put_away
    with bench.session({"clock": ("time",)}) as devices:
        assert devices["clock"].time(4.905) == pytest.approx(1.0)
        assert bench.describe(devices)["devices"]["clock"]["identity"] == {"serial": "clock-0001"}
    assert Clock.put_away == before + 1
    with bench.session({}) as devices:  # nothing asked for, nothing opened
        assert devices == {}
    assert Clock.put_away == before + 1


def test_a_session_refuses_what_the_bench_cannot_supply(repo: Path) -> None:
    bench, before = Bench.named("table"), Clock.put_away
    with pytest.raises(Refusal, match=r"no `radio`: add \[device.radio\]"), bench.session({"radio": ()}):
        pass
    with pytest.raises(Refusal, match=r"its `clock` .* has no \['chime'\]"), bench.session({"clock": ("chime",)}):
        pass
    assert Clock.put_away == before + 1  # opened, found wanting, and put away all the same
    with pytest.raises(Refusal, match="should read `module:attribute`"):
        Device("clock", "test_lab").open()
    with pytest.raises(Refusal, match="no driver `test_lab:Sundial`"):
        Device("clock", "test_lab:Sundial").open()


def test_the_simulated_bench_supplies_an_array_that_can_do_what_the_backend_asks() -> None:
    # The repository's own bench file, and the simulator as a driver: the path a hardware driver will take.
    bench = Bench.load(Path(__file__).resolve().parents[1] / "benches" / "simulated.toml")
    assert bench.name == "simulated"
    backend = ("controller", "tune", "capture", "cal_signal", "set_tx")  # array16/backend.py::Backend
    with bench.session({"array": backend}) as devices:
        assert devices["array"].kind == "quad"
        assert devices["array"].capture(64).samples.shape == (4, 64)
    assert bench.describe()["devices"]["array"]["options"] == {"radio": "quad", "seed": 0}


# ---- predict, then run ---------------------------------------------------------------------------------------------
def test_run_refuses_until_the_predictions_are_committed(repo: Path) -> None:
    drops = Drops(repo / "exp")
    with pytest.raises(Refusal, match="run `predict` first"):
        lab.run(drops, Bench())
    lab.predict(drops)
    with pytest.raises(Refusal, match=r"not committed.*exp/predictions.json, exp/PROTOCOL.md"):
        lab.run(drops, Bench())
    assert lab.measured(drops) == {}
    _register(repo)
    assert lab.run(drops, Bench()) == 0
    assert len(lab.measured(drops)) == 3


def test_a_record_says_what_made_it(repo: Path) -> None:
    drops = Drops(repo / "exp")
    lab.predict(drops)
    _register(repo)
    lab.run(drops, Bench("a laptop"))
    record = lab.measured(drops)["drop from 2 m"]
    assert record["shot"] == {"height": 2.0}
    assert record["measurement"]["seconds"] == pytest.approx(math.sqrt(4 / 9.81))
    made = record["provenance"]
    assert len(made["registered"]) == 40
    assert made["registered"] == made["commit"]
    assert made["bench"] == {"name": "a laptop", "workers": 1, "devices": {}}
    assert made["numpy"] == np.__version__
    assert len(made["predictions_sha256"]) == 64
    assert not made["dirty"]


def test_an_unregistered_run_is_allowed_and_marked(repo: Path) -> None:
    drops = Drops(repo / "exp")
    lab.predict(drops)
    assert lab.run(drops, Bench(), unregistered=True) == 0
    assert all(record["provenance"]["registered"] is None for record in lab.measured(drops).values())
    assert lab.report(drops)["integrity"]["records_made_before_registration"] == 3


def test_a_run_picks_up_where_it_stopped(repo: Path, capsys: pytest.CaptureFixture[str]) -> None:
    drops = Drops(repo / "exp", broken=2.0)
    lab.predict(drops)
    _register(repo)
    assert lab.run(drops, Bench()) == 1  # one batch fails, and the others are kept
    assert sorted(lab.measured(drops)) == ["drop from 1 m", "drop from 5 m"]
    assert "FAILED drop from 2 m" in capsys.readouterr().err
    mended = Drops(repo / "exp")
    assert lab.run(mended, Bench()) == 0
    assert "1 shots to measure in 1 batches (2 on file)" in capsys.readouterr().out
    assert len(lab.measured(mended)) == 3
    assert lab.run(mended, Bench()) == 0
    assert "0 shots to measure" in capsys.readouterr().out


def test_only_keeps_the_shots_it_names(repo: Path) -> None:
    drops = Drops(repo / "exp")
    lab.predict(drops)
    _register(repo)
    lab.run(drops, Bench(), only="from 5")
    assert list(lab.measured(drops)) == ["drop from 5 m"]


def test_many_workers_measure_what_one_does(repo: Path) -> None:
    drops = Drops(repo / "exp")
    lab.predict(drops)
    _register(repo)
    assert lab.run(drops, Bench(workers=3)) == 0
    several = {ident: record["measurement"] for ident, record in lab.measured(drops).items()}
    (repo / "exp" / "runs.jsonl").unlink()
    assert lab.run(drops, Bench(workers=1)) == 0
    assert {ident: record["measurement"] for ident, record in lab.measured(drops).items()} == several


def test_a_device_measures_and_its_identity_is_kept(repo: Path) -> None:
    drops = TimedDrops(repo / "exp")
    lab.predict(drops)
    _register(repo)
    with pytest.raises(Refusal, match="no `clock`"):
        lab.run(drops, Bench())  # the local bench has no clock
    assert lab.run(drops, Bench.named("table", workers=8)) == 0  # with a device, one process, whatever the workers
    record = lab.measured(drops)["drop from 1 m"]
    assert record["provenance"]["bench"]["devices"]["clock"] == {
        "driver": "test_lab:Clock",
        "options": {"gravity": 9.81},
        "identity": {"serial": "clock-0001"},
    }
    assert lab.report(drops)["integrity"]["benches"] == ["table"]


# ---- report --------------------------------------------------------------------------------------------------------
def test_report_tallies_the_checks_and_names_the_failures(repo: Path) -> None:
    drops = Drops(repo / "exp", gravity=9.0)  # a world that does not agree with the theory
    lab.predict(drops)
    _register(repo)
    lab.run(drops, Bench())
    result = lab.report(drops)
    assert result["summary"]["C1 time"]["judged"] == 3
    assert result["summary"]["C1 time"]["passed"] == 0
    assert result["integrity"]["predictions_match_the_code"]
    assert result["integrity"]["records_made_before_registration"] == 0
    text = (repo / "exp" / "results.md").read_text()
    assert "## Failed: C1 time" in text
    assert "| C1 time | drop from 1 m | **FAIL** |" in text


def test_a_theory_changed_after_registration_shows(repo: Path) -> None:
    drops = Drops(repo / "exp")
    lab.predict(drops)
    _register(repo)
    lab.run(drops, Bench())
    revised = Drops(repo / "exp", theory=9.0)
    assert not lab.report(revised)["integrity"]["predictions_match_the_code"]
    with pytest.raises(Refusal, match="not what the code predicts now"):
        lab.run(revised, Bench())
    lab.predict(revised)  # the predictions rewritten after the fact: the records carry the old file's hash
    assert lab.report(revised)["integrity"]["records_from_other_predictions"] == 3


def test_two_shots_with_one_name_are_refused(repo: Path) -> None:
    class Twice(Drops):
        def shots(self) -> Sequence[Drop]:
            return [Drop(1.0), Drop(1.0)]

    with pytest.raises(Refusal, match="two shots share a name"):
        lab.predict(Twice(repo / "exp"))


# ---- the command line ----------------------------------------------------------------------------------------------
def test_the_command_line_runs_the_three_steps(repo: Path, capsys: pytest.CaptureFixture[str]) -> None:
    drops = Drops(repo / "exp")
    assert lab.main(drops, ["run"]) == 2
    assert "refused:" in capsys.readouterr().err
    assert lab.main(drops, ["predict"]) == 0
    assert lab.main(drops, ["run"]) == 2  # not committed yet
    _register(repo)
    assert lab.main(drops, ["run", "--workers", "1"]) == 0
    assert lab.main(drops, ["report"]) == 0
    assert "C1 time: 3 of 3" in capsys.readouterr().out
    assert lab.main(drops, ["run", "--bench", "attic"]) == 2
