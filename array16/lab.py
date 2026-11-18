"""The experiment harness: every experiment predicts, runs and reports the same way, on whatever bench it is given.

An experiment is a set of shots. For each shot it has a prediction, made from theory alone, and a way to measure
it. What it measures with is a bench: the compute and the hardware it may use. The harness

- writes the predictions, and will not measure until they are committed, so the history shows which came first;
- measures the shots not yet on file, in parallel where no device is involved, and keeps each measurement with
  the commit, the bench and the versions that made it, so a run can be stopped and picked up again;
- judges the measurements against the predictions by the experiment's own criteria, and writes the verdicts.

  Experiment        what an experiment supplies: shots, predict, measure, judge
  Bench, Device     a resource setup: workers, and devices by role, each a simulator or a hardware driver
  Check             one measured number set against one predicted number
  predict, run, report, main     the three steps, and the command line every experiment script shares

A bench is a TOML file. A device's driver is any callable named as "module:attribute"; the harness calls it with
the device's options. A simulator and the hardware it stands in for are two drivers for one role:

    name = "simulated"
    workers = 8

    [device.array]
    driver = "array16.sim:SimBackend"
    kind = "quad"
"""

from __future__ import annotations

import argparse
import cmath
import dataclasses
import hashlib
import importlib
import json
import math
import os
import platform
import subprocess
import sys
import time
import tomllib
from abc import ABC, abstractmethod
from collections.abc import Iterator, Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Generic, TypeVar

import numpy as np

S = TypeVar("S")  # a shot
P = TypeVar("P")  # a prediction
M = TypeVar("M")  # a measurement

ROOT = Path(__file__).resolve().parents[1]
BENCHES = ROOT / "benches"


class Refusal(RuntimeError):
    """The harness will not go on, and the message says what would let it."""


# ---- plain data ----------------------------------------------------------------------------------------------------
def plain(value: Any) -> Any:
    """A value as JSON can hold it: a dataclass as a dict, a complex number as [re, im], a non-finite number as
    null, numpy scalars and arrays as Python numbers and lists."""
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {f.name: plain(getattr(value, f.name)) for f in dataclasses.fields(value)}
    if isinstance(value, np.ndarray):
        return plain(value.tolist())
    if isinstance(value, np.generic):
        return plain(value.item())
    if isinstance(value, complex):
        return [value.real, value.imag] if cmath.isfinite(value) else None
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, Mapping):
        return {str(key): plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    return value


def table(columns: Sequence[str], rows: Sequence[Sequence[str]]) -> list[str]:
    """A Markdown table, a line per list item."""
    return [
        "| " + " | ".join(columns) + " |",
        "|" + " --- |" * len(columns),
        *("| " + " | ".join(row) + " |" for row in rows),
    ]


# ---- the bench: compute and devices ---------------------------------------------------------------------------------
@dataclass(frozen=True)
class Device:
    """One device of a bench: the role an experiment knows it by, and the driver that plays it."""

    role: str
    driver: str  # "module:attribute", a callable taking the options
    options: Mapping[str, Any] = field(default_factory=dict)

    def open(self) -> Any:
        module, _, attribute = self.driver.partition(":")
        if not attribute:
            raise Refusal(f"device `{self.role}`: the driver `{self.driver}` should read `module:attribute`")
        try:
            factory = getattr(importlib.import_module(module), attribute)
        except (ImportError, AttributeError) as error:
            raise Refusal(f"device `{self.role}`: no driver `{self.driver}` ({error})") from error
        return factory(**self.options)


@dataclass(frozen=True)
class Bench:
    """A resource setup: how many workers compute may use, and the devices at hand, by role."""

    name: str = "local"
    workers: int = 1
    devices: tuple[Device, ...] = ()

    @classmethod
    def local(cls, workers: int | None = None) -> Bench:
        """This machine and no devices: all its cores but two, unless told otherwise."""
        return cls("local", workers or max(1, (os.cpu_count() or 1) - 2))

    @classmethod
    def load(cls, path: Path) -> Bench:
        spec = tomllib.loads(path.read_text())
        known = {"name", "workers", "device"}
        if set(spec) - known:
            raise Refusal(f"{path}: unknown keys {sorted(set(spec) - known)}; a bench has {sorted(known)}")
        devices = []
        for role, entry in spec.get("device", {}).items():
            if "driver" not in entry:
                raise Refusal(f"{path}: device `{role}` names no driver")
            devices.append(Device(role, entry["driver"], {k: v for k, v in entry.items() if k != "driver"}))
        workers = spec.get("workers") or cls.local().workers
        return cls(spec.get("name", path.stem), int(workers), tuple(devices))

    @classmethod
    def named(cls, name: str, workers: int | None = None) -> Bench:
        """The bench in benches/NAME.toml. `local` needs no file."""
        path = BENCHES / f"{name}.toml"
        if not path.exists():
            if name == "local":
                return cls.local(workers)
            have = sorted(p.stem for p in BENCHES.glob("*.toml"))
            raise Refusal(f"no bench `{name}`: {path} does not exist (benches: {', '.join(have) or 'none'})")
        bench = cls.load(path)
        return bench if workers is None else dataclasses.replace(bench, workers=workers)

    @contextmanager
    def session(self, needs: Mapping[str, Sequence[str]]) -> Iterator[dict[str, Any]]:
        """Open the devices an experiment needs, check that each can do what is asked of it, and close them after."""
        by_role = {device.role: device for device in self.devices}
        missing = sorted(set(needs) - set(by_role))
        if missing:
            raise Refusal(
                f"the bench `{self.name}` has no {', '.join(f'`{role}`' for role in missing)}: add "
                f"[device.{missing[0]}] with a driver to its file, or choose another bench with --bench"
            )
        opened: dict[str, Any] = {}
        try:
            for role, abilities in needs.items():
                opened[role] = by_role[role].open()
                lacking = [ability for ability in abilities if not hasattr(opened[role], ability)]
                if lacking:
                    raise Refusal(f"the bench `{self.name}`: its `{role}` ({by_role[role].driver}) has no {lacking}")
            yield opened
        finally:
            for device in opened.values():
                close = getattr(device, "close", None)
                if callable(close):
                    close()

    def describe(self, opened: Mapping[str, Any] | None = None) -> dict[str, Any]:
        """The bench as a record keeps it. A device that has `identity()` (a serial number, a firmware version)
        is asked for it."""
        devices: dict[str, Any] = {}
        for device in self.devices:
            entry: dict[str, Any] = {"driver": device.driver, "options": dict(device.options)}
            identity = getattr((opened or {}).get(device.role), "identity", None)
            if callable(identity):
                entry["identity"] = identity()
            devices[device.role] = entry
        return {"name": self.name, "workers": self.workers, "devices": devices}


# ---- what an experiment supplies ------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Check:
    """One measured number set against one predicted number."""

    criterion: str
    ident: str
    passed: bool
    said: str


@dataclass(frozen=True)
class Tally:
    """How one criterion fared."""

    passed: int
    judged: int
    failed: tuple[str, ...]  # each failure, in words

    @classmethod
    def of(cls, checks: Sequence[Check]) -> Tally:
        failed = tuple(f"{check.ident}: {check.said}" for check in checks if not check.passed)
        return cls(sum(check.passed for check in checks), len(checks), failed)


class Experiment(ABC, Generic[S, P, M]):
    """An experiment, as the harness sees it. Shots, predictions and measurements are dataclasses."""

    out: Path  # where its files go
    needs: Mapping[str, Sequence[str]] = {}  # device roles, and what each must be able to do
    registered: Sequence[Path] = ()  # files that must be committed before a run, besides the predictions

    @abstractmethod
    def shots(self) -> Sequence[S]:
        """Every shot, in a fixed order."""

    @abstractmethod
    def ident(self, shot: S) -> str:
        """A shot's name: unique, readable, and the same on every run."""

    @abstractmethod
    def predict(self, shot: S) -> P:
        """What theory says the shot does. No measurement may enter it."""

    @abstractmethod
    def measure(self, batch: Sequence[S], devices: Mapping[str, Any]) -> Mapping[str, M]:
        """Measure a batch of shots, on the devices asked for in `needs`: each shot's measurement by its ident. A
        shot that cannot be taken is left out."""

    @abstractmethod
    def load(self, record: Mapping[str, Any]) -> M:
        """A measurement back from the plain data it was stored as."""

    @abstractmethod
    def judge(self, rows: Sequence[tuple[S, P, M | None]]) -> Sequence[Check]:
        """Set the measurements against the predictions, by criteria fixed before the run."""

    def batches(self, shots: Sequence[S]) -> Sequence[Sequence[S]]:
        """How shots go together into one call of `measure`. By default, one at a time."""
        return [[shot] for shot in shots]

    def apparatus(self) -> Mapping[str, Any]:
        """What the predictions file should say about the apparatus."""
        return {}

    def summary(self, checks: Sequence[Check]) -> Mapping[str, Sequence[Check]]:
        """Extra groupings of the checks for the report: a title, and the checks under it."""
        return {}


# ---- files ----------------------------------------------------------------------------------------------------------
def _git(*args: str) -> str | None:
    done = subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True, check=False)
    return done.stdout.strip() if done.returncode == 0 else None


def committed(paths: Sequence[Path]) -> str | None:
    """The commit at HEAD if every path is tracked there and unchanged since, else None."""
    head = _git("rev-parse", "HEAD")
    if head is None:
        return None
    for path in paths:
        if not path.exists() or _git("ls-files", "--error-unmatch", str(path)) is None:
            return None
        if _git("status", "--porcelain", "--", str(path)):
            return None
    return head


def _predictions(experiment: Experiment[S, P, M]) -> dict[str, Any]:
    shots = [
        {"ident": experiment.ident(shot), "shot": plain(shot), "prediction": plain(experiment.predict(shot))}
        for shot in experiment.shots()
    ]
    idents = [entry["ident"] for entry in shots]
    if len(set(idents)) != len(idents):
        twice = sorted({ident for ident in idents if idents.count(ident) > 1})
        raise Refusal(f"two shots share a name: {twice[:3]}")
    return {"apparatus": plain(experiment.apparatus()), "shots": shots}


def _text(data: Any) -> str:
    return json.dumps(data, indent=1, sort_keys=True) + "\n"


def predict(experiment: Experiment[S, P, M]) -> Path:
    """Write the predictions file. Theory only."""
    experiment.out.mkdir(parents=True, exist_ok=True)
    path = experiment.out / "predictions.json"
    path.write_text(_text(_predictions(experiment)))
    return path


def measured(experiment: Experiment[S, P, M]) -> dict[str, dict[str, Any]]:
    """Every record on file, by its shot's ident. A later record of the same shot replaces an earlier one."""
    path = experiment.out / "runs.jsonl"
    lines = path.read_text().splitlines() if path.exists() else []
    return {record["ident"]: record for record in map(json.loads, lines)}


def _measure(experiment: Experiment[S, P, M], batch: Sequence[S], devices: Mapping[str, Any]) -> tuple[Any, float]:
    began = time.time()
    found = experiment.measure(batch, devices)
    return {ident: plain(measurement) for ident, measurement in found.items()}, time.time() - began


def run(
    experiment: Experiment[S, P, M],
    bench: Bench,
    *,
    only: str | None = None,
    unregistered: bool = False,
) -> int:
    """Measure every shot not yet on file. Returns how many batches failed.

    The predictions file and the experiment's `registered` files must be committed first. `unregistered` lifts
    that, for trying things out, and marks every record it makes. `only` keeps the shots whose ident contains it.
    """
    path = experiment.out / "predictions.json"
    if not path.exists() or path.read_text() != _text(_predictions(experiment)):
        raise Refusal(f"{path} is missing or not what the code predicts now: run `predict` first")
    at = committed([path, *experiment.registered])
    if at is None and not unregistered:
        raise Refusal(
            "the predictions are not committed, so nothing would show they came first: commit "
            f"{', '.join(str(p.relative_to(ROOT)) for p in [path, *experiment.registered])} and run again "
            "(or pass --unregistered to try things out; the records will say so)"
        )
    have = measured(experiment)
    todo = [s for s in experiment.shots() if experiment.ident(s) not in have and (only or "") in experiment.ident(s)]
    batches = [batch for batch in experiment.batches(todo) if batch]
    print(
        f"{len(todo)} shots to measure in {len(batches)} batches ({len(have)} on file); bench `{bench.name}`",
        flush=True,
    )
    stamp = {
        "commit": _git("rev-parse", "HEAD"),
        "dirty": bool(_git("status", "--porcelain", "--untracked-files=no")),
        "registered": at,
        "predictions_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "machine": f"{platform.system()} {platform.machine()}",
    }
    by_ident = {experiment.ident(shot): shot for shot in todo}
    failed = 0

    def keep(found: Mapping[str, Any], seconds: float, batch: Sequence[S], about: Mapping[str, Any]) -> None:
        made = {**stamp, "bench": about, "utc": datetime.now(UTC).isoformat(timespec="seconds")}
        with (experiment.out / "runs.jsonl").open("a") as out:
            for ident, measurement in found.items():
                record = {"ident": ident, "shot": plain(by_ident[ident]), "measurement": measurement}
                out.write(json.dumps({**record, "seconds": seconds / len(batch), "provenance": made}) + "\n")
        for shot in batch:
            ident = experiment.ident(shot)
            print(f"  {ident}: {'measured' if ident in found else 'not taken'}", flush=True)

    with bench.session(experiment.needs) as devices:
        about = bench.describe(devices)
        if devices or bench.workers <= 1 or len(batches) <= 1:  # a device is used by one process, in order
            for batch in batches:
                try:
                    found, seconds = _measure(experiment, batch, devices)
                except Exception as error:
                    failed += 1
                    print(f"  FAILED {experiment.ident(batch[0])}: {error!r}", file=sys.stderr, flush=True)
                    continue
                keep(found, seconds, batch, about)
        else:
            with ProcessPoolExecutor(max_workers=min(bench.workers, len(batches))) as pool:
                pending = {pool.submit(_measure, experiment, batch, {}): batch for batch in batches}
                for future in as_completed(pending):
                    batch = pending[future]
                    try:
                        found, seconds = future.result()
                    except Exception as error:
                        failed += 1
                        print(f"  FAILED {experiment.ident(batch[0])}: {error!r}", file=sys.stderr, flush=True)
                        continue
                    keep(found, seconds, batch, about)
    return failed


def rows(experiment: Experiment[S, P, M]) -> list[tuple[S, P, M | None]]:
    """Each shot with its prediction and, if it is on file, its measurement."""
    have = measured(experiment)
    out: list[tuple[S, P, M | None]] = []
    for shot in experiment.shots():
        record = have.get(experiment.ident(shot))
        out.append((shot, experiment.predict(shot), None if record is None else experiment.load(record["measurement"])))
    return out


def report(experiment: Experiment[S, P, M]) -> dict[str, Any]:
    """Judge what is on file, and write results.json and results.md."""
    have = measured(experiment)
    table_rows = rows(experiment)
    checks = list(experiment.judge(table_rows))
    groups: dict[str, list[Check]] = {}
    for check in checks:
        groups.setdefault(check.criterion, []).append(check)
    groups.update({title: list(members) for title, members in experiment.summary(checks).items()})
    path = experiment.out / "predictions.json"
    stamps = [record["provenance"] for record in have.values()]
    integrity = {
        "shots": len(table_rows),
        "measured": sum(m is not None for _, _, m in table_rows),
        "predictions_match_the_code": path.exists() and path.read_text() == _text(_predictions(experiment)),
        "records_made_before_registration": sum(stamp["registered"] is None for stamp in stamps),
        "records_from_other_predictions": sum(
            stamp["predictions_sha256"] != hashlib.sha256(path.read_bytes()).hexdigest() for stamp in stamps
        )
        if path.exists()
        else len(stamps),
        "records_from_uncommitted_code": sum(bool(stamp["dirty"]) for stamp in stamps),
        "commits": sorted({str(stamp["commit"])[:12] for stamp in stamps}),
        "benches": sorted({stamp["bench"]["name"] for stamp in stamps}),
    }
    tally = {title: Tally.of(members) for title, members in groups.items()}
    result = {"integrity": integrity, "summary": plain(tally), "checks": plain(checks)}
    (experiment.out / "results.json").write_text(_text(result))
    lines = ["<!-- Written by the `report` step of array16.lab. -->", "", "# The verdicts", ""]
    lines += [f"{integrity['measured']} of {integrity['shots']} shots are on file.", ""]
    lines += table(["Integrity", ""], [[key.replace("_", " "), str(value)] for key, value in integrity.items()])
    lines += [
        "",
        *table(["Criterion", "Passed", "Judged"], [[t, str(v.passed), str(v.judged)] for t, v in tally.items()]),
    ]
    for title, entry in tally.items():
        if entry.failed:
            lines += ["", f"## Failed: {title}", "", *(f"- {line}" for line in entry.failed)]
    said = [[c.criterion, c.ident, "pass" if c.passed else "**FAIL**", c.said] for c in checks]
    lines += ["", "## Every comparison", "", *table(["Criterion", "Shot", "Result", "Numbers"], said)]
    (experiment.out / "results.md").write_text("\n".join(lines) + "\n")
    return result


def main(experiment: Experiment[S, P, M], argv: Sequence[str] | None = None, about: str = "") -> int:
    """The command line every experiment script shares: predict, run, report."""
    parser = argparse.ArgumentParser(description=about)
    step = parser.add_subparsers(dest="step", required=True)
    step.add_parser("predict", help="write the predictions: theory only")
    go = step.add_parser("run", help="measure what is not yet on file")
    go.add_argument("--bench", default="local", help="the resource setup: benches/NAME.toml (default: local)")
    go.add_argument("--workers", type=int, help="processes for compute, instead of the bench's")
    go.add_argument("--only", help="only the shots whose name contains this")
    go.add_argument("--unregistered", action="store_true", help="measure before the predictions are committed")
    step.add_parser("report", help="judge the measurements against the predictions")
    args = parser.parse_args(argv)
    try:
        if args.step == "predict":
            print(f"predictions written: {predict(experiment)}")
        elif args.step == "run":
            failed = run(
                experiment, Bench.named(args.bench, args.workers), only=args.only, unregistered=args.unregistered
            )
            if failed:
                print(f"{failed} batches failed", file=sys.stderr)
                return 1
        else:
            result = report(experiment)
            for key, value in result["integrity"].items():
                print(f"{key.replace('_', ' ')}: {value}")
            for title, entry in result["summary"].items():
                print(f"{title}: {entry['passed']} of {entry['judged']}")
    except Refusal as refusal:
        print(f"refused: {refusal}", file=sys.stderr)
        return 2
    return 0
