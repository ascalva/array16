"""The detector's records of what a run loses (docs/KNOTTED_GUIDE.md section 8.3; H3 of experiments/twin_rung1/
REPORT.md section 6): one entry for each electron lost, never one a step. An entry holds the electron's row in the
launch, how it was lost, the step and the time, where along the path it was (its arc length, and the laps it had gone
from its start) and where in the section (u, v, the section grid's own coordinates), its kinetic energy and its sense
along the path then, and the point in space. A TwinDetector's watch writes them when asked (`losses`, from a bench
file or for one shot); array16.lab keeps them with a measurement as plain data, and Losses.load reads them back.

How an electron is lost:

  struck    it entered a solid cell of the section: the copper or a board
  left      it crossed the grid's edge: out past the aperture
  stopped   another observer stopped it: the gas

Where the field found an electron off the grid, its entry is where it crossed onto the copper, a board or the grid's
edge. Over a step the Boris push moves it in a straight line, from the last place the field found it on the grid to the
first place it found it off, so the crossing lies on that segment: the segment is looked at in sixteenths, so that the
first boundary the step crosses is the one taken, and that sixteenth is halved 30 times, to 2^-34 of the step (54 um a
step at 100 eV on R1 at 120 steps an RF period: 3e-15 m). The entry's time is the crossing's, its energy and sense the
step's. A stopped electron's entry is where the observer stopped it, at the end of its step.

  Losses    the entries, column by column, as array16.lab keeps them (plain) and reads them back (Losses.load)
  LossLog   what a watch keeps to write them: the places its last two calls were given, and the entries so far
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, fields
from typing import Any

import numpy as np
from numpy.typing import NDArray

from array16.particles import C0, ME, Arr, BoolArr, E, IntArr, SectionField
from array16.particles.sections import Placement

__all__ = ["HOW", "LEFT", "STOPPED", "STRUCK", "LossLog", "Losses", "StrArr"]

StrArr = NDArray[np.str_]
STRUCK, LEFT, STOPPED = "struck", "left", "stopped"
HOW = (STRUCK, LEFT, STOPPED)  # how an electron is lost, as its entry says it
_SCAN = 16  # the parts a step is first looked at in, so that the first boundary it crosses is the one taken
_HALVINGS = 30  # the halvings of that part: the crossing to 2^-34 of a step
_ACROSS = 1e-6  # of an electron's speed: less along the path and it moves square across it, sense 0
_KINDS: dict[str, type[np.generic]] = {"index": np.intp, "how": np.str_, "step": np.intp, "sense": np.intp}


def _floats(values: Any) -> Arr:
    """A column of numbers back from plain data, where a null stood for a number that was not finite."""
    out: Arr = np.array([math.nan if x is None else float(x) for x in values], dtype=np.float64)
    return out


@dataclass(frozen=True, eq=False)
class Losses:
    """The electrons a run lost, one entry each, in the order they were lost (by time, then row):

      index    the electron's row in the launch, by which an experiment joins the entry to the energy, sense and phase
               it launched the electron with
      how      how it was lost: struck, left or stopped (HOW)
      step     the step it was lost in, integrate's count: step dt <= time <= (step + 1) dt. The watch stopped it
               after the next step for one that struck or left (Detection.step is one more, or the same for one off
               the grid from its launch), after this one for one stopped
      time     when it was lost (s): when it crossed onto the copper, a board or the grid's edge, or the end of the
               step after which it was stopped
      s        its arc length along the path then (m, 0 to `length`)
      laps     the laps it had gone along the path from its start then, negative against the path's direction
      u, v     its place in the section (m, the section grid's coordinates: the null at the grid's null, u across the
               path in the boards' plane, v from the board toward the null)
      energy   its kinetic energy (eV), over the step in which it was lost
      sense    its sense along the path then: +1 along the path's tangent, -1 against it, 0 square across it (to a
               millionth of its speed)
      point    where it was (m, (K, 3))
      length   the lap (m)

    The columns are read-only arrays of one entry an electron; `point` is (K, 3)."""

    index: IntArr
    how: StrArr
    step: IntArr
    time: Arr
    s: Arr
    laps: Arr
    u: Arr
    v: Arr
    energy: Arr
    sense: IntArr
    point: Arr
    length: float

    def __post_init__(self) -> None:
        count = len(self.index)
        for item in fields(self):
            if item.name == "length":
                continue
            column = np.array(getattr(self, item.name), dtype=_KINDS.get(item.name, np.float64))
            shape: tuple[int, ...] = (count,)
            if item.name == "point":
                shape = (count, 3)
                column = column.reshape(0, 3) if column.size == 0 else column
            if column.shape != shape:
                raise ValueError(f"losses take one entry an electron: {item.name} is {column.shape}, not {shape}")
            column.setflags(write=False)
            object.__setattr__(self, item.name, column)
        unknown = sorted(set(self.how.tolist()) - set(HOW))
        if unknown:
            raise ValueError(f"an electron is lost as one of {list(HOW)}, not {unknown}")
        senses = sorted(set(self.sense.tolist()) - {-1, 0, 1})
        if senses:
            raise ValueError(f"a sense along the path is +1, -1 or 0, not {senses}")
        if not self.length > 0:
            raise ValueError(f"losses take the lap's length, a positive number of metres, not {self.length}")
        object.__setattr__(self, "length", float(self.length))

    def __len__(self) -> int:
        return len(self.index)

    @classmethod
    def load(cls, record: Mapping[str, Any]) -> Losses:
        """Losses back from the plain data array16.lab keeps of them (lab.plain, through JSON): a null where a number
        was not finite comes back as nan."""
        point = np.array(
            [[math.nan if c is None else float(c) for c in row] for row in record["point"]], dtype=np.float64
        )
        return cls(
            index=np.array(record["index"], dtype=np.intp),
            how=np.array(record["how"], dtype=np.str_),
            step=np.array(record["step"], dtype=np.intp),
            time=_floats(record["time"]),
            s=_floats(record["s"]),
            laps=_floats(record["laps"]),
            u=_floats(record["u"]),
            v=_floats(record["v"]),
            energy=_floats(record["energy"]),
            sense=np.array(record["sense"], dtype=np.intp),
            point=point.reshape(-1, 3),
            length=float(record["length"]),
        )


class LossLog:
    """What a watch keeps to write its Losses: the places its last two calls were given, the velocities of the last,
    each electron's arc length at its start, and the entries so far, nothing a step beyond that. The watch calls `note`
    after every step with what it found: the electrons the field found off the grid at the start of the step just
    taken (`hit`), those another observer stopped (`others`), the laps each has gone by the watch's count, and how far
    along the path each now is (m, from its start). An entry's laps are its arc length's advance on its start's, the
    whole laps added that bring them nearest the watch's count. With `start`, where the electrons were launched (N, 3),
    the first step's crossings are found as every other's and the laps counted from the launch; without it the laps
    are counted from where the first step left the electrons, and one lost in the first step is placed where the field
    first found it off the grid, or, off it from its launch, where the first step left it."""

    def __init__(self, swept: SectionField, count: int, start: Arr | None = None) -> None:
        self.swept, self.count = swept, count
        self._last: Arr | None = None if start is None else np.array(start, dtype=np.float64).reshape(count, 3)
        self._origin: Arr | None = None if self._last is None else swept.locate(self._last).s
        self._before: Arr | None = None
        self._pace: Arr | None = None  # gamma v over the step that ended at _last
        self._entries: list[dict[str, Any]] = []

    def note(self, step: int, t: float, x: Arr, u: Arr, hit: BoolArr, others: BoolArr, laps: Arr, gone: Arr) -> None:
        """Write the entries of the electrons lost at this call of the watch, after step `step` at time t."""
        if self._origin is None:
            self._origin = self.swept.locate(x).s
        dt = t / (step + 1)
        struck, stopped = np.flatnonzero(hit), np.flatnonzero(others)
        if len(struck):
            self._entries.append(self._crossed(struck, step, dt, x, u, laps))
        if len(stopped):
            where, _ = self._place(x[stopped])
            how = np.full(len(stopped), STOPPED)
            counted = gone[stopped] / self.swept.length
            self._entries.append(
                self._entry(stopped, how, step, np.full(len(stopped), t), where, counted, x[stopped], u[stopped])
            )
        self._before, self._last = self._last, np.array(x, dtype=np.float64)
        self._pace = np.array(u, dtype=np.float64)

    def losses(self) -> Losses:
        """The entries so far, in the order the electrons were lost: by time, then by row."""
        names = ("index", "how", "step", "time", "s", "laps", "u", "v", "energy", "sense", "point")
        if not self._entries:
            empty: dict[str, Any] = {name: np.zeros((0, 3) if name == "point" else 0) for name in names}
            return Losses(**empty, length=self.swept.length)
        joined: dict[str, Any] = {name: np.concatenate([entry[name] for entry in self._entries]) for name in names}
        order = np.lexsort((joined["index"], joined["time"]))
        return Losses(**{name: joined[name][order] for name in names}, length=self.swept.length)

    # ---- where an electron crossed ----------------------------------------------------------------------------------
    def _place(self, points: Arr) -> tuple[Placement, BoolArr]:
        """Where points are along the path and in the section (SectionField.locate, which searches every sample and
        keeps nothing of it), and which of them the field would find on its grid."""
        where = self.swept.locate(points)
        _, _, on = self.swept.grid.gradient_at(where.u, where.v)
        return where, on

    def _crossing(self, a: Arr, b: Arr) -> Arr:
        """The fraction of the way from a, which the field found on the grid, to b, which it found off it, at which
        each segment first leaves the grid: 1 where the segment's ends do not bear that out (b, the field's own
        finding)."""
        count, parts = len(a), np.arange(_SCAN + 1) / _SCAN
        span = b - a
        _, on = self._place((a[:, None, :] + parts[None, :, None] * span[:, None, :]).reshape(-1, 3))
        grid_on = on.reshape(count, _SCAN + 1)
        rows = np.flatnonzero(grid_on[:, 0] & ~grid_on.all(axis=1))
        first = np.argmin(grid_on[rows], axis=1)  # the first part found off the grid
        low, high = (first - 1) / _SCAN, first / _SCAN
        for _ in range(_HALVINGS if len(rows) else 0):
            mid = (low + high) / 2
            _, inside = self._place(a[rows] + mid[:, None] * span[rows])
            low, high = np.where(inside, mid, low), np.where(inside, high, mid)
        fraction: Arr = np.ones(count)
        fraction[rows] = high
        return fraction

    def _crossed(self, idx: IntArr, step: int, dt: float, x: Arr, u: Arr, laps: Arr) -> dict[str, Any]:
        """The entries of the electrons the field found off the grid at the start of step `step`: where each crossed,
        over the step before, between the last place it was found on the grid and the first it was found off."""
        count = len(idx)
        pace = u[idx] if self._pace is None else self._pace[idx]
        if self._last is None:  # off the grid from a launch the watch was not given: where the first step left it
            point, when, lost = x[idx], np.full(count, (step + 1) * dt), 0
        elif self._before is None:  # off the grid at its launch, or lost in the first step of a launch not given
            point, when, lost = self._last[idx], np.full(count, step * dt), max(step - 1, 0)
        else:
            a, b = self._before[idx], self._last[idx]
            fraction = self._crossing(a, b)
            point = a + fraction[:, None] * (b - a)
            when, lost = (step - 1 + fraction) * dt, step - 1
        where, _ = self._place(point)
        grid = self.swept.grid
        pu, pv = where.u, where.v
        inside = (
            np.isfinite(pu)
            & np.isfinite(pv)
            & (pu >= grid.u[0])
            & (pu <= grid.u[-1])
            & (pv >= grid.v[0])
            & (pv <= grid.v[-1])
        )
        how = np.where(inside, STRUCK, LEFT)
        return self._entry(idx, how, lost, when, where, laps[idx], point, pace)

    def _entry(
        self, idx: IntArr, how: StrArr, step: int, when: Arr, where: Placement, counted: Arr, point: Arr, pace: Arr
    ) -> dict[str, Any]:
        """One call's entries: the electrons' rows, how and when each was lost, its place, its laps (its arc length's
        advance on its start's, with the whole laps that bring them nearest the watch's count, `counted`), and its
        energy and sense from gamma v over its last step (`pace`)."""
        assert self._origin is not None
        advance = (where.s - self._origin[idx]) / self.swept.length
        laps = advance + np.round(counted - advance)
        tangent = np.cross(where.first, where.second)  # second = T x first
        along = (pace * tangent).sum(axis=1)
        w = (pace * pace).sum(axis=1) / C0**2  # (gamma v / c)^2 = gamma^2 - 1
        sense = np.where(np.abs(along) > _ACROSS * np.sqrt(w) * C0, np.sign(along), 0.0)
        return {
            "index": idx,
            "how": how,
            "step": np.full(len(idx), step),
            "time": when,
            "s": where.s,
            "laps": laps,
            "u": where.u,
            "v": where.v,
            "energy": w / (np.sqrt(1 + w) + 1) * (ME * C0**2 / E),  # (gamma - 1) m c^2 / e
            "sense": sense,
            "point": point,
        }
