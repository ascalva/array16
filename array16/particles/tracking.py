"""Following particles through fields: a general loop with observers, and electrons round a closed guide."""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import NamedTuple, Protocol

import numpy as np

from .base import C0, ME, QM_E, Arr, BoolArr, CArr, E, IntArr
from .curves import ClosedCurve
from .fields import BField, Fields, magnetic
from .kinematics import v_of
from .push import push


@dataclass(frozen=True, eq=False)
class Launch:
    """Electrons to send round a guide: kinetic energies (eV), start offsets from the curve's first point (N, 3),
    and senses (+1 along the field, -1 against it)."""

    w_ev: Arr
    offset: Arr
    sense: Arr

    @classmethod
    def of(cls, w_ev: Arr, offset: Arr, sense: Arr | None = None) -> Launch:
        """One electron per energy; a single offset (3,) is shared by all."""
        n = len(w_ev)
        return cls(w_ev, np.broadcast_to(offset, (n, 3)).copy(), np.ones(n) if sense is None else sense)


class Observer(Protocol):
    """Something that looks at the particles after every step. It may return the ones to stop from here on."""

    def __call__(self, step: int, t: float, x: Arr, u: Arr, alive: BoolArr) -> BoolArr | None: ...


class Flight(NamedTuple):
    """Where integrate left the particles."""

    x: Arr  # positions (N, 3); a stopped particle stays where it stopped
    u: Arr  # gamma v (N, 3)
    alive: BoolArr  # not stopped by any observer
    steps: int  # steps taken


def integrate(
    x0: Arr,
    u0: Arr,
    qm: Arr,
    dt: float,
    steps: int,
    fields: Fields,
    observers: Sequence[Observer] = (),
    *,
    drag: float = 0.0,
    acc: Arr | None = None,
) -> Flight:
    """Advance N particles by up to `steps` Boris steps of dt in the given fields.

    x0 and u0 = gamma v are (N, 3), qm is each particle's charge over mass (N,). The fields are evaluated at the
    particles' positions at the start of each step, at t = step x dt. After each step every observer is called in
    turn with the step just taken, the time, the state and who is still going; the particles it returns are
    stopped, and the next observer sees that. The run ends early once every particle has stopped. drag is a
    Stokes rate (1/s) and acc a uniform acceleration (3,), as in push.
    """
    x, u = np.array(x0, dtype=np.float64), np.array(u0, dtype=np.float64)
    if x.ndim != 2 or x.shape[1] != 3 or u.shape != x.shape or qm.shape != (len(x),):
        raise ValueError(
            f"expected x0 and u0 of shape (N, 3) and qm of shape (N,): got {x.shape}, {u.shape}, {qm.shape}"
        )
    alive: BoolArr = np.ones(len(x), dtype=np.bool_)
    none = np.zeros_like(x)
    taken = 0
    for step in range(steps):
        e, b = fields(x, step * dt)
        xn, un = push(x, u, qm, dt, none if e is None else e, none if b is None else b, drag, acc)
        x, u = np.where(alive[:, None], xn, x), np.where(alive[:, None], un, u)
        taken = step + 1
        for observer in observers:
            stopped = observer(step, taken * dt, x, u, alive)
            if stopped is not None:
                alive = alive & ~stopped
        if not alive.any():
            break
    return Flight(x, u, alive, taken)


class GuideWatch:
    """Watches particles go round a closed curve: the station each is at, its laps, its largest distance from the
    curve, and its offset on the curve's axes each time it completes a lap. One further than `wall` is lost."""

    def __init__(self, curve: ClosedCurve, wall: float, count: int) -> None:
        self.curve, self.wall = curve, wall
        self.station: IntArr = np.zeros(count, dtype=np.intp)
        self.advance: Arr = np.zeros(count)
        self.reach: Arr = np.zeros(count)
        self.crossings: list[list[complex]] = [[] for _ in range(count)]
        self._laps: IntArr = np.zeros(count, dtype=np.intp)

    def __call__(self, step: int, t: float, x: Arr, u: Arr, alive: BoolArr) -> BoolArr:
        curve = self.curve
        self.station, self.advance = curve.advance(x, self.station, self.advance, iters=4)
        across, far = curve.offset(x, self.station)
        self.reach = np.where(alive, np.maximum(self.reach, far), self.reach)
        inside = far < self.wall
        laps = (np.abs(self.advance) // len(curve)).astype(np.intp)
        e1, e2 = curve.axes
        for i in np.flatnonzero((laps > self._laps) & alive & inside):
            at = self.station[i]
            self.crossings[i].append(complex(float(across[i] @ e1[at]), float(across[i] @ e2[at])))
        self._laps = np.maximum(self._laps, laps)
        return ~inside

    @property
    def laps(self) -> Arr:
        """Laps completed, with their fraction; negative for a particle going against the curve's direction."""
        return self.advance / len(self.curve)


class PathRecorder:
    """Keeps one particle's position every `every` steps, while it is going."""

    def __init__(self, index: int, every: int = 1) -> None:
        self.index, self.every = index, every
        self.points: list[Arr] = []

    def __call__(self, step: int, t: float, x: Arr, u: Arr, alive: BoolArr) -> None:
        if step % self.every == 0 and alive[self.index]:
            self.points.append(x[self.index].copy())

    @property
    def path(self) -> Arr:
        return np.array(self.points).reshape(-1, 3)


class GuideTrack(NamedTuple):
    """What track_guide measured, one entry per electron."""

    held: BoolArr  # still inside the wall at the end
    laps: Arr  # laps completed (negative: against T)
    reach: Arr  # largest distance from the curve while held (m)
    path: Arr  # (m, 3) positions of the recorded electron; empty if none was recorded
    section: tuple[CArr, ...]  # per electron, its offset on the curve's axes each time it completes a lap


def track_guide(
    curve: ClosedCurve,
    field: BField,
    bt: float,
    launch: Launch,
    wall: float,
    laps: float,
    *,
    steps_per_turn: int = 20,
    record: int | None = None,
    record_every: int = 3,
) -> GuideTrack:
    """Electrons in a static magnetic field that guides them round a closed curve, by full-orbit integration.

    Each electron starts at the curve's first point plus its offset, moving along the local field (against it if
    its sense is -1), and counts as lost beyond `wall` metres from the curve. The run lasts until the slowest has
    had time for `laps` laps. bt is the guide field strength (T); it sets the time step, 1/steps_per_turn of the
    gyration period.
    """
    count = len(launch.w_ev)
    v = np.array([v_of(float(w)) for w in launch.w_ev])
    x0 = np.tile(curve.points[0], (count, 1)) + launch.offset
    b0 = field(x0)
    along = b0 / np.linalg.norm(b0, axis=1)[:, None]
    u0 = (v / np.sqrt(1 - (v / C0) ** 2))[:, None] * along * launch.sense[:, None]
    dt = 2 * math.pi / (E / ME * bt) / steps_per_turn
    watch = GuideWatch(curve, wall, count)
    recorder = PathRecorder(record, record_every) if record is not None else None
    observers: list[Observer] = [watch] if recorder is None else [watch, recorder]
    steps = int(laps * curve.length / v.min() / dt)
    flight = integrate(x0, u0, np.full(count, QM_E), dt, steps, magnetic(field), observers)
    section = tuple(np.array(c, dtype=np.complex128) for c in watch.crossings)
    path = recorder.path if recorder is not None else np.empty((0, 3))
    return GuideTrack(flight.alive, watch.laps, watch.reach, path, section)
