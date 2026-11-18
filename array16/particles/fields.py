"""Fields as functions of position and time, and ways of combining them."""

from __future__ import annotations

from collections.abc import Callable
from functools import reduce
from typing import Protocol

import numpy as np

from .base import Arr

BField = Callable[[Arr], Arr]  # positions (N, 3) in metres -> magnetic field (N, 3) in tesla


class Fields(Protocol):
    """The fields on a set of particles: (E, B), each (N, 3), at positions x (N, 3) and time t (s). Either is None
    where there is no such field, which saves evaluating and adding a zero."""

    def __call__(self, x: Arr, t: float) -> tuple[Arr | None, Arr | None]: ...


def uniform(b: Arr) -> BField:
    """The same field b (3,) everywhere."""
    return lambda x: np.broadcast_to(b, x.shape)


def scaled(field: BField, factor: float) -> BField:
    """A field multiplied by a number: a coil's field per ampere times its current."""
    return lambda x: field(x) * factor


def superpose(*fields: BField) -> BField:
    """The sum of several fields."""
    return lambda x: reduce(lambda a, b: a + b, (f(x) for f in fields))


def magnetic(field: BField) -> Fields:
    """A static magnetic field and no electric one."""
    return lambda x, t: (None, field(x))


def _sum(a: Arr | None, b: Arr | None) -> Arr | None:
    return b if a is None else a if b is None else a + b


def together(*fields: Fields) -> Fields:
    """Several sets of fields acting at once."""

    def total(x: Arr, t: float) -> tuple[Arr | None, Arr | None]:
        e: Arr | None = None
        b: Arr | None = None
        for field in fields:
            ei, bi = field(x, t)
            e, b = _sum(e, ei), _sum(b, bi)
        return e, b

    return total
