"""Electrodes printed on boards, in the gapless approximation: the gaps between electrodes are small against the
distance to the field's null, and the rest of the board's plane is grounded.

An electrode at V in a grounded plane has potential V / (2 pi) times the solid angle it subtends, V / pi times the
angle across long electrodes (Oliveira and Miranda 2001; Wesenberg 2008). Across long electrodes the field is
two-dimensional and E_x - i E_y is an analytic function of z = x + iy: a null is a zero, its gradient a
derivative, and the potential of one pattern follows from a simpler one by a conformal map.

  RailPair      one board: two rails at the RF volts either side of a grounded centre strip. The null stands
                above the strip (ELE-1), and the well is at most (5 sqrt 5 - 11) / (2 pi^2) = 0.0091 of a
                four-rod well's depth at the same volts and distance (ELE-3)
  FacingBoards  two boards across a gap, each split into halves at opposite volts: the null at the centre, and a
                quarter of a four-rod well's depth (ELE-2)

depth_share is the share that rf.electrode_depth_per_watt_ev and rf.electrode_gain take, and gradient_share the
null's gradient over V / h^2, which sets the null's Mathieu q (rf.electrode_drive, rf.electrode_drive_at_q). 0.5 mm
from the null, at q = 0.5 and through a 10 kOhm resonator, 1 eV takes 16 V and 13 mW on four rods at 1.07 GHz,
50 V and 126 mW on two facing boards at 1.68 GHz, and 542 V and 14.7 W on one board at its best at 3.46 GHz.
Coordinates: x across the electrodes, y away from the board; volts are peak RF amplitudes, across the resonator
that drives the pattern (half to half, or rails to the ground plane); SI units.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import NamedTuple

import numpy as np

from .base import Arr


class Saddle(NamedTuple):
    """Where a particle leaves a well most easily: the saddle's height (m) and the field there (V/m)."""

    height: float
    field: float


@dataclass(frozen=True)
class RailPair:
    """One board in the plane y = 0: a grounded strip of width a = gap down the centre line, an RF rail of width
    b = rail on either side of it, both rails at the same volts V, and the rest of the plane grounded. The rails run
    from c1 = a / 2 to c2 = a / 2 + b either side of the centre line.

    The complex field E_x - i E_y is (i V / pi) [1/(z - c2) - 1/(z - c1) + 1/(z + c1) - 1/(z + c2)]. It vanishes
    where z^2 = -c1 c2: the null stands at h = sqrt(c1 c2) = sqrt(a (a + 2b)) / 2 above the centre line, the
    geometric mean of the rails' edge distances, with gradient 4 V b / (pi h (a + b)^2) (ELE-1). Above the null
    the field on the axis rises to a saddle at h sqrt(1 + w), w = (a + b) / h, and falls away beyond it. The
    saddle's field squared is (4 / pi^2)(w - 2) / ((w + 2) w^2) of (V / h)^2, at most (5 sqrt 5 - 11) / (2 pi^2)
    = 0.0091, at w = 1 + sqrt 5 (rail / gap = 3.68) (ELE-3).
    """

    gap: float
    rail: float

    def __post_init__(self) -> None:
        if not (self.gap > 0 and self.rail > 0):
            raise ValueError(f"a rail pair needs a positive gap and rail width, not {self.gap} m and {self.rail} m")

    @property
    def edges(self) -> tuple[float, float]:
        """The distances (m) of a rail's inner and outer edges from the centre line: c1 = a / 2, c2 = a / 2 + b."""
        return self.gap / 2, self.gap / 2 + self.rail

    @property
    def null_height(self) -> float:
        """The height (m) of the field's null above the centre line: sqrt(a (a + 2b)) / 2 = sqrt(c1 c2) (ELE-1)."""
        return 0.5 * math.sqrt(self.gap * (self.gap + 2 * self.rail))

    @property
    def distance(self) -> float:
        """The distance (m) from the null to the board, null_height: the h of the four-rod well it is compared with."""
        return self.null_height

    def _rails(self) -> tuple[tuple[float, float], tuple[float, float]]:
        c1, c2 = self.edges
        return (c1, c2), (-c2, -c1)

    def potential(self, x: Arr | float, y: Arr | float, volts: float) -> Arr:
        """The potential (V) at (x, y), y > 0: for each rail [x1, x2], V / pi times the angle it subtends,
        atan((x2 - x) / y) - atan((x1 - x) / y)."""
        xs, ys = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
        angle = np.zeros(np.broadcast_shapes(xs.shape, ys.shape))
        for x1, x2 in self._rails():
            angle += np.arctan((x2 - xs) / ys) - np.arctan((x1 - xs) / ys)
        out: Arr = volts / math.pi * angle
        return out

    def field(self, x: Arr | float, y: Arr | float, volts: float) -> Arr:
        """The field (V/m) at (x, y), y > 0, with (E_x, E_y) on the last axis: minus the gradient of potential, for
        each rail (V / pi)(y / d2 - y / d1, (x2 - x) / d2 - (x1 - x) / d1) with d = y^2 + (x_i - x)^2."""
        xs, ys = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
        shape = np.broadcast_shapes(xs.shape, ys.shape)
        ex, ey = np.zeros(shape), np.zeros(shape)
        for x1, x2 in self._rails():
            u1, u2 = x1 - xs, x2 - xs
            d1, d2 = ys * ys + u1 * u1, ys * ys + u2 * u2
            ex += ys / d2 - ys / d1
            ey += u2 / d2 - u1 / d1
        out: Arr = volts / math.pi * np.stack([ex, ey], axis=-1)
        return out

    def gradient(self, volts: float) -> float:
        """The field's gradient (V/m^2) at the null, the derivative of the complex field there:
        4 V b / (pi h (a + b)^2). Near the null |E| is the gradient times the distance from it, in every direction
        (ELE-1), and rf.mathieu_q of it is the guide's Mathieu q."""
        return 4 * abs(volts) * self.rail / (math.pi * self.null_height * (self.gap + self.rail) ** 2)

    def escape(self, volts: float) -> Saddle:
        """The saddle of |E| on the axis above the null, where a particle leaves the well most easily (ELE-3).

        On the axis the field points along y, (2 V / pi)(c2 / (y^2 + c2^2) - c1 / (y^2 + c1^2)): zero at the null,
        largest where y^2 = h (h + a + b), and falling away beyond. Its largest value is
        2 V b / (pi (a + b)(a + b + 2h)); across the axis the field grows, so the point is a saddle of |E|."""
        h, width = self.null_height, self.gap + self.rail
        height = math.sqrt(h * (h + width))
        return Saddle(height, 2 * abs(volts) * self.rail / (math.pi * width * (width + 2 * h)))

    @property
    def depth_share(self) -> float:
        """The well's depth as a fraction of a four-rod well's at the same volts and distance: the escape field
        squared over (V / h)^2. With w = (a + b) / h it is (4 / pi^2)(w - 2) / ((w + 2) w^2), largest at
        w = 1 + sqrt 5 (rail / gap = 3.68), where it is (5 sqrt 5 - 11) / (2 pi^2) = 0.0091; for rail / gap from
        0.25 to 8 it is 0.0010 to 0.0091 (ELE-3)."""
        return (self.escape(1.0).field * self.null_height) ** 2

    @property
    def gradient_share(self) -> float:
        """The null's gradient as a fraction of four rods' V / h^2 at the same volts and distance:
        4 b h / (pi (a + b)^2), 0.309 at the deepest well (rail / gap 3.68). With depth_share it sets the volts and
        the frequency a depth takes, and the null's Mathieu q (rf.electrode_drive, rf.electrode_drive_at_q)."""
        return self.gradient(1.0) * self.null_height**2


@dataclass(frozen=True)
class FacingBoards:
    """Two boards facing across a gap 2h, at y = +h and y = -h, each split down x = 0 into two halves: the top
    board at +V/2 for x < 0 and -V/2 for x > 0, the bottom board the opposite.

    The potential is -(V / pi) atan(tan(pi y / 2h) tanh(pi x / 2h)), which is -(V / pi) arg cosh(pi z / 2h): the
    imaginary part of an analytic function, so harmonic, and on the boards it takes their volts. The field is
    E_x - i E_y = -i (V / 2h) tanh(pi z / 2h). It vanishes at the centre with gradient pi V / (4 h^2); on the
    midplane it runs across the gap, (V / 2h) tanh(pi x / 2h), rising toward V / 2h without reaching it, and off
    the midplane it is stronger. So V / 2h is the barrier, and the well is a quarter of a four-rod well's depth at
    the same volts and distance h (ELE-2).
    """

    gap: float

    def __post_init__(self) -> None:
        if not self.gap > 0:
            raise ValueError(f"facing boards need a positive gap, not {self.gap} m")

    @property
    def distance(self) -> float:
        """The distance (m) from the null to either board: h = gap / 2."""
        return self.gap / 2

    def potential(self, x: Arr | float, y: Arr | float, volts: float) -> Arr:
        """The potential (V) at (x, y), |y| < h: -(V / pi) atan(tan(pi y / 2h) tanh(pi x / 2h))."""
        xs, ys = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
        k = math.pi / (2 * self.distance)
        out: Arr = -volts / math.pi * np.arctan(np.tan(k * ys) * np.tanh(k * xs))
        return out

    def field(self, x: Arr | float, y: Arr | float, volts: float) -> Arr:
        """The field (V/m) at (x, y), |y| < h, with (E_x, E_y) on the last axis: from
        E_x - i E_y = -i (V / 2h) tanh(pi z / 2h), E_x = (V / 2h) Im tanh(pi z / 2h) and E_y = (V / 2h) Re of it."""
        xs, ys = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
        t = np.tanh(math.pi / (2 * self.distance) * (xs + 1j * ys))
        out: Arr = volts / (2 * self.distance) * np.stack([t.imag, t.real], axis=-1)
        return out

    def gradient(self, volts: float) -> float:
        """The field's gradient (V/m^2) at the centre: pi V / (4 h^2) (ELE-2)."""
        return math.pi * abs(volts) / (4 * self.distance**2)

    def barrier_field(self, volts: float) -> float:
        """The field (V/m) a particle crosses to leave: V / 2h. Far from the split the field across the gap tends to
        it from below; the boards are reached only through a stronger field, (V / 2h) coth(pi |x| / 2h) (ELE-2)."""
        return abs(volts) / (2 * self.distance)

    @property
    def depth_share(self) -> float:
        """The well's depth as a fraction of a four-rod well's at the same volts and distance: the barrier field
        squared over (V / h)^2, 1/4 (ELE-2)."""
        return (self.barrier_field(1.0) * self.distance) ** 2

    @property
    def gradient_share(self) -> float:
        """The null's gradient as a fraction of four rods' V / h^2 at the same volts and distance h: pi / 4 (ELE-2).
        With depth_share it sets the volts and the frequency a depth takes, and the null's Mathieu q
        (rf.electrode_drive, rf.electrode_drive_at_q)."""
        return self.gradient(1.0) * self.distance**2
