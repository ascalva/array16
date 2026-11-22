"""Cross-sections of guides made of boards, solved by finite volumes: the line their copper makes and the well it
holds (docs/KNOTTED_GUIDE.md sections 3.3, 3.6, 5.1, 5.4 and 8.1).

Along a guide whose bends are gentle against its gap the field is two-dimensional, one quasi-TEM mode of the copper
across the path. Three kinds of section, each a geometry in metres:

  FacingStrips  two boards across a gap 2h, each split into halves at +V/2 and -V/2, opposite on the two boards:
                section B (flex, strips 8h wide split 0.3h), section A (rigid slotted boards, strips 1.5h wide either
                side of a split h), and ELE-2 in the gapless limit
  SlotGuide     one double-sided board with a routed slot, its four slot-edge strips a quadrupole: top and bottom
                copper of each side at opposite volts (section 3.3, S2-C3's line)
  Rails         one board: two rails at V either side of a grounded strip, the rest of its face grounded (R1, section
                3.6; ELE-1 and ELE-3 in the gapless limit)

A Section puts a geometry on a laminate of some thickness, with copper of some thickness, at a band. solve(section)
returns Solved: the line (C' with and without the laminate, n_eff, Z0, the laminate's share of the field energy,
Q_c by Wheeler's rule on the copper's thickness, Q_d = 1/(share tan delta), Q), the well (the null and its gradient,
the escape saddle and the depth share, the Mathieu q and depth at given volts), and the potential per volt of drive on
the solver's grid for a swept field: u across the path in the board's plane, v normal to the board, from it toward
the null.

Method. Gauss's law by finite volumes on a graded tensor grid (nodes): the spacing is a fraction of the copper's
thickness at the copper's faces and corners and grows with distance from them, so 18 um copper and 25 um polyimide are
resolved inside a box of many h. Each section is solved on the quarter or half its symmetry allows; the linear system
is block-tridiagonal and is solved directly by block elimination (laplace2d), exact to rounding and six to eight times
faster than the study's conjugate gradients, which stay available. The copper's loss is Wheeler's
incremental-inductance rule (wheeler_qc): the air problem's capacitance with every copper face receded by an eighth of
its thickness gives dL'/dn, so R' = (R_s / mu0) dL'/dn, and Q_c = omega L' / R'. The faces against the laminate, the
copper's treated side, take their own surface resistance (Section.rough: Hammerstad's 1 + (2/pi) atan(1.4 (D/delta)^2)
for an rms roughness D, or a measured factor), from a second air solve with those faces alone receded.
stripline_wheeler holds the receded solve to stripline's closed form (Pozar's eq. 3.198), which is itself Wheeler's
rule applied to a closed-form Z0(t): it shows the numerical derivative right, so the copper's loss here is held to
Wheeler's closed form only. The check that does not use the rule is the full-wave script's face integral of the static
field, whose limit lies 3.5 to 5 % under the rule on R1's backed board (experiments/twin_openems_line); a
surface-impedance solve (Palace) is DR-0's item for it. The charge rule on copper without thickness (the surface
charge's square as the current) does not converge with the cell and is not used. The null's gradient is a
least-squares fit of the harmonic polynomials of the section's symmetry over a quarter of h round it; the escape is a
minimax flood from the null (escape_flood): the lowest |E| a particle must cross to reach copper, the board, or the
open space beyond the guide.

Held to: the gapless boards of electrodes.py (FacingBoards, RailPair: ELE-1 to ELE-3); section B at 915 MHz, n_eff
1.031, the polyimide's share 0.074, Z0 41.7 ohm, Q_c 228, Q 214 (the knotted study's line solve); the slot guide's Q
43 on FR-4 and 117 on RO4003C; the split as a lever, Q 210 at 0.3h and 290 at h; the share against d ln C'/d ln eps_r
(DR-0). SI units; volts are peak RF amplitudes across the resonator that drives the pattern (half to half, top to
bottom copper of one side, rails to ground).
"""

from __future__ import annotations

import heapq
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, replace
from typing import Literal

import numpy as np

from array16.particles import C0, EPS0, MU0, QM_E, Arr, BoolArr, SectionGrid
from array16.particles.rf import electrode_drive_at_q, mathieu_q, ponderomotive_ev

__all__ = [
    "COARSE",
    "COPPER_SIGMA",
    "FLEX_COPPER",
    "FR4",
    "POLYIMIDE",
    "RIGID_COPPER",
    "RO4003C",
    "Calibration",
    "FacingStrips",
    "Grid",
    "Laminate",
    "Rails",
    "Section",
    "SlotGuide",
    "Solved",
    "edge_weights",
    "energy2d",
    "escape_flood",
    "fit_well",
    "laplace2d",
    "nodes",
    "section_a",
    "section_b",
    "slot_guide",
    "solve",
    "stripline_wheeler",
    "wheeler_qc",
]

COPPER_SIGMA = 5.8e7  # S/m, bare copper
FLEX_COPPER = 18e-6  # m: the flex fab's copper
RIGID_COPPER = 35e-6  # m: 1 oz, the rigid boards' copper
DIRECT_LIMIT = 50_000_000  # doubles of inverses the block elimination keeps (400 MB); past it, it keeps checkpoints
_TOL = 1e-12  # coordinates (in h) closer than this are the same place
_WINDOW = 0.25  # h: the neighbourhood of the null its gradient is fitted over


@dataclass(frozen=True)
class Laminate:
    """A board's dielectric at the band: relative permittivity er and loss tangent tand."""

    name: str
    er: float
    tand: float

    def __post_init__(self) -> None:
        if not (self.er >= 1.0 and self.tand >= 0.0):
            raise ValueError(f"a laminate takes er >= 1 and tan delta >= 0, not {self.er} and {self.tand}")


POLYIMIDE = Laminate("polyimide", 3.3, 0.004)
FR4 = Laminate("FR-4", 4.4, 0.02)
RO4003C = Laminate("RO4003C", 3.55, 0.0027)


@dataclass(frozen=True)
class Grid:
    """The graded grid: at the copper's faces and corners (and at the null, with `null`) nodes `fine` of the copper's
    thickness apart, the spacing growing by `growth` per h of distance from them up to `coarse` of h; past `reach` h
    from every feature it grows by `stretch` more per h, so distant walls cost few nodes. Wheeler's rule recedes the
    copper by `recede` of its thickness. The default is the knotted study's converged grid; COARSE its quick grid.
    With `null` the grid is refined where the null sits (and, on rails, where the saddle sits): the well's numbers
    then hold to the grid's accuracy; without it the line's grid is the study's line grid."""

    fine: float = 1 / 16
    coarse: float = 0.05
    growth: float = 0.1
    recede: float = 1 / 8
    stretch: float = 0.25
    reach: float = 8.0
    null: bool = True

    def __post_init__(self) -> None:
        if not (self.fine > 0 and self.coarse > 0 and self.growth > 0 and 0 < self.recede < 0.5):
            raise ValueError(
                f"a grid takes fine, coarse and growth above 0 and a recession between 0 and half the copper, not "
                f"{self.fine}, {self.coarse}, {self.growth} and {self.recede}"
            )
        if not (self.stretch >= 0 and self.reach > 0):
            raise ValueError(f"a grid takes stretch >= 0 and reach > 0, not {self.stretch} and {self.reach}")


COARSE = Grid(fine=1 / 8, coarse=0.1)


@dataclass(frozen=True)
class _Box:
    """A rectangle in units of h: copper at `volts` per volt of drive (receding in Wheeler's rule when `lossy`), or
    laminate cells, or where a particle meets the board."""

    x0: float
    x1: float
    y0: float
    y1: float
    volts: float = 0.0
    lossy: bool = True


Mirror = Literal["odd", "even", "none"]


@dataclass(frozen=True)
class _Layout:
    """A section as the solver sees it, lengths in units of h: the domain [0, xmax] x [ylo, ymax], odd or even about
    x = 0 and about y = 0 (or the whole of y), the copper, the laminate's cells, the boards' bodies and the flight box
    beyond which a particle has escaped, the grid's features and finest spacings, v0 (v = (y - v0) h), and where the
    null is looked for: the origin, or on the axis above `above`."""

    h: float
    xmax: float
    ylo: float
    ymax: float
    odd_x: bool
    mirror_y: Mirror
    copper: tuple[_Box, ...]
    laminate: tuple[_Box, ...]
    solid: tuple[_Box, ...]
    flight: tuple[float, float]
    xfeat: tuple[float, ...]
    yfeat: tuple[float, ...]
    fine_x: float
    fine_y: float
    v0: float
    on_axis: bool
    above: float = 0.0


@dataclass(frozen=True)
class FacingStrips:
    """Two boards facing across a gap (m) between their copper faces, 2h, each with copper `width` wide (outer edge
    to outer edge) split down the centre line by `split`: the top board's halves at +V/2 for x < 0 and -V/2 for
    x > 0, the bottom board's the opposite. The laminate lies behind the copper, `margin` beyond its outer edges
    (math.inf: the board out to the walls), routed away under the split when `slotted`. Grounded walls stand `walls` h
    beyond the copper's outer edges and the laminate's back. Section B: width 8h, split 0.3h, margin 2h; section A:
    width 4h, split h, slotted. In the gapless limit (split and copper to zero, width to infinity) it is ELE-2's
    FacingBoards: gradient pi V/(4 h^2), a quarter of a four-rod well."""

    gap: float
    width: float
    split: float
    margin: float = math.inf
    slotted: bool = False
    walls: float = 8.0

    def __post_init__(self) -> None:
        if not (self.gap > 0 and 0 < self.split < self.width):
            raise ValueError(
                f"facing strips take a positive gap and a split narrower than the strips, not a gap of {self.gap} m "
                f"and a split of {self.split} m on {self.width} m"
            )
        if not (self.margin >= 0 and self.walls > 2):
            raise ValueError(
                f"facing strips take a margin >= 0 and walls more than 2h beyond the copper, not {self.margin} m and "
                f"{self.walls} h"
            )

    @property
    def distance(self) -> float:
        """h: the distance (m) from the null to either board's copper."""
        return self.gap / 2

    def scaled(self, factor: float) -> FacingStrips:
        """The same pattern `factor` times larger: gap, width, split and margin (the copper and the laminate keep their
        thickness)."""
        return replace(
            self,
            gap=self.gap * factor,
            width=self.width * factor,
            split=self.split * factor,
            margin=self.margin * factor,
        )

    def _layout(self, section: Section, grid: Grid) -> _Layout:
        h = self.gap / 2
        tc = section.copper / h
        tp = section.thickness / h if section.laminate is not None else 0.0
        half, inner, wide = self.width / (2 * h), self.split / (2 * h), self.margin / h
        rise = 1.0
        back, rc = rise + tc, grid.recede * tc
        xmax, ymax = half + self.walls, back + tp + self.walls
        xfeat = [inner, half, inner + rc, half - rc] + ([half + wide] if tp > 0 and math.isfinite(wide) else [])
        yfeat = [rise, back, rise + rc, back - rc] + ([back + tp] if tp > 0 else [])
        if grid.null:
            xfeat.append(0.0)
            yfeat.append(0.0)
        laminate = (_Box(inner if self.slotted else 0.0, half + wide, back, back + tp),) if tp > 0 else ()
        return _Layout(
            h=h,
            xmax=xmax,
            ylo=0.0,
            ymax=ymax,
            odd_x=True,
            mirror_y="odd",
            copper=(_Box(inner, half, rise, back, volts=-0.5),),
            laminate=laminate,
            solid=laminate,
            flight=(half + 2.0, back),
            xfeat=tuple(xfeat),
            yfeat=tuple(yfeat),
            fine_x=grid.fine * tc,
            fine_y=min(grid.fine * tc, tp / 3) if tp > 0 else grid.fine * tc,
            v0=-rise,
            on_axis=False,
        )


@dataclass(frozen=True)
class SlotGuide:
    """One double-sided board with a routed slot `slot` wide (m), its four slot-edge strips `strip` wide, their copper
    set back `pull` from the slot's wall, the top and bottom strips of each side at opposite volts, +V/2 and -V/2,
    and the other side the opposite: a quadrupole whose null sits in the slot on the board's mid-plane, at h = slot/2
    from either wall. Grounded planes stand `planes` beyond either face (another board's copper, or a shield) and
    grounded walls `wall` from the centre line (moved out to 4h beyond the strips where they would be nearer). The
    walls are ideal. With `planes_lose` the planes are copper of the section's thickness, receded in Wheeler's rule as
    the strips are, so Q_c counts their loss; without it they are ideal and Q_c is the strips' alone, the bound from
    above (on S2-C3's section, 1 mm off, the planes take 1.8 % of the loss: Q_c 154.4 against 157.2, the quadrupole's
    field held between the strips and in the board). The board is the section's laminate, its thickness the section's
    thickness; a particle that meets the slot's wall or leaves it half way to the planes is lost."""

    slot: float
    strip: float
    pull: float = 0.2e-3
    planes: float = 1.0e-3
    wall: float = 6.0e-3
    planes_lose: bool = False

    def __post_init__(self) -> None:
        if not (self.slot > 0 and self.strip > 0 and self.pull >= 0 and self.planes > 0 and self.wall > 0):
            raise ValueError(
                f"a slot guide takes a positive slot, strip, planes and wall and a set-back >= 0, not {self.slot}, "
                f"{self.strip}, {self.planes}, {self.wall} and {self.pull} m"
            )

    @property
    def distance(self) -> float:
        """h: the distance (m) from the null to the slot's wall, slot / 2."""
        return self.slot / 2

    def scaled(self, factor: float) -> SlotGuide:
        """The same pattern `factor` times larger: slot, strips, set-back, planes and walls (the board and the copper
        keep their thickness)."""
        return replace(
            self,
            slot=self.slot * factor,
            strip=self.strip * factor,
            pull=self.pull * factor,
            planes=self.planes * factor,
            wall=self.wall * factor,
        )

    def _layout(self, section: Section, grid: Grid) -> _Layout:
        h = self.slot / 2
        tc = section.copper / h
        if not section.thickness > 0:
            raise ValueError("a slot guide's board takes a thickness: the section's thickness, not 0")
        if not self.planes > section.copper:
            raise ValueError(
                f"a slot guide's planes stand beyond its copper: {self.planes} m planes against "
                f"{section.copper} m of copper"
            )
        bz, rc = section.thickness / (2 * h), grid.recede * tc
        x0 = 1.0 + self.pull / h
        x1 = x0 + self.strip / h
        xmax, top = max(self.wall / h, x1 + 4.0), bz + self.planes / h
        board = _Box(1.0, xmax, 0.0, bz)
        xfeat = [1.0, x0, x1, x0 + rc, x1 - rc] + ([0.0] if grid.null else [])
        yfeat = [bz, bz + tc, bz + rc, bz + tc - rc] + ([0.0] if grid.null else [])
        copper = [_Box(x0, x1, bz, bz + tc, volts=0.5)]
        ymax = top
        if self.planes_lose:  # the planes as copper boxes, their faces toward the board receded with the strips'
            ymax = top + tc
            copper.append(_Box(0.0, xmax, top, ymax))
            yfeat += [top, top + rc]
        return _Layout(
            h=h,
            xmax=xmax,
            ylo=0.0,
            ymax=ymax,
            odd_x=True,
            mirror_y="odd",
            copper=tuple(copper),
            laminate=(board,) if section.laminate is not None else (),
            solid=(board,),
            flight=(math.inf, bz + self.planes / (2 * h)),
            xfeat=tuple(xfeat),
            yfeat=tuple(yfeat),
            fine_x=grid.fine * tc,
            fine_y=min(grid.fine * tc, 2 * bz / 3),
            v0=-bz,
            on_axis=False,
        )


@dataclass(frozen=True)
class Rails:
    """One board: a grounded centre strip `gap` wide (m; RailPair's a), a rail `rail` wide (b) at V either side of it,
    and grounded copper beyond, out to the walls, the rails isolated from the grounded copper by `space`, centred on
    the rails' edges c1 = a/2 and c2 = a/2 + b; a ground plane on the board's back face when `backed`. The null
    stands about h = sqrt(c1 c2) above the copper (ELE-1), and the well escapes over the saddle about h sqrt(1 + w)
    up, w = (a + b)/h (ELE-3). Grounded walls stand `walls` null heights beyond the rails, above and below the board.
    In the gapless limit (space and copper to zero, no laminate) it is ELE-1's RailPair(gap, rail)."""

    gap: float
    rail: float
    space: float
    backed: bool = False
    walls: float = 100.0

    def __post_init__(self) -> None:
        if not (self.gap > 0 and self.rail > 0 and 0 < self.space < min(self.gap, self.rail)):
            raise ValueError(
                f"rails take a positive strip and rails and a space narrower than both, not a {self.gap} m strip, "
                f"{self.rail} m rails and {self.space} m spaces"
            )
        if not self.walls > 4:
            raise ValueError(f"rails take walls more than 4 null heights out, not {self.walls}")

    @property
    def distance(self) -> float:
        """The gapless null height sqrt(c1 c2) (m), ELE-1's; Solved.distance is the solved one."""
        return math.sqrt(self.gap / 2 * (self.gap / 2 + self.rail))

    def scaled(self, factor: float) -> Rails:
        """The same pattern `factor` times larger: strip, rails and spaces (the board and the copper keep their
        thickness)."""
        return replace(self, gap=self.gap * factor, rail=self.rail * factor, space=self.space * factor)

    def _layout(self, section: Section, grid: Grid) -> _Layout:
        h = self.distance
        c1, c2 = self.gap / (2 * h), (self.gap / 2 + self.rail) / h
        s, tc, tb = self.space / (2 * h), section.copper / h, section.thickness / h
        rc = grid.recede * tc
        xmax, ymax = c2 + self.walls, tc + self.walls
        ylo = -tb - tc if self.backed else -tb - self.walls
        copper = [
            _Box(0.0, c1 - s, 0.0, tc),
            _Box(c1 + s, c2 - s, 0.0, tc, volts=1.0),
            _Box(c2 + s, xmax, 0.0, tc),
        ]
        if self.backed:
            copper.append(_Box(0.0, xmax, ylo, -tb))
        xfeat = [c1 - s, c1 - s - rc, c1 + s, c1 + s + rc, c2 - s, c2 - s - rc, c2 + s, c2 + s + rc]
        yfeat = [0.0, tc, rc, tc - rc] + ([-tb] if tb > 0 else []) + ([-tb - rc] if self.backed else [])
        if grid.null:
            xfeat.append(0.0)
            yfeat += [tc + 1.0, tc + math.sqrt(1.0 + c1 + c2)]
        has_board = section.laminate is not None and tb > 0
        return _Layout(
            h=h,
            xmax=xmax,
            ylo=ylo,
            ymax=ymax,
            odd_x=False,
            mirror_y="none",
            copper=tuple(copper),
            laminate=(_Box(0.0, xmax, -tb, 0.0),) if has_board else (),
            solid=(_Box(0.0, xmax, ylo, 0.0),),
            flight=(math.inf, math.inf),
            xfeat=tuple(xfeat),
            yfeat=tuple(yfeat),
            fine_x=grid.fine * tc,
            fine_y=min(grid.fine * tc, tb / 3) if has_board else grid.fine * tc,
            v0=tc,
            on_axis=True,
            above=tc,
        )


Geometry = FacingStrips | SlotGuide | Rails


@dataclass(frozen=True)
class Section:
    """A cross-section to solve: a geometry, the laminate (None: the copper alone, in vacuum) and its thickness (m)
    behind the copper (facing strips), as the board (slot guide) or under the rails, the copper's thickness (m), the
    band f (Hz), the surface resistance over bare copper's (ENIG: 1.3 to 5 at 915 MHz), the grid, and the roughness of
    the copper's faces against the laminate, the treated side, as a further factor on their surface resistance
    (Hammerstad's 1 + (2/pi) atan(1.4 (D/delta)^2) for an rms roughness D against the skin depth delta, 2.2 um at 915
    MHz: 1.37 at 1.5 um, 1.55 at 2 um; or a measured factor)."""

    geometry: Geometry
    laminate: Laminate | None = None
    thickness: float = 0.0
    copper: float = FLEX_COPPER
    f: float = 915e6
    resistance: float = 1.0
    grid: Grid = Grid()
    rough: float = 1.0

    def __post_init__(self) -> None:
        if not (self.copper > 0 and self.f > 0 and self.resistance > 0 and self.rough > 0 and self.thickness >= 0):
            raise ValueError(
                f"a section takes copper, a band, a surface resistance and a roughness above 0 and a thickness >= 0, "
                f"not {self.copper} m, {self.f} Hz, {self.resistance}, {self.rough} and {self.thickness} m"
            )
        if self.laminate is not None and not self.thickness > 0:
            raise ValueError(f"a laminate takes a thickness, not {self.thickness} m of {self.laminate.name}")

    def scaled(self, factor: float) -> Section:
        """The same section with its pattern `factor` times larger: the copper, the laminate and the band kept."""
        return replace(self, geometry=self.geometry.scaled(factor))


@dataclass(frozen=True, eq=False)
class Solved:
    """A solved cross-section.

    The line: capacitance C' (F/m) with the laminate and air_capacitance without it, so n_eff = sqrt(C'/C'_air) and
    Z0 = 1/(c sqrt(C' C'_air)); share, the laminate's part of the field energy; qc by Wheeler's rule on the copper,
    qd = 1/(share tan delta), and quality, the line's Q, their parallel sum. The well: the null (u, v) (m) and the
    field's gradient there per volt of drive (1/m^2: |E| is gradient x volts x the distance from the null); the saddle
    (u, v) (m) where a particle leaves most easily and the field there per volt (escape, 1/m); distance, the h the
    shares are counted at: the null to the copper for facing strips and rails (on rails the solved null's height above
    the copper), to the slot's wall for the slot guide. The field: phi (len(u), len(v)), the potential per volt of
    drive at the nodes u (across the path in the board's plane, from the centre line) and v (normal to the board, from
    it toward the null): from the lower board's copper face for facing strips (the null at v = h), from the board's
    lower face for the slot guide (the null on the mid-plane), from the rails' copper face for rails."""

    section: Section
    capacitance: float
    air_capacitance: float
    share: float
    qc: float
    qd: float
    null: tuple[float, float]
    gradient: float
    saddle: tuple[float, float]
    escape: float
    distance: float
    u: Arr
    v: Arr
    phi: Arr

    @property
    def index(self) -> float:
        """n_eff = sqrt(C'/C'_air): the guided wave's index."""
        return math.sqrt(self.capacitance / self.air_capacitance)

    @property
    def z0(self) -> float:
        """Z0 (ohm) = 1/(c sqrt(C' C'_air))."""
        return 1 / (C0 * math.sqrt(self.capacitance * self.air_capacitance))

    @property
    def quality(self) -> float:
        """The line's Q: 1/(1/Q_c + 1/Q_d)."""
        return 1 / (1 / self.qc + 1 / self.qd)

    @property
    def wavelength(self) -> float:
        """The guided wavelength (m) at the band, c/(n_eff f): a ring of m guided waves is m times it round."""
        return C0 / (self.index * self.section.f)

    @property
    def depth_share(self) -> float:
        """The well's depth as a share of a four-rod well's at the same volts and distance h: (escape h)^2, as
        electrodes.FacingBoards and electrodes.RailPair count it."""
        return (self.escape * self.distance) ** 2

    @property
    def gradient_share(self) -> float:
        """The null's gradient over four rods' V/h^2 at the same volts and distance: gradient h^2."""
        return self.gradient * self.distance**2

    def shunt(self, m: int = 1) -> float:
        """The shunt resistance (ohm) of a ring of m guided waves of this line carrying a travelling wave: Z0 Q/(2 pi m)
        (section 5.2)."""
        return self.z0 * self.quality / (2 * math.pi * m)

    def mathieu_q(self, volts: float, qm: float = QM_E) -> float:
        """Mathieu's q at the null for `volts` of drive at the section's band."""
        return mathieu_q(self.gradient * abs(volts), self.section.f, qm)

    def depth_ev(self, volts: float, qm: float = QM_E) -> float:
        """The pseudo-potential well's depth (eV) at `volts` of drive at the section's band: the ponderomotive energy of
        the escape saddle's field."""
        return ponderomotive_ev(self.escape * abs(volts), self.section.f, qm)

    def volts_for(self, depth_ev: float, q: float) -> float:
        """The drive (V) a well depth_ev deep takes at Mathieu's q: 8 gamma depth/(share q) (RFG-3)."""
        return 8 * self.gradient_share * depth_ev / (self.depth_share * q)

    def power_for(self, depth_ev: float, q: float, m: int = 1) -> float:
        """The power (W) a ring of m guided waves takes for a well depth_ev deep at Mathieu's q: V^2/(2 R)."""
        return self.volts_for(depth_ev, q) ** 2 / (2 * self.shunt(m))

    def section_grid(self, reach: float | None = None) -> SectionGrid:
        """The potential per volt of drive as particles.sections.SectionGrid takes it, for SectionField to sweep along
        a path: the same nodes u (across the path) and v (from the board toward the null), the cells of copper and of
        the boards solid, and the null where the grid's own interpolated field vanishes (SectionGrid.find_null from
        the solved null, which it moves by the interpolation's error alone: on the axis for rails, onto a node for the
        odd sections). With `reach` (m) the grid keeps the nodes within reach of the null along u and along v, so that
        a swept grid stays inside its curve's thickness, where SectionField's window search is trusted."""
        section = self.section
        lay = section.geometry._layout(section, section.grid)
        u, v, phi = self.u, self.v, self.phi
        if reach is not None:
            near_u = np.abs(u - self.null[0]) <= reach * (1 + _TOL)
            near_v = np.abs(v - self.null[1]) <= reach * (1 + _TOL)
            u, v, phi = u[near_u], v[near_v], phi[np.ix_(near_u, near_v)]
        # each cell's centre in the solver's units of h, folded onto the part it solved (x >= 0, and y >= 0 when
        # mirrored), where the layout's copper and boards are drawn
        cx = np.abs((u[1:] + u[:-1]) / 2) / lay.h
        cy = (v[1:] + v[:-1]) / 2 / lay.h + lay.v0
        if lay.mirror_y != "none":
            cy = np.abs(cy)
        gx, gy = np.meshgrid(cx, cy, indexing="ij")
        solid: BoolArr = np.zeros(gx.shape, dtype=np.bool_)
        for box in (*lay.copper, *lay.solid):
            solid |= _inside(box, gx, gy, closed=True)
        first = SectionGrid(u, v, phi, self.null, solid)
        return SectionGrid(first.u, first.v, first.phi, first.find_null(), solid)


@dataclass(frozen=True)
class Calibration:
    """Wheeler's rule on a stripline against the closed form: Z0 (ohm), the conductor's attenuation by the rule and by
    the closed form (Np/m), and their ratio."""

    z0: float
    alpha: float
    closed_form: float
    ratio: float


# ---- the grid and the solver --------------------------------------------------------------------------------------
def nodes(
    lo: float,
    hi: float,
    features: Sequence[float],
    fine: float,
    coarse: float,
    growth: float = 0.3,
    stretch: float = 0.0,
    reach: float = math.inf,
) -> Arr:
    """Grid nodes from lo to hi, `fine` apart at the features, the spacing growing by `growth` per unit distance from
    the nearest feature up to `coarse`, and past `reach` from every feature by `stretch` more per unit distance. The
    ends and the features are nodes (a feature within fine/4 of an end is that end); the march's nearest node moves onto
    each feature, and nodes closer than fine/4 to the one before are dropped unless they are features or the end."""
    snapped = [lo if f - lo < fine / 4 else hi if hi - f < fine / 4 else f for f in features if lo <= f <= hi]
    feats = np.unique(np.array(snapped, dtype=np.float64))
    xs = [lo]
    while xs[-1] < hi:
        if len(feats):
            d = float(np.min(np.abs(feats - xs[-1])))
            step = min(coarse, fine + growth * d) + stretch * max(0.0, d - reach)
        else:
            step = coarse
        xs.append(xs[-1] + step)
    arr = np.array(xs, dtype=np.float64)
    arr[-1] = hi
    for f in feats:
        k = int(np.argmin(np.abs(arr - f)))
        if 0 < k < len(arr) - 1:
            arr[k] = f
    arr = np.union1d(arr, feats)
    pinned = np.isin(arr, feats) | (arr == hi)
    keep = np.concatenate([[True], np.diff(arr) > fine / 4]) | pinned
    last = int(np.flatnonzero(keep[:-1])[-1])
    if last > 0 and hi - arr[last] <= fine / 4 and not pinned[last]:
        keep[last] = False  # the end stays a node, not the march's node just before it
    out: Arr = arr[keep]
    return out


def edge_weights(x: Arr, y: Arr, eps: Arr) -> tuple[Arr, Arr]:
    """Finite-volume weights of a tensor grid's edges: the permittivity of the cells beside an edge times their share
    of its control face, over the edge's length. eps (nx - 1, ny - 1) holds each cell's relative permittivity."""
    hx, hy = np.diff(x), np.diff(y)
    ep = np.pad(eps, ((0, 0), (1, 1)))
    hp = np.pad(hy, (1, 1))
    wx: Arr = (ep[:, :-1] * hp[None, :-1] / 2 + ep[:, 1:] * hp[None, 1:] / 2) / hx[:, None]
    epx = np.pad(eps, ((1, 1), (0, 0)))
    hpx = np.pad(hx, (1, 1))
    wy: Arr = (epx[:-1, :] * hpx[:-1, None] / 2 + epx[1:, :] * hpx[1:, None] / 2) / hy[None, :]
    return wx, wy


def _apply(wx: Arr, wy: Arr, u: Arr) -> Arr:
    """The grid's Laplacian times u: at each node, the sum over its edges of weight x (u here - u there)."""
    r = np.zeros_like(u)
    fx = wx * (u[1:, :] - u[:-1, :])
    r[:-1, :] -= fx
    r[1:, :] += fx
    fy = wy * (u[:, 1:] - u[:, :-1])
    r[:, :-1] -= fy
    r[:, 1:] += fy
    return r


def _pcg(apply: Callable[[Arr], Arr], b: Arr, diag: Arr, free: BoolArr, tol: float, iters: int) -> Arr:
    """Jacobi-preconditioned conjugate gradients for A x = b on the free nodes."""
    x = np.zeros_like(b)
    r = np.where(free, b, 0.0)
    z = r / diag
    p = z.copy()
    rz = float((r * z).sum())
    norm_b = math.sqrt(float((r * r).sum())) + 1e-300
    for _ in range(iters):
        ap = np.where(free, apply(p), 0.0)
        alpha = rz / (float((p * ap).sum()) + 1e-300)
        x += alpha * p
        r -= alpha * ap
        if math.sqrt(float((r * r).sum())) < tol * norm_b:
            break
        z = r / diag
        rz_new = float((r * z).sum())
        p = z + (rz_new / (rz + 1e-300)) * p
        rz = rz_new
    return x


def _blocks(wx: Arr, wy: Arr, free: BoolArr, rhs: Arr, limit: int) -> Arr:
    """A x = rhs on the free nodes (x = 0 on the others) by block elimination: the grid's Laplacian is
    block-tridiagonal, a block per line of nodes along the first axis (the longer, so the blocks are the shorter
    lines), and the Schur complements S_i = D_i - C S_(i-1)^-1 C carry the elimination along it; back-substitution
    returns x. Exact to rounding. The inverses take n m^2 numbers for n lines of m nodes; past `limit` only every
    sqrt(n)-th is kept and the others are formed again during the back-substitution, twice the work in a fraction of
    the memory."""
    n, m = free.shape
    if m > n:
        out: Arr = _blocks(wy.T, wx.T, free.T, rhs.T, limit).T
        return out
    cx = np.where(free[:-1] & free[1:], wx, 0.0)
    cy = np.where(free[:, :-1] & free[:, 1:], wy, 0.0)
    deg = np.zeros((n, m))
    deg[:-1, :] += wx
    deg[1:, :] += wx
    deg[:, :-1] += wy
    deg[:, 1:] += wy
    deg = np.where(free, deg, 1.0)
    b = np.where(free, rhs, 0.0)
    k = np.arange(m)

    def schur(i: int, before: Arr | None) -> Arr:
        """S_i^-1 from S_(i-1)^-1."""
        s = np.zeros((m, m))
        s[k, k] = deg[i]
        s[k[:-1], k[1:]] = -cy[i]
        s[k[1:], k[:-1]] = -cy[i]
        if before is not None:
            s -= before * np.outer(cx[i - 1], cx[i - 1])
        inv: Arr = np.linalg.inv(s)
        return inv

    fits = n * m * m <= limit
    span = n if fits else max(1, math.isqrt(n))
    kept: dict[int, Arr] = {}
    g = np.empty((n, m))
    inv: Arr | None = None
    for i in range(n):
        g[i] = b[i] if inv is None else b[i] + cx[i - 1] * (inv @ g[i - 1])
        inv = schur(i, inv)
        if fits or i % span == 0:
            kept[i] = inv
    x = np.empty((n, m))
    for start in reversed(range(0, n, span)):
        stop = min(start + span, n)
        run = [kept[start]]
        for i in range(start + 1, stop):
            run.append(kept[i] if fits else schur(i, run[-1]))
        for i in reversed(range(start, stop)):
            x[i] = run[i - start] @ (g[i] if i == n - 1 else g[i] + cx[i] * x[i + 1])
    solved: Arr = np.where(free, x, 0.0)
    return solved


def laplace2d(
    x: Arr,
    y: Arr,
    eps: Arr,
    fixed: BoolArr,
    value: Arr,
    tol: float = 1e-12,
    iters: int = 60000,
    direct: bool = True,
    limit: int | None = None,
) -> Arr:
    """The potential on a tensor grid's nodes (nx, ny) with cells of relative permittivity eps (nx - 1, ny - 1):
    value where fixed, Gauss's law by finite volumes elsewhere, zero normal flux on the box's faces where they are not
    fixed. Solved directly by block elimination (`direct`, the default), keeping at most `limit` numbers of inverses
    (DIRECT_LIMIT), or by Jacobi-preconditioned conjugate gradients to a residual `tol` of the right-hand side's, as the
    knotted study solved it."""
    wx, wy = edge_weights(x, y, eps)
    free = ~fixed
    u0 = np.where(fixed, value, 0.0)
    rhs = -_apply(wx, wy, u0)
    if direct:
        xs = _blocks(wx, wy, free, rhs, DIRECT_LIMIT if limit is None else limit)
    else:
        diag = np.zeros_like(value)
        diag[:-1, :] += wx
        diag[1:, :] += wx
        diag[:, :-1] += wy
        diag[:, 1:] += wy
        xs = _pcg(lambda u: _apply(wx, wy, u), rhs, np.where(diag > 0, diag, 1.0), free, tol, iters)
    out: Arr = np.where(fixed, value, xs)
    return out


def energy2d(x: Arr, y: Arr, eps: Arr, phi: Arr) -> float:
    """The integral of eps |grad phi|^2 over the grid, by its edges: twice the energy per length in eps0 units."""
    wx, wy = edge_weights(x, y, eps)
    return float((wx * np.diff(phi, axis=0) ** 2).sum() + (wy * np.diff(phi, axis=1) ** 2).sum())


def escape_flood(emag: Arr, start: tuple[int, int], exits: BoolArr) -> tuple[float, tuple[int, int]]:
    """The lowest |E| a particle must cross to leave a well: a minimax flood from the null (Dijkstra on the largest
    field met), stopping at the first exit. Returns that field and the node where the path met it."""
    nx, ny = emag.shape
    e = emag.ravel().tolist()
    out = exits.ravel().tolist()
    seen = bytearray(nx * ny)
    first = start[0] * ny + start[1]
    heap = [(e[first], first)]
    level, at = 0.0, first
    while heap:
        val, k = heapq.heappop(heap)
        if seen[k]:
            continue
        seen[k] = 1
        if val > level:
            level, at = val, k
        if out[k]:
            break
        i, j = divmod(k, ny)
        for step, ok in ((ny, i + 1 < nx), (-ny, i > 0), (1, j + 1 < ny), (-1, j > 0)):
            if ok and not seen[k + step]:
                heapq.heappush(heap, (e[k + step], k + step))
    i, j = divmod(at, ny)
    return level, (i, j)


def wheeler_qc(
    c0: float,
    c1: float,
    recede: float,
    f: float,
    resistance: float = 1.0,
    rough: float = 1.0,
    c1_rough: float | None = None,
) -> float:
    """Q_c by Wheeler's incremental-inductance rule: the air problem's capacitance C0 and C1 with every copper face
    receded by `recede` (m), L' = mu0 eps0 / C', R' = (R_s / mu0) dL'/dn, Q_c = 2 recede / (delta (C0/C1 - 1)), for
    surface resistance `resistance` times bare copper's. With C1' (`c1_rough`), the capacitance with only the rough
    faces receded, those faces' part of dL'/dn, C0/C1' - 1, takes the further factor `rough`."""
    delta = math.sqrt(2 / (2 * math.pi * f * MU0 * COPPER_SIGMA))
    excess = c0 / c1 - 1
    if c1_rough is not None:
        excess += (rough - 1) * (c0 / c1_rough - 1)
    return 2 * recede / (delta * excess) / resistance


# ---- the section solve --------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class _Line:
    """A layout's line: the grid, C' with and without the laminate (eps0 per volt^2), the receded copper's, the share,
    the potential the well is read from, and C' with only the faces against the laminate receded (None when they
    take no further factor)."""

    x: Arr
    y: Arr
    c: float
    c0: float
    c1: float
    share: float
    phi: Arr
    c1_rough: float | None = None


def _grid(lay: _Layout, grid: Grid) -> tuple[Arr, Arr]:
    """The graded grid's nodes for a layout."""
    spread = {"coarse": grid.coarse, "growth": grid.growth, "stretch": grid.stretch, "reach": grid.reach}
    x = nodes(0.0, lay.xmax, lay.xfeat, lay.fine_x, **spread)
    y = nodes(lay.ylo, lay.ymax, lay.yfeat, lay.fine_y, **spread)
    return x, y


def _inside(box: _Box, a: Arr, b: Arr, closed: bool) -> BoolArr:
    """Where the points (a, b) lie in the box: on or inside it (closed), or strictly inside."""
    if closed:
        inside: BoolArr = (a >= box.x0 - _TOL) & (a <= box.x1 + _TOL) & (b >= box.y0 - _TOL) & (b <= box.y1 + _TOL)
    else:
        inside = (a > box.x0) & (a < box.x1) & (b > box.y0) & (b < box.y1)
    return inside


def _receded(box: _Box, lay: _Layout, cut: float, faces: frozenset[str] | None = None) -> _Box:
    """A copper box with its faces inside the domain moved in by `cut` (faces on the domain's edges stay); with
    `faces`, only those of "x0", "x1", "y0", "y1"."""
    if not box.lossy or cut == 0:
        return box
    chosen = frozenset(("x0", "x1", "y0", "y1")) if faces is None else faces

    def move(name: str, inside: bool) -> bool:
        return inside and name in chosen

    return replace(
        box,
        x0=box.x0 + cut if move("x0", box.x0 > _TOL) else box.x0,
        x1=box.x1 - cut if move("x1", box.x1 < lay.xmax - _TOL) else box.x1,
        y0=box.y0 + cut if move("y0", box.y0 > lay.ylo + _TOL) else box.y0,
        y1=box.y1 - cut if move("y1", box.y1 < lay.ymax - _TOL) else box.y1,
    )


def _against(box: _Box, laminate: Sequence[_Box]) -> frozenset[str]:
    """The faces of a copper box that lie on a laminate box's face, over a stretch of some length: the treated side."""
    out = set()
    for lam in laminate:
        across = min(box.x1, lam.x1) - max(box.x0, lam.x0) > _TOL
        along = min(box.y1, lam.y1) - max(box.y0, lam.y0) > _TOL
        if across and abs(box.y0 - lam.y1) <= _TOL:
            out.add("y0")
        if across and abs(box.y1 - lam.y0) <= _TOL:
            out.add("y1")
        if along and abs(box.x0 - lam.x1) <= _TOL:
            out.add("x0")
        if along and abs(box.x1 - lam.x0) <= _TOL:
            out.add("x1")
    return frozenset(out)


def _mirrors(lay: _Layout) -> int:
    """How many copies of the solved part make the whole section."""
    return 2 * (2 if lay.mirror_y != "none" else 1)


def _walls(lay: _Layout, gx: Arr, gy: Arr) -> tuple[BoolArr, BoolArr]:
    """The grounded walls, and every face held at 0 (the walls and the odd mirror planes)."""
    walls: BoolArr = (gx >= lay.xmax - _TOL) | (gy >= lay.ymax - _TOL)
    if lay.mirror_y == "none":
        walls |= gy <= lay.ylo + _TOL
    held = walls.copy()
    if lay.odd_x:
        held |= gx <= _TOL
    if lay.mirror_y == "odd":
        held |= gy <= _TOL
    return walls, held


def _line(lay: _Layout, grid: Grid, laminate: Laminate | None, copper: float, rough: float = 1.0) -> _Line:
    """Solve a layout's line: the air problem as built and with its copper receded (Wheeler), and with the laminate;
    with `rough` other than one, the air problem with only the copper's faces against the laminate receded too."""
    x, y = _grid(lay, grid)
    gx, gy = np.meshgrid(x, y, indexing="ij")
    _, held = _walls(lay, gx, gy)
    cx, cy = np.meshgrid((x[1:] + x[:-1]) / 2, (y[1:] + y[:-1]) / 2, indexing="ij")
    lam = np.zeros(cx.shape, dtype=bool)
    for box in lay.laminate:
        lam |= _inside(box, cx, cy, closed=False)
    ones = np.ones_like(cx)
    whole = _mirrors(lay)

    def solve(cut: float, cells: Arr, treated: bool = False) -> tuple[float, Arr]:
        fixed = held.copy()
        value = np.zeros_like(gx)
        for box in lay.copper:
            faces = _against(box, lay.laminate) if treated else None
            cu = _inside(_receded(box, lay, cut, faces), gx, gy, closed=True)
            fixed |= cu
            value = np.where(cu, box.volts, value)
        phi = laplace2d(x, y, cells, fixed, value)
        return whole * energy2d(x, y, cells, phi), phi

    cut = grid.recede * copper / lay.h
    c0, phi = solve(0.0, ones)
    c1, _ = solve(cut, ones)
    c1_rough = None
    if rough != 1.0 and any(_against(box, lay.laminate) for box in lay.copper):
        c1_rough, _ = solve(cut, ones, treated=True)
    c, share = c0, 0.0
    if laminate is not None and bool(lam.any()):
        eps = np.where(lam, laminate.er, 1.0)
        c, phi = solve(0.0, eps)
        share = whole * energy2d(x, y, np.where(lam, laminate.er, 0.0), phi) / c
    return _Line(x, y, c, c0, c1, share, phi, c1_rough)


def _unfold(a: Arr, axis: int, odd: bool) -> Arr:
    """Values on the half a >= 0 of an axis (index 0 on the mirror) continued to the other half, odd or even."""
    rest = np.flip(np.take(a, np.arange(1, a.shape[axis]), axis=axis), axis=axis)
    out: Arr = np.concatenate([-rest if odd else rest, a], axis=axis)
    return out


def _whole(lay: _Layout, line: _Line) -> tuple[Arr, Arr, Arr]:
    """The grid and the potential over the whole section: the solved part continued by its symmetry."""
    xf = np.concatenate([-line.x[:0:-1], line.x])
    phi = _unfold(line.phi, 0, lay.odd_x)
    if lay.mirror_y == "none":
        return xf, line.y, phi
    yf = np.concatenate([-line.y[:0:-1], line.y])
    return xf, yf, _unfold(phi, 1, lay.mirror_y == "odd")


def _exits(lay: _Layout, x: Arr, y: Arr) -> BoolArr:
    """Where a particle is lost or has escaped: at copper or a wall or beside them, on or in a board, or out of the
    flight box."""
    gx, gy = np.meshgrid(x, y, indexing="ij")
    walls, _ = _walls(lay, gx, gy)
    hit = walls.copy()
    for box in lay.copper:
        hit |= _inside(box, gx, gy, closed=True)
    near = hit.copy()
    near[1:, :] |= hit[:-1, :]
    near[:-1, :] |= hit[1:, :]
    near[:, 1:] |= hit[:, :-1]
    near[:, :-1] |= hit[:, 1:]
    for box in lay.solid:
        near |= _inside(box, gx, gy, closed=True)
    out: BoolArr = near | (gx >= lay.flight[0]) | (gy >= lay.flight[1])
    return out


def _null_odd(x: Arr, y: Arr, phi: Arr, radius: float) -> float:
    """The gradient (per volt, in units of h) of a null at the origin of a potential odd in x and in y: a least-squares
    fit over the nodes within `radius` of the harmonic polynomials of that symmetry, Im((z/radius)^(2k)) for k = 1
    to 3, whose first term c Im(z^2) has the gradient |2c|."""
    gx, gy = np.meshgrid(x, y, indexing="ij")
    near = np.hypot(gx, gy) <= radius
    z = (gx[near] + 1j * gy[near]) / radius
    basis = np.stack([(z ** (2 * k)).imag for k in (1, 2, 3)], axis=1)
    c = np.linalg.lstsq(basis, phi[near], rcond=None)[0]
    return abs(2 * float(c[0])) / radius**2


def _null_axis(x: Arr, y: Arr, phi: Arr, guess: float, radius: float) -> tuple[float, float]:
    """The height and gradient (in units of h) of a null on the axis x = 0 of a potential even in x, near `guess`: a
    least-squares fit over the nodes within `radius` of the harmonic polynomials even in x, Re(((y - guess) + i x)^k)
    for k = 0 to 5; the null is where the fit's slope up the axis vanishes, and the gradient its curvature there."""
    gx, gy = np.meshgrid(x, y, indexing="ij")
    near = np.hypot(gx, gy - guess) <= radius
    w = ((gy[near] - guess) + 1j * gx[near]) / radius
    basis = np.stack([(w**k).real for k in range(6)], axis=1)
    c = np.linalg.lstsq(basis, phi[near], rcond=None)[0]
    slope = np.polynomial.polynomial.polyder(c)
    bend = np.polynomial.polynomial.polyder(slope)
    t = 0.0
    for _ in range(50):
        step = float(np.polynomial.polynomial.polyval(t, slope)) / float(np.polynomial.polynomial.polyval(t, bend))
        t = min(max(t - step, -1.0), 1.0)
    return guess + t * radius, abs(float(np.polynomial.polynomial.polyval(t, bend))) / radius**2


def _saddle(xf: Arr, yf: Arr, emag: Arr, i: int, j: int) -> tuple[float, float]:
    """The saddle of |E| near node (i, j) of the whole grid: the stationary point of a quadratic fitted over the
    node's 3 x 3 neighbourhood, where it is a saddle within that neighbourhood; else the node."""
    if not (0 < i < len(xf) - 1 and 0 < j < len(yf) - 1):
        return float(xf[i]), float(yf[j])
    dx, dy = np.meshgrid(xf[i - 1 : i + 2] - xf[i], yf[j - 1 : j + 2] - yf[j], indexing="ij")
    a, b = dx.ravel(), dy.ravel()
    basis = np.stack([np.ones_like(a), a, b, a * a, a * b, b * b], axis=1)
    c = np.linalg.lstsq(basis, emag[i - 1 : i + 2, j - 1 : j + 2].ravel(), rcond=None)[0]
    hess = np.array([[2 * c[3], c[4]], [c[4], 2 * c[5]]])
    if np.linalg.det(hess) >= 0:
        return float(xf[i]), float(yf[j])
    step = np.linalg.solve(hess, -c[1:3])
    if abs(step[0]) > np.abs(a).max() or abs(step[1]) > np.abs(b).max():
        return float(xf[i]), float(yf[j])
    return float(xf[i] + step[0]), float(yf[j] + step[1])


def solve(section: Section) -> Solved:
    """Solve a cross-section: its line (C', n_eff, Z0, share, Q_c by Wheeler's rule, Q_d, Q), its well (null,
    gradient, saddle, escape) and its potential per volt over the whole section."""
    lay = section.geometry._layout(section, section.grid)
    line = _line(lay, section.grid, section.laminate, section.copper, section.rough)
    x, y, phi = line.x, line.y, line.phi
    xf, yf, whole = _whole(lay, line)
    emag_whole = np.hypot(np.gradient(whole, xf, axis=0), np.gradient(whole, yf, axis=1))
    i0, j0 = len(x) - 1, (len(y) - 1 if lay.mirror_y != "none" else 0)
    emag = emag_whole[i0:, j0:]
    if lay.on_axis:
        axis = np.gradient(phi[0], y)  # -E_y up the axis; E_x vanishes there by symmetry
        up = np.flatnonzero((y[:-1] > lay.above + _TOL) & (np.sign(axis[:-1]) * np.sign(axis[1:]) < 0))
        if not len(up):
            raise ValueError("no null stands above this section's copper: its field does not reverse on the axis")
        k = int(up[0])
        guess = float(y[k] - axis[k] * (y[k + 1] - y[k]) / (axis[k + 1] - axis[k]))
        null_y, gamma = _null_axis(x, y, phi, guess, min(_WINDOW, (guess - lay.above) / 2))
        start = (0, int(np.argmin(np.abs(y - null_y))))
    else:
        null_y = 0.0
        gamma = _null_odd(x, y, phi, _WINDOW)
        start = (0, 0)
    level, (si, sj) = escape_flood(emag, start, _exits(lay, x, y))
    su, sv = _saddle(xf, yf, emag_whole, si + i0, sj + j0)
    h = lay.h
    null_v = (null_y - lay.v0) * h
    distance = null_v if lay.on_axis else section.geometry.distance
    laminate = section.laminate
    qd = 1 / (line.share * laminate.tand) if laminate is not None and line.share * laminate.tand > 0 else math.inf
    u, v = xf * h, (yf - lay.v0) * h
    for a in (u, v, whole):
        a.setflags(write=False)
    return Solved(
        section=section,
        capacitance=line.c * EPS0,
        air_capacitance=line.c0 * EPS0,
        share=line.share,
        qc=wheeler_qc(
            line.c0,
            line.c1,
            section.grid.recede * section.copper,
            section.f,
            section.resistance,
            section.rough,
            line.c1_rough,
        ),
        qd=qd,
        null=(0.0, null_v),
        gradient=gamma / h**2,
        saddle=(abs(su) * h, (sv - lay.v0) * h),
        escape=level / h,
        distance=distance,
        u=u,
        v=v,
        phi=whole,
    )


def fit_well(section: Section, depth_ev: float = 1.0, q: float = 0.5, qm: float = QM_E) -> Solved:
    """The section scaled so that a well depth_ev deep sits at Mathieu's q at its band, and solved: solved at its size,
    its pattern scaled to the distance h the solved well's shares ask (f h from rf.electrode_drive_at_q), and solved
    again, the copper and laminate keeping their thickness (the split as a lever, section 5.4)."""
    first = solve(section)
    fh = electrode_drive_at_q(depth_ev, q, 1.0, 1.0, first.depth_share, first.gradient_share, qm).f
    return solve(section.scaled(fh / section.f / first.distance))


# ---- calibration ---------------------------------------------------------------------------------------------------
def stripline_wheeler(
    width: float, thick: float, spacing: float, f: float = 1e9, grid: Grid | None = None, planes_lose: bool = True
) -> Calibration:
    """The receded solve held to Wheeler's closed form on stripline (rule 5): a strip `width` wide and `thick` thick
    centred between grounded copper planes `spacing` apart (m), solved as the sections are, every copper face receded
    (the strip's and the planes'), against stripline's closed form (Pozar, Microwave Engineering, eq. 3.198, for
    sqrt(eps_r) Z0 < 120 ohm). That form is itself Wheeler's rule applied to a closed-form Z0(t), so agreement shows
    the numerical derivative right, not that the rule gives the loss. With planes_lose False the planes are ideal and
    only the strip recedes: the control, which must miss the closed form by the planes' part of the loss."""
    grid = grid if grid is not None else Grid()
    h = spacing / 2
    tc, half = thick / h, width / (2 * h)
    rc = grid.recede * tc
    top = 1.0 + tc
    lay = _Layout(
        h=h,
        xmax=half + 8.0,
        ylo=0.0,
        ymax=top,
        odd_x=False,
        mirror_y="even",
        copper=(_Box(0.0, half, 0.0, tc / 2, volts=1.0), _Box(0.0, half + 8.0, 1.0, top, lossy=planes_lose)),
        laminate=(),
        solid=(),
        flight=(math.inf, math.inf),
        xfeat=(half, half - rc),
        yfeat=(tc / 2, tc / 2 - rc, 1.0, 1.0 + rc),
        fine_x=grid.fine * tc,
        fine_y=grid.fine * tc,
        v0=0.0,
        on_axis=False,
    )
    line = _line(lay, grid, None, thick)
    qc = wheeler_qc(line.c0, line.c1, grid.recede * thick, f)
    z0 = 1 / (C0 * EPS0 * line.c0)
    alpha = math.pi * f / (C0 * qc)
    rs = math.sqrt(math.pi * f * MU0 / COPPER_SIGMA)
    b, t = spacing, thick
    a = 1 + 2 * width / (b - t) + (b + t) / (b - t) * math.log((2 * b - t) / t) / math.pi
    closed = 2.7e-3 * rs * z0 / (30 * math.pi * (b - t)) * a
    return Calibration(z0, alpha, closed, alpha / closed)


# ---- the lead designs' sections ------------------------------------------------------------------------------------
def section_b(h: float = 0.915e-3, f: float = 915e6, split: float = 0.3, grid: Grid | None = None) -> Section:
    """Section B, T915-S's flex strips (section 5.1): 18 um copper 8h wide split `split` h, on 25 um polyimide 2h wider
    than the copper; h = 0.915 mm at 915 MHz."""
    geometry = FacingStrips(gap=2 * h, width=8 * h, split=split * h, margin=2 * h)
    return Section(geometry, POLYIMIDE, 25e-6, FLEX_COPPER, f, grid=grid if grid is not None else Grid())


def section_a(
    laminate: Laminate = RO4003C,
    thickness: float = 0.51e-3,
    h: float = 0.92e-3,
    f: float = 915e6,
    grid: Grid | None = None,
) -> Section:
    """Section A (section 5.1): 35 um copper strips 1.5h wide either side of a split h on rigid boards routed as a slot
    under the split, the board out to the walls; RO4003C 0.51 mm at h = 0.92 mm and 915 MHz."""
    geometry = FacingStrips(gap=2 * h, width=4 * h, split=h, slotted=True)
    return Section(geometry, laminate, thickness, RIGID_COPPER, f, grid=grid if grid is not None else Grid())


def slot_guide(
    laminate: Laminate = RO4003C,
    thickness: float = 0.813e-3,
    strip: float = 0.5e-3,
    planes: float = 1.0e-3,
    f: float = 915e6,
    grid: Grid | None = None,
) -> Section:
    """The slot guide of section 3.3: a 1.0 mm slot routed in a double-sided board, 35 um copper strips `strip` wide
    set back 0.2 mm from its walls, the other board's copper `planes` off; RO4003C 0.813 mm (S2-C3) or FR-4 0.8 mm."""
    geometry = SlotGuide(slot=1.0e-3, strip=strip, pull=0.2e-3, planes=planes)
    return Section(geometry, laminate, thickness, RIGID_COPPER, f, grid=grid if grid is not None else Grid())
