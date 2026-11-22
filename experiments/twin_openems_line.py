#!/usr/bin/env python
"""R1's rail line in openEMS against its 2-D section: the twin bench's first full-wave check (stage 0, DR-0).

docs/KNOTTED_GUIDE.md sections 3.6 and 9 item 5 put R1's rail line on one board: a grounded centre strip between two
rails driven in common at 915 MHz, ELE-1's a = 0.40 mm and b = 1.46 mm with the null 0.575 mm up, on 1.6 mm FR-4 with
35 um copper. Section 8.1's Hypothesis and DR-0 ask the decomposed model, a 2-D section of the line closed into a ring,
to give the resonance within 0.5 % and Q within 15 % of full wave. This script computes both for one geometry:

  section    the decomposed model: the 2-D section by finite volumes with the knotted study's solver (nodes,
             laplace2d, energy2d, wheeler_qc): C', n_eff, Z0, the laminate's share of the field energy and
             Q_d = 1 / (share tan d), Q_c by Wheeler's incremental-inductance rule on the real 35 um copper (every
             face receded by t/8), on the study's converged grading (t/16 at the copper); the same on the FDTD's own
             grid, which parts the mesh's error from the full wave's; the section as written, without the plane
  straight   openEMS on half the section (a magnetic wall on the centre plane: the rails in phase) between two
             magnetic walls normal to the line. By images that is a ring of twice its length without curvature: its
             first resonance gives n_eff with no end correction, its ring-down the dielectric Q, probes across a gap
             and round a rail Z0, and the surface integral of R_s |H_t|^2 / 2 over the copper at the current
             antinode the copper's Q; on two meshes
  ring       openEMS on R1's ring of the same section at m = 1, in a cylindrical mesh so that every copper edge lies
             on a mesh line (no staircase): a quarter of it between a magnetic wall at alpha = 0 and an electric wall
             at alpha = pi/2, which hold m = 1's cos(alpha) standing waves; each mode's frequency, Q and rail volts

The board as modelled. The etched gaps, 0.15 mm, are centred on ELE-1's edges c1 = a/2 and c2 = a/2 + b, so the null
stays where the gapless model puts it and the copper is 0.25 mm (strip) and 1.31 mm (rails) wide. Left floating, the
centre strip takes 0.85 to 0.92 of the rails' volts through the gaps (with and without a plane under the board), which
removes the null; it is held at ground along the line by a via wall to a ground plane under the board, and the top
ground by a second one 1 mm beyond the outer gap, so all ground copper is one conductor. Via walls are continuous (a
dense via row) and carry no loss in either model; the line sits in a grounded shield 25 mm from it, the same in both
models. FR-4's loss in openEMS is a conductivity giving tan d = 0.02 at 915 MHz, so its Q is compared at that tan d.

The two rails are two conductors: straight, the even mode (rails in phase) and the odd one part by symmetry; closed
into the ring, the rails' lengths differ by -+(c1 + c2)/2R, and the ring resonates on each rail apart. The ring stage
reports both modes against a two-rail decomposed model (the rails' 2 x 2 line matrices, each rail on its own radius)
beside the single-mode one.

    .env/bin/python experiments/twin_openems_line.py            # every stage (resumes), about half an hour
    .env/bin/python experiments/twin_openems_line.py section    # the 2-D sections alone (also under uv)
    .env/bin/python experiments/twin_openems_line.py report     # run.log and the figure from results.json
    .env/bin/python experiments/twin_openems_line.py all --quick   # coarse meshes, a minute, into the scratch

Outputs: twin_openems_line/{run.log, results.json, twin_openems_line.png}. openEMS works in TWIN_OPENEMS_SCRATCH (the
system's temporary directory by default), each folder deleted once its numbers are stored; TWIN_OPENEMS_ENGINE picks
the engine (gpu by default, multithreaded for the CPU).
"""

from __future__ import annotations

import argparse
import importlib
import json
import math
import os
import shutil
import sys
import tempfile
import time
from collections.abc import Callable, Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

import knotted_guide_study as kg
import numpy as np

from array16 import lab
from array16 import particles as pt
from array16.particles import Arr, BoolArr, CArr

OUT = Path(__file__).resolve().parent / "twin_openems_line"
SCRATCH = Path(os.environ.get("TWIN_OPENEMS_SCRATCH", tempfile.gettempdir())) / "array16_twin_openems"
ENGINE = os.environ.get("TWIN_OPENEMS_ENGINE", "gpu")
LOG: list[str] = []

F915 = kg.F915
FC = 0.5e9  # Hz, the Gaussian excitation's half band: 0.415 to 1.415 GHz, below the shield's first mode near 3 GHz
SIM_NS = 30.0  # ns simulated: the 5.7 ns pulse, then about 18 periods of ring-down
SIGMA = kg.COPPER_SIGMA
R_S = math.sqrt(math.pi * F915 * pt.MU0 / SIGMA)  # ohm, bare copper's surface resistance at 915 MHz
DELTA = math.sqrt(2 / (2 * math.pi * F915 * pt.MU0 * SIGMA))  # m, its skin depth
# DR-0 (docs/KNOTTED_GUIDE.md sections 8.1 and 8.2), registered here before any run: the decomposed model's resonance
# within 0.5 % of full wave and its Q within 15 %
DR0_F, DR0_Q = 0.005, 0.15


def say(line: str = "") -> None:
    print(line, flush=True)
    LOG.append(line)


# ---- the line ------------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class RailLine:
    """R1's rail line on one board, lengths in mm, x across the line from its centre plane and y up from the board's
    top face. ELE-1's grounded strip `strip` (a) and rails `rail` (b) wide, etched apart by gaps `gap` wide centred on
    its edges c1 = a/2 and c2 = a/2 + b, copper `copper` thick on a board `board` thick (relative permittivity `er`,
    loss tangent `tand`). Backed: a ground plane under the board, a via wall on the centre plane tying the strip to it
    and one `fence` beyond the outer gap tying the top ground; else the board in air and the strip held at ground. A
    grounded shield `shield` from the centre plane, from the board's top face and (unbacked) below it."""

    strip: float = 0.40
    rail: float = 1.46
    gap: float = 0.15
    copper: float = kg.RIGID_COPPER_T / kg.MM
    board: float = 1.6
    er: float = kg.FR4[0]
    tand: float = kg.FR4[1]
    backed: bool = True
    fence: float = 1.0
    shield: float = 25.0

    def __post_init__(self) -> None:
        if not (0 < self.gap < min(self.strip, self.rail) and self.copper > 0 and self.board > 0 and self.fence > 0):
            raise ValueError(f"a rail line needs 0 < gap < strip, rail and positive copper, board, fence: {self}")
        if self.shield < self.fence_x + 2.0:
            raise ValueError(f"the shield at {self.shield} mm must clear the fence at {self.fence_x} mm by 2 mm")

    @property
    def edges(self) -> tuple[float, float]:
        """ELE-1's edges c1 = a/2 and c2 = a/2 + b (mm), the gaps' centres."""
        return self.strip / 2, self.strip / 2 + self.rail

    @property
    def spans(self) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
        """The copper across x >= 0 (mm): the centre strip from the centre plane, the rail, the top ground to the
        shield."""
        c1, c2 = self.edges
        g = self.gap / 2
        return (0.0, c1 - g), (c1 + g, c2 - g), (c2 + g, self.shield)

    @property
    def fence_x(self) -> float:
        """The top ground's via wall (mm from the centre plane)."""
        return self.edges[1] + self.gap / 2 + self.fence

    @property
    def floor(self) -> float:
        """The section's floor (mm): the plane's top face, or the shield below the board."""
        return -self.board if self.backed else -self.shield

    @property
    def null_gapless(self) -> float:
        """ELE-1's null height sqrt(c1 c2) (mm) over the gapless pattern."""
        return pt.RailPair(self.strip * kg.MM, self.rail * kg.MM).null_height / kg.MM


@dataclass(frozen=True)
class Grading:
    """A graded grid: `fine` (mm) at the copper's edges and faces, the spacing growing by `growth` per mm away from
    them up to `coarse` (the knotted study's nodes())."""

    fine: float
    coarse: float = 1.0
    growth: float = 0.1


R1 = RailLine()
AS_WRITTEN = replace(R1, backed=False)
CONVERGED = Grading(fine=R1.copper / 16)  # the study's section grid: t/16 at the copper, growth 0.1 per mm
FDTD_FINE = Grading(fine=0.005, growth=0.3)
FDTD_COARSE = Grading(fine=0.010, growth=0.3)
QUICK = Grading(fine=0.020, growth=0.4)


def section_lines(
    line: RailLine, grading: Grading, recede: float | None = None, extra_y: Sequence[float] = ()
) -> tuple[Arr, Arr]:
    """The section's grid nodes (mm): x from the centre plane to the shield, y from the floor (or, receded, the
    plane's copper below it) to the shield above, fine at the copper and at `extra_y`; with `recede`, lines that far
    inside every copper face too."""
    (_, s1), (r0, r1), (g0, _) = line.spans
    xs = [0.0, s1, r0, r1, g0, line.shield] + ([line.fence_x] if line.backed else [])
    ys = [line.floor, -line.board, 0.0, line.copper, line.shield, *extra_y]
    floor = line.floor
    if recede is not None:
        xs += [s1 - recede, r0 + recede, r1 - recede, g0 + recede]
        ys += [recede, line.copper - recede]
        if line.backed:
            ys.append(-line.board - recede)
            floor = -line.board - line.copper
    x = kg.nodes(0.0, line.shield, xs, grading.fine, grading.coarse, grading.growth)
    y = kg.nodes(floor, line.shield, ys, grading.fine, grading.coarse, grading.growth)
    return x, y


def full_lines(line: RailLine, grading: Grading) -> tuple[Arr, Arr]:
    """The FDTD's cross-section (mm): the half grid mirrored across the centre plane without a node on it, so that a
    magnetic wall on the first dual line (where openEMS puts it, half a cell inside the boundary) is the centre plane,
    and y as section_lines."""
    x, y = section_lines(line, grading)
    out: Arr = np.concatenate([-x[:0:-1], x[1:]])
    return out, y


def node_widths(v: Arr) -> Arr:
    """Each node's share of the grid along one axis (half a cell at the ends)."""
    h = np.diff(v)
    w = np.zeros(len(v))
    w[:-1] += h / 2
    w[1:] += h / 2
    return w


@dataclass(frozen=True)
class Masks:
    """Node masks of a section (any shape of nodes): the centre strip, the rail(s), the ground held at 0 V (top
    ground, plane, via walls, shield), the copper whose faces carry its loss (strip, rails, top ground, plane), and
    the free nodes."""

    centre: BoolArr
    rail: BoolArr
    ground: BoolArr
    lossy: BoolArr
    free: BoolArr


def section_masks(
    line: RailLine, ax: Arr, gy: Arr, cut: float = 0.0, tied: bool = True, half_wall: float = 0.0
) -> Masks:
    """The masks at nodes |x| = ax, y = gy (mm), every copper face receded by `cut`, the centre via wall `half_wall`
    either side of the centre plane (zero: a sheet on a node there); not tied: the centre strip without its via wall
    (the modal and floating analyses)."""
    e = 1e-9
    t = line.copper
    (_, s1), (r0, r1), (g0, _) = line.spans
    layer = (gy >= cut - e) & (gy <= t - cut + e)
    centre = layer & (ax <= s1 - cut + e)
    rail = layer & (ax >= r0 + cut - e) & (ax <= r1 - cut + e)
    top = layer & (ax >= g0 + cut - e)
    shield = (ax >= line.shield - e) | (gy >= line.shield - e)
    if not line.backed:  # the floor is the shield below the board; backed, it is the plane, whose face is lossy
        shield |= gy <= float(gy.min()) + e
    plane = (gy <= -line.board - cut + e) if line.backed else np.zeros_like(layer)
    walls = (gy >= -line.board - cut - e) & (gy <= cut + e)
    fence = walls & (np.abs(ax - line.fence_x) <= e) if line.backed else np.zeros_like(layer)
    middle = walls & (ax <= half_wall + e) if (line.backed and tied) else np.zeros_like(layer)
    ground = top | plane | fence | middle | shield
    lossy = (centre | rail | top | plane) & ~shield
    free = ~(centre | rail | ground)
    return Masks(centre, rail, ground, lossy, free)


def face_sum(
    x: Arr, y: Arr, m: Masks, fy: Arr | CArr, fx: Arr | CArr, wr: Arr | None = None, wrd: Arr | None = None
) -> float:
    """The sum over the lossy copper's faces of |F_t|^2 on the edge beside each face node times the node's width
    along the face (mm times F's unit squared): fy (nx, ny - 1) on the edges along y, fx (nx - 1, ny) on the edges
    along x, the Yee places of H_x and H_y (or E_y and E_x of the static field). For a ring, wr and wrd weight each
    node and each dual line by its radius."""
    wx, wy = node_widths(x), node_widths(y)
    if wr is not None and wrd is not None:
        wx = wx * wr
        wy_x = wrd[:, None] * wy[None, :]
    else:
        wy_x = np.broadcast_to(wy[None, :], fx.shape)
    across_y = (m.lossy[:, :-1] & m.free[:, 1:]) | (m.free[:, :-1] & m.lossy[:, 1:])
    across_x = (m.lossy[:-1, :] & m.free[1:, :]) | (m.free[:-1, :] & m.lossy[1:, :])
    sy = (np.abs(fy) ** 2 * wx[:, None])[across_y].sum()
    sx = (np.abs(fx) ** 2 * wy_x)[across_x].sum()
    return float(sy + sx)


def copper_q(energy: float, faces: float, f: float = F915) -> float:
    """Q_c = omega mu0 int|H|^2 dA / (R_s sum|H_t|^2 dl): the magnetic energy a length of line holds over the power
    its copper loses, from one cross-section of a TEM-like standing wave (energy in mm^2, faces in mm, the field's
    unit squared alike); R_s scales as sqrt(f)."""
    rs = R_S * math.sqrt(f / F915)
    return 2 * math.pi * f * pt.MU0 * energy * 1e-6 / (rs * faces * 1e-3)


# ---- the decomposed model: the 2-D section -------------------------------------------------------------------------
def solve_section(
    line: RailLine,
    grading: Grading,
    wheeler: bool = True,
    modal: bool = False,
    full: bool = False,
    at_null: bool = False,
) -> dict[str, Any]:
    """The 2-D section in its even mode (the rails at 1 V, as the guide drives it) and its odd mode (the rails at
    -1 and +1 V): on half the line (x >= 0, a node on the centre plane, which is free for the even mode and held at
    0 V for the odd), or `full`, on the FDTD's own cross-section (full_lines, the centre via wall between the two
    nodes beside the plane). Each mode: C' and the air problem's C' (pF/m, the whole line), n_eff, Z0 (ohm), the
    laminate's share, Q_d, Q_c by Wheeler's rule (receded t/8) and by the face integral of the static field on this
    grid (the 3-D integral's quasi-static twin), Q; the even mode's keys at the top, the odd's under "odd", and the
    rails' 2 x 2 matrices (inner, outer; eps0 units) under "rails". On the half line also the null's height above
    the copper's top face (mm) and its gradient over ELE-1's, read well only `at_null` (the grid fine there too);
    modal: the strip without its via wall, its floating volts and the two quasi-TEM modes' indices."""
    began = time.time()
    rc = line.copper / 8
    if full:
        x, y = full_lines(line, grading)
        wheeler = modal = False
    else:
        extra = [line.copper + line.null_gapless] if at_null else []
        x, y = section_lines(line, grading, rc if wheeler else None, extra)
    halves = 1.0 if full else 2.0
    half_wall = float(x[x > 0][0]) if full else 0.0
    gx, gy = np.meshgrid(x, y, indexing="ij")
    ax = np.abs(gx)
    _, cy = np.meshgrid((x[1:] + x[:-1]) / 2, (y[1:] + y[:-1]) / 2, indexing="ij")
    lam = (cy < 0) & (cy > -line.board)
    eps = np.where(lam, line.er, 1.0)
    ones = np.ones_like(eps)
    layer_only = np.where(lam, line.er, 0.0)
    plane = np.zeros_like(gx, dtype=bool) if full else gx <= 1e-9  # the centre plane's nodes, half line

    def solve(cells: Arr, m: Masks, odd: bool) -> Arr:
        rails = np.where(gx < 0, -1.0, 1.0) if (odd and full) else np.ones_like(gx)
        fixed = m.rail | m.ground | m.centre | (plane if odd else np.zeros_like(plane))
        phi: Arr = kg.laplace2d(x, y, cells, fixed, np.where(m.rail, rails, 0.0), tol=1e-12)
        return phi

    m = section_masks(line, ax, gy, half_wall=half_wall)
    out: dict[str, Any] = {"grading": asdict(grading), "full": full, "nodes": len(x) * len(y)}
    sums: dict[str, dict[str, float]] = {}
    phi_even = np.zeros_like(gx)
    for odd in (False, True):
        phi_air = solve(ones, m, odd)
        phi = solve(eps, m, odd)
        ca = halves * kg.energy2d(x, y, ones, phi_air)
        cl = halves * kg.energy2d(x, y, eps, phi)
        held = halves * kg.energy2d(x, y, layer_only, phi)
        fy = np.diff(phi_air, axis=1) / np.diff(y)[None, :]
        fx = np.diff(phi_air, axis=0) / np.diff(x)[:, None]
        faces = halves * face_sum(x, y, m, fy, fx)
        row: dict[str, Any] = {
            "C_pF": cl * pt.EPS0 * 1e12,
            "Ca_pF": ca * pt.EPS0 * 1e12,
            "index": math.sqrt(cl / ca),
            "z0": 1 / (pt.C0 * pt.EPS0 * math.sqrt(cl * ca)),
            "share": held / cl,
            "Qd": cl / (held * line.tand),
            "Qc_faces": copper_q(ca, faces),
        }
        sums[str(odd)] = {"C": cl, "Ca": ca, "S": held, "F": faces}
        if wheeler:
            ca1 = halves * kg.energy2d(x, y, ones, solve(ones, section_masks(line, ax, gy, rc), odd))
            row["Qc"] = kg.wheeler_qc(ca, ca1, rc * kg.MM, F915)
            sums[str(odd)]["Ca1"] = ca1
        row["Q"] = 1 / (1 / row.get("Qc", row["Qc_faces"]) + 1 / row["Qd"])
        if odd:
            out["odd"] = row
        else:
            out.update(row)
            phi_even = phi
    # the rails' matrices from the two modes: v M v = 2 (M11 + M12) even, 2 (M11 - M12) odd, M11 = M22 by symmetry
    out["rails"] = {
        key: [
            [(sums["False"][key] + sums["True"][key]) / 4, (sums["False"][key] - sums["True"][key]) / 4],
            [(sums["False"][key] - sums["True"][key]) / 4, (sums["False"][key] + sums["True"][key]) / 4],
        ]
        for key in sums["False"]
    }
    if not full:
        out.update(null_on_axis(line, x, y, phi_even))
    if modal:
        out.update(modes(line, x, y, gx, gy, eps, ones))
    out["seconds"] = time.time() - began
    return out


def two_rail_ring(rails: dict[str, Any], line: RailLine, radius: float, tand: float) -> list[dict[str, float]]:
    """The decomposed model of the ring with each rail its own conductor. The rails' 2 x 2 matrices per unit length
    (eps0 units: C with the laminate, the air problem's Ca, Wheeler's receded Ca1, the laminate's energy S, the face
    sum F) are taken per radian on each rail's centre radius R -+ (c1 + c2)/2, D M D with D = diag(sqrt R); m = 1's
    standing waves cos(alpha) have omega^2 the eigenvalues of (L C)^-1 per radian, L = mu0 Ca^-1. Each mode: its
    frequency, the rails' volts (inner over outer), Q_d = v C v / (tan d v S v), and Q_c by Wheeler's rule and by the
    face sum on the air problem that carries the mode's currents, u = Ca^-1 D C v; the lowest frequency first."""
    cbar = sum(line.edges) / 2
    root = np.diag(np.sqrt(np.array([radius - cbar, radius + cbar]) * kg.MM))
    m = {key: np.array(value, dtype=np.float64) for key, value in rails.items()}
    c_alpha = root @ m["C"] @ root * pt.EPS0
    l_alpha = root @ (pt.MU0 * np.linalg.inv(m["Ca"])) @ root
    s_alpha = root @ m["S"] @ root * pt.EPS0
    w2, vecs = np.linalg.eig(np.linalg.inv(l_alpha @ c_alpha))
    modes_out = []
    for k in np.argsort(w2.real):
        v = vecs[:, k].real
        v = v / v[np.argmax(np.abs(v))]
        u = np.linalg.solve(m["Ca"], root @ (c_alpha @ v))
        row = {
            "f": math.sqrt(w2[k].real) / (2 * math.pi),
            "inner_over_outer": float(v[0] / v[1]),
            "Qd": float(v @ c_alpha @ v / (tand * (v @ s_alpha @ v))),
            "Qc_faces": copper_q(float(u @ m["Ca"] @ u), float(u @ m["F"] @ u)),
        }
        if "Ca1" in m:
            ratio = float(u @ m["Ca"] @ u) / float(u @ m["Ca1"] @ u)
            row["Qc"] = 2 * line.copper / 8 * kg.MM / (DELTA * (ratio - 1))
        row["Q"] = 1 / (1 / row.get("Qc", row["Qc_faces"]) + 1 / row["Qd"])
        modes_out.append(row)
    return modes_out


def null_on_axis(line: RailLine, x: Arr, y: Arr, phi: Arr) -> dict[str, float]:
    """The null of E_y on the centre plane above the copper (mm over the top face) and the field's gradient there over
    ELE-1's 4 V b / (pi h (a + b)^2) at V = 1; NaN where the field has no null."""
    above = np.flatnonzero(y > line.copper)
    ys, column = y[above], phi[0, above]
    ey = -np.gradient(column, ys)
    zeros = kg.zero_crossings(ys, ey)
    if not zeros:
        return {"null_mm": math.nan, "gradient_ele1": math.nan}
    height = zeros[0][0]
    near = np.abs(ys - height) <= 0.15  # the potential near the null is a parabola, its curvature the gradient
    curvature = 2 * np.polyfit(ys[near] - height, column[near], 3)[-3]
    pair = pt.RailPair(line.strip * kg.MM, line.rail * kg.MM)
    return {"null_mm": height - line.copper, "gradient_ele1": abs(curvature) * 1e6 / pair.gradient(1.0)}


def modes(line: RailLine, x: Arr, y: Arr, gx: Arr, gy: Arr, eps: Arr, ones: Arr) -> dict[str, Any]:
    """The strip as its own conductor (no via wall): the capacitance matrices of rail and strip against ground with
    and without the laminate, the strip's floating volts (charge zero) over the rails', the two quasi-TEM modes'
    indices (eigenvalues of C_air^-1 C), and the null with the strip floating."""
    m = section_masks(line, gx, gy, tied=False)
    fixed = m.rail | m.centre | m.ground
    mats: dict[str, Arr] = {}
    fields: dict[str, list[Arr]] = {}
    for name, cells in (("air", ones), ("laminate", eps)):
        solved = []
        energies = []
        for vr, vc in ((1.0, 0.0), (0.0, 1.0), (1.0, 1.0)):
            value = np.where(m.rail, vr, 0.0) + np.where(m.centre, vc, 0.0)
            phi = kg.laplace2d(x, y, cells, fixed, value, tol=1e-12)
            solved.append(phi)
            energies.append(2 * kg.energy2d(x, y, cells, phi))
        c12 = (energies[2] - energies[0] - energies[1]) / 2
        mats[name] = np.array([[energies[0], c12], [c12, energies[1]]])
        fields[name] = solved
    lam = mats["laminate"]
    floating = -lam[0, 1] / lam[1, 1]
    phi_float = fields["laminate"][0] + floating * fields["laminate"][1]
    eig = np.linalg.eigvals(np.linalg.solve(mats["air"], lam))
    return {
        "floating": float(floating),
        "floating_air": float(-mats["air"][0, 1] / mats["air"][1, 1]),
        "mode_index": sorted(float(v) for v in np.sqrt(eig.real)),
        "null_floating_mm": null_on_axis(line, x, y, phi_float)["null_mm"],
        "C_matrix_eps0": lam.tolist(),
    }


def section_job(job: tuple[str, RailLine, Grading, bool, bool, bool]) -> dict[str, Any]:
    """One section solve, for a process pool; a name starting with "null" grades the grid fine at the null too."""
    name, line, grading, wheeler, modal, full = job
    return {"name": name, **solve_section(line, grading, wheeler, modal, full, name.startswith("null"))}


def stripline_job(job: tuple[float, float, float]) -> dict[str, Any]:
    """The study's calibration of Wheeler's rule on stripline (rule 5), for a process pool."""
    return {"name": "stripline", "case": list(job), **kg.stripline_wheeler(*job)}


def run_pool(func: Callable[[Any], dict[str, Any]], jobs: Sequence[Any]) -> list[dict[str, Any]]:
    """Each job through func in a process pool, in order."""
    workers = max(1, min(len(jobs), (os.cpu_count() or 2) - 1))
    if workers == 1:
        return [func(job) for job in jobs]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(func, jobs))


def stage_section(quick: bool) -> dict[str, Any]:
    """Every 2-D solve: the decomposed model of R1 as built and as written, on the converged grid and the FDTD's,
    the face integral's convergence, the gap and the shield as levers, and the stripline calibration."""
    converged = QUICK if quick else CONVERGED
    jobs: list[tuple[str, RailLine, Grading, bool, bool, bool]] = [
        ("built", R1, converged, True, True, False),
        ("written", AS_WRITTEN, converged, True, True, False),
        ("fdtd quick", R1, QUICK, False, False, True),  # the FDTD's own cross-sections: its static limit, cell for cell
        ("null built", R1, converged, False, False, False),
        ("null written", AS_WRITTEN, converged, False, False, False),
    ]
    if not quick:
        jobs += [
            ("fdtd fine", R1, FDTD_FINE, False, False, True),
            ("fdtd coarse", R1, FDTD_COARSE, False, False, True),
            ("gap 0.10", replace(R1, gap=0.10), converged, True, False, False),
            ("gap 0.20", replace(R1, gap=0.20), converged, True, False, False),
            ("shield 50", replace(R1, shield=50.0), converged, True, False, False),
            ("coarse 0.5", R1, replace(converged, coarse=0.5), True, False, False),
            (
                "RO4003C 1.524",
                replace(R1, board=1.524, er=kg.RO4003C[0], tand=kg.RO4003C[1]),
                converged,
                True,
                False,
                False,
            ),
        ]
        for div in (8, 16, 32, 64):  # the face integral's convergence on the half line, no receded lines
            jobs.append((f"faces t/{div}", R1, Grading(fine=R1.copper / div), False, False, False))
    began = time.time()
    rows = run_pool(section_job, jobs)
    cal = run_pool(stripline_job, kg.STRIPLINES[:1] if quick else list(kg.STRIPLINES))
    out: dict[str, Any] = {row.pop("name"): row for row in rows}
    out["stripline"] = cal
    out["line"] = asdict(R1)
    out["seconds"] = time.time() - began
    return out


# ---- the full wave: openEMS ---------------------------------------------------------------------------------------
def openems() -> tuple[Any, Any]:
    """openEMS and CSXCAD, imported when a full-wave stage runs (they live in .env, not in uv's environment)."""
    return importlib.import_module("openEMS"), importlib.import_module("CSXCAD")


def steps_for(grading: Grading, along: float, sim_ns: float) -> int:
    """A cap on the timesteps for `sim_ns` of simulated time (openEMS stops at the time itself): the Courant step
    of cells a quarter of `fine`, below any the grading makes."""
    dt = 1 / (pt.C0 * math.sqrt(2 / (grading.fine / 4 * kg.MM) ** 2 + 1 / (along * kg.MM) ** 2))
    return math.ceil(sim_ns * 1e-9 / dt)


def quantization(line: RailLine, dt: float) -> dict[str, float]:
    """openEMS holds each edge's update coefficient VV = (1 - x)/(1 + x), x = dt sigma / 2 eps, in single precision
    (FDTD_FLOAT is float), so a lossy cell's damping per step 1 - VV, about 1e-6 here, is rounded to a whole number of
    2^-24: the laminate's loss in the run is `factor` times the loss it was given (the bulk's edges; the few at its
    faces round apart), and a ring-down's Q_d times `factor` is the given tan d's. With the step in 2^-24 units."""
    x = math.pi * F915 * line.tand * dt  # dt kappa / (2 eps), kappa = 2 pi 915 MHz eps tan d
    exact = (1 - x) / (1 + x)
    held = float(np.float32(exact))
    return {"dt": dt, "damping_steps": (1 - exact) * 2**24, "factor": (1 - held) / (1 - exact)}


def nearest(lines: Arr, value: float) -> float:
    return float(lines[int(np.argmin(np.abs(lines - value)))])


def dual(lines: Arr) -> Arr:
    out: Arr = (lines[1:] + lines[:-1]) / 2
    return out


def straight_model(
    line: RailLine,
    grading: Grading,
    length: float,
    cells: int,
    steps: int,
    sim_ns: float,
    path: Path,
    f_dump: Sequence[float],
) -> dict[str, Any]:
    """Half the section along z between magnetic walls at z = 0 and z = `length` (mm), `cells` (even) cells between
    them, built and run in `path`. openEMS puts a magnetic wall on the first and last dual lines, so the x lines start
    one ghost line behind the centre plane (the half of full_lines) and the z lines half a cell beyond each end; the
    centre via wall fills the two lines beside the plane. A soft source across the outer gap near z = 0; voltage
    probes across both gaps on the first and last lines inside (half a cell from the walls) and the line nearest the
    middle; a current probe round the rail and an H dump (frequency domain, no interpolation) at the middle, a dual
    line. Returns the lines, the probes' places and the run's wall time."""
    oe, cs = openems()
    full, y = full_lines(line, grading)
    x = full[len(full) // 2 - 1 :]  # from the ghost line -x1 on
    x1 = float(x[1])
    dz = length / cells
    z: Arr = (np.arange(cells + 2, dtype=np.float64) - 0.5) * dz
    fdtd = oe.openEMS(NrTS=steps, EndCriteria=0, MaxTime=sim_ns * 1e-9)
    fdtd.SetGaussExcite(F915, FC)
    fdtd.SetBoundaryCond(["PMC", "PEC", "PEC", "PEC", "PMC", "PMC"])
    csx = cs.ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(kg.MM)
    grid.AddLine("x", x)
    grid.AddLine("y", y)
    grid.AddLine("z", z)
    z0, z1 = float(z[0]), float(z[-1])
    kappa = 2 * math.pi * F915 * pt.EPS0 * line.er * line.tand
    board = csx.AddMaterial("fr4", epsilon=line.er, kappa=kappa)
    board.AddBox([-x1, -line.board, z0], [line.shield, 0.0, z1], priority=1)
    copper = csx.AddMetal("copper")
    (_, s1), (r0, r1), (g0, _) = line.spans
    copper.AddBox([-x1, 0.0, z0], [s1, line.copper, z1], priority=10)
    copper.AddBox([r0, 0.0, z0], [r1, line.copper, z1], priority=10)
    copper.AddBox([g0, 0.0, z0], [line.shield, line.copper, z1], priority=10)
    if line.backed:
        copper.AddBox([-x1, -line.board, z0], [x1, 0.0, z1], priority=10)
        copper.AddBox([line.fence_x, -line.board, z0], [line.fence_x, 0.0, z1], priority=10)
    ym = nearest(y, line.copper / 2)
    source = csx.AddExcitation("source", exc_type=0, exc_val=[1.0, 0.0, 0.0])
    source.AddBox([r1, ym, float(z[2])], [g0, ym, float(z[2])])
    probes = {"v0": float(z[1]), "vm": nearest(z, length / 2 - dz / 2), "vl": float(z[-2])}
    for name, zp in probes.items():
        csx.AddProbe(name, p_type=0).AddBox([r1, ym, zp], [g0, ym, zp])  # rail to top ground
        csx.AddProbe(name + "_in", p_type=0).AddBox([r0, ym, zp], [s1, ym, zp])  # rail to strip
    xd, yd = dual(x), dual(y)
    zc = length / 2  # a dual line: z = (k - 1/2) dz, cells even
    loop = (
        [nearest(xd, line.edges[0]), nearest(yd, -line.board / 2), zc],
        [nearest(xd, line.edges[1]), nearest(yd, 3 * line.copper), zc],
    )
    csx.AddProbe("ic", p_type=1, norm_dir=2).AddBox(*loop)
    dump = csx.AddDump("hcut", dump_type=11, dump_mode=0, file_type=1, frequency=list(f_dump))
    dump.AddBox([float(x[0]), float(y[0]), zc], [float(x[-1]), float(y[-1]), zc])
    began = time.time()
    fdtd.Run(str(path), cleanup=True, verbose=0, engine=ENGINE)
    return {
        "x": x,
        "y": y,
        "z": z,
        "x1": x1,
        "z_current": zc,
        "probes": probes,
        "loop": loop,
        "wall_s": time.time() - began,
        "cells": (len(x) - 1) * (len(y) - 1) * (len(z) - 1),
        "steps": steps,
    }


def ring_model(
    line: RailLine,
    grading: Grading,
    radius: float,
    lines: int,
    steps: int,
    sim_ns: float,
    path: Path,
    f_dump: Sequence[float],
) -> dict[str, Any]:
    """R1's ring: the whole section (full_lines) round the centre-line radius `radius` (mm) in a cylindrical mesh,
    from a magnetic wall at alpha = 0 (openEMS's first dual line, so the alpha lines start half a cell behind it) to
    an electric wall at alpha = pi/2 (its last line), `lines` lines in alpha: soft sources across both outer gaps near
    alpha = 0 (the rails in phase), voltage probes across the four gaps on the first line inside, current probes
    round both rails and an H dump on the last dual line, half a cell from pi/2, built and run in `path`."""
    oe, cs = openems()
    x, y = full_lines(line, grading)
    x1 = float(x[x > 0][0])
    r = radius + x
    da = math.pi / (2 * lines - 3)
    a = (np.arange(lines) - 0.5) * da
    fdtd = oe.openEMS(NrTS=steps, EndCriteria=0, CoordSystem=1, MaxTime=sim_ns * 1e-9)
    fdtd.SetGaussExcite(F915, FC)
    fdtd.SetBoundaryCond(["PEC", "PEC", "PMC", "PEC", "PEC", "PEC"])
    csx = cs.ContinuousStructure(CoordSystem=1)
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(kg.MM)
    grid.SetLines("r", r)
    grid.SetLines("a", a)
    grid.SetLines("z", y)
    a0, a1 = float(a[0]), float(a[-1])
    kappa = 2 * math.pi * F915 * pt.EPS0 * line.er * line.tand
    board = csx.AddMaterial("fr4", epsilon=line.er, kappa=kappa)
    board.AddBox([float(r[0]), a0, -line.board], [float(r[-1]), a1, 0.0], priority=1)
    copper = csx.AddMetal("copper")
    (_, s1), (r0, r1), (g0, _) = line.spans
    copper.AddBox([radius - s1, a0, 0.0], [radius + s1, a1, line.copper], priority=10)
    for lo, hi in ((r0, r1), (g0, line.shield)):
        copper.AddBox([radius + lo, a0, 0.0], [radius + hi, a1, line.copper], priority=10)
        copper.AddBox([radius - hi, a0, 0.0], [radius - lo, a1, line.copper], priority=10)
    if line.backed:
        copper.AddBox([radius - x1, a0, -line.board], [radius + x1, a1, 0.0], priority=10)
        for rv in (radius - line.fence_x, radius + line.fence_x):
            copper.AddBox([rv, a0, -line.board], [rv, a1, 0.0], priority=10)
    zm = nearest(y, line.copper / 2)
    for name, sign, lo, hi in (
        ("source_out", 1.0, radius + r1, radius + g0),
        ("source_in", -1.0, radius - g0, radius - r1),
    ):
        source = csx.AddExcitation(name, exc_type=0, exc_val=[sign, 0.0, 0.0])
        source.AddBox([lo, float(a[2]), zm], [hi, float(a[2]), zm])
    gaps = {
        "v_out": (radius + r1, radius + g0),
        "v_out_in": (radius + r0, radius + s1),
        "v_in": (radius - r1, radius - g0),
        "v_in_in": (radius - r0, radius - s1),
    }
    av = float(a[1])
    for name, (start, stop) in gaps.items():  # each from its rail to the copper beside it
        csx.AddProbe(name, p_type=0).AddBox([start, av, zm], [stop, av, zm])
    rd, zd = dual(r), dual(y)
    ac = float(a[-1]) - da / 2
    loops = {}
    for name, sign in (("i_out", 1.0), ("i_in", -1.0)):
        lo, hi = sorted((radius + sign * line.edges[0], radius + sign * line.edges[1]))
        loops[name] = (
            [nearest(rd, lo), ac, nearest(zd, -line.board / 2)],
            [nearest(rd, hi), ac, nearest(zd, 3 * line.copper)],
        )
        csx.AddProbe(name, p_type=1, norm_dir=1).AddBox(*loops[name])
    dump = csx.AddDump("hcut", dump_type=11, dump_mode=0, file_type=1, frequency=list(f_dump))
    dump.AddBox([float(r[0]), ac, float(y[0])], [float(r[-1]), ac, float(y[-1])])
    began = time.time()
    fdtd.Run(str(path), cleanup=True, verbose=0, engine=ENGINE)
    return {
        "x": x,
        "r": r,
        "a": a,
        "y": y,
        "x1": x1,
        "alpha_voltage": av,
        "alpha_current": ac,
        "loops": loops,
        "wall_s": time.time() - began,
        "cells": (len(r) - 1) * (len(a) - 1) * (len(y) - 1),
        "steps": steps,
    }


def read_probe(path: Path, name: str) -> tuple[Arr, Arr]:
    data = np.loadtxt(path / name, comments="%")
    return data[:, 0], data[:, 1]


def pencil(t: Arr, y: Arr, order: int) -> tuple[CArr, CArr, float]:
    """Poles s (1/s) and complex amplitudes c of y(t) ~ sum c_k exp(s_k (t - t_0)) on uniform samples, by the matrix
    pencil (Hua and Sarkar 1990): the Hankel matrix's leading `order` right singular vectors, shifted, give the
    poles, and least squares the amplitudes; with the fit's rms residual over the signal's rms."""
    n = len(y)
    pencil_length = n // 3
    hankel = np.lib.stride_tricks.sliding_window_view(y, pencil_length + 1)
    _, _, vh = np.linalg.svd(hankel, full_matrices=False)
    v = vh[:order].conj().T
    zs = np.linalg.eigvals(np.linalg.pinv(v[:-1]) @ v[1:])
    dt = float(t[1] - t[0])
    s: CArr = np.asarray(np.log(zs.astype(np.complex128)) / dt, dtype=np.complex128)
    basis = zs[None, :] ** np.arange(n)[:, None]
    c: CArr = np.linalg.lstsq(basis, y.astype(np.complex128), rcond=None)[0]
    residual = float(np.sqrt(np.mean(np.abs(basis @ c - y) ** 2)) / np.sqrt(np.mean(y**2)))
    return s, c, residual


def ringdown(signals: dict[str, tuple[Arr, Arr]], reference: str, start: float, f_guess: float) -> dict[str, Any]:
    """The resonance nearest f_guess in the ring-down after `start` (s): its frequency and Q from the reference
    signal's poles, and every signal's complex amplitude at that pole (each signal fitted with the reference's poles
    on its own time samples)."""
    t, y = signals[reference]
    keep = t >= start
    s, c, residual = pencil(t[keep], y[keep], 8)
    upper = np.flatnonzero(s.imag > 0)
    k = int(upper[np.argmin(np.abs(s.imag[upper] / (2 * math.pi) - f_guess))])
    pole = complex(s[k])
    amplitudes: dict[str, list[float]] = {}
    for name, (ti, yi) in signals.items():
        sel = ti >= start
        tt = ti[sel] - ti[sel][0]
        basis = np.exp(np.outer(tt, s))
        ci = np.linalg.lstsq(basis, yi[sel].astype(np.complex128), rcond=None)[0]
        amp = complex(ci[k]) * complex(np.exp(pole * (t[keep][0] - ti[sel][0])))  # at the reference's first sample
        amplitudes[name] = [amp.real, amp.imag]
    others = [float(abs(s[j].imag) / (2 * math.pi)) for j in upper if j != k]
    return {
        "f": pole.imag / (2 * math.pi),
        "Q": pole.imag / (2 * abs(pole.real)),
        "residual": residual,
        "amplitudes": amplitudes,
        "other_poles_GHz": [v / 1e9 for v in others],
        "samples": int(keep.sum()),
        "amplitude": float(2 * abs(c[k])),
    }


def pencil_check(seed: int = 7) -> dict[str, float]:
    """Calibration of the ring-down fit (rule 5): a damped tone at 914.3 MHz, Q 60, with a weaker one at 1.8 GHz and
    white noise at 10^-3 of the tone, on the probes' 88 ps samples over 23 ns, recovered to these relative errors."""
    rng = np.random.default_rng(seed)
    t: Arr = np.arange(260, dtype=np.float64) * 88.3e-12
    f0, q0 = 914.3e6, 60.0
    y: Arr = np.exp(-math.pi * f0 * t / q0) * np.cos(2 * math.pi * f0 * t + 0.4)
    y = y + 0.05 * np.exp(-math.pi * 1.8e9 * t / 50) * np.cos(2 * math.pi * 1.8e9 * t) + 1e-3 * rng.standard_normal(260)
    got = ringdown({"y": (t, y)}, "y", 0.0, 915e6)
    return {"f_error": got["f"] / f0 - 1, "Q_error": got["Q"] / q0 - 1, "residual": got["residual"]}


def read_dump(path: Path, index: int = 0) -> CArr:
    """The H dump's field (3, n1, n2, n3) at its index-th frequency, complex, the native Yee components."""
    h5py = importlib.import_module("h5py")
    with h5py.File(path / "hcut.h5", "r") as f:
        data: CArr = np.asarray(f[f"FieldData/FD/f{index}"][()], dtype=np.complex128)
    return data


def straight_copper_q(line: RailLine, x: Arr, y: Arr, x1: float, h: CArr) -> dict[str, float]:
    """Q_c of the straight line from its H dump at the current antinode, read at the Yee places of H_x (x_i,
    y_j+1/2) and H_y (x_i+1/2, y_j) on the half model's lines (the ghost line behind the centre plane carries no
    field); and the share of |H|^2 the dump puts inside the copper, zero when the places are read right."""
    gx, gy = np.meshgrid(x, y, indexing="ij")
    m = section_masks(line, np.abs(gx), gy, half_wall=x1)
    hx, hy, hz = h[0, :, :-1, 0], h[1, :-1, :, 0], h[2, :-1, :-1, 0]
    wx, wy = node_widths(x), node_widths(y)
    dx, dy = np.diff(x), np.diff(y)
    ex = np.abs(hx) ** 2 * wx[:, None] * dy[None, :]
    ey = np.abs(hy) ** 2 * dx[:, None] * wy[None, :]
    ez = np.abs(hz) ** 2 * dx[:, None] * dy[None, :]
    energy = float(ex.sum() + ey.sum() + ez.sum())
    metal = ~m.free
    leak = float((ex[metal[:, :-1] & metal[:, 1:]].sum() + ey[metal[:-1, :] & metal[1:, :]].sum()) / energy)
    return {"Qc": copper_q(energy, face_sum(x, y, m, hx, hy)), "inside": leak}


def ring_copper_q(line: RailLine, x: Arr, r: Arr, y: Arr, x1: float, h: CArr) -> dict[str, float]:
    """Q_c of the ring from its H dump at the current antinode: H_r at (r_i, z_j+1/2), H_z at (r_i+1/2, z_j), every
    area and face weighted by its radius; and the share of |H|^2 inside the copper."""
    gx, gy = np.meshgrid(x, y, indexing="ij")
    m = section_masks(line, np.abs(gx), gy, half_wall=x1)
    hr, hz, ha = h[0, :, 0, :-1], h[2, :-1, 0, :], h[1, :-1, 0, :-1]
    rd = dual(r)
    wr, wz = node_widths(r), node_widths(y)
    dr, dz = np.diff(r), np.diff(y)
    er = np.abs(hr) ** 2 * (wr * r)[:, None] * dz[None, :]
    ez = np.abs(hz) ** 2 * (dr * rd)[:, None] * wz[None, :]
    ea = np.abs(ha) ** 2 * (dr * rd)[:, None] * dz[None, :]
    energy = float(er.sum() + ez.sum() + ea.sum())
    metal = ~m.free
    leak = float((er[metal[:, :-1] & metal[:, 1:]].sum() + ez[metal[:-1, :] & metal[1:, :]].sum()) / energy)
    return {"Qc": copper_q(energy, face_sum(x, y, m, hr, hz, wr=r, wrd=rd)), "inside": leak}


def amplitude(fit: dict[str, Any], name: str) -> complex:
    re, im = fit["amplitudes"][name]
    return complex(re, im)


def trace(signal: tuple[Arr, Arr], every: int = 2) -> dict[str, list[float]]:
    """A probe's signal for the figure, every other sample, over its largest value."""
    t, v = signal
    return {"t_ns": (t[::every] * 1e9).tolist(), "v": (v[::every] / float(np.max(np.abs(v)))).tolist()}


def stage_straight(
    line: RailLine, grading: Grading, section: dict[str, Any], sim_ns: float, keep: bool, tag: str
) -> dict[str, Any]:
    """The straight line between magnetic walls, half a guided wave long at 915 MHz by the decomposed index (the
    converged section's): its resonance, ring-down Q (the laminate's loss), Z0 from the probes, Q_c from the H dump
    at the resonance its own grid's 2-D twin predicts."""
    length = pt.C0 / (2 * F915 * section["built"]["index"]) / kg.MM
    f_twin = pt.C0 / (2 * length * kg.MM * section[f"fdtd {tag}"]["index"])
    cells = 2 * round(length / 3.0)  # even, about 1.5 mm a cell
    steps = steps_for(grading, length / cells, sim_ns)
    path = SCRATCH / f"straight_{tag}"
    shutil.rmtree(path, ignore_errors=True)
    path.mkdir(parents=True)
    built = straight_model(line, grading, length, cells, steps, sim_ns, path, [f_twin])
    signals = {name: read_probe(path, name) for name in ("v0", "v0_in", "vm", "vl", "vl_in", "ic")}
    step = quantization(line, 2 * float(signals["ic"][0][0]))  # a current probe reads half a step after the volts
    start = 2 * 9 / (2 * math.pi * FC) + 1e-9
    fit = ringdown(signals, "v0", start, F915)
    v0, vl, vin, ic = (amplitude(fit, k) for k in ("v0", "vl", "v0_in", "ic"))
    beta = math.pi / length
    zv = built["probes"]["v0"]
    # V0 cos(beta z) across a gap, (V0 / Z0) sin(beta z) in both rails together
    z0 = abs(v0) / math.cos(beta * zv) / (2 * abs(ic))
    copper = straight_copper_q(line, built["x"], built["y"], built["x1"], read_dump(path))
    out = {
        "grading": asdict(grading),
        "length_mm": length,
        "dz_mm": length / cells,
        "cells": built["cells"],
        "f_twin": f_twin,
        "steps": round(sim_ns * 1e-9 / step["dt"]),
        "quantization": step,
        "wall_s": built["wall_s"],
        "f": fit["f"],
        "index": pt.C0 / (2 * length * kg.MM * fit["f"]),
        "Qd": fit["Q"],
        "z0": z0,
        "Qc": copper["Qc"],
        "inside": copper["inside"],
        "Q": 1 / (1 / fit["Q"] + 1 / copper["Qc"]),
        "residual": fit["residual"],
        "other_poles_GHz": fit["other_poles_GHz"],
        "samples": fit["samples"],
        "end_over_start": abs(vl) / abs(v0),
        "phase_end_deg": math.degrees(float(np.angle(vl / v0))),
        "inner_over_outer": abs(vin) / abs(v0),
        "node_over_antinode": abs(amplitude(fit, "vm")) / abs(v0),
        "node_expected": abs(math.cos(beta * built["probes"]["vm"])) / math.cos(beta * zv),
        "trace": trace(signals["v0"]),
    }
    if not keep:
        shutil.rmtree(path, ignore_errors=True)
    return out


def stage_ring(
    line: RailLine, grading: Grading, section: dict[str, Any], sim_ns: float, keep: bool, tag: str
) -> dict[str, Any]:
    """R1's ring of the same section, one guided wave round at 915 MHz by the decomposed index (m = 1). Its two
    rails resonate apart (two_rail_ring): each mode's frequency, ring-down Q, the rails' volts (inner over outer),
    from the ring-down of the rail where it is strongest, and Q_c from the H dump at the frequency its own grid's
    2-D twin predicts for it."""
    radius = pt.C0 / (F915 * section["built"]["index"]) / kg.MM / (2 * math.pi)
    twin = two_rail_ring(section[f"fdtd {tag}"]["rails"], line, radius, line.tand)
    lines = round(math.pi / 2 * radius / 1.5) + 2  # about 1.5 mm of arc a cell on the centre line
    inner_arc = (radius - line.shield) * math.pi / (2 * lines - 3)
    steps = steps_for(grading, inner_arc, sim_ns)
    path = SCRATCH / f"ring_{tag}"
    shutil.rmtree(path, ignore_errors=True)
    path.mkdir(parents=True)
    built = ring_model(line, grading, radius, lines, steps, sim_ns, path, [mode["f"] for mode in twin])
    names = ("v_out", "v_out_in", "v_in", "v_in_in", "i_out", "i_in")
    signals = {name: read_probe(path, name) for name in names}
    step = quantization(line, 2 * float(signals["i_out"][0][0]))  # a current probe reads half a step after the volts
    start = 2 * 9 / (2 * math.pi * FC) + 1e-9
    found = []
    for k, (mode, reference) in enumerate(zip(twin, ("v_out", "v_in"), strict=True)):
        fit = ringdown(signals, reference, start, mode["f"])
        vo, vi = amplitude(fit, "v_out"), amplitude(fit, "v_in")
        copper = ring_copper_q(line, built["x"], built["r"], built["y"], built["x1"], read_dump(path, k))
        ratio = vi / vo
        found.append(
            {
                "f": fit["f"],
                "index_on_centre": pt.C0 / (2 * math.pi * radius * kg.MM * fit["f"]),
                "Qd": fit["Q"],
                "Qc": copper["Qc"],
                "inside": copper["inside"],
                "Q": 1 / (1 / fit["Q"] + 1 / copper["Qc"]),
                "inner_over_outer": abs(ratio) * (1.0 if ratio.real >= 0 else -1.0),
                "phase_inner_deg": math.degrees(float(np.angle(ratio))),
                "residual": fit["residual"],
                "reference": reference,
                "f_twin": mode["f"],
                "other_poles_GHz": fit["other_poles_GHz"],
                "outer_gaps_inner_over_outer": abs(amplitude(fit, "v_out_in")) / abs(vo),
            }
        )
    out = {
        "grading": asdict(grading),
        "radius_mm": radius,
        "length_mm": 2 * math.pi * radius,
        "alpha_lines": lines,
        "cells": built["cells"],
        "steps": round(sim_ns * 1e-9 / step["dt"]),
        "quantization": step,
        "wall_s": built["wall_s"],
        "modes": found,
        "twin": twin,
        "samples": int((signals["v_out"][0] >= start).sum()),
        "trace": trace(signals["v_out"]),
        "trace_inner": trace(signals["v_in"]),
    }
    if not keep:
        shutil.rmtree(path, ignore_errors=True)
    return out


# ---- the run --------------------------------------------------------------------------------------------------------
STAGES = ("section", "straight", "ring")


def load(out: Path) -> dict[str, Any]:
    try:
        loaded: dict[str, Any] = json.loads((out / "results.json").read_text())
        return loaded
    except (OSError, ValueError):
        return {}


def save(out: Path, results: dict[str, Any]) -> None:
    out.mkdir(parents=True, exist_ok=True)
    tmp = out / "results.json.tmp"
    tmp.write_text(json.dumps(lab.plain(results), indent=1) + "\n")
    tmp.replace(out / "results.json")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument("command", nargs="?", default="all", choices=("all", *STAGES, "report"))
    parser.add_argument("--quick", action="store_true", help="coarse meshes, a short ring-down, out to the scratch")
    parser.add_argument("--redo", action="store_true", help="recompute the stages asked for even when stored")
    parser.add_argument("--keep", action="store_true", help="keep openEMS's folders")
    args = parser.parse_args(argv)
    out = SCRATCH / "quick" if args.quick else OUT
    results = load(out)
    chosen = list(STAGES) if args.command == "all" else ([] if args.command == "report" else [args.command])
    sim_ns = 12.0 if args.quick else SIM_NS
    for stage in chosen:
        if stage in results and not args.redo and (stage != "section" or "built" in results[stage]):
            print(f"{stage}: stored, kept (--redo recomputes)", flush=True)
            continue
        began = time.time()
        print(f"{stage}: running", flush=True)
        if stage == "section":
            results["section"] = stage_section(args.quick)
            results["pencil"] = pencil_check()
        else:
            meshes = (("quick", QUICK),) if args.quick else (("fine", FDTD_FINE), ("coarse", FDTD_COARSE))
            if stage == "ring":
                meshes = meshes[:1]
            run = stage_straight if stage == "straight" else stage_ring
            results[stage] = {
                tag: run(R1, grading, results["section"], sim_ns, args.keep, tag) for tag, grading in meshes
            }
        results.setdefault("seconds", {})[stage] = time.time() - began
        save(out, results)
        print(f"{stage}: done in {time.time() - began:.0f} s", flush=True)
    if "section" in results:
        report(results, out, args.quick)
    return 0


def at_tand(row: dict[str, Any], factor: float = 1.0) -> float:
    """A full-wave ring-down Q at the laminate's constant tan d: openEMS's conductivity gives tan d (915 MHz / f),
    so the dielectric Q it reads at f is (f / 915 MHz) times the constant-tan-d one; and the run's loss is `factor`
    times the given one (quantization), so its Q times `factor` is the given loss's."""
    return float(row["Qd"]) * factor * F915 / float(row["f"])


def finite(value: Any) -> bool:
    """A stored number that is there: results.json keeps a NaN as null."""
    return isinstance(value, float | int) and math.isfinite(value)


def combine(*qs: float) -> float:
    return 1 / sum(1 / q for q in qs)


def verdict(deviation: float, band: float) -> str:
    return f"{deviation:+.2%}: {'pass' if abs(deviation) <= band else 'FAIL'} (band {band:.1%})"


def face_limit(qs: Sequence[float]) -> dict[str, float]:
    """The face integral's Q_c on cells halving at the copper, extrapolated to no cell: the corners' fields go as
    r^(-1/3), their squares as r^(-2/3), so the error falls about 2^(-1/3) per halving. Richardson at that ratio from
    the last two, and Aitken's delta-squared on the last three (the ratio the last three show, 2^(-p)); the one check
    of the copper's loss here that does not use Wheeler's rule."""
    if len(qs) < 3:
        raise ValueError(f"an extrapolation takes three cells or more, not {len(qs)}")
    a, b, c = qs[-3:]
    ratio = 2 ** (-1 / 3)
    seen = (c - b) / (b - a)
    return {
        "richardson": c + (c - b) * ratio / (1 - ratio),
        "aitken": c - (c - b) ** 2 / ((c - b) - (b - a)),
        "ratio": seen,
        "order": -math.log2(seen),
    }


def dr0(results: dict[str, Any]) -> dict[str, Any]:
    """DR-0's comparisons, decomposed against full wave, from the stored stages: the straight line on each mesh, and
    R1's ring against the single-mode and the two-rail decomposed models. The laminate's Q from the ring-down at the
    given tan d, its single-precision rounding (quantization) taken out; the copper's Q read raw from the H dump and
    calibrated by the face integral's own error on the same cells (the static twin against Wheeler's rule)."""
    sec = results["section"]
    built = sec["built"]
    out: dict[str, Any] = {"straight": {}, "ring": {}}
    faces = [sec[k]["Qc_faces"] for k in sorted((k for k in sec if k.startswith("faces t/")), key=lambda k: int(k[8:]))]
    if len(faces) >= 3:
        out["copper"] = face_limit(faces) | {"faces": faces, "wheeler": built["Qc"]}
    for tag, row in results.get("straight", {}).items():
        twin = sec[f"fdtd {tag}"]
        factor = row.get("quantization", {}).get("factor", 1.0)
        qd = at_tand(row, factor)
        calibration = built["Qc"] / twin["Qc_faces"]
        out["straight"][tag] = {
            "f_dev": row["f"] / F915 - 1,
            "Q": combine(qd, row["Qc"]),
            "Q_dev": combine(qd, row["Qc"]) / built["Q"] - 1,
            "Q_cal": combine(qd, row["Qc"] * calibration),
            "Q_cal_dev": combine(qd, row["Qc"] * calibration) / built["Q"] - 1,
            "Qd_tand": qd,
            "Qd_tand_rounded": at_tand(row),
            "quantization": factor,
            "calibration": calibration,
            "index_vs_twin": row["index"] / twin["index"] - 1,
            "Qd_vs_twin": qd / twin["Qd"] - 1,
            "z0_vs_twin": row["z0"] / twin["z0"] - 1,
            "Qc_vs_twin": row["Qc"] / twin["Qc_faces"] - 1,
            "index_vs_decomposed": row["index"] / built["index"] - 1,
            "z0_vs_decomposed": row["z0"] / built["z0"] - 1,
            "Qd_vs_decomposed": qd / built["Qd"] - 1,
            "Qc_vs_decomposed": row["Qc"] / built["Qc"] - 1,
        }
    for tag, ring in results.get("ring", {}).items():
        model = two_rail_ring(built["rails"], R1, ring["radius_mm"], R1.tand)
        factor = ring.get("quantization", {}).get("factor", 1.0)
        rows = []
        for mode, fit, twin in zip(model, ring["modes"], ring["twin"], strict=True):
            qd = at_tand(fit, factor)
            calibration = mode["Qc"] / twin["Qc_faces"]
            rows.append(
                {
                    "model": mode,
                    "f_dev": fit["f"] / mode["f"] - 1,
                    "f_dev_single": fit["f"] / F915 - 1,
                    "Q": combine(qd, fit["Qc"]),
                    "Q_dev": combine(qd, fit["Qc"]) / mode["Q"] - 1,
                    "Q_cal": combine(qd, fit["Qc"] * calibration),
                    "Q_cal_dev": combine(qd, fit["Qc"] * calibration) / mode["Q"] - 1,
                    "Qd_tand": qd,
                    "Qd_tand_rounded": at_tand(fit),
                    "quantization": factor,
                    "calibration": calibration,
                    "ratio_dev": fit["inner_over_outer"] - mode["inner_over_outer"],
                    "f_vs_twin": fit["f"] / twin["f"] - 1,
                    "Qd_vs_twin": qd / twin["Qd"] - 1,
                    "Qc_vs_twin": fit["Qc"] / twin["Qc_faces"] - 1,
                }
            )
        out["ring"][tag] = rows
    return out


def report(results: dict[str, Any], out: Path, quick: bool) -> None:
    """run.log, the DR-0 verdicts and the figure, from the stored results alone."""
    LOG.clear()
    sec = results["section"]
    verdicts = dr0(results)
    b, w = sec["built"], sec["written"]
    nb, nw = sec.get("null built", b), sec.get("null written", w)
    line = RailLine(**sec["line"])
    c1, c2 = line.edges
    lam_mm = pt.C0 / F915 / kg.MM
    say("twin bench, stage 0: R1's rail line in openEMS against its 2-D section (docs/KNOTTED_GUIDE.md 8.1, DR-0)")
    say(f"registered before any run (DR-0): the resonance within {DR0_F:.1%} and Q within {DR0_Q:.0%} of full wave")
    say(
        f"run: {'quick' if quick else 'full'}; engine {ENGINE}; seconds "
        + ", ".join(f"{k} {v:.0f}" for k, v in results.get("seconds", {}).items())
    )
    say("\n== the section as modelled ==")
    say(
        f"  one {line.board} mm FR-4 board (er {line.er}, tan d {line.tand}) with {line.copper * 1e3:.0f} um copper: a "
        f"centre strip and rails of {2 * line.spans[0][1]:.2f} and {line.spans[1][1] - line.spans[1][0]:.2f} mm of "
        f"copper between {line.gap:.2f} mm gaps centred on ELE-1's edges c1 = {c1:.2f} and c2 = {c2:.2f} mm (a = "
        f"{line.strip:.2f}, b = {line.rail:.2f} mm); a ground plane under the board, via walls on the centre plane and "
        f"{line.fence:.1f} mm beyond the outer gap, a grounded shield {line.shield:.0f} mm out"
    )
    say(
        f"  the null {nb['null_mm']:.3f} mm above the copper's top face (ELE-1's gapless sqrt(c1 c2) = "
        f"{line.null_gapless:.3f} mm over the plane), its gradient {nb['gradient_ele1']:.3f} of ELE-1's: q 0.5 takes "
        f"{50 / nb['gradient_ele1']:.1f} V where ELE-1 says 50 V"
    )
    say(
        f"  left floating, the strip takes {b['floating']:.3f} of the rails' volts on this board and "
        f"{w['floating']:.3f} on the board as written (no plane); the null then sits "
        + (f"{b['null_floating_mm']:.3f} mm" if finite(b.get("null_floating_mm")) else "nowhere")
        + " above the copper on this board and "
        + (f"{w['null_floating_mm']:.3f} mm" if finite(w.get("null_floating_mm")) else "nowhere")
        + " on that one: the via wall that holds the strip at 0 V is part of the design. Floating, the line has two "
        f"modes, n_eff {b['mode_index'][0]:.4f} and {b['mode_index'][1]:.4f}"
    )
    say("\n== calibration (rule 5) ==")
    for row in sec["stripline"]:
        say(
            f"  Wheeler's rule on stripline {row['case'][0]}, {row['case'][1]}, {row['case'][2]} mm (width, "
            f"thickness, spacing): Z0 {row['z0']:.1f} ohm, alpha_c {row['ratio']:.3f} of Pozar's closed form"
        )
    say(
        "  Pozar's eq. 3.198 is itself Wheeler's rule applied to a closed-form Z0(t): the stripline agreement shows "
        "the receded solve's numerical derivative right, so the copper's loss is held to Wheeler's closed form only"
    )
    pc = results["pencil"]
    say(
        f"  the ring-down fit (matrix pencil) on a damped tone, Q 60, with a second tone and 1e-3 noise: f to "
        f"{pc['f_error']:+.1e}, Q to {pc['Q_error']:+.1e}"
    )
    faces = sorted(
        ((v["grading"]["fine"], v["Qc_faces"]) for k, v in sec.items() if k.startswith("faces")), reverse=True
    )
    if len(faces) > 1:
        say(
            "  the face integral of the static field, R_s |E_n|^2 on the edge beside each copper face (the 3-D "
            "integral's quasi-static twin), against the cell at the copper: "
            + ", ".join(f"Q_c {q:.1f} at {cell * 1e3:.2f} um" for cell, q in faces)
            + f"; Wheeler's rule {b['Qc']:.1f}. The face integral converges slowly, the corners' fields going as "
            "r^(-1/3), so a copper loss read off a grid's faces is calibrated by its static twin on those cells"
        )
    copper = verdicts.get("copper")
    if copper:
        say(
            f"  extrapolated to no cell: its differences fall {copper['ratio']:.2f} per halving at the finest "
            f"(2^(-1/3) is 0.79, the corners' r^(-2/3)); Richardson at 2^(-1/3) puts its limit at "
            f"{copper['richardson']:.1f} and Aitken's on the last three at {copper['aitken']:.1f}, "
            f"{1 - copper['richardson'] / copper['wheeler']:.1%} and {1 - copper['aitken'] / copper['wheeler']:.1%} "
            f"under Wheeler's rule's {copper['wheeler']:.1f}: the one check of the copper's loss here that does not "
            "use the rule; a surface-impedance solve (Palace) is DR-0's item for it"
        )
    steady = [abs(sec[k]["Qc"] / b["Qc"] - 1) for k in ("coarse 0.5", "shield 50") if k in sec]
    if steady:
        say(f"  Wheeler's rule moves {max(steady):.2%} with far cells of 0.5 mm instead of 1 mm or the shield at 50 mm")
    say("\n== the decomposed model: the 2-D section on the study's converged grading ==")
    for name, row in (
        ("even mode (the rails in phase, as the guide drives it)", b),
        ("odd mode (the rails apart)", b["odd"]),
    ):
        say(
            f"  {name}: C' {row['C_pF']:.1f} pF/m, n_eff {row['index']:.4f}, Z0 {row['z0']:.2f} ohm, the laminate "
            f"{row['share']:.3f} of the field energy: Q_d {row['Qd']:.1f}; Q_c {row['Qc']:.1f} by Wheeler's rule; Q "
            f"{row['Q']:.1f}"
        )
    rr = b["rails"]
    say(
        f"  the rails as two conductors: C12/C11 {rr['C'][0][1] / rr['C'][0][0]:+.4f} with the laminate, "
        f"{rr['Ca'][0][1] / rr['Ca'][0][0]:+.4f} in air: the strip and its via wall screen one rail from the other"
    )
    shunt = b["z0"] * b["Q"] / (2 * math.pi)
    say(
        f"  m = 1 on this line is {lam_mm / b['index']:.1f} mm round; R = Z0 Q / 2 pi = {shunt:.1f} ohm, so 50 V takes "
        f"{50**2 / (2 * shunt):.1f} W (RFG-3; the laminate's Q_d {b['Qd']:.0f} and the low Z0 set it)"
    )
    say(
        f"  the board as written (no plane, the strip held at 0 V): n_eff {w['index']:.4f}, Z0 {w['z0']:.2f} "
        f"ohm, share "
        f"{w['share']:.3f}, Q_d {w['Qd']:.1f}, Q_c {w['Qc']:.1f}, Q {w['Q']:.1f}; m = 1 takes "
        f"{lam_mm / w['index']:.1f} "
        f"mm (section 3.6's estimate sqrt((1 + er)/2) = {math.sqrt((1 + line.er) / 2):.4f}: "
        f"{lam_mm / math.sqrt((1 + line.er) / 2):.1f} mm); its null {nw['null_mm']:.3f} mm up, gradient "
        f"{nw['gradient_ele1']:.3f} of ELE-1's"
    )
    for key in ("gap 0.10", "gap 0.20"):
        if key in sec:
            g = sec[key]
            say(
                f"  lever, the gaps at {key[4:]} mm: n_eff {g['index']:.4f}, Z0 {g['z0']:.2f} ohm, Q_c "
                f"{g['Qc']:.1f}, Q "
                f"{g['Q']:.1f}, the null {g['null_mm']:.3f} mm up"
            )
    if "RO4003C 1.524" in sec:
        ro = sec["RO4003C 1.524"]
        shunt_ro = ro["z0"] * ro["Q"] / (2 * math.pi)
        say(
            f"  lever, the board RO4003C 1.524 mm (er {kg.RO4003C[0]}, tan d {kg.RO4003C[1]}): n_eff "
            f"{ro['index']:.4f}, Z0 {ro['z0']:.2f} ohm, Q_d {ro['Qd']:.0f}, Q_c {ro['Qc']:.1f}, Q {ro['Q']:.1f}: R "
            f"{shunt_ro:.0f} ohm, 50 V takes {50**2 / (2 * shunt_ro):.1f} W"
        )
    if "shield 50" in sec:
        s = sec["shield 50"]
        say(
            f"  the shield at 50 mm instead of 25: n_eff {s['index'] / b['index'] - 1:+.4%}, Z0 "
            f"{s['z0'] / b['z0'] - 1:+.4%}, Q {s['Q'] / b['Q'] - 1:+.4%}"
        )
    for tag, row in results.get("straight", {}).items():
        v = verdicts["straight"][tag]
        twin = sec[f"fdtd {tag}"]
        say(f"\n== the straight line between magnetic walls, full wave ({tag} mesh) ==")
        say(
            f"  {row['grading']['fine'] * 1e3:.0f} um at the copper growing {row['grading']['growth']} per mm, "
            f"{row['length_mm']:.2f} mm between the walls in {row['dz_mm']:.2f} mm cells: {row['cells'] / 1e6:.2f} M "
            f"cells, {row['steps'] / 1e6:.2f} M steps, {row['wall_s'] / 60:.1f} min; {row['samples']} ring-down "
            f"samples "
            f"after the pulse, fit residual {row['residual']:.1e}"
        )
        say(
            f"  resonance {row['f'] / 1e6:.3f} MHz: n_eff {row['index']:.5f} (its own cells' 2-D twin "
            f"{twin['index']:.5f}, {v['index_vs_twin']:+.3%}; the converged section {b['index']:.5f}, "
            f"{v['index_vs_decomposed']:+.3%})"
        )
        quant = row.get("quantization")
        if quant:
            say(
                f"  dt {quant['dt']:.5e} s: the laminate's damping per step, 1 - VV, is {quant['damping_steps']:.3f} "
                f"units of 2^-24 and openEMS's single precision holds {round(quant['damping_steps'])}: the run's "
                f"laminate loses {quant['factor']:.4f} of the given loss"
            )
        say(
            f"  Q_d {row['Qd']:.2f} by the ring-down, {v['Qd_tand_rounded']:.2f} at tan d {line.tand} as the run "
            f"rounded it, {v['Qd_tand']:.2f} with the rounding taken out (twin {twin['Qd']:.2f}, "
            f"{v['Qd_vs_twin']:+.2%}; converged {b['Qd']:.2f}, {v['Qd_vs_decomposed']:+.2%}); Z0 {row['z0']:.3f} ohm "
            f"(twin {twin['z0']:.3f}, {v['z0_vs_twin']:+.2%}; converged {b['z0']:.3f}, {v['z0_vs_decomposed']:+.2%})"
        )
        say(
            f"  Q_c {row['Qc']:.1f} from the H dump's |H_t|^2 over the copper (twin {twin['Qc_faces']:.1f}, "
            f"{v['Qc_vs_twin']:+.2%}; Wheeler's rule converged {b['Qc']:.1f}, {v['Qc_vs_decomposed']:+.1%}); the dump "
            f"puts {row['inside']:.1e} of |H|^2 inside the copper"
        )
        say(
            f"  checks: the far wall's volts {row['end_over_start']:.4f} of the near wall's at "
            f"{row['phase_end_deg']:.1f} "
            f"deg (m = 1's cos); the strip's gap {row['inner_over_outer']:.4f} of the ground's (one mode); the middle "
            f"{row['node_over_antinode']:.4f} of the end where cos puts {row['node_expected']:.4f}; other poles "
            + ", ".join(f"{f:.3f}" for f in row["other_poles_GHz"])
            + " GHz"
        )
        say(
            f"  DR-0, straight ({tag}): resonance, decomposed 915.000 MHz, full wave {row['f'] / 1e6:.3f} MHz: "
            + verdict(v["f_dev"], DR0_F)
        )
        say(
            f"  DR-0, straight ({tag}): Q, decomposed {b['Q']:.2f}, full wave {v['Q']:.2f} with the copper read off "
            f"these cells: {verdict(v['Q_dev'], DR0_Q)}; with that reading calibrated by its twin ("
            f"x{v['calibration']:.3f}), {v['Q_cal']:.2f}: {verdict(v['Q_cal_dev'], DR0_Q)}. Q_d sets most of Q; the "
            "full wave's copper is PEC, its loss read off the cells, and the calibration is Wheeler's rule again: "
            "this holds n_eff, Z0 and Q_d to full wave, and the copper's loss to Wheeler's closed form only"
        )
    for tag, ring in results.get("ring", {}).items():
        rows = verdicts["ring"][tag]
        say(f"\n== R1's ring, full wave ({tag} mesh, cylindrical: every edge on a mesh line, no staircase) ==")
        say(
            f"  centre-line radius {ring['radius_mm']:.3f} mm, {ring['length_mm']:.2f} mm round (m = 1 at 915 "
            f"MHz by the "
            f"even mode's n_eff); a quarter between a magnetic wall at alpha 0 and an electric wall at pi/2, "
            f"{ring['alpha_lines']} alpha lines: {ring['cells'] / 1e6:.2f} M cells, {ring['steps'] / 1e6:.2f} M steps, "
            f"{ring['wall_s'] / 60:.1f} min"
        )
        quant = ring.get("quantization")
        if quant:
            say(
                f"  dt {quant['dt']:.5e} s: the laminate's damping per step {quant['damping_steps']:.3f} units of "
                f"2^-24, held as {round(quant['damping_steps'])}: the run's laminate loses {quant['factor']:.4f} of "
                "the given loss"
            )
        say(
            "  the single-mode decomposed model puts m = 1 at 915.000 MHz with the rails in phase; full wave finds no "
            "such mode:"
        )
        for k, (fit, row) in enumerate(zip(ring["modes"], rows, strict=True)):
            say(
                f"    mode {'AB'[k]}: {fit['f'] / 1e6:.3f} MHz ({row['f_dev_single']:+.2%} from 915), the inner "
                f"rail at "
                f"{fit['inner_over_outer']:+.3f} of the outer ({fit['phase_inner_deg']:+.1f} deg); Q_d {fit['Qd']:.2f} "
                f"({row['Qd_tand']:.2f} at tan d {line.tand}, the rounding out), Q_c {fit['Qc']:.1f} (H dump), Q "
                f"{row['Q']:.2f}; each "
                f"rail's two gaps {fit['outer_gaps_inner_over_outer']:.4f} apart; fit residual {fit['residual']:.1e}"
            )
        cbar = (c1 + c2) / 2
        say(
            f"  why: the rails' centre lines are {2 * math.pi * (ring['radius_mm'] - cbar):.1f} and "
            f"{2 * math.pi * (ring['radius_mm'] + cbar):.1f} mm round (-+{cbar / ring['radius_mm']:.2%}), while the "
            f"strip's via wall leaves them coupled by C12/C11 {rr['C'][0][1] / rr['C'][0][0]:+.3f} and their even and "
            f"odd indices {b['index']:.4f} and {b['odd']['index']:.4f} differ by "
            f"{(b['index'] - b['odd']['index']) / b['index']:.2%}: curvature mixes the even mode with the odd as "
            "strongly as the section separates them, so each ring mode lives mostly on one rail"
        )
        say(
            "  the two-rail decomposed model (the rails' 2 x 2 line matrices from the converged section, each rail per "
            "radian on its own radius):"
        )
        for k, row in enumerate(rows):
            mode = row["model"]
            say(
                f"    mode {'AB'[k]}: {mode['f'] / 1e6:.3f} MHz, the inner rail at "
                f"{mode['inner_over_outer']:+.3f} of the "
                f"outer, Q_d {mode['Qd']:.2f}, Q_c {mode['Qc']:.1f} (Wheeler), Q {mode['Q']:.2f}; full wave "
                f"against it: f {row['f_dev']:+.3%}, the rails' ratio {row['ratio_dev']:+.3f}, Q_d "
                f"{row['Qd_tand'] / mode['Qd'] - 1:+.2%}; on its own cells' twin: f {row['f_vs_twin']:+.3%}, Q_d "
                f"{row['Qd_vs_twin']:+.2%}, Q_c {row['Qc_vs_twin']:+.2%}"
            )
        for k, row in enumerate(rows):
            say(
                f"  DR-0, ring mode {'AB'[k]}, single-mode model: resonance "
                + verdict(row["f_dev_single"], DR0_F)
                + "; the rails not in phase"
            )
            say(
                f"  DR-0, ring mode {'AB'[k]}, two-rail model: resonance {verdict(row['f_dev'], DR0_F)}; Q, the copper "
                f"off these cells: {verdict(row['Q_dev'], DR0_Q)}; calibrated (x{row['calibration']:.3f}): "
                f"{verdict(row['Q_cal_dev'], DR0_Q)}"
            )
        say(
            "  Hypothesis (read from these results, not yet tested): the null needs the rails in phase, so R1's ring "
            "needs them held together against the curvature: bridges joining the rails across the strip at N >= 3 "
            "stations a lap (N not dividing 2m), the rails then slipping apart in phase by at most 2 pi m (c1 + c2) / "
            f"(2 R N) midway between bridges, {cbar / ring['radius_mm'] * 2 * math.pi / 6:.3f} rad at N = 6; or two "
            "compiled feeds that drive both ring modes into the in-phase pattern at 915 MHz (section 5.3), at the "
            "power the off-resonant mode costs"
        )
    say("\n== what limits the agreement ==")
    say(
        "  staircasing: none; the straight line's copper and the ring's (in r and z, the cylindrical mesh) lie on mesh "
        "lines, and the magnetic walls sit on dual lines placed at the centre plane and the line's ends "
        "(openEMS puts a "
        "magnetic wall half a cell inside the boundary)"
    )
    say(
        "  mesh: the FDTD's cells at the copper set n_eff, Z0 and the face integral; its own 2-D twin on the "
        "same cells "
        "carries the same error (the full wave adds the rest)"
    )
    say(
        "  ports: none; a soft source and probes, the line closed by mirror walls, so no end correction enters: the "
        "two walls' volts agree to 1e-6 and the middle reads what cos puts there"
    )
    say(
        "  precision: openEMS keeps its update coefficients in single precision, so a lossy cell's damping per step is "
        "a whole number of 2^-24 and the laminate's loss is off by up to 2^-25 / (dt sigma / eps) of itself: 2.7 % "
        "for FR-4 at the fine mesh's dt, growing as tan d times dt falls (low-loss laminates, finer cells); computed "
        "from dt and taken out above, it is a limit to register for every lossy openEMS run"
    )
    say(
        "  the copper: openEMS holds PEC copper of its real thickness; its loss is the surface integral of the "
        "full-wave H, which on 5 to 10 um cells reads the copper's Q high (its static twin shows by how much), and the "
        "via walls and the shield carry none in either model; the shield's PEC carries no radiation either. The "
        "copper's treated face, its roughness of the order of the 2.2 um skin depth, takes Hammerstad's factor on its "
        "share of the loss (array16.guide.sections, Section.rough), and the unbacked board as R1 draws it is the next "
        "full-wave ring"
    )
    save(out, {**results, "dr0": verdicts})
    figure(results, verdicts, out)
    (out / "run.log").write_text("\n".join(LOG) + "\n")


def even_field(line: RailLine, grading: Grading) -> tuple[Arr, Arr, Arr]:
    """The even mode's static potential (rails at 1 V, with the laminate) on the half line's grid, for the figure."""
    x, y = section_lines(line, grading)
    gx, gy = np.meshgrid(x, y, indexing="ij")
    cy = (y[1:] + y[:-1]) / 2
    eps = np.where(((cy < 0) & (cy > -line.board))[None, :], line.er, 1.0) * np.ones((len(x) - 1, 1))
    m = section_masks(line, gx, gy)
    phi: Arr = kg.laplace2d(x, y, eps, m.rail | m.ground | m.centre, np.where(m.rail, 1.0, 0.0), tol=1e-10)
    return x, y, phi


def figure(results: dict[str, Any], verdicts: dict[str, Any], out: Path) -> None:
    """One figure: the section and its field, the ring-downs, the ring's two modes, DR-0's two comparisons, and the
    copper's Q by the face integral against the cell."""
    plt = kg._pyplot()
    sec = results["section"]
    b = sec["built"]
    line = RailLine(**sec["line"])
    fig = plt.figure(figsize=(16.5, 9.6))
    grid = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.28, top=0.86)
    straight_c, ring_a, ring_b = kg.SERIES[0], kg.SERIES[1], kg.SERIES[2]

    # (a) the section, its copper, and the even mode's equipotentials
    ax = fig.add_subplot(grid[0, 0])
    x, y, phi = even_field(line, Grading(fine=line.copper / 8))
    xf = np.concatenate([-x[:0:-1], x])
    phif = np.concatenate([phi[:0:-1], phi])
    ax.add_patch(kg._rect(-3.4, -line.board, 6.8, line.board, kg.LAMINATE))
    for a0, a1 in line.spans:
        for sign in (1, -1):
            lo, hi = sorted((sign * a0, sign * min(a1, 3.4)))
            ax.add_patch(kg._rect(lo, 0.0, hi - lo, line.copper * 4, kg.COPPER))
    ax.add_patch(kg._rect(-3.4, -line.board - line.copper * 4, 6.8, line.copper * 4, kg.COPPER))
    for xv in (0.0, line.fence_x, -line.fence_x):
        ax.plot([xv, xv], [-line.board, 0.0], color=kg.INK2, lw=1.4)
    ax.contour(xf, y, phif.T, levels=np.linspace(0.05, 0.95, 10), colors=[kg.BLUES[5]], linewidths=0.7)
    null = sec.get("null built", b)["null_mm"] + line.copper
    ax.plot([0.0], [null], marker="x", color=kg.INK, ms=8, mew=1.8)
    ax.annotate(
        f"null {null - line.copper:.3f} mm up",
        (0.0, null),
        (0.35, null + 0.3),
        fontsize=8.2,
        arrowprops={"arrowstyle": "-", "color": kg.MUTED, "lw": 0.8},
    )
    ax.set_xlim(-3.4, 3.4)
    ax.set_ylim(-1.85, 1.6)
    ax.set_aspect("equal")
    ax.set_xlabel("across the line (mm)")
    ax.set_ylabel("up from the board's top face (mm)")
    ax.set_title("(a) the section: rails in phase, 0.1 V equipotentials\n(copper drawn 4x thick; via walls in grey)")

    # (b) ring-downs
    ax = fig.add_subplot(grid[0, 1])
    st = next(iter(results.get("straight", {}).values()), None)
    rg = next(iter(results.get("ring", {}).values()), None)
    rtag = next(iter(results.get("ring", {})), "")
    if st:
        ax.plot(
            st["trace"]["t_ns"], st["trace"]["v"], color=straight_c, lw=0.8, label="straight line, rail volts at a wall"
        )
    if rg:
        ax.plot(rg["trace"]["t_ns"], rg["trace"]["v"], color=ring_a, lw=0.8, label="ring, outer rail at alpha 0")
        ax.plot(
            rg["trace_inner"]["t_ns"], rg["trace_inner"]["v"], color=ring_b, lw=0.8, label="ring, inner rail at alpha 0"
        )
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("volts over each trace's largest")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=7.6)
    ax.set_title("(b) the ring-downs after the 5.7 ns pulse")

    # (c) the ring's modes: frequency against the rails' ratio
    ax = fig.add_subplot(grid[0, 2])
    ax.axhline(1.0, color=kg.GRID, lw=0.9)
    ax.axhline(0.0, color=kg.EDGE, lw=0.9)
    ax.plot([915.0], [1.0], "o", color=kg.INK, mfc=kg.SURFACE, ms=9, label="single-mode decomposed model")
    if rg:
        for k, (fit, row) in enumerate(zip(rg["modes"], verdicts["ring"][rtag], strict=True)):
            colour = (ring_a, ring_b)[k]
            ax.plot(
                [fit["f"] / 1e6], [fit["inner_over_outer"]], "o", color=colour, ms=9, label=f"full wave, mode {'AB'[k]}"
            )
            ax.plot(
                [row["model"]["f"] / 1e6],
                [row["model"]["inner_over_outer"]],
                "x",
                color=kg.INK,
                ms=9,
                mew=1.6,
                label="two-rail decomposed model" if k == 0 else None,
            )
            ax.annotate(
                f"{fit['f'] / 1e6:.1f} MHz",
                (fit["f"] / 1e6, fit["inner_over_outer"]),
                (8, -14 if k == 0 else 8),
                textcoords="offset points",
                fontsize=8.0,
            )
    ax.set_xlabel("m = 1 resonance (MHz)")
    ax.set_ylabel("inner rail's volts over the outer's")
    ax.legend(loc="lower left", fontsize=7.6)
    ax.set_title("(c) R1's ring: the rails resonate apart")

    # (d), (e) DR-0
    rows_f: list[tuple[str, float, str]] = []
    rows_q: list[tuple[str, float, float, str]] = []
    for tag, v in verdicts["straight"].items():
        rows_f.append((f"straight, {tag} mesh", v["f_dev"], straight_c))
        rows_q.append((f"straight, {tag} mesh", v["Q_dev"], v["Q_cal_dev"], straight_c))
    for tag, rows in verdicts["ring"].items():
        for k, row in enumerate(rows):
            colour = (ring_a, ring_b)[k]
            rows_f.append((f"ring {'AB'[k]} vs single-mode model", row["f_dev_single"], colour))
            rows_f.append((f"ring {'AB'[k]} vs two-rail model", row["f_dev"], colour))
            rows_q.append((f"ring {'AB'[k]} vs two-rail model ({tag})", row["Q_dev"], row["Q_cal_dev"], colour))
    ax = fig.add_subplot(grid[1, 0])
    ax.axvspan(-100 * DR0_F, 100 * DR0_F, color=kg.PALE, lw=0)
    for k, (_, dev, colour) in enumerate(rows_f):
        ax.plot([100 * dev], [k], "o", color=colour, ms=8)
        ax.text(100 * dev + 0.12, k, f"{100 * dev:+.2f} %", ha="left", va="center", fontsize=7.6, color=kg.INK2)
    ax.set_yticks(range(len(rows_f)))
    ax.set_yticklabels([r[0] for r in rows_f], fontsize=8.0)
    ax.set_ylim(len(rows_f) - 0.35, -0.6)
    ax.set_xlabel("full wave's resonance over the decomposed model's, minus 1 (%)")
    ax.set_title(f"(d) DR-0, resonance: within {100 * DR0_F:.1f} % (shaded)")
    ax = fig.add_subplot(grid[1, 1])
    ax.axvspan(-100 * DR0_Q, 100 * DR0_Q, color=kg.PALE, lw=0)
    for k, (_, raw, cal, colour) in enumerate(rows_q):
        ax.plot(
            [100 * raw],
            [k],
            "o",
            color=colour,
            mfc=kg.SURFACE,
            ms=8,
            mew=1.6,
            label="the copper read off the FDTD's cells" if k == 0 else None,
        )
        ax.plot([100 * cal], [k], "o", color=colour, ms=8, label="calibrated by its static twin" if k == 0 else None)
    ax.set_yticks(range(len(rows_q)))
    ax.set_yticklabels([r[0] for r in rows_q], fontsize=8.0)
    ax.set_ylim(len(rows_q) - 0.5, -0.5)
    ax.set_xlabel("full wave's Q over the decomposed model's, minus 1 (%)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=7.6)
    ax.set_title(f"(e) DR-0, Q: within {100 * DR0_Q:.0f} % (shaded)")

    # (f) the copper's Q against the cell at the copper
    ax = fig.add_subplot(grid[1, 2])
    pts = sorted((v["grading"]["fine"] * 1e3, v["Qc_faces"]) for k, v in sec.items() if k.startswith("faces"))
    ax.plot(
        [p_[0] for p_ in pts],
        [p_[1] for p_ in pts],
        "o-",
        color=kg.INK2,
        ms=6,
        lw=1.2,
        label="static field, the study's grading (0.1 per mm)",
    )
    twins = [
        (sec[f"fdtd {t}"]["grading"]["fine"] * 1e3, sec[f"fdtd {t}"]["Qc_faces"])
        for t in ("fine", "coarse")
        if f"fdtd {t}" in sec
    ]
    if twins:
        ax.plot(
            [p_[0] for p_ in twins],
            [p_[1] for p_ in twins],
            "s",
            color=kg.INK,
            mfc=kg.SURFACE,
            ms=9,
            mew=1.4,
            label="static field on the FDTD's cells (its twin)",
        )
    full = [(v["grading"]["fine"] * 1e3, v["Qc"]) for v in results.get("straight", {}).values()]
    if full:
        ax.plot(
            [p_[0] for p_ in full],
            [p_[1] for p_ in full],
            "o",
            color=straight_c,
            ms=7,
            label="full wave, the straight line's H dump",
        )
    ax.axhline(b["Qc"], color=kg.INK, lw=1.0)
    ax.text(0.55, b["Qc"] + 0.8, f"Wheeler's rule, converged: {b['Qc']:.1f}", fontsize=7.8, va="bottom", color=kg.INK2)
    ax.set_xscale("log")
    kg._log_ticks(ax.xaxis, [0.5, 1, 2, 5, 10])
    ax.set_xlabel("finest cell at the copper (um)")
    ax.set_ylabel("Q_c by the face integral, bare copper")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.08), fontsize=7.6)
    ax.set_title("(f) the copper's loss read off a grid's faces")

    title = "R1's rail line, 915 MHz, in openEMS against its 2-D section (DR-0)"
    if st:
        v = next(iter(verdicts["straight"].values()))
        title += (
            f": the straight line resonates {100 * v['f_dev']:+.2f} % from the decomposed model and its Q is "
            f"{100 * v['Q_cal_dev']:+.1f} % off ({100 * v['Q_dev']:+.1f} % with the copper read off the cells)"
        )
    if rg:
        a_, b_ = rg["modes"]
        title += (
            f"; closed into R1's ring the rails resonate apart, at {a_['f'] / 1e6:.1f} and {b_['f'] / 1e6:.1f} MHz "
            f"with the inner rail at {a_['inner_over_outer']:+.2f} and {b_['inner_over_outer']:+.2f} of the outer, "
            f"which a two-rail model of the section predicts to "
            f"{100 * max(abs(r['f_dev']) for r in verdicts['ring'][rtag]):.2f} %"
        )
    kg._suptitle(fig, title, width=170)
    fig.savefig(out / "twin_openems_line.png", bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    say("  wrote twin_openems_line.png")


if __name__ == "__main__":
    sys.exit(main())
