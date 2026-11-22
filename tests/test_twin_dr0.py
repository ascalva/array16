"""DR-0, the instruments' review of docs/KNOTTED_GUIDE.md section 8.2, held on the chain the twin bench runs: the
cross-section solver against the gapless boards, ELE-1 to ELE-3 to 0.5 %, through the SectionGrid the twin sweeps; the
laminate's share against d ln C'/d ln eps_r on the designs' own sections, 2 %; a tracked electron's secular frequency
in a solved section against Mathieu's exact one, 0.1 % at 120 steps an RF period, and missed at 30, the case that must
fail; the ring's splitting against twice its impedance harmonic on drawn realisations, 1 %; the multipactor model at
ECSS's held-out point, a consistency check (DR-0's blind check is tests/test_guide_multipactor.py's); the receded solve
against stripline's closed form, 1 %, the control with ideal planes missing it, which holds the copper's loss to
Wheeler's closed form only (Pozar's form is the rule on a closed-form Z0); and R1's line and ring in full wave
(experiments/twin_openems_line) against the decomposed model, 0.5 % in f and 15 % in Q, with the twin's own R1 section
against that run's. The full wave holds n_eff, Z0, the share rule through Q_d and the two-rail ring; its copper is PEC,
read off the cells and calibrated by the static twin, and the face integral's limit, 3.5 to 5 % under the rule, is the
check of the copper's loss that does not use it (tests/test_twin_openems_line.py)."""

from __future__ import annotations

import json
import math
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from array16 import particles as pt
from array16.guide import multipactor as mp
from array16.guide import ring
from array16.guide import sections as gs
from array16.particles import FacingBoards, InPhase, RailPair, SectionField
from array16.twin import Design, TwinGuide

ROOT = Path(__file__).resolve().parents[1]
DESIGNS = ROOT / "designs"
FULL_WAVE = ROOT / "experiments" / "twin_openems_line" / "results.json"
H = 1e-3  # m: the gapless sections' h


def _coarse(name: str, where: Path) -> Path:
    """A design file on the coarse grid, written beside the test: the instruments' checks need no converged grid."""
    text = (DESIGNS / f"{name}.toml").read_text()
    assert 'grid = "default"' in text
    out = where / f"{name}.toml"
    out.write_text(text.replace('grid = "default"', 'grid = "coarse"'))
    return out


def _ring_field(field: SectionField, x: pt.Arr) -> pt.Arr:
    e, _ = field(x, 0.0)
    assert e is not None
    return e


def _around(radius: float, count: int = 48) -> tuple[pt.Arr, pt.Arr]:
    angle = np.linspace(0, 2 * math.pi, count, endpoint=False)
    return radius * np.cos(angle), radius * np.sin(angle)


# ---- ELE-1 to ELE-3 to 0.5 % ----------------------------------------------------------------------------------------
def test_dr_0_ele_2_the_solved_strips_and_their_swept_grid_are_the_gapless_boards() -> None:
    # ELE-2: strips 24h wide split 0.02h on 2 um copper in vacuum hold the gradient share pi/4 and the depth share 1/4
    # to 0.5 %, and swept straight (a static volt) their field round the null, 0.1h and 0.3h out, is the gapless
    # boards' -i (V/2h) tanh(pi z/2h) to 0.5 % of its largest there
    boards = FacingBoards(2 * H)
    solved = gs.solve(gs.Section(gs.FacingStrips(gap=2 * H, width=24 * H, split=0.02 * H), copper=2e-6, grid=gs.COARSE))
    assert solved.gradient_share == pytest.approx(boards.gradient_share, rel=5e-3)
    assert solved.depth_share == pytest.approx(boards.depth_share, rel=5e-3)
    grid = solved.section_grid(3 * H)
    assert grid.null == (0.0, H)
    field = SectionField.straight(grid, 0.05, InPhase(1.0, 0.0))
    for radius in (0.1 * H, 0.3 * H):
        a, b = _around(radius)
        e = _ring_field(field, np.column_stack([a, b, np.linspace(0.0, 0.05, len(a))]))
        closed = boards.field(a, b, 1.0)
        assert np.abs(e[:, :2] - closed).max() <= 5e-3 * np.abs(closed).max()
        assert not e[:, 2].any()


def test_dr_0_ele_1_and_ele_3_the_solved_rails_and_their_swept_grid_are_the_gapless_pair() -> None:
    # ELE-1 and ELE-3 at the deepest rail pair, rail/gap 3.68: rails 0.01h from their ground on 0.01h copper in
    # vacuum, walls 100h out, put the null, its gradient, the saddle, its field and the depth share where RailPair does
    # to 0.5 %, and their swept grid's field 0.1h and 0.3h round the null at the pair's own there, to 0.5 % of its
    # largest
    pair = RailPair(0.4e-3, 3.68 * 0.4e-3)
    h, saddle = pair.null_height, pair.escape(1.0)
    solved = gs.solve(gs.Section(gs.Rails(pair.gap, pair.rail, 0.01 * h), copper=0.01 * h, grid=gs.COARSE))
    assert solved.null[1] == pytest.approx(h, rel=5e-3)
    assert solved.gradient == pytest.approx(pair.gradient(1.0), rel=5e-3)
    assert solved.saddle[1] == pytest.approx(saddle.height, rel=5e-3)
    assert solved.escape == pytest.approx(saddle.field, rel=5e-3)
    assert solved.depth_share == pytest.approx(pair.depth_share, rel=5e-3)
    grid = solved.section_grid(3 * h)
    assert grid.null[0] == 0.0
    assert grid.null[1] == pytest.approx(solved.null[1], abs=1e-3 * h)
    field = SectionField.straight(grid, 0.05, InPhase(1.0, 0.0))
    for radius in (0.1 * h, 0.3 * h):
        a, b = _around(radius)
        pu, pv = a, h + b  # the closed form's coordinates: v up from the board, the null at h
        x = np.column_stack([pu - grid.null[0], pv - grid.null[1], np.zeros(len(a))])
        e = _ring_field(field, x)
        closed = pair.field(pu, pv, 1.0)
        assert np.abs(e[:, :2] - closed).max() <= 5e-3 * np.abs(closed).max()


# ---- the laminate's share -------------------------------------------------------------------------------------------
@pytest.mark.parametrize("name", ["r1_rail_ring", "t915s"])
def test_dr_0_the_laminate_share_is_d_ln_c_over_d_ln_er_on_the_designs(name: str) -> None:
    # DR-0: on R1's rails over 1.6 mm FR-4 (the board holding three-quarters of the field) and section B's 25 um
    # polyimide (a fourteenth), the share the solve reports equals the slope of ln C' against ln eps_r, by central
    # differences 1 % either way, within 2 %
    section = replace(Design.load(DESIGNS / f"{name}.toml").section, grid=gs.COARSE)
    laminate = section.laminate
    assert laminate is not None
    share = gs.solve(section).share
    up, down = (gs.solve(replace(section, laminate=replace(laminate, er=laminate.er * k))) for k in (1.01, 0.99))
    slope = math.log(up.capacitance / down.capacitance) / math.log(1.01 / 0.99)
    assert slope == pytest.approx(share, rel=0.02)


# ---- the secular frequency to 0.1 % ---------------------------------------------------------------------------------
def _mathieu_beta(q: float, terms: int = 40) -> float:
    """Mathieu's characteristic exponent at a = 0, the beta of the Floquet solutions e^(i beta tau) P(tau) of
    y'' + 2q cos(2 tau) y = 0: the root in (0, 1) of beta^2 = q (c_2 + c_-2) / c_0, the two ratios continued fractions
    of (2j + beta)^2 c_2j = q (c_2j-2 + c_2j+2) (McLachlan; Abramowitz and Stegun 20.3), by bisection."""

    def balance(beta: float) -> float:
        up = down = 0.0
        for j in range(terms, 0, -1):
            up = q / ((2 * j + beta) ** 2 - q * up)
            down = q / ((beta - 2 * j) ** 2 - q * down)
        return beta * beta - q * (up + down)

    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if balance(mid) < 0 else (lo, mid)
    return (lo + hi) / 2


def _tracked_beta(field: SectionField, f: float, steps_per_period: int, periods: int, offset: float) -> pt.Arr:
    """Mathieu's beta as four electrons show it, started `offset` off the null along the two diagonals (the section's
    principal axes) and along u and v, 1 eV along the path: read once an RF period, x(k+1) + x(k-1) = 2 cos(mu) x(k)
    for a linear oscillation, cos mu a least-squares ratio, beta = mu / pi."""
    root = 1 / math.sqrt(2)
    starts = offset * np.array([[root, root], [root, -root], [1.0, 0.0], [0.0, 1.0]])
    v = pt.v_of(1.0)
    u0 = np.tile([0.0, 0.0, v / math.sqrt(1 - (v / pt.C0) ** 2)], (4, 1))
    seen = [starts.copy()]

    def strobe(step: int, t: float, x: pt.Arr, u: pt.Arr, alive: pt.BoolArr) -> None:
        if (step + 1) % steps_per_period == 0:
            seen.append(x[:, :2].copy())

    flight = pt.integrate(
        np.column_stack([starts, np.zeros(4)]),
        u0,
        np.full(4, pt.QM_E),
        1 / (f * steps_per_period),
        steps_per_period * periods,
        field,
        [strobe],
    )
    assert bool(flight.alive.all())
    assert not field.off_grid.any()
    path = np.array(seen)
    inner, outer = path[1:-1], path[2:] + path[:-2]
    cos_mu = (inner * outer).sum(axis=(0, 2)) / (2 * (inner * inner).sum(axis=(0, 2)))
    out: pt.Arr = np.arccos(cos_mu) / math.pi
    return out


def test_dr_0_a_tracked_electron_in_the_solved_section_keeps_mathieus_exact_secular_frequency() -> None:
    # DR-0: in T915-S's solved section B (the copper, the split and the polyimide resolved), swept straight and driven
    # at 915 MHz to q = 0.3 at the null, electrons 0.004h off it keep Mathieu's exact beta (0.2161, the first order
    # 1.8 % low) to 0.1 % at 120 steps an RF period; at 30 steps the integrator's own error, 0.19 %, misses the
    # criterion on the principal axes
    section = replace(Design.load(DESIGNS / "t915s.toml").section, grid=gs.COARSE)
    solved = gs.solve(section)
    h, f = solved.distance, section.f
    grid = solved.section_grid(2 * h)
    probe = 1e-3 * h
    du_up, _, _ = grid.gradient_at(np.array([0.0]), np.array([grid.null[1] + probe]))
    du_down, _, _ = grid.gradient_at(np.array([0.0]), np.array([grid.null[1] - probe]))
    assert abs(float(du_up[0] - du_down[0]) / (2 * probe)) == pytest.approx(solved.gradient, rel=1e-3)
    volts = 0.3 / pt.mathieu_q(solved.gradient, f)
    field = SectionField.straight(grid, 0.05, InPhase(volts, 2 * math.pi * f))
    exact = _mathieu_beta(0.3)
    assert exact / (0.3 / math.sqrt(2)) - 1 > 0.018
    fine = _tracked_beta(field, f, 120, 60, 4e-3 * h)
    assert float(np.abs(fine / exact - 1).max()) < 1e-3
    coarse = _tracked_beta(field, f, 30, 60, 4e-3 * h)
    assert float(np.abs(coarse[:2] / exact - 1).min()) > 1e-3


# ---- the splitting to 1 % -------------------------------------------------------------------------------------------
def test_dr_0_a_drawn_rings_splitting_is_twice_its_impedance_harmonic(tmp_path: Path) -> None:
    # DR-0: on three realisations of T915-S (spacers at 2 % rms, the closing joints, the split's etch, the laminate),
    # the twin's tuned +-1 pair splits by 2|c_2m| of each drawn impedance profile within 1 %; the design as drawn does
    # not split
    path = _coarse("t915s", tmp_path)
    for k in (1, 2, 3):
        guide = TwinGuide(path, k, cache=tmp_path / "cache")
        harmonic = ring.ring_harmonic(guide.drawn.z, 1, guide.build.z0)
        assert harmonic > 1e-3
        assert guide.tune().splitting == pytest.approx(2 * harmonic, rel=0.01)
    assert abs(TwinGuide(path, 0, cache=tmp_path / "cache").tune().splitting) < 1e-12


# ---- multipactor's held-out point, a consistency check ---------------------------------------------------------------
def test_the_multipactor_model_is_consistent_with_ecss_at_its_held_out_point() -> None:
    # calibrated by one gain on copper's lowest point (27.4 V at 0.70 GHz mm, met to 0.2 %), the statistical model reads
    # ECSS's 91.9 V at 1.68 GHz mm, the published point nearest T915-S's 1.674 GHz mm, within 10 %: computed while the
    # model was built, and ECSS's chart is itself a model of the same yield, so two models agree here
    assert mp.threshold(mp.CALIBRATION[0]) == pytest.approx(mp.CALIBRATION[1], rel=2e-3)
    assert mp.threshold(mp.HELD_OUT[0]) == pytest.approx(mp.HELD_OUT[1], rel=0.10)


# ---- the receded solve against Wheeler's closed form to 1 % -------------------------------------------------------
@pytest.mark.parametrize(
    ("width", "thick", "spacing"), [(1.0e-3, 35e-6, 1.6e-3), (2.0e-3, 18e-6, 1.0e-3), (1.0e-3, 18e-6, 1.0e-3)]
)
def test_dr_0_wheelers_rule_on_stripline(width: float, thick: float, spacing: float) -> None:
    # DR-0, the copper's loss: Wheeler's incremental-inductance rule on the converged grid, every copper face receded,
    # gives stripline's attenuation (Pozar eq. 3.198, the rule applied to a closed-form Z0(t)) within 1 %: the numerical
    # derivative is right, and the loss is held to Wheeler's closed form only; with the planes ideal it misses by their
    # share
    assert gs.stripline_wheeler(width, thick, spacing).ratio == pytest.approx(1.0, abs=0.01)
    assert gs.stripline_wheeler(width, thick, spacing, planes_lose=False).ratio < 0.9


# ---- full wave against the decomposed model -------------------------------------------------------------------------
def _full_wave() -> dict[str, object]:
    """The full-wave record experiments/twin_openems_line.py wrote: in the tree, so a missing one fails."""
    out: dict[str, object] = json.loads(FULL_WAVE.read_text())
    return out


def test_dr_0_full_wave_meets_the_decomposed_model_on_r1s_line_and_its_two_rail_ring() -> None:
    # DR-0 and section 8.1's Hypothesis, as experiments/twin_openems_line recorded them: R1's straight line between
    # magnetic walls resonates within 0.5 % of the 2-D section's 915 MHz on both meshes, its Q within 15 % read off the
    # cells and calibrated by their static twin, Q_d setting most of it (the copper's calibration is Wheeler's rule
    # again); closed into the ring, each of its two modes meets the two-rail model within 0.5 % in f and 15 % in Q,
    # while the single-mode ring misses both by more than 0.5 %: the rails resonate apart, 899 and 962 MHz, and DR-0
    # does not pass the single-mode ring E2-21's drive uses
    record = _full_wave()
    straight, closed = record["dr0"]["straight"], record["dr0"]["ring"]  # type: ignore[index]
    for mesh in ("fine", "coarse"):
        line = straight[mesh]
        assert abs(line["f_dev"]) < 5e-3
        assert abs(line["Q_dev"]) < 0.15
        assert abs(line["Q_cal_dev"]) < 0.15
        run = record["straight"][mesh]  # type: ignore[index]
        assert line["f_dev"] == pytest.approx(run["f"] / 915e6 - 1, abs=1e-9)
    for mode in closed["fine"]:
        assert abs(mode["f_dev"]) < 5e-3
        assert abs(mode["Q_dev"]) < 0.15
        assert abs(mode["Q_cal_dev"]) < 0.15
        assert abs(mode["f_dev_single"]) > 5e-3


def test_dr_0_the_twins_r1_section_is_the_full_wave_runs_own() -> None:
    # the full-wave script's 2-D solve of R1's board as written (no plane, the strip held at ground, 0.15 mm gaps, a
    # shield 25 mm out) is the section designs/r1_rail_ring.toml puts on the twin: array16.guide.sections gives its
    # n_eff, Z0, share, Q_d, Q_c and Q within 0.5 % on the coarse grid (0.16 % at most)
    written = _full_wave()["section"]["written"]  # type: ignore[index]
    solved = gs.solve(replace(Design.load(DESIGNS / "r1_rail_ring.toml").section, grid=gs.COARSE))
    for ours, theirs in (
        (solved.index, written["index"]),
        (solved.z0, written["z0"]),
        (solved.share, written["share"]),
        (solved.qd, written["Qd"]),
        (solved.qc, written["Qc"]),
        (solved.quality, written["Q"]),
    ):
        assert ours == pytest.approx(theirs, rel=5e-3)
