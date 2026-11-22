"""array16.guide.sections, the cross-section solver, held to the gapless boards (ELE-1 to ELE-3), the knotted study's
line solve (experiments/knotted_guide/run.log), stripline's closed form, DR-0's share rule, the split as a lever and the
slot guide's Q (docs/KNOTTED_GUIDE.md sections 3.3, 3.6, 5.1, 5.4, 8.1 and 8.2)."""

from __future__ import annotations

import math
from dataclasses import replace

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

from array16.guide import sections as gs
from array16.particles import FacingBoards, RailPair, mathieu_q, ponderomotive_ev

H = 1e-3  # m: the gapless boards' h
F915 = 915e6
FH = 0.839e6  # Hz m: f h of a 1 eV well at q = 0.5 on facing boards (50.3 V)
LINE_GRID = gs.Grid(fine=1 / 8, coarse=0.1, null=False)  # the study's coarse line grid (run.log's "coarse grid")


@pytest.fixture(scope="module")
def gapless_boards() -> gs.Solved:
    """ELE-2's boards as strips: 24h wide, split 0.02h, 2 um copper, in vacuum, on the coarse grid."""
    geometry = gs.FacingStrips(gap=2 * H, width=24 * H, split=0.02 * H)
    return gs.solve(gs.Section(geometry, copper=2e-6, grid=gs.COARSE))


@pytest.fixture(scope="module")
def slot_fr4() -> gs.Solved:
    """The slot guide on FR-4 0.8 mm: 0.5 mm strips of 35 um copper set back 0.2 mm, planes 1 mm off."""
    return gs.solve(gs.slot_guide(gs.FR4, 0.8e-3))


def test_ele_2_facing_strips_reach_the_gapless_boards(gapless_boards: gs.Solved) -> None:
    # ELE-2: the gapless boards' null has gradient pi V/(4 h^2), a well a quarter of four rods' deep at the same volts
    # and distance, and potential -(V/pi) atan(tan(pi y/2h) tanh(pi x/2h)); strips 24h wide split 0.02h hold each to
    # 0.5 %, and so the Mathieu q and the depth at 50 V
    boards, r = FacingBoards(2 * H), gapless_boards
    assert r.null == (0.0, H)
    assert r.distance == H
    assert r.gradient_share == pytest.approx(boards.gradient_share, rel=5e-3)
    assert r.depth_share == pytest.approx(boards.depth_share, rel=5e-3)
    assert r.saddle[1] == pytest.approx(H, abs=1e-3 * H)  # the barrier is on the midplane, far out along the gap
    assert r.mathieu_q(50.0) == pytest.approx(mathieu_q(boards.gradient(50.0), F915), rel=5e-3)
    assert r.depth_ev(50.0) == pytest.approx(ponderomotive_ev(boards.barrier_field(50.0), F915), rel=5e-3)
    u, v = np.meshgrid(r.u, r.v, indexing="ij")
    gap = (np.abs(u) < 6 * H) & (v > 0.05 * H) & (v < 1.95 * H)
    assert np.abs(r.phi[gap] - boards.potential(u[gap], v[gap] - H, 1.0)).max() < 5e-3


@pytest.mark.parametrize("ratio", [3.68, 1.0])
def test_ele_1_and_ele_3_rails_reach_the_gapless_pair(ratio: float) -> None:
    # ELE-1: the null at h = sqrt(c1 c2) with gradient 4 V b/(pi h (a + b)^2); ELE-3: the saddle h sqrt(1 + w) up the
    # axis with field 2 V b/(pi (a + b)(a + b + 2h)), the depth share (4/pi^2)(w - 2)/((w + 2) w^2), 0.0091 at rail/gap
    # 3.68; rails 0.01h from their ground on 0.01h copper in vacuum, walls 100h out, hold each to 0.5 %, and the
    # potential to 0.5 % of V above the board
    pair = RailPair(0.4e-3, ratio * 0.4e-3)
    h, saddle = pair.null_height, pair.escape(1.0)
    r = gs.solve(gs.Section(gs.Rails(pair.gap, pair.rail, 0.01 * h), copper=0.01 * h, grid=gs.COARSE))
    assert r.null[0] == 0.0
    assert r.null[1] == pytest.approx(h, rel=5e-3)
    assert r.distance == r.null[1]
    assert r.gradient == pytest.approx(pair.gradient(1.0), rel=5e-3)
    assert r.saddle[0] == pytest.approx(0.0, abs=1e-3 * h)
    assert r.saddle[1] == pytest.approx(saddle.height, rel=5e-3)
    assert r.escape == pytest.approx(saddle.field, rel=5e-3)
    assert r.depth_share == pytest.approx(pair.depth_share, rel=5e-3)
    u, v = np.meshgrid(r.u, r.v, indexing="ij")
    above = (np.abs(u) < 4 * h) & (v > 0.1 * h) & (v < 4 * h)
    assert np.abs(r.phi[above] - pair.potential(u[above], v[above], 1.0)).max() < 5e-3


@pytest.mark.parametrize(
    ("h", "f", "expected"),
    [
        (0.915e-3, F915, (82.62, 1.0312, 0.0739, 41.63, 226.0, 211.9)),
        (0.3425e-3, 2450e6, (91.09, 1.0659, 0.1533, 39.03, 145.3, 133.4)),
    ],
)
def test_section_b_on_the_study_line_grid(h: float, f: float, expected: tuple[float, ...]) -> None:
    # run.log, section B with the copper's 18 um and the polyimide's 25 um resolved on the coarse line grid: C' (pF/m),
    # n_eff, share, Z0 (ohm), Q_c by Wheeler's rule and Q, to the digits it prints; 915 MHz at h = 0.915 mm and
    # 2450 MHz at h = 0.3425 mm
    r = gs.solve(gs.section_b(h, f, grid=LINE_GRID))
    c, index, share, z0, qc, quality = expected
    assert r.capacitance * 1e12 == pytest.approx(c, abs=5e-3)
    assert r.index == pytest.approx(index, abs=5e-5)
    assert r.share == pytest.approx(share, abs=5e-5)
    assert r.z0 == pytest.approx(z0, abs=5e-3)
    assert r.qc == pytest.approx(qc, abs=0.05)
    assert r.quality == pytest.approx(quality, abs=0.05)
    assert r.qd == pytest.approx(1 / (r.share * gs.POLYIMIDE.tand))


def test_section_b_converged() -> None:
    # run.log's fine grid at 915 MHz: n_eff 1.0311, the polyimide 7.37 % of the field energy, Z0 41.65 ohm, Q_c 228.2
    # and Q 213.8 on bare copper, so T915-S's m = 1 ring takes 317.8 mm (section 5.1)
    r = gs.solve(gs.section_b(grid=gs.Grid(null=False)))
    assert r.capacitance * 1e12 == pytest.approx(82.57, abs=5e-3)
    assert r.index == pytest.approx(1.0311, abs=5e-5)
    assert r.share == pytest.approx(0.0737, abs=5e-5)
    assert r.z0 == pytest.approx(41.65, abs=5e-3)
    assert r.qc == pytest.approx(228.2, abs=0.05)
    assert r.quality == pytest.approx(213.8, abs=0.05)
    assert r.wavelength * 1e3 == pytest.approx(317.8, abs=0.05)


@pytest.mark.parametrize(
    ("width", "thick", "spacing"), [(1.0e-3, 35e-6, 1.6e-3), (2.0e-3, 18e-6, 1.0e-3), (1.0e-3, 18e-6, 1.0e-3)]
)
def test_wheeler_rule_on_stripline(width: float, thick: float, spacing: float) -> None:
    # section 3.3's calibration: Wheeler's rule on the converged grid, the strip's and the planes' copper receded,
    # against stripline's closed form (Pozar eq. 3.198), within 1 %; the control, with the planes ideal, misses it by
    # more than a tenth, the planes' part of the loss
    cal = gs.stripline_wheeler(width, thick, spacing)
    assert cal.ratio == pytest.approx(1.0, abs=0.01)
    assert cal.alpha == pytest.approx(cal.closed_form * cal.ratio)
    assert gs.stripline_wheeler(width, thick, spacing, planes_lose=False).ratio < 0.9


@pytest.mark.parametrize(
    "section",
    [
        gs.section_b(grid=gs.COARSE),
        gs.slot_guide(gs.FR4, 0.8e-3, grid=gs.COARSE),
        gs.Section(gs.Rails(0.40e-3, 1.46e-3, 0.15e-3, walls=10.0), gs.FR4, 1.6e-3, gs.RIGID_COPPER, grid=gs.COARSE),
    ],
    ids=["section B", "slot guide", "rails"],
)
def test_the_laminate_share_is_d_ln_c_over_d_ln_er(section: gs.Section) -> None:
    # DR-0: the laminate's share of the field energy equals d ln C'/d ln eps_r within 2 % (central differences at
    # +-1 % in eps_r), on a thin film, a laminate holding most of the field, and a board under rails
    laminate = section.laminate
    assert laminate is not None
    share = gs.solve(section).share
    up, down = (gs.solve(replace(section, laminate=replace(laminate, er=laminate.er * k))) for k in (1.01, 0.99))
    slope = math.log(up.capacitance / down.capacitance) / math.log(1.01 / 0.99)
    assert slope == pytest.approx(share, rel=0.02)


def test_the_split_as_a_lever() -> None:
    # section 5.4: at the same 1 eV well at q = 0.5 and 915 MHz, the well solved at each split, a split of h against
    # 0.3h raises Q from 210 to 290; the null's gradient falls from 0.77 to 0.61 V/h^2, so h closes from 0.90 to 0.72 mm
    # to keep q, the volts fall from 49.5 to 39.6 V and the power from 0.88 to 0.33 W
    narrow, wide = (gs.fit_well(gs.section_b(h=FH / F915, split=g, grid=gs.COARSE)) for g in (0.3, 1.0))
    assert (narrow.quality, wide.quality) == pytest.approx((210, 290), rel=0.01)
    assert wide.quality / narrow.quality == pytest.approx(1.38, abs=0.02)
    assert (narrow.gradient_share, wide.gradient_share) == pytest.approx((0.768, 0.614), abs=2e-3)
    assert (narrow.distance * 1e3, wide.distance * 1e3) == pytest.approx((0.899, 0.718), abs=2e-3)
    assert (narrow.volts_for(1.0, 0.5), wide.volts_for(1.0, 0.5)) == pytest.approx((49.5, 39.6), abs=0.1)
    assert (narrow.power_for(1.0, 0.5), wide.power_for(1.0, 0.5)) == pytest.approx((0.878, 0.326), rel=0.01)
    for r in (narrow, wide):
        volts = r.volts_for(1.0, 0.5)
        assert r.mathieu_q(volts) == pytest.approx(0.5, rel=0.01)
        assert r.depth_ev(volts) == pytest.approx(1.0, rel=0.01)


def test_the_slot_guide_q_on_fr4_and_ro4003c(slot_fr4: gs.Solved) -> None:
    # section 3.3: the quadrupole puts the laminate between opposite volts, so it holds about four-fifths of the field
    # energy; with 0.5 mm strips on 35 um copper Q_c is 158 by Wheeler's rule and Q 43 on FR-4 (its loss tangent sets
    # it) and 117 on RO4003C; the null on the mid-plane at 0.33 to 0.34 G0, the escape 0.10 to 0.11 D0
    ro = gs.solve(gs.slot_guide(gs.RO4003C, 0.813e-3))
    assert slot_fr4.quality == pytest.approx(43, abs=1)
    assert ro.quality == pytest.approx(117, abs=1)
    for r in (slot_fr4, ro):
        assert r.qc == pytest.approx(158, rel=0.02)
        assert 0.78 < r.share < 0.86
        assert r.null == pytest.approx((0.0, r.section.thickness / 2))
        assert 0.33 <= r.gradient_share / (math.pi / 4) <= 0.345
        assert 0.10 <= 4 * r.depth_share <= 0.115
    assert slot_fr4.qd == pytest.approx(1 / (slot_fr4.share * gs.FR4.tand))
    # the fab's ENIG raises the surface resistance 1.3 to 5 times; Q_c falls in proportion
    enig = gs.solve(replace(slot_fr4.section, resistance=1.3, grid=gs.COARSE))
    plain = gs.solve(replace(slot_fr4.section, grid=gs.COARSE))
    assert enig.qc == pytest.approx(plain.qc / 1.3, rel=1e-12)


def test_the_treated_faces_and_lossy_planes_take_their_share_of_the_copper_loss(slot_fr4: gs.Solved) -> None:
    # the copper's faces on the laminate, its treated side, take a roughness factor k on their part p of the loss:
    # Q_c falls as 1 / (1 + (k - 1) p), the same p at every k (Wheeler's rule recedes those faces alone for it), about
    # four-tenths of the strips' loss on the slot guide's board and a quarter on R1's rails over FR-4; k = 1 is the
    # solve without it, to the bit. The slot guide's planes made lossy copper take their own small part, so the ideal
    # planes' Q_c bounds it from above
    rails = gs.Section(gs.Rails(0.40e-3, 1.46e-3, 0.15e-3, walls=10.0), gs.FR4, 1.6e-3, gs.RIGID_COPPER, grid=gs.COARSE)
    for section, lo, hi in ((replace(slot_fr4.section, grid=gs.COARSE), 0.3, 0.5), (rails, 0.2, 0.35)):
        plain = gs.solve(section)
        assert gs.solve(replace(section, rough=1.0)).qc == plain.qc
        shares = [(plain.qc / gs.solve(replace(section, rough=k)).qc - 1) / (k - 1) for k in (1.37, 1.55)]
        assert shares[0] == pytest.approx(shares[1], rel=1e-9)
        assert lo < shares[0] < hi
    bare = gs.Section(gs.FacingStrips(gap=2 * H, width=8 * H, split=0.3 * H), copper=18e-6, grid=gs.COARSE)
    assert gs.solve(replace(bare, rough=1.5)).qc == gs.solve(bare).qc  # no laminate, no treated face
    section = replace(slot_fr4.section, grid=gs.COARSE)
    geometry = section.geometry
    assert isinstance(geometry, gs.SlotGuide)
    lossy = gs.solve(replace(section, geometry=replace(geometry, planes_lose=True)))
    plain = gs.solve(section)
    assert 0.9 * plain.qc < lossy.qc < plain.qc
    assert lossy.index == pytest.approx(plain.index, rel=1e-3)
    assert replace(geometry, planes_lose=True).scaled(2.0).planes_lose
    with pytest.raises(ValueError, match="roughness"):
        gs.Section(geometry, gs.FR4, 0.8e-3, rough=0.0)


def test_the_slot_guide_on_thin_copper_is_the_study_slot_line(slot_fr4: gs.Solved) -> None:
    # the study's slot line, zero-thickness copper on a uniform 25 um grid, gives n_eff 1.708 and share 0.845 on FR-4;
    # on 1 um copper the section solve reaches both within 0.2 %, so the 2 % the real 35 um copper takes off the index
    # (1.677) is the copper's
    thin = gs.solve(replace(slot_fr4.section, copper=1e-6, grid=gs.COARSE))
    assert thin.index == pytest.approx(1.708, rel=2e-3)
    assert thin.share == pytest.approx(0.845, rel=2e-3)
    assert slot_fr4.index == pytest.approx(1.677, abs=2e-3)
    assert slot_fr4.index < thin.index - 0.02


def test_r1_rails_on_their_board() -> None:
    # section 3.6: R1's rails, 1.46 mm either side of a 0.40 mm grounded strip, at 50 V and 915 MHz put the null 0.575
    # mm up at q = 0.5 in a 0.093 eV well (the gapless pair); on 1.6 mm FR-4 with 35 um copper and 0.1 mm isolation the
    # solve keeps all three to 1 %, while the laminate holds three-quarters of the line's field energy
    pair = RailPair(0.40e-3, 1.46e-3)
    section = gs.Section(gs.Rails(0.40e-3, 1.46e-3, 0.10e-3), gs.FR4, 1.6e-3, gs.RIGID_COPPER, F915, grid=gs.COARSE)
    r = gs.solve(section)
    assert r.null[1] == pytest.approx(pair.null_height, rel=0.01)
    assert r.mathieu_q(50.0) == pytest.approx(0.5, rel=0.01)
    assert r.depth_ev(50.0) == pytest.approx(0.093, rel=0.01)
    assert 0.7 < r.share < 0.8
    assert 1.5 < r.index < math.sqrt(gs.FR4.er)


def test_the_potential_arrays_carry_each_section_s_symmetry(gapless_boards: gs.Solved, slot_fr4: gs.Solved) -> None:
    # the swept field's arrays: u across the path from the centre line, v from the board toward the null; facing strips
    # and the slot guide odd in u and about the null's plane, rails even in u; read-only
    for r in (gapless_boards, slot_fr4):
        assert r.phi.shape == (len(r.u), len(r.v))
        assert np.array_equal(r.u, -r.u[::-1])
        assert np.allclose(r.v[::-1] - r.null[1], -(r.v - r.null[1]), atol=1e-15)
        assert np.array_equal(r.phi[::-1], -r.phi)
        assert np.array_equal(r.phi[:, ::-1], -r.phi)
    rails = gs.solve(gs.Section(gs.Rails(0.40e-3, 1.46e-3, 0.10e-3, walls=10.0), copper=35e-6, grid=gs.COARSE))
    assert np.array_equal(rails.phi[::-1], rails.phi)
    assert rails.phi[np.argmin(np.abs(rails.u)), np.argmin(np.abs(rails.v - rails.null[1]))] > 0
    with pytest.raises(ValueError, match="read-only"):
        gapless_boards.phi[0, 0] = 1.0


def test_sections_say_what_they_take() -> None:
    with pytest.raises(ValueError, match="split narrower"):
        gs.FacingStrips(gap=1e-3, width=1e-3, split=2e-3)
    with pytest.raises(ValueError, match="space narrower"):
        gs.Rails(0.40e-3, 1.46e-3, 0.5e-3)
    with pytest.raises(ValueError, match="takes a thickness"):
        gs.Section(gs.FacingStrips(gap=1e-3, width=4e-3, split=0.15e-3), gs.POLYIMIDE, 0.0)
    with pytest.raises(ValueError, match="board takes a thickness"):
        gs.solve(gs.Section(gs.SlotGuide(1e-3, 0.5e-3)))
    with pytest.raises(ValueError, match="er >= 1"):
        gs.Laminate("air gap", 0.5, 0.0)
    with pytest.raises(ValueError, match="recession"):
        gs.Grid(recede=0.6)


@settings(max_examples=40, deadline=None)
@given(st.data())
def test_block_elimination_solves_gauss_s_law(data: st.DataObject) -> None:
    # laplace2d's block elimination on random graded grids, cells and held nodes: Gauss's law holds on the free nodes,
    # the potential stays between the held values (the maximum principle), conjugate gradients agree, and so does the
    # elimination that keeps only checkpoints of its inverses
    nx, ny = data.draw(st.integers(3, 12)), data.draw(st.integers(3, 12))
    steps = st.floats(0.05, 2.0)
    x = np.concatenate([[0.0], np.cumsum(data.draw(arrays(np.float64, nx - 1, elements=steps)))])
    y = np.concatenate([[0.0], np.cumsum(data.draw(arrays(np.float64, ny - 1, elements=steps)))])
    eps = data.draw(arrays(np.float64, (nx - 1, ny - 1), elements=st.floats(1.0, 10.0)))
    fixed = data.draw(arrays(np.bool_, (nx, ny)))
    fixed[data.draw(st.integers(0, nx - 1)), data.draw(st.integers(0, ny - 1))] = True
    value = data.draw(arrays(np.float64, (nx, ny), elements=st.floats(-1.0, 1.0)))
    phi = gs.laplace2d(x, y, eps, fixed, value)
    wx, wy = gs.edge_weights(x, y, eps)
    flux = np.zeros_like(phi)
    fx, fy = wx * np.diff(phi, axis=0), wy * np.diff(phi, axis=1)
    flux[:-1] += fx
    flux[1:] -= fx
    flux[:, :-1] += fy
    flux[:, 1:] -= fy
    scale = float(wx.max() + wy.max())
    assert np.abs(flux[~fixed]).max(initial=0.0) < 1e-9 * scale
    assert phi.min() >= value[fixed].min() - 1e-9
    assert phi.max() <= value[fixed].max() + 1e-9
    assert np.allclose(gs.laplace2d(x, y, eps, fixed, value, tol=1e-14, direct=False), phi, atol=1e-8)
    assert np.allclose(gs.laplace2d(x, y, eps, fixed, value, limit=1), phi, atol=1e-12)


@settings(max_examples=60, deadline=None)
@given(
    st.floats(-5.0, 5.0),
    st.floats(0.5, 20.0),
    st.lists(st.floats(0.0, 1.0), max_size=5),
    st.floats(1e-3, 0.05),
    st.floats(0.05, 0.5),
    st.floats(0.05, 0.5),
)
def test_nodes_keep_their_features_and_spacing(
    lo: float, span: float, where: list[float], fine: float, coarse: float, growth: float
) -> None:
    # the graded grid: from lo to hi, increasing, every feature a node (one within fine/4 of an end is that end), no
    # spacing above coarse by more than a fine cell, and a stretched far field takes no more nodes than the capped one
    hi = lo + span
    features = [lo + w * span for w in where]
    expected = [lo if f - lo < fine / 4 else hi if hi - f < fine / 4 else f for f in features]
    out = gs.nodes(lo, hi, features, fine, coarse, growth)
    assert (out[0], out[-1]) == (lo, hi)
    assert np.all(np.diff(out) > 0)
    assert np.isin(expected, out).all()
    assert np.diff(out).max() <= coarse + fine
    stretched = gs.nodes(lo, hi, features, fine, coarse, growth, stretch=0.25, reach=1.0)
    assert (stretched[0], stretched[-1]) == (lo, hi)
    assert len(stretched) <= len(out)
    assert np.isin(expected, stretched).all()
