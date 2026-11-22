"""The knotted-guide study's fast parts, held to the numbers docs/KNOTTED_GUIDE.md states."""

from __future__ import annotations

import math
from pathlib import Path

import knotted_guide_study as study
import numpy as np
import pytest

from array16 import particles as pt


def test_the_lead_trefoil_and_its_framings() -> None:
    # sections 3.1 and 3.5: the r/R = 0.45 trefoil at 318 mm bends to 15.3 mm, turns through 5.37 pi, and its
    # twist-free frame comes back turned -153 degrees; k = 0 twists +0.48 deg/mm and keeps an in-plane radius of
    # 18.6 mm, k = -1 15.6 mm, k = -3 29.0 mm; k = 0's plates are 329.5 and 306.1 mm, k = -3's equal
    knot = study.torus(0.45, 0.318, 3000)
    assert knot.bend_radius * 1e3 == pytest.approx(15.3, abs=0.05)
    assert knot.total_curvature / math.pi == pytest.approx(5.37, abs=0.01)
    assert math.degrees(knot.holonomy) == pytest.approx(-153.0, abs=0.3)
    k0, k1, k3 = (study.framing(knot, k) for k in (0, -1, -3))
    assert math.degrees(k0.twist) * 1e-3 == pytest.approx(0.48, abs=0.005)
    assert (k0.inplane_radius * 1e3, k1.inplane_radius * 1e3, k3.inplane_radius * 1e3) == pytest.approx(
        (18.6, 15.6, 29.0), abs=0.1
    )
    plates = [float(study.plate(knot, k0.normal, offset)[0].sum()) * 1e3 for offset in (study.H915, -study.H915)]
    assert plates == pytest.approx([329.5, 306.1], abs=0.4)
    equal = [float(study.plate(knot, k3.normal, offset)[0].sum()) * 1e3 for offset in (study.H915, -study.H915)]
    assert equal[0] == pytest.approx(equal[1], abs=0.05)


def test_layouts_on_two_levels() -> None:
    # section 3.1's table: crossings, ramps and symmetry of the torus knots' layouts, and the knots' determinants
    cases = {
        (2, 3): (3, 3, 6, 3),
        (3, 2): (4, 3, 4, 2),
        (2, 5): (5, 5, 10, 5),
        (5, 2): (8, 5, 4, 2),
        (3, 4): (8, 3, 8, 4),
        (4, 3): (9, 3, 6, 3),
    }
    for (p, q), expected in cases.items():
        lay = study.layout_numbers(study.torus(0.45, 1.0, 1200, p, q))
        assert (lay["crossings"], lay["determinant"], lay["ramps"], lay["symmetry"]) == expected
        assert lay["least_ramps"] == expected[2]
    fig8 = study.layout_numbers(study.figure_eight_knot(1.0, 1200))
    assert (fig8["crossings"], fig8["determinant"]) == (4, 5)


def test_sticks_and_lattice() -> None:
    # section 3.4: the six-stick trefoil's corners turn 128.8 degrees, their planes 46.1 and 20.1 apart; the lattice
    # trefoil has 24 steps and 14 sticks in 8 planar runs, 582 mm round at 26.2 mm with corners of 7.9 mm
    turn, normals = study.corner_turns(study.six_stick_trefoil(1.0))
    assert np.degrees(turn) == pytest.approx(np.full(6, 128.8), abs=0.05)
    planes = {round(math.degrees(math.acos(abs(float(normals[i] @ normals[(i + 1) % 6])))), 1) for i in range(6)}
    assert planes == {20.1, 46.1}
    lattice = study.lattice_trefoil(1.0)
    assert float(np.abs(np.roll(lattice, -1, 0) - lattice).sum()) == pytest.approx(24.0)
    _, lnorm = study.corner_turns(lattice)
    assert (len(lattice), len(study.plane_runs(lnorm))) == (14, 8)
    arr = study.compute_arrangements(quick=True)
    assert arr["L26"]["length_mm"] == pytest.approx(582, abs=1.5)
    assert arr["sticks"]["path_mm"] == pytest.approx(291, abs=1.0)


def test_what_a_joint_and_a_bend_hold() -> None:
    # section 3.1: the trefoil is flat on its inner equator at r/R = 4/13; holding 16.6 D through sharp joints takes
    # 51 sticks. Section 6.1: the static force maxima, four rods 0.6997 D/r0 and facing boards 1.209 D/h in their
    # plane; toward a board q0 sec^2 reaches 0.908 at 0.47h for q0 = 0.5, where the force is 5.15 D/h
    ratios = np.linspace(0.29, 0.33, 41)
    least = [float(study.torus(float(r), 1.0, 1500).curvature.min()) for r in ratios]
    assert float(ratios[int(np.argmin(least))]) == pytest.approx(4 / 13, abs=0.002)
    assert next(m for m in range(6, 200) if 1 / math.sin(math.radians(720 / m)) ** 2 >= 16.6) == 51
    boards = study.boards_well()
    assert boards["rods_force_max"] == pytest.approx(0.6997, abs=2e-4)
    assert boards["force_max_in"] == pytest.approx(1.209, abs=2e-3)
    assert boards["edge"] == pytest.approx(0.47, abs=0.005)
    assert boards["force_at_edge"] == pytest.approx(5.15, abs=0.03)


def test_fields_between_two_planes() -> None:
    # sections 4.2 and 4.6: per volt on the axis, a segment's Q gives pi/4h^2, Y and X 1/2h, and a common-mode step
    # puts (1/2h) sech(pi s/2h) on the axis; odd disturbances fade as exp(-pi s/h), even ones as exp(-pi s/2h)
    nx, ns, step = 256, 1024, 1 / 16
    s = (np.arange(ns) - ns // 2) * step
    pats = study.segment_patterns(nx, ns, step, -4.0, 4.0)
    fields = {k: study.planes_axis(*v, step) for k, v in pats.items()}
    mid = ns // 2
    assert fields["Q"]["gradient"][mid] / (math.pi / 4) == pytest.approx(1.0, abs=0.005)
    assert abs(fields["Y"]["Ey"][mid]) == pytest.approx(0.5, abs=0.003)
    assert abs(fields["X"]["Ex"][mid]) == pytest.approx(0.5, abs=0.005)
    end = study.planes_axis(*study.segment_patterns(nx, ns, step, 0.0, ns * step)["C"], step)
    assert np.abs(end["Es"])[np.abs(s) < 3] == pytest.approx(0.5 / np.cosh(math.pi * s[np.abs(s) < 3] / 2), abs=2e-3)


def test_the_slot_guide_and_the_cut_out() -> None:
    # section 3.3: the slot guide in a 0.8 mm FR-4 board with a 1.0 mm slot holds 0.26 to 0.34 D0 and 0.50 to 0.59 G0
    # flush, 0.14 D0 and 0.37 G0 with its copper set back 0.2 mm; a cut-out of a gapless stack's middle board leaves
    # (h/(2h + t/2))^2 G0
    flush = study.slot_section(0.8, 1.0, 0.0, 4.4, 1.5, None, 0.02)
    assert 0.26 <= flush["depth_D0"] <= 0.34
    assert 0.50 <= flush["gradient_G0"] <= 0.59
    back = study.slot_section(0.8, 1.0, 0.2, 4.4, 1.5, None, 0.02)
    assert back["depth_D0"] == pytest.approx(0.14, abs=0.015)
    assert back["gradient_G0"] == pytest.approx(0.37, abs=0.02)
    cut = study.cutout(0.5, 0.4, 6.0, 0.1)
    assert cut["inside_gradient"] == pytest.approx(cut["predicted_inside"], rel=0.03)


def test_the_ring_drive() -> None:
    # section 5.3: section A's ring at 915 MHz (Q 251, Z0 73.5 ohm by the line study), R = Z0 Q / (2 pi m) = 2934
    # ohm; with six compensated feeds two cheap modes, 3632 V/A near sqrt(K) R / 2, over a floor 124 times lower; one
    # feed makes a standing wave, two compiled feeds a travelling one
    case = study.RingCase("compensated", study.F915, 1.13, 251.0, 73.5, 1, 6, compensated=True)
    ring = case.ring(6, None)
    assert ring.shunt == pytest.approx(2934, rel=2e-3)
    freqs = ring.resonance() * (1 + np.linspace(-0.01, 0.01, 101))
    sv = np.array([np.linalg.svd(g, compute_uv=False) for g in ring.operator(freqs)]) / math.sqrt(ring.samples)
    top = sv[int(np.argmax(sv[:, 0]))]
    assert top[0] == pytest.approx(3632, rel=0.01)
    assert top[1] == pytest.approx(top[0], rel=1e-3)
    assert top[0] / top[2] == pytest.approx(124, rel=0.05)
    one = study.best_drive(case.ring(1, None))
    two = study.best_drive(case.ring(2, None))
    assert one["ripple"] > 0.9
    assert two["ripple"] < 0.02


def test_multipactor_lines_and_the_chirp() -> None:
    # sections 5.5 and 6.4: 1 eV at q = 0.5 is 50.3 V and 5.2 dB below ECSS's copper boundary; the six-gap pi-mode
    # ramp on T915-S6 takes 1 to 8 eV in 117 laps, 32.7 us, on a chirp of 3 f_rev from 5.60 to 15.84 MHz, while the
    # tune falls from 92 to 32; six 1 mm dips of 15 % cost about 0.16 a band crossed at 60 meV a lap
    fd = math.sqrt(1.0 / (1.42 * 0.25))
    volts = 35.7 * 0.5 * fd**2
    ecss = float(np.interp(fd, *np.array(study.ECSS_CU).T))
    assert volts == pytest.approx(50.3, abs=0.1)
    assert 20 * math.log10(ecss / volts) == pytest.approx(5.2, abs=0.05)
    out = study.compute_acceleration(quick=True)
    ramp = out["T915-S6"]
    assert ramp["laps"] == 117
    assert ramp["time_us"] == pytest.approx(32.7, abs=0.1)
    assert (ramp["chirp_MHz"][0], ramp["chirp_MHz"][-1]) == pytest.approx((5.60, 15.84), abs=0.02)
    assert (ramp["tune"][0], ramp["tune"][-1]) == pytest.approx((92, 32), abs=0.5)
    assert out["design"]["G"] == pytest.approx(0.16, abs=0.02)
    assert study.mathieu_beta(0.5) * study.F915 / 2 / 1e6 == pytest.approx(171.0, abs=0.1)


def test_section_b_on_the_revision_grid_and_converged() -> None:
    # section 5.1: sectionB.py's uniform grid of 22.9 um cells holds the 25 um polyimide as 22.9 um, which gave its
    # n_eff 1.0312 and share 0.0730; with the copper's 18 um and the polyimide's 25 um resolved n_eff is 1.031, the
    # share 0.074 and Z0 41.7 ohm; Wheeler's rule gives Q_c 228 (section 5.4), the charge rule about twice that
    old = study.facing_line(0.915, 0.025)
    assert old["layer_um"] == pytest.approx(22.9, abs=0.05)
    assert (old["C_pF"], old["index"], old["share"], old["z0"]) == pytest.approx(
        (82.50, 1.0312, 0.0730, 41.7), rel=2e-3
    )
    new = study.facing_section(study.H915, study.F915, fine=1 / 8, coarse=0.1, growth=0.1)
    assert new["index"] == pytest.approx(1.031, abs=5e-4)
    assert new["share"] == pytest.approx(0.0737, abs=1e-3)
    assert new["z0"] == pytest.approx(41.7, abs=0.15)
    assert new["Qc"] == pytest.approx(228, rel=0.02)
    assert new["Q"] == pytest.approx(214, rel=0.02)
    assert old["Qc"] / new["Qc"] > 1.9
    assert pt_length(new["index"]) == pytest.approx(study.L915 / study.MM, abs=0.1)


def pt_length(index: float) -> float:
    """m = 1 at 915 MHz (mm) for an index."""
    return float(pt.C0 / (index * study.F915) / study.MM)


def test_wheeler_rule_on_stripline() -> None:
    # section 3.3's calibration: Wheeler's rule on the graded grid against stripline's closed form, within 2 %
    for width, thick, spacing in study.STRIPLINES:
        assert study.stripline_wheeler(width, thick, spacing)["ratio"] == pytest.approx(1.0, abs=0.02)


def test_the_slot_guide_as_a_line() -> None:
    # section 3.3: slotline.py's FR-4 row (n_eff 1.71, Z0 60, share 0.845, its charge rule's Q_c 292, lam_g 192 mm,
    # 0.335 G0 and 0.103 D0 on the axis); Wheeler's rule on 35 um copper gives Q_c 158, so Q 43 on FR-4 and 117 on
    # RO4003C, R(m = 1) 412 and 1243 ohm, and 40 V takes 1.9 and 0.64 W at m = 1, 1.29 W at m = 2 on RO4003C
    fr4 = study.slot_line(0.8, 4.4, 0.02, 0.5, 1.0)
    assert (fr4["C_pF"], fr4["index"], fr4["z0"], fr4["share"], fr4["Qc"], fr4["lam_g_mm"]) == pytest.approx(
        (94.8, 1.708, 60.1, 0.845, 292, 192), rel=2e-3
    )
    assert (fr4["gradient_G0"], fr4["barrier_D0"]) == pytest.approx((0.335, 0.103), abs=1e-3)
    ro = study.slot_line(0.813, 3.55, 0.0027, 0.5, 1.0)
    for row, board, q, shunt in ((fr4, 0.8, 43, 412), (ro, 0.813, 117, 1243)):
        qc = study.thick_slot_qc(board, 0.5, 1.0)
        assert qc == pytest.approx(158, rel=0.02)
        quality = 1 / (1 / qc + 1 / row["Qd"])
        assert quality == pytest.approx(q, abs=1)
        assert row["z0"] * quality / (2 * math.pi) == pytest.approx(shunt, rel=0.01)
    shunt_m2 = 1243 / 2
    assert 1600 / (2 * shunt_m2) == pytest.approx(1.29, abs=0.01)


REFERENCE_DRAWS = (  # lens5/ring.py's evaluate on ring915.py's draws (seed 7, Q 381, 2 %): 2Q|c_2m|, hybrid |b/a|,
    # depth per watt, re-solved |b/a|, depth per watt, second feed's amplitude and phase
    (1.1506357521841937, 0.18386959811125966, 0.6416204123378344, 0.010289723479643607, 0.9567970854499657, 0.625,
     -92.43243243243244),
    (0.7521990026695583, 0.12458477498807383, 0.7496002751500863, 0.00851173697262015, 0.948048373509738, 1.375,
     82.70270270270268),
    (1.5871757061862706, 0.2607672135196098, 0.51062172628541, 0.0016668356801840837, 0.9443285844754605, 1.875,
     -180.0),
)  # fmt: skip


def test_the_ring_is_the_buildability_lens_model() -> None:
    # section 5.3: the study's ring of 192 segments draws and scores as lens5/ring.py did for ring915.py
    length = pt.C0 / (1.031 * study.F915)
    rng = np.random.default_rng(7)
    flat = np.full(study.RING_NODES, 50.0), np.full(study.RING_NODES, 1.031), np.zeros(study.RING_NODES)
    ref = study.ring_score(length, 1, *flat, 381.0, False)
    for expected in REFERENCE_DRAWS:
        z, n, rel = study.ring_draw(rng, length, 0.02, 1.031, 0.01, 0.0)
        fixed = study.ring_score(length, 1, z, n, rel, 381.0, False)
        comp = study.ring_score(length, 1, z, n, rel, 381.0, True)
        got = (
            2 * 381.0 * study.ring_harmonic(z),
            fixed["purity"],
            fixed["score"] / ref["score"],
            comp["purity"],
            comp["score"] / ref["score"],
            comp["w_abs"],
            comp["w_deg"],
        )
        assert got == pytest.approx(expected, rel=1e-7, abs=1e-9)


def test_the_ring_at_quick_settings() -> None:
    # sections 5.2 and 5.3: at 2 % rms spacer error two re-solved feeds leave a backward wave near 0.01 and keep 0.94
    # of a perfect ring's depth per watt where the fixed hybrid leaves |b/a| above 0.05; at Q 214 the spacer rule is
    # 0.87 linewidths and one feed loses 2.9 % of the depth a lap
    out = study.compute_ring(quick=True)
    two = out["cases"]["bare copper, 2 %"]
    assert two["resolved_purity_median"] < 0.03
    assert two["resolved_depth_median"] > 0.9
    assert two["hybrid_purity_median"] > 0.05
    assert out["Q"]["bare"] == pytest.approx(214, rel=0.08)
    quality = 214
    assert 2 * quality * 0.5 * 0.02 / math.sqrt(24) == pytest.approx(0.87, abs=0.005)
    assert 1 - math.exp(-2 * math.pi / quality) == pytest.approx(0.029, abs=5e-4)


def test_the_bunch_and_the_controls() -> None:
    # section 6.6: a bunch's dipole per unit charge on T915-S's trefoil, 4.75, 15.5, 4.60, 0.28, 1.65 and 0.68 mm at
    # k = 1 to 6, axial at 3 and 6; magnetic 1046, 488 and 775 mm^2 at k = 1 to 3; vector areas 2971 and 8037 mm^2,
    # 0.37; the ring's dipole 35.8 mm at k = 1 alone. Section 3.2: 4_1 at 318 mm bends to 6.1 mm with strands 8.4 mm
    # apart and 6.61 pi; the (1, 3) curve at r/R 0.2 to 25.1 mm, 3.64 pi, holonomy -146 degrees
    knot = study.torus(study.RATIO915, study.L915, 3000)
    tref = study.bunch_multipoles(knot)
    dip = [r["dipole_mm"] for r in tref["harmonics"]]
    assert dip == pytest.approx([4.75, 15.49, 4.60, 0.28, 1.65, 0.68], abs=0.01)
    axial = [r["axial_mm"] > 0.1 for r in tref["harmonics"]]
    assert axial == [False, False, True, False, False, True]
    assert [r["magnetic_mm2"] for r in tref["harmonics"][:3]] == pytest.approx([1046, 488, 775], abs=1.5)
    ring = study.bunch_multipoles(study.scaled(pt.ClosedCurve.circle(1.0, 3000), study.L915))
    assert (tref["area_abs_mm2"], ring["area_abs_mm2"]) == pytest.approx((2971, 8037), abs=2)
    assert tref["area_abs_mm2"] / ring["area_abs_mm2"] == pytest.approx(0.37, abs=0.005)
    assert ring["harmonics"][0]["dipole_mm"] == pytest.approx(35.8, abs=0.05)
    assert max(r["dipole_mm"] for r in ring["harmonics"][1:]) < 1e-3
    assert study.strand_clearance(knot) / study.MM == pytest.approx(18.5, abs=0.05)
    controls = study.control_curves(study.L915, 3000)
    fig8 = controls["4₁, common form"]
    assert (fig8.bend_radius / study.MM, fig8.total_curvature / math.pi) == pytest.approx((6.10, 6.61), abs=0.01)
    assert study.strand_clearance(fig8) / study.MM == pytest.approx(8.44, abs=0.05)
    helix = controls["T(1, 3) at r/R 0.2"]
    assert helix.bend_radius / study.MM == pytest.approx(25.05, abs=0.05)
    assert helix.total_curvature / math.pi == pytest.approx(3.644, abs=0.005)
    assert math.degrees(helix.holonomy) == pytest.approx(-145.7, abs=0.3)


def test_the_odd_k_band() -> None:
    # section 3.5: closed with k = -3 the plates are one band of 2L, 639 mm, laid flat in about 44 x 570 mm, and the
    # k = 0 strip's compression is 4.6 times the wrinkling scale, k = -3's 29; k = -1's band overlaps itself
    out = study.compute_strips(quick=True)
    k3, k1, k0 = (out["framings"][k] for k in (-3, -1, 0))
    assert k3["band_mm"] == pytest.approx(639, abs=1.5)
    assert k3["band_box_mm"] == pytest.approx((44, 571), abs=(1.5))
    assert k1["band_clearance_mm"] < study.WIDTH915 / study.MM
    assert k3["band_clearance_mm"] > study.WIDTH915 / study.MM
    assert (k0["wrinkle_ratio"], k3["wrinkle_ratio"]) == pytest.approx((4.6, 29.2), abs=0.15)
    assert out["bending_strain"] == pytest.approx(0.0036, abs=1e-4)


def test_the_revision_arithmetic() -> None:
    # sections 3 to 9: 95 % at 95 % confidence takes 59 draws, 93 with one failure; T915-S6's chirp 5.60 to 15.84 MHz
    # in 117 laps, 32.7 us, segments 53.0 mm; the flat ring 2.7 to 3.3 times T915-S; S2-C3's lobes 30.0 and 42.2 mm
    # and its 418 mm layout at R0 = 33 mm; the critical feed capacitor 0.47 pF for 1.42 kOhm (0.35 pF for 2.54 kOhm)
    out = study.compute_numbers(quick=True)
    assert out["draws"] == {0: 59, 1: 93, 2: 124}
    s6 = out["T915-S6"]
    assert s6["chirp_MHz"] == pytest.approx((5.60, 15.84), abs=0.005)
    assert (s6["time_us"], s6["segment_mm"]) == pytest.approx((32.7, 53.0), abs=0.05)
    assert out["held"]["flat ring"]["over"] == pytest.approx((2.72, 3.31), abs=0.005)
    assert [x / study.MM for x in study.lobe_bends(33.2e-3, 4e-3)] == pytest.approx([30.0, 42.2], abs=0.05)
    assert study.layout_length(33e-3, 4e-3) / study.MM == pytest.approx(418.1, abs=0.05)
    assert study.feed_capacitor(1417) * 1e12 == pytest.approx(0.466, abs=0.001)
    assert study.feed_capacitor(2538) * 1e12 == pytest.approx(0.347, abs=0.001)
    assert study.half_db(30.24) == 30.0
    assert study.half_db(37.21) == 37.0
    t = out["T915-S"]
    assert (t["two_w_margin_mm"], t["free_space"]) == pytest.approx((0.65, 0.970), abs=0.005)
    assert t["pickups_mm"] == pytest.approx(53.0, abs=0.05)
    assert out["rails"][0.575e-3]["depth"] == pytest.approx(0.0927, abs=5e-5)
    # the power chain on section B at its converged Q_c 228, Q_d 3394 and Z0 41.65 ohm: Q 214, R 1.42 kOhm, 1 eV
    # (50.27 V) in 0.89 W, +30.0 dBm through the coupling; ENIG's Q_c / 1.3 and / 5, 1.14 to 4.23 W
    for factor, power in ((1.0, 0.891), (1.3, 1.142), (5.0, 4.232)):
        q = 1 / (factor / 228.2 + 1 / 3394)
        shunt = 41.65 * q / (2 * math.pi)
        assert 50.27**2 / (2 * shunt) == pytest.approx(power, abs=0.004)


def test_section_a_and_the_well() -> None:
    # section 5.4: section A on slotted RO4003C (35 um copper) has share 0.19, Z0 73.5 ohm and Q 251 by Wheeler's rule,
    # so 1 eV at its drawn 54.9 V takes 0.51 W at m = 1 and 1.03 W at m = 2; the well's solver reads ELE-2's gapless
    # boards, pi/4 and 1/4
    a = study.facing_section(
        study.HA,
        study.F915,
        half=2.0,
        split=1.0,
        fine=1 / 8,
        coarse=0.1,
        growth=0.1,
        **study.SECTION_A["RO4003C 0.51 mm, slotted"],
    )
    assert a["share"] == pytest.approx(0.187, abs=0.003)
    assert a["z0"] == pytest.approx(73.5, abs=0.2)
    assert a["Q"] == pytest.approx(251, rel=0.02)
    shunt = a["z0"] * a["Q"] / (2 * math.pi)
    assert (54.9**2 / (2 * shunt), 54.9**2 / shunt) == pytest.approx((0.51, 1.03), rel=0.03)
    gapless = study.facing_section(
        1e-3,
        study.F915,
        half=12.0,
        split=0.02,
        copper=2e-6,
        laminate=0.0,
        fine=1 / 8,
        coarse=0.1,
        growth=0.1,
        well=True,
    )
    assert gapless["gamma"] == pytest.approx(math.pi / 4, abs=0.003)
    assert gapless["depth_share"] == pytest.approx(0.25, abs=0.001)


def test_the_split_as_a_lever() -> None:
    # section 5.4: at the same 1 eV well at q = 0.5 and 915 MHz, the well solved at each split, the strips take 0.88 W
    # at a split of 0.3h and 0.33 W at h (the volts 49.5 and 39.6 V, h 0.90 and 0.72 mm)
    grid = {"fine": 1 / 8, "coarse": 0.1, "growth": 0.1}
    narrow, wide = study.lever_job((0.3, grid)), study.lever_job((1.0, grid))
    assert (narrow["P"], wide["P"]) == pytest.approx((0.88, 0.33), rel=0.05)
    assert (narrow["V"], wide["V"]) == pytest.approx((49.5, 39.6), abs=0.3)
    assert (narrow["h_mm"], wide["h_mm"]) == pytest.approx((0.90, 0.72), abs=0.01)
    assert wide["P"] < 0.5 * narrow["P"]


def test_the_rigid_and_high_bands() -> None:
    # sections 5.4 and 5.5: section B on rigid FR-4 0.8 mm at 2450 MHz, m = 2, has Q 52 and takes 10.8 W for 1 eV at
    # 50.3 V; on 25 um polyimide at 5800 MHz, m = 5, it takes about 13 W
    fr4 = study.facing_section(
        study.FH / 2.45e9, 2.45e9, fine=1 / 8, coarse=0.1, growth=0.1, **study.LAMINATES["FR-4 0.8 mm"]
    )
    assert fr4["Q"] == pytest.approx(52, rel=0.03)
    assert 50.27**2 / (2 * fr4["z0"] * fr4["Q"] / (4 * math.pi)) == pytest.approx(10.8, rel=0.05)
    high = study.facing_section(
        study.FH / 5.8e9, 5.8e9, fine=1 / 8, coarse=0.1, growth=0.1, **study.LAMINATES["polyimide 25 µm"]
    )
    assert 50.27**2 / (2 * high["z0"] * high["Q"] / (10 * math.pi)) == pytest.approx(12.8, rel=0.05)


def test_the_strips_at_2450_mhz_as_the_lens_ran_them() -> None:
    # section 5.5: lens5's 12-draw Monte Carlo of the 2.45 GHz strips at 10 um rms on h (2.9 %), Q 250: median
    # 2Q|c_2m| 0.993, depth per watt 0.934 re-solved against 0.654 fixed, as lens5/ring.py found; |c_2m| does not
    # depend on Q, so at Wheeler's Q 134.5 the median is 0.53 and 2Q|c_2m| <= 1 allows 5.5 %, 37 um of the gap
    lens = study.strips_2450(250.0, 1.067, 12)
    assert lens["split_median"] == pytest.approx(0.9934, abs=1e-4)
    assert (lens["resolved_depth_median"], lens["hybrid_depth_median"]) == pytest.approx((0.9342, 0.6543), abs=1e-4)
    split = lens["split_median"] * 134.5 / 250
    assert lens["sigma_pct"] / split == pytest.approx(5.5, abs=0.05)
    assert 10 / split * 2 == pytest.approx(37, abs=0.5)


def test_numbers_in_titles() -> None:
    assert study.sci(2.2e-4) == "2.2 × 10⁻⁴"
    assert study.sci(0.46) == "0.46"
    assert study.minus(-0.085, "+.3f") == "−0.085"


def test_a_quick_run_writes_its_figures(tmp_path: Path) -> None:
    assert study.main(["--quick", "--only", "knots,strips,acceleration", "--out", str(tmp_path)]) == 0
    assert sorted(p.name for p in tmp_path.glob("*.png")) == ["acceleration.png", "knots.png", "strips.png"]
