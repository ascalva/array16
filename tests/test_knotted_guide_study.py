"""The knotted-guide study's fast parts, held to the numbers docs/KNOTTED_GUIDE.md states."""

from __future__ import annotations

import math
from pathlib import Path

import knotted_guide_study as study
import numpy as np
import pytest


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
    # section 5.3: section A's ring at 915 MHz, R = Z0 Q / (2 pi m); with six feeds two cheap modes near
    # sqrt(K) R / 2 over a floor about a hundred times lower; one feed makes a standing wave, two compiled feeds a
    # travelling one
    case = study.DRIVE_CASES[1]
    ring = case.ring(6, None)
    assert ring.shunt == pytest.approx(73 * 193 / (2 * math.pi), rel=1e-12)
    freqs = ring.resonance() * (1 + np.linspace(-0.01, 0.01, 101))
    sv = np.array([np.linalg.svd(g, compute_uv=False) for g in ring.operator(freqs)]) / math.sqrt(ring.samples)
    top = sv[int(np.argmax(sv[:, 0]))]
    assert top[0] == pytest.approx(2778, rel=0.01)
    assert top[1] == pytest.approx(top[0], rel=1e-3)
    assert 80 < top[0] / top[2] < 120
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


def test_numbers_in_titles() -> None:
    assert study.sci(2.2e-4) == "2.2 × 10⁻⁴"
    assert study.sci(0.46) == "0.46"
    assert study.minus(-0.085, "+.3f") == "−0.085"


def test_a_quick_run_writes_its_figures(tmp_path: Path) -> None:
    assert study.main(["--quick", "--only", "knots,strips,acceleration", "--out", str(tmp_path)]) == 0
    assert sorted(p.name for p in tmp_path.glob("*.png")) == ["acceleration.png", "knots.png", "strips.png"]
