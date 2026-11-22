"""experiments/twin_openems_line.py, R1's line in full wave against its 2-D section (DR-0), on what needs no openEMS:
the board it models, the ring-down fit's calibration, openEMS's single-precision loss, the face integral's limit, the
2-D section stage on the QUICK grading against the record the script wrote, and DR-0's verdicts recomputed from that
record. The full-wave stages run under openEMS (the .env environment), and tests/test_twin_dr0.py holds the twin to
their record."""

from __future__ import annotations

import json
import math
from typing import Any

import pytest
import twin_openems_line as fw

RECORD = fw.OUT / "results.json"


def _record() -> dict[str, Any]:
    """The record the script wrote: in the tree, so a missing one fails."""
    out: dict[str, Any] = json.loads(RECORD.read_text())
    return out


def test_the_rail_line_is_r1s_board() -> None:
    # ELE-1's a = 0.40 and b = 1.46 mm, edges c1 = 0.20 and c2 = 1.66 mm, the gapless null sqrt(c1 c2) = 0.576 mm up;
    # 0.15 mm gaps centred on the edges leave 0.25 mm of strip and 1.31 mm of rail; 35 um copper on 1.6 mm FR-4; the
    # board as written has no plane
    line = fw.RailLine()
    assert (line.strip, line.rail, line.gap, line.board, line.er, line.tand) == (0.40, 1.46, 0.15, 1.6, 4.4, 0.02)
    assert line.copper == pytest.approx(0.035, rel=1e-12)
    assert line.edges == pytest.approx((0.2, 1.66), abs=1e-12)
    assert line.null_gapless == pytest.approx(0.5762, abs=1e-4)
    (s0, s1), (r0, r1), (g0, g1) = line.spans
    assert (s0, s1, r0, r1, g0, g1) == pytest.approx((0.0, 0.125, 0.275, 1.585, 1.735, 25.0), abs=1e-12)
    assert line.backed
    assert not fw.AS_WRITTEN.backed
    with pytest.raises(ValueError, match="gap < strip"):
        fw.RailLine(gap=0.5)


def test_the_ring_down_fit_recovers_a_damped_tone() -> None:
    # rule 5's calibration of the matrix pencil: a Q 60 tone at 914.3 MHz under a second tone and 1e-3 noise, on the
    # probes' 88 ps samples, recovered to 2.5e-6 in f and 6.2e-4 in Q, as the record has it
    got = fw.pencil_check()
    assert abs(got["f_error"]) < 3e-6
    assert abs(got["Q_error"]) < 1e-3
    assert got == pytest.approx(_record()["pencil"], rel=1e-9)


def test_openems_single_precision_rounds_the_laminates_loss() -> None:
    # at the fine mesh's dt the laminate's damping per step is 18.563 units of 2^-24, held as 19: the run's FR-4 loses
    # 1.0235 of the given loss, which the record takes out
    q = fw.quantization(fw.RailLine(), 9.62273e-15)
    assert q["damping_steps"] == pytest.approx(18.563, abs=1e-3)
    assert q["factor"] == pytest.approx((19 / q["damping_steps"]), rel=1e-6)
    assert q["factor"] == pytest.approx(1.0235, abs=1e-4)
    assert _record()["straight"]["fine"]["quantization"]["factor"] == pytest.approx(q["factor"], rel=1e-6)


def test_the_face_integrals_limit_is_the_check_that_does_not_use_wheelers_rule() -> None:
    # an error falling 2^(-1/3) per halving, the corners' r^(-2/3), extrapolates exactly by both rules; on R1's backed
    # board the face integral's Q_c at t/8 to t/64 cells, 119.9, 112.8, 107.8 and 103.9, extrapolates to 88.8
    # (Richardson at 2^(-1/3)) and 90.3 (Aitken's on the last three), 5.0 and 3.5 % under Wheeler's rule's 93.5
    exact = [90.0 + 30.0 * 2 ** (-k / 3) for k in range(4)]
    clean = fw.face_limit(exact)
    assert clean["richardson"] == pytest.approx(90.0, rel=1e-12)
    assert clean["aitken"] == pytest.approx(90.0, rel=1e-12)
    assert clean["order"] == pytest.approx(1 / 3, rel=1e-9)
    record = _record()
    section = record["section"]
    faces = [section[f"faces t/{div}"]["Qc_faces"] for div in (8, 16, 32, 64)]
    assert faces == pytest.approx([119.92, 112.82, 107.78, 103.87], abs=0.01)
    limit = fw.face_limit(faces)
    wheeler = section["built"]["Qc"]
    assert (limit["richardson"], limit["aitken"], wheeler) == pytest.approx((88.81, 90.25, 93.52), abs=0.01)
    assert 0.03 < 1 - limit["aitken"] / wheeler < 1 - limit["richardson"] / wheeler < 0.06
    stored = record["dr0"]["copper"]
    assert stored["faces"] == faces
    assert stored["wheeler"] == wheeler
    assert {key: stored[key] for key in limit} == pytest.approx(limit, rel=1e-12)
    with pytest.raises(ValueError, match="three cells"):
        fw.face_limit(faces[:2])


def test_the_section_stage_on_the_quick_grading_is_the_records() -> None:
    # the FDTD's own cross-section on the QUICK grading, solved again (0.2 s): C', n_eff, Z0, the share, Q_d, the face
    # integral's Q_c and Q as the record stored them, to rounding
    found = fw.section_job(("fdtd quick", fw.R1, fw.QUICK, False, False, True))
    stored = _record()["section"]["fdtd quick"]
    for key in ("C_pF", "Ca_pF", "index", "z0", "share", "Qd", "Qc_faces", "Q"):
        assert found[key] == pytest.approx(stored[key], rel=1e-9), key


def test_dr_0_s_verdicts_follow_from_the_record() -> None:
    # dr0() from the stored stages alone gives back the verdicts the record keeps: the straight line on both meshes,
    # the ring's two modes against the two-rail model, the copper's extrapolation
    record = _record()
    again = fw.dr0(record)
    assert set(again) == set(record["dr0"]) == {"straight", "ring", "copper"}
    assert json.loads(json.dumps(again)) == record["dr0"]
    for mode in again["ring"]["fine"]:
        assert abs(mode["f_dev"]) < 5e-3 < abs(mode["f_dev_single"])
    assert math.isfinite(again["straight"]["fine"]["Q_cal"])
