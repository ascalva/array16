"""The first run's fast parts, held to its numbers."""

from __future__ import annotations

import accelerator_study as study
import pytest


def test_radiated_field_is_the_far_field_formula() -> None:
    # 38.9 mW of EIRP, 1 m away: 1.08 V/m rms, 1.53 V/m peak
    assert study.radiated_field(15.9, 1.0) == pytest.approx(1.528, rel=1e-3)
    assert study.radiated_field(15.9, 0.3) == pytest.approx(5.09, rel=1e-3)


def test_scales_give_the_field_each_depth_takes() -> None:
    bands = study.scales()["bands"]
    assert bands["2.44 GHz"]["field_for_1eV_kVm"] == pytest.approx(73.1, rel=1e-3)
    assert bands["0.433 GHz"]["volts_across_1mm_for_1eV"] == pytest.approx(13.0, rel=5e-3)
    assert bands["0.433 GHz"]["depth_at_0.3m_eV"] == pytest.approx(1.5e-7, rel=0.05)


def test_four_gaps_on_one_chirp_reproduce_the_first_run() -> None:
    out = study.ring()
    assert out["harmonic"] == 13
    assert out["phase_step_deg"] == pytest.approx(90.0)
    assert (out["f_start_MHz"], out["f_end_MHz"]) == pytest.approx((58.5, 112.6), abs=0.1)
    assert out["laps"] == pytest.approx(135, abs=1)
    assert out["end_eV"] == pytest.approx(185, abs=1)
    assert out["captured"] == pytest.approx(0.45, abs=0.01)
    assert out["spread_eV"] == pytest.approx(1.6, abs=0.1)


def test_a_beams_own_numbers() -> None:
    own = study.own({"end_eV": 185.0, "length_mm": 932.0})
    assert own["space_charge"]["q 0.6"]["density_m3"] == pytest.approx(1.12e15, rel=0.01)
    assert own["pickup"]["1e+06"]["dBm"] == pytest.approx(-78, abs=1)
    assert own["zigzag_mm"] == pytest.approx(0.65, abs=0.01)
