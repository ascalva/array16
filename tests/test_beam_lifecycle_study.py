"""The beam-lifecycle study's fast parts, held to the numbers docs/BEAM_LIFECYCLE.md states."""

from __future__ import annotations

import math
from pathlib import Path

import beam_lifecycle_study as study
import numpy as np
import pytest

from array16 import particles as pt


def test_getting_electrons_in() -> None:
    # section 3.2 and 3.3: 4.3 nA per mW from magnesium as installed; one-turn fills of 24, 7.2 and 1.0 uA, the
    # last 0.24 W of 265 nm light as installed and 2.4 mW cleaned, the 932 mm ring at 50 eV 1.7 W and 17 mW
    out = study.compute_sources(quick=True)
    assert 2e-5 * out["photons_per_mW"] * pt.E * 1e9 == pytest.approx(4.3, abs=0.05)
    fills = out["fills"]
    assert [row["current_A"] * 1e6 for row in fills.values()] == pytest.approx([24, 7.2, 1.0], abs=0.35)
    assert fills["932 mm, 1 eV"]["lap_s"] == pytest.approx(1.57e-6, abs=0.01e-6)
    assert fills["932 mm, 1 eV"]["mW_magnesium"] == pytest.approx(240, abs=5)
    assert fills["932 mm, 1 eV"]["mW_cleaned"] == pytest.approx(2.4, abs=0.05)
    assert fills["932 mm, 50 eV"]["mW_magnesium"] == pytest.approx(1700, abs=20)


def test_cooling_and_heating_times() -> None:
    # section 3.4: resistive cooling 0.16 to 1.4 ms in a segment and 7.3 to 22 ms round a lap; gas heats 0.1 eV to
    # the 1 eV rim in 0.33 to 0.42 ms at 1e-5 mbar, a hundred times slower at 1e-7; cooling wins below 2e-6 to
    # 2e-5 mbar at 0.1 eV and 7e-7 to 6e-6 mbar at 1 eV
    out = study.compute_sources(quick=True)
    assert out["cooling_s"]["segment"] == pytest.approx((0.16e-3, 1.4e-3), rel=0.03)
    assert out["cooling_s"]["lap"] == pytest.approx((7.35e-3, 22e-3), rel=0.005)
    low, high = out["to_rim_s"]["1e-05"]
    assert (low, high) == pytest.approx((0.33e-3, 0.42e-3), rel=0.01)
    assert out["to_rim_s"]["1e-07"][0] == pytest.approx(100 * low, rel=1e-9)
    assert out["crossover_mbar"]["0.1 eV"] == pytest.approx((2.2e-6, 2.0e-5), rel=0.03)
    assert out["crossover_mbar"]["1 eV"] == pytest.approx((6.9e-7, 6.2e-6), rel=0.03)


def test_a_collision_adds_2u_to_first_order_and_more_in_exact_mathieu_motion() -> None:
    # section 3.4: from rest along the guide a collision adds 2.02, 2.25 and 2.88 U at q = 0.1, 0.3 and 0.5, U half
    # the secular energy; repeated collisions keep 82 to 85 % of the energy transverse. As q -> 0 the first tends to
    # the Proposition's 2U, and the steady state, from the first-order map [[5/3, 2/3], [1/3, 1/3]], to 2 sqrt(2/3)
    # with a transverse share sqrt(2/3)
    gains = [study.collision_gains(q) for q in (0.1, 0.3, 0.5)]
    assert [g["from_rest"] for g in gains] == pytest.approx([2.02, 2.25, 2.88], abs=0.005)
    assert [g["steady"] for g in gains] == pytest.approx([1.66, 1.88, 2.51], abs=0.005)
    assert [g["transverse_share"] for g in gains] == pytest.approx([0.818, 0.828, 0.852], abs=0.001)
    small = study.collision_gains(0.02)
    assert small["from_rest"] == pytest.approx(2.0, abs=2e-3)
    assert small["steady"] == pytest.approx(2 * math.sqrt(2 / 3), abs=2e-3)
    assert small["transverse_share"] == pytest.approx(math.sqrt(2 / 3), abs=1e-3)
    beta, _, _ = study.floquet(0.3)
    assert beta == pytest.approx(0.21606, abs=1e-5)  # Mathieu's characteristic exponent at a = 0, q = 0.3


def test_the_four_gap_ramps() -> None:
    # section 4: h = 13, 58.5 -> 112.5 MHz in 20.5 us over 135 laps, 2.635 MHz per us, 45 % captured at 1.6 eV rms;
    # h = 4, 18.0 -> 34.6 MHz, 49 % at 2.7 eV rms; buckets +-1.83 and +-3.52 eV; Q_s 0.189 and 0.098; f_s 0.852 MHz;
    # the ramp's 126 m of path against a 40 m free path at 1e-5 mbar leaves 4.5 % unscattered, 73 % and 97 % at 1e-6
    # and 1e-7 mbar
    out = study.compute_phase(quick=False)
    four, thirteen = out["h = 4"], out["h = 13"]
    assert (four["f_start_MHz"], four["f_end_MHz"]) == pytest.approx((18.0, 34.6), abs=0.05)
    assert four["time_us"] == pytest.approx(20.5, abs=0.05)
    assert four["laps"] == thirteen["laps"] == 135
    assert (four["captured"], four["spread_eV"]) == pytest.approx((0.49, 2.7), abs=0.05)
    assert (thirteen["f_start_MHz"], thirteen["f_end_MHz"]) == pytest.approx((58.5, 112.5), abs=0.05)
    assert (out["path_m"], out["free_path_m"]["1e-05"]) == pytest.approx((126, 40.5), abs=0.3)
    assert list(out["unscattered"].values()) == pytest.approx([0.045, 0.73, 0.97], abs=0.005)
    assert thirteen["chirp_MHz_per_us"] == pytest.approx(2.635, abs=0.002)
    assert (thirteen["captured"], thirteen["spread_eV"]) == pytest.approx((0.45, 1.6), abs=0.05)
    assert thirteen["captured_count"] == 1814
    assert (thirteen["bucket_50"], thirteen["bucket_185"]) == pytest.approx((1.83, 3.52), abs=0.005)
    assert (thirteen["tune_50"], thirteen["tune_185"]) == pytest.approx((0.189, 0.098), abs=5e-4)
    assert thirteen["f_s_MHz"] == pytest.approx(0.852, abs=5e-4)
    # the captured electrons' first periods, quartiles in laps, against 1/Q_s = 5.3 at 50 eV
    assert (thirteen["first_period_low"], thirteen["first_period_high"]) == pytest.approx((5.74, 6.74), abs=0.05)
    assert (four["spread_slope"], thirteen["spread_slope"]) == pytest.approx((0.5, 0.5), abs=0.02)


def test_speed_and_fields() -> None:
    # section 5 and 6: beta at 1 eV to 79 keV; q grows 1.03, 1.06, 1.24, 1.64 and 2 times; 62 V/m and 5.5 nT at
    # 0.5 mm from 1e7 electrons at 185 eV; 4 fC on a quarter, a 1.5 uA pulse of 2.7 ns, 73 uV on 50 ohm
    out = study.compute_speed(quick=True)
    assert list(out["beta"].values()) == pytest.approx([0.00198, 0.0063, 0.0198, 0.0269, 0.0625, 0.195, 0.5], rel=0.01)
    d = out["doppler"]
    assert (d["185 eV"]["n_eff 1"][0], d["185 eV"]["n_eff 2"][0]) == pytest.approx((1.03, 1.06), abs=0.005)
    assert (d["10000 eV"]["n_eff 1"][0], d["10000 eV"]["n_eff 2"][0]) == pytest.approx((1.24, 1.64), abs=0.005)
    assert d["79000 eV"]["n_eff 1"][0] == pytest.approx(2.0, abs=0.01)
    assert (out["E_at_half_mm"], out["B_at_half_mm_nT"]) == pytest.approx((62, 5.5), abs=0.5)
    assert out["pulse_fC"] == pytest.approx(4.0, abs=0.05)
    assert (out["pulse_uA"], out["pulse_ns"], out["pulse_uV"]) == pytest.approx((1.46, 2.73, 73), rel=0.01)


def test_each_quarter_takes_a_quarter_of_the_image() -> None:
    # INJ-2: the four weighting potentials are equal at the centre and sum to one
    assert study.quarter_shares(0.5e-3) == pytest.approx([0.25] * 4, abs=1e-6)


def test_resonances() -> None:
    # section 7: Madey's best detuning 2.606 and conductance 0.424 N I0 / V0; start currents 1.1 uA, 11 uA, 0.25 mA
    # and 3.9 mA; Smith-Purcell at 5.93 GHz for 100 eV over 1 mm, the band from 0.46 to 102 eV there; over the 135
    # laps of the ramp the tune sweeps 23.6 -> 12.3 in the 1 GHz guide and 135.6 -> 70.5 in the 5.75 GHz guide; one
    # 1 mm joint at 10 % grows the 5.75 GHz guide's action 1.75 times by 185 eV, at 5 % by 4 %, and the 1 GHz guide's
    # by about 2 %
    out = study.compute_resonances(quick=True)
    assert (out["madey_psi"], out["madey_constant"]) == pytest.approx((2.606, 0.424), abs=1e-3)
    assert list(out["start_A"].values()) == pytest.approx([1.1e-6, 11e-6, 0.25e-3, 3.9e-3], rel=0.03)
    assert out["ring_current_A"] == pytest.approx(10e-6, rel=0.03)
    assert out["band_eV"]["1 mm"][0] == pytest.approx(0.46, abs=0.01)
    assert out["band_eV"]["1 mm"][1] == pytest.approx(102, abs=0.5)
    assert out["band_eV"]["10 mm"] == pytest.approx((45.5, 10.6e3), rel=0.01)
    assert float(study.speed(100.0)) / 1e-3 / 1e9 == pytest.approx(5.93, abs=0.005)
    final = out["final_growth"]
    assert final["5.75 GHz facing boards, 610 MHz, 10 %"] == pytest.approx(1.75, abs=0.01)
    assert final["5.75 GHz facing boards, 610 MHz, 5 %"] == pytest.approx(1.04, abs=0.01)
    assert final["1 GHz guide, 106 MHz secular, 10 %"] == pytest.approx(1.02, abs=0.01)
    ends = out["tune_range"]
    assert ends["1 GHz guide, 106 MHz secular"] == pytest.approx((23.6, 12.3), abs=0.05)
    assert ends["5.75 GHz facing boards, 610 MHz"] == pytest.approx((135.6, 70.5), abs=0.05)
    assert len(out["_tunes"]["1 GHz guide, 106 MHz secular"]) == 135


def test_what_the_beam_can_tell_us() -> None:
    # section 9: 0.40 and 0.012 um per V/m; the Earth's 50 uT moves the beams 4.8 and 12 um; argon's minimum of
    # 0.125 pi a0^2 at 0.285 eV and 11 us at 10 eV; a Sagnac phase of 0.057 rad round the circle and 0.021 for the
    # trefoil; one radian of Aharonov-Bohm phase per 9.5 fT
    out = study.compute_world(quick=True)
    assert list(out["per_field_um"].values()) == pytest.approx([0.40, 0.012], rel=0.02)
    assert list(out["earth_shift_um"].values()) == pytest.approx([4.8, 12], rel=0.03)
    d_term, p_term, residual = study.fit_argon()
    assert residual < 1e-9
    energies = np.linspace(0.2, 0.4, 2001)
    sigma = study.argon_cross_section(energies, d_term, p_term)
    assert float(energies[np.argmin(sigma)]) == pytest.approx(0.285, abs=1e-3)
    assert float(sigma.min()) == pytest.approx(0.125 * math.pi * study.BOHR**2, rel=1e-3)
    assert out["argon"]["estimate_10eV_s"] == pytest.approx(11e-6, rel=0.03)
    assert out["argon"]["hold_s"]["0.285 eV"] == pytest.approx(11.6e-3, rel=0.02)
    assert (out["sagnac_circle_rad"], out["sagnac_trefoil_rad"]) == pytest.approx((0.057, 0.021), abs=0.001)
    assert out["areas_cm2"]["932 mm circle"] == pytest.approx(691, abs=1)
    assert out["fT_per_radian"] == pytest.approx(9.5, abs=0.05)
    assert out["quantum_ueV"] == pytest.approx(0.44, abs=0.005)


def test_a_quick_run_writes_its_figures(tmp_path: Path) -> None:
    assert study.main(["--quick", "--only", "sources,speed", "--out", str(tmp_path)]) == 0
    assert sorted(p.name for p in tmp_path.glob("*.png")) == ["sources.png", "speed.png"]
