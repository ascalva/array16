"""The virtual-antenna study's fast parts, held to the numbers docs/VIRTUAL_ANTENNA.md states."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import virtual_antenna_study as study

from array16 import particles as pt


def test_the_model_is_calibrated() -> None:
    # one dipole radiates Z0 k^2 |Il|^2 / 12 pi, and the closed-form radiation matrix is the far-field quadrature's
    out = study.selftest()
    assert out["dipole_power_closed"] < 1e-12
    assert max(v for key, v in out.items() if key.startswith("R ")) < 1e-6


def test_the_tile_radiates_through_few_drives_at_low_frequency() -> None:
    # section 4.2: eigenvalues of R within 10 dB of the largest, 3 / 6 / 14; condition 1.5e7 / 5.4e4 / 51
    counts, conditions = [], []
    for f in (2.3e9, 3.5e9, 5.8e9):
        ev = np.linalg.eigvalsh(study.radiation(study.panel(f)))
        counts.append(int((ev > 0.1 * ev[-1]).sum()))
        conditions.append(float(ev[-1] / ev[0]))
    assert counts == [3, 6, 14]
    assert conditions == pytest.approx([1.5e7, 5.4e4, 51], rel=0.03)


def test_the_tiles_price_list_in_a_ball() -> None:
    # section 4.1: a ball of radius 52 mm 150 mm in front at 5.8 GHz, 10 % loss: 0, 4.8, 5.7, 11.5, 14.9 dB
    s2 = study.price_list(study.panel(5.8e9), (0.0, 0.0, 0.15), study.C0 / 5.8e9, study.LOSS)
    prices = study.db(s2[0] / s2)
    assert prices[:5] == pytest.approx([0.0, 4.8, 5.7, 11.5, 14.9], abs=0.06)


def test_fine_structure_costs_a_power_law_in_the_spacing() -> None:
    # section 4.4: six alternating values at 0.15 lambda on a 64-element shell, relative to one
    shell = study.shell(5.8e9, 64, 0.25)
    lam = shell.wavelength
    power = np.array([study.alternating_cost(shell, (0, 0, 0), 0.15 * lam, m, 2, 0.0) for m in range(1, 7)])
    assert study.db(power / power[0])[1:] == pytest.approx([10.6, 23.1, 36.1, 48.7, 61.1], abs=0.06)


def test_a_focus_against_the_ceiling() -> None:
    # section 4.5: the whips at 433 MHz, 100 mm above the centre with 10 % loss, give 4.9 times RFG-2's ceiling,
    # 382 micro-eV per watt, and 8.9 times the ceiling over the ground, 0.55 of RFG-2's for a horizontal field there
    whips = study.whips(0.433e9)
    assert study.focus_share(whips, (0.0, 0.0, 0.1), study.LOSS, free_space=True) == pytest.approx(4.88, abs=0.01)
    assert study.focus_share(whips, (0.0, 0.0, 0.1), study.LOSS) == pytest.approx(8.85, abs=0.06)
    assert study.half_space_factors(whips.k, 0.1)[0] == pytest.approx(0.55, abs=0.01)
    assert study.focus_depth_ev(whips, (0.0, 0.0, 0.1), study.LOSS) == pytest.approx(382e-6, rel=0.01)
    assert pt.free_space_depth_ceiling_ev(1.0) == pytest.approx(78.2e-6, rel=1e-3)


def test_a_dense_shell_passes_the_ceiling_by_its_stored_energy() -> None:
    # section 4.5: sixteen elements on a shell stay near the ceiling (0.97 lossless and 0.72 with 10 % loss at
    # kb = pi, at most 1.57 and 1.04 at kb >= 3, reached at kb = 3); 256, dense enough to close, pass it beside the
    # cavity resonance at kb = 6.117 by 939 times lossless and 12.4 times with 10 % loss: the case that passes it
    def share(count: int, kb: float, loss: float) -> float:
        return study.focus_share(study.shell(kb * study.C0 / (2 * np.pi * 0.25), count, 0.25), (0, 0, 0), loss)

    assert (share(16, np.pi, 0.0), share(16, np.pi, study.LOSS)) == pytest.approx((0.97, 0.72), abs=0.005)
    assert (share(16, 3.0, 0.0), share(16, 3.0, study.LOSS)) == pytest.approx((1.57, 1.04), abs=0.005)
    assert share(16, 3.2, 0.0) < share(16, 3.0, 0.0)
    assert share(256, 6.1, 0.0) == pytest.approx(939, rel=0.01)
    assert share(256, 6.1, study.LOSS) == pytest.approx(12.4, abs=0.05)


def test_a_closed_shell_stores_the_dipole_wave_at_its_cavity_resonances() -> None:
    assert study.cavity_resonances(10.0) == pytest.approx([2.744, 6.117, 9.317], abs=1e-3)


def test_the_alphabets_line_ceilings() -> None:
    # section 6.2: rev A at 5.8 GHz 0.900 on line 1 and 0.707 on line 0; rev B at 2.3 GHz 0.981 and 0.808, gap 72.3
    rev_a = study.hull_numbers(study.alphabet_phases("rev A", 5.8e9))
    rev_b = study.hull_numbers(study.alphabet_phases("rev B", 2.3e9))
    assert (rev_a["line1"], rev_a["inradius"], rev_a["largest_gap_deg"]) == pytest.approx(
        (0.9003, 0.7071, 90.0), abs=1e-4
    )
    assert (rev_b["line1"], rev_b["inradius"]) == pytest.approx((0.981, 0.808), abs=1e-3)
    assert rev_b["largest_gap_deg"] == pytest.approx(72.3, abs=0.05)


def test_the_harmonic_weight_law() -> None:
    # section 6.2: a 2-bit staircase puts sinc^2(1/4) = 0.811 (-0.91 dB) on line 1 and -10.5 dB on line -3;
    # section 6.3: Frank shifts give the 16 lines -12.0 to -16.0 dB of the full-array gain, 0.773 in all
    out = study.compute_sequences(quick=True)
    lines = list(out["lines"])
    stair = study.db(out["staircase"]["power_per_line"])
    assert stair[lines.index(1)] == pytest.approx(-0.91, abs=0.01)
    assert stair[lines.index(-3)] == pytest.approx(-10.5, abs=0.05)
    frank = out["Frank shifts"]["best_beam_per_line"][:16]
    assert float(np.sum(frank)) == pytest.approx(0.773, abs=1e-3)
    assert (study.db(frank.max()), study.db(frank.min())) == pytest.approx((-12.0, -16.0), abs=0.05)


def test_the_guide_radiates_only_its_selected_orders() -> None:
    # section 7.2: the shares of every order sum to one; 3 stubs at m = 2 give Q_rad 1.1e3, 8 stubs at m = 5 give 20
    assert sum(float(study.order_power(order, [2.0])[0]) for order in range(-30, 31)) == pytest.approx(1.0, abs=1e-6)
    assert study.guide_q(2, 3) == pytest.approx(1.06e3, rel=0.01)
    assert study.guide_q(5, 8) == pytest.approx(20.4, rel=0.01)


def test_containment_through_electrodes() -> None:
    # section 7.3: 89.1 eV per watt on four rods and 22.3 on two facing boards at 1 GHz, 0.5 mm, 10 kOhm, up to
    # 2.9 and 1.2 eV, where the null's q reaches 0.908; 1 eV takes 15 V and 11 mW at q = 0.53 on four rods, and one
    # board would need q = 1.73, past the edge. At q = 0.5: 16 V and 13 mW at 1.07 GHz, 50 V and 126 mW at 1.68 GHz,
    # 542 V and 14.7 W at 3.46 GHz. One board keeps at most 0.0091 of a four-rod well, its null where it was asked.
    out = study.compute_containment(quick=True)
    four, faced, one = (out["at_1GHz"][k] for k in ("four rods", "two facing boards", "one board, two rails"))
    assert four["eV_per_W"] == pytest.approx(89.1, abs=0.05)
    assert faced["eV_per_W"] == pytest.approx(22.3, abs=0.05)
    assert (four["edge_eV"], faced["edge_eV"], one["edge_eV"]) == pytest.approx((2.89, 1.17, 0.276), abs=0.005)
    assert (four["V_for_1eV"], four["mW_for_1eV"], four["q_for_1eV"]) == pytest.approx((15.0, 11.2, 0.534), abs=0.05)
    assert four["q_for_1eV"] < 0.908 < one["q_for_1eV"]
    assert one["q_for_1eV"] == pytest.approx(1.73, abs=0.005)
    halves = [row[key] for row in (four, faced, one) for key in ("q_half_V", "q_half_mW", "q_half_GHz")]
    assert halves == pytest.approx([16.0, 12.8, 1.07, 50.3, 126.3, 1.68, 541.7, 14670, 3.46], rel=3e-3)
    assert out["rails"]["share"] == pytest.approx(0.0091, abs=5e-5)
    assert out["rails"]["null_mm"] == pytest.approx(0.5, abs=1e-9)


def test_charges_are_heard_on_their_electrodes() -> None:
    # section 7.4: one electron across 2 mm is 142 H; the column, at the limiting density in a radius r0/2 = 0.25 mm,
    # holds a quarter of N'_max = 2.2e8 per metre, 5.5e7 per metre and 1.5e7 on 277 mm (6.1e7 at q = 0.6, radiation
    # Q 2.0e6), and drives -38 dBm into 50 Ohm, 51 dB above its radiation
    out = study.compute_charges(quick=True)
    assert out["one_electron_H"] == pytest.approx(142, abs=0.5)
    assert out["column_radius_m"] == pytest.approx(0.25e-3)
    assert out["per_metre"] == pytest.approx(out["full_per_metre"] / 4, rel=1e-12)
    assert (out["full_per_metre"], out["per_metre"], out["count"]) == pytest.approx((2.19e8, 5.48e7, 1.52e7), rel=0.01)
    assert (out["count_q06"], out["radiation_Q_q06"]) == pytest.approx((6.07e7, 1.97e6), rel=0.01)
    signals = out["signals_dBm"]
    assert signals["column, coherent: induced into 50 Ohm"] == pytest.approx(-38, abs=0.5)
    gap = signals["column, coherent: induced into 50 Ohm"] - signals["column, coherent: radiated"]
    assert gap == pytest.approx(51, abs=0.5)


def test_the_v1_predictions() -> None:
    # section 10.1: the spectrum in R + rI is 0, -5.9, -22.3, -28.6 dB; the centre -31.9 dB, a wavelength across
    # -2.4 dB (the case that must fail); half a wavelength out, away beats across by 17 dB (V1-6 asks 10)
    out = study.compute_bench(quick=True)
    assert out["generalised_db"] == pytest.approx([0.0, -5.9, -22.3, -28.6], abs=0.06)
    table = out["table"]
    assert table["centre"]["held_db"] == pytest.approx(-31.9, abs=0.06)
    assert table["across 1"]["held_db"] == pytest.approx(-2.4, abs=0.06)
    assert table["across 0.5"]["held_db"] - table["away 0.5"]["held_db"] == pytest.approx(17, abs=0.5)


def test_numbers_in_titles() -> None:
    assert study.sci(9.4e-5) == "9.4 × 10⁻⁵"
    assert study.sci(1.14e6) == "1.1 × 10⁶"
    assert study.sci(0.0091) == "0.0091"
    assert study.sci(142.0) == "142"


def test_a_quick_run_writes_its_figures(tmp_path: Path) -> None:
    assert study.main(["--quick", "--only", "compiler,alphabets,selection", "--out", str(tmp_path)]) == 0
    assert sorted(p.name for p in tmp_path.glob("*.png")) == ["alphabets.png", "compiler.png", "selection.png"]
