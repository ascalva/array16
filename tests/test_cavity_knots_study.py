"""The cavity-knots study's fast parts, held to the numbers docs/CAVITY_KNOTS.md states."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import cavity_knots_study as study
import numpy as np
import pytest

from array16 import particles as pt


def test_bessel_functions_have_their_known_zeros() -> None:
    # zeros of J_m and J_m' (the cylinder's TM and TE modes), as tabulated to 15 digits
    assert study.bessel_zeros(0, 2) == pytest.approx([2.404825557695773, 5.520078110286311], abs=1e-12)
    assert study.bessel_zeros(1, 1) == pytest.approx([3.831705970207512], abs=1e-12)
    assert study.bessel_zeros(1, 2, derivative=True) == pytest.approx([1.841183781340659, 5.331442773525033], abs=1e-12)
    assert study.bessel_zeros(2, 2, derivative=True) == pytest.approx([3.054236928227140, 6.706133194158459], abs=1e-12)
    # J_0' = -J_1, so their zeros agree once the origin is left out
    assert study.bessel_zeros(0, 2, derivative=True) == pytest.approx(study.bessel_zeros(1, 2), abs=1e-12)


def test_spherical_bessel_functions_have_their_known_zeros() -> None:
    # j_1 vanishes where tan x = x; j_2 and j_3 at their tabulated first zeros; j_2's second zero
    assert study.sphere_zeros(1, 1) == pytest.approx([4.493409457909064], abs=1e-12)
    assert study.sphere_zeros(2, 2) == pytest.approx([5.763459196894550, 9.095011330476355], abs=1e-12)
    assert study.sphere_zeros(3, 1) == pytest.approx([6.987932000500520], abs=1e-12)
    assert study.sphere_zeros(0, 2) == pytest.approx([math.pi, 2 * math.pi], abs=1e-12)


def test_bessel_values_against_closed_forms_and_the_recurrence() -> None:
    x = np.linspace(0.05, 30, 400)
    j = study.bessel_ladder(list(range(-2, 9)), x)
    assert float(study.bessel_j(0, 1.0)) == pytest.approx(0.7651976865579666, abs=1e-15)
    assert float(study.bessel_j(1, 1.0)) == pytest.approx(0.4400505857449335, abs=1e-15)
    for n in range(-1, 8):  # J_(n+1) = (2n/x) J_n - J_(n-1)
        assert np.abs(j[n + 3] - (2 * n / x * j[n + 2] - j[n + 1])).max() < 1e-12
    spherical = (np.sin(x) / x, np.sin(x) / x**2 - np.cos(x) / x, (3 / x**2 - 1) * np.sin(x) / x - 3 * np.cos(x) / x**2)
    for order, closed in enumerate(spherical):
        tail = x > 2  # the closed forms lose digits near the origin; the series does not
        assert np.abs(study.spherical_j(order, x)[tail] - closed[tail]).max() < 1e-12
    assert float(study.spherical_j(3, 0.01)) == pytest.approx(0.01**3 / 105 * (1 - 0.01**2 / 18), rel=1e-12)


def test_harmonics_are_orthonormal_and_harmonic() -> None:
    for degree in (1, 2, 3, 4):
        basis = study.Harmonics.of(degree)
        points, weights = study.sphere_quadrature(2 * degree + 2)
        values = basis.coefficients @ study.monomial_values(basis.exponents, points)
        assert (values * weights) @ values.T == pytest.approx(np.eye(2 * degree + 1), abs=1e-12)
        _, _, hess = basis.values(np.ones(2 * degree + 1), np.random.default_rng(degree).normal(size=(5, 3)))
        assert np.abs(np.trace(hess, axis1=1, axis2=2)).max() < 1e-10


def test_a_sphere_family_vanishes_on_diameters_counted_by_morse() -> None:
    # section 4: within one family the nulls are diameters over the critical points of A; max - saddles + min = 2
    for degree, diameters in ((1, 1), (2, 3)):
        out = study.compute_sphere(degree, 1, 11, quick=True, feeds=False)
        assert out["morse_sum"] == 2
        assert out["diameters"] == diameters
    owner = study.compute_sphere(3, 1, 2026, quick=True, feeds=False)
    assert owner["morse_sum"] == 2


def test_the_sphere_families_at_2_45_ghz() -> None:
    # section 3's table: ka, radius and copper Q of the TE and TM families l = 1, 2, 3
    table = study.compute_charts(quick=True)["table"]
    expected = {
        "TE l=1": (4.4934, 8.75, 65300),
        "TE l=2": (5.7635, 11.22, 83800),
        "TE l=3": (6.9879, 13.61, 101600),
        "TM l=1": (2.7437, 5.34, 29300),
        "TM l=2": (3.8702, 7.54, 33700),
        "TM l=3": (4.9734, 9.69, 37200),
    }
    for name, (ka, radius, q) in expected.items():
        assert table[name]["ka"] == pytest.approx(ka, abs=1e-4)
        assert table[name]["radius_cm"] == pytest.approx(radius, abs=0.005)
        assert table[name]["Q"] == pytest.approx(q, rel=2e-3)


@pytest.fixture(scope="module")
def worked() -> dict[str, Any]:
    return study.compute_worked(quick=True)


def test_the_worked_cylinder(worked: dict[str, Any]) -> None:
    # section 7's table: the cavity, its Q's, its field at 1 W
    assert worked["d_over_a"] == pytest.approx(1.3780, abs=5e-5)
    assert worked["ka"] == pytest.approx(7.0831, abs=5e-5)
    assert (worked["radius_cm"], worked["height_cm"], worked["volume_l"]) == pytest.approx(
        (13.79, 19.01, 11.36), abs=6e-3
    )
    assert (worked["Q"]["TE113"], worked["Q"]["TE221"], worked["Q"]["mix"]) == pytest.approx(
        (55967, 89953, 81658), abs=1
    )
    assert worked["share_TE221"] == pytest.approx(0.83, abs=0.005)
    assert worked["linewidth_kHz"] == pytest.approx(30.0, abs=0.05)
    assert worked["E_rms_1W"] == pytest.approx(10.3e3, rel=5e-3)
    assert worked["peak_over_mean"] == pytest.approx(7.87, abs=0.005)


def test_the_worked_cylinders_null_loops(worked: dict[str, Any]) -> None:
    # three closed loops 34.15, 34.15 and 31.68 cm long, z from 0.171 to 0.829 d, at least 1.92 cm from any wall,
    # two folds each, unlinked; two open lines from side wall to side wall, 15.64 cm, and one from cap to cap, 19.96 cm
    loops = worked["loops"]
    assert [row["length_cm"] for row in loops] == pytest.approx([34.15, 34.15, 31.68], abs=0.01)
    assert all(row["z_over_d"] == pytest.approx((0.171, 0.829), abs=1e-3) for row in loops)
    assert min(row["wall_cm"] for row in loops) == pytest.approx(1.92, abs=0.005)
    assert [row["folds"] for row in loops] == [2, 2, 2]
    assert max(abs(v) for v in worked["links"].values()) < 1e-6
    lines = sorted((line["ends"], line["length_cm"]) for line in worked["open_lines"])
    assert [end for end, _ in lines] == ["cap to cap", "side wall to side wall", "side wall to side wall"]
    assert [length for _, length in lines] == pytest.approx([19.96, 15.64, 15.64], abs=0.06)


def test_the_worked_cylinders_well(worked: dict[str, Any]) -> None:
    # the escape saddle at s = 0.01919 of the mean |E|^2, rho = 0.545 a, z/d = 0.096 and 0.904; 375 micro-eV per
    # watt, 4.80 times RFG-2's ceiling, 25 meV for 67 W and 1 eV for 2.66 kW; Mathieu q 6.5e-4 and 5.3e-4 at 1 W,
    # 0.021 and 0.017 at 1 kW; B's quiver at most 369 beta^2 of the depth
    assert worked["saddle_share"] == pytest.approx(0.01919, abs=5e-6)
    saddle = worked["saddles"][0]
    assert math.hypot(saddle[0], saddle[1]) == pytest.approx(0.545, abs=1e-3)
    heights = sorted(min(p[2], 1.3780 - p[2]) / 1.3780 for p in worked["saddles"])
    assert heights == pytest.approx([0.096] * 4, abs=1e-3)
    assert worked["depth_ueV_per_W"] == pytest.approx(375.4, abs=0.2)
    assert worked["ceiling_ratio"] == pytest.approx(4.80, abs=0.005)
    assert worked["W_for_25meV"] == pytest.approx(67, abs=0.5)
    assert worked["kW_for_1eV"] == pytest.approx(2.66, abs=0.01)
    assert list(worked["q_1W"]) == pytest.approx([6.5e-4, 5.3e-4], abs=0.06e-4)
    assert list(worked["q_1kW"]) == pytest.approx([0.021, 0.017], abs=6e-4)
    assert worked["B_ratio"][0] == pytest.approx((0.09, 2.66), abs=0.01)
    assert worked["magnetic_over_beta2"] == pytest.approx(369, abs=1)


def test_folds_are_quadratic_across_and_quartic_along(worked: dict[str, Any]) -> None:
    # CAV-6 at the worked example's folds: exponents 2 in the normal plane and 4 along the tangent near the fold
    out = study.compute_folds(quick=True, worked=worked)
    for fold in out["folds"]:
        exponents = list(fold["near"].values())
        assert exponents == pytest.approx([2.0, 2.0, 4.0], abs=0.02)
    assert out["soft_range"] == pytest.approx((2.70, 9.08), abs=0.01)
    assert out["knot_turns"] == {"trefoil": 4, "Hopf link": 4}
    assert abs(out["hopf_linking"]) == pytest.approx(1.0, abs=1e-3)


def test_the_tolerance_table(worked: dict[str, Any]) -> None:
    # section 8.3: relative mode phase 47.8, 10.2, 5.1 and 2.0 degrees; 1019, 204, 102 and 41 nm on the radius
    out = study.compute_tolerance(quick=True, worked=worked)
    assert out["slope"] == pytest.approx(0.829, abs=5e-4)
    table = out["table"]
    assert [row["phase_deg"] for row in table] == pytest.approx([47.8, 10.2, 5.1, 2.0], abs=0.06)
    assert [row["radius_nm"] for row in table] == pytest.approx([1019, 204, 102, 41], abs=0.6)
    assert table[0]["floor_over_s"] == pytest.approx(14.8, rel=0.05)
    assert out["shared_phase_power"] == pytest.approx(4.3, abs=0.02)
    # a floor under a tenth of the depth: a split under 0.032 linewidths, 66 nm on the radius
    assert (out["tenth_split"], out["tenth_nm"]) == pytest.approx((0.032, 66), rel=0.03)


def test_depth_against_power_and_the_thresholds(worked: dict[str, Any]) -> None:
    # section 8.1 and 8.5: facing boards 3.7 eV and four rods 14.8 eV per watt at 2.45 GHz; the side wall at 0.942
    # of the mean, 18.4 meV per watt; air near 11 kW; a few Torr at 0.1 to 1 W
    out = study.compute_power(quick=True, worked=worked)
    routes = out["routes"]
    assert routes["facing boards, h = 0.5 mm, 10 kOhm (RFG-3)"]["eV_per_W"] == pytest.approx(3.7, abs=0.02)
    assert routes["four rods, h = 0.5 mm, 10 kOhm (RFG-3)"]["eV_per_W"] == pytest.approx(14.8, abs=0.06)
    assert routes["fields carried in from far away (RFG-2)"]["eV_per_W"] == pytest.approx(
        pt.free_space_depth_ceiling_ev(1.0), rel=1e-12
    )
    assert out["wall_peak"] == pytest.approx(0.942, abs=1e-3)
    assert out["wall_meV_per_W"] == pytest.approx(18.4, abs=0.05)
    assert out["air_kW"] == pytest.approx(10.8, abs=0.05)
    assert out["few_torr_W"] == pytest.approx((0.12, 1.09), abs=0.01)


def test_numbers_in_titles() -> None:
    assert study.sci(9.4e-5) == "9.4 × 10⁻⁵"
    assert study.sci(0.0091) == "0.0091"


def test_a_quick_run_writes_its_figures(tmp_path: Path) -> None:
    assert study.main(["--quick", "--only", "charts", "--out", str(tmp_path)]) == 0
    assert sorted(p.name for p in tmp_path.glob("*.png")) == ["charts.png"]
