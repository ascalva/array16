"""The secondary-emission instrument (array16/guide/emission.py): Furman and Pivi's model held to their Table I and
Table II and to its own algebra, its samplers to quadrature, the 2-D tracker held to the closed-form flights and to
multipactor's one-dimensional model in the parallel-plate limit, the chambers of the protocol's sections, and the
protocol's readings on a coarse grid (docs/KNOTTED_GUIDE.md sections 8.1, 10 E2-13 and 11 MPC-H1)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from array16.guide import emission as em
from array16.guide import multipactor as mp
from array16.guide import sections as gs
from array16.particles import Arr

ROOT = Path(__file__).resolve().parents[1]
QUICK = {"electrons": 300, "periods": 16, "settle": 6, "steps": 80}  # small populations, few periods


def _quadrature_mean(shape: float, limit: float, points: int = 200_001) -> tuple[float, float]:
    """The mean of Gamma(shape) held under `limit`, and its CDF at half the limit, by the trapezoid rule."""
    x = np.linspace(0.0, limit, points)
    log = (shape - 1) * np.log(np.maximum(x, 1e-300)) - x
    density = np.exp(log - log.max())
    density[0] = 0.0 if shape > 1 else density[0]
    total = np.trapezoid(density, x)
    half = x <= limit / 2
    return float(np.trapezoid(x * density, x) / total), float(np.trapezoid(density[half], x[half]) / total)


# ---- Furman and Pivi's copper --------------------------------------------------------------------------------------
def test_copper_is_furman_and_pivis_table() -> None:
    fp = em.FP_COPPER
    # Table I, copper: backscattered, rediffused, true secondaries, the emitted angular spectrum
    assert (fp.p1e_inf, fp.p1e_peak, fp.e_e, fp.w, fp.p, fp.sigma_e, fp.e1, fp.e2) == (
        0.02, 0.496, 0.0, 60.86, 1.0, 2.0, 0.26, 2.0,
    )  # fmt: skip
    assert (fp.p1r_inf, fp.e_r, fp.r, fp.q, fp.r1, fp.r2) == (0.2, 0.041, 0.104, 0.5, 0.26, 2.0)
    assert (fp.delta_ts, fp.e_ts, fp.s, fp.t1, fp.t2, fp.t3, fp.t4, fp.alpha) == (
        1.8848, 276.8, 1.54, 0.66, 0.8, 0.7, 1.0, 1.0,
    )  # fmt: skip
    # Table II, copper: p_n and eps_n (eV), n = 1 to 10
    assert fp.p_n == (2.5, 3.3, 2.5, 2.5, 2.8, 1.3, 1.5, 1.5, 1.5, 1.5)
    assert fp.eps_n == (1.5, 1.75, 1.0, 3.75, 8.5, 11.5, 2.5, 3.0, 2.5, 3.0)
    # the total their Table I lists: 2.1 at 271 eV
    at, top = fp.peak
    assert at == pytest.approx(271.0, abs=1.0)
    assert top == pytest.approx(2.1, abs=0.01)
    assert top == pytest.approx(2.0936, abs=1e-4)
    # at zero energy only backscattering is left, at P1e^; far above the peak each component meets its floor
    assert float(fp(0.0)) == pytest.approx(0.496, abs=1e-12)
    assert float(fp.rediffused(0.0)) == 0.0
    assert float(fp.true(0.0)) == 0.0
    assert float(fp.elastic(1e5)) == pytest.approx(0.02, abs=1e-9)
    assert float(fp.rediffused(1e5)) == pytest.approx(0.2, rel=1e-3)
    assert fp.crossover == pytest.approx(27.32, abs=0.01)
    assert float(fp(fp.crossover)) == pytest.approx(1.0, abs=1e-9)
    # against the chart's yield, ECSS's Table 9-1 average: one at 27.3 eV, not 19; 2.09, not 2.48; 0.496 at zero, not 0.8
    assert fp.crossover > mp.COPPER.e1 + 8
    assert top < mp.COPPER.dmax - 0.35
    assert float(fp(0.0)) < mp.COPPER.d0 - 0.3


def test_the_angle_raises_each_component_as_their_equations_47_and_48() -> None:
    fp = em.FP_COPPER
    for energy in (5.0, 40.0, 300.0, 2000.0):
        assert float(fp.elastic(energy, 0.0)) == pytest.approx(1.26 * float(fp.elastic(energy)), rel=1e-12)
        assert float(fp.rediffused(energy, 0.0)) == pytest.approx(1.26 * float(fp.rediffused(energy)), rel=1e-12)
        c = 0.5  # 60 degrees: dts^ (1 + 0.66 (1 - 0.5^0.8)) at Ets^ (1 + 0.7 (1 - 0.5))
        peak, where = 1.8848 * (1 + 0.66 * (1 - 0.5**0.8)), 276.8 * 1.35
        x = energy / where
        assert float(fp.true(energy, c)) == pytest.approx(peak * 1.54 * x / (0.54 + x**1.54), rel=1e-12)
    # grazing impacts free more at every energy above a few eV
    energies = np.geomspace(5.0, 5e3, 40)
    assert bool(np.all(fp(energies, 0.2) > fp(energies, 1.0)))


def test_conditioned_copper_scales_the_true_secondaries_alone() -> None:
    bare, cond = em.FP_COPPER, em.FP_CONDITIONED
    assert cond.peak[1] == pytest.approx(1.2, abs=1e-6)
    assert cond.delta_ts == pytest.approx(0.99097, abs=1e-4)  # 0.5258 of 1.8848
    energies = np.geomspace(0.5, 5e3, 50)
    assert np.array_equal(cond.elastic(energies), bare.elastic(energies))
    assert np.array_equal(cond.rediffused(energies), bare.rediffused(energies))
    ratio = cond.true(energies) / bare.true(energies)
    assert float(np.ptp(ratio)) < 1e-12
    assert cond.crossover == pytest.approx(90.43, abs=0.05)
    with pytest.raises(ValueError, match="alone reach"):
        bare.conditioned(0.3)
    with pytest.raises(ValueError, match="p_n and eps_n"):
        em.FurmanPivi(p_n=(1.0, 2.0))
    with pytest.raises(ValueError, match="penetrate"):
        em.FurmanPivi(p1e_peak=0.7)


@settings(max_examples=80, deadline=None)
@given(st.floats(0.0, 1e4), st.floats(0.0, 1.0))
def test_the_multiplicity_is_a_distribution_whose_mean_is_the_yield(energy: float, cos: float) -> None:
    for surface in (em.FP_COPPER, em.FP_CONDITIONED):
        p = surface.multiplicity(energy, cos)
        assert p.shape == (em.MULTIPLICITY + 1,)
        assert bool(np.all(p >= 0))
        assert float(p.sum()) == pytest.approx(1.0, abs=1e-12)
        assert float(p @ np.arange(em.MULTIPLICITY + 1)) == pytest.approx(float(surface(energy, cos)), abs=1e-12)


# ---- the samplers --------------------------------------------------------------------------------------------------
@pytest.mark.parametrize(("shape", "limit"), [(2.5, 0.3), (2.5, 5.0), (7.8, 3.0), (69.0, 20.0), (69.0, 68.5), (69.0, 90.0)])
def test_the_truncated_gamma_is_the_quadratures(shape: float, limit: float) -> None:
    rng = np.random.default_rng(7)
    count = 40_000
    draws = em._truncated_gamma(rng, np.full(count, shape), np.full(count, limit))
    assert bool(np.all((draws >= 0) & (draws <= limit)))
    mean, half = _quadrature_mean(shape, limit)
    assert float(draws.mean()) == pytest.approx(mean, abs=4 * float(draws.std()) / math.sqrt(count) + 1e-9)
    share = float(np.mean(draws <= limit / 2))
    assert share == pytest.approx(half, abs=4 * math.sqrt(max(half * (1 - half), 1e-6) / count) + 1e-4)


def test_the_backscattered_deficit_is_a_half_gaussian_under_the_impact() -> None:
    rng = np.random.default_rng(3)
    count = 40_000
    far = em._deficit(rng, np.full(count, 300.0), 2.0)
    assert float(far.mean()) == pytest.approx(2.0 * math.sqrt(2 / math.pi), rel=0.02)
    near = em._deficit(rng, np.full(count, 0.5), 2.0)  # under half the width: nearly uniform, exp(-d^2/8) on [0, 0.5]
    x = np.linspace(0, 0.5, 10_001)
    density = np.exp(-(x**2) / 8)
    assert float(near.mean()) == pytest.approx(float(np.trapezoid(x * density, x) / np.trapezoid(density, x)), rel=0.01)
    assert bool(np.all(near <= 0.5)) and float(em._deficit(rng, np.zeros(3), 2.0).max()) == 0.0


def test_one_emitted_electron_carries_the_components_spectra() -> None:
    fp = em.FP_COPPER
    rng = np.random.default_rng(11)
    count, e0 = 60_000, 300.0
    out = fp.emit(rng, np.full(count, e0), np.ones(count))
    assert np.allclose(out.yields, float(fp(e0)))  # each impact's weight is the yield
    de, dr, dts = float(fp.elastic(e0)), float(fp.rediffused(e0)), float(fp.true(e0))
    for code, part in ((em.BACKSCATTERED, de), (em.REDIFFUSED, dr), (em.TRUE, dts)):
        share = float(np.mean(out.component == code))
        assert share == pytest.approx(part / (de + dr + dts), abs=4 * math.sqrt(0.25 / count))
    assert bool(np.all((out.energy >= 0) & (out.energy <= e0)))
    back = out.energy[out.component == em.BACKSCATTERED]
    assert float(back.mean()) == pytest.approx(e0 - 2.0 * math.sqrt(2 / math.pi), abs=0.1)
    red = out.energy[out.component == em.REDIFFUSED]
    assert float(red.mean()) == pytest.approx(e0 * 1.5 / 2.5, rel=0.03)  # (q + 1) / (q + 2) of the impact
    # one true secondary: n - 1 binomial over M - 1 (the size-biased multiplicity), then p_n eps_n on average
    rest = 1 - de - dr
    p = dts / (rest * em.MULTIPLICITY)
    n = np.arange(1, em.MULTIPLICITY + 1)
    weights = np.array([math.comb(em.MULTIPLICITY - 1, k - 1) * p ** (k - 1) * (1 - p) ** (em.MULTIPLICITY - k) for k in n])
    expected = float(weights @ (np.array(fp.p_n) * np.array(fp.eps_n)))
    true = out.energy[out.component == em.TRUE]
    assert float(true.mean()) == pytest.approx(expected, rel=0.03)
    assert float(np.median(true)) < 6.0  # a few eV: the true secondaries' spectrum
    # Lambert's cos theta: the mean cosine 2/3, the mean square 1/2; the azimuth uniform
    assert float(out.cos.mean()) == pytest.approx(2 / 3, abs=0.005)
    assert float((out.cos**2).mean()) == pytest.approx(0.5, abs=0.005)
    assert float(np.cos(out.azimuth).mean()) == pytest.approx(0.0, abs=0.01)


def test_true_secondaries_keep_under_the_impacts_energy() -> None:
    fp = em.FP_COPPER
    rng = np.random.default_rng(5)
    for e0 in (0.05, 2.0, 8.0, 40.0):
        energies = fp.true_energy(rng, np.full(20_000, e0), np.ones(20_000))
        assert bool(np.all((energies >= 0) & (energies <= e0)))
    out = fp.emit(rng, np.zeros(1000), np.ones(1000))  # at zero energy only backscattering, with nothing to give
    assert bool(np.all(out.component == em.BACKSCATTERED))
    assert float(out.energy.max()) == 0.0


def test_the_charts_surface_is_multipactors_model() -> None:
    rng = np.random.default_rng(2)
    energies = np.linspace(0.0, 600.0, 20_000)
    out = em.CHART.emit(rng, energies, np.full(energies.size, 0.3))
    assert np.array_equal(out.yields, mp.COPPER(energies))
    assert bool(np.all(out.cos == 1.0))
    assert float(out.energy.mean()) == pytest.approx(3.0, rel=0.03)  # Anza et al.: exponential, mean 2 W_t


# ---- the chamber ---------------------------------------------------------------------------------------------------
@settings(max_examples=120, deadline=None)
@given(
    st.floats(-1.0, 1.0),
    st.floats(0.01, 1.0),
    st.floats(-1.0, 1.0),
    st.floats(0.01, 1.0),
    st.floats(0.05, 0.95),
    st.floats(0.05, 0.95),
    st.floats(0.0, 2 * math.pi),
    st.floats(1.0, 3.0),
)
def test_a_segment_meets_a_box_where_it_first_enters_it(
    u0: float, du: float, v0: float, dv: float, fu: float, fv: float, angle: float, stretch: float
) -> None:
    box = em.Box(u0, u0 + du, v0, v0 + dv)
    grid = gs.Grid()  # unused: a chamber needs a grid, the strike test does not read it
    del grid
    chamber = em.Chamber.plates(1e-3, 1e9)
    chamber = em.Chamber(chamber.grid, 1e9, (box,), (), em.Box(-10.0, 10.0, -10.0, 10.0))
    target = np.array([box.u0 + fu * du, box.v0 + fv * dv])
    start = target + 3.0 * np.array([math.cos(angle), math.sin(angle)])
    step = (target - start) * stretch
    strike = chamber.hit(start[None, :], step[None, :])
    assert int(strike.box[0]) == 0
    tau = float(strike.tau[0])

    def inside(t: float) -> bool:
        p = start + t * step
        return bool(box.u0 <= p[0] <= box.u1 and box.v0 <= p[1] <= box.v1)

    lo, hi = 0.0, 1.0 / stretch
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if inside(mid) else (mid, hi)
    assert tau == pytest.approx(hi, abs=1e-9)
    p = start + tau * step
    normal = strike.normal[0]
    on_face = abs(p[0] - (box.u0 if normal[0] < 0 else box.u1)) < 1e-9 if normal[0] else abs(
        p[1] - (box.v0 if normal[1] < 0 else box.v1)
    ) < 1e-9
    assert on_face
    assert float(normal @ step) < 0  # outward: against the motion
    # the same segment stopped short of the box misses it; one running away from it misses it
    short = chamber.hit(start[None, :], 0.999 * tau * step[None, :])
    away = chamber.hit(start[None, :], -step[None, :])
    assert int(short.box[0]) == -1 and int(away.box[0]) == -1
    assert math.isinf(float(away.tau[0]))


def test_parallel_plates_seed_their_inner_faces() -> None:
    chamber = em.plates(0.7)
    assert chamber.f == pytest.approx(0.7e9)
    points, normals, lengths = chamber.faces()
    lower, upper = normals[:, 1] > 0, normals[:, 1] < 0
    assert bool(np.all(lower | upper)) and lower.sum() == upper.sum() == 64
    assert np.allclose(points[lower, 1], em.EPS) and np.allclose(points[upper, 1], 1e-3 - em.EPS)
    assert float(lengths.sum()) == pytest.approx(4.0)
    gu, gv, on = chamber.grid.gradient_at(np.array([0.3, -0.7]), np.array([0.2e-3, 0.9e-3]))
    assert bool(on.all()) and np.allclose(gu, 0.0) and np.allclose(gv, 1e3)
    with pytest.raises(ValueError, match="band"):
        em.Chamber(chamber.grid, 0.0, chamber.copper, (), chamber.region)


def test_the_chambers_of_the_sections() -> None:
    b = em.Chamber.of(gs.solve(gs.section_b(grid=gs.COARSE)))
    h = 0.915e-3
    assert len(b.copper) == 4 and len(b.insulator) == 4
    edges = sorted({round(abs(x) / h, 6) for box in b.copper for x in (box.u0, box.u1)})
    assert edges == [0.15, 4.0]  # the split's edges and the strips' outer edges, 8h wide split 0.3h
    faces = sorted({round(v, 9) for box in b.copper for v in (box.v0, box.v1)})
    assert faces == [-18e-6, 0.0, 2 * h, 2 * h + 18e-6]
    assert b.region == em.Box(-6 * h, 6 * h, -18e-6, 2 * h + 18e-6)
    points, normals, _ = b.faces()
    assert not bool(b.solid(points).any()) and bool(b.inside(points).all())
    # no seed on a copper face against the polyimide: the strips' faces look across the gap and into the split
    toward = points[normals[:, 1] != 0]
    assert bool(np.all((np.abs(toward[:, 1] - 0.0) < 1e-8) | (np.abs(toward[:, 1] - 2 * h) < 1e-8)))
    null = np.array([[0.0, h]])
    assert bool(b.inside(null)[0]) and not bool(b.solid(null)[0])
    # R1: the strip, the rails and the outer ground, the FR-4 open in the isolation gaps
    r1 = em.Chamber.of(gs.solve(em.r1_section(gs.COARSE)))
    assert len(r1.copper) == 6 and len(r1.insulator) == 2
    gaps = np.array([[0.2e-3, -30e-6], [1.66e-3, -30e-6], [-0.2e-3, -30e-6]])
    assert not bool(r1.solid(gaps + [0.0, 30e-6 + 1e-7]).any())  # open above the FR-4 in each gap
    assert bool(r1.solid(gaps).all())
    assert r1.region.u1 == pytest.approx(2 * (0.2e-3 + 1.46e-3 + 0.075e-3) / (r1.region.u1 / 6.0) * r1.region.u1 / 6.0)


def test_r1_section_is_the_design_files() -> None:
    from array16.twin import Design

    assert em.r1_section() == Design.load(ROOT / "designs" / "r1_rail_ring.toml").section
    assert em.cases()[0].section == Design.load(ROOT / "designs" / "t915s.toml").section


# ---- flights -------------------------------------------------------------------------------------------------------
def _phase_flights(steps: int) -> tuple[Arr, Arr, Arr, Arr]:
    """Electrons leaving the lower plate at 1 GHz mm, 60 V, with 2 eV along the normal at eight drive phases: the
    tracker's first impacts against multipactor's closed form (whose drive is V sin(wt + phi) from leaving)."""
    fd, volts, w0 = 1.0, 60.0, 2.0
    chamber = em.plates(fd)
    k = mp.volts_scale(fd)
    phases = np.array([0.1, 0.4, 0.8, 1.2, 1.6, 2.4, 3.5, 5.0])
    fl = mp.flights(phases, math.sqrt(2 * w0 / k), volts / k)
    w = 2 * math.pi * chamber.f
    t0 = np.mod((phases - math.pi / 2) / w, 1 / chamber.f)  # cos(w t) = sin(w t + pi/2)
    start = np.column_stack([np.zeros(8), np.full(8, em.EPS)])
    velocity = np.column_stack([np.zeros(8), np.full(8, math.sqrt(2 * w0 * em.QM)), np.zeros(8)])
    imp = em.impacts(chamber, volts, start, velocity, t0, steps=steps, periods=10)
    assert bool(imp.struck.all())
    assert list(imp.point[:, 1] > 0.5e-3) == list(fl.far)
    return (imp.time - t0) * w, fl.tau, imp.energy, 0.5 * k * fl.speed**2


def test_impacts_meet_the_closed_form_flights_at_second_order() -> None:
    errors = []
    for steps in (200, 400):
        tau, exact_tau, energy, exact_energy = _phase_flights(steps)
        errors.append((float(np.max(np.abs(tau - exact_tau))), float(np.max(np.abs(energy / exact_energy - 1)))))
    assert errors[0][0] < 1e-3 and errors[0][1] < 2e-3
    assert errors[0][0] / errors[1][0] > 3.0 and errors[0][1] / errors[1][1] > 3.0  # twice the steps, four times less


# ---- the tracker in the parallel-plate limit -----------------------------------------------------------------------
def test_the_tracker_is_multipactors_model_between_plates() -> None:
    # Normal on the chart's yield through multipactor's gain is the one-dimensional model: the growth per RF period
    # is the operator's Euler-Lotka root, above and below copper's lowest point
    chamber = em.plates(0.70)
    for scale in (1.15, 0.85):
        volts = scale * mp.CALIBRATION[1]
        expected = mp.growth_rate(0.70, volts)
        found = em.growth(chamber, volts, em.CHART, gain=mp.GAIN, seeds=(0, 1), **QUICK)
        assert (found.rate > 0) == (scale > 1)
        assert found.rate == pytest.approx(expected, abs=4 * found.spread + 0.012)
        assert found.spread == pytest.approx(float(np.std([r.growth for r in found.runs], ddof=1)))
        run = found.runs[0]
        assert run.held == 0.0 and run.escaped == 0.0 and math.isnan(run.laminate) is False
        assert run.impacts > 1.5


def test_a_gap_without_volts_decays() -> None:
    # no field: secondaries cross on their own few eV, where copper frees about half an electron
    found = em.track(em.plates(1.0), 1e-9, em.FP_COPPER, **QUICK)
    assert found.growth < -0.3
    assert found.energy < 12.0


def test_bare_copper_multipacts_where_conditioned_copper_does_not() -> None:
    # the same field at 1.0 GHz mm and 57.6 V, 1.6 times ECSS's chart: bare copper's population grows; conditioned
    # copper, its true secondaries scaled to a peak of 1.2, decays (the case that must fail with the yield)
    chamber = em.plates(1.0)
    bare = em.track(chamber, 57.6, em.FP_COPPER, **QUICK)
    cond = em.track(chamber, 57.6, em.FP_CONDITIONED, **QUICK)
    assert bare.growth > 0.2
    assert cond.growth < -0.2
    with pytest.raises(ValueError, match="settle"):
        em.track(chamber, 57.6, periods=5, settle=5)
