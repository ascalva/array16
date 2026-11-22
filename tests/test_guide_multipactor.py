"""The multipactor instrument (array16/guide/multipactor.py), held to closed forms, to the core's tracker, to the
knotted study it was taken from, and to ECSS-E-ST-20-01C's copper chart, itself the standard's model on its yield:
MPC-1 to MPC-3 of docs/KNOTTED_GUIDE.md section 11, the consistency check at the held-out point and DR-0's blind check
on the twelve points registered before the model was run on them."""

from __future__ import annotations

import math

import knotted_guide_study as study
import numpy as np
import pytest
from hypothesis import assume, given, settings
from hypothesis import strategies as st

from array16 import particles as pt
from array16.guide import multipactor as mp
from array16.particles import Arr, BoolArr

FAST = 96, 16  # a coarser operator (phases, nodes), for searches whose answer is far from the line drawn


# ---- the yield and the emission ------------------------------------------------------------------------------------
def test_copper_yield_is_vaughans_curve_on_the_standards_parameters() -> None:
    y = mp.COPPER
    assert (y.e1, y.dmax, y.emax, y.d0) == (19.0, 2.48, 232.0, 0.8)  # ECSS-E-ST-20-01C Table 9-1
    assert float(y(19.0)) == pytest.approx(1.0, abs=1e-12)
    assert float(y(232.0)) == pytest.approx(2.48, abs=1e-12)
    slope = (float(y.curve(232.0 + 1e-3)) - float(y.curve(232.0 - 1e-3))) / 2e-3
    assert abs(slope) < 1e-6  # the peak
    assert y.e0 == pytest.approx(0.828, abs=1e-3)
    assert float(y.curve(5.0)) == pytest.approx(0.4538, abs=1e-4)  # Vaughan's curve alone
    assert [float(y(e)) for e in (0.0, 5.0, 12.0)] == [0.8, 0.8, 0.8]  # held up to the zero-energy value
    assert float(y(13.0)) > 0.8
    w = y.e0 + 3.6 * (y.emax - y.e0)  # where the 1993 tail takes over: continuous to Vaughan's rounding
    assert float(y.curve(w * (1 + 1e-9))) == pytest.approx(float(y.curve(w * (1 - 1e-9))), rel=1e-3)
    assert y.e2 == pytest.approx(4337, abs=1)
    assert float(y(y.e2)) == pytest.approx(1.0, abs=1e-9)


def test_copper_yield_is_the_studys_below_the_tail() -> None:
    # extracted, not rewritten: the knotted study's Yield, up to where Vaughan's 1993 tail takes over (833 eV)
    energies = np.linspace(0.0, 830.0, 4151)
    assert np.max(np.abs(mp.COPPER(energies) - study.Yield()(energies))) < 1e-9


@settings(max_examples=60, deadline=None)
@given(st.floats(0.0, 231.0), st.floats(0.0, 231.0), st.floats(233.0, 5e4), st.floats(233.0, 5e4))
def test_the_curve_rises_to_its_peak_and_falls_after(a: float, b: float, c: float, d: float) -> None:
    y = mp.COPPER
    lo, hi = sorted((a, b))
    assert float(y.curve(lo)) <= float(y.curve(hi)) + 1e-12
    assert float(y(lo)) >= y.d0
    lo, hi = sorted((c, d))
    assert float(y.curve(lo)) >= float(y.curve(hi)) - 1e-12


def test_emission_is_exponential_in_normal_energy() -> None:
    em = mp.Emission(1.5)  # Anza et al.: W_t = 1.5 eV, so a mean normal energy of 3 eV
    nodes = em.nodes(4000)
    assert float(nodes.mean()) == pytest.approx(3.0, rel=2e-3)
    assert float(np.median(nodes)) == pytest.approx(3.0 * math.log(2), rel=1e-3)
    draws = em.draw(np.random.default_rng(1), 200_000)
    assert float(draws.mean()) == pytest.approx(3.0, rel=0.01)
    assert float(np.mean(draws > 6.0)) == pytest.approx(math.exp(-2), rel=0.03)


# ---- flights, against closed forms and against the core's integrator ----------------------------------------------
@settings(max_examples=60, deadline=None)
@given(st.floats(0.3, 4.0), st.floats(0.0, 5.0), st.floats(0.0, 1.4))
def test_first_order_flights_land_where_the_closed_form_puts_them(fd: float, w0: float, phase: float) -> None:
    # order 1 from a phase in [0, pi/2): the field never pushes the electron back before it lands, at tau = pi
    k = mp.volts_scale(fd)
    beta = math.sqrt(2 * w0 / k)
    assume(1 - math.pi * beta > 0.05)
    volts = float(mp.resonant_volts(fd, 1, phase, w0))
    fl = mp.flights(phase, beta, volts / k)
    assert bool(fl.hit[0])
    assert bool(fl.far[0])
    assert float(fl.tau[0]) == pytest.approx(math.pi, abs=1e-9)
    energy = 0.5 * k * float(fl.speed[0]) ** 2
    assert energy == pytest.approx(float(mp.resonant_impact(fd, 1, phase, w0)), rel=1e-9)
    assert float(mp.impact_energy(fd, volts, phase, w0)) == pytest.approx(energy, rel=1e-12)


class _Plates:
    """Stops electrons that leave the gap 0 < x < d, keeping the step, the speed and the plate of each."""

    def __init__(self, gap: float, count: int) -> None:
        self.gap = gap
        self.step = np.full(count, -1)
        self.speed = np.zeros(count)
        self.far = np.zeros(count, dtype=np.bool_)

    def __call__(self, step: int, t: float, x: Arr, u: Arr, alive: BoolArr) -> BoolArr:
        far = x[:, 0] >= self.gap
        out: BoolArr = alive & (far | (x[:, 0] <= 0.0))
        self.step[out] = step + 1
        self.speed[out] = np.abs(u[out, 0])
        self.far[out] = far[out]
        return out


def test_flights_agree_with_the_core_integrator() -> None:
    # INT: Boris steps of pt.integrate in the same uniform field converge on the analytic flights, at first order
    # (the core samples the field at the start of each step); each electron's field has its own phase, set by the
    # y at which it travels
    f, d, volts, w0 = 1e9, 1e-3, 60.0, 2.0
    w = 2 * math.pi * f
    phases = np.array([0.1, 0.4, 0.8, 1.2, 1.6, 2.4, 3.5, 5.0])
    n = phases.size

    def fields(x: Arr, t: float) -> tuple[Arr | None, Arr | None]:
        e = np.zeros_like(x)
        e[:, 0] = -(volts / d) * np.sin(w * t + phases[np.rint(x[:, 1]).astype(int)])
        return e, None

    v0 = math.sqrt(2 * w0 * pt.E / pt.ME)
    x0 = np.column_stack([np.full(n, 1e-12), np.arange(n, dtype=float), np.zeros(n)])
    u0 = np.column_stack([np.full(n, v0), np.zeros(n), np.zeros(n)])
    fl = mp.flights(phases, v0 / (w * d), volts / mp.volts_scale(f * d / 1e6))
    assert bool(np.all(fl.hit))
    assert bool(np.any(~fl.far))  # some come back to their own plate
    errors = []
    for per_period in (2000, 8000):
        dt = 1 / f / per_period
        plates = _Plates(d, n)
        pt.integrate(x0, u0, np.full(n, pt.QM_E), dt, 4 * per_period, fields, [plates])
        assert bool(np.all(plates.step > 0))
        assert list(plates.far) == list(fl.far)
        error = float(np.max(np.abs(plates.step * dt * w - fl.tau)))
        assert error < 5 * dt * w  # within five steps: most where an electron turns back slowly to its own plate
        errors.append(error)
        kick = 2 * mp.QM * volts / d * dt  # what the field adds to a speed within a step, twice
        assert np.allclose(plates.speed, fl.speed * w * d, rtol=2e-3, atol=kick)
        assert bool(np.all(np.abs(plates.speed - v0) <= mp.speed_change_bound(volts / d, f) * (1 + 1e-9)))  # MPC-1
    assert errors[0] / errors[1] > 2.5  # a step four times shorter, an error near four times smaller


# ---- MPC-1: impacts carry at most (sqrt W0 + sqrt(8 D))^2 ----------------------------------------------------------
@settings(max_examples=80, deadline=None)
@given(st.floats(0.2, 10.0), st.floats(0.01, 2.0), st.floats(0.0, 2 * math.pi), st.floats(0.0, 20.0))
def test_mpc_1_impacts_never_pass_the_bound(fd: float, v: float, phase: float, w0: float) -> None:
    k = mp.volts_scale(fd)
    volts, beta = v * k, math.sqrt(2 * w0 / k)
    fl = mp.flights(phase, beta, v)
    assume(bool(fl.hit[0]))
    # the speed changes by at most 2 v (in units of w d): 2 e E0 / (m w) with E0 = V / d
    assert abs(float(fl.speed[0]) - beta) <= 2 * v * (1 + 1e-9) + 1e-12
    energy = 0.5 * k * float(fl.speed[0]) ** 2
    assert energy <= mp.impact_bound(w0, mp.gap_depth(fd, volts)) * (1 + 1e-9) + 1e-12
    # the same bound in SI: 2 e E0 / (m w) is 2 v w d
    f, d = 1e9, fd * 1e-3
    assert mp.speed_change_bound(volts / d, f) == pytest.approx(2 * v * 2 * math.pi * f * d, rel=1e-12)


def test_mpc_1_the_bound_is_reached() -> None:
    # an electron that leaves at the field's zero with the field rising, and lands half a period later, gains the
    # whole 2 e E0 / (m w): at w0 = 0 its impact energy is 8 D
    fd = 1.0
    volts = float(mp.resonant_volts(fd, 1, 0.0))
    energy = float(mp.impact_energy(fd, volts, 0.0))
    assert energy == pytest.approx(mp.impact_bound(0.0, mp.gap_depth(fd, volts)), rel=1e-9)


# ---- MPC-2: V = 35.7 q (fd)^2 and D = 1.42 q^2 (fd)^2 on facing boards ---------------------------------------------
@settings(max_examples=80, deadline=None)
@given(st.floats(0.05, 0.9), st.floats(0.01, 10.0), st.floats(0.1e-3, 5e-3))
def test_mpc_2_agrees_with_the_core_drive(q: float, depth: float, h: float) -> None:
    fb = pt.FacingBoards(2 * h)
    drive = pt.electrode_drive_at_q(depth, q, 1.0, h, fb.depth_share, fb.gradient_share)
    fd = drive.f * 2 * h / 1e6
    assert mp.facing_volts(q, fd) == pytest.approx(drive.volts, rel=1e-9)
    assert mp.facing_depth(q, fd) == pytest.approx(depth, rel=1e-9)
    assert mp.facing_fd(depth, q) == pytest.approx(fd, rel=1e-9)


def test_mpc_2_coefficients_and_the_one_ev_well() -> None:
    assert mp.facing_volts(1.0, 1.0) == pytest.approx(35.72, abs=0.005)
    assert mp.facing_depth(1.0, 1.0) == pytest.approx(1.4214, abs=5e-5)
    fd = mp.facing_fd(1.0, 0.5)
    assert (fd, mp.facing_volts(0.5, fd)) == pytest.approx((1.678, 50.27), abs=0.005)
    # the chart point follows from the depth and q alone: a 1 eV well at q = 0.5 keeps 5.2 dB in every band, its
    # gap scaled as 1 / f
    for f in (433e6, 915e6, 2450e6, 5800e6):
        gap = fd * 1e6 / f
        reading = mp.level1(mp.Gap("1 eV at q = 0.5", f, mp.facing_volts(0.5, fd), gap, gap))
        assert reading.fd == pytest.approx(fd)
        assert reading.margin == pytest.approx(5.22, abs=0.01)


# ---- the level-1 chart, classical ----------------------------------------------------------------------------------
def test_first_order_band_is_the_textbook_one() -> None:
    # docs/KNOTTED_GUIDE.md section 5.5: first order across g spans (wg)^2 / (|e/m| sqrt(pi^2 + 4)) to
    # (wg)^2 / (|e/m| pi), 42 to 50 V for g = h at 915 MHz (h = 0.917 mm) and 3.8 to 4.5 V for 0.3h
    anything = lambda e: np.full_like(e, 2.0)  # noqa: E731
    for g, (lo, hi) in ((0.917e-3, (42.4, 50.3)), (0.275e-3, (3.82, 4.53))):
        fd = 915e6 * g / 1e6
        k = mp.volts_scale(fd)
        band = mp.susceptibility(fd, anything, orders=(1,))[0]
        assert (band.low, band.high) == pytest.approx((k / math.sqrt(math.pi**2 + 4), k / math.pi), rel=5e-4)
        assert (band.low, band.high) == pytest.approx((lo, hi), abs=0.05)
    # copper cuts the band where impacts carry less than its first crossover
    band = mp.susceptibility(0.839)[0]
    assert band.energies[0] == pytest.approx(19.0, abs=0.1)
    assert band.low > mp.volts_scale(0.839) / math.sqrt(math.pi**2 + 4)


def test_phase_focusing_edges() -> None:
    # at w0 = 0, G = (n pi / 2) tan(phi): the resonance holds its phase for 0 <= tan(phi) < 2 / (n pi)
    for n in (1, 3, 5):
        edge = math.atan(2 / (n * math.pi))
        assert float(mp.focusing(1.0, n, edge)) == pytest.approx(1.0, abs=1e-12)
        assert float(mp.focusing(1.0, n, 0.3 * edge)) == pytest.approx(n * math.pi / 2 * math.tan(0.3 * edge))
    with pytest.raises(ValueError, match="odd orders"):
        mp.susceptibility(1.0, orders=(2,))


def test_classical_nose_is_pi_e1_over_two() -> None:
    # with no emission speed, order one first multiplies at phase 0 where 2K / pi^2 = E1: V = pi E1 / 2 at
    # K = pi^2 E1 / 2
    fds = np.linspace(0.63, 0.67, 81)
    low, _ = mp.chart(fds, orders=(1,))
    i = int(np.nanargmin(low))
    assert low[i] == pytest.approx(math.pi * 19.0 / 2, rel=3e-3)
    assert fds[i] == pytest.approx(math.sqrt(math.pi**2 * 19.0 / 2 / mp.volts_scale(1.0)), abs=1e-3)


# ---- the statistical theory: calibration and DR-0 ------------------------------------------------------------------
def test_ecss_table_is_the_studys_and_the_standards() -> None:
    assert mp.ECSS_COPPER == study.ECSS_CU
    table = np.array(mp.ECSS_COPPER)
    i = int(np.argmin(table[:, 1]))
    assert tuple(table[i]) == mp.CALIBRATION == (0.70, 27.4)
    assert tuple(table[list(table[:, 0]).index(1.68)]) == mp.HELD_OUT
    assert float(mp.ecss_copper(1.83 * 0.915)) == pytest.approx(91.37, abs=0.01)
    assert bool(np.isnan(mp.ecss_copper(0.3)))
    assert bool(np.isnan(mp.ecss_copper(11.0)))


def test_the_gain_is_calibrated_on_coppers_lowest_point() -> None:
    gain = mp.calibrate()
    assert gain == pytest.approx(mp.GAIN, abs=2e-4)
    assert mp.threshold(mp.CALIBRATION[0]) == pytest.approx(27.4, rel=2e-3)
    # the published inputs alone (gain 1) put copper's lowest point 2.0 % low
    assert mp.threshold(0.70, gain=1.0) == pytest.approx(26.84, abs=0.03)


def test_the_held_out_point_is_consistent() -> None:
    # 91.9 V at 1.68 GHz mm, the published point nearest T915-S's 1.674, computed while the model was built: the
    # calibrated instrument gives 91.2 V, two models of one yield agreeing within 10 %
    fd, published = mp.HELD_OUT
    found = mp.threshold(fd)
    assert abs(found / published - 1) < 0.10
    assert found == pytest.approx(91.2, abs=0.3)


def test_dr_0_the_blind_points_of_the_copper_chart() -> None:
    # DR-0's blind check: the twelve points of Table 5-1's copper column from 0.49 to 1.51 GHz mm that ECSS_COPPER
    # leaves out, registered before the model was run on any of them, each within 10 % of the calibrated instrument
    table = dict(mp.ECSS_COPPER)
    assert len(mp.BLIND) == 12
    assert not set(dict(mp.BLIND)) & set(table)
    assert all(0.47 < fd < 1.59 for fd, _ in mp.BLIND)
    for fd, published in mp.BLIND:
        assert abs(mp.threshold(fd) / published - 1) < mp.BLIND_BAND, f"{fd} GHz mm"


def test_the_published_inputs_follow_the_whole_column() -> None:
    # the model taken as published (gain 1) on points of every stretch of the column: the left branch, the
    # first-order branch, the order-3 dip and the high orders, each within the column's worst 5.4 %
    for fd, published in ((0.50, 32.5), (1.28, 54.2), (2.92, 149.2), (4.30, 211.2)):
        assert mp.threshold(fd, gain=1.0) == pytest.approx(published, rel=0.054)


def test_calibrations_that_must_fail() -> None:
    fd, volts = mp.CALIBRATION
    held, published = mp.HELD_OUT
    # the emission's energy calibrated on copper's lowest point instead of the gain: the minimum hardly moves
    # with it, so it lands far off and misses the held-out point by more than DR-0 allows
    lo, hi = 0.2, 1.5
    for _ in range(10):
        mid = 0.5 * (lo + hi)
        at = mp.threshold(fd, emission=mp.Emission(mid), gain=1.0, rtol=1e-3, phases=FAST[0], nodes=FAST[1])
        lo, hi = (mid, hi) if at > volts else (lo, mid)
    emission = mp.Emission(0.5 * (lo + hi))
    assert emission.thermal == pytest.approx(0.65, abs=0.05)
    assert mp.threshold(held, emission=emission, gain=1.0, phases=FAST[0], nodes=FAST[1]) / published - 1 > 0.10
    # the classical chart, its one emission energy calibrated on the same point
    lo, hi = 0.0, 1.0
    for _ in range(20):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if float(mp.chart([fd], w0=mid)[0][0]) > volts else (lo, mid)
    w0 = 0.5 * (lo + hi)
    assert float(mp.chart([fd], w0=w0)[0][0]) == pytest.approx(volts, abs=0.01)
    assert float(mp.chart([held], w0=w0)[0][0]) / published - 1 > 0.5


def test_the_tracker_agrees_with_the_operator() -> None:
    # the 1-D population tracker's growth per RF period against the operator's Euler-Lotka root, above and below
    # the threshold; and the operator's radius crosses one where its growth does
    for fd, scale in ((0.70, 1.15), (1.68, 0.9)):
        volts = scale * mp.threshold(fd)
        expected = mp.growth_rate(fd, volts)
        pop = mp.track(fd, volts, electrons=1500, periods=36)
        assert (expected > 0) == (scale > 1)
        assert math.isfinite(pop.growth)
        assert pop.spread < 0.004
        assert pop.growth == pytest.approx(expected, abs=4 * pop.spread + 0.005)
        assert 0 < pop.far < 1
        assert pop.impacts > 0
    at = mp.threshold(1.0)
    assert float(mp.effective_yield(1.0, at * 0.995)[0]) < 1 < float(mp.effective_yield(1.0, at * 1.005)[0])
    assert abs(mp.growth_rate(1.0, at)) < 0.01


def test_a_gap_without_volts_loses_its_electrons() -> None:
    # no field: secondaries cross on their emission speed alone, at a yield of d0 < 1 every crossing
    pop = mp.track(1.0, 1e-9, electrons=500, periods=20)
    assert pop.growth < 0
    assert pop.energy < 10
    assert mp.growth_rate(1.0, 1e-9) < 0


# ---- MPC-3: level-1 margins against ECSS copper --------------------------------------------------------------------
def test_mpc_3_level_1_margins() -> None:
    read = {gap.name: mp.level1(gap) for gap in mp.CASES}
    rods3, rods6 = read["first run's rods, q = 0.3"], read["first run's rods, q = 0.6"]
    # docs/PRIOR_ART.md: 6.5 dB at 13 V and 0.8 dB at 25 V, the rods' arcs reaching copper's lowest point
    assert (rods3.fd, rods3.threshold) == (0.70, 27.4)
    assert (rods3.margin, rods6.margin) == pytest.approx((6.48, 0.80), abs=0.01)
    # T915-S's 50 V across 1.83 mm at 915 MHz: 5.2 dB (docs/KNOTTED_GUIDE.md section 5.5)
    assert read["T915-S, 1 eV at q = 0.5"].margin == pytest.approx(5.24, abs=0.01)
    # R1's rails at 50 V: arcs of every length from the slit to the rails' width, copper's lowest point among them
    r1 = read["R1's rails at 50 V, the parallel-plate bound"]
    assert r1.fd == 0.70
    assert r1.margin == pytest.approx(-5.22, abs=0.01)
    # facing boards at q = 0.5: 2.0 dB at a 2.5 eV well, none left at 3 eV
    assert read["facing boards, 2.5 eV at q = 0.5"].margin == pytest.approx(2.04, abs=0.01)
    assert read["facing boards, 3 eV at q = 0.5"].margin == pytest.approx(-0.06, abs=0.01)


def test_mpc_3_the_instrument_reads_the_designs_alike() -> None:
    # the calibrated instrument at T915-S's own f d, and the growth it gives the designs' gaps at their volts
    t915 = next(g for g in mp.CASES if g.name.startswith("T915-S"))
    fd = t915.span[0]
    found = mp.level1(t915, mp.threshold)
    assert found.margin == pytest.approx(5.17, abs=0.05)
    assert mp.growth_rate(fd, t915.volts) < 0
    assert mp.growth_rate(0.70, 50.0) > 0.4  # R1's rails at 50 V on their 0.70 GHz mm arcs multiply
    assert mp.growth_rate(0.70, 25.0) < 0  # the rods at 25 V, under copper's lowest point


def test_gap_validation() -> None:
    assert math.isnan(mp.threshold(0.0))
    with pytest.raises(ValueError, match="positive f d"):
        mp.growth_rate(0.0, 10.0)
    with pytest.raises(ValueError, match="settle"):
        mp.track(1.0, 10.0, periods=5, settle=5)
    with pytest.raises(ValueError, match="dmax"):
        mp.Yield(dmax=0.9)
    with pytest.raises(ValueError, match="thermal"):
        mp.Emission(0.0)
