"""The ring resonator (array16/guide/ring.py), held to the theory of docs/KNOTTED_GUIDE.md sections 5.2 and 5.3, to
the knotted study's functions on the same seeds, and to the numbers the study and the buildability lens recorded."""

from __future__ import annotations

import math

import knotted_guide_study as study
import numpy as np
import pytest
from hypothesis import assume, given, settings
from hypothesis import strategies as st
from hypothesis.extra import numpy as hnp

from array16.guide import ring
from array16.particles import C0, Arr

# The line study's sections (experiments/knotted_guide/results.json): section B at 915 MHz on bare copper, Q by
# Wheeler's rule on its 18 um; section A on slotted RO4003C, its Q and Z0 on 35 um copper.
Q_BARE = 213.81496523444684
Q_A, Z0_A = 250.78042004090963, 73.51321999751003
L915 = C0 / (ring.N915 * ring.F915)  # T915-S: m = 1 at 915 MHz, 317.8 mm

LENS5 = (  # lens5/ring.py's evaluate on ring915.py's draws (seed 7, Q 381, 2 %): 2Q|c_2m|, the hybrid's |b/a| and
    # depth per watt, the re-solved feeds' |b/a| and depth per watt, the second feed's amplitude and phase
    (1.1506357521841937, 0.18386959811125966, 0.6416204123378344, 0.010289723479643607, 0.9567970854499657, 0.625,
     -92.43243243243244),
    (0.7521990026695583, 0.12458477498807383, 0.7496002751500863, 0.00851173697262015, 0.948048373509738, 1.375,
     82.70270270270268),
    (1.5871757061862706, 0.2607672135196098, 0.51062172628541, 0.0016668356801840837, 0.9443285844754605, 1.875,
     -180.0),
)  # fmt: skip

SINGULAR = {  # results.json: G's singular values (V/A) for six feeds on section A's ring, where the largest peaks
    "915 MHz, m = 1": (
        (3642.457186376325, 3642.449031033232, 28.40119891236767, 9.734635991642099, 9.734597807977972,
         5.126736711868821),
        877.5490548672465,
    ),
    "the same, feed susceptance compensated": (
        (3631.5165180737313, 3631.507652217766, 29.31509052525608, 10.01891145770291, 10.018868404326811,
         5.258936364787922),
        911.8890000000001,
    ),
}  # fmt: skip


@st.composite
def profiles(draw: st.DrawFn, sizes: tuple[int, ...] = (16, 32, 48)) -> tuple[Arr, Arr]:
    """A ring of segments: Z within 10 % of 50 ohm and the index between 1 and 1.2, segment by segment."""
    nodes = draw(st.sampled_from(sizes))
    lnz = draw(hnp.arrays(np.float64, nodes, elements=st.floats(-0.1, 0.1)))
    index = draw(hnp.arrays(np.float64, nodes, elements=st.floats(1.0, 1.2)))
    return 50.0 * np.exp(lnz), index


def test_a_perfect_ring_travels_under_quadrature_feeds_and_stands_under_one() -> None:
    # RNG-1 and section 5.2 on the ring of segments: two current feeds a quarter guided wave apart at +-90 degrees
    # leave no component of order -m (|b/a| = 0) at any frequency, one feed a standing wave (|b/a| = 1); the travelling
    # wave takes mean |V|^2 / (2R) watts, R = Z0 Q / (2 pi m) (RNG-2), and ring_score's weakest depth per watt sits
    # below 2R by less than a lap's loss, the wave fed at one place losing 1 - exp(-2 pi m / Q) (RNG-4)
    z, n, rel = ring.perfect_ring(ring.N915, z0=41.7)
    for m in (1, 2):
        length = m * L915
        f0 = m * C0 / (ring.N915 * length)
        quarter = ring.RING_NODES // (4 * m)
        for f in (0.8 * f0, f0 * (1 - 2 / Q_BARE), f0, f0 * (1 + 1 / Q_BARE)):
            g = ring.ring_response(f, length, z, n, Q_BARE, (0, quarter))
            for w in (-1j, 1j):
                assert ring.ring_purity(g @ np.array([1, w]), m) < 1e-10
            assert ring.ring_purity(g[:, 0], m) == pytest.approx(1.0, abs=1e-10)
        v = ring.ring_response(f0, length, z, n, Q_BARE, (0, quarter)) @ np.array([1, -1j])
        power = 0.5 * float(np.real(v[0] + v[quarter] * 1j))
        shunt = 41.7 * Q_BARE / (2 * math.pi * m)
        assert float(np.mean(np.abs(v) ** 2)) / power == pytest.approx(2 * shunt, rel=1e-6)
        best = ring.ring_score(length, m, z, n, rel, Q_BARE, False)
        assert best.purity < 1e-10
        assert math.exp(-2 * math.pi * m / Q_BARE) < best.score / (2 * shunt) < 1.0


@settings(max_examples=100, deadline=None)
@given(
    st.integers(1, 6),
    st.integers(1, 3),
    st.floats(-3.0, 3.0),
    st.floats(20.0, 1000.0),
    st.floats(20.0, 120.0),
    st.floats(0.0, 2 * math.pi),
)
def test_equal_feeds_reach_the_backward_wave_only_when_their_count_divides_2m(
    k: int, m: int, detune: float, q0: float, z0: float, phase: float
) -> None:
    # RNG-3: K equal feeds spaced L/K and phased for the travelling wave reach -m only when K divides 2m, and then
    # drive both senses alike (|b/a| = 1 on a perfect ring); otherwise none of -m reaches it (|b/a| = 0)
    nodes = 120
    length = m * L915
    f = m * C0 / (ring.N915 * length) * (1 + detune / q0)
    z, n, _ = ring.perfect_ring(ring.N915, nodes, z0)
    weights = np.exp(1j * phase - 2j * math.pi * m * np.arange(k) / k)
    g = ring.ring_response(f, length, z, n, q0, [j * nodes // k for j in range(k)])
    assert ring.ring_purity(g @ weights, m) == pytest.approx(1.0 if (2 * m) % k == 0 else 0.0, abs=1e-9)


def test_one_fourier_component_of_the_impedance_splits_the_pair_by_twice_itself() -> None:
    # section 5.2 and DR-0: ln Z = eps cos(2 pi 2m s / L + phi) splits +-m by Delta f / f = 2|c_2m| (Little, Laine and
    # Chu) to 1 %, and to 1e-5 once the segments' sampling of the profile, sin(x)/x at x = 2 pi m / N, is counted, in
    # linewidths 2Q|c_2m|; both keep Q0, and the lossy admittance is singular at both; drawn T915-S rings split by
    # their own 2|c_2m| to 1 %
    nodes = ring.RING_NODES
    k = np.arange(nodes)
    n = np.full(nodes, ring.N915)
    for m in (1, 2):
        length = m * L915
        z = 41.7 * np.exp(1e-3 * np.cos(2 * math.pi * 2 * m * k / nodes + 0.7))
        c = ring.ring_harmonic(z, m, 41.7)
        assert c == pytest.approx(5e-4, rel=1e-9)
        pair = ring.ring_resonances(length, z, n, Q_BARE, m)
        x = 2 * math.pi * m / nodes
        assert pair.splitting == pytest.approx(2 * c, rel=0.01)
        assert pair.splitting == pytest.approx(2 * c * math.sin(x) / x, rel=1e-5)
        assert pair.linewidths == pytest.approx(2 * Q_BARE * c, rel=0.01)
        assert pair.q == pytest.approx(Q_BARE, rel=1e-12)
        assert pair.f == pytest.approx(ring.F915 / (1 + 1 / (4 * Q_BARE**2)), rel=1e-6)
        for root in (pair.lower, pair.upper):
            s = np.linalg.svd(ring.ring_admittance(root, length, z, n, Q_BARE), compute_uv=False)
            assert s[-1] < 1e-12 * s[0]
        assert ring.ring_resonances(length, np.full(nodes, 41.7), n, Q_BARE, m).splitting < 1e-12
    rng = np.random.default_rng(7)
    for _ in range(3):
        z, n, _ = ring.ring_draw(rng, L915, 0.02, ring.N915, 0.01, 0.0)
        pair = ring.ring_resonances(L915, z, n, Q_BARE)
        assert pair.linewidths == pytest.approx(2 * Q_BARE * ring.ring_harmonic(z), rel=0.01)


def test_a_weak_reflection_builds_a_backward_wave_of_the_splitting_in_linewidths() -> None:
    # section 5.2's backscatter on the ring of segments: one segment of Z0 (1 + delta) reflects S = G (1 - e^(-2i beta
    # l)) / (1 - G^2 e^(-2i beta l)), G = delta / (2 + delta); driven by the quadrature pair at the pair's centre it
    # leaves |b/a| = |S| Q / (pi m) to 1 %, the splitting in linewidths 2Q|c_2m|: Kippenberg, Spillane and Vahala's
    # coupled modes, Q = omega tau with tau the whole lifetime (the full linewidth f / Q), |b/a| the backscatter rate
    # over the half linewidth. Section 5.2's Claim writes |S| Q / (2 pi m), Q counted on the half linewidth.
    nodes = ring.RING_NODES
    for m, q0, delta in ((1, Q_BARE, 1e-3), (1, 60.0, 1e-2), (2, Q_BARE, 1e-2)):
        length = m * L915
        z, n, _ = ring.perfect_ring(ring.N915)
        z[120] *= 1 + delta
        step = 2 * math.pi * m / nodes  # beta l of one segment at the resonance
        g = delta / (2 + delta)
        reflection = abs(g * (1 - np.exp(-2j * step)) / (1 - g**2 * np.exp(-2j * step)))
        pair = ring.ring_resonances(length, z, n, q0, m)
        v = ring.ring_response(pair.f, length, z, n, q0, (0, nodes // (4 * m))) @ np.array([1, -1j])
        assert ring.ring_purity(v, m) == pytest.approx(reflection * q0 / (math.pi * m), rel=0.01)
        assert ring.ring_purity(v, m) == pytest.approx(pair.linewidths, rel=0.01)


@settings(max_examples=50, deadline=None)
@given(
    st.sampled_from([1, 2]),
    st.floats(2e-4, 2e-3),
    st.floats(0.0, 2 * math.pi),
    st.floats(30.0, 1000.0),
    st.sampled_from([3, 4, 5, 6, 8]),
)
def test_a_profile_splits_the_pair_only_through_a_harmonic_dividing_2m(
    m: int, eps: float, phi: float, q0: float, other: int
) -> None:
    # section 5.2's selection rule: a profile of period L/h couples +-m only when h divides 2m. One component at 2m
    # splits the pair by 2|c_2m| sin(x)/x (x = 2 pi m / N) to 1e-4 whatever its phase and Q0; one at an h that does not
    # divide 2m leaves the pair one frequency
    assume((2 * m) % other)
    nodes = 120
    k = np.arange(nodes)
    n = np.full(nodes, ring.N915)
    length = m * L915
    z = 50.0 * np.exp(eps * np.cos(2 * math.pi * 2 * m * k / nodes + phi))
    pair = ring.ring_resonances(length, z, n, q0, m)
    x = 2 * math.pi * m / nodes
    assert pair.splitting == pytest.approx(2 * ring.ring_harmonic(z, m) * math.sin(x) / x, rel=1e-4)
    assert pair.q == pytest.approx(q0, rel=1e-9)
    whole = ring.ring_resonances(length, 50.0 * np.exp(eps * np.cos(2 * math.pi * other * k / nodes + phi)), n, q0, m)
    assert whole.splitting < 1e-10


@settings(max_examples=100, deadline=None)
@given(profiles(), st.floats(0.5, 1.5), st.floats(10.0, 1e4), st.data())
def test_the_ring_of_segments_is_reciprocal(
    profile: tuple[Arr, Arr], scale: float, q0: float, data: st.DataObject
) -> None:
    # the nodal admittance is symmetric, element for element, and so is its inverse: the volts at one node per ampere
    # into another are the same both ways
    z, n = profile
    f = scale * ring.F915
    y = ring.ring_admittance(f, L915, z, n, q0)
    assert np.array_equal(y, y.T)
    inverse = np.linalg.inv(y)
    assert np.abs(inverse - inverse.T).max() <= 1e-9 * np.abs(inverse).max()
    a, b = (data.draw(st.integers(0, len(z) - 1)) for _ in range(2))
    g = ring.ring_response(f, L915, z, n, q0, (a, b))
    assert abs(g[b, 0] - g[a, 1]) <= 1e-9 * np.abs(g).max()


@settings(max_examples=60, deadline=None)
@given(
    profiles(),
    st.integers(0, 47),
    st.sampled_from([1, 2]),
    st.floats(-2.0, 2.0),
    st.floats(30.0, 1000.0),
    st.complex_numbers(min_magnitude=0.2, max_magnitude=5.0),
)
def test_rotating_or_mirroring_the_ring_changes_none_of_its_numbers(
    profile: tuple[Arr, Arr], shift: int, m: int, detune: float, q0: float, w: complex
) -> None:
    # a ring's numbers are the ring's, not its labels': starting it at another node rolls its node volts and keeps
    # |b/a| and the pair; running it the other way keeps the pair
    z, n = profile
    nodes = len(z)
    r = shift % nodes
    length = m * L915
    f = m * C0 / (float(np.mean(n)) * length) * (1 + detune / q0)
    feeds = (0, nodes // (4 * m))
    weights = np.array([1.0, w])
    v = ring.ring_response(f, length, z, n, q0, feeds) @ weights
    turned = ring.ring_response(f, length, np.roll(z, r), np.roll(n, r), q0, tuple((p + r) % nodes for p in feeds))
    vr = turned @ weights
    assert np.abs(vr - np.roll(v, r)).max() <= 1e-9 * np.abs(v).max()
    assert ring.ring_purity(vr, m) == pytest.approx(ring.ring_purity(v, m), abs=1e-9)
    pair = ring.ring_resonances(length, z, n, q0, m)
    for other in (
        ring.ring_resonances(length, np.roll(z, r), np.roll(n, r), q0, m),
        ring.ring_resonances(length, z[::-1], n[::-1], q0, m),
    ):
        assert abs(other.lower - pair.lower) <= 1e-10 * abs(pair.lower)
        assert abs(other.upper - pair.upper) <= 1e-10 * abs(pair.upper)


@settings(max_examples=60, deadline=None)
@given(profiles(), st.sampled_from([1, 2]), st.floats(10.0, 1e4))
def test_an_index_profile_alone_leaves_the_pair_one_frequency(profile: tuple[Arr, Arr], m: int, q0: float) -> None:
    # with Z0 the same on every segment nothing reflects: whatever the index does segment by segment, +-m stay one
    # frequency, m c / (n L) with n the mean index, carried to f / (1 - i/(2 Q0)) by the loss
    _, n = profile
    length = m * L915
    pair = ring.ring_resonances(length, np.full(len(n), 41.7), n, q0, m)
    exact = m * C0 / (float(np.mean(n)) * length) / (1 - 0.5j / q0)
    assert abs(pair.lower - exact) <= 1e-12 * abs(exact)
    assert abs(pair.upper - exact) <= 1e-12 * abs(exact)


def test_the_scored_drive_turns_with_the_ring() -> None:
    # ring_score with its feeds given: a drawn T915-S ring started 37 nodes on, its feeds moved with it, scores the same
    # drive (frequency, weakest depth per watt, |b/a|, the second feed's weight) under the fixed hybrid and under two
    # re-solved feeds, its node volts turned with it
    z, n, rel = ring.ring_draw(np.random.default_rng(9), L915, 0.02, ring.N915, 0.01, 0.0)
    turned = np.roll(z, 37), np.roll(n, 37), np.roll(rel, 37)
    for compensate in (False, True):
        here = ring.ring_score(L915, 1, z, n, rel, Q_BARE, compensate)
        there = ring.ring_score(L915, 1, *turned, Q_BARE, compensate, feeds=(37, 37 + ring.RING_NODES // 4))
        numbers = (here.f, here.score, here.purity, here.w_abs, here.w_deg)
        assert (there.f, there.score, there.purity, there.w_abs, there.w_deg) == pytest.approx(numbers, rel=1e-9)
        assert np.abs(there.v - np.roll(here.v, 37)).max() <= 1e-9 * np.abs(here.v).max()


def test_the_t915s_monte_carlo_at_2_percent_is_the_study_s() -> None:
    # section 5.3: the first three draws of T915-S's 2 % case (seed 7) at the line study's Q 213.8, drawn and scored by
    # the module with the lens's figures, are the study's ring_draw, ring_score and ring_harmonic's on the same seed,
    # number for number; on them the re-solved feeds hold E2-8's band (|b/a| <= 0.03, 0.9 of a perfect ring's depth per
    # watt) and the fixed hybrid leaves more backward wave and less depth; T915-S's own laminate figure, section B's
    # solved one, moves the draws' bands and not their scores
    case = ring.RING_CASES[0]
    assert (case.label, case.sigma, case.seed, case.draws, case.n_eff) == ("bare copper, 2 %", 0.02, 7, 60, 1.031)
    run = ring.tolerance_run(case, Q_BARE, draws=3, build=ring.LENS)
    own = ring.tolerance_run(case, Q_BARE, draws=3)
    for name in ("split", "hybrid_purity", "hybrid_depth", "resolved_purity", "resolved_depth"):
        np.testing.assert_allclose(getattr(own, name), getattr(run, name), rtol=1e-6, atol=1e-12)
    rng = np.random.default_rng(case.seed)
    flat = np.full(study.RING_NODES, 50.0), np.full(study.RING_NODES, study.N915), np.zeros(study.RING_NODES)
    ref = study.ring_score(study.L915, 1, *flat, Q_BARE, False)
    rows = []
    for _ in range(3):
        z, n, rel = study.ring_draw(rng, study.L915, 0.02, study.N915, 0.01, 0.0)
        fixed, comp = (study.ring_score(study.L915, 1, z, n, rel, Q_BARE, c) for c in (False, True))
        rows.append(
            (
                2 * Q_BARE * study.ring_harmonic(z),
                fixed["purity"],
                fixed["score"] / ref["score"],
                comp["purity"],
                comp["score"] / ref["score"],
                comp["w_abs"],
                comp["w_deg"],
            )
        )
    columns = (run.split, run.hybrid_purity, run.hybrid_depth, run.resolved_purity, run.resolved_depth)
    mine = np.array([*columns, run.weight_abs, run.weight_deg]).T
    np.testing.assert_allclose(mine, np.array(rows), rtol=1e-12, atol=0)
    assert run.length == pytest.approx(study.L915, rel=1e-15)
    assert run.perfect_purity == pytest.approx(ref["purity"], rel=1e-9, abs=1e-15)
    assert run.summary()["split_median"] == pytest.approx(float(np.median(np.array(rows)[:, 0])), rel=1e-15)
    assert run.resolved_purity.max() < 0.03
    assert run.resolved_depth.min() > 0.9
    assert np.all(run.hybrid_purity > run.resolved_purity)
    assert np.all(run.hybrid_depth < run.resolved_depth)


def test_the_ring_is_the_buildability_lens_model() -> None:
    # section 5.3: the ring of 192 segments draws and scores as lens5/ring.py did for ring915.py, to 1e-10
    length = C0 / (1.031 * ring.F915)
    rng = np.random.default_rng(7)
    perfect = ring.ring_score(length, 1, *ring.perfect_ring(1.031), 381.0, False)
    for expected in LENS5:
        z, n, rel = ring.ring_draw(rng, length, 0.02, 1.031, 0.01, 0.0, build=ring.LENS)
        fixed = ring.ring_score(length, 1, z, n, rel, 381.0, False)
        comp = ring.ring_score(length, 1, z, n, rel, 381.0, True)
        got = (
            2 * 381.0 * ring.ring_harmonic(z),
            fixed.purity,
            fixed.score / perfect.score,
            comp.purity,
            comp.score / perfect.score,
            comp.w_abs,
            comp.w_deg,
        )
        assert got == pytest.approx(expected, rel=0, abs=1e-10)


def test_the_strips_at_2450_mhz_as_the_lens_ran_them() -> None:
    # section 5.5: lens5's 12-draw Monte Carlo of the 2.45 GHz strips at 10 um rms on h (2.9 %), Q 250: median 2Q|c_2m|
    # 0.993, depth per watt 0.934 re-solved against 0.654 fixed
    lens = ring.strips_2450(250.0, 1.067, 12)
    assert lens["sigma_pct"] == pytest.approx(2.92, abs=0.005)
    assert lens["split_median"] == pytest.approx(0.9934, abs=1e-4)
    assert (lens["resolved_depth_median"], lens["hybrid_depth_median"]) == pytest.approx((0.9342, 0.6543), abs=1e-4)


def test_the_drive_lens_ring() -> None:
    # section 5.3 and figure 6: section A's ring at 915 MHz (Q 250.8 and Z0 73.5 ohm by the line study, 290 mm, R =
    # Z0 Q / (2 pi m) = 2934 ohm); with six feeds two cheap modes near sqrt(K) R / 2 over a floor 124 times lower, as
    # results.json holds them; one tolerance draw (seed 3) driven by one feed, the untuned hybrid and two compiled feeds
    # leaves figure 6's ripples 0.983, 0.649 and 0.557; with identical joints one feed stands and two compiled travel
    cases = ring.drive_cases(Q_A, Z0_A)
    assert [c.label for c in cases[:2]] == list(SINGULAR)
    assert cases[0].length * 1e3 == pytest.approx(289.95, abs=0.005)
    for case in cases[:2]:
        six = case.ring(6, None)
        assert six.shunt == pytest.approx(2934.13, abs=0.01)
        freqs = six.resonance() * (1 + np.linspace(-0.02, 0.02, 801))
        sv = np.array([np.linalg.svd(g, compute_uv=False) for g in six.operator(freqs)]) / math.sqrt(six.samples)
        peak = int(np.argmax(sv[:, 0]))
        values, f_mhz = SINGULAR[case.label]
        assert sv[peak] == pytest.approx(values, rel=1e-9)
        assert freqs[peak] / 1e6 == pytest.approx(f_mhz, rel=1e-12)
        assert sv[peak][0] / sv[peak][2] == pytest.approx(126, rel=0.03)
        assert sv[peak][0] == pytest.approx(math.sqrt(6) * six.shunt / 2, rel=0.02)
    drawn = cases[0].ring(2, np.random.default_rng(3))
    for k, fixed, expected in (
        (1, None, 0.9830639467812782),
        (2, np.array([1, -1j]), 0.6485151679588262),
        (2, None, 0.5573694494097506),
    ):
        feeds = ring.feed_layout(k, drawn.length, 1, 6, False)
        built = ring.Ring(drawn.length, drawn.index, drawn.q, drawn.z0, drawn.m, drawn.joints, feeds)
        assert ring.best_drive(built, fixed=fixed).ripple == pytest.approx(expected, rel=1e-9)
    assert ring.best_drive(cases[1].ring(1, None)).ripple > 0.9
    assert ring.best_drive(cases[1].ring(2, None)).ripple < 0.02


def test_the_line_matrices_and_the_nodal_admittance_are_one_model() -> None:
    # the drive study's chain matrices and the buildability lens's nodal admittance give the same volts for the same
    # ring: section A's with six joints eight segments long at 1.40 to 1.65 Z0 and three tuned feeds, its 240 samples
    # the nodal ring's odd nodes, at five frequencies across the resonance, to 1e-8; and the chain matrices are
    # reciprocal between the feeds, with tuned feeds or coupling capacitors
    samples, nodes = 240, 480
    length = ring.drive_cases(Q_A, Z0_A)[0].length
    step = length / nodes
    joints = tuple(((40 + 80 * j) * step, 8 * step, Z0_A * (1.4 + 0.05 * j)) for j in range(6))
    at = (3, 63, 150)
    feeds = tuple((i + 0.5) * length / samples for i in at)
    line = ring.Ring(length, ring.A_INDEX, Q_A, Z0_A, 1, joints, feeds, samples=samples, compensated=True)
    z = np.full(nodes, Z0_A)
    for j in range(6):
        z[36 + 80 * j : 44 + 80 * j] = Z0_A * (1.4 + 0.05 * j)
    n = np.full(nodes, ring.A_INDEX)
    freqs = line.resonance() * np.array([0.97, 0.995, 1.0, 1.004, 1.03])
    chain = line.operator(freqs)
    for f, g in zip(freqs, chain, strict=True):
        nodal = ring.ring_response(float(f), length, z, n, Q_A, [2 * i + 1 for i in at], source=line.conductance)
        assert np.abs(nodal[1::2] - g).max() <= 1e-8 * np.abs(g).max()
    capacitors = ring.Ring(length, ring.A_INDEX, Q_A, Z0_A, 1, joints, feeds, samples=samples).operator(freqs)
    for g in (*chain, *capacitors):
        between = g[list(at)]
        assert np.abs(between - between.T).max() <= 1e-11 * np.abs(g).max()


def test_compiled_feeds_remove_the_backward_wave_the_hybrid_leaves() -> None:
    # section 5.3 on the ring of segments: on T915-S's 2 % draws two or three current feeds compiled from the ring's
    # own G at its pair's centre leave |b/a| under 0.003, against the fixed hybrid's 0.07 to 0.15, and keep 0.9 of a
    # perfect ring's depth per watt; on a perfect ring three equal feeds compile to the travelling wave's phases
    rng = np.random.default_rng(7)
    perfect = ring.ring_score(L915, 1, *ring.perfect_ring(ring.N915), Q_BARE, False)
    for _ in range(3):
        z, n, rel = ring.ring_draw(rng, L915, 0.02, ring.N915, 0.01, 0.0)
        hybrid = ring.ring_score(L915, 1, z, n, rel, Q_BARE, False)
        assert hybrid.purity > 0.05
        for feeds in ((0, 48), (0, 64, 128)):
            compiled = ring.compile_feeds(L915, 1, z, n, rel, Q_BARE, feeds)
            assert compiled.purity < 0.003
            assert compiled.score / perfect.score > 0.9
    z, n, rel = ring.perfect_ring(ring.N915)
    three = ring.compile_feeds(L915, 1, z, n, rel, Q_BARE, (0, 64, 128))
    assert three.purity < 1e-10
    assert three.weights == pytest.approx(np.exp(-2j * math.pi * np.arange(3) / 3), abs=1e-9)


def test_each_tolerance_draws_what_it_says() -> None:
    # ring_draw's parts one at a time on T915-S's ring: no tolerance is the ring as drawn; the closing joints alone are
    # two steps of joint_dz a scarf long; the spacers' errors give 2Q|c_2m| near 2Q S sigma / sqrt(N) in rms (section
    # 5.2's spacer rule, 0.87 linewidths at Q 214 and 2 %); the split's stations and the laminate draw their sd; one
    # seed, one ring
    rng = np.random.default_rng(5)
    z, n, rel = ring.ring_draw(rng, L915, 0.0, ring.N915, 0.0, 0.0, joint_dz=0.0, sigma_eps=0.0)
    assert np.array_equal(z, np.full(ring.RING_NODES, 50.0))
    assert np.array_equal(n, np.full(ring.RING_NODES, ring.N915))
    assert not rel.any()
    z, _, _ = ring.ring_draw(rng, L915, 0.0, ring.N915, 0.0, 0.0)
    steps = np.log(z / 50.0)
    assert set(np.round(steps, 12)) == {0.0, -0.02}
    inside = int(np.count_nonzero(steps))
    assert (
        2 * math.floor(ring.SCARF / (L915 / ring.RING_NODES))
        <= inside
        <= 2 * math.ceil(ring.SCARF / (L915 / ring.RING_NODES))
    )
    draws = [ring.ring_draw(rng, L915, 0.02, ring.N915, 0.0, 0.0, joint_dz=0.0) for _ in range(1000)]
    rms = math.sqrt(float(np.mean([(2 * Q_BARE * ring.ring_harmonic(z)) ** 2 for z, _, _ in draws])))
    rule = 2 * Q_BARE * ring.RING_SENSITIVITY * 0.02 / math.sqrt(ring.SPACERS)
    assert rule == pytest.approx(0.873, abs=5e-4)
    assert rms == pytest.approx(rule, rel=0.1)
    assert float(np.std([n[0] for _, n, _ in draws])) == pytest.approx(ring.T915S.index_per_eps * 0.1, rel=0.1)
    # section B's n share / (2 eps_r), DR-0's share rule, where the lens took 2.2 times it
    assert ring.T915S.index_per_eps == pytest.approx(1.0311 * 0.07366 / (2 * 3.3), rel=1e-3)
    assert ring.LENS.index_per_eps / ring.T915S.index_per_eps == pytest.approx(2.2, abs=0.01)
    split = [ring.ring_draw(rng, L915, 0.0, ring.N915, 0.0, 0.2, joint_dz=0.0)[0] for _ in range(1000)]
    assert float(np.std([np.log(z[0] / 50.0) for z in split])) == pytest.approx(0.03 * 0.2, rel=0.1)
    one, two = (ring.ring_draw(np.random.default_rng(11), L915, 0.02, ring.N915, 0.01, 0.2) for _ in range(2))
    assert all(np.array_equal(a, b) for a, b in zip(one, two, strict=True))
