"""The antenna model against textbook results, energy conservation and its own physics."""
import numpy as np
import pytest

from array16 import Array, antenna as A, board as B, dsp
from array16.sim import BoardModel, SimBackend, Source

F0, K0 = B.F0, 2 * np.pi * B.F0 / B.C0
FR4 = A.Substrate("FR-4, ground on L3", 4.4, 0.018, 1.275e-3, rough=1.8)
AIR = A.Substrate("air, lossless", 1.0, 0.0, 4e-3, sigma=1e30)


@pytest.fixture(scope="module")
def fr4_patch():
    return A.design(FR4, probe_r=0.15e-3)


@pytest.fixture(scope="module")
def air_array():
    return A.PatchArray(A.design(AIR, probe_r=0.5e-3, probe_comp=True))


def balanis_slot_integral(k, W, L=None):
    """Balanis 14-17 (L None) and 14-18a: the slot conductances by integration over angle, with 120 pi -> eta0."""
    t = np.linspace(1e-6, np.pi - 1e-6, 20001)
    c = np.cos(t)
    f = (np.sin(k * W * c / 2) / c) ** 2 * np.sin(t) ** 3
    if L is not None:
        tau = np.linspace(0, np.pi, 801)
        f = f * np.trapezoid(np.cos(k * L * np.sin(t)[:, None] * np.sin(tau)), tau, axis=1) / np.pi     # J0
    return np.trapezoid(f, t) / (np.pi * A.ETA0)


def test_closed_form_design_matches_balanis_example_14_1():
    W, L = A.textbook_patch(2.2, 1.588e-3, 10e9)                    # Balanis: W 1.186 cm, eps_eff 1.972, dL 0.081 cm, L 0.906 cm
    assert abs(W - 11.86e-3) < 0.02e-3 and abs(L - 9.06e-3) < 0.02e-3
    assert abs(A.eps_eff(2.2, 1.588e-3, W) - 1.972) < 0.002 and abs(A.open_end(2.2, 1.588e-3, W) - 0.81e-3) < 0.005e-3


@pytest.mark.parametrize("W,L", [(15.7e-3, 11.9e-3), (22e-3, 23e-3), (15.7e-3, 26e-3)])
def test_slot_conductances_match_balanis(W, L):
    assert A.slot_conductance(K0, W) == pytest.approx(balanis_slot_integral(K0, W), rel=2e-4)
    assert A.slot_mutual(K0, W, L, 0.0).real == pytest.approx(balanis_slot_integral(K0, W, L), rel=2e-3)


def test_energy_balance_and_reciprocity(air_array):
    r = air_array.solve(F0)
    S = r["S"]
    assert np.max(np.abs(S - S.T)) < 1e-12
    for w in (np.eye(16)[0], np.eye(16)[5], np.exp(-1j * K0 * (B.X * 0.4 + B.Y * 0.2))):
        p_in = (np.sum(np.abs(w) ** 2) - np.sum(np.abs(S @ w) ** 2)) / 2
        assert air_array.radiated_power(F0, w) == pytest.approx(p_in, rel=1e-4)    # far field = power into the ports


def test_embedded_gain_approaches_the_cell_aperture(air_array):
    # A lossless, matched element deep in an array at broadside has the gain of its cell: 4 pi A / lambda^2
    cell = 10 * np.log10(4 * np.pi * B.PITCH ** 2 / (B.C0 / F0) ** 2)
    g = 10 * np.log10(abs(air_array.eep(np.array(0.0), np.array(0.0), F0)[5]) ** 2)
    assert abs(g - cell) < 1.5


def test_fr4_patch_agrees_with_jackson_alexopoulos(fr4_patch):
    p, er, h = fr4_patch, FR4.eps_r, FR4.h
    iso = A.PatchArray(p, [0.0], [0.0])
    s11 = iso.solve(F0)["S"][0, 0]
    assert abs(s11) < 0.1                                           # designed for a match at 5.8 GHz
    eff = iso.radiated_power(F0, [1.0]) / ((1 - abs(s11) ** 2) / 2)
    ja = A.jackson_alexopoulos(FR4, p.W, p.L)
    assert eff == pytest.approx(ja["efficiency"], rel=0.05)
    assert 0.45 < eff < 0.65
    assert ja["q_total"] == pytest.approx(18.3, abs=0.5)               # by hand: Q_sp 34.1, Q_sw 161, Q_d 55.6, Q_c 816


def test_coupled_patches_in_the_simulator(fr4_patch):
    arr = A.PatchArray(fr4_patch)
    b = SimBackend("quad", seed=3, board=BoardModel(4, antenna=arr))
    a = Array(b)
    b.sources = [Source(0, 0, -75)]
    a.calibrate()
    e0 = arr.eep(np.array(0.0), np.array(0.0), F0)
    steer = lambda uu, vv: arr.eep(uu, vv, F0) / e0
    for d in ((20, 30), (40, 200)):
        b.sources = [Source(*d, power_dbm=-75)]
        z = a.snapshot()
        assert dsp.angle_between(B.angles(*dsp.locate(z, F0)), d) < 3          # boresight calibration: coupling bias
        assert dsp.angle_between(B.angles(*dsp.locate(z, F0, steer)), d) < 1   # exact array response: the hardware floor


def test_coupling_stays_near_the_slab_greens_function_reference(fr4_patch):
    # Reference: the model review's spectral-domain moment method with the grounded-slab Green's function
    # (surface waves included), centre element of the 4x4 on 1.28 mm FR-4: E -20.7, H -18.8, D -28.5 dB.
    # Known biases of the slot model at this pitch: E about right, H 1.5-2 dB weak, diagonal 2-3 dB weak.
    S = A.PatchArray(fr4_patch).solve(F0)["S"]
    e, h, d = (20 * np.log10(abs(S[5, j])) for j in (6, 9, 10))
    assert abs(e - (-20.7)) < 1.0
    assert -2.5 < h - (-18.8) < 0
    assert -3.5 < d - (-28.5) < 0


def test_close_facing_edges_are_flagged(air_array, fr4_patch):
    assert A.PatchArray(fr4_patch).coupling_reliable
    assert not air_array.coupling_reliable                          # 26 mm pitch, sheet patches about 20 mm long


def test_receive_coupling_fit_reproduces_the_embedded_patterns(fr4_patch):
    C, resid = A.fit_coupling(A.PatchArray(fr4_patch))
    assert resid < 0.02
    assert -23 < 20 * np.log10(abs(C[5, 6] / C[5, 5])) < -18          # E-plane neighbour about -21 dB


def test_coupling_model_with_no_coupling_is_the_default_element():
    ideal = A.CoupledIdeal()
    bm = BoardModel(1)
    for d in ((0, 0), (30, 45), (60, -120)):
        u, v = B.direction(*d)
        assert np.allclose(ideal.eep(np.array(u), np.array(v), F0), bm.elem_resp(u, v, F0))
