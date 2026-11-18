"""Topology and software against the classic physical limits of a small array: seven numerical demonstrations.

Companion to experiments/evolve_layout.py (whip layouts) for the rev B redesign (REDESIGN.md). Each demo
states its assumptions in its docstring; all are idealised (isotropic elements, no mutual coupling unless a
calibration residual is named), so they bound what the technique can do, not what a built array will do.

  1 coarray    16 sensors on a line: nested / coprime / minimum-redundancy layouts and co-array (spatially
               smoothed) root-MUSIC resolving more than 15 uncorrelated sources; what breaks it.
               Also the rank of the co-array manifold of the 16-whip plate layout per frequency.
  2 freqdiv    frequency diversity against the ambiguity of a sparse layout (4-whip ring far above its
               0.354-wavelength limit); Chinese-remainder two-baseline unwrapping and its noise tolerance.
  3 rotation   rotation synthesis: turntable under the 16-whip plate, and Earth rotation for a 10 m field
               layout at the 1420 MHz hydrogen line; uv coverage and dirty-beam sidelobes.
  4 mimo       MIMO-radar virtual arrays: panel (5.8 GHz) and whips (915 MHz); beamwidth, sidelobes, and
               the Doppler price of time-division transmit.
  5 switched   rev B's snapshot switching: which covariance entries 4 channels can measure, and how many
               sources a covariance fit then resolves; what a bursty (non-stationary) source does to it.
  6 threshold  MUSIC's two-source resolution threshold against SNR and calibration residual.
  7 astro      sensitivity at the radio-astronomy bands (408, 611, 1420.4, 1665.4, 2695, 4830, 4995 MHz):
               which sources 16 elements can detect, and what switching 4 of 16 costs.

    python experiments/topology_study.py            # all, about 2 minutes on 3 threads
    python experiments/topology_study.py mimo       # one demo
"""
from __future__ import annotations

import json
import sys
import time
from math import erf
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

C0 = 299_792_458.0
HERE = Path(__file__).resolve().parent
OUT = HERE / "topology"
LAYOUTS = json.loads((HERE / "evolve" / "results.json").read_text())["designs"]
WHIPS16 = np.array(LAYOUTS["rings16_225_112mm"]["positions_m"])      # two rings of 8 (225 / 112 mm)
EVOLVED16 = np.array(LAYOUTS["evolved16"]["positions_m"])
RING4 = np.array(LAYOUTS["ring4_133mm"]["positions_m"])


def nnls(A, b, maxiter=None, tol=1e-3):
    """Lawson-Hanson active-set non-negative least squares (the environment has no scipy)."""
    m, n = A.shape
    maxiter = maxiter or 3 * n
    x = np.zeros(n)
    P = np.zeros(n, bool)
    w = A.T @ (b - A @ x)
    w0 = np.abs(w).max()
    it = 0
    while (~P).any() and w[~P].max() > tol * w0 and it < maxiter:
        j = np.argmax(np.where(~P, w, -np.inf))
        P[j] = True
        while True:
            it += 1
            s = np.zeros(n)
            s[P] = np.linalg.lstsq(A[:, P], b, rcond=None)[0]
            if s[P].min() > 0 or it >= maxiter:
                break
            neg = P & (s <= 0)
            alpha = np.min(x[neg] / (x[neg] - s[neg]))
            x = x + alpha * (s - x)
            P &= x > 1e-14
        x = np.maximum(s, 0)
        w = A.T @ (b - A @ x)
    return x, float(np.linalg.norm(A @ x - b))


def cn(rng, shape):
    return (rng.normal(size=shape) + 1j * rng.normal(size=shape)) / np.sqrt(2)


def wrap_deg(a):
    return (a + 180.0) % 360.0 - 180.0


# =========================================================================================================
# 1  difference co-array
# =========================================================================================================
def nested(n1, n2):
    return np.r_[np.arange(1, n1 + 1), (n1 + 1) * np.arange(1, n2 + 1)] - 1


def coprime(m, n):
    """Extended coprime pair: n sensors at m*i and 2m sensors at n*i (share 0): 2m + n - 1 sensors."""
    return np.unique(np.r_[m * np.arange(n), n * np.arange(2 * m)])


def wichmann(r, s):
    """Wichmann's sparse ruler: 4r + s + 3 marks, every lag up to 4r(r+s+2) + 3(s+1) present."""
    gaps = [1] * r + [r + 1] + [2 * r + 1] * r + [4 * r + 3] * s + [2 * r + 2] * (r + 1) + [1] * r
    return np.r_[0, np.cumsum(gaps)]


def best_wichmann(n):
    best = None
    for r in range(0, n):
        s = n - 3 - 4 * r
        if s < 0:
            break
        p = wichmann(r, s)
        if best is None or contiguous(p) > contiguous(best):
            best = p
    return best if best is not None else np.arange(n)


def contiguous(pos):
    """Largest L such that every lag 0..L occurs among the pairwise differences."""
    lags = set(np.abs(pos[:, None] - pos[None, :]).ravel().tolist())
    L = 0
    while L + 1 in lags:
        L += 1
    return L


def coarray_root_music(R, pos, K):
    """Average the sample covariance over equal lags, build the (L+1)x(L+1) Hermitian Toeplitz matrix (equal
    to the spatially smoothed co-array covariance up to a square, Liu & Vaidyanathan 2015) and root-MUSIC it.
    Returns K estimates of u = sin(theta), sensors in half wavelengths."""
    L = contiguous(pos)
    d = pos[:, None] - pos[None, :]
    z = np.array([R[d == l].mean() for l in range(L + 1)])
    i = np.arange(L + 1)
    lag = i[:, None] - i[None, :]
    T = np.where(lag >= 0, z[np.abs(lag)], z[np.abs(lag)].conj())
    w, V = np.linalg.eigh(T)
    En = V[:, np.argsort(-np.abs(w))[K:]]
    Cm = En @ En.conj().T
    coef = np.array([np.trace(Cm, offset=k) for k in range(L, -L - 1, -1)])
    roots = np.roots(coef)
    roots = roots[np.abs(roots) < 1]
    roots = roots[np.argsort(1 - np.abs(roots))[:K]]
    return np.sort(np.angle(roots) / np.pi)


def element_root_music(R, K):
    """Plain root-MUSIC on a uniform line (needs K < N)."""
    n = len(R)
    w, V = np.linalg.eigh(R)
    En = V[:, : n - K]
    Cm = En @ En.conj().T
    coef = np.array([np.trace(Cm, offset=k) for k in range(n - 1, -n, -1)])
    roots = np.roots(coef)
    roots = roots[np.abs(roots) < 1]
    roots = roots[np.argsort(1 - np.abs(roots))[:K]]
    return np.sort(np.angle(roots) / np.pi)


def line_trial(rng, pos, u, snr_db, T, cal_deg=0.0, coherent_pairs=False):
    K = len(u)
    A = np.exp(1j * np.pi * pos[:, None] * u[None, :])
    S = cn(rng, (K, T))
    if coherent_pairs:                      # sources 2k and 2k+1 carry the same waveform (a direct path and an echo)
        S[1::2] = S[0::2] * np.exp(1j * rng.uniform(0, 2 * np.pi, (K // 2, 1)))
    g = np.exp(1j * np.radians(rng.normal(0, cal_deg, len(pos)))) * 10 ** (rng.normal(0, 0.5 if cal_deg else 0, len(pos)) / 20)
    X = g[:, None] * (A @ S) + cn(rng, (len(pos), T)) * 10 ** (-snr_db / 20)
    return X @ X.conj().T / T


def demo_coarray(quick=False):
    """Assumptions: narrowband far-field sources, mutually UNCORRELATED and stationary over the T snapshots,
    equal power, SNR per source per sensor; sensors on a half-wavelength grid; white noise; number of sources
    known. 'ok' = every source matched to within half the source spacing."""
    rng = np.random.default_rng(1)
    res = {}
    # degrees of freedom against the number of sensors
    tab = []
    for n in (4, 6, 8, 12, 16, 24, 32):
        nst = nested(n // 2, n - n // 2)
        cps = [(m, n + 1 - 2 * m) for m in range(2, n // 2) if n + 1 - 2 * m > 1 and np.gcd(m, n + 1 - 2 * m) == 1]
        cp = max((contiguous(coprime(m, k)), (m, k)) for m, k in cps) if cps else (n - 1, None)
        mra = best_wichmann(n)
        tab.append(dict(n=n, ula=n - 1, nested=contiguous(nst), nested_aperture=int(nst.max()),
                        coprime=cp[0], coprime_pair=cp[1], wichmann=contiguous(mra), wichmann_aperture=int(mra.max()),
                        pairs_bound=n * (n - 1) // 2))
    res["sources_vs_n"] = tab
    arrays = {"ula16": np.arange(16), "nested8+8": nested(8, 8), "coprime(5,7)": coprime(5, 7), "wichmann(2,5)": wichmann(2, 5)}
    res["arrays"] = {k: dict(positions=v.tolist(), contiguous_lags=contiguous(v), aperture=int(v.max())) for k, v in arrays.items()}

    trials = 20 if quick else 60
    pos = arrays["nested8+8"]
    K = 24
    u0 = np.linspace(-0.9, 0.9, K)
    half = (u0[1] - u0[0]) / 2

    def run(pos, u0, snr, T, **kw):
        errs, ok = [], 0
        for _ in range(trials):
            u = np.sort(u0 + rng.uniform(-0.2, 0.2, len(u0)) * (u0[1] - u0[0]))
            est = coarray_root_music(line_trial(rng, pos, u, snr, T, **kw), pos, len(u))
            if len(est) < len(u):
                errs.append(np.full(len(u), 1.0)); continue
            e = est - u
            errs.append(e)
            ok += np.max(np.abs(e)) < (u0[1] - u0[0]) / 2
        e = np.concatenate(errs)
        return dict(rmse_u=float(np.sqrt(np.mean(e ** 2))), ok=ok / trials)

    res["k24_vs_snapshots"] = [dict(snr_db=s, T=T, **run(pos, u0, s, T)) for s in (0, 10, 20) for T in (100, 300, 1000, 3000, 10000)]
    res["vs_k"] = [dict(K=k, **run(pos, np.linspace(-0.9, 0.9, k), 10, 3000)) for k in (8, 15, 24, 40, 56, 64, 71)]
    res["k24_cal3deg"] = run(pos, u0, 10, 3000, cal_deg=3.0)
    res["k24_coherent_pairs"] = run(pos, u0, 10, 3000, coherent_pairs=True)
    res["k24_wichmann"] = run(arrays["wichmann(2,5)"], u0, 10, 3000)
    res["k24_coprime"] = run(arrays["coprime(5,7)"], u0, 10, 3000)
    # the ordinary case, for scale: 8 sources, ULA element-space MUSIC against nested co-array MUSIC
    u8 = np.linspace(-0.8, 0.8, 8)
    e_ula, e_n = [], []
    for _ in range(trials):
        u = np.sort(u8 + rng.uniform(-0.2, 0.2, 8) * (u8[1] - u8[0]))
        est = element_root_music(line_trial(rng, arrays["ula16"], u, 10, 300), 8)
        e_ula.append(est - u if len(est) == 8 else np.ones(8))
        e_n.append(coarray_root_music(line_trial(rng, pos, u, 10, 300), pos, 8) - u)
    res["k8_T300_rmse_u"] = dict(ula16=float(np.sqrt(np.mean(np.concatenate(e_ula) ** 2))),
                                 nested=float(np.sqrt(np.mean(np.concatenate(e_n) ** 2))))

    # the plate: rank of the co-array manifold {vec(a a^H)} over azimuth, per frequency
    az = np.radians(np.arange(0, 360, 0.25))
    rk = []
    for f in (0.408e9, 0.611e9, 0.915e9, 1.4204e9, 1.6654e9, 2.3e9):
        row = dict(f_mhz=f / 1e6)
        for name, p in (("rings16", WHIPS16), ("evolved16", EVOLVED16)):
            k = 2 * np.pi * f / C0
            A = np.exp(1j * k * (np.cos(az)[:, None] * p[:, 0] + np.sin(az)[:, None] * p[:, 1]))     # (G, N)
            B = (A[:, :, None] * A[:, None, :].conj()).reshape(len(az), -1)
            sa = np.linalg.svd(A, compute_uv=False)
            sb = np.linalg.svd(B, compute_uv=False)
            dmax = np.max(np.linalg.norm(p[:, None] - p[None], axis=2))
            row[name] = dict(element_rank_1e2=int(np.sum(sa > 1e-2 * sa[0])), coarray_rank_1e2=int(np.sum(sb > 1e-2 * sb[0])),
                             coarray_rank_1e4=int(np.sum(sb > 1e-4 * sb[0])), modes_2kD_plus_1=float(2 * k * dmax + 1),
                             modes_2kR_plus_1=float(2 * k * np.max(np.linalg.norm(p, axis=1)) + 1))
        rk.append(row)
    res["plate_coarray_rank"] = rk

    # figure: one co-array MUSIC spectrum with 24 sources
    u = u0
    R = line_trial(np.random.default_rng(7), pos, u, 10, 3000)
    L = contiguous(pos)
    d = pos[:, None] - pos[None, :]
    z = np.array([R[d == l].mean() for l in range(L + 1)])
    i = np.arange(L + 1)
    lag = i[:, None] - i[None, :]
    Tm = np.where(lag >= 0, z[np.abs(lag)], z[np.abs(lag)].conj())
    w, V = np.linalg.eigh(Tm)
    En = V[:, np.argsort(-np.abs(w))[K:]]
    ug = np.linspace(-1, 1, 8001)
    sp = 1 / np.sum(np.abs(En.conj().T @ np.exp(1j * np.pi * i[:, None] * ug[None])) ** 2, axis=0)
    Ru = line_trial(np.random.default_rng(7), arrays["ula16"], u, 10, 3000)
    bart = np.real(np.einsum("ng,nm,mg->g", np.exp(1j * np.pi * np.arange(16)[:, None] * ug[None]).conj(), Ru,
                             np.exp(1j * np.pi * np.arange(16)[:, None] * ug[None])))
    fig, ax = plt.subplots(figsize=(9, 3.6))
    ax.plot(ug, 10 * np.log10(sp / sp.max()), lw=0.9, label="nested 8+8, co-array MUSIC")
    ax.plot(ug, 10 * np.log10(bart / bart.max()), lw=0.9, label="uniform 16, beamscan")
    for x in u:
        ax.axvline(x, color="0.8", lw=0.5, zorder=0)
    ax.set(xlabel="u = sin(angle)", ylabel="dB", ylim=(-60, 2), title="24 uncorrelated sources, 16 sensors, 10 dB, 3000 snapshots")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "coarray_music.png", dpi=130); plt.close(fig)
    return res


# =========================================================================================================
# 2  frequency diversity and Chinese-remainder unwrapping
# =========================================================================================================
def demo_freqdiv(quick=False):
    """Assumptions: one source 0-15 deg above the horizon, azimuth-only Bartlett beamscan (assumes 0 deg);
    3 deg / 0.5 dB calibration residual per whip, the SAME at every tone; 20 dB per-whip SNR and 64 samples
    in total, split evenly over J = 8 tones spread across the stated fractional bandwidth (an emitter that
    really occupies that band: a hopper, a wideband burst, or our own stepped radar/calibration signal);
    spectra added incoherently. Gross error = more than 10 deg."""
    rng = np.random.default_rng(2)
    trials = 400 if quick else 2500
    grid = np.radians(np.arange(0, 360, 0.5))
    res = {"ring": []}
    J = 8
    snr_eff = 20 + 10 * np.log10(64 / J)
    for name, p, fc in (("ring4_133mm", RING4, 1.4204e9), ("ring4_133mm", RING4, 2.3e9),
                        ("ring8_225mm", np.array(LAYOUTS["ring8_225mm"]["positions_m"]), 2.3e9)):
        n = len(p)
        pg = np.cos(grid)[:, None] * p[:, 0] + np.sin(grid)[:, None] * p[:, 1]
        for fb in (0.0, 0.01, 0.024, 0.05, 0.1, 0.2, 0.4, 0.67):
            f = fc * (1 + fb * np.linspace(-0.5, 0.5, J))
            az = rng.uniform(0, 2 * np.pi, trials)
            ce = np.cos(np.radians(rng.uniform(0, 15, trials)))
            g = 10 ** (rng.normal(0, 0.5, (trials, n)) / 20) * np.exp(1j * np.radians(rng.normal(0, 3.0, (trials, n))))
            pt = (np.cos(az)[:, None] * p[:, 0] + np.sin(az)[:, None] * p[:, 1]) * ce[:, None]
            pw = np.zeros((trials, len(grid)))
            for fj in f:
                k = 2 * np.pi * fj / C0
                y = np.exp(1j * k * pt) * g + cn(rng, (trials, n)) * 10 ** (-snr_eff / 20)
                pw += np.abs(y @ np.exp(-1j * k * pg).T) ** 2
            e = np.abs(wrap_deg(np.degrees(grid[np.argmax(pw, 1)] - az)))
            res["ring"].append(dict(layout=name, fc_mhz=fc / 1e6, r_over_lambda=float(np.linalg.norm(p[0]) * fc / C0),
                                    frac_bw=fb, gross_pct=float(np.mean(e > 10) * 100), median_deg=float(np.median(e)),
                                    p95_deg=float(np.percentile(e, 95))))
    # Chinese remainder: baselines p*g and q*g, g = half a wavelength, p and q coprime
    crt = []
    T = 20000 if not quick else 4000
    for p_, q_ in ((2, 3), (3, 4), (5, 7), (8, 9), (15, 16)):
        for sig in (2.0, 4.2, 8.0, 12.0, 20.0):
            x = rng.uniform(-0.5, 0.5, T)                       # g sin(theta) / lambda
            e1, e2 = np.radians(rng.normal(0, sig, (2, T)))
            ph1 = np.angle(np.exp(1j * (2 * np.pi * p_ * x + e1)))
            ph2 = np.angle(np.exp(1j * (2 * np.pi * q_ * x + e2)))
            # search the integer pair that makes the two baselines agree best
            n1 = np.arange(-(p_ // 2) - 1, p_ // 2 + 2)
            n2 = np.arange(-(q_ // 2) - 1, q_ // 2 + 2)
            x1 = (ph1[:, None] / (2 * np.pi) + n1[None]) / p_                    # (T, |n1|)
            x2 = (ph2[:, None] / (2 * np.pi) + n2[None]) / q_
            dd = np.abs(x1[:, :, None] - x2[:, None, :])
            dd[(np.abs(x1) > 0.5 + 0.5 / p_)[:, :, None].repeat(len(n2), 2)] = np.inf
            ij = dd.reshape(T, -1).argmin(1)
            i1, i2 = np.unravel_index(ij, (len(n1), len(n2)))
            xh = (p_ ** 2 * x1[np.arange(T), i1] + q_ ** 2 * x2[np.arange(T), i2]) / (p_ ** 2 + q_ ** 2)
            err = np.abs(xh - x)
            err = np.minimum(err, 1 - err)                                         # endfire wraps
            good = err < 0.25 / q_
            crt.append(dict(p=p_, q=q_, sigma_deg=sig, p_correct=float(good.mean()),
                            p_theory=float(erf(np.pi / (np.sqrt(2 * (p_ ** 2 + q_ ** 2)) * np.radians(sig)))),
                            rms_x_when_correct=float(np.sqrt(np.mean(err[good] ** 2))),
                            rms_x_single_halfwave=float(np.radians(sig) / (2 * np.pi))))
    res["crt"] = crt

    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    for key in sorted({(r["layout"], r["fc_mhz"]) for r in res["ring"]}):
        rows = [r for r in res["ring"] if (r["layout"], r["fc_mhz"]) == key]
        ax.plot([max(r["frac_bw"], 0.005) * 100 for r in rows], [r["gross_pct"] for r in rows], "o-",
                label=f"{key[0]} at {key[1]:.0f} MHz (r = {rows[0]['r_over_lambda']:.2f} wavelengths)")
    ax.set(xscale="log", xlabel="occupied bandwidth, % of centre (leftmost point = single tone)", ylabel="gross errors, %",
           title="Frequency diversity against sparse-ring ambiguity")
    ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(OUT / "freq_diversity.png", dpi=130); plt.close(fig)
    return res


# =========================================================================================================
# 3  rotation synthesis
# =========================================================================================================
def dirty_beam(uv, lim=1.0, npix=201):
    """Cross-correlation-only point-spread function over direction-cosine offsets (l, m) in [-lim, lim]."""
    ax = np.linspace(-lim, lim, npix)
    l, m = np.meshgrid(ax, ax)
    b = np.zeros_like(l)
    for i in range(0, len(uv), 2048):
        c = uv[i:i + 2048]
        b += np.cos(2 * np.pi * (c[:, 0, None, None] * l + c[:, 1, None, None] * m)).sum(0)
    return ax, b / len(uv)


def beam_stats(uv, lim=1.0):
    ax, b = dirty_beam(uv, lim)
    r = np.hypot(*np.meshgrid(ax, ax))
    umax = np.max(np.linalg.norm(uv, axis=1))
    # main lobe: out to the first minimum of the azimuthal average (rings two pixels wide)
    step = 2 * (ax[1] - ax[0])
    ib = np.rint(r / step).astype(int).ravel()
    prof = np.bincount(ib, b.ravel()) / np.maximum(np.bincount(ib), 1)
    rb = np.arange(len(prof)) * step
    i0 = next((i for i in range(1, len(prof) - 1) if prof[i] < prof[i - 1] and prof[i] <= prof[i + 1]), 3)
    out = (r > rb[i0]) & (r <= lim)
    ih = int(np.argmax(prof < 0.5))
    hp = rb[ih - 1] + step * (prof[ih - 1] - 0.5) / (prof[ih - 1] - prof[ih]) if ih > 0 else step
    cells = {(int(np.floor(a / 0.5)), int(np.floor(c / 0.5))) for a, c in uv}
    cells |= {(-a - 1, -c - 1) for a, c in cells}
    total = sum(1 for a in range(-int(umax * 2) - 1, int(umax * 2) + 1) for c in range(-int(umax * 2) - 1, int(umax * 2) + 1)
                if np.hypot((a + 0.5) * 0.5, (c + 0.5) * 0.5) <= umax)
    return dict(n_uv=int(len(uv)), umax_wavelengths=float(umax), hpbw_deg=float(np.degrees(2 * np.arcsin(min(1.0, hp)))),
                first_min_deg=float(np.degrees(np.arcsin(min(1.0, rb[i0])))), peak_sidelobe=float(np.max(np.abs(b[out]))),
                rms_sidelobe=float(np.sqrt(np.mean(b[out] ** 2))), uv_cells_filled=float(min(1.0, len(cells) / max(total, 1)))), (ax, b)


def baselines(p):
    i, j = np.triu_indices(len(p), 1)
    return p[i] - p[j]


def demo_rotation(quick=False):
    """Assumptions: the sky does not change while the aperture is synthesised (minutes for the turntable,
    hours for Earth rotation); gains are calibrated across the whole synthesis; a planar array, so the
    response depends only on the horizontal direction cosines (l, m) and the point-spread function is
    shift-invariant; cross-correlations only, natural weighting. Sidelobes are relative to the peak (1 =
    as strong as the main lobe); they measure how well several sources can be told apart in one image."""
    res = {"turntable": [], "earth": []}
    figs = {}
    for f in (1.4204e9, 2.3e9):
        lam = C0 / f
        for name, p, steps, span in (("rings16 static", WHIPS16, 1, 0), ("rings16 turntable 45 deg / 16 steps", WHIPS16, 16, 45),
                                     ("evolved16 static", EVOLVED16, 1, 0), ("evolved16 turntable 180 deg / 60 steps", EVOLVED16, 60, 180),
                                     ("ring4_133mm static", RING4, 1, 0), ("ring4_133mm turntable 90 deg / 30 steps", RING4, 30, 90),
                                     ("2 whips 450 mm apart, turntable 180 deg / 60 steps", np.array([[0.225, 0], [-0.225, 0]]), 60, 180)):
            b = baselines(p) / lam
            uv = []
            for s in range(steps):
                a = np.radians(span * s / steps)
                rot = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
                uv.append(b @ rot.T)
            uv = np.concatenate(uv)
            st, im = beam_stats(uv)
            res["turntable"].append(dict(f_mhz=f / 1e6, case=name, **st))
            figs[(f, name)] = (uv, im)
    # Earth rotation: 16 antennas on a 10 m field (the plate layout scaled 20x: rings of 4.5 m and 2.24 m), 1420.4 MHz
    lam = C0 / 1.4204e9
    lat = np.radians(40.0)                                      # ASSUMED site latitude
    p = EVOLVED16 * 20.0
    bl = baselines(p)                                           # east, north (up = 0)
    E, N = bl[:, 0], bl[:, 1]
    X, Y, Z = -N * np.sin(lat), E, N * np.cos(lat)              # equatorial components
    for src, dec in (("Cas A", 58.8), ("Cyg A", 40.7), ("Sun at equinox", 0.0)):
        d = np.radians(dec)
        for hours in (0.0, 2.0, 6.0):
            H = np.radians(15 * np.linspace(-hours / 2, hours / 2, max(1, int(hours * 30) + 1)))
            u = (np.sin(H)[:, None] * X + np.cos(H)[:, None] * Y) / lam
            v = (-np.sin(d) * np.cos(H)[:, None] * X + np.sin(d) * np.sin(H)[:, None] * Y + np.cos(d) * Z) / lam
            uv = np.column_stack([u.ravel(), v.ravel()])
            st, im = beam_stats(uv, lim=0.15)
            res["earth"].append(dict(source=src, dec_deg=dec, hours=hours, **st))
            if src == "Cas A":
                figs[("earth", hours)] = (uv, im)

    fig, axs = plt.subplots(2, 4, figsize=(13, 6.4))
    sel = [((2.3e9, "rings16 static"), "rings16 static, 2.3 GHz"), ((2.3e9, "rings16 turntable 45 deg / 16 steps"), "rings16 + turntable"),
           (("earth", 0.0), "10 m field, 1420 MHz, static"), (("earth", 6.0), "10 m field, 6 h Earth rotation")]
    for c, (key, title) in enumerate(sel):
        uv, (ax_, b) = figs[key]
        axs[0, c].plot(np.r_[uv[:, 0], -uv[:, 0]], np.r_[uv[:, 1], -uv[:, 1]], ".", ms=1.2)
        axs[0, c].set(title=title, xlabel="u (wavelengths)", ylabel="v", aspect="equal")
        axs[1, c].imshow(np.abs(b), extent=[ax_[0], ax_[-1], ax_[0], ax_[-1]], origin="lower", vmin=0, vmax=1, cmap="viridis")
        axs[1, c].set(xlabel="l offset", ylabel="m offset", title="|point-spread function|")
    fig.tight_layout(); fig.savefig(OUT / "rotation_uv.png", dpi=120); plt.close(fig)
    return res


# =========================================================================================================
# 4  MIMO virtual arrays
# =========================================================================================================
def cut_stats(x, w, lam, u0=0.0, elem=None):
    """Azimuth cut of sum_i w_i exp(j k x_i (u - u0)); HPBW (deg) and peak sidelobe (dB) over all real u."""
    u = np.linspace(-1, 1, 20001)
    af = np.abs((w[None, :] * np.exp(2j * np.pi / lam * x[None, :] * (u[:, None] - u0))).sum(1))
    if elem is not None:
        af = af * elem(u)
    ip = np.argmin(np.abs(u - u0))
    while 0 < ip < len(u) - 1 and (af[ip + 1] > af[ip] or af[ip - 1] > af[ip]):
        ip += 1 if af[ip + 1] > af[ip] else -1
    pk = af[ip]
    lo = hi = ip
    while lo > 0 and af[lo] > pk / np.sqrt(2):
        lo -= 1
    while hi < len(u) - 1 and af[hi] > pk / np.sqrt(2):
        hi += 1
    nl, nh = lo, hi
    while nl > 0 and af[nl - 1] < af[nl]:
        nl -= 1
    while nh < len(u) - 1 and af[nh + 1] < af[nh]:
        nh += 1
    side = np.r_[af[:nl], af[nh + 1:]]
    return dict(hpbw_deg=float(np.degrees(np.arcsin(min(1, u[hi])) - np.arcsin(max(-1, u[lo])))),
                psl_db=float(20 * np.log10(side.max() / pk)) if len(side) else None,
                peak_at_deg=float(np.degrees(np.arcsin(u[ip])))), (u, af / pk)


def demo_mimo(quick=False):
    """Assumptions: point targets in the far field, isotropic elements (so sidelobes at wide angles are
    pessimistic), perfectly orthogonal transmit waveforms unless 'TDM' is named, a target that does not move
    during the coherent interval unless a speed is named, calibration perfect. Patterns are two-way
    (transmit x receive) azimuth cuts of the panel at 5.8 GHz and of the whip plate at 915 MHz."""
    res = {"panel_5800": [], "whips_915": [], "tdm": []}
    lam = C0 / 5.8e9
    g = 0.025
    rx = (np.arange(4) - 1.5) * g                                  # the 4x4 collapses to 4 x-positions, weight 4 each
    one = np.ones(4)
    curves = {}

    def add(name, x, w, note, **kw):
        st0, cv = cut_stats(x, w, lam, 0.0, **kw)
        st30, _ = cut_stats(x, w, lam, 0.5, **kw)
        res["panel_5800"].append(dict(case=name, virtual_x_positions=int(len(np.unique(np.round(x, 6)))), note=note,
                                      hpbw_deg=st0["hpbw_deg"], psl_db=st0["psl_db"], hpbw_at_30deg=st30["hpbw_deg"], psl_db_at_30deg=st30["psl_db"]))
        curves[name] = cv

    def virt(tx, rxp, wt=None, wr=None):
        wt = np.ones(len(tx)) if wt is None else wt
        wr = np.ones(len(rxp)) if wr is None else wr
        return (tx[:, None] + rxp[None, :]).ravel(), (wt[:, None] * wr[None, :]).ravel()

    add("a: 1 TX horn, 4x4 RX (rev B radar)", rx, one, "receive pattern only")
    add("b: 4x4 phased TX x 4x4 RX, beam scanned", *virt(rx, rx), "needs one dwell per transmit beam")
    add("c: 4 quadrant TX (orthogonal) x 16 RX", *virt(np.array([-g, g]), rx), "existing hardware; quadrant = 2x2 subarray at fixed steer",
        elem=lambda u: np.abs(np.cos(np.pi * g / lam * u)))
    add("d: 2x2 TX at 100 mm pitch x 16 RX -> 8x8", *virt(np.array([-2 * g, 2 * g]), rx), "4 extra TX antennas beside the panel")
    add("e: 4 TX in a row at 100 mm x 16 RX -> 16x4", *virt((np.arange(4) - 1.5) * 4 * g, rx), "4 extra TX antennas in a row")
    add("f: 2x2 TX at 100 mm x 4 quadrant RX (no switching)", *virt(np.array([-2 * g, 2 * g]), np.array([-g, g])),
        "16 virtual subarrays at 50 mm = 0.97 wavelength: grating lobes", elem=lambda u: np.abs(np.cos(np.pi * g / lam * u)))
    # time-division transmit: Doppler phase per slot on case e (4 transmitters fire in turn)
    Tslot = 1e-3
    tx_e = (np.arange(4) - 1.5) * 4 * g
    res["tdm_v_unambiguous_mps"] = float(lam / (4 * 4 * Tslot))
    for v in (0.0, 0.3, 1.0, 2.0):
        ph = np.exp(1j * 4 * np.pi * v * Tslot * np.arange(4) / lam)             # extra phase of transmitter m
        x, w = virt(tx_e, rx, wt=ph)
        st, _ = cut_stats(x, w, lam, 0.0)
        u = np.linspace(-1, 1, 20001)
        af = np.abs((w[None] * np.exp(2j * np.pi / lam * x[None] * u[:, None])).sum(1))
        res["tdm"].append(dict(v_mps=v, phase_step_deg=float(np.degrees(4 * np.pi * v * Tslot / lam)), peak_at_deg=st["peak_at_deg"],
                               loss_at_true_dir_db=float(20 * np.log10(af[10000] / 16)), psl_db=st["psl_db"]))
    # whips at 915 MHz: X0..X3 transmit (four whips of the outer ring, 90 deg apart), all 16 receive
    lam9 = C0 / 915e6
    az = np.radians(np.arange(0, 360, 0.05))

    def az_pattern(pv, az0=0.0):
        k = 2 * np.pi / lam9
        ph = pv[:, 0, None] * (np.cos(az)[None] - np.cos(az0)) + pv[:, 1, None] * (np.sin(az)[None] - np.sin(az0))
        af = np.abs(np.exp(1j * k * ph).sum(0)) / len(pv)
        i0 = int(np.argmin(np.abs(az - az0)))
        lo = hi = i0
        n = len(az)
        while af[lo % n] > 1 / np.sqrt(2):
            lo -= 1
        while af[hi % n] > 1 / np.sqrt(2):
            hi += 1
        nl, nh = lo, hi
        while af[(nl - 1) % n] < af[nl % n]:
            nl -= 1
        while af[(nh + 1) % n] < af[nh % n]:
            nh += 1
        idx = np.arange(nh + 1, nl + n) % n
        return dict(hpbw_deg=float((hi - lo) * 0.05), psl_db=float(20 * np.log10(af[idx].max())))

    r_out = np.linalg.norm(WHIPS16, axis=1)
    outer = WHIPS16[r_out > 0.2]
    tx4 = outer[np.argsort(np.arctan2(outer[:, 1], outer[:, 0]))][::2]
    cases = {"16 RX only": WHIPS16,
             "4 TX (outer ring, 90 deg apart) x 16 RX = 64 virtual": (tx4[:, None, :] + WHIPS16[None]).reshape(-1, 2),
             "16 phased TX x 16 RX (not buildable: 4 TX ports)": (WHIPS16[:, None, :] + WHIPS16[None]).reshape(-1, 2)}
    for name, pv in cases.items():
        sts = [az_pattern(pv, a) for a in np.radians([0, 11.25, 22.5, 33.75])]
        res["whips_915"].append(dict(case=name, hpbw_deg=float(np.mean([s["hpbw_deg"] for s in sts])),
                                     worst_psl_db=float(max(s["psl_db"] for s in sts))))

    fig, ax = plt.subplots(figsize=(8, 3.8))
    for name in list(curves)[:5]:
        u, a = curves[name]
        ax.plot(np.degrees(np.arcsin(u)), 20 * np.log10(np.maximum(a, 1e-4)), lw=1, label=name)
    ax.set(xlim=(-90, 90), ylim=(-45, 1), xlabel="azimuth, deg", ylabel="two-way pattern, dB", title="Panel at 5.8 GHz: phased against MIMO virtual arrays")
    ax.legend(fontsize=7, loc="upper right"); fig.tight_layout(); fig.savefig(OUT / "mimo_virtual.png", dpi=130); plt.close(fig)
    return res


# =========================================================================================================
# 5  switched (snapshot) covariance
# =========================================================================================================
def switching_schemes():
    """Channel c reaches the four whips of group c only (REDESIGN.md, section 3: one SP5T per group), so two
    whips of one group are never sampled together. Groups here: whip index mod 4 on each ring."""
    groups = [[g, g + 4, g + 8, g + 12] for g in range(4)]
    mul2 = [0, 2, 3, 1]                                             # multiply by the generator in GF(4)
    oa = [(a, b, a ^ b, a ^ mul2[b]) for a in range(4) for b in range(4)]      # orthogonal array OA(16, 4, 4, 2)
    mk = lambda rows: [[groups[c][r[c]] for c in range(4)] for r in rows]
    return {"all 16 at once (16 radios, hypothetical)": [list(range(16))],
            "4 slots (each whip once)": mk([(s, s, s, s) for s in range(4)]),
            "16 slots (orthogonal array)": mk(oa)}


def cov_fit(slots_R, slots, A, n):
    """Non-negative least squares fit of source powers on an azimuth grid (plus white noise) to the measured
    covariance entries only. A: (G, N) steering matrix."""
    num = np.zeros((n, n), complex)
    cnt = np.zeros((n, n))
    for R, s in zip(slots_R, slots):
        ix = np.ix_(s, s)
        num[ix] += R
        cnt[ix] += 1
    i, j = np.nonzero(np.triu(cnt, 1))
    Rm = num[i, j] / cnt[i, j]
    M = A[:, i] * A[:, j].conj()                                       # (G, P)
    dg = np.real(np.diag(num) / np.maximum(np.diag(cnt), 1))
    rows = np.vstack([np.column_stack([M.real.T, np.zeros(len(i))]), np.column_stack([M.imag.T, np.zeros(len(i))]),
                      np.column_stack([np.ones((n, A.shape[0])), np.ones(n)])])
    rhs = np.r_[Rm.real, Rm.imag, dg]
    x, _ = nnls(rows, rhs, maxiter=120)
    return x[:-1], len(i)


def peaks_circular(p, k):
    """Centroids of the k strongest connected clusters of a circular non-negative spectrum (grid index units)."""
    g = len(p)
    on = p > 1e-3 * p.max()
    if on.all() or not on.any():
        return np.array([])
    start = int(np.argmin(on))
    order = (np.arange(g) + start) % g
    cl, cur = [], []
    for ix in order:
        if on[ix]:
            cur.append(ix)
        elif cur:
            cl.append(cur); cur = []
    if cur:
        cl.append(cur)
    out = []
    for c in cl:
        w = p[c]
        ang = np.angle(np.sum(w * np.exp(2j * np.pi * np.array(c) / g)))
        out.append((w.sum(), ang))
    out.sort(reverse=True)
    return np.array([a for _, a in out[:k]])


def demo_switched(quick=False):
    """Assumptions: K uncorrelated equal-power sources at the horizon, at least 0.8 beamwidth apart (30 deg
    at 915 MHz, 12 deg at 2300 MHz), 10 dB per-whip SNR, 3 deg / 0.5 dB calibration residual, 4096 samples
    IN TOTAL however they are split into slots (so a 16-slot scheme has 256 per slot); the sources are
    stationary unless 'bursty': the observation is cut into 16 equal intervals and each source is on or off
    for a whole interval with probability 1/2 (on at least once); every scheme sees the same 16 intervals
    (one slot of the 16-slot scheme = one interval, one slot of the 4-slot scheme = four). A trial is a
    success when all K sources are found within a third of the minimum separation. The number of sources
    is given to the peak picker. Port-to-group assignment: whip index mod 4 (an assumption; REDESIGN.md
    section 4 lists the groups by port, not by position on the plate)."""
    rng = np.random.default_rng(5)
    trials = 12 if quick else 100
    schemes = switching_schemes()
    p = WHIPS16
    n = 16
    grid = np.radians(np.arange(0, 360, 1.0))
    res = {"pairs_measured": {}, "runs": []}
    for name, slots in schemes.items():
        pr = {(min(a, b), max(a, b)) for s in slots for a in s for b in s if a != b}
        res["pairs_measured"][name] = len(pr)
    for f in (0.915e9, 2.3e9):
        k = 2 * np.pi * f / C0
        A = np.exp(1j * k * (np.cos(grid)[:, None] * p[:, 0] + np.sin(grid)[:, None] * p[:, 1]))
        sep = 30.0 if f < 1e9 else 12.0
        for bursty in (False, True):
            for K in ((1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 20) if not bursty else (2, 3, 4)):
                if sep * K > 300:
                    continue
                ok = {nm: 0 for nm in schemes}
                for _ in range(trials):
                    # K azimuths at least sep apart, uniformly over all such sets (gaps = sep + Dirichlet share)
                    gaps = sep + (360.0 - sep * K) * rng.dirichlet(np.ones(K))
                    th = np.sort((rng.uniform(0, 360) + np.cumsum(gaps)) % 360.0)
                    tr = np.radians(th)
                    At = np.exp(1j * k * (np.cos(tr)[None, :] * p[:, 0, None] + np.sin(tr)[None, :] * p[:, 1, None]))
                    gcal = 10 ** (rng.normal(0, 0.5, n) / 20) * np.exp(1j * np.radians(rng.normal(0, 3.0, n)))
                    env = np.ones((K, 16))
                    if bursty:
                        env = (rng.random((K, 16)) < 0.5).astype(float)
                        for q in np.nonzero(~env.any(1))[0]:
                            env[q, rng.integers(16)] = 1.0
                        env *= np.sqrt(2)
                    for nm, slots in schemes.items():
                        T = 4096 // len(slots)
                        per = 16 // len(slots)                          # burst intervals per slot
                        Rs = []
                        for js, s in enumerate(slots):
                            S = cn(rng, (K, T)) * np.repeat(env[:, js * per:(js + 1) * per], T // per, axis=1)
                            X = gcal[s, None] * (At[s] @ S) + cn(rng, (len(s), T)) * 10 ** (-10 / 20)
                            Rs.append(X @ X.conj().T / T)
                        pw, _ = cov_fit(Rs, slots, A, n)
                        est = np.degrees(peaks_circular(pw, K))
                        if len(est) == K:
                            d = np.abs(wrap_deg(est[:, None] - th[None, :]))
                            ok[nm] += bool(np.all(d.min(0) < sep / 3) and np.all(d.min(1) < sep / 3))
                for nm in schemes:
                    res["runs"].append(dict(f_mhz=f / 1e6, bursty=bursty, K=K, scheme=nm, success=ok[nm] / trials))
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
    for a, f in zip(axs, (915.0, 2300.0)):
        for nm in schemes:
            rows = [r for r in res["runs"] if r["f_mhz"] == f and not r["bursty"] and r["scheme"] == nm]
            a.plot([r["K"] for r in rows], [r["success"] for r in rows], "o-", label=f"{nm}: {res['pairs_measured'][nm]} pairs")
        a.set(title=f"{f:.0f} MHz, two rings of 8 whips", xlabel="number of sources")
    axs[0].set_ylabel("all sources found"); axs[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(OUT / "switched_cov.png", dpi=130); plt.close(fig)
    return res


# =========================================================================================================
# 6  MUSIC resolution threshold
# =========================================================================================================
def demo_threshold(quick=False):
    """Assumptions: two uncorrelated equal-power sources at the horizon, 915 MHz, the 16-whip two-ring layout,
    all 16 whips sampled together (256 samples), number of sources known; calibration residual per whip of
    0 / 1 / 3 deg (with 0 / 0.17 / 0.5 dB). Resolved = MUSIC shows two peaks, each within half the separation
    of a source."""
    rng = np.random.default_rng(6)
    trials = 60 if quick else 300
    p, f = WHIPS16, 915e6
    k = 2 * np.pi * f / C0
    az = np.radians(np.arange(-60, 60.001, 0.1))
    A = np.exp(1j * k * (np.cos(az)[:, None] * p[:, 0] + np.sin(az)[:, None] * p[:, 1]))
    bp = np.abs(A @ A[len(az) // 2].conj()) / 16
    hp = az[bp >= 1 / np.sqrt(2)]
    hpbw = float(np.degrees(hp.max() - hp.min()))
    res = {"hpbw_deg": hpbw, "runs": []}
    for frac in (0.1, 0.25, 0.5):
        sep = frac * hpbw
        th = np.radians([-sep / 2, sep / 2])
        At = np.exp(1j * k * (np.cos(th)[None] * p[:, 0, None] + np.sin(th)[None] * p[:, 1, None]))
        for cal in (0.0, 1.0, 3.0):
            for snr in (-10, 0, 10, 20, 30, 40):
                ok = 0
                for _ in range(trials):
                    g = 10 ** (rng.normal(0, cal / 6, 16) / 20) * np.exp(1j * np.radians(rng.normal(0, cal, 16)))
                    X = g[:, None] * (At @ cn(rng, (2, 256))) + cn(rng, (16, 256)) * 10 ** (-snr / 20)
                    w, V = np.linalg.eigh(X @ X.conj().T / 256)
                    sp = 1 / np.sum(np.abs(A.conj() @ V[:, :14]) ** 2, axis=1)
                    pk = np.nonzero((sp[1:-1] > sp[:-2]) & (sp[1:-1] > sp[2:]))[0] + 1
                    pk = pk[np.argsort(-sp[pk])[:2]]
                    if len(pk) == 2:
                        e = np.sort(np.degrees(az[pk])) - np.array([-sep / 2, sep / 2])
                        ok += bool(np.all(np.abs(e) < sep / 2))
                res["runs"].append(dict(sep_beamwidths=frac, sep_deg=sep, cal_deg=cal, snr_db=snr, p_resolved=ok / trials))
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
    for a, frac in zip(axs, (0.1, 0.25, 0.5)):
        for cal in (0.0, 1.0, 3.0):
            rows = [r for r in res["runs"] if r["sep_beamwidths"] == frac and r["cal_deg"] == cal]
            a.plot([r["snr_db"] for r in rows], [r["p_resolved"] for r in rows], "o-", label=f"calibration residual {cal:.0f} deg")
        a.set(title=f"separation {frac} beamwidth ({frac * hpbw:.1f} deg)", xlabel="SNR per whip, dB")
    axs[0].set_ylabel("P(two sources resolved)"); axs[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(OUT / "threshold.png", dpi=130); plt.close(fig)
    return res

# =========================================================================================================
# 7  sensitivity at the radio-astronomy frequencies
# =========================================================================================================
# Perley & Butler 2017 (arXiv:1609.05940, table of fitted spectra): log10 S[Jy] = sum a_i log10(f/GHz)^i
PB17 = {"Cyg A": ((3.3498, -1.0022, -0.225, 0.023, 0.043), (0.05, 12.0)),
        "Cas A (2016)": ((3.3584, -0.7518, -0.035, -0.071), (0.2, 4.0)),
        "Tau A": ((2.9516, -0.217, -0.047, -0.067), (0.05, 4.0)),
        "Vir A": ((2.4466, -0.8116, -0.048), (0.05, 3.0))}
# Learmonth quiet-Sun flux, 15 June 2021 (near solar minimum), solar flux units (1 sfu = 1e4 Jy):
# https://listserver.ips.gov.au/pipermail/ips-iflux-dailyvalues/2021-June/006577.html
QUIET_SUN_SFU = {410: 29, 610: 42, 1415: 61, 2695: 67, 4995: 134}
KB = 1.380649e-23


def demo_astro(quick=False):
    """Assumptions (ESTIMATES, not measurements): one polarisation; a whip over the plate has 2 dBi towards
    a source 30-60 deg up; a panel element has the lattice-cell gain 4 pi (25 mm)^2 / lambda^2; system
    temperature = receiver (REDESIGN.md noise figures) + 150 K of ground and sky for a whip (it sees the
    ground as well as it sees the sky) or + 50 K for the panel pointed up; 20 MHz of bandwidth, except where
    the protected band is narrower; cross-correlations only, point sources, gains stable for the whole
    integration. Interferometer noise: dS = SEFD / sqrt(Npairs_simultaneous * 2 * B * t), SEFD = 2 k Tsys / Ae,
    Ae = G lambda^2 / 4 pi. 'switched' = 4 of the 16 elements at a time (6 of the 120 pairs), any schedule."""
    rows = []
    bands = [  # f MHz, what, usable bandwidth MHz, front end, NF dB
        (408.0, "406.1-410 continuum", 3.9, "whip", 1.1), (611.0, "608-614 continuum", 6.0, "whip", 1.2),
        (1420.4, "H I, 1400-1427", 20.0, "whip", 1.2), (1665.4, "OH, 1660-1670", 10.0, "whip", 1.3),
        (2695.0, "2690-2700 continuum", 10.0, "panel", 2.0), (4830.0, "H2CO, 4800-5000", 20.0, "panel", 3.0),
        (4995.0, "4990-5000 continuum", 10.0, "panel", 3.0)]
    sun_f = np.array(sorted(QUIET_SUN_SFU))
    sun_s = np.array([QUIET_SUN_SFU[k] for k in sun_f]) * 1e4
    for f, what, bw, fe, nf in bands:
        lam = C0 / (f * 1e6)
        trx = 290 * (10 ** (nf / 10) - 1)
        if fe == "whip":
            g, tsys = 10 ** 0.2, trx + 150
        else:
            g, tsys = 4 * np.pi * 0.025 ** 2 / lam ** 2, trx + 50
        ae = g * lam ** 2 / (4 * np.pi)
        sefd = 2 * KB * tsys / ae / 1e-26                                    # Jy, one element
        ds16 = sefd / np.sqrt(120 * 2 * bw * 1e6)                             # 1 s, all 16 at once
        ds4 = sefd / np.sqrt(6 * 2 * bw * 1e6)                                # 1 s, 4 at a time
        row = dict(f_mhz=f, band=what, bw_mhz=bw, front_end=fe, tsys_k=float(tsys), ae_cm2=float(ae * 1e4),
                   sefd_element_jy=float(sefd), ds_1s_all16_jy=float(ds16), ds_1s_switched_jy=float(ds4), sources={})
        x = np.log10(f / 1e3)
        src = {k: (10 ** sum(a * x ** i for i, a in enumerate(c)), lo <= f / 1e3 <= hi) for k, (c, (lo, hi)) in PB17.items()}
        src["quiet Sun"] = (float(np.exp(np.interp(np.log(f), np.log(sun_f), np.log(sun_s)))), True)
        for k, (sj, valid) in src.items():
            row["sources"][k] = dict(flux_jy=float(sj), within_fit_range=bool(valid),
                                     t_10sigma_all16_s=float((10 * ds16 / sj) ** 2), t_10sigma_switched_s=float((10 * ds4 / sj) ** 2))
        # total-power line radiometry on one element: dT = Tsys / sqrt(dnu * t), 10 kHz channels
        row["dT_mK_10kHz_1000s_one_element"] = float(1e3 * tsys / np.sqrt(1e4 * 1000))
        row["dT_mK_10kHz_1000s_16_elements_added_incoherently"] = float(1e3 * tsys / np.sqrt(1e4 * 1000 * 16))
        rows.append(row)
    res = {"bands": rows}
    # what a longer baseline costs: delay across the band, fringe rate, field of view of one correlator channel
    geo = []
    for D in (0.45, 3.0, 10.0, 30.0, 100.0):
        for f, bw in ((408e6, 3.9e6), (1420.4e6, 20e6)):
            lam = C0 / f
            geo.append(dict(baseline_m=D, f_mhz=f / 1e6, resolution_deg=float(np.degrees(lam / D)),
                            delay_ns_at_horizon=float(D / C0 * 1e9), delay_times_bandwidth=float(D / C0 * bw),
                            channel_khz_for_1pct_loss=float(0.078 * C0 / D / 1e3),   # sinc(x) = 0.99 at x = 0.078
                            fringe_rate_mhz_max=float(7.2921e-5 * D / lam * 1e3),
                            sun_visibility=float(abs(np.sinc(0)) if D / lam * np.radians(0.53) < 1e-6 else
                                                 abs(2 * _j1(np.pi * D / lam * np.radians(0.53)) / (np.pi * D / lam * np.radians(0.53))))))
    res["baseline_costs"] = geo
    return res


def _j1(x):
    """Bessel J1 by its integral (no scipy)."""
    t = np.linspace(0, np.pi, 2001)
    y = np.cos(t - x * np.sin(t))
    return float((y.sum() - 0.5 * (y[0] + y[-1])) * (t[1] - t[0]) / np.pi)


DEMOS = dict(coarray=demo_coarray, freqdiv=demo_freqdiv, rotation=demo_rotation, mimo=demo_mimo,
             switched=demo_switched, threshold=demo_threshold, astro=demo_astro)


def main(argv):
    quick = "--quick" in argv
    names = [a for a in argv if a in DEMOS] or list(DEMOS)
    OUT.mkdir(exist_ok=True)
    path = OUT / ("results_quick.json" if quick else "results.json")
    res = json.loads(path.read_text()) if path.exists() else {}
    for nm in names:
        t0 = time.time()
        res[nm] = DEMOS[nm](quick)
        res[nm]["assumptions"] = " ".join(DEMOS[nm].__doc__.split())
        print(f"{nm}: {time.time() - t0:.0f} s", flush=True)
        path.write_text(json.dumps(res, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))


if __name__ == "__main__":
    main(sys.argv[1:])
