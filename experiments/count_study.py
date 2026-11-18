"""How many antennas? Direction finding, beams, resolvable sources, fix time and cost against the number of
elements N = 4 ... 64, for planar layouts on the 500 mm plate and 3-D layouts in a 500 x 500 x 300 mm volume.

Conventions are those of experiments/evolve_layout.py (imported, not copied): isotropic elements, a
calibration residual per element (phase sigma_cal, 0.5 dB amplitude), K = 64 samples, per-element SNR,
Bartlett beamscan, gross error = azimuth error over 10 deg, conditions "elevated" (3 deg, 20 dB, sources
0-30 deg up) and "stress" (5 deg, 0 dB, 0-30 deg up). New here:

  geometry    positions are (N, 3); the planar plate is z = 0; the volume is the plate's annulus
              (60 mm <= r <= 240 mm) times |z| <= 150 mm. Minimum spacing 75 mm, except planar N = 32, 48, 64
              (60, 48, 42 mm: more than 24 whips do not fit on the plate at 75 mm).
  estimators  "az"   = evolve_layout's azimuth-only beamscan, which assumes the source is on the horizon;
              "azel" = a two-dimensional beamscan over azimuth and elevation (software instead of geometry).
  sharing     "full"   = N coherent channels, one slot;
              "groups" = 4 channels, N/4 slots of disjoint groups, slots combined without phase (the power
                         spectra of the 4-element sub-arrays are added): valid for any signal;
              "ref"    = 4 channels, one stays on a reference element and three rove, ceil((N-1)/3) slots,
                         each element's phase taken against the reference and stitched into one vector.
  theory      calibration-limited Cramer-Rao-type bound, aperture directivity, degrees of freedom, slot counts.

    python experiments/count_study.py            # about 40 minutes with OMP_NUM_THREADS=3
    python experiments/count_study.py --quick    # a smoke run
    python experiments/count_study.py --report   # tables and figures from count/results.json
    python experiments/count_study.py --extras   # theory checks, equal-time experiment, fix time, astronomy, cost
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import evolve_layout as EL  # noqa: E402
from array16.evolve import climb, random_search, reselect  # noqa: E402

C0 = EL.C0
R_OUT, R_IN, Z_HALF = EL.R_OUT, EL.R_IN, 0.150
OUT = HERE / "count"
NS = [4, 6, 8, 12, 16, 24, 32, 48, 64]
ASTRO = {"408 MHz continuum": 0.408e9, "611 MHz continuum": 0.611e9, "H I 1420.4 MHz": 1.4204e9, "OH 1665.4 MHz": 1.6654e9}
FREQS = np.sort(np.concatenate([np.arange(0.4e9, 2.3001e9, 0.1e9), list(ASTRO.values())]))
CONDS = dict(elevated=dict(sigma_cal=3.0, snr_db=20.0, el_max=30.0, seed=4242),
             stress=dict(sigma_cal=5.0, snr_db=0.0, el_max=30.0, seed=777))
COST_PORT, COST_CELL = 11.25, 14.56       # REDESIGN.md section 10: front ends 180 / 16 ports, cells 233 / 16


def s_min(n: int, dim: int) -> float:
    return 0.075 if dim == 3 or n <= 24 else {32: 0.060, 48: 0.048, 64: 0.042}[n]


# ---- geometry ---------------------------------------------------------------------------------------------
def repair3(x: np.ndarray, dim: int, smin: float) -> np.ndarray:
    """evolve_layout.repair for (N, dim) designs: onto the annulus (and |z| <= Z_HALF), pairs pushed apart."""
    p = x.reshape(-1, dim).copy()
    for _ in range(80):
        r = np.hypot(p[:, 0], p[:, 1])
        p[:, :2] *= (np.clip(r, R_IN, R_OUT) / np.maximum(r, 1e-12))[:, None]
        if dim == 3:
            p[:, 2] = np.clip(p[:, 2], -Z_HALF, Z_HALF)
        d = p[:, None, :] - p[None, :, :]
        dist = np.linalg.norm(d, axis=2)
        np.fill_diagonal(dist, np.inf)
        i, j = np.nonzero(np.triu(dist < smin - 1e-9, 1))
        if not len(i):
            break
        for a, b in zip(i, j):
            v = p[a] - p[b]
            n = np.linalg.norm(v)
            v = v / n if n > 1e-12 else np.eye(dim)[0]
            push = (smin - n) / 2 + 1e-6
            p[a] += push * v
            p[b] -= push * v
    p = p[np.argsort(np.arctan2(p[:, 1], p[:, 0]), kind="stable")]
    return p.reshape(-1)


def violation3(p: np.ndarray, smin: float) -> float:
    r = np.hypot(p[:, 0], p[:, 1])
    d = np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)
    np.fill_diagonal(d, np.inf)
    v = np.sum(np.maximum(r - R_OUT, 0)) + np.sum(np.maximum(R_IN - r, 0)) + np.sum(np.triu(np.maximum(smin - d, 0), 1))
    if p.shape[1] == 3:
        v += np.sum(np.maximum(np.abs(p[:, 2]) - Z_HALF, 0))
    return float(v)


def to3(p: np.ndarray) -> np.ndarray:
    return p if p.shape[1] == 3 else np.column_stack([p, np.zeros(len(p))])


def min_spacing(p: np.ndarray) -> float:
    d = np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)
    np.fill_diagonal(d, np.inf)
    return float(d.min())


def rings(n: int, smin: float) -> np.ndarray:
    """Concentric rings filled from the rim inwards. N <= 16: rev B's two rings (225 and 112 mm, N/2 each).
    More: the fewest rings (radii spread from 240 to 80 mm, elements in proportion to radius) that keep the
    spacing."""
    if 8 <= n <= 16:
        return np.vstack([EL.ring(n // 2, 0.225), EL.ring(n - n // 2, 0.112, np.pi / (n // 2))])
    if n < 8:
        return EL.ring(n, 0.225)
    for m in range(2, 9):
        radii = np.linspace(0.240, 0.080, m)
        cnt = np.floor(n * radii / radii.sum()).astype(int)
        for i in range(n - cnt.sum()):
            cnt[i % m] += 1
        p = np.vstack([EL.ring(c, r, np.pi / c * (i % 2)) for i, (c, r) in enumerate(zip(cnt, radii))])
        if min_spacing(p) >= smin - 1e-9:
            return p
    return repair3(p.reshape(-1), 2, smin).reshape(-1, 2)


def stack(n: int) -> np.ndarray:
    """Rings of radius 225 mm stacked between z = -150 and +150 mm, each turned against the one below."""
    layers = max(2, math.ceil(n / 18))
    cnt = [n // layers + (1 if i < n % layers else 0) for i in range(layers)]
    z = np.linspace(-Z_HALF, Z_HALF, layers)
    return np.vstack([np.column_stack([EL.ring(c, 0.225, np.pi / c * (i % 2)), np.full(c, z[i])]) for i, c in enumerate(cnt)])


def cone(n: int, smin: float) -> np.ndarray:
    """The planar rings with height in proportion to the distance from the rim (rim at -150 mm, innermost ring
    at +150 mm): a frustum. For N < 8 a helix on the 225 mm cylinder."""
    if n < 8:
        b = 2 * np.pi * np.arange(n) / n
        return np.column_stack([0.225 * np.cos(b), 0.225 * np.sin(b), np.linspace(-Z_HALF, Z_HALF, n)])
    p = rings(n, smin)
    r = np.hypot(p[:, 0], p[:, 1])
    return np.column_stack([p, -Z_HALF + 2 * Z_HALF * (r.max() - r) / max(r.max() - r.min(), 1e-9)])


def helix(n: int) -> np.ndarray:
    """N elements on the 225 mm cylinder, heights spread evenly over 300 mm in golden-ratio order of azimuth."""
    b = 2 * np.pi * np.arange(n) / n
    order = np.argsort((np.arange(n) * 0.6180339887) % 1.0)
    z = np.empty(n)
    z[order] = np.linspace(-Z_HALF, Z_HALF, n)
    return np.column_stack([0.225 * np.cos(b), 0.225 * np.sin(b), z])


# ---- search model (evolve_layout.make_score for (N, dim) designs; the "az" mode with dim 2 is identical) ----
def _refine_lin(pw: np.ndarray, i: np.ndarray) -> np.ndarray:
    """Parabolic offset (in grid steps) along the last axis without wrap-around; 0 at the ends."""
    g = pw.shape[-1]
    y0 = np.take_along_axis(pw, np.clip(i - 1, 0, g - 1)[..., None], -1)[..., 0]
    y1 = np.take_along_axis(pw, i[..., None], -1)[..., 0]
    y2 = np.take_along_axis(pw, np.clip(i + 1, 0, g - 1)[..., None], -1)[..., 0]
    den = y0 - 2 * y1 + y2
    ok = (np.abs(den) > 1e-30) & (i > 0) & (i < g - 1)
    return np.where(ok, np.clip(0.5 * (y0 - y2) / np.where(ok, den, 1.0), -0.5, 0.5), 0.0)


def make_score3(n: int, dim: int, smin: float, mode: str, *, n_freq=16, trials=96, sigma_cal=4.0, snr_db=10.0,
                k_samples=64, el_max=EL.EL_MAX, seed=0):
    snr_eff = snr_db + 10 * np.log10(k_samples)
    edges = np.geomspace(*EL.BAND, n_freq + 1)
    az_g = np.radians(np.arange(0.0, 360.0, 3.0))
    el_g = np.radians(np.arange(-8.0 if dim == 3 else 0.0, 40.01, 4.0))
    ua = np.stack([np.cos(az_g), np.sin(az_g)], 1)                                           # (Ga, 2)

    def score(X: np.ndarray, gen: int) -> np.ndarray:
        rng = np.random.default_rng([seed, gen])
        f = edges[:-1] * (edges[1:] / edges[:-1]) ** rng.random(n_freq)
        k = 2 * np.pi * f / C0
        az = rng.uniform(0, 2 * np.pi, (n_freq, trials))
        el = np.radians(rng.uniform(0, el_max, (n_freq, trials)))
        ce, se = np.cos(el), np.sin(el)
        cal = (10 ** (rng.normal(0, 0.5, (n_freq, trials, n)) / 20)
               * np.exp(1j * np.radians(rng.normal(0, sigma_cal, (n_freq, trials, n)))))
        noise = ((rng.normal(size=(n_freq, trials, n)) + 1j * rng.normal(size=(n_freq, trials, n)))
                 / np.sqrt(2) * 10 ** (-snr_eff / 20))
        out = np.empty(len(X))
        for m, x in enumerate(X):
            p = x.reshape(n, dim)
            proj_t = (p[:, 0] * np.cos(az)[..., None] + p[:, 1] * np.sin(az)[..., None]) * ce[..., None]
            if dim == 3:
                proj_t = proj_t + p[:, 2] * se[..., None]
            y = np.exp(1j * k[:, None, None] * proj_t) * cal + noise                           # (F, T, N)
            if mode == "az":
                proj_g = p[:, 0] * np.cos(EL.GRID)[:, None] + p[:, 1] * np.sin(EL.GRID)[:, None]
                A = np.exp(-1j * k[:, None, None] * proj_g[None])
                pw = (np.abs(A @ y.transpose(0, 2, 1)) ** 2).transpose(0, 2, 1)
                est = EL._refine(pw, np.argmax(pw, axis=-1))
                err = np.degrees(np.angle(np.exp(1j * (est - az))))
                loss = np.sqrt(np.mean(np.minimum(np.abs(err), 30.0) ** 2, axis=1))
            else:
                h = ua @ p[:, :2].T                                                           # (Ga, N)
                proj_g = h[:, None, :] * np.cos(el_g)[None, :, None]
                if dim == 3:
                    proj_g = proj_g + (np.sin(el_g)[:, None] * p[:, 2])[None]
                A = np.exp(-1j * k[:, None, None] * proj_g.reshape(1, -1, n)).astype(np.complex64)
                pw = np.abs(A @ y.transpose(0, 2, 1).astype(np.complex64)) ** 2               # (F, G, T)
                pw = pw.transpose(0, 2, 1).reshape(n_freq, trials, len(az_g), len(el_g)).astype(float)
                idx = np.argmax(pw.reshape(n_freq, trials, -1), axis=-1)
                ia, ie = np.unravel_index(idx, (len(az_g), len(el_g)))
                pa = np.take_along_axis(pw, ie[..., None, None], 3)[..., 0]                    # (F, T, Ga)
                pe = np.take_along_axis(pw, ia[..., None, None], 2)[..., 0, :]                 # (F, T, Ge)
                g = len(az_g)
                y0 = np.take_along_axis(pa, ((ia - 1) % g)[..., None], -1)[..., 0]
                y1 = np.take_along_axis(pa, ia[..., None], -1)[..., 0]
                y2 = np.take_along_axis(pa, ((ia + 1) % g)[..., None], -1)[..., 0]
                den = y0 - 2 * y1 + y2
                off = np.where(np.abs(den) > 1e-30, 0.5 * (y0 - y2) / np.where(np.abs(den) > 1e-30, den, 1.0), 0.0)
                az_e = np.radians(3.0) * (ia + np.clip(off, -0.5, 0.5))
                el_e = el_g[0] + np.radians(4.0) * (ie + _refine_lin(pe, ie))
                e_az = np.degrees(np.angle(np.exp(1j * (az_e - az))))
                loss = np.sqrt(np.mean(np.minimum(np.abs(e_az), 30.0) ** 2, axis=1))
                if dim == 3:
                    e_el = np.degrees(el_e - el)
                    loss = loss + 0.5 * np.sqrt(np.mean(np.minimum(np.abs(e_el), 30.0) ** 2, axis=1))
            out[m] = loss.max() + 0.25 * loss.mean() + 100.0 * violation3(p, smin)
        return out

    return score


# ---- check model: evolve_layout.check's signal model with height, a 2-D estimator and channel sharing ------
def slots_for(n: int, sharing: str, m: int = 4):
    if sharing == "full":
        return [np.arange(n)]
    if sharing == "groups":                       # elements are sorted by azimuth: each group spans the plate
        g = n // m
        return [np.arange(n)[i::g] for i in range(g)]
    rov = np.arange(1, n)                         # "ref": element 0 on one channel all the time
    return [np.concatenate([[0], rov[i:i + m - 1]]) for i in range(0, n - 1, m - 1)]


def check3(p: np.ndarray, freqs, *, sharing="full", sigma_cal=3.0, snr_db=20.0, el_max=30.0, trials=600,
           k_samples=64, seed=12345, chunk=None):
    """Azimuth-elevation beamscan (coarse 3 x 4 deg, then 0.25 x 0.5 deg around the peak, then a parabola).
    Returns per frequency: median and 95th percentile azimuth error, azimuth gross errors (%), median and
    95th percentile elevation error, elevation gross errors (%, over 10 deg)."""
    rng = np.random.default_rng(seed)
    p = to3(np.asarray(p, float))
    n = len(p)
    flat = bool(np.allclose(p[:, 2], 0))
    az_g = np.radians(np.arange(0.0, 360.0, 3.0))
    el_g = np.radians(np.arange(0.0 if flat else -8.0, 44.01, 4.0))
    AZ, ELg = np.meshgrid(az_g, el_g, indexing="ij")
    U = np.stack([np.cos(ELg) * np.cos(AZ), np.cos(ELg) * np.sin(AZ), np.sin(ELg)], -1).reshape(-1, 3)
    d_az = np.radians(np.arange(-2.0, 2.001, 0.25))
    d_el = np.radians(np.arange(-3.0, 3.001, 0.5))
    slots = slots_for(n, sharing)
    chunk = chunk or max(8, int(4e7 // (len(U) * n)))
    res = []
    for f in freqs:
        k = 2 * np.pi * f / C0
        th = rng.uniform(0, 2 * np.pi, trials)
        el = np.radians(rng.uniform(0, el_max, trials))
        u = np.stack([np.cos(el) * np.cos(th), np.cos(el) * np.sin(th), np.sin(el)], 1)
        a = np.exp(1j * k * (u @ p.T))                                                        # (T, N)
        g = 10 ** (rng.normal(0, 0.5, (trials, n)) / 20) * np.exp(1j * np.radians(rng.normal(0, sigma_cal, (trials, n))))
        ag = a * g
        R = np.zeros((trials, n, n), complex) if sharing != "ref" else None
        yv = np.zeros((trials, n), complex)
        nref = 0
        for idx in slots:
            s = (rng.normal(size=(trials, k_samples)) + 1j * rng.normal(size=(trials, k_samples))) / np.sqrt(2)
            w = ((rng.normal(size=(trials, len(idx), k_samples)) + 1j * rng.normal(size=(trials, len(idx), k_samples)))
                 / np.sqrt(2) * 10 ** (-snr_db / 20))
            x = ag[:, idx, None] * s[:, None, :] + w
            if sharing == "ref":
                pref = np.sqrt(np.mean(np.abs(x[:, 0]) ** 2, axis=1))                         # (T,)
                yv[:, idx[1:]] = np.mean(x[:, 1:] * x[:, :1].conj(), axis=2) / pref[:, None]
                yv[:, 0] += pref
                nref += 1
            else:
                R[:, idx[:, None], idx[None, :]] = x @ x.conj().transpose(0, 2, 1) / k_samples
        if sharing == "ref":
            yv[:, 0] /= nref

        def power(A, sl):                              # A (G, N) or (t, G, N) -> (t, G)
            Ac = A.conj().astype(np.complex64)
            if sharing == "ref":
                return np.abs(np.einsum('...gn,tn->tg', Ac, yv[sl].astype(np.complex64)) if A.ndim == 2
                              else np.einsum('tgn,tn->tg', Ac, yv[sl].astype(np.complex64))) ** 2
            M = np.matmul(Ac if A.ndim == 3 else Ac[None], R[sl].astype(np.complex64))        # (t, G, N)
            return np.real(np.sum(M * (A if A.ndim == 3 else A[None]), axis=-1))

        Ac = np.exp(1j * k * (U @ p.T))
        az_e = np.empty(trials)
        el_e = np.empty(trials)
        for c0 in range(0, trials, chunk):
            sl = slice(c0, min(c0 + chunk, trials))
            t = sl.stop - sl.start
            i0 = np.argmax(power(Ac, sl), axis=1)
            a0, e0 = AZ.reshape(-1)[i0], ELg.reshape(-1)[i0]
            fa = (a0[:, None, None] + d_az[None, :, None]) + 0 * d_el[None, None, :]
            fe = (e0[:, None, None] + d_el[None, None, :]) + 0 * d_az[None, :, None]
            Uf = np.stack([np.cos(fe) * np.cos(fa), np.cos(fe) * np.sin(fa), np.sin(fe)], -1).reshape(t, -1, 3)
            pf = power(np.exp(1j * k * (Uf @ p.T)), sl).reshape(t, len(d_az), len(d_el)).astype(float)
            ii = np.argmax(pf.reshape(t, -1), axis=1)
            ia, ie = np.unravel_index(ii, (len(d_az), len(d_el)))
            pa = pf[np.arange(t), :, ie]
            pe = pf[np.arange(t), ia, :]
            az_e[sl] = a0 + d_az[0] + np.radians(0.25) * (ia + _refine_lin(pa, ia))
            el_e[sl] = e0 + d_el[0] + np.radians(0.5) * (ie + _refine_lin(pe, ie))
        ea = np.abs(np.degrees(np.angle(np.exp(1j * (az_e - th)))))
        ee = np.abs(np.degrees((np.abs(el_e) if flat else el_e) - el))
        res.append((np.median(ea), np.percentile(ea, 95), np.mean(ea > 10) * 100,
                    np.median(ee), np.percentile(ee, 95), np.mean(ee > 10) * 100))
    return np.array(res, float)


# ---- theory -----------------------------------------------------------------------------------------------
def crb(p: np.ndarray, f: float, *, sigma_cal=3.0, snr_db=20.0, k_samples=64, el_max=30.0, joint=True):
    """Calibration-limited bound. A per-element phase error phi_n (variance s2 = sigma_cal^2 + 1/(2 K rho),
    the second term the thermal phase noise of K samples at SNR rho) moves the estimate by
    (D^T D)^-1 D^T phi, with D the (N, 1 or 2) matrix of d(phase_n)/d(az, el), its column means removed (the
    common phase is not known). Covariance s2 (D^T D)^-1. Returns the RMS over azimuth and over elevation
    2 deg ... el_max of the azimuth and elevation standard deviations, in degrees."""
    p = to3(p)
    k = 2 * np.pi * f / C0
    s2 = np.radians(sigma_cal) ** 2 + 1 / (2 * k_samples * 10 ** (snr_db / 10))
    va, ve = [], []
    for el in np.radians(np.linspace(2.0, el_max, 8)):
        for az in np.radians(np.arange(0, 360, 15.0)):
            da = k * np.cos(el) * (-p[:, 0] * np.sin(az) + p[:, 1] * np.cos(az))
            de = k * (-np.sin(el) * (p[:, 0] * np.cos(az) + p[:, 1] * np.sin(az)) + np.cos(el) * p[:, 2])
            D = np.column_stack([da - da.mean(), de - de.mean()])
            if joint:
                Cv = s2 * np.linalg.inv(D.T @ D + 1e-12 * np.eye(2))
                va.append(Cv[0, 0]); ve.append(Cv[1, 1])
            else:
                va.append(s2 / (D[:, 0] @ D[:, 0])); ve.append(np.nan)
    return float(np.degrees(np.sqrt(np.mean(va)))), float(np.degrees(np.sqrt(np.mean(ve))))


def directivity(p: np.ndarray, f: float, n_az: int = 12) -> float:
    """Directivity (dBi) of the uniformly weighted beam of isotropic elements steered to the horizon, averaged
    over azimuth: D = |w^H a|^2 / (w^H S w), S_mn = sin(k d_mn) / (k d_mn) (the average of a a^H over the
    sphere)."""
    p = to3(p)
    k = 2 * np.pi * f / C0
    d = np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)
    S = np.sinc(k * d / np.pi)
    out = []
    for az in 2 * np.pi * (np.arange(n_az) + 0.37) / n_az:
        w = np.exp(1j * k * (p[:, 0] * np.cos(az) + p[:, 1] * np.sin(az)))
        out.append(len(p) ** 2 / np.real(w.conj() @ S @ w))
    return float(10 * np.log10(np.mean(out)))


def sidelobes(p: np.ndarray, f: float, n_az: int = 12):
    """Peak sidelobe (dB below the beam) of the uniformly weighted beam steered to the horizon: the highest
    local maximum other than the beam itself, over the whole sphere and along the horizon; mean and worst
    over n_az steering azimuths. Also the mean -3 dB azimuth beamwidth (deg)."""
    p = to3(p)
    n = len(p)
    k = 2 * np.pi * f / C0
    az = np.radians(np.arange(0, 360, 2.0))
    el = np.radians(np.arange(-88, 88.01, 2.0))
    AZ, ELg = np.meshgrid(az, el, indexing="ij")
    U = np.stack([np.cos(ELg) * np.cos(AZ), np.cos(ELg) * np.sin(AZ), np.sin(ELg)], -1).reshape(-1, 3)
    A = np.exp(1j * k * (U @ p.T))
    steer = np.radians(2.0 * np.round(360 / n_az * (np.arange(n_az) + 0.37) / 2.0))
    W = np.exp(1j * k * (np.cos(steer)[:, None] * p[:, 0] + np.sin(steer)[:, None] * p[:, 1]))  # (S, N)
    P = (np.abs(A @ W.conj().T) ** 2 / n ** 2).T.reshape(n_az, len(az), len(el))
    i_hor = len(el) // 2
    sph, hor, bw = [], [], []
    for s in range(n_az):
        q = P[s]
        pad = np.pad(q, ((0, 0), (1, 1)), constant_values=-1.0)
        nb = [np.roll(pad, (da, de), (0, 1)) for da in (-1, 0, 1) for de in (-1, 0, 1) if (da, de) != (0, 0)]
        peak = np.all(pad[None] >= np.array(nb), axis=0)[:, 1:-1]
        v = np.sort(q[peak])[::-1]
        sph.append(v[1] if len(v) > 1 else 0.0)
        h = q[:, i_hor]
        pk = (h >= np.roll(h, 1)) & (h >= np.roll(h, -1))
        vh = np.sort(h[pk])[::-1]
        hor.append(vh[1] if len(vh) > 1 else 0.0)
        i0 = int(np.argmax(h))
        hh = np.roll(h, len(az) // 2 - i0)
        c = len(az) // 2
        lo = c - np.argmax(hh[c::-1] < 0.5) if np.any(hh[c::-1] < 0.5) else 0
        hi = c + np.argmax(hh[c:] < 0.5) if np.any(hh[c:] < 0.5) else len(az)
        bw.append((hi - lo) * 2.0)
    db = lambda x: float(10 * np.log10(max(x, 1e-6)))
    return dict(sphere_mean=db(np.mean(sph)), sphere_worst=db(np.max(sph)), horizon_mean=db(np.mean(hor)),
                horizon_worst=db(np.max(hor)), beamwidth=float(np.mean(bw)))


def music_capacity(p: np.ndarray, f: float, counts, *, trials=60, k_samples=256, snr_db=20.0, sigma_cal=3.0, seed=99):
    """How many equal, uncorrelated sources on the horizon MUSIC separates with N coherent channels. Sources
    evenly spread in azimuth with a jitter of +/-0.2 of their spacing; success = S peaks, every one within
    min(5 deg, a third of the spacing) of its source. Returns {S: success rate}."""
    rng = np.random.default_rng(seed)
    p = to3(p)
    n = len(p)
    k = 2 * np.pi * f / C0
    grid = np.radians(np.arange(0, 360, 0.25))
    A = np.exp(1j * k * (np.cos(grid)[:, None] * p[:, 0] + np.sin(grid)[:, None] * p[:, 1]))     # (G, N)
    out = {}
    for S in counts:
        if S >= n:
            continue
        ok = 0
        tol = min(5.0, 360.0 / S / 3)
        for _ in range(trials):
            th = np.sort((2 * np.pi * (np.arange(S) + rng.uniform(-0.2, 0.2, S)) / S + rng.uniform(0, 2 * np.pi)) % (2 * np.pi))
            a = np.exp(1j * k * (np.cos(th)[:, None] * p[:, 0] + np.sin(th)[:, None] * p[:, 1])).T   # (N, S)
            g = 10 ** (rng.normal(0, 0.5, n) / 20) * np.exp(1j * np.radians(rng.normal(0, sigma_cal, n)))
            s = (rng.normal(size=(S, k_samples)) + 1j * rng.normal(size=(S, k_samples))) / np.sqrt(2)
            w = (rng.normal(size=(n, k_samples)) + 1j * rng.normal(size=(n, k_samples))) / np.sqrt(2) * 10 ** (-snr_db / 20)
            x = (a * g[:, None]) @ s + w
            R = x @ x.conj().T / k_samples
            _, V = np.linalg.eigh(R)
            En = V[:, : n - S]
            ps = 1.0 / np.maximum(np.sum(np.abs(A.conj() @ En) ** 2, axis=1), 1e-12)
            pk = np.nonzero((ps > np.roll(ps, 1)) & (ps >= np.roll(ps, -1)))[0]
            if len(pk) < S:
                continue
            est = np.sort(grid[pk[np.argsort(ps[pk])[::-1][:S]]])
            best = min(np.max(np.abs(np.degrees(np.angle(np.exp(1j * (np.roll(est, r) - th)))))) for r in range(S))
            ok += best < tol
        out[int(S)] = ok / trials
        if out[int(S)] < 0.5:
            break
    return out


def slot_counts(n: int, m: int = 4):
    """Slots needed when m coherent channels are shared by n elements."""
    if n <= m:
        return dict(groups=1, ref=1, pairs_bound=1, pairs_revb=1)
    return dict(groups=math.ceil(n / m), ref=math.ceil((n - 1) / (m - 1)),
                pairs_bound=math.ceil(n / m * math.ceil((n - 1) / (m - 1))),   # Schoenheim bound on C(n, 4, 2)
                pairs_revb=math.ceil(n / m) ** 2)                               # each channel owns n/m ports


def dof(f: float, r: float = R_OUT, h: float = Z_HALF):
    """Spatial degrees of freedom of the aperture: azimuth modes of a disc 2 k r + 1; modes of the sphere that
    encloses the volume (k a + 1)^2."""
    k = 2 * np.pi * f / C0
    return 2 * k * r + 1, (k * math.hypot(r, h) + 1) ** 2, (k * r + 1) ** 2


# ---- the study --------------------------------------------------------------------------------------------
def search(n: int, dim: int, canon: dict, quick: bool, log):
    smin = s_min(n, dim)
    mode = "az" if dim == 2 else "azel"
    score = make_score3(n, dim, smin, mode, n_freq=8 if quick else 12, trials=24 if quick else 64)
    rep = lambda x: repair3(x, dim, smin)
    lo = -R_OUT * np.ones(dim * n)
    hi = R_OUT * np.ones(dim * n)
    if dim == 3:
        lo[2::3], hi[2::3] = -Z_HALF, Z_HALF
    t0 = time.time()
    rs = random_search(score, lo, hi, batches=3 if quick else 12, batch=16 if quick else 40, seed=n + 10 * dim, repair=rep)
    cands = [rep(c[:, :dim].reshape(-1)) for c in canon.values()] + list(rs.finalists[:5])
    X, F = reselect(score, cands, draws=2 if quick else 6)
    hc = climb(score, X[0], lo, hi, steps=4 if quick else 50, lam=8, sigma=0.03, seed=n + 10 * dim, repair=rep)
    X2, F2 = reselect(score, np.vstack([X[:3], hc.x[None]]), draws=2 if quick else 10, first_gen=2_000_000)
    p = X2[0].reshape(n, dim)
    assert violation3(p, smin) < 1e-5, violation3(p, smin)
    log(f"  N={n} {dim}-D search: {rs.evaluations + hc.evaluations} evaluations, {time.time() - t0:.0f} s; "
        f"start {F[0]:.2f} -> {F2[0]:.2f} (climbed kept: {bool(np.allclose(X2[0], hc.x))})")
    return p, float(F2[0])


def selftest(log):
    """The generalised search model must reproduce evolve_layout's exactly for planar designs."""
    n = 8
    X = np.array([EL.repair(EL.ring(n, 0.112).reshape(-1)), EL.repair(EL.ring(n, 0.2).reshape(-1))])
    a = EL.make_score(n, n_freq=6, trials=24)(X, 3)
    b = make_score3(n, 2, EL.S_MIN, "az", n_freq=6, trials=24)(X, 3)
    assert np.allclose(a, b), (a, b)
    log(f"self-test: make_score3('az') == evolve_layout.make_score: {a.round(4)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--ns", type=str, default="")
    ap.add_argument("--extras", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    if args.extras:
        extras()
        return
    if args.report:
        report(json.loads((OUT / "results.json").read_text()))
        return
    quick = args.quick
    ns = [int(v) for v in args.ns.split(",")] if args.ns else ([4, 8, 16] if quick else NS)
    freqs = FREQS[::6] if quick else FREQS
    logf = open(OUT / ("run_quick.log" if quick else "run.log"), "a")

    def log(*a):
        print(*a, flush=True)
        print(*a, file=logf, flush=True)

    selftest(log)
    name = OUT / ("results_quick.json" if quick else "results.json")
    res = json.loads(name.read_text()) if name.exists() and args.ns else dict(designs={})
    res.update(freqs_hz=freqs.tolist(), conditions=CONDS, astro_hz=ASTRO,
               volume=dict(r_out_m=R_OUT, r_in_m=R_IN, z_half_m=Z_HALF))
    for n in ns:
        t0 = time.time()
        trials = 100 if quick else (600 if n <= 24 else 400)
        sm2 = s_min(n, 2)
        planar = {"ring": EL.ring(n, 0.225)}
        if n >= 8:
            planar["rings"] = rings(n, sm2)
        planar = {k: repair3(v.reshape(-1), 2, sm2).reshape(n, 2) for k, v in planar.items()}
        vol = {"stack": stack(n), "cone": cone(n, 0.075), "helix": helix(n)}
        vol = {k: repair3(v.reshape(-1), 3, 0.075).reshape(n, 3) for k, v in vol.items()}
        planar["opt"], sc2 = search(n, 2, planar, quick, log)
        vol["opt"], sc3 = search(n, 3, vol, quick, log)
        for dim, fam in ((2, planar), (3, vol)):
            for kname, p in fam.items():
                key = f"{'planar' if dim == 2 else 'vol'}_{kname}_{n}"
                d = dict(n=n, dim=dim, family=kname, positions_m=np.round(to3(p), 4).tolist(),
                         min_spacing_mm=round(min_spacing(p) * 1000, 1))
                line = f"  {key:18s}"
                for cname, c in CONDS.items():
                    if dim == 2:
                        r = EL.check(p, freqs, trials=trials, k_samples=64, **c)
                        d[f"az_{cname}"] = dict(median=r[:, 0].tolist(), p95=r[:, 1].tolist(), gross_pct=r[:, 2].tolist())
                        line += f" | az-only {cname}: p95 {r[:, 1].max():5.2f} gross {r[:, 2].max():4.1f}"
                    r = check3(p, freqs, trials=trials, **c)
                    d[f"azel_{cname}"] = dict(median=r[:, 0].tolist(), p95=r[:, 1].tolist(), gross_pct=r[:, 2].tolist(),
                                              el_median=r[:, 3].tolist(), el_p95=r[:, 4].tolist(), el_gross_pct=r[:, 5].tolist())
                    line += f" | az-el {cname}: p95 {r[:, 1].max():5.2f} gross {r[:, 2].max():4.1f} el p95 {r[:, 4].max():5.1f}"
                d["crb_az_el_deg"] = [crb(p, f) for f in freqs]
                d["directivity_dbi"] = [directivity(p, f) for f in freqs]
                d["sidelobes"] = [sidelobes(p, f) for f in freqs]
                res["designs"][key] = d
                log(line)
        # channel sharing and source counts, on the layouts that did best with the 2-D estimator
        for dim, fam in ((2, planar), (3, vol)):
            pre = "planar" if dim == 2 else "vol"
            best = min(fam, key=lambda q: (max(res["designs"][f"{pre}_{q}_{n}"]["azel_elevated"]["gross_pct"]),
                                           max(res["designs"][f"{pre}_{q}_{n}"]["azel_elevated"]["p95"])))
            d = res["designs"][f"{pre}_{best}_{n}"]
            res.setdefault("best", {})[f"{pre}_{n}"] = best
            p = np.array(d["positions_m"])
            for sharing in ("groups", "ref"):
                if n <= 4 or (sharing == "groups" and n % 4):
                    continue
                for cname in ("elevated", "stress"):
                    r = check3(p, freqs, sharing=sharing, trials=trials, **CONDS[cname])
                    d[f"{sharing}_{cname}"] = dict(median=r[:, 0].tolist(), p95=r[:, 1].tolist(), gross_pct=r[:, 2].tolist(),
                                                   el_median=r[:, 3].tolist(), el_p95=r[:, 4].tolist(), el_gross_pct=r[:, 5].tolist())
                    log(f"  {pre}_{best}_{n} shared 4 channels, {sharing:6s} {cname}: p95 {r[:, 1].max():5.2f} "
                        f"gross {r[:, 2].max():4.1f} median {np.median(r[:, 0]):.2f}")
            if dim == 2:
                counts = [1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 15, 19, 23, 27, 31, 39, 47, 55, 63]
                d["music"] = {f"{f / 1e9:.2f}": music_capacity(p, f, counts, trials=20 if quick else 60)
                              for f in (0.4e9, 1.0e9, 1.4204e9, 2.3e9)}
                log(f"  {pre}_{best}_{n} MUSIC success by source count: " + "; ".join(
                    f"{f} GHz {max([s for s, v in m.items() if v >= 0.9], default=0)}" for f, m in d["music"].items()))
        res.setdefault("slots", {})[str(n)] = slot_counts(n)
        name.write_text(json.dumps(res, indent=1))
        log(f"N={n} done in {time.time() - t0:.0f} s")
    report(res, quick)


# ---- report -----------------------------------------------------------------------------------------------
def report(res, quick=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    D = res["designs"]
    f = np.array(res["freqs_hz"])
    ns = sorted({d["n"] for d in D.values()})
    lines = []
    w = lines.append
    fi = lambda x: int(np.argmin(np.abs(f - x)))
    lowhalf = f <= 1.0e9

    def band(v, sel=None):
        v = np.array(v)
        return v if sel is None else v[sel]

    w("TABLE 1. Direction finding against N (N coherent channels; 3 deg calibration residual, 20 dB, K = 64, sources 0-30 deg up).")
    w("p95 = 95th percentile azimuth error: worst over 0.4-1 GHz / worst over 1-2.3 GHz; gross = worst over the band, % of fixes off by > 10 deg (stress: 5 deg, 0 dB).")
    w(f"{'layout':20s} {'min mm':>6s} | {'az-only p95':>13s} {'gross':>5s} {'stress':>6s} | {'az-el p95':>13s} {'gross':>5s} {'stress':>6s} | {'el p95':>11s} {'el gross':>8s} | bound az p95 (0.4 / 2.3 GHz)")
    for n in ns:
        for key, d in D.items():
            if d["n"] != n:
                continue
            a = d.get("az_elevated")
            az = (f"{band(a['p95'], lowhalf).max():5.2f} / {band(a['p95'], ~lowhalf).max():5.2f} {max(a['gross_pct']):5.1f} "
                  f"{max(d['az_stress']['gross_pct']):6.1f}") if a else f"{'':13s} {'':5s} {'':6s}"
            e = d["azel_elevated"]
            cr = np.array(d["crb_az_el_deg"])
            w(f"{key:20s} {d['min_spacing_mm']:6.0f} | {az} | {band(e['p95'], lowhalf).max():5.2f} / {band(e['p95'], ~lowhalf).max():5.2f} "
              f"{max(e['gross_pct']):5.1f} {max(d['azel_stress']['gross_pct']):6.1f} | {band(e['el_p95'], lowhalf).max():4.1f} / {band(e['el_p95'], ~lowhalf).max():4.1f} "
              f"{max(e['el_gross_pct']):8.1f} | {1.96 * cr[0, 0]:.2f} / {1.96 * cr[-1, 0]:.2f}")
        w("")
    w("TABLE 2. Beams against N (uniform weights, isotropic elements, steered to the horizon). Directivity dBi and peak sidelobe dB (whole sphere, mean over steering azimuth) at 0.408 / 1.0 / 1.42 / 2.3 GHz; 10 log N for reference.")
    for n in ns:
        for key, d in D.items():
            if d["n"] != n:
                continue
            ix = [fi(0.408e9), fi(1.0e9), fi(1.4204e9), fi(2.3e9)]
            w(f"{key:20s} 10logN {10 * np.log10(n):4.1f} | D " + " / ".join(f"{d['directivity_dbi'][i]:4.1f}" for i in ix)
              + " | PSL " + " / ".join(f"{d['sidelobes'][i]['sphere_mean']:5.1f}" for i in ix)
              + " | PSL horizon " + " / ".join(f"{d['sidelobes'][i]['horizon_mean']:5.1f}" for i in ix)
              + " | worst PSL in band " + f"{max(s['sphere_mean'] for s in d['sidelobes']):5.1f}"
              + " | beamwidth " + " / ".join(f"{d['sidelobes'][i]['beamwidth']:3.0f}" for i in ix))
        w("")
    w("TABLE 3. Four shared channels against N coherent channels (best layout by the az-el estimator). Azimuth median (band median) / p95 (band worst) / gross % (band worst), elevated; then gross % under stress.")
    for n in ns:
        for pre in ("planar", "vol"):
            b = res.get("best", {}).get(f"{pre}_{n}")
            if not b:
                continue
            d = D[f"{pre}_{b}_{n}"]
            row = f"{pre}_{b}_{n}".ljust(20)
            for mode, lab in (("azel", "N channels"), ("ref", "reference + 3 roving"), ("groups", "disjoint groups")):
                e = d.get(f"{mode}_elevated")
                if e:
                    row += f" | {lab}: {np.median(e['median']):.2f} / {max(e['p95']):5.2f} / {max(e['gross_pct']):4.1f}; stress {max(d[f'{mode}_stress']['gross_pct']):4.1f}"
            w(row)
    w("")
    w("TABLE 4. Slots per fix with 4 shared channels, sources MUSIC separates with N channels (success >= 90 %) and aperture degrees of freedom.")
    w(f"{'N':>3s} | slots: {'groups':>6s} {'ref':>4s} {'all pairs >=':>12s} {'cross-ch pairs':>15s} | MUSIC sources at 0.4 / 1.0 / 1.42 / 2.3 GHz (limit N-1) | port cost $ | 2-ch radios for N channels")
    for n in ns:
        s = res["slots"][str(n)]
        b = res.get("best", {}).get(f"planar_{n}")
        mu = D[f"planar_{b}_{n}"].get("music", {}) if b else {}
        cap = []
        for fk in ("0.40", "1.00", "1.42", "2.30"):
            m = {int(k2): v for k2, v in mu.get(fk, {}).items()}
            good = 0
            for sct in sorted(m):
                if m[sct] >= 0.9:
                    good = sct
                else:
                    break
            cap.append(good)
        w(f"{n:3d} | {s['groups']:13d} {s['ref']:4d} {s['pairs_bound']:12d} {s["pairs_revb"]:14d} | {' / '.join(f'{c:2d}' for c in cap)} (limit {n - 1}) "
          f"| {n * COST_PORT:6.0f} | {math.ceil(n / 2)} (${math.ceil(n / 2) * 540})")
    w("")
    w("Aperture degrees of freedom (azimuth modes of the 240 mm disc 2kr+1; sphere modes of the disc (kr+1)^2; of the 500 x 300 mm volume (ka+1)^2):")
    for fq in (0.408e9, 0.611e9, 1.0e9, 1.4204e9, 1.6654e9, 2.3e9):
        a, v, pl = dof(fq)
        w(f"  {fq / 1e9:.3f} GHz: azimuth {a:5.1f}; planar sphere modes {pl:6.1f}; volume {v:6.1f}")
    txt = "\n".join(lines)
    (OUT / ("tables_quick.txt" if quick else "tables.txt")).write_text(txt + "\n")
    print(txt)

    # figures: the best layout of each kind against N
    def best_series(pre, getter):
        xs, ys = [], []
        for n in ns:
            b = res.get("best", {}).get(f"{pre}_{n}")
            if b:
                xs.append(n); ys.append(getter(D[f"{pre}_{b}_{n}"]))
        return xs, ys

    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    for pre, col in (("planar", "C0"), ("vol", "C3")):
        x, y = best_series(pre, lambda d: np.max(np.array(d["azel_elevated"]["p95"])[lowhalf]))
        ax[0, 0].loglog(x, y, "o-", color=col, label=f"{pre}, 0.4-1 GHz")
        x, y = best_series(pre, lambda d: np.max(np.array(d["azel_elevated"]["p95"])[~lowhalf]))
        ax[0, 0].loglog(x, y, "s--", color=col, label=f"{pre}, 1-2.3 GHz")
        x, y = best_series(pre, lambda d: max(d["azel_stress"]["gross_pct"]))
        ax[0, 1].semilogx(x, y, "o-", color=col, label=f"{pre}: N channels, stress")
        x, y = best_series(pre, lambda d: max(d["groups_stress"]["gross_pct"]) if "groups_stress" in d else np.nan)
        ax[0, 1].semilogx(x, y, "^:", color=col, label=f"{pre}: 4 channels, disjoint groups, stress")
        x, y = best_series(pre, lambda d: max(s["sphere_mean"] for s in d["sidelobes"]))
        ax[1, 0].semilogx(x, y, "o-", color=col, label=f"{pre}: worst in band")
        x, y = best_series(pre, lambda d: d["directivity_dbi"][fi(0.408e9)])
        ax[1, 1].semilogx(x, y, "o-", color=col, label=f"{pre}, 408 MHz")
        x, y = best_series(pre, lambda d: d["directivity_dbi"][fi(2.3e9)])
        ax[1, 1].semilogx(x, y, "s--", color=col, label=f"{pre}, 2.3 GHz")
    xx = np.array(ns, float)
    ax[0, 0].loglog(xx, ax[0, 0].lines[0].get_ydata()[0] * np.sqrt(ns[0] / xx), "k:", label="1 / sqrt(N)")
    ax[1, 0].semilogx(xx, -10 * np.log10(xx), "k:", label="-10 log N (mean sidelobe of a random array)")
    ax[1, 1].semilogx(xx, 10 * np.log10(xx), "k:", label="10 log N")
    ax[0, 0].set_title("95th percentile azimuth error, az-el estimator, sources 0-30 deg up (deg)", fontsize=9)
    ax[0, 1].set_title("Gross errors (> 10 deg), worst frequency, 5 deg and 0 dB (%)", fontsize=9)
    ax[1, 0].set_title("Peak sidelobe over the sphere, beam on the horizon (dB)", fontsize=9)
    ax[1, 1].set_title("Directivity of the beam (dBi)", fontsize=9)
    for a in ax.ravel():
        a.set_xlabel("elements N"); a.grid(alpha=0.3, which="both"); a.legend(fontsize=7)
        a.set_xticks(ns); a.set_xticklabels([str(n) for n in ns]); a.minorticks_off()
    fig.tight_layout()
    fig.savefig(OUT / ("count_vs_n_quick.png" if quick else "count_vs_n.png"), dpi=130)
    plt.close(fig)


# ---- extras: theory against the stored results, and the fixed-time experiment -----------------------------
def extras():
    """python experiments/count_study.py --extras   (about 5 minutes; writes count/extras.txt and extras.json)

    E1  equal observing time: 4 shared channels and N elements, the samples per slot cut in proportion to the
        number of slots (K = 4096 / N), against N coherent channels with all the samples. Noise-limited
        (no phase residual, -10 dB) and calibration-limited (3 deg, 20 dB).
    E2  closed-form bounds for rings and stacked rings against the stored Monte Carlo.
    E3  peak sidelobe of the searched (irregular) layouts against ln(M) / N; rings against J0.
    E4  gross errors against the highest sidelobe.
    E5  time per fix; E6 radio-astronomy sensitivity; E7 cost."""
    res = json.loads((OUT / "results.json").read_text())
    D = res["designs"]
    f = np.array(res["freqs_hz"])
    fi = lambda x: int(np.argmin(np.abs(f - x)))
    L = []
    w = L.append
    out = {}

    # E1 ----------------------------------------------------------------------------------------------------
    w("E1. Equal observing time (4096 channel-samples per fix). 95th percentile azimuth error, deg, sources 0-30 deg up,")
    w("    at 0.408 / 1.0 / 1.42 / 2.3 GHz. 'shared' = 4 channels, disjoint groups, K = 4096/N per slot; 'full' = N channels, K = 1024.")
    fq = [0.408e9, 1.0e9, 1.4204e9, 2.3e9]
    e1 = {}
    for cname, c in (("noise-limited (0 deg residual, -10 dB)", dict(sigma_cal=0.0, snr_db=-10.0)),
                     ("calibration-limited (3 deg, 20 dB)", dict(sigma_cal=3.0, snr_db=20.0))):
        w(f"  {cname}")
        for n in (4, 8, 16, 32, 64):
            p = np.array(D[f"planar_{'opt' if n == 4 else 'ring'}_{n}"]["positions_m"])
            row = f"    N={n:2d} "
            if n > 4:
                r = check3(p, fq, sharing="groups", trials=300, k_samples=4096 // n, el_max=30.0, seed=5, **c)
                row += "shared p95 " + " / ".join(f"{v:5.2f}" for v in r[:, 1]) + f" gross {r[:, 2].max():4.1f} | "
                e1[f"{cname}|shared|{n}"] = r.tolist()
            else:
                row += " " * 56
            r = check3(p, fq, sharing="full", trials=300, k_samples=1024, el_max=30.0, seed=5, **c)
            row += "full p95 " + " / ".join(f"{v:5.2f}" for v in r[:, 1]) + f" gross {r[:, 2].max():4.1f}"
            e1[f"{cname}|full|{n}"] = r.tolist()
            w(row)
            print(row, flush=True)
    out["equal_time"] = e1

    # E2 ----------------------------------------------------------------------------------------------------
    w("")
    w("E2. Closed form: sigma_az = s / (k r cos(el) sqrt(N/2)) for a ring; sigma_el = s / (k sqrt(N) sqrt(r^2 sin^2(el)/2 + h^2 cos^2(el)))")
    w("    for rings at heights +/-h (h = 0 planar). s^2 = sigma_cal^2 + 1/(2 K rho). 95th percentile = 1.96 sigma, RMS over el 2-30 deg.")
    w("    N | ring az p95 closed form / Monte Carlo at 0.408 GHz | planar el p95 closed / MC | stack el p95 closed / MC at 0.408 GHz")
    s = math.sqrt(np.radians(3.0) ** 2 + 1 / (2 * 64 * 100))
    els = np.radians(np.linspace(2, 30, 8))
    k = 2 * np.pi * 0.408e9 / C0
    i408 = fi(0.408e9)
    e2 = {}
    for n in NS:
        if n < 8:
            continue
        r_, h_ = 0.225, Z_HALF
        az = np.degrees(1.96 * math.sqrt(np.mean((s / (k * r_ * np.cos(els) * math.sqrt(n / 2))) ** 2)))
        el_pl = np.degrees(1.96 * math.sqrt(np.mean((s / (k * math.sqrt(n) * r_ * np.sin(els) / math.sqrt(2))) ** 2)))
        el_st = np.degrees(1.96 * math.sqrt(np.mean(s ** 2 / (k ** 2 * n * (r_ ** 2 * np.sin(els) ** 2 / 2 + h_ ** 2 * np.cos(els) ** 2)))))
        dr, ds = D[f"planar_ring_{n}"]["azel_elevated"], D[f"vol_stack_{n}"]["azel_elevated"]
        w(f"   {n:2d} | {az:5.2f} / {dr['p95'][i408]:5.2f} | {el_pl:5.1f} / {dr['el_p95'][i408]:5.1f} | {el_st:5.2f} / {ds['el_p95'][i408]:5.2f}")
        e2[n] = dict(az_closed=az, az_mc=dr['p95'][i408], el_planar_closed=el_pl, el_planar_mc=dr['el_p95'][i408],
                     el_stack_closed=el_st, el_stack_mc=ds['el_p95'][i408])
    out["closed_form"] = e2
    # exponent of the 95th percentile against N (best layouts, N >= 8)
    for pre in ("planar", "vol"):
        xs = [n for n in NS if n >= 8]
        for lab, sel in (("0.4-1 GHz", f <= 1e9), ("1-2.3 GHz", f > 1e9)):
            ys = [np.max(np.array(D[f"{pre}_{res['best'][f'{pre}_{n}']}_{n}"]["azel_elevated"]["p95"])[sel]) for n in xs]
            b = np.polyfit(np.log(xs), np.log(ys), 1)[0]
            w(f"    fitted exponent of the worst 95th percentile against N, {pre}, {lab}: {b:+.2f} (theory -0.50)")
            out[f"exponent_{pre}_{lab}"] = float(b)

    # E3 ----------------------------------------------------------------------------------------------------
    w("")
    w("E3. Peak sidelobe (dB, whole sphere, beam on the horizon). Irregular layouts: theory 10 log(ln(M) / N), M = (k a + 1)^2 beam positions;")
    w("    rings: J0's first sidelobe, -7.9 dB, once N exceeds the ring's 2 k r + 1 azimuth modes; below that, grating lobes near 0 dB.")
    w("    N | planar searched, 1.0 GHz: theory / measured | 2.3 GHz: theory / measured | volume searched 2.3 GHz: theory / measured | ring 2.3 GHz measured (modes 22.7)")
    e3 = {}
    for n in NS:
        row = f"   {n:2d} |"
        for pre, fqs in (("planar", (1.0e9, 2.3e9)), ("vol", (2.3e9,))):
            for fx in fqs:
                az_m, vol_m, pl_m = dof(fx)
                M = pl_m if pre == "planar" else vol_m
                th = 10 * math.log10(min(1.0, math.log(M) / n))
                me = D[f"{pre}_opt_{n}"]["sidelobes"][fi(fx)]["sphere_mean"]
                row += f" {th:5.1f} / {me:5.1f} |"
                e3[f"{pre}_{n}_{fx / 1e9:.1f}"] = (th, me)
        row += f" {D[f'planar_ring_{n}']['sidelobes'][fi(2.3e9)]['sphere_mean']:5.1f}"
        w(row)
    out["sidelobes"] = e3

    # E4 ----------------------------------------------------------------------------------------------------
    w("")
    w("E4. Gross errors (az-el estimator, stress: 5 deg, 0 dB, K = 64) against the highest sidelobe in the band (worst steering azimuth).")
    w("    Theory: a sidelobe of power b wins with probability about Q((1 - sqrt(b)) sqrt(2 N) / (s sqrt(1 - b))), s = 0.124 rad here;")
    w("    below about -1.5 dB that is under 1e-3 for N >= 6. Planar and volume layouts, N <= 16:")
    e4 = []
    qf = lambda x: 0.5 * math.erfc(x / math.sqrt(2))
    for key, d in D.items():
        if d["n"] > 16:
            continue
        psl = max(sl["sphere_worst"] for sl in d["sidelobes"])
        b = 10 ** (psl / 10)
        s_ = math.sqrt(np.radians(5.0) ** 2 + 1 / (2 * 64 * 1.0))
        q = qf((1 - math.sqrt(b)) * math.sqrt(2 * d["n"]) / (s_ * math.sqrt(max(1 - b, 1e-9)))) if b < 1 else 0.5
        g = max(d["azel_stress"]["gross_pct"])
        e4.append((key, psl, q * 100, g))
        w(f"    {key:18s} worst sidelobe {psl:5.1f} dB  theory per-lobe {q * 100:6.2f} %  simulated gross {g:5.1f} %")
    out["gross_vs_psl"] = e4

    # E5 ----------------------------------------------------------------------------------------------------
    w("")
    w("E5. Time per fix, ms: slots x (15 us switch slot + K / B), K = 64 independent samples, B = occupied bandwidth (<= 20 MHz sample rate).")
    w("    N | N channels | 4 channels: disjoint groups | reference + 3 roving | every pair (full covariance, bound) || at B = 12.5 kHz / 200 kHz / 2 MHz / 20 MHz")
    e5 = {}
    for n in NS:
        sc = slot_counts(n)
        row = f"   {n:2d} |"
        for lab in ("one", "groups", "ref", "pairs_bound"):
            sl = 1 if lab == "one" else sc[lab]
            t = [sl * (15e-6 + 64 / B) * 1e3 for B in (12.5e3, 200e3, 2e6, 20e6)]
            e5[f"{n}_{lab}"] = t
            row += f" {sl:3d} slots: " + " / ".join(f"{v:7.2f}" for v in t) + " |"
        w(row)
    out["fix_time_ms"] = e5

    # E6 ----------------------------------------------------------------------------------------------------
    w("")
    w("E6. Radio astronomy (estimates). Element: whip over the plate, gain 2 dBi towards the source; Tsys = receiver (NF 1.2 dB -> 92 K) + sky and ground")
    w("    (408 MHz: 60 K; others 30 K). SEFD = 2 k Tsys / (m A_e), A_e = G lambda^2 / 4 pi, m = elements combined coherently at once.")
    w("    dS = SEFD / sqrt(B t), t = 600 s. Independent sky pixels = (k r + 1)^2 (planar plate); baselines N (N - 1) / 2.")
    kB = 1.380649e-23
    e6 = {}
    for name, fx, B, tsky in (("408 MHz", 0.408e9, 3.9e6, 60.0), ("611 MHz", 0.611e9, 6e6, 30.0), ("H I 1420 MHz", 1.4204e9, 27e6, 30.0),
                              ("OH 1665 MHz", 1.6654e9, 10e6, 30.0)):
        lam = C0 / fx
        ae = 10 ** 0.2 * lam ** 2 / (4 * np.pi)
        tsys = 92.0 + tsky
        row = f"    {name:13s} lambda {lam * 100:5.1f} cm, A_e {ae * 1e4:6.1f} cm^2, Tsys {tsys:.0f} K, B {B / 1e6:.1f} MHz, beam of the 480 mm plate {np.degrees(lam / 0.48):3.0f} deg, sky pixels {dof(fx)[2]:5.1f} | SEFD (MJy) and 1-sigma (Jy) for m ="
        for m in (1, 4, 16, 64):
            sefd = 2 * kB * tsys / (m * ae) / 1e-26
            row += f" {m}: {sefd / 1e6:6.1f}, {sefd / math.sqrt(B * 600):6.0f};"
            e6[f"{name}_{m}"] = (sefd, sefd / math.sqrt(B * 600))
        w(row)
    w("    baselines N(N-1)/2: " + ", ".join(f"N={n}: {n * (n - 1) // 2}" for n in NS))
    out["astronomy"] = e6

    # E7 ----------------------------------------------------------------------------------------------------
    w("")
    w("E7. Cost against N (external ports; rev B section 10: port front end $11.25, panel cell $14.56, channel chain $42.5, bladeRF $270 per channel;")
    w("    whip and cable $8 per element, estimate from the $125 ring kit). 'switched' = 4 channels kept; 'coherent' = N channels.")
    w("    N | switched: added to the 4-channel core $ | $ per element | coherent: channels + radios $ | 95th pct az 0.4-1 GHz (planar / volume) | elevation 95th (volume) | deg gained per $100 from the previous N")
    prev = None
    e7 = {}
    for n in NS:
        sw = n * (COST_PORT + 8.0)
        co = n * (COST_PORT + 8.0) + max(n - 4, 0) * (42.5 + 270.0)
        pp = np.max(np.array(D[f"planar_{res['best'][f'planar_{n}']}_{n}"]["azel_elevated"]["p95"])[f <= 1e9])
        pv = np.max(np.array(D[f"vol_{res['best'][f'vol_{n}']}_{n}"]["azel_elevated"]["p95"])[f <= 1e9])
        ev = np.max(np.array(D[f"vol_{res['best'][f'vol_{n}']}_{n}"]["azel_elevated"]["el_p95"]))
        gain = "" if prev is None else f"{(prev[1] - pp) / ((sw - prev[0]) / 100):6.2f}"
        w(f"   {n:2d} | {sw:6.0f} | {COST_PORT + 8.0:5.2f} | {co:7.0f} | {pp:5.2f} / {pv:5.2f} | {ev:5.2f} | {gain}")
        e7[n] = dict(switched=sw, coherent=co, p95_planar=pp, p95_vol=pv, el_p95_vol=ev)
        prev = (sw, pp)
    out["cost"] = e7
    (OUT / "extras.txt").write_text("\n".join(L) + "\n")
    (OUT / "extras.json").write_text(json.dumps(out, indent=1))
    print("\n".join(L))


if __name__ == "__main__":
    main()
