"""Band plan against geometry: which bands does each array serve, and what does choosing the bands first buy?

The owner's question (2026-11-15): "we can start restricting the frequency range if such a wide range seems
futile, we can be more selective, even cherry pick bands of interest, optimal frequency range sets for a given
geometry". Three parts, all by simulation, all with the models of experiments/evolve_layout.py:

1. Geometry -> bands. For rev B's geometries (the 25 mm 4x4 panel; the fixed whip layouts rings16_225_112mm
   and random8 of experiments/evolve/results.json; a half-wavelength 4x4 grid of whips) figures of merit over
   0.3-8 GHz: direction-finding error and gross-error rate (evolve_layout.check, sources 0-15 and 0-30 deg up),
   peak sidelobe (ambiguity, grating lobes), beamwidth and directivity. Windows where each geometry is good.
2. Bands -> geometry. For a few band sets the whip layout is optimised for those bands only (array16.evolve:
   random search, genetic algorithm, hill climb from the full-band layout; finalists re-scored with reselect;
   same plate and spacing rules) and compared with the full-band layouts at the same bands.
3. The panel's pitch if its top band were 2.5 GHz or 5.9 GHz only, against rev B's 25 mm.

Models. Whips: evolve_layout.check unchanged (azimuth-only Bartlett beamscan on a horizontal plate, isotropic
whips, 3 deg / 0.5 dB calibration residual, 20 dB per whip, 64 samples; gross = error above 10 deg). Panel:
the same Monte Carlo with a two-dimensional scan (check_panel below): the panel stands upright, sources are
within +/-60 deg of its normal in azimuth and 0-15 or 0-30 deg up, the estimator scans the whole front
hemisphere in direction cosines, and the error is the angle between the true and the estimated direction.
Elements are isotropic in the Monte Carlo (so are the whips); directivity uses a cos(theta) element for the
panel and isotropic elements for the whips (the whip's own 2-5 dBi is not included).

    python experiments/band_plan_study.py            # about 30 minutes with OMP_NUM_THREADS=4 on a busy machine;
                                                     # an interrupted run resumes part 2 from its finished searches
    python experiments/band_plan_study.py --quick    # a smoke run (about 1 minute), writes *_quick files
    python experiments/band_plan_study.py --parts 1,3
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import evolve_layout as el  # noqa: E402  (read-only: check, repair, violation, ring, _refine, GRID, plate rules)
from array16.evolve import climb, evolve, random_search, reselect  # noqa: E402

try:  # read-only helpers of the wideband study (grating-lobe scan limit, closest-pair ambiguity)
    import wideband_study as wb  # noqa: E402
except Exception:  # pragma: no cover - the study does not depend on it
    wb = None

C0 = el.C0
OUT = HERE / "band_plan"
PITCH_B = 0.025  # rev B panel pitch

CONDITIONS = dict(nominal=dict(sigma_cal=3.0, snr_db=20.0, el_max=15.0, seed=12345),
                  elevated=dict(sigma_cal=3.0, snr_db=20.0, el_max=30.0, seed=4242),
                  stress=dict(sigma_cal=5.0, snr_db=0.0, el_max=30.0, seed=777))

# Bands of interest, Hz (sources: see band_plan/README in the report; allocations are cited there).
BANDS = {
    "433": (433.05e6, 434.79e6),      # ISM (ITU Region 1) / EU short-range devices
    "868": (863e6, 870e6),            # EU short-range devices
    "915": (902e6, 928e6),            # ISM (ITU Region 2), FCC 15.247 / 15.249
    "1090": (1087e6, 1093e6),         # ADS-B / Mode S replies
    "1296": (1240e6, 1300e6),         # 23 cm amateur
    "1420": (1420.0e6, 1420.8e6),     # hydrogen line
    "1575": (1565.42e6, 1585.42e6),   # GPS L1 / Galileo E1
    "1621": (1616e6, 1626.5e6),       # Iridium
    "2440": (2400e6, 2483.5e6),       # ISM
    "3600": (3550e6, 3700e6),         # CBRS (US)
    "5500": (5150e6, 5850e6),         # 5 GHz Wi-Fi (U-NII-1 to U-NII-3)
    "5800": (5725e6, 5875e6),         # ISM
}
BAND_SETS = {
    "A_433_915_2440": ["433", "915", "2440"],
    "B_915_1575_2440": ["915", "1575", "2440"],
    "C_2440_5800": ["2440", "5800"],
}
SINGLE_BANDS = ["433", "2440", "5800"]


def centre(b: str) -> float:
    lo, hi = BANDS[b]
    return 0.5 * (lo + hi)


# ---- geometries -------------------------------------------------------------------------------------------
def grid(n_side: int, pitch: float) -> np.ndarray:
    c = (np.arange(n_side) - (n_side - 1) / 2) * pitch
    x, y = np.meshgrid(c, c)
    return np.column_stack([x.ravel(), y.ravel()])


def load_layouts() -> dict:
    d = json.loads((HERE / "evolve" / "results.json").read_text())["designs"]
    return {k: np.array(v["positions_m"]) for k, v in d.items()}


# ---- deterministic figures of merit (horizon manifold of a whip layout) ------------------------------------
def horizon_pattern_metrics(p: np.ndarray, f: float, el_test: float = 30.0) -> dict:
    """Beam patterns on the horizon for 72 steering azimuths: peak sidelobe (worst steer; 0 dB = a second
    direction that looks the same, i.e. a grating lobe / ambiguity), 3 dB beamwidth (mean), directivity of the
    array factor (isotropic elements, closed form, mean and worst steer), and what a source `el_test` deg up
    does to an estimator that assumes the horizon: the loss of the best match and the margin between the best
    match within 10 deg of the true azimuth and the best match elsewhere (negative = a wrong azimuth wins)."""
    k = 2 * np.pi * f / C0
    n = len(p)
    scan = np.radians(np.arange(0, 360, 0.25))
    steer_i = np.arange(0, len(scan), 20)                     # every 5 deg
    A = np.exp(1j * k * (np.cos(scan)[:, None] * p[:, 0] + np.sin(scan)[:, None] * p[:, 1]))   # (G, N)
    Bm = np.abs(A[steer_i].conj() @ A.T) ** 2 / n ** 2                                         # (S, G)
    G = len(scan)
    roll = (steer_i[:, None] + np.arange(G)[None, :]) % G
    R = np.take_along_axis(Bm, roll, 1)                       # peak at column 0
    psl, bw = [], []
    for r in R:
        up = np.nonzero((r[1:-1] <= r[:-2]) & (r[1:-1] <= r[2:]))[0] + 1    # local minima
        if len(up) < 2:
            psl.append(-99.0)
        else:
            rest = r[up[0]:up[-1] + 1]
            psl.append(10 * np.log10(max(rest.max(), 1e-10)))
        below = np.nonzero(r < 0.5)[0]
        bw.append(360.0 if not len(below) else 0.25 * (below[0] + (G - below[-1]) - 1))
    # directivity: N^2 / sum_mn cos(k s.d_mn) sinc(k |d_mn|)
    d = p[:, None, :] - p[None, :, :]
    dist = np.linalg.norm(d, axis=2)
    s = np.column_stack([np.cos(scan[steer_i]), np.sin(scan[steer_i])])
    ph = k * np.einsum('sc,mnc->smn', s, d)
    den = np.sum(np.cos(ph) * np.sinc(k * dist / np.pi)[None], axis=(1, 2))
    D = n ** 2 / den
    # elevated source against the horizon manifold
    ce = np.cos(np.radians(el_test))
    Ae = np.exp(1j * k * ce * (np.cos(scan[steer_i])[:, None] * p[:, 0] + np.sin(scan[steer_i])[:, None] * p[:, 1]))
    Me = np.take_along_axis(np.abs(Ae.conj() @ A.T) ** 2 / n ** 2, roll, 1)
    near = np.maximum(Me[:, :41].max(1), Me[:, -40:].max(1))
    far = Me[:, 41:-40].max(1)
    return dict(psl_db=float(np.max(psl)), beamwidth_deg=float(np.mean(bw)),
                directivity_dbi=float(10 * np.log10(np.mean(D))), directivity_min_dbi=float(10 * np.log10(D.min())),
                el_loss_db=float(-10 * np.log10(near.min())), el_margin_db=float(10 * np.log10((near / far).min())))


# ---- the panel ----------------------------------------------------------------------------------------------
def check_panel(p: np.ndarray, freqs, *, sigma_cal=3.0, snr_db=20.0, el_max=15.0, az_max=60.0, trials=600,
                k_samples=64, seed=12345, step=0.025):
    """evolve_layout.check for an upright planar array and a two-dimensional scan. p: element positions in the
    panel's plane (horizontal, vertical), metres. Returns per frequency the median and 95th percentile of the
    angle between the true and estimated directions (deg) and the gross-error rate (%, over 10 deg)."""
    rng = np.random.default_rng(seed)
    n = len(p)
    c = np.arange(-1, 1 + 1e-9, step)
    uu, vv = np.meshgrid(c, c)
    m = uu ** 2 + vv ** 2 <= 1.0
    ug, vg = uu[m], vv[m]
    fo = np.arange(-step, step + 1e-9, step / 25)
    fu, fv = [a.ravel() for a in np.meshgrid(fo, fo)]
    res = []
    for f in freqs:
        k = 2 * np.pi * f / C0
        az = np.radians(rng.uniform(-az_max, az_max, trials))
        elv = np.radians(rng.uniform(0, el_max, trials))
        u, v = np.cos(elv) * np.sin(az), np.sin(elv)
        a = np.exp(1j * k * (u[:, None] * p[:, 0] + v[:, None] * p[:, 1]))
        g = 10 ** (rng.normal(0, 0.5, (trials, n)) / 20) * np.exp(1j * np.radians(rng.normal(0, sigma_cal, (trials, n))))
        s = (rng.normal(size=(trials, k_samples)) + 1j * rng.normal(size=(trials, k_samples))) / np.sqrt(2)
        w = ((rng.normal(size=(trials, n, k_samples)) + 1j * rng.normal(size=(trials, n, k_samples)))
             / np.sqrt(2) * 10 ** (-snr_db / 20))
        x = (a * g)[:, :, None] * s[:, None, :] + w
        R = x @ x.conj().transpose(0, 2, 1) / k_samples
        Ac = np.exp(1j * k * (ug[:, None] * p[:, 0] + vg[:, None] * p[:, 1]))        # (G, N)
        ue, ve = np.empty(trials), np.empty(trials)
        for i0 in range(0, trials, 100):
            Rt = R[i0:i0 + 100]
            pc = np.real(np.sum((Ac.conj()[None] @ Rt) * Ac[None], axis=2))           # (t, G)
            j = np.argmax(pc, axis=1)
            fuu, fvv = ug[j][:, None] + fu[None], vg[j][:, None] + fv[None]            # (t, 2601)
            Af = np.exp(1j * k * (fuu[..., None] * p[:, 0] + fvv[..., None] * p[:, 1]))
            pf = np.real(np.sum((Af.conj() @ Rt) * Af, axis=2))
            jf = np.argmax(pf, axis=1)
            ue[i0:i0 + 100] = fuu[np.arange(len(j)), jf]
            ve[i0:i0 + 100] = fvv[np.arange(len(j)), jf]
        r2 = ue ** 2 + ve ** 2
        sc = np.where(r2 > 1, 1 / np.sqrt(np.maximum(r2, 1e-12)), 1.0)
        ue, ve = ue * sc, ve * sc
        dot = ue * u + ve * v + np.sqrt(np.maximum(1 - ue ** 2 - ve ** 2, 0)) * np.sqrt(np.maximum(1 - u ** 2 - v ** 2, 0))
        e = np.degrees(np.arccos(np.clip(dot, -1, 1)))
        res.append((float(np.median(e)), float(np.percentile(e, 95)), float(np.mean(e > 10) * 100)))
    return np.array(res)


def panel_beam(p: np.ndarray, f: float, scan_deg: float = 0.0) -> dict:
    """Beam of a planar array of cos(theta) elements (power pattern, front hemisphere only, lossless, no
    mismatch), steered `scan_deg` off the normal in the horizontal plane: directivity of the peak, 3 dB
    beamwidth in that plane and the highest lobe outside the main lobe in that plane."""
    k = 2 * np.pi * f / C0
    th = np.radians(np.arange(0.25, 90, 0.5))
    ph = np.radians(np.arange(0, 360, 2.0))
    T, P = np.meshgrid(th, ph, indexing="ij")
    u, v = np.sin(T) * np.cos(P), np.sin(T) * np.sin(P)
    u0 = np.sin(np.radians(scan_deg))
    af = np.abs(np.exp(1j * k * ((u - u0)[..., None] * p[:, 0] + v[..., None] * p[:, 1])).sum(-1)) ** 2
    pw = af * np.cos(T)
    total = np.sum(pw * np.sin(T)) * np.radians(0.5) * np.radians(2.0)
    cut = np.radians(np.arange(-90, 90.01, 0.1))
    pc = np.abs(np.exp(1j * k * (np.sin(cut) - u0)[:, None] * p[:, 0][None]).sum(-1)) ** 2 * np.cos(cut)
    i = int(np.argmax(pc))
    lo = i
    while lo > 0 and pc[lo - 1] < pc[lo]:
        lo -= 1
    hi = i
    while hi < len(pc) - 1 and pc[hi + 1] < pc[hi]:
        hi += 1
    rest = np.concatenate([pc[:lo], pc[hi + 1:]])
    half = np.nonzero(pc >= pc[i] / 2)[0]
    seg = half[(half >= lo) & (half <= hi)]
    return dict(directivity_dbi=float(10 * np.log10(4 * np.pi * pc[i] / total)),
                beamwidth_deg=float(0.1 * (seg[-1] - seg[0] + 1)),
                sidelobe_db=float(10 * np.log10(rest.max() / pc[i])) if len(rest) else -99.0,
                peak_deg=float(np.degrees(cut[i])))


def gl_free_scan_deg(f: float, d: float) -> float:
    """Largest scan angle of a square lattice of pitch d without a grating lobe in visible space."""
    if wb is not None:
        return float(wb.max_scan_without_grating_lobes(f, d))
    s = C0 / (f * d) - 1
    return 90.0 if s >= 1 else float(np.degrees(np.arcsin(max(s, 0.0))))


# ---- windows ------------------------------------------------------------------------------------------------
def windows(freqs: np.ndarray, ok: np.ndarray) -> list:
    out, start = [], None
    for i, good in enumerate(ok):
        if good and start is None:
            start = i
        if start is not None and (not good or i == len(ok) - 1):
            j = i if good else i - 1
            out.append((float(freqs[start]), float(freqs[j])))
            start = None
    return out


def fmt_windows(ws) -> str:
    return ", ".join(f"{a / 1e9:.2f}-{b / 1e9:.2f}" for a, b in ws) or "none"


# ---- part 1 -------------------------------------------------------------------------------------------------
def part1(quick: bool, log) -> dict:
    L = load_layouts()
    freqs = np.geomspace(0.3e9, 8.0e9, 20 if quick else 80)
    trials = 200 if quick else 1000
    geoms = {
        "rings16_225_112mm": dict(p=L["rings16_225_112mm"], kind="whips", note="two rings of 8, fixed (rev B, 16 whips)"),
        "random8": dict(p=L["random8"], kind="whips", note="fixed irregular 8 (rev B, 8 whips)"),
        "ring8_112mm": dict(p=el.ring(8, 0.112), kind="whips", note="rev B middle ring of 8 (for scale)"),
        "grid16_164mm": dict(p=grid(4, C0 / 915e6 / 2), kind="whips", note="half-wave 4x4 at 915 MHz (does not fit the 500 mm plate)"),
        "grid16_95mm": dict(p=grid(4, C0 / 1575.42e6 / 2), kind="whips", note="half-wave 4x4 at 1575 MHz"),
        "panel_25mm": dict(p=grid(4, PITCH_B), kind="panel", note="rev B panel, inner 4x4 at 25 mm"),
    }
    out = dict(freqs_hz=freqs.tolist(), geometries={})
    for name, g in geoms.items():
        t0 = time.time()
        p = g["p"]
        d = dict(kind=g["kind"], note=g["note"], positions_m=np.round(p, 4).tolist())
        for cname in ("nominal", "elevated"):
            c = CONDITIONS[cname]
            r = (el.check(p, freqs, trials=trials, **c) if g["kind"] == "whips"
                 else check_panel(p, freqs, trials=min(trials, 500), **c))
            d[cname] = dict(median=r[:, 0].tolist(), p95=r[:, 1].tolist(), gross_pct=r[:, 2].tolist())
        if g["kind"] == "whips":
            m = [horizon_pattern_metrics(p, f) for f in freqs]
            for key in m[0]:
                d[key] = [x[key] for x in m]
        else:
            b0 = [panel_beam(p, f, 0.0) for f in freqs]
            b45 = [panel_beam(p, f, 45.0) for f in freqs]
            d["directivity_dbi"] = [x["directivity_dbi"] for x in b0]
            d["beamwidth_deg"] = [x["beamwidth_deg"] for x in b0]
            d["psl_db"] = [max(x["sidelobe_db"], y["sidelobe_db"]) for x, y in zip(b0, b45)]   # broadside and 45 deg scan
            d["directivity_45_dbi"] = [x["directivity_dbi"] for x in b45]
            d["gl_free_scan_deg"] = [gl_free_scan_deg(f, PITCH_B) for f in freqs]
        p95n, grn = np.array(d["nominal"]["p95"]), np.array(d["nominal"]["gross_pct"])
        gre, p95e = np.array(d["elevated"]["gross_pct"]), np.array(d["elevated"]["p95"])
        robust = (gre <= 1.0) & (grn <= 1.0)
        d["windows"] = dict(
            df_2deg=windows(freqs, robust & (p95n <= 2.0) & (p95e <= 2.5)),
            df_1deg=windows(freqs, robust & (p95n <= 1.0) & (p95e <= 1.25)),
            df_5deg=windows(freqs, robust & (p95n <= 5.0)),
            no_gross=windows(freqs, robust),
            beam=windows(freqs, (np.array(d["psl_db"]) <= -6.0) & (np.array(d["beamwidth_deg"]) <= 60.0)),
        )
        out["geometries"][name] = d
        log(f"  {name:20s} DF <=2 deg, robust: {fmt_windows(d['windows']['df_2deg'])} GHz | <=1 deg: "
            f"{fmt_windows(d['windows']['df_1deg'])} | no gross errors: {fmt_windows(d['windows']['no_gross'])} | "
            f"beams (PSL<=-6 dB, <=60 deg): {fmt_windows(d['windows']['beam'])}   [{time.time() - t0:.0f} s]")
    # exactly at the bands of interest
    tb = 400 if quick else 3000
    fc = np.array([centre(b) for b in BANDS])
    out["bands"] = {b: dict(lo_hz=BANDS[b][0], hi_hz=BANDS[b][1]) for b in BANDS}
    out["at_bands"] = {}
    for name, g in geoms.items():
        row = {}
        for cname in ("nominal", "elevated", "stress"):
            c = CONDITIONS[cname]
            r = (el.check(g["p"], fc, trials=tb, **c) if g["kind"] == "whips"
                 else check_panel(g["p"], fc, trials=min(tb, 800), **c))
            row[cname] = {b: dict(median=float(r[i, 0]), p95=float(r[i, 1]), gross_pct=float(r[i, 2])) for i, b in enumerate(BANDS)}
        out["at_bands"][name] = row
    return out


# ---- part 2 -------------------------------------------------------------------------------------------------
def make_band_score(n: int, bands: list, *, n_freq=12, trials=96, sigma_cal=4.0, snr_db=10.0, k_samples=64,
                    el_max=el.EL_MAX, seed=0):
    """evolve_layout.make_score with the frequencies drawn inside the chosen bands instead of across 0.4-2.3 GHz
    (n_freq draws shared equally among the bands, uniform inside each); everything else is the same."""
    snr_eff = snr_db + 10 * np.log10(k_samples)
    per = max(1, n_freq // len(bands))
    lo = np.repeat([BANDS[b][0] for b in bands], per)
    hi = np.repeat([BANDS[b][1] for b in bands], per)
    nf = len(lo)

    def score(X: np.ndarray, gen: int) -> np.ndarray:
        rng = np.random.default_rng([seed, gen])
        f = lo + (hi - lo) * rng.random(nf)
        k = 2 * np.pi * f / C0
        az = rng.uniform(0, 2 * np.pi, (nf, trials))
        ce = np.cos(np.radians(rng.uniform(0, el_max, (nf, trials))))
        cal = (10 ** (rng.normal(0, 0.5, (nf, trials, n)) / 20)
               * np.exp(1j * np.radians(rng.normal(0, sigma_cal, (nf, trials, n)))))
        noise = ((rng.normal(size=(nf, trials, n)) + 1j * rng.normal(size=(nf, trials, n)))
                 / np.sqrt(2) * 10 ** (-snr_eff / 20))
        out = np.empty(len(X))
        for m, x in enumerate(X):
            p = x.reshape(n, 2)
            proj_t = (p[:, 0] * np.cos(az)[..., None] + p[:, 1] * np.sin(az)[..., None]) * ce[..., None]
            y = np.exp(1j * k[:, None, None] * proj_t) * cal + noise
            proj_g = p[:, 0] * np.cos(el.GRID)[:, None] + p[:, 1] * np.sin(el.GRID)[:, None]
            A = np.exp(-1j * k[:, None, None] * proj_g[None])
            pw = (np.abs(A @ y.transpose(0, 2, 1)) ** 2).transpose(0, 2, 1)
            est = el._refine(pw, np.argmax(pw, axis=-1))
            err = np.degrees(np.angle(np.exp(1j * (est - az))))
            loss = np.sqrt(np.mean(np.minimum(np.abs(err), 30.0) ** 2, axis=1))
            out[m] = loss.max() + 0.25 * loss.mean() + 100.0 * el.violation(p)
        return out

    return score


def search_bands(n: int, bands: list, ref: np.ndarray, budget: int, draws: int, seed: int,
                 methods=("random", "climbed", "polished")):
    """Searches for one band set; each method's finalists are re-scored on fresh draws (reselect).
    random: random search, `budget` evaluations. climbed: hill climb from the full-band layout `ref`,
    0.75 budget. polished: hill climb from the random search's re-selected winner, 0.25 budget (so
    random + polished is one search of 1.25 budget). evolved: the genetic algorithm, `budget` (optional; a first
    run with it found nothing the others did not, see the log of the first attempt in the report).
    Returns {method: (positions, fresh score)} and the score function."""
    score = make_band_score(n, bands)
    lo, hi = -el.R_OUT * np.ones(2 * n), el.R_OUT * np.ones(2 * n)
    pop = 40
    out = {}
    for name in methods:
        if name == "random":
            res = random_search(score, lo, hi, batches=budget // pop, batch=pop, seed=seed, repair=el.repair)
        elif name == "evolved":
            res = evolve(score, lo, hi, pop=pop, generations=budget // pop - 1, elite=2, sigma=0.15, crossover=0.5,
                         seed=seed, repair=el.repair)
        elif name == "climbed":
            res = climb(score, el.repair(ref.reshape(-1)), lo, hi, steps=int(0.75 * budget) // 9, lam=8, sigma=0.05,
                        seed=seed, repair=el.repair)
        else:
            res = climb(score, el.repair(out["random"][0].reshape(-1)), lo, hi, steps=int(0.25 * budget) // 9, lam=8,
                        sigma=0.03, seed=seed + 1, repair=el.repair)
        X, F = reselect(score, res.finalists[:10], draws=draws)
        out[name] = (X[0].reshape(n, 2), float(F[0]))
    return out, score


def band_check(p: np.ndarray, bands: list, trials: int) -> dict:
    """evolve_layout.check at the bottom, middle and top of each band; per band and condition the worst of the
    three: 95th percentile, gross-error rate, and the median at the centre."""
    f = np.concatenate([[BANDS[b][0], centre(b), BANDS[b][1]] for b in bands])
    out = {}
    for cname, c in CONDITIONS.items():
        r = el.check(p, f, trials=trials, **c).reshape(len(bands), 3, 3)
        out[cname] = {b: dict(median=float(r[i, 1, 0]), p95=float(r[i, :, 1].max()), gross_pct=float(r[i, :, 2].max()))
                      for i, b in enumerate(bands)}
    return out


def worst(bc: dict, cname: str, key: str) -> float:
    return max(v[key] for v in bc[cname].values())


def part2(quick: bool, log, prev: dict | None = None, save=None) -> dict:
    L = load_layouts()
    budget = 400 if quick else 6000
    draws = 3 if quick else 12
    trials = 300 if quick else 3000
    full = np.arange(0.40e9, 2.3001e9, 0.10e9)
    refs = {16: "rings16_225_112mm", 8: "random8"}
    baselines = {16: ["rings16_225_112mm", "evolved16", "random16", "ring16_225mm"],
                 8: ["random8", "evolved8", "climbed8", "ring8_225mm", "ring8_112mm", "ring8_78mm"]}
    extra = {16: {"ring16_240mm": el.ring(16, 0.240)}, 8: {"ring8_240mm": el.ring(8, 0.240)}}
    out = dict(budget=budget, sets={}, single={}, relaxed={})
    for key in ('sets', 'single', 'relaxed'):          # finished searches of an interrupted run are kept
        out[key].update((prev or {}).get(key, {}) if (prev or {}).get('budget') == budget else {})

    def full_band(p):
        r = {c: el.check(p, full, trials=300 if quick else 800, **CONDITIONS[c]) for c in ("nominal", "elevated")}
        return dict(p95_nominal_max=float(r["nominal"][:, 1].max()), gross_elevated_max=float(r["elevated"][:, 2].max()))

    def one(n, bands, tag, seed, methods=("random", "climbed", "polished"), with_baselines=True, scale=1.0):
        t0 = time.time()
        ref = L[refs[n]]
        found, score = search_bands(n, bands, ref, int(budget * scale), draws, seed, methods)
        cands = {}
        if with_baselines:
            cands.update({k: L[k] for k in baselines[n]})
            cands.update(extra[n])
        cands.update({f"cherry_{m}": p for m, (p, _) in found.items()})
        names = list(cands)
        X = np.array([cands[k].reshape(-1) for k in names])
        F = np.mean([score(X, 2_000_000 + j) for j in range(draws)], axis=0)       # same fresh draws for all
        rows = {}
        for k, fs in zip(names, F):
            p = cands[k]
            rows[k] = dict(search_score=float(fs), violation=float(el.violation(p)), positions_m=np.round(p, 4).tolist(),
                           check=band_check(p, bands, trials))
            if k.startswith("cherry_"):
                rows[k]["full_band_0p4_2p3"] = full_band(p)
        best = min((k for k in names if k.startswith("cherry_")), key=lambda k: rows[k]["search_score"])
        log(f"  {tag} N={n}: " + " | ".join(
            f"{k} score {rows[k]['search_score']:.2f}, stress p95 {worst(rows[k]['check'], 'stress', 'p95'):.2f}, "
            f"gross(el) {worst(rows[k]['check'], 'elevated', 'gross_pct'):.1f} %" for k in names)
            + f"   [{time.time() - t0:.0f} s]")
        return dict(bands=bands, n=n, reference=refs[n], best_cherry=best, layouts=rows)

    for tag, bands in BAND_SETS.items():
        for n in ((8,) if quick else (16, 8)):
            if f"{tag}_N{n}" not in out["sets"]:
                out["sets"][f"{tag}_N{n}"] = one(n, bands, tag, seed=100 + n + len(tag))
                save and save(out)
    if not quick:
        # one band at a time: the most any cherry-picked set could reach in that band (16 whips)
        for b in SINGLE_BANDS:
            if b not in out["single"]:
                out["single"][b] = one(16, [b], f"single_{b}", seed=300 + len(b) + int(b) % 97,
                                       methods=("random", "polished"), scale=0.5)
                save and save(out)
        # the spacing rule is 0.1 wavelength at the lowest frequency: relax it to the set's lowest band (floor 30 mm)
        s_keep = el.S_MIN
        try:
            for tag in ("B_915_1575_2440", "C_2440_5800"):
                if tag in out["relaxed"]:
                    continue
                bands = BAND_SETS[tag]
                el.S_MIN = max(0.030, 0.1 * C0 / min(BANDS[b][0] for b in bands))
                r = one(16, bands, f"relaxed_{tag}", seed=500 + len(tag), methods=("random", "polished"),
                        with_baselines=False, scale=0.75)
                r["min_spacing_m"] = el.S_MIN
                out["relaxed"][tag] = r
                save and save(out)
        finally:
            el.S_MIN = s_keep
    return out


# ---- part 3 -------------------------------------------------------------------------------------------------
def part3(quick: bool, log) -> dict:
    lam59, lam25 = C0 / 5.9e9, C0 / 2.5e9
    s60 = 1 + np.sin(np.radians(60))
    pitches = {
        "25.0 mm (rev B, half-wave at 6.0 GHz)": 0.025,
        "25.4 mm (half-wave at 5.9 GHz)": lam59 / 2,
        "27.2 mm (no grating lobe to +/-60 deg at 5.9 GHz)": lam59 / s60,
        "60.0 mm (half-wave at 2.5 GHz)": lam25 / 2,
        "64.3 mm (no grating lobe to +/-60 deg at 2.5 GHz)": lam25 / s60,
    }
    freqs = np.array([2.44e9, 3.5e9, 5.8e9])
    trials = 200 if quick else 800
    out = dict(freqs_hz=freqs.tolist(), pitches={})
    for name, d in pitches.items():
        p = grid(4, d)
        row = dict(pitch_m=float(d), active_aperture_m=float(4 * d), tile_6x6_m=float(6 * d),
                   unambiguous_to_hz=float(C0 / (2 * d)), gl_free_60deg_to_hz=float(C0 / (d * s60)))
        for cname in ("nominal", "elevated", "stress"):
            r = check_panel(p, freqs, trials=trials, **CONDITIONS[cname])
            row[cname] = dict(median=r[:, 0].tolist(), p95=r[:, 1].tolist(), gross_pct=r[:, 2].tolist())
        b0 = [panel_beam(p, f, 0.0) for f in freqs]
        b45 = [panel_beam(p, f, 45.0) for f in freqs]
        row["directivity_dbi"] = [x["directivity_dbi"] for x in b0]
        row["beamwidth_deg"] = [x["beamwidth_deg"] for x in b0]
        row["directivity_45_dbi"] = [x["directivity_dbi"] for x in b45]
        row["sidelobe_45_db"] = [x["sidelobe_db"] for x in b45]
        row["gl_free_scan_deg"] = [gl_free_scan_deg(f, d) for f in freqs]
        out["pitches"][name] = row
        log(f"  {name:52s} " + " | ".join(
            f"{f / 1e9:.2f} GHz: D {row['directivity_dbi'][i]:.1f} dBi, {row['beamwidth_deg'][i]:.0f} deg, "
            f"p95 {row['nominal']['p95'][i]:.2f}, gross(el) {row['elevated']['gross_pct'][i]:.1f} %" for i, f in enumerate(freqs)))
    return out


# ---- tables and figures ---------------------------------------------------------------------------------------
def tables(res: dict) -> str:
    T = []
    if "part1" in res:
        P = res["part1"]
        T.append("PART 1. Windows (GHz) where each geometry is good, 0.3-8 GHz. DF: 95th percentile (0-15 deg up) within the limit,\n"
                 "gross errors <= 1 % for sources 0-15 and 0-30 deg up. Beams: peak sidelobe <= -6 dB and beamwidth <= 60 deg.")
        T.append(f"{'geometry':20s} | {'DF <= 1 deg':24s} | {'DF <= 2 deg':24s} | {'DF <= 5 deg':24s} | {'no gross errors':24s} | beams")
        for k, g in P["geometries"].items():
            w = g["windows"]
            T.append(f"{k:20s} | {fmt_windows(w['df_1deg']):24s} | {fmt_windows(w['df_2deg']):24s} | {fmt_windows(w['df_5deg']):24s} | "
                     f"{fmt_windows(w['no_gross']):24s} | {fmt_windows(w['beam'])}")
        T.append("\nAt the bands of interest: 95th percentile error, deg (0-15 deg up) / gross errors, % (0-30 deg up) / stress 95th percentile")
        names = list(P["at_bands"])
        T.append(f"{'band':6s} | " + " | ".join(f"{k:22s}" for k in names))
        for b in P["bands"]:
            T.append(f"{b:6s} | " + " | ".join(
                f"{P['at_bands'][k]['nominal'][b]['p95']:6.2f} /{P['at_bands'][k]['elevated'][b]['gross_pct']:5.1f} /{P['at_bands'][k]['stress'][b]['p95']:7.2f}"
                for k in names))
        f = np.array(P["freqs_hz"])
        T.append("\nBeam figures (whips: horizon beam, isotropic elements; panel: broadside, cos element): directivity dBi / beamwidth deg / peak sidelobe dB")
        sel = [0.433e9, 0.915e9, 1.575e9, 2.44e9, 3.5e9, 5.8e9]
        T.append(f"{'geometry':20s} | " + " | ".join(f"{x / 1e9:5.3f} GHz        " for x in sel))
        for k, g in P["geometries"].items():
            cells = []
            for x in sel:
                i = int(np.argmin(np.abs(np.log(f / x))))
                cells.append(f"{g['directivity_dbi'][i]:5.1f} /{g['beamwidth_deg'][i]:6.1f} /{g['psl_db'][i]:6.1f}")
            T.append(f"{k:20s} | " + " | ".join(cells))
    if "part2" in res:
        P = res["part2"]
        T.append("\nPART 2. Whip layouts optimised for a band set against the full-band layouts, at those bands (worst of the bands).\n"
                 "score = the search's own score on fresh draws (lower is better); p95 in deg, gross in %; nominal / elevated / stress.")
        for group in ("sets", "single", "relaxed"):
            for tag, S in P.get(group, {}).items():
                T.append(f"\n[{group}] {tag}: bands {S['bands']}, N = {S['n']}, best cherry-picked = {S['best_cherry']}"
                         + (f", min spacing {S['min_spacing_m'] * 1000:.0f} mm" if 'min_spacing_m' in S else ""))
                T.append(f"  {'layout':20s} | score | p95 nom | p95 elev | p95 stress | gross nom | gross elev | gross stress | per band p95 nominal / stress")
                for k, r in S["layouts"].items():
                    c = r["check"]
                    per = "; ".join(f"{b}: {c['nominal'][b]['p95']:.2f}/{c['stress'][b]['p95']:.2f}" for b in S["bands"])
                    fb = r.get("full_band_0p4_2p3")
                    T.append(f"  {k:20s} | {r['search_score']:5.2f} | {worst(c, 'nominal', 'p95'):7.2f} | {worst(c, 'elevated', 'p95'):8.2f} | "
                             f"{worst(c, 'stress', 'p95'):10.2f} | {worst(c, 'nominal', 'gross_pct'):9.1f} | {worst(c, 'elevated', 'gross_pct'):10.1f} | "
                             f"{worst(c, 'stress', 'gross_pct'):12.1f} | {per}"
                             + (f" | over 0.4-2.3 GHz: p95 {fb['p95_nominal_max']:.2f}, gross(el) {fb['gross_elevated_max']:.1f} %" if fb else ""))
    if "part3" in res:
        P = res["part3"]
        fs = P["freqs_hz"]
        T.append("\nPART 3. Panel pitch (4x4 active). Per frequency: directivity dBi (broadside / 45 deg scan) | beamwidth deg | "
                 "DF 95th pct deg nominal / stress | gross % elevated | grating-lobe-free scan deg")
        for name, r in P["pitches"].items():
            T.append(f"\n  {name}: active aperture {r['active_aperture_m'] * 1000:.0f} mm, 6x6 tile {r['tile_6x6_m'] * 1000:.0f} mm, "
                     f"half-wave at {r['unambiguous_to_hz'] / 1e9:.2f} GHz, grating-lobe-free to +/-60 deg up to {r['gl_free_60deg_to_hz'] / 1e9:.2f} GHz")
            for i, f in enumerate(fs):
                T.append(f"    {f / 1e9:.2f} GHz: {r['directivity_dbi'][i]:5.1f} / {r['directivity_45_dbi'][i]:5.1f} | {r['beamwidth_deg'][i]:5.1f} | "
                         f"{r['nominal']['p95'][i]:6.2f} / {r['stress']['p95'][i]:6.2f} | {r['elevated']['gross_pct'][i]:5.1f} | {r['gl_free_scan_deg'][i]:4.0f}")
    return "\n".join(T) + "\n"


def plot(res: dict, sfx: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    marks = [centre(b) / 1e9 for b in ("433", "915", "1575", "2440", "3600", "5800")]
    if "part1" in res:
        P = res["part1"]
        f = np.array(P["freqs_hz"]) / 1e9
        fig, axes = plt.subplots(2, 3, figsize=(16, 8.5))
        panels = [("nominal", "p95", "95th percentile error, sources 0-15 deg up (deg)", True),
                  ("elevated", "gross_pct", "gross errors (>10 deg), sources 0-30 deg up (%)", False),
                  (None, "psl_db", "peak sidelobe, worst steer (dB; 0 = ambiguous)", False),
                  (None, "directivity_dbi", "directivity (dBi; whips isotropic, panel cos element)", False),
                  (None, "beamwidth_deg", "3 dB beamwidth (deg)", True)]
        for ax, (cond, key, title, logy) in zip(axes.ravel(), panels):
            for k, g in P["geometries"].items():
                y = np.array(g[cond][key] if cond else g[key], float)
                if key == "psl_db":
                    y = np.maximum(y, -30)
                (ax.loglog if logy else ax.semilogx)(f, y, label=k, lw=2 if k in ("rings16_225_112mm", "panel_25mm") else 1.1,
                                                     ls="--" if g["kind"] == "panel" else "-")
            for m in marks:
                ax.axvline(m, color="0.8", lw=0.7, zorder=0)
            ax.set_title(title, fontsize=9); ax.set_xlabel("GHz"); ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=6.5)
        ax = axes.ravel()[5]
        names = list(P["geometries"])
        for i, k in enumerate(names):
            w = P["geometries"][k]["windows"]
            for (a, b) in w["df_5deg"]:
                ax.barh(i, (b - a) / 1e9, left=a / 1e9, height=0.7, color="C0", alpha=0.25)
            for (a, b) in w["df_2deg"]:
                ax.barh(i, (b - a) / 1e9, left=a / 1e9, height=0.5, color="C0", alpha=0.55)
            for (a, b) in w["df_1deg"]:
                ax.barh(i, (b - a) / 1e9, left=a / 1e9, height=0.3, color="C0")
            for (a, b) in w["beam"]:
                ax.barh(i - 0.42, (b - a) / 1e9, left=a / 1e9, height=0.1, color="C3")
        for m in marks:
            ax.axvline(m, color="0.5", lw=0.7)
        ax.set_xscale("log"); ax.set_xlim(0.3, 8); ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=7)
        ax.set_title("windows: DF <=5 / <=2 / <=1 deg and robust (blue, light to dark); beams (red)", fontsize=9); ax.set_xlabel("GHz")
        fig.tight_layout(); fig.savefig(OUT / f"fom_vs_frequency{sfx}.png", dpi=120); plt.close(fig)
    if "part2" in res:
        S = res["part2"]["sets"]
        tags = list(S)
        fig, axes = plt.subplots(2, len(tags), figsize=(3.1 * len(tags), 6.4), squeeze=False)
        for j, tag in enumerate(tags):
            s = S[tag]
            ref, best = s["reference"], s["best_cherry"]
            for ax, k, col in ((axes[0, j], ref, "C0"), (axes[1, j], best, "C3")):
                p = np.array(s["layouts"][k]["positions_m"]) * 1000
                ax.add_patch(plt.Circle((0, 0), 250, fill=False, color="0.6"))
                ax.add_patch(plt.Circle((0, 0), el.R_IN * 1000, fill=False, ls=":", color="0.6"))
                ax.plot(p[:, 0], p[:, 1], "o", color=col, ms=4); ax.plot(0, 0, "+k")
                ax.set_aspect("equal"); ax.set_xlim(-260, 260); ax.set_ylim(-260, 260); ax.tick_params(labelsize=6)
                c = s["layouts"][k]["check"]
                ax.set_title(f"{tag}\n{k}: stress p95 {worst(c, 'stress', 'p95'):.2f} deg", fontsize=7)
        fig.suptitle("Top: full-band layout. Bottom: best layout optimised for the band set (mm on the 500 mm plate)", fontsize=9)
        fig.tight_layout(); fig.savefig(OUT / f"cherry_layouts{sfx}.png", dpi=120); plt.close(fig)
    if "part3" in res:
        P = res["part3"]
        fs = np.array(P["freqs_hz"]) / 1e9
        fig, axes = plt.subplots(1, 3, figsize=(14, 4))
        w = 0.8 / len(P["pitches"])
        for i, (name, r) in enumerate(P["pitches"].items()):
            x = np.arange(len(fs)) + i * w
            axes[0].bar(x, r["directivity_dbi"], w, label=name.split(" (")[0])
            axes[1].bar(x, r["nominal"]["p95"], w)
            axes[2].bar(x, r["elevated"]["gross_pct"], w)
        for ax, t in zip(axes, ("broadside directivity (dBi)", "DF 95th percentile, nominal (deg)", "gross errors, sources 0-30 deg up (%)")):
            ax.set_xticks(np.arange(len(fs)) + 0.4 - w / 2); ax.set_xticklabels([f"{x:.2f} GHz" for x in fs]); ax.set_title(t, fontsize=9)
            ax.grid(alpha=0.3, axis="y")
        axes[1].set_yscale("log"); axes[0].legend(fontsize=7)
        fig.tight_layout(); fig.savefig(OUT / f"panel_pitch{sfx}.png", dpi=120); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--parts", default="1,2,3")
    ap.add_argument("--plot-only", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    sfx = "_quick" if args.quick else ""
    path = OUT / f"results{sfx}.json"
    res = json.loads(path.read_text()) if path.exists() else {}
    if not args.plot_only:
        t0 = time.time()

        def log(s):
            print(s, flush=True)

        for part, fn in (("1", part1), ("2", part2), ("3", part3)):
            if part in args.parts.split(","):
                log(f"PART {part}")
                if part == "2":
                    def save(o):
                        res["part2"] = o
                        path.write_text(json.dumps(res, indent=1))
                    res["part2"] = fn(args.quick, log, res.get("part2"), save)
                else:
                    res[f"part{part}"] = fn(args.quick, log)
                path.write_text(json.dumps(res, indent=1))
        log(f"total {time.time() - t0:.0f} s")
    (OUT / f"tables{sfx}.txt").write_text(tables(res))
    plot(res, sfx)
    print(tables(res))


if __name__ == "__main__":
    main()
