"""Is a 3-D arrangement of the antennas worth it? A simulation study for rev B (REDESIGN.md, sections 4 and 9).

Part 1, external whips (0.4-2.3 GHz, 16 identical vertical whips/dipoles, amplitude pattern cos(elevation)):
  (a) the planar two rings of 8 (r = 225 / 112 mm) with the azimuth-only estimator (assumes elevation 0),
  (b) the same planar array with a joint azimuth + elevation estimator (for a planar array this is the
      search over the direction cosines u, v: its beam pattern depends on u, v only),
  (c) simple 3-D shapes inside 500 x 500 x 300 mm (inner ring raised, stacked rings, ring + mast, cone, helix),
  (d) 3-D layouts found by array16.evolve (evolve, climb, random_search at the same budget, finalists
      re-scored with reselect), minimum spacing 75 mm; and the best planar layout found the same way.
Sources 0-60 deg above the horizon, reported in elevation bins 0-10, 10-30, 30-60 deg.

Conventions are those of experiments/evolve_layout.py (its `check`): per-whip calibration residual 3 deg
and 0.5 dB, 20 dB SNR per whip (here: for a source on the horizon; the cos(elevation) pattern lowers it for
elevated sources, 6 dB at 60 deg), K = 64 samples, sample covariance, Bartlett beamscan, gross error =
azimuth off by more than 10 deg. The joint estimator scans a 2 deg azimuth x elevation grid, then refines
on 0.5, 0.1 and 0.025 deg grids. 3-D arrays are searched from -20 deg elevation (they can tell up from
down), planar ones from 0 (they cannot).

Model limits: ideal identical elements, no mutual coupling, no ground-plate image or edge diffraction, no
mast or feed-line scattering, one source, no multipath. A raised whip is no longer a monopole on the plate:
it has to be a sleeve dipole or carry its own small ground plane, and its feed line runs through the array.

Part 2, panel band (2.3-6 GHz): one 4x4 tile (25 mm pitch, element power pattern cos^1.3 of the angle off
boresight, nothing behind the ground plane) against 3 or 4 tiles on the faces of a vertical prism and one
tile on a rotator: gain and azimuth DF accuracy against azimuth, elements used per direction.

    python experiments/array3d_study.py            # about 25 minutes with OMP_NUM_THREADS=4
    python experiments/array3d_study.py --quick    # a smoke run
    python experiments/array3d_study.py --plot-only
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
from array16.evolve import climb, evolve, random_search, reselect  # noqa: E402
import evolve_layout as planar  # noqa: E402  (read-only: ring, check, constants)

C0 = planar.C0
R_OUT, S_MIN, HEIGHT = planar.R_OUT, planar.S_MIN, 0.300
BINS = [(0.0, 10.0), (10.0, 30.0), (30.0, 60.0)]
BANDS = {"0.4-1 GHz": (0.39e9, 1.001e9), "1-2.3 GHz": (1.001e9, 2.31e9)}
OUT = HERE / "array3d"


def direction(az, el):
    return np.stack([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)], axis=-1)


# ---- layouts -----------------------------------------------------------------------------------------------
def ring3(n, radius, z=0.0, rot=0.0):
    return np.column_stack([planar.ring(n, radius, rot), np.full(n, z)])


def hand_designs() -> dict:
    d = {}
    d["planar rings 225/112"] = np.vstack([ring3(8, 0.225), ring3(8, 0.112, 0, np.pi / 8)])
    d["inner ring +150 mm"] = np.vstack([ring3(8, 0.225), ring3(8, 0.112, 0.150, np.pi / 8)])
    d["inner ring +300 mm"] = np.vstack([ring3(8, 0.225), ring3(8, 0.112, 0.300, np.pi / 8)])
    d["stacked rings 225/225, 300 mm apart"] = np.vstack([ring3(8, 0.225), ring3(8, 0.225, 0.300, np.pi / 8)])
    d["ring of 12 + mast of 4"] = np.vstack([ring3(12, 0.225), np.column_stack(
        [np.zeros(4), np.zeros(4), [0.075, 0.150, 0.225, 0.300]])])
    d["cone 8/5/3"] = np.vstack([ring3(8, 0.225), ring3(5, 0.140, 0.150, 0.3), ring3(3, 0.060, 0.300, 0.7)])
    i = np.arange(16)
    b = 2 * np.pi * i * 2 / 16 + 0.1 * i                       # two turns, slightly irregular
    d["helix r 225, 2 turns"] = np.column_stack([0.225 * np.cos(b), 0.225 * np.sin(b), HEIGHT * i / 15])
    return d


def repair3(x: np.ndarray, flat: bool = False) -> np.ndarray:
    """Into the cylinder (r <= 240 mm, 0 <= z <= 300 mm), pairs pushed apart to S_MIN in 3-D, ordered by angle."""
    p = x.reshape(-1, 3).copy()
    for _ in range(60):
        r = np.hypot(p[:, 0], p[:, 1])
        p[:, :2] *= (np.minimum(r, R_OUT) / np.maximum(r, 1e-12))[:, None]
        p[:, 2] = 0.0 if flat else np.clip(p[:, 2], 0.0, HEIGHT)
        dist = np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)
        np.fill_diagonal(dist, np.inf)
        i, j = np.nonzero(np.triu(dist < S_MIN - 1e-9, 1))
        if not len(i):
            break
        for a, c in zip(i, j):
            v = p[a] - p[c]
            if flat:
                v[2] = 0.0
            n = np.linalg.norm(v)
            v = v / n if n > 1e-12 else np.array([1.0, 0.0, 0.0])
            push = (S_MIN - n) / 2 + 1e-6
            p[a] += push * v
            p[c] -= push * v
    p = p[np.argsort(np.arctan2(p[:, 1], p[:, 0]))]
    return p.reshape(-1)


def violation3(p: np.ndarray) -> float:
    r = np.hypot(p[:, 0], p[:, 1])
    d = np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)
    np.fill_diagonal(d, np.inf)
    return float(np.sum(np.maximum(r - R_OUT, 0)) + np.sum(np.maximum(-p[:, 2], 0))
                 + np.sum(np.maximum(p[:, 2] - HEIGHT, 0)) + np.sum(np.triu(np.maximum(S_MIN - d, 0), 1)))


# ---- covariance model (final numbers) ---------------------------------------------------------------------
def covariances(p, k, az, el, rng, *, sigma_cal=3.0, snr_db=20.0, k_samples=64, pattern=True, amp=None):
    """Sample covariances (T, N, N) for one source per trial; same noise model as evolve_layout.check."""
    t, n = len(az), len(p)
    a = np.exp(1j * k * (direction(az, el) @ p.T))
    if amp is not None:
        a = a * amp
    elif pattern:
        a = a * np.cos(el)[:, None]
    g = 10 ** (rng.normal(0, 0.5, (t, n)) / 20) * np.exp(1j * np.radians(rng.normal(0, sigma_cal, (t, n))))
    s = (rng.normal(size=(t, k_samples)) + 1j * rng.normal(size=(t, k_samples))) / np.sqrt(2)
    w = ((rng.normal(size=(t, n, k_samples)) + 1j * rng.normal(size=(t, n, k_samples)))
         / np.sqrt(2) * 10 ** (-snr_db / 20))
    x = (a * g)[:, :, None] * s[:, None, :] + w
    return x @ x.conj().transpose(0, 2, 1) / k_samples


def est_az_only(R, p, k):
    """Azimuth-only Bartlett, elevation assumed 0, heights ignored (as evolve_layout.check): 1 deg then 0.05 deg."""
    t = len(R)
    coarse = np.radians(np.arange(0, 360, 1.0))
    Ac = np.exp(1j * k * (np.cos(coarse)[:, None] * p[:, 0] + np.sin(coarse)[:, None] * p[:, 1]))
    pc = np.real(np.einsum('gi,tij,gj->tg', Ac.conj(), R, Ac))
    c0 = coarse[np.argmax(pc, axis=1)]
    fine = c0[:, None] + np.radians(np.arange(-1.5, 1.5001, 0.05))[None, :]
    Af = np.exp(1j * k * (np.cos(fine)[..., None] * p[:, 0] + np.sin(fine)[..., None] * p[:, 1]))
    pf = np.real(np.einsum('tgi,tij,tgj->tg', Af.conj(), R, Af))
    return fine[np.arange(t), np.argmax(pf, axis=1)]


def est_joint(R, p, k, el_lo=0.0, step=2.0, chunk=100):
    """Joint azimuth/elevation Bartlett: `step` deg grid, then 0.5, 0.1 and 0.025 deg local grids."""
    t = len(R)
    azg = np.radians(np.arange(0, 360, step))
    elg = np.radians(np.arange(el_lo, 88.001, step))
    AZ, EL = np.meshgrid(azg, elg, indexing="ij")
    A = np.exp(1j * k * (direction(AZ.ravel(), EL.ravel()) @ p.T)).astype(np.complex64)    # (G, N)
    Rc = R.astype(np.complex64)
    idx = np.empty(t, int)
    for s in range(0, t, chunk):
        B = A.conj()[None] @ Rc[s:s + chunk]                                                # (c, G, N)
        idx[s:s + chunk] = np.argmax(np.einsum('cgn,gn->cg', B, A).real, axis=1)
    az, el = AZ.ravel()[idx], EL.ravel()[idx]
    lo = np.radians(el_lo)
    for half, st in ((step, 0.5), (0.6, 0.1), (0.1, 0.025)):
        o = np.radians(np.arange(-half, half + 1e-9, st))
        DA, DE = np.meshgrid(o, o, indexing="ij")
        azf = az[:, None] + DA.ravel()[None]
        elf = np.clip(el[:, None] + DE.ravel()[None], lo, np.pi / 2)
        Af = np.exp(1j * k * (direction(azf, elf) @ p.T))                                   # (T, g, N)
        pf = np.real(np.einsum('tgi,tij,tgj->tg', Af.conj(), R, Af))
        j = np.argmax(pf, axis=1)
        az, el = azf[np.arange(t), j], elf[np.arange(t), j]
    return az, el


def wrap_deg(x):
    return np.abs(np.degrees(np.angle(np.exp(1j * x))))


def monte_carlo(p, freqs, *, trials=300, seed=2026, estimators=("az", "joint"), bins=BINS, pattern=True,
                sigma_cal=3.0, snr_db=20.0):
    """Errors (deg) for every frequency, elevation bin and trial. Returns dict of arrays (F, B, T)."""
    rng = np.random.default_rng(seed)
    flat = np.ptp(p[:, 2]) < 1e-9
    out = {k: np.full((len(freqs), len(bins), trials), np.nan) for k in ("az_only", "joint_az", "joint_el", "el_true")}
    for i, f in enumerate(freqs):
        k = 2 * np.pi * f / C0
        az = rng.uniform(0, 2 * np.pi, (len(bins), trials))
        el = np.radians(np.array([rng.uniform(lo, hi, trials) for lo, hi in bins]))
        R = covariances(p, k, az.ravel(), el.ravel(), rng, pattern=pattern, sigma_cal=sigma_cal, snr_db=snr_db)
        out["el_true"][i] = np.degrees(el)
        if "az" in estimators:
            out["az_only"][i] = wrap_deg(est_az_only(R, p, k) - az.ravel()).reshape(az.shape)
        if "joint" in estimators:
            a, e = est_joint(R, p, k, el_lo=0.0 if flat else -20.0)
            out["joint_az"][i] = wrap_deg(a - az.ravel()).reshape(az.shape)
            out["joint_el"][i] = np.abs(np.degrees(e) - np.degrees(el.ravel())).reshape(az.shape)
    return out


def stats(err, freqs):
    """Summaries per band and elevation bin, plus per-frequency curves."""
    freqs = np.asarray(freqs)
    res = {"per_freq": {}, "bands": {}}
    for key in ("az_only", "joint_az", "joint_el"):
        e = err[key]
        if np.isnan(e).all():
            continue
        res["per_freq"][key] = dict(median=np.median(e, axis=2).round(3).tolist(),
                                    p95=np.percentile(e, 95, axis=2).round(3).tolist(),
                                    gross_pct=(np.mean(e > 10, axis=2) * 100).round(2).tolist())
        for bname, (lo, hi) in {**BANDS, "0.4-2.3 GHz": (0, 1e12)}.items():
            m = (freqs >= lo) & (freqs <= hi)
            if not m.any():
                continue
            rows = []
            for b in range(e.shape[1]):
                x = e[m, b]
                rows.append(dict(median=float(np.median(x)), p95=float(np.percentile(x, 95)),
                                 gross_pct=float(np.mean(x > 10) * 100),
                                 gross_pct_worst_freq=float(np.mean(x > 10, axis=1).max() * 100)))
            res["bands"].setdefault(bname, {})[key] = rows
    return res


# ---- search model (fast, single snapshot, joint azimuth/elevation) ----------------------------------------
def make_score3(n, *, n_freq=6, trials=40, sigma_cal=4.0, snr_db=10.0, k_samples=64, el_max=60.0, seed=0, step=4.0):
    """Score of a batch of 3-D layouts: per frequency the RMS of the angle between the true and the estimated
    direction (capped at 30 deg), sources 0-el_max deg up; worst frequency plus a quarter of the mean."""
    snr_eff = snr_db + 10 * np.log10(k_samples)
    edges = np.geomspace(*planar.BAND, n_freq + 1)
    azg = np.radians(np.arange(0, 360, step))
    elg = np.radians(np.arange(0, 88.001, step))
    AZ, EL = np.meshgrid(azg, elg, indexing="ij")
    Dg = direction(AZ.ravel(), EL.ravel())                                # (G, 3)
    na, ne = len(azg), len(elg)
    st = np.radians(step)

    def parab(y0, y1, y2):
        den = y0 - 2 * y1 + y2
        ok = np.abs(den) > 1e-30
        return np.clip(np.where(ok, 0.5 * (y0 - y2) / np.where(ok, den, 1.0), 0.0), -0.5, 0.5)

    def score(X: np.ndarray, gen: int) -> np.ndarray:
        rng = np.random.default_rng([seed, gen])
        f = edges[:-1] * (edges[1:] / edges[:-1]) ** rng.random(n_freq)
        k = 2 * np.pi * f / C0
        az = rng.uniform(0, 2 * np.pi, (n_freq, trials))
        el = np.radians(rng.uniform(0, el_max, (n_freq, trials)))
        d = direction(az, el)                                             # (F, T, 3)
        cal = (10 ** (rng.normal(0, 0.5, (n_freq, trials, n)) / 20)
               * np.exp(1j * np.radians(rng.normal(0, sigma_cal, (n_freq, trials, n)))))
        noise = ((rng.normal(size=(n_freq, trials, n)) + 1j * rng.normal(size=(n_freq, trials, n)))
                 / np.sqrt(2) * 10 ** (-snr_eff / 20))
        fi, ti = np.meshgrid(np.arange(n_freq), np.arange(trials), indexing="ij")
        out = np.empty(len(X))
        for m, x in enumerate(X):
            p = x.reshape(n, 3)
            y = np.exp(1j * k[:, None, None] * (d @ p.T)) * np.cos(el)[..., None] * cal + noise
            A = np.exp(-1j * k[:, None, None] * (Dg @ p.T)[None])         # (F, G, N)
            pw = (np.abs(A @ y.transpose(0, 2, 1)) ** 2).transpose(0, 2, 1).reshape(n_freq, trials, na, ne)
            j = np.argmax(pw.reshape(n_freq, trials, -1), axis=-1)
            ia, ie = j // ne, j % ne
            oa = parab(pw[fi, ti, (ia - 1) % na, ie], pw[fi, ti, ia, ie], pw[fi, ti, (ia + 1) % na, ie])
            inner = (ie > 0) & (ie < ne - 1)
            oe = np.where(inner, parab(pw[fi, ti, ia, np.maximum(ie - 1, 0)], pw[fi, ti, ia, ie],
                                       pw[fi, ti, ia, np.minimum(ie + 1, ne - 1)]), 0.0)
            dh = direction(azg[ia] + oa * st, elg[ie] + oe * st)
            err = np.degrees(np.arccos(np.clip(np.sum(dh * d, axis=-1), -1, 1)))
            loss = np.sqrt(np.mean(np.minimum(err, 30.0) ** 2, axis=1))
            out[m] = loss.max() + 0.25 * loss.mean() + 100.0 * violation3(p)
        return out

    return score


def search(log, quick=False, n=16):
    """evolve, climb (from the rings with the inner ring raised 300 mm) and random search, same budget, in the
    500 x 500 x 300 mm volume; and evolve restricted to the plate (z = 0) with the same score."""
    score = make_score3(n, n_freq=4 if quick else 8, trials=16 if quick else 64)
    pop, gens = (8, 5) if quick else (40, 200)
    draws = 3 if quick else 10
    lo = np.tile([-R_OUT, -R_OUT, 0.0], n)
    hi = np.tile([R_OUT, R_OUT, HEIGHT], n)
    hi_flat = np.tile([R_OUT, R_OUT, 0.0], n)
    flat_repair = lambda x: repair3(x, flat=True)  # noqa: E731
    out, t0 = {}, time.time()
    ga = evolve(score, lo, hi, pop=pop, generations=gens, elite=2, sigma=0.15, crossover=0.5, seed=16, repair=repair3)
    start = repair3(hand_designs()["inner ring +300 mm"].reshape(-1))
    lam = 8
    hc = climb(score, start, lo, hi, steps=ga.evaluations // (lam + 1), lam=lam, sigma=0.05, seed=16, repair=repair3)
    rs = random_search(score, lo, hi, batches=ga.evaluations // pop, batch=pop, seed=16, repair=repair3)
    gp = evolve(score, lo, hi_flat, pop=pop, generations=gens, elite=2, sigma=0.15, crossover=0.5, seed=16,
                repair=flat_repair)
    for name, res in (("evolved 3-D", ga), ("climbed 3-D", hc), ("random 3-D", rs), ("evolved planar (joint score)", gp)):
        X, F = reselect(score, res.finalists[:10], draws=draws)
        out[name] = dict(x=X[0].reshape(n, 3), score=float(F[0]), evaluations=res.evaluations + len(X) * draws,
                         history=[dict(evals=h["evals"], best=h["best"]) for h in res.history])
        log(f"  search: {name:30s} {res.evaluations} evaluations + {len(X)} finalists re-scored: mean score {F[0]:.2f}")
    log(f"  searches took {time.time() - t0:.0f} s")
    hand = hand_designs()
    Xh = np.array([v.reshape(-1) for v in hand.values()])
    Fh = np.mean([score(Xh, 2_000_000 + g) for g in range(draws)], axis=0)
    scores = {k: float(v) for k, v in zip(hand, Fh)}
    for k, v in scores.items():
        log(f"  search score of {k:40s} {v:.2f}")
    return out, scores


# ---- cross-check against evolve_layout.check ---------------------------------------------------------------
def cross_check(log, quick=False):
    """This file's covariance model and azimuth-only estimator, run as evolve_layout.check runs (isotropic
    whips, sources 0-30 deg up), against check() itself and against the stored results.json."""
    freqs = np.array([0.4e9, 1.0e9, 1.6e9, 2.0e9, 2.3e9])
    trials = 300 if quick else 1500
    stored = json.loads((HERE / "evolve" / "results.json").read_text())["designs"]
    rows = {}
    for name, p2 in (("ring8_225mm", planar.ring(8, 0.225)),
                     ("rings16_225_112mm", np.vstack([planar.ring(8, 0.225), planar.ring(8, 0.112, np.pi / 8)]))):
        p3 = np.column_stack([p2, np.zeros(len(p2))])
        e = monte_carlo(p3, freqs, trials=trials, seed=99, estimators=("az",), bins=[(0.0, 30.0)], pattern=False)["az_only"][:, 0]
        mine = np.column_stack([np.median(e, 1), np.percentile(e, 95, 1), np.mean(e > 10, 1) * 100])
        ref = planar.check(p2, freqs, el_max=30.0, trials=trials, seed=4242)
        sf = np.array(stored[name]["freqs_hz"])
        ii = [int(np.argmin(np.abs(sf - f))) for f in freqs]
        st = np.column_stack([np.array(stored[name]["elevated"][k])[ii] for k in ("median", "p95", "gross_pct")])
        rows[name] = dict(freqs_hz=freqs.tolist(), this_study=mine.round(3).tolist(), evolve_layout_check=ref.round(3).tolist(),
                          stored_results_json=st.round(3).tolist())
        log(f"  cross-check {name} (isotropic, sources 0-30 deg up, azimuth-only): median / p95 / gross %")
        for f, a, b, c in zip(freqs, mine, ref, st):
            log(f"    {f / 1e9:.1f} GHz  this study {a[0]:.2f} / {a[1]:6.2f} / {a[2]:4.1f}   check() {b[0]:.2f} / {b[1]:6.2f} / {b[2]:4.1f}"
                f"   results.json {c[0]:.2f} / {c[1]:6.2f} / {c[2]:4.1f}")
    return rows


# ---- part 2: the panel ---------------------------------------------------------------------------------------
PITCH, FACE_W, PAT_EXP = 0.025, 0.170, 1.3          # 150 mm of cards plus frame = 170 mm face (estimate)
PAT_CUT = 80.0                                      # deg off boresight beyond which an element is not used or trusted
G0_DBI = {2.44e9: 6.4, 3.5e9: 10.7, 5.8e9: 15.5}   # REDESIGN.md section 4, realised boresight gain


def prism(n_faces):
    """Element positions (16 n, 3), outward normals, and tile index for 4x4 tiles on a vertical prism."""
    apothem = 0.0 if n_faces == 1 else FACE_W / (2 * np.tan(np.pi / n_faces))
    c = (np.arange(4) - 1.5) * PITCH
    Y, Z = np.meshgrid(c, c, indexing="ij")
    pos, nrm, tid = [], [], []
    for t in range(n_faces):
        ph = 2 * np.pi * t / n_faces
        nv, tv = np.array([np.cos(ph), np.sin(ph), 0.0]), np.array([-np.sin(ph), np.cos(ph), 0.0])
        pos.append(apothem * nv + Y.reshape(-1, 1) * tv + Z.reshape(-1, 1) * np.array([0, 0, 1.0]))
        nrm.append(np.tile(nv, (16, 1)))
        tid.append(np.full(16, t))
    return np.vstack(pos), np.vstack(nrm), np.concatenate(tid)


def panel_manifold(pos, nrm, k, az):
    d = direction(az, np.zeros_like(az))
    c = d @ nrm.T
    c = np.where(c > np.cos(np.radians(PAT_CUT)), c, 0.0)
    return c ** (PAT_EXP / 2) * np.exp(1j * k * (d @ pos.T))


def panel_df(pos, nrm, k, az_true, window, rng, trials):
    """Azimuth Bartlett with the pattern-weighted manifold; `window` = (centre, half width) of the search in rad."""
    az = az_true + np.radians(rng.uniform(-2.5, 2.5, trials))
    a = panel_manifold(pos, nrm, k, az)
    R = covariances(pos, k, az, np.zeros(trials), rng, amp=np.abs(a))
    def scan(grid):                                   # grid (T, g) or (g,)
        A = panel_manifold(pos, nrm, k, grid)
        A = A / np.maximum(np.linalg.norm(A, axis=-1, keepdims=True), 1e-9)
        if A.ndim == 2:
            return np.real(np.einsum('gi,tij,gj->tg', A.conj(), R, A))
        return np.real(np.einsum('tgi,tij,tgj->tg', A.conj(), R, A))
    coarse = window[0] + np.radians(np.arange(-window[1], window[1] + 1e-9, 0.5))
    c0 = coarse[np.argmax(scan(coarse), axis=1)]
    fine = c0[:, None] + np.radians(np.arange(-0.5, 0.5001, 0.02))[None]
    est = fine[np.arange(trials), np.argmax(scan(fine), axis=1)]
    e = wrap_deg(est - az)
    return float(np.median(e)), float(np.percentile(e, 95)), float(np.mean(e > 10) * 100)


def panel_study(log, quick=False):
    trials = 100 if quick else 500
    az_deg = np.arange(0, 181, 5.0)
    rng = np.random.default_rng(5)
    res = dict(az_deg=az_deg.tolist(), freqs_hz=list(G0_DBI), options={}, face_width_m=FACE_W)
    single_pos, single_nrm, _ = prism(1)
    for f, g0 in G0_DBI.items():
        k = 2 * np.pi * f / C0
        # one tile against the angle off its boresight
        off = np.arange(0, PAT_CUT - 4.9, 5.0)
        tile = np.array([panel_df(single_pos, single_nrm, k, np.radians(o), (0.0, PAT_CUT), rng, trials) for o in off])
        def tile_at(o):                               # o = angle off boresight, deg
            i = int(round(o / 5.0))
            return tile[i] if o <= PAT_CUT - 5 else np.array([np.nan, np.nan, 100.0])
        def gain(o):
            return g0 + 10 * PAT_EXP * np.log10(np.cos(np.radians(o))) if o < PAT_CUT else -np.inf
        opts = {}
        opts["one fixed tile"] = dict(elements_built=16, elements_per_direction=16,
                                      gain_dbi=[gain(a) for a in az_deg], df=[tile_at(a).tolist() for a in az_deg])
        opts["one tile on a rotator"] = dict(elements_built=16, elements_per_direction=16,
                                             gain_dbi=[g0 for _ in az_deg], df=[tile[0].tolist() for _ in az_deg])
        for nf in (3, 4):
            sector = 360.0 / nf
            fold = np.abs((az_deg + sector / 2) % sector - sector / 2)
            opts[f"{nf}-face prism, best tile"] = dict(
                elements_built=16 * nf, elements_per_direction=16,
                gain_dbi=[gain(o) for o in fold], df=[tile_at(o).tolist() for o in fold])
            pos, nrm, _ = prism(nf)
            uniq = np.arange(0, sector / 2 + 0.1, 5.0)
            joint = {float(u): panel_df(pos, nrm, k, np.radians(u), (np.radians(u), 180.0), rng, trials) for u in uniq}
            gj = []
            for a in az_deg:
                c = np.cos(np.radians(a) - 2 * np.pi * np.arange(nf) / nf)
                c = np.where(c > np.cos(np.radians(PAT_CUT)), c, 0.0)
                gj.append(g0 + 10 * np.log10(np.sum(c ** PAT_EXP)))
            lit = [int(16 * np.sum(np.cos(np.radians(a) - 2 * np.pi * np.arange(nf) / nf) > np.cos(np.radians(PAT_CUT)))) for a in az_deg]
            opts[f"{nf}-face prism, all visible tiles"] = dict(
                elements_built=16 * nf, elements_per_direction=lit, gain_dbi=gj,
                df=[list(joint[float(5.0 * round(o / 5.0))]) for o in fold])
        res["options"][str(f)] = opts
        log(f"  panel {f / 1e9:.2f} GHz (boresight {g0} dBi):")
        for name, o in opts.items():
            g = np.array(o["gain_dbi"], float)
            d = np.array(o["df"], float)
            cov = np.isfinite(g) & (g >= g0 - 3.0)
            ok = d[:, 2] < 1.0
            o["coverage_within_3dB_pct"] = float(np.mean(cov[:-1]) * 100)     # 0..175 deg mirrors to 360
            o["coverage_df_pct"] = float(np.mean(ok[:-1]) * 100)
            o["min_gain_dbi"] = float(g.min())
            o["df_p95_worst_covered"] = float(np.nanmax(d[ok, 1])) if ok.any() else None
            o["df_median_boresight"] = float(d[0, 0])
            log(f"    {name:34s} min gain {g.min():6.1f} dBi, azimuth within 3 dB of boresight {o['coverage_within_3dB_pct']:5.1f} %,"
                f" DF usable (<1 % gross) {o['coverage_df_pct']:5.1f} %, median / p95 at boresight {d[0, 0]:.2f} / {d[0, 1]:.2f}, worst covered {o['df_p95_worst_covered']:.2f} deg")
    return res


# ---- figures --------------------------------------------------------------------------------------------------
def plot(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    D = res["designs"]
    freqs = np.array(res["freqs_hz"]) / 1e9
    names = list(D)
    # layouts: top and side view
    cols = min(6, len(names))
    rows = -(-len(names) // cols)
    fig, axes = plt.subplots(2 * rows, cols, figsize=(2.7 * cols, 4.6 * rows), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    for i, k in enumerate(names):
        p = np.array(D[k]["positions_m"]) * 1000
        top, side = axes[2 * (i // cols), i % cols], axes[2 * (i // cols) + 1, i % cols]
        top.axis("on"); side.axis("on")
        top.add_patch(plt.Circle((0, 0), 250, fill=False, color="0.6"))
        sc = top.scatter(p[:, 0], p[:, 1], c=p[:, 2], cmap="viridis", vmin=0, vmax=300, s=28, edgecolor="k", lw=0.4)
        top.set_aspect("equal"); top.set_xlim(-265, 265); top.set_ylim(-265, 265)
        top.set_title(k, fontsize=7.5); top.tick_params(labelsize=6)
        side.scatter(p[:, 0], p[:, 2], c=p[:, 2], cmap="viridis", vmin=0, vmax=300, s=28, edgecolor="k", lw=0.4)
        side.axhline(0, color="0.6", lw=0.8)
        side.set_xlim(-265, 265); side.set_ylim(-30, 330); side.set_aspect("equal"); side.tick_params(labelsize=6)
        side.set_xlabel("x, mm", fontsize=6.5); side.set_ylabel("height, mm", fontsize=6.5)
    fig.suptitle("16 whips: top view (colour = height, 0-300 mm) and side view", fontsize=9)
    fig.tight_layout(); fig.savefig(OUT / "layouts.png", dpi=130); plt.close(fig)

    show = [k for k in ("planar rings 225/112", "inner ring +300 mm", "stacked rings 225/225, 300 mm apart",
                        "cone 8/5/3", "evolved 3-D", "evolved planar (joint score)") if k in D]
    binlab = [f"sources {int(a)}-{int(b)} deg up" for a, b in BINS]
    # azimuth p95 and elevation p95 against frequency, per elevation bin
    fig, axes = plt.subplots(2, 3, figsize=(15, 7.6), sharex=True)
    for b in range(3):
        ax = axes[0, b]
        pf = D["planar rings 225/112"]["per_freq"]
        ax.semilogy(freqs, np.array(pf["az_only"]["p95"])[:, b], "k--", lw=1.6, label="planar rings, azimuth-only (baseline)")
        for j, k in enumerate(show):
            ax.semilogy(freqs, np.array(D[k]["per_freq"]["joint_az"]["p95"])[:, b], lw=1.5, color=f"C{j}", label=f"{k}, joint")
        ax.set_title(f"azimuth error, 95th percentile; {binlab[b]}", fontsize=9); ax.grid(alpha=0.3, which="both")
        ax = axes[1, b]
        for j, k in enumerate(show):
            ax.semilogy(freqs, np.array(D[k]["per_freq"]["joint_el"]["p95"])[:, b], lw=1.5, color=f"C{j}", label=f"{k}, joint")
        ax.set_title(f"elevation error, 95th percentile; {binlab[b]}", fontsize=9); ax.grid(alpha=0.3, which="both")
        ax.set_xlabel("GHz")
    axes[0, 0].set_ylabel("deg"); axes[1, 0].set_ylabel("deg")
    axes[0, 0].legend(fontsize=6.5); axes[1, 0].legend(fontsize=6.5)
    fig.suptitle("16 vertical whips, 3 deg / 0.5 dB calibration residual, 20 dB SNR (horizon), 64 samples", fontsize=10)
    fig.tight_layout(); fig.savefig(OUT / "errors_vs_freq.png", dpi=130); plt.close(fig)

    # whole-band summary bars
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    y = np.arange(len(names))
    for ax, key, title in ((axes[0], "joint_az", "azimuth error, 95th percentile (joint estimator), deg"),
                           (axes[1], "joint_el", "elevation error, 95th percentile (joint estimator), deg"),
                           (axes[2], "az_only", "gross azimuth errors, azimuth-only estimator, % of trials")):
        for b in range(3):
            fld = "gross_pct" if key == "az_only" else "p95"
            v = [D[k]["bands"]["0.4-2.3 GHz"][key][b][fld] for k in names]
            ax.barh(y + (b - 1) * 0.27, v, 0.25, label=binlab[b], color=["#3b6fb6", "#e0a030", "#b6403b"][b])
        ax.set_yticks(y); ax.set_yticklabels(names if ax is axes[0] else [], fontsize=7.5); ax.invert_yaxis()
        ax.set_title(title, fontsize=9); ax.grid(alpha=0.3, axis="x")
        if key != "az_only":
            ax.set_xscale("log")
    axes[0].legend(fontsize=7.5)
    fig.suptitle("All of 0.4-2.3 GHz pooled", fontsize=10)
    fig.tight_layout(); fig.savefig(OUT / "summary.png", dpi=130); plt.close(fig)

    if res.get("search"):
        fig, ax = plt.subplots(figsize=(6, 3.8))
        for k, v in res["search"]["found"].items():
            ax.plot([h["evals"] for h in v["history"]], [h["best"] for h in v["history"]], lw=1, label=f"{k} (re-scored {v['score']:.2f})")
        ax.set_yscale("log"); ax.set_xlabel("evaluations"); ax.set_ylabel("best score of the generation (noisy)")
        ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=7)
        fig.tight_layout(); fig.savefig(OUT / "convergence.png", dpi=130); plt.close(fig)

    P = res.get("panel")
    if P:
        az = np.array(P["az_deg"])
        fig, axes = plt.subplots(2, 3, figsize=(15, 7.4), sharex=True)
        styles = {"one fixed tile": dict(color="k", lw=2), "one tile on a rotator": dict(color="C2", lw=2, ls="--"),
                  "3-face prism, best tile": dict(color="C0", lw=1.6), "3-face prism, all visible tiles": dict(color="C0", lw=1, ls=":"),
                  "4-face prism, best tile": dict(color="C3", lw=1.6), "4-face prism, all visible tiles": dict(color="C3", lw=1, ls=":")}
        for j, f in enumerate(P["freqs_hz"]):
            for name, o in P["options"][str(f)].items():
                g = np.array(o["gain_dbi"], float); g[~np.isfinite(g)] = np.nan
                d = np.array(o["df"], float)
                axes[0, j].plot(az, g, label=name, **styles[name])
                axes[1, j].semilogy(az, np.where(d[:, 2] < 50, d[:, 1], np.nan), label=name, **styles[name])
            axes[0, j].set_title(f"{float(f) / 1e9:.2f} GHz: gain toward the source, dBi", fontsize=9)
            axes[1, j].set_title("azimuth DF error, 95th percentile, deg", fontsize=9)
            axes[1, j].set_xlabel("source azimuth from the first tile's boresight, deg (mirror for 180-360)")
            for ax in axes[:, j]:
                ax.grid(alpha=0.3, which="both")
        axes[0, 0].legend(fontsize=7)
        fig.tight_layout(); fig.savefig(OUT / "panel.png", dpi=130); plt.close(fig)


# ---- main -----------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--plot-only", action="store_true")
    ap.add_argument("--skip-search", action="store_true")
    ap.add_argument("--panel-only", action="store_true", help="keep part 1 from results.json, redo part 2 and the figures")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    name = "results_quick.json" if args.quick else "results.json"
    if args.plot_only:
        plot(json.loads((OUT / name).read_text()))
        return
    t00 = time.time()

    def log(s):
        print(s, flush=True)

    if args.panel_only:
        res = json.loads((OUT / name).read_text())
        log("panel band: one tile, prisms, rotator")
        res["panel"] = panel_study(log, args.quick)
        (OUT / name).write_text(json.dumps(res, indent=1))
        plot(res)
        log(f"done in {time.time() - t00:.0f} s")
        return

    quick = args.quick
    freqs = np.arange(0.4e9, 2.3001e9, 0.4e9 if quick else 0.1e9)
    trials = 60 if quick else 600
    res = dict(freqs_hz=freqs.tolist(), elevation_bins_deg=BINS, trials_per_bin_per_freq=trials,
               stress=dict(sigma_cal_deg=5.0, snr_db_horizon=0.0, trials_per_bin_per_freq=max(trials // 2, 30)),
               conventions=dict(sigma_cal_deg=3.0, amp_db=0.5, snr_db_horizon=20.0, samples=64, gross_deg=10.0,
                                pattern="cos(elevation) amplitude", volume_m=[2 * R_OUT, 2 * R_OUT, HEIGHT], min_spacing_m=S_MIN),
               designs={})
    log("cross-check against experiments/evolve_layout.py")
    res["cross_check"] = cross_check(log, quick)
    designs = hand_designs()
    stored = json.loads((HERE / "evolve" / "results.json").read_text())["designs"]
    if "evolved16" in stored:
        p2 = np.array(stored["evolved16"]["positions_m"])
        designs["evolved16 planar (evolve_layout)"] = np.column_stack([p2, np.zeros(len(p2))])
    if not args.skip_search:
        log("search in the 500 x 500 x 300 mm volume (fast single-snapshot model, joint azimuth/elevation)")
        found, hand_scores = search(log, quick)
        res["search"] = dict(found={k: dict(score=v["score"], evaluations=v["evaluations"], history=v["history"]) for k, v in found.items()},
                             hand_scores=hand_scores)
        for k, v in found.items():
            assert violation3(v["x"]) < 1e-6, (k, violation3(v["x"]))
            designs[k] = v["x"]
    log(f"final numbers: covariance model, {len(freqs)} frequencies x {len(BINS)} elevation bins x {trials} trials")
    for i, (k, p) in enumerate(designs.items()):
        t0 = time.time()
        err = monte_carlo(p, freqs, trials=trials, seed=1000 + i)
        s = stats(err, freqs)
        s["stress"] = stats(monte_carlo(p, freqs, trials=max(trials // 2, 30), seed=5000 + i, sigma_cal=5.0, snr_db=0.0), freqs)["bands"]
        s["positions_m"] = np.round(p, 4).tolist()
        res["designs"][k] = s
        log(f"  {k}  ({time.time() - t0:.0f} s)")
        for band in s["bands"]:
            for b, (lo, hi) in enumerate(BINS):
                r = s["bands"][band]
                a, j, e = r["az_only"][b], r["joint_az"][b], r["joint_el"][b]
                log(f"    {band:12s} el {int(lo):2d}-{int(hi):2d}: az-only med {a['median']:6.2f} p95 {a['p95']:6.2f} gross {a['gross_pct']:5.1f} % (worst f {a['gross_pct_worst_freq']:5.1f})"
                    f" | joint az med {j['median']:5.2f} p95 {j['p95']:6.2f} gross {j['gross_pct']:4.1f} % (worst f {j['gross_pct_worst_freq']:4.1f})"
                    f" | joint el med {e['median']:5.2f} p95 {e['p95']:6.2f}")
        for b, (lo, hi) in enumerate(BINS):
            r = s["stress"]["0.4-2.3 GHz"]
            a, j, e = r["az_only"][b], r["joint_az"][b], r["joint_el"][b]
            log(f"    stress 5 deg/0 dB el {int(lo):2d}-{int(hi):2d}: az-only p95 {a['p95']:6.2f} gross {a['gross_pct']:5.1f} % | joint az med {j['median']:5.2f} p95 {j['p95']:6.2f}"
                f" gross {j['gross_pct']:4.1f} % (worst f {j['gross_pct_worst_freq']:4.1f}) | joint el med {e['median']:5.2f} p95 {e['p95']:6.2f}")
        (OUT / name).write_text(json.dumps(res, indent=1))
    log("panel band: one tile, prisms, rotator")
    res["panel"] = panel_study(log, quick)
    res["runtime_s"] = time.time() - t00
    (OUT / name).write_text(json.dumps(res, indent=1))
    plot(res)
    log(f"done in {time.time() - t00:.0f} s")


if __name__ == "__main__":
    main()
