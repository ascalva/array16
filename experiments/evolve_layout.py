"""One fixed layout of external whips for the whole 0.4-2.3 GHz band, found by evolutionary search.

Rev B (REDESIGN.md, section 4) mounts the external whips as rings sized for one band at a time
(8 whips: r = 225 / 112 / 78 mm for 0.4-0.8 / 0.8-1.6 / 1.15-2.3 GHz), so the ring is rebuilt for each
band. Question: can one fixed set of holes on the same 500 mm plate serve the whole band?

Plate: whips within 240 mm of the centre, none within 60 mm of it (the centre monopole is the in-place
calibration source), none closer than 75 mm to each other (0.1 wavelength at 0.4 GHz, to keep coupling
modest).

Search model (fast, inside the search): azimuth-only direction finding of one source 0-30 deg above the
horizon (the estimator assumes 0 deg; for big apertures this, not ambiguity, is what fails), isotropic whips; each whip keeps a calibration residual (phase sigma_cal, 0.5 dB amplitude); one coherent
snapshot at the per-whip SNR plus 10log10(K); Bartlett beamscan on a 1 deg grid with parabolic refinement.
Per frequency the capped RMS error sqrt(mean(min(|e|, 30 deg)^2)); a design's score is the worst
frequency plus a quarter of the mean. Frequencies (one per log-spaced stratum), directions and errors are
redrawn every generation, the same for every design of that generation.

Check model (separate code, run only on finished designs, fresh seeds): the design panel's ring Monte
Carlo (lead2/calc/ring_mc.py) for any layout, sources 0-15 deg (nominal) or 0-30 deg (elevated) up: K = 64 noisy samples, sample covariance, Bartlett beamscan
on 1 deg then 0.05 deg, gross error = more than 10 deg; median, 95th percentile and gross-error rate on a
50 MHz grid.

    python experiments/evolve_layout.py                  # about 10 minutes
    python experiments/evolve_layout.py --quick          # a smoke run
    python experiments/evolve_layout.py --replicates 5   # each method 5 times with new seeds (about 40 minutes)
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from array16.evolve import climb, evolve, random_search, reselect  # noqa: E402

C0 = 299_792_458.0
R_OUT, R_IN, S_MIN = 0.240, 0.060, 0.075
EL_MAX = 30.0      # search: sources 0-30 deg above the horizon (a drone, a rooftop); the estimator assumes 0
BAND = (0.4e9, 2.3e9)
OUT = Path(__file__).resolve().parent / "evolve"


# ---- the plate ------------------------------------------------------------------------------------------
def repair(x: np.ndarray) -> np.ndarray:
    """Move whips onto the plate's annulus and push apart any pair closer than S_MIN; order them by angle
    (so that crossover pairs up whips that sit in similar places)."""
    p = x.reshape(-1, 2).copy()
    for _ in range(60):
        r = np.hypot(p[:, 0], p[:, 1])
        p *= (np.clip(r, R_IN, R_OUT) / np.maximum(r, 1e-12))[:, None]
        d = p[:, None, :] - p[None, :, :]
        dist = np.linalg.norm(d, axis=2)
        np.fill_diagonal(dist, np.inf)
        i, j = np.nonzero(np.triu(dist < S_MIN - 1e-9, 1))
        if not len(i):
            break
        for a, b in zip(i, j):
            v = p[a] - p[b]
            n = np.linalg.norm(v)
            v = v / n if n > 1e-12 else np.array([1.0, 0.0])
            push = (S_MIN - n) / 2 + 1e-6
            p[a] += push * v
            p[b] -= push * v
    p = p[np.argsort(np.arctan2(p[:, 1], p[:, 0]))]
    return p.reshape(-1)


def violation(p: np.ndarray) -> float:
    r = np.hypot(p[:, 0], p[:, 1])
    d = np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)
    np.fill_diagonal(d, np.inf)
    return float(np.sum(np.maximum(r - R_OUT, 0)) + np.sum(np.maximum(R_IN - r, 0))
                 + np.sum(np.triu(np.maximum(S_MIN - d, 0), 1)))


def ring(n: int, radius: float, rot: float = 0.0) -> np.ndarray:
    b = 2 * np.pi * np.arange(n) / n + rot
    return np.column_stack([radius * np.cos(b), radius * np.sin(b)])


# ---- search model -----------------------------------------------------------------------------------------
GRID = np.radians(np.arange(360.0))


def _refine(pw: np.ndarray, i: np.ndarray) -> np.ndarray:
    """Parabolic peak refinement on a circular 1 deg grid; pw (..., G), i (...) -> angle in rad."""
    g = pw.shape[-1]
    y0 = np.take_along_axis(pw, ((i - 1) % g)[..., None], -1)[..., 0]
    y1 = np.take_along_axis(pw, i[..., None], -1)[..., 0]
    y2 = np.take_along_axis(pw, ((i + 1) % g)[..., None], -1)[..., 0]
    den = y0 - 2 * y1 + y2
    off = np.where(np.abs(den) > 1e-30, 0.5 * (y0 - y2) / np.where(np.abs(den) > 1e-30, den, 1.0), 0.0)
    return np.radians(i + np.clip(off, -0.5, 0.5))


def make_score(n: int, *, n_freq=16, trials=96, sigma_cal=4.0, snr_db=10.0, k_samples=64, el_max=EL_MAX, seed=0):
    snr_eff = snr_db + 10 * np.log10(k_samples)
    edges = np.geomspace(*BAND, n_freq + 1)

    def score(X: np.ndarray, gen: int) -> np.ndarray:
        rng = np.random.default_rng([seed, gen])
        f = edges[:-1] * (edges[1:] / edges[:-1]) ** rng.random(n_freq)
        k = 2 * np.pi * f / C0
        az = rng.uniform(0, 2 * np.pi, (n_freq, trials))
        ce = np.cos(np.radians(rng.uniform(0, el_max, (n_freq, trials))))
        cal = (10 ** (rng.normal(0, 0.5, (n_freq, trials, n)) / 20)
               * np.exp(1j * np.radians(rng.normal(0, sigma_cal, (n_freq, trials, n)))))
        noise = ((rng.normal(size=(n_freq, trials, n)) + 1j * rng.normal(size=(n_freq, trials, n)))
                 / np.sqrt(2) * 10 ** (-snr_eff / 20))
        out = np.empty(len(X))
        for m, x in enumerate(X):
            p = x.reshape(n, 2)
            proj_t = (p[:, 0] * np.cos(az)[..., None] + p[:, 1] * np.sin(az)[..., None]) * ce[..., None]   # (F, T, N)
            y = np.exp(1j * k[:, None, None] * proj_t) * cal + noise
            proj_g = p[:, 0] * np.cos(GRID)[:, None] + p[:, 1] * np.sin(GRID)[:, None]        # (G, N)
            A = np.exp(-1j * k[:, None, None] * proj_g[None])                                 # (F, G, N)
            pw = np.abs(A @ y.transpose(0, 2, 1)) ** 2                                         # (F, G, T)
            pw = pw.transpose(0, 2, 1)                                                         # (F, T, G)
            est = _refine(pw, np.argmax(pw, axis=-1))
            err = np.degrees(np.angle(np.exp(1j * (est - az))))
            loss = np.sqrt(np.mean(np.minimum(np.abs(err), 30.0) ** 2, axis=1))               # (F,)
            out[m] = loss.max() + 0.25 * loss.mean() + 100.0 * violation(p)
        return out

    return score


# ---- check model ------------------------------------------------------------------------------------------
def check(p: np.ndarray, freqs, *, sigma_cal=3.0, snr_db=20.0, el_max=15.0, trials=1500, k_samples=64, seed=12345):
    """Independent Monte Carlo (the design panel's ring_mc.py for any layout). Returns per-frequency
    median, 95th percentile (deg) and gross-error rate (%, errors over 10 deg)."""
    rng = np.random.default_rng(seed)
    n = len(p)
    coarse = np.radians(np.arange(0, 360, 1.0))
    res = []
    for f in freqs:
        k = 2 * np.pi * f / C0
        th = rng.uniform(0, 2 * np.pi, trials)
        ce = np.cos(np.radians(rng.uniform(0, el_max, trials)))[:, None]
        a = np.exp(1j * k * ce * (np.cos(th)[:, None] * p[:, 0] + np.sin(th)[:, None] * p[:, 1]))  # (T, N)
        g = 10 ** (rng.normal(0, 0.5, (trials, n)) / 20) * np.exp(1j * np.radians(rng.normal(0, sigma_cal, (trials, n))))
        s = (rng.normal(size=(trials, k_samples)) + 1j * rng.normal(size=(trials, k_samples))) / np.sqrt(2)
        w = ((rng.normal(size=(trials, n, k_samples)) + 1j * rng.normal(size=(trials, n, k_samples)))
             / np.sqrt(2) * 10 ** (-snr_db / 20))
        x = (a * g)[:, :, None] * s[:, None, :] + w                                            # (T, N, K)
        R = x @ x.conj().transpose(0, 2, 1) / k_samples                                         # (T, N, N)
        Ac = np.exp(1j * k * (np.cos(coarse)[:, None] * p[:, 0] + np.sin(coarse)[:, None] * p[:, 1]))
        pc = np.real(np.einsum('gi,tij,gj->tg', Ac.conj(), R, Ac))
        c0 = coarse[np.argmax(pc, axis=1)]
        fine = c0[:, None] + np.radians(np.arange(-1.5, 1.5001, 0.05))[None, :]                 # (T, 61)
        Af = np.exp(1j * k * (np.cos(fine)[..., None] * p[:, 0] + np.sin(fine)[..., None] * p[:, 1]))
        pf = np.real(np.einsum('tgi,tij,tgj->tg', Af.conj(), R, Af))
        est = fine[np.arange(trials), np.argmax(pf, axis=1)]
        e = np.abs(np.degrees(np.angle(np.exp(1j * (est - th)))))
        res.append((float(np.median(e)), float(np.percentile(e, 95)), float(np.mean(e > 10) * 100)))
    return np.array(res)


# ---- the study --------------------------------------------------------------------------------------------
def run(n: int, budget_gens: int, pop: int, quick: bool, log, seed_offset: int = 0):
    """GA, hill climb from rev B's middle ring, and random search, each with the same number of evaluations;
    then each method's finalists are re-scored on fresh draws and the best by mean is its answer."""
    score = make_score(n, n_freq=8 if quick else 16, trials=32 if quick else 96)
    lo, hi = -R_OUT * np.ones(2 * n), R_OUT * np.ones(2 * n)
    out = {}
    t0 = time.time()
    sd = n + 1000 * seed_offset
    ga = evolve(score, lo, hi, pop=pop, generations=budget_gens, elite=2, sigma=0.15, crossover=0.5, seed=sd, repair=repair)
    start = repair(ring(n, 0.112 if n >= 8 else 0.074).reshape(-1))      # rev B's middle ring: "a design we like"
    lam = 8
    hc = climb(score, start, lo, hi, steps=ga.evaluations // (lam + 1), lam=lam, sigma=0.05, seed=sd, repair=repair)
    rs = random_search(score, lo, hi, batches=ga.evaluations // pop, batch=pop, seed=sd, repair=repair)
    for name, res in (("evolved", ga), ("climbed", hc), ("random", rs)):
        X, F = reselect(score, res.finalists[:10], draws=4 if quick else 12)
        out[name] = dict(x=X[0], reselect_score=float(F[0]), evaluations=res.evaluations + len(X) * (4 if quick else 12),
                         history=res.history)
        log(f"  N={n}: {name:8s} {res.evaluations} evaluations + {len(X)} finalists re-scored: mean score {F[0]:.2f}")
    log(f"  N={n}: searches took {time.time() - t0:.0f} s")
    return out


def replicates(k: int, log=print):
    """Each method k times with new seeds: is one method really better, or was one run lucky?"""
    freqs = np.arange(0.40e9, 2.3001e9, 0.10e9)
    conds = dict(nominal=dict(sigma_cal=3.0, snr_db=20.0, el_max=15.0, seed=12345),
                 elevated=dict(sigma_cal=3.0, snr_db=20.0, el_max=30.0, seed=4242),
                 stress=dict(sigma_cal=5.0, snr_db=0.0, el_max=30.0, seed=777))
    out = {}
    for n, (gens, pop) in {4: (120, 32), 8: (200, 40), 16: (200, 40)}.items():
        for rep in range(1, k + 1):
            found = run(n, gens, pop, False, lambda *_: None, seed_offset=rep)
            for name, v in found.items():
                p = v["x"].reshape(n, 2)
                row = dict(fresh=v["reselect_score"], positions_m=np.round(p, 4).tolist())
                for c, kw in conds.items():
                    r = check(p, freqs, trials=800, **kw)
                    row[c] = dict(p95_max=float(r[:, 1].max()), gross_max=float(r[:, 2].max()), median=float(np.median(r[:, 0])))
                out.setdefault(str(n), {}).setdefault(name, []).append(row)
                log(f"N={n} rep {rep} {name:8s} fresh {row['fresh']:.2f} | " + " | ".join(
                    f"{c}: p95 {row[c]['p95_max']:.1f}, gross {row[c]['gross_max']:.1f} %" for c in conds))
    (OUT / "replicates.json").write_text(json.dumps(out, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--replicates", type=int, default=0)
    ap.add_argument("--plot-only", action="store_true", help="redraw the figures from results.json")
    args = ap.parse_args()
    if args.plot_only:
        res = json.loads((OUT / "results.json").read_text())
        hs = {int(n): v for n, v in res.get("histories", {}).items()}
        plot(res, hs, False)
        return
    if args.replicates:
        OUT.mkdir(exist_ok=True)
        replicates(args.replicates)
        return
    quick = args.quick
    OUT.mkdir(exist_ok=True)
    log = print
    plan = {4: (120, 32), 8: (200, 40), 16: (200, 40)} if not quick else {4: (6, 8), 8: (6, 8)}
    freqs = np.arange(0.40e9, 2.3001e9, 0.10e9 if quick else 0.05e9)
    rings = {4: [0.133, 0.074, 0.042], 8: [0.225, 0.112, 0.078], 16: [0.225]}
    bands = [(0.40e9, 0.80e9), (0.80e9, 1.60e9), (1.15e9, 2.30e9)]
    results = dict(plate=dict(r_out_m=R_OUT, r_in_m=R_IN, min_spacing_m=S_MIN), band_hz=BAND, designs={})
    trials = 300 if quick else 1500

    conditions = dict(nominal=dict(sigma_cal=3.0, snr_db=20.0, el_max=15.0, seed=12345),
                      elevated=dict(sigma_cal=3.0, snr_db=20.0, el_max=30.0, seed=4242),
                      stress=dict(sigma_cal=5.0, snr_db=0.0, el_max=30.0, seed=777))
    results["conditions"] = conditions

    def summarise(name, p, how):
        d = dict(how=how, positions_m=np.round(p, 4).tolist(), freqs_hz=freqs.tolist())
        line = f"  {name:22s}"
        for cname, c in conditions.items():
            r = check(p, freqs, trials=trials, **c)
            d[cname] = dict(median=r[:, 0].tolist(), p95=r[:, 1].tolist(), gross_pct=r[:, 2].tolist())
            line += f" | {cname}: p95 max {r[:, 1].max():6.2f}, gross max {r[:, 2].max():4.1f} %, med {np.median(r[:, 0]):.2f}"
        results["designs"][name] = d
        log(line)

    histories = {}
    for n, (gens, pop) in plan.items():
        log(f"N = {n} whips")
        for r in rings[n]:
            summarise(f"ring{n}_{int(round(r * 1000))}mm", ring(n, r), f"ring, radius {r * 1000:.0f} mm")
        if n == 16:
            summarise("rings16_225_112mm", np.vstack([ring(8, 0.225), ring(8, 0.112, np.pi / 8)]), "two concentric rings of 8")
        if n == 8:
            summarise("rings8_225_112mm", np.vstack([ring(4, 0.225), ring(4, 0.112, np.pi / 4)]), "two concentric rings of 4")
        found = run(n, gens, pop, quick, log)
        histories[n] = {k: v["history"] for k, v in found.items()}
        how = dict(evolved="genetic algorithm", climbed="hill climb from rev B's middle ring", random="best random layouts")
        for name, v in found.items():
            p = v["x"].reshape(n, 2)
            assert violation(p) < 1e-6, (name, violation(p))
            summarise(f"{name}{n}", p, f"{how[name]}, {v['evaluations']} evaluations")
            results["designs"][f"{name}{n}"]["search_score_fresh"] = v["reselect_score"]

    # rebuild-per-band envelope for rev B's rings (what section 4 asks the owner to do)
    for n in (4, 8):
        env = []
        for i, f in enumerate(freqs):
            cands = [results["designs"][f"ring{n}_{int(round(r * 1000))}mm"]["nominal"]["p95"][i]
                     for r, (fl, fh) in zip(rings[n], bands) if fl <= f <= fh]
            env.append(min(cands) if cands else None)
        results[f"ring{n}_rebuilt_per_band_p95"] = env

    results["histories"] = {str(n): {k: [dict(evals=h["evals"], best=h["best"]) for h in v] for k, v in hs.items()}
                            for n, hs in histories.items()}
    (OUT / ("results_quick.json" if quick else "results.json")).write_text(json.dumps(results, indent=1))
    plot(results, histories, quick)


def plot(results, histories, quick):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    D = results["designs"]
    sfx = "_quick" if quick else ""
    # layouts of the designs that matter
    names = [k for k in ("ring8_112mm", "evolved8", "climbed8", "random8", "rings16_225_112mm", "evolved16", "random16", "evolved4")
             if k in D]
    cols = 4
    rows = -(-len(names) // cols)
    fig, axes = plt.subplots(rows, cols, figsize=(3.0 * cols, 3.2 * rows))
    for ax in np.ravel(axes)[len(names):]:
        ax.axis("off")
    for ax, k in zip(np.ravel(axes), names):
        p = np.array(D[k]["positions_m"]) * 1000
        ax.add_patch(plt.Circle((0, 0), 250, fill=False, color="0.6"))
        ax.add_patch(plt.Circle((0, 0), R_IN * 1000, fill=False, ls=":", color="0.6"))
        ax.plot(p[:, 0], p[:, 1], "o", color="C0" if "ring" in k else "C3")
        ax.plot(0, 0, "+", color="k")
        ax.set_aspect("equal"); ax.set_xlim(-260, 260); ax.set_ylim(-260, 260)
        ax.set_title(k, fontsize=9); ax.tick_params(labelsize=7)
    fig.suptitle("Whip positions on the 500 mm plate (mm); + = centre calibration monopole, dotted = 60 mm keep-out", fontsize=9)
    fig.tight_layout(); fig.savefig(OUT / f"layouts{sfx}.png", dpi=130); plt.close(fig)
    # accuracy (nominal p95) and robustness (gross errors with elevated sources) against frequency
    ns = [n for n in (8, 16) if f"evolved{n}" in D] or [n for n in (4, 8) if f"evolved{n}" in D]
    fig, axes = plt.subplots(2, len(ns), figsize=(5.6 * len(ns), 7.4), squeeze=False)
    for j, n in enumerate(ns):
        keys = [k for k in D if (k.startswith(("ring", "rings")) and k[len("rings" if k.startswith("rings") else "ring"):].startswith(f"{n}_"))
                or k in (f"evolved{n}", f"climbed{n}", f"random{n}")]
        f = np.array(D[f"evolved{n}"]["freqs_hz"]) / 1e9
        for k in keys:
            style = dict(lw=2.2, color="C3") if k.startswith("evolved") else dict(lw=1.1)
            axes[0, j].semilogy(f, D[k]["nominal"]["p95"], label=k, **style)
            axes[1, j].plot(f, D[k]["elevated"]["gross_pct"], label=k, **style)
        env = results.get(f"ring{n}_rebuilt_per_band_p95")
        if env:
            axes[0, j].semilogy(f, [e if e is not None else np.nan for e in env], "k--", lw=1.4, label="rev B ring rebuilt per band")
        axes[0, j].set_title(f"{n} whips: 95th percentile error, sources 0-15 deg up (3 deg cal, 20 dB)", fontsize=8.5)
        axes[1, j].set_title(f"{n} whips: gross errors (>10 deg), sources 0-30 deg up", fontsize=8.5)
        axes[0, j].set_ylabel("deg"); axes[1, j].set_ylabel("% of trials"); axes[1, j].set_xlabel("GHz")
        for ax in axes[:, j]:
            ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=6.5)
        axes[1, j].set_ylim(-0.5, 12)
    fig.tight_layout(); fig.savefig(OUT / f"p95_vs_freq{sfx}.png", dpi=130); plt.close(fig)
    # convergence
    if histories:
        fig, axes = plt.subplots(1, len(histories), figsize=(4 * len(histories), 3.3))
        for ax, (n, hs) in zip(np.atleast_1d(axes), histories.items()):
            for k, h in hs.items():
                ax.plot([x["evals"] for x in h], [x["best"] for x in h], label=k, lw=1)
            ax.set_yscale("log"); ax.set_title(f"N={n}: best score per generation (noisy)", fontsize=9)
            ax.set_xlabel("evaluations"); ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=7)
        fig.tight_layout(); fig.savefig(OUT / f"convergence{sfx}.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
