"""Evolutionary search: change a design a little, keep what works, repeat.

`climb` is the plain version of that idea: one design, a few random small changes per step, and a change
is kept only if it scores better than the current design on the same random numbers. It is a (1+lambda)
evolution strategy; its step size grows after a success and shrinks after a failure (the 1/5 rule), so it
takes big steps far from an optimum and small ones near it.

`evolve` keeps a population instead: a real-coded genetic algorithm with tournament selection, blend
crossover, self-adaptive Gaussian mutation and elitism. Each generation costs more evaluations, but a
population is less likely to stall on the first local optimum, which matters when the score has many
(sparse array layouts, for instance).

Both minimise `score(X, gen)`, which scores a whole batch of designs (the rows of X) at once. A problem
can then vectorise, and draw its random numbers from `gen`: every design in one generation is judged on
the same noise (common random numbers), and the noise changes between generations so that no design
survives by being lucky on one draw. For the same reason elites are re-scored every generation instead of
keeping their old score.

What this cannot do: it optimises the model, not the hardware. A design that exploits a model error scores
well and builds badly. Re-check every winner with an independent model and fresh random numbers, and
finally with a measurement.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

Score = Callable[[np.ndarray, int], np.ndarray]


@dataclass
class Result:
    x: np.ndarray                 # best design of the last generation (or step)
    score: float                  # its score on that generation's random numbers
    history: list = field(default_factory=list)   # one dict per generation: gen, best, median, sigma, evals
    evaluations: int = 0
    finalists: np.ndarray | None = None   # the last population (or the best designs seen), best first


def reselect(score: Score, candidates, *, draws: int = 8, first_gen: int = 1_000_000) -> tuple[np.ndarray, np.ndarray]:
    """Re-score finalists on `draws` fresh generations of random numbers and rank them by their mean.

    With a noisy score the best design of the last generation is partly the luckiest one (the winner's
    curse); averaging over new draws, which no design has been selected on, removes that bias.
    Returns (candidates sorted best first, their mean scores)."""
    X = np.atleast_2d(np.asarray(candidates, float))
    F = np.mean([_batch_score(score, X, first_gen + k) for k in range(draws)], axis=0)
    order = np.argsort(F, kind="stable")
    return X[order], F[order]


def reflect(x, lo, hi):
    """Fold values that left [lo, hi] back inside, as a mirror would; clip anything still outside."""
    y = np.where(x < lo, 2 * lo - x, x)
    y = np.where(y > hi, 2 * hi - y, y)
    return np.clip(y, lo, hi)


def blend(a, b, rng, alpha=0.5):
    """BLX-alpha crossover: each coordinate drawn uniformly from the parents' interval widened by alpha."""
    u = rng.uniform(-alpha, 1 + alpha, a.shape)
    return a + u * (b - a)


def _batch_score(score: Score, X: np.ndarray, gen: int) -> np.ndarray:
    f = np.asarray(score(X, gen), dtype=float).reshape(-1)
    if f.shape[0] != X.shape[0]:
        raise ValueError(f"score returned {f.shape[0]} values for {X.shape[0]} designs")
    return np.where(np.isfinite(f), f, np.inf)


def evolve(score: Score, lo, hi, *, pop: int = 32, generations: int = 100, elite: int = 2,
           tournament: int = 3, crossover: float = 0.7, sigma: float = 0.1, seed: int = 0,
           init=None, repair: Callable | None = None, cross: Callable | None = None,
           mutate: Callable | None = None, callback: Callable | None = None) -> Result:
    """Minimise score over the box [lo, hi] with a population of `pop` designs.

    sigma: initial mutation step as a fraction of each coordinate's range; every design carries its own
    step, which is inherited, perturbed log-normally and so adapts to the landscape.
    init: designs to start from (one or several rows), e.g. a design you already like. The rest of the
    first population is half mutated copies of them, half random.
    repair(x) -> x: maps any design onto the feasible set (spacing rules, keep-outs); applied to every child.
    cross(a, b, rng) -> child and mutate(x, step, rng) -> x replace the default blend crossover and
    Gaussian mutation (step is in units of the range).
    """
    rng = np.random.default_rng(seed)
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    span = hi - lo
    d = lo.size
    if pop < 2 or not 0 <= elite < pop:
        raise ValueError("need pop >= 2 and 0 <= elite < pop")
    cross = cross or blend
    X = lo + rng.random((pop, d)) * span
    if init is not None:
        seeds = np.atleast_2d(np.asarray(init, float))[:pop]
        X[:len(seeds)] = seeds
        n_mut = (pop - len(seeds)) // 2
        for i in range(len(seeds), len(seeds) + n_mut):
            X[i] = reflect(seeds[i % len(seeds)] + sigma * span * rng.standard_normal(d), lo, hi)
    if repair is not None:
        X = np.array([repair(x) for x in X])
    S = np.full(pop, float(sigma))
    F = _batch_score(score, X, 0)
    evals = pop
    tau = 1 / np.sqrt(d)
    history = [dict(gen=0, best=float(F.min()), median=float(np.median(F)), sigma=float(sigma), evals=evals)]

    def pick():
        idx = rng.integers(0, pop, tournament)
        return idx[np.argmin(F[idx])]

    for gen in range(1, generations + 1):
        keep = np.argsort(F, kind="stable")[:elite]
        kids, ksig = [], []
        while len(kids) < pop - elite:
            a = pick()
            x, s = X[a].copy(), S[a]
            if rng.random() < crossover:
                b = pick()
                x = cross(X[a], X[b], rng)
                s = np.sqrt(S[a] * S[b])
            s = float(np.clip(s * np.exp(tau * rng.standard_normal()), 1e-7, 0.5))
            x = mutate(x, s, rng) if mutate else x + s * span * rng.standard_normal(d)
            x = reflect(x, lo, hi)
            if repair is not None:
                x = repair(x)
            kids.append(x)
            ksig.append(s)
        X = np.vstack([X[keep], np.array(kids)])
        S = np.concatenate([S[keep], ksig])
        F = _batch_score(score, X, gen)          # elites re-scored on this generation's noise
        evals += pop
        i = int(np.argmin(F))
        history.append(dict(gen=gen, best=float(F[i]), median=float(np.median(F[np.isfinite(F)])) if np.isfinite(F).any() else np.inf,
                            sigma=float(np.median(S)), evals=evals))
        if callback is not None:
            callback(gen, X[i], F[i])
    order = np.argsort(F, kind="stable")
    return Result(X[order[0]].copy(), float(F[order[0]]), history, evals, finalists=X[order].copy())


def climb(score: Score, x0, lo, hi, *, steps: int = 200, lam: int = 8, sigma: float = 0.05, seed: int = 0,
          repair: Callable | None = None, mutate: Callable | None = None,
          callback: Callable | None = None) -> Result:
    """Start from a design you like; each step try `lam` small random changes and keep the best of them only
    if it beats the current design on the same random numbers (the current design is re-scored every step).
    The step grows by e^0.2 after a success and shrinks by e^-0.05 after a failure: steps settle where about
    one in five succeeds."""
    rng = np.random.default_rng(seed)
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    span = hi - lo
    x = np.asarray(x0, float).copy()
    if repair is not None:
        x = repair(x)
    s = float(sigma)
    evals = 0
    history = []
    for step in range(steps):
        M = np.empty((lam, x.size))
        for j in range(lam):
            m = mutate(x, s, rng) if mutate else x + s * span * rng.standard_normal(x.size)
            m = reflect(m, lo, hi)
            M[j] = repair(m) if repair is not None else m
        F = _batch_score(score, np.vstack([x, M]), step)
        evals += lam + 1
        j = int(np.argmin(F[1:]))
        if F[1 + j] < F[0]:
            x, best, s = M[j], float(F[1 + j]), min(s * np.exp(0.2), 0.5)
        else:
            best, s = float(F[0]), max(s * np.exp(-0.05), 1e-7)
        history.append(dict(gen=step, best=best, median=float(np.median(F)), sigma=s, evals=evals))
        if callback is not None:
            callback(step, x, best)
    return Result(x.copy(), history[-1]["best"] if history else np.nan, history, evals, finalists=x[None].copy())


def random_search(score: Score, lo, hi, *, batches: int = 100, batch: int = 32, seed: int = 0,
                  repair: Callable | None = None) -> Result:
    """The baseline both methods must beat for the same number of evaluations: score random designs and keep
    the best (re-scored against each new batch so it is judged on the same noise)."""
    rng = np.random.default_rng(seed)
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    keep = max(1, batch // 4)
    pool, evals, history = None, 0, []          # the best `keep` designs so far, re-scored with each batch
    for b in range(batches):
        X = lo + rng.random((batch, lo.size)) * (hi - lo)
        if repair is not None:
            X = np.array([repair(x) for x in X])
        if pool is not None:
            X[:len(pool)] = pool
        F = _batch_score(score, X, b)
        evals += batch
        order = np.argsort(F, kind="stable")
        pool = X[order[:keep]].copy()
        history.append(dict(gen=b, best=float(F[order[0]]), median=float(np.median(F)), sigma=np.nan, evals=evals))
    return Result(pool[0].copy(), history[-1]["best"], history, evals, finalists=pool)
