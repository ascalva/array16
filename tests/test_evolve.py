"""The evolutionary optimiser on problems with known answers."""
import numpy as np

from array16.evolve import climb, evolve, random_search, reflect, reselect


def sphere(target):
    return lambda X, gen: np.sum((X - target) ** 2, axis=1)


def rastrigin(X, gen):
    return 10 * X.shape[1] + np.sum(X ** 2 - 10 * np.cos(2 * np.pi * X), axis=1)


def test_reflect_keeps_inside():
    lo, hi = np.zeros(3), np.ones(3)
    x = reflect(np.array([-0.2, 1.3, 0.5]), lo, hi)
    assert np.allclose(x, [0.2, 0.7, 0.5])
    assert np.all((reflect(np.array([-5.0, 9.0, 0.5]), lo, hi) >= 0) & (reflect(np.array([-5.0, 9.0, 0.5]), lo, hi) <= 1))


def test_evolve_finds_shifted_sphere():
    target = np.array([0.3, -1.2, 2.0, 0.7, -0.4, 1.1])
    r = evolve(sphere(target), -3 * np.ones(6), 3 * np.ones(6), pop=30, generations=150, seed=1)
    assert np.max(np.abs(r.x - target)) < 0.01
    assert r.evaluations == 30 * 151


def test_evolve_escapes_local_minima():
    # Rastrigin in 4-D has ~10^4 local minima; the global one is at 0 with score 0.
    r = evolve(rastrigin, -5.12 * np.ones(4), 5.12 * np.ones(4), pop=60, generations=200, seed=3)
    assert r.score < 1.0


def test_climb_never_gets_worse_on_a_fixed_score():
    target = np.array([1.0, -2.0, 0.5])
    r = climb(sphere(target), np.zeros(3), -4 * np.ones(3), 4 * np.ones(3), steps=300, lam=6, seed=2)
    best = [h["best"] for h in r.history]
    assert all(b2 <= b1 + 1e-12 for b1, b2 in zip(best, best[1:]))
    assert np.max(np.abs(r.x - target)) < 0.01


def test_same_seed_same_answer():
    f = sphere(np.array([0.5, 0.5]))
    a = evolve(f, -np.ones(2), np.ones(2), pop=10, generations=20, seed=7)
    b = evolve(f, -np.ones(2), np.ones(2), pop=10, generations=20, seed=7)
    assert np.array_equal(a.x, b.x) and a.score == b.score


def test_elites_are_rescored_on_new_noise():
    # A design that was lucky once must not keep its lucky score: every design is scored every generation.
    calls = []

    def noisy(X, gen):
        calls.append(len(X))
        rng = np.random.default_rng(gen)
        return np.sum(X ** 2, axis=1) + rng.normal(0, 0.1, len(X))

    evolve(noisy, -np.ones(2), np.ones(2), pop=8, generations=5, elite=2, seed=0)
    assert calls == [8] * 6


def test_repair_is_applied_to_every_design():
    seen = []

    def score(X, gen):
        seen.append(X.copy())
        return np.sum(X, axis=1)

    def repair(x):                       # feasible set: x0 >= x1
        return np.array([max(x[0], x[1]), min(x[0], x[1])])

    evolve(score, np.zeros(2), np.ones(2), pop=12, generations=10, repair=repair, seed=4)
    X = np.vstack(seen)
    assert np.all(X[:, 0] >= X[:, 1])
    assert np.all((X >= 0) & (X <= 1))


def test_init_seeds_the_population():
    liked = np.array([0.9, 0.9])
    first = []

    def score(X, gen):
        if gen == 0:
            first.append(X.copy())
        return np.sum((X - 1.0) ** 2, axis=1)

    evolve(score, np.zeros(2), 2 * np.ones(2), pop=10, generations=3, init=liked, seed=5)
    X0 = first[0]
    assert np.allclose(X0[0], liked)                                    # the liked design itself
    near = np.linalg.norm(X0[1:] - liked, axis=1) < 0.6                  # about half are small changes of it
    assert near.sum() >= 4


def test_both_methods_beat_random_search_for_the_same_budget():
    target = np.linspace(-1, 1, 8)
    f = sphere(target)
    lo, hi = -2 * np.ones(8), 2 * np.ones(8)
    ga = evolve(f, lo, hi, pop=20, generations=49, seed=6)
    hc = climb(f, np.zeros(8), lo, hi, steps=111, lam=8, seed=6)
    rs = random_search(f, lo, hi, batches=50, batch=20, seed=6)
    assert ga.evaluations == rs.evaluations == 1000 and hc.evaluations == 999
    assert ga.score < rs.score / 10 and hc.score < rs.score / 10


def test_reselect_removes_the_winners_curse():
    # Ten designs of equal true quality plus one truly better one; one noisy draw often crowns a lucky
    # equal, the mean over fresh draws finds the better one.
    X = np.vstack([np.zeros((10, 1)), [[1.0]]])
    true = np.where(X[:, 0] > 0.5, 0.0, 0.3)

    def noisy(Xb, gen):
        rng = np.random.default_rng(gen)
        return true[:len(Xb)] + rng.normal(0, 0.5, len(Xb))

    lucky = sum(int(np.argmin(noisy(X, g)) != 10) for g in range(200)) / 200
    ranked, means = reselect(noisy, X, draws=64)
    assert lucky > 0.5                     # a single draw usually picks a merely lucky design
    assert ranked[0, 0] == 1.0 and means[0] < means[1]


def test_finalists_are_returned_best_first():
    f = sphere(np.array([0.2, -0.3]))
    r = evolve(f, -np.ones(2), np.ones(2), pop=12, generations=15, seed=8)
    assert r.finalists.shape == (12, 2) and np.allclose(r.finalists[0], r.x)
    rs = random_search(f, -np.ones(2), np.ones(2), batches=5, batch=8, seed=8)
    assert len(rs.finalists) == 2 and np.allclose(rs.finalists[0], rs.x)
