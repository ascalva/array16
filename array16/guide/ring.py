"""The ring resonator: a closed guide as a ring of transmission line, its resonances, its feeds and its tolerances.

docs/KNOTTED_GUIDE.md sections 5.2 and 5.3, taken from the knotted study's typed ports (experiments/
knotted_guide_study.py, its drive and ring sections) and held to their numbers. Two forms of one line model, each
segment losing beta / (2 Q0):

  nodal admittance  N segments closed on themselves, segment k from node k to node k + 1 with its own Z0 and index,
                    fed by ideal currents at nodes (ring_admittance, ring_response), the buildability lens's model
                    (lens5/ring.py, ring915.py): the resonances of the +-m pair and their splitting 2Q|c_2m|
                    (ring_resonances, ring_harmonic), the backward wave |b/a| (ring_purity), the weakest depth per
                    watt under the fixed 90 degree hybrid or two re-solved feeds (ring_score), K feeds compiled for
                    a travelling wave (compile_feeds), and seeded tolerance draws: the spacers' gaps, the two closing
                    joints' impedance and position, the split's etch and the laminate (ring_draw), with T915-S's
                    Monte Carlo against a perfect ring (RING_CASES, tolerance_run) and the 2.45 GHz strips as the
                    lens ran them (ring_draw_stations, strips_2450). A build's figures are a Build: T915S's from
                    section B as solved, LENS the buildability lens's, whose dn / d eps_r, 0.0253, is 2.2 times the
                    solved section's n share / (2 eps_r).
  line matrices     Ring: chain matrices closed round the ring, joints as short sections of another impedance, feeds
                    as Norton currents through a coupling capacitor or a tuned conductance, V(s_i) = sum_k G_ik w_k
                    at samples (Ring.operator); K feeds compiled to a travelling wave by regularised least squares
                    with the target's phase relaxed (solve_weights, best_drive) against an untuned hybrid; the drive
                    study's rings (RingCase, drive_cases).

The backward wave, in coupled modes (T. J. Kippenberg, S. M. Spillane and K. J. Vahala, Opt. Lett. 27, 1669 (2002)):
with Q = omega tau, tau the mode's whole lifetime, its linewidth (full width) is f / Q; the +-m pair's coupling splits
it by 2|c_2m| of f, and driven at the pair's centre in one sense the ring carries |b/a| = 2Q|c_2m|, the splitting in
linewidths. A lumped reflection S splits it by |S| / (pi m) of f, so |b/a| = |S| Q / (pi m) (tests/test_guide_ring.py
holds it to 1 %); section 5.2's |S| Q / (2 pi m) is that with Q counted on the half linewidth.

SI units: frequencies in Hz, impedances in ohm, node volts per ampere of feed current.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any, NamedTuple

import numpy as np
from numpy.typing import NDArray

from array16.particles import C0, Arr, CArr

# T915-S (docs/KNOTTED_GUIDE.md section 9, item 1) and the drive lens's section A
F915 = 915e6  # Hz
N915 = 1.031  # section B's index at 915 MHz, its 18 um copper and 25 um polyimide resolved (the line study)
A_INDEX = 1.13  # section A's index as drawn (the drive lens's ring)
RING_NODES = 192
RING_SENSITIVITY = 0.50  # d ln Z0 / d ln h, the buildability lens's (section 5.1)
SPACERS = 24
SCARF = 29.4e-3  # m: three secular wavelengths at 8 eV (q = 0.5, 915 MHz), the closing joints' scarf
MUS = (0.0, 1e-3, 1e-2, 3e-2, 1e-1, 3e-1)  # the drive compile's regularisations


# ---- nodal admittance: the ring of segments -----------------------------------------------------------------------
def _closed(own: NDArray[Any], mutual: NDArray[Any]) -> NDArray[Any]:
    """The N x N matrix of a ring of two-ports, segment k (self term `own`, mutual term `mutual`) from node k to
    node k + 1."""
    count = len(own)
    out = np.zeros((count, count), dtype=own.dtype)
    idx = np.arange(count)
    nxt = (idx + 1) % count
    np.add.at(out, (idx, idx), own)
    np.add.at(out, (nxt, nxt), own)
    np.add.at(out, (idx, nxt), mutual)
    np.add.at(out, (nxt, idx), mutual)
    return out


def ring_admittance(f: complex, length: float, z: Arr, n: Arr, q0: float) -> CArr:
    """The nodal admittance (N x N, S) of a closed ring of N line segments, segment k from node k to node k + 1, of
    impedance z_k and index n_k, losing beta / (2 Q0) (the buildability lens's ring, lens5/ring.py). The frequency
    may be complex; Q0 = inf is the lossless ring."""
    count = len(z)
    beta = 2 * math.pi * f * n / C0
    g = (beta / (2 * q0) + 1j * beta) * length / count
    out: CArr = _closed(1 / (z * np.tanh(g)), -1 / (z * np.sinh(g)))
    return out


def ring_response(
    f: complex, length: float, z: Arr, n: Arr, q0: float, feeds: Sequence[int], source: complex = 0.0
) -> CArr:
    """G (N, K): the node volts per ampere injected at each feed node, ideal current feeds, or Norton feeds of source
    admittance `source` (S) at each feed node. Columns of the inverse nodal admittance."""
    y = ring_admittance(f, length, z, n, q0)
    at = np.asarray(feeds, dtype=np.intp)
    np.add.at(y, (at, at), source)
    rhs = np.zeros((len(z), len(at)), dtype=np.complex128)
    rhs[at, np.arange(len(at))] = 1.0
    out: CArr = np.linalg.solve(y, rhs)
    return out


def perfect_ring(index: float, nodes: int = RING_NODES, z0: float = 50.0) -> tuple[Arr, Arr, Arr]:
    """A ring as drawn: Z (ohm), n and dh/h at its nodes, uniform."""
    return np.full(nodes, z0), np.full(nodes, index), np.zeros(nodes)


def ring_purity(v: CArr, m: int) -> float:
    """|b/a|: the smaller over the larger of the node voltages' azimuthal components of order -m and +m."""
    k = np.arange(len(v))
    a = abs(np.sum(v * np.exp(2j * math.pi * m * k / len(v))))
    b = abs(np.sum(v * np.exp(-2j * math.pi * m * k / len(v))))
    return float(min(a, b) / max(a, b))


def ring_harmonic(z: Arr, m: int = 1, z0: float = 50.0) -> float:
    """|c_2m|: the lap harmonic 2m of ln Z, whose twice splits +-m by that fraction of the frequency."""
    k = np.arange(len(z))
    return float(abs(np.sum(np.log(z / z0) * np.exp(2j * math.pi * 2 * m * k / len(z))) / len(z)))


@dataclass(frozen=True)
class Pair:
    """The +-m pair's resonances as complex frequencies f (1 + i / (2Q)) (Hz), lower and upper in their real
    parts."""

    lower: complex
    upper: complex

    @property
    def f(self) -> float:
        """Their centre (Hz)."""
        return (self.lower.real + self.upper.real) / 2

    @property
    def q(self) -> float:
        """Their Q: the centre over twice their mean imaginary part."""
        return self.f / (self.lower.imag + self.upper.imag)

    @property
    def splitting(self) -> float:
        """Delta f / f: 2|c_2m| to first order (Little, Laine and Chu)."""
        return (self.upper.real - self.lower.real) / self.f

    @property
    def linewidths(self) -> float:
        """The splitting in linewidths f / Q: 2Q|c_2m| to first order."""
        return self.q * self.splitting


def _susceptance(f: float, length: float, z: Arr, n: Arr) -> tuple[Arr, Arr]:
    """B and dB/df (S, S/Hz) of the lossless ring, whose nodal admittance is -iB with B real and symmetric."""
    theta = 2 * math.pi * f * n / C0 * length / len(z)
    cot, csc = 1 / np.tan(theta), 1 / np.sin(theta)
    rate = theta / f
    return _closed(cot / z, -csc / z), _closed(-(csc**2) / z * rate, csc * cot / z * rate)


def _offsets(f: float, length: float, z: Arr, n: Arr) -> Arr:
    """The lossless ring's modes as offsets from f (Hz) to first order: B(f) v = d (-dB/df) v, a definite pencil
    (-dB/df is d Im Y / df, positive definite: each segment's susceptance rises with frequency, Foster's theorem)."""
    b, db = _susceptance(f, length, z, n)
    li = np.linalg.inv(np.linalg.cholesky(-db))
    out: Arr = np.linalg.eigvalsh(li @ b @ li.T)
    return out


def ring_resonances(length: float, z: Arr, n: Arr, q0: float, m: int = 1, iters: int = 8) -> Pair:
    """The +-m pair of a ring of segments: the complex frequencies where its nodal admittance is singular. The ring
    loses beta / (2 Q0) on every segment, so it is the lossless ring at f (1 - i / (2 Q0)): the pair is found on the
    lossless ring, from its two modes nearest m c / (n L) to first order, each refined by successive linear problems
    (Ruhe's method, quadratic) until a step is under 1e-11 of f, and carried over, its Q being Q0. The roots settle
    to about 1e-16 N^2 / m^2 of f (the nodal matrix's conditioning: 2e-14 at 192 nodes, 1e-12 at 960); one that has
    not settled to 1e-9 in `iters` steps raises ArithmeticError."""
    f0 = m * C0 / (float(np.mean(n)) * length)
    first = _offsets(f0, length, z, n)
    roots: list[float] = []
    for step in np.sort(first[np.argsort(np.abs(first))[:2]]):
        f, move = f0 + float(step), 0.0
        for _ in range(iters):
            near = _offsets(f, length, z, n)
            move = float(near[int(np.argmin(np.abs(near)))])
            f += move
            if abs(move) <= 1e-11 * f:
                break
        if abs(move) > 1e-9 * f:
            raise ArithmeticError(
                f"the +-{m} pair has not settled in {iters} steps: the last moved {move / f:.1e} of f"
            )
        roots.append(f)
    lower, upper = sorted(roots)
    carry = 1 / (1 - 0.5j / q0)
    return Pair(lower * carry, upper * carry)


@dataclass(frozen=True, eq=False)
class Score:
    """A drive of the ring of segments as scored: its frequency (Hz) and the ring's nominal m c / (n L); the weakest
    depth per watt (V^2 / W, the depth going as |V|^2 / (1 + dh/h)^2); |b/a|; the weakest depth over the strongest
    (DR-6's min/max); the feeds' currents (A per ampere of the first), the node volts and the depths."""

    f: float
    f0: float
    score: float
    purity: float
    depth_ratio: float
    weights: CArr
    v: CArr
    depth: Arr

    @property
    def w_abs(self) -> float:
        """The second feed's current over the first's, in amplitude."""
        return float(abs(self.weights[1] / self.weights[0]))

    @property
    def w_deg(self) -> float:
        """The second feed's current over the first's, in phase (deg)."""
        w = complex(self.weights[1] / self.weights[0])
        return math.degrees(math.atan2(w.imag, w.real))


def ring_score(
    length: float,
    m: int,
    z: Arr,
    n: Arr,
    rel: Arr,
    q0: float,
    compensate: bool,
    feeds: tuple[int, int] | None = None,
) -> Score:
    """The ring driven by two current feeds a quarter guided wave apart (nodes 0 and N / 4m unless given), the
    second at w times the first: over 49 frequencies within 3/Q0 of the resonance and over w = +-i (the fixed 90
    degree hybrid) or, compensated, 13 amplitudes from 0.5 to 2 by 37 phases (two re-solved feeds), the best weakest
    depth per watt, the depth going as |V|^2 / (1 + dh/h)^2 (lens5/ring.py's evaluate, vectorised over w)."""
    count = len(z)
    at = (0, count // (4 * m)) if feeds is None else feeds
    f0 = m * C0 / (float(np.mean(n)) * length)
    if compensate:
        ws = np.array(
            [
                a * np.exp(1j * p)
                for a in np.linspace(0.5, 2.0, 13)
                for p in np.linspace(-math.pi, math.pi, 37, endpoint=False)
            ]
        )
    else:
        ws = np.array([1j, -1j])
    best: tuple[float, float, complex, CArr, Arr] | None = None
    for f in np.linspace(f0 * (1 - 3 / q0), f0 * (1 + 3 / q0), 49):
        inv = np.linalg.inv(ring_admittance(float(f), length, z, n, q0))
        v = inv[:, [at[0]]] + inv[:, [at[1]]] * ws[None, :]
        power = 0.5 * np.real(v[at[0]] + np.conj(ws) * v[at[1]])
        depth = np.abs(v) ** 2 / (1 + rel[:, None]) ** 2
        score = depth.min(axis=0) / power
        k = int(np.argmax(score))
        if best is None or float(score[k]) > best[0]:
            best = (float(score[k]), float(f), complex(ws[k]), v[:, k], depth[:, k])
    assert best is not None
    score_best, f_best, w, v_best, depth_best = best
    return Score(
        f=f_best,
        f0=f0,
        score=score_best,
        purity=ring_purity(v_best, m),
        depth_ratio=float(depth_best.min() / depth_best.max()),
        weights=np.array([1.0, w], dtype=np.complex128),
        v=v_best,
        depth=depth_best,
    )


def compile_feeds(
    length: float,
    m: int,
    z: Arr,
    n: Arr,
    rel: Arr,
    q0: float,
    feeds: Sequence[int],
    f: float | None = None,
    mu: float = 1e-3,
) -> Score:
    """K ideal current feeds at nodes, their currents re-solved for a travelling wave of order m from the ring's own
    G at frequency f (the pair's centre when None): the drive study's compile (solve_weights) on the ring of
    segments, scored as ring_score scores a drive."""
    at = tuple(feeds)
    if f is None:
        f = ring_resonances(length, z, n, q0, m).f
    g = ring_response(f, length, z, n, q0, at)
    target = np.exp(-2j * math.pi * m * np.arange(len(z)) / len(z))
    w = solve_weights(g, target, mu)
    w = w / w[0]
    v: CArr = g @ w
    power = 0.5 * float(np.real(np.sum(v[list(at)] * np.conj(w))))
    depth: Arr = np.abs(v) ** 2 / (1 + rel) ** 2
    return Score(
        f=f,
        f0=m * C0 / (float(np.mean(n)) * length),
        score=float(depth.min()) / power,
        purity=ring_purity(v, m),
        depth_ratio=float(depth.min() / depth.max()),
        weights=w,
        v=v,
        depth=depth,
    )


# ---- tolerance draws and the Monte Carlo --------------------------------------------------------------------------
@dataclass(frozen=True)
class Build:
    """What a built ring's tolerances act through, T915-S's by default: the ring's nodes; the spacer stations and
    d ln Z0 / d ln h; the two closing joints' scarf and the sd of their spacing (m); the split's d ln Z / d ln g and
    its stations; the laminate's dn / d eps_r, section B's solved n share / (2 eps_r) = 1.0311 x 0.0737 / 6.6 (DR-0's
    share rule; array16.twin.TwinGuide.build takes it from each design's solve); the nominal Z0 (ohm)."""

    nodes: int = RING_NODES
    spacers: int = SPACERS
    sensitivity: float = RING_SENSITIVITY
    scarf: float = SCARF
    joint_jitter: float = 0.5e-3
    split_sensitivity: float = 0.03
    split_stations: int = 6
    index_per_eps: float = 0.01151
    z0: float = 50.0


T915S = Build()
LENS = Build(index_per_eps=0.0253)  # the buildability lens's figures (lens5/ring.py), its laminate's dn / d eps_r


def ring_draw(
    rng: np.random.Generator,
    length: float,
    sigma: float,
    index: float,
    joint_sd: float,
    split_sd: float,
    joint_dz: float = -0.02,
    sigma_eps: float = 0.1,
    build: Build = T915S,
) -> tuple[Arr, Arr, Arr]:
    """One built ring (ring915.py's realise), by the build's figures (T915-S's by default): the spacer stations at a
    random offset, each with a relative gap error of sd sigma, linear between them; ln Z = S dh/h (S = 0.5); a slow
    split-width error of sd split_sd at the split's stations, d ln Z / d ln g about 0.03; two closing joints a scarf
    long, a quarter of the length (a quarter guided wave at m = 1, +-0.5 mm) apart at a random place, each dz =
    joint_dz + N(0, joint_sd) in ln Z; the laminate's permittivity off by N(0, sigma_eps), dn / d eps_r the build's
    (T915-S's 0.0115; the lens's 0.0253). Returns Z (ohm), n and dh/h at the nodes."""
    s = np.arange(build.nodes) * length / build.nodes
    st = np.arange(build.spacers) * length / build.spacers + rng.uniform(0, length / build.spacers)
    st = np.sort(st % length)
    dh = rng.normal(0, sigma, build.spacers)
    rel: Arr = np.interp(s, np.concatenate([st - length, st, st + length]), np.concatenate([dh, dh, dh]))
    lnz = build.sensitivity * rel
    if split_sd > 0:
        k = rng.normal(0, split_sd, build.split_stations)
        stations = np.linspace(0, length, build.split_stations + 1)
        lnz = lnz + build.split_sensitivity * np.interp(s, stations, np.concatenate([k, k[:1]]))
    s0 = rng.uniform(0, length)
    for base in (s0, s0 + length / 4 + rng.normal(0, build.joint_jitter)):
        dz = joint_dz + rng.normal(0, joint_sd)
        lnz = lnz + np.where(np.abs(((s - base + length / 2) % length) - length / 2) < build.scarf / 2, dz, 0.0)
    z: Arr = build.z0 * np.exp(lnz)
    n: Arr = np.full(build.nodes, index * (1 + build.index_per_eps * rng.normal(0, sigma_eps) / index))
    return z, n, rel


def ring_draw_stations(
    rng: np.random.Generator, length: float, sigma: float, index: float, sigma_eps: float = 0.1, build: Build = LENS
) -> tuple[Arr, Arr, Arr]:
    """lens5/ring.py's realise for the 2.45 GHz strips, without joints: the spacer stations at fixed places (the
    first at s = 0), each with a relative gap error of sd sigma, linear between them; ln Z = S dh/h; the laminate's
    permittivity off by N(0, sigma_eps), dn / d eps_r the build's (the lens's 0.0253 by default, as it ran them)."""
    s = np.arange(build.nodes) * length / build.nodes
    st = np.arange(build.spacers) * length / build.spacers
    dh = rng.normal(0, sigma, build.spacers)
    rel: Arr = np.interp(s, np.concatenate([st, [length]]), np.concatenate([dh, dh[:1]]))
    z: Arr = build.z0 * np.exp(build.sensitivity * rel)
    n: Arr = np.full(build.nodes, index * (1 + build.index_per_eps * rng.normal(0, sigma_eps) / index))
    return z, n, rel


class ToleranceCase(NamedTuple):
    """A case of T915-S's Monte Carlo (section 5.3): its label; its finish, a key of the line study's Q ("bare",
    "enig 1.3", "enig 5") or a Q; the spacers' rms relative gap error; draws; seed; the joints' sd in ln Z; the
    split's etch sd; the index it is drawn with (None: N915). It unpacks as the study's tuples did."""

    label: str
    finish: str | float
    sigma: float
    draws: int
    seed: int
    joint_sd: float
    split_sd: float
    n_drawn: float | None

    @property
    def n_eff(self) -> float:
        return N915 if self.n_drawn is None else self.n_drawn

    @property
    def length(self) -> float:
        """m = 1 at 915 MHz (m)."""
        return C0 / (self.n_eff * F915)


RING_CASES = (
    ToleranceCase("bare copper, 2 %", "bare", 0.02, 60, 7, 0.01, 0.0, None),
    ToleranceCase("bare copper, 1 %", "bare", 0.01, 30, 8, 0.01, 0.0, None),
    ToleranceCase("bare copper, 3 %", "bare", 0.03, 30, 9, 0.01, 0.0, None),
    ToleranceCase("ENIG high, 2 %", "enig 1.3", 0.02, 30, 10, 0.01, 0.0, None),
    ToleranceCase("ENIG low, 2 %", "enig 5", 0.02, 30, 11, 0.01, 0.0, None),
    ToleranceCase("bare copper, 2 %, split etched 20 %", "bare", 0.02, 30, 12, 0.01, 0.2, None),
    ToleranceCase("joints alone, bare copper", "bare", 0.0, 30, 14, 0.01, 0.0, None),
    ToleranceCase("as drawn before the revision: n 1.067, Q 338, 2 %", 338.0, 0.02, 30, 13, 0.01, 0.0, 1.067),
)


@dataclass(frozen=True, eq=False)
class ToleranceRun:
    """One case's draws, each scored against a perfect ring of the same Q and index: the splitting 2Q|c_2m|
    (linewidths); the fixed hybrid's |b/a| and depth per watt over the perfect ring's; the two re-solved feeds' the
    same; the re-solved second feed's amplitude and phase (deg); the perfect ring's |b/a| under the hybrid."""

    label: str
    q0: float
    length: float
    split: Arr
    hybrid_purity: Arr
    hybrid_depth: Arr
    resolved_purity: Arr
    resolved_depth: Arr
    weight_abs: Arr
    weight_deg: Arr
    perfect_purity: float

    def summary(self) -> dict[str, float]:
        """Medians, 90th percentiles of the splitting and the backward waves, 10th of the depths, the weights'
        ranges: the rows of section 5.3's table."""
        out = {"Q0": self.q0, "draws": float(len(self.split)), "length_mm": self.length * 1e3}
        for name in ("split", "hybrid_purity", "resolved_purity"):
            col = getattr(self, name)
            out[f"{name}_median"], out[f"{name}_p90"] = float(np.median(col)), float(np.percentile(col, 90))
        for name in ("hybrid_depth", "resolved_depth"):
            col = getattr(self, name)
            out[f"{name}_median"], out[f"{name}_p10"] = float(np.median(col)), float(np.percentile(col, 10))
        for name in ("weight_abs", "weight_deg"):
            col = getattr(self, name)
            out[f"{name}_min"], out[f"{name}_max"] = float(col.min()), float(col.max())
        out["perfect_purity"] = self.perfect_purity
        return out


def tolerance_run(case: ToleranceCase, q0: float, draws: int | None = None, build: Build = T915S) -> ToleranceRun:
    """The case's Monte Carlo at Q0: its draws (the first `draws` of them when given) from its seed, each driven by
    the fixed hybrid and by two re-solved feeds (ring_score), against the perfect ring; by T915-S's figures, or LENS's
    to repeat the knotted study's draws. Each ring is scored about its own resonance, so the laminate's figure moves
    the draws' tuned bands and not their scores."""
    rng = np.random.default_rng(case.seed)
    count = case.draws if draws is None else draws
    length = case.length
    built = [
        ring_draw(rng, length, case.sigma, case.n_eff, case.joint_sd, case.split_sd, build=build) for _ in range(count)
    ]
    perfect = ring_score(length, 1, *perfect_ring(case.n_eff, build.nodes, build.z0), q0, False)
    rows = []
    for z, n, rel in built:
        fixed = ring_score(length, 1, z, n, rel, q0, False)
        comp = ring_score(length, 1, z, n, rel, q0, True)
        rows.append(
            (
                2 * q0 * ring_harmonic(z, 1, build.z0),
                fixed.purity,
                fixed.score / perfect.score,
                comp.purity,
                comp.score / perfect.score,
                comp.w_abs,
                comp.w_deg,
            )
        )
    cols = np.array(rows, dtype=np.float64).reshape(count, 7).T
    return ToleranceRun(
        case.label,
        q0,
        length,
        split=cols[0],
        hybrid_purity=cols[1],
        hybrid_depth=cols[2],
        resolved_purity=cols[3],
        resolved_depth=cols[4],
        weight_abs=cols[5],
        weight_deg=cols[6],
        perfect_purity=perfect.purity,
    )


def strips_2450(q0: float, index: float, draws: int) -> dict[str, float]:
    """The 2.45 GHz strips' tolerance as the buildability lens ran it (lens5/ring.py): m = 2, 24 spacers at fixed
    stations with 10 um rms gap error on its h = 0.8391e6 / 2.45e9 m, no joints, seed 1, the 12 draws at 5 um drawn
    first; the median 2Q|c_2m|, and the fixed hybrid's and two re-solved feeds' depth per watt over a perfect ring's."""
    h = 0.8391e6 / 2.45e9
    length = 2 * C0 / (index * 2450e6)
    rng = np.random.default_rng(1)
    for _ in range(12):
        ring_draw_stations(rng, length, 5e-6 / h, index)
    built = [ring_draw_stations(rng, length, 10e-6 / h, index) for _ in range(draws)]
    ref = ring_score(length, 2, *perfect_ring(index), q0, False)
    runs = [
        (ring_score(length, 2, z, n, rel, q0, False), ring_score(length, 2, z, n, rel, q0, True)) for z, n, rel in built
    ]
    split = float(np.median([2 * q0 * ring_harmonic(z, 2) for z, _, _ in built]))
    return {
        "Q0": q0,
        "index": index,
        "sigma_pct": 100 * 10e-6 / h,
        "draws": float(draws),
        "split_median": split,
        "allowed_pct": 100 * 10e-6 / h / split,
        "allowed_um": 10e-6 / split * 2 / 1e-6,
        "hybrid_depth_median": float(np.median([fixed.score for fixed, _ in runs])) / ref.score,
        "resolved_depth_median": float(np.median([comp.score for _, comp in runs])) / ref.score,
        "hybrid_purity_median": float(np.median([fixed.purity for fixed, _ in runs])),
        "resolved_purity_median": float(np.median([comp.purity for _, comp in runs])),
    }


# ---- line matrices: the drive study's ring ------------------------------------------------------------------------
@dataclass(frozen=True)
class Ring:
    """A closed guide as a ring of transmission line (section 5.3): length (m), index, Q, impedance z0 (ohm) and
    order m; joints as short sections (centre, length, impedance); feeds at positions, each a Norton current with a
    source admittance, a series capacitor into 50 ohm sized for critical coupling of the set (compensated: its
    susceptance tuned out, a pure conductance). V(s_i) = sum_k G_ik w_k over `samples` positions."""

    length: float
    index: float
    q: float
    z0: float
    m: int
    joints: tuple[tuple[float, float, float], ...]
    feeds: tuple[float, ...]
    samples: int = 240
    compensated: bool = False

    @property
    def shunt(self) -> float:
        """R = Z0 Q / (2 pi m) (RNG-2)."""
        return self.z0 * self.q / (2 * math.pi * self.m)

    @property
    def conductance(self) -> float:
        """Each feed's source conductance for critical coupling of the set."""
        return 1 / (len(self.feeds) * self.shunt)

    @property
    def positions(self) -> Arr:
        out: Arr = (np.arange(self.samples) + 0.5) * self.length / self.samples
        return out

    def resonance(self) -> float:
        """f_m = m c / (n L), pulled down by the feeds' coupling capacitors when they are not compensated."""
        f0 = self.m * C0 / (self.index * self.length)
        if self.compensated:
            return f0
        x = math.sqrt(self.conductance / (50 * (1 - 50 * self.conductance)))
        return f0 * (1 - len(self.feeds) * x * self.z0 / (4 * math.pi * self.m))

    def operator(self, freqs: Arr) -> CArr:
        """G (F, samples, K): volts at the samples per ampere of each feed's Norton current, at each frequency."""
        nf, k = len(freqs), len(self.feeds)
        beta = 2 * math.pi * freqs * self.index / C0
        gamma = beta / (2 * self.q) + 1j * beta
        if self.compensated:
            yf = np.full(nf, self.conductance, dtype=np.complex128)
        else:
            x = math.sqrt(self.conductance / (50 * (1 - 50 * self.conductance)))
            yf = 1 / (50 + 1 / (1j * x * freqs / freqs[nf // 2]))
        # position, kind (0 joint start, 1 joint end, 2 feed, 3 sample), index, impedance
        events: list[tuple[float, int, int, float]] = []
        for c, ell, zj in self.joints:
            events.append(((c - ell / 2) % self.length, 0, 0, zj))
            events.append(((c + ell / 2) % self.length, 1, 0, 0.0))
        for j, p in enumerate(self.feeds):
            events.append((p % self.length, 2, j, 0.0))
        for i, p in enumerate(self.positions):
            events.append((float(p), 3, i, 0.0))
        events.sort(key=lambda e: e[0])
        t = np.tile(np.eye(2, dtype=np.complex128), (nf, 1, 1))
        b = np.zeros((nf, 2, k), dtype=np.complex128)
        rec_t = np.zeros((self.samples, nf, 2), dtype=np.complex128)
        rec_b = np.zeros((self.samples, nf, k), dtype=np.complex128)
        pos, z = 0.0, self.z0
        inside = sum(1 for c, ell, _ in self.joints if (c - ell / 2) % self.length > (c + ell / 2) % self.length)
        if inside:
            z = next(zj for c, ell, zj in self.joints if (c - ell / 2) % self.length > (c + ell / 2) % self.length)
        for p, kind, idx, zj in [*events, (self.length, 4, 0, 0.0)]:
            if p > pos:
                gl = gamma * (p - pos)
                ch, sh = np.cosh(gl), np.sinh(gl)
                mat = np.empty((nf, 2, 2), dtype=np.complex128)
                mat[:, 0, 0], mat[:, 0, 1], mat[:, 1, 0], mat[:, 1, 1] = ch, -z * sh, -sh / z, ch
                t, b = mat @ t, mat @ b
                pos = p
            if kind == 0:
                z = zj
            elif kind == 1:
                z = self.z0
            elif kind == 2:
                shunt = np.tile(np.eye(2, dtype=np.complex128), (nf, 1, 1))
                shunt[:, 1, 0] = -yf
                t, b = shunt @ t, shunt @ b
                b[:, 1, idx] += 1.0
            elif kind == 3:
                rec_t[idx], rec_b[idx] = t[:, 0, :], b[:, 0, :]
        x0 = np.linalg.solve(np.eye(2)[None] - t, b)  # the state at s = 0 per feed current
        g: CArr = np.einsum("ifj,fjk->fik", rec_t, x0) + np.transpose(rec_b, (1, 0, 2))
        return g


def ripple(v: CArr) -> float:
    """(max |V| - min |V|) / (max |V| + min |V|): zero for a travelling wave (RNG-1), one for a standing wave."""
    a = np.abs(v)
    return float((a.max() - a.min()) / (a.max() + a.min()))


def solve_weights(g: CArr, target: CArr, mu: float, iters: int = 8) -> CArr:
    """Regularised least squares to the travelling-wave target, then the target's phase relaxed to the field's (a
    uniform |V| with any phase), repeated."""
    u, s, vh = np.linalg.svd(g, full_matrices=False)
    filt = s / (s * s + mu * float(s[0]) ** 2)
    w: CArr = vh.conj().T @ (filt * (u.conj().T @ target))
    for _ in range(iters):
        w = vh.conj().T @ (filt * (u.conj().T @ np.exp(1j * np.angle(g @ w))))
    return w


@dataclass(frozen=True, eq=False)
class Drive:
    """A compiled drive of the line-matrix ring: its ripple, the share of the available power the ring takes, its
    frequency (Hz), and the feeds' weights and the sample volts, both over the volts' rms."""

    ripple: float
    efficiency: float
    f: float
    w: CArr
    v: CArr


def best_drive(ring: Ring, fixed: CArr | None = None, span: float = 0.03, count: int = 41, least: float = 0.8) -> Drive:
    """The least ripple over a frequency scan about the resonance and over the regularisations, among drives that
    deliver at least `least` of the available power to the ring (when none does, among those within 0.9 of the most
    efficient), or of the fixed weights alone."""
    freqs = ring.resonance() * (1 + np.linspace(-span, span, count))
    gs = ring.operator(freqs)
    target = np.exp(-2j * math.pi * ring.m * ring.positions / ring.length)
    cands: list[Drive] = []
    for i in range(count):
        ws = [fixed] if fixed is not None else [solve_weights(gs[i], target, mu) for mu in MUS]
        for w in ws:
            v = gs[i] @ w
            rms = float(np.sqrt(np.mean(np.abs(v) ** 2)))
            eta = (rms**2 / (2 * ring.shunt)) / (float(np.sum(np.abs(w) ** 2)) / (8 * ring.conductance))
            cands.append(Drive(ripple(v), eta, float(freqs[i]), w / rms, v / rms))
    ok = [c for c in cands if c.efficiency >= least]
    if not ok:
        top = max(c.efficiency for c in cands)
        ok = [c for c in cands if c.efficiency >= 0.9 * top]
    return min(ok, key=lambda c: c.ripple)


def feed_layout(k: int, length: float, m: int, joints: int, at_joints: bool) -> tuple[float, ...]:
    """K feeds: two a quarter guided wave apart; K not dividing 2m equally spaced; K dividing 2m equally spaced with
    alternate feeds moved a quarter wave (so they reach -m); or at the joints."""
    lam = length / m
    if at_joints:
        return tuple((j + 0.5) * length / joints for j in range(k))
    if k == 1:
        return (0.03 * length,)
    if k == 2:
        return (0.03 * length, 0.03 * length + lam / 4)
    if (2 * m) % k:
        return tuple(0.03 * length + j * length / k for j in range(k))
    return tuple(0.03 * length + j * length / k + (lam / 4 if j % 2 else 0.0) for j in range(k))


@dataclass(frozen=True)
class RingCase:
    """A ring of the drive study: band, index, Q, impedance, order, joint count and length."""

    label: str
    f: float
    index: float
    q: float
    z0: float
    m: int
    joints: int
    joint_length: float = 3e-3
    joint_z: float = 1.4
    compensated: bool = False

    @property
    def length(self) -> float:
        return self.m * C0 / (self.index * self.f)

    def ring(
        self,
        k: int,
        rng: np.random.Generator | None,
        at_joints: bool = False,
        tol: float = 0.10,
        jitter: float = 0.5e-3,
    ) -> Ring:
        """The ring with K feeds, its joints identical or drawn: impedance +-tol of Z0 and position +-jitter (sd)."""
        joints = tuple(
            (
                (j + 0.5) * self.length / self.joints + (0.0 if rng is None else float(rng.normal(0, jitter))),
                self.joint_length,
                self.z0 * (self.joint_z + (0.0 if rng is None else float(rng.normal(0, tol)))),
            )
            for j in range(self.joints)
        )
        return Ring(
            self.length,
            self.index,
            self.q,
            self.z0,
            self.m,
            joints,
            feed_layout(k, self.length, self.m, self.joints, at_joints),
            compensated=self.compensated,
        )


DRIVE_2450 = (  # the drive lens's 2450 MHz rings
    RingCase("2450 MHz, m = 2", 2450e6, 1.20, 128, 66, 2, 6),
    RingCase("must fail: m = 3, three joints, 2450 MHz", 2450e6, 1.20, 128, 66, 3, 3),
)


def drive_cases(q: float, z0: float) -> tuple[RingCase, ...]:
    """The drive study's rings: section A's at 915 MHz with its Q and Z0 from the line study (Wheeler's rule on 35
    um copper: Q 250.8, Z0 73.5 ohm; the index as drawn), the same with its feeds' susceptance compensated, and the
    drive lens's rings at 2450 MHz."""
    return (
        RingCase("915 MHz, m = 1", F915, A_INDEX, q, z0, 1, 6),
        RingCase("the same, feed susceptance compensated", F915, A_INDEX, q, z0, 1, 6, compensated=True),
        *DRIVE_2450,
    )
