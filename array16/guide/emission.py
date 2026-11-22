"""Secondary emission by Furman and Pivi's model, and a 2-D tracker of an electron population in a solved section's RF
field: multipactor in the boards' own field (docs/KNOTTED_GUIDE.md sections 5.5, 8.1, 8.2 DR-8, 10 E2-13, 11 MPC-H1).

array16.guide.multipactor reads a gap as parallel plates: one field line's length, a uniform field along it, copper's
yield at normal incidence and secondaries leaving along the normal. Here the field is the section's own, the potential
per volt that array16.guide.sections.solve returns, at the design's volts and band; electrons fly in it in two
dimensions, strike the copper's faces at any angle, and free secondaries by Furman and Pivi's three components, each
with its own yield and spectrum.

  FurmanPivi     the model (M. A. Furman and M. T. F. Pivi, "Probabilistic model for the simulation of secondary
                 electron emission", Phys. Rev. ST Accel. Beams 5, 124404 (2002)): backscattered electrons (eqs. 25
                 and 47a; a half Gaussian under the impact's energy, eq. 26), rediffused ones (eqs. 28 and 47b; E^q
                 under it, eq. 29) and true secondaries (eqs. 31, 32 and 48; n of them up to M = 10, binomial per
                 penetrating electron, eqs. 43 and 46, each E^(p_n - 1) exp(-E / eps_n), eq. 33, their sum held under
                 the impact's energy, eq. 52), all leaving by cos^alpha of the angle from the normal; the defaults
                 copper's, their Tables I and II
  FP_COPPER      bare copper as they fitted it (N. Hilleret's CERN sample, chemically cleaned, not baked in situ):
                 2.09 at 271 eV at normal incidence (their 2.1 at 271 eV), one at 27.3 eV, 0.496 at zero energy
  FP_CONDITIONED conditioned copper: the true secondaries alone scaled to a peak of 1.2, where V. Baglin et al.
                 (Proc. EPAC 2000, 217, section 4.3) found unbaked copper's maximum settling under an electron dose
                 above 10^-3 C/mm^2; one at 90.4 eV
  Normal, CHART  the level-1 chart's surface: multipactor's yield at the impact's energy (ECSS-E-ST-20-01C Table 9-1
                 copper), every secondary along the normal with Anza et al.'s exponential normal energy
  Chamber        where electrons fly and what they strike: the field per volt as a SectionGrid, the copper's boxes and
                 the insulators' (laminate, a board's body), and the region past which an electron has escaped;
                 Chamber.of(solved) for a solved section, Chamber.plates and plates(fd) for parallel plates
  impacts        electrons flown to their first impact: when, where, how fast and at what angle, on which face
  track, growth  the population in time: seeds leave the copper's open faces at random times of the first period,
                 flown by kick-drift-kick steps; each strike multiplies an electron's weight by the yield and sends it
                 off as one emitted electron of the component drawn (an unbiased estimator: the mean emission per
                 impact, the spectrum of one emitted electron), splitting and Russian roulette at the strikes holding
                 the number; each period's growth the weight's, and the growth in nepers per RF period over seeds
  generations    the population generation by generation: the effective yield, the spectral radius of the operator
                 that carries one generation's emissions to the next (multipactor.effective_yield in two dimensions)
  read, Reading  a case read both ways on a surface, with the laminate absorbing or emitting as the copper
  plate_threshold, calibrate
                 the parallel-plate limit: the volts where the tracked growth crosses zero at f d, and the gain that
                 would put it on copper's lowest chart point (27.4 V at 0.70 GHz mm, multipactor.CALIBRATION)
  Case, cases, r1_section, facing_case
                 the protocol's readings: T915-S's section B at 50.3 V, section B with its split widened to h and
                 refitted for 1 eV at q 0.5 (h 0.718 mm, 39.6 V: the split as a lever), R1's rails at 50 V; and facing
                 strips fitted to a depth D at q 0.5, the gap scaled as the root of D (E2-13's family)

**The model's limits.** The field is the section's quasi-static field, uniform along the path: no axial field, and
neither a travelling wave's magnetic field nor its Doppler shift, each v_e / v_phase of what is kept (under 1 % for
strikes under 100 eV on T915-S). Electrons are non-relativistic and carry no space charge: the growth is the
population's small-signal rate, before its own field limits it. The laminate charges nothing: an electron that strikes
it is absorbed (laminate None), or the laminate emits as the copper does; the two bracket an uncharged dielectric, and
E2-13 counts where they agree. Copper's parameters come from one sample (bare) and one scaled case (conditioned);
conditioning also raises the reflection at low energy, which R. Cimino et al. measured near one at zero energy on
scrubbed copper (Phys. Rev. Lett. 93, 014801 (2004)) and the conditioned case keeps at 0.496. Furman and Pivi set
alpha = 1 for the backscattered and rediffused electrons for expediency, noting that their angular spectra differ, and
took the angular dependence from aluminium. A flight that strikes nothing within the horizon (48 RF periods, as
multipactor's) is dropped as held, one that leaves the region or the grid has escaped, and the region is a box: for
facing strips 2h past the copper's outer edges and its back faces, for R1's rails 3.5 mm either way and 1.7 mm up.

**Calibration (DR-0), in the parallel-plate limit** (Claim). On the level-1 chart's own yield (Normal, multipactor's
GAIN) the tracker is multipactor's one-dimensional model in two dimensions: its growth per RF period meets
multipactor.growth_rate within 0.003 Np at 1.15 and 0.85 times copper's lowest point and at 0.9 and 1.1 times the
held-out point, its effective yield multipactor.effective_yield within 0.6 %, and its threshold (plate_threshold, to
the bisection's 0.7 %, three seeds, 500 electrons) lies at 27.3 V at 0.70 GHz mm (-0.3 % of ECSS's 27.4 V), 91.6 V at
the held-out 1.68 GHz mm (-0.3 %) and within -1.7 to +1.8 % of the twelve blind points. Furman and Pivi's bare copper as
published (gain one) multipacts later than the chart on its left branch and sooner on its right: +32.9 % at 0.66 GHz
mm, +19.5 % at 0.70 (32.7 V), +4.6 % at 0.87, +0.3 % at 1.08, -3.1 % at 1.35 and -8.3 % at 1.68 (84.2 V), and not at
all from 0.49 to 0.59 GHz mm within three times ECSS's volts: its crossover, 27.3 eV against the chart's 19 eV, moves
the chart's nose to higher f d. One gain set at 0.70 GHz mm (0.837) would leave the held-out point 23 % low, outside
DR-0's 10 %, so its copper is read as published, with the chart's yield beside it. Conditioned copper multipacts at
1.6 to 2.1 times ECSS's volts from 1.51 to 2.62 GHz mm, and not at 1.20 within three times them. At T915-S's 1.674 GHz
mm, 50.3 V keeps 4.5 dB on bare copper as published (5.2 dB on the chart) and 9.4 dB conditioned.

**The protocol's readings** (Claims; read on the default grid, 1000 electrons, 30 RF periods of which 10 settle and 10
generations of which 3 settle, 200 steps a period, seeds 0 to 3: the growth in Np per RF period and the effective
yield, each the mean over the seeds with the standard deviation of the growth). Every reading decays, in both laminate
models and on every yield:

  case                          yield          laminate absorbs        laminate emits as the copper
  T915-S section B, 50.3 V      bare           -0.165 +- 0.005, 0.798  -0.142 +- 0.004, 0.844
                                conditioned    -0.176 +- 0.013, 0.679  -0.192 +- 0.012, 0.717
                                the chart's    -0.094 +- 0.004, 0.866  -0.060 +- 0.006, 0.896
  split h refitted, 39.6 V      bare           -0.180 +- 0.043, 0.713  -0.175 +- 0.006, 0.833
                                conditioned    -0.186 +- 0.051, 0.626  -0.218 +- 0.021, 0.713
                                the chart's    -0.136 +- 0.002, 0.827  -0.134 +- 0.002, 0.845
  R1's rails, 50 V              bare           dies; -1.06, 0.644      dies; -1.12, 0.660
                                conditioned    dies; -1.09, 0.637      dies; -1.09, 0.644
                                the chart's    dies; -1.13, 0.735      dies; -1.11, 0.736

The effective yields' spreads over the seeds are 0.002 to 0.017. On T915-S's section the strikes carry 16 eV on
average, under bare copper's crossover, 9 to 17 % of them on the polyimide of the split and the margins, and 4 to 7 %
of the population leaves the gap each period; both 1 eV wells decay, as MPC-H1 states for 1 eV (its falsifier, E2-13,
is the family facing_case builds). Each R1 population is gone within 4 to 9 RF periods (Run.lived), and its effective
yield, 0.64 to 0.66 on bare copper and 0.74 on the chart's, is under one: no generation frees as many electrons as
the one before. How fast it falls while it lives, 1.06 to 1.13 Np a period (Run.decay), is how fast its electrons
leave the region: in a region half as large again it lives 16 to 23 periods and falls 0.58 to 0.69 Np a period, its
effective yield 0.640 against 0.648 in the default region. Of the seeds, the secondaries of a 50 eV strike, 78 % leave
the rails' field before they strike again, those that strike carry 4.5 eV at the median, and the lineage settles on
electrons that fall back onto the face they left at 0.5 to 3.5 eV. **Hypothesis (read from these numbers):** MPC-3's
parallel-plate bound, -5.2 dB at 50 V, is conservative on R1's rails: in their own field the electrons leave before
they multiply; DR-8's measured yield and a 3-D field test it. Settings: 100, 200 and 400 steps a period give -0.175,
-0.162 and -0.162 on T915-S (bare, absorbing, two seeds), the effective yield 0.80 at each; 3000 electrons give -0.172,
the coarse grid -0.166.
SI units, except energies in eV and f d in GHz mm.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, replace
from functools import cached_property
from typing import NamedTuple, Protocol

import numpy as np

from array16.guide import multipactor as mp
from array16.guide.sections import (
    FR4,
    RIGID_COPPER,
    Grid,
    Rails,
    Section,
    Solved,
    fit_well,
    section_b,
    solve,
)
from array16.particles import ME, Arr, BoolArr, E, IntArr, SectionGrid

__all__ = [
    "BACKSCATTERED",
    "CHART",
    "EPS",
    "FP_CONDITIONED",
    "FP_COPPER",
    "HORIZON",
    "MULTIPLICITY",
    "QM",
    "REDIFFUSED",
    "ROUNDS",
    "SEED_ENERGY",
    "SPLIT",
    "TRUE",
    "Box",
    "Case",
    "Chamber",
    "Emitted",
    "FurmanPivi",
    "Growth",
    "Impacts",
    "Lineage",
    "Normal",
    "Reading",
    "Run",
    "Strike",
    "Surface",
    "calibrate",
    "cases",
    "facing_case",
    "generations",
    "growth",
    "impacts",
    "plate_threshold",
    "plates",
    "r1_section",
    "read",
    "track",
]

QM = E / ME  # |e/m| (C/kg)
MULTIPLICITY = 10  # M: the most electrons one impact frees ("typically we set M = 10", Furman and Pivi's note 3)
BACKSCATTERED, REDIFFUSED, TRUE = 0, 1, 2  # the components, as Emitted.component names them
EPS = 1e-9  # m: an emitted electron starts this far off the face it leaves
SEED_ENERGY = 50.0  # eV: the seeds leave as the secondaries of an impact at this energy
ROUNDS = 32  # strikes one electron may make within a step; past them it is left at the face, its weight dropped
HORIZON = mp.HORIZON / (2 * math.pi)  # RF periods: a flight that strikes nothing within 48 is dropped, as multipactor's
_DRAWS = 200  # rounds of a rejection sampler; each keeps 0.15 or more of what it draws


# ---- Furman and Pivi -----------------------------------------------------------------------------------------------
class Emitted(NamedTuple):
    """What one impact each sends off: the yield (the electron's weight is multiplied by it), and one emitted
    electron's energy (eV), the cosine of its angle from the surface's normal, its azimuth about the normal (rad) and
    the component it was drawn from (BACKSCATTERED, REDIFFUSED, TRUE)."""

    yields: Arr
    energy: Arr
    cos: Arr
    azimuth: Arr
    component: IntArr


class Surface(Protocol):
    """A surface's secondary emission: for impacts at energies (eV) and cosines of incidence from the normal, the yield
    and one emitted electron each."""

    def emit(self, rng: np.random.Generator, energy: Arr, cos: Arr) -> Emitted: ...


def _deficit(rng: np.random.Generator, top: Arr, sigma: float) -> Arr:
    """|N(0, sigma)| held under `top` (each >= 0): how far under the impact energy a backscattered electron leaves
    (Furman and Pivi's eq. 26, a half Gaussian ending at the impact energy). Half-normal draws kept under the top where
    it is at least sigma / 2 (0.38 of them or more kept), uniform ones on [0, top] kept with exp(-d^2 / 2 sigma^2)
    (0.88 or more) where it is less."""
    out = np.zeros(top.shape)
    todo: IntArr = np.flatnonzero(top > 0)
    for _ in range(_DRAWS):
        if not todo.size:
            break
        tt = top[todo]
        near = tt < 0.5 * sigma
        d = np.where(near, tt * rng.uniform(size=todo.size), sigma * np.abs(rng.standard_normal(todo.size)))
        keep = np.where(near, rng.uniform(size=todo.size) <= np.exp(-0.5 * (d / sigma) ** 2), d <= tt)
        out[todo[keep]] = d[keep]
        todo = todo[~keep]
    out[todo] = 0.0
    return out


def _truncated_gamma(rng: np.random.Generator, shape: Arr, limit: Arr) -> Arr:
    """Gamma(shape, 1) held under `limit` (Furman and Pivi's eq. 52: n true secondaries' summed energy under the
    impact's, in units of eps_n). Where the limit is at least shape - 1, the gamma's draws kept under it (0.2 or more of
    them); under that, x = limit (1 - y) with y from Beta(1, shape - limit) kept with ((1 - y) e^y)^limit, which is
    exact, the target being (1 - y)^(shape - 1) e^(limit y) in y (0.15 or more kept)."""
    out = np.zeros(shape.shape)
    todo: IntArr = np.flatnonzero(limit > 0)
    for _ in range(_DRAWS):
        if not todo.size:
            break
        a, cut = shape[todo], limit[todo]
        wide = cut >= a - 1
        x = np.empty(todo.size)
        keep = np.zeros(todo.size, dtype=np.bool_)
        if wide.any():
            g = rng.gamma(a[wide])
            x[wide] = g
            keep[wide] = g <= cut[wide]
        narrow = ~wide
        if narrow.any():
            b, c = a[narrow] - cut[narrow], cut[narrow]
            y = 1 - rng.uniform(size=b.size) ** (1 / b)
            x[narrow] = c * (1 - y)
            keep[narrow] = rng.uniform(size=b.size) <= ((1 - y) * np.exp(y)) ** c
        out[todo[keep]] = x[keep]
        todo = todo[~keep]
    out[todo] = np.minimum(shape[todo] - 1, limit[todo])  # the mode: not reached in practice
    return out


@dataclass(frozen=True)
class FurmanPivi:
    """Furman and Pivi's secondary emission (Phys. Rev. ST Accel. Beams 5, 124404 (2002)), the yield per impact at
    energy E (eV) and angle theta from the normal split three ways:

        backscattered  delta_e = [P1e(inf) + (P1e^ - P1e(inf)) exp(-(|E - Ee^| / W)^p / p)] [1 + e1 (1 - cos^e2)]
        rediffused     delta_r = P1r(inf) [1 - exp(-(E / Er)^r)] [1 + r1 (1 - cos^r2)]
        true           delta_ts = d^ D(E / E^),  D(x) = s x / (s - 1 + x^s),
                       d^ = dts^ [1 + t1 (1 - cos^t2)],  E^ = Ets^ [1 + t3 (1 - cos^t4)]

    (their eqs. 25, 28, 31, 32, 47 and 48). An impact frees one backscattered electron with probability delta_e, one
    rediffused with delta_r, and otherwise n true secondaries, n binomial over M = 10 with mean
    delta_ts / (1 - delta_e - delta_r) (eqs. 43 and 46). A backscattered electron leaves with a half Gaussian of width
    sigma_e under E (eq. 26), a rediffused one with (q + 1) E'^q / E^(q + 1) (eq. 29); n true secondaries with energies
    each E'^(p_n - 1) exp(-E' / eps_n) (eq. 33), their sum held under E (eq. 52, through the normalised incomplete
    gamma). Every emitted electron leaves by cos^alpha of its angle from the normal, its azimuth uniform.

    The defaults are copper's (their Table I and Table II, fitted to N. Hilleret's CERN sample, chemically cleaned but
    not baked in situ, at normal incidence; the angular parameters from aluminium data): 2.09 at 271 eV at normal
    incidence (their 2.1 at 271 eV), one first at 27.3 eV, 0.496 at zero energy, where only backscattering is left."""

    p1e_inf: float = 0.02
    p1e_peak: float = 0.496
    e_e: float = 0.0  # eV
    w: float = 60.86  # eV
    p: float = 1.0
    sigma_e: float = 2.0  # eV
    e1: float = 0.26
    e2: float = 2.0
    p1r_inf: float = 0.2
    e_r: float = 0.041  # eV
    r: float = 0.104
    q: float = 0.5
    r1: float = 0.26
    r2: float = 2.0
    delta_ts: float = 1.8848
    e_ts: float = 276.8  # eV
    s: float = 1.54
    t1: float = 0.66
    t2: float = 0.8
    t3: float = 0.7
    t4: float = 1.0
    p_n: tuple[float, ...] = (2.5, 3.3, 2.5, 2.5, 2.8, 1.3, 1.5, 1.5, 1.5, 1.5)
    eps_n: tuple[float, ...] = (1.5, 1.75, 1.0, 3.75, 8.5, 11.5, 2.5, 3.0, 2.5, 3.0)  # eV
    alpha: float = 1.0
    name: str = "copper, bare (Furman and Pivi)"

    def __post_init__(self) -> None:
        if len(self.p_n) != MULTIPLICITY or len(self.eps_n) != MULTIPLICITY:
            raise ValueError(
                f"p_n and eps_n take {MULTIPLICITY} values each, not {len(self.p_n)} and {len(self.eps_n)}"
            )
        positive = (self.w, self.p, self.sigma_e, self.e_r, self.r, self.e_ts, self.s, *self.p_n, *self.eps_n)
        if not all(x > 0 for x in positive) or self.s <= 1 or self.q <= -1 or self.alpha < 0 or self.delta_ts < 0:
            raise ValueError(
                "a Furman and Pivi surface takes positive widths, scales, powers and p_n, eps_n, s > 1, q > -1, "
                f"alpha >= 0 and dts^ >= 0: {self}"
            )
        grazing = max(self.p1e_peak, self.p1e_inf) * (1 + self.e1) + self.p1r_inf * (1 + self.r1)
        if min(self.p1e_inf, self.p1e_peak, self.p1r_inf) < 0 or not grazing < 1:
            raise ValueError(
                f"backscattering and rediffusion must leave some electrons to penetrate at every angle: "
                f"P1e and P1r >= 0 with (P1e^ (1 + e1) + P1r(inf) (1 + r1)) < 1, not {grazing:.3f}"
            )

    # ---- the yields --------------------------------------------------------------------------------------------
    def elastic(self, energy: Arr | float, cos: Arr | float = 1.0) -> Arr:
        """delta_e: backscattered electrons per impact (eqs. 25 and 47a)."""
        e, c = np.asarray(energy, dtype=np.float64), np.clip(np.asarray(cos, dtype=np.float64), 0.0, 1.0)
        base = self.p1e_inf + (self.p1e_peak - self.p1e_inf) * np.exp(
            -((np.abs(e - self.e_e) / self.w) ** self.p) / self.p
        )
        out: Arr = base * (1 + self.e1 * (1 - c**self.e2))
        return out

    def rediffused(self, energy: Arr | float, cos: Arr | float = 1.0) -> Arr:
        """delta_r: rediffused electrons per impact (eqs. 28 and 47b)."""
        e, c = np.asarray(energy, dtype=np.float64), np.clip(np.asarray(cos, dtype=np.float64), 0.0, 1.0)
        out: Arr = (
            self.p1r_inf * -np.expm1(-((np.maximum(e, 0.0) / self.e_r) ** self.r)) * (1 + self.r1 * (1 - c**self.r2))
        )
        return out

    def true(self, energy: Arr | float, cos: Arr | float = 1.0) -> Arr:
        """delta_ts: true secondaries per impact (eqs. 31, 32 and 48)."""
        e, c = np.asarray(energy, dtype=np.float64), np.clip(np.asarray(cos, dtype=np.float64), 0.0, 1.0)
        peak = self.delta_ts * (1 + self.t1 * (1 - c**self.t2))
        x = np.maximum(e, 0.0) / (self.e_ts * (1 + self.t3 * (1 - c**self.t4)))
        out: Arr = peak * self.s * x / (self.s - 1 + x**self.s)
        return out

    def __call__(self, energy: Arr | float, cos: Arr | float = 1.0) -> Arr:
        """The yield: delta_e + delta_r + delta_ts. At normal incidence it is a multipactor.Spectrum."""
        out: Arr = self.elastic(energy, cos) + self.rediffused(energy, cos) + self.true(energy, cos)
        return out

    @cached_property
    def peak(self) -> tuple[float, float]:
        """The yield's maximum at normal incidence: (energy in eV, yield), by a scan and golden sections."""
        grid = np.geomspace(1.0, 2e4, 2001)
        i = int(np.argmax(self(grid)))
        lo, hi = float(grid[max(i - 1, 0)]), float(grid[min(i + 1, grid.size - 1)])
        golden = (math.sqrt(5) - 1) / 2
        for _ in range(80):
            a, b = hi - golden * (hi - lo), lo + golden * (hi - lo)
            lo, hi = (lo, b) if float(self(a)) > float(self(b)) else (a, hi)
        at = 0.5 * (lo + hi)
        return at, float(self(at))

    @cached_property
    def crossover(self) -> float:
        """The first energy (eV) at which the yield at normal incidence reaches one; nan if it never does."""
        at, top = self.peak
        if top < 1 or float(self(0.0)) >= 1:
            return math.nan
        lo, hi = 0.0, at
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if float(self(mid)) < 1 else (lo, mid)
        return 0.5 * (lo + hi)

    def conditioned(self, peak: float, name: str | None = None) -> FurmanPivi:
        """The same surface with its true secondaries scaled (dts^ alone) so that the yield's maximum at normal
        incidence is `peak`: how electron-cloud work models a conditioned surface."""
        floor = max(float(self.elastic(e) + self.rediffused(e)) for e in np.geomspace(1.0, 2e4, 401))
        if not peak > floor:
            raise ValueError(
                f"the backscattered and rediffused electrons alone reach {floor:.3f}; a peak of {peak} is under it"
            )
        lo, hi = 0.0, 1.0
        while replace(self, delta_ts=self.delta_ts * hi).peak[1] < peak:
            hi *= 2
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if replace(self, delta_ts=self.delta_ts * mid).peak[1] < peak else (lo, mid)
        label = name if name is not None else f"{self.name}, conditioned to {peak:g}"
        return replace(self, delta_ts=self.delta_ts * 0.5 * (lo + hi), name=label)

    # ---- the emission ------------------------------------------------------------------------------------------
    def multiplicity(self, energy: Arr | float, cos: Arr | float = 1.0) -> Arr:
        """P_n, the probability that an impact frees n electrons, n = 0 to M ((..., M + 1)): eq. 43 with the binomial of
        eq. 46, P'_n = C(M, n) p^n (1 - p)^(M - n), p = delta_ts / ((1 - delta_e - delta_r) M):
        P_0 = (1 - de - dr) P'_0, P_1 = (1 - de - dr) P'_1 + de + dr, P_n = (1 - de - dr) P'_n for n >= 2. Its sum is
        one and its mean the yield."""
        de, dr, dts = self.elastic(energy, cos), self.rediffused(energy, cos), self.true(energy, cos)
        rest = 1 - de - dr
        p = np.clip(dts / (rest * MULTIPLICITY), 0.0, 1.0)[..., None]
        n = np.arange(MULTIPLICITY + 1)
        coefficients = np.array([math.comb(MULTIPLICITY, k) for k in n], dtype=np.float64)
        out: Arr = rest[..., None] * coefficients * p**n * (1 - p) ** (MULTIPLICITY - n)
        out[..., 1] += de + dr
        return out

    def true_energy(self, rng: np.random.Generator, energy: Arr, cos: Arr) -> Arr:
        """One true secondary's energy (eV) for each impact: from an event of n true secondaries, n drawn as one emitted
        electron sees it (n - 1 binomial over M - 1, the size-biased multiplicity of eq. 46), the event's summed energy
        eps_n Gamma(n p_n) held under the impact's (eq. 52), and this electron's share of it Beta(p_n, (n - 1) p_n), as
        the joint spectrum of n independent gammas held to their sum makes it."""
        e, c = np.asarray(energy, dtype=np.float64), np.asarray(cos, dtype=np.float64)
        rest = 1 - self.elastic(e, c) - self.rediffused(e, c)
        p = np.clip(self.true(e, c) / (rest * MULTIPLICITY), 0.0, 1.0)
        n = 1 + rng.binomial(MULTIPLICITY - 1, p)
        shape, eps = np.asarray(self.p_n)[n - 1], np.asarray(self.eps_n)[n - 1]
        total = eps * _truncated_gamma(rng, n * shape, np.maximum(e, 0.0) / eps)
        share = np.where(n > 1, rng.beta(shape, np.where(n > 1, (n - 1) * shape, 1.0)), 1.0)
        out: Arr = total * share
        return out

    def emit(self, rng: np.random.Generator, energy: Arr, cos: Arr) -> Emitted:
        """One emitted electron for each impact, weighted by the yield: its component drawn in proportion to delta_e,
        delta_r and delta_ts, its energy from that component's spectrum, its angle from cos^alpha."""
        e = np.asarray(energy, dtype=np.float64).ravel()
        c = np.clip(np.asarray(cos, dtype=np.float64).ravel(), 0.0, 1.0)
        de, dr, dts = self.elastic(e, c), self.rediffused(e, c), self.true(e, c)
        total = de + dr + dts
        pick = rng.uniform(size=e.size) * total
        component: IntArr = np.where(pick < de, BACKSCATTERED, np.where(pick < de + dr, REDIFFUSED, TRUE)).astype(
            np.intp
        )
        out = np.zeros(e.size)
        back, red, tru = component == BACKSCATTERED, component == REDIFFUSED, component == TRUE
        out[back] = np.maximum(e[back] - _deficit(rng, np.maximum(e[back], 0.0), self.sigma_e), 0.0)
        out[red] = np.maximum(e[red], 0.0) * rng.uniform(size=int(red.sum())) ** (1 / (self.q + 1))
        out[tru] = self.true_energy(rng, e[tru], c[tru])
        angle = rng.uniform(size=e.size) ** (1 / (self.alpha + 1))
        return Emitted(total, out, angle, rng.uniform(0.0, 2 * math.pi, e.size), component)


FP_COPPER = FurmanPivi()  # Table I and Table II as published: copper chemically cleaned, not baked in situ
FP_CONDITIONED = FP_COPPER.conditioned(1.2, "copper, conditioned to a peak of 1.2 (Baglin et al.)")


@dataclass(frozen=True)
class Normal:
    """The surface the level-1 chart is computed on: a yield against the impact's energy (multipactor's Vaughan curve on
    ECSS-E-ST-20-01C Table 9-1 copper by default), every secondary leaving along the normal with a normal energy from
    the emission (Anza et al.'s exponential, mean 2 W_t). Between parallel plates it is multipactor's one-dimensional
    model."""

    sey: mp.Spectrum = mp.COPPER
    emission: mp.Emission = mp.EMISSION
    name: str = "copper, ECSS Table 9-1 (the level-1 chart's)"

    def emit(self, rng: np.random.Generator, energy: Arr, cos: Arr) -> Emitted:
        e = np.asarray(energy, dtype=np.float64).ravel()
        n = e.size
        yields = np.asarray(self.sey(e), dtype=np.float64)
        return Emitted(yields, self.emission.draw(rng, n), np.ones(n), np.zeros(n), np.full(n, TRUE, dtype=np.intp))


CHART = Normal()


# ---- the chamber ---------------------------------------------------------------------------------------------------
class Box(NamedTuple):
    """An axis-aligned rectangle in the section (m): u0 <= u <= u1, v0 <= v <= v1."""

    u0: float
    u1: float
    v0: float
    v1: float


class Strike(NamedTuple):
    """Where segments first meet a box: the fraction of the segment (inf where none), the box (-1 where none), the
    outward normal of the face struck (N, 2), and whether the box is copper."""

    tau: Arr
    box: IntArr
    normal: Arr
    copper: BoolArr


def _slab(s: Arr, d: Arr, a: Arr, b: Arr) -> tuple[Arr, Arr]:
    """The fractions of the way along s + tau d at which a coordinate enters and leaves [a, b]: (-inf, inf) for a
    segment that runs inside the slab without crossing it, (inf, -inf) for one that runs outside it."""
    with np.errstate(divide="ignore", invalid="ignore"):
        t0 = (a - s) / d
        t1 = (b - s) / d
    flat = d == 0
    inside = (s >= a) & (s <= b)
    lo: Arr = np.where(flat, np.where(inside, -np.inf, np.inf), np.minimum(t0, t1))
    hi: Arr = np.where(flat, np.where(inside, np.inf, -np.inf), np.maximum(t0, t1))
    return lo, hi


@dataclass(frozen=True, eq=False)
class Chamber:
    """Where electrons fly and what they strike: the field per volt of drive as a SectionGrid (E = -V(t) grad phi), the
    band f (Hz), the copper's boxes and the insulators' (laminate, a board's body), and the region beyond which an
    electron has escaped (m)."""

    grid: SectionGrid
    f: float
    copper: tuple[Box, ...]
    insulator: tuple[Box, ...]
    region: Box
    name: str = ""

    def __post_init__(self) -> None:
        if not self.f > 0 or not self.copper:
            raise ValueError(f"a chamber takes a band f > 0 and copper to strike, not {self.f} Hz and {self.copper}")
        r = self.region
        if not (r.u1 > r.u0 and r.v1 > r.v0):
            raise ValueError(f"a chamber's region must have extent: {r}")

    @classmethod
    def plates(cls, gap: float, f: float, thickness: float | None = None) -> Chamber:
        """Parallel copper plates `gap` (m) apart, the lower at v <= 0 and the upper at v >= gap, a metre either way
        along u: the field per volt 1 / gap across them, the volts V across the gap, as multipactor counts them."""
        if not (gap > 0 and f > 0):
            raise ValueError(f"plates take a gap > 0 and a band f > 0, not {gap} m and {f} Hz")
        t = 0.1 * gap if thickness is None else thickness
        u = np.array([-1.0, 0.0, 1.0])
        v = np.array([-t, 0.5 * gap, gap + t])
        grid = SectionGrid.sampled(lambda uu, vv: vv / gap, u, v, (0.0, 0.5 * gap))
        copper = (Box(-1.0, 1.0, -t, 0.0), Box(-1.0, 1.0, gap, gap + t))
        return cls(
            grid, f, copper, (), Box(-1.0, 1.0, -t, gap + t), f"plates {gap * 1e3:g} mm apart at {f / 1e9:g} GHz"
        )

    @classmethod
    def of(cls, solved: Solved, region: Box | None = None) -> Chamber:
        """A solved section's chamber: its potential per volt on the solver's grid (Solved.section_grid), its copper's
        boxes and its boards' (the laminate, a board's body) unfolded over the whole section, at its band. The region
        is the layout's flight box (facing strips: 2h beyond the copper's outer edges and the copper's back faces);
        along an axis the layout leaves open, twice (across) or once (up) the farthest copper feature from the centre
        line: for R1's rails 3.5 mm either way and 1.7 mm up."""
        section = solved.section
        lay = section.geometry._layout(section, section.grid)
        h, v0 = lay.h, lay.v0
        mirrored = lay.mirror_y != "none"

        def unfold(boxes: Sequence[object]) -> tuple[Box, ...]:
            out = []
            for b in boxes:
                x0, x1, y0, y1 = (float(getattr(b, k)) for k in ("x0", "x1", "y0", "y1"))
                ys = [(y0, y1), (-y1, -y0)] if mirrored else [(y0, y1)]
                for a0, a1 in ((x0, x1), (-x1, -x0)):
                    for c0, c1 in ys:
                        out.append(Box(a0 * h, a1 * h, (c0 - v0) * h, (c1 - v0) * h))
            return tuple(out)

        reach = max(abs(x) for x in lay.xfeat)
        fx = lay.flight[0] if math.isfinite(lay.flight[0]) else 2 * reach
        fy = lay.flight[1] if math.isfinite(lay.flight[1]) else reach
        ylo = -fy if mirrored else lay.ylo
        box = region if region is not None else Box(-fx * h, fx * h, (ylo - v0) * h, (fy - v0) * h)
        kind = type(section.geometry).__name__
        return cls(
            solved.section_grid(),
            section.f,
            unfold(lay.copper),
            unfold(lay.solid),
            box,
            f"{kind} at {section.f / 1e6:g} MHz, h {solved.distance * 1e3:.3f} mm",
        )

    @cached_property
    def _boxes(self) -> tuple[Arr, BoolArr]:
        boxes = np.array([*self.copper, *self.insulator], dtype=np.float64).reshape(-1, 4)
        copper = np.zeros(len(boxes), dtype=np.bool_)
        copper[: len(self.copper)] = True
        return boxes, copper

    def solid(self, points: Arr) -> BoolArr:
        """Which points (N, 2) lie in or on a box."""
        boxes, _ = self._boxes
        u, v = points[:, 0, None], points[:, 1, None]
        out: BoolArr = ((u >= boxes[:, 0]) & (u <= boxes[:, 1]) & (v >= boxes[:, 2]) & (v <= boxes[:, 3])).any(axis=1)
        return out

    def inside(self, points: Arr) -> BoolArr:
        """Which points (N, 2) lie strictly inside the region."""
        r = self.region
        u, v = points[:, 0], points[:, 1]
        out: BoolArr = (u > r.u0) & (u < r.u1) & (v > r.v0) & (v < r.v1)
        return out

    def hit(self, start: Arr, step: Arr) -> Strike:
        """Where segments from `start` (N, 2) along `step` (N, 2) first meet a box: the slab test against each box, the
        face struck being the one whose slab the segment entered last, its normal pointing back along the segment."""
        boxes, copper = self._boxes
        n = len(start)
        lo_u, hi_u = _slab(start[:, 0, None], step[:, 0, None], boxes[None, :, 0], boxes[None, :, 1])
        lo_v, hi_v = _slab(start[:, 1, None], step[:, 1, None], boxes[None, :, 2], boxes[None, :, 3])
        enter = np.maximum(lo_u, lo_v)
        meets = (enter <= np.minimum(hi_u, hi_v)) & (enter <= 1.0) & (np.minimum(hi_u, hi_v) >= 0.0)
        tau = np.where(meets, np.maximum(enter, 0.0), np.inf)
        which: IntArr = np.argmin(tau, axis=1)
        rows = np.arange(n)
        first = tau[rows, which]
        struck = np.isfinite(first)
        across = lo_u[rows, which] >= lo_v[rows, which]
        normal = np.zeros((n, 2))
        normal[:, 0] = np.where(across, -np.sign(step[:, 0]), 0.0)
        normal[:, 1] = np.where(across, 0.0, -np.sign(step[:, 1]))
        box: IntArr = np.where(struck, which, -1)
        return Strike(first, box, normal, np.where(struck, copper[which], False))

    def faces(self, count: int = 64) -> tuple[Arr, Arr, Arr]:
        """Points on the copper's faces open to the vacuum inside the region, EPS off their faces, with the faces'
        outward normals (N, 2) and the length (m) each point stands for: where seed electrons leave."""
        points, normals, lengths = [], [], []
        r = self.region
        for b in self.copper:
            for axis, at, sign in ((0, b.u0, -1.0), (0, b.u1, 1.0), (1, b.v0, -1.0), (1, b.v1, 1.0)):
                lo, hi = (max(b.v0, r.v0), min(b.v1, r.v1)) if axis == 0 else (max(b.u0, r.u0), min(b.u1, r.u1))
                if not hi > lo:
                    continue
                along = lo + (np.arange(count) + 0.5) * ((hi - lo) / count)
                p = np.empty((count, 2))
                p[:, axis] = at + sign * EPS
                p[:, 1 - axis] = along
                _, _, on = self.grid.gradient_at(p[:, 0], p[:, 1])
                ok = self.inside(p) & ~self.solid(p) & on
                if ok.any():
                    nrm = np.zeros((int(ok.sum()), 2))
                    nrm[:, axis] = sign
                    points.append(p[ok])
                    normals.append(nrm)
                    lengths.append(np.full(int(ok.sum()), (hi - lo) / count))
        if not points:
            raise ValueError(f"no copper face of {self.name or 'the chamber'} lies open inside its region")
        return np.concatenate(points), np.concatenate(normals), np.concatenate(lengths)


def plates(fd: float, gap: float = 1e-3) -> Chamber:
    """Parallel plates at f d (GHz mm): `gap` apart (1 mm) at fd / gap."""
    return Chamber.plates(gap, fd * mp.GHZ_MM / gap)


# ---- flights -------------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class _Drive:
    """The chamber at volts V (read through a gain) and its band: the acceleration a = |e/m| V cos(omega t) grad phi."""

    chamber: Chamber
    scale: float  # |e/m| V / gain
    omega: float

    def accel(self, points: Arr, t: Arr | float) -> tuple[Arr, BoolArr]:
        gu, gv, on = self.chamber.grid.gradient_at(points[:, 0], points[:, 1])
        k = self.scale * np.cos(self.omega * np.asarray(t, dtype=np.float64))
        out = np.empty((len(points), 2))
        out[:, 0] = k * gu
        out[:, 1] = k * gv
        return out, on


def _velocity(energy: Arr, cos: Arr, azimuth: Arr, normal: Arr) -> tuple[Arr, Arr]:
    """Emitted electrons' velocities from their energies (eV), cosines from the faces' normals (N, 2) and azimuths
    about them: in the section (N, 2), along the normal and the face, and along the path (N,)."""
    speed = np.sqrt(2 * QM * np.maximum(energy, 0.0))
    sin = np.sqrt(np.maximum(1 - cos**2, 0.0))
    tangent = np.stack([-normal[:, 1], normal[:, 0]], axis=1)
    inplane = speed[:, None] * (cos[:, None] * normal + (sin * np.cos(azimuth))[:, None] * tangent)
    out: tuple[Arr, Arr] = (inplane, speed * sin * np.sin(azimuth))
    return out


class Impacts(NamedTuple):
    """Electrons flown to their first impact: the time (s; nan where none within the flight), the point (N, 2), the
    velocity there (N, 3: across, toward and along the path), its energy (eV), the cosine of incidence, the outward
    normal of the face struck (N, 2), whether the box struck is copper, which struck at all, and which left the region
    or the grid (the rest of those that struck nothing flew the whole span)."""

    time: Arr
    point: Arr
    velocity: Arr
    energy: Arr
    cos: Arr
    normal: Arr
    copper: BoolArr
    struck: BoolArr
    escaped: BoolArr


def impacts(
    chamber: Chamber,
    volts: float,
    start: Arr,
    velocity: Arr,
    t0: Arr | float = 0.0,
    *,
    gain: float = 1.0,
    steps: int = 200,
    periods: float = 8.0,
) -> Impacts:
    """Electrons from `start` (N, 2) with `velocity` (N, 3) at times t0, flown in the chamber driven at `volts` by
    kick-drift-kick steps of 1/steps of an RF period on a clock common to all, to their first impact within `periods`
    RF periods; an electron that leaves the region or the grid has none. Each drift is a segment tested against the
    boxes; the impact's velocity is the drift's corrected by the acceleration there over the part of the step before
    it, so that the impact's energy is good to second order in the step."""
    x = np.array(start, dtype=np.float64).reshape(-1, 2)
    n = len(x)
    vel = np.array(velocity, dtype=np.float64).reshape(-1, 3)
    born = np.broadcast_to(np.asarray(t0, dtype=np.float64), (n,)).copy()
    drive = _Drive(chamber, QM * volts / gain, 2 * math.pi * chamber.f)
    dt = 1 / (chamber.f * steps)
    span = periods / chamber.f
    time, point = np.full(n, np.nan), np.full((n, 2), np.nan)
    hitv, energy, cosine = np.full((n, 3), np.nan), np.full(n, np.nan), np.full(n, np.nan)
    normal = np.zeros((n, 2))
    copper, struck, escaped = (np.zeros(n, dtype=np.bool_) for _ in range(3))
    state = np.zeros(n, dtype=np.intp)  # 0 waiting, 1 flying, 2 done
    half = vel[:, :2].copy()  # the in-plane velocity, pre-kicked for the step under way
    first = math.floor(float(born.min()) / dt) if n else 0
    for k in range(first, math.ceil((float(born.max()) + span) / dt) + 1 if n else first):
        t, t_end = k * dt, (k + 1) * dt
        fresh = np.flatnonzero((state == 0) & (born < t_end))
        if fresh.size:
            a0, _ = drive.accel(x[fresh], born[fresh])
            half[fresh] = vel[fresh, :2] + a0 * (0.5 * (t_end - born[fresh]))[:, None]
            state[fresh] = 1
        go = np.flatnonzero(state == 1)
        if not go.size:
            if not (state == 0).any():
                break
            continue
        t_start = np.maximum(born[go], t)
        dur = t_end - t_start
        disp = half[go] * dur[:, None]
        hit = chamber.hit(x[go], disp)
        s = hit.box >= 0
        i, tau, d = go[s], hit.tau[s], dur[s]
        if i.size:
            t_hit = t_start[s] + tau * d
            p = x[i] + tau[:, None] * disp[s]
            a, _ = drive.accel(p, t_hit)
            v = half[i] + a * ((tau - 0.5) * d)[:, None]
            vz = vel[i, 2]
            speed2 = (v**2).sum(axis=1) + vz**2
            time[i], point[i] = t_hit, p
            hitv[i, :2], hitv[i, 2] = v, vz
            energy[i] = 0.5 * speed2 / QM
            cosine[i] = np.abs((v * hit.normal[s]).sum(axis=1)) / np.sqrt(np.maximum(speed2, 1e-300))
            normal[i], copper[i], struck[i], state[i] = hit.normal[s], hit.copper[s], True, 2
        rest = go[~s]
        x[rest] = x[rest] + disp[~s]
        a1, on = drive.accel(x[rest], t_end)
        out = ~on | ~chamber.inside(x[rest])
        gone = out | ((t_end - born[rest]) > span)
        state[rest[gone]] = 2  # left the chamber, or flew the whole span: no impact
        escaped[rest[out]] = True
        keep = ~gone
        half[rest[keep]] = half[rest[keep]] + a1[keep] * (0.5 * dur[~s][keep] + 0.5 * dt)[:, None]
        if not (state < 2).any():
            break
    return Impacts(time, point, hitv, energy, cosine, normal, copper, struck, escaped)


class Run(NamedTuple):
    """What track measured, over the RF periods after the first `settle`."""

    growth: float  # Np per RF period: the mean of the log of each period's growth; -inf if the population died
    spread: float  # its standard error over the periods
    factors: Arr  # each period's growth, from the first
    impacts: float  # strikes per period over the electrons at its start
    energy: float  # the mean energy of a strike (eV), weighted by the electrons' weights
    laminate: float  # the share of the strikes' weight on insulators
    escaped: float  # the share of the population's weight leaving the region each period
    held: float  # the share dropped each period for striking nothing within the horizon

    @property
    def lived(self) -> int:
        """The RF periods the population lived through: all of them, unless it died."""
        zero = np.flatnonzero(self.factors <= 0)
        return int(zero[0]) if zero.size else len(self.factors)

    @property
    def decay(self) -> float:
        """The mean log of the periods' growths over those it lived through after the first (Np per RF period): the
        rate a population that dies fell at."""
        lived = self.factors[1 : self.lived]
        return float(np.mean(np.log(lived))) if lived.size else math.nan


class Growth(NamedTuple):
    """The growth over seeds: the mean of their growths (Np per RF period), its spread (their standard deviation), and
    each run."""

    rate: float
    spread: float
    runs: tuple[Run, ...]

    @property
    def grows(self) -> bool:
        """Whether the population grows: multipactor."""
        return self.rate > 0

    @property
    def factor(self) -> float:
        """The population's growth over one RF period: e^rate."""
        return math.exp(self.rate) if self.rate > -700 else 0.0


_SLOTS = ("x", "half", "vz", "weight", "last", "born", "dur", "waiting")


class _Swarm:
    """The electrons a run follows, in slots that splitting adds and each period's compaction removes: the position
    (m), the in-plane velocity pre-kicked for the step under way, the axial velocity, the weight, the time each last
    left a surface, its birth, the duration of the last segment it flew this step, and whether it waits to be born."""

    def __init__(self, x: Arr, velocity: Arr, vz: Arr, born: Arr) -> None:
        n = len(born)
        cap = max(4 * n, 64)
        self.x, self.half = np.zeros((cap, 2)), np.zeros((cap, 2))
        self.vz, self.weight, self.last, self.born, self.dur = (np.zeros(cap) for _ in range(5))
        self.waiting = np.zeros(cap, dtype=np.bool_)
        self.x[:n], self.half[:n], self.vz[:n] = x, velocity, vz
        self.weight[:n], self.last[:n], self.born[:n], self.waiting[:n] = 1.0, born, born, True
        self.n = n

    def add(self, count: int) -> IntArr:
        """`count` new slots after the last, the arrays grown as needed: their indices."""
        need = self.n + count
        if need > len(self.weight):
            cap = max(2 * len(self.weight), need)
            for name in _SLOTS:
                old = getattr(self, name)
                new = np.zeros((cap, *old.shape[1:]), dtype=old.dtype)
                new[: self.n] = old[: self.n]
                setattr(self, name, new)
        out: IntArr = np.arange(self.n, need)
        self.n = need
        return out

    def compact(self) -> None:
        """The slots of weight zero removed, the rest kept in order."""
        keep = np.flatnonzero(self.weight[: self.n] > 0)
        for name in _SLOTS:
            arr = getattr(self, name)
            arr[: keep.size] = arr[keep]
        self.weight[keep.size : self.n] = 0.0
        self.n = keep.size


def _emit(
    rng: np.random.Generator, surface: Surface, laminate: Surface | None, energy: Arr, cos: Arr, copper: BoolArr
) -> tuple[Arr, Arr, Arr, Arr]:
    """Each impact's yield (zero where an insulator absorbs) and one emitted electron's energy, cosine and azimuth."""
    k = energy.size
    y, e, c, a = np.zeros(k), np.zeros(k), np.ones(k), np.zeros(k)
    for where, surf in ((copper, surface), (~copper, laminate)):
        if surf is None or not where.any():
            continue
        out = surf.emit(rng, energy[where], cos[where])
        y[where], e[where], c[where], a[where] = out.yields, out.energy, out.cos, out.azimuth
    return y, e, c, a


SPLIT = 64  # the most electrons one emission is split into


def track(
    chamber: Chamber,
    volts: float,
    surface: Surface = FP_COPPER,
    *,
    laminate: Surface | None = None,
    gain: float = 1.0,
    electrons: int = 1000,
    periods: int = 30,
    settle: int = 10,
    steps: int = 200,
    horizon: float = HORIZON,
    seed: int = 0,
) -> Run:
    """A population of about `electrons` in the chamber driven at `volts` (peak, read as volts / gain), followed for
    `periods` RF periods by kick-drift-kick steps of 1/steps of a period. They leave the copper's open faces
    (Chamber.faces, by length) at random times of the first period as the secondaries of a SEED_ENERGY impact. Each
    strike on copper multiplies an electron's weight by the surface's yield and sends it off from the strike as one
    emitted electron; a strike on an insulator does the same through `laminate`, or absorbs it (None: the laminate
    charges nothing either way). An electron that leaves the region or the grid has escaped, one that strikes nothing
    within `horizon` RF periods of leaving a surface is dropped as held, and one that strikes more than ROUNDS times
    within a step is left at its face, its weight dropped.

    The number is held near `electrons` by splitting and Russian roulette at the strikes, as transport Monte Carlo
    keeps its weights: with w the population's weight over `electrons` at the step's start, an emission of weight W'
    leaves as m electrons of weight w, m = floor(W'/w + U) for U uniform (none at roulette; past SPLIT, SPLIT
    electrons sharing W'), each drawn from the surface on its own. That keeps the expected weight, and every copy leaves
    the strike with its own draw, so that a population that mostly escapes keeps its number and does not collapse onto
    copies of one flight. Each period's growth is the population's weight at its end over its weight at its start; the
    growth is the mean of their logs over the periods after the first `settle`, and -inf if the population dies."""
    if not (0 <= settle < periods and electrons > 0 and steps >= 8 and horizon > 0):
        raise ValueError(
            f"a run takes 0 <= settle < periods, electrons, steps >= 8 and a horizon > 0, not {settle}, {periods}, "
            f"{electrons}, {steps} and {horizon}"
        )
    rng = np.random.default_rng(seed)
    period = 1 / chamber.f
    dt = period / steps
    drive = _Drive(chamber, QM * volts / gain, 2 * math.pi * chamber.f)
    points, normals, lengths = chamber.faces()
    pick = rng.choice(len(lengths), size=electrons, p=lengths / lengths.sum())
    seeded = surface.emit(rng, np.full(electrons, SEED_ENERGY), np.ones(electrons))
    v_in, vz = _velocity(seeded.energy, seeded.cos, seeded.azimuth, normals[pick])
    sw = _Swarm(points[pick], v_in, vz, rng.uniform(0.0, period, electrons))
    start_weight = float(electrons)
    target = 1.0  # the weight each emission leaves with
    factors: list[float] = []
    impacts_per, escaped_per, held_per = [], [], []
    energy = weighed = on_laminate = 0.0
    for period_index in range(periods):
        counting = period_index >= settle
        count0 = sw.n
        hits = 0
        escaped = held = 0.0
        for k in range(steps):
            t = (period_index * steps + k) * dt
            t_end = t + dt
            n = sw.n
            fresh = np.flatnonzero(sw.waiting[:n] & (sw.born[:n] < t_end) & (sw.weight[:n] > 0))
            if fresh.size:
                sw.waiting[fresh] = False
                a0, _ = drive.accel(sw.x[fresh], sw.born[fresh])
                sw.half[fresh] += a0 * (0.5 * (t_end - sw.born[fresh]))[:, None]
            go = np.flatnonzero(~sw.waiting[:n] & (sw.weight[:n] > 0))
            if not go.size:
                continue
            touched = [go]
            idx, pos, vel = go, sw.x[go], sw.half[go]
            t_start = np.maximum(sw.born[go], t)
            dur = t_end - t_start
            for _ in range(ROUNDS):
                disp = vel * dur[:, None]
                strike = chamber.hit(pos, disp)
                s = strike.box >= 0
                done = idx[~s]
                sw.x[done], sw.half[done], sw.dur[done] = pos[~s] + disp[~s], vel[~s], dur[~s]
                if not s.any():
                    idx = idx[:0]
                    break
                i, tau, d = idx[s], strike.tau[s], dur[s]
                t_hit = t_start[s] + tau * d
                p = pos[s] + tau[:, None] * disp[s]
                a, _ = drive.accel(p, t_hit)
                v = vel[s] + a * ((tau - 0.5) * d)[:, None]
                speed2 = (v**2).sum(axis=1) + sw.vz[i] ** 2
                w_ev = 0.5 * speed2 / QM
                nrm, copper = strike.normal[s], strike.copper[s]
                cos_in = np.abs((v * nrm).sum(axis=1)) / np.sqrt(np.maximum(speed2, 1e-300))
                before = sw.weight[i].copy()
                hits += i.size
                if counting:
                    energy += float(np.sum(before * w_ev))
                    weighed += float(np.sum(before))
                    on_laminate += float(np.sum(before[~copper]))
                y, e_out, c_out, az = _emit(rng, surface, laminate, w_ev, cos_in, copper)
                # each emission leaves as m electrons of the target weight, m its weight over the target rounded at
                # random (none, at roulette); past SPLIT, SPLIT electrons share its weight
                w = before * y
                m = np.floor(w / target + rng.uniform(size=i.size)).astype(np.intp)
                w = np.where(m > SPLIT, w / SPLIT, np.where(m > 0, target, 0.0))
                m = np.minimum(m, SPLIT)
                sw.weight[i], sw.x[i], sw.last[i] = w, p, t_hit
                live = w > 0
                extra = np.where(live, m - 1, 0)
                rep = np.repeat(np.arange(i.size), extra)
                new = sw.add(rep.size)
                if rep.size:
                    _, e_new, c_new, az_new = _emit(rng, surface, laminate, w_ev[rep], cos_in[rep], copper[rep])
                    sw.weight[new], sw.x[new], sw.last[new], sw.born[new] = w[rep], p[rep], t_hit[rep], t_hit[rep]
                    sw.waiting[new] = False
                    touched.append(new)
                seg = np.concatenate([np.flatnonzero(live), rep])
                slots = np.concatenate([i[live], new])
                energies = np.concatenate([e_out[live], e_new]) if rep.size else e_out[live]
                cosines = np.concatenate([c_out[live], c_new]) if rep.size else c_out[live]
                azimuths = np.concatenate([az[live], az_new]) if rep.size else az[live]
                v_new, vz_new = _velocity(energies, cosines, azimuths, nrm[seg])
                rest = (1 - tau[seg]) * d[seg]
                pos = p[seg] + EPS * nrm[seg]
                a0, _ = drive.accel(pos, t_hit[seg])
                sw.vz[slots] = vz_new
                idx, vel, t_start, dur = slots, v_new + a0 * (0.5 * rest)[:, None], t_hit[seg], rest
            if idx.size:  # struck ROUNDS times within the step: left at the face
                sw.weight[idx] = 0.0
            moved = np.concatenate(touched)
            moved = moved[sw.weight[moved] > 0]
            a1, on = drive.accel(sw.x[moved], t_end)
            gone = ~on | ~chamber.inside(sw.x[moved])
            escaped += float(sw.weight[moved[gone]].sum())
            sw.weight[moved[gone]] = 0.0
            ok = moved[~gone]
            sw.half[ok] += a1[~gone] * (0.5 * sw.dur[ok] + 0.5 * dt)[:, None]
            stale = ok[(t_end - sw.last[ok]) > horizon * period]
            held += float(sw.weight[stale].sum())
            sw.weight[stale] = 0.0
        total = float(sw.weight[: sw.n].sum())
        factors.append(total / start_weight)
        if counting:
            impacts_per.append(hits / max(count0, 1))
            escaped_per.append(escaped / start_weight)
            held_per.append(held / start_weight)
        if total <= 0:
            break
        sw.compact()
        # the target follows the population's growth, steered toward `electrons` alive
        target *= (total / start_weight) * min(max(sw.n / electrons, 0.25), 4.0)
        start_weight = total
    kept = np.array(factors[settle:])
    if kept.size and bool(np.all(kept > 0)):
        logs = np.log(kept)
        growth_rate = float(np.mean(logs))
        spread = float(np.std(logs) / math.sqrt(logs.size)) if logs.size > 1 else math.inf
    else:
        growth_rate, spread = -math.inf, math.inf
    return Run(
        growth_rate,
        spread,
        np.array(factors),
        float(np.mean(impacts_per)) if impacts_per else math.nan,
        energy / weighed if weighed > 0 else math.nan,
        on_laminate / weighed if weighed > 0 else math.nan,
        float(np.mean(escaped_per)) if escaped_per else math.nan,
        float(np.mean(held_per)) if held_per else math.nan,
    )


def growth(
    chamber: Chamber,
    volts: float,
    surface: Surface = FP_COPPER,
    *,
    seeds: Sequence[int] = (0, 1, 2, 3),
    laminate: Surface | None = None,
    gain: float = 1.0,
    electrons: int = 1000,
    periods: int = 30,
    settle: int = 10,
    steps: int = 200,
    horizon: float = HORIZON,
) -> Growth:
    """track over seeds: the mean of their growths (Np per RF period) and its spread, their standard deviation."""
    runs = tuple(
        track(
            chamber,
            volts,
            surface,
            laminate=laminate,
            gain=gain,
            electrons=electrons,
            periods=periods,
            settle=settle,
            steps=steps,
            horizon=horizon,
            seed=s,
        )
        for s in seeds
    )
    rates = np.array([r.growth for r in runs])
    if not np.all(np.isfinite(rates)):
        return Growth(-math.inf if np.any(rates == -math.inf) else math.nan, math.nan, runs)
    spread = float(np.std(rates, ddof=1)) if len(rates) > 1 else math.nan
    return Growth(float(np.mean(rates)), spread, runs)


# ---- generations ---------------------------------------------------------------------------------------------------
class Lineage(NamedTuple):
    """What generations measured, over the generations after the first `settle`: the secondaries per generation, each
    flight weighted by exp(-rate T) for its flight of T RF periods (at rate 0 the effective yield), as the geometric
    mean of the generations' with the standard error of its log, and each generation's; and of the flights, the share
    that strike, their mean flight (RF periods), the mean energy of a strike (eV), the share of the strikes on
    insulators and the share of the flights that leave the region."""

    secondaries: float
    spread: float
    multipliers: Arr
    struck: float
    flight: float
    energy: float
    laminate: float
    escaped: float


def _systematic(weight: Arr, count: int, rng: np.random.Generator) -> IntArr:
    """Systematic resampling: count indices drawn in proportion to weight, with the least spread."""
    cumulative = np.cumsum(weight)
    cumulative /= cumulative[-1]
    out: IntArr = np.minimum(
        np.searchsorted(cumulative, (rng.uniform() + np.arange(count)) / count, side="right"), weight.size - 1
    )
    return out


def generations(
    chamber: Chamber,
    volts: float,
    surface: Surface = FP_COPPER,
    *,
    laminate: Surface | None = None,
    gain: float = 1.0,
    rate: float = 0.0,
    electrons: int = 1000,
    count: int = 10,
    settle: int = 3,
    steps: int = 200,
    horizon: float = HORIZON,
    seed: int = 0,
) -> Lineage:
    """The population generation by generation, as multipactor's statistical theory counts it: `electrons` emitted
    electrons (generation 0 the seeds of track) flown to their first impacts (impacts, within `horizon` RF periods),
    each impact weighted by the surface's yield times exp(-rate T), T its flight in RF periods, zero for a flight that
    strikes nothing, escapes, or meets an insulator that absorbs; the next generation `electrons` emissions drawn from
    the impacts in proportion to weight (systematic resampling), each with its own draw from the surface, at its
    impact's time within the RF period. A generation's multiplier is its mean weight. At rate 0 it is the effective
    yield, the spectral radius of the operator that carries one generation's emissions to the next
    (multipactor.effective_yield in two dimensions): above one the population grows, below one it dies."""
    if not (0 <= settle < count and electrons > 0 and steps >= 8 and horizon > 0):
        raise ValueError(
            f"generations take 0 <= settle < count, electrons, steps >= 8 and a horizon > 0, not {settle}, {count}, "
            f"{electrons}, {steps} and {horizon}"
        )
    rng = np.random.default_rng(seed)
    period = 1 / chamber.f
    points, normals, lengths = chamber.faces()
    pick = rng.choice(len(lengths), size=electrons, p=lengths / lengths.sum())
    seeded = surface.emit(rng, np.full(electrons, SEED_ENERGY), np.ones(electrons))
    v, vz = _velocity(seeded.energy, seeded.cos, seeded.azimuth, normals[pick])
    start, velocity = points[pick], np.column_stack([v, vz])
    t0 = rng.uniform(0.0, period, electrons)
    multipliers: list[float] = []
    struck = flight = energy = laminate_share = escaped = 0.0
    strikes = 0
    for generation in range(count):
        imp = impacts(chamber, volts, start, velocity, t0, gain=gain, steps=steps, periods=horizon)
        hit = np.flatnonzero(imp.struck)
        y = np.zeros(electrons)
        if hit.size:
            y[hit], _, _, _ = _emit(rng, surface, laminate, imp.energy[hit], imp.cos[hit], imp.copper[hit])
        took = (imp.time - t0) / period
        w = np.zeros(electrons)
        w[hit] = y[hit] * np.exp(-rate * took[hit])
        multipliers.append(float(w.mean()))
        if generation >= settle:
            struck += hit.size / electrons
            escaped += float(imp.escaped.mean())
            if hit.size:
                strikes += hit.size
                flight += float(took[hit].sum())
                energy += float(imp.energy[hit].sum())
                laminate_share += float((~imp.copper[hit]).sum())
        if not w.sum() > 0:
            break
        pick = _systematic(w, electrons, rng)
        _, e, c, a = _emit(rng, surface, laminate, imp.energy[pick], imp.cos[pick], imp.copper[pick])
        v, vz = _velocity(e, c, a, imp.normal[pick])
        start = imp.point[pick] + EPS * imp.normal[pick]
        velocity = np.column_stack([v, vz])
        t0 = np.mod(imp.time[pick], period)
    kept = np.array(multipliers[settle:])
    measured = max(len(multipliers) - settle, 1)
    if kept.size and bool(np.all(kept > 0)):
        logs = np.log(kept)
        value = float(np.exp(np.mean(logs)))
        spread = float(np.std(logs) / math.sqrt(logs.size)) if logs.size > 1 else math.inf
    else:
        value, spread = 0.0, math.inf
    return Lineage(
        value,
        spread,
        np.array(multipliers),
        struck / measured,
        flight / strikes if strikes else math.nan,
        energy / strikes if strikes else math.nan,
        laminate_share / strikes if strikes else math.nan,
        escaped / measured,
    )


# ---- the parallel-plate limit --------------------------------------------------------------------------------------
def plate_threshold(
    fd: float,
    surface: Surface = FP_COPPER,
    *,
    gain: float = 1.0,
    guess: float | None = None,
    span: tuple[float, float] = (0.6, 1.8),
    iterations: int = 8,
    seeds: Sequence[int] = (0, 1, 2),
    electrons: int = 500,
    periods: int = 24,
    settle: int = 8,
    steps: int = 120,
) -> float:
    """The volts (peak) across parallel plates at f d (GHz mm) where the tracked growth crosses zero: bisection in the
    log of the volts on the sign of the mean growth over `seeds`, the same seeds at every volts, from a bracket of
    `span` times the guess (ECSS's copper chart at f d, or the calibrated statistical threshold off its table), widened
    by those factors until it holds a sign change; nan if none within four widenings."""
    chamber = plates(fd)
    start = guess if guess is not None else float(mp.ecss_copper(fd))
    if not math.isfinite(start):
        start = mp.threshold(fd)

    def rate(volts: float) -> float:
        return growth(
            chamber,
            volts,
            surface,
            seeds=seeds,
            gain=gain,
            electrons=electrons,
            periods=periods,
            settle=settle,
            steps=steps,
        ).rate

    lo, hi = start * span[0], start * span[1]
    r_lo, r_hi = rate(lo), rate(hi)
    for _ in range(4):
        if r_lo < 0 <= r_hi:
            break
        if r_lo >= 0:
            lo, r_lo = lo * span[0], rate(lo * span[0])
        if r_hi < 0:
            hi, r_hi = hi * span[1], rate(hi * span[1])
    else:
        if not r_lo < 0 <= r_hi:
            return math.nan
    for _ in range(iterations):
        mid = math.sqrt(lo * hi)
        if rate(mid) >= 0:
            hi = mid
        else:
            lo = mid
    return math.sqrt(lo * hi)


def calibrate(
    surface: Surface = FP_COPPER,
    target: tuple[float, float] = mp.CALIBRATION,
    *,
    seeds: Sequence[int] = (0, 1, 2),
    electrons: int = 500,
    periods: int = 24,
    settle: int = 8,
    steps: int = 120,
    iterations: int = 8,
) -> float:
    """The volts gain that puts the tracked parallel-plate threshold on the target (copper's lowest chart point, 27.4 V
    at 0.70 GHz mm): the target's volts over the threshold at gain one, as multipactor.calibrate finds GAIN."""
    fd, volts = target
    found = plate_threshold(
        fd,
        surface,
        gain=1.0,
        guess=volts,
        seeds=seeds,
        electrons=electrons,
        periods=periods,
        settle=settle,
        steps=steps,
        iterations=iterations,
    )
    return volts / found


# ---- the protocol's readings ---------------------------------------------------------------------------------------
def r1_section(grid: Grid | None = None) -> Section:
    """R1's rails as designs/r1_rail_ring.toml states them: a 0.40 mm grounded strip, 1.46 mm rails, 0.15 mm spaces, no
    plane under the board, a shield 43.4 null heights out, on 1.6 mm FR-4 in 35 um copper at 915 MHz."""
    geometry = Rails(gap=0.40e-3, rail=1.46e-3, space=0.15e-3, backed=False, walls=43.4)
    return Section(geometry, FR4, 1.6e-3, RIGID_COPPER, 915e6, grid=grid if grid is not None else Grid())


class Case(NamedTuple):
    """A reading the protocol takes: a section and its volts (peak)."""

    name: str
    section: Section
    volts: float

    def chamber(self) -> Chamber:
        """The case's chamber, its section solved."""
        return Chamber.of(solve(self.section))


def facing_case(depth: float, q: float = 0.5, grid: Grid | None = None) -> Case:
    """Section B's facing strips fitted to a well `depth` (eV) deep at Mathieu's q at 915 MHz (the gap scaled as the
    root of the depth, E2-13's family), at the volts that well takes."""
    fitted = fit_well(section_b(grid=grid), depth, q)
    return Case(f"facing strips, {depth:g} eV at q = {q:g}", fitted.section, fitted.volts_for(depth, q))


def cases(grid: Grid | None = None) -> tuple[Case, ...]:
    """The three readings of the brief: T915-S's section B at 50.3 V (1 eV at q 0.5 by MPC-2's facing boards; 0.998 eV
    and q 0.491 as solved), section B with its split widened to h and refitted for 1 eV at q 0.5 (h 0.718 mm, 39.6 V),
    and R1's rails at 50 V."""
    lever = fit_well(section_b(split=1.0, grid=grid))
    return (
        Case("T915-S section B at 50.3 V", section_b(grid=grid), 50.3),
        Case("section B, split h, 1 eV at q 0.5", lever.section, lever.volts_for(1.0, 0.5)),
        Case("R1's rails at 50 V", r1_section(grid), 50.0),
    )


class Reading(NamedTuple):
    """A case read on a surface: the tracked growth per RF period over seeds (growth), the rate a population that dies
    fell at while it lived (its runs' mean decay; nan while it lives), and the effective yield over the same seeds,
    their mean and spread (standard deviation), with each lineage."""

    case: str
    surface: str
    laminate: str
    volts: float
    growth: Growth
    decay: float
    secondaries: float
    secondaries_spread: float
    lineages: tuple[Lineage, ...]

    @property
    def grows(self) -> bool:
        """Multipactor: the population grows, its strikes freeing more than one electron a generation."""
        return self.growth.rate > 0 and self.secondaries > 1


def read(
    case: Case,
    surface: Surface = FP_COPPER,
    *,
    laminate: Surface | None = None,
    gain: float = 1.0,
    chamber: Chamber | None = None,
    seeds: Sequence[int] = (0, 1, 2, 3),
    electrons: int = 1000,
    periods: int = 30,
    settle: int = 10,
    count: int = 10,
    generation_settle: int = 3,
    steps: int = 200,
    horizon: float = HORIZON,
) -> Reading:
    """A case read both ways on a surface, the laminate absorbing (None) or emitting through `laminate`: track over
    `periods` RF periods and generations over `count` generations, on each of `seeds`."""
    ch = chamber if chamber is not None else case.chamber()
    found = growth(
        ch,
        case.volts,
        surface,
        seeds=seeds,
        laminate=laminate,
        gain=gain,
        electrons=electrons,
        periods=periods,
        settle=settle,
        steps=steps,
        horizon=horizon,
    )
    lineages = tuple(
        generations(
            ch,
            case.volts,
            surface,
            laminate=laminate,
            gain=gain,
            electrons=electrons,
            count=count,
            settle=generation_settle,
            steps=steps,
            horizon=horizon,
            seed=s,
        )
        for s in seeds
    )
    yields = np.array([lineage.secondaries for lineage in lineages])
    died = [run.decay for run in found.runs if run.lived < len(run.factors) or run.factors[-1] <= 0]
    return Reading(
        case.name,
        str(getattr(surface, "name", type(surface).__name__)),
        "absorbs" if laminate is None else f"emits as {getattr(laminate, 'name', type(laminate).__name__)}",
        case.volts,
        found,
        float(np.mean(died)) if died else math.nan,
        float(np.mean(yields)),
        float(np.std(yields, ddof=1)) if len(yields) > 1 else math.nan,
        lineages,
    )
