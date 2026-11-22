"""Two-surface multipactor between electrodes driven at f across a gap d: resonances, impact energies, the level-1
chart, copper's secondary-emission yield, the statistical threshold and a tracker of the electron population.

Between plates d apart with V sin(wt) across them, an electron that leaves a plate at drive phase phi with speed
beta w d is, in units of d and of the drive's radians tau since it left (with v = V / K),

    y(tau) = beta tau + v (tau cos phi - sin(tau + phi) + sin phi),    K = (w d)^2 / |e/m| = 224.5 (f d)^2 volts,

f d in GHz mm. It strikes the far plate where y = 1 or its own where y = 0, with (K / 2) y'^2 electron-volts, and
the secondaries it frees leave at phase phi + tau, plus pi on the far plate. Volts and energies scale with K, so a
gap is f d alone: the axis of the published charts.

  volts_scale, flights, impact_energy
                 K; first impacts of many electrons at once, by a scan and bisection, exact in the uniform field
  resonant_volts, resonant_impact, focusing, susceptibility, chart
                 the classical theory (Hatch and Williams; Vaughan 1988): order n (odd) strikes the far plate at
                 tau = n pi, at K (1 - n pi beta) / (n pi cos phi + 2 sin phi) volts, and holds its phase where
                 |beta + n pi v sin phi| < beta + 2 v cos phi; the level-1 chart is the band of volts where such a
                 flight strikes with a yield above one
  Yield          Vaughan's yield (1989, revised 1993) with copper's parameters from ECSS-E-ST-20-01C Table 9-1
  Emission       secondaries leave normal to the plate with a Rayleigh speed (Anza et al. 2010): normal energy
                 exponential with mean 2 W_t
  effective_yield, growth_rate, threshold, calibrate
                 the statistical theory the standard's handbook names for level 1 (S. Anza et al., Phys. Plasmas
                 17, 062110 (2010)): the operator carrying one generation's emission phases to the next, weighted
                 by the yield; its spectral radius is the secondaries per generation once the phases settle, its
                 Euler-Lotka root the growth per RF period, and the lowest volts where the radius is one the
                 chart's lower boundary
  track          the 1-D tracker: a population between the plates, flight by flight, its weight multiplied by the
                 yield at each impact and resampled every RF period; the growth per period
  impact_bound, speed_change_bound, gap_depth                          MPC-1
  facing_volts, facing_depth, facing_fd                                MPC-2
  ECSS_COPPER, ecss_copper, margin_db, Gap, Reading, level1, CASES     MPC-3

**Calibration (DR-0).** The statistical model takes only published inputs: copper's yield from the standard's
Table 9-1 and Anza et al.'s emission at W_t = 1.5 eV. The reference is ECSS's computed chart: Table 5-1 is the
standard's own model run on that Table 9-1 yield, not a measurement, so agreement shows two models of one yield
agree. So taken, the boundary here lies 2.0 % under the standard's 27.4 V at 0.70 GHz mm and within -3.5 to +4.0 % of
every point of ECSS_COPPER from 0.50 to 10.41 GHz mm (-5.4 % at 0.47; 1.9 % rms). Calibration sets one gain on volts
at copper's lowest point (`CALIBRATION`): GAIN = 1.0207, the instrument reading a gap at V as the model reads V /
GAIN; it sits on the low-energy floor d0, which alone moves that point's threshold from 32.0 to 26.8 V. At the
table's 91.9 V at 1.68 GHz mm (`HELD_OUT`), the published point nearest T915-S's 1.674 GHz mm, where the lead design's
margin is read, the instrument gives 91.2 V, -0.8 %: a consistency check, that point having been computed while the
model's structure was chosen. Calibrated, the first-order branch from 0.62 to 2.10 GHz mm is within 1.6 % and
ECSS_COPPER within -3.5 to +6.1 % (2.5 % rms), the most at the order-3 transition near 2.6 GHz mm. **The blind check
(DR-0)** is `BLIND`: the twelve copper points of Table 5-1 between 0.49 and 1.51 GHz mm that ECSS_COPPER leaves out,
registered with DR-0's 10 % band before the model was run on any of them; run once registered, they read -1.7 to
+1.5 % (1.1 % rms).
**Cases that must fail:** the emission's energy calibrated on the same point instead (0.65 eV; the minimum hardly
moves with it) reads 106.4 V there, +16 %, and the classical chart, its one emission energy calibrated the same way
(0.50 eV), 148.9 V, +62 %: the spread of the secondaries' speeds sets the first-order branch, and the low-energy
yield the minimum.

**Furman and Pivi** (M. A. Furman and M. T. F. Pivi, Phys. Rev. ST Accel. Beams 5, 124404 (2002)) is the next
step: true secondaries, rediffused and backscattered electrons, each with its own spectrum, for the 2-D tracker in
the strips' own field (E2-13, MPC-H1) and for ENIG's measured yield (DR-8). It changes the low-energy yield, on
which the minimum rests (the floor d0 alone moves the threshold at 0.70 GHz mm from 32.0 to 26.8 V), and asks for
the calibration again.

The model is one-dimensional: a uniform field between infinite plates, secondaries leaving normal to the plate, a
flight that strikes neither plate within 48 RF periods lost, one electron's flight at a time. The 2-D tracker in the
strips' own field (E2-13, MPC-H1) adds the field's shape, Furman and Pivi the emission's spectra, and a space-charge
model the population's own field. R1's rails read against this chart, -5.2 dB at 50 V, are the parallel-plate bound on
their arcs until that tracker runs in the rail field. SI units, except f d in GHz mm, energies in eV and phases in
radians of the drive.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from functools import cached_property
from typing import NamedTuple

import numpy as np

from array16.particles.base import ME, Arr, BoolArr, E, IntArr

QM = E / ME  # |e/m| (C/kg)
GHZ_MM = 1e6  # m/s: one GHz times one mm
TWO_PI = 2 * math.pi
HORIZON = 96 * math.pi  # rad: a flight that strikes neither plate within 48 RF periods is lost (four periods
# would lose the high orders' flights: 922 V at 10.41 GHz mm against the standard's 521 V)
STEP = 0.05  # rad: the scan for a flight's first impact, then refined by bisection
REFINE = 40  # bisections: the impact's time to STEP / 2^40
BLOCK = 4096  # flights scanned at once
K_LOW, K_HIGH = 0.56, 0.25  # Vaughan's exponents below and above the peak (1993 revision)
PHASES, NODES = 192, 32  # the operator's emission phases and emission energies

Spectrum = Callable[[Arr], Arr]  # impact energies (eV) -> secondaries per impact


# ---- the yield and the emission ------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Yield:
    """Secondary electrons per impact against impact energy (eV), at normal incidence: Vaughan's formula,
    dmax (w e^(1 - w))^k with w = (W - e0) / (emax - e0), k = 0.56 below the peak and 0.25 above it, and
    dmax 1.125 / w^0.35 beyond w = 3.6 (J. R. M. Vaughan, IEEE Trans. Electron Devices 36, 1963 (1989), and its
    revision 40, 830 (1993)). e0 is set so that the curve passes one at the first crossover e1. Below the peak
    the yield is at least d0, the value at zero energy, where reflection holds it up; the curve meets d0 near
    12.7 eV for copper.

    The defaults are copper's: e1 = 19 eV, dmax = 2.48 at emax = 232 eV, d0 = 0.8, ECSS-E-ST-20-01C (15 June
    2020) Table 9-1, the average of more than 20 air-exposed samples at normal incidence and room temperature, from
    which the standard's copper chart is computed. The knotted study's Yield, ported with the 1993 tail."""

    e1: float = 19.0
    dmax: float = 2.48
    emax: float = 232.0
    d0: float = 0.8

    def __post_init__(self) -> None:
        if not (0 < self.e1 < self.emax and self.dmax > 1 and 0 <= self.d0 < 1):
            raise ValueError(
                f"a yield needs 0 < e1 < emax, dmax > 1 and 0 <= d0 < 1, not e1 {self.e1}, emax {self.emax}, "
                f"dmax {self.dmax}, d0 {self.d0}"
            )

    @cached_property
    def e0(self) -> float:
        """Where Vaughan's curve starts (eV): w1 e^(1 - w1) = dmax^(-1/0.56) below the peak puts one at e1, and
        e0 = (e1 - w1 emax) / (1 - w1); 0.828 eV for copper."""
        target = self.dmax ** (-1 / K_LOW)
        lo, hi = 0.0, 1.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if mid * math.exp(1 - mid) < target else (lo, mid)
        w1 = 0.5 * (lo + hi)
        return (self.e1 - w1 * self.emax) / (1 - w1)

    def curve(self, energy: Arr | float) -> Arr:
        """Vaughan's curve alone (zero below e0)."""
        e = np.asarray(energy, dtype=np.float64)
        w = (e - self.e0) / (self.emax - self.e0)
        safe = np.maximum(w, 1e-300)
        body = self.dmax * (safe * np.exp(1 - safe)) ** np.where(w < 1, K_LOW, K_HIGH)
        tail = self.dmax * 1.125 / np.maximum(safe, 3.6) ** 0.35
        out: Arr = np.where(w <= 0, 0.0, np.where(w > 3.6, tail, body))
        return out

    def __call__(self, energy: Arr | float) -> Arr:
        """The yield: the curve, held up to d0 below the peak."""
        e = np.asarray(energy, dtype=np.float64)
        d = self.curve(e)
        out: Arr = np.where(e < self.emax, np.maximum(d, self.d0), d)
        return out

    @cached_property
    def e2(self) -> float:
        """The second crossover (eV), where the curve falls back to one above the peak: 4.3 keV for copper."""
        lo, hi = self.emax, self.emax * 1e4
        for _ in range(200):
            mid = math.sqrt(lo * hi)
            lo, hi = (mid, hi) if float(self.curve(mid)) > 1 else (lo, mid)
        return math.sqrt(lo * hi)


COPPER = Yield()  # ECSS-E-ST-20-01C Table 9-1, as published


@dataclass(frozen=True)
class Emission:
    """How secondaries leave a plate: normal to it, with speeds of Rayleigh's distribution, (v / v_t^2)
    exp(-v^2 / 2 v_t^2), thermal energy W_t = m v_t^2 / 2e (S. Anza et al., "Nonstationary statistical theory for
    multipactor", Phys. Plasmas 17, 062110 (2010), their eq. 12, with W_t = 1.5 eV). The normal energy is then
    exponential with mean 2 W_t. One dimension keeps only the normal motion, which the field changes."""

    thermal: float = 1.5  # eV

    def __post_init__(self) -> None:
        if not self.thermal > 0:
            raise ValueError(f"an emission needs a positive thermal energy, not {self.thermal} eV")

    @property
    def mean(self) -> float:
        """The mean normal energy (eV): 2 W_t."""
        return 2 * self.thermal

    def draw(self, rng: np.random.Generator, count: int) -> Arr:
        """count normal energies (eV) at random."""
        out: Arr = rng.exponential(self.mean, count)
        return out

    def nodes(self, count: int) -> Arr:
        """count normal energies (eV) of equal weight: the midpoints of count equal slices of probability."""
        out: Arr = -self.mean * np.log1p(-(np.arange(count) + 0.5) / count)
        return out


EMISSION = Emission()  # Anza et al.'s W_t = 1.5 eV, kept as published


# ---- flights in the gap --------------------------------------------------------------------------------------------
def volts_scale(fd: float) -> float:
    """K = (2 pi f d)^2 / |e/m| (V), f d in GHz mm: 224.5 (f d)^2. A gap's volts and energies scale with it."""
    return (TWO_PI * fd * GHZ_MM) ** 2 / QM


def _positive(fd: float) -> float:
    """K for a gap that has one: f d must be positive."""
    if not fd > 0:
        raise ValueError(f"a gap needs a positive f d, not {fd} GHz mm")
    return volts_scale(fd)


class Flights(NamedTuple):
    """Where electrons strike, one entry each."""

    tau: Arr  # rad of the drive from leaving to striking; nan where lost
    speed: Arr  # |dy/dtau| at the strike, in units of w d; 0 where lost
    far: BoolArr  # struck the far plate, else its own
    hit: BoolArr  # struck a plate within the horizon


def flights(
    phase: Arr | float, beta: Arr | float, v: Arr | float, horizon: float = HORIZON, step: float = STEP
) -> Flights:
    """The first impact of each electron that leaves a plate at drive phase `phase` with speed beta w d, across a
    gap of v = V / K: its position, in units of d, is beta tau + v (tau cos phi - sin(tau + phi) + sin phi). The
    arguments broadcast together and are flattened. A scan in steps of `step` finds the first sample past either
    plate, bisection the crossing to rounding; a crossing that grazes a plate and returns within one step, by at
    most v step^2 / 8 of the gap, goes unseen. An electron that leaves with no speed into a field that pushes it
    back (beta = 0, v sin phi < 0) never leaves: it strikes its own plate at once, with no speed."""
    p, b, vv = (a.ravel() for a in np.broadcast_arrays(*(np.asarray(x, dtype=np.float64) for x in (phase, beta, v))))
    n = p.size
    c, s = np.cos(p), np.sin(p)
    lo, hi = np.zeros(n), np.zeros(n)
    far: BoolArr = np.zeros(n, dtype=np.bool_)
    stuck = (b <= 0) & (vv * s < 0)
    hit: BoolArr = stuck.copy()
    grid = step * np.arange(1, math.ceil(TWO_PI / step) + 1)
    for first in range(0, n, BLOCK):
        idx: IntArr = np.arange(first, min(n, first + BLOCK))
        idx = idx[~stuck[idx]]
        start = 0.0
        while idx.size and start < horizon:
            t = start + grid
            y = b[idx, None] * t + vv[idx, None] * (t * c[idx, None] - np.sin(t + p[idx, None]) + s[idx, None])
            out = (y >= 1.0) | (y <= 0.0)
            got = out.any(axis=1)
            j = out.argmax(axis=1)[got]
            k = idx[got]
            hi[k] = start + grid[j]
            lo[k] = np.where(j > 0, start + grid[np.maximum(j - 1, 0)], start)
            far[k] = y[got, j] >= 1.0
            hit[k] = True
            idx = idx[~got]
            start += grid[-1]
    k = np.flatnonzero(hit & ~stuck)
    a, z, pk, bk, vk, ck, sk, fk = lo[k], hi[k], p[k], b[k], vv[k], c[k], s[k], far[k]
    for _ in range(REFINE):
        mid = 0.5 * (a + z)
        y = bk * mid + vk * (mid * ck - np.sin(mid + pk) + sk)
        past = np.where(fk, y >= 1.0, y <= 0.0)
        a, z = np.where(past, a, mid), np.where(past, mid, z)
    tau, speed = np.full(n, np.nan), np.zeros(n)
    tau[stuck] = 0.0
    tau[k] = z
    speed[k] = np.abs(bk + vk * (ck - np.cos(z + pk)))
    return Flights(tau, speed, far, hit)


def impact_energy(fd: float, volts: Arr | float, phase: Arr | float, w0: float = 0.0) -> Arr:
    """The energy (eV) of an electron's first impact, leaving at `phase` with normal energy w0 (eV) across a gap
    f d (GHz mm) with `volts` peak across it; nan where it strikes neither plate within the horizon."""
    k = _positive(fd)
    fl = flights(phase, math.sqrt(2 * w0 / k), np.asarray(volts, dtype=np.float64) / k)
    out: Arr = np.where(fl.hit, 0.5 * k * fl.speed**2, np.nan)
    shape = np.broadcast_shapes(np.shape(volts), np.shape(phase))
    return out.reshape(shape)


# ---- the classical theory and the level-1 chart --------------------------------------------------------------------
def _beta(fd: float, w0: float) -> float:
    return math.sqrt(2 * w0 / _positive(fd))


def resonant_volts(fd: float, order: int, phase: Arr | float, w0: float = 0.0) -> Arr:
    """The volts (peak) at which an electron leaving at `phase` with normal energy w0 (eV) strikes the far plate
    after `order` half periods (odd): K (1 - n pi beta) / (n pi cos phi + 2 sin phi), K = volts_scale(fd)."""
    k, beta, ph = volts_scale(fd), _beta(fd, w0), np.asarray(phase, dtype=np.float64)
    out: Arr = k * (1 - order * math.pi * beta) / (order * math.pi * np.cos(ph) + 2 * np.sin(ph))
    return out


def resonant_impact(fd: float, order: int, phase: Arr | float, w0: float = 0.0) -> Arr:
    """The impact energy (eV) of that resonant flight: (K / 2)(beta + 2 v cos phi)^2, v = V / K."""
    k, beta, ph = volts_scale(fd), _beta(fd, w0), np.asarray(phase, dtype=np.float64)
    v = resonant_volts(fd, order, ph, w0) / k
    out: Arr = 0.5 * k * (beta + 2 * v * np.cos(ph)) ** 2
    return out


def focusing(fd: float, order: int, phase: Arr | float, w0: float = 0.0) -> Arr:
    """How a small error in the emission phase grows from one resonant flight to the next:
    G = (beta + n pi v sin phi) / (beta + 2 v cos phi). The resonance holds its phase where |G| < 1 (phase
    focusing; for w0 = 0, 0 <= tan phi < 2 / (n pi))."""
    k, beta, ph = volts_scale(fd), _beta(fd, w0), np.asarray(phase, dtype=np.float64)
    v = resonant_volts(fd, order, ph, w0) / k
    out: Arr = (beta + order * math.pi * v * np.sin(ph)) / (beta + 2 * v * np.cos(ph))
    return out


class Band(NamedTuple):
    """One order's band of the level-1 chart at one f d."""

    order: int
    low: float  # V: the lowest volts at which a focused resonant flight strikes with a yield above one
    high: float  # V: the highest
    phases: tuple[float, float]  # rad: the emission phases at low and at high
    energies: tuple[float, float]  # eV: the impact energies at low and at high


def susceptibility(
    fd: float,
    sey: Spectrum = COPPER,
    orders: Sequence[int] = (1, 3, 5, 7, 9, 11),
    w0: float = 0.0,
    samples: int = 1201,
) -> tuple[Band, ...]:
    """The level-1 chart at f d (GHz mm): for each order (odd), the band of volts where an electron leaving with
    normal energy w0 (eV) at a phase that phase focusing holds (|G| < 1) strikes the far plate after exactly that
    many half periods, its flight staying inside the gap, with a yield above one. The phases run over `samples`
    points of (-pi/2, atan(2 / (n pi))), each flight checked by `flights`. Orders with no such flight are left out."""
    k, beta = volts_scale(fd), _beta(fd, w0)
    bands = []
    for n in orders:
        if n < 1 or n % 2 == 0:
            raise ValueError(f"two-surface resonances have odd orders, not {n}")
        if n * math.pi * beta >= 1:
            continue
        phase = np.linspace(-math.pi / 2, math.atan2(2.0, n * math.pi), samples)[1:-1]
        denom = n * math.pi * np.cos(phase) + 2 * np.sin(phase)
        v = (1 - n * math.pi * beta) / np.where(denom > 0, denom, np.inf)
        speed = beta + 2 * v * np.cos(phase)
        focus = (beta + n * math.pi * v * np.sin(phase)) / np.where(speed > 0, speed, np.inf)
        ok = (denom > 0) & (speed > 0) & (np.abs(focus) < 1)
        fl = flights(phase[ok], beta, v[ok])
        ok[ok] = fl.far & (np.abs(fl.tau - n * math.pi) < 1e-6)
        energy = 0.5 * k * speed**2
        ok &= np.asarray(sey(energy)) > 1
        if not ok.any():
            continue
        volts = v * k
        i, j = int(np.flatnonzero(ok)[np.argmin(volts[ok])]), int(np.flatnonzero(ok)[np.argmax(volts[ok])])
        bands.append(
            Band(
                n,
                float(volts[i]),
                float(volts[j]),
                (float(phase[i]), float(phase[j])),
                (float(energy[i]), float(energy[j])),
            )
        )
    return tuple(bands)


def chart(
    fds: Arr | Sequence[float], sey: Spectrum = COPPER, w0: float = 0.0, orders: Sequence[int] = (1, 3, 5, 7, 9, 11)
) -> tuple[Arr, Arr]:
    """The level-1 chart's lowest and highest volts over the orders at each f d (nan where no order multipacts)."""
    low, high = [], []
    for fd in np.asarray(fds, dtype=np.float64):
        bands = susceptibility(float(fd), sey, orders, w0)
        low.append(min((b.low for b in bands), default=math.nan))
        high.append(max((b.high for b in bands), default=math.nan))
    return np.array(low), np.array(high)


# ---- the statistical theory: the calibrated instrument -------------------------------------------------------------
@dataclass(frozen=True, eq=False)
class _Transfer:
    """One generation of flights from `phases` emission phases and `nodes` emission energies, at several volts:
    each flight's source phase bin, its target bins and weights (linear in the target phase), its yield and its
    time. The operator weights them by yield times exp(-rate tau / 2 pi)."""

    phases: int
    count: int  # volts
    source: IntArr  # (count, phases * nodes)
    target: IntArr  # (count, phases * nodes, 2)
    share: Arr  # (count, phases * nodes, 2): the weight of each target bin, over the nodes
    yields: Arr  # (count, phases * nodes): the yield at the strike over the nodes, zero where lost
    tau: Arr  # (count, phases * nodes): rad

    def radius(self, which: int, rate: float = 0.0) -> float:
        """The spectral radius of the operator at volts number `which`, each flight weighted by
        exp(-rate tau / 2 pi), rate being a growth per RF period in nepers."""
        m = self.phases
        weight = self.yields[which] * np.exp(-rate * self.tau[which] / TWO_PI)
        op = np.zeros((m, m))
        for side in (0, 1):
            np.add.at(op, (self.target[which, :, side], self.source[which]), weight * self.share[which, :, side])
        return float(np.max(np.abs(np.linalg.eigvals(op))))


def _transfer(fd: float, volts: Arr, sey: Spectrum, emission: Emission, phases: int, nodes: int) -> _Transfer:
    k = _positive(fd)
    phi = (np.arange(phases) + 0.5) * TWO_PI / phases
    beta = np.sqrt(2 * emission.nodes(nodes) / k)
    count = volts.size
    p = np.broadcast_to(np.repeat(phi, nodes), (count, phases * nodes))
    b = np.broadcast_to(np.tile(beta, phases), (count, phases * nodes))
    fl = flights(p, b, (volts / k)[:, None] * np.ones((1, phases * nodes)))
    shape = (count, phases * nodes)
    tau = np.nan_to_num(fl.tau).reshape(shape)
    far = fl.far.reshape(shape)
    yields = np.where(fl.hit.reshape(shape), sey(0.5 * k * fl.speed.reshape(shape) ** 2), 0.0) / nodes
    where = np.mod(p + tau + np.where(far, math.pi, 0.0), TWO_PI) * phases / TWO_PI - 0.5
    below = np.floor(where)
    frac = where - below
    target = np.stack([below.astype(np.intp) % phases, (below.astype(np.intp) + 1) % phases], axis=-1)
    share = np.stack([1 - frac, frac], axis=-1)
    source = np.broadcast_to(np.repeat(np.arange(phases), nodes), shape).astype(np.intp)
    return _Transfer(phases, count, source, target, share, yields, tau)


# ECSS-E-ST-20-01C (15 June 2020) Table 5-1, copper: the lowest two-surface breakdown threshold, f x d (GHz mm)
# against peak volts, computed by the standard from its Table 9-1 yield; the knotted study's ECSS_CU, checked
# against the standard's table (which lists more points, to 100 GHz mm).
ECSS_COPPER: tuple[tuple[float, float], ...] = (
    (0.47, 36.6), (0.50, 32.5), (0.56, 29.7), (0.62, 28.1), (0.70, 27.4), (0.82, 28.7), (0.92, 31.9), (1.02, 37.0),
    (1.14, 44.4), (1.28, 54.2), (1.43, 67.0), (1.59, 83.3), (1.68, 91.9), (1.78, 101.4), (1.88, 110.8), (1.99, 120.9),
    (2.10, 130.1), (2.22, 139.7), (2.34, 148.2), (2.48, 156.7), (2.62, 160.2), (2.77, 154.3), (2.92, 149.2),
    (3.09, 150.6), (3.27, 155.2), (3.45, 163.0), (3.65, 171.7), (3.85, 183.4), (4.07, 195.8), (4.30, 211.2),
    (4.55, 227.6), (4.81, 244.6), (5.08, 262.5), (5.37, 277.3), (5.67, 290.7), (5.99, 302.9), (6.33, 314.4),
    (6.69, 329.6), (7.07, 350.3), (7.47, 372.1), (7.90, 394.5), (8.34, 418.0), (8.82, 441.6), (9.32, 466.6),
    (9.85, 493.1), (10.41, 521.0),
)  # fmt: skip
CALIBRATION = (0.70, 27.4)  # copper's lowest two-surface point: the calibration
HELD_OUT = (1.68, 91.9)  # the published point nearest T915-S's 1.674 GHz mm: a consistency check, at 10 %
# ECSS-E-ST-20-01C Table 5-1, copper: the twelve points between 0.47 and 1.59 GHz mm that ECSS_COPPER leaves out, every
# other point of the table there, on the first-order branch where T915-S's margin is read. DR-0's blind check: no run
# of this model had computed them when they were registered here with BLIND_BAND (2026-11-22).
BLIND: tuple[tuple[float, float], ...] = (
    (0.49, 33.4), (0.53, 30.8), (0.59, 28.8), (0.66, 27.7), (0.74, 27.5), (0.78, 27.7),
    (0.87, 29.9), (0.97, 34.0), (1.08, 40.3), (1.21, 49.0), (1.35, 60.3), (1.51, 74.8),
)  # fmt: skip
BLIND_BAND = 0.10  # DR-0: ECSS to 10 %
GAIN = 1.0207  # the volts gain `calibrate` finds at CALIBRATION: the instrument reads V as the model reads V / GAIN


def _radii(fd: float, volts: Arr, sey: Spectrum, emission: Emission, phases: int, nodes: int, rate: float) -> Arr:
    op = _transfer(fd, volts, sey, emission, phases, nodes)
    return np.array([op.radius(i, rate) for i in range(volts.size)])


def effective_yield(
    fd: float,
    volts: Arr | float,
    sey: Spectrum = COPPER,
    emission: Emission = EMISSION,
    *,
    gain: float = GAIN,
    phases: int = PHASES,
    nodes: int = NODES,
) -> Arr:
    """The secondaries per generation once the emission phases have settled, at each of `volts` across a gap f d
    (GHz mm), read through the volts gain: the spectral radius of the operator T(phi' <- phi), the mean yield of
    the flights from phase phi that free secondaries at phi', over the emission energies (the stationary limit of
    Anza et al.'s theory). Above one the population grows, below one it dies. The operator is discretised on
    `phases` emission phases and `nodes` emission energies."""
    vs = np.atleast_1d(np.asarray(volts, dtype=np.float64)) / gain
    return _radii(fd, vs, sey, emission, phases, nodes, 0.0)


def growth_rate(
    fd: float,
    volts: float,
    sey: Spectrum = COPPER,
    emission: Emission = EMISSION,
    *,
    gain: float = GAIN,
    phases: int = PHASES,
    nodes: int = NODES,
) -> float:
    """The population's growth per RF period, in nepers, once its phases have settled: the rate r at which the
    operator with each flight weighted by exp(-r tau / 2 pi) has spectral radius one (the Euler-Lotka equation of
    a branching process typed by phase). Zero where the effective yield is one; -inf where every flight is lost."""
    op = _transfer(fd, np.array([volts / gain]), sey, emission, phases, nodes)
    if not op.yields.any():
        return -math.inf
    lo, hi = -1.0, 1.0
    for _ in range(60):
        if op.radius(0, lo) >= 1:
            break
        lo *= 2
    else:
        return -math.inf
    for _ in range(60):
        if op.radius(0, hi) <= 1:
            break
        hi *= 2
    else:
        return math.inf
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if op.radius(0, mid) > 1 else (lo, mid)
    return 0.5 * (lo + hi)


def threshold(
    fd: float,
    sey: Spectrum = COPPER,
    emission: Emission = EMISSION,
    *,
    growth: float = 0.0,
    gain: float = GAIN,
    lowest: float = 5.0,
    highest: float = 4000.0,
    ratio: float = 1.12,
    rtol: float = 1e-3,
    phases: int = PHASES,
    nodes: int = NODES,
) -> float:
    """The lowest volts (peak) across a gap f d (GHz mm) at which the population grows by `growth` nepers per RF
    period or more: where the operator, each flight weighted by exp(-growth tau / 2 pi), reaches spectral radius
    one. With growth = 0 it is where the effective yield is one, the boundary ECSS-E-ST-20-01C Table 5-1 lists. A
    coarse operator scans the model's volts from `lowest` in steps of `ratio`, the full one checks its bracket and
    bisects it to `rtol`; nan if none to `highest`, or if f d is not positive (no flight resonates in a gap of
    none); the result is read through the gain."""
    if not fd > 0:
        return math.nan
    grid = np.geomspace(lowest, highest, math.ceil(math.log(highest / lowest) / math.log(ratio)) + 1)

    def full(volts: Arr) -> Arr:
        return _radii(fd, volts, sey, emission, phases, nodes, growth)

    coarse = _radii(fd, grid, sey, emission, max(32, phases // 3), max(8, nodes // 3), growth)
    above = np.flatnonzero(coarse >= 1)
    found = False
    if above.size:
        i = int(above[0])
        lo, hi = float(grid[max(i - 2, 0)]), float(grid[min(i + 1, grid.size - 1)])
        ends = full(np.array([lo, hi]))
        found = bool(ends[0] < 1 <= ends[1])
    if not found:  # the coarse scan misled or found nothing: scan at full resolution
        hits = np.flatnonzero(full(grid) >= 1)
        if not hits.size:
            return math.nan
        i = int(hits[0])
        if i == 0:
            return gain * float(grid[0])
        lo, hi = float(grid[i - 1]), float(grid[i])
    while hi / lo > 1 + rtol:
        mid = math.sqrt(lo * hi)
        lo, hi = (lo, mid) if float(full(np.array([mid]))[0]) >= 1 else (mid, hi)
    return gain * math.sqrt(lo * hi)


def calibrate(
    target: tuple[float, float] = CALIBRATION,
    sey: Spectrum = COPPER,
    emission: Emission = EMISSION,
    *,
    phases: int = PHASES,
    nodes: int = NODES,
) -> float:
    """The volts gain that puts the model's threshold on the target, the target volts over the model's own
    threshold (gain 1) at the target f d (GHz mm). On ECSS copper's lowest point it is GAIN."""
    fd, volts = target
    return volts / threshold(fd, sey, emission, gain=1.0, rtol=1e-6, phases=phases, nodes=nodes)


# ---- the tracker ---------------------------------------------------------------------------------------------------
class Population(NamedTuple):
    """What track measured."""

    growth: float  # nepers per RF period: the mean over the periods after the first `settle`; -inf if all lost
    spread: float  # the standard error of that mean
    factors: Arr  # the population's growth over each period
    impacts: float  # impacts per electron per period
    energy: float  # mean impact energy (eV), weighted by the electrons' weights
    far: float  # the share of impacts on the far plate


def _resample(weight: Arr, count: int, rng: np.random.Generator) -> IntArr:
    """Systematic resampling: count indices drawn in proportion to weight, with the least spread."""
    cumulative = np.cumsum(weight)
    cumulative /= cumulative[-1]
    out: IntArr = np.minimum(
        np.searchsorted(cumulative, (rng.uniform() + np.arange(count)) / count, side="right"), weight.size - 1
    )
    return out


def track(
    fd: float,
    volts: float,
    sey: Spectrum = COPPER,
    emission: Emission = EMISSION,
    *,
    gain: float = GAIN,
    electrons: int = 2000,
    periods: int = 40,
    settle: int = 8,
    seed: int = 0,
) -> Population:
    """A population of electrons between plates f d (GHz mm) apart with `volts` (peak) across them, followed flight
    by flight. They leave either plate at random phases of the first period; each impact multiplies an electron's
    weight by the yield at its energy and sends it off again from that plate with a fresh normal energy drawn from
    the emission; a flight that strikes neither plate within the horizon is lost. At the end of every RF period the
    population's weight over its number is that period's growth, and it is resampled in proportion to weight
    (systematic resampling), which keeps its number and spends it where the weight is. The growth is the mean of
    the log over the periods after the first `settle`, in nepers per period: above zero, multipactor. The volts
    are read through the gain, as by `threshold`."""
    if not (0 <= settle < periods and electrons > 0):
        raise ValueError(f"a run needs 0 <= settle < periods and electrons, not {settle}, {periods}, {electrons}")
    em = emission
    rng = np.random.default_rng(seed)
    k = _positive(fd)
    v = volts / gain / k
    n = electrons

    def launch(phase: Arr) -> Flights:
        return flights(phase, np.sqrt(2 * em.draw(rng, phase.size) / k), v)

    left = rng.uniform(0.0, TWO_PI, n)  # drive phase (rad) at which each left a plate
    phase = np.mod(left + math.pi * rng.integers(0, 2, n), TWO_PI)  # the field's phase in its own sense
    fl = launch(phase)
    strike = left + np.where(fl.hit, fl.tau, np.inf)
    speed, far = fl.speed, fl.far
    weight = np.where(fl.hit, 1.0, 0.0)
    factors: list[float] = []
    hits = 0
    energy = far_hits = weighed = 0.0
    for period in range(periods):
        end = TWO_PI * (period + 1)
        while True:
            due = np.flatnonzero((weight > 0) & (strike < end))
            if not due.size:
                break
            w = 0.5 * k * speed[due] ** 2
            hits += due.size
            far_hits += float(far[due].sum())
            energy += float(np.sum(weight[due] * w))
            weighed += float(np.sum(weight[due]))
            weight[due] *= sey(w)
            phase[due] = np.mod(phase[due] + (strike[due] - left[due]) + np.where(far[due], math.pi, 0.0), TWO_PI)
            left[due] = strike[due]
            fl = launch(phase[due])
            strike[due] = left[due] + np.where(fl.hit, fl.tau, np.inf)
            speed[due], far[due] = fl.speed, fl.far
            weight[due] = np.where(fl.hit, weight[due], 0.0)
        total = float(weight.sum())
        factors.append(total / n)
        if total <= 0:
            break
        pick = _resample(weight, n, rng)
        left, phase, strike, speed, far = left[pick], phase[pick], strike[pick], speed[pick], far[pick]
        weight = np.ones(n)
    kept = np.array(factors[settle:])
    if kept.size and bool(np.all(kept > 0)):
        logs = np.log(kept)
        growth = float(np.mean(logs))
        spread = float(np.std(logs) / math.sqrt(logs.size)) if logs.size > 1 else math.inf
    else:  # every electron lost
        growth, spread = -math.inf, math.inf
    return Population(
        growth,
        spread,
        np.array(factors),
        hits / (n * len(factors)),
        energy / weighed if weighed > 0 else math.nan,
        far_hits / hits if hits else math.nan,
    )


# ---- MPC-1: the impact bound ---------------------------------------------------------------------------------------
def speed_change_bound(field: float, f: float) -> float:
    """The most an electron's speed changes in a uniform field E0 sin(wt + phi): 2 e E0 / (m w) (m/s), since
    v(t) - v(0) = (e E0 / m w)(cos phi - cos(wt + phi)) (MPC-1)."""
    return 2 * QM * field / (TWO_PI * f)


def impact_bound(w0: float, depth: float) -> float:
    """The most energy (eV) an electron that left with w0 (eV) brings to a plate through a field whose
    ponderomotive energy is `depth` (eV): (sqrt w0 + sqrt(8 D))^2, since 2 e E0 / (m w) is sqrt(8 D) in energy
    (MPC-1). Between facing strips D is the well's depth, the barrier field's ponderomotive energy."""
    return (math.sqrt(w0) + math.sqrt(8 * depth)) ** 2


def gap_depth(fd: float, volts: float) -> float:
    """The ponderomotive energy (eV) of a gap's field V / d: |e/m| (V / d)^2 / (4 w^2) = V^2 / (4 K)."""
    return volts**2 / (4 * _positive(fd))


# ---- MPC-2: facing boards at fixed q -------------------------------------------------------------------------------
def facing_volts(q: float, fd: float) -> float:
    """The volts (peak, half to half) that give Mathieu's q on facing boards a gap d = 2h apart at f: the null's
    gradient pi V / d^2 (ELE-2) and q = 2 |e/m| G / w^2 give V = 2 pi q (f d)^2 / |e/m|, 35.72 q (f d)^2 with
    f d in GHz mm (MPC-2)."""
    return TWO_PI * q * (fd * GHZ_MM) ** 2 / QM


def facing_depth(q: float, fd: float) -> float:
    """The well's depth (eV) on facing boards at Mathieu's q: q V / (8 pi) (ELE-2), so q^2 (f d)^2 / (4 |e/m|),
    1.4214 q^2 (f d)^2 with f d in GHz mm (MPC-2). At fixed q the chart point (f d, V) follows from the depth."""
    return q * q * (fd * GHZ_MM) ** 2 / (4 * QM)


def facing_fd(depth: float, q: float) -> float:
    """The f d (GHz mm) at which facing boards hold a well `depth` (eV) deep at Mathieu's q. Inverse of
    facing_depth."""
    return math.sqrt(4 * QM * depth) / (q * GHZ_MM)


# ---- MPC-3: level-1 margins against ECSS copper --------------------------------------------------------------------
def ecss_copper(fd: Arr | float) -> Arr:
    """ECSS copper's lowest threshold (V, peak) at f d (GHz mm), linear between the table's points; nan outside
    0.47 to 10.41 GHz mm, the stretch copied (the standard lists no copper point below 0.47 GHz mm, and runs on to
    100 GHz mm)."""
    x = np.asarray(fd, dtype=np.float64)
    table = np.array(ECSS_COPPER)
    out: Arr = np.where((x >= table[0, 0]) & (x <= table[-1, 0]), np.interp(x, table[:, 0], table[:, 1]), np.nan)
    return out


def margin_db(volts: float, threshold: float) -> float:
    """The level-1 margin in dB of power: 20 log10(threshold / volts). ECSS asks 6 dB for metal-only parts on its
    own yield (level 1), 7 to 8 dB at level 2."""
    return 20 * math.log10(threshold / volts)


class Gap(NamedTuple):
    """A place where electrons may multipact: two surfaces at `volts` (peak) apart at frequency f, joined by field
    lines from `shortest` to `longest` metres. Facing plates have one length; thin rods and electrodes side by side
    on one board are joined by arcs from the gap upward, every one of them across the full volts."""

    name: str
    f: float  # Hz
    volts: float  # V, peak
    shortest: float  # m
    longest: float  # m

    @property
    def span(self) -> tuple[float, float]:
        """The f d (GHz mm) of the shortest and the longest field line."""
        return self.f * self.shortest / 1e6, self.f * self.longest / 1e6


class Reading(NamedTuple):
    """A gap read against a chart."""

    gap: Gap
    fd: float  # GHz mm: the field line on which the threshold is lowest
    threshold: float  # V
    margin: float  # dB


def level1(gap: Gap, chart: Callable[[float], float] | None = None) -> Reading:
    """The level-1 reading of a gap: the lowest threshold over its field lines' f d, at the span's ends and at the
    table's points inside it, and the margin of its volts below that (ECSS copper unless another chart is given,
    such as a calibrated statistical threshold)."""
    lo, hi = gap.span
    where = [lo, hi, *(x for x, _ in ECSS_COPPER if lo < x < hi)]
    look = chart or (lambda x: float(ecss_copper(x)))
    values = [(look(x), x) for x in where]
    found = [(t, x) for t, x in values if math.isfinite(t)]
    if not found:
        return Reading(gap, math.nan, math.nan, math.inf)
    t, x = min(found)
    return Reading(gap, x, t, margin_db(gap.volts, t))


F915 = 915e6
ROD_R0, ROD_WIRE = 0.5e-3, 25e-6  # m: the first run's rods, 0.5 mm from the axis, 25 um in radius
ROD_GAP = ROD_R0 * math.sqrt(2) - 2 * ROD_WIRE  # m: neighbouring rods' surfaces, 0.657 mm apart
R1_RAIL = 1.46e-3  # m: R1's rails either side of its 0.40 mm grounded strip (docs/KNOTTED_GUIDE.md section 3.6)


def _facing(depth: float, q: float = 0.5, f: float = F915) -> Gap:
    fd = facing_fd(depth, q)
    gap = fd * 1e6 / f
    return Gap(f"facing boards, {depth:g} eV at q = {q:g}", f, facing_volts(q, fd), gap, gap)


# The gaps of section 10's designs read at level 1. The rods' neighbours are joined by circular arcs through both,
# from the chord to the semicircle on it and beyond; R1's rails meet the grounded strip and the outer ground along
# arcs from the slit to the rails' width, pi times the radius long (the gapless model's slit is zero), each read as a
# pair of parallel plates that far apart: a bound, the rails' field not being uniform along an arc.
CASES: tuple[Gap, ...] = (
    Gap("first run's rods, q = 0.3", 1e9, 13.0, ROD_GAP, ROD_GAP * math.pi / 2),
    Gap("first run's rods, q = 0.6", 1e9, 25.0, ROD_GAP, ROD_GAP * math.pi / 2),
    Gap("T915-S, 1 eV at q = 0.5", F915, 50.0, 1.83e-3, 1.83e-3),
    Gap("R1's rails at 50 V, the parallel-plate bound", F915, 50.0, 0.0, math.pi * R1_RAIL),
    _facing(2.5),
    _facing(3.0),
)

__all__ = [
    "BLIND",
    "BLIND_BAND",
    "BLOCK",
    "CALIBRATION",
    "CASES",
    "COPPER",
    "ECSS_COPPER",
    "EMISSION",
    "F915",
    "GAIN",
    "GHZ_MM",
    "HELD_OUT",
    "HORIZON",
    "NODES",
    "PHASES",
    "QM",
    "ROD_GAP",
    "STEP",
    "Band",
    "Emission",
    "Flights",
    "Gap",
    "Population",
    "Reading",
    "Spectrum",
    "Yield",
    "calibrate",
    "chart",
    "ecss_copper",
    "effective_yield",
    "facing_depth",
    "facing_fd",
    "facing_volts",
    "flights",
    "focusing",
    "gap_depth",
    "growth_rate",
    "impact_bound",
    "impact_energy",
    "level1",
    "margin_db",
    "resonant_impact",
    "resonant_volts",
    "speed_change_bound",
    "susceptibility",
    "threshold",
    "track",
    "volts_scale",
]
