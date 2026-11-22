"""The figures of docs/BEAM_LIFECYCLE.md: the life of a beam in array16's ring guide, from the source on the board,
round the ring as it is accelerated, to the detector at the end of the exit guide, and what it measures on the way.

The document answers the owner's questions of 2026-11-21: how electrons get into the guide, what happens as they
accelerate, which fields they make, whether they are relativistic, which resonances they meet, what ejecting and
sensing them means, and what they can tell us. This script computes its figures (section of the document; figure
number):

  ring3d      the 932 mm ring on two facing boards in 3-D: the LED, its pad and feeder, the gate at the junction,
              four gaps on one chirp, four bunches of the h = 4 ramp in flight, split pickups, a Smith-Purcell comb,
              the ejection junction, exit guide, retarding barrier and Faraday pad; beside it the beam's field lines
              across the guide ending on the four quarter electrodes (2, 3, 8; figure 1)
  sources     electrons per second against 265 nm power for each pad, the one-turn fills and the filament; the
              transverse energy of a held electron against time under gas collisions, beside the cooling times
              (3.2 to 3.4; 2)
  phase       the h = 4 ramp's longitudinal phase space in 3-D: phase, energy offset and lap, the separatrix a
              surface widening as W^1/2 round 4,000 tracked electrons (4; 3)
  speed       beta against energy; the travelling wave's factor 1/(1 -+ beta n_eff); the beam's |E| and |B|; the
              current a passing bunch induces on a quarter electrode (5, 6; 4)
  resonances  the board orotron's start current; Smith-Purcell frequencies; the transverse tune, 2 nu mod 1 and the
              action's growth at a joint along the ramp (7; 5)
  world       displacement by stray E and B; hold time in 1e-5 mbar of argon; Sagnac and Aharonov-Bohm phases (9; 6)

    uv run python experiments/beam_lifecycle_study.py              # every figure, a few seconds
    uv run python experiments/beam_lifecycle_study.py --quick      # fewer electrons, coarser grids: .../quick/
    uv run python experiments/beam_lifecycle_study.py --only phase,speed

Models: the four-gap ring of experiments/accelerator_study.py (the 932 mm trefoil's length, four 0.5 V gaps on one
chirp, a continuous 50 eV beam, the motion along the path alone), run at harmonic 4 and 13; facing boards 2h = 1 mm
apart in the gapless model of array16.particles.electrodes (ELE-2), the beam's own field a line charge between
grounded planes; a collision's heating from the Floquet solutions of Mathieu's equation, averaged over the drive phase
and the scattering direction; a joint as a thin gradient error once per lap in the smooth-focusing approximation;
argon's cross-section from modified effective-range theory (O'Malley 1963) with Golden and Bandel's scattering length,
fitted to their Ramsauer-Townsend minimum. Outputs: beam_lifecycle/{*.png, results.json, run.log}.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from array16 import lab
from array16 import particles as pt
from array16.particles import Arr, BoolArr, CArr, ClosedCurve

OUT = Path(__file__).resolve().parent / "beam_lifecycle"
LOG: list[str] = []
QM = pt.E / pt.ME
REST_EV = pt.ME * pt.C0**2 / pt.E
HBAR = 1.054571817e-34
PLANCK = 2 * math.pi * HBAR
RING_L = ClosedCurve.torus_knot(0.060, 0.028, 600).length  # the four-gap model's 932 mm
SHORT_L = ClosedCurve.torus_knot(0.020, 0.006, 2001).length  # the first run's 277 mm guide
HALF_GAP = 0.5e-3  # facing boards 1 mm apart: h = 0.5 mm from the null to each board
GAS_T = 293.0  # K
SIGMA_GAS = 1e-19  # m^2, the repository's estimate for residual gas (docs/PROGRESS.md)

# Ink, surfaces and the categorical order of the reference palette (fixed order, never cycled).
INK, INK2, MUTED, GRID, EDGE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
BLUES = ("#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
SURFACE, COPPER, BOARD = "#fcfcfb", "#c8874a", "#7fa37a"


def say(line: str = "") -> None:
    print(line, flush=True)
    LOG.append(line)


def speed(w_ev: Any) -> Any:
    """Speed (m/s) of an electron of kinetic energy w (eV), relativistic; arrays or numbers."""
    return pt.C0 * np.sqrt(1 - 1 / (1 + np.asarray(w_ev, dtype=np.float64) / REST_EV) ** 2)


def beta_of(w_ev: Any) -> Any:
    return speed(w_ev) / pt.C0


def gas_density(p_mbar: float) -> float:
    """Molecules per cubic metre at p (mbar) and room temperature."""
    return p_mbar * 100 / (pt.KB * GAS_T)


# ---- the four-gap ramp ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Ramp:
    """A tracked ramp of the four-gap ring: the drive's chirp, what it captured, and every electron's phase offset
    from the synchronous phase (rad) and energy offset (eV) as it crosses the first gap on each lap. `laps` is the
    last lap sampled at the first gap; the ramp takes `turns` laps, its last gap crossing raising the synchronous
    energy to `final`."""

    harmonic: int
    length: float
    f_start: float
    f_end: float
    time: float
    laps: int
    synchronous: Arr  # (laps + 1,) eV at the first gap
    phase: Arr  # (laps + 1, count)
    offset: Arr  # (laps + 1, count)
    energy: Arr  # (count,) at the end
    captured: BoolArr
    phase_s: float
    gap_volts: float
    gaps: int
    turns: int
    final: float

    @property
    def captured_share(self) -> float:
        return float(self.captured.mean())

    @property
    def spread(self) -> float:
        return float(self.energy[self.captured].std())


def ramp(
    harmonic: int,
    count: int = 4000,
    start: float = 50.0,
    end: float = 185.0,
    gap_volts: float = 0.5,
    phase_s: float = math.radians(30),
    gaps: int = 4,
    seed: int = 1,
) -> Ramp:
    """The four-gap ring of experiments/accelerator_study.py at a chosen harmonic: one chirp drives four gaps a
    quarter lap apart (in phase at h = 4, 90 degrees apart at h = 13); the synchronous electron gains gap_volts
    sin(phase_s) at each, and the drive's phase is tabulated at its passages, so the chirp is the one that keeps it
    synchronous. A continuous beam of `count` electrons at start +- 1 eV is tracked by its arrival times."""
    length = RING_L
    f0 = harmonic * float(speed(start)) / length
    times, sync = [-(length / gaps) / float(speed(start)), 0.0], [start, start + gap_volts * math.sin(phase_s)]
    while sync[-1] < end - 1e-9:
        times.append(times[-1] + (length / gaps) / float(speed(sync[-1])))
        sync.append(sync[-1] + gap_volts * math.sin(phase_s))
    t_sync, w_sync = np.array(times), np.array(sync)
    table = (np.arange(len(t_sync)) - 1) * harmonic / gaps
    rng = np.random.default_rng(seed)
    t = rng.uniform(-0.5, 0.5, count) / f0
    w = start + rng.uniform(-1, 1, count)
    phases, offsets, synchronous = [], [], []
    for j in range(len(t_sync) - 1):
        seen = 2 * math.pi * (np.interp(t, t_sync, table) - j * harmonic / gaps) + phase_s
        if j % gaps == 0:
            phases.append(np.angle(np.exp(1j * (seen - phase_s))))
            offsets.append(w - w_sync[j])
            synchronous.append(w_sync[j])
        w = np.maximum(w + gap_volts * np.sin(seen), 0.5)
        t = t + (length / gaps) / speed(w)
    captured = np.abs(w - w_sync[-1]) < 0.05 * w_sync[-1]
    return Ramp(
        harmonic,
        length,
        f0,
        harmonic * float(speed(w_sync[-1])) / length,
        float(t_sync[-1]),
        len(synchronous) - 1,
        np.array(synchronous),
        np.array(phases),
        np.array(offsets),
        w,
        captured,
        phase_s,
        gap_volts,
        gaps,
        (len(t_sync) - 1) // gaps,
        float(w_sync[-1]),
    )


def separatrix(energy: float, harmonic: int, gap_volts: float, gaps: int, phase_s: float, phi: Arr) -> Arr:
    """The bucket's half-height (eV) at phase offsets phi (rad) from the synchronous phase, for a non-relativistic
    beam of energy W crossing n_g gaps of peak volts V per lap: from H = -(pi h / 2W) dW^2 + n_g e V (cos phi +
    phi sin phi_s), whose separatrix passes through the unstable fixed point pi - phi_s (McMillan); NaN outside."""
    absolute = phi + phase_s
    level = math.cos(math.pi - phase_s) + (math.pi - phase_s) * math.sin(phase_s)
    inside = np.cos(absolute) + absolute * math.sin(phase_s) - level
    height = np.sqrt(np.maximum(2 * energy * gaps * gap_volts * inside / (math.pi * harmonic), 0.0))
    out: Arr = np.where((inside >= 0) & (absolute <= math.pi - phase_s), height, np.nan)
    return out


def bucket_height(
    energy: float, harmonic: int, gap_volts: float = 0.5, gaps: int = 4, phase_s: float = math.radians(30)
) -> float:
    """The bucket's half-height (eV) at the synchronous phase: sqrt((2 W n_g e V / pi h)(2 cos phi_s - (pi - 2 phi_s)
    sin phi_s))."""
    return math.sqrt(
        2
        * energy
        * gaps
        * gap_volts
        / (math.pi * harmonic)
        * (2 * math.cos(phase_s) - (math.pi - 2 * phase_s) * math.sin(phase_s))
    )


def synchrotron_tune(
    energy: float, harmonic: int, gap_volts: float = 0.5, gaps: int = 4, phase_s: float = math.radians(30)
) -> float:
    """Q_s = sqrt(h n_g e V cos phi_s / (4 pi W)), oscillations per lap (Proposition: phase stability)."""
    return math.sqrt(harmonic * gaps * gap_volts * math.cos(phase_s) / (4 * math.pi * energy))


def synchrotron_frequency(
    harmonic: int, length: float = RING_L, gap_volts: float = 0.5, gaps: int = 4, phase_s: float = math.radians(30)
) -> float:
    """f_s = (1/L) sqrt(h n_g e V cos phi_s / (2 pi m)), the same at every energy (Hz)."""
    return math.sqrt(harmonic * gaps * gap_volts * pt.E * math.cos(phase_s) / (2 * math.pi * pt.ME)) / length


# ---- the ring and its parts ------------------------------------------------------------------------------------------
# Where the parts sit round the 932 mm ring, in degrees of azimuth (the beam runs counter-clockwise).
GAPS_DEG = (0.0, 90.0, 180.0, 270.0)
INJECT_DEG, PICKUP_DEG, COMB_DEG, EJECT_DEG = 135.0, (200.0, 215.5), (238.0, 257.3), 315.0


def line_charge_field_lines(h: float, levels: int) -> tuple[Arr, Arr, Arr, Arr]:
    """The beam's own field across the guide: a line charge midway between two grounded boards 2h apart. Its complex
    potential is -(lambda / 2 pi eps0) log tanh(pi z / 4h), |tanh| = 1 on both boards, so field lines are the level
    sets of arg tanh(pi z / 4h); returns the grid, the potential's level function and the field lines' function."""
    xs = np.linspace(-4 * h, 4 * h, 801)
    ys = np.linspace(-h, h, 241)
    gx, gy = np.meshgrid(xs, ys)
    w = np.tanh(math.pi * (gx + 1j * gy) / (4 * h))
    potential = -np.log(np.maximum(np.abs(w), 1e-300))
    flux = np.angle(w)
    return xs, ys, potential, flux


def quarter_shares(h: float) -> list[float]:
    """The share of a centred line charge's image on each quarter electrode, top-left, top-right, bottom-left,
    bottom-right: the change of arg tanh(pi z / 4h) along each, over 2 pi (INJ-2)."""
    far = 40 * h
    shares = []
    for y in (h, -h):
        for x0, x1 in ((-far, 0.0), (0.0, far)):
            xs = np.linspace(x0, x1, 4001)
            arg = np.unwrap(np.angle(np.tanh(math.pi * (xs + 1j * y) / (4 * h))))
            shares.append(float(abs(arg[-1] - arg[0]) / (2 * math.pi)))
    return shares


def compute_ring3d(quick: bool) -> dict[str, Any]:
    say("\n== ring3d: the 932 mm ring on two facing boards, and the beam's field across the guide (figure 1) ==")
    run = ramp(4, count=1500 if quick else 4000)
    lap = run.laps // 2
    radius = RING_L / (2 * math.pi)
    shares = quarter_shares(HALF_GAP)
    say(
        f"  ring {RING_L * 1e3:.1f} mm round, radius {radius * 1e3:.1f} mm; the h = 4 ramp drawn at lap {lap}, "
        f"{run.synchronous[lap]:.1f} eV synchronous: four bunches a quarter lap apart; the image of a centred line "
        "charge on the four quarters: " + ", ".join(f"{s:.4f}" for s in shares)
    )
    return {
        "radius_mm": radius * 1e3,
        "lap": lap,
        "energy_eV": float(run.synchronous[lap]),
        "quarter_shares": shares,
        "_ramp": run,
    }


# ---- getting electrons in --------------------------------------------------------------------------------------------
PHOTONS_PER_MW = 1e-3 / (PLANCK * pt.C0 / 265e-9)  # 265 nm photons per second in a milliwatt
PADS = {  # quantum efficiency at 265 nm (section 3.2)
    "magnesium as installed, 2 × 10⁻⁵": 2e-5,
    "copper, sputtered and annealed, 9 × 10⁻⁵": 9e-5,
    "silver, sputtered and annealed, 2 × 10⁻⁴": 2e-4,
    "magnesium laser-cleaned in vacuum, 2 × 10⁻³": 2e-3,
}
FILAMENT_A = {"2000 K": 1.4e-6, "2200 K": 18e-6}  # the tungsten filament of section 3.2, Richardson-Dushman
GAIN_Q = (0.1, 0.3, 0.5)  # the Mathieu q of section 3.4's heating band


def floquet(q: float, steps: int = 1000) -> tuple[float, CArr, CArr]:
    """Mathieu's equation x'' + 2 q cos(2 tau) x = 0 (a = 0, tau = w t / 2) over one period: the Floquet exponent beta
    and the Floquet solution phi = e^(i beta tau) P(tau), P of period pi, with phi', on `steps` points of the period.
    Two solutions are integrated by fourth-order Runge-Kutta; phi is the monodromy matrix's eigenvector of eigenvalue
    e^(i pi beta), scaled so that the mean of P, the secular harmonic, is 1."""
    h = math.pi / steps

    def slope(t: float, y: Arr) -> Arr:
        c = -2 * q * math.cos(2 * t)
        out: Arr = np.array([y[1], c * y[0], y[3], c * y[2]])
        return out

    y: Arr = np.array([1.0, 0.0, 0.0, 1.0])
    rows = [y]
    for n in range(steps):
        t = n * h
        k1 = slope(t, y)
        k2 = slope(t + h / 2, y + h / 2 * k1)
        k3 = slope(t + h / 2, y + h / 2 * k2)
        k4 = slope(t + h, y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        rows.append(y)
    path = np.array(rows)
    values, vectors = np.linalg.eig(np.array([[path[-1, 0], path[-1, 2]], [path[-1, 1], path[-1, 3]]]))
    k = int(np.argmax(np.angle(values)))
    beta = float(np.angle(values[k])) / math.pi
    phi: CArr = vectors[0, k] * path[:-1, 0] + vectors[1, k] * path[:-1, 2]
    dphi: CArr = vectors[0, k] * path[:-1, 1] + vectors[1, k] * path[:-1, 3]
    scale = complex(np.mean(phi * np.exp(-1j * beta * h * np.arange(steps))))
    return beta, phi / scale, dphi / scale


def collision_gains(q: float) -> dict[str, float]:
    """What a collision with a heavy scatterer at a random drive phase does to a held electron in exact Mathieu motion
    (section 3.4): the electron keeps its position and speed, its velocity is redirected isotropically in three
    dimensions, both transverse directions are held (the y direction sees the drive half a period later) and the
    motion along the guide is free. With x = Re[A phi] in each transverse direction the secular energy is
    (m/2) w_s^2 |A|^2, w_s = beta w / 2, and a collision sets the new A from the Wronskian of x and phi; averaged over
    the drive phase, the secular phases and the direction, the mean energies after a collision, transverse and along
    the guide, are linear in those before. Returns c = 2 <dE>/E for an electron at rest along the guide (2 to first
    order in q, the Proposition), c in the steady state of repeated collisions (from that linear map's leading
    eigenvalue), and the steady state's transverse share of the energy."""
    beta, phi, dphi = floquet(q)
    wronskian = abs(complex(np.mean(np.conj(phi) * dphi - np.conj(dphi) * phi)))
    size, rate = np.abs(phi) ** 2, np.abs(dphi) ** 2
    same = float(np.mean(size * rate))
    cross = float(np.mean(size * np.roll(rate, -len(rate) // 2)))
    # (transverse, along the guide) energies after against before, in units of m (w/2)^2
    step = np.array(
        [
            [
                4 / wronskian**2 * (same / 2 + (same + cross) / 6),
                8 * beta**2 * float(np.mean(size)) / (3 * wronskian**2),
            ],
            [float(np.mean(rate)) / (6 * beta**2), 1 / 3],
        ]
    )
    values, vectors = np.linalg.eig(step)
    k = int(np.argmax(values.real))
    mode = vectors[:, k].real / vectors[:, k].real.sum()
    return {
        "from_rest": 2 * (float(step[0, 0] + step[1, 0]) - 1),
        "steady": 2 * (float(values[k].real) - 1),
        "transverse_share": float(mode[0]),
    }


def fills() -> dict[str, dict[str, float]]:
    """What holding 1e7 electrons in one turn takes, for a gate that closes within a lap: I = N e / T_rev."""
    out = {}
    for name, length, energy in (
        ("277 mm, 50 eV", SHORT_L, 50.0),
        ("932 mm, 50 eV", RING_L, 50.0),
        ("932 mm, 1 eV", RING_L, 1.0),
    ):
        lap = length / float(speed(energy))
        current = 1e7 * pt.E / lap
        rate = current / pt.E
        out[name] = {
            "lap_s": lap,
            "current_A": current,
            "mW_magnesium": rate / (PADS["magnesium as installed, 2 × 10⁻⁵"] * PHOTONS_PER_MW),
            "mW_cleaned": rate / (PADS["magnesium laser-cleaned in vacuum, 2 × 10⁻³"] * PHOTONS_PER_MW),
        }
    return out


def heating(energy0: float, energy1: float, pressure: float, gain: float) -> float:
    """Time (s) for a held electron's transverse energy to grow from energy0 to energy1 (eV) by collisions with
    residual gas: each collision multiplies it by 1 + gain/2 on average (the mean gain is gain times the mean
    pseudo-potential, half the secular energy), at the rate n sigma v(E), so dE/dt = E n sigma v / kappa with
    kappa = 1 / ln(1 + gain/2) collisions per e-fold; integrated, 2 kappa (tau_c(E0) - tau_c(E1))."""
    kappa = 1 / math.log(1 + gain / 2)
    n = gas_density(pressure)
    tau = [1 / (n * SIGMA_GAS * float(speed(e))) for e in (energy0, energy1)]
    return 2 * kappa * (tau[0] - tau[1])


def heating_curve(energy0: float, pressure: float, gain: float, times: Arr) -> Arr:
    """E(t) for the same law: 1/v(E) = 1/v0 - n sigma t / (2 kappa), until it diverges."""
    kappa = 1 / math.log(1 + gain / 2)
    n = gas_density(pressure)
    inverse = 1 / float(speed(energy0)) - n * SIGMA_GAS * times / (2 * kappa)
    v = np.where(inverse > 0, 1 / np.maximum(inverse, 1e-300), np.inf)
    out: Arr = 0.5 * pt.ME * v**2 / pt.E
    return out


def cooling_times() -> dict[str, tuple[float, float]]:
    """Resistive cooling (Wineland and Dehmelt): a pickup of resistance R at D = 1 mm damps the secular energy at
    e^2 R / (m D^2); a board tank at 106 MHz of 1 to 3 pF and Q 50 to 150 is 25 to 225 kOhm. Round a lap, with
    50 pF per metre of guide, tau = m D^2 w C' L / (e^2 Q)."""
    omega, distance = 2 * math.pi * pt.secular_frequency(0.3, 1e9), 1e-3
    resist = (50 / (omega * 3e-12), 150 / (omega * 1e-12))
    segment = tuple(pt.ME * distance**2 / (pt.E**2 * r) for r in resist)
    lap = tuple(pt.ME * distance**2 * omega * 50e-12 * RING_L / (pt.E**2 * q) for q in (150, 50))
    return {"segment": (min(segment), max(segment)), "lap": (min(lap), max(lap)), "ohms": resist}


def compute_sources(quick: bool) -> dict[str, Any]:
    say("\n== sources: getting electrons in, and what gas does to them (figure 2) ==")
    rows = fills()
    say(
        f"  265 nm: {PHOTONS_PER_MW:.3g} photons per second per milliwatt; magnesium as installed gives "
        f"{2e-5 * PHOTONS_PER_MW * pt.E * 1e9:.2f} nA per milliwatt"
    )
    for name, row in rows.items():
        say(
            f"  1e7 electrons in one turn of the {name} ring: lap {row['lap_s'] * 1e9:.0f} ns, "
            f"{row['current_A'] * 1e6:.1f} uA; "
            f"{row['mW_magnesium']:.0f} mW of 265 nm light on magnesium as installed, {row['mW_cleaned']:.1f} mW "
            "cleaned"
        )
    cool = cooling_times()
    say(
        f"  resistive cooling at D = 1 mm: {cool['segment'][0] * 1e3:.2f} to {cool['segment'][1] * 1e3:.2f} ms in a "
        f"segment ({cool['ohms'][0] / 1e3:.0f} to {cool['ohms'][1] / 1e3:.0f} kOhm), {cool['lap'][0] * 1e3:.2f} to "
        f"{cool['lap'][1] * 1e3:.1f} ms round a lap"
    )
    gains = {q: collision_gains(q) for q in GAIN_Q}
    first = [row["from_rest"] for row in gains.values()]
    say(
        "  a collision in exact Mathieu motion, the velocity redirected isotropically: c = 2 <dE>/E "
        + ", ".join(f"{row['from_rest']:.3f}" for row in gains.values())
        + f" at q = {', '.join(f'{q:g}' for q in GAIN_Q)} from rest along the guide (e-fold every "
        + ", ".join(f"{1 / math.log(1 + c / 2):.2f}" for c in first)
        + " collisions); in the steady state "
        + ", ".join(f"{row['steady']:.3f}" for row in gains.values())
        + ", the energy "
        + ", ".join(f"{row['transverse_share'] * 100:.1f}" for row in gains.values())
        + " % transverse (e-fold every "
        + ", ".join(f"{1 / math.log(1 + row['steady'] / 2):.2f}" for row in gains.values())
        + " collisions)"
    )
    to_rim = {}
    for pressure in (1e-5, 1e-7, 1e-9):
        values = [heating(0.1, 1.0, pressure, g) for g in first]
        to_rim[f"{pressure:g}"] = (min(values), max(values))
        say(f"  {pressure:g} mbar: 0.1 eV to the 1 eV rim in {min(values) * 1e3:.3g} to {max(values) * 1e3:.3g} ms")
    crossover = {}
    for energy in (0.1, 1.0):
        kappa = 1 / math.log(1 + gains[0.1]["from_rest"] / 2)
        bounds = [kappa * pt.KB * GAS_T / (SIGMA_GAS * float(speed(energy)) * tau) / 100 for tau in cool["segment"]]
        crossover[f"{energy:g} eV"] = (min(bounds), max(bounds))
        say(
            f"  at {energy:g} eV resistive cooling overtakes the heating below {min(bounds):.1e} to "
            f"{max(bounds):.1e} mbar"
        )
    return {
        "photons_per_mW": PHOTONS_PER_MW,
        "fills": rows,
        "cooling_s": cool,
        "collision_gain": {f"q {q:g}": row for q, row in gains.items()},
        "to_rim_s": to_rim,
        "crossover_mbar": crossover,
        "_gains": first,
    }


# ---- acceleration --------------------------------------------------------------------------------------------------
def first_periods(run: Ramp) -> tuple[float, float]:
    """The quartiles of the captured electrons' first synchrotron periods in laps: the interval between their phase's
    first two upward crossings of the synchronous phase, each placed between laps by linear interpolation."""
    phases = run.phase[:, run.captured]
    periods = []
    for k in range(phases.shape[1]):
        x = phases[:, k]
        up = np.flatnonzero((x[:-1] < 0) & (x[1:] >= 0))[:2]
        if len(up) == 2:
            at = up - x[up] / (x[up + 1] - x[up])
            periods.append(float(at[1] - at[0]))
    low, high = np.percentile(periods, [25, 75])
    return float(low), float(high)


def compute_phase(quick: bool) -> dict[str, Any]:
    say("\n== phase: the h = 4 ramp's longitudinal phase space, and the h = 13 ramp's numbers (figure 3) ==")
    count = 1500 if quick else 4000
    four, thirteen = ramp(4, count), ramp(13, count)
    out: dict[str, Any] = {}
    for run in (four, thirteen):
        chirp = (run.f_end - run.f_start) / run.time
        heights = [bucket_height(e, run.harmonic) for e in (50.0, 185.0)]
        row = {
            "f_start_MHz": run.f_start / 1e6,
            "f_end_MHz": run.f_end / 1e6,
            "time_us": run.time * 1e6,
            "laps": run.turns,
            "captured": run.captured_share,
            "spread_eV": run.spread,
            "chirp_MHz_per_us": chirp / 1e12,
            "chirp_formula": run.harmonic * pt.E / (pt.ME * run.length**2) * 1.0 / 1e12,
            "tune_50": synchrotron_tune(50.0, run.harmonic),
            "tune_185": synchrotron_tune(185.0, run.harmonic),
            "f_s_MHz": synchrotron_frequency(run.harmonic) / 1e6,
            "bucket_50": heights[0],
            "bucket_185": heights[1],
        }
        kept = run.captured
        early = int(np.argmin(np.abs(run.synchronous - 58.1)))
        late = int(np.argmin(np.abs(run.synchronous - 175.1)))
        growth = float(run.offset[late, kept].std() / run.offset[early, kept].std())
        row["spread_growth"] = growth
        row["spread_expected"] = math.sqrt(run.synchronous[late] / run.synchronous[early])
        row["captured_count"] = int(kept.sum())
        settled = np.arange(len(run.synchronous)) >= 30
        spreads = run.offset[:, kept].std(axis=1)
        row["spread_slope"] = float(np.polyfit(np.log(run.synchronous[settled]), np.log(spreads[settled]), 1)[0])
        row["first_period_low"], row["first_period_high"] = first_periods(run)
        out[f"h = {run.harmonic}"] = row
        say(
            f"  h = {run.harmonic}: chirp {row['f_start_MHz']:.1f} -> {row['f_end_MHz']:.1f} MHz in "
            f"{row['time_us']:.1f} "
            f"us ({row['laps']} laps), {row['chirp_MHz_per_us']:.3f} MHz per us (h dW e / m L^2 = "
            f"{row['chirp_formula']:.3f}); {row['captured'] * 100:.1f} % captured at {row['spread_eV']:.2f} eV rms; "
            f"bucket +-{heights[0]:.2f} eV at 50 eV, +-{heights[1]:.2f} eV at 185 eV; Q_s {row['tune_50']:.3f} and "
            f"{row['tune_185']:.3f}, f_s {row['f_s_MHz']:.3f} MHz at every energy; the spread of the "
            f"{row['captured_count']} captured grows {growth:.2f} times from {run.synchronous[early]:.1f} to "
            f"{run.synchronous[late]:.1f} eV, against sqrt of the energies' ratio {row['spread_expected']:.2f}; "
            f"from lap 30 on it grows as W^{row['spread_slope']:.3f}; the captured electrons' first synchrotron "
            f"periods run {row['first_period_low']:.2f} to {row['first_period_high']:.2f} laps (quartiles), against "
            f"1/Q_s = {1 / row['tune_50']:.2f} at 50 eV"
        )
    path = thirteen.turns * thirteen.length
    free = {f"{p:g}": 1 / (gas_density(p) * SIGMA_GAS) for p in (1e-5, 1e-6, 1e-7)}
    survive = {k: math.exp(-path / v) for k, v in free.items()}
    say(
        f"  the ramp's {path:.1f} m of path against a free path of {free['1e-05']:.1f} m at 1e-5 mbar (sigma "
        f"{SIGMA_GAS:g} m^2, {GAS_T:.0f} K): unscattered "
        + ", ".join(f"{v * 100:.1f} % at {k} mbar" for k, v in survive.items())
    )
    out["path_m"] = path
    out["free_path_m"] = free
    out["unscattered"] = survive
    out["_four"] = four
    return out


# ---- the beam's fields and its speed -------------------------------------------------------------------------------
def doppler_factors(energy: float, n_eff: float) -> tuple[float, float]:
    """How the null's Mathieu q of a travelling-wave well changes for an electron riding with the wave and against
    it: q / (1 - beta n_eff) and q / (1 + beta n_eff) (Proposition ACC-3)."""
    b = float(beta_of(energy))
    return 1 / (1 - b * n_eff) if b * n_eff < 1 else math.inf, 1 / (1 + b * n_eff)


def strip_kernel_integral(t: Arr, h: float) -> Arr:
    """The midplane potential of two grounded boards 2h apart under a unit potential on the top board from -inf to
    t: (1/pi)(atan(tanh(pi t / 4h)) + pi/4), the integral of the strip's Poisson kernel 1/(4h cosh(pi t / 2h))."""
    out: Arr = (np.arctan(np.tanh(math.pi * t / (4 * h))) + math.pi / 4) / math.pi
    return out


def induced_current(count: float, energy: float, bunch: float, pickup: float, h: float, samples: int) -> dict[str, Arr]:
    """The charge and current a uniform bunch of `count` electrons, `bunch` metres long at `energy` eV, induces on one
    quarter of a pickup segment `pickup` metres long on facing boards 2h apart (Shockley-Ramo, i = q v E_w): the
    quarter's weighting potential on the axis is half the strip's for the whole width of the top board."""
    v = float(speed(energy))
    span = pickup + 2 * bunch + 10 * h
    times = np.linspace(-span / v, span / v, samples)
    along = np.linspace(-bunch / 2, bunch / 2, 401)

    def weight(z: Arr) -> Arr:
        out: Arr = 0.5 * (strip_kernel_integral(z + pickup / 2, h) - strip_kernel_integral(z - pickup / 2, h))
        return out

    charge = np.array([count * pt.E * float(weight(along + v * t).mean()) for t in times])
    current = np.gradient(charge, times)
    return {"times": times, "charge": charge, "current": current}


def compute_speed(quick: bool) -> dict[str, Any]:
    say("\n== speed: relativity, the travelling wave, the beam's fields and its image current (figure 4) ==")
    energies = (1.0, 10.0, 100.0, 185.0, 1e3, 1e4, 79e3)
    betas = {f"{e:g} eV": float(beta_of(e)) for e in energies}
    say("  beta: " + ", ".join(f"{k} {v:.4g}" for k, v in betas.items()))
    factors = {f"{e:g} eV": {f"n_eff {n:g}": doppler_factors(e, n) for n in (1.0, 2.0)} for e in (185.0, 1e4, 79e3)}
    say(
        "  riding with the wave q grows "
        + "; ".join(f"{k}: " + ", ".join(f"{n} {v[0]:.3g}" for n, v in row.items()) for k, row in factors.items())
    )
    count, energy = 1e7, 185.0
    current = count * pt.E * float(speed(energy)) / RING_L
    line = count * pt.E / RING_L
    e_half = line / (2 * math.pi * pt.EPS0 * HALF_GAP)
    b_half = pt.MU0 * current / (2 * math.pi * HALF_GAP)
    say(
        f"  1e7 electrons at 185 eV round {RING_L * 1e3:.0f} mm: {current * 1e6:.1f} uA, {e_half:.1f} V/m and "
        f"{b_half * 1e9:.2f} nT at 0.5 mm"
    )
    pulse = induced_current(1e5, energy, 22e-3, 40e-3, HALF_GAP, 1200 if quick else 4000)
    peak = float(pulse["current"].max())
    above = pulse["times"][pulse["current"] > peak / 2]
    width = float(above.max() - above.min()) if len(above) else math.nan
    rising = pulse["times"][(pulse["current"] > peak / 2) & (pulse["times"] < 0)]
    width = float(rising.max() - rising.min()) if len(rising) else width
    charge = float(pulse["charge"].max())
    say(
        f"  10^5 electrons at 185 eV in a 22 mm bunch through a 40 mm pickup: {charge * 1e15:.2f} fC on a quarter, a "
        f"{peak * 1e6:.2f} uA pulse {width * 1e9:.2f} ns wide, {peak * 50 * 1e6:.0f} uV on 50 ohm"
    )
    return {
        "beta": betas,
        "doppler": factors,
        "current_uA": current * 1e6,
        "E_at_half_mm": e_half,
        "B_at_half_mm_nT": b_half * 1e9,
        "pulse_fC": charge * 1e15,
        "pulse_uA": peak * 1e6,
        "pulse_ns": width * 1e9,
        "pulse_uV": peak * 50 * 1e6,
        "_pulse": pulse,
        "_line": line,
        "_current": current,
    }


# ---- resonances ------------------------------------------------------------------------------------------------------
def madey_constant() -> tuple[float, float]:
    """The best detuning and the cold-beam conductance constant of section 7.1: the largest value of
    -d/dpsi sinc^2(psi / 2) (Madey's gain function), times pi / 2, and where it is reached."""
    psi = np.linspace(0.5, 6.0, 200001)
    x = psi / 2
    gain = -0.5 * (np.sin(2 * x) / x**2 - 2 * np.sin(x) ** 2 / x**3)
    best = int(np.argmax(gain))
    return float(psi[best]), float(gain[best] * math.pi / 2)


def start_current(volts: float, periods: Any, reach: float, resistance: float, constant: float) -> Any:
    """The board orotron's start-oscillation current (Hypothesis RES-H1): V0 e^(4 pi h / d) / (4 G N^3 R), G the
    conductance constant; reach is h/d."""
    return (
        volts * math.exp(4 * math.pi * reach) / (4 * constant * np.asarray(periods, dtype=np.float64) ** 3 * resistance)
    )


def energy_for_speed(v: float) -> float:
    """Kinetic energy (eV) of an electron moving at v (m/s)."""
    return REST_EV * (1 / math.sqrt(1 - (v / pt.C0) ** 2) - 1)


def joint_growth(
    run: Ramp, secular: float, depth_error: float, joint: float, count: int, seed: int = 0
) -> tuple[Arr, Arr]:
    """The transverse tune nu = f_sec / f_rev on each lap of a ramp, at the synchronous energy with which the lap
    begins, and the mean transverse action over `count` electrons of random betatron phase, relative to its start,
    when a joint `joint` metres long with a gradient error `depth_error` kicks them once per lap, on every lap of the
    ramp: in the smooth-focusing frame, a turn by 2 pi nu, then P -= (2 delta l_j / beta) X with beta = v / w_sec
    (the focusing strength goes as the gradient squared)."""
    rng = np.random.default_rng(seed)
    omega = 2 * math.pi * secular
    angle = rng.uniform(0, 2 * math.pi, count)
    x, p = np.sqrt(2) * np.cos(angle), np.sqrt(2) * np.sin(angle)
    tunes, actions = [], [1.0]
    for energy in run.synchronous:
        v = float(speed(energy))
        nu = secular * run.length / v
        kick = 2 * depth_error * joint * omega / v
        c, s = math.cos(2 * math.pi * nu), math.sin(2 * math.pi * nu)
        x, p = c * x + s * p, -s * x + c * p
        p = p - kick * x
        tunes.append(nu)
        actions.append(float(np.mean(0.5 * (x**2 + p**2))))
    return np.array(tunes), np.array(actions)


def compute_resonances(quick: bool) -> dict[str, Any]:
    say("\n== resonances: a board orotron, Smith-Purcell lines, and the tune along the ramp (figure 5) ==")
    psi, constant = madey_constant()
    say(
        f"  Madey: -d/dpsi sinc^2(psi/2) peaks at psi = {psi:.3f}; times pi/2, {constant:.4f}; 4 G = {4 * constant:.3f}"
    )
    cases = {}
    for periods, reach, ohms in ((50, 0.25, 1e4), (50, 0.25, 1e3), (50, 0.5, 1e3), (20, 0.5, 1e3)):
        cases[f"N {periods}, h/d {reach}, {ohms / 1e3:g} kOhm"] = float(
            start_current(100.0, periods, reach, ohms, constant)
        )
    ring_current = 1e7 * pt.E * float(speed(100.0)) / RING_L
    say(
        "  start currents at 100 eV over 1 mm: "
        + "; ".join(f"{k}: {sci(v)} A" for k, v in cases.items())
        + f"; 1e7 electrons round 932 mm at 100 eV carry {ring_current * 1e6:.1f} uA"
    )
    band = {}
    for period in (0.5e-3, 1e-3, 3e-3, 10e-3):
        band[f"{period * 1e3:g} mm"] = (energy_for_speed(0.4e9 * period), energy_for_speed(6e9 * period))
    say(
        "  the 0.4-6 GHz band at normal emission: "
        + "; ".join(f"d {k}: {v[0]:.3g} to {v[1]:.3g} eV" for k, v in band.items())
        + f"; 100 eV over 1 mm at {float(speed(100.0)) / 1e-3 / 1e9:.2f} GHz"
    )
    run = ramp(4, count=200)
    guides = {
        "1 GHz guide, 106 MHz secular": pt.secular_frequency(0.3, 1e9),
        "5.75 GHz facing boards, 610 MHz": pt.secular_frequency(0.3, 5.75e9),
    }
    count = 1000 if quick else 4000
    tunes, ends, growth = {}, {}, {}
    for name, secular in guides.items():
        nu, _ = joint_growth(run, secular, 0.0, 1e-3, 10)
        tunes[name] = nu
        ends[name] = (float(nu[0]), secular * run.length / float(speed(run.final)))
        crossings = int(abs(math.floor(2 * ends[name][0]) - math.floor(2 * ends[name][1])))
        for delta in (0.1, 0.05):
            _, action = joint_growth(run, secular, delta, 1e-3, count)
            growth[f"{name}, {delta * 100:g} %"] = action
        rate = np.diff(2 * nu)
        lingering = [int(i) for i in np.flatnonzero(np.abs(rate - np.round(rate)) < 0.05) if round(float(rate[i])) != 0]
        say(
            f"  {name}: over the {run.turns} laps of the ramp nu {ends[name][0]:.1f} -> {ends[name][1]:.1f} (50 -> "
            f"{run.final:.0f} eV), {crossings} half-integer resonances crossed; d(2 nu)/dlap {rate[0]:.2f} -> "
            f"{rate[-1]:.2f}, near a whole number at laps {lingering[:3]}...{lingering[-3:]}; mean action after the "
            f"ramp with one 1 mm joint at 10 % {growth[f'{name}, 10 %'][-1]:.2f} (over the last 12 laps "
            f"{growth[f'{name}, 10 %'][-12:].min():.2f} to {growth[f'{name}, 10 %'][-12:].max():.2f}), at 5 % "
            f"{growth[f'{name}, 5 %'][-1]:.3f} times"
        )
    secular = guides["5.75 GHz facing boards, 610 MHz"]
    stopband = np.array([2 * 0.1 * 1e-3 * 2 * math.pi * secular / (math.pi * float(speed(w))) for w in run.synchronous])
    return {
        "madey_psi": psi,
        "madey_constant": constant,
        "start_A": cases,
        "ring_current_A": ring_current,
        "band_eV": band,
        "tune_range": ends,
        "final_growth": {k: float(v[-1]) for k, v in growth.items()},
        "_tunes": tunes,
        "_stopband": stopband,
        "_growth": growth,
        "_constant": constant,
    }


# ---- what we can learn ---------------------------------------------------------------------------------------------
BOHR = 5.29177210903e-11
RYDBERG_EV = 13.605693122994
ARGON = {"A": -1.65, "alpha": 11.08}  # Golden and Bandel's scattering length (a0); argon's polarisability (a0^3)
ARGON_MINIMUM = (0.285, 0.125)  # eV, pi a0^2: Golden and Bandel's Ramsauer-Townsend minimum


def argon_cross_section(energy: Any, d_term: float, p_term: float) -> Arr:
    """e-Ar total cross-section (m^2) by modified effective-range theory (O'Malley 1963): tan eta_0 = -A k [1 +
    (4 alpha / 3) k^2 ln k] - (pi alpha / 3) k^2 + D k^3, tan eta_1 = pi alpha k^2 / 15 + A_1 k^3, and the
    polarisation phase pi alpha k^2 / ((2l + 3)(2l + 1)(2l - 1)) for l >= 2, in atomic units."""
    k = np.sqrt(np.asarray(energy, dtype=np.float64) / RYDBERG_EV)
    a, alpha = ARGON["A"], ARGON["alpha"]
    tan0 = -a * k * (1 + 4 * alpha / 3 * k**2 * np.log(k)) - math.pi * alpha / 3 * k**2 + d_term * k**3
    total = np.sin(np.arctan(tan0)) ** 2
    for order in range(1, 12):
        tangent = math.pi * alpha * k**2 / ((2 * order + 3) * (2 * order + 1) * (2 * order - 1))
        if order == 1:
            tangent = tangent + p_term * k**3
        total = total + (2 * order + 1) * np.sin(np.arctan(tangent)) ** 2
    out: Arr = 4 * math.pi / k**2 * total * BOHR**2
    return out


def fit_argon() -> tuple[float, float, float]:
    """D and A_1 such that the cross-section's minimum sits at Golden and Bandel's 0.285 eV with their 0.125 pi a0^2:
    Newton's method on the slope there and the value, from D = 47.5, A_1 = -12.5 a0^3; returns D, A_1 and the
    largest relative residual."""
    target, value = ARGON_MINIMUM
    unit = math.pi * BOHR**2

    def conditions(params: Arr) -> Arr:
        h = 1e-4 * target
        lo, mid, hi = (float(argon_cross_section(e, *params)) / unit for e in (target - h, target, target + h))
        return np.array([(hi - lo) / (2 * h) * target / value, mid / value - 1])

    params = np.array([47.5, -12.5])
    for _ in range(30):
        f = conditions(params)
        jac = np.column_stack([(conditions(params + e) - f) / 1e-4 for e in np.eye(2) * 1e-4])
        params = params - np.linalg.solve(jac, f)
    return float(params[0]), float(params[1]), float(np.abs(conditions(params)).max())


def compute_world(quick: bool) -> dict[str, Any]:
    say("\n== world: stray fields, argon's Ramsauer-Townsend window, Sagnac and Aharonov-Bohm (figure 6) ==")
    guides = {"106 MHz": pt.secular_frequency(0.3, 1e9), "610 MHz": pt.secular_frequency(0.3, 5.75e9)}
    per_field = {k: QM / (2 * math.pi * f) ** 2 for k, f in guides.items()}  # m per V/m
    earth = 50e-6
    shifts = {
        "185 eV, 610 MHz": QM * float(speed(185.0)) * earth / (2 * math.pi * guides["610 MHz"]) ** 2,
        "1 eV, 106 MHz": QM * float(speed(1.0)) * earth / (2 * math.pi * guides["106 MHz"]) ** 2,
    }
    say(
        "  a stray field E moves a held electron by eE/(m w_sec^2): "
        + ", ".join(f"{v * 1e6:.3g} um per V/m at {k}" for k, v in per_field.items())
        + "; the Earth's 50 uT moves a beam by evB/(m w_sec^2): "
        + ", ".join(f"{v * 1e6:.2g} um for {k}" for k, v in shifts.items())
    )
    d_term, p_term, residual = fit_argon()
    energies = np.logspace(math.log10(0.03), 0, 300)
    sigma = argon_cross_section(energies, d_term, p_term)
    density = gas_density(1e-5)
    hold = 1 / (density * sigma * speed(energies))
    at = {
        f"{e:g} eV": float(1 / (density * argon_cross_section(e, d_term, p_term) * speed(e)))
        for e in (0.1, 0.25, 0.285, 1.0)
    }
    estimate_10 = 1 / (density * 2e-19 * float(speed(10.0)))
    best = int(np.argmax(hold))
    say(
        f"  argon by MERT, A = {ARGON['A']} a0, alpha = {ARGON['alpha']} a0^3, fitted D = {d_term:.2f} a0^3 and A_1 = "
        f"{p_term:.2f} a0^3 (residual {residual:.0e}): sigma {float(argon_cross_section(0.285, d_term, p_term)):.3g} "
        "m^2 "
        f"at 0.285 eV; in 1e-5 mbar an electron lives "
        + ", ".join(f"{v * 1e3:.3g} ms at {k}" for k, v in at.items())
        + f", longest at {energies[best]:.3f} eV; at 10 eV with the section's 2e-19 m^2, {estimate_10 * 1e6:.0f} us"
    )
    omega_earth = 7.2921159e-5 * math.sin(math.radians(40.7))
    sagnac_per_area = 2 * pt.ME * omega_earth / HBAR
    circle = RING_L**2 / (4 * math.pi)
    points = ClosedCurve.torus_knot(0.060, 0.028, 600).points
    trefoil = 0.5 * float(np.sum(points[:, 0] * np.roll(points[:, 1], -1) - np.roll(points[:, 0], -1) * points[:, 1]))
    small = SHORT_L**2 / (4 * math.pi)
    tesla_per_radian = HBAR / (pt.E * circle)
    quantum = HBAR * 2 * math.pi * guides["106 MHz"]
    say(
        f"  Sagnac at 40.7 N: {sagnac_per_area:.3f} rad per m^2; {sagnac_per_area * circle:.4f} rad round the 932 mm "
        "circle "
        f"({circle * 1e4:.0f} cm^2), {sagnac_per_area * trefoil:.4f} rad for the trefoil's projected "
        f"{trefoil * 1e4:.0f} "
        f"cm^2; Aharonov-Bohm: one radian per {tesla_per_radian * 1e15:.2f} fT through the circle; hbar w_sec at 106 "
        "MHz "
        f"{quantum / pt.E * 1e6:.2f} ueV ({quantum / pt.KB * 1e3:.2f} mK), {0.1 * pt.E / quantum:.2g} quanta in 0.1 eV"
    )
    return {
        "per_field_um": {k: v * 1e6 for k, v in per_field.items()},
        "earth_shift_um": {k: v * 1e6 for k, v in shifts.items()},
        "argon": {
            "D": d_term,
            "A1": p_term,
            "residual": residual,
            "hold_s": at,
            "longest_eV": float(energies[best]),
            "estimate_10eV_s": estimate_10,
        },
        "sagnac_rad_per_m2": sagnac_per_area,
        "sagnac_circle_rad": sagnac_per_area * circle,
        "sagnac_trefoil_rad": sagnac_per_area * trefoil,
        "areas_cm2": {"932 mm circle": circle * 1e4, "trefoil, projected": trefoil * 1e4, "277 mm circle": small * 1e4},
        "fT_per_radian": tesla_per_radian * 1e15,
        "quantum_ueV": quantum / pt.E * 1e6,
        "_argon": (energies, sigma, hold),
        "_guides": guides,
    }


# ---- figures ---------------------------------------------------------------------------------------------------------
SUPERSCRIPTS = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def sci(value: float, digits: int = 2) -> str:
    """A number for a title: plain between 0.001 and 9999, else m × 10ⁿ."""
    if value == 0 or not math.isfinite(value):
        return f"{value:g}"
    exponent = math.floor(math.log10(abs(value)))
    if -3 <= exponent <= 3:
        return f"{value:.{max(digits - 1 - exponent, 0)}f}"
    mantissa = f"{value / 10**exponent:.{digits - 1}f}".rstrip("0").rstrip(".")
    power = "10" + str(exponent).translate(SUPERSCRIPTS)
    return power if mantissa == "1" else f"{mantissa} × {power}"


def _pyplot() -> Any:
    """matplotlib's pyplot on the Agg backend, in the house style of experiments/virtual_antenna_study.py."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from cycler import cycler

    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 9.5,
            "axes.labelsize": 9.5,
            "figure.titlesize": 11,
            "legend.fontsize": 8.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "axes.edgecolor": EDGE,
            "axes.labelcolor": INK,
            "axes.titlecolor": INK,
            "text.color": INK,
            "xtick.color": INK2,
            "ytick.color": INK2,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.axisbelow": True,
            "axes.titlelocation": "left",
            "lines.linewidth": 1.8,
            "lines.markersize": 5,
            "legend.frameon": False,
            "axes.prop_cycle": cycler(color=list(SERIES)),
            "savefig.dpi": 140,
            "mathtext.default": "regular",
        }
    )
    return plt


def _ramp_colours() -> Any:
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("energy", ["#6da7ec", "#256abf", "#6250d6", "#e87ba4", "#eb6834"])


def _suptitle(fig: Any, text: str, width: int = 150) -> None:
    """A figure's title, made from its data, wrapped so that it never sets the figure's width."""
    import re
    import textwrap

    glued = re.sub(r"(\d) (?=\S)", "\\1\u00a0", text)
    fig.suptitle(textwrap.fill(glued, width), x=0.01, ha="left", fontsize=10.5)


def _save(fig: Any, out: Path, name: str) -> None:
    import matplotlib.pyplot as plt

    fig.savefig(out / f"{name}.png", bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    say(f"  wrote {name}.png")


def _annulus(
    ax: Any, r0: float, r1: float, z: float, colour: str, alpha: float, arc: tuple[float, float] = (0.0, 360.0)
) -> None:
    theta = np.radians(np.linspace(arc[0], arc[1], 241))
    radii = np.array([r0, r1])
    t, r = np.meshgrid(theta, radii)
    ax.plot_surface(r * np.cos(t), r * np.sin(t), np.full_like(t, z), color=colour, alpha=alpha, lw=0, shade=False)


def _strip(
    ax: Any, start: Arr, direction: Arr, length: float, half: float, z: float, colour: str, alpha: float
) -> None:
    across = np.array([-direction[1], direction[0], 0.0])
    s = np.array([0.0, length])
    w = np.array([-half, half])
    gs, gw = np.meshgrid(s, w)
    x = start[0] + gs * direction[0] + gw * across[0]
    y = start[1] + gs * direction[1] + gw * across[1]
    ax.plot_surface(x, y, np.full_like(x, z), color=colour, alpha=alpha, lw=0, shade=False)
    for side in (-half, half):
        edge = start[:2, None] + np.outer(direction[:2], s) + side * across[:2, None]
        ax.plot(edge[0], edge[1], [z, z], color=EDGE, lw=0.8)


def figure_ring3d(res: dict[str, Any], out: Path) -> None:
    """Figure 1: the 932 mm ring on two facing boards in 3-D, and the beam's field lines across the guide."""
    plt = _pyplot()
    fig = plt.figure(figsize=(17.0, 8.2), layout="constrained")
    grid = fig.add_gridspec(1, 2, width_ratios=[1.75, 1.0])
    ax: Any = fig.add_subplot(grid[0, 0], projection="3d")
    radius = res["radius_mm"]
    scale = 7.0  # vertical scale
    top, bottom = HALF_GAP * 1e3 * scale, -HALF_GAP * 1e3 * scale
    half = 6.0  # each board half, mm
    plus, minus = "#cde2fb", "#fbe3d6"
    for z, alpha in ((bottom, 0.55), (top, 0.22)):
        inner, outer = (plus, minus) if z > 0 else (minus, plus)
        _annulus(ax, radius - half, radius - 0.15, z, inner, alpha)
        _annulus(ax, radius + 0.15, radius + half, z, outer, alpha)
    theta = np.radians(np.linspace(0, 360, 721))
    ax.plot(radius * np.cos(theta), radius * np.sin(theta), 0 * theta, color=INK2, lw=0.6, ls=(0, (2, 3)))
    for k, angle in enumerate(GAPS_DEG):
        a = math.radians(angle)
        for z in (top, bottom):
            ax.plot(
                [(radius - half) * math.cos(a), (radius + half) * math.cos(a)],
                [(radius - half) * math.sin(a), (radius + half) * math.sin(a)],
                [z, z],
                color=INK,
                lw=2.6,
            )
        ax.text(
            (radius + 22) * math.cos(a), (radius + 22) * math.sin(a), top + 4, f"gap {k + 1}", fontsize=8.5, ha="center"
        )
    run: Ramp = res["_ramp"]
    lap = res["lap"]
    bucket = bucket_height(float(run.synchronous[lap]), run.harmonic)
    keep = np.abs(run.offset[lap]) < 1.5 * bucket
    shift = -run.phase[lap][keep] / (2 * math.pi) * (RING_L * 1e3 / run.harmonic)
    energies = run.synchronous[lap] + run.offset[lap][keep]
    rng = np.random.default_rng(4)
    norm = plt.Normalize(float(energies.min()), float(energies.max()))
    for k in range(run.harmonic):
        centre = math.radians(52.0 + 360.0 * k / run.harmonic)
        angle = centre + shift / radius
        r = radius + rng.normal(0, 0.25, len(angle))
        lift = rng.normal(0, 0.12, len(angle)) * scale
        ax.scatter(
            r * np.cos(angle),
            r * np.sin(angle),
            lift,
            c=energies,
            cmap=_ramp_colours(),
            norm=norm,
            s=3.0,
            depthshade=False,
        )
    junction = math.radians(INJECT_DEG)
    along = np.array([-math.sin(junction), math.cos(junction), 0.0])
    joint = np.array([radius * math.cos(junction), radius * math.sin(junction), 0.0])
    feeder = joint - 80 * along
    for z, alpha in ((bottom, 0.55), (top, 0.22)):
        _strip(ax, feeder, along, 80.0, half, z, "#d9d8d0", alpha)
    ax.plot(*np.c_[feeder, joint], color=INK2, lw=0.8, ls=(0, (2, 3)))
    gate = joint - 18 * along
    across = np.array([-along[1], along[0], 0.0])
    ax.plot(*np.c_[gate - half * across, gate + half * across], color=SERIES[7], lw=4.0)
    ax.text(*(gate + np.array([0, 0, top + 4])), "gate", fontsize=8, color=SERIES[7], ha="center")
    pad = feeder + 4 * along
    corners = np.array(
        [pad + sx * 3 * along + sy * 3 * across for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1), (-1, -1))]
    )
    ax.plot(corners[:, 0], corners[:, 1], bottom + 0.05, color="#898781", lw=2.2)
    led = pad + np.array([0, 0, 30.0])
    ax.scatter(*led, marker="v", s=110, color=SERIES[6], depthshade=False)
    for corner in corners[:4]:
        ax.plot([led[0], corner[0]], [led[1], corner[1]], [led[2], bottom], color=SERIES[6], lw=0.8, alpha=0.7)
    ax.text(*(led + np.array([-10, 0, 6])), "265 nm LED over a\nmagnesium pad; the feeder", fontsize=8, ha="right")
    ax.text(*(joint + np.array([-8, 10, top + 8])), "junction", fontsize=8)
    lo, hi = (math.radians(v) for v in PICKUP_DEG)
    arc = np.linspace(lo, hi, 60)
    for z in (top, bottom):
        for r in (radius - 3, radius + 3):
            ax.plot(r * np.cos(arc), r * np.sin(arc), z + 0 * arc, color=SERIES[2], lw=3.0)
    mid = 0.5 * (lo + hi)
    ax.text(
        (radius + 40) * math.cos(mid),
        (radius + 40) * math.sin(mid),
        top + 6,
        "split pickups:\nfour quarters, 40 mm",
        fontsize=8,
        ha="center",
    )
    lo, hi = (math.radians(v) for v in COMB_DEG)
    for a in np.linspace(lo, hi, 50):
        ax.plot(
            [(radius - 3) * math.cos(a), (radius + 3) * math.cos(a)],
            [(radius - 3) * math.sin(a), (radius + 3) * math.sin(a)],
            [top, top],
            color=SERIES[6],
            lw=0.8,
        )
    mid = 0.5 * (lo + hi)
    ax.text(
        (radius - 34) * math.cos(mid),
        (radius - 34) * math.sin(mid),
        top + 2,
        "Smith-Purcell comb:\n50 periods of 1 mm",
        fontsize=8,
        ha="center",
    )
    exit_at = math.radians(EJECT_DEG)
    forward = np.array([-math.sin(exit_at), math.cos(exit_at), 0.0])
    start = np.array([radius * math.cos(exit_at), radius * math.sin(exit_at), 0.0])
    outward = np.array([math.cos(exit_at), math.sin(exit_at), 0.0])
    direction = forward * math.cos(0.35) + outward * math.sin(0.35)
    for z, alpha in ((bottom, 0.55), (top, 0.22)):
        _strip(ax, start, direction, 95.0, half, z, "#d9d8d0", alpha)
    ax.plot(*np.c_[start, start + 95 * direction], color=INK2, lw=0.8, ls=(0, (2, 3)))
    side = np.array([-direction[1], direction[0], 0.0])
    barrier = start + 62 * direction
    ax.plot(*np.c_[barrier - half * side, barrier + half * side], color=SERIES[3], lw=4.0)
    cup = start + 92 * direction
    square = np.array(
        [cup + sx * 4 * direction + sy * 4 * side for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1), (-1, -1))]
    )
    ax.plot(square[:, 0], square[:, 1], bottom + 0.05, color=SERIES[7], lw=2.6)
    ax.text(*(barrier + np.array([0, 0, top + 5])), "retarding\nbarrier", fontsize=8, ha="center")
    ax.text(*(cup + np.array([12, 0, top + 2])), "Faraday pad", fontsize=8)
    ax.text(*(start + np.array([30, -38, top + 2])), "ejection junction,\nexit guide", fontsize=8, ha="center")
    ax.set_xlim(-radius - 30, radius + 30)
    ax.set_ylim(-radius - 30, radius + 30)
    ax.set_zlim(-12, 40)
    ax.set_box_aspect((2 * radius + 60, 2 * radius + 60, 52), zoom=1.45)
    ax.set_axis_off()
    ax.view_init(elev=40, azim=-62)
    ax.text2D(
        0.02,
        0.95,
        f"ring {RING_L * 1e3:.0f} mm round (radius {radius:.1f} mm); facing boards 1 mm apart, drawn 7 mm apart\n"
        "(vertical scale ×7), each split along the path into halves at ±V/2 (blue, orange)",
        transform=ax.transAxes,
        fontsize=8.6,
    )
    bar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=_ramp_colours()), ax=ax, shrink=0.45, pad=0.0, location="bottom"
    )
    bar.set_label(f"electron energy (eV), the h = 4 ramp at lap {lap}: four bunches a quarter lap apart")
    side_ax = fig.add_subplot(grid[0, 1])
    h = HALF_GAP
    xs, ys, potential, flux = line_charge_field_lines(h, 24)
    side_ax.contour(
        xs * 1e3,
        ys * 1e3,
        potential,
        levels=np.linspace(0.05, 3.0, 10),
        colors=[BLUES[3]],
        linewidths=0.6,
        linestyles="solid",
    )
    upper, lower = ys > 0, ys < 0
    side_ax.contour(
        xs * 1e3,
        ys[upper] * 1e3,
        flux[upper],
        levels=np.linspace(0, math.pi, 21)[1:-1],
        colors=[INK2],
        linewidths=0.8,
        linestyles="solid",
    )
    side_ax.contour(
        xs * 1e3,
        ys[lower] * 1e3,
        flux[lower],
        levels=np.linspace(-math.pi, 0, 21)[1:-1],
        colors=[INK2],
        linewidths=0.8,
        linestyles="solid",
    )
    shares = res["quarter_shares"]
    labels = (("top, x < 0", -1, 1), ("top, x > 0", 1, 1), ("bottom, x < 0", -1, -1), ("bottom, x > 0", 1, -1))
    for (name, sx, sy), share in zip(labels, shares, strict=True):
        x0, x1 = (-2.0, -0.03) if sx < 0 else (0.03, 2.0)
        colour = plus if sx * sy < 0 else minus
        edge = SERIES[0] if sx * sy < 0 else SERIES[1]
        side_ax.fill_between([x0, x1], sy * 0.5, sy * 0.58, color=colour, ec=edge, lw=1.2)
        side_ax.text(
            0.5 * (x0 + x1),
            sy * 0.66,
            f"{name}: {'+' if sx * sy < 0 else '−'}V/2; {share:.4f} of the image",
            ha="center",
            va="center",
            fontsize=8,
        )
    side_ax.scatter([0], [0], s=70, color=SERIES[1], zorder=5, edgecolor=INK, lw=0.8)
    side_ax.set_xlim(-2, 2)
    side_ax.set_ylim(-0.75, 0.75)
    side_ax.set_aspect("equal")
    side_ax.set_xlabel("across the guide (mm)")
    side_ax.set_ylabel("between the boards (mm)")
    side_ax.set_title(
        "the cross-section, true scale: the beam on the RF null (dot); its field lines (grey) end on\nthe four "
        "quarters, "
        "a quarter of its image on each, so the quarters are its pickups (INJ-2);\nblue, its equipotentials",
        fontsize=9,
    )
    _suptitle(
        fig,
        "A three-dimensional model of the beam's life: electrons made at the pad merge through the gated junction, are "
        "bunched and accelerated by four gaps on one chirp, heard by the pickups and the comb, and leave through the "
        "ejection junction to the retarding barrier and the Faraday pad",
    )
    _save(fig, out, "ring3d")


def figure_sources(res: dict[str, Any], out: Path) -> None:
    """Figure 2: electrons per second against 265 nm power; transverse energy against time under gas collisions."""
    plt = _pyplot()
    fig, (left, right) = plt.subplots(1, 2, figsize=(15.5, 6.4), layout="constrained")
    power = np.logspace(-2, 3.5, 200)
    for (name, qe), colour in zip(PADS.items(), (SERIES[1], SERIES[3], SERIES[2], SERIES[0]), strict=True):
        left.loglog(power, qe * PHOTONS_PER_MW * power, color=colour, label=name)
    for (name, row), ls in zip(res["fills"].items(), ("--", "-.", ":"), strict=True):
        rate = row["current_A"] / pt.E
        left.axhline(rate, color=INK2, lw=1.0, ls=ls)
        left.text(
            1.2e-2, rate * 1.15, f"10⁷ in one turn, {name}: {row['current_A'] * 1e6:.1f} µA", fontsize=7.8, color=INK2
        )
        left.scatter([row["mW_magnesium"]], [rate], color=SERIES[1], s=28, zorder=4)
        left.scatter([row["mW_cleaned"]], [rate], color=SERIES[0], s=28, zorder=4)
    for (name, amps), place in zip(FILAMENT_A.items(), (0.62, 0.62), strict=True):
        left.axhline(amps / pt.E, color=SERIES[7], lw=1.0, alpha=0.8)
        left.text(
            2500,
            amps / pt.E * place,
            f"W filament, {name}: {amps * 1e6:.1f} µA",
            fontsize=7.8,
            color=SERIES[7],
            ha="right",
        )
    left.set_xlabel("265 nm optical power on the pad (mW)")
    left.set_ylabel("electrons per second")
    left.legend(loc="lower right", fontsize=8, title="quantum efficiency at 265 nm", title_fontsize=8)
    one = res["fills"]["932 mm, 1 eV"]
    left.set_title(
        f"getting electrons in: {sci(PHOTONS_PER_MW, 3)} photons per second per mW; 10⁷ electrons in one 1 eV turn "
        "of the 932 "
        f"mm ring\ntake {one['current_A'] * 1e6:.2f} µA, {one['mW_magnesium'] / 1e3:.2f} W on magnesium as installed "
        "or "
        f"{one['mW_cleaned']:.1f} mW cleaned (dots), or the 2000 K filament",
        fontsize=9,
    )
    times = np.logspace(-6, 1, 600)
    for pressure, colour in ((1e-5, SERIES[1]), (1e-7, SERIES[3]), (1e-9, SERIES[0])):
        curves = np.array([heating_curve(0.1, pressure, g, times) for g in res["_gains"]])
        low, high = curves.min(axis=0), curves.max(axis=0)
        right.fill_between(times, np.minimum(low, 3), np.minimum(high, 3), color=colour, alpha=0.35, lw=0)
        right.loglog(times, np.minimum(curves[1], 3), color=colour, lw=1.6, label=f"{pressure:g} mbar")
    cool = res["cooling_s"]
    right.axvspan(*cool["segment"], color=SERIES[2], alpha=0.15, lw=0)
    right.text(
        math.sqrt(cool["segment"][0] * cool["segment"][1]),
        0.012,
        "resistive\ncooling in a\nsegment",
        ha="center",
        fontsize=7.8,
        color=INK2,
    )
    right.axvspan(*cool["lap"], color=SERIES[6], alpha=0.12, lw=0)
    right.text(math.sqrt(cool["lap"][0] * cool["lap"][1]), 0.012, "round\na lap", ha="center", fontsize=7.8, color=INK2)
    right.axhline(1.0, color=INK, lw=1.0, ls="--")
    right.text(1.3e-6, 1.1, "the 1 eV rim", fontsize=8)
    right.axhline(pt.KB * 300 / pt.E, color=INK2, lw=0.9, ls=":")
    right.text(1.3e-6, pt.KB * 300 / pt.E * 1.1, "300 K, where cooling ends", fontsize=8, color=INK2)
    right.set_xlim(1e-6, 10)
    right.set_ylim(0.01, 3)
    right.set_xlabel("time (s)")
    right.set_ylabel("transverse (secular) energy of a held electron (eV)")
    right.legend(loc="upper left", fontsize=8, title="band: q 0.1 to 0.5, line q 0.3", title_fontsize=8)
    rim = res["to_rim_s"]
    right.set_title(
        f"what gas does: each collision adds about the secular energy ({min(res['_gains']):.2f} to "
        f"{max(res['_gains']):.2f} U), so 0.1 eV reaches\nthe rim in "
        + ", ".join(f"{v[0] * 1e3:.2g}-{v[1] * 1e3:.2g} ms at {k} mbar" for k, v in list(rim.items())[:2])
        + f", {rim['1e-09'][0]:.2g}-{rim['1e-09'][1]:.2g} s at 1e-09 mbar",
        fontsize=9,
    )
    cross = res["crossover_mbar"]
    _suptitle(
        fig,
        "Getting electrons in and keeping them: a gate that closes within a lap holds I T_rev / e electrons; gas heats "
        f"what it holds, and resistive cooling overtakes it below {cross['0.1 eV'][0]:.0e} to {cross['0.1 eV'][1]:.0e} "
        f"mbar at 0.1 eV and {cross['1 eV'][0]:.0e} to {cross['1 eV'][1]:.0e} mbar at 1 eV",
    )
    _save(fig, out, "sources")


def figure_phase(res: dict[str, Any], out: Path) -> None:
    """Figure 3: the h = 4 ramp's phase space in 3-D, the separatrix widening as W^1/2."""
    plt = _pyplot()
    run: Ramp = res["_four"]
    fig = plt.figure(figsize=(16.0, 7.6), layout="constrained")
    grid = fig.add_gridspec(1, 2, width_ratios=[1.6, 1.0])
    ax: Any = fig.add_subplot(grid[0, 0], projection="3d")
    laps = np.arange(run.laps + 1)
    phi = np.radians(np.linspace(-75, 125, 160))
    heights = np.array(
        [separatrix(float(w), run.harmonic, run.gap_volts, run.gaps, run.phase_s, phi) for w in run.synchronous]
    )
    gl, gp = np.meshgrid(laps, np.degrees(phi), indexing="ij")
    for sign in (1.0, -1.0):
        ax.plot_surface(gl, gp, sign * heights, color=SERIES[0], alpha=0.12, lw=0, shade=False)
    for at in range(0, run.laps + 1, 15):
        ring = np.r_[heights[at], -heights[at][::-1]]
        ax.plot(at + 0 * ring, np.r_[np.degrees(phi), np.degrees(phi)[::-1]], ring, color=SERIES[0], lw=1.0)
    kept = np.flatnonzero(run.captured)
    rng = np.random.default_rng(3)
    shown = rng.choice(kept, size=min(450, len(kept)), replace=False)
    every = 9
    rows = np.arange(0, run.laps + 1, every)
    xs = np.repeat(rows[:, None], len(shown), axis=1)
    ys = np.degrees(run.phase[np.ix_(rows, shown)])
    zs = run.offset[np.ix_(rows, shown)]
    inside = np.abs(zs) < 9
    ax.scatter(xs[inside], ys[inside], zs[inside], s=1.2, color=SERIES[1], alpha=0.7, depthshade=False)
    for k in shown[:3]:
        ax.plot(laps, np.degrees(run.phase[:, k]), run.offset[:, k], color=SERIES[7], lw=1.0)
    lost = np.flatnonzero(~run.captured)[:300]
    early = rows[rows < 30]
    lx, ly = np.degrees(run.phase[np.ix_(early, lost)]), run.offset[np.ix_(early, lost)]
    ok = np.abs(ly) < 9
    ax.scatter(
        np.repeat(early[:, None], len(lost), axis=1)[ok],
        lx[ok],
        ly[ok],
        s=0.5,
        color=MUTED,
        alpha=0.4,
        depthshade=False,
    )
    ax.set_ylim(-80, 130)
    ax.set_zlim(-9, 9)
    ax.set_xlabel("lap", fontsize=8)
    ax.set_ylabel("phase from the synchronous (deg)", fontsize=8)
    ax.set_zlabel("energy offset (eV)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_box_aspect((1.9, 1.0, 0.9))
    ax.view_init(elev=22, azim=-62)
    row = res["h = 4"]
    ax.set_title(
        f"h = 4: the chirp {row['f_start_MHz']:.1f} → {row['f_end_MHz']:.1f} MHz in {row['time_us']:.1f} µs; "
        f"{row['captured'] * 100:.0f} % of {len(run.energy)} electrons captured (orange, every\nninth lap; three "
        "traced in red; grey, "
        f"some of those lost), the separatrix (blue) ±{row['bucket_50']:.2f} eV at 50 eV and ±{row['bucket_185']:.2f} "
        "eV at 185 eV",
        fontsize=9,
    )
    side = fig.add_subplot(grid[0, 1])
    spread = run.offset[:, run.captured].std(axis=1)
    side.plot(run.synchronous, spread, color=SERIES[1], label="rms energy spread of the captured electrons")
    reference = spread[30] * np.sqrt(run.synchronous / run.synchronous[30])
    side.plot(run.synchronous, reference, color=INK2, ls="--", lw=1.2, label="∝ W^½ through lap 30")
    bucket = np.array([bucket_height(float(w), run.harmonic) for w in run.synchronous])
    side.plot(run.synchronous, bucket, color=SERIES[0], label="the bucket's half-height, ∝ W^½")
    side.set_xscale("log")
    side.set_yscale("log")
    side.set_xlabel("synchronous energy W (eV)")
    side.set_ylabel("eV")
    side.legend(loc="upper left", fontsize=8)
    thirteen = res["h = 13"]
    side.set_title(
        f"the invariant is ∮ΔW dt: from lap 30 on the captured electrons' spread grows as "
        f"W^{row['spread_slope']:.2f}\n(h = 13: W^{thirteen['spread_slope']:.2f}, {thirteen['captured_count']} "
        "captured), the bucket's half-height as W^½",
        fontsize=9,
    )
    _suptitle(
        fig,
        f"Longitudinal phase space along the h = 4 ramp: bunches turn about the synchronous phase at "
        f"{row['f_s_MHz']:.3f} MHz at every energy (h = 13: {thirteen['f_s_MHz']:.3f} MHz, Q_s "
        f"{thirteen['tune_50']:.3f} "
        f"at 50 eV and {thirteen['tune_185']:.3f} at 185 eV); the drive chirps linearly, {row['chirp_MHz_per_us']:.3f} "
        f"MHz per µs ({thirteen['chirp_MHz_per_us']:.2f} at h = 13), inside the HackRF's 20 MHz",
    )
    _save(fig, out, "phase")


def figure_speed(res: dict[str, Any], out: Path) -> None:
    """Figure 4: beta against energy, the travelling-wave factors, the beam's fields, the induced current."""
    plt = _pyplot()
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.0), layout="constrained")
    energy = np.logspace(-1, 6, 400)
    beta = beta_of(energy)
    a = axes[0, 0]
    a.loglog(energy, beta, color=SERIES[0])
    for key, value in res["beta"].items():
        e = float(key.split()[0])
        a.scatter([e], [value], color=SERIES[1], s=26, zorder=4)
        a.annotate(f"{key}: {value:.3g}", (e, value), xytext=(6, -12), textcoords="offset points", fontsize=7.8)
    a.set_xlabel("kinetic energy (eV)")
    a.set_ylabel("β = v / c")
    a.set_title(
        "relativity: β = 0.0020 at 1 eV, 0.027 at 185 eV, 0.195 at 10 keV; \u03b3 − 1 = W / 511 keV", fontsize=9
    )
    b = axes[0, 1]
    for n_eff, ls in ((1.0, "-"), (2.0, "--")):
        with np.errstate(divide="ignore"):
            riding = np.where(beta * n_eff < 1, 1 / (1 - beta * n_eff), np.nan)
        b.semilogx(energy, riding, color=SERIES[1], ls=ls, label=f"riding with the wave, n_eff = {n_eff:g}: 1/(1 − βn)")
        b.semilogx(
            energy, 1 / (1 + beta * n_eff), color=SERIES[0], ls=ls, label=f"against it, n_eff = {n_eff:g}: 1/(1 + βn)"
        )
    for key, row in res["doppler"].items():
        e = float(key.split()[0])
        b.scatter([e, e], [row["n_eff 1"][0], row["n_eff 2"][0]], color=INK, s=18, zorder=4)
    b.axhline(1, color=INK2, lw=0.8)
    b.set_ylim(0.3, 4)
    b.set_xlim(1, 1e6)
    b.set_xlabel("kinetic energy (eV)")
    b.set_ylabel("the null's Mathieu q over its value at rest")
    b.legend(loc="upper left", fontsize=8)
    d = res["doppler"]
    b.set_title(
        f"the travelling wave at speed: q grows {d['185 eV']['n_eff 1'][0]:.2f} and {d['185 eV']['n_eff 2'][0]:.2f} "
        "times "
        f"at 185 eV,\n{d['10000 eV']['n_eff 1'][0]:.2f} and {d['10000 eV']['n_eff 2'][0]:.2f} at 10 keV; at 79 keV it "
        "doubles in air and the electron rides unconfined on laminate",
        fontsize=9,
    )
    c = axes[1, 0]
    r = np.logspace(-4.7, math.log10(2.5e-3), 300)
    line, current = res["_line"], res["_current"]
    free = line / (2 * math.pi * pt.EPS0 * r)
    u = math.pi * r / (4 * HALF_GAP)
    boards = line / (2 * math.pi * pt.EPS0) * (math.pi / (2 * HALF_GAP)) / np.sinh(2 * u)
    c.loglog(r * 1e3, free, color=SERIES[0], label="|E|, the line charge alone")
    c.loglog(r * 1e3, boards, color=SERIES[0], ls="--", label="|E| on the midplane between grounded boards 1 mm apart")
    twin = c.twinx()
    twin.loglog(r * 1e3, pt.MU0 * current / (2 * math.pi * r) * 1e9, color=SERIES[1], label="|B| (right axis)")
    twin.set_ylabel("|B| (nT)")
    twin.spines["right"].set_visible(True)
    c.axvline(0.5, color=INK2, lw=0.8, ls=":")
    c.set_ylim(0.1, 3e4)
    c.set_xlabel("distance from the beam (mm)")
    c.set_ylabel("|E| (V/m)")
    handles = c.get_legend_handles_labels()[0] + twin.get_legend_handles_labels()[0]
    c.legend(handles=handles, loc="lower left", fontsize=8)
    c.set_title(
        f"the fields of 10⁷ electrons at 185 eV round 932 mm ({res['current_uA']:.1f} µA): {res['E_at_half_mm']:.0f} "
        "V/m "
        f"and {res['B_at_half_mm_nT']:.1f} nT\nat 0.5 mm, B/E = v/c²; between the boards the field dies within a gap",
        fontsize=9,
    )
    e = axes[1, 1]
    pulse = res["_pulse"]
    e.plot(pulse["times"] * 1e9, pulse["current"] * 1e6, color=SERIES[2], label="induced current")
    volt = e.twinx()
    volt.plot(
        pulse["times"] * 1e9, pulse["charge"] * 1e15, color=SERIES[6], lw=1.2, ls="--", label="induced charge (right)"
    )
    volt.set_ylabel("induced charge on the quarter (fC)")
    volt.spines["right"].set_visible(True)
    e.set_xlabel("time (ns), the bunch's centre at the pickup's centre at 0")
    e.set_ylabel("current into the quarter (µA); ×50 Ω for µV")
    handles = e.get_legend_handles_labels()[0] + volt.get_legend_handles_labels()[0]
    e.legend(handles=handles, loc="upper right", fontsize=8)
    e.set_title(
        f"10⁵ electrons at 185 eV in a 22 mm bunch passing a 40 mm pickup: {res['pulse_fC']:.2f} fC on a quarter, a\n"
        f"{res['pulse_uA']:.2f} µA pulse {res['pulse_ns']:.2f} ns wide as it enters and its mirror as it leaves, "
        f"{res['pulse_uV']:.0f} µV on 50 Ω (Shockley-Ramo)",
        fontsize=9,
    )
    _suptitle(
        fig,
        "Speed and fields: non-relativistic electrons whose speed shows first in the drive's Doppler shift, and "
        "whose own field the electrodes hear",
    )
    _save(fig, out, "speed")


def figure_resonances(res: dict[str, Any], out: Path) -> None:
    """Figure 5: the orotron's start current, Smith-Purcell frequencies, the tune and the action along the ramp."""
    plt = _pyplot()
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.4), layout="constrained")
    a = axes[0, 0]
    periods = np.arange(10, 201)
    constant = res["_constant"]
    for (reach, ohms), colour, ls in (
        ((0.25, 1e4), SERIES[0], "-"),
        ((0.25, 1e3), SERIES[0], "--"),
        ((0.5, 1e4), SERIES[1], "-"),
        ((0.5, 1e3), SERIES[1], "--"),
    ):
        a.semilogy(
            periods,
            start_current(100.0, periods, reach, ohms, constant) * 1e6,
            color=colour,
            ls=ls,
            label=f"h/d = {reach}, {ohms / 1e3:g} kΩ",
        )
    a.axhline(res["ring_current_A"] * 1e6, color=INK, lw=1.0, ls="-.")
    a.text(
        12,
        res["ring_current_A"] * 1e6 * 1.2,
        f"10⁷ electrons round 932 mm at 100 eV: {res['ring_current_A'] * 1e6:.1f} µA",
        fontsize=8,
    )
    for name, value in res["start_A"].items():
        n = int(name.split(",")[0].split()[1])
        a.scatter([n], [value * 1e6], color=INK, s=20, zorder=4)
    a.set_xlabel("periods N of the comb (d = 1 mm, 100 eV, so 5.93 GHz)")
    a.set_ylabel("start-oscillation current (µA)")
    a.legend(loc="upper right", fontsize=8)
    cases = list(res["start_A"].values())
    a.set_title(
        f"a board orotron (hypothesis): I_st = V₀ e^(4πh/d) / (4 G N³ R), G = {constant:.3f} (Madey, ψ = "
        f"{res['madey_psi']:.3f});\n{sci(cases[0] * 1e6)} µA at N = 50, h/d = 0.25, 10 kΩ, where a ring carries "
        f"{res['ring_current_A'] * 1e6:.0f} µA; {sci(cases[2] * 1e3)} mA at h/d = 0.5, 1 kΩ",
        fontsize=9,
    )
    b = axes[0, 1]
    energy = np.logspace(-1, 4.5, 300)
    v = speed(energy)
    b.axhspan(0.4, 6.0, color=SERIES[2], alpha=0.12, lw=0)
    b.text(0.13, 2.0, "the array's 0.4 to 6 GHz", fontsize=8, color=INK2)
    for period, colour in zip((0.5e-3, 1e-3, 3e-3, 10e-3), SERIES, strict=False):
        window = res["band_eV"][f"{period * 1e3:g} mm"]
        b.loglog(
            energy,
            v / period / 1e9,
            color=colour,
            label=f"d = {period * 1e3:g} mm: {window[0]:.3g} to {window[1]:.3g} eV in the band",
        )
    b.scatter([100], [float(speed(100.0)) / 1e-3 / 1e9], color=INK, s=24, zorder=4)
    b.set_xlabel("electron energy (eV)")
    b.set_ylabel("Smith-Purcell frequency at normal emission, v/d (GHz)")
    b.set_ylim(0.05, 300)
    b.legend(loc="upper left", fontsize=8)
    b.set_title(
        "Smith-Purcell lines: passing electrodes of period d, an electron emits at v/d; 100 eV over 1 mm\n"
        "is 5.93 GHz (dot); gaps, pickup segments and joints set d",
        fontsize=9,
    )
    c = axes[1, 0]
    tunes = res["_tunes"]
    name = "5.75 GHz facing boards, 610 MHz"
    nu = tunes[name]
    laps = np.arange(len(nu))
    width = res["_stopband"]
    c.fill_between(laps, 0, width, color=SERIES[7], alpha=0.18, lw=0, label="the stopband of one 1 mm joint at 10 %")
    c.fill_between(laps, 1 - width, 1, color=SERIES[7], alpha=0.18, lw=0)
    c.plot(laps, np.mod(2 * nu, 1), color=SERIES[1], lw=0.5, alpha=0.5)
    c.scatter(laps, np.mod(2 * nu, 1), color=SERIES[1], s=12, zorder=3, label="2\u03bd mod 1, 5.75 GHz facing boards")
    rate = np.diff(2 * nu)
    for target in (-2.0, -1.0):
        at = int(np.argmin(np.abs(rate - target)))
        c.axvline(at, color=INK2, lw=0.9, ls=":")
        c.text(at + 1, 1.02, f"d(2\u03bd)/dlap = {target:g}", fontsize=7.8, color=INK2)
    c.set_ylim(0, 1.1)
    c.set_xlabel("lap of the ramp (50 → 185 eV)")
    c.set_ylabel("2\u03bd mod 1")
    twin = c.twinx()
    twin.plot(laps, nu, color=SERIES[1], lw=1.2)
    twin.plot(laps, tunes["1 GHz guide, 106 MHz secular"], color=SERIES[0], lw=1.2, label="\u03bd, 1 GHz guide")
    twin.set_ylabel("tune \u03bd = f_sec / f_rev (lines)")
    twin.spines["right"].set_visible(True)
    c.legend(loc="lower left", fontsize=8)
    ends = res["tune_range"]
    c.set_title(
        f"the tune sweeps {ends[name][0]:.1f} → {ends[name][1]:.1f} (5.75 GHz) and "
        f"{ends['1 GHz guide, 106 MHz secular'][0]:.1f} → {ends['1 GHz guide, 106 MHz secular'][1]:.1f} (1 GHz); with "
        "one joint per lap only 2\u03bd mod 1\nmatters, and it lingers where its change per lap passes a whole number",
        fontsize=9,
    )
    d = axes[1, 1]
    growth = res["_growth"]
    styles = {
        "5.75 GHz facing boards, 610 MHz, 10 %": (SERIES[1], "-"),
        "5.75 GHz facing boards, 610 MHz, 5 %": (SERIES[1], "--"),
        "1 GHz guide, 106 MHz secular, 10 %": (SERIES[0], "-"),
        "1 GHz guide, 106 MHz secular, 5 %": (SERIES[0], "--"),
    }
    for key, values in growth.items():
        colour, ls = styles[key]
        d.plot(np.arange(len(values)), values, color=colour, ls=ls, label=f"{key}: ×{values[-1]:.3g}")
    for target in (-2.0, -1.0):
        d.axvline(int(np.argmin(np.abs(rate - target))), color=INK2, lw=0.9, ls=":")
    d.set_xlabel("lap of the ramp")
    d.set_ylabel("mean transverse action over its start")
    d.legend(loc="upper left", fontsize=8, title="one 1 mm joint, gradient error", title_fontsize=8)
    final = res["final_growth"]
    d.set_title(
        f"a joint as a stopband: the action grows {final['5.75 GHz facing boards, 610 MHz, 10 %']:.2f} times by 185 "
        "eV in "
        f"the 5.75 GHz guide at 10 %,\n{(final['5.75 GHz facing boards, 610 MHz, 5 %'] - 1) * 100:.0f} % at 5 %, and "
        f"{(final['1 GHz guide, 106 MHz secular, 10 %'] - 1) * 100:.0f} % in the 1 GHz guide: δ·\u2113_j ≲ 5 × 10⁻⁵ m",
        fontsize=9,
    )
    _suptitle(
        fig,
        "Resonances: a comb the beam can drive, the lines it radiates passing periodic electrodes, and the "
        "half-integer resonances its ramp crosses",
    )
    _save(fig, out, "resonances")


def figure_world(res: dict[str, Any], out: Path) -> None:
    """Figure 6: displacements by stray fields, argon's hold time, Sagnac and Aharonov-Bohm phases."""
    plt = _pyplot()
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.0), layout="constrained")
    guides = res["_guides"]
    a = axes[0, 0]
    field = np.logspace(-3, 3, 100)
    for (name, secular), colour in zip(guides.items(), (SERIES[0], SERIES[1]), strict=True):
        a.loglog(
            field,
            QM * field / (2 * math.pi * secular) ** 2 * 1e6,
            color=colour,
            label=f"{name} secular: {res['per_field_um'][name]:.3g} µm per V/m",
        )
    a.set_xlabel("stray electric field (V/m)")
    a.set_ylabel("displacement of a held electron (µm)")
    a.legend(loc="upper left", fontsize=8)
    a.set_title(
        "an electric field map: E moves a held electron by eE/(mω_sec²); the volts that null the\nexcess "
        "micromotion read E",
        fontsize=9,
    )
    b = axes[0, 1]
    tesla = np.logspace(-9, -3, 100)
    for (label, energy, name), colour, ls in (
        (("185 eV, 610 MHz", 185.0, "610 MHz"), SERIES[1], "-"),
        (("47 eV, 106 MHz", 47.0, "106 MHz"), SERIES[0], "--"),
        (("1 eV, 106 MHz", 1.0, "106 MHz"), SERIES[0], "-"),
    ):
        b.loglog(
            tesla * 1e6,
            QM * float(speed(energy)) * tesla / (2 * math.pi * guides[name]) ** 2 * 1e6,
            color=colour,
            ls=ls,
            label=label,
        )
    b.axvline(50, color=INK2, lw=0.9, ls=":")
    b.text(55, 2e-3, "the Earth's\n50 µT", fontsize=8, color=INK2)
    shifts = res["earth_shift_um"]
    b.set_xlabel("stray magnetic field across the path (µT)")
    b.set_ylabel("displacement of the beam (µm)")
    b.legend(loc="upper left", fontsize=8)
    b.set_title(
        "a magnetometer: B across the path moves the beam by evB/(mω_sec²): the Earth's field moves\n185 eV in the "
        "610 "
        f"MHz guide {shifts['185 eV, 610 MHz']:.1f} µm, 1 eV in the 106 MHz guide {shifts['1 eV, 106 MHz']:.0f} µm",
        fontsize=9,
    )
    c = axes[1, 0]
    energies, sigma, hold = res["_argon"]
    c.loglog(energies, hold * 1e3, color=SERIES[0], label="hold time in 10⁻⁵ mbar of argon")
    twin = c.twinx()
    twin.loglog(energies, sigma, color=SERIES[1], ls="--", lw=1.3, label="\u03c3 (right axis)")
    twin.set_ylabel("e-Ar total cross-section (m²)")
    twin.spines["right"].set_visible(True)
    argon = res["argon"]
    c.scatter([argon["longest_eV"]], [max(hold) * 1e3], color=INK, s=24, zorder=4)
    c.set_xlabel("electron energy (eV)")
    c.set_ylabel("hold time (ms)")
    handles = c.get_legend_handles_labels()[0] + twin.get_legend_handles_labels()[0]
    c.legend(handles=handles, loc="lower left", fontsize=8)
    c.set_title(
        f"hold-time spectroscopy: argon's Ramsauer-Townsend window keeps a {argon['longest_eV']:.2f} eV electron "
        f"{max(hold) * 1e3:.1f} ms\n(0.25 eV: {argon['hold_s']['0.25 eV'] * 1e3:.1f} ms) against "
        f"{argon['estimate_10eV_s'] * 1e6:.0f} µs at 10 eV; MERT with Golden and Bandel's A = −1.65 a₀",
        fontsize=9,
    )
    d = axes[1, 1]
    area = np.logspace(0, 4, 100)  # cm^2
    d.loglog(
        area,
        res["sagnac_rad_per_m2"] * area * 1e-4,
        color=SERIES[0],
        lw=2.2,
        label="Sagnac, the Earth's rotation at 40.7° N",
    )
    for flux_density, ls, label in ((1e-15, "--", "1 fT"), (1e-12, "-.", "1 pT"), (1e-9, ":", "1 nT")):
        d.loglog(
            area,
            pt.E * flux_density * area * 1e-4 / HBAR,
            color=SERIES[1],
            ls=ls,
            label=f"Aharonov-Bohm, {label} through the area",
        )
    for name, value in res["areas_cm2"].items():
        d.axvline(value, color=INK2, lw=0.8, ls=":")
        d.text(value * 1.05, 2e-4, name, rotation=90, fontsize=7.8, color=INK2)
    d.set_ylim(1e-4, 1e3)
    d.set_xlabel("area the electron's paths enclose (cm²)")
    d.set_ylabel("phase (rad)")
    d.legend(loc="upper left", fontsize=8)
    d.set_title(
        f"as a wave: {res['sagnac_circle_rad']:.3f} rad of Sagnac phase round the 932 mm circle "
        f"({res['areas_cm2']['932 mm circle']:.0f} cm²), {res['sagnac_trefoil_rad']:.3f} for the\ntrefoil; one "
        "radian of "
        f"Aharonov-Bohm phase per {res['fT_per_radian']:.1f} fT through the circle",
        fontsize=9,
    )
    _suptitle(
        fig,
        "What the beam can tell us: stray electric and magnetic fields, a gas and its cross-sections, the Earth's "
        "rotation and magnetic flux",
    )
    _save(fig, out, "world")


# ---- the run ---------------------------------------------------------------------------------------------------------
SECTIONS: dict[str, tuple[Callable[[bool], dict[str, Any]], Callable[[dict[str, Any], Path], None]]] = {
    "ring3d": (compute_ring3d, figure_ring3d),
    "sources": (compute_sources, figure_sources),
    "phase": (compute_phase, figure_phase),
    "speed": (compute_speed, figure_speed),
    "resonances": (compute_resonances, figure_resonances),
    "world": (compute_world, figure_world),
}
ORDER = tuple(SECTIONS)
KEEP = 64  # results.json keeps arrays and lists of at most this many numbers; the figures hold the rest


def slim(value: Any) -> Any:
    """The results without private entries (keys starting with '_') and large arrays, for results.json."""
    if isinstance(value, np.ndarray):
        return value if value.size <= KEEP else None
    if isinstance(value, dict):
        kept = {str(k): slim(v) for k, v in value.items() if not str(k).startswith("_")}
        return {k: v for k, v in kept.items() if v is not None}
    if isinstance(value, (list, tuple)):
        if len(value) > KEEP:
            return None
        return [slim(v) for v in value]
    return value


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument(
        "--quick", action="store_true", help="fewer electrons and coarser grids, into beam_lifecycle/quick/"
    )
    parser.add_argument("--only", default="", help="comma-separated sections, of: " + ", ".join(ORDER))
    parser.add_argument("--out", type=Path, default=None, help="where the figures, results.json and run.log go")
    args = parser.parse_args(argv)
    chosen = [name for name in args.only.split(",") if name] or list(ORDER)
    unknown = [name for name in chosen if name not in ORDER]
    if unknown:
        print(f"no section {unknown}; there are {list(ORDER)}", file=sys.stderr)
        return 2
    out: Path = args.out if args.out is not None else (OUT / "quick" if args.quick else OUT)
    out.mkdir(parents=True, exist_ok=True)
    LOG.clear()
    began = time.time()
    say(f"beam lifecycle study: {'quick' if args.quick else 'full'} run of {', '.join(chosen)}")
    results: dict[str, Any] = {}
    for name in chosen:
        start = time.time()
        results[name] = SECTIONS[name][0](args.quick)
        say(f"  [{name}: {time.time() - start:.1f} s]")
    say("\n== figures ==")
    for name in chosen:
        SECTIONS[name][1](results[name], out)
    say(f"\ndone in {time.time() - began:.0f} s")
    if chosen == list(ORDER):
        (out / "results.json").write_text(json.dumps(lab.plain(slim(results)), indent=1) + "\n")
        (out / "run.log").write_text("\n".join(LOG) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
