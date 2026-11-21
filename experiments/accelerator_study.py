"""The first run, on the core: what radio-frequency fields take to hold a charge on a closed path and speed it up.

The owner's idea (2026-11-16): a field that acts as a small particle accelerator, its shape holding the particle's
path, with the particles' own fields refining it. This script puts numbers on the radio-frequency parts of it.

  scales     the pseudo-potential: the field a given depth takes, beside what the array radiates and what its
             elements make close by
  damping    calibration: a charged grain in air in an ideal quadrupole, against the damped pseudo-potential
  dust       a quadrupole guide on a trefoil knot at 60 Hz: charged glass grains in air, carried by a travelling wave
  electrons  the same shape at 1 GHz, for electrons
  ring       four gaps on a closed path, one chirp, fixed phase steps: acceleration along the path
  own        the particles' own fields: how many a guide holds, what a beam induces, when a line of grains buckles
  resonance  cyclotron resonance, which adds a static field (not in the default run)

    uv run python experiments/accelerator_study.py              # about 70 s: accelerator/
    uv run python experiments/accelerator_study.py resonance    # chosen sections only

Models: point charges and thin rods (line charges), no image charges, no gas collisions; one particle at a time
except where a section says otherwise. The magnetic guide of the first run is experiment 1, guide_study.py.
"""

from __future__ import annotations

import json
import math
import sys
import time
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from array16 import lab
from array16 import particles as pt
from array16.particles import Arr, BoolArr, ClosedCurve, GuideWatch, PathRecorder, QuadRods

OUT = Path(__file__).resolve().parent / "accelerator"
LOG: list[str] = []
G = 9.81

# The array's radiated power at full drive, as effective isotropic radiated power (REDESIGN.md, section 7): in its
# four-radio mode, and through the panel with the HackRF alone, the radio of the owner's decision of 2026-11-20
# (computed there for a HackRF One's output).
EIRP_DBM = {0.433e9: 15.9, 2.44e9: 8.1, 5.8e9: 6.2}
HACKRF_EIRP_DBM = {2.44e9: 1.5, 5.8e9: -17.8}

# A resonant gap's shunt resistance (ohm), two estimates: a metal cavity, and a resonator printed on the board with Q
# of 50 to 200 and sub-picofarad electrodes, R = Q / (w C).
SHUNT_OHM = {"metal cavity, 0.25-1 Mohm": (0.25e6, 1e6), "board resonator, 10-100 kohm": (10e3, 100e3)}

# A 10 um borosilicate sphere in air (the first run's grain, after Vinitsky et al. 2014).
GRAIN_MASS, GRAIN_QM, GRAIN_DRAG = 1.17e-12, 0.0176, 1450.0  # kg, C/kg, Stokes rate 1/s
DUST_F, DUST_R0, DUST_WIRE = 60.0, 6e-3, 0.75e-3  # Hz, m, m


def say(line: str = "") -> None:
    print(line, flush=True)
    LOG.append(line)


# ---- scales ----------------------------------------------------------------------------------------------------
def radiated_field(eirp_dbm: float, distance: float) -> float:
    """Peak field (V/m) at a distance from a source of the given EIRP, by the far-field formula sqrt(60 P) / r."""
    return math.sqrt(60 * 1e-3 * 10 ** (eirp_dbm / 10)) / distance


def scales() -> dict[str, Any]:
    say("\n== scales: the field a pseudo-potential of a given depth takes (electron) ==")
    rows = {}
    for f, eirp in EIRP_DBM.items():
        far = radiated_field(eirp, 0.3)
        row = {
            "eirp_dBm": eirp,
            "field_for_25meV_kVm": pt.field_for_depth(0.025, f) / 1e3,
            "field_for_1eV_kVm": pt.field_for_depth(1.0, f) / 1e3,
            "field_for_100eV_kVm": pt.field_for_depth(100.0, f) / 1e3,
            "radiated_at_0.3m_Vm": far,
            "depth_at_0.3m_eV": pt.ponderomotive_ev(far, f),
            "volts_across_1mm_for_1eV": pt.field_for_depth(1.0, f) * 1e-3,
            "distance_for_1eV_by_far_field_mm": radiated_field(eirp, 1.0) / pt.field_for_depth(1.0, f) * 1e3,
        }
        rows[f"{f / 1e9:g} GHz"] = row
        say(
            f"  {f / 1e9:5.3f} GHz: 25 meV takes {row['field_for_25meV_kVm']:5.1f} kV/m, 1 eV "
            f"{row['field_for_1eV_kVm']:5.1f}, "
            f"100 eV {row['field_for_100eV_kVm']:6.0f}. At {eirp:+.1f} dBm EIRP the field 0.3 m away is {far:.1f} V/m "
            f"({row['depth_at_0.3m_eV']:.1e} eV); 1 eV is {row['volts_across_1mm_for_1eV']:.0f} V across a millimetre"
        )
        if f in HACKRF_EIRP_DBM:
            alone = HACKRF_EIRP_DBM[f]
            row["hackrf_eirp_dBm"] = alone
            row["hackrf_radiated_at_0.3m_Vm"] = radiated_field(alone, 0.3)
            row["hackrf_depth_at_0.3m_eV"] = pt.ponderomotive_ev(row["hackrf_radiated_at_0.3m_Vm"], f)
            say(
                f"             through the panel with the HackRF alone, {alone:+.1f} dBm EIRP: "
                f"{row['hackrf_radiated_at_0.3m_Vm']:.2f} V/m 0.3 m away ({row['hackrf_depth_at_0.3m_eV']:.1e} eV)"
            )
    powers = (3.2e-3, 1.0)
    gap = {
        name: {f"{p * 1e3:g} mW": [pt.resonator_volts(p, r) for r in ohms] for p in powers}
        for name, ohms in SHUNT_OHM.items()
    }
    say("  a resonant gap makes sqrt(2 P R) volts (RFG-3) from its shunt resistance R, estimated:")
    for name, volts in gap.items():
        say(f"     {name}: " + "; ".join(f"{v[0]:.0f}-{v[1]:.0f} V from {k}" for k, v in volts.items()))
    path = {p: pt.KB * 293 / (p * 100 * 1e-19) for p in (1013.0, 1e-3, 1e-5, 1e-7)}  # pressure in mbar
    say(
        "  an electron's mean free path (cross-section 1e-19 m^2, an estimate): "
        + ", ".join(f"{v:.2g} m at {p:g} mbar" for p, v in path.items())
    )
    return {"bands": rows, "gap_volts": gap, "mean_free_path_m": {f"{p:g} mbar": v for p, v in path.items()}}


# ---- calibration: a damped grain in an ideal quadrupole -----------------------------------------------------------
class Tail:
    """Over the steps from `start` on: the mean of each particle's y, and the largest of its coordinates."""

    def __init__(self, start: int, count: int) -> None:
        self.start, self.total, self.seen = start, np.zeros(count), 0
        self.largest: Arr = np.zeros(count)

    def __call__(self, step: int, t: float, x: Arr, u: Arr, alive: BoolArr) -> None:
        if step >= self.start:
            self.total += x[:, 1]
            self.seen += 1
            self.largest = np.maximum(self.largest, np.abs(x).max(axis=1))

    @property
    def mean(self) -> Arr:
        return self.total / max(self.seen, 1)


def rod_gradient(volts: Arr | float, r0: float, wire: float) -> Any:
    """Field gradient on the axis of four thin rods with `volts` between neighbours (ROD-1)."""
    return 2 * volts / (math.log(r0 / wire) * r0**2)


def damping() -> dict[str, Any]:
    say("\n== damping: a 10 um glass grain in air, ideal quadrupole at 60 Hz, r0 6 mm ==")
    omega = 2 * math.pi * DUST_F
    volts = np.array([250.0, 500, 1000, 2000, 4000, 8000, 16000, 32000])
    grad = rod_gradient(volts, DUST_R0, DUST_WIRE)
    count, per, periods = len(volts), 200, 300

    def field(x: Arr, t: float) -> tuple[Arr, None]:
        return np.stack([-grad * x[:, 0], grad * x[:, 1], 0 * grad], 1) * math.cos(omega * t), None

    tail = Tail(per * (periods - 50), count)
    x0 = np.tile([2e-4, 2e-4, 0.0], (count, 1))
    pt.integrate(
        x0,
        np.zeros((count, 3)),
        np.full(count, GRAIN_QM),
        1 / (DUST_F * per),
        per * periods,
        field,
        [tail],
        drag=GRAIN_DRAG,
        acc=np.array([0.0, -G, 0.0]),
    )
    q = 2 * GRAIN_QM * grad / omega**2
    predicted = -G / ((GRAIN_QM * grad) ** 2 / (2 * (omega**2 + GRAIN_DRAG**2)))  # the damped pseudo-potential
    held = tail.largest < 0.8 * DUST_R0
    say("     volts   Mathieu q   held    sag tracked / predicted (mm)")
    for v, q_, h, s, p in zip(volts, q, held, tail.mean, predicted, strict=True):
        say(f"     {v:6.0f}   {q_:8.1f}   {'held' if h else 'LOST'}   {s * 1e3:8.3f} / {p * 1e3:8.3f}")
    at = int(np.argmin(np.abs(volts - 1000)))
    assert bool(held[at]), "the calibration case is lost"
    assert abs(tail.mean[at] / predicted[at] - 1) < 0.5, "the sag is not the damped pseudo-potential's"
    return {"volts": volts, "q": q, "held": held, "sag_mm": tail.mean * 1e3, "sag_predicted_mm": predicted * 1e3}


# ---- watching a transverse section --------------------------------------------------------------------------------
class Section:
    """One particle's offset on the rods' axes, every `every` steps, and every particle's mean height above the
    curve."""

    def __init__(self, rods: QuadRods, watch: GuideWatch, index: int, every: int) -> None:
        self.rods, self.watch, self.index, self.every = rods, watch, index, every
        self.offsets: list[tuple[float, float]] = []
        self.height: Arr = np.zeros(len(watch.station))
        self.seen = 0

    def __call__(self, step: int, t: float, x: Arr, u: Arr, alive: BoolArr) -> None:
        if step % self.every:
            return
        across, _ = self.rods.curve.offset(x, self.watch.station)
        self.height += across[:, 2]
        self.seen += 1
        if alive[self.index]:
            at = self.watch.station[self.index]
            frame = self.rods.frame
            self.offsets.append((float(across[self.index] @ frame.a[at]), float(across[self.index] @ frame.b[at])))


# ---- dust: a knotted quadrupole guide in air ----------------------------------------------------------------------
def dust() -> dict[str, Any]:
    say("\n== dust: quadrupole guide on a trefoil knot, 10 um glass grains in air ==")
    curve = ClosedCurve.torus_knot(0.060, 0.016, 601)
    rods = QuadRods.along(curve, DUST_R0)
    T, n, length = curve.tangents, len(curve), curve.length
    volts, speed = 1000.0, 0.05
    lam = math.pi * pt.EPS0 * volts / math.log(DUST_R0 / DUST_WIRE)
    omega = 2 * math.pi * DUST_F
    period = length / round(length / 0.03)  # the travelling wave's, along the guide
    slope = float(np.abs(T[:, 2]).max())
    wave = 1.6 * (GRAIN_DRAG * speed + G * slope) / GRAIN_QM  # V/m: 1.6 times what drag and the steepest climb take
    qm = GRAIN_QM * np.array([0.5, 1.0, 2.0])
    count, per = len(qm), 200
    dt, ds = 1 / (DUST_F * per), length / n
    k, w = 2 * math.pi / period, 2 * math.pi * speed / period
    watch = GuideWatch(curve, 0.8 * DUST_R0, count)

    def field(x: Arr, t: float) -> tuple[Arr, None]:
        at = watch.station
        along = at * ds + np.sum((x - curve.points[at]) * T[at], 1)
        e = rods.field_all(x, lam) * math.cos(omega * t)
        return e + T[at] * (wave * np.sin(k * along - w * t))[:, None], None

    section, path = Section(rods, watch, 1, 20), PathRecorder(1, 400)
    began = time.time()
    flight = pt.integrate(
        np.tile(curve.points[0], (count, 1)) + 3e-4,
        np.zeros((count, 3)),
        qm,
        dt,
        int(1.4 * length / speed / dt),
        field,
        [watch, section, path],
        drag=GRAIN_DRAG,
        acc=np.array([0.0, 0.0, -G]),
    )
    lap = [times[0] if times else math.nan for times in watch.times]
    say(
        f"  knot {length * 1e3:.0f} mm long, tightest bend {curve.bend_radius * 1e3:.0f} mm, steepest slope "
        f"{slope:.2f}; "
        f"rods +/-{volts / 2:.0f} V at {DUST_F:.0f} Hz"
    )
    say(
        f"  travelling wave {wave / 1e3:.1f} kV/m ({wave * period / (2 * math.pi):.0f} V on a "
        f"{period * 1e3:.0f} mm period) "
        f"moving at {speed * 100:.0f} cm/s"
    )
    for i in range(count):
        say(
            f"     charge x{qm[i] / GRAIN_QM:.1f}: {'held' if flight.alive[i] else 'LOST'}, "
            f"{watch.laps[i]:.2f} laps, first lap "
            f"in {lap[i]:.1f} s, reach {watch.reach[i] * 1e3:.2f} mm of {DUST_R0 * 1e3:.0f}, mean height "
            f"{section.height[i] / section.seen * 1e3:+.3f} mm"
        )
    say(f"  ({time.time() - began:.0f} s)")
    return {
        "length_mm": length * 1e3,
        "bend_mm": curve.bend_radius * 1e3,
        "volts": volts,
        "wave_kVm": wave / 1e3,
        "wave_volts": wave * period / (2 * math.pi),
        "held": flight.alive,
        "laps": watch.laps,
        "lap_s": lap,
        "reach_mm": watch.reach * 1e3,
        "mean_height_mm": section.height / section.seen * 1e3,
        "_curve": curve,
        "_path": path.path,
    }


# ---- electrons: the same shape at 1 GHz ----------------------------------------------------------------------------
def electrons() -> dict[str, Any]:
    say("\n== electrons: quadrupole guide on a trefoil knot at 1 GHz ==")
    f, r0, wire = 1e9, 0.5e-3, 25e-6
    curve = ClosedCurve.torus_knot(0.020, 0.006, 2001)
    rods = QuadRods.along(curve, r0)
    energy = np.array([0.5, 1, 2, 4, 8, 0.5, 1, 2, 4, 8, 16, 32.0])
    q = np.array([0.3] * 5 + [0.6] * 7)
    lam = np.array([rods.line_charge(pt.gradient_for_q(float(qi), f)) for qi in q])
    count = len(energy)
    v = np.array([pt.v_of(float(w)) for w in energy])
    gamma = 1 / np.sqrt(1 - (v / pt.C0) ** 2)
    omega, near = 2 * math.pi * f, np.arange(-60, 61, dtype=np.intp)
    watch = GuideWatch(curve, 0.8 * r0, count)

    def field(x: Arr, t: float) -> tuple[Arr, None]:
        return rods.field(x, watch.station, lam, near) * math.cos(omega * t), None

    shown = 8  # the 4 eV electron at q = 0.6
    section, path = Section(rods, watch, shown, 5), PathRecorder(shown, 100)
    dt = 1 / (f * 30)
    began = time.time()
    flight = pt.integrate(
        np.tile(curve.points[0], (count, 1)) + 3e-5 * rods.frame.a[0],
        (gamma * v)[:, None] * curve.tangents[0][None, :],
        np.full(count, pt.QM_E),
        dt,
        int(5 * curve.length / v.min() / dt),
        field,
        [watch, section, path],
    )
    say(
        f"  knot {curve.length * 1e3:.0f} mm long, tightest bend {curve.bend_radius * 1e3:.1f} mm, rods "
        f"{r0 * 1e3:.1f} mm from the axis"
    )
    cases = {}
    for qi in (0.3, 0.6):
        mine = np.flatnonzero(q == qi)
        depth = pt.quadrupole_depth_ev(qi, f, r0)
        limit = pt.bend_limit_ev(depth, curve.bend_radius, r0)
        volts = rods.volts(float(lam[mine[0]]), wire)
        say(
            f"   q {qi}: {volts:.0f} V between rods, secular {pt.secular_frequency(qi, f) / 1e6:.0f} MHz, "
            f"depth {depth:.2f} eV, "
            f"bend estimate of the most held {limit:.1f} eV"
        )
        for i in mine:
            say(
                f"      {energy[i]:4.1f} eV: {'held' if flight.alive[i] else 'LOST'} after {watch.laps[i]:5.1f} laps, "
                f"reach {watch.reach[i] * 1e6:5.0f} um"
            )
        cases[f"q {qi}"] = {
            "rod_volts": volts,
            "secular_MHz": pt.secular_frequency(qi, f) / 1e6,
            "depth_eV": depth,
            "bend_limit_eV": limit,
            "energy_eV": energy[mine],
            "held": flight.alive[mine],
            "laps": watch.laps[mine],
            "reach_um": watch.reach[mine] * 1e6,
        }
    say(f"  ({time.time() - began:.0f} s)")
    return {
        "length_mm": curve.length * 1e3,
        "bend_mm": curve.bend_radius * 1e3,
        "r0_mm": r0 * 1e3,
        "cases": cases,
        "_rods": rods,
        "_path": path.path,
        "_section": np.array(section.offsets),
        "_shown": f"{energy[shown]:g} eV at q = {q[shown]:g}",
    }


# ---- ring: four gaps, one chirp ------------------------------------------------------------------------------------
def ring() -> dict[str, Any]:
    say("\n== ring: four gaps on a closed path, one chirp, fixed phase steps ==")
    length = ClosedCurve.torus_knot(0.060, 0.028, 600).length
    start, gap_volts, phase, gaps, band = 50.0, 0.5, math.radians(30), 4, 54e6
    rest = pt.ME * pt.C0**2 / pt.E

    def speed(w: Any) -> Any:
        return pt.C0 * np.sqrt(1 - 1 / (1 + w / rest) ** 2)

    harmonic = round(60e6 / (float(speed(start)) / length))
    f0 = harmonic * float(speed(start)) / length
    # the synchronous particle gains gap_volts sin(phase) at each gap; the drive's phase is tabulated at its passages
    times, synchronous = [-(length / gaps) / float(speed(start)), 0.0], [start, start + gap_volts * math.sin(phase)]
    while harmonic * float(speed(synchronous[-1])) / length < f0 + band:
        times.append(times[-1] + (length / gaps) / float(speed(synchronous[-1])))
        synchronous.append(synchronous[-1] + gap_volts * math.sin(phase))
    t_sync, w_sync = np.array(times), np.array(synchronous)
    table = (np.arange(len(t_sync)) - 1) * harmonic / gaps
    rng = np.random.default_rng(1)
    count = 4000
    t = rng.uniform(-0.5, 0.5, count) / f0  # a continuous beam: every phase of the drive
    w = start + rng.uniform(-1, 1, count)
    sample = []
    for j in range(len(t_sync) - 1):
        seen = 2 * math.pi * (np.interp(t, t_sync, table) - j * harmonic / gaps) + phase
        w = np.maximum(w + gap_volts * np.sin(seen), 0.5)
        t = t + (length / gaps) / speed(w)
        if j % gaps == 0:
            sample.append(w[:40].copy())
    captured = np.abs(w - w_sync[-1]) < 0.05 * w_sync[-1]
    out = {
        "length_mm": length * 1e3,
        "harmonic": harmonic,
        "f_start_MHz": f0 / 1e6,
        "f_end_MHz": harmonic * float(speed(w_sync[-1])) / length / 1e6,
        "laps": (len(t_sync) - 2) / gaps,
        "time_us": float(t_sync[-1] * 1e6),
        "start_eV": start,
        "end_eV": float(w_sync[-1]),
        "captured": float(captured.mean()),
        "spread_eV": float(w[captured].std()),
        "gap_volts": gap_volts,
        "phase_step_deg": 360 * (harmonic / gaps % 1),
    }
    say(
        f"  a continuous {start:.0f} eV beam, {gaps} gaps of {gap_volts} V, harmonic {harmonic}, "
        f"{out['phase_step_deg']:.0f} deg "
        f"between gaps: chirp {out['f_start_MHz']:.1f} -> {out['f_end_MHz']:.1f} MHz in {out['time_us']:.0f} us "
        f"({out['laps']:.0f} laps)"
    )
    say(
        f"  reaches {out['end_eV']:.0f} eV; {out['captured'] * 100:.0f} % captured, {out['spread_eV']:.1f} eV "
        f"rms spread"
    )
    return {**out, "_sample": np.array(sample), "_synchronous": w_sync[1::gaps]}


# ---- the particles' own fields ---------------------------------------------------------------------------------
def own(ring_out: dict[str, Any]) -> dict[str, Any]:
    say("\n== own: the particles' own fields ==")
    f, beam = 1e9, 50e-6
    out: dict[str, Any] = {}
    for q in (0.3, 0.6):
        density = pt.space_charge_limit(q, f)
        per_metre = density * math.pi * beam**2
        current = per_metre * pt.E * pt.v_of(8.0)
        out[f"q {q}"] = {"density_m3": density, "per_metre": per_metre, "uA_at_8eV": current * 1e6}
        say(
            f"  the 1 GHz guide at q {q} holds up to {density:.1e} electrons/m^3: {per_metre:.1e} per metre in a beam "
            f"{beam * 1e6:.0f} um in radius, {current * 1e6:.1f} uA at 8 eV"
        )
    speed = pt.v_of(ring_out["end_eV"])
    length = ring_out["length_mm"] / 1e3
    signal = {}
    for electrons_ in (1e5, 1e6, 1e7):
        current = electrons_ * pt.E * speed / length
        dbm = 20 * math.log10(0.3 * 2 * current * 50 / math.sqrt(2) / 0.2236)
        signal[f"{electrons_:.0e}"] = {"uA": current * 1e6, "dBm": dbm}
        say(
            f"  {electrons_:.0e} electrons circulating at {ring_out['end_eV']:.0f} eV: {current * 1e6:.2f} uA, "
            f"about {dbm:.0f} dBm "
            "on a 50 ohm pickup coupling 30 % (an estimate)"
        )
    omega = 2 * math.pi * DUST_F
    grad = rod_gradient(1000.0, DUST_R0, DUST_WIRE)
    stiffness = GRAIN_MASS * (GRAIN_QM * grad) ** 2 / (2 * (omega**2 + GRAIN_DRAG**2))
    charge = GRAIN_QM * GRAIN_MASS
    zigzag = (4.207 * charge**2 / (4 * math.pi * pt.EPS0) / stiffness) ** (1 / 3)
    say(f"  grains closer than {zigzag * 1e3:.2f} mm along the 60 Hz guide buckle from a line into a zigzag")
    return {"space_charge": out, "pickup": signal, "zigzag_mm": zigzag * 1e3}


# ---- cyclotron resonance: needs a static field -------------------------------------------------------------------
def resonance_energy(e0: float, detune: float, f: float = 2.45e9, steps: int = 40000) -> tuple[Arr, Arr, float]:
    """Rotating-wave equations for an electron in a uniform field B_res (1 + detune) and a linearly polarised
    field e0 (V/m, peak) at the cyclotron frequency: times (s), energies (eV), and the time scale of the first peak."""
    w = 2 * math.pi * f
    kick = pt.E * e0 / (2 * pt.ME * pt.C0)
    scale = (12 * math.pi * pt.ME**2 * pt.C0**2 / (w * pt.E**2 * e0**2)) ** (1 / 3)
    dt, a = 6 * scale / steps, 0j
    times, energies = [], []

    def rate(z: complex) -> complex:
        return kick - 1j * (w * (1 + detune) / math.sqrt(1 + abs(z) ** 2) - w) * z

    for i in range(steps):
        k1 = rate(a)
        k2 = rate(a + dt / 2 * k1)
        k3 = rate(a + dt / 2 * k2)
        k4 = rate(a + dt * k3)
        a += dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if i % 40 == 0:
            times.append((i + 1) * dt)
            energies.append((math.sqrt(1 + abs(a) ** 2) - 1) * pt.ME * pt.C0**2 / pt.E)
    return np.array(times), np.array(energies), scale


def resonance() -> dict[str, Any]:
    say("\n== resonance: cyclotron resonance at 2.45 GHz (adds a static field) ==")
    f = 2.45e9
    w = 2 * math.pi * f
    b0 = w * pt.ME / pt.E
    rows = {}
    for e0 in (5.0, 250.0, 1e3, 1e4):
        peak = {f"{d:g}": float(resonance_energy(e0, d)[1].max()) for d in (-1e-3, -1e-4, 0.0, 1e-4, 1e-3, 1e-2)}
        times, energies, _ = resonance_energy(e0, 0.0)
        rows[f"{e0:g} V/m"] = {"peak_eV": peak, "time_us": float(times[int(np.argmax(energies))] * 1e6)}
        say(
            f"  {e0:7.0f} V/m: {peak['0']:7.1f} eV after {rows[f'{e0:g} V/m']['time_us']:.3f} us in a uniform "
            f"{b0 * 1e3:.1f} mT; "
            f"with the field off by 1e-4 {peak['0.0001']:.1f}, by 1e-3 {peak['0.001']:.1f}, by 1e-2 {peak['0.01']:.1f}"
        )
    e0 = 3e4  # full orbits against the rotating-wave equations
    _, energies, scale = resonance_energy(e0, 0.0)
    dt, top = 1 / (f * 100), 0.0
    x, u = np.zeros((1, 3)), np.zeros((1, 3))
    for it in range(int(6 * scale / dt)):
        e = np.array([[e0 * math.cos(w * (it + 0.5) * dt), 0, 0]])
        x, u = pt.push(x, u, np.array([pt.QM_E]), dt, e, np.array([[0, 0, b0]]))
        top = max(top, (math.sqrt(1 + float(np.sum(u * u)) / pt.C0**2) - 1) * pt.ME * pt.C0**2 / pt.E)
    assert abs(top / energies.max() - 1) < 0.05
    say(f"  check at 30 kV/m: full orbits {top:.0f} eV, rotating-wave {energies.max():.0f} eV")
    loop = pt.MU0 / (2 * 0.01) * 1e3  # mT per ampere-turn at the centre of a loop of 10 mm radius
    say(
        f"  a winding printed on a board gives {loop:.3f} mT per ampere-turn at the centre of a 10 mm loop, where the "
        f"resonance sits at {pt.E / pt.ME / (2 * math.pi) * loop * 1e-3 / 1e6:.2f} MHz per ampere-turn"
    )
    return {
        "b_res_mT": b0 * 1e3,
        "uniform": rows,
        "check": {"full_orbit_eV": top, "rotating_wave_eV": float(energies.max())},
    }


# ---- figures ----------------------------------------------------------------------------------------------------
def figures(results: dict[str, Any]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if "scales" in results:
        fig, ax = plt.subplots(figsize=(10, 5))
        f = np.logspace(math.log10(0.05e9), math.log10(10e9), 200)
        for depth, label in ((0.025, "25 meV (room temperature)"), (1.0, "1 eV"), (100.0, "100 eV"), (1e4, "10 keV")):
            ax.loglog(f / 1e9, [pt.field_for_depth(depth, float(fi)) for fi in f], label=label)
        for volts in (1, 10, 100, 1000):
            ax.axhline(volts / 1e-3, color="0.6", lw=0.6, ls=":")
            ax.text(0.052, volts / 1e-3 * 1.12, f"{volts} V across 1 mm", fontsize=7, color="0.35")
        for band, eirp in EIRP_DBM.items():
            ax.plot(band / 1e9, radiated_field(eirp, 0.3), "kv", ms=6)
        ax.plot([], [], "kv", label="the array's radiated field 0.3 m away, four radios")
        for band, eirp in HACKRF_EIRP_DBM.items():
            ax.plot(band / 1e9, radiated_field(eirp, 0.3), "v", color="k", mfc="none", ms=6)
        ax.plot([], [], "v", color="k", mfc="none", label="the same through the panel, the HackRF alone")
        drive = math.sqrt(2 * 50 * 1e-3 * 10**0.5)  # volts, peak: +5 dBm on 50 ohm
        ax.axhline(
            drive / 1e-3,
            color="k",
            lw=1,
            ls="--",
            label=f"{drive:.2f} V across 1 mm: a radio channel's +5 dBm, as it comes",
        )
        ax.axvspan(0.4, 6.0, color="C0", alpha=0.07)
        ax.set_xlabel("frequency (GHz); shaded, the array's band")
        ax.set_ylabel("peak field (V/m)")
        ax.set_title("The field a pseudo-potential of a given depth takes, for an electron", fontsize=10)
        ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
        ax.grid(alpha=0.3, which="both")
        fig.tight_layout()
        fig.savefig(OUT / "scales.png", dpi=130)
        plt.close(fig)

    if "dust" in results and "electrons" in results:
        dust_out, el = results["dust"], results["electrons"]
        fig = plt.figure(figsize=(15, 5))
        left: Any = fig.add_subplot(1, 3, 1, projection="3d")
        curve: ClosedCurve = dust_out["_curve"]
        c, path = curve.points * 1e3, dust_out["_path"] * 1e3
        left.plot(c[:, 0], c[:, 1], c[:, 2], color="0.75", lw=3)
        left.plot(path[:, 0], path[:, 1], path[:, 2], color="C3", lw=0.9)
        left.set_title(
            f"A 10 um grain in air, 60 Hz guide on a {dust_out['length_mm']:.0f} mm trefoil:\n"
            f"{dust_out['laps'][1]:.2f} laps, first lap in {dust_out['lap_s'][1]:.1f} s",
            fontsize=9,
        )
        middle: Any = fig.add_subplot(1, 3, 2, projection="3d")
        rods: QuadRods = el["_rods"]
        for j in range(4):
            rod = rods.mid[j] * 1e3
            middle.plot(rod[:, 0], rod[:, 1], rod[:, 2], color="C0" if rods.sign[j] > 0 else "C1", lw=0.5)
        path = el["_path"] * 1e3
        middle.plot(path[:, 0], path[:, 1], path[:, 2], color="C3", lw=0.6)
        middle.set_title(
            f"An electron of {el['_shown']}, 1 GHz guide on a {el['length_mm']:.0f} mm trefoil:\n"
            "rods of the two signs in blue and orange, the path in red",
            fontsize=9,
        )
        for axis in (left, middle):
            axis.set_xlabel("mm", fontsize=7)
            axis.set_ylabel("mm", fontsize=7)
            axis.set_zlabel("mm", fontsize=7)
            axis.tick_params(labelsize=6)
            axis.set_box_aspect((1, 1, 0.5))
        right = fig.add_subplot(1, 3, 3)
        offsets = el["_section"] * 1e6
        r0 = el["r0_mm"] * 1e3
        grid = np.linspace(-r0, r0, 121)
        gx, gy = np.meshgrid(grid, grid)
        depth = el["cases"]["q 0.6"]["depth_eV"]
        levels = [0.01, 0.05, 0.1, 0.25, 0.5, 1.0]
        contours = right.contour(gx, gy, depth * (gx**2 + gy**2) / r0**2, levels=levels, colors="0.6", linewidths=0.6)
        right.clabel(contours, fmt="%g eV", fontsize=6)
        right.plot(offsets[:, 0], offsets[:, 1], ".", color="C3", ms=1.5)
        for j in range(4):
            angle = j * math.pi / 2
            right.plot(r0 * math.cos(angle), r0 * math.sin(angle), "o", color="C0" if j % 2 == 0 else "C1", ms=8)
        right.set_aspect("equal")
        right.set_xlabel("offset along the first rod (um)")
        right.set_ylabel("offset along the second (um)")
        right.set_title(
            "The same electron across the guide: the pseudo-potential's\ncontours, the four rods, and where it went",
            fontsize=9,
        )
        fig.tight_layout()
        fig.savefig(OUT / "guides3d.png", dpi=130)
        plt.close(fig)

    if "ring" in results:
        ring_out = results["ring"]
        fig, ax = plt.subplots(figsize=(6.5, 4))
        sample = ring_out["_sample"]
        ax.plot(sample, color="C0", lw=0.4, alpha=0.6)
        ax.plot(ring_out["_synchronous"][: len(sample)], "k--", lw=1.2, label="the synchronous energy")
        ax.set_xlabel("lap")
        ax.set_ylabel("eV")
        ax.set_title(
            f"Four {ring_out['gap_volts']} V gaps, one chirp {ring_out['f_start_MHz']:.0f} to "
            f"{ring_out['f_end_MHz']:.0f} MHz: "
            f"40 of 4000 electrons, {ring_out['captured'] * 100:.0f} % captured",
            fontsize=9,
        )
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(OUT / "ring.png", dpi=130)
        plt.close(fig)


SECTIONS: dict[str, Callable[[], dict[str, Any]]] = {
    "scales": scales,
    "damping": damping,
    "dust": dust,
    "electrons": electrons,
    "ring": ring,
    "resonance": resonance,
}
DEFAULT = ("scales", "damping", "dust", "electrons", "ring", "own")


def main(argv: Sequence[str] | None = None) -> int:
    chosen = list(sys.argv[1:] if argv is None else argv) or list(DEFAULT)
    unknown = [name for name in chosen if name not in SECTIONS and name != "own"]
    if unknown:
        print(f"no section {unknown}; there are {[*SECTIONS, 'own']}", file=sys.stderr)
        return 2
    OUT.mkdir(exist_ok=True)
    began = time.time()
    results: dict[str, Any] = {}
    for name in chosen:
        if name == "own":
            results["own"] = own(results["ring"] if "ring" in results else ring())
        else:
            results[name] = SECTIONS[name]()
    say(f"\ndone in {time.time() - began:.0f} s")
    figures(results)
    if chosen == list(DEFAULT):
        kept = {name: {k: v for k, v in out.items() if not k.startswith("_")} for name, out in results.items()}
        (OUT / "results.json").write_text(json.dumps(lab.plain(kept), indent=1) + "\n")
        (OUT / "run.log").write_text("\n".join(LOG) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
