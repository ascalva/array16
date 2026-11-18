"""Radio-frequency confinement: the pseudo-potential of an oscillating field, and the quadrupole guide's numbers.

A charge in a field E(x) cos(wt) that varies little over its quiver motion is pushed toward where the field is
weak, as if by a potential q^2 E^2 / (4 m w^2). Along a line where the field vanishes, a quadrupole
E = g cos(wt) (-x, +y) holds it: its motion is Mathieu's equation with parameter q = 2 |q/m| g / w^2, stable
below 0.908, and for small q a slow oscillation at w_sec = q w / (2 sqrt 2).

Energies are in eV for a particle of one elementary charge; qm is its charge over mass (C/kg).
"""

from __future__ import annotations

import math

from .base import EPS0, QM_E, E


def _w(f: float) -> float:
    return 2 * math.pi * f


def ponderomotive_ev(e_peak: float, f: float, qm: float = QM_E) -> float:
    """Depth (eV) of the pseudo-potential where the field's peak amplitude is e_peak (V/m), at frequency f (Hz),
    against a place where it vanishes: |q/m| E^2 / (4 w^2)."""
    return abs(qm) * e_peak**2 / (4 * _w(f) ** 2)


def field_for_depth(depth_ev: float, f: float, qm: float = QM_E) -> float:
    """The peak field (V/m) whose pseudo-potential is depth_ev deep at frequency f. Inverse of ponderomotive_ev."""
    return math.sqrt(4 * _w(f) ** 2 * depth_ev / abs(qm))


def mathieu_q(gradient: float, f: float, qm: float = QM_E) -> float:
    """Mathieu's q for a quadrupole of field gradient `gradient` (V/m^2, peak) at frequency f."""
    return 2 * abs(qm) * gradient / _w(f) ** 2


def gradient_for_q(q: float, f: float, qm: float = QM_E) -> float:
    """The field gradient (V/m^2, peak) that gives Mathieu's q at frequency f. Inverse of mathieu_q."""
    return q * _w(f) ** 2 / (2 * abs(qm))


def secular_frequency(q: float, f: float) -> float:
    """The slow oscillation's frequency (Hz) in a quadrupole of Mathieu's q driven at f: q f / (2 sqrt 2), to
    first order in q (2 % high at q = 0.3, 10 % at q = 0.6)."""
    return q * f / (2 * math.sqrt(2))


def quadrupole_depth_ev(q: float, f: float, r0: float, qm: float = QM_E) -> float:
    """Depth (eV) of a quadrupole's pseudo-potential at distance r0 (m) from its axis: w_sec^2 r0^2 / (2 |q/m|),
    which is q g r0^2 / 8 for gradient g."""
    return _w(secular_frequency(q, f)) ** 2 * r0**2 / (2 * abs(qm))


def volts_for_depth(depth_ev: float, q: float) -> float:
    """The drive V0 (volts, peak; the quadrupole's potential is V0 (a^2 - b^2) / (2 r0^2)) that makes a well
    depth_ev deep at Mathieu's q: 8 depth / q, whatever the frequency and the size (RFG-1)."""
    return 8 * depth_ev / q


def frequency_times_distance(volts: float, q: float, qm: float = QM_E) -> float:
    """The product f r0 (Hz m) at which the drive V0 gives Mathieu's q: sqrt(2 |q/m| V0 / q) / (2 pi) (RFG-1)."""
    return math.sqrt(2 * abs(qm) * volts / q) / (2 * math.pi)


def bend_limit_ev(depth_ev: float, bend_radius: float, r0: float) -> float:
    """The kinetic energy (eV) at which a guide's bend pushes a particle half way to r0: the centrifugal force
    2 W / R against the guide's restoring force at r0 / 2, depth x R / (2 r0). An estimate of the most a guide
    holds round a bend of radius bend_radius; tracking has held about half of it."""
    return depth_ev * bend_radius / (2 * r0)


def space_charge_limit(q: float, f: float, qm: float = QM_E) -> float:
    """The largest density (1/m^3) of singly charged particles a quadrupole guide holds: where the beam's own
    field cancels the guide's restoring force, 2 eps0 w_sec^2 / (|q/m| e)."""
    return 2 * EPS0 * _w(secular_frequency(q, f)) ** 2 / (abs(qm) * E)
