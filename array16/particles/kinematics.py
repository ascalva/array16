"""Energy, speed and rigidity of a particle."""

from __future__ import annotations

import math

from .base import C0, ME, E


def v_of(w_ev: float, m: float = ME) -> float:
    """Speed (m/s) of a particle of mass m (kg) with kinetic energy w_ev (eV)."""
    g = 1 + w_ev * E / (m * C0**2)
    return C0 * math.sqrt(1 - 1 / g**2)


def rigidity(w_ev: float, m: float = ME) -> float:
    """Momentum per elementary charge (T m): the field-times-radius that bends the particle."""
    v = v_of(w_ev, m)
    return m * v / math.sqrt(1 - (v / C0) ** 2) / E


def energy_of(b_rho: float, m: float = ME) -> float:
    """Kinetic energy (eV) of a singly charged particle of mass m with rigidity b_rho (T m). Inverse of rigidity."""
    rest = m * C0**2
    return (math.hypot(E * b_rho * C0, rest) - rest) / E
