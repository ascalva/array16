"""Radio-frequency confinement: the pseudo-potential of an oscillating field, the quadrupole guide's numbers, and
the depth a watt buys.

A charge in a field E(x) cos(wt) that varies little over its quiver motion is pushed toward where the field is
weak, as if by a potential q^2 E^2 / (4 m w^2). Along a line where the field vanishes, a quadrupole
E = g cos(wt) (-x, +y) holds it: its motion is Mathieu's equation with parameter q = 2 |q/m| g / w^2, stable
below 0.908, and for small q a slow oscillation at w_sec = q w / (2 sqrt 2). The pseudo-potential is first order
in q: its depth holds for q below about 0.5, and past 0.908 there is no well.

Power sets the depth. A field that carries P in from far away, with no sources near a point, makes its strongest
field there when it converges as an electric dipole wave, sqrt(4 Z0 k^2 P / (3 pi)), so a well at most
|q/m| Z0 P / (3 pi c^2) deep at any frequency: 78.2 micro-eV per watt for an electron, so 1 eV takes 12.8 kW
(dipole_wave_peak_field, free_space_depth_ceiling_ev; RFG-2). A structure that returns the wave passes that ceiling
by its stored energy: a cavity, a surround dense enough to close, an element whose own stored field reaches the
point. Fed through a resonator of shunt resistance R, the same power puts sqrt(2 P R) volts on electrodes a
distance h from the null and stores the field there instead of radiating it: the well is (3 pi / 2)(R / Z0) /
(k h)^2 times deeper, times the electrodes' share of a four-rod well (resonator_volts, electrode_depth_per_watt_ev,
electrode_gain; RFG-3). Its depth is share V q / (8 gamma), gamma being the null's gradient over V / h^2, so a
watt buys that depth while the null's q stays under 0.908 (electrode_drive, electrode_drive_at_q). At 1 GHz,
h = 0.5 mm and R = 10 kOhm four rods make 89.1 eV per watt up to 2.9 eV, and 1 eV takes 15 V and 11 mW at
q = 0.53; electrodes.py gives the share and the gradient of patterns printed on boards. A guide holds particles up
to the density at which their plasma frequency reaches sqrt(2) times the secular frequency, q f / 2
(plasma_frequency; RFG-4).

Energies are in eV for a particle of one elementary charge; qm is its charge over mass (C/kg).
"""

from __future__ import annotations

import math
from typing import NamedTuple

from .base import C0, EPS0, QM_E, Z0, E


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


def plasma_frequency(density: float, qm: float = QM_E) -> float:
    """The plasma frequency (Hz) of singly charged particles at a density (1/m^3): sqrt(n e |q/m| / eps0) / (2 pi),
    8.98 sqrt(n) for electrons. At a guide's space-charge limit it is sqrt(2) times the secular frequency, q f / 2,
    whatever the particle (RFG-4): 300 MHz in a 1 GHz guide at q = 0.6."""
    return math.sqrt(density * E * abs(qm) / EPS0) / (2 * math.pi)


def dipole_wave_peak_field(power: float, f: float) -> float:
    """The peak field (V/m) at the centre of a converging electric dipole wave that carries power (W) at
    frequency f (Hz): sqrt(4 Z0 k^2 P / (3 pi)), k = 2 pi f / c; 265 V/m for a watt at 1 GHz.

    At a point of a region without sources only the three electric dipole modes of a field are nonzero, so a wave
    that is all electric dipole puts all its power into the field there. Over fields without sources near the point
    that carry P in from far away, with nothing returning the wave, it is the strongest field P makes at a point
    (Bassett 1986, cited; RFG-2). Spread over a spot (lambda / 2)^2 wide, the same power gives a field squared
    2 pi / 3 times weaker. A structure that returns the wave to the point passes it by its stored energy, by the
    order of its Q in field squared and so in depth: a cavity, a surround of elements dense enough to close, an
    element whose own stored field reaches the point; electrodes fed through a resonator do it by electrode_gain."""
    k = _w(f) / C0
    return math.sqrt(4 * Z0 * k**2 * power / (3 * math.pi))


def free_space_depth_ceiling_ev(power: float, qm: float = QM_E) -> float:
    """The deepest pseudo-potential well (eV) that power (W) carried in from far away makes at a point without
    sources near it: |q/m| Z0 P / (3 pi c^2), the depth at the centre of the dipole wave, the same at every
    frequency (RFG-2). For an electron 78.2 micro-eV per watt: 1 eV takes 12.8 kW in free space. Stored field passes
    it, as dipole_wave_peak_field says."""
    return abs(qm) * Z0 * power / (3 * math.pi * C0**2)


def resonator_volts(power: float, shunt_resistance: float) -> float:
    """The peak volts across a resonator of shunt resistance R (ohm) that takes power P (W): sqrt(2 P R). 3.2 mW
    into 0.25 MOhm is 40 V; 1 W into 1 MOhm is 1.4 kV."""
    return math.sqrt(2 * power * shunt_resistance)


def electrode_depth_per_watt_ev(
    shunt_resistance: float, distance: float, f: float, share: float = 1.0, qm: float = QM_E
) -> float:
    """The depth (eV) per watt of a well made by electrodes a distance h (m) from the null, fed at frequency f (Hz)
    through a resonator of shunt resistance R (ohm): share |q/m| R / (2 w^2 h^2).

    The resonator puts V = sqrt(2 P R) volts, peak, across itself, and so between the two electrode sets it
    drives: rod pair to rod pair, half to half on facing boards, the rails to the ground plane on one board. Four
    rods at r0 = h driven at RFG-1's V0 = V make a field V / h at the rods and a well |q/m| V^2 / (4 w^2 h^2) deep;
    share is the geometry's fraction of that four-rod well: 1 for four rods, 1/4 for two facing boards, at most
    0.0091 for two rails on one board (the depth_share of electrodes.FacingBoards and electrodes.RailPair).

    The depth is the pseudo-potential's, first order in the null's Mathieu q, which grows with the volts
    (electrode_drive): the slope holds up to q = 0.908, past which there is no well. At 1 GHz, h = 0.5 mm and
    R = 10 kOhm: 89.1 eV per watt on four rods up to 2.9 eV at 32 mW, 22.3 on two facing boards up to 1.2 eV at
    53 mW, 0.81 on one board up to 0.28 eV at 0.34 W (RFG-3)."""
    return share * abs(qm) * shunt_resistance / (2 * _w(f) ** 2 * distance**2)


def electrode_gain(shunt_resistance: float, distance: float, f: float, share: float = 1.0) -> float:
    """How many times deeper, per watt, electrodes a distance h (m) from the null, fed through a resonator of shunt
    resistance R (ohm) at frequency f (Hz), make a well than the free-space ceiling: (3 pi / 2)(R / Z0) share /
    (k h)^2, k = 2 pi f / c, whatever the particle. The resonator stores the power in its near field instead of
    radiating it, which gains R / Z0, and the electrodes gather the field onto a distance h instead of a
    wavelength, which gains 1 / (k h)^2. 1.14e6 on four rods at 1 GHz, h = 0.5 mm and R = 10 kOhm; 2.85e5 on two
    facing boards (RFG-3). Like the depth per watt, it holds while the null's Mathieu q stays under 0.908."""
    k = _w(f) / C0
    return 1.5 * math.pi * (shunt_resistance / Z0) * share / (k * distance) ** 2


class ElectrodeDrive(NamedTuple):
    """What a well takes from electrodes fed through a resonator: the power (W) into the resonator, the peak volts
    across it, Mathieu's q at the null, and the drive frequency (Hz). The well is there only for q under 0.908."""

    power: float
    volts: float
    q: float
    f: float


def electrode_drive(
    depth_ev: float,
    shunt_resistance: float,
    distance: float,
    f: float,
    share: float = 1.0,
    gradient_share: float = 1.0,
    qm: float = QM_E,
) -> ElectrodeDrive:
    """What a well depth_ev deep takes at frequency f (Hz) from electrodes a distance h (m) from the null, fed
    through a resonator of shunt resistance R (ohm): the volts V = 2 w h sqrt(depth / (share |q/m|)), the power
    V^2 / (2 R), and the null's Mathieu q = 2 |q/m| gamma V / (w h)^2, gamma being gradient_share, the null's
    gradient over V / h^2 (1 for four rods; the gradient_share of electrodes.FacingBoards and
    electrodes.RailPair). The depth is share V q / (8 gamma) (RFG-3).

    q comes out as the volts make it: the depth, first order in q, holds for q below about 0.5, and past 0.908
    there is no well. At 1 GHz, h = 0.5 mm and 10 kOhm, 1 eV takes 15 V and 11 mW at q = 0.53 on four rods, 30 V at
    q = 0.84 on two facing boards, and 157 V at q = 1.73 on one board at its best, where it takes a higher frequency
    (electrode_drive_at_q)."""
    w = _w(f)
    volts = 2 * w * distance * math.sqrt(depth_ev / (share * abs(qm)))
    q = 2 * abs(qm) * gradient_share * volts / (w * distance) ** 2
    return ElectrodeDrive(volts**2 / (2 * shunt_resistance), volts, q, f)


def electrode_drive_at_q(
    depth_ev: float,
    q: float,
    shunt_resistance: float,
    distance: float,
    share: float = 1.0,
    gradient_share: float = 1.0,
    qm: float = QM_E,
) -> ElectrodeDrive:
    """What a well depth_ev deep takes at Mathieu's q from electrodes a distance h (m) from the null, fed through a
    resonator of shunt resistance R (ohm): the volts V = 8 gamma depth / (share q), the power V^2 / (2 R), and the
    frequency at which those volts give that q, f h = sqrt(2 |q/m| gamma V / q) / (2 pi), gamma being
    gradient_share as in electrode_drive. For four rods (share and gamma 1) these are RFG-1's drive and frequency
    times distance (RFG-3).

    At q = 0.5 and 10 kOhm, 1 eV takes 16 V and 13 mW on four rods, at 1.07 GHz for h = 0.5 mm; 50 V and 126 mW on
    two facing boards, at 1.68 GHz; 542 V and 14.7 W on one board at its best, at 3.46 GHz."""
    volts = 8 * gradient_share * depth_ev / (share * q)
    f = math.sqrt(2 * abs(qm) * gradient_share * volts / q) / (2 * math.pi * distance)
    return ElectrodeDrive(volts**2 / (2 * shunt_resistance), volts, q, f)
