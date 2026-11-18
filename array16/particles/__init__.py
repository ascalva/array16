"""Charged particles in prescribed fields: the core the particle experiments are built on.

  base         array types, physical constants
  kinematics   v_of, rigidity, energy_of
  fields       Fields (E and B at positions and a time), and ways of combining them
  push         push: one relativistic Boris step
  curves       ClosedCurve: tangents, twist-free axes, holonomy, writhe
  rods         QuadRods: four thin rods along a curve
  rf           the pseudo-potential, and a quadrupole guide's Mathieu q, secular frequency, depth and limits
  coils        Coil: a filament winding round a curve, and its Biot-Savart field
  optics       GuideOptics: the first-order theory of a guide field; orbit_circle, window
  tracking     integrate: the general loop with observers; track_guide: electrons round a closed guide

Objects are immutable and functions pure, except the observers, which collect what they watch. Models: point
charges, thin wires, filament coils, no image charges, no collisions, no space charge. SI units throughout.
"""

from .base import C0, EPS0, KB, ME, MU0, QM_E, Z0, Arr, BoolArr, CArr, E, IntArr, Winding
from .coils import Coil
from .curves import ClosedCurve, Frame
from .fields import BField, Fields, magnetic, scaled, superpose, together, uniform
from .kinematics import energy_of, rigidity, v_of
from .optics import GuideOptics, OrbitFit, orbit_circle, window
from .push import cross, push
from .rf import (
    bend_limit_ev,
    field_for_depth,
    frequency_times_distance,
    gradient_for_q,
    mathieu_q,
    ponderomotive_ev,
    quadrupole_depth_ev,
    secular_frequency,
    space_charge_limit,
    volts_for_depth,
)
from .rods import QuadRods
from .tracking import Flight, GuideTrack, GuideWatch, Launch, Observer, PathRecorder, integrate, track_guide

__all__ = [
    "C0",
    "EPS0",
    "KB",
    "ME",
    "MU0",
    "QM_E",
    "Z0",
    "Arr",
    "BField",
    "BoolArr",
    "CArr",
    "ClosedCurve",
    "Coil",
    "E",
    "Fields",
    "Flight",
    "Frame",
    "GuideOptics",
    "GuideTrack",
    "GuideWatch",
    "IntArr",
    "Launch",
    "Observer",
    "OrbitFit",
    "PathRecorder",
    "QuadRods",
    "Winding",
    "bend_limit_ev",
    "cross",
    "energy_of",
    "field_for_depth",
    "frequency_times_distance",
    "gradient_for_q",
    "integrate",
    "magnetic",
    "mathieu_q",
    "orbit_circle",
    "ponderomotive_ev",
    "push",
    "quadrupole_depth_ev",
    "rigidity",
    "scaled",
    "secular_frequency",
    "space_charge_limit",
    "superpose",
    "together",
    "track_guide",
    "uniform",
    "v_of",
    "volts_for_depth",
    "window",
]
