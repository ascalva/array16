"""Charged particles in prescribed fields: the core the particle experiments are built on.

  base         array types, physical constants
  kinematics   v_of, rigidity, energy_of
  fields       Fields (E and B at positions and a time), and ways of combining them
  push         push: one relativistic Boris step
  curves       ClosedCurve: tangents, twist-free axes, holonomy, writhe
  rods         QuadRods: four thin rods along a curve
  sections     SectionGrid, SectionField: a solved cross-section swept along a path with a frame, driven as a
               TravellingWave, a StandingWave or InPhase; BendShift, the bend's first-order shift of the facing
               boards' null toward each bend's centre (BEND_IN_PLANE, BEND_TOWARD), off unless asked for
  rf           the pseudo-potential; a quadrupole guide's Mathieu q, secular frequency, depth and limits; the depth
               a watt buys, radiated freely or stored in a resonator on electrodes, and the volts and q it takes
  electrodes   RailPair, FacingBoards: electrodes printed on boards, their nulls, gradients and share of a
               four-rod well
  coils        Coil: a filament winding round a curve, and its Biot-Savart field
  optics       GuideOptics: the first-order theory of a guide field; orbit_circle, window
  tracking     integrate: the general loop with observers; track_guide: electrons round a closed guide

Objects are immutable and functions pure, except the observers, which collect what they watch. Models: point
charges, thin wires, filament coils, board electrodes without gaps; no image charges, no collisions, no space
charge. SI units throughout.
"""

from .base import C0, EPS0, KB, ME, MU0, QM_E, Z0, Arr, BoolArr, CArr, E, IntArr, Winding
from .coils import Coil
from .curves import ClosedCurve, Frame
from .electrodes import FacingBoards, RailPair, Saddle
from .fields import BField, Fields, magnetic, scaled, superpose, together, uniform
from .kinematics import energy_of, rigidity, v_of
from .optics import GuideOptics, OrbitFit, orbit_circle, window
from .push import cross, push
from .rf import (
    ElectrodeDrive,
    bend_limit_ev,
    dipole_wave_peak_field,
    electrode_depth_per_watt_ev,
    electrode_drive,
    electrode_drive_at_q,
    electrode_gain,
    field_for_depth,
    free_space_depth_ceiling_ev,
    frequency_times_distance,
    gradient_for_q,
    mathieu_q,
    plasma_frequency,
    ponderomotive_ev,
    quadrupole_depth_ev,
    resonator_volts,
    secular_frequency,
    space_charge_limit,
    volts_for_depth,
)
from .rods import QuadRods
from .sections import (
    BEND_IN_PLANE,
    BEND_TOWARD,
    BendShift,
    InPhase,
    SectionField,
    SectionGrid,
    StandingWave,
    TravellingWave,
)
from .tracking import Flight, GuideTrack, GuideWatch, Launch, Observer, PathRecorder, integrate, track_guide

__all__ = [
    "BEND_IN_PLANE",
    "BEND_TOWARD",
    "C0",
    "EPS0",
    "KB",
    "ME",
    "MU0",
    "QM_E",
    "Z0",
    "Arr",
    "BField",
    "BendShift",
    "BoolArr",
    "CArr",
    "ClosedCurve",
    "Coil",
    "E",
    "ElectrodeDrive",
    "FacingBoards",
    "Fields",
    "Flight",
    "Frame",
    "GuideOptics",
    "GuideTrack",
    "GuideWatch",
    "InPhase",
    "IntArr",
    "Launch",
    "Observer",
    "OrbitFit",
    "PathRecorder",
    "QuadRods",
    "RailPair",
    "Saddle",
    "SectionField",
    "SectionGrid",
    "StandingWave",
    "TravellingWave",
    "Winding",
    "bend_limit_ev",
    "cross",
    "dipole_wave_peak_field",
    "electrode_depth_per_watt_ev",
    "electrode_drive",
    "electrode_drive_at_q",
    "electrode_gain",
    "energy_of",
    "field_for_depth",
    "free_space_depth_ceiling_ev",
    "frequency_times_distance",
    "gradient_for_q",
    "integrate",
    "magnetic",
    "mathieu_q",
    "orbit_circle",
    "plasma_frequency",
    "ponderomotive_ev",
    "push",
    "quadrupole_depth_ev",
    "resonator_volts",
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
