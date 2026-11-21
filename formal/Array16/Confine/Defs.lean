module

public import Mathlib.Analysis.SpecialFunctions.Sqrt
public import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic

/-!
# Wells made by radio-frequency fields: RFG-2 to RFG-4, the definitions

The definitions the laws of `Array16/Confine.lean` are stated over. They mirror the functions of
`array16/particles/rf.py`, and a reader checks them against that module and against `docs/LAWS.md`.

Depths are in electronvolts, as in the code: a depth in joules divided by the charge. So a particle
enters through `qm`, the size of its charge over its mass. Frequencies here are angular, `ω = 2πf`,
where the code takes `f` and forms `ω` itself. The code returns a field and a voltage as square
roots; the laws use only their squares, and the definitions here are the squares.
-/

@[expose] public section

namespace Array16.Confine

open Real

/-- The depth in eV of the well a field of peak size squared `Esq`, oscillating at `ω`, makes for a
particle of charge over mass `qm` (`rf.ponderomotive_ev`). -/
noncomputable def ponderomotive (qm ω Esq : ℝ) : ℝ :=
  qm * Esq / (4 * ω ^ 2)

/-- The square of the peak field at the centre of an electric dipole wave carrying power `P` at
`ω`, in a medium of impedance `Z0` and wave speed `c` (the square of
`rf.dipole_wave_peak_field`). -/
noncomputable def dipoleWaveFieldSq (Z0 c ω P : ℝ) : ℝ :=
  4 * Z0 * (ω / c) ^ 2 * P / (3 * π)

/-- The square of the peak field of power `P` spread evenly over a spot half a wavelength on a
side, `2 Z0 P / (λ/2)²` with `λ = 2π c / ω`. -/
noncomputable def spotFieldSq (Z0 c ω P : ℝ) : ℝ :=
  2 * Z0 * (ω / c) ^ 2 * P / π ^ 2

/-- The deepest well in eV that freely radiated power `P` makes
(`rf.free_space_depth_ceiling_ev`). -/
noncomputable def ceiling (qm Z0 c P : ℝ) : ℝ :=
  qm * Z0 * P / (3 * π * c ^ 2)

/-- The square of the peak volts that power `P` puts across a resonator of shunt resistance `R`
(the square of `rf.resonator_volts`). -/
def resonatorVoltsSq (P R : ℝ) : ℝ :=
  2 * P * R

/-- The depth in eV per watt of the well of electrodes of share `share` a distance `h` from the
null, fed through a resonator of shunt resistance `R` (`rf.electrode_depth_per_watt_ev`). -/
noncomputable def electrodeDepthPerWatt (share qm R ω h : ℝ) : ℝ :=
  share * qm * R / (2 * ω ^ 2 * h ^ 2)

/-- Mathieu's `q` at a null whose field gradient is `γ V / h²`, driven at `ω`, for a particle of
charge over mass `qm` (`rf.mathieu_q` of that gradient, as `rf.electrode_drive` forms it; `γ` is
the geometry's `gradient_share`, 1 for four rods). -/
noncomputable def mathieuQ (qm γ V ω h : ℝ) : ℝ :=
  2 * qm * (γ * V / h ^ 2) / ω ^ 2

/-- The electrodes' depth per watt over the free-space ceiling per watt (`rf.electrode_gain`). -/
noncomputable def electrodeGain (R Z0 c ω h share : ℝ) : ℝ :=
  3 * π / 2 * (R / Z0) * share / (ω / c * h) ^ 2

/-- The angular secular frequency of a quadrupole guide at Mathieu's `q`, here `qM`, driven at `ω`
(`rf.secular_frequency`, times `2π`). -/
noncomputable def secular (qM ω : ℝ) : ℝ :=
  qM * ω / (2 * √2)

/-- The density at which a uniform column's own field cancels a guide's well of angular secular
frequency `ωs`, for particles of charge `e` and charge over mass `qm` (`rf.space_charge_limit`). -/
noncomputable def spaceChargeLimit (eps0 qm e ωs : ℝ) : ℝ :=
  2 * eps0 * ωs ^ 2 / (qm * e)

/-- The square of the angular plasma frequency at density `n` (the square of
`rf.plasma_frequency`, times `(2π)²`). -/
noncomputable def plasmaFrequencySq (eps0 qm e n : ℝ) : ℝ :=
  n * e * qm / eps0

end Array16.Confine
