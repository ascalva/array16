module

public import Array16.Confine.Defs
public import Mathlib.Analysis.Real.Pi.Bounds
public import Mathlib.Tactic.FieldSimp
public import Mathlib.Tactic.Linarith
public import Mathlib.Tactic.LinearCombination
public import Mathlib.Tactic.Ring

/-!
# Wells made by radio-frequency fields: RFG-2 to RFG-4

How deep a well a watt makes, by three routes that `array16/particles/rf.py` computes.

Carried in from far away and converging on a point as an electric dipole wave, a power `P` makes
a field whose square goes as `ω²`. The well a field makes goes as its square over `ω²`. So the
depth is `qm Z0 P / (3π c²)` at every frequency (RFG-2): 78.2 µeV per watt for an electron.

Fed through a resonator of shunt resistance `R` onto electrodes a distance `h` from the null, the
same power makes a well `(3π/2) (R/Z0) share / ((ω/c) h)²` times deeper (RFG-3): 1.14e6 times at
1 GHz, 0.5 mm and 10 kΩ on four rods. At volts `V` that well is `share V q / (8 γ)`, with `q` the
null's Mathieu `q` and `γ` its gradient over `V / h²`.

A guide filled to its space-charge limit has a plasma frequency `√2` times its secular frequency,
which is `qM ω / 2` (RFG-4).

Scope. These are identities between the formulas of `Array16/Confine/Defs.lean` over exact real
numbers, proved by algebra. What they rest on is stated elsewhere and not proved here. That the
electric dipole wave is the strongest field at a point per watt, over fields without sources near
the point that carry the power in from far away, is a theorem of the multipole expansion (Bassett
1986), cited by the ledger: RFG-2 here is the depth of that wave, and it is a ceiling of those
fields by that theorem. A structure that returns the wave, a cavity or a dense surround, stores
field and passes it. That a well's depth is `qm E² / (4ω²)` is the pseudo-potential, first order in
Mathieu's `q`, and that a well exists only below `q = 0.908` is Mathieu's stability (INT-3). That a
resonator holds `V² = 2PR` and that
electrodes of share `share` make that fraction of a four-rod well are the definitions of shunt
resistance and of share. That a uniform column cancels a harmonic well at `n = 2 ε₀ ω_sec² / (qm e)`
is Gauss's law for a cylinder. And the code computes in binary64, not in `ℝ`: the tests named
`test_rfg_N` attack the float implementation of each law.
-/

@[expose] public section

namespace Array16.Confine

open Real

/-! ### Freely radiated power: the same depth at every frequency -/

/-- **RFG-2** The well of an electric dipole wave carrying power `P` is `qm Z0 P / (3π c²)` deep,
whatever the frequency. -/
theorem ponderomotive_dipoleWaveFieldSq {ω c : ℝ} (hω : ω ≠ 0) (hc : c ≠ 0) (qm Z0 P : ℝ) :
    ponderomotive qm ω (dipoleWaveFieldSq Z0 c ω P) = ceiling qm Z0 c P := by
  unfold ponderomotive dipoleWaveFieldSq ceiling
  field_simp

/-- **RFG-2** So two frequencies give the same depth. -/
theorem ponderomotive_dipoleWaveFieldSq_eq {ω ω' c : ℝ} (hω : ω ≠ 0) (hω' : ω' ≠ 0) (hc : c ≠ 0)
    (qm Z0 P : ℝ) :
    ponderomotive qm ω (dipoleWaveFieldSq Z0 c ω P) =
      ponderomotive qm ω' (dipoleWaveFieldSq Z0 c ω' P) := by
  rw [ponderomotive_dipoleWaveFieldSq hω hc, ponderomotive_dipoleWaveFieldSq hω' hc]

/-- **RFG-2** The dipole wave's field squared is `2π/3` times that of the same power spread over a
spot half a wavelength on a side. -/
theorem dipoleWaveFieldSq_eq_spotFieldSq (Z0 c ω P : ℝ) :
    dipoleWaveFieldSq Z0 c ω P = 2 * π / 3 * spotFieldSq Z0 c ω P := by
  unfold dipoleWaveFieldSq spotFieldSq
  field_simp
  ring

/-- **RFG-2** control: the half-wave spot does not reach the ceiling. Its well is `3/(2π)` of it,
which is less than half. -/
theorem ponderomotive_spotFieldSq_ne {ω c : ℝ} (hω : ω ≠ 0) (hc : c ≠ 0) {qm Z0 P : ℝ}
    (h : ceiling qm Z0 c P ≠ 0) :
    ponderomotive qm ω (spotFieldSq Z0 c ω P) ≠ ceiling qm Z0 c P := by
  have e : ponderomotive qm ω (spotFieldSq Z0 c ω P) = 3 / (2 * π) * ceiling qm Z0 c P := by
    unfold ponderomotive spotFieldSq ceiling
    field_simp
    ring
  have h3 : 3 / (2 * π) < 1 := by
    rw [div_lt_one (by positivity)]
    linarith [pi_gt_three]
  intro heq
  rw [e] at heq
  have h0 : (3 / (2 * π) - 1) * ceiling qm Z0 c P = 0 := by linear_combination heq
  rcases mul_eq_zero.mp h0 with h1 | h1
  · linarith
  · exact h h1

/-! ### Through a resonator onto electrodes: the gain over the ceiling -/

/-- **RFG-3** Electrodes of share `share` at `h`, at the volts power `P` puts across a resonator
of shunt resistance `R`, make a well `P` times the depth per watt deep. -/
theorem share_mul_ponderomotive_resonatorVoltsSq (share qm R ω h P : ℝ) :
    share * ponderomotive qm ω (resonatorVoltsSq P R / h ^ 2) =
      P * electrodeDepthPerWatt share qm R ω h := by
  unfold ponderomotive resonatorVoltsSq electrodeDepthPerWatt
  by_cases hh : h = 0
  · simp [hh]
  by_cases hω : ω = 0
  · simp [hω]
  field_simp
  ring

/-- **RFG-3** That depth per watt is the gain `(3π/2) (R/Z0) share / ((ω/c) h)²` times the
free-space ceiling per watt, whatever the particle. -/
theorem electrodeDepthPerWatt_eq_gain_mul_ceiling {Z0 c ω h : ℝ} (hZ : Z0 ≠ 0) (hc : c ≠ 0)
    (hω : ω ≠ 0) (hh : h ≠ 0) (share qm R : ℝ) :
    electrodeDepthPerWatt share qm R ω h = electrodeGain R Z0 c ω h share * ceiling qm Z0 c 1 := by
  unfold electrodeDepthPerWatt electrodeGain ceiling
  field_simp

/-- **RFG-3** So at power `P` the electrodes' well is the gain times the well of RFG-2. -/
theorem share_mul_ponderomotive_eq_gain_mul {Z0 c ω h : ℝ} (hZ : Z0 ≠ 0) (hc : c ≠ 0)
    (hω : ω ≠ 0) (hh : h ≠ 0) (share qm R P : ℝ) :
    share * ponderomotive qm ω (resonatorVoltsSq P R / h ^ 2) =
      electrodeGain R Z0 c ω h share * ponderomotive qm ω (dipoleWaveFieldSq Z0 c ω P) := by
  rw [share_mul_ponderomotive_resonatorVoltsSq,
    electrodeDepthPerWatt_eq_gain_mul_ceiling hZ hc hω hh, ponderomotive_dipoleWaveFieldSq hω hc]
  unfold ceiling
  ring

/-- **RFG-3** At volts `V` the electrodes' well is `share V q / (8 γ)`, with `q` the null's
Mathieu `q` and `γ` its gradient over `V / h²`. The `q` grows with the volts, so the depth a watt
buys holds only up to the `q` at which Mathieu's equation stays stable (0.908, INT-3). -/
theorem share_mul_ponderomotive_eq_mathieuQ {γ ω h : ℝ} (hγ : γ ≠ 0) (hω : ω ≠ 0) (hh : h ≠ 0)
    (share qm V : ℝ) :
    share * ponderomotive qm ω (V ^ 2 / h ^ 2) = share * V * mathieuQ qm γ V ω h / (8 * γ) := by
  unfold ponderomotive mathieuQ
  field_simp
  ring

/-- **RFG-3** control: unlike the ceiling, the gain does depend on the frequency. At twice the
frequency it is a quarter. -/
theorem electrodeGain_two_mul_ne {R Z0 c ω h share : ℝ} (hg : electrodeGain R Z0 c ω h share ≠ 0) :
    electrodeGain R Z0 c (2 * ω) h share ≠ electrodeGain R Z0 c ω h share := by
  have e : electrodeGain R Z0 c (2 * ω) h share = electrodeGain R Z0 c ω h share / 4 := by
    unfold electrodeGain
    rw [show 2 * ω / c * h = 2 * (ω / c * h) by ring, mul_pow, ← div_div]
    norm_num
    ring
  intro heq
  rw [e] at heq
  exact hg (by linarith)

/-! ### The space-charge limit: plasma frequency and secular frequency -/

/-- **RFG-4** At the space-charge limit the plasma frequency squared is twice the secular frequency
squared. -/
theorem plasmaFrequencySq_spaceChargeLimit {eps0 qm e : ℝ} (heps : eps0 ≠ 0) (hqm : qm ≠ 0)
    (he : e ≠ 0) (ωs : ℝ) :
    plasmaFrequencySq eps0 qm e (spaceChargeLimit eps0 qm e ωs) = 2 * ωs ^ 2 := by
  unfold plasmaFrequencySq spaceChargeLimit
  field_simp

/-- The secular frequency squared is `(qM ω)² / 8`. -/
theorem secular_sq (qM ω : ℝ) : secular qM ω ^ 2 = (qM * ω) ^ 2 / 8 := by
  unfold secular
  rw [div_pow, mul_pow, mul_pow, sq_sqrt (by norm_num : (0 : ℝ) ≤ 2)]
  norm_num

/-- **RFG-4** So in a guide at Mathieu's `qM` driven at `ω` the plasma frequency at the limit is
`qM ω / 2`, whatever the particle. -/
theorem plasmaFrequencySq_spaceChargeLimit_secular {eps0 qm e : ℝ} (heps : eps0 ≠ 0)
    (hqm : qm ≠ 0) (he : e ≠ 0) (qM ω : ℝ) :
    plasmaFrequencySq eps0 qm e (spaceChargeLimit eps0 qm e (secular qM ω)) = (qM * ω / 2) ^ 2 := by
  rw [plasmaFrequencySq_spaceChargeLimit heps hqm he, secular_sq]
  ring

/-- **RFG-4** control: at half the limiting density the plasma frequency is the secular frequency,
not `√2` times it. -/
theorem plasmaFrequencySq_half_spaceChargeLimit_ne {eps0 qm e ωs : ℝ} (heps : eps0 ≠ 0)
    (hqm : qm ≠ 0) (he : e ≠ 0) (hωs : ωs ≠ 0) :
    plasmaFrequencySq eps0 qm e (spaceChargeLimit eps0 qm e ωs / 2) ≠ 2 * ωs ^ 2 := by
  have e1 : plasmaFrequencySq eps0 qm e (spaceChargeLimit eps0 qm e ωs / 2) = ωs ^ 2 := by
    unfold plasmaFrequencySq spaceChargeLimit
    field_simp
  rw [e1]
  have : 0 < ωs ^ 2 := by positivity
  linarith

end Array16.Confine
