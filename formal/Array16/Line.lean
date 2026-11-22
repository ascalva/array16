module

public import Array16.Line.Defs
public import Mathlib.Analysis.Calculus.Deriv.Comp
public import Mathlib.Analysis.Calculus.Deriv.Inv
public import Mathlib.Analysis.Calculus.Deriv.Mul
public import Mathlib.Analysis.Calculus.Deriv.Pow
public import Mathlib.Analysis.SpecialFunctions.ExpDeriv
public import Mathlib.Analysis.SpecialFunctions.Log.Deriv
public import Mathlib.Tactic.FieldSimp
public import Mathlib.Tactic.Linarith
public import Mathlib.Tactic.Ring

/-!
# The line of a cross-section: LIN-1

`array16/guide/sections.py` solves a cross-section's potential and reports the laminate's share of
its field energy, from which the dielectric's `Q_d = 1 / (share tan δ)` and the index's response to
the laminate's permittivity follow. DR-0 holds the solve to a rule: the share is
`d ln C' / d ln ε_r`. This is why the rule is exact.

The energy of a potential is `A u + ε B u`, and the solved potential `u ε` makes it stationary. As
`ε` changes the potential moves, and that movement changes the energy by nothing to first order,
because the energy is stationary in the potential. What is left is the explicit `ε`: the
capacitance changes at the rate `B (u ε)`, the laminate's energy per unit permittivity. In
logarithms that is `ε B / C`, the share.

The smallest solve, two capacitors in series with one free node, is worked through: its divider is
stationary, its energy there is the series capacitance, and its share is `ca / (ca + ε cl)`, the
logarithmic slope of that capacitance.

Scope. These are laws of an energy `A u + ε B u` over exact real numbers, proved by the rules of
differentiation. What they rest on is stated elsewhere and not proved here. That a section's
capacitance is its least field energy and that the laminate enters linearly in `ε` is
electrostatics (Thomson's principle). That the finite-volume solve reaches the stationary point of
its own grid's energy is the solver's residual, which its tests hold. And the code computes in
binary64 and takes the slope by central differences at `±1 %`: the test named `test_lin_1` attacks
the float implementation.
-/

@[expose] public section

namespace Array16.Line

/-- **LIN-1** The capacitance of a solved section changes with the laminate's permittivity at the
rate of the laminate's energy per unit permittivity. `u` is the solved potential as `ε` varies,
and `hstat` says the energy is stationary in the potential at the solve. -/
theorem hasDerivAt_energy_of_stationary {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    {A B : E → ℝ} {A' B' : E →L[ℝ] ℝ} {u : ℝ → E} {u' : E} {ε : ℝ}
    (hA : HasFDerivAt A A' (u ε)) (hB : HasFDerivAt B B' (u ε)) (hu : HasDerivAt u u' ε)
    (hstat : A' + ε • B' = 0) :
    HasDerivAt (fun e ↦ energy A B e (u e)) (B (u ε)) ε := by
  have h1 : HasDerivAt (fun e ↦ A (u e)) (A' u') ε := hA.comp_hasDerivAt ε hu
  have h2 : HasDerivAt (fun e ↦ B (u e)) (B' u') ε := hB.comp_hasDerivAt ε hu
  have h3 := h1.add ((hasDerivAt_id ε).mul h2)
  have h0 : A' u' + ε * B' u' = 0 := by
    have := congrArg (fun l : E →L[ℝ] ℝ ↦ l u') hstat
    simpa using this
  have e : A' u' + (1 * B (u ε) + id ε * B' u') = B (u ε) := by
    simp only [id]
    linarith
  rw [e] at h3
  exact h3

/-- The slope of `ln C` against `ln ε` is `ε C' / C`. -/
theorem hasDerivAt_log_comp_exp {C : ℝ → ℝ} {C' ε : ℝ} (hε : 0 < ε) (hC : HasDerivAt C C' ε)
    (h0 : C ε ≠ 0) :
    HasDerivAt (fun t ↦ Real.log (C (Real.exp t))) (ε * C' / C ε) (Real.log ε) := by
  have hexp : HasDerivAt Real.exp ε (Real.log ε) := by
    simpa [Real.exp_log hε] using Real.hasDerivAt_exp (Real.log ε)
  have hC' : HasDerivAt C C' (Real.exp (Real.log ε)) := by rwa [Real.exp_log hε]
  have h1 : HasDerivAt (fun t ↦ C (Real.exp t)) (C' * ε) (Real.log ε) := hC'.comp _ hexp
  have h2 := h1.log (by rwa [Real.exp_log hε])
  rw [Real.exp_log hε] at h2
  have e : C' * ε / C ε = ε * C' / C ε := by ring
  rw [e] at h2
  exact h2

/-- **LIN-1** So the laminate's share of a solved section's field energy is the slope of the
logarithm of its capacitance against the logarithm of the laminate's permittivity:
`share = d ln C' / d ln ε_r`. -/
theorem hasDerivAt_log_energy_of_stationary {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    {A B : E → ℝ} {A' B' : E →L[ℝ] ℝ} {u : ℝ → E} {u' : E} {ε : ℝ} (hε : 0 < ε)
    (hA : HasFDerivAt A A' (u ε)) (hB : HasFDerivAt B B' (u ε)) (hu : HasDerivAt u u' ε)
    (hstat : A' + ε • B' = 0) (h0 : energy A B ε (u ε) ≠ 0) :
    HasDerivAt (fun t ↦ Real.log (energy A B (Real.exp t) (u (Real.exp t)))) (share A B ε (u ε))
      (Real.log ε) :=
  hasDerivAt_log_comp_exp (C := fun e ↦ energy A B e (u e)) hε
    (hasDerivAt_energy_of_stationary hA hB hu hstat) h0

/-! ### The smallest solve: two capacitors in series -/

/-- **LIN-1** The divider makes the series pair's energy stationary in its free node. -/
theorem hasDerivAt_stack_energy_divider {ca cl ε : ℝ} (h : ca + ε * cl ≠ 0) :
    HasDerivAt (fun x ↦ energy (stackAir ca) (stackLaminate cl) ε x) 0 (divider ca cl ε) := by
  have h1 : HasDerivAt (fun x : ℝ ↦ ca * x ^ 2) (ca * (2 * divider ca cl ε)) (divider ca cl ε) := by
    simpa using (hasDerivAt_pow 2 (divider ca cl ε)).const_mul ca
  have h2 : HasDerivAt (fun x : ℝ ↦ ε * (cl * (1 - x) ^ 2))
      (ε * (cl * (2 * (1 - divider ca cl ε) * -1))) (divider ca cl ε) := by
    have := ((hasDerivAt_id (divider ca cl ε)).const_sub 1).pow 2
    simpa using (this.const_mul cl).const_mul ε
  have h3 := h1.add h2
  have e : ca * (2 * divider ca cl ε) + ε * (cl * (2 * (1 - divider ca cl ε) * -1)) = 0 := by
    unfold divider
    field_simp
    ring
  rw [e] at h3
  exact h3

/-- **LIN-1** Its energy there is the series capacitance. -/
theorem energy_stack_divider {ca cl ε : ℝ} (h : ca + ε * cl ≠ 0) :
    energy (stackAir ca) (stackLaminate cl) ε (divider ca cl ε) = stackCapacitance ca cl ε := by
  unfold energy stackAir stackLaminate divider stackCapacitance
  field_simp
  ring

/-- **LIN-1** The laminate's share of it is `ca / (ca + ε cl)`: the laminate holds most of the
energy when it is the smaller capacitor. -/
theorem share_stack_divider {ca cl ε : ℝ} (h : ca + ε * cl ≠ 0) (hca : ca ≠ 0) (hcl : cl ≠ 0)
    (hε : ε ≠ 0) :
    share (stackAir ca) (stackLaminate cl) ε (divider ca cl ε) = ca / (ca + ε * cl) := by
  have key : 1 - ε * cl / (ca + ε * cl) = ca / (ca + ε * cl) := by
    field_simp
    ring
  unfold share
  rw [energy_stack_divider h]
  unfold stackLaminate divider stackCapacitance
  rw [key]
  field_simp

/-- **LIN-1** And the series capacitance changes with `ε` at the rate of the laminate's energy per
unit permittivity at the divider, so its logarithmic slope `ε C' / C` is that share. -/
theorem hasDerivAt_stackCapacitance {ca cl ε : ℝ} (h : ca + ε * cl ≠ 0) :
    HasDerivAt (stackCapacitance ca cl) (stackLaminate cl (divider ca cl ε)) ε := by
  have hn : HasDerivAt (fun e : ℝ ↦ ca * (e * cl)) (ca * cl) ε := by
    simpa using ((hasDerivAt_id ε).mul_const cl).const_mul ca
  have hd : HasDerivAt (fun e : ℝ ↦ ca + e * cl) cl ε := by
    simpa using ((hasDerivAt_id ε).mul_const cl).const_add ca
  have h3 := hn.div hd h
  have e : (ca * cl * (ca + ε * cl) - ca * (ε * cl) * cl) / (ca + ε * cl) ^ 2 =
      stackLaminate cl (divider ca cl ε) := by
    have key : 1 - ε * cl / (ca + ε * cl) = ca / (ca + ε * cl) := by
      field_simp
      ring
    unfold stackLaminate divider
    rw [key]
    field_simp
    ring
  rw [e] at h3
  exact h3

/-- **LIN-1** control: a potential that moves with `ε` without being the solve does not obey the
rule. With both capacitors 1 and the free node put at `x = ε`, at `ε = 1` the laminate holds no
energy and the energy still changes with `ε`, at the rate 2. -/
theorem not_hasDerivAt_energy_unsolved :
    ¬HasDerivAt (fun e ↦ energy (stackAir 1) (stackLaminate 1) e e) (stackLaminate 1 1) 1 := by
  have h : HasDerivAt (fun e : ℝ ↦ energy (stackAir 1) (stackLaminate 1) e e) 2 1 := by
    have h1 := ((hasDerivAt_pow 3 (1 : ℝ)).sub (hasDerivAt_pow 2 (1 : ℝ))).add (hasDerivAt_id 1)
    have e : (fun e : ℝ ↦ energy (stackAir 1) (stackLaminate 1) e e) =
        (fun x : ℝ ↦ x ^ 3 - x ^ 2) + id := by
      funext x
      simp only [energy, stackAir, stackLaminate, Pi.add_apply, id]
      ring
    rw [e]
    convert h1 using 1
    norm_num
  intro h'
  have := h.unique h'
  norm_num [stackLaminate] at this

end Array16.Line
