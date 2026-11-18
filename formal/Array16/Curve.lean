module

public import Array16.Curve.Defs
public import Mathlib.Algebra.Order.Chebyshev
public import Mathlib.Tactic.GCongr
public import Mathlib.Tactic.Positivity
public import Mathlib.Tactic.Ring

/-!
# Bending of a closed curve: CUR-2

A guide's axis bends, and what the bends cost an electron goes as the squared curvature summed round
the curve, `ClosedCurve.bending` in `array16/particles/curves.py`. This is the inequality that ties
it to the total curvature: `bending ≥ total curvature ^ 2 / length`, by Cauchy and Schwarz.

It turns two theorems of the geometry of curves into floors on the bending. Every closed curve has
total curvature at least `2π` (Fenchel), and every knotted one more than `4π` (Fáry, Milnor). So a
knot of length `L` bends by more than `16π² / L`, four times what a circle of that length does,
whatever its shape. Neither theorem is in Mathlib, and neither is proved here: each enters as the
hypothesis `c ≤ totalCurvature h κ` of `le_bending_of_le_totalCurvature`, and the ledger cites it
(CUR-3).
-/

@[expose] public section

namespace Array16.Curve

/-- **CUR-2** The total curvature squared is at most the length times the bending. -/
theorem sq_totalCurvature_le {n : ℕ} (h : ℝ) (κ : Fin n → ℝ) :
    totalCurvature h κ ^ 2 ≤ n * h * bending h κ := by
  have key : (∑ i, κ i) ^ 2 ≤ n * ∑ i, κ i ^ 2 := by
    simpa using sq_sum_le_card_mul_sum_sq (s := Finset.univ) (f := κ)
  have hh2 : 0 ≤ h ^ 2 := by positivity
  unfold totalCurvature bending
  calc (h * ∑ i, κ i) ^ 2 = h ^ 2 * (∑ i, κ i) ^ 2 := by ring
    _ ≤ h ^ 2 * (n * ∑ i, κ i ^ 2) := mul_le_mul_of_nonneg_left key hh2
    _ = n * h * (h * ∑ i, κ i ^ 2) := by ring

/-- **CUR-2** So a curve whose total curvature is at least `c` bends by at least `c ^ 2` over its
length: with `c = 2π` for any closed curve, and `c = 4π` for a knotted one. -/
theorem le_bending_of_le_totalCurvature {n : ℕ} {h c : ℝ} (hh : 0 < h) (hn : 0 < n) (hc : 0 ≤ c)
    (κ : Fin n → ℝ) (hκ : c ≤ totalCurvature h κ) : c ^ 2 / (n * h) ≤ bending h κ := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have hL : (0 : ℝ) < n * h := mul_pos hn' hh
  rw [div_le_iff₀ hL]
  calc c ^ 2 ≤ totalCurvature h κ ^ 2 := pow_le_pow_left₀ hc hκ 2
    _ ≤ n * h * bending h κ := sq_totalCurvature_le h κ
    _ = bending h κ * (n * h) := by ring

/-- **CUR-2** control: the inequality is sharp. A curve of constant curvature, a circle, attains
it. -/
theorem sq_totalCurvature_const {n : ℕ} (h k : ℝ) :
    totalCurvature h (fun _ : Fin n ↦ k) ^ 2 = n * h * bending h (fun _ : Fin n ↦ k) := by
  simp only [totalCurvature, bending, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
    nsmul_eq_mul]
  ring

end Array16.Curve
