module

public import Mathlib.Basic.Real.Basic
public import Mathlib.Algebra.BigOperators.Fin

/-!
# Bending of a closed curve: CUR-2, the definitions

The definitions the law of `Array16/Curve.lean` is stated over. They mirror `ClosedCurve` in
`array16/particles/curves.py`: a curve is `n` points a step `h` apart in arc length, with the
curvature `κ i` at each, so its length is `n * h`.
-/

@[expose] public section

namespace Array16.Curve

/-- The curvature summed round the curve (`ClosedCurve.total_curvature`). -/
def totalCurvature {n : ℕ} (h : ℝ) (κ : Fin n → ℝ) : ℝ :=
  h * ∑ i, κ i

/-- The squared curvature summed round the curve (`ClosedCurve.bending`). -/
def bending {n : ℕ} (h : ℝ) (κ : Fin n → ℝ) : ℝ :=
  h * ∑ i, κ i ^ 2

end Array16.Curve
