module

public import Array16.Electrode.Defs
public import Mathlib.Analysis.Calculus.Deriv.Add
public import Mathlib.Analysis.Calculus.Deriv.Inv
public import Mathlib.Analysis.SpecialFunctions.Sqrt
public import Mathlib.Tactic.FieldSimp
public import Mathlib.Tactic.LinearCombination
public import Mathlib.Tactic.Ring

/-!
# Electrodes printed on a board: ELE-1

Two rails on one board, either side of a grounded centre strip, hold a field null above the centre
line. `RailPair` in `array16/particles/electrodes.py` puts it at the geometric mean of the rails'
edge distances, `h = √(c1 c2)`. This is why: the field is a rational function of `z = x + i y`
whose numerator is `2 (c2 - c1) (z² + c1 c2)`, so it vanishes exactly where `z² = -(c1 c2)`, which
above the board is `z = i √(c1 c2)`. There its derivative is `i` times
`4 (c2 - c1) / (√(c1 c2) (c1 + c2)²)`, which is `4 b / (h (a + b)²)` for a gap `a` and rails `b`:
the gradient of `RailPair.gradient` over `V / π`.

Scope. These are laws of the rational function `railField` over exact complex numbers, proved by
algebra and the rules of differentiation. What they rest on is stated elsewhere and not proved
here. That the rails' field is `(i V / π) * railField` is the gapless surface-electrode model
(Oliveira and Miranda 2001; Wesenberg 2008), in which an electrode at `V` in a grounded plane has
potential `V / π` times the angle it subtends: the tests check the potential harmonic, the volts on
the rails, and the field against its gradient. And the code computes in binary64, not in `ℂ`: the
tests named `test_ele_1` attack the float implementation.
-/

@[expose] public section

namespace Array16.Electrode

open Complex

/-- Off the rails' edges the field is `2 (c2 - c1) (z² + c1 c2)` over the product of the four
distances to the edges. -/
theorem railField_eq {c1 c2 : ℝ} {z : ℂ} (h1 : z ≠ c1) (h1' : z ≠ -c1) (h2 : z ≠ c2)
    (h2' : z ≠ -c2) :
    railField c1 c2 z =
      2 * (c2 - c1) * (z ^ 2 + c1 * c2) / ((z - c1) * (z + c1) * (z - c2) * (z + c2)) := by
  have a1 : z - c1 ≠ 0 := sub_ne_zero.mpr h1
  have a2 : z + c1 ≠ 0 := fun h ↦ h1' (eq_neg_of_add_eq_zero_left h)
  have a3 : z - c2 ≠ 0 := sub_ne_zero.mpr h2
  have a4 : z + c2 ≠ 0 := fun h ↦ h2' (eq_neg_of_add_eq_zero_left h)
  unfold railField
  field_simp
  ring

/-- **ELE-1** Off the rails' edges, the field of two rails of some width vanishes exactly where
`z² = -(c1 c2)`. -/
theorem railField_eq_zero_iff {c1 c2 : ℝ} (hc : c1 ≠ c2) {z : ℂ} (h1 : z ≠ c1) (h1' : z ≠ -c1)
    (h2 : z ≠ c2) (h2' : z ≠ -c2) :
    railField c1 c2 z = 0 ↔ z ^ 2 = -(c1 * c2) := by
  have a1 : z - c1 ≠ 0 := sub_ne_zero.mpr h1
  have a2 : z + c1 ≠ 0 := fun h ↦ h1' (eq_neg_of_add_eq_zero_left h)
  have a3 : z - c2 ≠ 0 := sub_ne_zero.mpr h2
  have a4 : z + c2 ≠ 0 := fun h ↦ h2' (eq_neg_of_add_eq_zero_left h)
  have hc' : (c2 : ℂ) - c1 ≠ 0 := sub_ne_zero.mpr (fun h ↦ hc (ofReal_injective h).symm)
  rw [railField_eq h1 h1' h2 h2', div_eq_zero_iff]
  constructor
  · rintro (h | h)
    · rcases mul_eq_zero.mp h with h | h
      · rcases mul_eq_zero.mp h with h | h
        · norm_num at h
        · exact absurd h hc'
      · linear_combination h
    · exact absurd h (by simp [a1, a2, a3, a4])
  · intro h
    left
    rw [h]
    ring

/-- A point above the board is not on it: `i s` with `s ≠ 0` is no real number. -/
theorem I_mul_ne_ofReal {s : ℝ} (hs : s ≠ 0) (c : ℝ) : I * s ≠ c := by
  intro h
  have him : (I * (s : ℂ)).im = s := by simp
  rw [h, ofReal_im] at him
  exact hs him.symm

/-- **ELE-1** The null stands above the centre line at the geometric mean of the rails' edge
distances: the field vanishes at `z = i √(c1 c2)`. -/
theorem railField_I_mul_sqrt {c1 c2 : ℝ} (h1 : 0 < c1) (h2 : 0 < c2) :
    railField c1 c2 (I * (√(c1 * c2) : ℝ)) = 0 := by
  have hs : √(c1 * c2) ≠ 0 := (Real.sqrt_pos.mpr (mul_pos h1 h2)).ne'
  have e : (I * (√(c1 * c2) : ℝ)) ^ 2 = -((c1 : ℂ) * c2) := by
    rw [mul_pow, I_sq, ← ofReal_pow, Real.sq_sqrt (mul_pos h1 h2).le]
    push_cast
    ring
  rw [railField_eq (I_mul_ne_ofReal hs c1) (by simpa using I_mul_ne_ofReal hs (-c1))
    (I_mul_ne_ofReal hs c2) (by simpa using I_mul_ne_ofReal hs (-c2)), e]
  simp

/-- **ELE-1** For a gap `a` and rails `b` wide, edges at `a/2` and `a/2 + b`, that height is
`½ √(a (a + 2b))` (`RailPair.null_height`). -/
theorem sqrt_edges (a b : ℝ) : √(a / 2 * (a / 2 + b)) = 1 / 2 * √(a * (a + 2 * b)) := by
  rw [show a / 2 * (a / 2 + b) = (1 / 2) ^ 2 * (a * (a + 2 * b)) by ring,
    Real.sqrt_mul (by positivity), Real.sqrt_sq (by norm_num)]

/-- Off the rails' edges the field's derivative is the sum of the four edges' terms. -/
theorem hasDerivAt_railField {c1 c2 : ℝ} {z : ℂ} (h1 : z ≠ c1) (h1' : z ≠ -c1) (h2 : z ≠ c2)
    (h2' : z ≠ -c2) :
    HasDerivAt (railField c1 c2)
      (-1 / (z - c2) ^ 2 + 1 / (z - c1) ^ 2 - 1 / (z + c1) ^ 2 + 1 / (z + c2) ^ 2) z := by
  have a1 : z - c1 ≠ 0 := sub_ne_zero.mpr h1
  have a2 : z + c1 ≠ 0 := fun h ↦ h1' (eq_neg_of_add_eq_zero_left h)
  have a3 : z - c2 ≠ 0 := sub_ne_zero.mpr h2
  have a4 : z + c2 ≠ 0 := fun h ↦ h2' (eq_neg_of_add_eq_zero_left h)
  have d1 := (HasDerivAt.sub_const (c2 : ℂ) (hasDerivAt_id' z)).inv a3
  have d2 := (HasDerivAt.sub_const (c1 : ℂ) (hasDerivAt_id' z)).inv a1
  have d3 := (HasDerivAt.add_const (c1 : ℂ) (hasDerivAt_id' z)).inv a2
  have d4 := (HasDerivAt.add_const (c2 : ℂ) (hasDerivAt_id' z)).inv a4
  convert ((d1.sub d2).add d3).sub d4 using 1
  · funext y
    simp only [railField, one_div, Pi.sub_apply, Pi.add_apply, Pi.inv_apply]
  · ring

/-- Where `z² = -(c1 c2)`, that derivative is `4 z (c2 - c1) / (c1 c2 (c1 + c2)²)`. -/
theorem railField_deriv_of_sq {c1 c2 : ℝ} (h1 : 0 < c1) (h2 : 0 < c2) {z : ℂ}
    (hz : z ^ 2 = -((c1 : ℂ) * c2)) (e1 : z ≠ c1) (e1' : z ≠ -c1) (e2 : z ≠ c2) (e2' : z ≠ -c2) :
    -1 / (z - c2) ^ 2 + 1 / (z - c1) ^ 2 - 1 / (z + c1) ^ 2 + 1 / (z + c2) ^ 2 =
      4 * z * (c2 - c1) / (c1 * c2 * (c1 + c2) ^ 2) := by
  have a1 : z - c1 ≠ 0 := sub_ne_zero.mpr e1
  have a2 : z + c1 ≠ 0 := fun h ↦ e1' (eq_neg_of_add_eq_zero_left h)
  have a3 : z - c2 ≠ 0 := sub_ne_zero.mpr e2
  have a4 : z + c2 ≠ 0 := fun h ↦ e2' (eq_neg_of_add_eq_zero_left h)
  have b1 : (c1 : ℂ) ≠ 0 := ofReal_ne_zero.mpr h1.ne'
  have b2 : (c2 : ℂ) ≠ 0 := ofReal_ne_zero.mpr h2.ne'
  have b3 : (c1 : ℂ) + c2 ≠ 0 := by exact_mod_cast (add_pos h1 h2).ne'
  have s1 : z ^ 2 - (c1 : ℂ) ^ 2 ≠ 0 := by
    rw [show z ^ 2 - (c1 : ℂ) ^ 2 = (z - c1) * (z + c1) by ring]
    exact mul_ne_zero a1 a2
  have s2 : z ^ 2 - (c2 : ℂ) ^ 2 ≠ 0 := by
    rw [show z ^ 2 - (c2 : ℂ) ^ 2 = (z - c2) * (z + c2) by ring]
    exact mul_ne_zero a3 a4
  have e : -1 / (z - c2) ^ 2 + 1 / (z - c1) ^ 2 - 1 / (z + c1) ^ 2 + 1 / (z + c2) ^ 2 =
      4 * z * c1 / (z ^ 2 - (c1 : ℂ) ^ 2) ^ 2 - 4 * z * c2 / (z ^ 2 - (c2 : ℂ) ^ 2) ^ 2 := by
    field_simp
    ring
  rw [e, hz, show -((c1 : ℂ) * c2) - (c1 : ℂ) ^ 2 = -((c1 : ℂ) * (c1 + c2)) by ring,
    show -((c1 : ℂ) * c2) - (c2 : ℂ) ^ 2 = -((c2 : ℂ) * (c1 + c2)) by ring]
  field_simp

/-- **ELE-1** At the null `i √(c1 c2)` the field's derivative is `i` times the real number
`4 (c2 - c1) / (√(c1 c2) (c1 + c2)²)`. Times `V / π` its size is the gradient of
`RailPair.gradient`, the same whichever way the null is left, as a derivative is. -/
theorem hasDerivAt_railField_I_mul_sqrt {c1 c2 : ℝ} (h1 : 0 < c1) (h2 : 0 < c2) :
    HasDerivAt (railField c1 c2)
      (I * ((4 * (c2 - c1) / (√(c1 * c2) * (c1 + c2) ^ 2) : ℝ) : ℂ)) (I * (√(c1 * c2) : ℝ)) := by
  have hs : √(c1 * c2) ≠ 0 := (Real.sqrt_pos.mpr (mul_pos h1 h2)).ne'
  have hsq : (√(c1 * c2) : ℂ) ^ 2 = (c1 : ℂ) * c2 := by
    rw [← ofReal_pow, Real.sq_sqrt (mul_pos h1 h2).le]
    push_cast
    ring
  have e1 := I_mul_ne_ofReal hs c1
  have e1' : I * (√(c1 * c2) : ℂ) ≠ -(c1 : ℂ) := by simpa using I_mul_ne_ofReal hs (-c1)
  have e2 := I_mul_ne_ofReal hs c2
  have e2' : I * (√(c1 * c2) : ℂ) ≠ -(c2 : ℂ) := by simpa using I_mul_ne_ofReal hs (-c2)
  have hz : (I * (√(c1 * c2) : ℂ)) ^ 2 = -((c1 : ℂ) * c2) := by
    rw [mul_pow, I_sq, hsq]
    ring
  convert hasDerivAt_railField e1 e1' e2 e2' using 1
  rw [railField_deriv_of_sq h1 h2 hz e1 e1' e2 e2', ← hsq]
  have hs' : (√(c1 * c2) : ℂ) ≠ 0 := ofReal_ne_zero.mpr hs
  have b3 : (c1 : ℂ) + c2 ≠ 0 := by exact_mod_cast (add_pos h1 h2).ne'
  push_cast
  field_simp

/-- **ELE-1** For a gap `a` and rails `b` wide, edges at `a/2` and `a/2 + b`, that size is
`4 b / (h (a + b)²)` with `h = ½ √(a (a + 2b))` the null's height: `RailPair.gradient` over
`V / π`. -/
theorem gradient_edges (a b : ℝ) :
    4 * ((a / 2 + b) - a / 2) / (√(a / 2 * (a / 2 + b)) * (a / 2 + (a / 2 + b)) ^ 2) =
      4 * b / (1 / 2 * √(a * (a + 2 * b)) * (a + b) ^ 2) := by
  rw [sqrt_edges]
  ring

/-- **ELE-1** control: twice as high as the null the field is not zero. -/
theorem railField_two_mul_I_mul_sqrt_ne {c1 c2 : ℝ} (h1 : 0 < c1) (h2 : 0 < c2) (hc : c1 ≠ c2) :
    railField c1 c2 (I * ((2 * √(c1 * c2) : ℝ) : ℂ)) ≠ 0 := by
  have hs : 2 * √(c1 * c2) ≠ 0 :=
    mul_ne_zero two_ne_zero (Real.sqrt_pos.mpr (mul_pos h1 h2)).ne'
  have e : (I * ((2 * √(c1 * c2) : ℝ) : ℂ)) ^ 2 = -(4 * ((c1 : ℂ) * c2)) := by
    rw [mul_pow, I_sq, ← ofReal_pow, mul_pow, Real.sq_sqrt (mul_pos h1 h2).le]
    push_cast
    ring
  rw [Ne, railField_eq_zero_iff hc (I_mul_ne_ofReal hs c1) (by simpa using I_mul_ne_ofReal hs (-c1))
    (I_mul_ne_ofReal hs c2) (by simpa using I_mul_ne_ofReal hs (-c2)), e]
  intro h
  have h0 : (c1 : ℂ) * c2 = 0 := by linear_combination (-1 / 3 : ℂ) * h
  have : c1 * c2 = 0 := by exact_mod_cast h0
  exact (mul_pos h1 h2).ne' this

end Array16.Electrode
