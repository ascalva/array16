module

public import Array16.Ring.Defs
public import Mathlib.Algebra.Field.GeomSum
public import Mathlib.Analysis.SpecialFunctions.Complex.Arg
public import Mathlib.Analysis.SpecialFunctions.Complex.Log
public import Mathlib.Tactic.FieldSimp
public import Mathlib.Tactic.Linarith
public import Mathlib.Tactic.LinearCombination
public import Mathlib.Tactic.Ring

/-!
# The ring resonator: RNG-1 to RNG-4

A closed guide is a ring of transmission line, and `array16/guide/ring.py` drives it so that the
volts travel round it one way: a travelling wave has the same size everywhere, so the well it makes
is as deep at every place on the path.

Two waves running opposite ways, `a e^(-iθ) + b e^(iθ)`, have a size between `|‖a‖ - ‖b‖|` and
`‖a‖ + ‖b‖`, and reach both; the size is the same everywhere exactly when one of them is absent,
and its ripple is the smaller amplitude over the larger (RNG-1). A travelling wave of peak volts
`V` on a ring of `m` guided waves stores `C V² L / 2` and so dissipates `π f C V² L / Q`, which is
`V² / (2 R)` with `R = Z0 Q / (2π m)` (RNG-2). `K` equal feeds phased for the wave of order `m`
launch the backward wave by the sum of the `K`-th roots of unity to the power `2m`: all of it when
`K` divides `2m` and none otherwise; two feeds a quarter guided wave apart at 1 and `-i` launch
none (RNG-3). A wave that runs a lap on a line losing `β / (2 Q)` keeps `e^(-2π m / Q)` of its
size squared (RNG-4).

Scope. These are laws of the formulas of `Array16/Ring/Defs.lean` over exact real and complex
numbers, proved by algebra. What they rest on is stated elsewhere and not proved here. That a
closed guide is a ring of uniform line whose volts are two such waves, that a feed launches both
with its current times the wave's phase at its place, that a travelling wave stores `C V² / 2` per
metre and that `Z0 = 1 / (v C)` are transmission-line theory, the model `ring.py` computes with;
whether a board's ring is that line is a matter for the full-wave solve and the bench. And the code
computes in binary64: the tests named `test_rng_N` attack the float implementation of each law.
-/

@[expose] public section

namespace Array16.Ring

open Complex

/-- A phase factor running the other way has size one. -/
theorem norm_exp_neg_mul_I (θ : ℝ) : ‖exp (-(θ * I))‖ = 1 := by
  rw [show -((θ : ℂ) * I) = ((-θ : ℝ) : ℂ) * I by push_cast; ring]
  exact norm_exp_ofReal_mul_I _

/-! ### Uniform means travelling -/

/-- **RNG-1** Two waves running opposite ways are nowhere larger than the sum of their
amplitudes. -/
theorem norm_wave_le (a b : ℂ) (θ : ℝ) : ‖wave a b θ‖ ≤ ‖a‖ + ‖b‖ := by
  unfold wave
  refine (norm_add_le _ _).trans (le_of_eq ?_)
  rw [norm_mul, norm_mul, norm_exp_neg_mul_I, norm_exp_ofReal_mul_I, mul_one, mul_one]

/-- **RNG-1** And nowhere smaller than their difference. -/
theorem abs_sub_le_norm_wave (a b : ℂ) (θ : ℝ) : |‖a‖ - ‖b‖| ≤ ‖wave a b θ‖ := by
  have h := abs_norm_sub_norm_le (a * exp (-(θ * I))) (-(b * exp (θ * I)))
  rw [norm_neg, norm_mul, norm_mul, norm_exp_neg_mul_I, norm_exp_ofReal_mul_I, mul_one, mul_one,
    sub_neg_eq_add] at h
  exact h

/-- Where the two waves are in phase they add. -/
theorem wave_polar (r s α β : ℝ) :
    wave (r * exp (α * I)) (s * exp (β * I)) ((α - β) / 2) =
      ((r : ℂ) + s) * exp (((α + β) / 2 : ℝ) * I) := by
  unfold wave
  rw [mul_assoc, mul_assoc, ← exp_add, ← exp_add,
    show (α : ℂ) * I + -(((α - β) / 2 : ℝ) * I) = ((α + β) / 2 : ℝ) * I by push_cast; ring,
    show (β : ℂ) * I + ((α - β) / 2 : ℝ) * I = ((α + β) / 2 : ℝ) * I by push_cast; ring]
  ring

/-- A quarter wave on they are opposed and subtract. -/
theorem wave_polar_quarter (r s α β : ℝ) :
    wave (r * exp (α * I)) (s * exp (β * I)) ((α - β) / 2 + Real.pi / 2) =
      ((r : ℂ) - s) * -I * exp (((α + β) / 2 : ℝ) * I) := by
  unfold wave
  rw [mul_assoc, mul_assoc, ← exp_add, ← exp_add,
    show (α : ℂ) * I + -(((α - β) / 2 + Real.pi / 2 : ℝ) * I) =
      ((α + β) / 2 : ℝ) * I + -(Real.pi / 2 * I) by push_cast; ring,
    show (β : ℂ) * I + ((α - β) / 2 + Real.pi / 2 : ℝ) * I =
      ((α + β) / 2 : ℝ) * I + Real.pi / 2 * I by push_cast; ring,
    exp_add, exp_add, exp_neg (Real.pi / 2 * I), exp_pi_div_two_mul_I, inv_I]
  ring

/-- **RNG-1** Somewhere on the ring they reach the sum. -/
theorem exists_norm_wave_eq_add (a b : ℂ) : ∃ θ, ‖wave a b θ‖ = ‖a‖ + ‖b‖ := by
  refine ⟨(arg a - arg b) / 2, ?_⟩
  have h := wave_polar ‖a‖ ‖b‖ (arg a) (arg b)
  rw [norm_mul_exp_arg_mul_I, norm_mul_exp_arg_mul_I] at h
  rw [h, norm_mul, norm_exp_ofReal_mul_I, mul_one, ← ofReal_add, norm_real,
    Real.norm_of_nonneg (by positivity)]

/-- **RNG-1** And somewhere the difference. -/
theorem exists_norm_wave_eq_abs_sub (a b : ℂ) : ∃ θ, ‖wave a b θ‖ = |‖a‖ - ‖b‖| := by
  refine ⟨(arg a - arg b) / 2 + Real.pi / 2, ?_⟩
  have h := wave_polar_quarter ‖a‖ ‖b‖ (arg a) (arg b)
  rw [norm_mul_exp_arg_mul_I, norm_mul_exp_arg_mul_I] at h
  rw [h, norm_mul, norm_mul, norm_exp_ofReal_mul_I, mul_one, norm_neg, norm_I, mul_one,
    ← ofReal_sub, norm_real, Real.norm_eq_abs]

/-- **RNG-1** So the ripple of the size, `(max - min) / (max + min)`, is the smaller amplitude over
the larger: the `|b/a|` of `ring.ring_purity`. -/
theorem ripple_extremes {x y : ℝ} (hx : 0 < x) (hyx : y ≤ x) : ripple (x + y) |x - y| = y / x := by
  unfold ripple
  rw [abs_of_nonneg (sub_nonneg.mpr hyx), div_eq_div_iff (by linarith) hx.ne']
  ring

/-- **RNG-1** A wave running one way has the same size everywhere. -/
theorem norm_wave_zero_right (a : ℂ) (θ : ℝ) : ‖wave a 0 θ‖ = ‖a‖ := by
  unfold wave
  rw [zero_mul, add_zero, norm_mul, norm_exp_neg_mul_I, mul_one]

/-- **RNG-1** And so has a wave running the other way. -/
theorem norm_wave_zero_left (b : ℂ) (θ : ℝ) : ‖wave 0 b θ‖ = ‖b‖ := by
  unfold wave
  rw [zero_mul, zero_add, norm_mul, norm_exp_ofReal_mul_I, mul_one]

/-- The size squared of two waves running opposite ways, at phase `θ`. -/
theorem sq_norm_wave (a b : ℂ) (θ : ℝ) :
    ‖wave a b θ‖ ^ 2 =
      ‖a‖ ^ 2 + ‖b‖ ^ 2 +
        2 * ((a * (starRingEnd ℂ) b).re * Real.cos (2 * θ) +
          (a * (starRingEnd ℂ) b).im * Real.sin (2 * θ)) := by
  have hw : wave a b θ =
      a * (Real.cos θ - Real.sin θ * I) + b * (Real.cos θ + Real.sin θ * I) := by
    unfold wave
    rw [show -((θ : ℂ) * I) = ((-θ : ℝ) : ℂ) * I by push_cast; ring, exp_mul_I, exp_mul_I,
      ← ofReal_cos, ← ofReal_sin, ← ofReal_cos, ← ofReal_sin, Real.cos_neg, Real.sin_neg]
    push_cast
    ring
  rw [hw, Complex.sq_norm, Complex.sq_norm, Complex.sq_norm, normSq_apply, normSq_apply,
    normSq_apply, Real.cos_two_mul, Real.sin_two_mul]
  simp only [add_re, add_im, mul_re, mul_im, sub_re, sub_im, ofReal_re, ofReal_im, I_re, I_im,
    conj_re, conj_im]
  linear_combination
    (a.re ^ 2 + a.im ^ 2 + b.re ^ 2 + b.im ^ 2 - 2 * (a.re * b.re + a.im * b.im)) *
      Real.cos_sq_add_sin_sq θ

/-- **RNG-1** And only such a wave: if the size is the same everywhere, one of the two waves is
absent. So a uniform well is a travelling wave, and the target of a drive is one wave. -/
theorem eq_zero_or_eq_zero_of_norm_wave_eq {a b : ℂ} (h : ∀ θ, ‖wave a b θ‖ = ‖wave a b 0‖) :
    a = 0 ∨ b = 0 := by
  have h0 := sq_norm_wave a b 0
  have h1 := sq_norm_wave a b (Real.pi / 2)
  have h2 := sq_norm_wave a b (Real.pi / 4)
  rw [h (Real.pi / 2)] at h1
  rw [h (Real.pi / 4)] at h2
  rw [show (2 : ℝ) * 0 = 0 by ring, Real.cos_zero, Real.sin_zero] at h0
  rw [show 2 * (Real.pi / 2) = Real.pi by ring, Real.cos_pi, Real.sin_pi] at h1
  rw [show 2 * (Real.pi / 4) = Real.pi / 2 by ring, Real.cos_pi_div_two, Real.sin_pi_div_two] at h2
  have hre : (a * (starRingEnd ℂ) b).re = 0 := by linarith
  have him : (a * (starRingEnd ℂ) b).im = 0 := by linarith
  have hc : a * (starRingEnd ℂ) b = 0 := Complex.ext (by simpa using hre) (by simpa using him)
  rcases mul_eq_zero.mp hc with ha | hb
  · exact Or.inl ha
  · exact Or.inr (by simpa using hb)

/-- **RNG-1** control: two equal waves stand. A quarter wave from their crest the volts vanish, a
ripple of one. -/
theorem wave_self_pi_div_two (a : ℂ) : wave a a (Real.pi / 2) = 0 := by
  unfold wave
  rw [exp_neg, show ((Real.pi / 2 : ℝ) : ℂ) = Real.pi / 2 by push_cast; ring,
    exp_pi_div_two_mul_I, inv_I]
  ring

/-! ### The shunt resistance of a travelling wave -/

/-- **RNG-2** A ring carrying a travelling wave of peak volts `V` dissipates `π f C V² L / Q`. -/
theorem dissipated_stored (f Q C V L : ℝ) :
    dissipated f Q (stored C V L) = Real.pi * f * C * V ^ 2 * L / Q := by
  unfold dissipated stored
  ring

/-- **RNG-2** That is `V² / (2 R)` with `R = Z0 Q / (2π m)`, for a ring of `m` guided waves
(`L f = m v`) of a line whose impedance is `Z0 = 1 / (v C)`. -/
theorem sq_div_dissipated_eq_shunt {f Q C V L m v Z0 : ℝ} (hres : L * f = m * v)
    (hZ : Z0 * (v * C) = 1) (hf : f ≠ 0) (hm : m ≠ 0) (hQ : Q ≠ 0) (hV : V ≠ 0) :
    V ^ 2 / (2 * dissipated f Q (stored C V L)) = shunt Z0 Q m := by
  have hv : v ≠ 0 := by
    rintro rfl
    simp at hZ
  have hC : C ≠ 0 := by
    rintro rfl
    simp at hZ
  have hL : L = m * v / f := by
    field_simp
    linarith
  have hZ0 : Z0 = 1 / (v * C) := by
    field_simp
    linarith
  rw [dissipated_stored, hL, hZ0]
  unfold shunt
  have hpi := Real.pi_ne_zero
  field_simp

/-- **RNG-2** control: the order matters. A ring of twice as many guided waves of the same line
has half the shunt resistance. -/
theorem shunt_two_mul_ne {Z0 Q m : ℝ} (h : shunt Z0 Q m ≠ 0) :
    shunt Z0 Q (2 * m) ≠ shunt Z0 Q m := by
  have e : shunt Z0 Q (2 * m) = shunt Z0 Q m / 2 := by
    unfold shunt
    rw [show 2 * Real.pi * (2 * m) = 2 * Real.pi * m * 2 by ring, ← div_div]
  intro heq
  rw [e] at heq
  exact h (by linarith)

/-! ### Which feeds reach the backward wave -/

/-- The powers of a `K`-th root of unity sum to `K` when it is one and to nothing otherwise. -/
theorem sum_exp_eq_ite {K : ℕ} (hK : K ≠ 0) (n : ℤ) :
    ∑ j ∈ Finset.range K, exp (2 * Real.pi * I * (n * j / K)) =
      if (K : ℤ) ∣ n then (K : ℂ) else 0 := by
  have hK' : (K : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hK
  obtain ⟨ζ, hζ⟩ : ∃ ζ : ℂ, ζ = exp (2 * Real.pi * I * (n / K)) := ⟨_, rfl⟩
  have hpow : ∀ j : ℕ, exp (2 * Real.pi * I * (n * j / K)) = ζ ^ j := fun j ↦ by
    rw [hζ, ← exp_nat_mul]
    congr 1
    ring
  have hone : ζ ^ K = 1 := by
    rw [hζ, ← exp_nat_mul,
      show (K : ℂ) * (2 * Real.pi * I * (n / K)) = n * (2 * Real.pi * I) by field_simp]
    exact exp_int_mul_two_pi_mul_I n
  simp_rw [hpow]
  split_ifs with hd
  · obtain ⟨k, rfl⟩ := hd
    have h1 : ζ = 1 := by
      rw [hζ, show 2 * Real.pi * I * (((K * k : ℤ) : ℂ) / K) = k * (2 * Real.pi * I) by
        push_cast
        field_simp]
      exact exp_int_mul_two_pi_mul_I k
    simp [h1]
  · have hne : ζ ≠ 1 := by
      intro h1
      rw [hζ] at h1
      obtain ⟨k, hk⟩ := exp_eq_one_iff.mp h1
      have hpi : (2 * Real.pi * I : ℂ) ≠ 0 := by
        simp [Real.pi_ne_zero]
      have h2 : (n : ℂ) / K = k :=
        mul_left_cancel₀ hpi (by rw [hk]; ring)
      have h3 : (n : ℂ) = ((K * k : ℤ) : ℂ) := by
        rw [(div_eq_iff hK').mp h2]
        push_cast
        ring
      exact hd ⟨k, by exact_mod_cast h3⟩
    rw [geom_sum_eq hne, hone, sub_self, zero_div]

/-- **RNG-3** `K` equal feeds phased for the travelling wave of order `m` all add in it. -/
theorem launch_travelling (K : ℕ) (m : ℤ) : launch K (travelling K m) m = K := by
  unfold launch travelling
  rw [Finset.sum_congr rfl fun j _ ↦ by rw [← exp_add, neg_add_cancel, exp_zero]]
  simp

/-- **RNG-3** They reach the backward wave, of order `-m`, only when their count divides `2m`, and
then as strongly as the wanted one; otherwise not at all. -/
theorem launch_travelling_neg {K : ℕ} (hK : K ≠ 0) (m : ℤ) :
    launch K (travelling K m) (-m) = if (K : ℤ) ∣ 2 * m then (K : ℂ) else 0 := by
  have hterm : ∀ j : ℕ,
      exp (-(2 * Real.pi * I * (m * j / K))) * exp (2 * Real.pi * I * ((-m : ℤ) * j / K)) =
        exp (2 * Real.pi * I * ((-(2 * m) : ℤ) * j / K)) := fun j ↦ by
    rw [← exp_add]
    congr 1
    push_cast
    ring
  unfold launch travelling
  rw [Finset.sum_congr rfl fun j _ ↦ hterm j, sum_exp_eq_ite hK]
  simp only [dvd_neg]

/-- **RNG-3** So a count that does not divide `2m` leaves the backward wave alone: three feeds at
`m = 1` or `2`, or any count above `2m`. -/
theorem launch_travelling_neg_of_not_dvd {K : ℕ} (hK : K ≠ 0) {m : ℤ} (h : ¬(K : ℤ) ∣ 2 * m) :
    launch K (travelling K m) (-m) = 0 := by
  simp [launch_travelling_neg hK, h]

/-- **RNG-3** Two feeds a quarter guided wave apart, the second at `-i` times the first (the 90
degree hybrid), launch nothing backward. -/
theorem pairBackward_neg_I : pairBackward (-I) = 0 := by
  unfold pairBackward
  rw [exp_neg, exp_pi_div_two_mul_I, inv_I]
  linear_combination I_sq

/-- **RNG-3** And add in the wanted wave. -/
theorem pairForward_neg_I : pairForward (-I) = 2 := by
  unfold pairForward
  rw [exp_pi_div_two_mul_I]
  linear_combination -I_sq

/-- **RNG-3** control: a count that divides `2m` drives both senses alike, a standing wave: two
equal feeds half a ring apart at `m = 1`. -/
theorem launch_travelling_neg_two : launch 2 (travelling 2 1) (-1) = 2 := by
  rw [launch_travelling_neg (by norm_num)]
  norm_num

/-- **RNG-3** control: the same pair of feeds a quarter wave apart, fed in phase, does launch the
backward wave. -/
theorem pairBackward_one_ne : pairBackward 1 ≠ 0 := by
  unfold pairBackward
  rw [exp_neg, exp_pi_div_two_mul_I, inv_I]
  intro h
  have := congrArg Complex.re h
  simp at this

/-! ### The loss of a lap -/

/-- The size squared of a wave after a length `s` of lossy line. -/
theorem sq_norm_decayed (a : ℂ) (α β s : ℝ) :
    ‖decayed a α β s‖ ^ 2 = ‖a‖ ^ 2 * Real.exp (-(2 * α * s)) := by
  unfold decayed
  have hre : (-((α + β * I) * s) : ℂ).re = -(α * s) := by simp
  rw [norm_mul, mul_pow, norm_exp, hre, sq (Real.exp _), ← Real.exp_add,
    show -(α * s) + -(α * s) = -(2 * α * s) by ring]

/-- **RNG-4** On a line of quality `Q`, losing `β / (2 Q)`, a wave that runs a lap of `m` guided
waves (`β L = 2π m`) keeps `e^(-2π m / Q)` of its size squared, and so of the depth of its well. -/
theorem sq_norm_decayed_lap {β Q L m : ℝ} (hQ : Q ≠ 0) (hlap : β * L = 2 * Real.pi * m) (a : ℂ) :
    ‖decayed a (β / (2 * Q)) β L‖ ^ 2 = ‖a‖ ^ 2 * Real.exp (-(2 * Real.pi * m / Q)) := by
  have e : 2 * (β / (2 * Q)) * L = 2 * Real.pi * m / Q := by
    rw [← hlap]
    field_simp
  rw [sq_norm_decayed, e]

/-- **RNG-4** So fed at one place the depth falls by `1 - e^(-2π m / Q)` round the lap. -/
theorem one_sub_sq_norm_decayed_lap_div {β Q L m : ℝ} (hQ : Q ≠ 0)
    (hlap : β * L = 2 * Real.pi * m) {a : ℂ} (ha : a ≠ 0) :
    1 - ‖decayed a (β / (2 * Q)) β L‖ ^ 2 / ‖a‖ ^ 2 = 1 - Real.exp (-(2 * Real.pi * m / Q)) := by
  have h : ‖a‖ ^ 2 ≠ 0 := pow_ne_zero 2 (norm_ne_zero_iff.mpr ha)
  rw [sq_norm_decayed_lap hQ hlap, mul_div_cancel_left₀ _ h]

/-- **RNG-4** control: at a finite `Q` the lap does cost depth; the wave fed at one place is not
uniform. -/
theorem sq_norm_decayed_lap_lt {β Q L m : ℝ} (hQ : 0 < Q) (hm : 0 < m)
    (hlap : β * L = 2 * Real.pi * m) {a : ℂ} (ha : a ≠ 0) :
    ‖decayed a (β / (2 * Q)) β L‖ ^ 2 < ‖a‖ ^ 2 := by
  have h : 0 < ‖a‖ ^ 2 := by positivity
  have hx : Real.exp (-(2 * Real.pi * m / Q)) < 1 := by
    have hneg : -(2 * Real.pi * m / Q) < 0 := by
      have := Real.pi_pos
      rw [neg_lt_zero]
      positivity
    have := Real.exp_lt_exp.mpr hneg
    rwa [Real.exp_zero] at this
  rw [sq_norm_decayed_lap hQ.ne' hlap]
  nlinarith

end Array16.Ring
