module

public import Array16.Guide.Defs
public import Mathlib.Analysis.Convex.Function
public import Mathlib.Tactic.FieldSimp
public import Mathlib.Tactic.LinearCombination
public import Mathlib.Tactic.NormNum
public import Mathlib.Tactic.Ring

/-!
# The lap map of a guide field: GDE-1 to GDE-7

`GuideOptics` in `array16/particles.py` predicts where an electron goes in a magnetic guide field
that follows a closed curve, a ring or a knot, from three numbers per curve: the turn `α` of the
transverse plane per lap, the drift `G` per unit gyroradius and the wander `H` of the field line.
One lap is `z ↦ exp (i α) * (z + ρ G + H)`.

With a turn (`r ≠ 1`, `‖r‖ = 1`) there is exactly one closed orbit, the offsets at the start of
each lap lie on a circle about it, and the orbit is bounded for every number of laps. Without a
turn (`r = 1`, a planar ring) every lap adds the same step and the orbit leaves every bounded
region. The reach is convex in `ρ`, so the values of `ρ` that stay inside a wall form one unbroken
range, which `particles.window` reports by its two ends.

The proofs are algebra. The map is conjugate to a turn about its fixed point, the offset at the
start of each lap is affine in `ρ`, and the length of an affine function is convex.

Scope. These are laws of the lap map over exact complex numbers. Three things they rest on are
stated elsewhere and not proved here. That a guiding centre obeys the lap map is first-order drift
theory, a hypothesis that full-orbit tracking tests (`experiments/guide/`). That the turn is `2π`
times the curve's writhe modulo `2π` is the Călugăreanu-White-Fuller theorem, checked numerically
by two independent computations (`tests/test_particles.py`). And the code computes in binary64, not
in `ℂ`: the tests named `test_gde_N` attack the float implementation of each law.
-/

@[expose] public section

namespace Array16.Guide

open Function

/-! ### With a turn: one closed orbit, and every orbit a circle about it -/

/-- **GDE-1** The closed orbit is a fixed point of the lap map. -/
theorem lap_closedOrbit {r : ℂ} (hr : r ≠ 1) (d : ℂ) :
    lap r d (closedOrbit r d) = closedOrbit r d := by
  have h : (1 : ℂ) - r ≠ 0 := sub_ne_zero.mpr (Ne.symm hr)
  unfold lap closedOrbit
  field_simp
  ring

/-- **GDE-2** It is the only one: if the plane turns, an orbit that closes after one lap is the
closed orbit. -/
theorem eq_closedOrbit_of_lap_eq {r : ℂ} (hr : r ≠ 1) {d z : ℂ} (hz : lap r d z = z) :
    z = closedOrbit r d := by
  have h : (1 : ℂ) - r ≠ 0 := sub_ne_zero.mpr (Ne.symm hr)
  unfold lap at hz
  rw [closedOrbit, ← mul_div_assoc, eq_div_iff h]
  linear_combination -hz

/-- One lap turns the offset from the closed orbit by `r`. -/
theorem lap_sub_closedOrbit {r : ℂ} (hr : r ≠ 1) (d z : ℂ) :
    lap r d z - closedOrbit r d = r * (z - closedOrbit r d) := by
  have h := lap_closedOrbit hr d
  unfold lap at h ⊢
  linear_combination h

/-- **GDE-2** What a lap moves an offset by is `r - 1` times its offset from the closed orbit, so
it is nothing only there. -/
theorem lap_sub_self {r : ℂ} (hr : r ≠ 1) (d z : ℂ) :
    lap r d z - z = (r - 1) * (z - closedOrbit r d) := by
  linear_combination lap_sub_closedOrbit hr d z

/-- **GDE-3** After `n` laps the offset from the closed orbit has turned by `r ^ n`. -/
theorem iterate_lap_sub_closedOrbit {r : ℂ} (hr : r ≠ 1) (d z : ℂ) (n : ℕ) :
    (lap r d)^[n] z - closedOrbit r d = r ^ n * (z - closedOrbit r d) := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [iterate_succ_apply', lap_sub_closedOrbit hr, ih]
    ring

/-- **GDE-4** A pure turn keeps the distance from the closed orbit: the offsets at the start of
each lap lie on a circle about it. -/
theorem norm_iterate_lap_sub_closedOrbit {r : ℂ} (hr : r ≠ 1) (hr1 : ‖r‖ = 1) (d z : ℂ) (n : ℕ) :
    ‖(lap r d)^[n] z - closedOrbit r d‖ = ‖z - closedOrbit r d‖ := by
  rw [iterate_lap_sub_closedOrbit hr, norm_mul, norm_pow, hr1, one_pow, one_mul]

/-- **GDE-5** So the orbit is bounded for every number of laps, also with an offset `w` reached
within the lap added: the bound does not depend on `n`. -/
theorem norm_iterate_lap_add_le {r : ℂ} (hr : r ≠ 1) (hr1 : ‖r‖ = 1) (d z w : ℂ) (n : ℕ) :
    ‖(lap r d)^[n] z + w‖ ≤ ‖closedOrbit r d‖ + ‖z - closedOrbit r d‖ + ‖w‖ := by
  have e : (lap r d)^[n] z + w = closedOrbit r d + ((lap r d)^[n] z - closedOrbit r d) + w := by
    ring
  calc ‖(lap r d)^[n] z + w‖
      = ‖closedOrbit r d + ((lap r d)^[n] z - closedOrbit r d) + w‖ := by rw [← e]
    _ ≤ ‖closedOrbit r d + ((lap r d)^[n] z - closedOrbit r d)‖ + ‖w‖ := norm_add_le _ _
    _ ≤ ‖closedOrbit r d‖ + ‖(lap r d)^[n] z - closedOrbit r d‖ + ‖w‖ := by
      gcongr
      exact norm_add_le _ _
    _ = ‖closedOrbit r d‖ + ‖z - closedOrbit r d‖ + ‖w‖ := by
      rw [norm_iterate_lap_sub_closedOrbit hr hr1]

/-! ### Without a turn: the same step every lap, and no orbit stays -/

/-- **GDE-6** Without a turn every lap adds the same step. -/
theorem iterate_lap_one (d z : ℂ) (n : ℕ) : (lap 1 d)^[n] z = z + n * d := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [iterate_succ_apply', ih]
    unfold lap
    push_cast
    ring

/-- **GDE-6** So no orbit closes unless the step vanishes. -/
theorem lap_one_ne {d : ℂ} (hd : d ≠ 0) (z : ℂ) : lap 1 d z ≠ z := by
  intro h
  unfold lap at h
  exact hd (by linear_combination h)

/-- **GDE-6** And every orbit leaves every bounded region. -/
theorem exists_lt_norm_iterate_lap_one {d : ℂ} (hd : d ≠ 0) (z : ℂ) (B : ℝ) :
    ∃ n : ℕ, B < ‖(lap 1 d)^[n] z‖ := by
  have hd' : 0 < ‖d‖ := norm_pos_iff.mpr hd
  obtain ⟨n, hn⟩ := exists_nat_gt ((B + ‖z‖) / ‖d‖)
  refine ⟨n, ?_⟩
  rw [iterate_lap_one]
  have h1 : B + ‖z‖ < n * ‖d‖ := (div_lt_iff₀ hd').mp hn
  have h2 : ‖(n : ℂ) * d‖ = n * ‖d‖ := by rw [norm_mul, Complex.norm_natCast]
  have h3 : ‖(n : ℂ) * d‖ ≤ ‖z + n * d‖ + ‖z‖ := by
    calc ‖(n : ℂ) * d‖ = ‖z + n * d - z‖ := by ring_nf
      _ ≤ ‖z + n * d‖ + ‖z‖ := norm_sub_le _ _
  linarith

/-! ### The reach is convex in `ρ` -/

/-- The offset at the start of each lap is an affine function of `ρ`. -/
theorem iterate_lap_affine (r G H z : ℂ) (n : ℕ) :
    ∃ a b : ℂ, ∀ ρ : ℝ, (lap r (ρ * G + H))^[n] z = a + ρ * b := by
  induction n with
  | zero => exact ⟨z, 0, fun ρ ↦ by simp⟩
  | succ n ih =>
    obtain ⟨a, b, h⟩ := ih
    refine ⟨r * (a + H), r * (b + G), fun ρ ↦ ?_⟩
    rw [iterate_succ_apply', h]
    unfold lap
    ring

/-- The length of an affine function of `ρ` is convex in `ρ`. -/
theorem convexOn_norm_affine (a b : ℂ) : ConvexOn ℝ Set.univ fun ρ : ℝ ↦ ‖a + ρ * b‖ := by
  refine ⟨convex_univ, fun x _ y _ p q hp hq hpq ↦ ?_⟩
  have hpq' : (p : ℂ) + q = 1 := by exact_mod_cast hpq
  have e : a + ((p * x + q * y : ℝ) : ℂ) * b = p * (a + x * b) + q * (a + y * b) := by
    push_cast
    linear_combination (-a) * hpq'
  simp only [smul_eq_mul]
  calc ‖a + ((p * x + q * y : ℝ) : ℂ) * b‖
      = ‖(p : ℂ) * (a + x * b) + q * (a + y * b)‖ := by rw [e]
    _ ≤ ‖(p : ℂ) * (a + x * b)‖ + ‖(q : ℂ) * (a + y * b)‖ := norm_add_le _ _
    _ = p * ‖a + x * b‖ + q * ‖a + y * b‖ := by
      rw [norm_mul, norm_mul, Complex.norm_of_nonneg hp, Complex.norm_of_nonneg hq]

/-- The largest of finitely many lengths of affine functions of `ρ` is convex in `ρ`. -/
theorem convexOn_sup'_norm_affine {ι : Type*} (s : Finset ι) (hs : s.Nonempty) (a b : ι → ℂ) :
    ConvexOn ℝ Set.univ fun ρ : ℝ ↦ s.sup' hs fun i ↦ ‖a i + ρ * b i‖ := by
  refine ⟨convex_univ, fun x _ y _ p q hp hq hpq ↦ ?_⟩
  simp only [smul_eq_mul]
  refine Finset.sup'_le hs _ fun i hi ↦ ?_
  have h := (convexOn_norm_affine (a i) (b i)).2 (Set.mem_univ x) (Set.mem_univ y) hp hq hpq
  simp only [smul_eq_mul] at h
  have hx : ‖a i + x * b i‖ ≤ s.sup' hs fun i ↦ ‖a i + x * b i‖ :=
    Finset.le_sup' (fun i ↦ ‖a i + x * b i‖) hi
  have hy : ‖a i + y * b i‖ ≤ s.sup' hs fun i ↦ ‖a i + y * b i‖ :=
    Finset.le_sup' (fun i ↦ ‖a i + y * b i‖) hi
  calc ‖a i + ((p * x + q * y : ℝ) : ℂ) * b i‖
      ≤ p * ‖a i + x * b i‖ + q * ‖a i + y * b i‖ := h
    _ ≤ p * (s.sup' hs fun i ↦ ‖a i + x * b i‖) + q * (s.sup' hs fun i ↦ ‖a i + y * b i‖) := by
      gcongr

/-- **GDE-7** The reach of a guiding centre launched at a fixed offset is convex in `ρ`. -/
theorem reach_convexOn {ι : Type*} (s : Finset ι) (hs : s.Nonempty) (laps : ℕ) (r G H z : ℂ)
    (g h : ι → ℂ) : ConvexOn ℝ Set.univ (reach s hs laps r G H z g h) := by
  choose A B hAB using iterate_lap_affine r G H z
  have e : reach s hs laps r G H z g h = fun ρ : ℝ ↦
      (Finset.range (laps + 1) ×ˢ s).sup'
        (Finset.Nonempty.product (Finset.nonempty_range_iff.mpr (Nat.succ_ne_zero laps)) hs)
        fun p ↦ ‖(A p.1 + h p.2) + ρ * (B p.1 + g p.2)‖ := by
    funext ρ
    unfold reach
    refine Finset.sup'_congr _ rfl fun p _ ↦ ?_
    rw [hAB]
    congr 1
    ring
  rw [e]
  exact convexOn_sup'_norm_affine _ _ _ _

/-- **GDE-7** So is the reach of the closed orbit. -/
theorem matchedReach_convexOn {ι : Type*} (s : Finset ι) (hs : s.Nonempty) (r G H : ℂ)
    (g h : ι → ℂ) : ConvexOn ℝ Set.univ (matchedReach s hs r G H g h) := by
  have e : matchedReach s hs r G H g h = fun ρ : ℝ ↦
      s.sup' hs fun i ↦ ‖(H * (r / (1 - r)) + h i) + ρ * (G * (r / (1 - r)) + g i)‖ := by
    funext ρ
    unfold matchedReach closedOrbit
    refine Finset.sup'_congr _ rfl fun i _ ↦ ?_
    congr 1
    ring
  rw [e]
  exact convexOn_sup'_norm_affine _ _ _ _

/-- **GDE-7** Hence the values of `ρ` whose reach stays inside a wall form one unbroken range: what
`particles.window` reports by its two ends. -/
theorem reach_lt_ordConnected {ι : Type*} (s : Finset ι) (hs : s.Nonempty) (laps : ℕ)
    (r G H z : ℂ) (g h : ι → ℂ) (wall : ℝ) :
    Set.OrdConnected {ρ : ℝ | reach s hs laps r G H z g h ρ < wall} := by
  have c := (reach_convexOn s hs laps r G H z g h).convex_lt wall
  simp only [Set.mem_univ, true_and] at c
  exact c.ordConnected

/-- **GDE-7** The same for the closed orbit. -/
theorem matchedReach_lt_ordConnected {ι : Type*} (s : Finset ι) (hs : s.Nonempty) (r G H : ℂ)
    (g h : ι → ℂ) (wall : ℝ) :
    Set.OrdConnected {ρ : ℝ | matchedReach s hs r G H g h ρ < wall} := by
  have c := (matchedReach_convexOn s hs r G H g h).convex_lt wall
  simp only [Set.mem_univ, true_and] at c
  exact c.ordConnected

/-! ### Controls: each hypothesis is needed

Nearby statements that are false, each refuted by a proof of its negation. -/

/-- **GDE-1** control: without a turn, the formula for the closed orbit does not give a fixed
point. -/
theorem lap_one_closedOrbit_ne {d : ℂ} (hd : d ≠ 0) :
    lap 1 d (closedOrbit 1 d) ≠ closedOrbit 1 d :=
  lap_one_ne hd _

/-- **GDE-2** control: without a turn and without a step every offset is a fixed point, so the
closed orbit is not the only one. -/
theorem lap_one_zero (z : ℂ) : lap 1 0 z = z := by
  unfold lap
  ring

/-- **GDE-4** control: a map that scales as it turns does not keep the distance, so `‖r‖ = 1` is
needed. -/
theorem norm_lap_two_sub_closedOrbit_ne :
    ‖lap 2 0 1 - closedOrbit 2 0‖ ≠ ‖(1 : ℂ) - closedOrbit 2 0‖ := by
  norm_num [lap, closedOrbit]

end Array16.Guide
