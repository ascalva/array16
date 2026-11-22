module

public import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic

/-!
# The ring resonator: RNG-1 to RNG-4, the definitions

The definitions the laws of `Array16/Ring.lean` are stated over. They mirror
`array16/guide/ring.py`, and a reader checks them against that module and against `docs/LAWS.md`.

A closed guide is a ring of transmission line. Along it the volts are two waves running opposite
ways, `a e^(-iβs) + b e^(iβs)`, written here in the phase `θ = βs`. A lap of `m` guided waves is
`θ = 2πm`. Feeds are currents injected at places on the ring: a feed of current `w` at phase `θ`
launches the wave of order `+m` with `w e^(iθ)` and the wave of order `-m` with `w e^(-iθ)`, each
times the same factor, which the laws here leave out.
-/

@[expose] public section

namespace Array16.Ring

open Complex

/-- The volts at phase `θ = βs` of two waves running opposite ways, of amplitudes `a` (order `+m`)
and `b` (order `-m`): what `ring.ring_purity` takes apart and `ring.ripple` measures. -/
noncomputable def wave (a b : ℂ) (θ : ℝ) : ℂ :=
  a * exp (-(θ * I)) + b * exp (θ * I)

/-- The ripple of a size that runs between `lo` and `hi` (`ring.ripple`). -/
noncomputable def ripple (hi lo : ℝ) : ℝ :=
  (hi - lo) / (hi + lo)

/-- The energy a travelling wave of peak volts `V` stores on a length `L` of line of capacitance
`C` per metre: `C V² / 4` per metre in the electric field and as much in the magnetic. -/
noncomputable def stored (C V L : ℝ) : ℝ :=
  C * V ^ 2 * L / 2

/-- The power a resonance at `f` of quality `Q` dissipates while it stores `W`: `2πf W / Q`. -/
noncomputable def dissipated (f Q W : ℝ) : ℝ :=
  2 * Real.pi * f * W / Q

/-- The shunt resistance of a ring of `m` guided waves of a line of impedance `Z0` and quality `Q`
(`Ring.shunt` in `ring.py`, `Solved.shunt` in `sections.py`). -/
noncomputable def shunt (Z0 Q m : ℝ) : ℝ :=
  Z0 * Q / (2 * Real.pi * m)

/-- What `K` feeds at the places `s_j = j L / K` of a ring, of currents `w j`, launch into the wave
of order `n`: the sum of `w j e^(2πi n j / K)`. Order `m` is the wanted wave and `-m` the backward
one. -/
noncomputable def launch (K : ℕ) (w : ℕ → ℂ) (n : ℤ) : ℂ :=
  ∑ j ∈ Finset.range K, w j * exp (2 * Real.pi * I * (n * j / K))

/-- The currents of `K` equal feeds phased for the travelling wave of order `m`: feed `j` lags by
the wave's own phase at its place, `e^(-2πi m j / K)`. -/
noncomputable def travelling (K : ℕ) (m : ℤ) (j : ℕ) : ℂ :=
  exp (-(2 * Real.pi * I * (m * j / K)))

/-- What two feeds a quarter guided wave apart, of currents 1 and `w`, launch into the wanted
wave: the second feed stands at phase `π/2` (`ring.ring_score`'s pair, nodes 0 and `N / 4m`). -/
noncomputable def pairForward (w : ℂ) : ℂ :=
  1 + w * exp (Real.pi / 2 * I)

/-- What the same two feeds launch into the backward wave. -/
noncomputable def pairBackward (w : ℂ) : ℂ :=
  1 + w * exp (-(Real.pi / 2 * I))

/-- A wave of amplitude `a` after a length `s` of line that attenuates by `α` and turns the phase
by `β` per metre (`ring.ring_admittance`'s `g`, with `α = β / (2 Q)`). -/
noncomputable def decayed (a : ℂ) (α β s : ℝ) : ℂ :=
  a * exp (-((α + β * I) * s))

end Array16.Ring
