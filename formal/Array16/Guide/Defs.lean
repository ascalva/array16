module

public import Mathlib.Analysis.Complex.Basic

/-!
# The lap map of a guide field: GDE-1 to GDE-7, the definitions

The definitions the laws of `Array16/Guide.lean` are stated over. They mirror `GuideOptics` in
`array16/particles/optics.py`, and a reader checks them against that class and against
`docs/LAWS.md`.

A magnetic guide field follows a closed curve. To first order in gyroradius over bend radius, a
guiding centre's offset from the curve is a complex number `z` on axes carried along the curve
without twisting. One lap adds a step `d` (the curvature drift plus the wander of the field line)
and then turns the plane by `r = exp (i α)`, where `α` is the curve's holonomy. The step is
`ρ * G + H`, with `ρ` the gyroradius, signed by the direction of travel.
-/

@[expose] public section

namespace Array16.Guide

/-- One lap: add the step `d`, then turn the transverse plane by `r` (`GuideOptics.lap`). -/
def lap (r d z : ℂ) : ℂ :=
  r * (z + d)

/-- The offset of the orbit that closes on itself after one lap (`GuideOptics.closed_orbit`). -/
noncomputable def closedOrbit (r d : ℂ) : ℂ :=
  d * (r / (1 - r))

/-- The largest offset over `laps + 1` laps and the stations `s` within each lap, for a guiding
centre launched at `z` (`GuideOptics.excursion`). `g i` and `h i` are the drift per unit `ρ` and the
wander accumulated at station `i`, and `ρ * G + H` is the step of a whole lap. -/
noncomputable def reach {ι : Type*} (s : Finset ι) (hs : s.Nonempty) (laps : ℕ) (r G H z : ℂ)
    (g h : ι → ℂ) (ρ : ℝ) : ℝ :=
  (Finset.range (laps + 1) ×ˢ s).sup'
    (Finset.Nonempty.product (Finset.nonempty_range_iff.mpr (Nat.succ_ne_zero laps)) hs)
    fun p ↦ ‖(lap r (ρ * G + H))^[p.1] z + ρ * g p.2 + h p.2‖

/-- The largest offset of the closed orbit itself over the stations `s`
(`GuideOptics.matched_excursion`). -/
noncomputable def matchedReach {ι : Type*} (s : Finset ι) (hs : s.Nonempty) (r G H : ℂ)
    (g h : ι → ℂ) (ρ : ℝ) : ℝ :=
  s.sup' hs fun i ↦ ‖closedOrbit r (ρ * G + H) + ρ * g i + h i‖

end Array16.Guide
