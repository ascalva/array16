module

public import Mathlib.Analysis.SpecialFunctions.Pow.Real

/-!
# The line of a cross-section: LIN-1, the definitions

The definitions the laws of `Array16/Line.lean` are stated over. They mirror the line solve of
`array16/guide/sections.py` (`_line`, `energy2d`), and a reader checks them against that module and
against `docs/LAWS.md`.

A cross-section's capacitance per metre, at one volt of drive, is twice the field energy of its
potential, and the potential is the one of least energy among those that hold the copper's volts
(Thomson's principle; the finite-volume solve is that minimisation on a grid). With a laminate of
relative permittivity `ε` in part of the section, the energy of a potential `u` is `A u + ε B u`:
`A` the part in air and `B` the part in the laminate per unit permittivity, both independent of
`ε`. The laminate's share is its part of the whole.

The smallest such solve has one free node: two capacitors in series across one volt, air (`ca`)
and laminate (`ε cl`), the free node at `x` volts across the air.
-/

@[expose] public section

namespace Array16.Line

/-- The field energy of a potential `u`, in units where the capacitance is the energy itself:
`A u` in air and `ε B u` in the laminate (`sections.energy2d` over the cells' permittivity). -/
def energy {E : Type*} (A B : E → ℝ) (ε : ℝ) (u : E) : ℝ :=
  A u + ε * B u

/-- The laminate's share of the field energy of the potential `u` (`Solved.share`). -/
noncomputable def share {E : Type*} (A B : E → ℝ) (ε : ℝ) (u : E) : ℝ :=
  ε * B u / energy A B ε u

/-- The energy in an air capacitor `c` with `x` volts across it. -/
def stackAir (c x : ℝ) : ℝ :=
  c * x ^ 2

/-- The energy per unit permittivity in a laminate capacitor `c` in series with it, the two across
one volt. -/
def stackLaminate (c x : ℝ) : ℝ :=
  c * (1 - x) ^ 2

/-- The volts the solve puts across the air of the series pair: the divider's. -/
noncomputable def divider (ca cl ε : ℝ) : ℝ :=
  ε * cl / (ca + ε * cl)

/-- The capacitance of the series pair. -/
noncomputable def stackCapacitance (ca cl ε : ℝ) : ℝ :=
  ca * (ε * cl) / (ca + ε * cl)

end Array16.Line
