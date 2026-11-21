module

public import Mathlib.Analysis.Complex.Basic

/-!
# Electrodes printed on a board: ELE-1, the definitions

The definition the laws of `Array16/Electrode.lean` are stated over. It mirrors `RailPair.field` in
`array16/particles/electrodes.py`, and a reader checks it against that class and against
`docs/LAWS.md`.

Across long electrodes the field is two-dimensional, and `E_x - i E_y` is an analytic function of
`z = x + i y`, with `y` the height above the board. Two rails at `V` in a grounded plane, from `c1`
to `c2` either side of the centre line, have `E_x - i E_y = (i V / π) * railField c1 c2 z`.
-/

@[expose] public section

namespace Array16.Electrode

/-- The field of a rail pair with edges at `±c1` and `±c2`, as `E_x - i E_y` without its factor
`i V / π` (`RailPair.field`). -/
noncomputable def railField (c1 c2 : ℝ) (z : ℂ) : ℂ :=
  1 / (z - c2) - 1 / (z - c1) + 1 / (z + c1) - 1 / (z + c2)

end Array16.Electrode
