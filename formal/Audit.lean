import Array16

/-!
# The audit of what the proofs rest on

`lake env lean Audit.lean` lists what stands behind every declaration in the `Array16` modules and
fails if any rests on more than Lean's three standard ones (`propext`, `Classical.choice`,
`Quot.sound`), or on an unfinished proof: one a tactic or an option produced is caught here, through
the environment, where a search of the text would miss it. It reports the count it audited.
(After `formal/Audit.lean` in the owner's mind-palace.)
-/

open Lean Elab Command

elab "#audit_array16" : command => do
  let env ← getEnv
  let standard := #[``propext, ``Classical.choice, ``Quot.sound]
  let mut checked : Nat := 0
  let mut unfinished : Array Name := #[]
  let mut escapes : Array (Name × Name) := #[]
  for (name, _) in env.constants.map₁.toList do
    let some idx := env.getModuleIdxFor? name | continue
    let some mod := env.header.moduleNames[idx.toNat]? | continue
    unless (`Array16).isPrefixOf mod do continue
    checked := checked + 1
    for ax in ← liftCoreM (collectAxioms name) do
      if ax == ``sorryAx then unfinished := unfinished.push name
      else unless standard.contains ax do escapes := escapes.push (name, ax)
  unless escapes.isEmpty do
    throwError m!"rests on more than the standard three: {escapes.toList}"
  unless unfinished.isEmpty do
    throwError m!"{unfinished.size} declarations rest on an unfinished proof: {unfinished.toList}"
  if checked == 0 then
    throwError "no declaration audited: is the library built?"
  logInfo m!"{checked} declarations audited; all rest on the standard three only"

#audit_array16
