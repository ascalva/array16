<!-- 2026-11-18. The owner's question of 2026-11-16, answered from experiments/vivaldi_openems.py and the openEMS
source it runs on (openEMS 0.37). The counts below are for the design H 39.6, w_t 1.5, z_t 7.7 mm. -->

# Can the openEMS mesh be finer where it matters and coarser elsewhere?

Yes, within one rule, and the study's mesh already is.

**What openEMS offers.** It is a finite-difference time-domain solver on a rectilinear grid: three lists of mesh
lines, one per axis. The spacing along each list is free, so a mesh can be graded: fine at a feature, growing
away from it. A mesh line runs through the whole domain, so refining at one feature refines a whole slab across
the other two axes. The mesh is fixed for a run: there is no refinement that follows the fields in time. The
only sub-grid in the solver is for cylindrical coordinates (`MultiGrid`, "a cylindrical sub-grid radius").

**What the study's mesh does** (`mesh_lines`, `grade`).

| Region | Step |
| --- | --- |
| In the plane of a card, over the whole card | 0.357 mm (35 steps per half pitch) |
| Across the 0.8 mm laminate | 3 slices of 0.267 mm |
| Beside a card, to the next one | growing by 1.3 per step from 0.267 mm to 1.5 mm |
| Air beyond the cards, to the absorbing boundary | growing by 1.3 per step from 0.357 mm to 2.5 mm |

That gives 1.71 million cells for one row of six elements and 15.2 million for the 6×6 tile. 74 % of each slab
of the row model, and 63 % of the tile's, is at the finest step, because the card is meshed evenly from edge to
edge. The check on the finer mesh moved the score by 0.09 dB (`summary.md`), so the mesh is fine enough for the
search that used it.

**Where it can be coarser.**

- **Along the flare's height.** The fields vary fastest at the slot's throat and the feed, 7.7 mm up, and slowly
  at the mouth. Keeping 0.357 mm within 3 mm of the throat and growing gently to three times that at the mouth
  takes 56 steps over the card's height instead of 111: about half the cells.
- **Across the card it cannot be much coarser.** The flare's edge sweeps the whole width of a cell as it rises,
  and a rectilinear grid follows a curved edge in steps of its own size.
- **The time step is set by the smallest cell.** The three 0.267 mm slices of laminate set it at 0.89 of what
  0.357 mm cubes would allow. Finer cells anywhere cost twice: more cells, and more steps for the same simulated
  time.
- **Edges of metal.** A mesh line a third of a cell inside a metal edge and two thirds outside takes up the
  field's singularity there, and costs no cells.
- **Symmetry.** A magnetic or electric wall on a plane of symmetry halves the domain.

**For refinement in a small volume only**, the tool is a finite-element solver on an unstructured mesh, whose
cells can shrink at the throat alone. `REDESIGN.md`, section 12, already plans one, EMerge, for the single
periodic cell; openEMS stays the tool for the whole tile.
