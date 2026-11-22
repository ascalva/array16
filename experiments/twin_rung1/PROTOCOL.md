# Experiment 2, stage 1: the calibration rung on the twin bench. The protocol

Written on 2026-11-22, before any of its shots was tracked, and revised the same day on the reviews of the run's
physics and gates, still before any was tracked; section 6 lists everything that was run or known. This file,
`predictions.json`, `predictions.md` and `slot_straight.toml` are committed together with `experiments/twin_rung1.py`,
and the commit that adds the tracking's output comes after them: `array16.lab` refuses to measure until they are
committed, and stores with every measurement the commit it was registered at; `twin_rung1.py run` also refuses while
`predictions.md` is not what the code predicts. The numbers that decide the verdicts are constants at the top of
`experiments/twin_rung1.py` (`DISCOUNT`, `CONTROL_BAND`, `RATIO_TOL`, `WELL_BAND`, `WELL_STEPS`, `BETA_TOL`,
`SENSE_STEPS`, `LAMINATE_SDS`, `WINDOW`), in the same commit. Anything done differently from this protocol is listed
in the report.

## 1. The question

Stage 1 of experiment 2 (`docs/KNOTTED_GUIDE.md` section 10) is the first run that tracks electrons through the computed
field of a design rather than an ideal one: the design's cross-section solved by finite volumes
(`array16.guide.sections`), swept along its path in its frame and driven by its realised ring
(`array16.twin.TwinGuide.field()`), with `array16.particles.integrate` stepping the electrons, on the twin bench of
`benches/twin_r1.toml`. Before the knots' predictions of stage 2 rest on that instrument, it is calibrated on cases the
theory answers:

- **E2-21.** R1's rails' section swept round a ring with an ideal in-phase drive, the single-mode ring's, at 915 MHz:
  does its well lose electrons where the exact motion in the section says, at the design's 50 V (q 0.5) and at q 0.3,
  the control where the pseudopotential's depth should hold; and how much energy does it hold round its bend, at 50 V
  as drawn and in the bench's realisation (FR-4's permittivity and Earth's field drawn), and at q 0.3, the control
  where the static balance should apply? E2-21 is a statement about that model of R1, not about R1 as built: full wave
  finds R1's two rails resonating apart, and DR-0 does not pass the single-mode ring for R1 (section 7).
- **E2-22.** At the same section, volts and q, what does the knot cost? The flat facing-board ring of T915-S's 318 mm
  and the (1, 3) helix at r/R = 0.2 against T915-S itself.
- **The straight slot guide** (stage 1's instrument check): tracked in three dimensions through its swept section, does
  it give back its 2-D section's secular frequency and depth?
- **The case that must fail** (section 8.2): one board's rails at 915 MHz driven for a 1 eV well, past Mathieu's q of
  0.908. Nothing may be held, and the same launch at the design's 50 V, the control, holds everything.

## 2. The theory, and what it claims

Lengths in h, the distance from the null to the copper (to the slot's wall for the slot guide); D the well's depth at
the drive; volts are peak RF amplitudes across the resonator; the electron's |q/m| is 1.759 × 10¹¹ C/kg. Every number
below comes from the solved sections of the design files, the ring model (`array16.guide.ring`) and the functions of
`experiments/twin_rung1.py`, and `predictions.md` lists each prediction with what it rests on.

### 2.1 The well

**Claim (the pseudopotential).** An electron in the field E(x) cos ωt feels the potential |q/m|·E²/(4ω²), first order in
the null's Mathieu q = 2|q/m|·G·V/ω², G the null's gradient per volt (`array16.particles.rf`). The well's depth D is its
value at the escape saddle the section's minimax flood finds (`Solved.escape`): the lowest field an electron crosses to
reach copper, a board, or the open space beyond. On R1's solved section the saddle stands on the axis 1.173 mm up,
0.609 mm (1.08h) above the null, and at 50 V mean (48.95 V at the start station, node 0, where the two compiled feeds'
ripple leaves 0.979 of the mean) the pseudopotential there is 0.0909 eV deep: 0.0949 eV at the mean volts, 0.0904 eV
at the weakest node. The slot guide's saddle is on its axis too, 0.526 mm above the null, out of the slot into the gap
between the boards; at 40 V (39.74 V at the start) its depth there is 0.229 eV, 0.109 of gapless boards' at the same
volts and h, inside section 3.3's 0.10 to 0.14.

**Proposition (the rails' well is bounded under the null).** On rails the field's gradient grows toward the grounded
strip, so the null's local Mathieu q, 2|q/m|·V·|∇E|/ω², reaches 0.908 below the null: at the start station's 48.95 V
0.20h under it (0.19h at the mean 50 V), where the pseudopotential is only 0.75 D, and at q 0.3 0.36h under it, where
it is 4.3 D. An electron sent up toward the saddle
swings back down through the null; where its swing reaches that zone it is pumped out, whatever D the saddle sets. So
at q 0.5 the rails keep less than the pseudopotential's depth, and the escape is not D with a first-order correction
(a well stiffer near the null, Mathieu's exact exponent being 5.7 % above the first order's at q 0.5) but what the
exact motion gives; in the slot guide's well, which has no such zone, the exact motion escapes a little below D.

**Claim (the exact motion in the section).** The theory of the escape is the motion itself: on rails the saddle stands
on the axis of the section's mirror, where the field has no part across, so an electron sent from the null at the
saddle stays on the axis and moves as y″ = |q/m|·V(t)·∂φ/∂v, the field the swept grid gives there;
`twin_rung1.section_flights` integrates it by RK4 at 480 steps an RF period, at the start station's volts and the tuned
band, for each shot's own energies, its two zero crossings and its 40 RF periods, and loses an electron where the grid
does. In the slot guide the field on the axis points across it, and the same function integrates the motion in the
section's two dimensions. Its readings, registered as the predictions:

| Shot | Depth D at the start | First loss | Robust threshold | The zone where q passes 0.908 |
| --- | --- | --- | --- | --- |
| R1 at 50 V (q 0.51), as drawn | 0.0909 eV | 0.0118 eV, 0.130 D | 0.0337 eV, 0.37 D | 0.20h under the null, 0.75 D there |
| R1 at 50 V, the laminate's draw of the band over ±3 sds | 0.086 to 0.096 eV | 0.0115 to 0.0213 eV | | |
| R1 at q 0.3 (29.68 V mean), as drawn | 0.0320 eV | 0.0314 eV, 0.979 D | 0.0326 eV, 1.016 D | 0.36h under the null, 4.3 D there |
| the straight slot guide at 40 V (q 0.44) | 0.2294 eV | 0.938 D | the edges' threshold 0.927 D | |

At q 0.5 the losses come in islands of energy: from 0.13 D up, electrons of one launch phase are lost at some energies
and held at the next, and the robust threshold, the energy above which fewer than half of a window's electrons are held
(section 5), sits near 0.37 D; the electron's well at 50 V is 0.01 to 0.03 eV, against the pseudopotential's 0.091 eV.
At q 0.3 the zone lies beyond the barrier and the escape comes at D, within 2 %. In the slot guide, at q 0.44, the
escape comes at about 0.93 D. The same motion integrated by the tracker's own scheme, leapfrog with the field at each
step's start, at 120 and at 240 steps an RF period, loses the first electron at the same energy of the scan in both of
R1's as-drawn wells.

**What is measured.** Electrons start on the null with no motion along the path and a velocity straight at the saddle,
at the drive's two zero crossings, where the micromotion leaves a launch velocity equal to the slow one (ẋ = ẋ_s(1 −
(q/2) cos ωt) on the null); on R1 the first energy at which one is lost, against the exact motion's, and in the slot
guide the edges' threshold, against D.

### 2.2 The energy a bend holds: the static balance on the solved section

**Claim (JNT-2c, section 6.1).** Particles hold to the static balance: an electron of kinetic energy W on a bend of
curvature κ needs the force 2W·κ toward the bend's centre, and sits where the pseudopotential's gradient supplies it. On
the tilted well U + 2W κ·r an equilibrium exists until the tilt reaches the fold, where the local minimum meets a
saddle, or until the local Mathieu q at the equilibrium (from the field's gradient there) reaches 0.908 first. The
section's static strength in a direction is that largest force; the static limit of a path is the least, over its
samples, of the strength along κ there, times the drive's (|V|/V̄)² there, over 2|κ|.

**Method.** `experiments/twin_rung1.py` computes the strength on the section grid the swept field uses (the solved φ,
its gradient by second-order differences and bilinear between nodes, the copper and boards solid), on a patch of nodes
h/400 apart round the null smoothed over five, by walking the tilted well's valley out from the null: at each distance
along the tilt the least of U across it is the equilibrium there, and F = −dU/dt the force it holds. The local q comes
from φ's second derivatives at the grid's nodes, bilinear between them, so it is second order everywhere; at the cap
the limit is the last equilibrium under it, within the walk's step of the crossing. Its tests hold it to closed forms:
ELE-2's gapless boards hold section 6.1's 1.209·D/h across their plane, and toward a board the local q follows q₀
sec²(πy/2h) to 0.2 % and reaches 0.908 at 0.4675h (the closed form's 0.4677h), where the walk holds 5.135·D/h against
the closed form's 5.154; on ELE-1's rails (R1's 0.40 and 1.46 mm) the grid's strength across the board's plane meets
the valley walked on `RailPair`'s exact field, 3.546·D/h, within 0.5 % (0.06 % achieved).

**The sections' strengths** (the directions are where the electron is pushed):

| Section | Across the boards' plane | Toward the board below | Away from it |
| --- | --- | --- | --- |
| R1's rails on FR-4 at 50 V | 3.467 D/h, fold at (−0.283h, +0.126h), local q 0.35 | 9.21 D/h, q limit at 0.19h | 1.731 D/h, fold at 0.23h |
| the same at q 0.3 | 3.467 D/h, the same fold, local q 0.21 | 38.9 D/h, q limit at 0.35h | 1.731 D/h, the same fold |
| section B (T915-S, the flat ring, the helix) | 1.199 D/h, fold at 0.44h, local q 0.32 | 5.42 D/h, q limit at 0.495h | 5.42 D/h, the same |

The rails' well is shallow toward its saddle and stiff across the board's plane: the field across rises without a saddle
within reach (8.3 D at 2h), and R1 bends across. Section B's in-plane strength is ELE-2's 1.209 less 0.8 %, the split
and the strips' finite width.

**The paths' static limits**, with the drive as drawn (`TwinGuide.drive()` at realisation 0):

| Design | Path | Static limit | Where | Per eV of the weakest depth |
| --- | --- | --- | --- | --- |
| R1 | circle, 210.67 mm, R = 33.53 mm | 9.307 eV (9.77 at the mean volts); 3.280 eV at q 0.3 | the weakest node, 0.80 lap | 103 |
| T915-S | r/R = 0.45 trefoil framed k = 0, 317.76 mm | 11.64 eV | a tightest-bend lobe, 0.80 lap | 11.9 |
| flat ring 318 | circle, R = 50.57 mm | 32.43 eV | the weakest node | 33.1 |
| helix (1, 3) | r/R = 0.2 torus curve framed k = 0 | 18.45 eV | 0.78 lap | 18.9 |

T915-S's limit lies between the two bounds the guide document took: 9.80 eV were every bend across the boards' plane
(its tightest, 15.29 mm), 11.91 eV were the bends toward the boards free (its in-plane radius, 18.57 mm). The helix's
bounds are 16.08 and 18.69 eV.

**Claim (the tracked fraction).** Tracking holds less than the static balance: four rods on the first run's trefoil held
0.96 of their static limit at q 0.3 and lost 1.03 of it, and held 0.61 at q 0.6 (section 6.1), and T915-S is designed on
0.61 to 0.95 of its static limit until E2-9 tracks facing boards. Each facing-board path here is predicted to hold
**0.61 to 0.95 of its static limit**.

**Hypothesis (R1 at q 0.5: four rods' discount on rails).** R1 holds 0.61 to 0.95 of its static limit at 50 V. Four rods
have one q everywhere; the rails at q 0.5 have the zone 0.19h under the null where the local q passes 0.908, and a fold
0.28h out, so any secular swing that reaches the zone ejects the electron (section 2.1), and the dynamic aperture the
zone sets may hold well under 0.61 of the static limit. The band is the four-rod discount's, registered as a hypothesis
with that competing cause named; 0.61 to 0.95 is not carried to T915-R until this shot and its control are in.

**Claim (the control: the static balance at q 0.3).** At q 0.3 the zone lies 0.35h under the null, beyond the fold's
reach, as four rods' cap lay beyond theirs at q 0.3, where they held 0.96 and lost 1.03 of the static limit: R1 at q
0.3 holds **0.85 to 1.05 of its static limit**, 3.280 eV, the four rods' bracket widened by a tenth below for the rails'
own launch and by the exact exponent's stiffening (1.8 % at q 0.3) above.

**Hypothesis (the price of the knot).** Paths that share a section, its volts and its q hold the same fraction of their
static limits, so the tracked ratio of two paths' limits is the ratio of their static limits: **the flat ring holds
2.786 times what T915-S holds, and the helix 1.585 times**, each within 15 %.

### 2.3 Earth's field on a flat ring

**Proposition (the senses' mean; algebra on the static balance).** A uniform field B_z along a flat ring's axis adds
s·e·v·B_z to the bend's force for an electron going round in sense s, so the limit W_s solves 2W/R − s·e·B_z·√(2W/m) =
F: a quadratic in √W whose roots for s = ±1 multiply to F·R/2. The two senses' limits move apart, √W₊ − √W₋ =
B_z·R·√(e/2m), and their geometric mean stays at the limit without the field, whatever B_z. The tests hold `ring_limit`
to both.

**Claim (Earth's field in any orientation).** The part of the field in the ring's plane pushes the electron up and down
once a lap, e·v·B_h, and up is where the rails are soft (1.73 D/h). Over 400 orientations of 50 µT the static balance on
R1 holds 7.37 to 11.19 eV in either sense, the two senses' geometric mean 8.09 to 9.31 eV (0.869 to 1.000 of the limit
without the field), and |√W₋ − √W₊| at most 0.586 √eV. The bench draws the orientation; the prediction does not know it.

### 2.4 The secular frequency, and the case that must fail

**Proposition (Mathieu).** Near the null the motion is Mathieu's y″ + 2q cos(2τ)·y = 0 along each principal axis, its
slow part at the exponent β, the root of the continued fractions (McLachlan; Abramowitz and Stegun 20.3): read once an
RF period, x(k + 1) + x(k − 1) = 2 cos(πβ)·x(k). On the straight slot guide at the start station, 39.74 V and the tuned
914.983 MHz put the null at q 0.4421, so β = 0.326034, 4.3 % above the first order's 0.3126. The slot section's grid
gives the solved gradient at the null to 2 × 10⁻⁶.

**Proposition (no well).** Past q = 0.908 Mathieu's equation at a = 0 has no stable solution, so no well holds: on R1's
section a 1 eV well at 915 MHz takes 162.4 V and puts the null at q 1.64, where a population grows by 2.35 Np an RF
period, and q 1.89 takes 187.0 V (a 1.33 eV well), 2.66 Np. Section 8.2's q 1.89 for the 1 eV well is the rails 0.5 mm
up; R1's solved null is 0.564 mm up, so both are run.

**Claim (the control).** The must-fail shots' own launch, the same starts, energies, phases and 200 RF periods on the
bench's guide, at the design's 50 V (q 0.51), holds every electron: its highest energy, 4 eV, is 0.43 of the static
limit, under the 0.61 a hold keeps and under the 0.75 a cubic fold loses an electron started at rest at its untilted
null, and the starts within 0.05h of the null swing far short of the zone 0.19h under it. "None held" at q 1.64 is
read only beside it, since a sign error, a field scaled wrong or a grid that loses everything would also hold none.

### 2.5 The guide document's numbers against the theory

E2-21 in `docs/KNOTTED_GUIDE.md` section 10 reads: "R1's rail ring at 915 MHz (null 0.575 mm, 50 V, q 0.5) holds its
0.093 eV well and electrons to 0.58 to 0.86 of the bend rule's 2.6 eV round its 32 mm radius". E2-22 reads: "the flat
facing-board ring of 318 mm holds 2.7 to 3.3 times what T915-S holds, and the (1, 3) helix at r/R = 0.2 1.4 to 1.6
times". Each number checked against the design files as solved:

| Number | Guide document | Now | Why |
| --- | --- | --- | --- |
| R1's null | 0.575 mm | 0.564 mm | the solved section: 0.15 mm isolation and 35 µm copper; the gapless √(c₁c₂) is 0.576 mm |
| its q at 50 V | 0.5 | 0.506 | the solved gradient at the tuned 914.83 MHz |
| its well | 0.093 eV | the pseudopotential's 0.0949 eV at the mean volts, 0.0909 at the start; the electron's, first lost at 0.0118 eV | the solved saddle and the two feeds' ripple at Q 36; the zone under the null where q passes 0.908 |
| its ring | about 200 mm, 32 mm in radius | 210.67 mm, 33.53 mm | one guided wave on the solved index 1.555; the document estimated n_eff from √((1 + ε_r)/2) |
| the bend rule D·R/2h | 2.6 eV | 2.82 eV | the same rule on the solved well, null and radius |
| what R1 holds | 0.58 to 0.86 of the rule: 1.5 to 2.2 eV | Hypothesis: 0.61 to 0.95 of the static limit, 9.31 eV: 5.68 to 8.84 eV | see below |
| flat ring over T915-S | 2.7 to 3.3 | 2.786, band 2.37 to 3.20 | see below |
| helix over T915-S | 1.4 to 1.6 | 1.585, band 1.35 to 1.82 | see below |

**What R1 holds.** The bend rule D·R/2h takes the well to be harmonic with its depth reached at h, so that the bend's
force at h/2 is D/h. The rails' depth is reached at the saddle above the null, where the well is soft (1.73 D/h), while
across the board's plane, where R1 bends, the static balance holds 3.467 D/h: 3.47 times the rule, 9.31 eV at the
weakest node. The document's 0.58 to 0.86 were four rods' rings tracked as fractions of the rule (section 6.1); with the
static balance in place of the rule its discount is the one section 6.1 designs T915-S on, 0.61 to 0.95, and on rails
at q 0.5 that discount is the hypothesis of section 2.2. The band is therefore 5.68 to 8.84 eV as drawn, and the
document's 1.5 to 2.2 eV lies below it; at q 0.3 the control's band is 2.79 to 3.44 eV.

**The price of the knot.** The document's 2.7 to 3.3 is run.log's 2.72 to 3.31: the flat ring's radius over T915-S's
in-plane radius (the bends toward the boards free) and over its tightest bend (every bend across the plane), two bounds
on the boards' stiffness toward the knot's bends. The solved section fixes that stiffness, 5.42 D/h toward the boards at
the local q limit against 1.199 D/h across, and the static balance along the k = 0 frame puts T915-S at 11.64 eV,
between its bounds, and the ratio at 2.786, inside the document's band near its foot. The helix's 1.4 to 1.6 rounds
run.log's 1.35 to 1.64, whose lower end divides the helix's tightest bend by T915-S's in-plane radius: the two bounds at
once, which one section cannot have. Taken together the bounds give 1.569 to 1.640, and the solved stiffness 1.585, at
the top of the document's band. With the stiffness known, what the static balance leaves open is how much of its static
limit each path holds when tracked; the bands give that fraction 15 % of room between paths, and the theory's ratios sit
at their centres.

## 3. The apparatus

All numerical, on the bench `twin-r1` (`benches/twin_r1.toml`):

| Role | Driver | Options |
| --- | --- | --- |
| guide | `array16.twin:TwinGuide` | `designs/r1_rail_ring.toml`, realisation 1 |
| detector | `array16.twin:TwinDetector` | pickups at 0 and 0.25 lap, a 1 pF collector |
| vacuum | `array16.twin:TwinVacuum` | 1 × 10⁻⁶ mbar: a mean free path of 400 m |
| array | `array16.sim:SimBackend` | the HackRF Pro alone; not used here |

The experiment needs the guide's `tune`, `drive` and `field`, the detector's `watch`, the vacuum's `scatter` and
`free_path`, and each one's `identity`. Shots on R1 "on the bench" use the bench's guide; the other shots open their
design at realisation 0 beside it, sharing its cache of solved sections, and their records carry each guide's identity
(file, sha256, realisation, section key, solver).

| Design | File | Section | Path and frame | Drive |
| --- | --- | --- | --- | --- |
| R1 | `designs/r1_rail_ring.toml` | rails 0.40 / 1.46 mm, 0.15 mm isolation, 1.6 mm FR-4, 35 µm copper | circle in the board's plane, u inward | 50 V, two feeds |
| T915-S | `designs/t915s.toml` | section B: h 0.915 mm, 8h wide, split 0.3h, 25 µm polyimide, 18 µm copper | r/R 0.45 trefoil, k = 0 | 50 V, two feeds |
| flat ring 318 | `designs/flat_ring_318.toml` | section B | circle bent in the boards' plane | 50 V, two feeds |
| helix (1, 3) | `designs/helix_13.toml` | section B | r/R 0.2 (1, 3) torus curve, k = 0 | 50 V, two feeds |
| straight slot guide | `experiments/twin_rung1/slot_straight.toml` | S2-C3's: 1.0 mm slot, 0.5 mm strips set back 0.2 mm, 0.813 mm RO4003C, planes 1.0 mm off | straight, one guided wave (213.49 mm) repeating | 40 V, two feeds |

Each ring is m = 1 at 915 MHz, one guided wave on its solved index, and is driven at its tuned centre through two
current feeds compiled for the travelling wave (`array16.guide.ring.compile_feeds`). The bench's realisation of R1 draws
FR-4's permittivity (sd 0.1, the index's sd 0.0133, 0.86 % in the band) and Earth's 50 µT in a uniform orientation; R1
draws no spacers, so its impedance and drive pattern are as drawn. A shot at q 0.3 drives the same ring at the volts
that put the null's q there on the solved section at the design's band, 29.68 V mean.

- **Field.** The section's φ per volt swept along the path (`array16.particles.SectionField`), E = −V(s, t)·∇φ in the
  frame and B = T × E/v at the wave's phase velocity, V(s, t) the ring's node volts continued between the nodes
  (`array16.twin.RingWave`), Earth's field added on the bench. An electron off the section's grid (4h either way of the
  null) or in its copper or boards is lost.
- **Tracker.** `array16.particles.integrate`, the relativistic Boris step, at 120 steps an RF period (section 10), held
  to Mathieu's exact secular frequency on section B to 0.1 % (DR-0).
- **Observers.** For every scan, the bench's gas (`TwinVacuum.scatter`, seeded) and then the detector's watch
  (`array16.twin.Watch` with the detector's weight and collector), which counts each electron held, lost to the guide
  or stopped by the gas, and its laps. The secular frequency's four electrons run with the watch alone: at 10⁻⁶ mbar
  the gas would stop one in 10⁵ in its 60 RF periods. The pickups' currents bear on no criterion and are not recorded.

## 4. The shots

Fourteen (`twin_rung1.py::shots`), seeded 20261122 plus each one's place in the list; the energies of each hold and each
well run a factor 1.025 apart, each range covering its band with room either side.

| Shot | Design, realisation | What | Energies | Electrons | Run |
| --- | --- | --- | --- | --- | --- |
| E2-21 well, R1, realisation 0 | R1, 0 | escape from the null toward the saddle | 0.0045 to 0.1393 eV across, 140 | 280 | 40 RF periods |
| E2-21 well, R1, on the bench | R1, bench | the same | the same | 280 | 40 RF periods |
| E2-21 well, R1, realisation 0 at q 0.3 | R1, 0, 29.68 V | the same | 0.006 to 0.0489 eV across, 86 | 172 | 40 RF periods |
| E2-21 hold, R1, realisation 0 | R1, 0 | held round the bend, both senses, eight phases | 1.0 to 12.72 eV along, 104 | 1664 | 20 laps at 1 eV: 795 486 steps |
| E2-21 hold, R1, on the bench | R1, bench | the same, Earth's field drawn | the same | 1664 | the same |
| E2-21 hold, R1, realisation 0 at q 0.3 | R1, 0, 29.68 V | the same, eight phases | 0.8 to 4.97 eV along, 75 | 1200 | 20 laps at 0.8 eV: 889 380 steps |
| E2-22 hold, T915-S, realisation 0 | T915-S, 0 | held round the knot, two phases | 3.5 to 15.78 eV, 62 | 248 | 20 laps at 3.5 eV: 641 462 steps |
| E2-22 hold, flat ring 318, realisation 0 | flat ring, 0 | held round the ring | 9.5 to 43.91 eV, 63 | 252 | 20 laps at 9.5 eV: 389 356 steps |
| E2-22 hold, helix (1, 3), realisation 0 | helix, 0 | held round the helix | 5.5 to 25.42 eV, 63 | 252 | 20 laps at 5.5 eV: 511 711 steps |
| slot beta, straight slot guide, realisation 0 | slot guide, 0 | secular frequency | at rest, 0.004h off the null on both diagonals and both axes | 4 | 60 RF periods |
| slot well, straight slot guide, realisation 0 | slot guide, 0 | escape toward the saddle | 0.16 to 0.327 eV across, 30 | 60 | 40 RF periods |
| must fail, R1, on the bench at a 1 eV well | R1, bench, 162.4 V | held near the null | 0.5, 1, 2, 4 eV along | 32 | 200 RF periods |
| must fail, R1, on the bench at q 1.89 | R1, bench, 187.0 V | the same | the same | 32 | 200 RF periods |
| must fail control, R1, on the bench | R1, bench, 50 V | the same launch | the same | 32 | 200 RF periods |

- **A hold.** Every energy in both senses, each RF phase one batch: the drive's phase at the start station drawn
  uniformly within each of the period's equal parts, one part a phase (eight on R1, two on E2-22's paths), from the
  shot's seed. Each electron starts at the path's first sample on the static equilibrium its energy takes there
  (`Valley.at` on the device's own section and drive, the bend's force alone), with the micromotion of that place and
  phase added to its position and velocity, −(|q/m|/ω²)·Re(V e^(iωt₀))·∇φ and (|q/m|/ω)·Im(V e^(iωt))·∇φ, the
  velocity's taken half a step before t₀ as the Boris step's leapfrog takes it, so that it starts on its slow motion (in
  a uniform field such a launch keeps its place to 20 nm over ten RF periods, where one at rest drifts a third of a
  millimetre). The run lasts until the slowest has gone 20 laps at its launch speed (2 % over, `run_steps`); an
  electron counts as held once it has made 20 laps, whatever came after, lost when the guide lost it sooner, stopped
  when the gas did, and short when the run ended before its 20 laps; neither of the last two is judged.
- **A well.** Every energy at the two zero crossings, one batch each, on the null with the velocity at the saddle and
  nothing along the path; held or lost after 40 RF periods, seven slow periods.
- **The secular frequency.** Four electrons at rest, read once an RF period, β from the least-squares ratio of section
  2.4 over both coordinates, with no gas.
- **A failure and its control.** One launch for the three shots, drawn from the master seed: electrons within 0.05h of
  the null, the four energies along the path in both senses at four phases, no micromotion matched; held or lost after
  200 RF periods.

## 5. What counts as agreement

Fixed before the run. Three readings of a scan's row (one sense, or a well's one row):

- the **first loss**, the lowest energy at which an electron was lost;
- the **robust threshold**, the lowest energy at which fewer than half of the electrons judged in the window of five
  energies centred on it were held, its geometric mean with the energy below (below the scan if that is the lowest,
  above it if there is none); a sharp edge gives that edge back, and a lone resonant loss below it does not move it;
- the **edges' threshold**, the geometric mean of the highest energy held below the first loss and the first loss.

A hold's limit is the two senses' robust thresholds' geometric mean, and the first loss is reported beside it.

| | Criterion | Passes when |
| --- | --- | --- |
| E2-21 well | R1 at 50 V as drawn | the first loss lies within two scan steps, and half a step for the scan's own, of the exact motion's 0.01179 eV: 0.01108 to 0.01254 eV |
| E2-21 well | R1 at 50 V on the bench | the first loss lies within the exact motion's first losses over the laminate's ±3 sds of the band, 13 draws, 0.01150 to 0.02132 eV, widened the same: 0.01081 to 0.02268 eV |
| E2-21 well control | R1 at q 0.3 as drawn | the first loss lies within two and a half steps of the exact motion's 0.03138 eV (0.979 D): 0.02950 to 0.03338 eV |
| E2-21 hold | what R1 holds at 50 V (Hypothesis) | as drawn the limit lies within 0.61 to 0.95 of the static limit, 5.68 to 8.84 eV; on the bench within 0.61 of Earth's least mean and 0.95 of its most, widened by three sds of the laminate's draw, (1.026)² either way: 4.69 to 9.30 eV |
| E2-21 hold control | what R1 holds at q 0.3 | the limit lies within 0.85 to 1.05 of the static limit, 3.280 eV: 2.79 to 3.44 eV |
| E2-21 senses | the senses as theory says | as drawn (both drives) the two senses' robust thresholds lie within two scan steps of each other; on the bench \|√W₋ − √W₊\| ≤ 0.670 √eV, Earth's most (0.586) and two scan steps in √W at 11.19 eV |
| E2-22 static limit | each path against its own | the limit lies within 0.61 to 0.95 of the path's static limit: T915-S 7.10 to 11.06 eV, the flat ring 19.78 to 30.81 eV, the helix 11.26 to 17.53 eV |
| E2-22 price of the knot | the ratios | the flat ring's limit over T915-S's lies within 15 % of 2.786 (2.368 to 3.204), the helix's within 15 % of 1.585 (1.347 to 1.822) |
| slot secular frequency | 3-D against 2-D | every electron is held and its β lies within 0.2 % of 0.326034 (0.325382 to 0.326686) |
| slot depth | 3-D against 2-D | the edges' threshold lies within 0.80 to 1.25 of the depth at the start, 0.2294 eV (0.1835 to 0.2867 eV) |
| must fail control | the same launch at 50 V | no electron is lost to the guide, and one is held |
| must fail | no well past q 0.908 | no electron of either shot is held, beside a control that passed |

**Why these sizes.** The discount is section 6.1's, the range four rods tracked and T915-S is designed on; it is wide
because tracking at q near 0.5 has held as little as 0.61 of the static limit and as much as 0.96. The control's 0.85
to 1.05 is section 2.2's. The ratios' 15 % is the room the static balance leaves to that fraction between paths that
share a section, volts and q: the room experiment 1 gave a first-order theory against tracking
(`experiments/guide/PROTOCOL.md` section 5), half again the document's own band's half width (10 %). It separates the
bounds the document's band spanned: the flat ring's price with every bend across the plane at full weight, 3.31, lies
above it, and with the bends toward the boards at no weight, 2.72, inside it. The wells' two steps: the tracked motion
is the exact motion on the same grid, so the room is for the tracker's integrator (its leapfrog at 120 and 240 steps
loses the first electron at the same energy as the RK4) and for the 3-D sweep's own small terms, the wave's magnetic
force (v_e/v_phase of the electric, 10⁻³ here) and its drift along the path; a half step more puts the band's edges
between the scan's energies. The slot depth's band keeps the pseudopotential's room at q near 0.5 and the scan's step,
with the exact motion at 0.93 D inside it. The secular frequency's 0.2 % is twice DR-0's 0.1 %, for a section DR-0 did
not test, and ten times what the integrator misses by at 120 steps (1.5 to 1.9 × 10⁻⁴ on section B). The senses' two
steps let each sense's robust threshold sit up to a step from its truth, half a step of the scan's resolution and half
a step of its window's sampling, so their difference two. The laminate's three sds move the tuned band by 2.6 % and
the pseudopotential's depth, which goes as V²/f², by 5.2 %; the bench's well takes the exact motion over that draw
instead, since at q 0.5 its first loss moves with q by more than the depth does.

**Earth on the bench.** The hold's launch takes the bend's force alone, so on the bench Earth's force, up to 89 eV/m at
9 eV against the bend's 540 eV/m there, starts each electron on a secular swing, and the bench's band, a factor 1.98
wide because the prediction does not know the orientation, can hardly fail: the as-drawn shots are E2-21's sharp tests
of the hold, and the bench's shot tests the senses' split. The bench's well band is the laminate's alone; Earth's force
on an electron at the first loss, 3 eV/m at 0.012 eV, is 1 % of the rails' strength toward the saddle.

**What a failure would mean.**

- *E2-21 well control* off: the tracker does not feel the pseudopotential's depth at q 0.3, or the swept field is not
  the section's; the 50 V well, R1's holds and E2-22's then stand only as ratios.
- *E2-21 well* off at 50 V with the control in: the tracker does not follow the exact motion where it is resonant, the
  islands at q 0.5 (the integrator at 120 steps, the 3-D sweep's frame or its magnetic term), and the q 0.5 holds read
  through that.
- *E2-21 hold control* off: the static balance does not hold on rails even at q 0.3, the tracker or the static
  computation is at fault, and the 50 V hold says nothing alone.
- *E2-21 hold* low at 50 V, with the control in: the local q zone 0.19h under the null sets the rails' dynamic aperture
  at q 0.5, not the static balance; R1's and T915-R's bends are then designed at q near 0.3 or on the tracked fraction.
  Low and near 1.5 to 2.2 eV with the control low too: the rails hold to the bend rule, as the document read them, not
  to their static balance across the board's plane. High: the static balance understates the hold.
- *E2-21 senses* off as drawn: something other than Earth tells the senses apart (the wave's own fields, the launch, the
  integrator): the instrument is at fault, and nothing else here stands until it is not. On the bench: Earth's term is
  not what the bench draws.
- *E2-22 static limit* off on one path: JNT-2c's discount does not hold for facing boards on that path; E2-9 is then
  registered on the tracked fraction.
- *E2-22 price of the knot* above the band: the knot costs more than its static price, the turning of its curvature
  through the frame, not only its tightest point, setting what it holds; stage 2 then designs T915-S on the tracked
  price and the flat ring's 33 eV per eV is its ceiling. Below: the ring loses more of its limit than the knot does.
- *Slot secular frequency* off: the swept slot section is not its solve, or the integrator is off on it, and stage 3's
  predictions on S2-C3 cannot rest on this instrument. *Slot depth* off: the same, for the escape.
- *Must fail control* off: the tracker, the field or the launch loses electrons at a drive that holds them, and the
  must-fail shots' "none held" says nothing. *Must fail* held: the tracker or the field is wrong, and nothing else here
  stands.

## 6. What was known when this was written

- **The twin's own tests**, run before this was written (`tests/test_twin.py`, `tests/test_twin_dr0.py`): electrons at
  rest within 0.05h of R1's null held for 5 RF periods; 5 eV electrons sent up leaving the aperture and 100 eV ones sent
  down striking the copper; electrons aimed at R1's copper turned back near the rail edges at 5 eV and by the rails' 50
  V at 20 eV; one 1 eV electron on a straight guide passing a pickup; DR-0's secular frequency in section B at q 0.3,
  met to 0.1 % at 120 steps an RF period and missed at 30. None follows an electron along a bend or for a lap.
- **The four-rod tracking of section 6.1**, the discount's source, and the full-wave run of R1's line and ring
  (`experiments/twin_openems_line`).
- **This experiment's tests** (`tests/test_twin_rung1.py`) and the checks made while writing it, all on coarse-grid
  copies of R1, T915-S and the slot guide, whose sections this protocol does not predict: holds of two energies, 1.0 and
  1.025 eV on R1 for 0.004 of a lap (on the bench's realisation, 159 steps) and 3.5 and 3.59 eV on T915-S for 0.003 of a
  lap (97 steps), every electron held; R1's well at 0.02 and 0.0205 eV for 2 RF periods, at 50 V and at q 0.3, held; the
  slot guide's four electrons for 4 RF periods, reading β = 0.3266 on the coarse section; a must-fail at a 1 eV well and
  its control for 2 RF periods, every electron still on the grid; a hold cut to 20 steps, every electron short of its
  laps; timing runs of 300 steps, 124 and 208 electrons started on the null of R1 and of T915-S at 0.7 eV (300 to 330
  µs a step); single electrons launched in a uniform oscillating field, which is no design, for ten RF periods; and the
  section model on a harmonic well, held at q 0.8 and lost at q 1.0. No registered shot was run, no electron went round
  more than 0.004 of a lap, and none was started near a registered threshold.
- **The theory** of section 2 was computed while writing this, from the solved sections and the ring model only: the
  static balance, and the exact motion in the section (`section_flights`) on the registered wells' own energies and
  phases, by RK4 at 480 steps and by the tracker's leapfrog at 120 and 240 steps, with the slot guide's in two
  dimensions on its section and, before that, on a wider scan.
- **The reviews of this protocol**, of its physics and of its gates, before any shot was tracked: their 1-D
  integrations on the axis of `RailPair`'s gapless field (the first loss at 0.13 D at q 0.51, escape near D at q 0.3,
  islands of energy over 200 RF periods), their Floquet growth past the edge (2.35 Np an RF period at q 1.64, 2.66 at
  1.89), and the independent check of the rails' strengths by Newton continuation (3.546 D/h across, 1.77 D/h up on the
  gapless pair). This revision registers the exact motion in place of the pseudopotential's depth for R1's wells, the
  controls, the robust threshold, the eight phases on R1's holds and the hold's laps counted as made.
- **The multipactor instrument on R1** (`array16.guide.multipactor`, MPC-3): R1's rails at 50 V read −5.22 dB at level
  1 against copper's chart, the parallel-plate bound on their arcs, and the calibrated population grows by more than
  0.4 Np an RF period on their 0.70 GHz mm arcs (`tests/test_guide_multipactor.py::test_mpc_3_level_1_margins`,
  `::test_mpc_3_the_instrument_reads_the_designs_alike`; `docs/LAWS.md`, MPC-3). The twin tracks without secondary
  emission, so E2-21 tests holding at a drive that MPC-3 says multipacts on bare copper; 13.7 V or less keeps level 1's
  6 dB there, and the 2-D tracking in the rail field (DR-8, E2-13) settles it. The owner chose to register as it
  stands (2026-11-22): this rung calibrates the tracking chain against the theory at 50 V, and stage 2's tracking with
  secondary emission settles what a built board holds.
- Nothing else. No electron had been tracked round R1, the flat ring, the helix or T915-S for a lap, nor sent at a
  well's saddle near its first loss on the default grids.

## 7. The model

- **The field** (`array16.twin`): the section quasi-static and swept along the path, its null where the 2-D solve puts
  it, |V| varying round the ring as the realised drive gives it (2.4 % round R1, 0.4 % round the facing rings) and the
  field along the path from the bend and the ripple coming with the 3-D solve, the boundary-element solve
  `array16/guide/bem.py` will hold (section 4.5), with the bend's shift of the null (7.7 µm on T915-S's tightest bend,
  section 4.4) and the joints. The three dimensions of the straight slot guide are the swept section's field in space,
  with the travelling wave's E and B, tracked in three dimensions; the boundary-element solve of the straight slot guide
  is the field this check is repeated in once it is built.
- **The magnetic field** is T × E/v at the wave's phase velocity from the loaded section's E, where a quasi-TEM line on
  a laminate carries the air line's pattern, (n/c)·T × E_air: the magnetic force is v_e/v_phase of the electric, about
  1 % at 9 eV, and with the forward wave's sign the Doppler shift and the magnetic term cancel in the pseudopotential to
  that order, so the senses agree as drawn; with the wrong sign they would differ by about 8 %, which the senses'
  criterion would see. The ring's backward part, |b/a| under 10⁻³, is given the forward wave's slowness. A design that
  needs the air pattern sweeps the air solve's φ for B.
- **R1's ring** is the single-mode ring of `array16.guide.ring`, both rails driven in phase. DR-0 does not pass that
  model for R1: on the backed board with via walls full wave finds the closed ring's two rails resonating apart, at
  −1.71 % and +5.15 % of 915 MHz (899 and 962 MHz), and only the two-rail model, whose line matrices live in
  `experiments/twin_openems_line.py`, meets them (within 0.15 %). The rails' centre lines differ by ±0.93 mm in radius,
  ±2.8 % in length, against a linewidth of 1/Q = 2.7 %, and their coupling through the grounded strip is weak (C₁₂/C₁₁
  −0.02 to −0.055), so this holds on the unbacked board as drawn too. E2-21 is therefore stated for the ideal in-phase
  drive. Before any R1 prediction is used for hardware, the two-rail line matrices with the bridges' admittances join
  `array16.guide.ring` and `TwinGuide`, the bridges (or their absence) the design file, the unbacked R1 is run in full
  wave, and the residual odd-mode field at the null is registered as its own prediction.
- **The electrons** go one at a time: space charge, image charges and beam loading come with WarpX or the 2.5-D model
  of section 6.5 (E2-26).
- **The surfaces** are copper with no charge of their own: contact potentials and patches come as the twin's DC draws
  (section 6.1), the charge on exposed dielectric with the 3-D solve, and secondary emission with Furman and Pivi's
  model in the 2-D multipactor tracker (E2-13, MPC-H1), which is where MPC-3's reading of R1 at 50 V is settled.
- **The gas** scatters with one cross-section, 1.02 × 10⁻¹⁹ m², at every energy: at 10⁻⁶ mbar it stops 1 to 4 % of a
  hold's electrons, which are counted apart. LXCat's cross-sections by energy will take its place for a named gas.

## 8. Running it

    uv run python experiments/twin_rung1.py predict                  # the predictions; about 50 s once solved
    uv run python experiments/twin_rung1.py run --bench twin_r1      # the tracking; resumable, one shot with --only
    uv run python experiments/twin_rung1.py report                   # the verdicts

The bench's devices run in one process, shot after shot: about two hours on this machine, R1's three holds about 100
minutes of it (eight batches of 150 to 208 electrons at about 300 µs a step), then T915-S's, the helix's and the flat
ring's holds, about a quarter of an hour; the wells, the secular frequency and the failures take seconds to a minute
each. `run` and `report` recompute the predictions to check them against the files, about 50 s.

## Sources

1. N. W. McLachlan, *Theory and Application of Mathieu Functions* (Oxford, Clarendon Press, 1947).
2. M. Abramowitz and I. A. Stegun, *Handbook of Mathematical Functions*, chapter 20 (National Bureau of Standards,
   1964).
3. `docs/KNOTTED_GUIDE.md` sections 3.6, 6.1, 8.2, 8.3 and 10, and `experiments/knotted_guide/run.log`, the numbers of
   section 2.5.
