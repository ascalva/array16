<!-- Draft pending an adversarial check, 2026-11-15. Written by the lead engineer from two investigator dossiers
(structured fields; shielding). Every derived number was recomputed in a check script (lead_calc.py, in the session
scratch folder, not kept in the repo), and the external facts were checked against the sources in section 7.
experiments/local_fields.md and experiments/fields_study.md did not exist when this was written; where this note
leans on the local-fields simulation it quotes that script's run log of 2026-11-15 21:34, which may change. -->

**Status: draft, pending an adversarial check.**

# Field uses: knots, shields, and what a milliwatt field can do

The owner asked: "i can see a knot almost acting like a stable field over a volume, what could such field geometries
be used for? another interesting question: the earth magnetic field shields us from radiation, can a local field
operate at some level close to that? what could such a field produce or shield from?"

This note answers both for rev B (`REDESIGN.md`: 0.4-6 GHz, a 16-element Vivaldi panel at 2.3-6 GHz, 16 external
whips, four coherent channels, transmit at milliwatts). Two parallel notes cover neighbouring ground and are not
repeated here: `local_fields.md` (focal spots, quiet zones, vortex lines and the optimiser's knot attempts near the
array) and `fields_study.md` (magnets, ferrites/YIG, the Sun, vortex beams, superdirectivity). "(est.)" marks an
estimate; all other numbers are derived here or cited.

## 1. Short answers

- **A knot in the array's field is stable in shape, not in substance.** The knotted curves are lines of exact
  darkness (zero field) with the phase winding around them. Small disturbances bend them but cannot untie them;
  untying takes two lines touching and reconnecting [1, 2]. Nothing is held there, though: at full drive a 27-litre
  cube 0.3 m from the whips contains about 3 pJ of field energy, roughly 5×10^7 times less than the heat in 1 mm³ of
  air. All of it is replaced every nanosecond, and the knot is gone within nanoseconds of switching off. A knot is
  also a curve, not a volume: it threads a dim region (17-33 dB below the focus in the parallel simulation) and
  occupies no volume itself.
- **Rev B cannot form one.** In the parallel local-fields simulation the optimiser formed a trefoil and two linked
  rings only with 64 dipoles on a 250 mm sphere around the volume. Rev B's panel and its 16 whips managed straight
  vortex lines, but not rings, links or the trefoil.
- **What such geometries are used for elsewhere.** In optics, dark lines and doughnut beams trap and spin micron-sized
  particles, beat the diffraction limit (STED microscopy), trap cold atoms, sense rotation and carry information
  [11-23]. In plasma, twisted and knotted magnetic fields hold themselves together [7-9]. At this project's
  milliwatts and centimetre wavelengths, forces stay below 100 pN and torques below 10^-11 N·m, so nothing moves.
  What survives is measurement:
  - nulls, for interference rejection and precise direction finding;
  - rotational Doppler, which detects spin and its sense;
  - mapping the structured field itself.
- **Can a local field shield like Earth's? No.** What matters is field times depth. Turning a charged particle back
  needs a B·L of about its momentum per charge: 5.7 T·m for a 1 GeV cosmic-ray proton. Earth supplies about 96 T·m
  from only 30 µT because its field extends thousands of kilometres. Over 1 m that takes about 6 T. The array's field
  is 17 nT at 0.3 m and reverses 10^9 times a second, which averages its push to about 10^-7 eV. Even superconducting
  spacecraft shields studied at 8 T·m (35-96 t) cut the cosmic-ray dose by only 8-24 % [33].
- **No field shields light, X-rays, gamma rays, neutrons or radio waves.** They carry no charge, and in air
  electromagnetic waves pass through each other. Matter does that job: the atmosphere (10 t/m²) for us, lead for
  gamma rays, water for neutrons, metal boxes for radio [42].
- **On the ground, Earth's field does less than people think.** Cosmic rays give about 0.4 mSv a year of a roughly
  2.4 mSv natural background [41]. Their sea-level dose varies by only about 10 % between latitudes [40]. The
  atmosphere is the main shield; the magnetosphere matters most for satellites, aircrews, astronauts, power grids and
  GPS [39].
- **The array's real shield is a receive null.** Four coherent channels give three nulls. `REDESIGN.md` §9 expects
  20-30 dB; adaptive weights should reach about 40 dB on a narrowband interferer (est.). Nothing is transmitted.
- **What its field can produce: signals, sensing, and microwatts.** At full drive the field 20 cm away is at most 3 %
  of the FCC public exposure limit. All the power the four channels can produce (≤25 mW) would warm 1 g of water by
  at most 0.006 K/s. At 1 m a receiver picks up microwatts.
- **What it can watch: Earth's shield at work.** It can follow the Sun's 10.7 cm flux (2.8 GHz is in the panel's
  band), solar radio bursts that jam GPS (+19 dB on a whip for a 10^6 sfu burst, est.), and the ionosphere's electron
  content through two-frequency GNSS.

## 2. Knotted and structured fields

### 2.1 In what sense a knot is stable

Near the array, any field component is a sum of 16 (or 32) fixed patterns weighted by complex excitations. Where its
complex amplitude is exactly zero (two real conditions in three dimensions) the zeros form lines. The phase turns by
2π times an integer around each line, and the lines can close, link and knot. Optical beams with isolated knotted
dark lines were built this way [1].

| Sense of "stable" | Does the array's knot have it? | Numbers |
| --- | --- | --- |
| Topological | Yes | Perturbations move the lines; the knot type changes only when two lines touch and reconnect [1, 2] |
| Against calibration error | Yes, if the knot can be formed | A perturbing field δ moves a zero line by about \|δ\|/\|∇ψ\| ≈ (\|δ\|/\|ψ_local\|)·λ/2π (est.). With rev B's 3°/0.26 dB per element, the error field is -36 dB (16 elements) or -42 dB (64) below the focus. The knot's dark region is 17-33 dB below the focus, so its lines move 0.01-0.08 λ, less than the spacing of its strands (est.). In the local-fields run the 64-element trefoil kept its topology with these errors |
| In time | Only while driven | A CW pattern is fixed in space but its energy flows through at c: 0.3 m in 1.0 ns. Switch off and it is gone in nanoseconds |
| In energy | No | At 0.3 m from the four whips at full drive the field holds 1.1×10^-10 J/m³: 3.1 pJ in 0.027 m³, against 1.5×10^-4 J of thermal energy in 1 mm³ of air |
| Self-sustaining | No; that needs conducting matter | See below |

Knots that hold themselves together:
- **Free space.** Maxwell's equations have solutions whose field lines are all closed and linked (Rañada's Hopf
  fields). They propagate, expanding "like an opening umbrella" while keeping their Hopf structure [5]. A family of
  null fields keeps torus knots and links in its field lines as it flows along the Poynting vector [6]. Either way the
  knot flies away at the speed of light; it does not stay over a volume.
- **Plasma.** In a slightly resistive plasma the linkage of field lines (magnetic helicity) decays much more slowly
  than the field energy, so the plasma relaxes to a twisted minimum-energy state. This explains the reversed-field
  pinch [7]. Measurements on the Madison Symmetric Torus found helicity falling less than energy during relaxation
  [8]. Knotted flux tubes in resistive MHD simulations relax into self-organised nested-torus structures that last
  much longer than the Alfvén time [9].
- **Fluids.** Knotted vortices in water, about 10 cm across, stretch, reconnect and end as unlinked rings after about
  0.5 s [4].
- **Light in turbulence.** An optical knot's invariant survives weak but not strong turbulence [2]. Persistence of
  the most robust knots fell below 26 % in moderate turbulence and to 4-8 % in the strongest case tested [3].
- **Magnets.** Skyrmions, small stable whirls of magnetisation, are proposed as bits in racetrack memories [10].

So the array's knot is stable the way a standing wave's shape is stable, not the way an object is.

### 2.2 What structured fields are used for

| Use | Structure | Where it works | Why it works there | Source |
| --- | --- | --- | --- | --- |
| Optical spanner | doughnut beam with a phase singularity | absorbing micron particles | light's orbital angular momentum spins trapped particles; the sense reverses with the singularity's sign | [11] |
| STED microscopy (Nobel 2014) | doughnut-shaped depletion beam | fluorescent molecules | only the dark centre fluoresces, which beats the diffraction limit | [12] |
| Dark-core atom trap | blue-detuned doughnut beam | laser-cooled Rb, 10^8 atoms, 150 ms | the dipole force pushes atoms into darkness | [13] |
| Rotational Doppler | orbital-angular-momentum (OAM) modes | mm waves, light, a 5.8 GHz radar | the shift is lΩ/2π, independent of carrier frequency | [14-18] |
| Vortex metrology | phase singularities as markers | optics, nanometre displacement | a zero is a sharper marker than a peak | [19] |
| OAM multiplexing | several vortex beams | 32 Gbit/s over 2.5 m (mm-wave) | it is a special case of MIMO: no extra capacity | [20, 21] |
| Topological coding | framed optical knots; skyrmion textures in microwave pulses | optics; a broadband horn | the invariants survive smooth deformation | [22, 23] |
| Acoustic levitation | twin, vortex and cage traps | 64-element ultrasonic arrays, 40 kHz | per watt, sound carries c/c_sound ≈ 9×10^5 times more momentum than light | [24, 25] |
| Plasma confinement | twisted and linked field lines | tokamaks, RFPs, spheromaks | field lines guide charged particles; helicity is conserved | [7-9] |
| Ion trap (Paul, Nobel 1989) | oscillating quadrupole saddle | RF on mm electrodes | time-averaged confinement, up to several eV deep | [26] |
| Magnetic bottle | Earth's dipole | Van Allen belts, 0.2-10 Earth radii | particles bounce between the poles | [27] |

### 2.3 Radio versions at this project's scale

**Forces and torques: no.**

| Quantity | Bound for rev B | Comparison |
| --- | --- | --- |
| Power the four channels deliver, an upper bound on radiated power (REDESIGN §7: +6/+5/+3 dBm each at ≤1.3/2.44/5.8 GHz) | 16 / 12.6 / 8 mW; at most 25 mW at Nuand's +8 dBm | |
| Total force on a perfect absorber, P/c | 53 / 42 / 27 pN; at most 84 pN | the weight of a 25 µm water droplet |
| Radiation pressure, 0.3 m from the four whips at 433 MHz | 1.1×10^-10 Pa | at the FCC limit of 10 W/m² [29]: 3.3×10^-8 Pa |
| Gradient force on a 1 mm water drop, full-contrast standing pattern at 3.6 V/m | 4×10^-18 N | its weight is 4.1×10^-5 N |
| OAM torque lP/ω per unit l | 5.9×10^-12 N·m (433 MHz), 8.3×10^-13 (2.44 GHz), 2.2×10^-13 (5.8 GHz) | measured on a 30 cm copper ring with a GHz turnstile: about 10^-8 N·m [28], 1,700 times more |
| First turn of a frictionless 10 g, 30 cm ring absorbing all of it | 6 h (433 MHz) to 32 h (5.8 GHz) | a 1 g, 5 cm disc: 27 min to 2.4 h, but it is far smaller than the wavelength and cannot absorb the beam |

Optics manages because P/c does not depend on wavelength, but the size of object a beam can grip does. A
wavelength-sized object is (12 cm / 1 µm)³ ≈ 1.7×10^15 times heavier at 2.44 GHz than in an infrared trap, where a
1 µm bead weighs 5 fN and piconewtons dominate. Sound carries the same intensity with 8.7×10^5 times more momentum
(the ratio of the speeds of light and sound). That is why ultrasonic arrays levitate beads [24, 25] and radio arrays
cannot.

**Rotational Doppler: yes. This is the best structured-field demonstration this array can do.**
- Physics: a component with azimuthal phase e^(ilφ), scattered by a body spinning at Ω, comes back shifted by lΩ/2π.
  The shift does not depend on the carrier, and polarisation adds ±1 [14, 15]. It was first measured with mm waves
  [15], and a spinning object was detected with light [16].
- Radio-frequency demonstrations:
  - 24.83 Hz measured at 50π rad/s with l = 1, error at most 0.67 % [17];
  - a commercial 5.8 GHz radar on a spinning OAM metasurface: ±24 Hz for l = ±1 and ±48 Hz for l = ±2 at about 24
    revolutions per second [18].
- Numbers: 1200 rpm gives 20 Hz per unit of l, 60 Hz at l = 3; 3000 rpm gives 50 Hz per l. For comparison, linear
  Doppler is 16.3 Hz per m/s at 2.44 GHz and 38.7 Hz per m/s at 5.8 GHz; a blade tip at 15 cm and 1200 rpm, seen
  edge-on, gives ±307 Hz and ±729 Hz.
- **No vortex transmitter is needed.** Illuminate along the spin axis, receive on a ring, and split the ring's signals
  into azimuthal modes (a spatial DFT). Mode m carries its line at -mΩ/2π, so the sign gives the sense of rotation.

The check script models this (est.): scalar point scatterers, forward scatter, the panel facing a three-blade fan.
- Fan: 15 cm blades and a foil tab on one blade, at 0.4 m, 1200 rpm.
- Illuminator: a horn on the fan's axis 1 m from the panel, shining through the fan.
- Errors: rev B's 3°/0.5 dB per channel.
- Assumed scattering lengths: 5 cm per blade and 1 cm for the tab.

| Receiver | Mode bin | Line, counter-clockwise | Line, clockwise | Level re direct signal |
| --- | --- | --- | --- | --- |
| 4 quadrant channels | +1 | -20 Hz | +20 Hz | -34 dB |
| | -1 | +20 Hz | -20 Hz | -35 dB |
| | ±2 (one bin with four channels) | ±40 Hz, both | ±40 Hz, both | -43 dB |
| 12 perimeter elements (Hadamard snapshots) | +1 / +2 / +3 | -20 / -40 / -60 Hz | +20 / +40 / +60 Hz | -33 / -40 / -39 dB |
| one element alone | | ±20, ±40, ±60 Hz, equal on both sides | the same | -33 to -37 dB |

What the model shows:
- The static direct signal leaks into the m ≠ 0 bins at -21 to -35 dB, but at 0 Hz. In the 12-element square ring,
  the m = ±4 bins carry a -17 dB static term from the geometry alone.
- Aliases separate by frequency. With four channels, the m = -3 component lands in the +1 bin but at +60 Hz.
- Misalignment does not matter much. With the fan axis 10 mm and the horn 20 mm off axis, the line frequencies and
  levels are unchanged (within 0.1 dB); the static leakage rises to about -21 dB.
- At 2.44 GHz with the fan at 0.25 m, the ±1 lines are again -34 dB.
- **A single element cannot tell the sense of rotation; the ring can.**

Signal budget (est.):
- With a 15.249-level horn at 1 m, the direct signal at a quadrant is about -38 dBm, so the lines are at -71 to
  -81 dBm. The receiver noise is -172 dBm/Hz.
- So the limits are the phase-noise skirt of the direct signal 20-200 Hz from the carrier, and fan vibration, not
  thermal noise.
- Phase noise from the common 10 MHz reference cancels at nanosecond delays. If the AD9361's uncorrelated
  synthesizer noise there is -80 dBc/Hz (assumed; to be measured), the lines clear the skirt by 37-47 dB in a 1 Hz
  bin, and the mode split removes another 20-30 dB of it.

The whip plate is a weaker host:
- Its vertical whips have their null along the ring axis (`fields_study.py` table 1b).
- A fan above the plate, lit by the centre monopole (CAL-ANT), still forms an axisymmetric set-up, so the mode split
  holds (est.).
- At the ≤-41 dBm EIRP of MANIFOLD mode, 0.01 m² at 1 m returns about -112 dBm. After the whips' weak response near
  the axis, that leaves 30-60 dB of SNR in 1 Hz (est.).

**Nulls as precision tools: yes.**
- The radio counterpart of STED's dark centre is the monopulse or steered null: a zero can be located much more
  finely than a peak.
- The thermal-noise accuracy is σθ ≈ θ3dB/(k_m√(2·SNR)), with k_m ≈ 1.6 [30]. For the panel's 26° beam at 5.8 GHz
  that is 1.1° at 20 dB SNR and 0.11° at 40 dB.
- Calibration then sets the floor: 0.5-1° with coupling and the measured manifold (`REDESIGN.md` §9).

**Near-field focusing: real but weak.** `local_fields.md` gives spot sizes and V/m per mW. With at most 25 mW
radiated, a focused spot delivers microwatts per square centimetre: enough to signal to, or trickle-power, a tag
within tens of centimetres, nothing more (est.).

**OAM multiplexing: no gain.** OAM channels are a subset of what MIMO already provides [20]. `fields_study.py`
(table 1d) shows mode l fading with distance as D^(-2(|l|+1)), so every mode but l = 0 falls faster than the usual
D^-2.

**Knots as information: weak.** Optical framed knots carry information in their invariants [22]. Reading a radio knot
needs a 3-D field scan, and rev B cannot form one.

**Polarisation textures: blocked on rev B.** The panel is single-polarisation (`fields_study.py` section 3).

### 2.4 Ranking

| Rank | Use at this scale | Real on rev B? | Key number | Value |
| --- | --- | --- | --- | --- |
| 1 | Receive nulls: interference rejection and null-based DF | yes, daily | 3 nulls, 20-40 dB (est.); DF 0.1-1° | serves all four jobs |
| 2 | Rotational Doppler: spin rate and sense | yes | 20 Hz per l at 1200 rpm; lines 33-43 dB below the direct signal (model) | best demonstration of OAM physics |
| 3 | Mapping the array's own structured field with a probe | yes | 4,096 points in a 10 cm cube at λ/8, about 35 min (est.) | turns the local-fields simulation into a measurement |
| 4 | Near-field focusing of power or data | partly | µW/cm² at 20 cm (est.) | small |
| 5 | Knots as markers or keys | no: needs about 64 elements (simulation) | rev B makes straight vortex lines only | none yet |
| 6 | OAM multiplexing | real, but it is MIMO | no capacity gain [20] | none |
| 7 | Torque, trapping, levitation | no | ≤84 pN; ≤6×10^-12 N·m | none (use ultrasound, section 6) |
| 8 | Shielding particles | no | ≤10^-7 eV | none |

## 3. Shielding

### 3.1 How Earth's field shields, and why scale decides

A magnetic field bends a charged particle onto a circle of radius ρ = p/(qB). To turn a particle back, the field
region must be at least one such radius deep: B·L ≥ p/q. When B·L is smaller, it deflects the particle by about
B·L/(p/q) radians. A check of the formula: a 7 TeV proton in the LHC's 8.33 T dipoles has ρ = 2,803 m; the LHC
design gives 2,804 m [31].

B·L needed (derived):

| Particle | p·c | B·ρ needed | Gyroradius in 30 µT |
| --- | --- | --- | --- |
| Solar-wind proton, 1 keV | 1.37 MeV | 0.0046 T·m | 150 m |
| Beta electron, 1 MeV | 1.42 MeV | 0.0047 T·m | 160 m |
| Alpha particle, 5 MeV (charge 2) | 193 MeV | 0.32 T·m | 11 km |
| Solar energetic proton, 10 MeV | 137 MeV | 0.46 T·m | 15 km |
| Solar energetic proton, 100 MeV | 445 MeV | 1.5 T·m | 49 km |
| Cosmic-ray proton, 1 GeV | 1.70 GeV | 5.7 T·m | 190 km |
| Cosmic-ray proton, 10 GeV | 10.9 GeV | 36 T·m | 1,200 km |
| Sea-level muon, 4 GeV mean energy [38] | 4.0 GeV | 13 T·m | 440 km |

B·L available:

| Source | B·L | Notes |
| --- | --- | --- |
| Earth, equatorial, from the surface outward (B0·RE/2) | 96 T·m | B0 ≈ 30 µT; the surface field is 25-65 µT [32] |
| One LHC dipole, 8.33 T over 14.3 m | 119 T·m | inside a 56 mm bore at 1.9 K [31] |
| Spacecraft toroid study (SR2S) | 8, 11.5, 23 T·m | 96 t at 8 T·m cuts the cosmic-ray blood-forming-organ dose from 21.9 to 16.7 cSv/yr (24 %); 23 T·m (about 140 t) leaves about 10 cSv/yr [33] |
| Spacecraft solenoid study (NIAC) | 8 T·m | 35 t; the same dose from 25.6 to 23.6 cSv/yr (8 %) [33] |
| 3 T MRI magnet across 1 m | 3 T·m | est. |
| NdFeB magnet, 0.5 T over 5 cm | 0.025 T·m | est. |
| The array's RF field at 1 m, if it were static | 5×10^-9 T·m | 5.1 nT peak (four whips, 433 MHz, full drive) |

How Earth does it:
- **Cosmic rays.** For a dipole, the vertical cutoff is (B0·RE·c/4)·cos⁴(latitude). That gives 14.3 GV at the
  magnetic equator (B·ρ = 48 T·m, half the 96 T·m integral), 8.1 GV at 30°, 3.6 GV at 45° and 0.9 GV at 60° (derived,
  B0 = 30 µT). The ECSS standard fits observed cutoffs with 16 GV for the constant [34]. Near the poles the field lines
  lead in and the cutoff falls towards zero.
- **The solar wind.** Slow wind moves at about 400 km/s [36]; take about 5 protons/cm³ (est.). Its ram pressure of
  1.3 nPa balances 58 nT at the boundary, which puts the magnetopause at about 10 Earth radii (derived; typically
  6-15 [35]). This is a collective plasma effect. It needs a field region larger than the ion gyroradius (about 76 km
  in 60 nT) and the ion inertial length (about 100 km) (derived). Lunar magnetic anomalies a few hundred kilometres
  across, with nanotesla-order fields, do form mini-magnetospheres: an electric field set up by the electrons
  deflects the ions [37].
- **The atmosphere.** It weighs 101,325 Pa / g = 10.3 t/m², the same as 10 m of water (derived), and absorbs the
  primary cosmic rays. What reaches the ground is mostly muons, about 1 per cm² per minute [38], plus neutrons. NASA
  credits the magnetosphere "together with Earth's atmosphere" [39].
- **What the field does for people on the ground.** The sea-level cosmic dose (270 µSv/yr, plus 48 µSv/yr from
  neutrons) varies by about 10 % with latitude [40]. Cosmic rays are about 0.4 mSv/yr of an average natural
  background of about 2.4 mSv/yr; radon in homes gives more than half of the 2 mSv/yr terrestrial part [41]. At
  ground level the field's share is of the order of a few hundredths of a millisievert a year (est.).
- **Earth's field also traps.** The Van Allen belts are a magnetic bottle in which particles bounce between the poles
  [27].

### 3.2 What a local field can and cannot shield

| Radiation | Charged? | Local field needed | Practical shield [42] | Verdict |
| --- | --- | --- | --- | --- |
| Solar wind | yes | 0.005 T·m, but over ≥100 km of plasma | none needed on the ground | space only |
| Beta electrons, 1 MeV | yes | 0.005 T·m: a 0.5 T magnet bends them on a 9 mm radius | thin metal or plastic | a magnet works; matter is simpler |
| Alpha particles, 5 MeV | yes | 0.32 T·m | paper, skin | not needed |
| Solar protons, 10-100 MeV | yes | 0.5-1.5 T·m | spacecraft walls | superconducting magnets in space |
| Cosmic-ray protons, muons (GeV) | yes | 6-36 T·m | the atmosphere, rock | no local field |
| X-rays, gamma rays | no | none | lead, concrete | no field |
| Neutrons | no | none | water | no field |
| Radio and light | no | none: fields superpose | a closed metal box | no field; cancellation only at points |

**Static fields.** A permanent magnet offers about 0.025 T·m, enough for electrons up to a few MeV. Even an MRI-class
3 T·m would turn back only protons below about 360 MeV, and nothing neutral (est.). An electric field stops a
particle of kinetic energy T only with a potential of T/q: 1 GV for a 1 GeV proton (derived).

**RF fields.** The push reverses every half cycle; what remains on average is the ponderomotive energy
U_p = e²E0²/(4mω²) [43]. At 0.3 m from the four whips at full drive (433 MHz) that is 1.5×10^-7 eV for electrons and
8×10^-11 eV for protons, against 0.025 eV of thermal energy (derived). A 1 eV barrier at 433 MHz would need 13 kV/m for
electrons and 0.56 MV/m for protons. A 1 keV barrier for protons would need 18 MV/m, six times the 3 MV/m at which air
breaks down [44]. Paul traps do hold ions with RF, but with hundreds of volts on millimetre electrodes and depths of a
few eV [26].

**Radio waves.** In air, fields add linearly, so no field blocks another. Active cancellation makes quiet pockets
only. In a diffuse field, a remote canceller gives a 10 dB quiet zone about λ/10 across [45]. That is an acoustics
result; applying it to radio is my inference. It means 6.9 cm at 433 MHz, 1.2 cm at 2.44 GHz and 0.5 cm at 5.8 GHz.
It needs a coherent copy of the unwanted signal, and the energy moves into brighter fringes elsewhere.

**Passive shields.** These work. Copper's skin depth is 2.1 µm at 1 GHz [46], so a closed metal box stops radio; its
seams and slots set the leakage. Mu-metal (relative permeability 80,000-100,000) diverts static and slow magnetic
fields; a five-layer box cut Earth's field 1,500-fold [47]. That shields fields, not particles.

### 3.3 The array's real shield: receive-side nulling

Weights on the four coherent channels make a pattern with zeros toward up to three interferers. Nothing is
transmitted, so no licence or exposure question arises.

**Steered nulls** point toward a known direction, and calibration limits their depth. The residual relative to one
element is about σ², the mean-square complex weight error (derived):

| Per-channel error | Residual |
| --- | --- |
| 0.15° / 0.02 dB | -49 dB |
| 1° / 0.1 dB | -34 dB |
| 3° / 0.5 dB (rev B's antenna-side residual) | -22 dB |

The local-fields run log gives -27 to -35 dB below the focus for a point null with 3°/0.26 dB errors.

**Adaptive nulls** are computed from the received data. They follow the interferer's actual array response, so
calibration error barely limits their depth. The limits are:
1. **Aperture and bandwidth.** One complex weight per channel nulls exactly at one frequency. Across a band B, with a
   path difference τ across the array, the residual is about (2πτB)²/12 (derived):
   - 42 mm four-whip ring at L-band: -60 dB at 2 MHz, -40 dB at 20 MHz;
   - 225 mm ring: -45 dB at 2 MHz, -25 dB at 20 MHz.
2. **Channel-to-channel filter ripple.** At 0.1-0.3 dB and 1-3° it leaves -25 to -35 dB unless each channel is
   equalised (est.).
3. **Dynamic range.** The interferer must not compress the front end: use the robust gain state; the QAMP compresses
   at -20 dBm (`REDESIGN.md` §5).
4. **Multipath.** Every strong echo of the interferer is another direction and uses up a null.

For comparison:
- A commercial 7-element GNSS anti-jam antenna claims up to 40 dB of extra protection [48].
- A fixed combiner null designed for 40 dB measured only 5-7 dB because of combiner errors and coupling [49]. That is
  the case for adaptive weights.

What it achieves for rev B's jobs:
- DF and passive radar keep working with up to three strong emitters nulled.
- GNSS reception keeps working near a jammer, using four whips.
- The radar table (`REDESIGN.md` §7) already counts on it: a 20 dB null on the transmit horn lifts the 2.44 GHz range
  from 240 m to 570 m, and the 915 MHz range from 100 m to 320 m.
- Real-time nulling uses the four simultaneous channels. The time-switched 16-whip layouts do not add real-time nulls.

## 4. What local fields can produce, and the limits

Far-field values (S = EIRP/4πr², E = √(377·S)). Inside about a wavelength they are only indicative (est.).

| Case | EIRP | At 0.3 m | At 1 m | Share of FCC public limit at 1 m [29] | FCC limit reached at |
| --- | --- | --- | --- | --- | --- |
| 4 whips, 433 MHz, full drive | 38.9 mW | 3.6 V/m, 17 nT peak | 1.08 V/m, 5.1 nT | 1.1×10^-3 | 33 mm |
| 4 whips, 915 MHz, full drive | 28.8 mW | 3.1 V/m | 0.93 V/m | 3.8×10^-4 | 19 mm |
| Panel, 2.44 GHz, full drive | 6.5 mW | 1.5 V/m | 0.44 V/m | 5×10^-5 | 7 mm |
| Panel, 5.8 GHz, full drive | 4.2 mW | 1.2 V/m | 0.35 V/m | 3×10^-5 | 6 mm |
| 15.249 level | 0.75 mW | 0.5 V/m | 0.15 V/m | 6×10^-6 | 2 mm |
| 15.245 level (hardware needs external gain) | 75 mW | 5.0 V/m | 1.5 V/m | 6×10^-4 | 24 mm |

The FCC public limit is f/1500 mW/cm² (f in MHz) for 300-1,500 MHz and 1 mW/cm² (10 W/m², 61 V/m) above 1.5 GHz,
averaged over 30 minutes [29]. For comparison, Earth's static field is 25,000-65,000 nT [32].

What these fields can produce:
- **Heat: negligible.**
  - At most 25 mW is conducted in total. Fully absorbed in 1 g of water it gives 0.006 K/s.
  - Fully absorbed in 1 g of tissue it would be 25 W/kg against the 1.6 W/kg limit [29]. That needs tissue touching
    an element, so keep hands about 20 cm from transmitting elements (est.); there the far-field estimate is at most
    2.7 % of the limit.
  - Below 1 mW no evaluation is needed at any distance [50]. An amateur station must evaluate exposure before
    transmitting [55].
  - A microwave oven may leak up to 5 mW/cm² at 5 cm in service [51]; its own power, roughly 1 kW (est.), is about
    40,000 times the array's 25 mW.
- **Sparks or plasma: no.** Air breaks down at about 3 MV/m [44]. The array makes 3.6 V/m at 0.3 m: 10^6 short in
  field, 10^12 in power.
- **Wireless power: microwatts.** At 915 MHz into a 0 dBi receiver, power at 1 m is:
  - -33 dBm (0.5 µW) at the 15.249 level;
  - -17 dBm (20 µW) at full drive, which needs 15.247-type modulation or an amateur licence (est.);
  - +4 dBm from a 4 W EIRP reader under 15.247 [54], or -5 dBm at 3 m.
- **Force: at most 84 pN** (section 2.3).
- **Signals, beams, nulls and radar echoes: yes.** That is what the array is for.

Transmit rules (`REDESIGN.md` §7 has the full table):

| Rule | Bands | Limit | Use here |
| --- | --- | --- | --- |
| 47 CFR 15.249 [52] | 902-928, 2400-2483.5, 5725-5875 MHz | 50 mV/m at 3 m = -1.25 dBm EIRP average | beacons, rotational Doppler, box tests |
| 47 CFR 15.245 [53] | 902-928, 2435-2465, 5785-5815 MHz | 500 mV/m at 3 m = 18.75 dBm EIRP | field-disturbance sensors (a Doppler on a fan qualifies), not perimeter protection |
| 47 CFR 15.247 [54] | 902-928, 2400-2483.5, 5725-5850 MHz | 1 W conducted, 6 dBi antenna | digital modulation or frequency hopping only |
| Amateur (with a licence) | 70, 33, 23, 13, 5 cm | high | exposure evaluation required [55] |
| Board calibration (`REDESIGN.md` §6) | outside 15.205 bands | ≤-60 dBm EIRP, or ≤-41 dBm in MANIFOLD mode | the probe atlas (section 6) |

## 5. The array as a space-weather instrument

Rather than imitate Earth's shield, the array can watch it work. `fields_study.md` covers the Sun's radio emission;
this section lists only the space-weather observables and the expected signals.

| Observable | Band and hardware | Expected signal | Reference values |
| --- | --- | --- | --- |
| Solar flux at 10.7 cm (F10.7) | 2.8 GHz, panel beam on the Sun | 100 sfu adds 0.43 % (0.019 dB) to the system noise (est.; G/T interpolated to -18.8 dB/K from `REDESIGN.md` §5); 19σ in 1 s at 20 MHz if the gain were stable, so it needs switching | 97 sfu at Penticton on 2 November 2026, 90-day mean 109 [57]; from below 50 to above 300 over a solar cycle [56] |
| Solar radio bursts | L-band, any whip channel | 10^4 sfu: +2.6 dB; 10^6 sfu: +19 dB (est.: 2 dBi whip, 200 K system) | 6 December 2006: above 10^6 sfu, and many GPS receivers lost track across the sunlit side [58, 59]; bursts are broadband, 10 MHz-10 GHz [60] |
| Ionospheric electron content (TEC) | GPS L1 + L2 (or L5) [61], one bladeRF per band (348 MHz apart, more than one AD9361 covers) | 0.162 m of delay per TECU at L1; L1-L2 difference 0.105 m (0.35 ns) per TECU; L1-L5 0.129 m | storms have exceeded 200 TECU over the US [60]: 32 m at L1 |
| Scintillation (S4) | L1 signal-to-noise at 50 Hz or more | S4 = standard deviation over mean of signal power, about 1 minute [60] | equatorial fades up to 40 dB [60] |
| Geomagnetic storms | not the array: an RM3100 magnetometer, as in HamSCI's Personal Space Weather Station [65] | K = 9 at Boulder means fluctuations of 500 nT or more [64] | about 1 % of Earth's field |

Notes:
- **A built-in Dicke switch.** Each panel element's SW_CAL arm can switch its LNA from the flare to the muted
  calibration tree, a matched load at room temperature. That gives a radiometer reference without new parts (est.;
  check the tree's match and that the muted synthesizer adds no noise).
- **GNSS software.** GNSS-SDR tracks GPS L1 C/A, L2C and L5 and Galileo E1 and E5a [62].
- **A preselector check.** L2 (1227.60 MHz) and L5 (1176.45 MHz) [61] sit at the B2/B3 crossover (LFCN-1200+ /
  HFCN-1200+); check which band passes them (not verified here).
- **Hardware delay.** The inter-radio delay that `REDESIGN.md` §6 measures (<0.1 ns) corresponds to about 0.3 TECU
  of hardware bias (est.).
- **What the array adds.** The array can confirm that a burst comes from the Sun's direction, and can null
  multipath for cleaner TEC (est.).
- **Where the two questions meet.** Coronal mass ejections carry twisted magnetic flux ropes. Storms are strongest
  when the field they bring points south for a long time [63].

## 6. Experiments worth doing here, legal and safe

1. **See the knot instead of imagining it: a field atlas by reciprocity.**
   - Set-up: a small probe antenna on a 3-D-printer gantry moves through a 10 cm cube in front of the panel. It sends
     the PN calibration code at ≤-41 dBm EIRP in MANIFOLD mode, at 5.8 GHz (outside 15.205 bands); the array only
     receives.
   - Reciprocity turns this into each element's pattern on the grid. Any weighting (a focus, a quiet zone, a vortex
     line, a knot attempt) can then be computed offline and its zero lines traced in measured data.
   - Size and time: a λ/8 grid is 6.5 mm, 4,096 points, about 35 minutes at 0.5 s per point (est.).
   - Signal: about -73 dBm per element at 0.3 m against -98 dBm of noise in 20 MHz, plus 30 dB from the 1023-chip
     code (est.).
   - Risk: the probe's cable scatters; run it along the axis with ferrites.
   - Expected result: straight vortex lines with rev B; knots need about 64 elements around the volume.
2. **Rotational Doppler of a desk fan.**
   - Set-up: the panel faces along the fan's axis at 0.3-0.5 m. A horn on the axis behind the fan sends CW at the
     15.249 level (or under 15.245 at 5785-5815 MHz). Put a foil tab on one blade.
   - Processing: a spatial DFT over the four quadrant channels, then 1 s FFTs.
   - Expected result: at 1200 rpm, lines at -20 Hz in the +1 bin and +20 Hz in the -1 bin, swapping when the fan
     reverses, while a single element shows both signs equally. Hadamard snapshots of all 16 elements then resolve
     |m| ≤ 5, with some crosstalk from the square ring's uneven spacing (the model's m = +5 bin shows the m = +1 line
     at -50 dB).
   - Safety: keep the fan guard on.
3. **The adaptive-null shield.**
   - Set-up: one 15.249-level beacon or a Wi-Fi access point as the interferer, a second beacon as the wanted signal.
     Null the interferer with four whips or the panel quadrants.
   - Measure: suppression against bandwidth (2-20 MHz), number of snapshots, and calibration state. Compare steered
     nulls (about 22-34 dB) with adaptive ones (20-40 dB on one narrowband interferer, less with wide bands or
     multipath; est.).
4. **What really shields radio: a Faraday-box test.**
   - Set-up: a battery-powered 15.249-level beacon (a cable into the box would itself leak) inside a cookie tin, a
     steel box and a mesh cage, then with a 1 mm slot cut. The array outside measures the drop at 915, 2440 and
     5800 MHz.
   - Expected result: seams and slots, not wall thickness, set the result, since the skin depth is about 2 µm.
5. **Watch Earth's shield (receive-only).**
   - Daily 2.8 GHz drift scans with the Dicke switch, compared with Penticton's F10.7.
   - The L-band noise floor on all channels during flare alerts.
   - L1/L2 or L1/L5 TEC with GNSS-SDR.
   - Optionally, an RM3100 magnetometer for storms.
6. **See why scale decides, with a cloud chamber and a magnet.**
   - A dry-ice and isopropanol cloud chamber shows mostly muon tracks [66].
   - A 0.5 T NdFeB magnet under it curls low-energy electron tracks: a 100 keV electron circles on 2.2 mm. Muon tracks
     stay straight (27 m radius at 0.5 T) (derived).
   - Safety: dry ice at -78 °C needs gloves and goggles, and the CO₂ needs a ventilated room [66]; isopropanol is
     flammable.
7. **Tweezers done right: ultrasound.**
   - Set-up: a 40 kHz, 64-transducer phased array driven by an Arduino Mega (Ultraino) [25] makes twin, vortex and
     cage traps that levitate beads [24]. The repo's optimiser (`array16.evolve`) could drive it with an acoustic
     propagation model in place of the radio one: the same structured-field ideas, in a medium where they produce
     force.
   - Safety: holding a 2 mm polystyrene bead needs about 140 dB at the trap (est.). The 1984 guideline for airborne
     ultrasound above 20 kHz is 110 dB; 140-150 dB vibrates ear-canal hairs, and above about 155 dB the body heats
     [67]. Keep ears away, keep sessions short, and wear hearing protection.

Not worth doing: torque or levitation with radio (piconewtons and 10^-12 N·m), any "shield" made of radio waves, and
anything above the legal transmit levels.

Open points:
- the AD9361's uncorrelated synthesizer noise 20-200 Hz from the carrier;
- the whips' patterns at 2.44 and 5.8 GHz for ring experiments;
- the B2/B3 preselector edge near 1.2 GHz;
- the λ/10 quiet-zone size for a vector field and this array's geometry;
- the knot results in `local_fields.md` once it is published.

## 7. Sources

1. Dennis, King, Jack, O'Holleran, Padgett, Isolated optical vortex knots, Nature Physics 6, 118 (2010): https://eprints.gla.ac.uk/41285
2. Pires et al., Stability of optical knots in atmospheric turbulence, Nature Communications 16 (2025): https://ideas.repec.org/a/nat/natcom/v16y2025i1d10.1038_s41467-025-57827-1.html
3. Tsvetkov, Pires, Litchinitser, vortex-line robustness in turbulence (preprint, 2025): https://arxiv.org/abs/2511.08836
4. Kleckner, Scheeler, Irvine, The life of a vortex knot: https://arxiv.org/abs/1310.3321
5. Irvine and Bouwmeester, Linked and knotted beams of light, Nature Physics 4, 716 (2008): https://irvinelab.uchicago.edu/papers/nphys1056.pdf
6. Kedia, Bialynicki-Birula, Peralta-Salas, Irvine, Tying knots in light fields, PRL 111, 150404 (2013): https://arxiv.org/abs/1302.0342
7. Taylor, Relaxation of toroidal plasma and generation of reverse magnetic fields, PRL 33, 1139 (1974): https://ui.adsabs.harvard.edu/abs/1974PhRvL..33.1139T/abstract
8. Ji, Prager, Sarff, Conservation of magnetic helicity during plasma relaxation, PRL 74, 2945 (1995): https://w3.pppl.gov/~hji/paper/helicity_abs.html
9. Smiet et al., Self-organizing knotted magnetic structures in plasma, PRL 115, 095001 (2015): https://arxiv.org/abs/1507.08780
10. Skyrmion racetrack memory (building on Fert, Cros, Sampaio, Nature Nanotechnology 8, 152, 2013): https://arxiv.org/pdf/1907.03449
11. He, Friese, Heckenberg, Rubinsztein-Dunlop, PRL 75, 826 (1995): https://doi.org/10.1103/PhysRevLett.75.826
12. Nobel Prize in Chemistry 2014, press release: https://www.nobelprize.org/prizes/chemistry/2014/press-release/
13. Kuga et al., Novel optical trap of atoms with a doughnut beam, PRL 78, 4713 (1997): https://link.aps.org/doi/10.1103/PhysRevLett.78.4713
14. Padgett, A new twist on the Doppler shift, Physics Today 67(2) (2014): https://physicstoday.aip.org/quick-study/a-new-twist-on-the-doppler-shift
15. Courtial, Robertson, Dholakia, Allen, Padgett, Rotational frequency shift of a light beam, PRL 81, 4828 (1998): https://link.aps.org/doi/10.1103/PhysRevLett.81.4828
16. Lavery, Speirits, Barnett, Padgett, Detection of a spinning object using light's orbital angular momentum, Science 341, 537 (2013): https://eprints.gla.ac.uk/85025
17. Zhao et al., rotational Doppler shift with a radio-frequency OAM beam, Optics Letters 41, 2549 (2016): https://opg.optica.org/abstract.cfm?URI=ol-41-11-2549
18. Liu, Chu, Giddens, Li, Hao, Linear and rotational Doppler shifts from designer surfaces, Scientific Reports (2019): https://pmc.ncbi.nlm.nih.gov/articles/PMC6586861/
19. Wang et al., Optical vortex metrology for nanometric speckle displacement measurement, Optics Express 14, 120 (2006): https://doi.org/10.1364/OPEX.14.000120
20. Edfors and Johansson, Is OAM based radio communication an unexploited area?, IEEE TAP 60, 1126 (2012): https://lup.lub.lu.se/record/2062936
21. Yan et al., High-capacity millimetre-wave communications with OAM multiplexing, Nature Communications 5, 4876 (2014): https://pmc.ncbi.nlm.nih.gov/articles/PMC4175588
22. Larocque et al., Optical framed knots as information carriers, Nature Communications 11, 5119 (2020): https://pmc.ncbi.nlm.nih.gov/articles/PMC7547656
23. Wang et al., Free-space propagation and skyrmion topology of toroidal electromagnetic pulses (2023): https://arxiv.org/abs/2311.01765
24. Marzo et al., Holographic acoustic elements for manipulation of levitated objects, Nature Communications 6, 8661 (2015): https://research-information.bris.ac.uk/en/publications/holographic-acoustic-elements-for-manipulation-of-levitated-objec/
25. Marzo, Corkett, Drinkwater, Ultraino: an open phased-array system for airborne ultrasound, IEEE TUFFC 65 (2018): https://research-information.bris.ac.uk/en/publications/ultraino-an-open-phased-array-system-for-narrowband-airborne-ultr/
26. Quadrupole ion trap: https://en.wikipedia.org/wiki/Quadrupole_ion_trap
27. Van Allen radiation belt: https://en.wikipedia.org/wiki/Van_Allen_radiation_belt
28. Emile et al., torque on a 30 cm ring in the microwave range, PRL 112, 053902 (2014): https://hal-univ-rennes1.archives-ouvertes.fr/hal-00998554
29. 47 CFR 1.1310, RF exposure limits: https://www.law.cornell.edu/cfr/text/47/1.1310
30. Monopulse angle accuracy: https://rfessentials.com/rf-knowledge-base/how-do-i-design-a-monopulse-tracking-radar-for-precision-angle-measurement/
31. LHC Machine, JINST 3 S08001, chapter 3 (magnets): https://jinst.sissa.it/LHC/LHCmachine/ch03.pdf
32. NOAA NCEI geomagnetism FAQ: https://www.ncei.noaa.gov/products/geomagnetism-frequently-asked-questions
33. Ambroglini, Battiston, Burger, superconducting magnet shields for long missions, Frontiers in Oncology 6:97 (2016): https://pmc.ncbi.nlm.nih.gov/articles/PMC4896949
34. ECSS-E-ST-10-04C via SPENVIS, geomagnetic cutoff: https://www.spenvis.oma.be/ecss/e_st_10_04c/ecss_13_08_00.htm
35. Magnetopause: https://en.wikipedia.org/wiki/Magnetopause
36. NOAA SWPC, solar wind: https://www.spaceweather.gov/phenomena/solar-wind
37. CERN Courier on Bamford et al., PRL 109, 081101 (2012), lunar mini-magnetospheres: https://cerncourier.com/a/deflector-shields-protect-the-lunar-surface/
38. Particle Data Group, cosmic rays review (2020): https://pdg.lbl.gov/2020/reviews/rpp2020-rev-cosmic-rays.pdf
39. NASA, Earth's magnetosphere: protecting our planet from harmful space energy: https://science.nasa.gov/science-research/earth-science/earths-magnetosphere-protecting-our-planet-from-harmful-space-energy/
40. Ferrari and Szuszkiewicz, Cosmic ray recipes (2006): https://arxiv.org/abs/astro-ph/0601158
41. ICRPaedia, sources of radiation exposure (UNSCEAR 2008): https://icrpaedia.org/Sources_of_Radiation_Exposure
42. US NRC, types of radiation: https://www.nrc.gov/reading-rm/basic-ref/students/science-101/what-are-different-types-of-radiation.html
43. Ponderomotive force: https://en.wikipedia.org/wiki/Ponderomotive_force
44. Dielectric strength: https://en.wikipedia.org/wiki/Dielectric_strength
45. Rafaely, Zones of quiet in a broadband diffuse sound field (2023): https://arxiv.org/abs/2310.04191
46. Skin effect: https://en.wikipedia.org/wiki/Skin_effect
47. Mu-metal: https://en.wikipedia.org/wiki/Mu-metal
48. NovAtel GAJT-710ML: https://www.gnss.ca/brochures/gajt-710ml-ps.html
49. Kunysz, Advanced Pinwheel compact CRPA, ION 2001: https://hexagondownloads.blob.core.windows.net/public/Novatel/assets/Documents/Papers/ion_2001/ion_2001.pdf
50. 47 CFR 1.1307: https://www.law.cornell.edu/cfr/text/47/1.1307
51. 21 CFR 1030.10, microwave ovens: https://www.law.cornell.edu/cfr/text/21/1030.10
52. 47 CFR 15.249: https://www.law.cornell.edu/cfr/text/47/15.249
53. 47 CFR 15.245: https://www.law.cornell.edu/cfr/text/47/15.245
54. 47 CFR 15.247: https://www.law.cornell.edu/cfr/text/47/15.247
55. 47 CFR 97.13: https://www.law.cornell.edu/cfr/text/47/97.13
56. NOAA SWPC, F10.7 radio emissions: https://www.spaceweather.gov/phenomena/f107-cm-radio-emissions
57. SpaceWeatherLive, solar activity report for 2 November 2026: https://www.spaceweatherlive.com/en/archive/2026/11/02/rsga.html
58. Cerruti et al., Effect of intense December 2006 solar radio bursts on GPS receivers, Space Weather 6 (2008): https://digitalcommons.njit.edu/fac_pubs/12574
59. NOAA via ScienceDaily (2007): https://www.sciencedaily.com/releases/2007/04/070407145619.htm
60. Kintner, A beginner's guide to space weather and GPS (Cornell): https://gps.ece.cornell.edu/SpaceWeatherIntro_ed2_10-31-06_ed.pdf
61. GPS signals: https://en.wikipedia.org/wiki/GPS_signals
62. GNSS-SDR: https://github.com/gnss-sdr/gnss-sdr
63. NOAA SWPC, coronal mass ejections: https://www.spaceweather.gov/phenomena/coronal-mass-ejections
64. NOAA SWPC, station K and A indices: https://www.spaceweather.gov/products/station-k-and-indices
65. HamSCI, RM3100 magnetometer for Personal Space Weather Stations (2026): https://hamsci.org/bibcite/reference/548 (kit: https://tapr.org/product/tangerine-sdr-magnetometer)
66. CERN S'Cool LAB, cloud chamber: https://scoollab.web.cern.ch/cloud-chamber
67. ICNIRP statement on the 1984 interim guidelines on airborne ultrasound (2024): https://icnirp.org/cms/upload/publications/ICNIRPUltrasoundStatement2024.pdf
