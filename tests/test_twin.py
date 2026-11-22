"""array16.twin, the twin bench's drivers (docs/KNOTTED_GUIDE.md section 8.3): design files read as the guide document
draws R1, T915-S, the flat ring and the (1, 3) helix; solves cached by the hash of their content and solver; a
realisation's tune, drive and field() from the solved section, the ring of segments and the swept section; the
detector's counts of held and lost electrons, its pickups' currents and its dump; the vacuum's pressure, free path and
collisions; and the bench that names them."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import knotted_guide_study as study
import numpy as np
import pytest

import array16.twin as twin
from array16 import particles as pt
from array16.guide import sections as gs
from array16.lab import Bench
from array16.particles import FacingBoards, RailPair, TravellingWave
from array16.twin import (
    Design,
    RingWave,
    TwinDetector,
    TwinGuide,
    TwinVacuum,
    framing,
    section_key,
)

ROOT = Path(__file__).resolve().parents[1]
DESIGNS = ROOT / "designs"
F915 = 915e6
H915 = 0.915e-3  # T915-S's h (m)
STRAIGHT = """
name = "straight facing strips"
seed = 1
[band]
f = 915e6
m = 1
[section]
kind = "facing"
gap = 1.83e-3
width = 7.32e-3
split = 0.2745e-3
margin = 1.83e-3
laminate = "polyimide"
thickness = 25e-6
copper = 18e-6
grid = "coarse"
[path]
kind = "straight"
length = 0.1
samples = 64
[drive]
volts = 50.0
"""


def _quick(name: str, where: Path) -> Path:
    """A design file on the coarse grid, beside the test."""
    text = (DESIGNS / f"{name}.toml").read_text()
    assert 'grid = "default"' in text
    out = where / f"{name}.toml"
    out.write_text(text.replace('grid = "default"', 'grid = "coarse"'))
    return out


@pytest.fixture(scope="module")
def quick(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A folder of the four designs on the coarse grid, with one cache for the module's guides."""
    where = tmp_path_factory.mktemp("designs")
    for name in ("r1_rail_ring", "t915s", "flat_ring_318", "helix_13"):
        _quick(name, where)
    (where / "straight.toml").write_text(STRAIGHT)
    return where


def _guide(quick: Path, name: str, realisation: int = 0) -> TwinGuide:
    return TwinGuide(quick / f"{name}.toml", realisation, cache=quick / "cache")


# ---- the design files -----------------------------------------------------------------------------------------------
def test_the_designs_are_read_as_the_guide_document_draws_them() -> None:
    # R1 (section 3.6): rails 1.46 mm either side of a 0.40 mm grounded strip on 1.6 mm FR-4 with 35 um copper, at
    # 50 V and m = 1 at 915 MHz; T915-S, the flat ring and the helix (section 9): section B's facing flex strips, 8h
    # wide split 0.3h at h = 0.915 mm on 25 um polyimide, at the same volts, on the r/R = 0.45 trefoil framed k = 0, a
    # circle bent in the boards' plane, and the (1, 3) torus curve at r/R = 0.2 framed k = 0
    r1 = Design.load(DESIGNS / "r1_rail_ring.toml")
    assert r1.section.geometry == gs.Rails(0.40e-3, 1.46e-3, 0.15e-3, backed=False, walls=43.4)
    assert (r1.section.laminate, r1.section.thickness, r1.section.copper) == (gs.FR4, 1.6e-3, gs.RIGID_COPPER)
    assert (r1.f, r1.m, r1.volts, r1.q, r1.feeds) == (F915, 1, 50.0, None, 2)
    assert (r1.path.kind, r1.path.length, r1.frame.theta, r1.frame.k) == ("circle", None, math.pi, None)
    assert (r1.tolerance.eps, r1.tolerance.earth, r1.tolerance.spacers, r1.tolerance.joint_dz) == (0.1, 50e-6, 0, 0)
    assert r1.source == "designs/r1_rail_ring.toml"
    assert len(r1.digest) == 64
    facing = [Design.load(DESIGNS / f"{name}.toml") for name in ("t915s", "flat_ring_318", "helix_13")]
    assert facing[0].section == facing[1].section == facing[2].section
    assert len({section_key(d.section) for d in facing}) == 1  # one solve serves the three
    b = gs.section_b().geometry
    assert isinstance(b, gs.FacingStrips)
    geometry = facing[0].section.geometry
    assert isinstance(geometry, gs.FacingStrips)
    for name in ("gap", "width", "split", "margin"):
        assert getattr(geometry, name) == pytest.approx(getattr(b, name), rel=1e-12)
    assert facing[0].section.laminate == gs.POLYIMIDE
    assert (facing[0].section.thickness, facing[0].section.copper) == (25e-6, gs.FLEX_COPPER)
    assert {d.volts for d in facing} == {50.0}
    assert (facing[0].path.kind, facing[0].path.ratio, facing[0].path.p, facing[0].path.q) == ("torus", 0.45, 2, 3)
    assert (facing[2].path.kind, facing[2].path.ratio, facing[2].path.p, facing[2].path.q) == ("torus", 0.2, 1, 3)
    assert facing[0].frame.k == facing[2].frame.k == 0
    assert (facing[1].path.kind, facing[1].frame.theta) == ("circle", math.pi)
    assert facing[0].tolerance.spacers == 0.02
    assert facing[0].tolerance.stations == 24


def test_a_design_file_says_what_it_takes_when_it_is_refused() -> None:
    base = (DESIGNS / "t915s.toml").read_text()
    for change, message in (
        (base + "\ncolour = 'gold'\n", r"unknown keys \['colour'\]"),
        (base.replace("volts = 50.0", "volts = 50.0\nq = 0.5"), "the volts or the null's q"),
        (base.replace('laminate = "polyimide"', 'laminate = "kapton"'), "laminate is one of"),
        (base.replace('kind = "facing"', 'kind = "coax"'), "kind is one of"),
        (base.replace('kind = "torus"', 'kind = "straight"'), "a straight path has no holonomy"),
        (base.replace("ratio = 0.45", "ratio = 1.5"), "0 < ratio"),
        (base.replace("split = 0.2745e-3", "split = 9e-3"), "a split narrower than the strips"),
        (base.replace("[drive]", "[drive]\nphase = 0"), r"\[drive\]: unknown keys"),
    ):
        with pytest.raises(ValueError, match=message):
            Design.parse(change, "t915s.toml")


def test_the_framings_are_the_knotted_studys() -> None:
    # the closing framing is the study's (section 3.5), arithmetic and all: on T915-S's trefoil the k = 0 start keeps
    # the boards' in-plane radius at 18.57 mm, and on the (1, 3) helix at 29.14 mm (run.log)
    for ratio, p, q, n, radius in ((0.45, 2, 3, 4000, 18.57e-3), (0.2, 1, 3, 3000, 29.14e-3)):
        curve = study.torus(ratio, study.L915, n, p, q)
        ours, theirs = framing(curve, 0), study.framing(curve, 0)
        assert ours.theta0 == theirs.theta0
        assert np.array_equal(ours.across, theirs.across)
        assert np.array_equal(ours.kappa_g, theirs.kappa_g)
        assert ours.inplane_radius == pytest.approx(radius, abs=5e-6)


# ---- solves cached by content ---------------------------------------------------------------------------------------
def test_solves_are_cached_by_their_content_and_solver(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cache = tmp_path / "cache"
    path = _quick("t915s", tmp_path)
    first = TwinGuide(path, cache=cache).solved
    key = section_key(first.section)
    assert sorted(p.name for p in cache.iterdir()) == [f"{key}.npz"]

    def refuse(section: gs.Section) -> gs.Solved:
        raise AssertionError("solved again")

    with monkeypatch.context() as patch:
        patch.setattr(twin, "solve", refuse)
        again = TwinGuide(path, 3, cache=cache).solved  # another realisation, the same section
        flat = TwinGuide(_quick("flat_ring_318", tmp_path), cache=cache).solved  # another design, the same section
    for read in (again, flat):
        assert all(np.array_equal(getattr(read, a), getattr(first, a)) for a in ("u", "v", "phi"))
        assert (read.index, read.z0, read.qc, read.qd, read.share) == (
            first.index,
            first.z0,
            first.qc,
            first.qd,
            first.share,
        )
        assert (read.null, read.saddle, read.gradient, read.escape) == (
            first.null,
            first.saddle,
            first.gradient,
            first.escape,
        )
        assert not read.phi.flags.writeable
    other = path.read_text().replace("split = 0.2745e-3", "split = 0.30e-3")
    (tmp_path / "wider.toml").write_text(other)
    wider = TwinGuide(tmp_path / "wider.toml", cache=cache).solved
    assert section_key(wider.section) != key
    assert len(list(cache.iterdir())) == 2
    (cache / f"{key}.npz").write_bytes(b"not a solve")  # a damaged file is solved again and written anew
    assert TwinGuide(path, cache=cache).solved.index == first.index
    with monkeypatch.context() as patch:
        patch.setattr(twin, "solver_version", lambda: "another solver")
        assert section_key(first.section) != key


# ---- tune, drive and field ------------------------------------------------------------------------------------------
def test_the_design_as_drawn_resonates_at_its_band_and_takes_the_power_of_its_shunt(quick: Path) -> None:
    # realisation 0 is the ring as drawn: m guided waves round, so its +-1 pair sits on the band, f/(1 + 1/4Q^2) with
    # the line's own Q, unsplit; two compiled feeds a quarter wave apart drive it at 50 V mean for mean |V|^2/(2R),
    # R = Z0 Q/(2 pi m) (RNG-2), the null at q 0.51 and a 0.09 eV well (R1, section 3.6), the backward wave under 1e-3
    guide = _guide(quick, "r1_rail_ring")
    solved, tuning = guide.solved, guide.tune()
    q0 = solved.quality
    assert guide.length == pytest.approx(solved.wavelength, rel=1e-12)
    assert tuning.f == pytest.approx(F915 / (1 + 1 / (4 * q0 * q0)), rel=1e-9)
    assert tuning.q == pytest.approx(q0, rel=1e-9)
    assert abs(tuning.splitting) < 1e-12
    assert tuning.shunt == pytest.approx(solved.z0 * q0 / (2 * math.pi), rel=1e-12)
    driven = guide.drive()
    assert driven.feeds == (0, 48)
    assert driven.volts == pytest.approx(50.0, rel=1e-12)
    assert driven.watts == pytest.approx(float(np.mean(np.abs(driven.profile) ** 2)) / (2 * tuning.shunt), rel=1e-3)
    assert 7.5 < driven.watts < 9.0
    assert driven.q == pytest.approx(pt.mathieu_q(solved.gradient * 50.0, tuning.f), rel=1e-12)
    assert driven.q == pytest.approx(0.505, abs=0.005)
    assert 0.085 < driven.depth_ev < solved.depth_ev(50.0) < 0.10
    assert driven.purity < 1e-3
    assert 0.85 < driven.depth_ratio < 0.95  # Q 36 loses 9 % of the depth between feeds a quarter wave apart
    same = guide.drive(watts=driven.watts)
    assert np.allclose(same.profile, driven.profile, rtol=1e-12, atol=0)
    assert guide.drive(volts=100.0).watts == pytest.approx(4 * driven.watts, rel=1e-12)
    assert guide.drive(q=0.5).q == pytest.approx(0.5, rel=1e-12)
    three = guide.drive(feeds=3)
    assert three.feeds == (0, 64, 128)
    assert three.purity < 1e-3


def test_a_realisation_draws_what_its_design_says(quick: Path) -> None:
    # T915-S: the spacers' 2 % rms in dh/h, the joints, the etch and the laminate split its pair; the same number draws
    # the same ring; Earth's field is 50 uT in a drawn orientation. R1 on one board draws only its laminate and Earth's
    # field: its impedance stays uniform and its resonance moves as its index, by FR-4's 0.1 in eps_r over the share
    drawn = [_guide(quick, "t915s", k).drawn for k in (1, 1, 2)]
    for name in ("z", "n", "rel", "earth"):
        assert np.array_equal(getattr(drawn[0], name), getattr(drawn[1], name))
    assert not np.array_equal(drawn[0].z, drawn[2].z)
    for one in drawn:
        assert 0.01 < float(np.std(one.rel)) < 0.03
        assert float(np.linalg.norm(one.earth)) == pytest.approx(50e-6, rel=1e-12)
    knot = _guide(quick, "t915s", 1)
    assert knot.tune().linewidths > 0.1
    driven = knot.drive()
    # the weakest well: the barrier field of the drawn gap h (1 + dh/h) at each node's volts, V escape / (1 + dh/h)
    barrier = float(np.min(np.abs(driven.profile) / (1 + knot.drawn.rel)))
    assert driven.depth_ev == pytest.approx(pt.ponderomotive_ev(knot.solved.escape * barrier, driven.f), rel=1e-12)
    assert float(np.std(knot.drawn.rel)) > 0.01  # the gaps' draw moves it
    wave = knot.field(driven).swept.drive
    assert isinstance(wave, RingWave)
    nodes = len(driven.profile)
    at_nodes = wave.phasor(2 * math.pi * np.arange(nodes) / nodes)  # the gap folded in as (1 + dh/h)^-2
    assert np.allclose(at_nodes, driven.profile / (1 + knot.drawn.rel) ** 2, rtol=1e-9, atol=0)
    nominal = _guide(quick, "r1_rail_ring")
    assert _guide(quick, "r1_rail_ring", 0).field().earth is None
    guide = _guide(quick, "r1_rail_ring", 1)
    real = guide.drawn
    assert np.array_equal(real.z, np.full(len(real.z), nominal.solved.z0))
    assert not real.rel.any()
    sd = 0.1 * nominal.solved.index * nominal.solved.share / (2 * gs.FR4.er)
    assert guide.build.index_per_eps * 0.1 == pytest.approx(sd, rel=1e-12)
    shift = float(np.mean(real.n)) / nominal.solved.index - 1
    assert 0 < abs(shift) < 5 * sd / nominal.solved.index
    assert guide.tune().f / nominal.tune().f - 1 == pytest.approx(-shift / (1 + shift), rel=1e-6)


def test_the_field_is_the_swept_section_driven_by_the_realised_ring(quick: Path) -> None:
    # R1 as drawn: on the path, its null, nothing; 0.1h above it, E = -V(s, t) grad phi turned into the frame (u inward,
    # v up), V the ring's node volts continued round the lap; B = T x E / v at the travelling wave's phase velocity;
    # a realisation adds its uniform Earth's field to B
    guide = _guide(quick, "r1_rail_ring")
    driven = guide.drive()
    field = guide.field(driven)
    curve, grid = guide.curve, guide.grid()
    assert curve is not None
    at = np.arange(0, len(curve), 97)
    e, b = field(curve.points[at], 1.3e-10)
    assert e is not None
    assert b is not None
    h = guide.solved.distance
    up = np.array([0.0, 0.0, 1.0])
    x = curve.points[at] + 0.1 * h * up
    t = 2.1e-10
    e, b = field(x, t)
    assert e is not None
    assert b is not None
    wave = RingWave.of(driven.profile, 2 * math.pi * driven.f, 1)
    volts = wave(2 * math.pi * at / len(curve), t)
    du, dv, on = grid.gradient_at(np.full(len(at), grid.null[0]), np.full(len(at), grid.null[1] + 0.1 * h))
    assert on.all()
    inward = -curve.points[at] / np.linalg.norm(curve.points[at], axis=1)[:, None]
    expected = -volts[:, None] * (du[:, None] * inward + dv[:, None] * up)
    scale = float(np.abs(expected).max())
    assert np.abs(e - expected).max() < 1e-6 * scale
    assert np.abs(field(curve.points[at], t)[0]).max() < 1e-9 * scale  # type: ignore[arg-type]
    speed = driven.f * guide.length
    tangent = curve.tangents[at]
    assert np.allclose(b, np.cross(tangent, e) / speed, rtol=0, atol=1e-12 * scale / speed)
    earthy = _guide(quick, "r1_rail_ring", 2)
    shaken = earthy.field()
    assert shaken.earth is not None
    _, b_earth = shaken(x, t)
    _, b_bare = earthy.sweep(earthy.grid(), RingWave.of(earthy.drive().profile, 2 * math.pi * earthy.tune().f, 1))(x, t)
    assert b_earth is not None
    assert b_bare is not None
    assert np.allclose(b_earth - b_bare, shaken.earth, rtol=1e-12, atol=1e-18)


def test_a_ring_of_m_guided_waves_is_m_wavelengths_round(quick: Path) -> None:
    # T915-S's section at m = 2: twice the guided wavelength round, its pair on the band, its two feeds an eighth of the
    # lap apart (a quarter guided wave), compiled for the travelling wave of order 2
    text = (quick / "t915s.toml").read_text().replace("m = 1", "m = 2")
    (quick / "t915s_m2.toml").write_text(text)
    guide = _guide(quick, "t915s_m2")
    assert guide.length == pytest.approx(2 * guide.solved.wavelength, rel=1e-12)
    q0 = guide.solved.quality
    assert guide.tune().f == pytest.approx(F915 / (1 + 1 / (4 * q0 * q0)), rel=1e-9)
    driven = guide.drive()
    assert driven.feeds == (0, 24)
    assert driven.purity < 1e-3
    assert guide.tune().shunt == pytest.approx(guide.solved.z0 * q0 / (4 * math.pi), rel=1e-12)


def test_the_section_grid_marks_the_copper_and_the_boards_solid(quick: Path) -> None:
    # R1's rails: below the copper's top face every cell is copper or FR-4 but the etched gaps, 0.15 mm centred on the
    # rails' edges 0.20 and 1.66 mm out, and above it nothing; section B: both boards' copper and polyimide, the split
    # open, the gap between them open
    rails = _guide(quick, "r1_rail_ring").grid()
    uc, vc = (rails.u[1:] + rails.u[:-1]) / 2, (rails.v[1:] + rails.v[:-1]) / 2
    assert rails.solid is not None
    assert rails.solid[:, vc < -35e-6].all()
    assert not rails.solid[:, vc > 0].any()
    layer = rails.solid[:, (vc > -35e-6) & (vc < 0)]
    gaps = (np.abs(np.abs(uc) - 0.20e-3) < 0.075e-3) | (np.abs(np.abs(uc) - 1.66e-3) < 0.075e-3)
    assert np.array_equal(layer, np.repeat(~gaps[:, None], layer.shape[1], axis=1))
    strips = _guide(quick, "t915s").grid()
    uc, vc = (strips.u[1:] + strips.u[:-1]) / 2, (strips.v[1:] + strips.v[:-1]) / 2
    assert strips.solid is not None
    gap, copper, film = 2 * H915, 18e-6, 25e-6
    assert not strips.solid[:, (vc > 0) & (vc < gap)].any()
    for lo, hi in ((-copper - film, -copper), (gap + copper, gap + copper + film)):
        assert strips.solid[:, (vc > lo) & (vc < hi)].all()
    for lo, hi in ((-copper, 0.0), (gap, gap + copper)):
        layer = strips.solid[:, (vc > lo) & (vc < hi)]
        split = np.abs(uc) < 0.3 * H915 / 2
        assert np.array_equal(layer, np.repeat(~split[:, None], layer.shape[1], axis=1))


def test_the_ring_wave_is_the_travelling_wave_its_phasors_make() -> None:
    # node volts V0 e^(i (phase0 - m angle)) in the ring model's convention, e^(i omega t), are the travelling wave
    # V0 cos(m angle - omega t - phase0), between the nodes too, and carry its slowness
    nodes, omega = 192, 2 * math.pi * F915
    angle = 2 * math.pi * np.arange(nodes) / nodes
    probe = np.random.default_rng(5).uniform(0, 2 * math.pi, 500)
    for m, phase in ((1, 0.3), (2, -1.1)):
        wave = RingWave.of(np.asarray(50.0 * np.exp(1j * (phase - m * angle)), dtype=np.complex128), omega, m)
        exact = TravellingWave(50.0, omega, m, -phase)
        for t in (0.0, 1.3e-10, 7.7e-10):
            assert np.abs(wave(probe, t) - exact(probe, t)).max() < 1e-5 * 50.0
        assert wave.slowness(0.2) == pytest.approx(exact.slowness(0.2), rel=1e-12)


def test_the_gapless_model_is_the_closed_form_swept(quick: Path) -> None:
    # E2-21 asks the ELE-1 field swept round the ring: the gapless model samples RailPair h/40 apart with its null at
    # sqrt(c1 c2), and its swept field 0.2h round the null is RailPair's to the interpolation (0.3 % of the largest
    # there), within ELE-1's 0.5 %; the flat ring's gapless boards are ELE-2's FacingBoards with the null midway
    flat = _guide(quick, "flat_ring_318")
    boards = FacingBoards(2 * H915)
    assert flat.grid("gapless").null == (0.0, H915)
    sweep = flat.sweep(flat.grid("gapless"), pt.InPhase(1.0, 0.0))
    ring = flat.curve
    assert ring is not None
    angle = np.linspace(0, 2 * math.pi, 24, endpoint=False)
    a, b = 0.3 * H915 * np.cos(angle), 0.3 * H915 * np.sin(angle)
    inward = -ring.points[7] / np.linalg.norm(ring.points[7])
    up = np.array([0.0, 0.0, 1.0])
    e, _ = sweep(ring.points[7] + a[:, None] * inward + b[:, None] * up, 0.0)
    assert e is not None
    closed = boards.field(a, b, 1.0)
    assert np.abs(np.column_stack([e @ inward, e @ up]) - closed).max() < 5e-3 * np.abs(closed).max()
    guide = _guide(quick, "r1_rail_ring")
    pair = RailPair(0.40e-3, 1.46e-3)
    h = pair.null_height
    grid = guide.grid("gapless")
    assert grid.null == (0.0, h)
    field = guide.sweep(grid, pt.InPhase(1.0, 0.0))
    curve = guide.curve
    assert curve is not None
    a, b = 0.2 * h * np.cos(angle), 0.2 * h * np.sin(angle)
    inward = -curve.points[0] / np.linalg.norm(curve.points[0])
    e, _ = field(curve.points[0] + a[:, None] * inward + b[:, None] * up, 0.0)
    assert e is not None
    closed = pair.field(a, h + b, 1.0)
    swept = np.column_stack([e @ inward, e @ up])
    assert np.abs(swept - closed).max() < 5e-3 * np.abs(closed).max()
    with pytest.raises(ValueError, match="'solved' or 'gapless'"):
        guide.grid("ideal")  # type: ignore[arg-type]


def test_a_knot_is_swept_in_its_closing_frame(quick: Path) -> None:
    # T915-S framed k = 0: u runs across the boards' plane and v along minus the board normal, the frame closes on
    # itself after a lap (seam 0), the null stays on the path at every sample, and the in-plane radius the frame keeps
    # is the study's 18.6 mm
    guide = _guide(quick, "t915s")
    assert guide.framing is not None
    assert guide.framing.inplane_radius == pytest.approx(18.57e-3, rel=2e-3)
    field = guide.field()
    assert np.allclose(field.swept.first, guide.framing.across, rtol=0, atol=1e-12)
    assert np.allclose(field.swept.second, -guide.framing.normal, rtol=0, atol=1e-12)
    assert abs(field.swept.seam) < 1e-9
    curve = guide.curve
    assert curve is not None
    e_on, _ = field(curve.points[::50], 3e-10)
    off = curve.points[::50] + 0.05 * H915 * guide.framing.normal[::50]
    e_off, _ = field(off, 3e-10)
    assert e_on is not None
    assert e_off is not None
    assert np.abs(e_on).max() < 1e-9 * np.abs(e_off).max()
    assert not field.off_grid.any()


# ---- the detector ---------------------------------------------------------------------------------------------------
def test_the_detector_counts_the_held_and_the_lost(quick: Path) -> None:
    # R1 as drawn, 5 RF periods at 120 steps: electrons at rest within 0.05h of the null are held; electrons sent up at
    # 5 eV, fifty times the well, leave the aperture, and those sent down at 100 eV, past the 50 V between the rails and
    # the strip that turns slower ones back within a fifth of a period, strike the copper; each is lost with the step,
    # time and arc length it left at; the dump takes the held charge, 20 e on 1 pF, and the loss current is the lost
    # charge over the run
    guide = _guide(quick, "r1_rail_ring")
    field = guide.field()
    curve = guide.curve
    assert curve is not None
    h, top = guide.solved.distance, guide.solved.null[1]
    rng = np.random.default_rng(3)
    count = 60
    x0 = curve.points[rng.integers(0, len(curve), count)] + rng.normal(0, 0.03 * h, (count, 3))
    u0 = np.zeros((count, 3))
    u0[20:40, 2], u0[40:, 2] = pt.v_of(5.0), -pt.v_of(100.0)
    detector = TwinDetector(collector=1e-12)
    watch = detector.watch(field, count)
    steps, dt = 600, 1 / (120 * guide.tune().f)
    flight = pt.integrate(x0, u0, np.full(count, pt.QM_E), dt, steps, field, [watch])
    seen = watch.reading()
    assert (seen.held, seen.lost, seen.stopped) == (20, 40, 0)
    assert np.array_equal(seen.fate == 0, flight.alive)
    assert (seen.fate[:20] == 0).all()
    assert (seen.fate[20:] == 1).all()
    lost = seen.fate == 1
    assert ((seen.step[lost] > 0) & (seen.step[lost] < steps)).all()
    assert np.allclose(seen.time[lost], (seen.step[lost] + 1) * dt, rtol=1e-12)
    radius = guide.length / (2 * math.pi)
    around = np.mod(np.arctan2(flight.x[:, 1], flight.x[:, 0]), 2 * math.pi) * radius
    miss = np.mod(seen.at - around + guide.length / 2, guide.length) - guide.length / 2
    assert np.abs(miss[lost]).max() < 1e-6
    struck = flight.x[40:, 2]  # the copper's top face is the null's height below the path, a step is 54 um
    assert (struck < -top + 1e-12).all()
    assert (struck > -top - 3 * pt.v_of(100.0) * dt).all()
    grid = guide.grid()  # the aperture's top: the last node within 4h of the null
    assert (flight.x[20:40, 2] > grid.v[-1] - grid.null[1]).all()
    assert grid.v[-1] - grid.null[1] == pytest.approx(guide.design.aperture * h, rel=0.05)
    assert np.isnan(seen.time[~lost]).all()
    assert seen.dump_charge == pytest.approx(20 * pt.E, rel=1e-12)
    assert seen.dump_volts == pytest.approx(20 * pt.E / 1e-12, rel=1e-12)
    assert seen.loss_current == pytest.approx(40 * pt.E / (steps * dt), rel=1e-12)


def test_the_detector_notes_what_another_observer_stops(quick: Path) -> None:
    # an observer before the watch (the vacuum's collisions, say) stops five of twelve held electrons at step 30: the
    # watch counts them stopped, not lost and not held, at that step's time and place
    guide = _guide(quick, "r1_rail_ring")
    field = guide.field()
    curve = guide.curve
    assert curve is not None
    count = 12
    x0 = curve.points[np.arange(count) * 150]
    taken = np.zeros(count, dtype=np.bool_)
    taken[[1, 4, 6, 7, 11]] = True

    def gas(step: int, t: float, x: pt.Arr, u: pt.Arr, alive: pt.BoolArr) -> pt.BoolArr:
        return taken & (step == 30)

    watch = TwinDetector().watch(field, count)
    dt = 1 / (120 * guide.tune().f)
    flight = pt.integrate(x0, np.zeros((count, 3)), np.full(count, pt.QM_E), dt, 60, field, [gas, watch])
    seen = watch.reading()
    assert (seen.held, seen.lost, seen.stopped) == (7, 0, 5)
    assert np.array_equal(seen.fate == 2, taken)
    assert (seen.step[taken] == 30).all()
    assert np.allclose(seen.time[taken], 31 * dt, rtol=1e-12)
    assert np.allclose(seen.at[taken], np.arange(count)[taken] * 150 * guide.length / len(curve), atol=1e-6)
    assert np.array_equal(flight.alive, ~taken)


def test_a_pickup_carries_the_image_of_a_passing_charge(quick: Path) -> None:
    # one electron along the null of a straight guide at 1 eV feels nothing and passes a pickup 10 mm long with 1 mm
    # edges, half a lap on: the image on it rises to 30 % of e and falls back, the current peaking at
    # 0.3 e v pi/(2 edge) as it enters and leaves, and nothing flows while it is beyond the edges
    guide = _guide(quick, "straight")
    field = guide.field()
    assert guide.curve is None
    v = pt.v_of(1.0)
    gamma = 1 / math.sqrt(1 - (v / pt.C0) ** 2)
    detector = TwinDetector(pickups=[0.5], length=10e-3, coupling=0.3, edge=1e-3)
    watch = detector.watch(field, 1, start=np.zeros((1, 3)))
    dt = 0.02e-3 / v  # 20 um a step
    steps = int(0.09 / (v * dt))
    pt.integrate(np.zeros((1, 3)), np.array([[0.0, 0.0, gamma * v]]), np.array([pt.QM_E]), dt, steps, field, [watch])
    seen = watch.reading()
    assert seen.held == 1
    assert seen.laps[0] == pytest.approx(steps * dt * v / 0.1, rel=1e-9)
    current = seen.currents[:, 0]
    where = v * seen.times  # the electron's place at each reading
    peak = 0.3 * pt.E * v * math.pi / (2 * 1e-3)
    assert float(current.max()) == pytest.approx(peak, rel=5e-3)
    assert float(current.min()) == pytest.approx(-peak, rel=5e-3)
    assert where[int(np.argmax(current))] == pytest.approx(0.05 - 5e-3, abs=2e-5)
    beyond = (np.abs(where - 0.05) > 5.5e-3 + 2e-5) & (np.abs(where - 0.05) < 0.045)
    assert float(np.abs(current[beyond]).max()) < 1e-12 * peak
    assert float(np.sum(current[where < 0.05]) * dt) == pytest.approx(0.3 * pt.E, rel=1e-9)


# ---- the vacuum -----------------------------------------------------------------------------------------------------
def test_the_vacuum_gives_its_pressure_free_path_and_collisions() -> None:
    # sections 6.5 and 8.2: 0.4 m at 1e-3 mbar; at 6.3e-5 mbar 6.35 m, 20 laps of T915-S, and 11 us between a 1 eV
    # electron's collisions; at 1e-6 mbar 0.67 ms. Over one free path the gas scatters 1 - 1/e of the electrons
    assert float(TwinVacuum(1e-3).free_path(1.0)) == pytest.approx(0.4, rel=1e-12)
    turbo = TwinVacuum(6.3e-5)
    assert turbo.pressure() == 6.3e-5
    assert float(turbo.free_path()) == pytest.approx(20 * 0.3178, rel=0.01)
    assert turbo.lifetime(1.0) == pytest.approx(11e-6, rel=0.03)
    assert TwinVacuum(1e-6).lifetime(1.0) == pytest.approx(0.67e-3, rel=0.03)
    assert TwinVacuum(1e-6).density() == pytest.approx(1e-4 / (pt.KB * 295.0), rel=1e-12)
    vacuum = TwinVacuum(1e-3)
    count, steps = 20000, 200
    v = pt.v_of(1.0)
    gamma = 1 / math.sqrt(1 - (v / pt.C0) ** 2)
    dt = float(vacuum.free_path(1.0)) / (v * steps)
    counts = []
    for _ in range(2):
        gas = vacuum.scatter(seed=11)
        nothing: pt.Fields = lambda x, t: (None, None)  # noqa: E731
        flight = pt.integrate(
            np.zeros((count, 3)),
            np.tile([gamma * v, 0.0, 0.0], (count, 1)),
            np.full(count, pt.QM_E),
            dt,
            steps,
            nothing,
            [gas],
        )
        assert gas.count == int((~flight.alive).sum())
        counts.append(gas.count)
    assert counts[0] == counts[1]
    sd = math.sqrt((1 - 1 / math.e) / math.e / count)
    assert counts[0] / count == pytest.approx(1 - 1 / math.e, abs=4 * sd)
    # each step's chance, exactly: integrate calls the gas after step k at t = (k + 1) dt, over a step of dt
    gas = vacuum.scatter(seed=2)
    u = np.tile([gamma * v, 0.0, 0.0], (3, 1))
    for k, step in ((0, 1e-9), (7, 2.5e-9), (199, dt)):
        chance = gas.chance(k, (k + 1) * step, u)
        assert np.allclose(chance, -math.expm1(-v * step / float(vacuum.free_path(1.0))), rtol=1e-9, atol=0)
    with pytest.raises(ValueError, match="positive pressure"):
        TwinVacuum(0.0)
    with pytest.raises(ValueError, match=r"-1\.0 K"):
        TwinVacuum(1e-6, temperature=-1.0)
    with pytest.raises(ValueError, match="'residual'"):
        TwinVacuum(1e-6, gas="argon")


# ---- the bench ------------------------------------------------------------------------------------------------------
def _strings(value: Any) -> list[str]:
    if isinstance(value, dict):
        return [s for item in value.values() for s in _strings(item)]
    if isinstance(value, list):
        return [s for item in value for s in _strings(item)]
    return [value] if isinstance(value, str) else []


def test_the_twin_bench_opens_its_devices_and_stamps_each() -> None:
    # benches/twin_r1.toml: the HackRF Pro alone as the array, R1's guide in realisation 1, a detector with two pickups
    # and the vacuum at 1e-6 mbar; a session opens what an experiment needs, nothing is solved to open or describe it,
    # and each record keeps the design's file and sha256 and no local path
    bench = Bench.named("twin_r1")
    assert bench.name == "twin-r1"
    assert {device.role for device in bench.devices} == {"array", "guide", "detector", "vacuum"}
    needs = {"guide": ("tune", "drive", "field"), "detector": ("watch",), "vacuum": ("pressure", "free_path")}
    with bench.session(needs) as devices:
        guide = devices["guide"]
        assert isinstance(guide, TwinGuide)
        assert guide.realisation == 1
        about = bench.describe(devices)
        assert "solved" not in vars(guide)
        assert devices["vacuum"].pressure() == 1e-6
        assert len(devices["detector"].pickups) == 2
    stamp = about["devices"]["guide"]["identity"]
    assert stamp["design"] == "designs/r1_rail_ring.toml"
    assert stamp["sha256"] == Design.load(DESIGNS / "r1_rail_ring.toml").digest
    assert stamp["realisation"] == 1
    text = json.dumps(about)
    assert not any(s.startswith(("/", "~")) for s in _strings(json.loads(text)))


def test_the_drivers_refuse_what_they_cannot_be() -> None:
    with pytest.raises(ValueError, match="a whole number from 0"):
        TwinGuide("designs/r1_rail_ring.toml", -1)
    with pytest.raises(ValueError, match="a whole number from 0"):
        TwinGuide("designs/r1_rail_ring.toml", True)
    guide = TwinGuide("designs/r1_rail_ring.toml")
    with pytest.raises(ValueError, match="one of them"):
        guide.drive(volts=50.0, watts=1.0)
    with pytest.raises(ValueError, match="longer than their edge"):
        TwinDetector(length=1e-3, edge=2e-3)
