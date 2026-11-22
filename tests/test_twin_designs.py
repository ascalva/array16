"""The twin's stage 2a designs and options (array16/twin; docs/KNOTTED_GUIDE.md sections 8.3 and 10): the paths H1 sets
against stage 1's (experiments/twin_rung1/REPORT.md section 6) read as their files state them, each its stage 1
counterpart with one thing changed, drawn by array16.guide.geometry and swept in the frame its file names; the bend's
shift of the null, off in every design file, put in by a design's [bend] or a sweep's bend=True, moving the swept null
by h^2 (inplane kappa_u, toward kappa_v) along the twin's own framing, and refused where it cannot be; the loss records
with the null moved, each strike on the face of the copper as the field placed it; and benches/twin_r2.toml, T915-S's
bench with its detector's entries on."""

from __future__ import annotations

import json
import math
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import array16.twin as twin
from array16 import particles as pt
from array16.guide import geometry
from array16.lab import Bench
from array16.particles import BEND_IN_PLANE, BEND_TOWARD, Arr, BendShift, InPhase
from array16.twin import BendSpec, Design, TwinDetector, TwinField, TwinGuide, section_key

ROOT = Path(__file__).resolve().parents[1]
DESIGNS = ROOT / "designs"
H915 = 0.915e-3  # T915-S's h (m)
STAGE_1 = ("r1_rail_ring", "t915s", "flat_ring_318", "helix_13")
STAGE_2A = ("ellipse_318", "flat_ring_318_stiff", "helix_13_r010", "helix_13_r030")
STEPS = 120  # steps an RF period, stage 1's

Guides = Callable[[str], TwinGuide]


@pytest.fixture(scope="module")
def guides(tmp_path_factory: pytest.TempPathFactory) -> Guides:
    """The designs as drawn (realisation 0) on their own grid, sharing one cache of solves for the module: section B is
    solved once, in about a second."""
    cache = tmp_path_factory.mktemp("cache")
    made: dict[str, TwinGuide] = {}

    def guide(name: str) -> TwinGuide:
        if name not in made:
            made[name] = TwinGuide(DESIGNS / f"{name}.toml", 0, cache=cache)
        return made[name]

    return guide


def _design(name: str) -> Design:
    return Design.load(DESIGNS / f"{name}.toml")


def _named(design: Design, like: Design) -> Design:
    """A design under another's name, line, seed, file and digest: what is left is what it draws."""
    return replace(design, name=like.name, about=like.about, seed=like.seed, source=like.source, digest=like.digest)


def _kappa(guide: TwinGuide, swept: pt.SectionField) -> tuple[Arr, Arr]:
    """The path's curvature (1/m) along the swept section's first and second normals at every sample."""
    assert guide.curve is not None
    kv = geometry.curvature_vector(guide.curve)
    out: tuple[Arr, Arr] = ((kv * swept.first).sum(1), (kv * swept.second).sum(1))
    return out


# ---- the designs ----------------------------------------------------------------------------------------------------
def test_each_stage_2a_design_is_its_stage_1_counterpart_with_one_thing_changed() -> None:
    # T915-S's section B, its 50 V through two feeds and its tolerances on every path; the ellipse is the flat ring with
    # its path changed, the stiff ring the flat ring with its frame turned from 180 to 90 degrees, and the two helices
    # the r/R = 0.2 helix at r/R = 0.1 and 0.3; no design file asks for the bend's shift of the null, and the seeds
    # their realisations draw from are all different
    flat, helix = _design("flat_ring_318"), _design("helix_13")
    ellipse, stiff = _design("ellipse_318"), _design("flat_ring_318_stiff")
    low, high = _design("helix_13_r010"), _design("helix_13_r030")
    assert replace(_named(ellipse, flat), path=flat.path) == flat
    assert replace(_named(stiff, flat), frame=flat.frame) == flat
    assert (ellipse.path.kind, ellipse.path.tightest, ellipse.path.samples) == ("ellipse", 18.56998e-3, 4000)
    assert (ellipse.frame.theta, stiff.frame.theta, flat.frame.theta) == (math.pi, math.pi / 2, math.pi)
    for helical, ratio in ((low, 0.1), (high, 0.3)):
        assert replace(_named(helical, helix), path=replace(helical.path, ratio=0.2)) == helix
        assert (helical.path.kind, helical.path.ratio, helical.path.p, helical.path.q) == ("torus", ratio, 1, 3)
        assert helical.frame.k == 0
    t915s = _design("t915s")
    every = [_design(name) for name in STAGE_1 + STAGE_2A]
    assert len({section_key(d.section) for d in every if d.section == t915s.section}) == 1
    assert sum(d.section == t915s.section for d in every) == 7
    assert all(d.tolerance == t915s.tolerance and (d.volts, d.feeds) == (50.0, 2) for d in every[1:])
    assert all(d.bend is None for d in every)
    assert len({d.seed for d in every}) == len(every)


def test_the_ellipse_keeps_t915s_inplane_radius_at_t915s_length(guides: Guides) -> None:
    # one section, band and order give one length, 317.76 mm; the ellipse's points lie on it, its tightest radius b^2/a
    # is T915-S's in its boards' plane, 18.570 mm, its widest a^2/b 119.7 mm, and its frame is the flat ring's: u
    # inward and v up, geometry.flat(curve, "in-plane"), every bend wholly in the boards' plane and toward the centre
    ellipse, t915s = guides("ellipse_318"), guides("t915s")
    assert ellipse.length == t915s.length
    assert ellipse.length == pytest.approx(0.31776, abs=1e-5)
    curve, frame = ellipse.curve, t915s.framing
    assert curve is not None
    assert frame is not None
    a, b = float(curve.points[:, 0].max()), float(curve.points[:, 1].max())
    assert b * b / a == pytest.approx(frame.inplane_radius, rel=1e-6)
    assert (a, b) == pytest.approx((64.330e-3, 34.563e-3), abs=1e-6)
    assert a * a / b == pytest.approx(119.73e-3, abs=0.01e-3)
    assert curve.bend_radius == pytest.approx(b * b / a, rel=1e-5)
    assert curve.length == pytest.approx(ellipse.length, rel=1e-14)
    x, y, z = curve.points.T
    assert float(np.abs((x / a) ** 2 + (y / b) ** 2 - 1).max()) < 1e-12
    assert not z.any()
    assert (x[0], y[0]) == pytest.approx((a, 0.0), abs=1e-15)  # the start: the end of the major axis
    swept = ellipse.sweep(ellipse.grid(), InPhase(1.0, 0.0))
    flat = geometry.flat(curve, "in-plane")
    assert float(np.abs(swept.first - flat.across).max()) < 1e-12
    assert float(np.abs(swept.second + flat.normal).max()) < 1e-12
    across, toward = _kappa(ellipse, swept)
    assert float(np.abs(toward).max()) < 1e-9 * float(across.max())
    assert float(across.min()) > 0
    assert 1 / float(across.max()) == pytest.approx(b * b / a, rel=1e-5)


def test_the_stiff_ring_is_the_flat_ring_bent_toward_its_boards(guides: Guides) -> None:
    # the same 50.57 mm circle; the flat ring's u runs inward with the whole bend across its boards' plane, the stiff
    # ring's down the ring's axis with v inward and the whole bend toward its boards, geometry.flat(curve, "toward")
    stiff, flat = guides("flat_ring_318_stiff"), guides("flat_ring_318")
    assert stiff.curve is not None
    assert flat.curve is not None
    assert np.array_equal(stiff.curve.points, flat.curve.points)
    radius = stiff.length / (2 * math.pi)
    assert radius == pytest.approx(50.573e-3, abs=1e-6)
    bent, plain = stiff.sweep(stiff.grid(), InPhase(1.0, 0.0)), flat.sweep(flat.grid(), InPhase(1.0, 0.0))
    toward = geometry.flat(stiff.curve, "toward")
    assert float(np.abs(bent.first - toward.across).max()) < 1e-12
    assert float(np.abs(bent.second + toward.normal).max()) < 1e-12
    assert float(np.abs(np.abs(bent.first[:, 2]) - 1).max()) < 1e-12  # u along the ring's axis
    assert toward.inplane_radius == math.inf
    assert toward.toward_radius == pytest.approx(radius, rel=1e-6)
    for swept, bend_u, bend_v in ((bent, 0.0, 1 / radius), (plain, 1 / radius, 0.0)):
        ku, kv = _kappa(stiff, swept)
        assert ku == pytest.approx(np.full(len(ku), bend_u), abs=1e-6 / radius)
        assert kv == pytest.approx(np.full(len(kv), bend_v), abs=1e-6 / radius)


@pytest.mark.parametrize(
    ("name", "ratio", "tightest", "inplane"),
    [
        ("helix_13_r010", 0.1, 31.47e-3, 60.08e-3),
        ("helix_13", 0.2, 25.05e-3, 29.14e-3),
        ("helix_13_r030", 0.3, 20.15e-3, 21.06e-3),
    ],
)
def test_the_helices_are_the_13_torus_curve_framed_k_0(
    guides: Guides, name: str, ratio: float, tightest: float, inplane: float
) -> None:
    # the (1, 3) torus curve at r/R = 0.1, 0.2 and 0.3 of T915-S's length, geometry.helix's to the bit: the tighter
    # the winding, the tighter its bend, and framed k = 0 the tighter its least radius in the boards' plane, each above
    # T915-S's 18.57 mm
    guide = guides(name)
    assert guide.curve is not None
    assert guide.framing is not None
    assert np.array_equal(guide.curve.points, geometry.helix(ratio, guide.length, 3000).points)
    assert guide.curve.bend_radius == pytest.approx(tightest, abs=0.01e-3)
    assert guide.framing.k == 0
    assert guide.framing.inplane_radius == pytest.approx(inplane, abs=0.01e-3)
    assert guide.framing.inplane_radius > 18.57e-3


def test_the_twin_frames_with_the_geometry_modules_framing() -> None:
    # one framing for the twin and the study's ports: array16.twin's names are array16.guide.geometry's
    assert twin.framing is geometry.framing
    assert twin.Framing is geometry.Framing
    assert twin.curvature_vector is geometry.curvature_vector


# ---- the bend's shift of the null -----------------------------------------------------------------------------------
def test_the_bend_shift_is_off_in_every_design_and_on_where_a_design_or_a_sweep_asks(guides: Guides) -> None:
    # no design file asks for it, so every guide sweeps its null on the path and records what stage 1 recorded; a
    # design's [bend] shift = true or a sweep's bend=True moves it, with T915-S's h and the facing boards' 2 ln 2/pi^2
    # and 1/pi (or the design's own), and bend=False takes it out again
    for name in STAGE_1 + STAGE_2A:
        guide = TwinGuide(DESIGNS / f"{name}.toml")
        assert guide.bend_shift() is None
        assert set(guide.identity()) == {"design", "name", "sha256", "realisation", "section", "solver"}
    t915s = guides("t915s")
    assert t915s.sweep(t915s.grid(), InPhase(1.0, 0.0)).bend is None
    assert t915s.bend_shift(True) == BendShift(H915, BEND_IN_PLANE, BEND_TOWARD)
    assert t915s.sweep(t915s.grid(), InPhase(1.0, 0.0), bend=True).bend == BendShift(H915)
    text = (DESIGNS / "t915s.toml").read_text()
    stated = Design.parse(text + "\n[bend]\nshift = true\n", "t915s_bent.toml")
    assert stated.bend == BendSpec(True, BEND_IN_PLANE, BEND_TOWARD)
    own = Design.parse(text + "\n[bend]\nshift = false\ninplane = 0.2\ntoward = 0.3\n", "t915s_bent.toml")
    assert own.bend == BendSpec(False, 0.2, 0.3)
    assert replace(own, bend=None) == replace(stated, bend=None)


def test_the_bent_null_moves_by_its_first_order_shift_along_the_twins_framing(guides: Guides) -> None:
    # T915-S framed k = 0 with bend=True: at every sample the swept field is the field with the null on the path moved
    # rigidly by delta = h^2 (inplane kappa_g, -toward kappa_n) on the twin's framing (u across the boards' plane, v
    # along minus the board normal), and at the tightest bend, 15.29 mm, the shift lies between that bend's 7.7 um in
    # the boards' plane and its 17.4 um toward them
    guide = guides("t915s")
    frame = guide.framing
    assert frame is not None
    grid = guide.grid()
    bent = guide.sweep(grid, InPhase(1.0, 0.0), bend=True)
    flat = guide.sweep(grid, InPhase(1.0, 0.0))
    samples = np.arange(0, 4000, 250)
    du, dv = bent.null_shift(bent.locate(bent.points[samples]))
    h2 = H915 * H915
    assert du == pytest.approx(BEND_IN_PLANE * h2 * frame.kappa_g[samples], rel=1e-9, abs=1e-15)
    assert dv == pytest.approx(-BEND_TOWARD * h2 * frame.kappa_n[samples], rel=1e-9, abs=1e-15)
    rng = np.random.default_rng(7)
    for k, i in enumerate(samples):
        offsets = rng.uniform(-0.5, 0.5, size=(6, 2)) * H915
        delta = du[k] * bent.first[i] + dv[k] * bent.second[i]
        x = bent.points[i] + offsets[:, :1] * bent.first[i] + offsets[:, 1:] * bent.second[i]
        e_bent, _ = bent(x, 0.0)
        e_flat, _ = flat(x - delta, 0.0)
        assert e_bent is not None
        assert e_flat is not None
        assert float(np.abs(e_bent - e_flat).max()) < 1e-9 * float(np.abs(e_flat).max())
    curve = guide.curve
    assert curve is not None
    tightest = int(np.argmax(curve.curvature))
    shift = math.hypot(*(float(c[0]) for c in bent.null_shift(bent.locate(bent.points[tightest : tightest + 1]))))
    assert 7.7e-6 * 0.99 < shift < 17.4e-6 * 1.01


def test_a_design_says_what_the_bend_and_the_ellipse_take_when_refused() -> None:
    t915s, r1, flat = ((DESIGNS / f"{name}.toml").read_text() for name in ("t915s", "r1_rail_ring", "flat_ring_318"))
    straight = t915s.replace('kind = "torus"', 'kind = "straight"\nlength = 0.3').replace("[frame]\nk = 0", "")
    for change, message in (
        (t915s + "\n[bend]\nshift = 'yes'\n", "takes shift = true or false"),
        (t915s + "\n[bend]\nshift = true\ncolour = 1\n", r"\[bend\]: unknown keys \['colour'\]"),
        (t915s + "\n[bend]\nshift = true\ninplane = nan\n", "finite coefficients"),
        (r1 + "\n[bend]\nshift = true\n", "the rails section's bend correction takes its own coefficients"),
        (r1 + "\n[bend]\nshift = true\ninplane = 0.1\n", "takes its own coefficients"),
        (straight + "\n[bend]\nshift = true\n", "a straight path has no bend"),
        (flat.replace('kind = "circle"', 'kind = "ellipse"'), "an ellipse takes \\[path\\] tightest"),
        (flat.replace("samples = 4000", "samples = 4000\ntightest = 0.02"), "no other path does"),
        (
            flat.replace('kind = "circle"', 'kind = "ellipse"\nlength = 0.3\ntightest = 0.05'),
            "under that of the circle of its length",
        ),
        (flat.replace('kind = "circle"', 'kind = "spiral"'), "kind is circle, ellipse, torus or straight"),
    ):
        with pytest.raises(ValueError, match=message):
            Design.parse(change, "refused.toml")
    assert Design.parse(r1 + "\n[bend]\nshift = true\ninplane = 0.1\ntoward = 0.2\n").bend == BendSpec(True, 0.1, 0.2)
    rails = TwinGuide(DESIGNS / "r1_rail_ring.toml")
    with pytest.raises(ValueError, match="Rails section's bend correction takes its own coefficients"):
        rails.bend_shift(True)
    with pytest.raises(ValueError, match="bend true or false"):
        rails.field(bend="yes")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="bend true or false"):
        rails.bend_shift(1)  # type: ignore[arg-type]
    assert "solved" not in vars(rails)  # refused before anything is solved


# ---- the loss records with the null moved ---------------------------------------------------------------------------
def _momentum(energy: float, headings: Arr) -> Arr:
    speed = pt.v_of(energy)
    out: Arr = speed / math.sqrt(1 - (speed / pt.C0) ** 2) * headings
    return out


def test_with_the_null_moved_a_strike_is_recorded_on_the_copper_as_the_field_placed_it(guides: Guides) -> None:
    # the stiff ring with no volts, so that every electron flies straight: two electrons sent from 1 mm off the split,
    # outward and inward, strike the outer and the inner board's copper at v = 0 and v = 2h in the grid's coordinates
    # with the bend in or out; with it in, the boards stand h^2/(pi R) = 5.27 um nearer the centre in space, where the
    # field moved them, and the entries say so: 'struck', on the face, the step before the watch stopped them. (Sent
    # down the split itself they strike the polyimide behind it, 18 um further, at v = -18 um and 2h + 18 um.)
    guide = guides("flat_ring_318_stiff")
    grid = guide.grid()
    dt = 1 / (STEPS * guide.tune().f)
    radius = guide.length / (2 * math.pi)
    for bend, moved in ((False, 0.0), (True, BEND_TOWARD * H915 * H915 / radius)):
        swept = guide.sweep(grid, InPhase(0.0, 0.0), bend=bend)
        field = TwinField(swept, guide.curve)
        x0 = swept.points[[0, 1000]] + 1e-3 * swept.first[[0, 1000]]  # 1 mm off the split, under the copper
        u0 = _momentum(30.0, np.stack([-swept.second[0], swept.second[1000]]))
        watch = TwinDetector(losses=True).watch(field, 2, start=x0)
        pt.integrate(x0, u0, np.full(2, pt.QM_E), dt, 60, field, [watch])
        seen = watch.reading()
        losses = seen.losses
        assert losses is not None
        assert seen.lost == 2
        assert list(losses.how) == ["struck", "struck"]
        assert np.array_equal(losses.step, seen.step[losses.index] - 1)
        order = np.argsort(losses.index)
        faces = np.array([0.0, 2 * H915])
        assert losses.v[order] == pytest.approx(faces, abs=1e-11)
        where = swept.locate(losses.point[order])
        assert where.v == pytest.approx(faces + moved, abs=1e-11)
        assert where.u == pytest.approx(np.full(2, 1e-3), abs=1e-12)
        assert losses.u == pytest.approx(np.full(2, 1e-3), abs=1e-12)  # the bend toward the boards moves no u


# ---- the stage 2 bench ----------------------------------------------------------------------------------------------
def _strings(value: Any) -> list[str]:
    if isinstance(value, dict):
        return [s for item in value.values() for s in _strings(item)]
    if isinstance(value, list):
        return [s for item in value for s in _strings(item)]
    return [value] if isinstance(value, str) else []


def test_the_stage_2_bench_is_t915s_built_once_with_its_detectors_entries_on() -> None:
    # benches/twin_r2.toml: the HackRF Pro alone as the array, T915-S's guide in realisation 1, a detector with two
    # pickups that writes an entry for each electron lost, and the vacuum at 1e-6 mbar; nothing is solved to open or
    # describe it, its stamp keeps the design's file and sha256, the detector's entries, and no local path, and its
    # workers are this machine's cores but two
    bench = Bench.named("twin_r2")
    assert bench.name == "twin-r2"
    assert bench.workers == Bench.local().workers
    assert {device.role for device in bench.devices} == {"array", "guide", "detector", "vacuum"}
    needs = {"guide": ("tune", "drive", "field", "identity"), "detector": ("watch", "identity"), "vacuum": ("scatter",)}
    with bench.session(needs) as devices:
        guide, detector = devices["guide"], devices["detector"]
        assert isinstance(guide, TwinGuide)
        assert isinstance(detector, TwinDetector)
        assert (guide.design.name, guide.realisation, guide.design.bend) == ("T915-S", 1, None)
        assert detector.losses is True
        assert len(detector.pickups) == 2
        assert devices["vacuum"].pressure() == 1e-6
        about = bench.describe(devices)
        assert "solved" not in vars(guide)
    stamp = about["devices"]
    assert stamp["guide"]["identity"]["design"] == "designs/t915s.toml"
    assert stamp["guide"]["identity"]["sha256"] == _design("t915s").digest
    assert stamp["detector"]["identity"]["losses"] is True
    assert stamp["detector"]["options"]["losses"] is True
    assert not any(s.startswith(("/", "~")) for s in _strings(json.loads(json.dumps(about))))
