"""The twin bench: the guide, its detector and its vacuum, simulated, as devices of an array16.lab bench
(docs/KNOTTED_GUIDE.md section 8.3). A bench file names them as it names hardware, and hardware takes their place at
the build stage as a change of bench file:

    [device.guide]
    driver = "array16.twin:TwinGuide"
    design = "designs/r1_rail_ring.toml"
    realisation = 1

  Design        a design file (designs/*.toml): the band and the ring's order, the cross-section, the path and its
                frame, the drive, and the tolerances a realisation draws
  TwinGuide     one realisation of a design: tune (the ring's +-m pair, as a sweep finds it), drive (compiled current
                feeds: the volts round the ring, the power, the backward wave, the well) and field() (a Fields for
                integrate), from the solved cross-section (array16.guide.sections, each solve cached under .twin_cache/
                by the hash of its content and the solver's source), the ring of segments (array16.guide.ring) and the
                section swept along the path (array16.particles.sections)
  TwinDetector  observers for integrate: the particles held and lost, where and when each was lost, the currents the
                pickups carry, and the charge the held particles bring to the dump; and, asked for (`losses`), an
                entry for each electron lost: where along the path and in the section it left, when, with what
                energy and sense, and how (array16.twin.records)
  TwinVacuum    the pressure, an electron's mean free path and lifetime between collisions, and an observer that stops
                the electrons the gas scatters

A realisation is the design built once: realisation 0 is the design as drawn, and realisation k > 0 draws, from the
design's seed and k, what a build departs by and the predictions do not know: the spacers' gaps, the closing joints'
impedance and place, the split's etch and the laminate's permittivity (array16.guide.ring.ring_draw), and Earth's field
in an orientation drawn uniformly. The calibration finds them blind, as on V1's simulated array.

The model. The section's field is quasi-static and swept along the path, its null where the 2-D solve puts it; the
bend's shift of the null and the field along the path come with the 3-D solve array16/guide/bem.py will hold (section
4.5). Its magnetic field is T x E / v at the wave's phase velocity from the loaded section's E, where a line on a
laminate carries the air line's pattern, (n/c) T x E_air: the magnetic force is v_e / v_phase of the electric, 1 % at
9 eV on R1, and the Doppler shift and the magnetic term cancel in the pseudopotential to that order; RingWave gives
the backward part of the wave, |b/a| under 1e-3 on the designs, the forward wave's slowness. The ring is one line mode,
both of R1's rails driven alike: full wave finds R1's two rails resonating apart (experiments/twin_openems_line), and
the two-rail line matrices, with the bridges between the rails, will join array16.guide.ring for the drive R1 is built
with. The gas scatters with one total cross-section at every energy; LXCat's cross-sections by energy will take its
place for a named gas. SI units, pressures in mbar.
"""

from __future__ import annotations

import hashlib
import json
import math
import operator
import os
import tomllib
import zipfile
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, fields
from functools import cached_property
from pathlib import Path
from typing import Any, Literal

import numpy as np

from array16.guide import sections as _sections
from array16.guide.ring import (
    RING_NODES,
    RING_SENSITIVITY,
    SCARF,
    SPACERS,
    Build,
    compile_feeds,
    feed_layout,
    perfect_ring,
    ring_draw,
    ring_resonances,
)
from array16.guide.sections import (
    COARSE,
    FR4,
    POLYIMIDE,
    RO4003C,
    FacingStrips,
    Grid,
    Laminate,
    Rails,
    Section,
    SlotGuide,
    Solved,
    solve,
)
from array16.lab import ROOT, plain
from array16.particles import (
    C0,
    KB,
    ME,
    Arr,
    BoolArr,
    CArr,
    ClosedCurve,
    E,
    FacingBoards,
    IntArr,
    RailPair,
    SectionField,
    SectionGrid,
    mathieu_q,
    ponderomotive_ev,
)
from array16.particles.sections import PathDrive

from .records import Losses, LossLog

__all__ = [
    "CACHE",
    "Collisions",
    "Design",
    "Detection",
    "Driven",
    "FrameSpec",
    "Framing",
    "Losses",
    "PathSpec",
    "Pickup",
    "Realisation",
    "RingWave",
    "Tolerance",
    "Tuning",
    "TwinDetector",
    "TwinField",
    "TwinGuide",
    "TwinVacuum",
    "Watch",
    "curvature_vector",
    "framing",
    "section_key",
    "solve_cached",
    "solver_version",
]

CACHE = ROOT / ".twin_cache"  # solved sections by content hash; git ignores it
FORMAT = 1  # the cache files' layout
LAMINATES = {laminate.name: laminate for laminate in (FR4, RO4003C, POLYIMIDE)}
GRIDS = {"default": Grid(), "coarse": COARSE}
GEOMETRIES: dict[str, type[FacingStrips] | type[SlotGuide] | type[Rails]] = {
    "facing": FacingStrips,
    "slot": SlotGuide,
    "rails": Rails,
}
GAPLESS_STEP = 1 / 40  # of h: the node spacing of a closed-form section sampled for the swept field
FREE_PATH = (0.4, 1e-3, 295.0)  # an electron's mean free path (m) at 1e-3 mbar and 295 K: the estimate of sections
# 6.5 and 8.2, 11 us between collisions at 1 eV and 6.3e-5 mbar, which sets TwinVacuum's cross-section


# ---- the design file ------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class PathSpec:
    """The path the section is swept along: a circle in the xy plane, a (p, q) torus knot on a torus of radii r/R =
    ratio (the trefoil is (2, 3)), or a straight line repeating every lap (a ring unrolled); its length (m), or None for
    m guided waves at the band on the solved section's index, the ring resonant there; and its samples."""

    kind: Literal["circle", "torus", "straight"]
    length: float | None
    samples: int
    ratio: float = 0.0
    p: int = 2
    q: int = 3


@dataclass(frozen=True)
class FrameSpec:
    """How the section turns along the path. With k None, the curve's twist-free axes turned by theta (rad): the first
    normal, u's direction, is cos theta e1 + sin theta e2. With k, the board normal closes after a lap with k half
    turns against the twist-free frame, turning (k pi - holonomy)/L a metre from theta0 (rad; None: the start that
    keeps the largest radius of curvature in the boards' plane, the knotted study's rule), and u runs across the
    boards' plane."""

    theta: float = 0.0
    k: int | None = None
    theta0: float | None = None


@dataclass(frozen=True)
class Tolerance:
    """What a realisation draws (array16.guide.ring.ring_draw), each zero in the design as drawn: the spacers' rms
    relative gap error at `stations` stations, ln Z0 moving `sensitivity` times dh/h; the two closing joints' step in
    ln Z, its sd and the sd of their spacing (m) over a scarf (m); the split's etch sd at `split_stations`, ln Z moving
    `split_sensitivity` times it; the laminate permittivity's sd; and Earth's field (T), in an orientation drawn
    uniformly."""

    spacers: float = 0.0
    stations: int = SPACERS
    sensitivity: float = RING_SENSITIVITY
    joint_dz: float = 0.0
    joint_sd: float = 0.0
    scarf: float = SCARF
    jitter: float = 0.5e-3
    split_sd: float = 0.0
    split_sensitivity: float = 0.03
    split_stations: int = 6
    eps: float = 0.0
    earth: float = 0.0


_TOP = {"name", "about", "seed", "band", "section", "path", "frame", "drive", "tolerance"}
_TABLES = {
    "band": {"f", "m"},
    "path": {"kind", "length", "samples", "ratio", "p", "q"},
    "frame": {"theta_deg", "k", "theta0_deg"},
    "drive": {"volts", "q", "feeds"},
    "tolerance": {f.name for f in fields(Tolerance)},
}
_SECTION = {"kind", "laminate", "er", "tand", "thickness", "copper", "resistance", "rough", "grid", "aperture"}


def _known(where: str, given: Mapping[str, Any], allowed: set[str]) -> None:
    extra = sorted(set(given) - allowed)
    if extra:
        raise ValueError(f"{where}: unknown keys {extra}; it takes {sorted(allowed)}")


def _number(value: Any) -> Any:
    """TOML's integers as floats where a float is meant, so that 100 and 100.0 hash alike."""
    return float(value) if isinstance(value, int) and not isinstance(value, bool) else value


def _section(spec: Mapping[str, Any], f: float, where: str) -> tuple[Section, float]:
    kind = spec.get("kind")
    if kind not in GEOMETRIES:
        raise ValueError(f"{where}: [section] kind is one of {sorted(GEOMETRIES)}, not {kind!r}")
    cls = GEOMETRIES[kind]
    shape = {field.name for field in fields(cls)}
    _known(f"{where} [section] ({kind})", spec, _SECTION | shape)
    geometry = cls(**{key: _number(value) for key, value in spec.items() if key in shape})
    name = spec.get("laminate", "none")
    laminate: Laminate | None
    if name == "none":
        laminate = None
    elif name in LAMINATES:
        base = LAMINATES[name]
        laminate = Laminate(base.name, float(spec.get("er", base.er)), float(spec.get("tand", base.tand)))
    else:
        raise ValueError(f"{where}: [section] laminate is one of {sorted(LAMINATES)} or 'none', not {name!r}")
    grid = spec.get("grid", "default")
    if grid not in GRIDS:
        raise ValueError(f"{where}: [section] grid is one of {sorted(GRIDS)}, not {grid!r}")
    section = Section(
        geometry,
        laminate,
        float(spec.get("thickness", 0.0)),
        float(spec.get("copper", _sections.FLEX_COPPER)),
        f,
        float(spec.get("resistance", 1.0)),
        GRIDS[grid],
        float(spec.get("rough", 1.0)),
    )
    aperture = float(spec.get("aperture", 4.0))
    if not aperture > 0:
        raise ValueError(f"{where}: the aperture is a positive number of null distances, not {aperture}")
    return section, aperture


@dataclass(frozen=True)
class Design:
    """A design as its file states it: name, a line about it and the seed its realisations are drawn from; the band
    f (Hz) and the ring's order m; the cross-section and the aperture (the half-width round the null, in null
    distances h, that the swept field covers: a particle beyond it has left the guide); the path and its frame; the
    drive, as peak volts (across the resonator, as array16.guide.sections counts them) or as the null's Mathieu q,
    through `feeds` compiled current feeds; and the tolerances. `source` is where it was read from, relative to the
    repository, and `digest` its file's sha256."""

    name: str
    about: str
    seed: int
    f: float
    m: int
    section: Section
    aperture: float
    path: PathSpec
    frame: FrameSpec
    volts: float | None
    q: float | None
    feeds: int
    tolerance: Tolerance
    source: str = ""
    digest: str = ""

    @classmethod
    def load(cls, path: Path) -> Design:
        """A design file. Unknown keys are refused, with the keys each table takes."""
        data = path.read_bytes()
        try:
            source = path.resolve().relative_to(ROOT).as_posix()
        except ValueError:
            source = path.name
        return cls.parse(data.decode(), source, hashlib.sha256(data).hexdigest())

    @classmethod
    def parse(cls, text: str, source: str = "", digest: str = "") -> Design:
        spec = tomllib.loads(text)
        where = source or "a design"
        _known(where, spec, _TOP)
        for table, allowed in _TABLES.items():
            _known(f"{where} [{table}]", spec.get(table, {}), allowed)
        band = spec.get("band", {})
        f, m = float(band.get("f", 0.0)), int(band.get("m", 1))
        if not (f > 0 and m >= 1):
            raise ValueError(f"{where}: [band] takes a frequency f > 0 (Hz) and an order m >= 1, not {f} and {m}")
        section, aperture = _section(spec.get("section", {}), f, where)
        path = spec.get("path", {})
        kind = path.get("kind", "circle")
        if kind not in ("circle", "torus", "straight"):
            raise ValueError(f"{where}: [path] kind is circle, torus or straight, not {kind!r}")
        length = path.get("length")
        route = PathSpec(
            kind,
            None if length is None else float(length),
            int(path.get("samples", 2000)),
            float(path.get("ratio", 0.0)),
            int(path.get("p", 2)),
            int(path.get("q", 3)),
        )
        if route.samples < 3 or (route.length is not None and not route.length > 0):
            raise ValueError(f"{where}: a path takes 3 samples or more and a positive length")
        if kind == "torus" and not 0 < route.ratio < 1:
            raise ValueError(f"{where}: a torus knot takes 0 < ratio = r/R < 1, not {route.ratio}")
        frame = spec.get("frame", {})
        k = frame.get("k")
        if k is not None and kind == "straight":
            raise ValueError(f"{where}: a straight path has no holonomy to close: give theta_deg, not k")
        theta0 = frame.get("theta0_deg")
        turn = FrameSpec(
            math.radians(float(frame.get("theta_deg", 0.0))),
            None if k is None else int(k),
            None if theta0 is None else math.radians(float(theta0)),
        )
        drive = spec.get("drive", {})
        volts, q = drive.get("volts"), drive.get("q")
        if (volts is None) == (q is None):
            raise ValueError(f"{where}: [drive] gives the volts or the null's q, one of them")
        feeds = int(drive.get("feeds", 2))
        if feeds < 1:
            raise ValueError(f"{where}: a drive takes one feed or more, not {feeds}")
        counts = {"stations", "split_stations"}
        tolerance = Tolerance(
            **{key: int(value) if key in counts else _number(value) for key, value in spec.get("tolerance", {}).items()}
        )
        return cls(
            name=str(spec.get("name", source)),
            about=str(spec.get("about", "")),
            seed=int(spec.get("seed", 0)),
            f=f,
            m=m,
            section=section,
            aperture=aperture,
            path=route,
            frame=turn,
            volts=None if volts is None else float(volts),
            q=None if q is None else float(q),
            feeds=feeds,
            tolerance=tolerance,
            source=source,
            digest=digest,
        )


# ---- solves cached by content hash ----------------------------------------------------------------------------------
def solver_version() -> str:
    """The sha256 of array16/guide/sections.py's source, the first 16 digits: a solve is reused only by the solver
    that made it."""
    return hashlib.sha256(Path(_sections.__file__).read_bytes()).hexdigest()[:16]


def section_key(section: Section) -> str:
    """The content hash of a section's solve: its geometry's kind and every number of the section, the grid's
    included, with the solver's version and the cache's format."""
    spec = {
        "format": FORMAT,
        "solver": solver_version(),
        "geometry": type(section.geometry).__name__,
        "section": plain(section),
    }
    return hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()


_SCALARS = ("capacitance", "air_capacitance", "share", "qc", "qd", "gradient", "escape", "distance")


def _write(path: Path, key: str, solved: Solved) -> None:
    meta = {name: getattr(solved, name) for name in _SCALARS}
    meta |= {"key": key, "null": list(solved.null), "saddle": list(solved.saddle)}
    path.parent.mkdir(parents=True, exist_ok=True)
    scratch = path.with_name(f".{path.stem}.{os.getpid()}.tmp")
    with scratch.open("wb") as out:
        np.savez_compressed(out, u=solved.u, v=solved.v, phi=solved.phi, meta=np.array(json.dumps(meta)))
    os.replace(scratch, path)  # whole or not at all, whichever process finishes first


def _read(path: Path, key: str, section: Section) -> Solved:
    with np.load(path, allow_pickle=False) as data:
        meta = json.loads(str(data["meta"]))
        if meta["key"] != key:
            raise ValueError(f"{path.name} holds the solve {meta['key']}, not {key}")
        arrays = [np.array(data[name]) for name in ("u", "v", "phi")]
    for array in arrays:
        array.setflags(write=False)
    scalars = {name: float(meta[name]) for name in _SCALARS}
    null, saddle = meta["null"], meta["saddle"]
    return Solved(
        section=section,
        null=(float(null[0]), float(null[1])),
        saddle=(float(saddle[0]), float(saddle[1])),
        u=arrays[0],
        v=arrays[1],
        phi=arrays[2],
        **scalars,
    )


def solve_cached(section: Section, cache: Path | None = CACHE) -> Solved:
    """array16.guide.sections.solve, kept under `cache` by section_key: a section solved once, by this solver, is read
    back rather than solved again (None: always solved). A file that cannot be read is solved and written anew."""
    if cache is None:
        return solve(section)
    key = section_key(section)
    path = cache / f"{key}.npz"
    if path.exists():
        try:
            return _read(path, key, section)
        except (OSError, EOFError, ValueError, KeyError, zipfile.BadZipFile):
            pass
    solved = solve(section)
    _write(path, key, solved)
    return solved


# ---- frames: the knotted study's framing ----------------------------------------------------------------------------
def curvature_vector(curve: ClosedCurve) -> Arr:
    """dT/ds at every point (n, 3), by second differences, without its part along the tangent (the knotted study's)."""
    pts, ds = curve.points, curve.length / len(curve)
    k = (np.roll(pts, -1, 0) - 2 * pts + np.roll(pts, 1, 0)) / ds**2
    t = curve.tangents
    out: Arr = k - (k * t).sum(1)[:, None] * t
    return out


@dataclass(frozen=True, eq=False)
class Framing:
    """Board normals along a curve (experiments/knotted_guide_study.py's Framing): its twist-free axes turned by
    theta0 + twist s, so that after a lap they come back turned by k pi. normal is the board normal n_b, across the
    direction T x n_b in the boards' plane; kappa_g and kappa_n are the curvature in the boards' plane and toward them
    (1/m)."""

    k: int
    theta0: float
    twist: float
    normal: Arr
    across: Arr
    kappa_g: Arr
    kappa_n: Arr

    @property
    def inplane_radius(self) -> float:
        """The smallest radius of curvature in the boards' plane (m)."""
        return float(1 / np.abs(self.kappa_g).max())


def framing(curve: ClosedCurve, k: int, theta0: float | None = None, starts: int = 720) -> Framing:
    """The framing that closes with k half turns against the twist-free frame, (k pi - holonomy) / L a metre; without
    theta0, the start that keeps the largest in-plane radius (section 3.5). The knotted study's, the same arithmetic."""
    n, length = len(curve), curve.length
    e1, e2 = curve.axes
    kv = curvature_vector(curve)
    a, b = (kv * e1).sum(1), (kv * e2).sum(1)
    twist = (k * math.pi - curve.holonomy) / length
    s = np.arange(n) * length / n
    if theta0 is None:
        trial = np.linspace(0, math.pi, starts, endpoint=False)
        th = trial[:, None] + twist * s[None, :]
        kg = np.cos(th) * b[None, :] - np.sin(th) * a[None, :]  # kappa . (T x n_b), n_b = cos th e1 + sin th e2
        theta0 = float(trial[int(np.argmin(np.abs(kg).max(1)))])
    th = theta0 + twist * s
    normal = np.cos(th)[:, None] * e1 + np.sin(th)[:, None] * e2
    across = np.cross(curve.tangents, normal)
    return Framing(k, theta0, twist, normal, across, (kv * across).sum(1), (kv * normal).sum(1))


# ---- the ring's volts along the path --------------------------------------------------------------------------------
@dataclass(frozen=True, eq=False)
class RingWave:
    """The volts of a realised ring along its path, a PathDrive: V(s, t) = Re(V(s) e^(i omega t)), the ring model's
    phasors, with V(s) the node volts of its ring of segments continued between the nodes by their Fourier series
    (a table `oversample` times finer than the nodes, then linearly); and the magnetic field of the travelling wave of
    order m the feeds were compiled for, T x E / v. On a perfect ring it is TravellingWave(|V|, omega, m) to the
    table's interpolation, (2 pi m / table)^2 / 8 of the volts."""

    table: CArr
    omega: float
    m: int

    @classmethod
    def of(cls, volts: CArr, omega: float, m: int, oversample: int = 16) -> RingWave:
        """From node volts (N,) equally spaced round the lap, node k at s = k L / N."""
        count = len(volts)
        spectrum = np.fft.fft(np.asarray(volts, dtype=np.complex128))
        fine = count * oversample
        padded = np.zeros(fine, dtype=np.complex128)
        half = (count + 1) // 2
        padded[:half] = spectrum[:half]
        padded[fine - (count - half) :] = spectrum[half:]
        table: CArr = np.fft.ifft(padded) * oversample
        return cls(table, omega, m)

    def phasor(self, angle: Arr) -> CArr:
        """V(s) at angles 2 pi s / L (complex volts)."""
        size = len(self.table)
        pos = np.mod(np.asarray(angle, dtype=np.float64) * (size / (2 * math.pi)), size)
        pos = np.where(np.isfinite(pos), pos, 0.0)
        i = np.minimum(pos.astype(np.intp), size - 1)
        frac = pos - i
        out: CArr = (1 - frac) * self.table[i] + frac * self.table[(i + 1) % size]
        return out

    def __call__(self, angle: Arr, t: float) -> Arr:
        out: Arr = np.real(self.phasor(angle) * np.exp(1j * self.omega * t))
        return out

    def slowness(self, length: float) -> float | None:
        """The forward wave's 1/v: the backward part, |b/a| of it, takes the same sign here, not the opposite."""
        return 2 * math.pi * self.m / (self.omega * length)


@dataclass(frozen=True, eq=False)
class TwinField:
    """A realisation's field as integrate takes it (a Fields): the swept section, and Earth's field (T, uniform) where
    the realisation draws it. `curve` is the path (None on a straight one); the swept field's flags and its observer
    pass through."""

    swept: SectionField
    curve: ClosedCurve | None
    earth: Arr | None = None

    def __call__(self, x: Arr, t: float) -> tuple[Arr | None, Arr | None]:
        e, b = self.swept(x, t)
        if self.earth is None:
            return e, b
        background: Arr = np.tile(self.earth, (len(x), 1))
        return e, background if b is None else b + background

    @property
    def length(self) -> float:
        return self.swept.length

    @property
    def off_grid(self) -> BoolArr:
        """The particles the last call found off the section's grid or in its copper or boards."""
        return self.swept.off_grid

    def stop_off_grid(self, step: int, t: float, x: Arr, u: Arr, alive: BoolArr) -> BoolArr:
        """An Observer for integrate: stops what the last call found off the grid (SectionField.stop_off_grid)."""
        return self.swept.stop_off_grid(step, t, x, u, alive)


# ---- the guide ------------------------------------------------------------------------------------------------------
@dataclass(frozen=True, eq=False)
class Realisation:
    """What one realisation of a design draws: its number; the ring's impedance (ohm), index and gap error dh/h at its
    nodes; and Earth's field (T, (3,): zero when the design draws none)."""

    number: int
    z: Arr
    n: Arr
    rel: Arr
    earth: Arr


@dataclass(frozen=True)
class Tuning:
    """A realisation's resonance as a frequency sweep finds it: the +-m pair's centre (Hz) and its Q, both roots
    (complex Hz, f (1 + i / 2Q)), the splitting delta f / f and in linewidths; the band it was designed for (Hz), the
    order, the ring's length (m), its mean index, its line's Z0 (ohm) and the ring's shunt resistance Z0 Q / (2 pi m)
    (ohm)."""

    f: float
    q: float
    lower: complex
    upper: complex
    splitting: float
    linewidths: float
    band: float
    m: int
    length: float
    index: float
    z0: float
    shunt: float


@dataclass(frozen=True, eq=False)
class Driven:
    """A realisation driven through compiled current feeds at f (Hz): the power it takes (W); the peak volts round the
    ring, their mean, least and most (V); the backward wave |b/a|; the weakest well over the strongest (DR-6's
    min/max); the weakest well's depth (eV) and the null's Mathieu q at the mean volts; the feeds' nodes and currents
    (A, complex); and the node volts (V, complex phasors, e^(i omega t)), the pattern field() sweeps."""

    f: float
    watts: float
    volts: float
    least: float
    most: float
    purity: float
    depth_ratio: float
    depth_ev: float
    q: float
    feeds: tuple[int, ...]
    currents: CArr
    profile: CArr


def _resolve(path: str | Path) -> Path:
    out = Path(path)
    return out if out.is_absolute() else ROOT / out


def _span(lo: float, hi: float, step: float, centre: float) -> Arr:
    """Nodes `step` apart from about lo to about hi, one of them on `centre`."""
    below, above = math.floor((centre - lo) / step + 1e-9), math.floor((hi - centre) / step + 1e-9)
    out: Arr = centre + step * np.arange(-below, above + 1, dtype=np.float64)
    return out


class TwinGuide:
    """One realisation of a design, as a device of the twin bench (docs/KNOTTED_GUIDE.md section 8.3): what the
    resonance (tune), the drive (drive) and the field (field) of the guide built once would be. Options, as a bench
    file gives them: the design file (relative to the repository), the realisation (0, the default, is the design as
    drawn), and the cache folder for solved sections (.twin_cache/ by default). Nothing is solved until it is asked
    for; a solve is read back from the cache when its content has been solved before."""

    def __init__(self, design: str | Path, realisation: int = 0, cache: str | Path | None = None) -> None:
        try:
            number = -1 if isinstance(realisation, bool) else operator.index(realisation)
        except TypeError:
            number = -1
        if number < 0:
            raise ValueError(f"a realisation is a whole number from 0 (the design as drawn), not {realisation!r}")
        self.design = Design.load(_resolve(design))
        self.realisation = number
        self.cache = CACHE if cache is None else _resolve(cache)
        self._grids: dict[str, SectionGrid] = {}
        self._tunings: dict[int, Tuning] = {}

    # ---- what the design takes ---------------------------------------------------------------------------------
    @cached_property
    def solved(self) -> Solved:
        """The design's cross-section, solved (or read back from the cache)."""
        return solve_cached(self.design.section, self.cache)

    @cached_property
    def length(self) -> float:
        """The ring's length (m): the design's, or m guided waves at the band on the solved index."""
        given = self.design.path.length
        return given if given is not None else self.design.m * self.solved.wavelength

    @cached_property
    def curve(self) -> ClosedCurve | None:
        """The path as a closed curve of the design's samples (None for a straight path)."""
        path, length = self.design.path, self.length
        if path.kind == "circle":
            return ClosedCurve.circle(length / (2 * math.pi), path.samples)
        if path.kind == "torus":
            knot = ClosedCurve.torus_knot(1.0, path.ratio, path.samples, path.p, path.q)
            return ClosedCurve(knot.points * (length / knot.length), length)
        return None

    @cached_property
    def framing(self) -> Framing | None:
        """The closing framing, where the design asks for one."""
        frame = self.design.frame
        if frame.k is None or self.curve is None:
            return None
        return framing(self.curve, frame.k, frame.theta0)

    @cached_property
    def build(self) -> Build:
        """What the ring's tolerances act through: the design's stations and sensitivities, the solved Z0, and the
        index's change with the laminate's permittivity from the solved share, n share / (2 er) (d ln n / d ln er =
        share / 2, DR-0's share rule)."""
        tol, solved, laminate = self.design.tolerance, self.solved, self.design.section.laminate
        per_eps = solved.index * solved.share / (2 * laminate.er) if laminate is not None else 0.0
        return Build(
            nodes=RING_NODES,
            spacers=tol.stations,
            sensitivity=tol.sensitivity,
            scarf=tol.scarf,
            joint_jitter=tol.jitter,
            split_sensitivity=tol.split_sensitivity,
            split_stations=tol.split_stations,
            index_per_eps=per_eps,
            z0=solved.z0,
        )

    @cached_property
    def drawn(self) -> Realisation:
        """This realisation's draws: none for realisation 0; else from the seed sequence (design seed, realisation)."""
        build, index = self.build, self.solved.index
        if self.realisation == 0:
            z, n, rel = perfect_ring(index, build.nodes, build.z0)
            return Realisation(0, z, n, rel, np.zeros(3))
        tol = self.design.tolerance
        rng = np.random.default_rng(np.random.SeedSequence([self.design.seed, self.realisation]))
        z, n, rel = ring_draw(
            rng, self.length, tol.spacers, index, tol.joint_sd, tol.split_sd, tol.joint_dz, tol.eps, build
        )
        way = rng.normal(size=3)
        return Realisation(self.realisation, z, n, rel, tol.earth * way / np.linalg.norm(way))

    def grid(self, model: Literal["solved", "gapless"] = "solved") -> SectionGrid:
        """The section's potential per volt round the null, the aperture either way: the solved section (copper and
        boards solid), or the gapless closed form of electrodes.py (RailPair, FacingBoards) sampled h/40 apart."""
        if model not in self._grids:
            if model == "solved":
                solved = self.solved
                self._grids[model] = solved.section_grid(self.design.aperture * solved.distance)
            elif model == "gapless":
                self._grids[model] = self._gapless()
            else:
                raise ValueError(f"a section model is 'solved' or 'gapless', not {model!r}")
        return self._grids[model]

    def _gapless(self) -> SectionGrid:
        geometry, a = self.design.section.geometry, self.design.aperture
        if isinstance(geometry, Rails):
            pair = RailPair(geometry.gap, geometry.rail)
            h = pair.null_height
            u = _span(-a * h, a * h, GAPLESS_STEP * h, 0.0)
            v = _span(max(0.02 * h, (1 - a) * h), (1 + a) * h, GAPLESS_STEP * h, h)
            return SectionGrid.sampled(lambda pu, pv: pair.potential(pu, pv, 1.0), u, v, (0.0, h))
        if isinstance(geometry, FacingStrips):
            boards = FacingBoards(geometry.gap)
            h = boards.distance
            u = _span(-a * h, a * h, GAPLESS_STEP * h, 0.0)
            v = _span(0.0, 2 * h, GAPLESS_STEP * h, h)
            return SectionGrid.sampled(lambda pu, pv: boards.potential(pu, pv - h, 1.0), u, v, (0.0, h))
        raise ValueError("the slot guide has no gapless closed form: its field is the solved section's")

    def sweep(self, grid: SectionGrid, drive: PathDrive) -> SectionField:
        """A section grid swept along the design's path in its frame, driven by any PathDrive: field() with the
        realised ring's volts, or a StandingWave, an InPhase or a static drive of the experiment's own."""
        path, frame = self.design.path, self.design.frame
        if self.curve is None:
            return SectionField.straight(grid, self.length, drive, samples=path.samples, theta=frame.theta)
        if self.framing is not None:
            return SectionField.along(grid, self.curve, drive, normals=self.framing.across)
        return SectionField.along(grid, self.curve, drive, theta=frame.theta)

    # ---- what the device answers -------------------------------------------------------------------------------
    def tune(self, m: int | None = None) -> Tuning:
        """The realised ring's +-m pair (the design's m by default), from its ring of segments losing the solved
        line's Q (array16.guide.ring.ring_resonances)."""
        order = self.design.m if m is None else m
        if order not in self._tunings:
            drawn, q0 = self.drawn, self.solved.quality
            pair = ring_resonances(self.length, drawn.z, drawn.n, q0, order)
            z0 = self.solved.z0
            self._tunings[order] = Tuning(
                f=pair.f,
                q=pair.q,
                lower=pair.lower,
                upper=pair.upper,
                splitting=pair.splitting,
                linewidths=pair.linewidths,
                band=self.design.f,
                m=order,
                length=self.length,
                index=float(np.mean(drawn.n)),
                z0=z0,
                shunt=z0 * q0 / (2 * math.pi * order),
            )
        return self._tunings[order]

    def feed_nodes(self, count: int | None = None) -> tuple[int, ...]:
        """The ring nodes the feeds drive (the design's count by default), as array16.guide.ring.feed_layout places
        them from node 0: two a quarter guided wave apart, more equally spaced, alternate ones moved a quarter wave
        where their count divides 2m."""
        k, nodes, m = self.design.feeds if count is None else count, self.build.nodes, self.design.m
        places = feed_layout(k, float(nodes), m, 0, False)
        return tuple(round(p - 0.03 * nodes) % nodes for p in places)

    def drive(
        self,
        volts: float | None = None,
        *,
        watts: float | None = None,
        q: float | None = None,
        feeds: int | None = None,
    ) -> Driven:
        """The realised ring driven at its tuned centre through compiled current feeds (array16.guide.ring.
        compile_feeds): to `volts` peak on average round the ring, or with `watts`, or to the null's Mathieu `q` at
        the mean volts; with none of them, as the design says."""
        if sum(given is not None for given in (volts, watts, q)) > 1:
            raise ValueError("a drive is set by its volts, its watts or its q, one of them")
        design, solved, drawn = self.design, self.solved, self.drawn
        tuning = self.tune()
        f = tuning.f
        nodes = self.feed_nodes(feeds)
        score = compile_feeds(self.length, design.m, drawn.z, drawn.n, drawn.rel, solved.quality, nodes, f)
        power = float(score.depth.min()) / score.score  # W per (A of the first feed)^2
        mean = float(np.mean(np.abs(score.v)))
        if watts is not None:
            scale = math.sqrt(watts / power)
        else:
            if volts is None and q is None:
                volts, q = design.volts, design.q
            if volts is None:
                assert q is not None
                volts = q / mathieu_q(solved.gradient, f)
            scale = volts / mean
        profile: CArr = scale * score.v
        size = np.abs(profile)
        level = float(np.mean(size))
        barrier = float(np.min(size / (1 + drawn.rel)))
        return Driven(
            f=f,
            watts=scale * scale * power,
            volts=level,
            least=float(size.min()),
            most=float(size.max()),
            purity=score.purity,
            depth_ratio=score.depth_ratio,
            depth_ev=ponderomotive_ev(solved.escape * barrier, f),
            q=mathieu_q(solved.gradient * level, f),
            feeds=nodes,
            currents=scale * score.weights,
            profile=profile,
        )

    def field(self, driven: Driven | None = None, *, model: Literal["solved", "gapless"] = "solved") -> TwinField:
        """The realisation's field, a Fields for integrate: the section (`model`) swept along the path and driven by
        the realised ring's volts (`driven`, the design's drive by default) as a RingWave, with Earth's field where the
        realisation draws it. A gap error dh/h at a node scales the null's gradient there by (1 + dh/h)^-2, which the
        volts carry, so that the null's q and secular frequency follow the drawn gap; the section keeps the shape the
        design draws, so its barrier field is scaled as the gradient is, (1 + dh/h)^-1 more than the gap gives."""
        driven = self.drive() if driven is None else driven
        drawn = self.drawn
        wave = RingWave.of(driven.profile / (1 + drawn.rel) ** 2, 2 * math.pi * driven.f, self.design.m)
        earth = drawn.earth if bool(np.any(drawn.earth)) else None
        return TwinField(self.sweep(self.grid(model), wave), self.curve, earth)

    def identity(self) -> dict[str, Any]:
        """What a record keeps of this device: the design file and its sha256, the realisation, and the section's solve
        key and solver. Nothing is solved to say it."""
        return {
            "design": self.design.source,
            "name": self.design.name,
            "sha256": self.design.digest,
            "realisation": self.realisation,
            "section": section_key(self.design.section)[:16],
            "solver": solver_version(),
        }


# ---- the detector ---------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Pickup:
    """An electrode along the path that a passing charge induces its image on: centred a fraction `at` of a lap round,
    `length` long (m), taking `coupling` of a charge's image while the charge is beside it, the share falling off as
    cos^2 over `edge` (m) at either end."""

    at: float
    length: float
    coupling: float
    edge: float

    def weight(self, s: Arr, lap: float) -> Arr:
        """The share of a charge's image the electrode takes, at arc lengths s (m) round a lap of `lap` m."""
        d = np.abs(np.mod(s - self.at * lap + lap / 2, lap) - lap / 2)
        inner, outer = self.length / 2 - self.edge / 2, self.length / 2 + self.edge / 2
        ramp = np.cos(math.pi / 2 * np.clip((d - inner) / self.edge, 0.0, 1.0)) ** 2
        out: Arr = self.coupling * np.where(d <= inner, 1.0, np.where(d >= outer, 0.0, ramp))
        return out


@dataclass(frozen=True, eq=False)
class Detection:
    """What the detector saw: the particles held, lost to the guide (off its grid: copper, a board, or out of the
    aperture) and stopped by another observer (the gas); each particle's fate (0 held, 1 lost, 2 stopped), the step,
    time (s) and arc length (m) it ended at (-1 and nan for the held), and the laps it went round; the pickups'
    currents (A, (T, P)) at the times recorded (s); the charge the held particles bring to the dump when the drive
    stops (C) and the volts it puts on the collector's charge amplifier (V); the mean current lost to the guide
    over the run (A); and, where the watch was asked for them, the entries of the electrons lost and stopped
    (array16.twin.records.Losses: where each crossed onto the copper, a board or the grid's edge, or was stopped),
    else None."""

    held: int
    lost: int
    stopped: int
    fate: IntArr
    step: IntArr
    time: Arr
    at: Arr
    laps: Arr
    times: Arr
    currents: Arr
    dump_charge: float
    dump_volts: float
    loss_current: float
    losses: Losses | None = None


class Watch:
    """An observer for integrate, made by TwinDetector.watch for one run: after each step it stops the particles the
    field found off its grid (a step after they left it, as SectionField.stop_off_grid), notes those another observer
    stopped, follows each particle's arc length and laps along the path (from `start`, where the particles were
    launched, or else from where the first step left them), and records the current each pickup carries, the change
    of the image charge induced on it over the step. With `losses` it also writes an entry for each electron lost or
    stopped (array16.twin.records.LossLog), one an electron and nothing a step: everything else it records is the
    same either way."""

    def __init__(
        self,
        field: TwinField,
        count: int,
        pickups: Sequence[Pickup],
        charge: float,
        weight: float,
        collector: float,
        every: int = 1,
        start: Arr | None = None,
        *,
        losses: bool = False,
    ) -> None:
        self.field, self.count, self.pickups = field, count, tuple(pickups)
        self.charge, self.weight, self.collector, self.every = charge, weight, collector, every
        self.fate: IntArr = np.zeros(count, dtype=np.intp)
        self.step: IntArr = np.full(count, -1, dtype=np.intp)
        self.time: Arr = np.full(count, math.nan)
        self.at: Arr = np.full(count, math.nan)
        self.laps: Arr = np.zeros(count)
        self.times: list[float] = []
        self.currents: list[Arr] = []
        self._station: IntArr | None = None
        self._advance: Arr = np.zeros(count)
        self._start: Arr | None = None
        self._induced: Arr | None = None
        self._then = 0.0
        self._last = 0.0
        self._log = LossLog(field.swept, count, start) if losses else None
        if start is not None:
            s, _ = self._arc(np.asarray(start, dtype=np.float64))
            if self.pickups:
                self._induced = self._image(s, np.ones(count, dtype=np.bool_))

    def _arc(self, x: Arr) -> tuple[Arr, Arr]:
        """Each particle's arc length (m, 0 to L) and how far it has gone along the path since the first call (m)."""
        swept, lap = self.field.swept, self.field.length
        clean = np.where(np.isfinite(x), x, 0.0)
        if swept.axis is not None:
            along = (clean - swept.points[0]) @ swept.axis
            if self._start is None:
                self._start = along.copy()
            return np.mod(along, lap), along - self._start
        curve = self.field.curve
        assert curve is not None
        if self._station is None:  # the swept field's own search of every sample: the segment each one is on
            self._station = swept.locate(clean).segment.copy()
        else:  # a particle moves at most four samples a step: 0.4 mm on 2000 samples of R1, keV electrons
            self._station, self._advance = curve.advance(clean, self._station, self._advance, iters=4)
        st, ds = self._station, lap / len(curve)
        tangential = ((clean - curve.points[st]) * curve.tangents[st]).sum(1)
        gone = self._advance * ds + tangential
        if self._start is None:
            self._start = gone.copy()
        return np.mod(st * ds + tangential, lap), gone - self._start

    def __call__(self, step: int, t: float, x: Arr, u: Arr, alive: BoolArr) -> BoolArr:
        off = self.field.off_grid
        hit: BoolArr = np.zeros(self.count, dtype=np.bool_)
        if off.shape == alive.shape:
            hit = off & alive & (self.fate == 0)
        others = ~alive & (self.fate == 0)
        s, gone = self._arc(x)
        for mask, fate in ((hit, 1), (others, 2)):
            idx = np.flatnonzero(mask)
            self.fate[idx], self.step[idx], self.time[idx], self.at[idx] = fate, step, t, s[idx]
        going = self.fate == 0
        self.laps = np.where(going, gone / self.field.length, self.laps)
        if self._log is not None:
            self._log.note(step, t, x, u, hit, others, self.laps, gone)
        self._last = t
        if self.pickups and (step + 1) % self.every == 0:
            induced = self._image(s, going)
            if self._induced is not None and t > self._then:
                self.times.append(t)
                self.currents.append((induced - self._induced) / (t - self._then))
            self._induced, self._then = induced, t
        return hit

    def _image(self, s: Arr, going: BoolArr) -> Arr:
        """The image charge (C) the particles still going induce on each pickup."""
        lap, scale = self.field.length, -self.charge * self.weight
        return np.array([scale * float(p.weight(s[going], lap).sum()) for p in self.pickups])

    def reading(self) -> Detection:
        """What the watch has seen so far."""
        held = int((self.fate == 0).sum())
        lost = int((self.fate == 1).sum())
        dump = held * abs(self.charge) * self.weight
        return Detection(
            held=held,
            lost=lost,
            stopped=int((self.fate == 2).sum()),
            fate=self.fate.copy(),
            step=self.step.copy(),
            time=self.time.copy(),
            at=self.at.copy(),
            laps=self.laps.copy(),
            times=np.array(self.times),
            currents=np.array(self.currents).reshape(len(self.times), len(self.pickups)),
            dump_charge=dump,
            dump_volts=dump / self.collector,
            loss_current=lost * abs(self.charge) * self.weight / self._last if self._last > 0 else 0.0,
            losses=None if self._log is None else self._log.losses(),
        )


class TwinDetector:
    """The detector of the twin bench (docs/KNOTTED_GUIDE.md sections 6.7 and 8.3): pickups along the path, the dump's
    collector with its charge amplifier, and the count of what is lost. Options, as a bench file gives them: the
    pickups' places (fractions of a lap), their length (m), coupling (the share of a passing charge's image each
    takes; 30 % in section 6.5) and edge (m); the collector's capacitance (F); how many particles each tracked one
    stands for; and whether each watch writes an entry for every electron lost (`losses`: where along the path and in
    the section it left, when, with what energy and sense, and how; array16.twin.records), which a shot may also ask
    of its own watch. `watch` makes the observer for one run."""

    def __init__(
        self,
        pickups: Sequence[float] = (0.0, 0.25),
        length: float = 10e-3,
        coupling: float = 0.30,
        edge: float = 1e-3,
        collector: float = 1e-12,
        weight: float = 1.0,
        every: int = 1,
        losses: bool = False,
    ) -> None:
        if not (length > edge > 0 and 0 <= coupling <= 1 and collector > 0 and weight > 0 and every >= 1):
            raise ValueError(
                "a detector takes pickups longer than their edge, a coupling from 0 to 1, a positive collector, a "
                f"positive weight and a recording interval of a step or more, not {length}, {edge}, {coupling}, "
                f"{collector}, {weight} and {every}"
            )
        if not isinstance(losses, bool):
            raise ValueError(
                f"a detector writes the entries of the electrons lost or not, losses true or false, not {losses!r}"
            )
        self.pickups = tuple(Pickup(float(at), length, coupling, edge) for at in pickups)
        self.collector, self.weight, self.every, self.losses = collector, weight, every, losses

    def watch(
        self,
        field: TwinField,
        count: int,
        *,
        charge: float = -E,
        weight: float | None = None,
        start: Arr | None = None,
        losses: bool | None = None,
    ) -> Watch:
        """The observer for one run of `count` particles of `charge` (C, an electron's by default) in `field`, counting
        their laps from `start` (N, 3), where they were launched, when it is given; writing an entry for each electron
        lost when `losses` asks it (the detector's own setting by default), from the launch on when `start` is given."""
        heft = self.weight if weight is None else weight
        keep = self.losses if losses is None else losses
        return Watch(field, count, self.pickups, charge, heft, self.collector, self.every, start, losses=keep)

    def identity(self) -> dict[str, Any]:
        """What a record keeps of the detector: its pickups, collector and weight, and that it writes the entries of the
        electrons lost where it does."""
        out: dict[str, Any] = {"pickups": plain(self.pickups), "collector_F": self.collector, "weight": self.weight}
        if self.losses:
            out["losses"] = True
        return out


# ---- the vacuum -----------------------------------------------------------------------------------------------------
class Collisions:
    """An observer for integrate that stops each electron the gas scatters: over a step it goes |v| dt and is scattered
    with probability 1 - exp(-|v| dt / free path). Seeded; `count` is how many it has stopped."""

    def __init__(self, free_path: Callable[[Arr], Arr], seed: int) -> None:
        self.free_path = free_path
        self.rng = np.random.default_rng(seed)
        self.count = 0

    def chance(self, step: int, t: float, u: Arr) -> Arr:
        """Each electron's chance of being scattered over the step just taken: integrate calls its observers after
        step `step` at t = (step + 1) dt, so the step lasted t / (step + 1); 1 - exp(-|v| dt / free path)."""
        dt = t / (step + 1)
        gamma = np.sqrt(1 + (u * u).sum(1) / C0**2)
        speed = np.sqrt((u * u).sum(1)) / gamma
        energy = (gamma - 1) * ME * C0**2 / E
        out: Arr = -np.expm1(-speed * dt / self.free_path(energy))
        return out

    def __call__(self, step: int, t: float, x: Arr, u: Arr, alive: BoolArr) -> BoolArr:
        hit: BoolArr = alive & (self.rng.random(len(x)) < self.chance(step, t, u))
        self.count += int(hit.sum())
        return hit


class TwinVacuum:
    """The vacuum round the guide (docs/KNOTTED_GUIDE.md sections 6.5, 8.2 and 8.3): a residual gas at `pressure`
    (mbar) and `temperature` (K), scattering electrons with one total cross-section (m^2) at every energy. The default
    cross-section is the one that puts the mean free path at 0.4 m at 1e-3 mbar and 295 K, the estimate those sections
    use (1.02e-19 m^2: 11 us between collisions at 1 eV and 6.3e-5 mbar, 20 laps of T915-S). LXCat's cross-sections
    for a named gas will take its place by energy; until they are read in, the gas is the residual gas, named
    "residual"."""

    def __init__(
        self,
        pressure: float = 1e-6,
        temperature: float = 295.0,
        cross_section: float | None = None,
        gas: str = "residual",
    ) -> None:
        path, at, warm = FREE_PATH
        sigma = KB * warm / (at * 100.0 * path) if cross_section is None else cross_section
        if not (pressure > 0 and temperature > 0 and sigma > 0):
            raise ValueError(
                f"a vacuum takes a positive pressure, temperature and cross-section, not {pressure} mbar, "
                f"{temperature} K and {sigma} m^2"
            )
        if gas != "residual":
            raise ValueError(
                f"the vacuum's gas is the residual gas, 'residual', until LXCat's are read in, not {gas!r}"
            )
        self._mbar, self.temperature, self.cross_section, self.gas = pressure, temperature, sigma, gas

    def pressure(self) -> float:
        """The gauge's reading (mbar)."""
        return self._mbar

    def density(self) -> float:
        """The gas's number density (1/m^3): p / (k T)."""
        return self._mbar * 100.0 / (KB * self.temperature)

    def free_path(self, energy_ev: Arr | float = 1.0) -> Arr:
        """An electron's mean free path (m) at energies (eV): 1 / (n sigma), the same at every energy."""
        out: Arr = np.full(np.shape(energy_ev), 1 / (self.density() * self.cross_section))
        return out

    def lifetime(self, energy_ev: float = 1.0) -> float:
        """The mean time (s) between an electron's collisions at an energy (eV): the free path over its speed."""
        gamma = 1 + energy_ev * E / (ME * C0**2)
        return float(self.free_path(energy_ev)) / (C0 * math.sqrt(1 - 1 / gamma**2))

    def scatter(self, seed: int = 0) -> Collisions:
        """The observer that stops the electrons the gas scatters, seeded."""
        return Collisions(self.free_path, seed)

    def identity(self) -> dict[str, Any]:
        return {
            "pressure_mbar": self._mbar,
            "temperature_K": self.temperature,
            "cross_section_m2": self.cross_section,
            "gas": self.gas,
        }
