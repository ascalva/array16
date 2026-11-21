"""The figures of docs/VIRTUAL_ANTENNA.md: what the field compiler's theory looks like on array16's own elements.

The document states a generic process that makes the array act as any antenna whose field over a region lies in
the span of its sixteen element fields. This script computes the figures that carry that theory, each from the
same Hertzian model the document's numbers came from (section of the document; figure number):

  compiler     the field compiler as a flow, each step's formula in its box (3.4; figure 1)
  price        the price list 1/sigma_i^2 of the tile, the whips and a reference shell (4.1, 4.3; 2)
  radiation    the eigenvalues of the radiation matrix R across 0.4-6 GHz (4.2; 3)
  fine         the power that holds alternating values at sub-wavelength spacing (4.4; 4)
  ceiling      a focus against the free-space ceiling RFG-2, for shells and for the board's arrays (4.5; 5)
  band         how far the operator moves across frequency (4.8; 6)
  behind       a virtual dipole behind the tile: error map, window, cone residuals (5.3; 7)
  focused      a virtual source focused in front of the tile, and the tile's on-axis maxima (5.4; 8)
  whips        the whip rings as an interior array: orders, their cost, Bessel weights, virtual whips (5.5; 9)
  alphabets    the weights the states reach: hulls for line 0, discs for line 1 (6.2; 10)
  sequences    state rasters, line weights and line spectra of a staircase, Frank shifts, Hadamard words (6.3; 11)
  lines        recovery of every element's signal by slot means and by comb lines (6.4; 12)
  nulls        line-1 null depth against sequence length, and a realised pattern against its target (6.5; 13)
  dark         the tile's drive space: radiating and reactive drives, and a dark drive's field (7.1; 14)
  selection    the ring selection rule, F_l(kR), and the radiate-or-hold trade (7.2; 15)
  containment  the power a 1 eV well takes, free space against electrodes at the null's Mathieu q; the facing boards
               and one board (7.3; 16)
  charges      a held column's signals and its equivalent circuit (7.4; 17)
  modules      the compiler's modules, types and bench roles (9; 18)
  bench        the V1 bench in 3-D, its signal chains, and its predicted held-out residuals (10.1; 19)

    uv run python experiments/virtual_antenna_study.py              # every figure, about half a minute
    uv run python experiments/virtual_antenna_study.py --quick      # coarse grids and few draws: virtual_antenna/quick/
    uv run python experiments/virtual_antenna_study.py --only price,behind

Model (experiments/local_fields.py, ported here typed): every element is a set of Hertzian current elements with
the exact dipole field, every near-field term kept; an infinite perfect ground by images where the array has one;
no mutual coupling, so the drives are the element currents. The tile is sixteen y-directed dipoles 13 mm over
ground at 25 mm pitch; the whips are two rings of eight monopoles (225 mm, and 112 mm turned 22.5 degrees) with a
sinusoidal current; a reference shell is dipoles on a Fibonacci sphere, tangential, alternating in polarisation.
The radiation matrix is the closed form of the dipoles' mutual radiation resistance, checked against the
far-field quadrature in selftest. Outputs: virtual_antenna/{*.png, results.json, run.log}.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from array16 import lab
from array16 import particles as pt
from array16.control import H16
from array16.dsp import slot_means
from array16.particles import Arr, CArr, IntArr

OUT = Path(__file__).resolve().parent / "virtual_antenna"
LOG: list[str] = []
C0, Z0 = pt.C0, pt.Z0

# Rev B's delay bits (REDESIGN.md section 3) and their sixteen sums; a delay tau multiplies the phasor by
# exp(+i 2 pi f tau) under the exp(-i omega t) convention of the field.
DELAY_BITS = (38.5e-12, 52.5e-12, 117.5e-12, 139.0e-12)
DELAYS = np.array(sorted(sum(c) for r in range(5) for c in itertools.combinations(DELAY_BITS, r)))
LOSS = 0.1  # an ohmic loss of a tenth of an element's radiated power: R + r I with r = 0.1 mean(diag R)

# Ink, surfaces and the categorical order of the reference palette (fixed order, never cycled).
INK, INK2, MUTED, GRID, EDGE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
BLUES = ("#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")


def say(line: str = "") -> None:
    print(line, flush=True)
    LOG.append(line)


def db(x: Any) -> Any:
    """10 log10 of a power ratio, floored at 1e-300."""
    return 10 * np.log10(np.maximum(np.asarray(x, dtype=np.float64), 1e-300))


# ---- the element model ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Elements:
    """N elements, each a weighted set of Hertzian current elements: the field is field(el, points) @ w.

    positions and moments (A m per unit drive) list every current element, images included when the array stands
    on a ground plane; transfer (S, N) maps the N drives onto them; feeds (N, 3) are the elements' feed points."""

    name: str
    f: float
    positions: Arr
    moments: Arr
    transfer: Arr
    feeds: Arr
    ground: bool

    @property
    def k(self) -> float:
        return 2 * math.pi * self.f / C0

    @property
    def wavelength(self) -> float:
        return C0 / self.f

    @property
    def count(self) -> int:
        return int(self.transfer.shape[1])


def _elements(
    name: str, f: float, rows: Sequence[tuple[Sequence[float], Sequence[float], int]], feeds: Arr, ground: bool
) -> Elements:
    """Elements from (position, moment, drive index) rows."""
    count = len(feeds)
    transfer = np.zeros((len(rows), count))
    for s, (_, _, n) in enumerate(rows):
        transfer[s, n] = 1.0
    return Elements(
        name,
        f,
        np.array([r[0] for r in rows], dtype=np.float64),
        np.array([r[1] for r in rows], dtype=np.float64),
        transfer,
        np.asarray(feeds, dtype=np.float64),
        ground,
    )


def panel(
    f: float, nx: int = 4, ny: int = 4, pitch: float = 0.025, height: float = 0.013, ground: bool = True
) -> Elements:
    """The tile: nx by ny y-directed dipoles at height over an infinite ground (image of opposite sign), element
    p = iy nx + ix. Without ground, the same dipoles alone in free space."""
    xs = (np.arange(nx) - (nx - 1) / 2) * pitch
    ys = (np.arange(ny) - (ny - 1) / 2) * pitch
    gx, gy = np.meshgrid(xs, ys)
    xy = np.c_[gx.ravel(), gy.ravel()]
    layers = ((1.0, height), (-1.0, -height)) if ground else ((1.0, height),)
    rows = [((x, y, z), (0.0, s, 0.0), e) for e, (x, y) in enumerate(xy) for s, z in layers]
    label = f"tile {f / 1e9:g} GHz" if ground else f"tile in free space {f / 1e9:g} GHz"
    return _elements(label, f, rows, np.c_[xy, np.full(len(xy), height)], ground)


def monopoles(f: float, xy: Any, length: float, segments: int = 8, name: str = "monopoles") -> Elements:
    """Vertical monopoles on an infinite ground at the points xy, each length long, carrying the thin-wire current
    sin(k (length - z)) in segments, with their images."""
    places = np.atleast_2d(np.asarray(xy, dtype=np.float64))
    k = 2 * math.pi * f / C0
    zs = (np.arange(segments) + 0.5) * length / segments
    current = np.sin(k * (length - zs)) * length / segments
    rows = [
        ((x, y, sz), (0.0, 0.0, c), e)
        for e, (x, y) in enumerate(places)
        for sz, c in zip(np.r_[zs, -zs], np.r_[current, current], strict=True)
    ]
    return _elements(name, f, rows, np.c_[places, np.zeros(len(places))], True)


def whip_places(radii: tuple[float, float] = (0.225, 0.112)) -> Arr:
    """The sixteen whip ports: an outer ring of eight and an inner ring of eight turned 22.5 degrees."""
    angle = np.arange(8) * math.pi / 4
    outer = np.c_[radii[0] * np.cos(angle), radii[0] * np.sin(angle)]
    inner = np.c_[radii[1] * np.cos(angle + math.pi / 8), radii[1] * np.sin(angle + math.pi / 8)]
    places: Arr = np.r_[outer, inner]
    return places


def whips(f: float, length: float | None = None, select: Sequence[int] | None = None, segments: int = 4) -> Elements:
    """The whip rings, each whip a quarter wave at f unless a fixed length is given."""
    places = whip_places()
    if select is not None:
        places = places[list(select)]
    cut = C0 / f / 4 if length is None else length
    label = f"{len(places)} whips {f / 1e9:g} GHz" + ("" if length is None else f", {cut * 1e3:.0f} mm")
    return monopoles(f, places, cut, segments, label)


def fibonacci(
    n: int, ct_lo: float, ct_hi: float, radius: float = 1.0, origin: Any = (0.0, 0.0, 0.0), ascending: bool = False
) -> Arr:
    """n points spread evenly over the zone ct_lo <= cos(theta) <= ct_hi of a sphere (a Fibonacci lattice)."""
    i = np.arange(n) + 0.5
    ct = ct_lo + (ct_hi - ct_lo) * i / n if ascending else ct_hi - (ct_hi - ct_lo) * i / n
    phi = math.pi * (1 + math.sqrt(5)) * i
    st = np.sqrt(1 - ct**2)
    points: Arr = radius * np.c_[st * np.cos(phi), st * np.sin(phi), ct] + np.asarray(origin, dtype=np.float64)
    return points


def shell(f: float, count: int, radius: float) -> Elements:
    """A reference shell: count Hertzian dipoles on a Fibonacci sphere of the given radius in free space, tangential,
    alternately along theta-hat and phi-hat."""
    p = fibonacci(count, -1.0, 1.0, radius)
    ct = p[:, 2] / radius
    st = np.sqrt(1 - ct**2)
    phi = np.arctan2(p[:, 1], p[:, 0])
    theta_hat = np.c_[ct * np.cos(phi), ct * np.sin(phi), -st]
    phi_hat = np.c_[-np.sin(phi), np.cos(phi), np.zeros(count)]
    moments = np.where((np.arange(count) % 2 == 0)[:, None], theta_hat, phi_hat)
    rows = [(tuple(p[n]), tuple(moments[n]), n) for n in range(count)]
    return _elements(f"{count}-element shell r {radius * 1e3:.0f} mm {f / 1e9:g} GHz", f, rows, p, False)


def dipole(f: float, position: Sequence[float], moment: Sequence[float]) -> Elements:
    """One Hertzian dipole in free space: a target field."""
    return _elements("dipole", f, [(tuple(position), tuple(moment), 0)], np.atleast_2d(np.asarray(position)), False)


def field(el: Elements, points: Any, component: int | None = None, chunk: int = 4096) -> CArr:
    """The field at points of each unit drive: (P, 3, N), or (P, N) for one Cartesian component. The exact field
    of a current element m: E = (i Z0 k / 4 pi) exp(ikr)/r {m - n(n.m) + [3 n(n.m) - m](1/(kr)^2 - i/(kr))}."""
    pts = np.atleast_2d(np.asarray(points, dtype=np.float64))
    out = np.empty((len(pts), 3, el.count), dtype=np.complex128)
    k = el.k
    pref = 1j * Z0 * k / (4 * math.pi)
    m = el.moments[None]
    for a in range(0, len(pts), chunk):
        d = pts[a : a + chunk, None, :] - el.positions[None]
        r = np.linalg.norm(d, axis=2)
        n = d / r[..., None]
        kr = k * r
        nm = np.einsum("psc,sc->ps", n, el.moments)[..., None]
        g = (pref * np.exp(1j * kr) / r)[..., None]
        near = (1 / kr**2 - 1j / kr)[..., None]
        e = g * ((m - n * nm) + (3 * n * nm - m) * near)
        out[a : a + chunk] = np.einsum("psc,sn->pcn", e, el.transfer)
    if component is None:
        return out
    one: CArr = out[:, component, :]
    return one


def stacked(el: Elements, points: Any) -> CArr:
    """field with the three components of every point stacked: (3P, N)."""
    g: CArr = field(el, points).reshape(-1, el.count)
    return g


def radiation(el: Elements) -> Arr:
    """The radiation matrix: P = w^H R w watts for peak-amplitude drives w. In closed form, the mutual radiation
    resistance of current elements a, b a distance r apart is (Z0 k^2 / 8 pi)[A m_a.m_b + B (m_a.n)(m_b.n)],
    A = sin x/x + cos x/x^2 - sin x/x^3, B = -sin x/x - 3 cos x/x^2 + 3 sin x/x^3, x = kr; Z0 k^2 |m|^2 / 12 pi
    alone. Over a ground the images make a field mirror-symmetric about it, so the half space takes half the power
    of sources and images together."""
    k = el.k
    d = el.positions[:, None, :] - el.positions[None, :, :]
    r = np.linalg.norm(d, axis=2)
    x = k * r
    safe = np.where(x > 1e-3, x, 1.0)
    a_far = np.sin(safe) / safe + np.cos(safe) / safe**2 - np.sin(safe) / safe**3
    b_far = -np.sin(safe) / safe - 3 * np.cos(safe) / safe**2 + 3 * np.sin(safe) / safe**3
    a = np.where(x > 1e-3, a_far, 2 / 3 - 2 * x**2 / 15)
    b = np.where(x > 1e-3, b_far, x**2 / 15)
    n = d / np.where(r > 0, r, 1.0)[..., None]
    mm = el.moments @ el.moments.T
    mn = np.einsum("abc,ac->ab", n, el.moments)
    nm = np.einsum("abc,bc->ab", n, el.moments)
    sources = Z0 * k**2 / (8 * math.pi) * (a * mm + b * mn * nm)
    out: Arr = el.transfer.T @ sources @ el.transfer * (0.5 if el.ground else 1.0)
    return 0.5 * (out + out.T)


def radiation_by_quadrature(el: Elements, nth: int = 48, nph: int = 96) -> Arr:
    """The radiation matrix by integrating the far field over the sphere (the upper half over a ground), the way
    experiments/local_fields.py computes it: the check on the closed form."""
    x, wg = np.polynomial.legendre.leggauss(nth)
    ct, wt = (0.5 * (x + 1), 0.5 * wg) if el.ground else (x, wg)
    ph = (np.arange(nph) + 0.5) * 2 * math.pi / nph
    cg, pg = np.meshgrid(ct, ph, indexing="ij")
    sg = np.sqrt(1 - cg**2)
    n = np.stack([sg * np.cos(pg), sg * np.sin(pg), cg], -1).reshape(-1, 3)
    weight = np.repeat(wt, nph) * 2 * math.pi / nph
    nm = n @ el.moments.T
    fs = (el.moments[None] - n[:, None, :] * nm[..., None]) * np.exp(-1j * el.k * (n @ el.positions.T))[..., None]
    far = 1j * Z0 * el.k / (4 * math.pi) * np.einsum("dsc,sn->dcn", fs, el.transfer)
    out: Arr = np.einsum("d,dcn,dcm->nm", weight, far.conj(), far).real / (2 * Z0)
    return out


def lossy(r: Arr, loss: float = LOSS) -> Arr:
    """R + r I with r the given fraction of the mean diagonal: an ohmic loss per element."""
    out: Arr = r + loss * float(np.mean(np.diag(r))) * np.eye(len(r))
    return out


# ---- linear algebra -------------------------------------------------------------------------------------------------
def generalised_eigh(a: Any, b: Any) -> tuple[Arr, CArr]:
    """A v = lambda B v for Hermitian A and positive definite B: eigenvalues ascending, B-orthonormal vectors."""
    inv = np.linalg.inv(np.linalg.cholesky(b))
    m = inv @ a @ inv.conj().T
    val, u = np.linalg.eigh(0.5 * (m + m.conj().T))
    vec: CArr = inv.conj().T @ u
    return np.asarray(val, dtype=np.float64), vec


def loaded_fit(a: CArr, b: CArr, loading: float) -> CArr:
    """argmin |a w - b|^2 + mu |w|^2 with mu = loading times the largest eigenvalue of a^H a (Tikhonov in the
    identity: a cap on white-noise gain)."""
    gram = a.conj().T @ a
    mu = loading * float(np.linalg.eigvalsh(gram)[-1])
    w: CArr = np.linalg.solve(gram + mu * np.eye(len(gram)), a.conj().T @ b)
    return w


def residual_db(a: CArr, b: CArr, w: CArr) -> float:
    """|a w - b|^2 / |b|^2 in dB."""
    return float(db(np.linalg.norm(a @ w - b) ** 2 / np.linalg.norm(b) ** 2))


def scaled_residual_db(a: CArr, b: CArr, w: CArr) -> float:
    """The residual of a w against b after the best complex scale of a w: a drive with one free scale."""
    v = a @ w
    scale = (v.conj() @ b) / (v.conj() @ v)
    return float(db(np.linalg.norm(scale * v - b) ** 2 / np.linalg.norm(b) ** 2))


def ball_cells(centre: Any, radius: float, step: float) -> tuple[Arr, float]:
    """Cell centres of a cubic grid of pitch about step inside a ball, and the volume of a cell."""
    n = max(math.ceil(2 * radius / step), 6)
    g = (np.arange(n) + 0.5) / n * 2 * radius - radius
    gx, gy, gz = np.meshgrid(g, g, g, indexing="ij")
    p = np.c_[gx.ravel(), gy.ravel(), gz.ravel()]
    inside: Arr = p[np.linalg.norm(p, axis=1) <= radius] + np.asarray(centre, dtype=np.float64)
    return inside, (2 * radius / n) ** 3


def square(z: float, half: float, step: float) -> Arr:
    """A square of points 2 half wide on the plane at height z, at the given pitch."""
    xs = np.arange(-half, half + 1e-12, step)
    gx, gy = np.meshgrid(xs, xs)
    out: Arr = np.c_[gx.ravel(), gy.ravel(), np.full(gx.size, z)]
    return out


def bessel_j(m: int, x: Any) -> Arr:
    """J_m(x) = (1/pi) integral over [0, pi] of cos(m t - x sin t), by the trapezoid rule (exact to rounding for
    |x| up to about 40 with 257 nodes, the integrand being smooth and periodic)."""
    xs = np.atleast_1d(np.asarray(x, dtype=np.float64))
    t = np.linspace(0.0, math.pi, 257)
    vals = np.cos(m * t[None, :] - xs[:, None] * np.sin(t)[None, :])
    out: Arr = np.trapezoid(vals, t, axis=1) / math.pi
    return out


def ceiling_per_watt(f: float) -> float:
    """|E|^2 per radiated watt at the centre of a converging electric dipole wave, 4 Z0 k^2 / (3 pi) (RFG-2)."""
    return float(pt.dipole_wave_peak_field(1.0, f) ** 2)


def half_space_factors(k: float, height: float) -> Arr:
    """The ceiling over a perfect ground at a height, per Cartesian component, relative to free space: the image
    takes away 3A/2 of a tangential component's and adds 3(sin x - x cos x)/x^3 to the normal one, x = 2 k h."""
    x = 2 * k * height
    a = math.sin(x) / x + math.cos(x) / x**2 - math.sin(x) / x**3
    tangential = 1 - 1.5 * a
    normal = 1 + 3 * (math.sin(x) - x * math.cos(x)) / x**3
    return np.array([tangential, tangential, normal])


def focus_share(el: Elements, point: Any, loss: float, component: int | None = None, free_space: bool = False) -> float:
    """The most |E|^2 per radiated watt at a point, as a share of the ceiling there: RFG-2's free-space ceiling
    for an array without ground or when free_space is set, else the ceiling over the ground with the image term.
    The largest eigenvalue of G_f (R + rI)^-1 G_f^H, polarisation free, or g^H (R + rI)^-1 g for one component
    (SYN-4)."""
    gf = field(el, [point])[0]
    r = radiation(el)
    rt = lossy(r, loss) if loss > 0 else r + 1e-12 * float(np.mean(np.diag(r))) * np.eye(el.count)
    m = gf @ np.linalg.solve(rt, gf.conj().T)
    ceiling = ceiling_per_watt(el.f)
    scale = half_space_factors(el.k, float(point[2])) if el.ground and not free_space else np.ones(3)
    if component is not None:
        return float(m[component, component].real) / (ceiling * float(scale[component]))
    s = 1 / np.sqrt(scale)
    best = float(np.linalg.eigvalsh(0.5 * (m + m.conj().T) * s[:, None] * s[None, :])[-1])
    return best / ceiling


def selftest() -> dict[str, float]:
    """The model's calibrations: one dipole's radiated power, its near and far fields, and the closed-form radiation
    matrix against the far-field quadrature for the tile and the whips."""
    one = dipole(1e9, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0))
    exact = Z0 * one.k**2 / (12 * math.pi)
    p_closed = float(radiation(one)[0, 0])
    p_quad = float(radiation_by_quadrature(one)[0, 0])
    far = field(one, [[100.0, 0.0, 0.0]])[0, :, 0]
    near = field(one, [[0.0, 0.0, 0.002]])[0, :, 0]
    far_err = abs(abs(far[2]) / (Z0 * one.k / (4 * math.pi * 100)) - 1)
    near_err = abs(abs(near[2]) / (Z0 * 2 / (4 * math.pi * one.k * 0.002**3)) - 1)
    out = {"dipole_power_closed": abs(p_closed / exact - 1), "dipole_power_quadrature": abs(p_quad / exact - 1)}
    out |= {"far_field": far_err, "near_field": near_err}
    for el in (panel(5.8e9), panel(2.3e9), whips(0.9e9)):
        rc, rq = radiation(el), radiation_by_quadrature(el)
        out[f"R {el.name}"] = float(np.abs(rc - rq).max() / np.abs(rq).max())
    assert out["dipole_power_closed"] < 1e-12, out
    assert out["dipole_power_quadrature"] < 1e-6, out
    assert far_err < 1e-3, out
    assert near_err < 1e-2, out
    assert all(v < 1e-6 for key, v in out.items() if key.startswith("R ")), out
    say(
        "selftest: one dipole radiates Z0 k^2 |Il|^2 / 12 pi to "
        f"{out['dipole_power_closed']:.0e} (closed form) and {out['dipole_power_quadrature']:.0e} (quadrature); "
        f"far and near fields to {far_err:.0e} and {near_err:.0e}; closed-form R against quadrature to "
        f"{max(v for key, v in out.items() if key.startswith('R ')):.0e}"
    )
    return out


# ---- 4.1 and 4.3: the price list ------------------------------------------------------------------------------------
def price_list(el: Elements, centre: Any, radius: float, loss: float) -> Arr:
    """sigma_i^2 of G^H G v = sigma^2 (R + rI) v for the vector field in a ball sampled at lambda/6, descending."""
    pts, cell = ball_cells(centre, radius, el.wavelength / 6)
    g = stacked(el, pts) * math.sqrt(cell)
    r = radiation(el)
    rt = lossy(r, loss) if loss > 0 else r
    a = g @ np.linalg.inv(np.linalg.cholesky(rt)).conj().T
    s = np.linalg.svd(a, compute_uv=False)
    out: Arr = s**2
    return out


def region_modes(ka: float) -> int:
    """2L(L+2), L = ceil(ka): the TE and TM spherical waves of degree 1 to L (SYN-5)."""
    order = math.ceil(ka)
    return 2 * order * (order + 2)


def compute_price(quick: bool) -> dict[str, Any]:
    say("\n== price: the price list 1/sigma_i^2 (dB above the cheapest mode) of the vector field in a ball ==")
    cases: list[tuple[str, Elements, tuple[float, float, float], float, float]] = []
    tile_ball = C0 / 5.8e9  # a wavelength at 5.8 GHz, 150 mm in front
    for f in (2.3e9, 3.5e9, 5.8e9):
        el = panel(f)
        cases.append(("tile", el, (0.0, 0.0, 0.15), tile_ball, math.pi * 0.01 / el.wavelength**2))
    for f in (0.433e9, 0.9e9, 2.3e9):
        el = whips(f)
        cases.append(("whips", el, (0.0, 0.0, 0.10), 0.10, 2 * (el.k * 0.225) ** 2))
    el = shell(2.4e9, 16, 0.25)
    cases.append(("shell", el, (0.0, 0.0, 0.0), el.wavelength / 2, 2 * (el.k * 0.25) ** 2))
    rows = []
    for group, el, centre, radius, source in cases:
        ka = el.k * radius
        knee = min(el.count, source, region_modes(ka))
        row: dict[str, Any] = {"group": group, "array": el.name, "f_GHz": el.f / 1e9, "ball_mm": radius * 1e3}
        row |= {"centre_mm": [c * 1e3 for c in centre], "ka": ka, "region_modes": region_modes(ka)}
        row |= {"source_count": source, "knee": knee}
        for loss, key in ((0.0, "lossless"), (LOSS, "loss")):
            s2 = price_list(el, centre, radius, loss)
            price = db(s2[0] / s2)
            row[f"price_db_{key}"] = price
            row[f"within_10dB_{key}"] = int((price < 10).sum())
            row[f"within_20dB_{key}"] = int((price < 20).sum())
        rows.append(row)
        say(
            f"  {el.name:34s} ball {radius * 1e3:5.1f} mm ka {ka:5.2f} 2L(L+2) {row['region_modes']:4d} "
            f"source {source:6.1f} knee {knee:5.1f}: within 10/20 dB {row['within_10dB_loss']}/"
            f"{row['within_20dB_loss']} (10 % loss), {row['within_10dB_lossless']}/{row['within_20dB_lossless']} "
            f"lossless; prices (10 % loss) {np.round(row['price_db_loss'][:10], 1).tolist()}"
        )
    return {"rows": rows}


# ---- 4.2 and 7.1: the radiation matrix across the band -------------------------------------------------------------
def radiation_eigen(el: Elements) -> Arr:
    """The eigenvalues of R, descending, relative to the largest."""
    ev = np.linalg.eigvalsh(radiation(el))[::-1]
    out: Arr = np.maximum(ev / ev[0], 1e-30)
    return out


def compute_radiation(quick: bool) -> dict[str, Any]:
    say("\n== radiation: eigenvalues of R (relative to the largest) from 0.4 to 6 GHz ==")
    freqs = np.linspace(0.4e9, 6.0e9, 29 if quick else 113)
    whip_cut = C0 / 0.9e9 / 4
    tile = np.array([radiation_eigen(panel(float(f))) for f in freqs])
    whip = np.array([radiation_eigen(whips(float(f), length=whip_cut, segments=16)) for f in freqs])
    table = {}
    for f in (2.3e9, 3.5e9, 5.8e9):
        ev = np.linalg.eigvalsh(radiation(panel(f)))[::-1]
        table[f"tile {f / 1e9:g} GHz"] = {
            "within_10dB": int((ev > 0.1 * ev[0]).sum()),
            "condition": float(ev[0] / ev[-1]),
            "pi_A_over_lambda2": math.pi * 0.01 / (C0 / f) ** 2,
            "smallest_below_mean_dB": float(db(ev[-1] / np.mean(ev))),
        }
    for f in (0.433e9, 0.9e9, 2.3e9):
        ev = np.linalg.eigvalsh(radiation(whips(f)))[::-1]
        table[f"whips cut per frequency {f / 1e9:g} GHz"] = {"within_10dB": int((ev > 0.1 * ev[0]).sum())}
        ev = np.linalg.eigvalsh(radiation(whips(f, length=whip_cut, segments=16)))[::-1]
        table[f"whips cut for 0.9 GHz {f / 1e9:g} GHz"] = {"within_10dB": int((ev > 0.1 * ev[0]).sum())}
    for key, row in table.items():
        say(f"  {key:34s} " + ", ".join(f"{k} {v:.3g}" for k, v in row.items()))
    return {"freqs_GHz": freqs / 1e9, "tile_db": db(tile), "whips_db": db(whip), "whip_cut_mm": whip_cut * 1e3} | {
        "table": table
    }


# ---- 4.4: the price of fine structure -------------------------------------------------------------------------------
def alternating_cost(el: Elements, centre: Any, spacing: float, count: int, component: int, loss: float) -> float:
    """The least radiated power (W) holding +1, -1, +1, ... V/m of one component at count points spaced spacing
    along x: e^H (G R^-1 G^H)^-1 e (SYN-2)."""
    x = (np.arange(count) - (count - 1) / 2) * spacing
    pts = np.asarray(centre, dtype=np.float64) + np.c_[x, 0 * x, 0 * x]
    g = field(el, pts, component)
    r = radiation(el)
    rt = lossy(r, loss) if loss > 0 else r
    k = g @ np.linalg.solve(rt, g.conj().T)
    e = (-1.0) ** np.arange(count)
    return float(np.real(e @ np.linalg.solve(k, e)))


def compute_fine(quick: bool) -> dict[str, Any]:
    say("\n== fine: power to hold alternating values at M points, relative to one point ==")
    ref = shell(5.8e9, 64, 0.25)
    lam = ref.wavelength
    spacings = (1.0, 0.5, 0.35, 0.25, 0.15)
    shell_rows = {}
    for d in spacings:
        p = np.array([alternating_cost(ref, (0, 0, 0), d * lam, m, 2, 0.0) for m in range(1, 7)])
        shell_rows[d] = db(p / p[0])
        say(f"  64-element shell, delta {d:4.2f} lambda: {np.round(shell_rows[d], 1).tolist()} dB for M = 1 to 6")
    tile = panel(5.8e9)
    tile_rows = {}
    for d in (1.5, 0.5):
        p = np.array([alternating_cost(tile, (0, 0, 0.15), d * lam, m, 1, LOSS) for m in range(1, 7)])
        tile_rows[d] = db(p / p[0])
        say(f"  tile 150 mm in front (10 % loss), delta {d:3.1f} lambda: {np.round(tile_rows[d], 1).tolist()} dB")
    fifth = float(tile_rows[1.5][4] - tile_rows[1.5][3])
    say(f"  tile, delta 1.5 lambda: four columns hold four values, and the fifth costs {fifth:.1f} dB more")
    dense = np.geomspace(0.1, 1.2, 9 if quick else 23)
    per_point = []
    for d in dense:
        p1 = alternating_cost(ref, (0, 0, 0), d * lam, 1, 2, 0.0)
        p6 = alternating_cost(ref, (0, 0, 0), d * lam, 6, 2, 0.0)
        per_point.append(float(db(p6 / p1)) / 5)
    return {
        "shell": {f"{d:g}": v for d, v in shell_rows.items()},
        "tile": {f"{d:g}": v for d, v in tile_rows.items()},
        "dense_spacing_lambda": dense,
        "dB_per_point": np.array(per_point),
    }


# ---- 4.5: a focus against the ceiling -------------------------------------------------------------------------------
def cavity_resonances(top: float) -> Arr:
    """The x below top where (x j1(x))' = 0, x cos x - sin x + x^2 sin x = 0: the TM resonances of order one of a
    sphere, at which a closed shell of currents holds the electric dipole wave at its centre without radiating it
    (2.744, 6.117, 9.317, ...)."""

    def g(x: Any) -> Any:
        return x * np.cos(x) - np.sin(x) + x**2 * np.sin(x)

    grid = np.linspace(1.0, top, int(200 * top))
    values = g(grid)
    roots = []
    for i in np.flatnonzero(np.sign(values[:-1]) != np.sign(values[1:])):
        lo, hi = float(grid[i]), float(grid[i + 1])
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if np.sign(g(mid)) == np.sign(g(lo)) else (lo, mid)
        roots.append(0.5 * (lo + hi))
    return np.array(roots)


def compute_ceiling(quick: bool) -> dict[str, Any]:
    say("\n== ceiling: the most |E|^2 per radiated watt at a point, as a share of the ceiling RFG-2 ==")
    kb = np.geomspace(0.8, 60.0, 14 if quick else 120)
    radius = 0.25
    shells: dict[int, dict[str, Arr]] = {}
    for count in (16, 64, 256):
        lossless, loss = [], []
        for x in kb:
            el = shell(float(x) * C0 / (2 * math.pi * radius), count, radius)
            lossless.append(focus_share(el, (0, 0, 0), 0.0))
            loss.append(focus_share(el, (0, 0, 0), LOSS))
        shells[count] = {"lossless": np.array(lossless), "loss": np.array(loss)}
    at_pi = {}
    for count in (16, 64):
        el = shell(math.pi * C0 / (2 * math.pi * radius), count, radius)
        at_pi[count] = (focus_share(el, (0, 0, 0), 0.0), focus_share(el, (0, 0, 0), LOSS))
        say(
            f"  {count} on a shell of radius lambda/2: {at_pi[count][0]:.2f} of the ceiling lossless, "
            f"{at_pi[count][1]:.2f} with 10 % loss; 3N/(16 (kb)^2) = {3 * count / (16 * math.pi**2):.2f}"
        )
    for count in (16, 64, 256):
        lim = math.sqrt(count / 0.1)
        i = int(np.argmin(np.abs(kb - 2 * lim)))
        say(
            f"  N = {count}: the line 3N/(16 (kb)^2) holds while N < 0.05 x 2 (kb)^2, kb > {lim:.1f}; at kb "
            f"{kb[i]:.1f} the model gives {shells[count]['lossless'][i]:.4f} against "
            f"{3 * count / (16 * kb[i] ** 2):.4f}"
        )
    sparse = shell(3.0 * C0 / (2 * math.pi * radius), 16, radius)
    dense = shell(6.1 * C0 / (2 * math.pi * radius), 256, radius)
    passing = {
        "16 at kb = 3": (focus_share(sparse, (0, 0, 0), 0.0), focus_share(sparse, (0, 0, 0), LOSS)),
        "256 at kb = 6.1": (focus_share(dense, (0, 0, 0), 0.0), focus_share(dense, (0, 0, 0), LOSS)),
    }
    beyond = kb >= 3
    first = (float(shells[16]["lossless"][beyond][0]), float(shells[16]["loss"][beyond][0]), float(kb[beyond][0]))
    most16 = (max(first[0], passing["16 at kb = 3"][0]), max(first[1], passing["16 at kb = 3"][1]))
    say(
        f"  16 on a shell at kb >= 3: at most {most16[0]:.2f} of the ceiling lossless and {most16[1]:.2f} with 10 % "
        f"loss, at kb = 3 ({first[0]:.2f} and {first[1]:.2f} by kb = {first[2]:.2f})"
    )
    say(
        f"  256 on a shell at kb = 6.1, beside the cavity resonance at 6.117: {passing['256 at kb = 6.1'][0]:.0f} "
        f"times the ceiling lossless, {passing['256 at kb = 6.1'][1]:.1f} times with 10 % loss"
    )
    resonances = cavity_resonances(float(kb[-1]))
    at_resonance = {}
    for x in resonances[resonances < 16]:
        el = shell(float(x) * C0 / (2 * math.pi * radius), 256, radius)
        at_resonance[float(x)] = focus_share(el, (0, 0, 0), LOSS)
    between = (kb > 2) & (kb < 16)
    say(
        "  256 on a shell with 10 % loss at the cavity resonances: "
        + ", ".join(f"kb {x:.3f} {v:.3g}" for x, v in at_resonance.items())
        + f" times the ceiling, and down to {float(shells[256]['loss'][between].min()):.2f} between them"
    )
    rings = []
    for f in (0.433e9, 0.9e9, 2.3e9):
        el = whips(f)
        rf = (0.0, 0.0, 0.1)
        kr = el.k * float(np.min(np.linalg.norm(el.feeds - np.array(rf), axis=1)))
        rings.append(
            {
                "f_GHz": f / 1e9,
                "kr": kr,
                "loss": focus_share(el, rf, LOSS),
                "loss_free": focus_share(el, rf, LOSS, free_space=True),
            }
        )
        say(
            f"  whips {f / 1e9:5.3f} GHz, 100 mm above the centre (kr {kr:4.2f} from the nearest feed): "
            f"{rings[-1]['loss']:.3g} of the ceiling over the ground, {rings[-1]['loss_free']:.3g} of RFG-2's, "
            "10 % loss"
        )
    heights = np.array([0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5])
    board: dict[str, dict[str, Arr]] = {}
    for label, el, component in (
        ("whips 0.915 GHz", whips(0.915e9), None),
        ("tile 2.44 GHz, E_y", panel(2.44e9), 1),
        ("tile 5.8 GHz, E_y", panel(5.8e9), 1),
    ):
        hs = heights[heights >= 0.03] if label.startswith("whips") else heights
        krs = np.array([el.k * float(np.min(np.linalg.norm(el.feeds - np.array([0, 0, z]), axis=1))) for z in hs])
        share_lossless = np.array([focus_share(el, (0, 0, z), 0.0, component) for z in hs])
        share_loss = np.array([focus_share(el, (0, 0, z), LOSS, component) for z in hs])
        free_lossless = np.array([focus_share(el, (0, 0, z), 0.0, component, free_space=True) for z in hs])
        free_loss = np.array([focus_share(el, (0, 0, z), LOSS, component, free_space=True) for z in hs])
        board[label] = {"height_mm": hs * 1e3, "kr": krs, "lossless": share_lossless, "loss": share_loss}
        board[label] |= {"lossless_free": free_lossless, "loss_free": free_loss}
        say(
            f"  {label:18s} on the axis, "
            + ", ".join(f"{z * 1e3:.0f} mm {v:.3g}" for z, v in zip(hs, share_lossless, strict=True))
            + " lossless, of the ceiling over the ground"
        )
        say(
            f"  {label:18s} on the axis, "
            + ", ".join(f"{z * 1e3:.0f} mm {v:.3g}" for z, v in zip(hs, free_lossless, strict=True))
            + " lossless, of RFG-2's ceiling"
        )
    say(
        "  a closed shell stores the electric dipole wave without radiating it where (x j1(x))' = 0: kb = "
        + ", ".join(f"{x:.3f}" for x in resonances[:6])
    )
    return {
        "kb": kb,
        "radius_m": radius,
        "shells": shells,
        "at_kb_pi": at_pi,
        "passing": passing,
        "dense_at_resonance": at_resonance,
        "rings": rings,
        "board": board,
        "cavity_kb": resonances,
    }


# ---- 4.8: across a band ---------------------------------------------------------------------------------------------
def compute_band(quick: bool) -> dict[str, Any]:
    say("\n== band: 1 - |corr|^2 between the operator at f and at f + df ==")
    offsets = np.array([2e7, 1e8, 5e8]) if quick else np.array([1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]) * 1e6
    whip_cut = C0 / 0.9e9 / 4
    cases: list[tuple[str, Callable[[float], Elements], float, tuple[float, float, float], float]] = [
        ("tile 5.8 GHz", panel, 5.8e9, (0.0, 0.0, 0.15), 0.05),
        ("tile 2.4 GHz", panel, 2.4e9, (0.0, 0.0, 0.15), 0.05),
        ("whips cut for 0.9 GHz", lambda f: whips(f, length=whip_cut, segments=8), 0.9e9, (0.0, 0.0, 0.10), 0.10),
    ]
    rows = {}
    for name, make, f0, centre, radius in cases:
        pts, _ = ball_cells(centre, radius, C0 / f0 / 9)
        g0 = stacked(make(f0), pts)
        out = []
        for df in offsets:
            g1 = stacked(make(f0 + float(df)), pts)
            corr = abs(np.vdot(g0, g1)) / (np.linalg.norm(g0) * np.linalg.norm(g1))
            out.append(1 - corr**2)
        rows[name] = np.array(out)
        say(f"  {name:22s} " + ", ".join(f"{d / 1e6:g} MHz: {v:.2g}" for d, v in zip(offsets, out, strict=True)))
    turn = 360 * DELAYS[-1] * 20e6
    say(f"  rev B's longest delay, {DELAYS[-1] * 1e12:.1f} ps, turns {turn:.1f} degrees across 20 MHz")
    return {"offsets_Hz": offsets, "rows": rows, "longest_delay_turn_deg": turn}


# ---- 5.3: virtual sources behind the board --------------------------------------------------------------------------
TILE_SIDE = 0.1  # the tile's aperture, 4 x 25 mm


def cone_shells(half_angle_deg: float, per_shell: int, radii: Sequence[float] = (0.3, 0.6, 2.0)) -> list[Arr]:
    """Fibonacci points in a cone about boresight on spheres about the tile's centre."""
    ct = math.cos(math.radians(half_angle_deg))
    return [fibonacci(per_shell, ct, 1.0, r) for r in radii]


def shell_system(el: Elements, target: Elements, shells: Sequence[Arr]) -> tuple[CArr, CArr]:
    """The fit's matrix and target with every shell scaled by the target's norm on it, so that each shell counts
    alike although the field falls as 1/r."""
    blocks, rhs = [], []
    for pts in shells:
        t = stacked(target, pts)[:, 0]
        scale = float(np.linalg.norm(t))
        blocks.append(stacked(el, pts) / scale)
        rhs.append(t / scale)
    a: CArr = np.vstack(blocks)
    b: CArr = np.concatenate(rhs)
    return a, b


def closed_form_drive(el: Elements, source: Sequence[float]) -> CArr:
    """The closed-form drive of a virtual point source behind the board: a delay rho/c to each element, amplitude
    1/rho, obliquity cos(theta) and the filter (ik - 1/rho), the Rayleigh integral sampled at the elements."""
    rel = el.feeds - np.asarray(source, dtype=np.float64)
    rho = np.linalg.norm(rel, axis=1)
    cos_theta = rel[:, 2] / rho
    w: CArr = np.asarray(cos_theta * (1j * el.k - 1 / rho) * np.exp(1j * el.k * rho) / rho, dtype=np.complex128)
    return w


def compute_behind(quick: bool) -> dict[str, Any]:
    say("\n== behind: a co-polarised virtual dipole behind the tile, fitted on shells at 0.3, 0.6 and 2 m ==")
    per = 120 if quick else 300
    cones = (20.0, 30.0, 45.0, 60.0, 85.0)
    sources = [(5.8e9, 0.025), (5.8e9, 0.05), (5.8e9, 0.1), (5.8e9, 0.3), (2.4e9, 0.05)]
    rows = {}
    for f, depth in sources:
        el = panel(f)
        src = (0.0, 0.0, -depth)
        target = dipole(f, src, (0.0, 1.0, 0.0))
        fit, closed, sampled, sens = [], [], [], []
        r = radiation(el)
        for cone in cones:
            a, b = shell_system(el, target, cone_shells(cone, per))
            w = loaded_fit(a, b, 1e-3)
            fit.append(residual_db(a, b, w))
            closed.append(scaled_residual_db(a, b, closed_form_drive(el, src)))
            sampled.append(scaled_residual_db(a, b, field(target, el.feeds)[:, 1, 0]))
            uniform = np.ones(el.count) / 4
            sens.append(float(db((np.vdot(w, w).real / np.real(w.conj() @ r @ w)) / (1 / (uniform @ r @ uniform)))))
        window = math.degrees(math.atan(TILE_SIDE / (2 * depth)))
        key = f"{depth * 1e3:.0f} mm behind, {f / 1e9:g} GHz"
        rows[key] = {"f_GHz": f / 1e9, "depth_mm": depth * 1e3, "window_deg": window, "fit_db": np.array(fit)}
        rows[key] |= {"closed_form_db": np.array(closed), "sampled_db": np.array(sampled), "S_db": np.array(sens)}
        say(
            f"  {key:24s} window {window:4.1f} deg; loaded fit "
            + " / ".join(f"{v:.1f}" for v in fit)
            + " dB in cones of "
            + "/".join(f"{c:g}" for c in cones)
            + " deg; closed form "
            + " / ".join(f"{v:.1f}" for v in closed)
            + "; S "
            + " / ".join(f"{v:.1f}" for v in sens)
            + " dB"
        )
    maps = {}
    xs = np.linspace(-0.35, 0.35, 57 if quick else 141)
    zs = np.linspace(0.02, 0.6, 47 if quick else 117)
    gx, gz = np.meshgrid(xs, zs)
    grid = np.c_[gx.ravel(), 0 * gx.ravel(), gz.ravel()]
    depth = 0.05
    along = gx * depth / (gz + depth)  # where the ray from the source to the point crosses the aperture plane
    margin = TILE_SIDE / 2 - np.abs(along)
    s1 = np.hypot(along, depth)
    s2 = np.hypot(gx - along, gz)
    edges = {}
    for f in (5.8e9, 2.4e9):
        el = panel(f)
        target = dipole(f, (0.0, 0.0, -depth), (0.0, 1.0, 0.0))
        a, b = shell_system(el, target, cone_shells(45.0, per))
        w = loaded_fit(a, b, 1e-3)
        made, want = field(el, grid) @ w, field(target, grid)[..., 0]
        err = np.linalg.norm(made - want, axis=1) / np.linalg.norm(want, axis=1)
        key = f"{f / 1e9:g} GHz"
        maps[key] = 20 * np.log10(np.maximum(err, 1e-6)).reshape(gx.shape)
        # the lit side of a Fresnel edge wave, 0.225 / v, v = t sqrt(2 / (lambda L)), L = d r / (d + r)
        v = np.where(margin > 0, margin, np.nan) * np.sqrt(2 / (C0 / f * s1 * s2 / (s1 + s2)))
        edges[key] = 20 * np.log10(0.225 / v)
        inside = np.abs(gx) < TILE_SIDE / 2 * (gz + depth) / depth
        say(
            f"  map at {key}: median error {float(np.median(maps[key][inside & (gz > 0.1)])):.1f} dB inside the "
            f"window beyond 100 mm, {float(np.median(maps[key][~inside & (gz > 0.1)])):.1f} dB outside it"
        )
    return {
        "cones_deg": np.array(cones),
        "rows": rows,
        "map_x": xs,
        "map_z": zs,
        "maps": maps,
        "edge_db": edges,
        "map_depth_mm": depth * 1e3,
    }


# ---- 5.4: focused sources in front ---------------------------------------------------------------------------------
def conjugate_focus(el: Elements, point: Sequence[float]) -> CArr:
    """The unit-norm drive that maximises |E|^2 at a point: the first right singular vector of G there."""
    _, _, vh = np.linalg.svd(field(el, [point])[0])
    w: CArr = vh[0].conj()
    return w


def compute_focused(quick: bool) -> dict[str, Any]:
    say("\n== focused: a virtual source 60 mm in front of the tile, fitted beyond it on a 200 mm square at 0.2 m ==")
    focus = (0.0, 0.0, 0.06)
    rows = {}
    drive = None
    for f in (2.44e9, 3.5e9, 5.8e9):
        el = panel(f)
        target = dipole(f, focus, (0.0, 1.0, 0.0))
        pts = square(0.2, 0.1, el.wavelength / (4 if quick else 6))
        a, b = stacked(el, pts), stacked(target, pts)[:, 0]
        w = loaded_fit(a, b, 1e-6)
        v = a @ w
        w = w * (v.conj() @ b) / (v.conj() @ v)
        rows[f"{f / 1e9:g} GHz"] = residual_db(a, b, w)
        say(f"  {f / 1e9:4.2f} GHz: reproduced to {rows[f'{f / 1e9:g} GHz']:.1f} dB over the square")
        if f == 5.8e9:
            drive = w
    assert drive is not None
    el = panel(5.8e9)
    target = dipole(5.8e9, focus, (0.0, 1.0, 0.0))
    xs = np.linspace(-0.15, 0.15, 60 if quick else 180)  # an even count keeps the grid off the source itself
    zs = np.linspace(0.02, 0.32, 60 if quick else 180)
    gx, gz = np.meshgrid(xs, zs)
    grid = np.c_[gx.ravel(), 0 * gx.ravel(), gz.ravel()]
    made = field(el, grid) @ drive
    want = field(target, grid)[..., 0]
    at_focus = float(np.linalg.norm(field(el, [focus])[0] @ drive))
    intensity = db((np.abs(made) ** 2).sum(1) / at_focus**2).reshape(gx.shape)
    phase = np.angle(made[:, 1]).reshape(gx.shape)
    err = 20 * np.log10(np.maximum(np.linalg.norm(made - want, axis=1) / np.linalg.norm(want, axis=1), 1e-6))
    cone = math.degrees(math.atan(TILE_SIDE / (2 * focus[2])))
    beyond = (gz > focus[2]) & (np.abs(gx) < (gz - focus[2]) * TILE_SIDE / (2 * focus[2]))
    regions = {
        "inside the cone, beyond 100 mm": beyond & (gz > 0.1),
        "outside the cone, beyond 100 mm": ~beyond & (gz > 0.1),
        "before the focus": gz < focus[2],
    }
    median_error = {name: float(np.median(err.reshape(gx.shape)[mask])) for name, mask in regions.items()}
    say(
        "  error against a source at the focus, median: "
        + ", ".join(f"{k} {v:.1f} dB" for k, v in median_error.items())
    )
    axis_z = np.linspace(0.016, 0.6, 600 if quick else 1200)
    ray = np.c_[0 * axis_z, 0 * axis_z, axis_z]
    profiles: dict[str, dict[str, Any]] = {}
    for f, asked in ((5.8e9, (0.03, 0.05, 0.1, 0.2)), (2.44e9, (0.05, 0.1))):
        tile = panel(f)
        g_axis = field(tile, ray)
        for zf in asked:
            w = conjugate_focus(tile, (0.0, 0.0, zf))
            prof = (np.abs(g_axis @ w) ** 2).sum(1)
            inner = np.flatnonzero((prof[1:-1] > prof[:-2]) & (prof[1:-1] > prof[2:])) + 1
            peak = float(axis_z[inner[np.argmax(prof[inner])]]) if len(inner) else None
            profiles[f"{f / 1e9:g} GHz, {zf * 1e3:.0f} mm"] = {"f_GHz": f / 1e9, "asked_mm": zf * 1e3}
            profiles[f"{f / 1e9:g} GHz, {zf * 1e3:.0f} mm"] |= {
                "peak_mm": None if peak is None else peak * 1e3,
                "profile": prof / prof.max(),
            }
            say(
                f"  conjugate focus asked at {zf * 1e3:3.0f} mm, {f / 1e9:4.2f} GHz: on-axis maximum "
                + ("none" if peak is None else f"at {peak * 1e3:.1f} mm")
            )
    return {
        "residual_db": rows,
        "focus_mm": [c * 1e3 for c in focus],
        "cone_deg": cone,
        "median_error_db": median_error,
        "map_x": xs,
        "map_z": zs,
        "intensity_db": intensity,
        "phase": phase,
        "error_db": err.reshape(gx.shape),
        "axis_mm": axis_z * 1e3,
        "profiles": profiles,
    }


# ---- 5.5: the whip rings ---------------------------------------------------------------------------------------------
def band_points(n: int, radius: float, max_elevation_deg: float = 45.0) -> Arr:
    """Points over all azimuths from the horizon up to an elevation, on a sphere about the plate's centre."""
    return fibonacci(n, 0.0, math.sin(math.radians(max_elevation_deg)), radius, ascending=True)


def compute_whips(quick: bool) -> dict[str, Any]:
    say("\n== whips: the rings at 0.9 GHz as an interior array ==")
    f = 0.9e9
    lam = C0 / f
    configs = {"both rings": None, "outer ring": tuple(range(8)), "inner ring": tuple(range(8, 16))}
    arrays = {name: whips(f, select=sel) for name, sel in configs.items()}
    far = band_points(500 if quick else 1500, 6.0)
    centre_whip = monopoles(f, [(0.0, 0.0)], lam / 4, 4)
    e0 = field(centre_whip, far)[..., 0]
    phi = np.arctan2(far[:, 1], far[:, 0])
    orders = np.arange(9)
    patterns: dict[str, dict[str, Arr]] = {}
    for name, el in arrays.items():
        a = stacked(el, far)
        r = radiation(el)
        res, sens = [], []
        for m in orders:
            b = (e0 * np.exp(1j * m * phi)[:, None]).reshape(-1)
            w = loaded_fit(a, b, 1e-4)
            res.append(residual_db(a, b, w))
            sens.append(float(db(np.vdot(w, w).real * r[0, 0] / np.real(w.conj() @ r @ w))))
        patterns[name] = {"residual_db": np.array(res), "S_db": np.array(sens)}
        say(
            f"  {name:10s} exp(i m phi), m = 0..8: "
            + ", ".join(f"{v:.1f}" for v in res)
            + " dB; S "
            + ", ".join(f"{v:.1f}" for v in sens)
            + " dB"
        )
    shells = [band_points(200 if quick else 800, rr) for rr in (1.0, 2.0, 6.0)]
    both = arrays["both rings"]
    blocks = [stacked(both, pts) for pts in shells]
    step = 0.1 if quick else 0.025
    xs = np.arange(-0.35, 0.35 + 1e-9, step)
    residual = np.empty((len(xs), len(xs)))
    for i, y in enumerate(xs):
        for j, x in enumerate(xs):
            virtual = monopoles(f, [(float(x), float(y))], lam / 4, 4)
            rows, rhs = [], []
            for block, pts in zip(blocks, shells, strict=True):
                t = stacked(virtual, pts)[:, 0]
                scale = float(np.linalg.norm(t))
                rows.append(block / scale)
                rhs.append(t / scale)
            a, b = np.vstack(rows), np.concatenate(rhs)
            residual[i, j] = residual_db(a, b, loaded_fit(a, b, 1e-3))
    radius = np.hypot(*np.meshgrid(xs, xs))
    worst = {f"{r * 1e3:.0f} mm": float(np.max(residual[radius <= r + 1e-9])) for r in (0.1, 0.15, 0.2)}
    say(
        "  a virtual whip anywhere within a radius of the centre, worst residual: "
        + ", ".join(f"{k} {v:.1f} dB" for k, v in worst.items())
    )
    on_line = {}
    for x in (0.0, 0.05, 0.1, 0.15, 0.2, 0.3):
        virtual = monopoles(f, [(x, 0.0)], lam / 4, 4)
        rows, rhs = [], []
        for block, pts in zip(blocks, shells, strict=True):
            t = stacked(virtual, pts)[:, 0]
            rows.append(block / float(np.linalg.norm(t)))
            rhs.append(t / float(np.linalg.norm(t)))
        a, b = np.vstack(rows), np.concatenate(rhs)
        w = loaded_fit(a, b, 1e-3)
        r = radiation(both)
        on_line[f"{x * 1e3:.0f} mm"] = {"residual_db": residual_db(a, b, w)}
        on_line[f"{x * 1e3:.0f} mm"]["S_db"] = float(db(np.vdot(w, w).real * r[0, 0] / np.real(w.conj() @ r @ w)))
    say(
        "  a virtual whip at (x, 0): "
        + ", ".join(f"{k} {v['residual_db']:.1f} dB (S {v['S_db']:.1f})" for k, v in on_line.items())
    )
    elevation = np.linspace(0.0, 45.0, 91)
    bessel = {}
    for name, radius in (("outer ring, 225 mm", 0.225), ("inner ring, 112 mm", 0.112)):
        ka = 2 * math.pi / lam * radius
        bessel[name] = {
            "ka": ka,
            "J": np.array([np.abs(bessel_j(m, ka * np.cos(np.radians(elevation)))) for m in range(9)]),
        }
    zero = math.degrees(math.acos(3.8317 / (2 * math.pi / lam * 0.225)))
    say(f"  J_1(ka cos el) of the outer ring (ka {2 * math.pi / lam * 0.225:.2f}) changes sign at {zero:.1f} deg")
    return {
        "orders": orders,
        "patterns": patterns,
        "map_x": xs,
        "map_residual_db": residual,
        "worst_within_db": worst,
        "on_line": on_line,
        "elevation_deg": elevation,
        "bessel": bessel,
        "places": whip_places(),
        "J1_zero_deg": zero,
    }


# ---- 6.2: the reachable weights of an alphabet ----------------------------------------------------------------------
def alphabet_phases(board: str, f: float) -> Arr:
    """Rev A's 2-bit delay lines (quarter turns at 5.8 GHz, growing with f) or rev B's sixteen delays, radians."""
    if board == "rev A":
        out: Arr = np.arange(4, dtype=np.float64) * (math.pi / 2) * (f / 5.8e9)
        return out
    out = 2 * math.pi * f * DELAYS
    return out


def hull_numbers(phases: Arr) -> dict[str, float]:
    """For unit phasors at the given phases: the largest gap, the line-1 ceiling Perimeter(conv A)/2 pi
    (sum of sin(gap/2) over pi) and the line-0 inradius cos(largest gap/2) (TMA-3, TMA-4)."""
    p = np.sort(np.mod(phases, 2 * math.pi))
    gaps = np.diff(np.r_[p, p[0] + 2 * math.pi])
    gaps = gaps[gaps > 1e-12]
    largest = float(gaps.max())
    return {
        "largest_gap_deg": math.degrees(largest),
        "line1": float(np.sum(np.sin(gaps / 2)) / math.pi),
        "inradius": math.cos(largest / 2) if largest < math.pi else 0.0,
    }


def compute_alphabets(quick: bool) -> dict[str, Any]:
    say("\n== alphabets: the hull (line 0) and the line-1 disc of each alphabet ==")
    cases = [("rev A", 5.8e9), ("rev B", 2.3e9), ("rev B", 3.5e9), ("rev B", 5.8e9)]
    panels = []
    for board, f in cases:
        ph = alphabet_phases(board, f)
        nums = hull_numbers(ph)
        panels.append({"board": board, "f_GHz": f / 1e9, "phases": ph} | nums)
        say(
            f"  {board} at {f / 1e9:g} GHz: largest gap {nums['largest_gap_deg']:.1f} deg, line-1 ceiling "
            f"{nums['line1']:.4f}, line-0 inradius {nums['inradius']:.3f}"
        )
    for board, f in (("rev A", 5.0e9), ("rev A", 6.0e9), ("rev B", 2.44e9)):
        nums = hull_numbers(alphabet_phases(board, f))
        say(f"  {board} at {f / 1e9:g} GHz: " + ", ".join(f"{k} {v:.4g}" for k, v in nums.items()))
    fb = np.linspace(2.3e9, 6.0e9, 371 if quick else 1481)
    sweep_b = [hull_numbers(alphabet_phases("rev B", float(f))) for f in fb]
    sweep_bp = [
        hull_numbers(np.r_[alphabet_phases("rev B", float(f)), alphabet_phases("rev B", float(f)) + math.pi])
        for f in fb
    ]
    fa = np.linspace(5.0e9, 6.0e9, 101)
    sweep_a = [hull_numbers(alphabet_phases("rev A", float(f))) for f in fa]
    line1 = np.array([s["line1"] for s in sweep_b])
    say(
        f"  rev B 2.3-6 GHz: line-1 ceiling {line1.min():.4f} to {line1.max():.4f}, least at "
        f"{fb[np.argmin(line1)] / 1e9:.2f} GHz; inradius {min(s['inradius'] for s in sweep_b):.2f} to "
        f"{max(s['inradius'] for s in sweep_b):.2f}; with polarity {min(s['line1'] for s in sweep_bp):.3f} to "
        f"{max(s['line1'] for s in sweep_bp):.3f}"
    )
    return {
        "panels": panels,
        "rev_b_GHz": fb / 1e9,
        "rev_b_line1": line1,
        "rev_b_inradius": np.array([s["inradius"] for s in sweep_b]),
        "rev_b_polarity_line1": np.array([s["line1"] for s in sweep_bp]),
        "rev_a_GHz": fa / 1e9,
        "rev_a_line1": np.array([s["line1"] for s in sweep_a]),
        "rev_a_inradius": np.array([s["inradius"] for s in sweep_a]),
    }


# ---- 6.3: sequences and their comb lines -----------------------------------------------------------------------------
def line_weights(g: CArr, lines: Any) -> CArr:
    """The harmonic weight law (TMA-1): element n holding g[n, l] in slot l of L has on comb line h the weight
    sinc(h/L) exp(-i pi h/L) (1/L) sum_l g[n, l] exp(-2 pi i h l / L). Returns (lines, N)."""
    hs = np.atleast_1d(np.asarray(lines, dtype=np.float64))
    length = g.shape[1]
    slots = np.arange(length)
    kernel = np.exp(-2j * math.pi * np.outer(hs, slots) / length) / length
    c: CArr = (np.sinc(hs / length) * np.exp(-1j * math.pi * hs / length))[:, None] * (kernel @ g.T)
    return c


def rev_a_positions(pitch: float = 0.026) -> Arr:
    """The rev A panel's element centres as the time-modulation studies place them (element i: column i mod 4)."""
    i = np.arange(16)
    out: Arr = np.c_[(i % 4 - 1.5) * pitch, (1.5 - i // 4) * pitch, np.zeros(16)]
    return out


def steering(positions: Arr, f: float, theta: float, phi: float) -> CArr:
    """exp(i k r.u): the far-field response of isotropic elements toward (theta, phi)."""
    k = 2 * math.pi * f / C0
    u = np.array([math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)])
    a: CArr = np.exp(1j * k * positions @ u)
    return a


def compute_sequences(quick: bool) -> dict[str, Any]:
    say("\n== sequences: what each schedule puts on each comb line (L = 16 slots, no marker) ==")
    length, count = 16, 16
    quarter = np.array([1, 1j, -1, -1j], dtype=np.complex128)
    pos = rev_a_positions()
    a = steering(pos, 5.8e9, math.radians(20.0), 0.0)
    shifts = np.round(-np.angle(np.conj(a)) / (2 * math.pi) * length).astype(int) % length
    stair = np.array([np.roll(np.repeat(np.arange(4), 4), s) for s in shifts])
    frank_base = np.array([((slot // 4) * (slot % 4)) % 4 for slot in range(length)])
    frank = np.array([np.roll(frank_base, n) for n in range(count)])
    hadamard = np.where(np.asarray(H16).T > 0, 0, 2)
    lines = np.arange(-8, 9)
    out: dict[str, Any] = {"lines": lines}
    for name, states in (("staircase", stair), ("Frank shifts", frank), ("Hadamard words", hadamard)):
        g = quarter[states]
        c = line_weights(g, lines)
        power = np.mean(np.abs(c) ** 2, axis=1)
        beams = np.sum(np.abs(c), axis=1) ** 2 / count**2
        out[name] = {"states": states, "weights": c, "power_per_line": power, "best_beam_per_line": beams}
        shown = np.maximum(db(power), -99.0)
        say(
            f"  {name:15s} per-element power on lines -8..8 (dB, floored at -99): "
            + " ".join(f"{v:5.1f}" for v in shown)
        )
    c1 = line_weights(quarter[stair], [1])[0]
    beam = float(db(abs(c1 @ a) ** 2 / count**2))
    frank_lines = out["Frank shifts"]["best_beam_per_line"][:16]
    say(
        f"  staircase steered to 20 deg: line-1 beam {beam:.2f} dB of 16 ideal elements; Frank lines carry "
        f"{db(frank_lines.max()):.1f} to {db(frank_lines.min()):.1f} dB of the full-array gain, "
        f"{frank_lines.sum():.3f} in all"
    )
    sv = np.linalg.svd(line_weights(quarter[hadamard], np.arange(-8, 8)), compute_uv=False)
    say(f"  Hadamard words on lines -8..7: decode matrix condition number {sv[0] / sv[-1]:.3f}")
    out["staircase_beam_db"] = beam
    out["hadamard_condition"] = float(sv[0] / sv[-1])
    return out


# ---- 6.4: recovering every element's signal from one channel ---------------------------------------------------------
def recovery_errors(
    bandwidth: float, rng: np.random.Generator, periods: int, per_slot: int = 8, slot: float = 8e-6, snr_db: float = 40
) -> tuple[float, float]:
    """One simulated capture: sixteen elements under the Hadamard words behind a blank marker (17 slots), a random
    array vector, a signal of the given two-sided bandwidth (0 for a tone) and white noise. Returns the error of the
    array vector (1 - |cos|^2 of its direction, dB) decoded by slot means and by comb lines."""
    words = np.vstack([np.zeros((1, 16)), np.asarray(H16, dtype=np.float64)])
    slots = len(words)
    period = slots * per_slot
    n = period * periods
    fs = per_slot / slot
    a = (rng.normal(size=16) + 1j * rng.normal(size=16)) / math.sqrt(2)
    if bandwidth > 0:
        spectrum = np.zeros(n, dtype=np.complex128)
        band = np.abs(np.fft.fftfreq(n, 1 / fs)) <= bandwidth / 2
        spectrum[band] = rng.normal(size=band.sum()) + 1j * rng.normal(size=band.sum())
        s = np.fft.ifft(spectrum)
        s = s / math.sqrt(float(np.mean(np.abs(s) ** 2)))
    else:
        s = np.ones(n, dtype=np.complex128)
    g = np.tile(np.repeat(words, per_slot, axis=0), (periods, 1))
    noise = (rng.normal(size=n) + 1j * rng.normal(size=n)) * 10 ** (-snr_db / 20) / math.sqrt(2) * 4
    y = (g * a[None, :]).sum(1) * s + noise

    def error(x: Any) -> float:
        _, vec = np.linalg.eigh(x.conj().T @ x)
        u = vec[:, -1].conj()  # the principal direction of x^H x is the conjugate of the array vector
        cos2 = abs(np.vdot(a, u)) ** 2 / (np.vdot(a, a).real * np.vdot(u, u).real)
        return float(db(1 - cos2))

    segments = y.reshape(periods, slots, 1, per_slot)
    starts = (np.arange(periods)[:, None] * slots + np.arange(slots)[None, :]) * per_slot
    means = slot_means(segments, starts, fs, 0.0)[..., 0]
    by_slots = means @ np.linalg.pinv(words).T
    hs = np.arange(-24, 25)
    c = np.stack([(g[:period] * np.exp(-2j * math.pi * h * np.arange(period) / period)[:, None]).mean(0) for h in hs])
    yf = np.fft.fft(y)
    keep = periods // 2
    z = np.array([yf[(h * periods + np.arange(-keep + 1, keep)) % n] for h in hs])
    by_lines = (np.linalg.pinv(c) @ z).T
    return error(by_slots), error(by_lines)


def compute_lines(quick: bool) -> dict[str, Any]:
    say("\n== lines: sixteen elements, Hadamard words and a marker, 8 us slots (f_p 7.35 kHz), 40 dB SNR ==")
    rng = np.random.default_rng(7)
    bands = np.array([0, 1000, 6000]) if quick else np.array([0, 30, 100, 300, 1000, 2000, 3000, 4500, 6000, 7000])
    draws, periods = (1, 256) if quick else (3, 1024)
    slots, lines = [], []
    for bw in bands:
        trial = np.array([recovery_errors(float(bw), rng, periods) for _ in range(draws)])
        slots.append(float(np.median(trial[:, 0])))
        lines.append(float(np.median(trial[:, 1])))
        say(f"  B {bw:5.0f} Hz: slot means {slots[-1]:6.1f} dB, comb lines {lines[-1]:6.1f} dB")
    comb = {
        f"{board} L={length}": 1 / (length * tau)
        for board, tau in (("rev A", 10e-6), ("rev B", 15e-6))
        for length in (4, 8, 16)
    }
    return {"bandwidth_Hz": bands, "slot_means_db": np.array(slots), "lines_db": np.array(lines)} | {
        "f_p_Hz": comb,
        "f_p_simulated_Hz": 1 / (17 * 8e-6),
    }


# ---- 6.5: what limits a line -----------------------------------------------------------------------------------------
def line_kernel(length: int) -> CArr:
    """The line-1 weight of each slot: sinc(1/L) exp(-i pi/L) exp(-2 pi i l/L) / L."""
    slots = np.arange(length)
    k: CArr = np.sinc(1 / length) * np.exp(-1j * math.pi / length) * np.exp(-2j * math.pi * slots / length) / length
    return k


def synthesis_feedback(target: CArr, alphabet: CArr, length: int, sweeps: int = 3) -> CArr:
    """Each element's sequence over the alphabet whose line-1 weight approaches its target: error feedback over
    the slots from four starting slots, then coordinate passes; the best start is kept. Returns (N, L) gains."""
    ker = line_kernel(length)
    count = len(target)
    rows = np.arange(count)
    best = np.zeros((count, length), dtype=np.intp)
    best_err = np.full(count, np.inf)
    for start in range(0, length, max(1, length // 4)):
        g = np.zeros((count, length), dtype=np.intp)
        acc = np.zeros(count, dtype=np.complex128)
        for j in range(length):
            slot = (start + j) % length
            rem = target * (j + 1) / length - acc
            cand = alphabet * ker[slot]
            s = np.argmin(np.abs(rem[:, None] - cand[None, :]), axis=1)
            g[:, slot] = s
            acc = acc + cand[s]
        for _ in range(sweeps):
            for slot in range(length):
                base = acc - alphabet[g[:, slot]] * ker[slot]
                cand2 = base[:, None] + alphabet[None, :] * ker[slot]
                s = np.argmin(np.abs(cand2 - target[:, None]), axis=1)
                g[:, slot] = s
                acc = cand2[rows, s]
        err = np.abs(acc - target)
        better = err < best_err
        best[better] = g[better]
        best_err = np.where(better, err, best_err)
    out: CArr = alphabet[best]
    return out


def synthesis_search(
    target: CArr, alphabet: CArr, length: int, rng: np.random.Generator, trials: int = 6, sweeps: int = 30
) -> CArr:
    """The same aim by search: from the rounded ramp and from random sequences, coordinate passes in random slot
    order until no slot changes; the best of the trials is kept. Returns (N, L) gains."""
    ker = line_kernel(length)
    count = len(target)
    rows = np.arange(count)
    slots = np.arange(length)
    best = np.zeros((count, length), dtype=np.intp)
    best_err = np.full(count, np.inf)
    for trial in range(trials):
        if trial == 0:
            ramp = (target / np.abs(target))[:, None] * np.exp(2j * math.pi * (slots + 0.5) / length)[None, :]
            g = np.argmax(np.real(np.conj(alphabet)[None, None, :] * ramp[:, :, None]), axis=2)
        else:
            g = rng.integers(len(alphabet), size=(count, length))
        acc = (alphabet[g] * ker[None, :]).sum(1)
        for _ in range(sweeps):
            changed = np.zeros(count, dtype=bool)
            for slot in rng.permutation(length):
                cand = (acc - alphabet[g[:, slot]] * ker[slot])[:, None] + alphabet[None, :] * ker[slot]
                s = np.argmin(np.abs(cand - target[:, None]), axis=1)
                changed |= s != g[:, slot]
                g[:, slot] = s
                acc = cand[rows, s]
            if not changed.any():
                break
        err = np.abs(acc - target)
        better = err < best_err
        best[better] = g[better]
        best_err = np.where(better, err, best_err)
    out: CArr = alphabet[best]
    return out


def lcmv(a0: CArr, a1: CArr) -> CArr:
    """The least-norm weights (applied as sum w_n a_n) with unit response toward a0 and none toward a1."""
    cm = np.stack([a0, a1])
    w: CArr = cm.conj().T @ np.linalg.solve(cm @ cm.conj().T, np.array([1.0, 0.0]))
    return w


def static_null(a0: CArr, a1: CArr, alphabet: CArr, passes: int = 6) -> CArr:
    """One static state per element: from the beam's nearest states, coordinate passes that deepen the null while
    the beam keeps at least 0.4 of sixteen ideal elements' gain."""
    q = alphabet[np.argmin(np.abs(np.conj(a0)[:, None] - alphabet[None, :]), axis=1)].copy()

    def cost(v: CArr) -> float:
        main = float(abs(v @ a0) ** 2)
        return float(abs(v @ a1) ** 2 / max(main, 1e-12) + (0 if main > 0.5 * 256 * 0.81 else 10))

    current = cost(q)
    for _ in range(passes):
        for n in range(len(q)):
            for s in alphabet:
                old = q[n]
                q[n] = s
                c = cost(q)
                if c < current:
                    current = c
                else:
                    q[n] = old
    return q


def compute_nulls(quick: bool) -> dict[str, Any]:
    say("\n== nulls: a line-1 beam with a null, rev A's panel at 5.8 GHz, sequences of L slots ==")
    rng = np.random.default_rng(11)
    pos = rev_a_positions()
    f = 5.8e9
    pairs: list[tuple[CArr, CArr]] = []
    while len(pairs) < (6 if quick else 40):
        a0 = steering(pos, f, math.acos(rng.uniform(math.cos(math.pi / 4), 1)), rng.uniform(0, 2 * math.pi))
        a1 = steering(pos, f, math.acos(rng.uniform(math.cos(math.pi / 3), 1)), rng.uniform(0, 2 * math.pi))
        if abs(np.vdot(a0, a1)) / 16 < 0.5:
            pairs.append((a0, a1))
    lengths = [16, 64] if quick else [4, 8, 16, 32, 64, 128, 256]
    curves: dict[str, list[float]] = {}
    gains: dict[str, list[float]] = {}
    static: dict[str, float] = {}
    for bits in (2, 4):
        alphabet = np.exp(2j * math.pi * np.arange(2**bits) / 2**bits)
        depth = [abs(q @ a1) ** 2 / abs(q @ a0) ** 2 for a0, a1 in pairs for q in [static_null(a0, a1, alphabet)]]
        static[f"{bits}-bit"] = float(db(np.median(depth)))
        for method in ("feedback", "search"):
            key = f"{bits}-bit {method}"
            curves[key], gains[key] = [], []
            for length in lengths:
                d, gn = [], []
                for a0, a1 in pairs:
                    w = lcmv(a0, a1)
                    target = 0.97 * np.sinc(1 / 2**bits) * w / np.abs(w).max()
                    if method == "feedback":
                        g = synthesis_feedback(target, alphabet, length)
                    else:
                        g = synthesis_search(target, alphabet, length, rng)
                    c = (g * line_kernel(length)[None, :]).sum(1)
                    d.append(abs(c @ a1) ** 2 / abs(c @ a0) ** 2)
                    gn.append(abs(c @ a0) ** 2 / 256)
                curves[key].append(float(db(np.median(d))))
                gains[key].append(float(db(np.mean(gn))))
            say(
                f"  {key:16s} median null "
                + ", ".join(f"L {n}: {v:.1f}" for n, v in zip(lengths, curves[key], strict=True))
                + f" dB; beam {np.mean(gains[key]):.2f} dB of 16 ideal elements"
            )
        say(f"  {bits}-bit static states: median null {static[f'{bits}-bit']:.1f} dB")
    floors = {f"{s:g} deg": float(db(math.radians(s) ** 2 / 16)) for s in (1.0, 3.0, 11.4)}
    say("  calibration floor sigma^2/N: " + ", ".join(f"{k} {v:.1f} dB" for k, v in floors.items()))
    a0 = steering(pos, f, 0.0, 0.0)
    a1 = steering(pos, f, math.radians(30.0), 0.0)
    w = lcmv(a0, a1)
    alphabet = np.array([1, 1j, -1, -1j], dtype=np.complex128)
    target = 0.97 * np.sinc(0.25) * w / np.abs(w).max()
    realised = {
        "target": target,
        "one static state": static_null(a0, a1, alphabet),
        "L 16, feedback": (synthesis_feedback(target, alphabet, 16) * line_kernel(16)[None, :]).sum(1),
        "L 64, search": (synthesis_search(target, alphabet, 64, rng) * line_kernel(64)[None, :]).sum(1),
    }
    theta = np.radians(np.linspace(-90, 90, 721))
    cut = np.stack([steering(pos, f, float(t), 0.0) for t in theta])
    pattern = {name: db(np.abs(cut @ c) ** 2 / abs(c @ a0) ** 2) for name, c in realised.items()}
    depth_at = {name: float(db(abs(c @ a1) ** 2 / abs(c @ a0) ** 2)) for name, c in realised.items()}
    say("  broadside beam, null at 30 deg: " + ", ".join(f"{k} {v:.1f} dB" for k, v in depth_at.items()))
    return {"lengths": lengths, "curves": curves, "gains": gains, "static": static, "floors": floors} | {
        "theta_deg": np.degrees(theta),
        "pattern": pattern,
        "pattern_null_db": depth_at,
        "pairs": len(pairs),
    }


# ---- 7.1: the tile's drive space -------------------------------------------------------------------------------------
def darkest_states(r: Arr, phases: Arr, rng: np.random.Generator, start: IntArr, starts: int) -> tuple[float, IntArr]:
    """Unit-amplitude drives over a phase alphabet that radiate least: coordinate passes from a given start and from
    random ones. With |w_n| = 1 the power moves with w_n as 2 Re(conj(w_n) b_n), b_n = sum over m != n of R_nm w_m,
    so each pass picks per element the state that minimises that term."""
    best: tuple[float, IntArr] = (math.inf, start)
    states = np.exp(1j * phases)
    for trial in range(starts):
        idx = start.copy() if trial == 0 else rng.integers(0, len(phases), len(start))
        w = states[idx]
        for _ in range(30):
            changed = False
            for n in range(len(w)):
                b = r[n] @ w - r[n, n] * w[n]
                s = int(np.argmin(np.real(np.conj(states) * b)))
                changed |= s != idx[n]
                idx[n] = s
                w[n] = states[s]
            if not changed:
                break
        power = float(np.real(w.conj() @ r @ w))
        if power < best[0]:
            best = (power, idx.copy())
    return best


def continuous_darkest(r: Arr, rng: np.random.Generator, starts: int, steps: int, seeds: Sequence[Arr] = ()) -> float:
    """The least power of a unit-amplitude drive with continuous phases, by gradient descent from random starts and
    from the given phase vectors (so that it never ends above the discrete drives it is compared with)."""
    ph = np.vstack(
        [
            *(np.asarray(seed, dtype=np.float64)[None, :] for seed in seeds),
            rng.uniform(0, 2 * math.pi, (starts, len(r))),
        ]
    )
    for _ in range(steps):
        w = np.exp(1j * ph)
        ph = ph - 0.05 * 2 * np.imag(np.conj(w) * (w @ r.T))
    w = np.exp(1j * ph)
    return float(np.min(np.real(np.einsum("sn,nm,sm->s", w.conj(), r, w))))


def compute_dark(quick: bool) -> dict[str, Any]:
    say("\n== dark: radiating and reactive drives of the tile, and how dark rev B's delay states make it ==")
    rng = np.random.default_rng(3)
    freqs = np.geomspace(0.4e9, 6.0e9, 25 if quick else 97)
    over = np.array([radiation_eigen(panel(float(f))) for f in freqs])
    free = np.array([radiation_eigen(panel(float(f), ground=False)) for f in freqs])
    signs = np.array(list(itertools.product((1.0, -1.0), repeat=15)))
    signs = np.c_[np.ones(len(signs)), signs]
    table = {}
    keep: dict[str, Any] = {}
    for model, ground in (("over ground", True), ("free space", False)):
        for f in (2.44e9, 3.5e9, 5.8e9):
            r = radiation(panel(f, ground=ground))
            r = r / float(np.mean(np.diag(r)))
            ev = np.linalg.eigvalsh(r)[::-1]
            uniform = float(np.sum(r))
            power = np.einsum("pi,ij,pj->p", signs, r, signs)
            pattern = signs[int(np.argmin(power))]
            ph = 2 * math.pi * f * DELAYS
            gap = np.angle(np.exp(1j * (ph[None, :] - ph[:, None])))
            pa, pb = np.unravel_index(int(np.argmin(np.abs(np.abs(gap) - math.pi))), gap.shape)
            start = np.where(pattern > 0, pa, pb).astype(np.intp)
            dark, idx = darkest_states(r, ph, rng, start, 6 if quick else 40)
            seeds = (np.where(pattern > 0, 0.0, math.pi), ph[idx])
            cont = continuous_darkest(r, rng, 6 if quick else 30, 600 if quick else 3000, seeds)
            key = f"{model} {f / 1e9:g} GHz"
            table[key] = {
                "above_10pct": int((ev > 0.1 * ev[0]).sum()),
                "below_1pct": int((ev < 0.01 * ev[0]).sum()),
                "smallest_over_mean": float(ev[-1]),
                "pair_deg": math.degrees(abs(float(gap[pa, pb]))),
                "sign_pattern_db": float(db(power.min() / uniform)),
                "delay_states_db": float(db(dark / uniform)),
                "continuous_db": float(db(cont / uniform)),
            }
            if ground and f == 2.44e9:
                keep = {"pattern": pattern, "states": idx, "phases": ph}
            say(f"  {key:22s} " + ", ".join(f"{k} {v:.3g}" for k, v in table[key].items()))
    for f in (2.4e9, 5.8e9):
        r = radiation(panel(f, ground=False))
        ev = np.linalg.eigvalsh(r / float(np.mean(np.diag(r))))
        say(f"  free space {f / 1e9:g} GHz: smallest eigenvalue {ev[0]:.2g} of the mean")
    el = panel(2.44e9)
    xs = np.linspace(-0.1, 0.1, 81 if quick else 201)
    zs = np.linspace(0.002, 0.15, 60 if quick else 149)
    gx, gz = np.meshgrid(xs, zs)
    g = field(el, np.c_[gx.ravel(), 0 * gx.ravel(), gz.ravel()])
    drives = {
        "uniform drive": np.ones(16, dtype=np.complex128),
        "darkest delay states": np.exp(1j * keep["phases"][keep["states"]]),
    }
    cuts = {name: (np.abs(g @ w) ** 2).sum(1).reshape(gx.shape) for name, w in drives.items()}
    ref = float(cuts["uniform drive"].max())
    heights = np.geomspace(1e-3, 0.15, 60)
    above = np.c_[np.full(60, 0.0125), np.full(60, 0.0125), 0.013 + heights]
    ga = field(el, above)
    drives["best sign pattern"] = keep["pattern"].astype(np.complex128)
    column, row = np.arange(16) % 4, np.arange(16) // 4
    drives["checkerboard"] = ((-1.0) ** (column + row)).astype(np.complex128)
    kappa = math.sqrt(2 * math.pi**2 / 0.025**2 - el.k**2)
    profiles = {name: (np.abs(ga @ w) ** 2).sum(1) for name, w in drives.items()}
    r = radiation(el)
    radiated = {
        name: round(float(db(np.real(w.conj() @ r @ w) / np.real(np.ones(16) @ r @ np.ones(16)))), 6) + 0.0
        for name, w in drives.items()
    }
    say(
        "  at 2.44 GHz over ground, radiated power against uniform drive: "
        + ", ".join(f"{k} {v:.1f} dB" for k, v in radiated.items())
    )
    ratio = profiles["darkest delay states"] / profiles["uniform drive"]
    checker = profiles["checkerboard"] / profiles["uniform drive"]
    for hmm in (2, 5, 10, 20, 50):
        i = int(np.argmin(np.abs(heights - hmm * 1e-3)))
        say(
            f"    {hmm:2d} mm above an element: |E|^2 of the darkest drive {db(ratio[i]):+.1f} dB, of the checkerboard "
            f"{db(checker[i]):+.1f} dB against uniform drive"
        )
    say(
        f"  the checkerboard's carpet decays as exp(-kappa h), kappa = sqrt(2 pi^2 / p^2 - k^2) = {kappa:.0f} per metre"
    )
    return {
        "freqs_GHz": freqs / 1e9,
        "kappa_per_m": kappa,
        "over_db": db(over),
        "free_db": db(free),
        "table": table,
        "cut_x": xs,
        "cut_z": zs,
        "cuts_db": {name: db(c / ref) for name, c in cuts.items()},
        "heights_mm": heights * 1e3,
        "profiles": profiles,
        "radiated_db": radiated,
        "states": keep["states"],
    }


# ---- 7.2: the guide as an antenna with a dial ------------------------------------------------------------------------
def order_power(order: int, x: Any) -> Arr:
    """F_l(x) = (3/4) integral over [0, pi] of sin^3(t) J_l(x sin t)^2: the share of a ring's radiation that goes
    into orbital order l, the rings' sum over every l being 1 (NFD-5)."""
    xs = np.atleast_1d(np.asarray(x, dtype=np.float64))
    t, w = np.polynomial.legendre.leggauss(96)
    theta = 0.5 * math.pi * (t + 1)
    weight = 0.5 * math.pi * w
    j = np.array([bessel_j(abs(order), xs * math.sin(float(th))) for th in theta])  # (nodes, x)
    out: Arr = 0.75 * np.einsum("n,n,nx->x", weight, np.sin(theta) ** 3, j**2)
    return out


GUIDE_LENGTH, GUIDE_VOLTS, GUIDE_PF_PER_M, Q_OHM = 0.277, 13.0, 50e-12, 150.0


def guide_q(mode: int, stubs: int) -> float:
    """The radiation Q of the 277 mm ring guide at its resonance of order m with N switched stubs of 0.5 pF at
    10 mm (dipole moment p0 = C_s V s / 2 each): Q = omega W / (P1 N^2 sum_p F_|m - pN|(kR)), W = C' V^2 L / 2."""
    omega = 2 * math.pi * mode * C0 / GUIDE_LENGTH
    x = omega / C0 * GUIDE_LENGTH / (2 * math.pi)
    stored = GUIDE_PF_PER_M * GUIDE_VOLTS**2 * GUIDE_LENGTH / 2
    p0 = 0.5e-12 * GUIDE_VOLTS * 0.01 / 2
    single = omega**4 * p0**2 / (12 * math.pi * pt.EPS0 * C0**3)
    share = sum(float(order_power(abs(mode - p * stubs), x)[0]) for p in range(-30, 31))
    return omega * stored / (single * stubs**2 * share)


def compute_selection(quick: bool) -> dict[str, Any]:
    say("\n== selection: what a ring of N perturbations carrying mode m radiates, and what it costs the well ==")
    x = np.linspace(0.02, 12.0, 120 if quick else 400)
    f_orders = np.array([order_power(order, x) for order in range(7)])
    total = sum(float(order_power(order, [2.0])[0]) for order in range(-30, 31))
    say(f"  sum over every order of F_l(2.0) = {total:.6f}")
    radius = GUIDE_LENGTH / (2 * math.pi)
    marks = {f"guide m = {m}, {m * C0 / GUIDE_LENGTH / 1e9:.2f} GHz": m * 1.0 for m in (1, 2, 5)}
    for f, r_ring in ((0.433e9, 0.225), (0.9e9, 0.225), (2.3e9, 0.225), (1.5e9, 0.112)):
        marks[f"whips {r_ring * 1e3:.0f} mm, {f / 1e9:g} GHz"] = 2 * math.pi * f / C0 * r_ring
    stubs = {}
    for mode, count in ((2, 3), (2, 8), (5, 3), (5, 8)):
        q = guide_q(mode, count)
        stubs[f"{count} stubs at m = {mode}"] = {"Q_rad": q, "eta": Q_OHM / (Q_OHM + q)}
    omega = 2 * math.pi * C0 / GUIDE_LENGTH
    p = GUIDE_PF_PER_M * GUIDE_VOLTS * 1e-3 * GUIDE_LENGTH
    two_wire = (
        omega
        * GUIDE_PF_PER_M
        * GUIDE_VOLTS**2
        * GUIDE_LENGTH
        / 2
        / (omega**4 * p**2 / (12 * math.pi * pt.EPS0 * C0**3) * float(order_power(1, omega / C0 * radius)[0]))
    )
    stubs["two-wire ring, 1 mm, m = 1"] = {"Q_rad": two_wire, "eta": Q_OHM / (Q_OHM + two_wire)}
    for key, row in stubs.items():
        say(f"  {key:28s} Q_rad {row['Q_rad']:9.3g}, radiates {row['eta'] * 100:5.1f} % at Q_ohm {Q_OHM:.0f}")
    chart = {}
    for mode in range(7):
        for count in range(1, 9):
            orders = sorted({mode - p * count for p in range(-12, 13) if abs(mode - p * count) <= 6}, key=abs)
            chart[(mode, count)] = {"orders": orders, "split": mode > 0 and (2 * mode) % count == 0}
    return {"x": x, "F": f_orders, "sum_F": total, "marks": marks, "stubs": stubs, "chart": chart}


# ---- 7.3: what containment takes -------------------------------------------------------------------------------------
def best_rail_pair(height: float) -> pt.RailPair:
    """The rail pair whose null stands at the given height with the deepest well: rail/gap 3.68 (ELE-3)."""
    s = ((1 + math.sqrt(5)) + math.sqrt((1 + math.sqrt(5)) ** 2 - 4)) / 2  # sqrt(1 + 2 rail/gap) at the optimum
    gap = 2 * height / s
    return pt.RailPair(gap=gap, rail=(s**2 - 1) / 2 * gap)


def focus_depth_ev(el: Elements, point: Sequence[float], loss: float) -> float:
    """The deepest electron well per radiated watt at a point: the pseudo-potential of the largest |E|^2."""
    gf = field(el, [point])[0]
    m = gf @ np.linalg.solve(lossy(radiation(el), loss), gf.conj().T)
    e2 = float(np.linalg.eigvalsh(0.5 * (m + m.conj().T))[-1])
    return pt.ponderomotive_ev(math.sqrt(e2), el.f)


def compute_containment(quick: bool) -> dict[str, Any]:
    say("\n== containment: depth per watt, free space against electrodes behind a resonator ==")
    freqs = np.geomspace(0.4e9, 6.0e9, 60)
    height = 0.5e-3
    boards = pt.FacingBoards(gap=2 * height)
    rails = best_rail_pair(height)
    shares = {"four rods": 1.0, "two facing boards": boards.depth_share, "one board, two rails": rails.depth_share}
    gammas = {
        "four rods": 1.0,
        "two facing boards": boards.gradient_share,
        "one board, two rails": rails.gradient_share,
    }
    ceiling = pt.free_space_depth_ceiling_ev(1.0)
    ladder = {
        f"{name}, {r / 1e3:g} kOhm": np.array(
            [pt.electrode_depth_per_watt_ev(r, height, float(f), share) for f in freqs]
        )
        for name, share in shares.items()
        for r in (1e3, 1e4)
    }
    # the null's Mathieu q when the well is 1 eV deep, against frequency: the same through any resistance
    q_for_1ev = {
        name: np.array([pt.electrode_drive(1.0, 1e4, height, float(f), shares[name], gammas[name]).q for f in freqs])
        for name in shares
    }
    at_1ghz: dict[str, dict[str, float]] = {}
    for name, share in shares.items():
        drive = pt.electrode_drive(1.0, 1e4, height, 1e9, share, gammas[name])
        per_watt = pt.electrode_depth_per_watt_ev(1e4, height, 1e9, share)
        edge = 0.908 / drive.q * drive.volts  # at one frequency q grows as the volts
        half = pt.electrode_drive_at_q(1.0, 0.5, 1e4, height, share, gammas[name])
        at_1ghz[name] = {
            "eV_per_W": per_watt,
            "gain_over_ceiling": pt.electrode_gain(1e4, height, 1e9, share),
            "mW_for_1eV": 1e3 * drive.power,
            "V_for_1eV": drive.volts,
            "q_for_1eV": drive.q,
            "edge_V": edge,
            "edge_mW": 1e3 * edge**2 / 2e4,
            "edge_eV": per_watt * edge**2 / 2e4,
            "gradient_share": gammas[name],
            "q_half_V": half.volts,
            "q_half_mW": 1e3 * half.power,
            "q_half_GHz": half.f / 1e9,
            "q_half_mW_1kOhm": 1e3 * pt.electrode_drive_at_q(1.0, 0.5, 1e3, height, share, gammas[name]).power,
        }
    for name, row in at_1ghz.items():
        say(
            f"  {name:22s} at 1 GHz, {height * 1e3:.1f} mm, 10 kOhm: {row['eV_per_W']:.3g} eV per watt, "
            f"{row['gain_over_ceiling']:.3g} times the ceiling; 1 eV takes {row['V_for_1eV']:.3g} V and "
            f"{row['mW_for_1eV']:.3g} mW at q = {row['q_for_1eV']:.3g}"
            + (", past 0.908: no well" if row["q_for_1eV"] > 0.908 else "")
        )
        say(
            f"  {'':22s} the slope holds to q = 0.908: {row['edge_eV']:.3g} eV at {row['edge_V']:.3g} V and "
            f"{row['edge_mW']:.3g} mW; at q = 0.5 1 eV takes {row['q_half_V']:.3g} V and {row['q_half_mW']:.3g} mW "
            f"at {row['q_half_GHz']:.3g} GHz ({row['q_half_mW_1kOhm']:.3g} mW through 1 kOhm)"
        )
    say(
        f"  the ceiling: {ceiling * 1e6:.1f} micro-eV per radiated watt at any frequency; "
        f"{pt.dipole_wave_peak_field(1.0, 1e9):.0f} V/m at 1 W and 1 GHz; 1 eV takes {1 / ceiling / 1e3:.1f} kW"
    )
    foci = {
        "16-element shell, kb = pi": (shell(2.4e9, 16, C0 / 2.4e9 / 2), (0.0, 0.0, 0.0)),
        "whips 0.433 GHz, 100 mm up": (whips(0.433e9), (0.0, 0.0, 0.1)),
        "whips 0.9 GHz, 100 mm up": (whips(0.9e9), (0.0, 0.0, 0.1)),
        "tile 2.44 GHz, 50 mm out": (panel(2.44e9), (0.0, 0.0, 0.05)),
        "tile 5.8 GHz, 150 mm out": (panel(5.8e9), (0.0, 0.0, 0.15)),
    }
    points = {}
    for name, (el, rf) in foci.items():
        points[name] = {"f_GHz": el.f / 1e9, "eV_per_W": focus_depth_ev(el, rf, LOSS)}
        say(f"  focus of the {name:28s} (10 % loss): {points[name]['eV_per_W'] * 1e6:8.2f} micro-eV per watt")
    bx = np.linspace(-4 * height, 4 * height, 121 if quick else 241)
    by = np.linspace(-0.98 * height, 0.98 * height, 41 if quick else 99)
    gbx, gby = np.meshgrid(bx, by)
    e_boards = boards.field(gbx, gby, 1.0)
    rx = np.linspace(-7 * height, 7 * height, 121 if quick else 281)
    ry = np.linspace(0.02 * height, 7 * height, 81 if quick else 211)
    grx, gry = np.meshgrid(rx, ry)
    e_rails = rails.field(grx, gry, 1.0)
    saddle = rails.escape(1.0)
    say(
        f"  one board at its best: gap {rails.gap * 1e3:.3f} mm, rails {rails.rail * 1e3:.3f} mm, null "
        f"{rails.null_height * 1e3:.2f} mm up, saddle {saddle.height * 1e3:.2f} mm up, share {rails.depth_share:.4f}"
    )
    return {
        "freqs_GHz": freqs / 1e9,
        "height_mm": height * 1e3,
        "ceiling_eV_per_W": ceiling,
        "ladder": ladder,
        "q_for_1eV": q_for_1ev,
        "at_1GHz": at_1ghz,
        "foci": points,
        "boards_x": bx,
        "boards_y": by,
        "boards_field": e_boards,
        "rails_x": rx,
        "rails_y": ry,
        "rails_field": e_rails,
        "rails": {"gap_mm": rails.gap * 1e3, "rail_mm": rails.rail * 1e3, "null_mm": rails.null_height * 1e3}
        | {"saddle_mm": saddle.height * 1e3, "share": rails.depth_share},
        "boards_share": boards.depth_share,
    }


# ---- 7.4: charges as antennas ----------------------------------------------------------------------------------------
def dbm(power: float) -> float:
    return 10 * math.log10(power / 1e-3)


def compute_charges(quick: bool) -> dict[str, Any]:
    say(
        "\n== charges: the 1 GHz, q = 0.3 column on the 277 mm guide, at the limiting density in a radius r0/2, heard "
        "by radiation and by its electrodes =="
    )
    q, f_rf, r0, length = 0.3, 1e9, 0.5e-3, GUIDE_LENGTH
    f_sec = pt.secular_frequency(q, f_rf)
    density = pt.space_charge_limit(q, f_rf)
    f_plasma = pt.plasma_frequency(density)
    depth = pt.quadrupole_depth_ev(q, f_rf, r0)
    full = 4 * math.pi * pt.EPS0 * depth / pt.E  # per metre, a column filling the guide: N'_max = 4 pi eps0 D / e^2
    per_metre = density * math.pi * (r0 / 2) ** 2  # a column of radius r0 / 2 holds a quarter of it
    count = per_metre * length
    omega = 2 * math.pi * f_sec
    amplitude, gap = 1e-4, 2e-3
    one_current = pt.E * omega * amplitude / gap
    larmor = (pt.E * amplitude) ** 2 * omega**4 / (12 * math.pi * pt.EPS0 * C0**3)
    signals = {
        "column, coherent: induced into 50 Ohm": dbm(0.5 * (count * one_current) ** 2 * 50),
        "column, coherent: radiated": dbm(count**2 * larmor),
        "column, Schottky (random phases): induced into 50 Ohm": dbm(count * 0.5 * one_current**2 * 50),
        "column, Schottky (random phases): radiated": dbm(count * larmor),
        "one electron: induced into a 100 kOhm resonant pickup": dbm(0.5 * one_current**2 * 1e5),
        "one electron: induced into 50 Ohm": dbm(0.5 * one_current**2 * 50),
        "one electron: radiated": dbm(larmor),
    }
    thermal = dbm(pt.KB * 290)
    for key, v in signals.items():
        say(f"  {key:54s} {v:7.1f} dBm")
    say(
        f"  k_B T in 1 Hz {thermal:.1f} dBm; the HackRF's floor in 1 Hz {thermal + 16:.0f} to {thermal + 22:.0f} dBm "
        "(noise figure 16 to 22 dB with its amplifier off, REDESIGN.md section 5)"
    )
    l_one = pt.ME * gap**2 / pt.E**2
    damping = pt.E**2 / (6 * math.pi * pt.EPS0 * pt.ME * C0**3) * omega**2
    offsets = np.geomspace(10.0, 1e7, 400)
    columns: dict[str, dict[str, Any]] = {}
    for n in (1e5, 1e6, count):
        ind = l_one / n
        cap = 1 / (omega**2 * ind)
        r_int = n * damping * ind
        w = 2 * math.pi * (f_sec + offsets)
        z_col = 1j * w * ind + 1 / (1j * w * cap) + r_int
        z = 50 * z_col / (50 + z_col)
        columns[f"{n:.2g}"] = {"N": n, "L_H": ind, "C_F": cap, "gamma": 50 / ind, "Q": omega * ind / 50}
        columns[f"{n:.2g}"] |= {"radiation_Q": omega / (n * damping), "Z_over_R": np.abs(z) / 50}
        row = columns[f"{n:.2g}"]
        say(
            f"  {n:9.3g} electrons: {ind * 1e6:9.3g} uH, {cap * 1e12:8.3g} pF, damped into 50 Ohm at "
            f"{row['gamma']:.3g} /s (Q {row['Q']:.3g}), radiation Q {row['radiation_Q']:.3g}"
        )
    say(
        f"  one electron across {gap * 1e3:.0f} mm is {l_one:.0f} H; f_sec {f_sec / 1e6:.0f} MHz, plasma frequency "
        f"at the limit {f_plasma / 1e6:.0f} MHz (density {density:.2g} per m^3; a column filling the guide, "
        f"N'_max = 4 pi eps0 D / e^2 at D = {depth:.2f} eV, {full:.2g} per metre; of radius r0/2 = "
        f"{r0 / 2 * 1e3:.2f} mm, {per_metre:.2g} per metre, {count:.3g} on the guide); coherent current "
        f"{count * one_current * 1e6:.0f} uA; kinetic inductance {pt.ME / (pt.E**2 * per_metre):.2f} H/m, Budker "
        f"parameter {per_metre * 2.8179403e-15:.2g}"
    )
    q6 = 0.6
    omega6 = 2 * math.pi * pt.secular_frequency(q6, f_rf)
    count6 = pt.space_charge_limit(q6, f_rf) * math.pi * (r0 / 2) ** 2 * length
    radiation_q6 = omega6 / (count6 * pt.E**2 / (6 * math.pi * pt.EPS0 * pt.ME * C0**3) * omega6**2)
    say(
        f"  at q = {q6}: the same column holds {count6:.3g} electrons on the guide, f_sec "
        f"{omega6 / (2 * math.pi) / 1e6:.0f} MHz, radiation Q {radiation_q6:.2g}"
    )
    return {
        "signals_dBm": signals,
        "thermal_dBm": thermal,
        "hackrf_floor_dBm": (thermal + 16, thermal + 22),
        "f_sec_Hz": f_sec,
        "f_plasma_Hz": f_plasma,
        "density_m3": density,
        "count": count,
        "per_metre": per_metre,
        "full_per_metre": full,
        "column_radius_m": r0 / 2,
        "count_q06": count6,
        "radiation_Q_q06": radiation_q6,
        "one_electron_H": l_one,
        "offsets_Hz": offsets,
        "columns": columns,
        "coherent_uA": count * one_current * 1e6,
    }


# ---- 10.1: V1 --------------------------------------------------------------------------------------------------------
V1_F = 2.44e9


def v1_grid(ys: Any, xs: Any, zs: Any) -> Arr:
    gx, gy, gz = np.meshgrid(xs, ys, zs, indexing="ij")
    out: Arr = np.c_[gx.ravel(), gy.ravel(), gz.ravel()]
    return out


def probe_response(el: Elements, points: Arr, half: float = 0.02, taps: int = 17) -> CArr:
    """What a vertical dipole of half-length half centred at each point reads: E_z weighted along it by a
    triangular current (weights summing to 1)."""
    zs = np.linspace(-half, half, taps)
    weight = 1 - np.abs(zs) / half
    weight = weight / weight.sum()
    out = np.zeros((len(points), el.count), dtype=np.complex128)
    for z, c in zip(zs, weight, strict=True):
        out += c * field(el, points + np.array([0.0, 0.0, z]), 2)
    return out


def capped_fit(g: CArr, t: CArr, rt: Arr, cap: float) -> CArr:
    """argmin |G w - t|^2 + alpha w^H R~ w with the least alpha that keeps w^H R~ w within the cap (bisection)."""
    scale = float(np.trace(g.conj().T @ g).real / np.trace(rt))
    gram, rhs = g.conj().T @ g, g.conj().T @ t

    def solve(alpha: float) -> CArr:
        w: CArr = np.linalg.solve(gram + alpha * scale * rt, rhs)
        return w

    lo, hi = 1e-14, 1e6
    w = solve(lo)
    if float(np.real(w.conj() @ rt @ w)) <= cap:
        return w
    for _ in range(120):
        mid = math.sqrt(lo * hi)
        if float(np.real(solve(mid).conj() @ rt @ solve(mid))) > cap:
            lo = mid
        else:
            hi = mid
    return solve(hi)


def compute_bench(quick: bool) -> dict[str, Any]:
    say("\n== bench: V1, four quarter-wave monopoles at 2.44 GHz on a lambda/4 square, fitted to a real one ==")
    lam = C0 / V1_F
    length, half_side = lam / 4, lam / 8
    fit_pts = v1_grid(np.arange(100, 221, 30) * 1e-3, np.arange(-90, 91, 30) * 1e-3, np.arange(30, 181, 30) * 1e-3)
    held_pts = v1_grid(np.arange(115, 206, 30) * 1e-3, np.arange(-75, 76, 30) * 1e-3, np.arange(45, 166, 30) * 1e-3)
    corners = np.array(
        [[-half_side, -half_side], [half_side, -half_side], [-half_side, half_side], [half_side, half_side]]
    )
    array = monopoles(V1_F, corners, length, 8, "V1 array")
    r = radiation(array)
    rt = lossy(r)
    one = float(radiation(monopoles(V1_F, [(0.0, 0.0)], length, 8))[0, 0])
    g_fit, g_held = probe_response(array, fit_pts), probe_response(array, held_pts)
    val, _ = generalised_eigh(g_fit.conj().T @ g_fit, rt)
    general = db(np.sort(val)[::-1] / val.max())
    sv = np.linalg.svd(g_fit, compute_uv=False)
    plain = db(sv**2 / sv[0] ** 2)
    say(
        f"  {len(fit_pts)} fit and {len(held_pts)} held-out positions; generalised spectrum in R + rI "
        + ", ".join(f"{v:.1f}" for v in general)
        + " dB; plain singular values "
        + ", ".join(f"{v:.1f}" for v in plain)
        + " dB"
    )

    def predict(xy: Sequence[float]) -> tuple[float, float, float]:
        target = monopoles(V1_F, [tuple(xy)], length, 8)
        t_fit, t_held = probe_response(target, fit_pts)[:, 0], probe_response(target, held_pts)[:, 0]
        w = capped_fit(g_fit, t_fit, rt, 10 * one)
        power = float(np.real(w.conj() @ r @ w)) / one
        return residual_db(g_fit, t_fit, w), residual_db(g_held, t_held, w), power

    table: dict[str, dict[str, Any]] = {}
    named = [("centre", (0.0, 0.0))]
    for label, direction in (("across", (1, 0)), ("away", (0, -1)), ("toward", (0, 1))):
        for frac in (0.125, 0.25, 0.5, 1.0):
            if label == "toward" and frac > 0.5:
                continue
            named.append(
                (
                    f"{label} {frac:g}",
                    (direction[0] * (half_side + frac * lam), direction[1] * (half_side + frac * lam)),
                )
            )
    for name, xy in named:
        fit, held, power = predict(xy)
        table[name] = {"xy_mm": [v * 1e3 for v in xy], "fit_db": fit, "held_db": held, "power_ratio": power}
        say(
            f"  {name:13s} ({xy[0] * 1e3:6.1f}, {xy[1] * 1e3:6.1f}) mm: fit {fit:6.1f} dB, held out {held:6.1f} dB, "
            f"power {power:4.2f} P1"
        )
    curves = {}
    steps = 8 if quick else 36
    for label, direction, reach in (("across", (1, 0), 1.1), ("away", (0, -1), 1.1), ("toward", (0, 1), 0.55)):
        dist = np.linspace(0.0, half_side + reach * lam, steps)
        held_curve = np.array([predict((direction[0] * d, direction[1] * d))[1] for d in dist])
        curves[label] = {"distance_mm": dist * 1e3, "held_db": held_curve}
    gap = table["away 0.5"]["held_db"], table["across 0.5"]["held_db"]
    say(f"  half a wavelength beyond the edge, away beats across by {gap[1] - gap[0]:.1f} dB (criterion V1-6: 10 dB)")
    return {
        "wavelength_mm": lam * 1e3,
        "half_side_mm": half_side * 1e3,
        "whip_mm": length * 1e3,
        "fit_points": fit_pts,
        "held_points": held_pts,
        "corners": corners,
        "generalised_db": general,
        "plain_db": plain,
        "table": table,
        "curves": curves,
    }


# ---- figures ---------------------------------------------------------------------------------------------------------
SURFACE, WARM, PALE = "#fcfcfb", "#fbe3d6", "#ecebe6"
SUPERSCRIPTS = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def sci(value: float, digits: int = 2) -> str:
    """A number for a title: plain between 0.001 and 9999, else m × 10ⁿ."""
    if value == 0 or not math.isfinite(value):
        return f"{value:g}"
    exponent = math.floor(math.log10(abs(value)))
    if -3 <= exponent <= 3:
        return f"{value:.{max(digits - 1 - exponent, 0)}f}"
    mantissa = f"{value / 10**exponent:.{digits - 1}f}".rstrip("0").rstrip(".")
    power = "10" + str(exponent).translate(SUPERSCRIPTS)
    return power if mantissa == "1" else f"{mantissa} × {power}"


def _pyplot() -> Any:
    """matplotlib's pyplot on the Agg backend, set to the house style: one sans face, recessive grid and axes, the
    categorical colours in their fixed order."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from cycler import cycler

    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 9.5,
            "axes.labelsize": 9.5,
            "figure.titlesize": 11,
            "legend.fontsize": 8.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "axes.edgecolor": EDGE,
            "axes.labelcolor": INK,
            "axes.titlecolor": INK,
            "text.color": INK,
            "xtick.color": INK2,
            "ytick.color": INK2,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.axisbelow": True,
            "axes.titlelocation": "left",
            "lines.linewidth": 1.8,
            "lines.markersize": 5,
            "legend.frameon": False,
            "axes.prop_cycle": cycler(color=list(SERIES)),
            "savefig.dpi": 140,
            "mathtext.default": "regular",
        }
    )
    return plt


def _ramp(reverse: bool = False) -> Any:
    """The sequential blue ramp, light to dark (or dark to light)."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("ramp", list(BLUES[::-1] if reverse else BLUES))


def _suptitle(fig: Any, text: str, width: int = 150) -> None:
    """A figure's title, made from its data, wrapped so that it never sets the figure's width."""
    import re
    import textwrap

    glued = re.sub(r"(\d) (?=\S)", "\\1\u00a0", text)  # a number keeps its unit on its line
    fig.suptitle(textwrap.fill(glued, width), x=0.01, ha="left", fontsize=10.5)


def _log_ticks(axis: Any, ticks: Sequence[float], labels: Sequence[str] | None = None) -> None:
    """Fixed ticks with plain labels on a logarithmic axis, and no labels on its minor ticks."""
    from matplotlib.ticker import NullFormatter

    axis.set_ticks(list(ticks))
    axis.set_ticklabels(list(labels) if labels is not None else [f"{t:g}" for t in ticks])
    axis.set_minor_formatter(NullFormatter())


def _save(fig: Any, out: Path, name: str) -> None:
    import matplotlib.pyplot as plt

    fig.savefig(out / f"{name}.png", bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    say(f"  wrote {name}.png")


def _box(
    ax: Any, x: float, y: float, w: float, h: float, title: str, body: str, tone: str = BLUES[1], board: str = ""
) -> None:
    """A rounded box centred at (x, y): a bold title, a body, and beneath a rule what stands for it on the board."""
    from matplotlib.patches import FancyBboxPatch

    ax.add_patch(
        FancyBboxPatch(
            (x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.004,rounding_size=0.025", fc=tone, ec=EDGE, lw=0.8
        )
    )
    top = y + h / 2
    ax.text(x, top - 0.06 * h, title, ha="center", va="top", fontsize=9.5, fontweight="bold", color=INK)
    ax.text(x, top - 0.23 * h, body, ha="center", va="top", fontsize=8.3, color=INK, linespacing=1.3)
    if board:
        ax.plot([x - 0.45 * w, x + 0.45 * w], [y - 0.09 * h, y - 0.09 * h], color=EDGE, lw=0.7)
        ax.text(x, y - 0.13 * h, board, ha="center", va="top", fontsize=7.9, color=INK2, linespacing=1.25)


def _arrow(ax: Any, start: tuple[float, float], end: tuple[float, float], style: str = "arc3") -> None:
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops={
            "arrowstyle": "-|>",
            "color": INK2,
            "lw": 1.1,
            "connectionstyle": style,
            "shrinkA": 1,
            "shrinkB": 1,
        },
    )


def figure_compiler(results: dict[str, Any], out: Path) -> None:
    """Figure 1: the field compiler as a flow, each step's formula in its box and the board's objects beneath it."""
    plt = _pyplot()
    price = "0, 4.8, 5.7, 11.5 dB"
    if "price" in results:
        row = next(r for r in results["price"]["rows"] if r["group"] == "tile" and abs(r["f_GHz"] - 5.8) < 1e-9)
        price = ", ".join(f"{v:.1f}" for v in row["price_db_loss"][:4]) + " dB"
    line1 = "0.980 to 0.992"
    if "alphabets" in results:
        sweep = results["alphabets"]["rev_b_line1"]
        line1 = f"{float(np.min(sweep)):.3f} to {float(np.max(sweep)):.3f}"
    steps = (
        (
            "1  region",
            "T sampled at λ/2 or finer,\nD = min(N, 2L(L+2)) points,\nL = ⌈ka⌉",
            "probe positions\non a printer gantry",
        ),
        (
            "2  operator",
            "G_mn = g_n(r_m), stacked;\nR from the model, a scan\nor the S matrix; R + rI",
            "switched-state frames:\nh = Hᵀy / N, H16",
        ),
        ("3  spectrum", "GᴴG v = σ²(R + rI) v;\nthe price list 1/σᵢ²", f"price list of the tile,\n5.8 GHz: {price}"),
        (
            "4  target, budget",
            "a field, a focus, a quiet zone,\n|E|², a null curve;\ncaps P₀, EIRP, s²",
            "weight error s² = 0.0037\n(3° and 0.26 dB)",
        ),
        (
            "5  solve",
            "w = (GᴴG + αR)⁻¹ Gᴴe;\nfocus w = R⁻¹g;\nα set by the power cap",
            "milliseconds\nper linear solve",
        ),
        (
            "6  map to the board",
            "nearest static state;\ndwell fractions (line 0);\nstaircase on line h",
            f"rev A 4 phases, rev B 16 delays,\nline-1 ceiling {line1};\nsequencer ≤ 256 words",
        ),
        (
            "7  verify",
            "ε = ‖G_held w − t‖ / ‖t‖\non held-out positions,\nby criteria registered first",
            "held-out grid;\narray16.lab",
        ),
        (
            "8  update",
            "h ← h + s̄(y − sᵀh) / ‖s‖²,\none projection per sample",
            "one Kaczmarz step per slot,\n10⁵ per second on rev A",
        ),
    )
    fig, ax = plt.subplots(figsize=(12.5, 5.0), layout="constrained")
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 2.3)
    ax.axis("off")
    w, h, top, bottom = 0.9, 0.82, 1.79, 0.63
    centres = [
        (0.5, top),
        (1.5, top),
        (2.5, top),
        (3.5, top),
        (3.5, bottom),
        (2.5, bottom),
        (1.5, bottom),
        (0.5, bottom),
    ]
    for (title, body, board), (x, y) in zip(steps, centres, strict=True):
        _box(ax, x, y, w, h, title, body, board=board)
    for i in (0, 1, 2):
        _arrow(ax, (centres[i][0] + w / 2, top), (centres[i + 1][0] - w / 2, top))
    _arrow(ax, (3.5, top - h / 2), (3.5, bottom + h / 2))
    for i in (4, 5, 6):
        _arrow(ax, (centres[i][0] - w / 2, bottom), (centres[i + 1][0] + w / 2, bottom))
    _arrow(ax, (0.75, bottom + h / 2), (2.35, top - h / 2))
    ax.text(
        0.03, 1.27, "back to 3 while the held-out\nresidual fails its criterion", fontsize=8.3, color=INK2, va="center"
    )
    ax.text(
        0.0,
        0.07,
        "9  band: steps 2 to 6 on a frequency grid of about 500 MHz for the tile and 150 MHz for the whips, kept as a "
        "table;"
        " stop when no σᵢ moves by more than the calibration residual",
        fontsize=8.3,
    )
    ax.set_title(
        "The field compiler: one operator and one loop, from a region of space to a drive on the board and back",
        fontsize=11,
    )
    _save(fig, out, "compiler")


def figure_price(res: dict[str, Any], out: Path) -> None:
    """Figure 2: the price list 1/sigma_i^2 against mode index, with the knee of the degree-of-freedom count."""
    plt = _pyplot()
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.3), sharey=True, layout="constrained")
    names = {"tile": "the tile", "whips": "the 16 whips", "shell": "a 16-element reference shell, radius 250 mm"}
    for ax, group in zip(axes, ("tile", "whips", "shell"), strict=True):
        rows = [r for r in res["rows"] if r["group"] == group]
        for i, row in enumerate(rows):
            idx = np.arange(1, len(row["price_db_loss"]) + 1)
            ax.plot(idx, row["price_db_loss"], "-o", color=SERIES[i], ms=3.5, label=f"{row['f_GHz']:g} GHz")
            ax.plot(idx, row["price_db_lossless"], "--", color=SERIES[i], lw=1.1)
            ax.plot([row["knee"]], [-2.2], "^", color=SERIES[i], ms=9, clip_on=False, zorder=5)
        ax.axhline(10, color=MUTED, lw=0.8, ls=":")
        ax.text(16.4, 10.8, "10 dB", fontsize=7.8, color=INK2, ha="right")
        first = rows[0]
        where = "at its centre" if group == "shell" else f"{first['centre_mm'][2]:.0f} mm up"
        counts = " / ".join(str(r["within_10dB_loss"]) for r in rows)
        ax.set_title(
            f"{names[group]}\na ball of radius {first['ball_mm']:.0f} mm, {where}; within 10 dB: {counts}",
            fontsize=9,
        )
        ax.set_xlabel("mode index i")
        ax.set_xticks([1, 4, 8, 12, 16])
        ax.set_ylim(-2.5, 72)
        ax.set_xlim(0.5, 16.5)
        ax.legend(loc="upper left")
    axes[0].set_ylabel("price 1/σᵢ², dB above the cheapest mode")
    tile = [r["within_10dB_loss"] for r in res["rows"] if r["group"] == "tile"]
    whips = [r["within_10dB_loss"] for r in res["rows"] if r["group"] == "whips"]
    _suptitle(
        fig,
        f"The price list 1/σᵢ² of the vector field in a ball (solid with 10 % loss, dashed lossless; triangles the "
        "knee "
        f"min(N, source count, 2L(L+2))): the tile makes {tile[0]}, {tile[1]} and {tile[2]} field shapes within 10 dB "
        f"of its cheapest at 2.3, 3.5 and 5.8 GHz, the whips "
        + (f"{whips[0]} at each frequency" if len(set(whips)) == 1 else ", ".join(str(n) for n in whips))
        + ", a closed shell all 16",
    )
    _save(fig, out, "price")


def figure_radiation(res: dict[str, Any], out: Path) -> None:
    """Figure 3: the eigenvalues of R across 0.4-6 GHz for the tile and one fixed set of whips."""
    plt = _pyplot()
    freqs = np.asarray(res["freqs_GHz"])
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4), sharey=True, layout="constrained")
    mesh = None
    for ax, key, title in (
        (axes[0], "tile_db", "the tile over ground, 25 mm pitch"),
        (axes[1], "whips_db", f"the 16 whips, each {res['whip_cut_mm']:.0f} mm (a quarter wave at 0.9 GHz)"),
    ):
        data = np.asarray(res[key])
        mesh = ax.pcolormesh(
            freqs, np.arange(1, data.shape[1] + 1), data.T, cmap=_ramp(), vmin=-60, vmax=0, shading="nearest"
        )
        count = (data > -10).sum(1)
        ax.step(freqs, count + 0.5, where="mid", color=SERIES[1], lw=2.2)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("frequency (GHz)")
        ax.set_xlim(freqs[0], freqs[-1])
        ax.grid(False)
    lam = C0 / (freqs * 1e9)
    axes[0].plot(freqs, np.clip(math.pi * 0.01 / lam**2, 0, 16.5), color=INK, lw=1.3, ls="--")
    axes[0].text(4.6, 6.0, "πA/λ², A = 0.01 m²", fontsize=8.3, color=SURFACE, rotation=38)
    axes[0].text(0.45, 1.85, "10 dB line", fontsize=8.3, color=SERIES[1])
    axes[1].axvline(0.9, color=SURFACE, lw=1.0, ls=":")
    axes[1].text(0.95, 15.6, "0.9 GHz", fontsize=8, color=SURFACE)
    axes[0].set_ylabel("eigenvalue index, the largest at the bottom")
    axes[0].set_ylim(0.5, 16.5)
    cbar = fig.colorbar(mesh, ax=axes, fraction=0.03, pad=0.01)
    cbar.set_label("eigenvalue of R, dB of the largest")
    t = res["table"]
    keys = ("tile 2.3 GHz", "tile 3.5 GHz", "tile 5.8 GHz")
    ratios = [t[k]["within_10dB"] / t[k]["pi_A_over_lambda2"] for k in keys]
    _suptitle(
        fig,
        f"How many drives radiate: the tile has {t[keys[0]]['within_10dB']}, {t[keys[1]]['within_10dB']} and "
        f"{t[keys[2]]['within_10dB']} eigenvalues of R within 10 dB of its largest (orange) at 2.3, 3.5 and 5.8 GHz, "
        f"{min(ratios):.1f} to {max(ratios):.1f} times πA/λ²; its smallest lies "
        f"{-t[keys[0]]['smallest_below_mean_dB']:.0f} dB below the mean at 2.3 GHz: the superdirective drives",
    )
    _save(fig, out, "radiation")


def figure_fine(res: dict[str, Any], out: Path) -> None:
    """Figure 4: the power that holds M alternating values at a spacing, against M, with the law's slope."""
    plt = _pyplot()
    fig, (left, right) = plt.subplots(
        1, 2, figsize=(12.5, 4.4), gridspec_kw={"width_ratios": [1.15, 1]}, layout="constrained"
    )
    m = np.arange(1, 7)
    for i, (key, row) in enumerate(res["shell"].items()):
        d = float(key)
        left.plot(m, row, "-o", color=SERIES[i], ms=4, label=f"δ = {d:g} λ")
        if d < 0.36:
            left.plot(m, 20 * (m - 1) * math.log10(0.63 / d), ":", color=SERIES[i], lw=1.4)
    tile = np.asarray(res["tile"]["1.5"])
    left.plot(m, tile, "--s", color=INK2, ms=4, lw=1.2, label="the tile, δ = 1.5 λ, 150 mm out")
    left.set_xlabel("number of alternating values M")
    left.set_ylabel("least radiated power, dB above one value")
    left.set_title("a 64-element shell at 5.8 GHz; dotted, the law (0.63 λ/δ)^(2(M−1))", fontsize=9)
    left.legend(fontsize=8, loc="upper left")
    spacing = np.asarray(res["dense_spacing_lambda"])
    per = np.asarray(res["dB_per_point"])
    right.semilogx(spacing, per, "-o", color=SERIES[0], ms=4, label="64-element shell, mean over M = 2 to 6")
    law = np.geomspace(spacing[0], 0.35, 50)
    right.semilogx(law, 20 * np.log10(0.63 / law), ":", color=INK, lw=1.5, label="20 log₁₀(0.63 λ/δ) per value")
    right.axvline(0.35, color=MUTED, lw=0.8, ls="--")
    right.text(0.36, 0.5, "0.35 λ", fontsize=8, color=INK2)
    _log_ticks(right.xaxis, [0.1, 0.15, 0.25, 0.35, 0.5, 1.0])
    right.set_xlabel("spacing δ (wavelengths)")
    right.set_ylabel("dB per added value")
    right.set_title("the cost per value: a power law below 0.35 λ, about 1 dB from λ/2 up", fontsize=9)
    right.legend(fontsize=8, loc="upper right")
    six = res["shell"]["0.15"][-1]
    half = res["shell"]["0.5"][-1]
    _suptitle(
        fig,
        f"The price of fine structure: six alternating values cost {six:.0f} dB more than one at 0.15 λ and "
        f"{half:.1f} dB at λ/2, exponential in the number of sub-wavelength features; along a line 150 mm before the "
        f"tile, four columns give four values and the fifth costs {tile[4] - tile[3]:.1f} dB more",
    )
    _save(fig, out, "fine")


def figure_ceiling(res: dict[str, Any], out: Path) -> None:
    """Figure 5: a focus as a share of the free-space ceiling RFG-2, for shells and for the board's arrays."""
    plt = _pyplot()
    fig, (left, right) = plt.subplots(1, 2, figsize=(12.5, 4.8), layout="constrained")
    kb = np.asarray(res["kb"])
    for x in np.asarray(res["cavity_kb"]):
        if x < 20:
            left.axvline(x, color=MUTED, lw=0.7, ls=(0, (1, 3)))
    for i, (count, row) in enumerate(res["shells"].items()):
        n = int(count)
        left.loglog(kb, row["lossless"], color=SERIES[i], lw=1.5, label=f"N = {n}")
        left.loglog(kb, row["loss"], "--", color=SERIES[i], lw=1.1)
        left.loglog(kb, 3 * n / (16 * kb**2), ":", color=SERIES[i], lw=1.0)
        limit = math.sqrt(10 * n)
        left.plot([limit], [3 * n / (16 * limit**2)], "o", mfc=SURFACE, mec=SERIES[i], ms=7, mew=1.6)
    for ring, offset in zip(res["rings"], ((6, 2), (6, -11), (6, 2)), strict=True):
        x = 2 * math.pi * ring["f_GHz"] * 1e9 / C0 * 0.225
        left.plot([x], [ring["loss_free"]], "D", color=INK, ms=5)
        left.annotate(
            f"whips, {ring['f_GHz']:g} GHz", (x, ring["loss_free"]), offset, textcoords="offset points", fontsize=7.8
        )
    left.axhline(1, color=INK, lw=1.0)
    left.text(kb[-1], 1.18, "the ceiling 4Z₀k²P/(3π)", ha="right", fontsize=8)
    left.set_xlabel("shell radius kb; dotted verticals, the cavity resonances d(x j₁(x))/dx = 0")
    left.set_ylabel("|E|² per watt at the centre, share of the ceiling")
    left.set_title(
        "N tangential dipoles on a shell: lossless solid, 10 % loss dashed,\n3N/(16(kb)²) dotted, circles where "
        "N = 0.05·2(kb)²",
        fontsize=9,
    )
    left.legend(fontsize=8, loc="lower left")
    left.set_ylim(1e-3, 100)
    _log_ticks(left.xaxis, [1, 3, 10, 30, 60])
    right.axvspan(0.2, 3, color=GRID, alpha=0.8, lw=0)
    right.text(0.55, 30, "kr < 3", fontsize=8.3, color=INK2)
    for i, (label, row) in enumerate(res["board"].items()):
        right.loglog(row["kr"], row["lossless"], "-o", color=SERIES[i], ms=3.5, label=label)
        right.loglog(row["kr"], row["loss"], "--", color=SERIES[i], lw=1.1)
        for j, offset in ((0, ((5, 2), (6, -2), (-38, -12))[i]), (len(row["kr"]) - 1, (4, 3))):
            right.annotate(
                f"{row['height_mm'][j]:.0f} mm",
                (row["kr"][j], row["lossless"][j]),
                offset,
                textcoords="offset points",
                fontsize=7.4,
                color=INK2,
            )
    right.axhline(1, color=INK, lw=1.0)
    right.set_xlim(0.5, 80)
    _log_ticks(right.xaxis, [1, 3, 10, 30])
    right.set_xlabel("kr from the focus to the nearest element")
    right.set_ylabel("share of the ceiling over the ground")
    right.set_title("foci on the board's axis, 20 to 500 mm up; dashed with 10 % loss", fontsize=9)
    right.legend(fontsize=8, loc="lower left")
    sixteen = res["at_kb_pi"][16] if 16 in res["at_kb_pi"] else res["at_kb_pi"]["16"]
    dense = res["passing"]["256 at kb = 6.1"]
    top = max(float(np.max(row["lossless"])) for row in res["board"].values())
    _suptitle(
        fig,
        f"A focus against the free-space ceiling (RFG-2): 16 elements on a shell of radius λ/2 reach {sixteen[0]:.2f} "
        f"of it ({sixteen[1]:.2f} with 10 % loss). Passing it takes stored field: a shell dense enough to close, at "
        f"its cavity resonances (256 elements at kb = 6.1, {dense[0]:.0f} times; {dense[1]:.0f} with 10 % loss), or a "
        f"focus in an element's own near field ({top:.0f} times the ceiling over the ground, 20 mm before the tile "
        "at 2.44 GHz)",
    )
    _save(fig, out, "ceiling")


def figure_band(res: dict[str, Any], out: Path) -> None:
    """Figure 6: how far the operator moves across frequency."""
    plt = _pyplot()
    fig, ax = plt.subplots(figsize=(8.5, 4.6), layout="constrained")
    offsets = np.asarray(res["offsets_Hz"]) / 1e6
    for i, (name, row) in enumerate(res["rows"].items()):
        ax.loglog(offsets, row, "-o", color=SERIES[i], ms=4, label=name)
    ax.axvline(20, color=INK, lw=1.2)
    ax.text(21.5, 2.3e-7, "the HackRF's\n20 MHz", fontsize=8, va="bottom")
    for x, text in ((500, "table grid, tile"), (150, "table grid, whips")):
        ax.axvline(x, color=MUTED, lw=0.9, ls="--")
        ax.text(x * 1.06, 2.3e-7, text, rotation=90, fontsize=7.8, color=INK2, va="bottom")
    ax.axhline(0.01, color=MUTED, lw=0.8, ls=":")
    ax.text(1.0, 0.012, "1 %", fontsize=8, color=INK2)
    ax.set_xlim(0.8, 1300)
    ax.set_ylim(1e-7, 1.5)
    _log_ticks(ax.xaxis, [1, 10, 100, 1000])
    ax.set_xlabel("frequency offset Δf (MHz)")
    ax.set_ylabel("1 − |corr|² between G(f) and G(f + Δf)")
    ax.legend(fontsize=8.5, loc="upper left")
    i20 = int(np.argmin(np.abs(offsets - 20)))
    rows = list(res["rows"].values())
    _suptitle(
        fig,
        f"Across the HackRF's 20 MHz the operator moves by {sci(rows[0][i20])} (tile, 5.8 GHz) and "
        f"{sci(rows[2][i20])} (whips), and rev B's longest delay turns {res['longest_delay_turn_deg']:.1f}°: one "
        "solution per capture, and a table on a grid of about 500 MHz (tile) and 150 MHz (whips) across the band",
        width=105,
    )
    _save(fig, out, "band")


def _tile_and_window(ax: Any, depth_mm: float, top_mm: float) -> None:
    """The tile seen edge-on at z = 0 and the rays from a source depth_mm behind it through the aperture's edges."""
    half = TILE_SIDE / 2 * 1e3
    for s in (-1.0, 1.0):
        ax.plot([s * half, s * half * (1 + top_mm / depth_mm)], [0, top_mm], color=INK, lw=1.0, ls=":")
    ax.plot([-half, half], [0, 0], color=INK, lw=4, solid_capstyle="butt")


def figure_behind(res: dict[str, Any], out: Path) -> None:
    """Figure 7: a virtual dipole behind the tile: error maps with the window and the edge term, and the residual
    against cone half-angle for sources 25 to 300 mm behind."""
    plt = _pyplot()
    fig = plt.figure(figsize=(12.5, 5.4), layout="constrained")
    grid = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.15])
    xs, zs = np.asarray(res["map_x"]) * 1e3, np.asarray(res["map_z"]) * 1e3
    depth = float(res["map_depth_mm"])
    maps = []
    mesh = None
    for j, key in enumerate(("5.8 GHz", "2.4 GHz")):
        ax = fig.add_subplot(grid[0, j])
        maps.append(ax)
        mesh = ax.pcolormesh(xs, zs, res["maps"][key], cmap=_ramp(reverse=True), vmin=-40, vmax=0, shading="auto")
        edge = np.asarray(res["edge_db"][key])
        lines = ax.contour(xs, zs, edge, levels=[-10], colors=[SERIES[1]], linestyles=["--"], linewidths=1.4)
        ax.clabel(lines, fmt={-10.0: "edge term −10 dB"}, fontsize=7.5)
        _tile_and_window(ax, depth, float(zs[-1]))
        ax.plot([0], [-depth], "*", color=SERIES[1], ms=13, clip_on=False)
        ax.text(14, -depth, f"virtual dipole {depth:.0f} mm behind", fontsize=7.8, va="center")
        ax.set_xlim(xs[0], xs[-1])
        ax.set_ylim(-depth - 15, zs[-1])
        ax.set_aspect("equal")
        ax.grid(False)
        ax.set_xlabel("x (mm)")
        inside = np.abs(xs[None, :]) < TILE_SIDE / 2 * 1e3 * (zs[:, None] + depth) / depth
        far = zs[:, None] > 100
        ax.set_title(
            f"{key}: median {float(np.median(res['maps'][key][inside & far])):.0f} dB inside the window,\n"
            f"{float(np.median(res['maps'][key][~inside & far])):.0f} dB outside it (beyond 100 mm)",
            fontsize=9,
        )
    maps[0].set_ylabel("z, in front of the tile (mm)")
    cbar = fig.colorbar(mesh, ax=maps, location="bottom", fraction=0.05, shrink=0.7)
    cbar.set_label("|E − E_dipole| / |E_dipole| (dB); dotted, the window of rays through the aperture")
    ax = fig.add_subplot(grid[0, 2])
    cones = np.asarray(res["cones_deg"])
    for i, (key, row) in enumerate(res["rows"].items()):
        ax.plot(cones, row["fit_db"], "-o", color=SERIES[i], ms=4, label=key)
        ax.plot(cones, row["closed_form_db"], "--", color=SERIES[i], lw=1.1)
        ax.plot([row["window_deg"]], [1.2], "v", color=SERIES[i], ms=7, clip_on=False)
    ax.set_xlabel("cone half-angle about boresight (degrees)")
    ax.set_ylabel("residual over shells at 0.3, 0.6 and 2 m (dB)")
    ax.set_ylim(-38, 1.5)
    ax.set_title(
        "fitted weights (solid) against the closed form with one\nfree scale (dashed); triangles, each window",
        fontsize=9,
    )
    ax.legend(fontsize=7.8, loc="lower right")
    row = res["rows"]["50 mm behind, 5.8 GHz"]
    _suptitle(
        fig,
        f"A virtual dipole behind the tile is reproduced inside the window of rays through the aperture: 50 mm behind "
        f"at 5.8 GHz, {row['fit_db'][0]:.0f} dB in a 20° cone, {row['fit_db'][1]:.0f} dB in 30° and "
        f"{row['fit_db'][3]:.0f} dB in 60°, its window being 45°; the fit is the process and the closed form its start",
    )
    _save(fig, out, "behind")


def _width_3db(xs: Arr, row_db: Arr) -> float:
    """The width of the region about the middle of a cut where it stays within 3 dB of its value there."""
    mid = len(xs) // 2
    level = float(row_db[mid]) - 3.0
    lo = mid
    while lo > 0 and row_db[lo - 1] >= level:
        lo -= 1
    hi = mid
    while hi < len(xs) - 1 and row_db[hi + 1] >= level:
        hi += 1
    return float(xs[hi] - xs[lo])


def figure_focused(res: dict[str, Any], out: Path) -> None:
    """Figure 8: a virtual source focused 60 mm in front of the tile: intensity, phase, its error against a true
    source at the focus, and the tile's on-axis maxima under conjugate weights."""
    plt = _pyplot()
    fig = plt.figure(figsize=(12.5, 9.0), layout="constrained")
    grid = fig.add_gridspec(2, 3, height_ratios=[1.3, 0.75])
    xs, zs = np.asarray(res["map_x"]) * 1e3, np.asarray(res["map_z"]) * 1e3
    zf = float(res["focus_mm"][2])
    half = TILE_SIDE / 2 * 1e3
    reach = math.tan(math.radians(float(res["cone_deg"])))
    panels = (
        ("intensity_db", "|E|², dB of its value at the focus", _ramp(), -25.0, 5.0),
        ("phase", "phase of E_y (radians)", "twilight", -math.pi, math.pi),
        ("error_db", "|E − E_source| / |E_source| (dB)", _ramp(reverse=True), -40.0, 0.0),
    )
    row = int(np.argmin(np.abs(zs - zf)))
    width = _width_3db(xs, np.asarray(res["intensity_db"])[row])
    lam = C0 / 5.8e9 * 1e3
    titles = (
        f"intensity: {width:.0f} mm wide (−3 dB) across the focus, λ/2 = {lam / 2:.0f} mm",
        "phase of E_y: beyond the focus, circles about it",
        "error against a true source at the focus",
    )
    maps = []
    for j, ((key, label, cmap, lo, hi), title) in enumerate(zip(panels, titles, strict=True)):
        ax = fig.add_subplot(grid[0, j])
        maps.append(ax)
        mesh = ax.pcolormesh(xs, zs, res[key], cmap=cmap, vmin=lo, vmax=hi, shading="auto")
        ray = SURFACE if j == 1 else INK
        for s in (-1.0, 1.0):
            ax.plot([s * half, 0], [0, zf], color=ray, lw=0.9, ls=":")
            ax.plot([0, s * (zs[-1] - zf) * reach], [zf, zs[-1]], color=ray, lw=1.0, ls="--")
        ax.plot([-100, 100], [200, 200], color=SERIES[1], lw=2.5)
        ax.plot([0], [zf], "+", color=SERIES[1], ms=11, mew=2)
        ax.plot([-half, half], [0, 0], color=INK, lw=4, solid_capstyle="butt")
        ax.set_aspect("equal")
        ax.set_xlim(xs[0], xs[-1])
        ax.set_ylim(0, zs[-1])
        ax.grid(False)
        ax.set_xlabel("x (mm)")
        ax.set_title(title, fontsize=9)
        cbar = fig.colorbar(mesh, ax=ax, location="bottom", fraction=0.06, pad=0.02)
        cbar.set_label(label, fontsize=8.5)
    maps[0].set_ylabel("z, in front of the tile (mm)")
    axis = np.asarray(res["axis_mm"])
    profiles = []
    for j, freq in enumerate((5.8, 2.44)):
        ax = fig.add_subplot(grid[1, 0:2] if j == 0 else grid[1, 2])
        profiles.append(ax)
        rows = [row for row in res["profiles"].values() if abs(row["f_GHz"] - freq) < 1e-9]
        for i, row in enumerate(rows):
            ax.plot(axis, db(row["profile"]), color=SERIES[i], lw=1.5, label=f"asked at {row['asked_mm']:.0f} mm")
            ax.axvline(row["asked_mm"], color=SERIES[i], lw=0.8, ls=":")
            if row["peak_mm"] is not None:
                k = int(np.argmin(np.abs(axis - row["peak_mm"])))
                ax.plot([row["peak_mm"]], [float(db(row["profile"][k]))], "o", color=SERIES[i], ms=6)
        ax.set_xscale("log")
        ax.set_xlim(axis[0], axis[-1])
        _log_ticks(ax.xaxis, [20, 50, 100, 200, 500])
        ax.set_ylim(-30, 1.5)
        ax.set_xlabel("distance on the axis (mm)")
        peaks = [f"{row['peak_mm']:.0f}" for row in rows if row["peak_mm"] is not None]
        ax.set_title(
            f"{freq:g} GHz, conjugate weights: "
            + (f"on-axis maxima at {', '.join(peaks)} mm (dots)" if peaks else "no on-axis maximum"),
            fontsize=9,
        )
        ax.legend(fontsize=7.8, loc="lower left")
    profiles[0].set_ylabel("|E|² on the axis (dB of its peak)")
    r = res["residual_db"]
    _suptitle(
        fig,
        f"A field-only virtual antenna in empty space is a focus: a source fitted 60 mm in front of the tile at "
        "5.8 GHz "
        f"(cross; orange bar, the 200 mm square at 0.2 m it is fitted on) equals the field of a true source there only "
        f"beyond the focus, inside the {res['cone_deg']:.0f}° cone: {r['2.44 GHz']:.1f}, {r['3.5 GHz']:.1f} and "
        f"{r['5.8 GHz']:.1f} dB on the square at 2.44, 3.5 and 5.8 GHz",
    )
    _save(fig, out, "focused")


def figure_whips(res: dict[str, Any], out: Path) -> None:
    """Figure 9: the whip rings at 0.9 GHz: residual and cost against azimuthal order, the weight each ring gives
    an order, and the residual of a virtual whip against its position over the plate."""
    plt = _pyplot()
    fig = plt.figure(figsize=(12.5, 9.2), layout="constrained")
    grid = fig.add_gridspec(3, 2, width_ratios=[1.05, 1], height_ratios=[1, 0.8, 0.8])
    orders = np.asarray(res["orders"])
    colours = {"both rings": SERIES[0], "outer ring": SERIES[1], "inner ring": SERIES[2]}
    top = fig.add_subplot(grid[0, 0])
    weight = fig.add_subplot(grid[1, 0], sharex=top)
    cost = fig.add_subplot(grid[2, 0], sharex=top)
    for name, colour in colours.items():
        row = res["patterns"][name]
        top.plot(orders, row["residual_db"], "-o", color=colour, ms=4.5, label=name)
        made = np.asarray(row["residual_db"]) < -6
        cost.plot(orders, np.where(made, row["S_db"], np.nan), "-o", color=colour, ms=4.5, label=name)
    top.set_ylabel("residual (dB)")
    top.set_title("exp(imφ) times a whip's pattern, all azimuths, elevations 0 to 45°, 1 to 6 m", fontsize=9)
    top.legend(fontsize=8, loc="lower right")
    for (name, row), colour in zip(res["bessel"].items(), (SERIES[1], SERIES[2]), strict=True):
        j = np.abs(np.asarray(row["J"]))
        weight.fill_between(orders, j.min(1), j.max(1), color=colour, alpha=0.25, lw=0)
        weight.plot(orders, j.min(1), "-", color=colour, lw=1.4, label=f"{name}, ka = {row['ka']:.2f}")
        weight.plot(orders, j.max(1), "-", color=colour, lw=0.8)
    weight.set_ylabel("|J_m(ka cos el)|")
    weight.set_title(
        f"the weight |J_m(ka cos el)| of order m, least to largest over 0 to 45°; outer J₁ = 0 at "
        f"{res['J1_zero_deg']:.1f}°",
        fontsize=9,
    )
    weight.legend(fontsize=8, loc="upper right")
    cost.set_ylabel("S against one whip (dB)")
    cost.set_xlabel("azimuthal order m")
    cost.set_title("what each order costs, where it is made (residual under −6 dB): S = ‖w‖² / (wᴴRw)", fontsize=9)
    cost.set_xticks(orders)
    plate = fig.add_subplot(grid[:, 1])
    xs = np.asarray(res["map_x"]) * 1e3
    mesh = plate.pcolormesh(xs, xs, res["map_residual_db"], cmap=_ramp(reverse=True), vmin=-45, vmax=0, shading="auto")
    places = np.asarray(res["places"]) * 1e3
    plate.plot(places[:8, 0], places[:8, 1], "o", mfc=SERIES[1], mec=INK, ms=8, label="outer ring, 225 mm")
    plate.plot(places[8:, 0], places[8:, 1], "o", mfc=SERIES[2], mec=INK, ms=8, label="inner ring, 112 mm")
    turn = np.linspace(0, 2 * math.pi, 200)
    plate.plot(200 * np.cos(turn), 200 * np.sin(turn), color=SERIES[3], lw=1.4, ls="--", label="200 mm from the centre")
    plate.set_aspect("equal")
    plate.grid(False)
    plate.set_xlabel("x of the virtual whip (mm)")
    plate.set_ylabel("y of the virtual whip (mm)")
    plate.legend(fontsize=7.8, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=3)
    cbar = fig.colorbar(mesh, ax=plate, location="right", fraction=0.05, shrink=0.6)
    cbar.set_label("residual over 1 to 6 m, elevations 0 to 45° (dB)")
    radius = np.hypot(*np.meshgrid(xs, xs))
    worst = {r: float(np.max(np.asarray(res["map_residual_db"])[radius <= r + 1e-9])) for r in (150, 200)}
    plate.set_title(
        f"a virtual quarter-wave whip at each place: {worst[150]:.0f} dB or better within 150 mm of the centre, "
        f"{worst[200]:.0f} dB within 200 mm",
        fontsize=9,
    )
    both = res["patterns"]["both rings"]
    _suptitle(
        fig,
        f"The whip rings at 0.9 GHz as an interior array of order about 3: exp(imφ) to "
        f"{max(both['residual_db'][:3]):.0f}"
        f" dB for m ≤ 2, {both['residual_db'][3]:.0f} dB at m = 3 for S = {both['S_db'][3]:.0f} dB and "
        f"{both['residual_db'][4]:.0f} dB at m = 4 for S = {both['S_db'][4]:.0f} dB: an order beyond ka costs "
        "sensitivity, the size asymmetry in one ring",
    )
    _save(fig, out, "whips")


def figure_alphabets(res: dict[str, Any], out: Path) -> None:
    """Figure 10: the weights each alphabet reaches: its hull on line 0, the disc of radius perimeter / 2 pi on
    line 1; and both across the band."""
    plt = _pyplot()
    fig = plt.figure(figsize=(12.5, 7.6), layout="constrained")
    grid = fig.add_gridspec(2, 4, height_ratios=[1, 0.75])
    turn = np.linspace(0, 2 * math.pi, 361)
    handles: list[Any] = []
    for j, panel_row in enumerate(res["panels"]):
        ax = fig.add_subplot(grid[0, j])
        ph = np.sort(np.mod(np.asarray(panel_row["phases"]), 2 * math.pi))
        pts = np.exp(1j * ph)
        ax.plot(np.cos(turn), np.sin(turn), color=EDGE, lw=0.9)
        hull = ax.fill(pts.real, pts.imag, fc=BLUES[1], ec=SERIES[0], lw=1.3, label="line 0: the hull")
        line1, inner = panel_row["line1"], panel_row["inradius"]
        disc = ax.plot(line1 * np.cos(turn), line1 * np.sin(turn), "--", color=SERIES[1], lw=1.6)
        ring = ax.plot(inner * np.cos(turn), inner * np.sin(turn), ":", color=SERIES[0], lw=1.4)
        dots = ax.plot(pts.real, pts.imag, "o", color=INK, ms=4)
        if j == 0:
            handles = [dots[0], hull[0], ring[0], disc[0]]
        ax.set_aspect("equal")
        ax.set_xlim(-1.15, 1.15)
        ax.set_ylim(-1.15, 1.15)
        ax.set_xticks([-1, 0, 1])
        ax.set_yticks([-1, 0, 1])
        ax.set_xlabel("Re")
        ax.set_title(
            f"{panel_row['board']} at {panel_row['f_GHz']:g} GHz, {len(ph)} states:\n"
            f"line 1 {line1:.3f}, line 0 {inner:.3f}; gap {panel_row['largest_gap_deg']:.0f}°",
            fontsize=9,
        )
    fig.legend(
        handles,
        ["the states", "line 0: their hull", "the hull's inradius", "line 1: the disc of radius perimeter/2π"],
        loc="outside lower center",
        ncol=4,
        fontsize=8.3,
    )
    first = fig.add_subplot(grid[1, :2])
    first.plot(res["rev_b_GHz"], res["rev_b_line1"], color=SERIES[0], label="rev B, sixteen delays")
    first.plot(res["rev_b_GHz"], res["rev_b_polarity_line1"], color=SERIES[2], label="rev B with its polarity")
    first.plot(res["rev_a_GHz"], res["rev_a_line1"], color=SERIES[1], label="rev A, 2-bit delay lines")
    first.axhline(float(np.sinc(1 / 16)), color=MUTED, lw=0.9, ls=":")
    first.text(4.95, float(np.sinc(1 / 16)) - 0.006, "a 4-bit phase grid, 0.9936", fontsize=7.8, color=INK2)
    first.set_ylim(0.86, 1.0)
    first.set_ylabel("line-1 ceiling")
    first.set_xlabel("frequency (GHz)")
    least = int(np.argmin(res["rev_b_line1"]))
    first.set_title(
        f"line 1: rev B keeps {float(np.min(res['rev_b_line1'])):.4f} to {float(np.max(res['rev_b_line1'])):.4f}, "
        f"least at {res['rev_b_GHz'][least]:.2f} GHz",
        fontsize=9,
    )
    first.legend(fontsize=8, loc="lower left")
    second = fig.add_subplot(grid[1, 2:], sharex=first)
    second.plot(res["rev_b_GHz"], res["rev_b_inradius"], color=SERIES[0], label="rev B, sixteen delays")
    second.plot(res["rev_a_GHz"], res["rev_a_inradius"], color=SERIES[1], label="rev A, 2-bit delay lines")
    second.set_ylim(0.4, 1.0)
    second.set_ylabel("line-0 inradius")
    second.set_xlabel("frequency (GHz)")
    second.set_xlim(2.3, 6.0)
    second.set_title("line 0: the largest disc inside the hull, every phase at that amplitude", fontsize=9)
    second.legend(fontsize=8, loc="lower left")
    a = res["panels"][0]
    _suptitle(
        fig,
        f"Continuous weights from discrete states: a periodic sequence reaches the hull of its states on line 0 and, "
        f"on line 1, a disc of radius perimeter/2π, {a['line1']:.3f} for rev A's four phases and within 0.17 dB of a "
        "perfect"
        " phase shifter for rev B's sixteen delays",
    )
    _save(fig, out, "alphabets")


def figure_sequences(res: dict[str, Any], out: Path) -> None:
    """Figure 11: state rasters, line weights and line spectra of a staircase, Frank shifts and the Hadamard words."""
    plt = _pyplot()
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    states_map = ListedColormap(list(SERIES[:4]))
    fig, axes = plt.subplots(
        3, 3, figsize=(12.5, 10.6), gridspec_kw={"width_ratios": [1, 1.1, 1.15]}, layout="constrained"
    )
    lines = np.asarray(res["lines"])
    names = ("staircase", "Frank shifts", "Hadamard words")
    notes = (
        "four steps, shifted to steer line 1 to 20°",
        "i^((l div 4)(l mod 4)), delayed n slots",
        "column n of H16, as 0 and 180°",
    )
    for r, (name, note) in enumerate(zip(names, notes, strict=True)):
        row = res[name]
        raster, weights, spectrum = axes[r]
        raster.imshow(row["states"], cmap=states_map, vmin=-0.5, vmax=3.5, aspect="auto", interpolation="nearest")
        raster.set_ylabel("element n")
        raster.set_title(f"{name}: {note}", fontsize=8.6)
        raster.set_xticks([0, 4, 8, 12, 15])
        raster.grid(False)
        c = np.asarray(row["weights"]).T
        scale = 0.9 / float(np.abs(c).max())
        hh, nn = np.meshgrid(lines, np.arange(c.shape[0]))
        weights.quiver(
            hh,
            nn,
            (c * scale).real,
            -(c * scale).imag,  # the element axis runs downward: minus keeps +Im pointing up
            angles="xy",
            scale_units="xy",
            scale=1,
            pivot="middle",
            color=SERIES[0],
            width=0.006,
            headwidth=3.5,
            headlength=3.5,
            headaxislength=3,
        )
        weights.set_xlim(lines[0] - 0.7, lines[-1] + 0.7)
        weights.set_ylim(c.shape[0] - 0.4, -0.6)
        weights.set_ylabel("element n")
        weights.set_xticks([-8, -4, 0, 4, 8])
        weights.set_title(f"c_(n,h) as arrows, the longest {float(np.abs(c).max()):.3f}", fontsize=8.6)
        power = np.asarray(row["power_per_line"])
        shown = np.maximum(db(power), -40)
        spectrum.bar(lines, shown + 40, bottom=-40, width=0.7, color=SERIES[0])
        spectrum.set_ylim(-40, 2)
        spectrum.set_ylabel("power per element (dB)")
        spectrum.set_xticks([-8, -4, 0, 4, 8])
        distinct = float(np.sum(power[(lines >= -8) & (lines <= 7)]))
        best = int(np.argmax(shown))
        spectrum.set_title(
            f"line {lines[best]:+d}: {shown[best]:.2f} dB; lines −8 to 7 hold {distinct:.3f} in all", fontsize=8.6
        )
    axes[-1, 0].set_xlabel("slot l")
    axes[-1, 1].set_xlabel("comb line h")
    axes[-1, 2].set_xlabel("comb line h")
    handles = [Patch(color=SERIES[s], label=f"state {90 * s}°") for s in range(4)]
    fig.legend(handles=handles, loc="outside lower left", ncol=4, fontsize=8.3)
    _suptitle(
        fig,
        f"On comb line h the board is a phased array with weights c_(n,h) = sinc(h/L) e^(−iπh/L) DFT_h(g_n): the "
        f"staircase steers line 1 at {res['staircase_beam_db']:.2f} dB of 16 ideal elements; Frank shifts give each "
        f"line a linear phase across the elements, one DFT beam per line; the Hadamard words decode with condition "
        f"number {res['hadamard_condition']:.2f}",
    )
    _save(fig, out, "sequences")


def figure_lines(res: dict[str, Any], out: Path) -> None:
    """Figure 12: recovery of every element's signal by slot means and by comb lines, against bandwidth."""
    plt = _pyplot()
    fig, ax = plt.subplots(figsize=(10, 4.9), layout="constrained")
    band = np.asarray(res["bandwidth_Hz"], dtype=np.float64)
    x = np.where(band > 0, band, 10.0)
    ax.semilogx(x, res["slot_means_db"], "-o", color=SERIES[1], ms=5, label="slot means")
    ax.semilogx(x, res["lines_db"], "-o", color=SERIES[0], ms=5, label="comb lines")
    for key, fp in res["f_p_Hz"].items():
        board, length = key.split(" L=")
        style = "-" if board == "rev A" else "--"
        ax.axvline(fp, color=MUTED, lw=0.9, ls=style)
        ax.text(fp * 1.025, -73, f"{board}, L = {length}: {fp / 1e3:.2f} kHz", rotation=90, fontsize=7.2, color=INK2)
    ax.axvline(res["f_p_simulated_Hz"], color=INK, lw=1.4)
    ax.text(
        res["f_p_simulated_Hz"] / 1.03, -9, "f_p of the capture,\n8 µs × 17 slots", fontsize=7.8, ha="right", va="top"
    )
    _log_ticks(ax.xaxis, [10, 100, 1e3, 1e4], ["a tone", "100 Hz", "1 kHz", "10 kHz"])
    ax.set_xlim(7, 3.2e4)
    ax.set_ylim(-90, -5)
    ax.set_xlabel(
        "two-sided signal bandwidth B; grey lines, the comb spacing f_p = 1/(Lτ), rev A (10 µs) solid, rev B (15 µs) "
        "dashed"
    )
    ax.set_ylabel("error of the recovered array vector (dB)")
    ax.legend(fontsize=8.5, loc="center left")
    i = int(np.argmin(np.abs(band - 6000)))
    _suptitle(
        fig,
        f"One channel, sixteen elements, Hadamard words behind a marker at 40 dB SNR: at {band[i] / 1e3:.0f} kHz the "
        f"slot means fall to {res['slot_means_db'][i]:.1f} dB while the comb lines hold {res['lines_db'][i]:.1f} dB; "
        "one-channel beamforming needs B < f_p, not a steady tone",
        width=125,
    )
    _save(fig, out, "lines")


def figure_nulls(res: dict[str, Any], out: Path) -> None:
    """Figure 13: line-1 null depth against sequence length, and a realised pattern against its target."""
    plt = _pyplot()
    fig = plt.figure(figsize=(12.5, 5.0), layout="constrained")
    grid = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 0.6])
    left = fig.add_subplot(grid[0, 0])
    lengths = np.asarray(res["lengths"])
    for key, values in res["curves"].items():
        colour = SERIES[0] if key.startswith("2") else SERIES[1]
        style = "-o" if key.endswith("feedback") else "--s"
        left.plot(lengths, values, style, color=colour, ms=4.5, label=key)
    for bits, colour, (x, dy) in (("2-bit", SERIES[0], (11.0, 0.8)), ("4-bit", SERIES[1], (90.0, -2.6))):
        left.axhline(res["static"][bits], color=colour, lw=1.0, ls=":")
        left.text(x, res["static"][bits] + dy, f"{bits}, one static state", fontsize=7.6, color=INK2)
    for key, value in res["floors"].items():
        left.axhspan(value - 0.7, value + 0.7, color=GRID, lw=0)
        left.text(lengths[-1], value + 1.0, f"σ²/N at {key}", fontsize=7.6, color=INK2, ha="right")
    left.axvline(256, color=INK, lw=1.0)
    left.text(256 / 1.06, -85, "256 words: what plays today", fontsize=7.6, ha="right", va="bottom")
    left.set_xscale("log", base=2)
    _log_ticks(left.xaxis, [float(v) for v in lengths], [str(int(v)) for v in lengths])
    left.set_xlabel("slots per period L")
    left.set_ylabel(f"median null depth, {res['pairs']} beam-and-null pairs (dB)")
    left.set_ylim(-86, -18)
    left.legend(fontsize=7.8, loc="lower left")
    left.set_title("rev A's panel at 5.8 GHz: error feedback (solid) and search (dashed)", fontsize=9)
    theta = np.asarray(res["theta_deg"])
    whole = fig.add_subplot(grid[0, 1])
    zoom = fig.add_subplot(grid[0, 2])
    for i, (name, pattern) in enumerate(res["pattern"].items()):
        target = name == "target"
        colour = INK if target else SERIES[i - 1]
        label = "target" if target else f"{name}: {res['pattern_null_db'][name]:.1f} dB"
        for ax in (whole, zoom):
            ax.plot(theta, pattern, ":" if target else "-", color=colour, lw=1.6 if target else 1.3, label=label)
    whole.axvspan(25, 35, color=GRID, lw=0)
    whole.set_xlim(-90, 90)
    whole.set_ylim(-60, 3)
    whole.set_xlabel("angle in the x-z plane (degrees)")
    whole.set_ylabel("pattern, dB of the beam's peak")
    whole.set_title("a broadside beam with a null at 30°, realised on line 1", fontsize=9)
    whole.legend(fontsize=7.6, loc="lower left")
    zoom.set_xlim(25, 35)
    zoom.set_ylim(-80, -10)
    zoom.set_xlabel("degrees")
    zoom.set_title("the null, 25 to 35°", fontsize=9)
    curve = res["curves"]["2-bit search"]
    at = {int(n): v for n, v in zip(lengths, curve, strict=True)}
    _suptitle(
        fig,
        f"Nulls from 2-bit states deepen with the sequence: search reaches {at.get(16, curve[0]):.1f} dB at L = 16 and "
        f"{at.get(64, curve[-1]):.1f} dB at L = 64, where one static state gives {res['static']['2-bit']:.1f} dB; a "
        f"state table known to 3° sets {res['floors']['3 deg']:.1f} dB, so the calibration, not the sequence, sets the "
        "null on hardware",
    )
    _save(fig, out, "nulls")


def figure_dark(res: dict[str, Any], out: Path) -> None:
    """Figure 14: the tile's drive space: eigenvalues of R across the band, field cuts of uniform drive and of the
    darkest delay-state drive on one scale, the carpet above an element, and how dark each alphabet gets."""
    plt = _pyplot()
    fig = plt.figure(figsize=(12.5, 9.0), layout="constrained")
    grid = fig.add_gridspec(2, 3, height_ratios=[1, 0.9])
    eig = fig.add_subplot(grid[0, 0])
    freqs = np.asarray(res["freqs_GHz"])
    for n in range(16):
        eig.semilogx(freqs, np.asarray(res["free_db"])[:, n], color=MUTED, lw=0.7, alpha=0.8)
        eig.semilogx(freqs, np.asarray(res["over_db"])[:, n], color=SERIES[0], lw=1.0)
    eig.plot([], [], color=SERIES[0], label="over ground")
    eig.plot([], [], color=MUTED, label="free space")
    for level, text in ((-10, "10 %"), (-20, "1 %")):
        eig.axhline(level, color=SERIES[1], lw=1.0, ls="--")
        eig.text(0.42, level + 1.5, text, fontsize=7.8, color=INK2)
    eig.set_ylim(-90, 3)
    eig.set_xlim(freqs[0], freqs[-1])
    _log_ticks(eig.xaxis, [0.4, 1, 2, 3, 6])
    eig.set_xlabel("frequency (GHz)")
    eig.set_ylabel("eigenvalues of R, dB of the largest")
    eig.set_title("sixteen drives: radiating above the lines, reactive below", fontsize=9)
    eig.legend(fontsize=7.8, loc="lower right")
    xs, zs = np.asarray(res["cut_x"]) * 1e3, np.asarray(res["cut_z"]) * 1e3
    radiated = res["radiated_db"]
    cuts = []
    mesh = None
    for j, name in enumerate(("uniform drive", "darkest delay states")):
        ax = fig.add_subplot(grid[0, 1 + j])
        cuts.append(ax)
        mesh = ax.pcolormesh(xs, zs, res["cuts_db"][name], cmap=_ramp(), vmin=-50, vmax=0, shading="auto")
        ax.plot((np.arange(4) - 1.5) * 25, np.full(4, 13), "s", mfc="none", mec=SERIES[1], ms=6, mew=1.4)
        ax.set_aspect("equal")
        ax.grid(False)
        ax.set_xlabel("x (mm), midway between two rows")
        told = "the reference" if name == "uniform drive" else f"radiates {radiated[name]:.1f} dB against it"
        ax.set_title(f"{name} at 2.44 GHz: {told}", fontsize=9)
    cuts[0].set_ylabel("height above the ground (mm)")
    cbar = fig.colorbar(mesh, ax=cuts, location="bottom", fraction=0.05, shrink=0.7)
    cbar.set_label("|E|², dB of the uniform drive's largest in the cut; squares, the element columns 12.5 mm away")
    above = fig.add_subplot(grid[1, 0])
    heights = np.asarray(res["heights_mm"])
    ref = np.asarray(res["profiles"]["uniform drive"])
    for i, name in enumerate(("darkest delay states", "best sign pattern", "checkerboard")):
        above.plot(heights, db(np.asarray(res["profiles"][name]) / ref), color=SERIES[i + 1], label=name)
    kappa = float(res["kappa_per_m"])
    near = (heights >= 12) & (heights <= 42)
    checker = db(np.asarray(res["profiles"]["checkerboard"]) / ref)
    start = int(np.argmin(np.abs(heights - 20)))
    above.plot(
        heights[near],
        checker[start] - 20 * kappa * (heights[near] - heights[start]) * 1e-3 / math.log(10),
        ":",
        color=INK,
        lw=1.3,
        label=f"e^(−2κh), κ = {kappa:.0f} per metre",
    )
    above.axhline(0, color=MUTED, lw=0.8)
    above.set_xlim(0, 60)
    above.set_ylim(-45, 8)
    above.set_xlabel("height above an element (mm)")
    above.set_ylabel("|E|² against uniform drive (dB)")
    above.set_title("the carpet: the same field near the board, falling away above it", fontsize=9)
    above.legend(fontsize=7.6, loc="lower left")
    table = fig.add_subplot(grid[1, 1:])
    groups = [k for k in res["table"] if k.startswith("over ground")] + [
        k for k in res["table"] if k.startswith("free space")
    ]
    positions = np.arange(len(groups))
    for i, (key, label) in enumerate(
        (
            ("sign_pattern_db", "best sign pattern"),
            ("delay_states_db", "rev B's delay states"),
            ("continuous_db", "continuous phases"),
        )
    ):
        values = [res["table"][g][key] for g in groups]
        table.plot(positions + (i - 1) * 0.16, values, "o", color=SERIES[i], ms=7, label=label)
    table.set_xticks(positions)
    table.set_xticklabels([g.replace(" ", "\n", 2).replace("\n", " ", 1) for g in groups], fontsize=8)
    table.set_ylabel("radiated power against uniform drive (dB)")
    table.set_title("unit-amplitude drives that radiate least, tile model", fontsize=9)
    table.legend(fontsize=7.8, loc="lower right")
    dark = res["table"]["over ground 2.44 GHz"]
    _suptitle(
        fig,
        f"The tile's reactive drives: at 2.44 GHz over ground {dark['above_10pct']} of its 16 drives radiate above "
        f"10 % of the strongest and {dark['below_1pct']} below 1 %; rev B's delay states reach "
        f"{dark['delay_states_db']:.1f} "
        f"dB against uniform drive (the best sign pattern {dark['sign_pattern_db']:.1f} dB) and keep the field of "
        "uniform drive within a few millimetres of the board",
    )
    _save(fig, out, "dark")


def figure_selection(res: dict[str, Any], out: Path) -> None:
    """Figure 15: the ring selection rule, the share F_l(kR) of each orbital order, and the radiate-or-hold trade."""
    plt = _pyplot()
    fig = plt.figure(figsize=(12.5, 8.6), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 0.95], width_ratios=[1.35, 1])
    chart = fig.add_subplot(grid[0, :])
    split = np.array([[1.0 if res["chart"][(m, n)]["split"] else 0.0 for n in range(1, 9)] for m in range(7)])
    chart.imshow(split, cmap=_ramp(), vmin=0, vmax=3.2, aspect="auto", extent=(0.5, 8.5, 6.5, -0.5))
    for m in range(7):
        for n in range(1, 9):
            orders = [o for o in res["chart"][(m, n)]["orders"] if abs(o) <= 4]
            text = ", ".join(f"{o:+d}" if o else "0" for o in orders)
            if n == 1:
                text = "every l"
            elif n == 2:
                text = "even l" if m % 2 == 0 else "odd l"
            chart.text(n, m, text, ha="center", va="center", fontsize=7.8)
    chart.set_xticks(range(1, 9))
    chart.set_yticks(range(7))
    chart.set_xlabel("perturbations N on the ring")
    chart.set_ylabel("guide mode m")
    chart.grid(False)
    chart.set_title(
        "orbital orders l = m − pN radiated (those with |l| ≤ 4); shaded where N divides 2m, where +m and −m couple "
        "and split into standing waves",
        fontsize=9,
    )
    share = fig.add_subplot(grid[1, 0])
    x = np.asarray(res["x"])
    for order in range(7):
        share.plot(x, np.asarray(res["F"])[order], color=SERIES[order], lw=1.5, label=f"l = {order}")
    for k, (name, kr) in enumerate(sorted(res["marks"].items(), key=lambda item: item[1])):
        share.axvline(kr, color=MUTED, lw=0.7, ls=":")
        right = kr > 9
        share.text(
            kr - 0.08 if right else kr + 0.06,
            0.97 - 0.07 * (k % 7),
            name,
            fontsize=6.8,
            color=INK2,
            va="top",
            ha="right" if right else "left",
        )
    share.set_xlim(0, 12)
    share.set_ylim(0, 1.0)
    share.set_xlabel("kR")
    share.set_ylabel("F_l(kR): the share of order l")
    share.set_title("each order radiates once kR passes about |l|; dotted, array16's rings", fontsize=9)
    share.legend(fontsize=7.6, ncol=2, loc="upper right")
    trade = fig.add_subplot(grid[1, 1])
    eta = np.linspace(0, 1, 50)
    trade.plot(eta, 1 - eta, color=INK, lw=1.2, label="U/U₀ = 1 − η")
    for i, (name, row) in enumerate(res["stubs"].items()):
        trade.plot(
            [row["eta"]], [1 - row["eta"]], "o", color=SERIES[i], ms=8, label=f"{name}: Q_rad {sci(row['Q_rad'])}"
        )
    trade.set_xlabel("fraction radiated η = Q_ohm / (Q_ohm + Q_rad)")
    trade.set_ylabel("well kept, U/U₀")
    trade.set_title(f"radiate or hold, at Q_ohm = {Q_OHM:.0f}: the 277 mm ring at 13 V", fontsize=9)
    trade.legend(fontsize=7.4, loc="upper right")
    stubs = res["stubs"]
    _suptitle(
        fig,
        f"The guide as an antenna with a dial: N stubs on a ring carrying mode m radiate only the orders l = m − pN, "
        f"each with the share F_l(kR); 3 stubs at m = 2 radiate {stubs['3 stubs at m = 2']['eta'] * 100:.0f} % of the "
        f"drive and 8 stubs at m = 5 {stubs['8 stubs at m = 5']['eta'] * 100:.0f} %, and the well keeps what is not "
        "radiated",
    )
    _save(fig, out, "selection")


def figure_containment(res: dict[str, Any], out: Path) -> None:
    """Figure 16: the power an electron well of 1 eV takes, radiated freely against stored behind electrodes, and
    where the null's Mathieu q lets it stand; the facing boards' field and one board's."""
    plt = _pyplot()
    fig = plt.figure(figsize=(12.5, 7.6), layout="constrained")
    grid = fig.add_gridspec(2, 2, width_ratios=[1.05, 1], height_ratios=[0.8, 1])
    ladder = fig.add_subplot(grid[:, 0])
    freqs = np.asarray(res["freqs_GHz"])
    colours = {"four rods": SERIES[0], "two facing boards": SERIES[1], "one board, two rails": SERIES[2]}
    for name, values in res["ladder"].items():
        geometry, resistance = name.rsplit(", ", 1)
        style = "-" if resistance.startswith("10") else "--"
        power = 1.0 / np.asarray(values)
        cut = int(np.argmax(np.asarray(res["q_for_1eV"][geometry]) <= 0.908))  # q falls as the frequency rises
        ladder.loglog(freqs[: cut + 1], power[: cut + 1], style, color=colours[geometry], lw=0.9, alpha=0.35)
        ladder.loglog(
            freqs[cut:], power[cut:], style, color=colours[geometry], lw=1.7, label=geometry if style == "-" else None
        )
    offsets = {"four rods": (8, -13), "two facing boards": (8, -13), "one board, two rails": (-108, 7)}
    for geometry, row in res["at_1GHz"].items():
        watts = row["q_half_mW"] / 1e3
        ladder.plot([row["q_half_GHz"]], [watts], "o", color=colours[geometry], ms=6)
        ladder.annotate(
            f"q = 0.5: {row['q_half_V']:.0f} V, " + (f"{watts * 1e3:.0f} mW" if watts < 1 else f"{watts:.1f} W"),
            (row["q_half_GHz"], watts),
            offsets[geometry],
            textcoords="offset points",
            fontsize=7.4,
            color=INK2,
        )
    ladder.axhline(1 / res["ceiling_eV_per_W"], color=INK, lw=1.4, label="carried in from far away (RFG-2)")
    short = {
        "16-element shell, kb = pi": ("shell of 16, λ/2", (-66, 5)),
        "whips 0.433 GHz, 100 mm up": ("whips, 100 mm up", (6, -10)),
        "whips 0.9 GHz, 100 mm up": ("whips, 100 mm up", (6, 4)),
        "tile 2.44 GHz, 50 mm out": ("tile, 50 mm", (6, 4)),
        "tile 5.8 GHz, 150 mm out": ("tile, 150 mm", (-60, 5)),
    }
    for name, row in res["foci"].items():
        label, offset = short.get(name, (name, (4, 4)))
        ladder.plot([row["f_GHz"]], [1 / row["eV_per_W"]], "D", color=INK2, ms=5)
        ladder.annotate(label, (row["f_GHz"], 1 / row["eV_per_W"]), offset, textcoords="offset points", fontsize=7.4)
    ladder.set_xlim(0.4, 6.0)
    _log_ticks(ladder.xaxis, [0.4, 1, 2, 4, 6])
    ladder.set_ylim(1e-3, 3e6)
    ladder.set_xlabel("frequency (GHz)")
    ladder.set_ylabel("power for an electron well 1 eV deep (W)")
    ladder.set_title(
        f"electrodes {res['height_mm']:.1f} mm from the null behind 10 kΩ (solid) or 1 kΩ (dashed); faint where the "
        "null's q\npasses 0.908 and no well stands; dots at q = 0.5; diamonds, foci of the shell, whips and tile",
        fontsize=8.8,
    )
    ladder.legend(fontsize=7.8, loc="upper left")
    h = float(res["height_mm"])
    boards = fig.add_subplot(grid[0, 1])
    bx, by = np.asarray(res["boards_x"]) * 1e3, np.asarray(res["boards_y"]) * 1e3
    e = np.asarray(res["boards_field"])
    magnitude = np.hypot(e[..., 0], e[..., 1]) * (2 * h * 1e-3)
    mesh = boards.pcolormesh(bx, by, magnitude, cmap=_ramp(), vmin=0, vmax=1.6, shading="auto")
    boards.streamplot(bx, by, e[..., 0], e[..., 1], color=INK2, linewidth=0.6, density=0.8, arrowsize=0.6)
    for y, signs in ((h, (SERIES[1], SERIES[0])), (-h, (SERIES[0], SERIES[1]))):
        boards.plot([bx[0], 0], [y, y], color=signs[0], lw=6, solid_capstyle="butt", clip_on=False)
        boards.plot([0, bx[-1]], [y, y], color=signs[1], lw=6, solid_capstyle="butt", clip_on=False)
    boards.plot([0], [0], "+", color=SERIES[1], ms=11, mew=2)
    boards.set_aspect("equal")
    boards.set_xlim(bx[0], bx[-1])
    boards.set_ylim(-h, h)
    boards.grid(False)
    boards.set_xlabel("x (mm)")
    boards.set_ylabel("y (mm)")
    boards.set_title(
        f"two facing boards {2 * h:.1f} mm apart, halves at +V/2 (orange) and −V/2 (blue):\na null at the centre, "
        f"a barrier V/2h, {res['boards_share']:.2f} of a four-rod well",
        fontsize=8.8,
    )
    cbar = fig.colorbar(mesh, ax=boards, location="bottom", fraction=0.08, shrink=0.7)
    cbar.set_label("|E| / (V/2h)")
    rails = fig.add_subplot(grid[1, 1])
    rx, ry = np.asarray(res["rails_x"]) * 1e3, np.asarray(res["rails_y"]) * 1e3
    er = np.asarray(res["rails_field"])
    rel = db((er[..., 0] ** 2 + er[..., 1] ** 2) * (h * 1e-3) ** 2)
    mesh = rails.pcolormesh(rx, ry, rel, cmap=_ramp(), vmin=-40, vmax=0, shading="auto")
    rails.streamplot(rx, ry, er[..., 0], er[..., 1], color=INK2, linewidth=0.5, density=0.8, arrowsize=0.5)
    geo = res["rails"]
    c1, c2 = geo["gap_mm"] / 2, geo["gap_mm"] / 2 + geo["rail_mm"]
    for s in (-1, 1):
        rails.plot([s * c1, s * c2], [0, 0], color=SERIES[1], lw=7, solid_capstyle="butt", clip_on=False)
    rails.plot([0], [geo["null_mm"]], "+", color=SERIES[1], ms=11, mew=2)
    rails.plot([0], [geo["saddle_mm"]], "x", color=INK, ms=9, mew=2)
    rails.set_aspect("equal")
    rails.set_xlim(rx[0], rx[-1])
    rails.set_ylim(0, ry[-1])
    rails.grid(False)
    rails.set_xlabel("x (mm)")
    rails.set_ylabel("height above the board (mm)")
    rails.set_title(
        f"one board: rails {geo['rail_mm']:.2f} mm wide either side of a {geo['gap_mm']:.2f} mm grounded strip;\n"
        f"the null (+) {geo['null_mm']:.2f} mm up, the saddle (×) {geo['saddle_mm']:.2f} mm up; {geo['share']:.4f} of "
        "a four-rod well",
        fontsize=8.8,
    )
    cbar = fig.colorbar(mesh, ax=rails, location="bottom", fraction=0.06, shrink=0.7)
    cbar.set_label("|E|², dB of (V/h)²")
    four, faced, one = (res["at_1GHz"][key] for key in ("four rods", "two facing boards", "one board, two rails"))
    _suptitle(
        fig,
        f"What containment takes: power carried in from far away makes at most {res['ceiling_eV_per_W'] * 1e6:.1f} µeV "
        f"per watt for an electron at any frequency, so 1 eV takes {1e-3 / res['ceiling_eV_per_W']:.1f} kW; behind "
        f"10 kΩ, electrodes {h:.1f} mm from the null make 1 eV at q = 0.5 with {four['q_half_V']:.0f} V and "
        f"{four['q_half_mW']:.0f} mW on four rods at {four['q_half_GHz']:.2f} GHz, {faced['q_half_V']:.0f} V and "
        f"{faced['q_half_mW']:.0f} mW on two facing boards at {faced['q_half_GHz']:.2f} GHz, {one['q_half_V']:.0f} V "
        f"and {one['q_half_mW'] / 1e3:.1f} W on one board at {one['q_half_GHz']:.2f} GHz",
    )
    _save(fig, out, "containment")


def _resistor(ax: Any, x: float, y0: float, y1: float, label: str) -> None:
    """A vertical zigzag resistor from y0 to y1 at x, labelled on its right."""
    ys = np.linspace(y0 + 0.25 * (y1 - y0), y0 + 0.75 * (y1 - y0), 13)
    xs = x + 0.07 * np.array([0, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, 0])
    ax.plot([x, x], [y0, ys[0]], color=INK, lw=1.2)
    ax.plot(xs, ys, color=INK, lw=1.2)
    ax.plot([x, x], [ys[-1], y1], color=INK, lw=1.2)
    ax.text(x + 0.14, (y0 + y1) / 2, label, va="center", fontsize=8)


def _inductor(ax: Any, x: float, y0: float, y1: float, label: str) -> None:
    """A vertical coil of four turns from y0 to y1 at x, drawn as half circles, labelled on its right."""
    lo, hi = y0 + 0.2 * (y1 - y0), y0 + 0.8 * (y1 - y0)
    r = (hi - lo) / 8
    for k in range(4):
        t = np.linspace(-math.pi / 2, math.pi / 2, 40)
        centre = lo + r * (2 * k + 1)
        ax.plot(x + 0.6 * r * np.cos(t) * 1.4, centre + r * np.sin(t), color=INK, lw=1.2)
    ax.plot([x, x], [y0, lo], color=INK, lw=1.2)
    ax.plot([x, x], [hi, y1], color=INK, lw=1.2)
    ax.text(x + 0.18, (y0 + y1) / 2, label, va="center", fontsize=8)


def _capacitor(ax: Any, x: float, y0: float, y1: float, label: str) -> None:
    mid = (y0 + y1) / 2
    ax.plot([x, x], [y0, mid - 0.035], color=INK, lw=1.2)
    ax.plot([x, x], [mid + 0.035, y1], color=INK, lw=1.2)
    for y in (mid - 0.035, mid + 0.035):
        ax.plot([x - 0.1, x + 0.1], [y, y], color=INK, lw=1.8)
    ax.text(x + 0.14, mid, label, va="center", fontsize=8)


def figure_charges(res: dict[str, Any], out: Path) -> None:
    """Figure 17: a held column's signals, radiated against induced; its equivalent circuit across its pickup; and
    the impedance dip it makes there."""
    plt = _pyplot()
    fig = plt.figure(figsize=(12.5, 5.4), layout="constrained")
    grid = fig.add_gridspec(1, 3, width_ratios=[1.4, 0.8, 1.05])
    ladder = fig.add_subplot(grid[0, 0])
    items = sorted(res["signals_dBm"].items(), key=lambda kv: kv[1])
    lo, hi = res["hackrf_floor_dBm"]
    ladder.axvspan(lo, hi, color=GRID, lw=0)
    ladder.axvline(res["thermal_dBm"], color=INK, lw=1.0, ls="--")
    for i, (name, value) in enumerate(items):
        colour = SERIES[1] if "radiated" in name else SERIES[0]
        ladder.plot([value], [i], "o", color=colour, ms=8)
        ladder.text(value + 6, i, f"{value:.0f} dBm", va="center", fontsize=7.8)
    ladder.set_yticks(range(len(items)))
    ladder.set_yticklabels([name.replace(": ", ":\n").replace("Ohm", "Ω") for name, _ in items], fontsize=7.4)
    ladder.text(res["thermal_dBm"] - 2, -0.9, "k_BT in 1 Hz", ha="right", fontsize=7.6)
    ladder.text(hi + 2, -0.9, "the HackRF's floor in 1 Hz", fontsize=7.6, color=INK2)
    ladder.plot([], [], "o", color=SERIES[0], label="induced on the electrodes")
    ladder.plot([], [], "o", color=SERIES[1], label="radiated")
    ladder.legend(fontsize=7.8, loc="lower right")
    ladder.set_xlim(-245, -10)
    ladder.set_ylim(-1.3, len(items) - 0.4)
    ladder.set_xlabel("power (dBm)")
    ladder.set_title(
        f"{sci(res['count'])} electrons of the 1 GHz, q = 0.3 column\n(the limiting density in a radius of "
        f"{res['column_radius_m'] * 1e3:.2f} mm)\nmoving 0.1 mm at {res['f_sec_Hz'] / 1e6:.0f} MHz, a 2 mm electrode "
        "pair",
        fontsize=8.8,
    )
    circuit = fig.add_subplot(grid[0, 1])
    circuit.axis("off")
    circuit.set_xlim(0, 2)
    circuit.set_ylim(0, 3)
    column = res["columns"][f"{res['count']:.2g}"]
    top, bottom = 2.6, 0.45
    circuit.plot([0.45, 1.45], [top, top], color=INK, lw=1.2)
    circuit.plot([0.45, 1.45], [bottom, bottom], color=INK, lw=1.2)
    _inductor(circuit, 0.45, 1.85, top, f"l_N = md²/(Ne²)\n= {column['L_H'] * 1e6:.1f} µH")
    _capacitor(circuit, 0.45, 1.15, 1.85, f"c_N = {column['C_F'] * 1e12:.2f} pF")
    _resistor(circuit, 0.45, bottom, 1.15, f"radiation:\nQ {sci(column['radiation_Q'])}")
    _resistor(circuit, 1.45, bottom, top, "50 Ω\npickup")
    circuit.text(
        1.0,
        0.12,
        f"resonant at {res['f_sec_Hz'] / 1e6:.0f} MHz; into 50 Ω it rings\ndown at {sci(column['gamma'])} per second, "
        f"Q {column['Q']:.0f}",
        ha="center",
        fontsize=7.8,
    )
    circuit.set_title("the column across its electrodes", fontsize=8.8)
    dip = fig.add_subplot(grid[0, 2])
    offsets = np.asarray(res["offsets_Hz"])
    for i, row in enumerate(res["columns"].values()):
        dip.semilogx(offsets, row["Z_over_R"], color=SERIES[i], label=f"{sci(row['N'])} electrons, Q {sci(row['Q'])}")
    _log_ticks(dip.xaxis, [10, 1e3, 1e5, 1e7], ["10 Hz", "1 kHz", "100 kHz", "10 MHz"])
    dip.set_xlabel(f"offset above f_sec = {res['f_sec_Hz'] / 1e6:.0f} MHz")
    dip.set_ylabel("|Z| / 50 Ω at the pickup")
    dip.set_ylim(0, 1.05)
    dip.legend(fontsize=7.8, loc="upper left")
    dip.set_title("the impedance dip: its width counts the electrons", fontsize=8.8)
    sig = res["signals_dBm"]
    gap = sig["column, coherent: induced into 50 Ohm"] - sig["column, coherent: radiated"]
    _suptitle(
        fig,
        f"Charges are antennas whose port is the electrode: the column's coherent motion is heard on its electrodes "
        f"{gap:.0f} dB above what it radiates, one electron reaches a 100 kΩ resonant pickup "
        f"{sig['one electron: induced into a 100 kOhm resonant pickup'] - res['thermal_dBm']:.0f} dB above k_BT in "
        "1 Hz, "
        f"and seen from a 2 mm pair one electron is {res['one_electron_H']:.0f} H",
    )
    _save(fig, out, "charges")


def figure_modules(results: dict[str, Any], out: Path) -> None:
    """Figure 18: the compiler's modules and types with the data between them, and the three bench roles."""
    plt = _pyplot()
    from matplotlib.patches import FancyBboxPatch

    fig, ax = plt.subplots(figsize=(12.5, 7.4), layout="constrained")
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 7.4)
    ax.axis("off")
    w, h = 2.0, 1.35
    row1, row2, row3 = 6.3, 3.95, 1.45
    boxes = {
        "region": (
            1.2,
            row1,
            "Region(points (M,3),\ncomponents, weights (M,));\nslab, ball, cone, shell,\ncurve; interleave",
        ),
        "elements": (
            3.6,
            row1,
            "Elements(positions, moments,\ntransfer (S,N), feeds,\nground); field → (M,3,N);\nradiation → (N,N)",
        ),
        "targets": (
            8.4,
            row1,
            "FieldValues | VirtualAntenna\n| Focus | Quiet | Intensity\n| NullCurve | Pattern;\n"
            "Budget(power_cap, loss, …)",
        ),
        "states": (10.8, row1, "Alphabet(gains (N,S), f,\nshared); rev_a, rev_b;\nceiling, perimeter / 2π;\nnearest"),
        "measure": (
            1.2,
            row2,
            "Design(states (K,N));\nhadamard, frank,\ncommutation; decode → h;\nkaczmarz; fit_gains → D",
        ),
        "operator": (
            3.6,
            row2,
            "Operator(G (M·C,N),\nR (N,N), f, region,\nprovenance); model,\nmeasured, radiation_from_*",
        ),
        "spectrum": (6.0, row2, "Spectrum(σ², V, U)\n.count(within_db);\nspectrum(op, loss=0.1)"),
        "solve": (
            8.4,
            row2,
            "Solution(w, residual_db,\npower, sensitivity_db, α);\nlinear, focus, intensity,\nnull_curve, eirp",
        ),
        "sequences": (
            10.8,
            row2,
            "Schedule(states (L,N),\nslot_s, marker);\nline_weights, compile_line,\ncompile_average, to_protocol",
        ),
        "compile": (
            1.2,
            row3,
            "Problem → Plan: the\nloop's steps 1 to 9 as\none pure pipeline;\nroute(direction, B, …)",
        ),
        "verify": (3.6, row3, "Score(fit_db, held_db,\npower_ratio,\nsensitivity_db);\nscore(…)"),
    }
    for name, (x, y, body) in boxes.items():
        _box(ax, x, y, w, h, name, body, tone=PALE if name in ("compile", "verify") else BLUES[1])

    def link(a: str, b: str, label: str, start_dx: float = 0.0, end_dx: float = 0.0) -> None:
        xa, ya, _ = boxes[a]
        xb, yb, _ = boxes[b]
        if abs(ya - yb) < 1e-9:
            sign = 1.0 if xb > xa else -1.0
            start, end = (xa + sign * w / 2, ya), (xb - sign * w / 2, yb)
            ax.text((xa + xb) / 2, ya + 0.08, label, fontsize=7.6, color=INK2, ha="center", va="bottom")
        else:
            start, end = (xa + start_dx, ya - h / 2), (xb + end_dx, yb + h / 2)
            ax.text((start[0] + end[0]) / 2 + 0.08, (start[1] + end[1]) / 2, label, fontsize=7.6, color=INK2)
        _arrow(ax, start, end)

    link("region", "operator", "points", start_dx=0.5, end_dx=-0.55)
    link("elements", "operator", "g_n(r), R")
    link("measure", "operator", "h, D")
    link("operator", "spectrum", "G,\nR + rI")
    link("spectrum", "solve", "σᵢ, vᵢ")
    link("targets", "solve", "Target, Budget")
    link("states", "sequences", "Alphabet")
    link("solve", "sequences", "w")
    link("operator", "verify", "held-out G")
    link("verify", "compile", "Score")
    start, end = (boxes["solve"][0] - 0.6, row2 - h / 2), (boxes["verify"][0] + 0.55, row3 + h / 2)
    _arrow(ax, start, end)
    ax.text((start[0] + end[0]) / 2 + 0.05, (start[1] + end[1]) / 2 + 0.2, "w", fontsize=7.6, color=INK2)
    left, bottom, width, height = 6.0, 0.15, 5.95, 2.05
    ax.add_patch(
        FancyBboxPatch(
            (left, bottom), width, height, boxstyle="round,pad=0.004,rounding_size=0.04", fc=PALE, ec=EDGE, lw=0.8
        )
    )
    ax.text(
        left + 0.12,
        bottom + height - 0.1,
        "devices: an experiment names roles, a bench file gives a driver per role",
        fontsize=9.2,
        fontweight="bold",
        va="top",
    )
    rows = (
        ("role", "must be able to", "simulated driver", "hardware driver"),
        (
            "array",
            "tune, play, capture\nor transmit, identity",
            "simbench:SimArray",
            "HackRF Pro (libhackrf);\nthe controller (protocol)",
        ),
        ("probe", "emit and off,\nor read", "simbench:SimBeacon", "a synthesiser on SPI,\nor a log detector"),
        ("positioner", "move_to, where, home", "simbench:SimGantry", "a 3-D printer, G-code"),
    )
    columns = (left + 0.12, left + 1.2, left + 2.75, left + 4.25)
    for r, row in enumerate(rows):
        for c, text in enumerate(row):
            ax.text(
                columns[c],
                bottom + height - 0.42 - r * 0.38,
                text,
                fontsize=7.6,
                va="top",
                fontweight="bold" if r == 0 else "normal",
                linespacing=1.05,
            )
    ax.text(
        left + 0.12,
        bottom + 0.08,
        "the array's captures y return to measure; schedules leave through to_protocol",
        fontsize=7.6,
        color=INK2,
    )
    _arrow(ax, (boxes["sequences"][0], row2 - h / 2), (boxes["sequences"][0], bottom + height))
    ax.text(boxes["sequences"][0] + 0.07, 2.75, "protocol.Sequence", fontsize=7.6, color=INK2)
    ax.set_title(
        "The generic process as code: array16/synthesis/, frozen dataclasses and pure functions, with state only in "
        "the "
        "devices",
        fontsize=11,
    )
    _save(fig, out, "modules")


def figure_bench(res: dict[str, Any], out: Path) -> None:
    """Figure 19: the V1 bench in 3-D, its two signal chains, and the predicted held-out residuals."""
    plt = _pyplot()
    from matplotlib.patches import FancyBboxPatch

    fig = plt.figure(figsize=(12.5, 10.8), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[1.1, 1.0], width_ratios=[1.0, 1.0])
    scene: Any = fig.add_subplot(grid[0, 0], projection="3d")
    lam, half, whip = float(res["wavelength_mm"]), float(res["half_side_mm"]), float(res["whip_mm"])
    plate_x, plate_y = (-200.0, 200.0), (-160.0, 250.0)
    scene.plot(
        [plate_x[0], plate_x[1], plate_x[1], plate_x[0], plate_x[0]],
        [plate_y[0], plate_y[0], plate_y[1], plate_y[1], plate_y[0]],
        [0, 0, 0, 0, 0],
        color=EDGE,
        lw=1.2,
    )
    for cx, cy in np.asarray(res["corners"]) * 1e3:
        scene.plot([cx, cx], [cy, cy], [0, whip], color=INK, lw=2.6)
    colours = {"centre": INK2, "across": SERIES[0], "away": SERIES[1], "toward": SERIES[2]}
    for name, row in res["table"].items():
        x, y = row["xy_mm"]
        kind = name.split()[0]
        scene.plot([x, x], [y, y], [0, whip], color=colours[kind], lw=1.2, alpha=0.85)
        scene.scatter([x], [y], [0], color=colours[kind], s=12)
    fit, held = np.asarray(res["fit_points"]) * 1e3, np.asarray(res["held_points"]) * 1e3
    scene.scatter(fit[:, 0], fit[:, 1], fit[:, 2], s=3, color=MUTED, alpha=0.6, label=f"{len(fit)} fit positions")
    scene.scatter(held[:, 0], held[:, 1], held[:, 2], s=7, marker="x", color=SERIES[1], label=f"{len(held)} held out")
    box_x, box_y, box_z = (-90, 90), (100, 220), (30, 180)
    for xx in box_x:
        for yy in box_y:
            scene.plot([xx, xx], [yy, yy], list(box_z), color=INK2, lw=0.7)
    for zz in box_z:
        scene.plot(
            [box_x[0], box_x[1], box_x[1], box_x[0], box_x[0]],
            [box_y[0], box_y[0], box_y[1], box_y[1], box_y[0]],
            [zz] * 5,
            color=INK2,
            lw=0.7,
        )
    scene.plot([0, 0], [160, 160], [85, 125], color=SERIES[3], lw=4, label="the probe, a 40 mm dipole")
    scene.set_xlabel("x (mm)", fontsize=8)
    scene.set_ylabel("y (mm)", fontsize=8)
    scene.set_zlabel("z (mm)", fontsize=8)
    scene.tick_params(labelsize=7)
    scene.set_box_aspect((400, 410, 190))
    scene.view_init(elev=26, azim=-62)
    scene.legend(fontsize=7.6, loc="upper left")
    scene.set_title(
        f"four {whip:.1f} mm monopoles (black) on a {2 * half:.1f} mm square, λ/4 at 2.44 GHz;\nthe target whips "
        "across (blue), away (orange), toward (aqua); the scanned slab",
        fontsize=8.8,
    )
    curves = fig.add_subplot(grid[0, 1])
    for name in ("across", "away", "toward"):
        row = res["curves"][name]
        beyond = (np.asarray(row["distance_mm"]) - half) / lam
        held_db = np.asarray(row["held_db"])
        curves.plot(beyond, held_db, color=colours[name], label=name)
        curves.fill_between(beyond, held_db - 3, held_db + 3, color=colours[name], alpha=0.12, lw=0)
    for name, row in res["table"].items():
        kind = name.split()[0]
        x = -0.125 if kind == "centre" else float(name.split()[1])
        curves.plot([x], [row["held_db"]], "o", color=colours[kind], ms=6, mec=SURFACE)
    curves.axhline(-24, color=INK, lw=0.9, ls="--")
    curves.text(1.1, -23.3, "V1-4: the centre at −24 dB or better", ha="right", fontsize=7.6)
    curves.axhline(-6, color=INK, lw=0.9, ls=":")
    curves.text(1.1, -5.3, "V1-8: a wavelength across at −6 dB or worse", ha="right", fontsize=7.6)
    curves.axvline(-0.125, color=MUTED, lw=0.8)
    curves.text(-0.115, -48.5, "centre", fontsize=7.6, color=INK2)
    curves.set_xlabel("the target whip's distance beyond the square's edge (wavelengths)")
    curves.set_ylabel("predicted held-out residual (dB)")
    curves.set_ylim(-50, 0)
    curves.legend(fontsize=8, loc="lower right", title="shaded ±3 dB: V1-5", title_fontsize=7.6)
    t = res["table"]
    curves.set_title(
        f"the prediction: the centre {t['centre']['held_db']:.1f} dB; λ/2 beyond the edge away "
        f"{t['away 0.5']['held_db']:.1f} dB, across {t['across 0.5']['held_db']:.1f} dB (V1-6)",
        fontsize=8.8,
    )
    chain = fig.add_subplot(grid[1, :])
    chain.axis("off")
    chain.set_xlim(0, 12.5)
    chain.set_ylim(-0.08, 4.5)

    def block(x: float, y: float, text: str, tone: str = BLUES[1], width: float = 1.45) -> None:
        chain.add_patch(
            FancyBboxPatch(
                (x - width / 2, y - 0.3), width, 0.6, boxstyle="round,pad=0.004,rounding_size=0.05", fc=tone, ec=EDGE
            )
        )
        chain.text(x, y, text, ha="center", va="center", fontsize=7.6, linespacing=1.15)

    def wire(a: tuple[float, float], b: tuple[float, float], dashed: bool = False, style: str = "arc3") -> None:
        chain.annotate(
            "",
            xy=b,
            xytext=a,
            arrowprops={
                "arrowstyle": "-|>",
                "color": INK2,
                "lw": 1.0,
                "ls": "--" if dashed else "-",
                "connectionstyle": style,
                "shrinkA": 0,
                "shrinkB": 0,
            },
        )

    xs = (0.85, 2.5, 4.15, 6.0, 7.75, 9.45, 11.3)
    control, y1, y1b, y2, y2b = 3.85, 2.95, 2.1, 1.05, 0.3
    chain.text(0.1, 4.45, "receive: sessions A, A\u2032, E1 and the targets", fontsize=8.8, fontweight="bold", va="top")
    block(7.75, control, "controller:\n74HC595, trigger", tone=PALE)
    block(11.3, control, "10 MHz out to the\nsynthesiser and controller;\ntrigger in", tone=PALE, width=2.0)
    wire((7.75, control - 0.3), (7.75, y1 + 0.3))
    wire((8.48, control + 0.12), (10.3, control + 0.12))
    wire((10.3, control - 0.12), (8.48, control - 0.12))
    chain.text(9.39, control + 0.17, "trigger", fontsize=7.4, color=INK2, ha="center")
    chain.text(9.39, control - 0.07, "10 MHz", fontsize=7.4, color=INK2, ha="center", va="top")
    for x, text, tone in zip(
        xs,
        (
            "synthesiser on\nthe 10 MHz,\n−20 dBm",
            "two-way\nsplitter",
            "probe dipole,\n40 mm, on the\nprinter head",
            "four whips\non the plate",
            "4-way switch\nwith all-off\n(the marker)",
            "combiner",
            "HackRF Pro\nreceives at\n2.44 GHz",
        ),
        (BLUES[1], BLUES[1], WARM, WARM, BLUES[1], BLUES[1], PALE),
        strict=True,
    ):
        block(x, y1, text, tone)
    block(5.95, y1b, "attenuator; a thin cable in a printed drag loop: the cabled reference", width=4.6)
    for a, b in itertools.pairwise(xs):
        wire((a + 0.73, y1), (b - 0.73, y1), dashed=a == xs[2])
    wire((2.5, y1 - 0.3), (3.65, y1b), style="arc3,rad=0.25")
    wire((8.25, y1b), (9.45, y1 - 0.3), style="arc3,rad=0.3")
    chain.text(0.1, 1.68, "transmit: sessions T, the HackRF the only source", fontsize=8.8, fontweight="bold", va="top")
    for x, text, tone in zip(
        xs,
        (
            "HackRF Pro\ntransmits CW",
            "coupler,\n10 to 20 dB",
            "4-way switch",
            "four whips",
            "probe dipole\non the head",
            "combiner",
            "log detector\n(AD8317) to the\ncontroller's ADC",
        ),
        (PALE, BLUES[1], BLUES[1], WARM, WARM, BLUES[1], BLUES[1]),
        strict=True,
    ):
        block(x, y2, text, tone)
    block(5.95, y2b, "the reference at 0°, 90° or off, up the drag loop", width=4.0)
    for a, b in itertools.pairwise(xs):
        wire((a + 0.73, y2), (b - 0.73, y2), dashed=a == xs[3])
    wire((2.5, y2 - 0.3), (3.95, y2b), style="arc3,rad=0.25")
    wire((7.95, y2b), (9.45, y2 - 0.3), style="arc3,rad=0.3")
    _suptitle(
        fig,
        f"V1, a virtual monopole from four switched monopoles, fitted on {len(fit)} positions and checked on "
        f"{len(held)} others, in both directions: the spectrum in R + rI is "
        + ", ".join(f"{v:.1f}" for v in res["generalised_db"])
        + f" dB, and a target a wavelength across the region ({t['across 1']['held_db']:.1f} dB) is the case that must "
        "fail",
    )
    _save(fig, out, "bench")


# ---- the run ---------------------------------------------------------------------------------------------------------
SECTIONS: dict[str, tuple[Callable[[bool], dict[str, Any]], Callable[[dict[str, Any], Path], None]]] = {
    "price": (compute_price, figure_price),
    "radiation": (compute_radiation, figure_radiation),
    "fine": (compute_fine, figure_fine),
    "ceiling": (compute_ceiling, figure_ceiling),
    "band": (compute_band, figure_band),
    "behind": (compute_behind, figure_behind),
    "focused": (compute_focused, figure_focused),
    "whips": (compute_whips, figure_whips),
    "alphabets": (compute_alphabets, figure_alphabets),
    "sequences": (compute_sequences, figure_sequences),
    "lines": (compute_lines, figure_lines),
    "nulls": (compute_nulls, figure_nulls),
    "dark": (compute_dark, figure_dark),
    "selection": (compute_selection, figure_selection),
    "containment": (compute_containment, figure_containment),
    "charges": (compute_charges, figure_charges),
    "bench": (compute_bench, figure_bench),
}
DIAGRAMS: dict[str, Callable[[dict[str, Any], Path], None]] = {"compiler": figure_compiler, "modules": figure_modules}
ORDER = ("compiler", *SECTIONS, "modules")
KEEP = 64  # results.json keeps arrays and lists of at most this many numbers; the figures hold the rest


def slim(value: Any) -> Any:
    """The results without the large arrays behind the maps and curves, for results.json."""
    if isinstance(value, np.ndarray):
        return value if value.size <= KEEP else None
    if isinstance(value, dict):
        kept = {str(k): slim(v) for k, v in value.items()}
        return {k: v for k, v in kept.items() if v is not None}
    if isinstance(value, (list, tuple)):
        if len(value) > KEEP:
            return None
        return [slim(v) for v in value]
    return value


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument("--quick", action="store_true", help="coarse grids and few draws, into virtual_antenna/quick/")
    parser.add_argument("--only", default="", help="comma-separated sections, of: " + ", ".join(ORDER))
    parser.add_argument("--out", type=Path, default=None, help="where the figures, results.json and run.log go")
    args = parser.parse_args(argv)
    chosen = [name for name in args.only.split(",") if name] or list(ORDER)
    unknown = [name for name in chosen if name not in ORDER]
    if unknown:
        print(f"no section {unknown}; there are {list(ORDER)}", file=sys.stderr)
        return 2
    out: Path = args.out if args.out is not None else (OUT / "quick" if args.quick else OUT)
    out.mkdir(parents=True, exist_ok=True)
    LOG.clear()
    began = time.time()
    say(f"virtual antenna study: {'quick' if args.quick else 'full'} run of {', '.join(chosen)}")
    results: dict[str, Any] = {"selftest": selftest()}
    for name in chosen:
        if name in SECTIONS:
            start = time.time()
            results[name] = SECTIONS[name][0](args.quick)
            say(f"  [{name}: {time.time() - start:.1f} s]")
    say("\n== figures ==")
    for name in chosen:
        if name in DIAGRAMS:
            DIAGRAMS[name](results, out)
        else:
            SECTIONS[name][1](results[name], out)
    say(f"\ndone in {time.time() - began:.0f} s")
    if chosen == list(ORDER):
        (out / "results.json").write_text(json.dumps(lab.plain(slim(results)), indent=1) + "\n")
        (out / "run.log").write_text("\n".join(LOG) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
