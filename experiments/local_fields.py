"""Local 3-D field structures near the array: what 16 (or 32) complex excitations can and cannot shape.

The owner's question: "can you produce fields that twist into each other to produce exotic field shapes that can
operate at a local level?"  In vacuum, fields do not act on each other: Maxwell's equations are linear, so the field
of 16 elements is the SUM  E(r) = sum_n w_n G_n(r)  of 16 fixed patterns G_n, and the only thing we control is the 16
complex numbers w_n.  Everything "exotic" here is interference of those 16 patterns.  This file computes the patterns
exactly (all near-field terms) and asks what the 16 numbers can do within about half a metre:

  1  dof       how many independent field shapes the array can make in a box (singular values of G)
  2  focus     a focal spot at 0.2 / 0.3 / 0.5 m: lateral size, depth, gain over an ordinary beam, V/m per mW
  3  spots     two and three simultaneous spots
  4  quiet     a quiet zone (null) next to a bright point: depth, size, tolerance; rev B's 4-channel constraint
  5  superosc  a spot narrower than the diffraction spot, and what it costs
  6  vortex    a steerable vortex line (phase winding in slices), charge 2, rev B's quantised delays
  7  pol       polarisation singularities (C points) if the panel were dual-polarised
  8  curves    the repo's optimiser (climb / evolve / random_search, finalists re-scored with reselect) placing
               vortex lines on a target curve: tilted line, ring, two linked rings, trefoil knot; 16 and 32 elements;
               checked by tracing the vortex lines on a 3-D grid (separate code from the score)
  9  room      a toy multipath room (image sources): time-reversal focusing against free-space focusing

Models (nothing here is a full-wave result):
  * every element is a set of Hertzian current elements; the field of each is the exact dipole field
        E = (i eta k / 4 pi) e^{ikr}/r { m - n(n.m) + [3n(n.m) - m] (1/(kr)^2 - i/(kr)) },  m = current moment (A m),
    time convention e^{-i omega t}.  Derivation: the textbook field of a dipole p with p = i m / omega.
  * panel: 4x4 at 25 mm pitch in the plane z = 0 (boresight +z), each element a y-directed Hertzian dipole 13 mm
    above an infinite perfect ground (image of opposite sign).  A Vivaldi flare is not a dipole; this keeps the
    polarisation, the position and a broad forward pattern, which is what decides the answers here.
  * whips: 16 quarter-wave monopoles on an infinite ground (z = 0), outer ring of 8 at r = 225 mm, inner ring of 8 at
    112 mm turned 22.5 deg; cosine current in 4 segments plus images.  "In front of the array" = above the plate.
  * no mutual coupling (the excitations are the element currents), no finite-ground diffraction.
  * errors: 3 deg and 0.26 dB rms per element (rev B's calibration residual), redrawn for every generation.

    OMP_NUM_THREADS=2 .env/bin/python experiments/local_fields.py            # about 15 minutes
    OMP_NUM_THREADS=2 .env/bin/python experiments/local_fields.py --quick    # smoke run
    ... --only focus,vortex                                                  # some sections

Outputs: experiments/local_fields/{results.json, run.log, *.png}
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from array16.evolve import climb, evolve, random_search, reselect  # noqa: E402

C0 = 299_792_458.0
ETA = 376.730313668
OUT = Path(__file__).resolve().parent / "local_fields"
SIG_PH, SIG_A = np.deg2rad(3.0), 0.03          # 3 deg, 0.26 dB rms per element
TTD_PS = np.array([38.5, 52.5, 117.5, 139.0])  # rev B delay bits (REDESIGN.md section 3)
RES: dict = {}
LOG = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


def db(x):
    return 10 * np.log10(np.maximum(np.asarray(x, float), 1e-300))


# ---- arrays ---------------------------------------------------------------------------------------------
class Arr:
    """N elements, each a weighted set of Hertzian current elements.  E(points) = G(points) @ w."""

    def __init__(self, name, f, pos, mom, T, epos, ground=True, comp=1):
        self.name, self.f = name, f
        self.k, self.lam = 2 * np.pi * f / C0, C0 / f
        self.pos, self.mom, self.T = np.asarray(pos, float), np.asarray(mom, float), np.asarray(T, float)
        self.N = self.T.shape[1]
        self.epos = np.asarray(epos, float)     # element positions (N, 3)
        self.ground, self.comp = ground, comp   # comp: the field component used as the scalar psi

    def G(self, pts, comp=None, chunk=12000):
        """(P, 3, N) vector transfer matrix, or (P, N) for one Cartesian component."""
        pts = np.atleast_2d(np.asarray(pts, float))
        P = len(pts)
        out = np.empty((P, self.N) if comp is not None else (P, 3, self.N), complex)
        pref = 1j * ETA * self.k / (4 * np.pi)
        for a in range(0, P, chunk):
            d = pts[a:a + chunk, None, :] - self.pos[None]
            r = np.linalg.norm(d, axis=2)
            n = d / r[..., None]
            kr = self.k * r
            nm = np.einsum("psc,sc->ps", n, self.mom)
            g = pref * np.exp(1j * kr) / r
            near = 1 / kr**2 - 1j / kr
            if comp is None:
                E = g[..., None] * ((self.mom[None] - n * nm[..., None])
                                    + (3 * n * nm[..., None] - self.mom[None]) * near[..., None])
                out[a:a + chunk] = np.einsum("psc,sn->pcn", E, self.T)
            else:
                mc, nc = self.mom[None, :, comp], n[..., comp]
                out[a:a + chunk] = (g * ((mc - nc * nm) + (3 * nc * nm - mc) * near)) @ self.T
        return out

    def psi(self, pts):
        return self.G(pts, comp=self.comp)

    def R(self, nth=40, nph=80):
        """Radiated power = w^H R w (W, peak-amplitude phasors), by integrating the far field."""
        x, wg = np.polynomial.legendre.leggauss(nth)
        ct, wt = (0.5 * (x + 1), 0.5 * wg) if self.ground else (x, wg)
        ph = (np.arange(nph) + 0.5) * 2 * np.pi / nph
        CT, PH = np.meshgrid(ct, ph, indexing="ij")
        ST = np.sqrt(1 - CT**2)
        n = np.stack([ST * np.cos(PH), ST * np.sin(PH), CT], -1).reshape(-1, 3)
        w = np.repeat(wt, nph) * 2 * np.pi / nph
        nm = n @ self.mom.T
        Fs = (self.mom[None] - n[:, None, :] * nm[..., None]) * np.exp(-1j * self.k * (n @ self.pos.T))[..., None]
        F = 1j * ETA * self.k / (4 * np.pi) * np.einsum("dsc,sn->dcn", Fs, self.T)
        return np.einsum("d,dcn,dcm->nm", w, F.conj(), F) / (2 * ETA)


def panel(f, nx=4, ny=4, pitch=0.025, h=0.013, pol="y", name=None):
    xs = (np.arange(nx) - (nx - 1) / 2) * pitch
    ys = (np.arange(ny) - (ny - 1) / 2) * pitch
    X, Y = np.meshgrid(xs, ys)
    xy = np.c_[X.ravel(), Y.ravel()]                 # element p = iy * nx + ix
    axes = {"y": [1], "x": [0], "xy": [0, 1]}[pol]
    n = len(xy)
    pos, mom, T = [], [], []
    for q, ax in enumerate(axes):
        for e, (x, y) in enumerate(xy):
            for sgn, z in ((1.0, h), (-1.0, -h)):
                pos.append([x, y, z])
                m = np.zeros(3)
                m[ax] = sgn
                mom.append(m)
                t = np.zeros(n * len(axes))
                t[q * n + e] = 1
                T.append(t)
    epos = np.tile(np.c_[xy, np.full(n, h)], (len(axes), 1))
    a = Arr(name or f"panel {nx}x{ny} {f / 1e9:g} GHz", f, pos, mom, T, epos, ground=True, comp=1)
    a.nx, a.ny = nx, ny
    return a


def whips(f, radii=(0.225, 0.112), nseg=4, sel=None, name=None):
    ang = np.arange(8) * np.pi / 4
    xy = np.r_[np.c_[radii[0] * np.cos(ang), radii[0] * np.sin(ang)],
               np.c_[radii[1] * np.cos(ang + np.pi / 8), radii[1] * np.sin(ang + np.pi / 8)]]
    if sel is not None:
        xy = xy[sel]
    L, k = C0 / f / 4, 2 * np.pi * f / C0
    zs = (np.arange(nseg) + 0.5) * L / nseg
    cur = np.cos(k * zs) * L / nseg
    n = len(xy)
    pos, mom, T = [], [], []
    for e, (x, y) in enumerate(xy):
        for z, c in zip(np.r_[zs, -zs], np.r_[cur, cur]):
            pos.append([x, y, z])
            mom.append([0, 0, c])
            t = np.zeros(n)
            t[e] = 1
            T.append(t)
    return Arr(name or f"{n} whips {f / 1e9:g} GHz", f, pos, mom, T, np.c_[xy, np.zeros(n)], ground=True, comp=2)


def surround(f, n=32, radius=0.25, name=None):
    """Hypothetical: n z-directed Hertzian dipoles spread over a sphere around the origin (free space)."""
    i = np.arange(n) + 0.5
    ct = 1 - 2 * i / n
    ph = np.pi * (1 + 5**0.5) * i
    st = np.sqrt(1 - ct**2)
    p = radius * np.c_[st * np.cos(ph), st * np.sin(ph), ct]
    return Arr(name or f"{n} dipoles on a {radius * 1e3:.0f} mm sphere {f / 1e9:g} GHz", f, p,
               np.tile([0, 0, 1.0], (n, 1)), np.eye(n), p, ground=False, comp=2)


def selftest():
    a = Arr("one dipole", 1e9, [[0, 0, 0]], [[0, 0, 1.0]], [[1.0]], [[0, 0, 0]], ground=False)
    p = a.R()[0, 0].real
    exact = ETA * a.k**2 / (12 * np.pi)
    assert abs(p / exact - 1) < 1e-6, (p, exact)
    # far field on the equator: |E| = eta k / (4 pi r)
    e = a.G([[100.0, 0, 0]])[0, :, 0]
    assert abs(abs(e[2]) / (ETA * a.k / (4 * np.pi * 100)) - 1) < 1e-3 and abs(e[0]) < 1e-6 * abs(e[2])
    # near field on the axis: quasi-static dipole, E_z = 2 p / (4 pi eps0 r^3) = eta m * 2 / (4 pi k r^3)
    r = 0.002
    e = a.G([[0, 0, r]])[0, :, 0]
    assert abs(abs(e[2]) / (ETA * 2 / (4 * np.pi * a.k * r**3)) - 1) < 1e-2
    say("selftest ok: radiated power of one dipole within 1e-6 of eta k^2 |Il|^2 / 12 pi; far and near limits ok")
    # vortex tracer on a field whose zeros are known: two linked rings, psi = ring1 * ring2
    g = np.linspace(-2, 2, 81) + 0.013
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    psi = (((X + 0.5)**2 + Y**2 - 1) + 1j * Z) * (((X - 0.5)**2 + Z**2 - 1) + 1j * Y)
    lines = trace_lines(vortex_faces(psi), [g, g, g])
    assert len(lines) == 2 and all(c for _, c in lines), [(len(q), c) for q, c in lines]
    lk = linking_number(lines[0][0], lines[1][0])
    assert abs(abs(lk) - 1) < 0.01, lk
    say(f"selftest ok: tracer finds the two closed loops of a synthetic Hopf link, linking number {lk:+.3f}")


# ---- helpers --------------------------------------------------------------------------------------------
def unit(w):
    return w / np.linalg.norm(w)


def focus_w(arr, rf):
    """Unit-norm excitation that maximises |E|^2 at rf (first right singular vector of the 3 x N matrix)."""
    Gf = arr.G(np.asarray(rf)[None])[0]
    _, s, Vh = np.linalg.svd(Gf)
    return Vh[0].conj(), s[0]**2


def errors(rng, N, n=1):
    return (1 + SIG_A * rng.standard_normal((n, N))) * np.exp(1j * SIG_PH * rng.standard_normal((n, N)))


def fwhm(x, I, i0, win):
    """Half-power full width of the local maximum of I nearest index i0 (searching +-win samples).
    Returns (width, x_peak) or (None, None) if I has no local maximum there."""
    lo, hi = max(i0 - win, 0), min(i0 + win + 1, len(I))
    j = lo + int(np.argmax(I[lo:hi]))
    if j in (lo, hi - 1) and not (j == i0):
        if j == 0 or j == len(I) - 1 or not (I[j] >= I[j - 1] and I[j] >= I[j + 1]):
            return None, None
    half = I[j] / 2
    a = j
    while a > 0 and I[a] > half:
        a -= 1
    b = j
    while b < len(I) - 1 and I[b] > half:
        b += 1
    if I[a] > half or I[b] > half:
        return None, float(x[j])
    xa = x[a] + (half - I[a]) / (I[a + 1] - I[a]) * (x[a + 1] - x[a])
    xb = x[b - 1] + (half - I[b - 1]) / (I[b] - I[b - 1]) * (x[b] - x[b - 1])
    return float(xb - xa), float(x[j])


def eigh(A, B):
    """Generalised Hermitian eigenproblem A v = lambda B v (B positive definite); eigenvalues ascending.
    (numpy only: the environment has no scipy.)"""
    Lc = np.linalg.cholesky(B)
    Li = np.linalg.inv(Lc)
    M = Li @ A @ Li.conj().T
    val, U = np.linalg.eigh(0.5 * (M + M.conj().T))
    return val, Li.conj().T @ U


class cKDTree:
    """Brute-force nearest neighbour with the scipy cKDTree.query interface (k = 1)."""

    def __init__(self, pts):
        self.p = np.asarray(pts, float)

    def query(self, q, chunk=2000):
        q = np.atleast_2d(np.asarray(q, float))
        d, i = np.empty(len(q)), np.empty(len(q), int)
        for a in range(0, len(q), chunk):
            D = np.linalg.norm(q[a:a + chunk, None, :] - self.p[None], axis=2)
            i[a:a + chunk] = D.argmin(1)
            d[a:a + chunk] = D.min(1)
        return d, i


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def winding2d(F):
    """Topological charge on every plaquette of a 2-D complex array F[i, j]; positive about e_i x e_j."""
    ph = np.angle(F)
    di, dj = wrap(ph[1:] - ph[:-1]), wrap(ph[:, 1:] - ph[:, :-1])
    return np.rint((di[:, :-1] + dj[1:, :] - di[:, 1:] - dj[:-1, :]) / (2 * np.pi)).astype(int)


def intensity(arr, pts, w):
    E = arr.G(pts) @ w
    return (np.abs(E)**2).sum(-1)


# ---- 1. degrees of freedom ------------------------------------------------------------------------------
def sec_dof(cases):
    say("\n== 1. degrees of freedom: singular values of the map (excitations -> scalar field in a 4-wavelength cube)")
    rows = []
    for arr, centre in cases:
        g = (np.arange(17) - 8) * arr.lam / 4
        X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
        pts = np.c_[X.ravel(), Y.ravel(), Z.ravel()] + centre
        s = np.linalg.svd(arr.psi(pts), compute_uv=False)
        p = db(s**2 / s[0]**2)
        n10, n20, n40 = int((p > -10).sum()), int((p > -20).sum()), int((p > -40).sum())
        part = float((s**2).sum()**2 / (s**4).sum())
        rows.append(dict(array=arr.name, N=arr.N, centre=list(map(float, centre)), sv_db=[round(float(v), 1) for v in p],
                         within_10dB=n10, within_20dB=n20, within_40dB=n40, participation=round(part, 2)))
        say(f"  {arr.name:46s} N={arr.N:2d}  modes within 10/20/40 dB of the strongest: {n10}/{n20}/{n40}"
            f"   participation number {part:.1f}")
    RES["dof"] = rows


# ---- 2. focusing ----------------------------------------------------------------------------------------
def sec_focus(arrs):
    say("\n== 2. focal spot (unit-norm excitation maximising |E|^2 at the point; 'beam' = same aimed at 100 m)")
    rows = []
    fig, ax = plt.subplots(2, len(arrs), figsize=(4.2 * len(arrs), 7))
    for c, arr in enumerate(arrs):
        Rm = arr.R()
        for zf in (0.2, 0.3, 0.5):
            rf = np.array([0, 0, zf])
            w, _ = focus_w(arr, rf)
            wb, _ = focus_w(arr, rf * 100 / zf)
            s = np.linspace(0.04, 1.2, 1161)
            ray = np.c_[np.zeros_like(s), np.zeros_like(s), s]
            Iz, Izb = intensity(arr, ray, w), intensity(arr, ray, wb)
            i0 = int(np.argmin(abs(s - zf)))
            depth, zpk = fwhm(s, Iz, i0, int(0.5 * zf / (s[1] - s[0])))
            lat = np.linspace(-0.6, 0.6, 1201)
            Ix = intensity(arr, np.c_[lat, 0 * lat, 0 * lat + zf], w)
            Iy = intensity(arr, np.c_[0 * lat, lat, 0 * lat + zf], w)
            wx, _ = fwhm(lat, Ix, 600, 3)
            wy, _ = fwhm(lat, Iy, 600, 3)
            If = Iz[i0]
            Ef = arr.G(rf[None])[0] @ w
            zshare = float(abs(Ef[2])**2 / (abs(Ef)**2).sum())
            gain = float(db(If / Izb[i0]))
            prad = float(np.real(w.conj() @ Rm @ w))
            e1mw = float(np.sqrt(If * 1e-3 / prad))               # V/m (peak) at 1 mW radiated
            s1mw = e1mw**2 / (2 * ETA)                            # W/m^2, plane-wave equivalent
            D = np.ptp(arr.epos[:, 0])
            rows.append(dict(array=arr.name, focus_m=zf, aperture_mm=round(D * 1e3), fraunhofer_m=round(2 * D**2 / arr.lam, 3),
                             fwhm_x_mm=None if wx is None else round(wx * 1e3, 1), fwhm_y_mm=None if wy is None else round(wy * 1e3, 1),
                             depth_mm=None if depth is None else round(depth * 1e3, 1),
                             axial_peak_m=None if zpk is None else round(zpk, 3),
                             gain_over_beam_db=round(gain, 2), peak_on_axis_db_above_focus=round(float(db(Iz.max() / If)), 1),
                             e_v_per_m_at_1mW=round(e1mw, 3), s_uW_per_cm2_at_1mW=round(s1mw * 100, 4),
                             lambda_mm=round(arr.lam * 1e3, 1), ez_share_of_power_at_focus=round(zshare, 3)))
            r = rows[-1]
            say(f"  {arr.name:26s} focus {zf:.1f} m: FWHM x/y {r['fwhm_x_mm']}/{r['fwhm_y_mm']} mm, depth {r['depth_mm']} mm "
                f"(axial peak at {r['axial_peak_m']} m), gain over beam {gain:+.2f} dB, "
                f"{e1mw:.2f} V/m = {s1mw * 100:.3f} uW/cm^2 per mW radiated, E_z share {zshare:.2f}  [2D^2/lambda = {r['fraunhofer_m']} m]")
            if zf == 0.3:
                gx, gz = np.linspace(-0.4, 0.4, 161), np.linspace(0.04, 0.8, 153)
                XX, ZZ = np.meshgrid(gx, gz)
                I = intensity(arr, np.c_[XX.ravel(), 0 * XX.ravel(), ZZ.ravel()], w).reshape(XX.shape)
                im = ax[0, c].pcolormesh(gx, gz, db(I / If), vmin=-25, vmax=10, cmap="magma", shading="auto")
                ax[0, c].plot(0, zf, "c+", ms=10)
                ax[0, c].set(title=f"{arr.name}\nfocus at 0.3 m, |E|^2 dB re focus", xlabel="x (m)", ylabel="z (m)")
                plt.colorbar(im, ax=ax[0, c])
                ax[1, c].plot(s, db(Iz / If), label="focused at 0.3 m")
                ax[1, c].plot(s, db(Izb / If), "--", label="ordinary beam")
                ax[1, c].axvline(zf, color="k", lw=0.5)
                ax[1, c].set(xlabel="distance along the axis (m)", ylabel="|E|^2 dB re focus", ylim=(-20, 25))
                ax[1, c].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "focus.png", dpi=110)
    plt.close(fig)
    RES["focus"] = rows


# ---- 3. several spots -----------------------------------------------------------------------------------
def multi_w(arr, pts, iters=50):
    """Equal field at several points with the least excitation norm; the spot phases are free and iterated."""
    g = []
    for p in pts:
        Gp = arr.G(np.asarray(p)[None])[0]
        U, s, Vh = np.linalg.svd(Gp)
        g.append(U[:, 0].conj() @ Gp)             # field along each point's own best polarisation
    g = np.array(g)
    d = np.ones(len(pts), complex)
    for _ in range(iters):
        w = np.linalg.pinv(g) @ d
        d = np.exp(1j * np.angle(g @ w))
    return unit(np.linalg.pinv(g) @ d)


def sec_spots(cases):
    say("\n== 3. simultaneous spots (dB relative to a single focus at the same point with the same excitation norm)")
    rows = []
    fig, ax = plt.subplots(1, len(cases), figsize=(5 * len(cases), 4.3))
    for c, (arr, pts, zc) in enumerate(cases):
        for n in (2, 3):
            P = np.array(pts[:n])
            w = multi_w(arr, P)
            rel = [float(db(intensity(arr, p[None], w)[0] / focus_w(arr, p)[1])) for p in P]
            g = np.linspace(-0.45, 0.45, 181)
            XX, YY = np.meshgrid(g, g)
            I = intensity(arr, np.c_[XX.ravel(), YY.ravel(), 0 * XX.ravel() + zc], w).reshape(XX.shape)
            Isp = np.mean([intensity(arr, p[None], w)[0] for p in P])
            far = np.min(np.hypot(XX[..., None] - P[:, 0], YY[..., None] - P[:, 1]), axis=-1) > 0.75 * arr.lam * zc / max(np.ptp(arr.epos[:, 0]), arr.lam)
            side = float(db(I[far].max() / Isp)) if far.any() else None
            rows.append(dict(array=arr.name, spots=P.tolist(), rel_single_focus_db=[round(v, 2) for v in rel],
                             ideal_split_db=round(float(-db(n)), 2), worst_elsewhere_db_re_spots=None if side is None else round(side, 1)))
            say(f"  {arr.name:26s} {n} spots: {['%.1f' % v for v in rel]} dB (ideal power split {-db(n):.1f} dB); "
                f"brightest point elsewhere in the plane {rows[-1]['worst_elsewhere_db_re_spots']} dB re spots")
            if n == 3:
                im = ax[c].pcolormesh(g, g, db(I / Isp), vmin=-25, vmax=3, cmap="magma", shading="auto")
                ax[c].plot(P[:, 0], P[:, 1], "c+", ms=10)
                ax[c].set(title=f"{arr.name}: 3 spots, plane z = {zc} m", xlabel="x (m)", ylabel="y (m)", aspect="equal")
                plt.colorbar(im, ax=ax[c])
    fig.tight_layout()
    fig.savefig(OUT / "spots.png", dpi=110)
    plt.close(fig)
    RES["spots"] = rows


# ---- rev B's constraint: 4 channels, per-element 4-bit true-time delay ------------------------------------
def ttd_phases(f):
    tau = np.array([sum(c) for r in range(5) for c in itertools.combinations(TTD_PS, r)]) * 1e-12
    return 2 * np.pi * f * tau                      # e^{-i w t}: a delay tau multiplies the phasor by e^{+i w tau}


def channel_of(arr):
    ix, iy = np.arange(arr.N) % arr.nx, np.arange(arr.N) // arr.nx
    return (ix // 2) + 2 * (iy // 2)


def revb_project(arr, w):
    """Nearest rev B excitation: one complex weight per 2x2 quadrant, one of 16 delays per element."""
    ps, ch = ttd_phases(arr.f), channel_of(arr)
    out = np.zeros(arr.N, complex)
    for c in range(4):
        idx = np.nonzero(ch == c)[0]
        best = (-np.inf, None)
        for th in np.linspace(0, 2 * np.pi, 180, endpoint=False):
            err = wrap(np.angle(w[idx])[:, None] - th - ps[None])
            s = np.argmin(abs(err), axis=1)
            val = float((abs(w[idx]) * np.cos(err[np.arange(len(idx)), s])).sum())
            if val > best[0]:
                best = (val, np.exp(1j * (th + ps[s])))
        out[idx] = best[0] / len(idx) * best[1]
    return out


def revb_decode(arr, X):
    """Design vector -> excitations: 8 numbers = 4 complex channel weights, 16 numbers in [0, 16) = delay states."""
    X = np.atleast_2d(X)
    ps, ch = ttd_phases(arr.f), channel_of(arr)
    cw = X[:, :4] + 1j * X[:, 4:8]
    st = np.clip(np.floor(X[:, 8:24]).astype(int), 0, 15)
    W = cw[:, ch] * np.exp(1j * ps[st])
    return W / np.maximum(np.linalg.norm(W, axis=1, keepdims=True), 1e-12)


# ---- 4. quiet zone --------------------------------------------------------------------------------------
def ball(centre, a, n=2):
    g = np.arange(-n, n + 1) / n * a
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    p = np.c_[X.ravel(), Y.ravel(), Z.ravel()]
    return p[np.linalg.norm(p, axis=1) <= a * 1.0001] + centre


def sec_quiet(cases, quick):
    say("\n== 4. quiet zone at B while a bright point is kept at A")
    rows = []
    rng = np.random.default_rng(5)
    for arr, rA, rB in cases:
        w0, IA0 = focus_w(arr, rA)
        GB = arr.G(rB[None])[0]
        IB0 = float(intensity(arr, rB[None], w0)[0])
        # (a) exact null of all three components at the point B: remove the part of w0 that the 3 rows of GB see
        w = unit(w0 - np.linalg.pinv(GB) @ (GB @ w0))
        IAn, IBn = float(intensity(arr, rA[None], w)[0]), float(intensity(arr, rB[None], w)[0])
        e = errors(rng, arr.N, 300)
        IBerr = float(np.mean((np.abs(np.einsum("cn,dn->dc", GB, w * e))**2).sum(1)))
        radii = np.array([0.02, 0.05, 0.1, 0.2, 0.35]) * arr.lam
        sizes = {}
        for rad in radii:
            sh = ball(rB, rad, 3)
            sh = sh[np.linalg.norm(sh - rB, axis=1) > 0.8 * rad]
            sizes[f"{rad / arr.lam:.2f} lambda"] = round(float(db(intensity(arr, sh, w).max() / IAn)), 1)
        row = dict(array=arr.name, A=rA.tolist(), B=rB.tolist(), B_before_db_re_A=round(float(db(IB0 / IA0)), 1),
                   point_null=dict(cost_at_A_db=round(float(db(IAn / IA0)), 2), exact_db_re_A=round(float(db(IBn / IAn)), 0),
                                   with_errors_db_re_A=round(float(db(IBerr / IAn)), 1), worst_on_shell_db_re_A=sizes))
        say(f"  {arr.name}: A={rA.tolist()} B={rB.tolist()}  B before: {row['B_before_db_re_A']} dB re A")
        say(f"    point null: costs {row['point_null']['cost_at_A_db']} dB at A; exact {row['point_null']['exact_db_re_A']} dB; "
            f"with 3 deg / 0.26 dB errors {row['point_null']['with_errors_db_re_A']} dB; worst on a shell around B: {sizes}")
        # (b) zone: minimise the mean |E|^2 in a ball of radius a for unit field at A (diagonally loaded)
        gA = np.linalg.svd(arr.G(rA[None])[0])[0][:, 0].conj() @ arr.G(rA[None])[0]
        zone = {}
        for a_l in (0.125, 0.25, 0.5):
            pts = ball(rB, a_l * arr.lam, 3)
            Gz = arr.G(pts).reshape(-1, arr.N)
            Z = Gz.conj().T @ Gz / len(pts)
            for load in (1e-4, 1e-2):
                wz = unit(np.linalg.solve(Z + load * np.trace(Z).real / arr.N * np.eye(arr.N), gA.conj()))
                IAz = float(intensity(arr, rA[None], wz)[0])
                Iz = intensity(arr, pts, wz)
                Ie = np.mean([(np.abs(Gz @ (wz * ee))**2).sum() / len(pts) for ee in e[:100]])
                zone[f"a={a_l} lambda, loading {load:g}"] = dict(
                    cost_at_A_db=round(float(db(IAz / IA0)), 2), mean_db_re_A=round(float(db(Iz.mean() / IAz)), 1),
                    worst_db_re_A=round(float(db(Iz.max() / IAz)), 1), mean_with_errors_db_re_A=round(float(db(Ie / IAz)), 1))
                z = zone[f"a={a_l} lambda, loading {load:g}"]
                say(f"    zone radius {a_l} lambda ({a_l * arr.lam * 1e3:.0f} mm), loading {load:g}: mean {z['mean_db_re_A']} dB, worst "
                    f"{z['worst_db_re_A']} dB, with errors {z['mean_with_errors_db_re_A']} dB re A; costs {z['cost_at_A_db']} dB at A")
        row["zone"] = zone
        # (c) the same zone (a = lambda/4) under rev B's constraint, searched with the repo's optimiser
        if hasattr(arr, "nx") and arr.N == 16:
            pts = ball(rB, 0.25 * arr.lam, 3)
            Gz = arr.G(pts).reshape(-1, arr.N)
            GA = arr.G(rA[None])[0]

            def score(X, gen, decode=lambda X: revb_decode(arr, X)):
                W = decode(X)
                ee = errors(np.random.default_rng(777 + gen), arr.N, 6)
                tot = 0
                for k in range(len(ee)):
                    We = W * ee[k]
                    IA = (np.abs(We @ GA.T)**2).sum(1)
                    IZ = (np.abs(We @ Gz.T)**2).reshape(len(W), -1, 3).sum(2).mean(1)
                    # quiet zone relative to A, and do not give away more than 6 dB at A
                    tot = tot + db(IZ / IA) + 2.0 * np.maximum(0, -6 - db(IA / IA0))
                return tot / len(ee)

            lo, hi = np.r_[-np.ones(8), np.zeros(16)], np.r_[np.ones(8), 16 * np.ones(16) - 1e-9]
            gens = 60 if quick else 400
            res = evolve(score, lo, hi, pop=64, generations=gens, seed=1)
            res2 = climb(score, res.x, lo, hi, steps=gens, lam=16, sigma=0.03, seed=2)
            rs = random_search(score, lo, hi, batches=gens, batch=64, seed=3)
            cand, sc = reselect(score, np.vstack([res.finalists[:8], res2.x, rs.finalists[:2]]), draws=24)
            wr = revb_decode(arr, cand[0])[0]
            IAr = float(intensity(arr, rA[None], wr)[0])
            # unconstrained reference on the same score: encode by hand is not possible, evaluate directly
            wz = unit(np.linalg.solve(Gz.conj().T @ Gz / len(pts) + 1e-2 * np.trace(Gz.conj().T @ Gz).real / len(pts) / 16 * np.eye(16), gA.conj()))
            ref = float(np.mean(score(wz[None], 1_000_000 + 0, decode=lambda W: np.atleast_2d(W))))
            # what the delays alone give (channels equal): project the unconstrained solution
            wp = unit(revb_project(arr, wz))
            proj = float(np.mean([score(wp[None], 1_000_000 + k, decode=lambda W: np.atleast_2d(W)) for k in range(24)]))
            row["revB_zone_quarter_lambda"] = dict(
                evolve_then_climb_db=round(float(sc[0]), 1), cost_at_A_db=round(float(db(IAr / IA0)), 2),
                random_search_db=round(float(reselect(score, rs.finalists[:1], draws=24)[1][0]), 1),
                unconstrained_db=round(ref, 1), projected_unconstrained_db=round(proj, 1), evaluations=res.evaluations + res2.evaluations)
            say(f"    rev B (4 channel weights x 16 delay states), zone lambda/4, with errors: {sc[0]:.1f} dB re A after evolve+climb "
                f"(random search {row['revB_zone_quarter_lambda']['random_search_db']} dB; unconstrained {ref:.1f} dB; "
                f"unconstrained solution snapped to rev B {proj:.1f} dB); costs {row['revB_zone_quarter_lambda']['cost_at_A_db']} dB at A")
        rows.append(row)
    # figure: the first case
    arr, rA, rB = cases[0]
    w0, IA0 = focus_w(arr, rA)
    pts = ball(rB, 0.25 * arr.lam, 3)
    Gz = arr.G(pts).reshape(-1, arr.N)
    gA = np.linalg.svd(arr.G(rA[None])[0])[0][:, 0].conj() @ arr.G(rA[None])[0]
    Z = Gz.conj().T @ Gz / len(pts)
    wz = unit(np.linalg.solve(Z + 1e-2 * np.trace(Z).real / arr.N * np.eye(arr.N), gA.conj()))
    g = np.linspace(-0.25, 0.25, 201)
    XX, YY = np.meshgrid(g, g)
    P = np.c_[XX.ravel(), YY.ravel(), 0 * XX.ravel() + rA[2]]
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.3))
    for a_, w_, t in ((ax[0], w0, "focus at A only"), (ax[1], wz, "focus at A, quiet ball (lambda/4) at B")):
        I = intensity(arr, P, w_).reshape(XX.shape)
        im = a_.pcolormesh(g, g, db(I / IA0), vmin=-45, vmax=0, cmap="magma", shading="auto")
        a_.plot(rA[0], rA[1], "c+", ms=10)
        a_.add_patch(plt.Circle(rB[:2], 0.25 * arr.lam, fill=False, color="w"))
        a_.set(title=f"{arr.name}: {t}", xlabel="x (m)", ylabel="y (m)", aspect="equal")
        plt.colorbar(im, ax=a_)
    fig.tight_layout()
    fig.savefig(OUT / "quiet_zone.png", dpi=110)
    plt.close(fig)
    RES["quiet"] = rows


# ---- 5. superoscillation --------------------------------------------------------------------------------
def sec_superosc(cases):
    say("\n== 5. a spot narrower than the diffraction spot (zeros forced at +-d around the focus, least-norm excitation)")
    rows = []
    for arr, rf in cases:
        g0 = arr.psi(rf[None])[0]                 # scalar component (E_y panel, E_z whips), so the whole row is one field
        w0 = unit(g0.conj())
        I0 = float(abs(g0 @ w0)**2)
        lat = np.linspace(-0.6, 0.6, 2401)
        line = np.c_[lat, 0 * lat, 0 * lat + rf[2]]
        Gl = arr.psi(line)
        I = abs(Gl @ w0)**2
        f0, _ = fwhm(lat, I, 1200, 3)
        for frac in (1.0, 0.6, 0.4, 0.25):
            if frac == 1.0:
                w = w0
            else:
                d = frac * f0
                pts = np.array([rf, rf + [d, 0, 0], rf - [d, 0, 0]])
                g = arr.psi(pts)
                w = unit(np.linalg.pinv(g) @ np.array([1, 0, 0]))
            I = abs(Gl @ w)**2
            fw, _ = fwhm(lat, I, 1200, 3)
            e = errors(np.random.default_rng(3), arr.N, 200)
            Ie = np.mean(abs((w * e) @ Gl.T)**2, axis=0)
            fwe, _ = fwhm(lat, Ie, 1200, 3)
            rows.append(dict(array=arr.name, focus=rf.tolist(), zeros_at_fraction_of_fwhm=frac, fwhm_mm=round(fw * 1e3, 1),
                             fwhm_lambda=round(fw / arr.lam, 3), spot_db_re_plain_focus=round(float(db(I[1200] / I0)), 1),
                             brightest_sidelobe_db_re_spot=round(float(db(I.max() / I[1200])), 1),
                             fwhm_with_errors_mm=None if fwe is None else round(fwe * 1e3, 1)))
            r = rows[-1]
            say(f"  {arr.name:26s} zeros at {frac:4.2f} FWHM: spot {r['fwhm_mm']} mm ({r['fwhm_lambda']} lambda), spot power "
                f"{r['spot_db_re_plain_focus']} dB re plain focus, brightest sidelobe {r['brightest_sidelobe_db_re_spot']:+} dB re spot, "
                f"FWHM with errors {r['fwhm_with_errors_mm']} mm")
    RES["superosc"] = rows


# ---- 6. vortex line -------------------------------------------------------------------------------------
def vortex_w(arr, l, theta_deg, phi_deg=0.0):
    th, ph = np.deg2rad(theta_deg), np.deg2rad(phi_deg)
    n = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
    az = np.arctan2(arr.epos[:, 1], arr.epos[:, 0])
    # e^{-i w t}: the far field of element x goes as e^{ikr} e^{-ik n.x}, so steering towards n needs e^{+ik n.x}
    return unit(np.exp(1j * l * az) * np.exp(1j * arr.k * (arr.epos @ n)))


def find_vortices(arr, w, z, half=0.3, n=241):
    g = np.linspace(-half, half, n)
    XX, YY = np.meshgrid(g, g, indexing="ij")
    psi = (arr.psi(np.c_[XX.ravel(), YY.ravel(), 0 * XX.ravel() + z]) @ w).reshape(XX.shape)
    q = winding2d(psi)                               # positive about +z
    i, j = np.nonzero(q)
    c = 0.5 * (g[1:] + g[:-1])
    return g, psi, np.c_[c[i], c[j]], q[i, j]


def sec_vortex(arrs):
    say("\n== 6. a steerable vortex line (elements phased exp(i l azimuth) plus a steering phase)")
    rows = []
    zs = (0.1, 0.2, 0.3, 0.4)
    for arr in arrs:
        for label, proj in (("ideal", False), ("rev B delays", True)):
            if proj and not (hasattr(arr, "nx") and arr.N == 16):
                continue
            for l, th in ((1, 0), (1, 15), (1, 30), (2, 0)):
                w = vortex_w(arr, l, th)
                if proj:
                    w = unit(revb_project(arr, w))
                pos = []
                for z in zs:
                    g, psi, xy, q = find_vortices(arr, w, z)
                    exp = np.array([z * np.tan(np.deg2rad(th)), 0.0])
                    if len(xy) == 0:
                        pos.append(None)
                        continue
                    d = np.linalg.norm(xy - exp, axis=1)
                    near = np.argsort(d)[:abs(l)]
                    I0 = focus_w(arr, np.array([exp[0], 0, z]))[1]
                    ring = max(abs(psi.ravel())**2) / I0
                    pos.append(dict(z=z, xy_mm=[[round(float(v) * 1e3, 1) for v in xy[k]] for k in near],
                                    charge=[int(q[k]) for k in near], peak_db_re_focus=round(float(db(ring)), 1)))
                tilt = [None if p is None else round(float(np.degrees(np.arctan2(np.mean([v[0] for v in p["xy_mm"]]) / 1e3, p["z"]))), 1) for p in pos]
                rows.append(dict(array=arr.name, excitation=label, l=l, steer_deg=th, slices=pos, line_tilt_deg=tilt))
                say(f"  {arr.name:22s} {label:12s} l={l} steer {th:2d} deg: vortex at "
                    + "; ".join("none" if p is None else f"z={p['z']}: {p['xy_mm']} mm q={p['charge']}" for p in pos)
                    + f"  -> tilt {tilt} deg")
    # figure: phase slices for the 5.8 GHz panel
    arr = arrs[0]
    fig, ax = plt.subplots(2, 4, figsize=(17, 8))
    for r, th in enumerate((0, 20)):
        w = vortex_w(arr, 1, th)
        for c, z in enumerate(zs):
            g, psi, xy, q = find_vortices(arr, w, z, half=0.25, n=201)
            a = ax[r, c]
            im = a.pcolormesh(g, g, np.angle(psi).T, cmap="twilight", vmin=-np.pi, vmax=np.pi, shading="auto")
            a.contour(g, g, db(abs(psi.T)**2 / (abs(psi)**2).max()), levels=[-20, -10, -3], colors="w", linewidths=0.6)
            for (x, y), qq in zip(xy, q):
                a.plot(x, y, "o", mfc="none", mec="lime" if qq > 0 else "red", ms=9)
            a.set(title=f"steer {th} deg, z = {z} m: phase of E_y", xlabel="x (m)", ylabel="y (m)", aspect="equal")
    fig.suptitle(f"{arr.name}: phase of E_y in four planes; circles = vortices (green +1, red -1); white = -3/-10/-20 dB")
    fig.tight_layout()
    fig.savefig(OUT / "vortex_slices.png", dpi=100)
    plt.close(fig)
    RES["vortex"] = rows


# ---- 7. polarisation singularities (dual-polarised panel, hypothetical) -----------------------------------
def sec_pol(f):
    say("\n== 7. polarisation structure if the panel had both polarisations (32 excitations; rev B has one)")
    arr = panel(f, pol="xy", name=f"dual-pol panel {f / 1e9:g} GHz")
    one = panel(f)
    a = vortex_w(one, 0, 0)            # l = 0 in left circular
    b = vortex_w(one, 1, 0) * 1.5      # l = 1 in right circular
    w = unit(np.r_[(a + b) / np.sqrt(2), 1j * (a - b) / np.sqrt(2)])     # E = A e_L + B e_R, e_L,R = (x +- i y)/sqrt 2
    rows = []
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    for c, z in enumerate((0.1, 0.2, 0.3)):
        g = np.linspace(-0.2, 0.2, 161)
        XX, YY = np.meshgrid(g, g, indexing="ij")
        E = (arr.G(np.c_[XX.ravel(), YY.ravel(), 0 * XX.ravel() + z]) @ w).reshape(*XX.shape, 3)
        Ex, Ey, Ez = E[..., 0], E[..., 1], E[..., 2]
        S0 = abs(Ex)**2 + abs(Ey)**2
        S1, S2, S3 = abs(Ex)**2 - abs(Ey)**2, 2 * np.real(Ex.conj() * Ey), 2 * np.imag(Ex.conj() * Ey)
        sig = S1 + 1j * S2
        q = winding2d(sig)
        i, j = np.nonzero(q)
        cc = 0.5 * (g[1:] + g[:-1])
        keep = np.hypot(cc[i], cc[j]) < 0.12
        cpts = [dict(xy_mm=[round(float(cc[a_]) * 1e3, 1), round(float(cc[b_]) * 1e3, 1)], index=float(q[a_, b_]) / 2,
                     s3_over_s0=round(float(S3[a_, b_] / S0[a_, b_]), 2)) for a_, b_ in zip(i[keep], j[keep])]
        ezr = float(np.sqrt((abs(Ez)**2)[np.hypot(XX, YY) < 0.05].mean() / S0[np.hypot(XX, YY) < 0.05].mean()))
        rows.append(dict(z=z, c_points=cpts, ez_over_et_rms_within_50mm=round(ezr, 3)))
        say(f"  z = {z} m: C points (circular polarisation) within 120 mm of the axis: {cpts};  |Ez|/|Et| rms within 50 mm: {ezr:.2f}")
        a_ = ax[c]
        im = a_.pcolormesh(g, g, (S3 / S0).T, cmap="coolwarm", vmin=-1, vmax=1, shading="auto")
        psi_ = 0.5 * np.angle(sig)
        st = 8
        a_.quiver(XX[::st, ::st], YY[::st, ::st], np.cos(psi_[::st, ::st]), np.sin(psi_[::st, ::st]), headwidth=0, headlength=0,
                  headaxislength=0, pivot="mid", scale=40, width=0.003)
        for p in cpts:
            a_.plot(p["xy_mm"][0] / 1e3, p["xy_mm"][1] / 1e3, "o", mfc="none", mec="k", ms=10)
        a_.set(title=f"z = {z} m: S3/S0 (colour), major axis (ticks), C points (circles)", xlabel="x (m)", ylabel="y (m)", aspect="equal")
        plt.colorbar(im, ax=a_)
    fig.tight_layout()
    fig.savefig(OUT / "polarisation.png", dpi=100)
    plt.close(fig)
    RES["pol"] = rows


# ---- 8. vortex lines on a target curve ------------------------------------------------------------------
def targets(lam, centre):
    """Target curves, sized in wavelengths.  Each: list of components (dense points, closed?)."""
    t = np.linspace(0, 2 * np.pi, 400, endpoint=False)
    th = np.deg2rad(30)
    out = {}
    s = np.linspace(-1.5, 1.5, 200)[:, None] * lam
    out["line"] = [((s * np.array([np.sin(th), 0, np.cos(th)])) + centre, False)]
    e1, e2 = np.array([np.cos(th), 0, -np.sin(th)]), np.array([0, 1.0, 0])
    out["ring"] = [(lam * (np.cos(t)[:, None] * e1 + np.sin(t)[:, None] * e2) + centre, True)]
    r = 0.9 * lam
    out["link"] = [(np.c_[r * np.cos(t) - r / 2, r * np.sin(t), 0 * t] + centre, True),
                   (np.c_[r * np.cos(t) + r / 2, 0 * t, r * np.sin(t)] + centre, True)]
    sc = 0.45 * lam
    out["trefoil"] = [(sc * np.c_[np.sin(t) + 2 * np.sin(2 * t), np.cos(t) - 2 * np.cos(2 * t), -np.sin(3 * t)] + centre, True)]
    return out


def sample_curve(comp, m):
    pts, closed = comp
    idx = np.linspace(0, len(pts), m, endpoint=False).astype(int) if closed else np.linspace(0, len(pts) - 1, m).astype(int)
    tan = np.gradient(pts, axis=0) if not closed else (np.roll(pts, -1, 0) - np.roll(pts, 1, 0))
    tan = tan / np.linalg.norm(tan, axis=1, keepdims=True)
    return pts[idx], tan[idx]


def circles(c, t, rho, K):
    """K points on a circle of radius rho around each point c, in the plane normal to the tangent t,
    ordered so that the angle increases right-handedly about t."""
    a = np.where(abs(t[:, 2:3]) < 0.9, [[0, 0, 1.0]], [[1.0, 0, 0]])
    e1 = np.cross(t, a)
    e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
    e2 = np.cross(t, e1)
    th = 2 * np.pi * np.arange(K) / K
    return c[:, None, :] + rho * (np.cos(th)[None, :, None] * e1[:, None, :] + np.sin(th)[None, :, None] * e2[:, None, :])


def enc(w):
    return np.r_[w.real, w.imag] / max(np.abs(w.real).max(), np.abs(w.imag).max())


def dec(X):
    X = np.atleast_2d(X)
    n = X.shape[1] // 2
    W = X[:, :n] + 1j * X[:, n:]
    return W / np.maximum(np.linalg.norm(W, axis=1, keepdims=True), 1e-12)


class CurveProblem:
    """Score = 1 - (share of the field on small circles around the curve that is the first azimuthal harmonic
    e^{+-i theta}), averaged over the curve and over element errors.  0 = a clean unit vortex all along the curve."""

    def __init__(self, arr, comps, m_per=32, K=8, rho=0.15, draws=4):
        self.arr, self.K, self.N = arr, K, arr.N
        self.c, self.t, self.sl = [], [], []
        a = 0
        for comp in comps:
            c, t = sample_curve(comp, m_per)
            self.c.append(c)
            self.t.append(t)
            self.sl.append(slice(a, a + len(c)))
            a += len(c)
        self.c, self.t = np.vstack(self.c), np.vstack(self.t)
        self.M = len(self.c)
        self.rho = rho * arr.lam
        self.Gc = arr.psi(circles(self.c, self.t, self.rho, K).reshape(-1, 3))        # (M K, N)
        self.h = np.exp(-1j * 2 * np.pi * np.arange(K) / K) / K
        self.draws = draws

    def quality(self, W, e=None):
        """(designs,) mean over the curve of the best-handed first-harmonic share, per component sign free."""
        W = W if e is None else W * e
        psi = (W @ self.Gc.T).reshape(len(W), self.M, self.K)
        pw = (abs(psi)**2).mean(2) + 1e-300
        qp, qm = abs(psi @ self.h)**2 / pw, abs(psi @ self.h.conj())**2 / pw
        tot = 0
        for sl in self.sl:
            tot = tot + np.maximum(qp[:, sl].sum(1), qm[:, sl].sum(1))
        return tot / self.M

    def score(self, X, gen):
        W = dec(X)
        e = errors(np.random.default_rng(10_000 + gen), self.N, self.draws)
        return 1 - np.mean([self.quality(W, e[k]) for k in range(self.draws)], axis=0)

    def worst_score(self, X, gen):
        """Second-stage score: half the mean, half the WORST point of the curve.  One sample with the wrong
        winding is where the dark line leaves the target, so the worst point decides the topology."""
        W = dec(X)
        e = errors(np.random.default_rng(20_000 + gen), self.N, self.draws)
        tot = 0
        for k in range(self.draws):
            psi = ((W * e[k]) @ self.Gc.T).reshape(len(W), self.M, self.K)
            pw = (abs(psi)**2).mean(2) + 1e-300
            qp, qm = abs(psi @ self.h)**2 / pw, abs(psi @ self.h.conj())**2 / pw
            q = np.zeros((len(W), self.M))
            for sl in self.sl:
                use = qp[:, sl].sum(1) >= qm[:, sl].sum(1)
                q[:, sl] = np.where(use[:, None], qp[:, sl], qm[:, sl])
            tot = tot + 0.5 * (1 - q.mean(1)) + 0.5 * (1 - q.min(1))
        return tot / self.draws

    def direct(self, load=1e-3):
        """Linear-algebra answer for the pooled ratio: a generalised eigenvector (the problem is quadratic in w)."""
        G3 = self.Gc.reshape(self.M, self.K, self.N)
        B = np.einsum("mkn,mkp->np", G3.conj(), G3) / self.K
        B = B + load * np.trace(B).real / self.N * np.eye(self.N)
        best = (np.inf, None)
        for signs in itertools.product((1, -1), repeat=len(self.sl)):
            A = np.zeros((self.N, self.N), complex)
            for s, sl in zip(signs, self.sl):
                g1 = np.einsum("k,mkn->mn", self.h if s > 0 else self.h.conj(), G3[sl])
                A += g1.conj().T @ g1
            _, V = eigh(A, B)
            w = unit(V[:, -1])
            sc = float(np.mean([self.score(enc(w)[None], 500_000 + k)[0] for k in range(8)]))
            if sc < best[0]:
                best = (sc, w)
        return best[1]


def vortex_faces(psi):
    """Integer vortex flux through every face of a 3-D grid: (through x-faces, y-faces, z-faces)."""
    ph = np.angle(psi)
    dx, dy, dz = wrap(ph[1:] - ph[:-1]), wrap(ph[:, 1:] - ph[:, :-1]), wrap(ph[:, :, 1:] - ph[:, :, :-1])
    wz = dx[:, :-1, :] + dy[1:, :, :] - dx[:, 1:, :] - dy[:-1, :, :]
    wx = dy[:, :, :-1] + dz[:, 1:, :] - dy[:, :, 1:] - dz[:, :-1, :]
    wy = dz[:-1, :, :] + dx[:, :, 1:] - dz[1:, :, :] - dx[:, :, :-1]
    return [np.rint(w / (2 * np.pi)).astype(int) for w in (wx, wy, wz)]


def trace_lines(W, axes):
    """Follow vortex lines from face to face through the cells.  Returns [(points (n, 3), closed?)]."""
    faces = {}
    for a in range(3):
        for idx in np.argwhere(W[a] != 0):
            faces[(a, *map(int, idx))] = int(W[a][tuple(idx)])
    ncell = [len(ax) - 1 for ax in axes]
    half = [0.5 * (ax[1] - ax[0]) for ax in axes]

    def centre(fc):
        a, i, j, k = fc
        p = [axes[0][i], axes[1][j], axes[2][k]]
        for b in range(3):
            if b != a:
                p[b] += half[b]
        return np.array(p)

    used, lines = set(), []
    for start in faces:
        if start in used:
            continue
        line, cur, closed = [start], start, False
        used.add(start)
        while True:
            a, c = cur[0], list(cur[1:])
            if faces[cur] < 0:
                c[a] -= 1
            if any(c[b] < 0 or c[b] >= ncell[b] for b in range(3)):
                break
            cands = []
            for b in range(3):
                up = list(c)
                up[b] += 1
                fu, fl = (b, *up), (b, *c)
                if faces.get(fu, 0) > 0 and (fu not in used or fu == start):
                    cands.append(fu)
                if faces.get(fl, 0) < 0 and (fl not in used or fl == start):
                    cands.append(fl)
            if not cands:
                break
            if start in cands and len(line) >= 3:
                closed = True
                break
            cands = [f for f in cands if f != start] or cands
            pc = centre(cur)
            cur = min(cands, key=lambda f: np.linalg.norm(centre(f) - pc))
            used.add(cur)
            line.append(cur)
        lines.append((np.array([centre(f) for f in line]), closed))
    return lines


def linking_number(a, b):
    """Gauss linking integral of two closed polylines."""
    da, db_ = np.roll(a, -1, 0) - a, np.roll(b, -1, 0) - b
    ma, mb = a + da / 2, b + db_ / 2
    r = ma[:, None, :] - mb[None, :, :]
    cr = np.cross(da[:, None, :], db_[None, :, :])
    return float((np.einsum("ijk,ijk->ij", cr, r) / np.linalg.norm(r, axis=2)**3).sum() / (4 * np.pi))


def min_self_distance(comps):
    P = np.vstack([c[0] for c in comps])
    cid = np.concatenate([np.full(len(c[0]), i) for i, c in enumerate(comps)])
    par = np.concatenate([np.arange(len(c[0])) / len(c[0]) for c in comps])
    d = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)
    dp = abs(par[:, None] - par[None, :])
    dp = np.minimum(dp, 1 - dp)
    far = (cid[:, None] != cid[None, :]) | (dp > 0.12)
    return float(d[far].min()) if far.any() else np.inf


def verify(arr, w, comps, prob, kind, step=1 / 20, margin=0.75):
    """Independent check: trace the vortex lines of psi on a 3-D grid and compare them with the target."""
    lam = arr.lam
    P = np.vstack([c[0] for c in comps])
    lo, hi = P.min(0) - margin * lam, P.max(0) + margin * lam
    # shift the grid by an odd fraction of a cell: a target that lies exactly in a grid plane (a ring at z = 0) makes
    # its vortex line run along cell faces, where rounding decides which cells it visits and the tracer fragments it
    lo = lo - step * lam * np.array([0.37, 0.23, 0.41])
    axes = [np.arange(lo[a], hi[a] + 1e-9, step * lam) for a in range(3)]
    X, Y, Z = np.meshgrid(*axes, indexing="ij")
    psi = (arr.psi(np.c_[X.ravel(), Y.ravel(), Z.ravel()]) @ w).reshape(X.shape)
    W = vortex_faces(psi)
    lines = trace_lines(W, axes)
    vp = np.vstack([l[0] for l in lines]) if lines else np.zeros((0, 3))
    out = dict(vortex_faces_in_box=int(len(vp)), lines=len(lines), closed_loops=int(sum(l[1] for l in lines)))
    # distance from each target sample to the nearest vortex line point
    if len(vp):
        d = cKDTree(vp).query(P)[0]
    else:
        d = np.full(len(P), np.inf)
    out.update(dist_median_lambda=round(float(np.median(d)) / lam, 3), dist_p90_lambda=round(float(np.quantile(d, 0.9)) / lam, 3),
               dist_max_lambda=round(float(d.max()) / lam, 3), dist_max_mm=round(float(d.max()) * 1e3, 1),
               target_within_0p1_lambda=round(float((d < 0.1 * lam).mean()), 3))
    # share of all vortex-line length in the box that is NOT the target (collateral structure)
    tt = cKDTree(P)
    if len(vp):
        dv = tt.query(vp)[0]
        out["vortex_length_off_target"] = round(float((dv > 0.15 * lam).mean()), 3)
    # charge on circles around the target (24 points, radius as in the score), exact and with errors
    cir = circles(prob.c, prob.t, prob.rho, 24).reshape(-1, 3)
    Gc = arr.psi(cir)

    def charges(wv):
        ph = np.angle((Gc @ wv).reshape(prob.M, 24))
        return np.rint(wrap(np.roll(ph, -1, 1) - ph).sum(1) / (2 * np.pi)).astype(int)

    q = charges(w)
    ok = 0
    signs = []
    for sl in prob.sl:
        s = 1 if (q[sl] == 1).sum() >= (q[sl] == -1).sum() else -1
        signs.append(s)
        ok += (q[sl] == s).sum()
    out["charge_correct"] = round(float(ok / prob.M), 3)
    rng = np.random.default_rng(99)
    okn = []
    for e in errors(rng, arr.N, 40):
        qe = charges(w * e)
        okn.append(sum((qe[sl] == s).sum() for sl, s in zip(prob.sl, signs)) / prob.M)
    out["charge_correct_with_errors"] = round(float(np.mean(okn)), 3)
    out["all_correct_with_errors_share_of_draws"] = round(float(np.mean(np.array(okn) == 1.0)), 2)
    # brightness of the field next to the line, relative to what the same excitation norm gives when simply focused there
    cen = P.mean(0)
    bright = (abs(prob.Gc @ w)**2).mean() / np.linalg.norm(arr.psi(cen[None]))**2
    out["field_on_circles_db_re_focus"] = round(float(db(bright)), 1)
    # topology: closed loops that follow one target component once around inside a tube
    tube = min(0.45 * min_self_distance(comps), 0.3 * lam) if kind != "line" else 0.1 * lam
    out["tube_radius_lambda"] = round(tube / lam, 3)
    cid = np.concatenate([np.full(len(c[0]), i) for i, c in enumerate(comps)])
    par = np.concatenate([np.arange(len(c[0])) / len(c[0]) for c in comps])
    realised = {}
    for pts, closed in lines:
        if not closed or len(pts) < 12:
            continue
        dd, ii = tt.query(pts)
        if dd.max() > tube or len(set(cid[ii])) != 1:
            continue
        dp = np.diff(np.r_[par[ii], par[ii][0]])
        dp = (dp + 0.5) % 1 - 0.5
        if abs(round(dp.sum())) == 1:
            realised[int(cid[ii][0])] = pts
    out["components_realised"] = len(realised)
    # geometry: a unit vortex of one sign on at least 95 % of the target, and a vortex line within 0.1 lambda of
    # at least 95 % of it.  topology: the traced lines are closed loops that each follow one component once around
    # inside the tube (and, for the link, have linking number +-1).  success needs both.
    out["geometry_ok"] = bool(out["charge_correct"] >= 0.95 and out["target_within_0p1_lambda"] >= 0.95)
    if kind == "line":
        out["geometry_ok"] = bool(out["charge_correct"] == 1.0 and d.max() < 0.1 * lam)
        out["topology_ok"] = out["geometry_ok"]
    elif kind == "link":
        lk = linking_number(realised[0], realised[1]) if len(realised) == 2 else None
        out["linking_number"] = None if lk is None else round(lk, 2)
        out["topology_ok"] = bool(lk is not None and abs(round(lk)) == 1)
    else:
        out["topology_ok"] = bool(len(realised) == len(comps))
    out["success"] = bool(out["geometry_ok"] and out["topology_ok"])
    return out, vp


def sec_curves(cases, quick, kinds):
    say("\n== 8. vortex lines on a target curve (score: 1 - first-harmonic share on circles of radius 0.15 lambda,")
    say("      with 3 deg / 0.26 dB element errors redrawn each generation; check: lines traced on a lambda/20 grid;")
    say("      success = geometry (unit charge on >= 95 % of the target and a vortex line within 0.1 lambda of >= 95 % of it)")
    say("      and topology (closed loops that follow each component once; linking number +-1 for the link))")
    gens = 40 if quick else 800
    rows = []
    fig = plt.figure(figsize=(4.2 * len(kinds), 3.9 * len(cases)))
    for ci, (arr, centre) in enumerate(cases):
        tg = targets(arr.lam, np.asarray(centre, float))
        for ki, kind in enumerate(kinds):
            t0 = time.time()
            comps = tg[kind]
            prob = CurveProblem(arr, comps, m_per=24 if kind == "line" else 48 if kind != "trefoil" else 96, draws=4 if quick else 8)
            lo, hi = -np.ones(2 * arr.N), np.ones(2 * arr.N)
            wd = prob.direct()
            ev = evolve(prob.score, lo, hi, pop=64, generations=gens, seed=11)
            cl_ev = climb(prob.score, ev.x, lo, hi, steps=gens, lam=16, sigma=0.02, seed=12)
            cl_d = climb(prob.score, enc(wd), lo, hi, steps=gens, lam=16, sigma=0.02, seed=13)
            rs = random_search(prob.score, lo, hi, batches=gens, batch=64, seed=14)
            names = ["direct (eigenvector)", "direct + climb", "evolve", "evolve + climb", "random search"]
            cand = np.vstack([enc(wd), cl_d.x, ev.finalists[0], cl_ev.x, rs.finalists[0]])
            # finalists re-scored on fresh draws; the pool also holds the rest of the last population
            pool = np.vstack([cand, ev.finalists[1:8]])
            ranked, sc = reselect(prob.score, pool, draws=16)
            per = {n: round(float(reselect(prob.score, c[None], draws=16)[1][0]), 4) for n, c in zip(names, cand)}
            w = dec(ranked[0])[0]
            winner = next((n for n, c in zip(names, cand) if np.allclose(c, ranked[0])), "evolve finalist")
            chk, vp = verify(arr, w, comps, prob, kind, step=1 / 12 if quick else 1 / 20)
            # second stage: climb on the worst-point score from the winner; keep it only if it wins on fresh draws
            pol = climb(prob.worst_score, ranked[0], lo, hi, steps=5 * gens, lam=16, sigma=0.01, seed=15)
            pr, psc = reselect(prob.worst_score, np.vstack([ranked[0], pol.x]), draws=16)
            chk2, vp2 = verify(arr, dec(pr[0])[0], comps, prob, kind, step=1 / 12 if quick else 1 / 20)
            keys = ("charge_correct", "charge_correct_with_errors", "target_within_0p1_lambda", "dist_median_lambda", "dist_max_lambda",
                    "field_on_circles_db_re_focus", "components_realised", "linking_number", "geometry_ok", "topology_ok", "success")
            polish = dict(worst_point_score_before_after=[round(float(v), 4) for v in
                                                          reselect(prob.worst_score, np.vstack([ranked[0], pol.x]), draws=16, first_gen=2_000_000)[1][::-1]]
                          if not np.allclose(pr[0], ranked[0]) else None,
                          kept=bool(not np.allclose(pr[0], ranked[0])), evaluations=pol.evaluations, **{k: chk2.get(k) for k in keys})
            row = dict(array=arr.name, N=arr.N, target=kind, scores=per, winner=winner, score=round(float(sc[0]), 4),
                       evaluations=dict(evolve=ev.evaluations, climb=cl_ev.evaluations, random=rs.evaluations), **chk, polish=polish)
            rows.append(row)
            say(f"  {arr.name:44s} {kind:8s} scores {per}")
            say(f"      winner: {winner}; check: charge correct {chk['charge_correct']:.2f} (with errors {chk['charge_correct_with_errors']:.2f}), "
                f"target within 0.1 lambda of a vortex line {chk['target_within_0p1_lambda']:.2f}, distance median/p90/max "
                f"{chk['dist_median_lambda']}/{chk['dist_p90_lambda']}/{chk['dist_max_lambda']} lambda, field near the line "
                f"{chk['field_on_circles_db_re_focus']} dB re focus, closed loops on target {chk['components_realised']}/{len(comps)}"
                + (f", linking number {chk.get('linking_number')}" if kind == "link" else "")
                + f"; geometry {'ok' if chk['geometry_ok'] else 'no'}, topology {'ok' if chk['topology_ok'] else 'no'}"
                + f"  => {'SUCCESS' if chk['success'] else 'FAIL'}   [{time.time() - t0:.0f} s]")
            say(f"      after the worst-point climb ({'kept' if polish['kept'] else 'not better on fresh draws, first stage kept'}): charge correct "
                f"{chk2['charge_correct']:.2f} (with errors {chk2['charge_correct_with_errors']:.2f}), within 0.1 lambda {chk2['target_within_0p1_lambda']:.2f}, "
                f"max miss {chk2['dist_max_lambda']} lambda, closed loops on target {chk2['components_realised']}/{len(comps)}"
                + (f", linking number {chk2.get('linking_number')}" if kind == "link" else "")
                + f"; geometry {'ok' if chk2['geometry_ok'] else 'no'}, topology {'ok' if chk2['topology_ok'] else 'no'}"
                + f"  => {'SUCCESS' if chk2['success'] else 'FAIL'}   [{time.time() - t0:.0f} s]")
            if (chk2["success"], chk2["topology_ok"], chk2["geometry_ok"]) > (chk["success"], chk["topology_ok"], chk["geometry_ok"]):
                chk, vp = chk2, vp2          # draw the better of the two stages
            ax = fig.add_subplot(len(cases), len(kinds), ci * len(kinds) + ki + 1, projection="3d")
            if len(vp):
                ax.scatter(*(vp.T * 1e3), s=1.5, c="0.55", alpha=0.5)
            for pts, _ in comps:
                ax.plot(*(pts.T * 1e3), "r-", lw=1.2)
            ax.set_title(f"{arr.name}\n{kind}: {'ok' if chk['success'] else 'topology only' if chk['topology_ok'] else 'fail'}, charge {chk['charge_correct']:.2f}, "
                         f"max miss {chk['dist_max_lambda']} lam", fontsize=7)
            ax.tick_params(labelsize=5)
    fig.suptitle("red = target curve; grey = vortex lines of the optimised field (all of them, in the checked box); axes in mm", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "curves.png", dpi=110)
    plt.close(fig)
    RES["curves"] = rows


# ---- 9. time reversal in a toy room ---------------------------------------------------------------------
def room_H(src, pts, k, L, gamma, order):
    """Scalar transfer from sources to points in a box room by image sources (Allen and Berkley 1979):
    image at (1 - 2q) x_s + 2 n L per axis, amplitude gamma^(number of wall hits)."""
    n = np.arange(-order, order + 1)
    img, amp = [], []
    for nx, ny, nz in itertools.product(n, n, n):
        for q in itertools.product((0, 1), repeat=3):
            nn, qq = np.array([nx, ny, nz]), np.array(q)
            hits = int((abs(nn - qq) + abs(nn)).sum())
            img.append((1 - 2 * qq, 2 * nn * L))
            amp.append(gamma**hits)
    amp = np.array(amp)
    sg = np.array([i[0] for i in img], float)
    sh = np.array([i[1] for i in img], float)
    H = np.zeros((len(pts), len(src)), complex)
    for a, s in enumerate(src):
        ip = sg * s + sh                                  # (I, 3)
        for b in range(0, len(pts), 400):
            r = np.linalg.norm(pts[b:b + 400, None, :] - ip[None], axis=2)
            H[b:b + 400, a] = (np.exp(1j * k * r) / (4 * np.pi * r)) @ amp
    return H


def sec_room(quick):
    say("\n== 9. toy multipath room (image sources, scalar field): time reversal against free-space focusing")
    f0 = 2.3e9
    L = np.array([5.0, 4.0, 2.6])
    gamma = 0.7
    order = 2 if quick else 3
    lam = C0 / f0
    wh = whips(f0)
    centre = np.array([1.0, 2.1, 1.0])
    target = np.array([3.7, 1.35, 1.25])
    g = np.linspace(-0.3, 0.3, 61 if quick else 101)
    XX, YY = np.meshgrid(g, g, indexing="ij")
    near = np.c_[XX.ravel(), YY.ravel(), 0 * XX.ravel()] + target
    rng = np.random.default_rng(4)
    bg = np.c_[rng.uniform(0.3, 4.7, 600), rng.uniform(0.3, 3.7, 600), rng.uniform(0.3, 2.3, 600)]
    bg = bg[(np.linalg.norm(bg - target, axis=1) > 0.5) & (np.linalg.norm(bg - centre, axis=1) > 0.8)]
    pts = np.vstack([target[None], near, bg])
    rows = []
    fig, ax = plt.subplots(2, 3, figsize=(15, 9))
    for ri, (label, sel) in enumerate((("16 whips", np.arange(16)), ("4 whips (X0-X3 transmit)", np.arange(0, 8, 2)))):
        src = wh.epos[sel] + centre
        freqs = {"single frequency": np.array([f0]), "20 MHz (11 tones)": f0 + np.linspace(-10e6, 10e6, 11)}
        if quick:
            freqs.pop("20 MHz (11 tones)")
        for bw, fs in freqs.items():
            acc = {m: np.zeros(len(pts), complex) for m in ("fs_in_fs", "fs_in_room", "tr_in_room")}
            inc = {m: np.zeros(len(pts)) for m in acc}
            for f in fs:
                k = 2 * np.pi * f / C0
                Hr = room_H(src, pts, k, L, gamma, order)
                Hf = room_H(src, pts, k, L, 0.0, 0)
                wf = unit(Hf[0].conj())
                wt = unit(Hr[0].conj())
                for m, (H, w) in dict(fs_in_fs=(Hf, wf), fs_in_room=(Hr, wf), tr_in_room=(Hr, wt)).items():
                    y = H @ w
                    acc[m] += y            # coherent over frequency: the field at the instant of the focus
                    inc[m] += abs(y)**2
            ref = abs(acc["fs_in_fs"][0])**2
            for m in acc:
                I = abs(acc[m])**2
                Imap = I[1:1 + len(near)].reshape(XX.shape)
                i0 = len(g) // 2
                wx, _ = fwhm(g, Imap[:, i0], i0, 2)
                wy, _ = fwhm(g, Imap[i0, :], i0, 2)
                Ibg = I[1 + len(near):]
                pk = np.unravel_index(np.argmax(Imap), Imap.shape)
                rows.append(dict(array=label, band=bw, method=m, at_target_db_re_free_space=round(float(db(I[0] / ref)), 1),
                                 fwhm_x_mm=None if wx is None else round(wx * 1e3), fwhm_y_mm=None if wy is None else round(wy * 1e3),
                                 fwhm_x_lambda=None if wx is None else round(wx / lam, 2),
                                 target_over_mean_background_db=round(float(db(I[0] / Ibg.mean())), 1),
                                 target_over_p99_background_db=round(float(db(I[0] / np.quantile(Ibg, 0.99))), 1),
                                 brightest_in_patch_offset_mm=round(float(np.hypot(g[pk[0]], g[pk[1]])) * 1e3),
                                 brightest_in_patch_db_re_target=round(float(db(Imap.max() / I[0])), 1)))
                r = rows[-1]
                say(f"  {label:26s} {bw:18s} {m:11s}: at target {r['at_target_db_re_free_space']:+.1f} dB re free space, spot "
                    f"{r['fwhm_x_mm']} x {r['fwhm_y_mm']} mm, target/background mean {r['target_over_mean_background_db']} dB "
                    f"(p99 {r['target_over_p99_background_db']} dB), brightest point of the 0.6 m patch is {r['brightest_in_patch_offset_mm']} mm "
                    f"from the target at {r['brightest_in_patch_db_re_target']:+.1f} dB")
                if bw == "single frequency":
                    a = ax[ri, list(acc).index(m)]
                    im = a.pcolormesh(g, g, db(Imap.T / ref), vmin=-20, vmax=8, cmap="magma", shading="auto")
                    a.plot(0, 0, "c+", ms=10)
                    a.set(title=f"{label}, {m}\n|field|^2 dB re free-space focus", xlabel="x - target (m)", ylabel="y - target (m)", aspect="equal")
                    plt.colorbar(im, ax=a)
    fig.tight_layout()
    fig.savefig(OUT / "time_reversal.png", dpi=100)
    plt.close(fig)
    RES["room"] = dict(room_m=L.tolist(), wall_reflection=gamma, image_order=order, f0=f0, array_centre=centre.tolist(),
                       target=target.tolist(), distance_m=round(float(np.linalg.norm(target - centre)), 2), rows=rows)


# ---- main -----------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    only = set(a.only.split(",")) - {""}
    OUT.mkdir(exist_ok=True)
    t0 = time.time()
    selftest()
    p24, p58, w09, w23 = panel(2.4e9), panel(5.8e9), whips(0.9e9), whips(2.3e9)
    p58x2 = panel(5.8e9, nx=8, ny=4, name="two tiles side by side (8x4) 5.8 GHz")
    s32, s64 = surround(5.8e9, 32), surround(5.8e9, 64)
    cp, cw, c0 = np.array([0, 0, 0.15]), np.array([0, 0, 0.35]), np.zeros(3)

    def run(name, fn, *args):
        if only and name not in only:
            return
        t = time.time()
        fn(*args)
        say(f"  [{name}: {time.time() - t:.0f} s]")
        (OUT / ("results_quick.json" if a.quick else "results.json")).write_text(json.dumps(RES, indent=1))

    run("dof", sec_dof, [(p24, cp), (p58, cp), (p58x2, cp), (w09, cw), (w23, cw), (s32, c0), (s64, c0)])
    run("focus", sec_focus, [p24, p58, w09, w23])
    run("spots", sec_spots, [(p58, [[-0.12, 0, 0.3], [0.12, 0, 0.3], [0, 0.2, 0.3]], 0.3),
                             (w23, [[-0.1, 0, 0.3], [0.1, 0, 0.3], [0, 0.15, 0.3]], 0.3)])
    run("quiet", sec_quiet, [(p58, np.array([0, 0, 0.3]), np.array([0.08, 0, 0.3])),
                             (p24, np.array([0, 0, 0.3]), np.array([0.08, 0, 0.3])),
                             (w23, np.array([0, 0, 0.3]), np.array([0.08, 0, 0.3]))], a.quick)
    run("superosc", sec_superosc, [(p58, np.array([0, 0, 0.3])), (w23, np.array([0, 0, 0.3]))])
    run("vortex", sec_vortex, [p58, p24])
    run("pol", sec_pol, 5.8e9)
    kinds = ["line", "ring", "link", "trefoil"]
    run("curves", sec_curves, [(p58, cp), (p58x2, cp), (w23, cw), (s32, c0), (s64, c0)], a.quick, kinds)
    run("room", sec_room, a.quick)
    say(f"\ntotal {time.time() - t0:.0f} s")
    (OUT / ("run_quick.log" if a.quick else "run.log")).write_text("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
