"""Is there a benefit to a 2-D or 3-D lattice? Square versus triangular panels, lattices versus rings and
irregular layouts for the 16 whips, and what periodic 3-D lattices (cubic, body-/face-centred, hexagonal
close-packed) do to the array factor.

A periodic array has an array factor that repeats in wavenumber space: it returns to full height at every
point G of the reciprocal lattice. Two consequences, both set by the shortest G (in cycles per metre):

  beam steering   a grating lobe enters real space when |G| <= (1 + sin(scan)) / wavelength
  direction find  two directions inside a field of view of half-angle t are indistinguishable when
                  |G| <= 2 sin(t) / wavelength (their wave vectors differ by exactly G)

so the lattice whose reciprocal lattice has the longest shortest vector, per element, wins: the equilateral
triangle in 2-D (Sharp 1961), the body-centred cubic in 3-D (Petersen and Middleton 1962).

Parts (all outputs under experiments/lattice/):
  A  panel, 16 active elements: square 25 mm, the same cards staggered by half a pitch, equilateral
     triangular lattices; grating-lobe-free frequency against scan, sidelobes, what the 4-bit delays lose,
     and a 2-D direction-finding Monte Carlo (same error model as evolve_layout.check).
  B  16 whips on the 500 mm plate, 0.4-2.3 GHz, azimuth DF: square and triangular lattices against the
     rings and searched layouts of experiments/evolve/results.json (evolve_layout.check, same conditions).
  C  3-D lattices: shortest reciprocal vector and elements per volume, and the ambiguity level of
     16-element pieces of each lattice against spacing in wavelengths.

    python experiments/lattice_study.py            # about 3 minutes
    python experiments/lattice_study.py --quick
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from evolve_layout import check, ring  # noqa: E402

C0 = 299_792_458.0
OUT = Path(__file__).resolve().parent / "lattice"
TTD_PS = (38.5, 52.5, 117.5, 139.0)      # rev B's four delay sections (REDESIGN.md, section 3)


# ---- lattice theory -----------------------------------------------------------------------------------------
def reciprocal(prim: np.ndarray, basis=None, nmax: int = 4):
    """Reciprocal-lattice points (cycles per unit length) of a lattice with primitive vectors `prim` (rows)
    and an optional basis; returns (|G|, relative height of the array factor there), shortest first."""
    prim = np.asarray(prim, float)
    dim = len(prim)
    basis = np.zeros((1, dim)) if basis is None else np.asarray(basis, float)
    B = np.linalg.inv(prim).T
    rows = []
    for m in itertools.product(range(-nmax, nmax + 1), repeat=dim):
        if not any(m):
            continue
        G = np.asarray(m) @ B
        s = abs(np.exp(2j * np.pi * basis @ G).sum()) / len(basis)
        rows.append((float(np.linalg.norm(G)), float(s)))
    rows.sort()
    return rows, abs(np.linalg.det(prim)) / len(basis)


def lattice_table():
    s3 = np.sqrt(3)
    lat = {
        "2-D square": (np.eye(2), None),
        "2-D staggered cards (pitch = card spacing)": ([[1, 0], [0.5, 1]], None),
        "2-D equilateral triangular": ([[1, 0], [0.5, s3 / 2]], None),
        "3-D simple cubic": (np.eye(3), None),
        "3-D body-centred cubic": (np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]) / s3, None),
        "3-D face-centred cubic": (np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]) / np.sqrt(2), None),
        "3-D hexagonal close-packed": ([[1, 0, 0], [0.5, s3 / 2, 0], [0, 0, np.sqrt(8 / 3)]],
                                       [[0, 0, 0], [0.5, s3 / 6, np.sqrt(2 / 3)]]),
    }
    out = {}
    for name, (prim, basis) in lat.items():
        rows, cell = reciprocal(np.asarray(prim, float), basis)
        dim = len(prim)
        nn = 1.0                                             # all are written with nearest-neighbour distance 1
        g_full = min(g for g, s in rows if s > 0.999)
        partial = [(g, s) for g, s in rows if 0.05 < s <= 0.999 and g < g_full - 1e-9]
        # elements needed per unit area/volume when the lattice is scaled to the same shortest full-height G
        dens = (1 / cell) / g_full ** dim                    # relative to square/cubic = 1
        out[name] = dict(shortest_G_times_spacing=g_full, cell_per_element=cell, nn=nn,
                         elements_vs_square_or_cubic=dens,
                         spacing_gain_pct=100 * (g_full - 1),
                         partial_lobes=[dict(G_times_spacing=g, height_db=20 * np.log10(s)) for g, s in partial[:1]])
    return out


def f_grating(prim_m: np.ndarray, scan_deg: float, plane_deg: float | None = None) -> float:
    """Highest frequency with no grating-lobe peak in real space while the beam scans to scan_deg: anywhere
    on that cone (plane_deg None) or in the plane at azimuth plane_deg."""
    rowsB = np.linalg.inv(np.asarray(prim_m, float)).T
    G = np.array([np.asarray(m) @ rowsB for m in itertools.product(range(-3, 4), repeat=2) if any(m)])
    s = np.sin(np.radians(scan_deg))
    az = np.radians(np.arange(0, 360, 1.0)) if plane_deg is None else np.radians([plane_deg, plane_deg + 180.0])
    u0 = s * np.column_stack([np.cos(az), np.sin(az)])                       # (A, 2)
    d = u0 @ G.T                                                              # (A, G)
    g2 = (G ** 2).sum(1)[None, :]
    if s > 0.99999:
        k = np.where(d < -1e-12, -g2 / (2 * np.minimum(d, -1e-12)), np.inf)
    else:
        k = (d + np.sqrt(d ** 2 + (1 - s ** 2) * g2)) / (1 - s ** 2)         # |k u0 + G| = k
    return float(C0 * k.min())


# ---- array factor on rays ---------------------------------------------------------------------------------
def beyond_first_null(vals: np.ndarray) -> np.ndarray:
    """vals (rays, R) sampled outward from a peak; True from the first local minimum of each ray onward."""
    d = np.diff(vals, axis=1)
    rising = np.concatenate([np.zeros((len(vals), 1), bool), d > 0], axis=1)
    return np.cumsum(rising, axis=1) > 0


def af_rays(p: np.ndarray, rmax: float, dr: float = 0.004, ndir: int = 360):
    """|array factor| / N against the wavenumber offset q (cycles per metre): rays out from the origin.
    p (N, 2 or 3) in metres. Returns q radii (R,), values (rays, R)."""
    r = np.arange(0, rmax + dr, dr)
    if p.shape[1] == 2:
        a = np.radians(np.arange(ndir) * 360.0 / ndir)
        dirs = np.column_stack([np.cos(a), np.sin(a)])
    else:                                                    # Fibonacci sphere (half is enough: |AF| is even)
        i = np.arange(ndir) + 0.5
        z = 1 - i / ndir
        a = np.pi * (1 + 5 ** 0.5) * i
        dirs = np.column_stack([np.sqrt(1 - z * z) * np.cos(a), np.sqrt(1 - z * z) * np.sin(a), z])
    proj = dirs @ p.T                                                         # (rays, N)
    vals = np.abs(np.exp(2j * np.pi * r[None, :, None] * proj[:, None, :]).sum(-1)) / len(p)
    return r, vals


def level_within(r, vals, radius):
    """Highest array-factor value outside the main lobe for offsets up to `radius` (0 when the main lobe
    fills that region)."""
    m = beyond_first_null(vals) & (r[None, :] <= radius)
    return float(vals[m].max()) if m.any() else 0.0


# ---- part A: the panel ------------------------------------------------------------------------------------
def panel_layouts():
    d = 0.025
    a = d * 2 / np.sqrt(3)                   # 28.87 mm: equilateral lattice whose card spacing is still 25 mm
    i, j = np.meshgrid(np.arange(4), np.arange(4), indexing="ij")
    i, j = i.ravel(), j.ravel()

    def cards(dx, dy, stagger):              # card i at x = i dx; flare j on it at y = j dy (+ stagger on odd cards)
        p = np.column_stack([i * dx, j * dy + stagger * (i % 2)])
        return p - p.mean(0)

    return {
        "square 25 mm (rev B)": dict(p=cards(d, d, 0), prim=[[d, 0], [0, d]]),
        "staggered cards 25 x 25 mm": dict(p=cards(d, d, d / 2), prim=[[0, d], [d, d / 2]]),
        "triangular 28.9 mm (cards 25 mm apart)": dict(p=cards(d, a, a / 2), prim=[[0, a], [d, a / 2]]),
        "triangular 25 mm (cards 21.7 mm apart)": dict(p=cards(d * np.sqrt(3) / 2, d, d / 2),
                                                       prim=[[0, d], [d * np.sqrt(3) / 2, d / 2]]),
        "square 26.9 mm (same area as triangular 28.9)": dict(p=cards(d * np.sqrt(a / d), d * np.sqrt(a / d), 0),
                                                              prim=[[d * np.sqrt(a / d), 0], [0, d * np.sqrt(a / d)]]),
    }


def ttd_loss(p: np.ndarray, f: float, scan_deg=(15, 30, 45, 60)):
    """Gain lost by steering with the 16 delay states instead of exact delays, over a set of directions; the
    delay span the worst direction needs."""
    states = np.array(sorted({sum(c) for n in range(5) for c in itertools.combinations(TTD_PS, n)})) * 1e-12
    res = {}
    for th in scan_deg:
        loss, span = [], []
        for az in np.radians(np.arange(0, 360, 7.5)):
            u0 = np.sin(np.radians(th)) * np.array([np.cos(az), np.sin(az)])
            tau = -(p @ u0) / C0
            tau = tau - tau.min()
            span.append(tau.max())
            best = 0.0
            for off in np.linspace(-20e-12, max(states[-1] - tau.max(), 0) + 20e-12, 80):
                q = states[np.argmin(np.abs(states[None, :] - (tau + off)[:, None]), axis=1)]
                best = max(best, abs(np.exp(2j * np.pi * f * (q - tau)).sum()) ** 2 / len(p) ** 2)
            loss.append(-10 * np.log10(best))
        res[str(th)] = dict(mean_loss_db=float(np.mean(loss)), worst_loss_db=float(np.max(loss)),
                            span_needed_ps=float(np.max(span) * 1e12))
    res["span_available_ps"] = float(states[-1] * 1e12)
    return res


def panel_df(p, f, *, fov_deg=60.0, sigma_cal=3.0, snr_db=20.0, k_samples=64, trials=600, seed=2024):
    """Two-angle direction finding with a planar array: one source inside the cone fov_deg of broadside,
    evolve_layout.check's error model (calibration residual, K noisy samples, sample covariance, Bartlett
    beamscan over the whole front hemisphere, coarse then fine). Error = angle between true and estimated
    direction; gross = more than 10 deg. Elements are isotropic over the hemisphere."""
    rng = np.random.default_rng(seed)
    n = len(p)
    k = 2 * np.pi * f / C0
    g = np.arange(-1, 1.0001, 0.025)
    U, V = np.meshgrid(g, g, indexing="ij")
    keep = U ** 2 + V ** 2 <= 1.0
    grid = np.column_stack([U[keep], V[keep]])
    Ac = np.exp(1j * k * (grid @ p.T)).astype(np.complex64)                   # (G, N)
    fo = np.arange(-0.025, 0.02501, 0.0025)
    FU, FV = np.meshgrid(fo, fo, indexing="ij")
    fine = np.column_stack([FU.ravel(), FV.ravel()])
    errs = []
    for start in range(0, trials, 100):
        t = min(100, trials - start)
        ct = rng.uniform(np.cos(np.radians(fov_deg)), 1.0, t)
        st, az = np.sqrt(1 - ct ** 2), rng.uniform(0, 2 * np.pi, t)
        u = np.column_stack([st * np.cos(az), st * np.sin(az)])
        a = np.exp(1j * k * (u @ p.T))
        gg = 10 ** (rng.normal(0, 0.5, (t, n)) / 20) * np.exp(1j * np.radians(rng.normal(0, sigma_cal, (t, n))))
        s = (rng.normal(size=(t, k_samples)) + 1j * rng.normal(size=(t, k_samples))) / np.sqrt(2)
        w = ((rng.normal(size=(t, n, k_samples)) + 1j * rng.normal(size=(t, n, k_samples)))
             / np.sqrt(2) * 10 ** (-snr_db / 20))
        x = (a * gg)[:, :, None] * s[:, None, :] + w
        R = x @ x.conj().transpose(0, 2, 1) / k_samples
        pc = np.real(np.einsum("gi,tij,gj->tg", Ac.conj(), R.astype(np.complex64), Ac, optimize=True))
        c0 = grid[np.argmax(pc, axis=1)]                                      # (t, 2)
        fu = c0[:, None, :] + fine[None, :, :]                                # (t, F, 2)
        Af = np.exp(1j * k * (fu @ p.T))
        pf = np.real(np.einsum("tgi,tij,tgj->tg", Af.conj(), R, Af, optimize=True))
        est = fu[np.arange(t), np.argmax(pf, axis=1)]
        est = est / np.maximum(1.0, np.linalg.norm(est, axis=1))[:, None]
        ze = np.sqrt(np.maximum(0, 1 - (est ** 2).sum(1)))
        cosang = (est * u).sum(1) + ze * ct
        errs.append(np.degrees(np.arccos(np.clip(cosang, -1, 1))))
    e = np.concatenate(errs)
    return float(np.median(e)), float(np.percentile(e, 95)), float(np.mean(e > 10) * 100)


def part_a(quick: bool, log):
    out = {}
    freqs = [2.4e9, 3.5e9, 5.0e9, 5.8e9, 6.0e9, 6.5e9, 7.0e9]
    for name, L in panel_layouts().items():
        p, prim = L["p"], np.asarray(L["prim"])
        rows, cell = reciprocal(prim)
        r, vals = af_rays(p, rmax=2 * 7.0e9 / C0, dr=0.1)
        row = dict(positions_m=np.round(p, 5).tolist(), cell_mm2=cell * 1e6,
                   aperture_vs_square25_db=10 * np.log10(cell / 0.025 ** 2),
                   shortest_G_per_m=rows[0][0],
                   f_grating_free_GHz={str(s): f_grating(prim, s) / 1e9 for s in (0, 30, 45, 60, 90)},
                   f_grating_free_60_along_cards_GHz=f_grating(prim, 60, 90.0) / 1e9,
                   f_grating_free_60_across_cards_GHz=f_grating(prim, 60, 0.0) / 1e9,
                   f_df_unambiguous_GHz={str(t): C0 * rows[0][0] / (2 * np.sin(np.radians(t))) / 1e9 for t in (45, 60, 90)},
                   extent_mm=[float(np.ptp(p[:, 0]) * 1e3), float(np.ptp(p[:, 1]) * 1e3)],
                   ttd_5p8GHz=ttd_loss(p, 5.8e9), by_freq={})
        for f in freqs:
            lam = C0 / f
            sl_broad = level_within(r, vals, 1 / lam)                         # beam at broadside, all real space
            sl_scan = level_within(r, vals, (1 + np.sin(np.radians(60))) / lam)   # beam anywhere within 60 deg
            amb = level_within(r, vals, 2 * np.sin(np.radians(60)) / lam)     # DF, both directions within 60 deg
            med, p95, gross = panel_df(p, f, trials=200 if quick else 600)
            hard = panel_df(p, f, trials=200 if quick else 600, sigma_cal=5.0, snr_db=0.0, seed=99)
            row["by_freq"][f"{f / 1e9:.1f}"] = dict(
                sidelobe_broadside_db=20 * np.log10(max(sl_broad, 1e-6)),
                sidelobe_scan60_db=20 * np.log10(max(sl_scan, 1e-6)),
                df_ambiguity_fov60=amb, df_median_deg=med, df_p95_deg=p95, df_gross_pct=gross,
                df_stress_p95_deg=hard[1], df_stress_gross_pct=hard[2])
        out[name] = row
        log(f"A {name}: grating-free to {row['f_grating_free_GHz']['60']:.2f} GHz at 60 deg, "
            f"{row['f_grating_free_GHz']['90']:.2f} at 90; aperture {row['aperture_vs_square25_db']:+.2f} dB")
    return out


# ---- part B: 16 whips ------------------------------------------------------------------------------------
def whip_layouts():
    evo = json.loads((OUT.parent / "evolve" / "results.json").read_text())["designs"]
    g = (np.arange(4) - 1.5)
    X, Y = np.meshgrid(g, g, indexing="ij")
    sq = np.column_stack([X.ravel(), Y.ravel()])

    def hexpatch(a):                         # triangular lattice centred on the (empty) calibration monopole
        pts = [a * np.array([i + 0.5 * j, j * np.sqrt(3) / 2]) for i in range(-4, 5) for j in range(-4, 5)]
        pts = np.array([q for q in pts if 1e-6 < np.hypot(*q) <= 2 * a + 1e-9])   # rings of 6 at a, 1.73a, 2a
        rad = np.hypot(pts[:, 0], pts[:, 1])
        ang = np.degrees(np.arctan2(pts[:, 1], pts[:, 0])) % 360
        drop = (rad > 1.9 * a) & np.isin(np.round(ang), [120, 300])          # 18 points: drop two opposite ones
        return pts[~drop]

    L = {
        "square lattice 4x4, 113 mm": sq * 0.113,
        "square lattice 4x4, 85 mm": sq * 0.085,
        "triangular lattice, 120 mm": hexpatch(0.120),
        "triangular lattice, 75 mm": hexpatch(0.075),
        "ring of 16, r 225 mm": ring(16, 0.225),
        "two rings of 8, r 225 / 112 mm (rev B)": np.array(evo["rings16_225_112mm"]["positions_m"]),
        "searched, genetic (evolved16)": np.array(evo["evolved16"]["positions_m"]),
        "searched, random (random16)": np.array(evo["random16"]["positions_m"]),
    }
    return L


def az_ambiguity(p, f):
    """Highest correlation between the responses to two horizon azimuths outside each other's main lobe."""
    k = 2 * np.pi * f / C0
    az = np.radians(np.arange(0, 360, 1.0))
    A = np.exp(1j * k * (np.column_stack([np.cos(az), np.sin(az)]) @ p.T))
    Cm = np.abs(A @ A.conj().T) / len(p)
    worst = 0.0
    for i in range(0, 360, 2):
        row = np.roll(Cm[i], -i)                                               # peak at index 0
        for half in (row[:181], row[::-1][:180]):
            d = np.diff(half)
            rise = np.nonzero(d > 0)[0]
            if len(rise):
                worst = max(worst, float(half[rise[0]:].max()))
    return worst


def part_b(quick: bool, log):
    freqs = np.arange(0.40e9, 2.3001e9, 0.10e9 if quick else 0.05e9)
    conds = dict(nominal=dict(sigma_cal=3.0, snr_db=20.0, el_max=15.0, seed=12345),
                 elevated=dict(sigma_cal=3.0, snr_db=20.0, el_max=30.0, seed=4242),
                 stress=dict(sigma_cal=5.0, snr_db=0.0, el_max=30.0, seed=777))
    out = dict(freqs_hz=freqs.tolist(), designs={})
    for name, p in whip_layouts().items():
        d = np.linalg.norm(p[:, None] - p[None], axis=2)
        np.fill_diagonal(d, np.inf)
        rad = np.hypot(p[:, 0], p[:, 1])
        row = dict(positions_m=np.round(p, 4).tolist(), min_spacing_mm=float(d.min() * 1e3),
                   r_min_mm=float(rad.min() * 1e3), r_max_mm=float(rad.max() * 1e3),
                   rms_radius_mm=float(np.sqrt((rad ** 2).mean()) * 1e3),
                   ambiguity=[az_ambiguity(p, f) for f in freqs])
        for c, kw in conds.items():
            r = check(p, freqs, trials=300 if quick else 1000, **kw)
            row[c] = dict(median=r[:, 0].tolist(), p95=r[:, 1].tolist(), gross_pct=r[:, 2].tolist())
        amb = np.array(row["ambiguity"])
        bad = freqs[amb > 0.9]
        row["first_freq_ambiguity_over_0p9_GHz"] = float(bad[0] / 1e9) if len(bad) else None
        out["designs"][name] = row
        log(f"B {name}: nominal p95 {min(row['nominal']['p95']):.2f}-{max(row['nominal']['p95']):.2f} deg, "
            f"gross max {max(row['nominal']['gross_pct']):.1f} / {max(row['elevated']['gross_pct']):.1f} / "
            f"{max(row['stress']['gross_pct']):.1f} %, ambiguity max {amb.max():.2f}"
            f" (first > 0.9 at {row['first_freq_ambiguity_over_0p9_GHz']} GHz)")
    return out


# ---- part C: 3-D lattices ---------------------------------------------------------------------------------
def chunk(prim, basis=None, n=16, centre=(0.13, 0.07, 0.03)):
    """The n lattice points nearest a (generic) centre: a compact piece of the lattice."""
    prim = np.asarray(prim, float)
    basis = np.zeros((1, len(prim))) if basis is None else np.asarray(basis, float)
    pts = np.array([np.asarray(m) @ prim + b for m in itertools.product(range(-3, 4), repeat=len(prim)) for b in basis])
    c = np.asarray(centre[:len(prim)])
    return pts[np.argsort(np.linalg.norm(pts - c, axis=1), kind="stable")[:n]]


def part_c(log):
    s3 = np.sqrt(3)
    g = np.arange(4) - 1.5
    X, Y = np.meshgrid(g, g, indexing="ij")
    flat = np.column_stack([X.ravel(), Y.ravel(), np.zeros(16)])
    tri = np.column_stack([X.ravel() * s3 / 2, Y.ravel() + 0.5 * (np.arange(16) // 4 % 2), np.zeros(16)])
    pieces = {
        "planar square 4x4": flat,
        "planar triangular 4x4": tri,
        "simple cubic 2x2x4": np.array(list(itertools.product(range(2), range(2), range(4))), float),
        "simple cubic, 16 nearest": chunk(np.eye(3)),
        "body-centred cubic, 16 nearest": chunk(np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]) / s3),
        "face-centred cubic, 16 nearest": chunk(np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]) / np.sqrt(2)),
        "hexagonal close-packed, 16 nearest": chunk([[1, 0, 0], [0.5, s3 / 2, 0], [0, 0, np.sqrt(8 / 3)]],
                                                    [[0, 0, 0], [0.5, s3 / 6, np.sqrt(2 / 3)]]),
    }
    s_over_lam = [0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 1.0]
    out = dict(spacing_over_wavelength=s_over_lam, pieces={})
    for name, p in pieces.items():
        p = p - p.mean(0)
        r, vals = af_rays(p, rmax=2.0, dr=0.004, ndir=4000)                  # nearest-neighbour spacing = 1
        lev = [level_within(r, vals, 2 * x) for x in s_over_lam]             # any two directions on the sphere
        ext = np.ptp(p, axis=0)
        first = next((x for x in np.arange(0.3, 1.0001, 0.01) if level_within(r, vals, 2 * x) > 0.9), None)
        out["pieces"][name] = dict(extent_in_spacings=ext.tolist(), ambiguity_level=lev,
                                   spacing_over_wavelength_where_level_passes_0p9=first)
        log(f"C {name}: extent {np.round(ext, 2).tolist()}, level " + " ".join(f"{v:.2f}" for v in lev)
            + f" | >0.9 from s/lambda = {first}")
    return out


# ---- figures ----------------------------------------------------------------------------------------------
def figures(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    A, B, Cc = res["panel"], res["whips"], res["lattices_3d"]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for name, row in A.items():
        f = [float(k) for k in row["by_freq"]]
        ax[0].plot(f, [v["df_ambiguity_fov60"] for v in row["by_freq"].values()], marker="o", label=name)
        ax[1].plot(f, [v["df_stress_gross_pct"] for v in row["by_freq"].values()], marker="o", label=name)
        ax[2].plot(f, [v["sidelobe_scan60_db"] for v in row["by_freq"].values()], marker="o", label=name)
    ax[0].set(xlabel="GHz", ylabel="highest false peak / true peak", title="Panel: DF ambiguity level, sources within 60 deg")
    ax[1].set(xlabel="GHz", ylabel="% of errors over 10 deg", title="Panel DF, 0 dB and 5 deg residual: gross errors")
    ax[2].set(xlabel="GHz", ylabel="dB below the beam", title="Panel: worst sidelobe, beam anywhere within 60 deg")
    ax[0].legend(fontsize=7)
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "panel.png", dpi=130)
    plt.close(fig)

    f = np.array(B["freqs_hz"]) / 1e9
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for name, row in B["designs"].items():
        ls = "--" if "lattice" in name else "-"
        ax[0].semilogy(f, row["nominal"]["p95"], ls, label=name)
        ax[1].plot(f, row["elevated"]["gross_pct"], ls, label=name)
        ax[2].plot(f, row["ambiguity"], ls, label=name)
    ax[0].set(xlabel="GHz", ylabel="deg", title="16 whips: 95th percentile error (20 dB, 3 deg residual)")
    ax[1].set(xlabel="GHz", ylabel="%", title="Gross errors, sources 0-30 deg above the horizon")
    ax[2].set(xlabel="GHz", ylabel="highest false peak / true peak", title="Azimuth ambiguity level (horizon)")
    ax[2].legend(fontsize=7)
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "whips.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.2))
    for name, row in Cc["pieces"].items():
        ax.plot(Cc["spacing_over_wavelength"], row["ambiguity_level"], marker="o", ls="--" if "planar" in name else "-", label=name)
    ax.set(xlabel="nearest-neighbour spacing / wavelength", ylabel="highest false peak / true peak",
           title="16-element pieces of lattices: ambiguity level, any two directions")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "lattices_3d.png", dpi=130)
    plt.close(fig)

    L = panel_layouts()
    fig, ax = plt.subplots(1, len(L), figsize=(3.2 * len(L), 3.4))
    for a, (name, v) in zip(ax, L.items()):
        a.plot(v["p"][:, 0] * 1e3, v["p"][:, 1] * 1e3, "o")
        a.set_title(name, fontsize=7)
        a.set_aspect("equal")
        a.set(xlim=(-60, 60), ylim=(-60, 60), xlabel="mm (across cards)")
        a.grid(alpha=0.3)
    ax[0].set_ylabel("mm (along cards)")
    fig.tight_layout()
    fig.savefig(OUT / "panel_layouts.png", dpi=130)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    t0 = time.time()
    res = dict(theory=lattice_table())
    for name, v in res["theory"].items():
        log(f"T {name}: shortest G x spacing {v['shortest_G_times_spacing']:.4f} "
            f"(spacing {v['spacing_gain_pct']:+.1f} %), elements for the same shortest G {v['elements_vs_square_or_cubic']:.3f}"
            + "".join(f", partial lobe {q['height_db']:.1f} dB at {q['G_times_spacing']:.4f}" for q in v["partial_lobes"]))
    # pitch a restricted top frequency and scan range allow (square; triangular is 2/sqrt(3) times this)
    res["max_square_pitch_mm"] = {f"{f:.2f} GHz": {str(s): 1e3 * C0 / (f * 1e9) / (1 + np.sin(np.radians(s)))
                                                   for s in (30, 45, 60, 90)} for f in (2.4835, 3.8, 5.0, 5.85, 6.0)}
    res["panel"] = part_a(args.quick, log)
    res["whips"] = part_b(args.quick, log)
    res["lattices_3d"] = part_c(log)
    log(f"done in {time.time() - t0:.0f} s")
    (OUT / "results.json").write_text(json.dumps(res, indent=1))
    (OUT / "run.log").write_text("\n".join(lines) + "\n")
    figures(res)


if __name__ == "__main__":
    main()
