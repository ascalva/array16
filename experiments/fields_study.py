"""Fields at the extremes, and the electronic equivalents of moving the array.

The owner's thought: currents make fields that "construct and destruct" and form new field structures, and designs
pushed to their extremes do interesting things. This file tests, for the rev B array (REDESIGN.md), the four
"extreme" field structures that an array of currents can make without moving, and says what each costs:

  1  oam      vortex beams (orbital angular momentum) from a ring of N elements: helical phase, the modes a ring
              of 8 or 16 can make, their purity, cone angle (divergence), on-axis null, and the ring-to-ring link
              (its singular values are the OAM modes: ordinary MIMO, and all but l = 0 die with distance).
  2  superdir superdirective end-fire rows of 2-4 whips at 0.4 GHz: gain against spacing, radiation resistance,
              efficiency, Q / bandwidth and excitation-error tolerance; the receive version (SNR gain against
              external and receiver noise); and the 4x4 panel pushed past 4 pi A / lambda^2.
  3  pol      electronic polarisation rotation with dual-polarised elements: what one V-pol panel loses, the axial
              ratio the rev B delay line can make with no new phase shifter, a field vector spinning at kHz
              (= two circular lines), and a Monte Carlo of polarisation measurement with two V and two H channels.
  4  commute  pseudo-Doppler: one receiver commutated around a ring. Its first harmonic is phase mode l = 1 of the
              ring taken in time instead of space; bias against N and radius.

Models (all stated where used; nothing here is a full-wave result):
  * whips: thin centre-fed dipoles in free space, sinusoidal current, self and mutual impedance by the induced-EMF
    integral (Balanis, Antenna Theory, ch. 8). A monopole row over an ideal ground has half the impedances and
    +3 dB directivity. Loss is a series resistance per element (dipole-equivalent ohms).
  * panel flare: element power pattern cos^1.5 in the front half space, nothing behind (rotation_study's estimate),
    no mutual coupling. Rev B's calibration residual: 3 deg, 0.5 dB rms per element.

Run:  OMP_NUM_THREADS=3 .env/bin/python experiments/fields_study.py            (about a minute)
Outputs: experiments/fields/{tables.md, results.json, *.png}
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C0 = 299_792_458.0
ETA = 376.730313668
OUT = Path(__file__).resolve().parent / "fields"
OUT.mkdir(exist_ok=True)
RNG = np.random.default_rng(16)
_LOG: list[str] = []
RES: dict = {}


def log(s=""):
    print(s)
    _LOG.append(s)


def table(title, head, rows):
    log()
    log(title)
    log("| " + " | ".join(head) + " |")
    log("| " + " | ".join("---" for _ in head) + " |")
    for r in rows:
        log("| " + " | ".join(str(x) for x in r) + " |")


def jv(n, x):
    """Bessel J_n(x) by its integral (scipy is not in this environment)."""
    t = (np.arange(4096) + 0.5) * np.pi / 4096
    return float(np.mean(np.cos(n * t - x * np.sin(t))))


def db(x):
    return 10 * np.log10(np.maximum(x, 1e-300))


# =====================================================================================================================
# 1. OAM from a ring
# =====================================================================================================================
def ring_xy(N, r, rot=0.0):
    a = 2 * np.pi * np.arange(N) / N + rot
    return np.c_[r * np.cos(a), r * np.sin(a)], a


def panel_xy(pitch=0.025, n=4):
    c = (np.arange(n) - (n - 1) / 2) * pitch
    X, Y = np.meshgrid(c, c)
    return np.c_[X.ravel(), Y.ravel()]


def af(xy, w, k, th, ph):
    """Array factor for isotropic elements at xy (m), weights w, on a (theta, phi) grid; axis = z (ring normal)."""
    ux = np.sin(th) * np.cos(ph)
    uy = np.sin(th) * np.sin(ph)
    return np.tensordot(np.exp(1j * k * (ux[..., None] * xy[:, 0] + uy[..., None] * xy[:, 1])), w, axes=([-1], [0]))


def oam_stats(xy, phi_el, l, k, elem="iso", nth=361, nph=256):
    """Vortex beam from elements at xy with weights exp(j l phi_el). Returns purity (power fraction of azimuthal
    harmonic l over the front half space, and inside the main cone), peak cone half-angle, peak directivity, the
    directivity of the same elements fed in phase, and the on-axis level."""
    w = np.exp(1j * l * phi_el) / np.sqrt(len(xy))
    th = np.linspace(0, np.pi / 2, nth)
    ph = 2 * np.pi * np.arange(nph) / nph
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    g = {"iso": np.ones_like(th), "flare": np.cos(th) ** 1.5, "whip": np.sin(th) ** 2}[elem]
    F = af(xy, w, k, TH, PH)
    P = np.abs(F) ** 2 * g[:, None]
    cm = np.fft.fft(F, axis=1) / nph                      # azimuthal harmonics per theta
    m = np.fft.fftfreq(nph, 1 / nph).astype(int)
    pl = np.abs(cm[:, m == (l % nph if l >= 0 else l)][:, 0]) ** 2 * g
    pt = np.sum(np.abs(cm) ** 2, axis=1) * g
    wth = np.sin(th)
    dth = th[1] - th[0]
    prad = np.sum(pt * wth) * dth * 2 * np.pi             # front half space only
    # elements radiate into the half space only for 'flare'; iso/whip are symmetric about the plane
    full = prad * (1 if elem == "flare" else 2)
    ring = P.mean(axis=1)
    ipk = int(np.argmax(ring))
    # main cone: out to the first minimum of the phi-averaged pattern beyond the peak
    j = ipk
    while j + 1 < nth and ring[j + 1] < ring[j]:
        j += 1
    sl = slice(0, j + 1)
    F0 = af(xy, np.ones(len(xy)) / np.sqrt(len(xy)), k, TH, PH)
    P0 = np.abs(F0) ** 2 * g[:, None]
    full0 = np.sum(P0.mean(axis=1) * wth) * dth * 2 * np.pi * (1 if elem == "flare" else 2)
    return dict(purity_half=float(np.sum(pl * wth) / np.sum(pt * wth)),
                purity_cone=float(np.sum((pl * wth)[sl]) / np.sum((pt * wth)[sl])),
                cone_deg=float(np.degrees(th[ipk])), null_deg=float(np.degrees(th[j])),
                D_peak_dBi=float(db(4 * np.pi * P.max() / full)),
                D_l0_dBi=float(db(4 * np.pi * P0.max() / full0)),
                axis_dB=float(db(P[0].mean() / P.max())),
                ripple_dB=float(db(P[ipk].max() / max(P[ipk].min(), 1e-30))))


def ring_link_svd(N, r, f, dists):
    """Two identical coaxial rings of N isotropic elements facing each other at distance D: H_mn = exp(-jkd)/d.
    H is circulant, so its singular vectors are the DFT vectors = OAM modes l; returns |sigma_l| (free-space
    amplitude gain per mode, normalised so that one element to one element at D is lambda/(4 pi D))."""
    k = 2 * np.pi * f / C0
    xy, _ = ring_xy(N, r)
    out = {}
    for D in dists:
        d = np.sqrt(D ** 2 + np.sum((xy[:, None, :] - xy[None, :, :]) ** 2, axis=2))
        H = np.exp(-1j * k * d) / d * (C0 / f) / (4 * np.pi)
        lam = np.fft.fft(H[0])                             # eigenvalues of a circulant matrix; index = mode l
        out[D] = np.abs(lam)                               # coherent N-to-N amplitude gain of mode l
    return out


def study_oam():
    log("## 1. Vortex (OAM) beams from a ring")
    res = {}
    # 1a analytic picture: weights exp(j l phi_n) on N elements at radius a give
    #    AF = N sum_q j^(l+qN) J_(l+qN)(k a sin th) exp(j (l+qN) phi): the wanted mode plus aliases l +- N, l +- 2N ...
    rows = []
    for N in (8, 16):
        for ka in (1.0, 2.0, 4.0, 8.0, 11.0):
            r = []
            for l in range(0, N // 2 + 1):
                xy, a = ring_xy(N, ka)                     # k = 1, radius = ka
                s = oam_stats(xy, a, l, 1.0)
                r.append(f"{100 * s['purity_half']:.0f}")
                res[f"iso_N{N}_ka{ka}_l{l}"] = s
            rows.append([N, ka] + r + [""] * (9 - len(r)))
    table("1a. Mode purity over the whole front half space (% of radiated power in azimuthal harmonic l), isotropic "
          "elements, ring of N at radius a. Aliases are l - N and l + N, so |l| = N/2 is never pure (it is a "
          "standing pattern cos(N phi / 2)), and purity falls once k a exceeds about N - |l|",
          ["N", "k a"] + [f"l = {l}" for l in range(9)], rows)

    # 1b the project's rings
    cases = [("plate outer ring, 8 whips, 225 mm", 8, 0.225, "whip", (0.4e9, 0.915e9, 1.42e9, 2.3e9)),
             ("16 whips on one 225 mm ring (hypothetical)", 16, 0.225, "whip", (0.4e9, 0.915e9, 1.42e9, 2.3e9)),
             ("ring of 8 transverse elements (dipoles/patches facing up), 225 mm", 8, 0.225, "iso",
              (0.4e9, 0.915e9, 1.42e9, 2.3e9)),
             ("ring of 16 transverse elements, 225 mm", 16, 0.225, "iso", (0.915e9, 1.42e9, 2.3e9))]
    rows = []
    for nm, N, r, el, fs in cases:
        for f in fs:
            k = 2 * np.pi * f / C0
            xy, a = ring_xy(N, r)
            for l in (1, 2, 3):
                s = oam_stats(xy, a, l, k, el)
                res[f"{nm}|{f/1e9:.3f}|l{l}"] = s
                rows.append([nm, f"{f/1e9:.3f}", f"{k*r:.2f}", l, f"{100*s['purity_cone']:.0f} / {100*s['purity_half']:.0f}",
                             f"{s['cone_deg']:.0f}", f"{s['D_peak_dBi']:.1f}", f"{s['D_l0_dBi']:.1f}",
                             f"{s['ripple_dB']:.1f}"])
    table("1b. Rings on the 500 mm plate. 'whip' = element pattern sin^2 (null along the ring axis, where the vortex "
          "is), 'transverse' = isotropic stand-in for elements polarised in the plate. Purity = main cone / whole "
          "half space; cone = half-angle of the brightest ring (divergence); D = peak directivity of the vortex "
          "beam against the same ring fed in phase (its own peak, which for whips is not on the axis); ripple = "
          "max/min around the brightest ring",
          ["ring", "GHz", "k a", "l", "purity % cone / half", "cone half-angle deg", "D vortex dBi", "D in-phase dBi",
           "ripple dB"], rows)

    # 1c the panel: its 4x4 is not a ring. Use the 12 perimeter elements, the 8 at r = 39.5 mm, the 4 quadrant sums
    P = panel_xy()
    rad = np.hypot(P[:, 0], P[:, 1])
    ang = np.arctan2(P[:, 1], P[:, 0])
    sets = {"all 16, phase l x azimuth of each element": np.ones(16, bool),
            "12 perimeter elements": rad > 0.02,
            "8 at r 39.5 mm": np.abs(rad - 0.0395) < 1e-3,
            "4 quadrant sums (what four receive channels see)": None}
    rows = []
    for f in (2.44e9, 3.5e9, 5.8e9):
        k = 2 * np.pi * f / C0
        for nm, msk in sets.items():
            for l in (1, 2):
                if msk is None:
                    # each quadrant fed uniformly, quadrant q gets phase l * (45 + 90 q) deg
                    q = (np.round((ang - np.pi / 4) / (np.pi / 2)).astype(int)) % 4
                    w_phi = np.pi / 4 + q * np.pi / 2
                    s = oam_stats(P, w_phi, l, k, "flare")
                else:
                    s = oam_stats(P[msk], ang[msk], l, k, "flare")
                res[f"panel|{nm}|{f/1e9:.2f}|l{l}"] = s
                rows.append([f"{f/1e9:.2f}", nm, l, f"{100*s['purity_cone']:.0f} / {100*s['purity_half']:.0f}",
                             f"{s['cone_deg']:.0f}", f"{s['D_peak_dBi']:.1f}", f"{s['D_l0_dBi']:.1f}",
                             "< -100" if s['axis_dB'] < -100 else f"{s['axis_dB']:.0f}"])
    table("1c. The 4x4 panel as a vortex source (flare pattern cos^1.5). l = 1 from the four quadrants is exactly the "
          "monopulse difference pair (A - C) + j (B - D); l = 2 from four quadrants is the 'quadrupole' pattern "
          "A - B + C - D and is not a vortex (purity 50 %: l = +2 and -2 together)",
          ["GHz", "elements used", "l", "purity % cone / half", "cone half-angle deg", "D vortex dBi",
           "D in-phase dBi", "on-axis dB"], rows)

    # 1d link: ring to ring
    N, r, f = 8, 0.225, 2.3e9
    dists = [0.3, 1.0, 3.0, 10.0, 30.0, 100.0]
    sv = ring_link_svd(N, r, f, dists)
    lam = C0 / f
    rows = []
    for D in dists:
        s = sv[D]
        ref = N * lam / (4 * np.pi * D)                    # all-in-phase far-field value: N^2 power gain over 1x1
        rows.append([D, f"{D / (2 * (2 * r) ** 2 / lam):.2f}"] +
                    [f"{2 * db(s[l] / ref):.1f}" for l in range(0, 5)])
    RES["oam_link"] = {str(D): (2 * db(sv[D][:5] / (N * lam / (4 * np.pi * D)))).tolist() for D in dists}
    # fitted slopes between 30 and 100 m
    sl = [float((2 * db(sv[100.0][l]) - 2 * db(sv[30.0][l])) / np.log10(100 / 30) / 10) for l in range(5)]
    RES["oam_link_slope"] = sl
    table(f"1d. Two coaxial rings of {N} at r = {r*1e3:.0f} mm, {f/1e9:.1f} GHz, facing each other (Rayleigh distance "
          f"2 (2r)^2 / lambda = {2*(2*r)**2/lam:.1f} m). The channel matrix is circulant, so its singular vectors ARE "
          "the OAM modes; power gain of each mode in dB relative to the in-phase far-field value (N^2 x Friis). "
          f"Far-field slopes fitted 30-100 m: power ~ D^({', '.join(f'{x:.1f}' for x in sl)}) for l = 0..4 "
          "(theory -2 (|l| + 1))",
          ["D m", "D / Rayleigh", "l = 0", "l = 1", "l = 2", "l = 3", "l = 4"], rows)

    # figure
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    th = np.linspace(0, np.pi / 2, 361)
    k = 2 * np.pi * 2.3e9 / C0
    for N, ls in ((8, "-"), (16, "--")):
        xy, a = ring_xy(N, 0.225)
        for l, c in zip((0, 1, 2, 3), ("k", "C0", "C1", "C2")):
            w = np.exp(1j * l * a) / np.sqrt(N)
            F = af(xy, w, k, th[:, None], (2 * np.pi * np.arange(128) / 128)[None, :])
            ax[0].plot(np.degrees(th), db(np.mean(np.abs(F) ** 2, axis=1) / N), ls, color=c,
                       label=f"N = {N}, l = {l}")
    ax[0].set_ylim(-40, 3)
    ax[0].set_xlabel("angle from ring axis (deg)")
    ax[0].set_ylabel("azimuth-averaged power, dB re N-element peak")
    ax[0].set_title("225 mm ring, 2.3 GHz (k a = 10.8): vortex cones and on-axis null", fontsize=9)
    ax[0].legend(fontsize=7, ncol=2)
    kas = np.linspace(0.3, 16, 120)
    for N, ls in ((8, "-"), (16, "--")):
        for l, c in zip((1, 2, 3), ("C0", "C1", "C2")):
            pur = []
            for ka in kas:
                xy, a = ring_xy(N, ka)
                pur.append(oam_stats(xy, a, l, 1.0, nth=91, nph=64)["purity_half"])
            ax[1].plot(kas, 100 * np.array(pur), ls, color=c, label=f"N = {N}, l = {l}")
    for f_, nm in ((0.4e9, "0.4"), (0.915e9, "0.915"), (1.42e9, "1.42"), (2.3e9, "2.3 GHz")):
        x = 2 * np.pi * f_ / C0 * 0.225
        ax[1].axvline(x, color="0.8", lw=0.8)
        ax[1].text(x, 3, nm, fontsize=7, rotation=90, va="bottom")
    ax[1].set_xlabel("k a  (225 mm ring: marked frequencies)")
    ax[1].set_ylabel("power in harmonic l, whole half space (%)")
    ax[1].set_title("Mode purity: fails once k a > N - |l|", fontsize=9)
    ax[1].legend(fontsize=7)
    Ds = np.geomspace(0.2, 200, 60)
    sv = ring_link_svd(8, 0.225, 2.3e9, Ds)
    for l in range(5):
        ax[2].loglog(Ds, [(sv[D][l] / (8 * lam / (4 * np.pi * D))) ** 2 for D in Ds], label=f"l = {l}")
    ax[2].axvline(2 * 0.45 ** 2 / lam, color="0.6", ls=":")
    ax[2].set_ylim(1e-12, 10)
    ax[2].set_xlabel("distance between rings (m)")
    ax[2].set_ylabel("mode power gain re in-phase far field")
    ax[2].set_title("Ring-to-ring link, 8 + 8 at 225 mm, 2.3 GHz: OAM modes = MIMO singular modes", fontsize=9)
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "oam.png", dpi=130)
    plt.close(fig)
    RES["oam"] = res


# =====================================================================================================================
# 2. Superdirectivity
# =====================================================================================================================
_GL = np.polynomial.legendre.leggauss(400)


def z_mutual(k, L, d):
    """Mutual impedance (referred to the current maxima) of two parallel side-by-side thin dipoles of length L,
    sinusoidal current, axes d apart: induced-EMF integral with the closed-form E_z of the sinusoidal dipole
    (Balanis eq. 8-58..8-68). d = wire radius gives the self impedance."""
    x, w = _GL
    z = 0.25 * L * (x + 1)                                  # 0 .. L/2 (the integrand is even and has a kink at 0)
    wz = 0.5 * L * w                                        # = 2 x (L/4) w
    R1 = np.sqrt(d * d + (z - L / 2) ** 2)
    R2 = np.sqrt(d * d + (z + L / 2) ** 2)
    r = np.sqrt(d * d + z * z)
    Ez = -1j * ETA / (4 * np.pi) * (np.exp(-1j * k * R1) / R1 + np.exp(-1j * k * R2) / R2
                                    - 2 * np.cos(k * L / 2) * np.exp(-1j * k * r) / r)
    return -np.sum(Ez * np.sin(k * (L / 2 - np.abs(z))) * wz)


def z_matrix(f, L, xs, a=2e-3):
    """Input-referred impedance matrix of parallel dipoles at positions xs along a line."""
    k = 2 * np.pi * f / C0
    n = len(xs)
    Z = np.empty((n, n), complex)
    for i in range(n):
        for j in range(n):
            if i == j:
                # reactance at the wire radius; resistance at zero separation, so that Re Z is exactly the
                # pattern-overlap matrix (positive semi-definite), which the supergain optimum depends on
                Z[i, j] = z_mutual(k, L, 1e-7).real + 1j * z_mutual(k, L, a).imag
            else:
                Z[i, j] = z_mutual(k, L, abs(xs[i] - xs[j]))
    return Z / np.sin(k * L / 2) ** 2


def endfire(f, L, xs, Rloss=0.0, a=2e-3, I=None, df=1e-3):
    """Maximum end-fire gain of a row of parallel dipoles (direction +x, in the plane broadside to the wires).
    Directivity D = eta (1 - cos(kL/2))^2 |s^H I_m|^2 / (pi I_m^H R_m I_m) with I_m the current maxima; here all in
    input-referred quantities. Returns gain, directivity, efficiency, array Q, currents."""
    k = 2 * np.pi * f / C0
    Z = z_matrix(f, L, xs, a)
    R = Z.real
    s = np.exp(-1j * k * np.asarray(xs))                    # field of element n toward +x carries exp(+jk x_n)
    scale = ETA * (1 - np.cos(k * L / 2)) ** 2 / (np.pi * np.sin(k * L / 2) ** 2)
    if I is None:
        I = np.linalg.solve(R + Rloss * np.eye(len(xs)), s)
    num = np.abs(np.vdot(s, I)) ** 2
    prad = np.real(np.vdot(I, R @ I))
    ploss = Rloss * np.real(np.vdot(I, I))
    # Q of the array with these currents held: omega * sum |I|^2 (X' + |X|/omega) / (2 sum |I|^2 R), each element
    # tuned by a lossless series reactance (Yaghjian & Best 2005 form, applied to the active impedances)
    Zp = z_matrix(f * (1 + df), L, xs, a)
    Zm = z_matrix(f * (1 - df), L, xs, a)
    Za = (Z @ I) / I
    dXa = ((Zp @ I) / I - (Zm @ I) / I).imag / (2 * df)     # omega dX/domega
    p = np.abs(I) ** 2
    Q = float(np.sum(p * (dXa + np.abs(Za.imag))) / (2 * np.sum(p * (Za.real + Rloss))))
    return dict(G=scale * num / (prad + ploss), D=scale * num / prad, eff=prad / (prad + ploss), Q=Q, I=I,
                Reff=prad / np.real(np.vdot(I, I)), Za=Za, R=R, s=s, scale=scale)


def tolerance(f, L, xs, I, Rloss, sig_db, sig_deg, trials=4000):
    """Mean and 10th-percentile gain when each element current has independent amplitude and phase errors."""
    k = 2 * np.pi * f / C0
    Z = z_matrix(f, L, xs)
    R = Z.real
    s = np.exp(-1j * k * np.asarray(xs))
    scale = ETA * (1 - np.cos(k * L / 2)) ** 2 / (np.pi * np.sin(k * L / 2) ** 2)
    n = len(xs)
    e = 10 ** (RNG.normal(0, sig_db, (trials, n)) / 20) * np.exp(1j * np.radians(RNG.normal(0, sig_deg, (trials, n))))
    Ie = I[None, :] * e
    num = np.abs(Ie @ s.conj()) ** 2
    den = np.real(np.einsum("ti,ij,tj->t", Ie.conj(), R + Rloss * np.eye(n), Ie))
    g = scale * num / den
    return float(db(np.mean(g))), float(db(np.percentile(g, 10)))


def rx_snr_gain(f, L, xs, t_ext_over_t_rx):
    """Receive: each whip has its own LNA; noise = external (isotropic, coherence = R / R11) + receiver (independent).
    SNR gain of the best weights over one element, end-fire source. Open-circuit model: an estimate."""
    k = 2 * np.pi * f / C0
    R = z_matrix(f, L, xs).real
    Rn = R / R[0, 0]
    s = np.exp(-1j * k * np.asarray(xs))
    t = t_ext_over_t_rx
    M = t * Rn + np.eye(len(xs))
    return float(np.real(np.vdot(s, np.linalg.solve(M, s))) * (t + 1))


def study_superdir():
    log()
    log("## 2. Superdirectivity")
    f = 0.4e9
    lam = C0 / f
    L = 0.475 * lam                                         # near-resonant dipole (monopole 178 mm over ground)
    Z11 = z_matrix(f, L, [0.0])[0, 0]
    Zhw = z_matrix(f, 0.5 * lam, [0.0], a=1e-5 * lam)[0, 0]
    Z12 = z_matrix(f, 0.5 * lam, [0.0, 0.5 * lam], a=1e-5 * lam)[0, 1]
    log(f"check: half-wave thin dipole Z11 = {Zhw.real:.1f} {Zhw.imag:+.1f}j ohm (textbook 73.1 + 42.5j); "
        f"side-by-side at 0.5 wl Z12 = {Z12.real:.1f} {Z12.imag:+.1f}j (textbook -12.5 - 29.9j); "
        f"single-element directivity {db(endfire(f, 0.5*lam, [0.0], a=1e-5*lam)['D']):.2f} dBi (2.15)")
    log(f"whip model at 0.4 GHz: L = 0.475 wl = {L*1e3:.0f} mm dipole (monopole {L*500:.0f} mm), radius 2 mm, "
        f"Z11 = {Z11.real:.1f} {Z11.imag:+.1f}j ohm, single-element Q = {endfire(f, L, [0.0])['Q']:.1f}")
    RES["sd_check"] = dict(Z11=[Zhw.real, Zhw.imag], Z12=[Z12.real, Z12.imag])

    # 2a gain against spacing, lossless and lossy, N = 2, 3, 4
    ds = [0.5, 0.35, 0.25, 0.2, 0.15, 0.1, 0.07, 0.05, 0.03]
    losses = [0.0, 0.5, 2.0, 5.0]
    res = {}
    rows = []
    for N in (2, 3, 4):
        for d in ds:
            xs = np.arange(N) * d * lam
            e0 = endfire(f, L, xs, 0.0)
            # uniform-amplitude ordinary end-fire (phase = -k d per element) for comparison
            Iord = np.exp(-1j * 2 * np.pi * d * np.arange(N))
            eo = endfire(f, L, xs, 0.0, I=Iord)
            g = [endfire(f, L, xs, rl) for rl in losses]
            t1 = tolerance(f, L, xs, e0["I"], 0.0, 0.5, 3.0)
            t2 = tolerance(f, L, xs, e0["I"], 0.0, 0.1, 1.0)
            # design that knows about the errors: loading equal to the error variance times R11
            sig2 = (np.radians(3.0)) ** 2 + (0.5 / 8.686) ** 2
            er = endfire(f, L, xs, sig2 * Z11.real)
            t3 = tolerance(f, L, xs, er["I"], 0.0, 0.5, 3.0)
            res[f"N{N}_d{d}"] = dict(D=float(db(e0["D"])), D_ord=float(db(eo["D"])), Q=e0["Q"], Reff=float(e0["Reff"]),
                                     G_loss={str(rl): float(db(x["G"])) for rl, x in zip(losses, g)},
                                     eff_loss={str(rl): float(x["eff"]) for rl, x in zip(losses, g)},
                                     tol_3deg=t1, tol_1deg=t2, tol_robust=t3,
                                     Imax_over_Imin=float(np.max(np.abs(e0["I"])) / np.min(np.abs(e0["I"]))))
            rows.append([N, d, f"{d*lam*1e3:.0f}", f"{db(eo['D']):.1f}", f"{db(e0['D']):.1f}", f"{e0['Reff']:.2f}",
                         f"{e0['Q']:.0f}", f"{100*0.707/max(e0['Q'], 1e-9):.1f}",
                         " / ".join(f"{db(x['G']):.1f}" for x in g[1:]),
                         " / ".join(f"{100*endfire(f, L, xs, 0.0, I=e0['I'])['eff']:.0f}" if rl == 0 else
                                    f"{100*_eff(e0, rl):.0f}" for rl in losses[1:]),
                         f"{t2[0]:.1f} ({t2[1]:.1f})", f"{t1[0]:.1f} ({t1[1]:.1f})", f"{t3[0]:.1f} ({t3[1]:.1f})"])
    table("2a. End-fire row of N whips at 0.4 GHz (dipole model in free space; over an ideal ground add 3 dB and halve "
          "the ohms). 'ordinary' = equal amplitudes, phase -k d per element. 'max' = lossless optimum currents. "
          "R_eff = radiation resistance seen by the optimum current set (ohm, per unit current; a lone whip has 63). "
          "Q and the 2:1-VSWR bandwidth 0.707/Q hold the optimum currents fixed and tune each element. Gains with "
          "loss re-optimise the currents for that loss; efficiency is for the LOSSLESS-optimum currents. Tolerance: "
          "mean (10th percentile) gain of the lossless-optimum currents with random element errors; 'robust' = "
          "currents designed for the 3 deg / 0.5 dB errors, then given them",
          ["N", "d / wl", "d mm", "ordinary dBi", "max D dBi", "R_eff ohm", "Q", "BW %",
           "gain dBi at 0.5 / 2 / 5 ohm loss", "efficiency % at 0.5 / 2 / 5 ohm", "1 deg, 0.1 dB", "3 deg, 0.5 dB",
           "robust design, 3 deg"], rows)
    RES["superdir_whips"] = res

    # 2b receive
    rows = []
    for N in (2, 3, 4):
        for d in (0.5, 0.25, 0.15, 0.1, 0.05):
            xs = np.arange(N) * d * lam
            rows.append([N, d] + [f"{db(rx_snr_gain(f, L, xs, t)):.1f}" for t in (0.0, 0.3, 1.0, 3.0, 10.0, 1e3)])
    RES["superdir_rx"] = rows
    table("2b. Receive version (each whip has its own LNA, as in rev B): SNR gain in dB of the best weights over ONE "
          "whip, end-fire source, against the ratio of external (diffuse) noise to receiver noise. Receiver-noise "
          "limited (ratio 0) the gain is 10 log N whatever the spacing; external-noise limited it tends to the "
          "directivity gain. At 0.4 GHz the ratio is about 1-10 (estimate: LNA 84 K + losses, sky and man-made "
          "100-900 K); at 2.3 GHz and up it is below 1",
          ["N", "d / wl", "T_ext/T_rx = 0", "0.3", "1", "3", "10", "1000"], rows)

    # 2c the 4x4 panel past 4 pi A / lambda^2
    P = panel_xy()
    th = np.linspace(0, np.pi / 2, 181)
    ph = np.linspace(0, 2 * np.pi, 361)[:-1]
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    g = np.cos(TH) ** 1.5
    dO = np.sin(TH) * (th[1] - th[0]) * (ph[1] - ph[0])
    rows = []
    res = {}
    for fz in (2.3e9, 2.44e9, 3.5e9, 5.8e9):
        k = 2 * np.pi * fz / C0
        ux, uy = np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH)
        A = np.exp(1j * k * (ux[..., None] * P[:, 0] + uy[..., None] * P[:, 1]))        # (th, ph, 16)
        Bm = np.einsum("tp,tpi,tpj->ij", g * dO, A, A.conj()) / (4 * np.pi)
        for scan in (0.0, 45.0):
            t0 = np.radians(scan)
            s = np.exp(1j * k * np.sin(t0) * P[:, 0])
            g0 = np.cos(t0) ** 1.5

            def D_of(w):
                return g0 * np.abs(np.vdot(w, s)) ** 2 / np.real(np.vdot(w, Bm @ w))

            def wng(w):
                return np.abs(np.vdot(w, s)) ** 2 / np.real(np.vdot(w, w)) / 16.0

            wu = s.copy()
            wm = np.linalg.solve(Bm + 1e-9 * np.eye(16), s)
            # robust: loading for 3 deg / 0.5 dB
            sig2 = np.radians(3.0) ** 2 + (0.5 / 8.686) ** 2
            wr = np.linalg.solve(Bm + sig2 * np.trace(Bm).real / 16 * np.eye(16), s)

            def mc(w, sd=0.5, sp=3.0, T=3000):
                e = 10 ** (RNG.normal(0, sd, (T, 16)) / 20) * np.exp(1j * np.radians(RNG.normal(0, sp, (T, 16))))
                we = w[None, :].conj() * e                    # errors sit on the element responses
                num = np.abs(we @ s) ** 2
                den = np.real(np.einsum("ti,ij,tj->t", we, Bm.T, we.conj()))
                return float(db(np.mean(g0 * num / den)))

            area = 4 * np.pi * 0.1 * 0.1 / (C0 / fz) ** 2 * np.cos(t0)
            res[f"{fz/1e9:.2f}_{scan:.0f}"] = dict(bound=float(db(area)), uni=float(db(D_of(wu))),
                                                   max=float(db(D_of(wm))), wng_max=float(db(wng(wm))),
                                                   max_err=mc(wm), rob=float(db(D_of(wr))),
                                                   wng_rob=float(db(wng(wr))), rob_err=mc(wr), uni_err=mc(wu))
            rows.append([f"{fz/1e9:.2f}", f"{scan:.0f}", f"{0.025/(C0/fz):.2f}", f"{db(area):.1f}",
                         f"{db(D_of(wu)):.1f}", f"{mc(wu):.1f}", f"{db(D_of(wm)):.1f}", f"{db(wng(wm)):.0f}",
                         f"{mc(wm):.1f}", f"{db(D_of(wr)):.1f}", f"{db(wng(wr)):.1f}", f"{mc(wr):.1f}"])
    RES["superdir_panel"] = res
    table("2c. The 4x4 panel (25 mm pitch, element pattern cos^1.5, no coupling: an idealised upper bound on what "
          "weights alone can do). Bound = 4 pi A cos(scan) / lambda^2 with A = 100 x 100 mm. WNG = white-noise gain "
          "relative to uniform weights (0 dB = full 12 dB array gain against receiver noise). 'with errors' = mean "
          "directivity with 3 deg / 0.5 dB element errors",
          ["GHz", "scan deg", "pitch / wl", "bound dBi", "uniform D dBi", "uniform with errors", "max D dBi",
           "its WNG dB", "max with errors", "robust D dBi", "its WNG dB", "robust with errors"], rows)

    # 2d Chu / Harrington numbers for the project
    rows = []
    for nm, a_, fz in (("one 178 mm whip + image (sphere radius 178 mm)", 0.178, 0.4e9),
                       ("4-whip ring r 133 mm + whips (radius 222 mm)", 0.222, 0.4e9),
                       ("2 whips 75 mm apart + image (radius 182 mm)", 0.182, 0.4e9),
                       ("500 mm plate + whips (radius 290 mm)", 0.290, 0.4e9),
                       ("500 mm plate + whips (radius 290 mm)", 0.290, 2.3e9),
                       ("panel active 4x4, 100 mm square (radius 71 mm)", 0.0707, 2.3e9),
                       ("panel active 4x4 (radius 71 mm)", 0.0707, 5.8e9)):
        ka = 2 * np.pi * fz / C0 * a_
        gh = ka ** 2 + 2 * ka
        qchu = 1 / ka ** 3 + 1 / ka
        rows.append([nm, f"{fz/1e9:.1f}", f"{ka:.2f}", f"{db(gh):.1f}", f"{qchu:.2f}"])
    table("2d. Harrington's 'normal gain' (ka)^2 + 2 ka for the smallest sphere that holds the antenna, and the Chu "
          "lower bound on Q of a small single-mode antenna, 1/(ka)^3 + 1/(ka) (McLean 1996). Gain above the normal "
          "gain is superdirective by definition; the Chu Q only binds for k a below about 1",
          ["antenna", "GHz", "k a", "normal gain dBi", "Chu Q bound"], rows)

    # figure
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    dd = np.linspace(0.02, 0.6, 60)
    for N, c in ((2, "C0"), (3, "C1"), (4, "C2")):
        Dm, G2, Qs, Go = [], [], [], []
        for d in dd:
            xs = np.arange(N) * d * lam
            e0 = endfire(f, L, xs, 0.0)
            Dm.append(db(e0["D"]))
            Qs.append(e0["Q"])
            G2.append(db(endfire(f, L, xs, 2.0)["G"]))
            Go.append(db(endfire(f, L, xs, 0.0, I=np.exp(-1j * 2 * np.pi * d * np.arange(N)))["D"]))
        ax[0].plot(dd, Dm, color=c, label=f"N = {N} lossless optimum")
        ax[0].plot(dd, G2, "--", color=c, label=f"N = {N}, 2 ohm loss per element")
        ax[0].plot(dd, Go, ":", color=c, label=f"N = {N} ordinary end-fire")
        ax[1].semilogy(dd, np.maximum(Qs, 0.1), color=c, label=f"N = {N}")
    ax[0].set_xlabel("spacing / wavelength")
    ax[0].set_ylabel("end-fire gain (dBi, dipole model)")
    ax[0].set_title("Superdirective whip rows at 0.4 GHz", fontsize=9)
    ax[0].legend(fontsize=6.5)
    ax[1].set_xlabel("spacing / wavelength")
    ax[1].set_ylabel("Q with the optimum currents held")
    ax[1].set_title("Stored energy: Q (2:1 bandwidth = 0.707 / Q)", fontsize=9)
    ax[1].legend(fontsize=7)
    sig = np.linspace(0, 6, 25)
    for N, c in ((2, "C0"), (3, "C1"), (4, "C2")):
        for d, ls in ((0.1, "-"), (0.25, "--")):
            xs = np.arange(N) * d * lam
            I = endfire(f, L, xs, 0.0)["I"]
            ax[2].plot(sig, [tolerance(f, L, xs, I, 0.0, s_ / 6.0, s_, 1500)[0] for s_ in sig], ls, color=c,
                       label=f"N = {N}, d = {d} wl")
    ax[2].axvline(3.0, color="0.7", ls=":")
    ax[2].set_xlabel("rms phase error per element (deg), with dB error = deg / 6")
    ax[2].set_ylabel("mean end-fire gain (dBi)")
    ax[2].set_title("Excitation-error tolerance (rev B residual: 3 deg, 0.5 dB)", fontsize=9)
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "superdirective.png", dpi=130)
    plt.close(fig)


def _eff(e0, rl):
    I = e0["I"]
    prad = np.real(np.vdot(I, e0["R"] @ I))
    return prad / (prad + rl * np.real(np.vdot(I, I)))


# =====================================================================================================================
# 3. Electronic polarisation rotation
# =====================================================================================================================
TTD_BITS = np.array([38.5, 52.5, 117.5, 139.0]) * 1e-12     # REDESIGN section 3


def axial_ratio_db(ev, eh):
    """Axial ratio (dB, >= 0) of the field with complex components ev, eh."""
    s0 = np.abs(ev) ** 2 + np.abs(eh) ** 2
    s3 = 2 * np.imag(np.conj(ev) * eh)
    x = np.clip(np.abs(s3) / s0, 0, 1)
    chi = 0.5 * np.arcsin(x)                                # ellipticity angle
    return 20 * np.log10(1 / np.maximum(np.tan(chi), 1e-9))


def stokes(ev, eh):
    return np.array([np.mean(np.abs(ev) ** 2 + np.abs(eh) ** 2), np.mean(np.abs(ev) ** 2 - np.abs(eh) ** 2),
                     2 * np.mean(np.real(np.conj(ev) * eh)), 2 * np.mean(np.imag(np.conj(ev) * eh))])


def study_pol():
    log()
    log("## 3. Electronic polarisation rotation with dual-polarised elements")
    # 3a what the single V-pol panel loses
    psi = np.array([0, 15, 30, 45, 60, 75, 85, 90])
    rows = [["V-pol only (rev B)"] + [f"{-db(np.cos(np.radians(p))**2):.1f}" if p < 90 else "inf" for p in psi] +
            ["3.0", "3.0 (and half the noise-free information)"],
            ["V + H, weights (cos psi, sin psi) or (1, -+j)"] + ["0.0"] * len(psi) + ["0.0", "0.0 (sqrt 2 better "
                                                                                     "sensitivity: two noise samples)"]]
    table("3a. Polarisation loss in dB against the tilt psi of a linearly polarised wave from vertical",
          ["element"] + [f"{p} deg" for p in psi] + ["circular", "unpolarised (Sun, HI line)"], rows)

    # 3b axial ratio and cross-polar level against V/H amplitude and phase error
    rows = []
    for d_db, d_ph in ((0.0, 0.0), (0.1, 1.0), (0.5, 3.0), (0.7, 4.2), (1.0, 10.0), (0.0, 13.0), (3.0, 0.0)):
        a = 10 ** (d_db / 20)
        ev, eh = 1.0, a * np.exp(1j * np.radians(90 + d_ph))
        ar = axial_ratio_db(ev, eh)
        # cross-circular discrimination
        co = abs(ev - 1j * eh) / np.sqrt(2)
        cx = abs(ev + 1j * eh) / np.sqrt(2)
        xpd = 20 * np.log10(max(co, cx) / max(min(co, cx), 1e-12))
        # linear synthesis at 45 deg: tilt error and residual ellipticity
        ev2, eh2 = 1.0, a * np.exp(1j * np.radians(d_ph))
        S = stokes(np.array([ev2]), np.array([eh2]))
        tilt = 0.5 * np.degrees(np.arctan2(S[2], S[1])) - 45
        rows.append([d_db, d_ph, f"{ar:.2f}", f"{xpd:.0f}", f"{tilt:+.1f}", f"{axial_ratio_db(ev2, eh2):.0f}"])
    table("3b. Synthesised polarisation against the V/H channel imbalance left after calibration. (0.7 dB, 4.2 deg) is "
          "the difference of two channels each with rev B's 0.5 dB / 3 deg residual",
          ["amplitude dB", "phase deg", "circular: axial ratio dB", "circular: cross-hand rejection dB",
           "45 deg linear: tilt error deg", "45 deg linear: axial ratio dB (inf = pure)"], rows)

    # 3c circular polarisation from the existing 4-bit delay line (analogue, before the combiner)
    states = np.array([np.sum(TTD_BITS[[b for b in range(4) if (s >> b) & 1]]) for s in range(16)])
    dts = (states[:, None] - states[None, :]).ravel()
    fs = np.linspace(2.3e9, 6.0e9, 371)
    best_ar, best_err = [], []
    for f in fs:
        ph = (2 * np.pi * f * dts) % (2 * np.pi)
        e = np.minimum(np.abs(ph - np.pi / 2), np.abs(ph - 3 * np.pi / 2))
        i = np.argmin(e)
        best_err.append(np.degrees(e[i]))
        best_ar.append(axial_ratio_db(1.0, np.exp(1j * ph[i])))
    best_ar, best_err = np.array(best_ar), np.array(best_err)
    # a dedicated delay step: quarter-wave at band centre only
    fixed = [axial_ratio_db(1.0, np.exp(1j * 2 * np.pi * f * 62.5e-12)) for f in fs]
    RES["pol_ttd"] = dict(worst_phase_err_deg=float(best_err.max()), worst_ar_db=float(best_ar.max()),
                          median_ar_db=float(np.median(best_ar)), frac_below_3db=float(np.mean(best_ar < 3)),
                          frac_below_1db=float(np.mean(best_ar < 1)))
    rows = [[f"{f/1e9:.2f}", f"{best_err[np.argmin(np.abs(fs - f))]:.1f}", f"{best_ar[np.argmin(np.abs(fs - f))]:.2f}"]
            for f in (2.3e9, 2.44e9, 3.0e9, 3.5e9, 4.0e9, 5.0e9, 5.8e9, 6.0e9)]
    table("3c. If a V flare and an H flare were summed in one quadrant, the existing delay states (38.5 / 52.5 / 117.5 / "
          f"139 ps, 16 states each) already give quadrature: best state pair at each frequency. Over 2.3-6 GHz the "
          f"worst phase error is {best_err.max():.1f} deg (axial ratio {best_ar.max():.2f} dB), median axial ratio "
          f"{np.median(best_ar):.2f} dB; {100*np.mean(best_ar < 1):.0f} % of the band is below 1 dB. The same pairs "
          "minus/plus 90 deg give the other hand; state pairs near 0 / 180 deg give slant linear",
          ["GHz", "phase error from 90 deg", "axial ratio dB"], rows)

    # 3d a linear polarisation spinning at f_rot is two circular lines at +- f_rot
    fsr, T, frot = 200e3, 0.05, 8.0e3
    t = np.arange(int(fsr * T)) / fsr
    ev, eh = np.cos(2 * np.pi * frot * t), np.sin(2 * np.pi * frot * t)
    r = (ev - 1j * eh) / np.sqrt(2)
    l = (ev + 1j * eh) / np.sqrt(2)
    fr = np.fft.fftfreq(len(t), 1 / fsr)
    R, Lh, V = np.fft.fft(r) / len(t), np.fft.fft(l) / len(t), np.fft.fft(ev) / len(t)
    RES["pol_spin"] = dict(rhcp_line_hz=float(fr[np.argmax(np.abs(R))]), lhcp_line_hz=float(fr[np.argmax(np.abs(Lh))]),
                           v_lines_db=float(db(np.sort(np.abs(V) ** 2)[-1])))
    # with switches only (0/180 on each of V and H = 4 states per turn): square-wave commutation
    sv = np.sign(np.cos(2 * np.pi * frot * t) + 1e-12)
    sh = np.sign(np.sin(2 * np.pi * frot * t) + 1e-12)
    r4 = (sv - 1j * sh) / 2
    R4 = np.abs(np.fft.fft(r4) / len(t)) ** 2
    lines = {int(round(fr[i] / frot)): float(db(R4[i])) for i in np.argsort(R4)[-6:]}
    RES["pol_spin"]["four_state_lines_db"] = lines
    log()
    log(f"3d. A linear polarisation turned electronically at {frot/1e3:.0f} kHz (weights cos, sin on V, H): seen by a "
        f"circular antenna it is ONE line at {fr[np.argmax(np.abs(R))]/1e3:+.0f} kHz (one hand) and one at "
        f"{fr[np.argmax(np.abs(Lh))]/1e3:+.0f} kHz (the other hand), each -3.0 dB: a spinning field vector = a "
        "frequency shift of each circular component (the rotational Doppler shift, with no moving part). With only "
        "0/180 switches on V and H (four states per turn) the wanted line is "
        f"{lines.get(-1, lines.get(1, float('nan'))):.1f} dB and the strongest images are at "
        + ", ".join(f"{m:+d} f_rot: {v:.1f} dB" for m, v in sorted(lines.items()) if abs(m) != 1) + ".")

    # 3e polarisation measurement: two V and two H quadrant channels against a V-only panel
    def mc(mode, snr_db=20.0, nsamp=64, T=4000, xpol_db=-20.0, sd=0.5, sp=3.0):
        """Wave: fully polarised, random tilt psi and ellipticity chi (or pure linear). Each of four channels is a
        quadrant sum (4 elements: +6 dB SNR against the element SNR snr_db). Channel errors sd dB / sp deg each,
        plus a fixed cross-polar leak of xpol_db with random phase. Returns tilt error (deg), ellipticity error
        (deg), and recovered power (dB re the true power, combining what the mode has)."""
        psi = RNG.uniform(-90, 90, T)
        chi = RNG.uniform(-45, 45, T) if mode.endswith("any") else np.zeros(T)
        p, c = np.radians(psi), np.radians(chi)
        # Jones vector (ev, eh) of an ellipse with tilt psi from vertical and ellipticity chi
        ev = np.cos(p) * np.cos(c) - 1j * np.sin(p) * np.sin(c)
        eh = np.sin(p) * np.cos(c) + 1j * np.cos(p) * np.sin(c)
        sig = RNG.normal(size=(T, nsamp)) + 1j * RNG.normal(size=(T, nsamp))
        sig /= np.sqrt(2)
        amp = 10 ** ((snr_db + 6.0) / 20)
        x = 10 ** (xpol_db / 20)
        ch_pol = "VVHH" if mode.startswith("dual") else "VVVV"
        Y = []
        for q in ch_pol:
            g = 10 ** (RNG.normal(0, sd, T) / 20) * np.exp(1j * np.radians(RNG.normal(0, sp, T)))
            leak = x * np.exp(1j * RNG.uniform(0, 2 * np.pi, T))
            h = (ev + leak * eh) if q == "V" else (eh + leak * ev)
            n = (RNG.normal(size=(T, nsamp)) + 1j * RNG.normal(size=(T, nsamp))) / np.sqrt(2)
            Y.append(amp * (g * h)[:, None] * sig + n)
        if mode.startswith("dual"):
            yv, yh = (Y[0] + Y[1]) / 2, (Y[2] + Y[3]) / 2
            S0 = np.mean(np.abs(yv) ** 2 + np.abs(yh) ** 2, axis=1) - 2 * 0.5
            S1 = np.mean(np.abs(yv) ** 2 - np.abs(yh) ** 2, axis=1)
            S2 = 2 * np.mean(np.real(np.conj(yv) * yh), axis=1)
            S3 = 2 * np.mean(np.imag(np.conj(yv) * yh), axis=1)
            tilt = 0.5 * np.degrees(np.arctan2(S2, S1))
            ell = 0.5 * np.degrees(np.arcsin(np.clip(S3 / np.sqrt(S1 ** 2 + S2 ** 2 + S3 ** 2), -1, 1)))
            dt = (tilt - psi + 90) % 180 - 90
            # tilt is undefined near circular: only score |chi| < 30 deg
            ok = np.abs(chi) < 30
            # beam SNR if V and H are combined with the estimated polarisation (maximum ratio): signal power
            # amp^2 (|ev|^2 + |eh|^2) = amp^2 against noise 1/2 per combined pair
            p_rec = db(S0 / amp ** 2)
            return dict(tilt=[float(np.median(np.abs(dt[ok]))), float(np.percentile(np.abs(dt[ok]), 95))],
                        ell=[float(np.median(np.abs(ell - chi))), float(np.percentile(np.abs(ell - chi), 95))],
                        power=[float(np.median(p_rec)), float(np.percentile(p_rec, 5))])
        y = (Y[0] + Y[1] + Y[2] + Y[3]) / 4
        P = np.mean(np.abs(y) ** 2, axis=1) - 0.25
        p_rec = db(np.maximum(P, 1e-12) / amp ** 2)
        return dict(tilt=[float("nan")] * 2, ell=[float("nan")] * 2,
                    power=[float(np.median(p_rec)), float(np.percentile(p_rec, 5))])

    rows = []
    res = {}
    for nm, mode, kw in (("V-only panel (rev B), linear wave of unknown tilt", "single_lin", {}),
                         ("V-only panel, any polarisation", "single_any", {}),
                         ("2 V + 2 H quadrants, linear wave, 20 dB", "dual_lin", {}),
                         ("2 V + 2 H quadrants, any polarisation, 20 dB", "dual_any", {}),
                         ("same, element cross-polar -15 dB", "dual_any", dict(xpol_db=-15.0)),
                         ("same, 0 dB element SNR", "dual_any", dict(snr_db=0.0)),
                         ("same, calibrated to 0.1 dB / 1 deg, cross-polar -30 dB", "dual_any",
                          dict(sd=0.1, sp=1.0, xpol_db=-30.0))):
        o = mc(mode, **kw)
        res[nm] = o
        rows.append([nm, f"{o['tilt'][0]:.1f} / {o['tilt'][1]:.1f}", f"{o['ell'][0]:.1f} / {o['ell'][1]:.1f}",
                     f"{o['power'][0]:+.1f} / {o['power'][1]:+.1f}"])
    RES["pol_mc"] = res
    table("3e. One capture (64 samples, 20 dB per element and sample unless noted, quadrant sums, 0.5 dB / 3 deg channel "
          "residual, element cross-polar -20 dB with unknown phase). Polarisation errors: median / 95th percentile; "
          "recovered power: median / 5th percentile, dB against the wave's full power. The V-only panel cannot "
          "measure polarisation in one capture; rotation_study.py measures it by rolling the panel",
          ["configuration", "tilt error deg (|ellipticity| < 30)", "ellipticity-angle error deg",
           "recovered power dB"], rows)

    # figure
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    ax[0].plot(fs / 1e9, best_ar, label="best pair of existing delay states")
    ax[0].plot(fs / 1e9, fixed, "--", label="one fixed 62.5 ps step (quarter wave at 4 GHz)")
    ax[0].axhline(3, color="0.7", ls=":")
    ax[0].set_ylim(0, 12)
    ax[0].set_xlabel("GHz")
    ax[0].set_ylabel("axial ratio (dB)")
    ax[0].set_title("Circular polarisation from V + H flares with rev B's 4-bit delay line", fontsize=9)
    ax[0].legend(fontsize=7)
    pp = np.linspace(0, 90, 91)
    ax[1].plot(pp, -db(np.cos(np.radians(pp)) ** 2 + 1e-4), label="V-only: linear wave at tilt psi")
    ax[1].axhline(3.0, color="C1", ls="--", label="V-only: circular or unpolarised")
    ax[1].axhline(0.0, color="C2", label="V + H combined")
    ax[1].set_ylim(-1, 30)
    ax[1].set_xlabel("tilt of the wave from vertical (deg)")
    ax[1].set_ylabel("polarisation loss (dB)")
    ax[1].set_title("What one polarisation loses", fontsize=9)
    ax[1].legend(fontsize=7)
    m = np.arange(-9, 10)
    R4s = np.fft.fftshift(R4)
    frs = np.fft.fftshift(fr)
    ax[2].stem(m, [db(R4s[np.argmin(np.abs(frs - mm * frot))]) for mm in m], basefmt=" ")
    ax[2].set_ylim(-40, 2)
    ax[2].set_xlabel("harmonic of the rotation rate")
    ax[2].set_ylabel("dB re total")
    ax[2].set_title("Field spun with 0/180 switches only (4 states per turn): lines seen in one circular hand",
                    fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "polarisation.png", dpi=130)
    plt.close(fig)


# =====================================================================================================================
# 4. Pseudo-Doppler commutation = phase mode 1 in time
# =====================================================================================================================
def study_commute():
    log()
    log("## 4. Pseudo-Doppler: one receiver commutated around a ring")
    # dwell on element n: x_n = exp(j k r cos(phi0 - phi_n)). First harmonic of the commutation:
    # c1 = sum_n x_n exp(+j phi_n) = N sum_q j^(1+qN) J_(1+qN)(kr) exp(j (1+qN) phi0); bearing = arg(c1) - 90 deg.
    rows = []
    res = {}
    phi0 = np.radians(np.arange(0, 360, 1.0))
    for N, r_wl in ((4, 0.18), (4, 0.25), (4, 0.32), (8, 0.25), (8, 0.5), (8, 0.69), (8, 1.0), (8, 1.73),
                    (16, 0.5), (16, 1.0), (16, 1.73)):
        kr = 2 * np.pi * r_wl
        an = 2 * np.pi * np.arange(N) / N
        x = np.exp(1j * kr * np.cos(phi0[:, None] - an[None, :]))
        c1 = x @ np.exp(1j * an)                             # = N j J1(kr) exp(j phi0) + aliases
        est = np.angle(c1 * np.exp(-1j * np.pi / 2))
        err = np.degrees(np.angle(np.exp(1j * (est - phi0))))
        # sign of J1 flips the bearing by 180 deg beyond its first zero; fold that out (a known constant)
        if jv(1, kr) < 0:
            err = np.degrees(np.angle(np.exp(1j * (est - phi0 + np.pi))))
        amp = np.abs(c1) / N
        dev = kr                                             # peak phase deviation (rad)
        res[f"N{N}_r{r_wl}"] = dict(max_bias=float(np.max(np.abs(err))), J1=float(jv(1, kr)),
                                    amp_min=float(amp.min()))
        rows.append([N, r_wl, f"{kr:.2f}", f"{jv(1, kr):+.3f}", f"{np.max(np.abs(err)):.1f}",
                     f"{db(amp.min()**2):.1f}", f"{db(1.0 / N):.1f}"])
    RES["commute"] = res
    table("4a. Bearing from the first commutation harmonic (the classical pseudo-Doppler output), noise-free: worst "
          "bias over azimuth from the aliases J_(N-1), before any lookup correction; |c1|^2 = level of the harmonic "
          "against the carrier; last column = SNR against using all N elements at once (one element is connected at "
          "a time). r = 0.18-0.32 wl are rev B's 4-whip rings, 0.69 and 1.73 wl are the 225 mm ring at 0.915 and "
          "2.3 GHz",
          ["N", "r / wl", "k r", "J1(k r)", "max bias deg", "harmonic level dB (worst azimuth)",
           "SNR against N simultaneous, dB"], rows)
    # commutation numbers for rev B: slot >= 15 us
    for N in (4, 8, 16):
        frot = 1 / (N * 15e-6)
        log(f"   rev B minimum slot 15 us: a ring of {N} is commutated at up to {frot/1e3:.1f} kHz; at 0.915 GHz on "
            f"the 225 mm ring the equivalent rim speed is {2*np.pi*0.225*frot:.0f} m/s and the peak 'Doppler' "
            f"deviation k r f_rot = {2*np.pi*0.225/(C0/0.915e9)*frot/1e3:.0f} kHz")


def main():
    t0 = time.time()
    study_oam()
    log(f"[oam: {time.time() - t0:.0f} s]")
    t1 = time.time()
    study_superdir()
    log(f"[superdirectivity: {time.time() - t1:.0f} s]")
    t1 = time.time()
    study_pol()
    study_commute()
    log(f"[polarisation, commutation: {time.time() - t1:.0f} s]")
    (OUT / "tables.md").write_text("\n".join(_LOG) + "\n")
    (OUT / "results.json").write_text(json.dumps(RES, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
    log(f"[total {time.time() - t0:.0f} s] wrote {OUT}")


if __name__ == "__main__":
    sys.exit(main())
