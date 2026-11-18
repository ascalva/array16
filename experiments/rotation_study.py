"""Moving the array itself: what does rotating it about one, two or three axes buy, and what does it cost?

The owner's questions: "does rotation along one to many axis have any consequences?" and "i was even thinking of
using a gyroscopic movement such that we can consistently rotate it along any direction, or some complex path".
This file answers with Monte Carlo. Two arrays from rev B (REDESIGN.md):

  plate   16 whips on the 500 mm plate, two rings of 8 at 225 / 112 mm (inner ring turned 22.5 deg), 0.4-2.3 GHz
  panel   the active 4x4 of the Vivaldi tile, 25 mm pitch, 2.3-6 GHz, one linear polarisation, front hemisphere

Conventions follow experiments/evolve_layout.py (its check model): 3 deg phase / 0.5 dB calibration residual per
element (fixed to the element, so it turns with the array), 20 dB SNR per element and sample, 64 samples, Bartlett
beamscan, gross error = more than 10 deg. Section 0 cross-checks this code against evolve_layout.check().

Two estimators:

  az-only (sections a, b)   exactly evolve_layout.check(): sample covariance per orientation, Bartlett spectra of
      the K orientations added (the K dwells are not phase-coherent with each other), 1 deg then 0.05 deg search.
      Whips are isotropic and the estimator assumes the source is on the horizon.

  sphere (sections c-f)   direction on the whole sphere plus linear polarisation. Per orientation k the beam power
      b_k(d) = |s_k(d)^H y_k|^2; the element's response to a field with coherency matrix C (2x2, real part) is
      w_k = h_k^T C h_k with h_k the element's real effective-height vector in the source's (vertical, horizontal)
      basis at orientation k. w_k is linear in c = (C11, C12, C22): w_k = q_k . c, q_k = (hv^2, 2 hv hh, hh^2).
      Least squares over c gives the statistic J(d) = m^T (Q + ridge)^-1 m, m = sum_k b_k q_k, Q = sum_k q_k q_k^T,
      maximised over d (3 deg sphere grid, then three local refinements down to 0.04 deg). The polarisation angle
      is 0.5 atan2(2 C12, C11 - C22). For one orientation Q has rank 1 and J reduces to the plain Bartlett power.
      Signal model: one coherent snapshot at SNR + 10 log10(64) (evolve_layout's search model; section 0 checks it
      against the 64-sample covariance).
      Elements: a whip is a short dipole along the plate normal (response = axis . E, free space: no ground
      plane shadow); a panel flare responds to the field component along its polarisation axis times
      cos^0.75(angle from boresight) in amplitude (power cos^1.5, -4.5 dB at 60 deg: estimate), and behind the
      panel the truth is a -20 dB response with a random phase per element (estimate; the estimator's manifold has
      nothing there).

Orientation: R = Rz(alpha) Ry(beta) Rz(gamma) (degrees): alpha = azimuth axis (world vertical), beta = tilt of the
body z axis from the zenith, gamma = spin about the body z axis. Plate: body z = plate normal (home 0, 0, 0).
Panel: body z = boresight, body x = polarisation axis (home 0, 90, 0 = boresight on the horizon, V-pol; gamma
rolls the polarisation about boresight).

Path length L is in degrees of axis travel with all axes allowed to move at once at the same maximum rate, so
time = L / rate (plus dwell if the mount stops for each capture). K captures are spread evenly in time.

    OMP_NUM_THREADS=3 python experiments/rotation_study.py            # about 15-25 minutes
    OMP_NUM_THREADS=3 python experiments/rotation_study.py --quick    # smoke run
    python experiments/rotation_study.py --only a,b                   # some sections
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import evolve_layout as EL  # noqa: E402  (read-only import: ring(), check())
from array16.evolve import evolve, random_search, reselect  # noqa: E402

C0 = 299_792_458.0
OUT = HERE / "rotation"
LINES: list[str] = []


def log(s=""):
    print(s, flush=True)
    LINES.append(s)


def table(title, header, rows):
    log()
    log(title)
    log("| " + " | ".join(header) + " |")
    log("|" + "|".join(" --- " for _ in header) + "|")
    for r in rows:
        log("| " + " | ".join(str(x) for x in r) + " |")


# =========================================================================================================
# az-only model (evolve_layout.check generalised to K orientations about the vertical axis)
# =========================================================================================================
def az_mc(p, f, steps_deg, *, el_max=15.0, trials=1500, k_each=64, sigma_cal=3.0, snr_db=20.0,
          orient_err=0.0, flex=0.0, seed=12345, est2d=False):
    """p (N, 2) m; steps_deg: the K turntable angles. The first angle is the surveyed home position (exact);
    later ones carry an independent error of orient_err deg (1 sigma) unknown to the estimator. flex (deg): each
    element's phase drifts as flex * (sin(theta + psi_n) - sin(psi_n)), psi_n random: zero at home, rms flex over a
    turn. With K = 1 and the defaults the random stream is the one evolve_layout.check() draws."""
    rng = np.random.default_rng(seed)
    n = len(p)
    r, b = np.hypot(p[:, 0], p[:, 1]), np.arctan2(p[:, 1], p[:, 0])
    k = 2 * np.pi * f / C0
    steps = np.radians(np.atleast_1d(np.asarray(steps_deg, float)))
    K = len(steps)
    th = rng.uniform(0, 2 * np.pi, trials)
    el = np.radians(rng.uniform(0, el_max, trials))
    ce = np.cos(el)[:, None]
    g = 10 ** (rng.normal(0, 0.5, (trials, n)) / 20) * np.exp(1j * np.radians(rng.normal(0, sigma_cal, (trials, n))))
    Rs = []
    for i, t0 in enumerate(steps):
        s = (rng.normal(size=(trials, k_each)) + 1j * rng.normal(size=(trials, k_each))) / np.sqrt(2)
        w = ((rng.normal(size=(trials, n, k_each)) + 1j * rng.normal(size=(trials, n, k_each)))
             / np.sqrt(2) * 10 ** (-snr_db / 20))
        tt = t0 + (np.radians(rng.normal(0, orient_err, trials)) if (orient_err and i) else 0.0) * np.ones(trials)
        a = np.exp(1j * k * ce * r * np.cos(th[:, None] - b[None, :] - tt[:, None]))
        if flex:
            psi = rng.uniform(0, 2 * np.pi, (trials, n))
            a = a * np.exp(1j * np.radians(flex) * (np.sin(t0 + psi) - np.sin(psi)))
        x = (a * g)[:, :, None] * s[:, None, :] + w
        Rs.append(x @ x.conj().transpose(0, 2, 1) / k_each)

    def spec(az, cel=1.0):                     # az (..., G) -> summed Bartlett power, same shape
        pw = 0.0
        for t0, R in zip(steps, Rs):
            A = np.exp(1j * k * cel * r * np.cos(az[..., None] - b - t0))
            if A.ndim == 2:
                pw = pw + np.real(np.einsum('gi,tij,gj->tg', A.conj(), R, A))
            else:
                pw = pw + np.real(np.einsum('tgi,tij,tgj->tg', A.conj(), R, A))
        return pw

    coarse = np.radians(np.arange(0, 360, 1.0))
    if not est2d:
        c0 = coarse[np.argmax(spec(coarse), axis=1)]
        fine = c0[:, None] + np.radians(np.arange(-1.5, 1.5001, 0.05))[None, :]
        est = fine[np.arange(trials), np.argmax(spec(fine), axis=1)]
    else:                                      # the same data with a joint azimuth-elevation search (el 0-60 deg)
        els = np.radians(np.arange(0, 60.01, 4.0))
        pc = np.stack([spec(coarse, np.cos(e)) for e in els], axis=1)            # (T, E, G)
        idx = pc.reshape(trials, -1).argmax(axis=1)
        e0, c0 = els[idx // len(coarse)], coarse[idx % len(coarse)]
        for da, de in ((np.arange(-1.5, 1.501, 0.25), np.arange(-4, 4.01, 1.0)),
                       (np.arange(-0.25, 0.2501, 0.05), np.arange(-1, 1.01, 0.25))):
            A_, E_ = np.meshgrid(np.radians(da), np.radians(de), indexing="ij")
            az = c0[:, None] + A_.ravel()[None, :]
            ee = np.clip(e0[:, None] + E_.ravel()[None, :], 0, np.pi / 2)
            pw = 0.0
            for t0, R in zip(steps, Rs):
                A = np.exp(1j * k * np.cos(ee)[..., None] * r * np.cos(az[..., None] - b - t0))
                pw = pw + np.real(np.einsum('tgi,tij,tgj->tg', A.conj(), R, A))
            j = pw.argmax(axis=1)
            c0, e0 = az[np.arange(trials), j], ee[np.arange(trials), j]
        est = c0
    e = np.abs(np.degrees(np.angle(np.exp(1j * (est - th)))))
    return float(np.median(e)), float(np.percentile(e, 95)), float(np.mean(e > 10) * 100)


def fmt3(r):
    return f"{r[0]:.2f} / {r[1]:.1f} / {r[2]:.1f} %"


# =========================================================================================================
# sphere model
# =========================================================================================================
def rotm(ang_deg):
    """(..., 3) degrees (alpha, beta, gamma) -> (..., 3, 3) = Rz(alpha) Ry(beta) Rz(gamma)."""
    a, b, g = (np.radians(ang_deg[..., i]) for i in range(3))
    ca, sa, cb, sb, cg, sg = np.cos(a), np.sin(a), np.cos(b), np.sin(b), np.cos(g), np.sin(g)
    R = np.empty(ang_deg.shape[:-1] + (3, 3))
    R[..., 0, 0] = ca * cb * cg - sa * sg
    R[..., 0, 1] = -ca * cb * sg - sa * cg
    R[..., 0, 2] = ca * sb
    R[..., 1, 0] = sa * cb * cg + ca * sg
    R[..., 1, 1] = -sa * cb * sg + ca * cg
    R[..., 1, 2] = sa * sb
    R[..., 2, 0] = -sb * cg
    R[..., 2, 1] = sb * sg
    R[..., 2, 2] = cb
    return R


def dirs(az_deg, el_deg):
    az, el = np.radians(az_deg), np.radians(el_deg)
    return np.stack([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)], axis=-1)


def polbasis(d):
    """Unit vectors (vertical-ish v = d(d)/d(el), horizontal h) transverse to d (..., 3)."""
    h = np.stack([-d[..., 1], d[..., 0], np.zeros_like(d[..., 0])], axis=-1)
    nrm = np.linalg.norm(h, axis=-1, keepdims=True)
    h = np.where(nrm > 1e-9, h / np.maximum(nrm, 1e-9), np.array([0.0, 1.0, 0.0]))
    v = np.cross(d, h)
    return v, h


def fib_sphere(spacing_deg):
    n = int(round(41253.0 / spacing_deg ** 2))
    i = np.arange(n) + 0.5
    z = 1 - 2 * i / n
    phi = i * np.pi * (3 - np.sqrt(5))
    rr = np.sqrt(1 - z * z)
    return np.stack([rr * np.cos(phi), rr * np.sin(phi), z], axis=-1)


def make_plate():
    p = np.vstack([EL.ring(8, 0.225), EL.ring(8, 0.112, np.pi / 8)])
    return dict(kind="plate", pos=np.column_stack([p, np.zeros(len(p))]), axis=np.array([0.0, 0.0, 1.0]),
                home=np.array([0.0, 0.0, 0.0]))


def make_panel():
    c = (np.arange(4) - 1.5) * 0.025
    X, Y = np.meshgrid(c, c, indexing="ij")
    return dict(kind="panel", pos=np.column_stack([X.ravel(), Y.ravel(), np.zeros(16)]),
                axis=np.array([1.0, 0.0, 0.0]), home=np.array([0.0, 90.0, 0.0]))


BACK = 0.1          # panel: amplitude of the incoherent response behind the ground plane (-20 dB, estimate)
PEXP = 0.75         # panel: element amplitude pattern cos^PEXP (power cos^1.5, estimate)
RIDGE = 0.01


def heights(arr, R, D, V, H):
    """Real effective-height vectors: R (K, 3, 3) nominal, D/V/H (..., 3) -> (..., K, 2)."""
    ax = R @ arr["axis"]                                   # (K, 3)
    hv, hh = V @ ax.T, H @ ax.T                            # (..., K)
    if arr["kind"] == "panel":
        g = np.clip(D @ R[:, :, 2].T, 0, None) ** PEXP
        hv, hh = hv * g, hh * g
    return np.stack([hv, hh], axis=-1)


def qvec(h):
    return np.stack([h[..., 0] ** 2, 2 * h[..., 0] * h[..., 1], h[..., 1] ** 2], axis=-1)


def solve_J(m, Q):
    """m (..., 3), Q (..., 3, 3) -> J = m^T (Q + ridge)^-1 m and c = (Q + ridge)^-1 m."""
    tr = np.trace(Q, axis1=-2, axis2=-1)
    Qr = Q + (RIDGE * tr + 1e-12)[..., None, None] * np.eye(3)
    c = np.linalg.solve(Qr, m[..., None])[..., 0]
    return np.sum(m * c, axis=-1), c


_GRIDS: dict = {}


def grid(spacing):
    if spacing not in _GRIDS:
        D = fib_sphere(spacing)
        _GRIDS[spacing] = (D,) + polbasis(D)
    return _GRIDS[spacing]


STAGES = ((2.4, 0.6), (0.6, 0.15), (0.15, 0.0375))


def sphere_mc(arr, f, ang, az, el, psi, *, sigma_cal=3.0, amp_db=0.5, snr_db=20.0, k_samples=64, jitter=0.0,
              bias=0.0, flex=0.0, pfluct_db=0.0, seed=0, spacing=3.0, stages=STAGES, covariance=False,
              chunk=100, upper=False):
    """One source per trial at (az, el) deg with linear polarisation psi (deg from vertical); ang (K, 3) deg is the
    orientation sequence the estimator believes. Truth: each axis angle = believed + a constant offset per trial
    (bias deg, 1 sigma: home-switch / mount error) + an independent error per capture (jitter deg, 1 sigma).
    flex (deg): per-element phase drift flex * (sin(Theta_k + psi_n) - sin(psi_n)), Theta_k = summed axis travel
    from the first capture. pfluct_db: source power change between captures (dB, 1 sigma).
    upper: the estimator is told the source is above the horizon (coarse search limited to el >= -2 deg).
    Returns arrays: ang (great-circle error, deg), daz, del, pol (deg, angle between estimated and true E lines)."""
    rng = np.random.default_rng(seed)
    T, K, N = len(az), len(ang), len(arr["pos"])
    k0 = 2 * np.pi * f / C0
    d = dirs(az, el)
    v, h = polbasis(d)
    E = np.cos(np.radians(psi))[:, None] * v + np.sin(np.radians(psi))[:, None] * h
    ang_t = ang[None] + rng.normal(0, 1, (T, 1, 3)) * bias + rng.normal(0, 1, (T, K, 3)) * jitter
    Rt = rotm(ang_t)                                                           # (T, K, 3, 3)
    ph = k0 * np.einsum('tkij,nj,ti->tkn', Rt, arr["pos"], d)
    amp = np.einsum('tkij,j,ti->tk', Rt, arr["axis"], E)
    sig = np.exp(1j * ph)
    if arr["kind"] == "panel":
        c = np.einsum('tki,ti->tk', Rt[..., 2], d)
        front = c > 0
        coh = amp * np.clip(c, 0, None) ** PEXP
        back = BACK * amp[..., None] * np.exp(1j * rng.uniform(0, 2 * np.pi, (T, K, N)))
        sig = np.where(front[..., None], coh[..., None] * sig, back)
    else:
        rng.uniform(0, 1, (T, K, N))                                           # keep streams aligned between arrays
        sig = amp[..., None] * sig
    cal = 10 ** (rng.normal(0, amp_db, (T, N)) / 20) * np.exp(1j * np.radians(rng.normal(0, sigma_cal, (T, N))))
    sig = sig * cal[:, None, :]
    psi_n = rng.uniform(0, 2 * np.pi, (T, 1, N))
    if flex:
        Th = np.radians(np.sum(ang - ang[0], axis=1))[None, :, None]
        sig = sig * np.exp(1j * np.radians(flex) * (np.sin(Th + psi_n) - np.sin(psi_n)))
    pf = rng.normal(0, 1, (T, K))
    if pfluct_db:
        sig = sig * 10 ** (pf * pfluct_db / 20)[..., None]
    if covariance:        # 64 samples and a sample covariance per capture, reduced to its principal component
        s = (rng.normal(size=(T, K, 1, k_samples)) + 1j * rng.normal(size=(T, K, 1, k_samples))) / np.sqrt(2)
        w = ((rng.normal(size=(T, K, N, k_samples)) + 1j * rng.normal(size=(T, K, N, k_samples)))
             / np.sqrt(2) * 10 ** (-snr_db / 20))
        x = sig[..., None] * s + w
        Rc = x @ x.conj().swapaxes(-1, -2) / k_samples
    else:
        se = 10 ** (-(snr_db + 10 * np.log10(k_samples)) / 20)
        y = sig + (rng.normal(size=(T, K, N)) + 1j * rng.normal(size=(T, K, N))) / np.sqrt(2) * se
        y = y.astype(np.complex64)

    def beam(S, tsl):                       # S (..., K, N) conj-steering against the data of trials tsl
        if covariance:
            return np.real(np.einsum('...kn,tknm,...km->t...k' if S.ndim == 3 else 'tfkn,tknm,tfkm->tfk',
                                     S, Rc[tsl], S.conj()))
        return np.abs(np.einsum('gkn,tkn->tgk' if S.ndim == 3 else 'tfkn,tkn->tfk', S, y[tsl])) ** 2

    Rn = rotm(ang)
    pw = np.einsum('kij,nj->kni', Rn, arr["pos"]).astype(np.float32)           # (K, N, 3)
    D, V, H = grid(spacing)
    q = qvec(heights(arr, Rn, D, V, H))                                        # (G, K, 3)
    Qg = np.einsum('gkc,gkd->gcd', q, q)
    tr = np.trace(Qg, axis1=-2, axis2=-1)
    Qi = np.linalg.inv(Qg + (RIDGE * tr + 1e-12)[:, None, None] * np.eye(3))
    S = np.exp(-1j * k0 * np.einsum('gc,knc->gkn', D.astype(np.float32), pw)).astype(np.complex64)
    out = dict(ang=np.empty(T), daz=np.empty(T), del_=np.empty(T), pol=np.empty(T), est_el=np.empty(T))
    for t0 in range(0, T, chunk):
        ts = slice(t0, min(T, t0 + chunk))
        b = beam(S, ts)                                                        # (t, G, K)
        m = np.einsum('tgk,gkc->tgc', b, q)
        J = np.einsum('tgc,gcd,tgd->tg', m, Qi, m)
        if upper:
            J = np.where(D[:, 2] >= -0.035, J, -np.inf)
        d0 = D[np.argmax(J, axis=1)]                                           # (t, 3)
        cbest = None
        for half, step in stages:
            o = np.radians(np.arange(-half, half + 1e-9, step))
            X, Y = np.meshgrid(o, o, indexing="ij")
            t1, t2 = polbasis(d0)
            Df = d0[:, None, :] + X.ravel()[None, :, None] * t1[:, None, :] + Y.ravel()[None, :, None] * t2[:, None, :]
            Df /= np.linalg.norm(Df, axis=-1, keepdims=True)
            Vf, Hf = polbasis(Df)
            Sf = np.exp(-1j * k0 * np.einsum('tfc,knc->tfkn', Df.astype(np.float32), pw)).astype(np.complex64)
            bf = beam(Sf, ts)                                                  # (t, F, K)
            qf = qvec(heights(arr, Rn, Df, Vf, Hf))                            # (t, F, K, 3)
            mf = np.einsum('tfk,tfkc->tfc', bf, qf)
            Jf, cf = solve_J(mf, np.einsum('tfkc,tfkd->tfcd', qf, qf))
            j = np.argmax(Jf, axis=1)
            ar = np.arange(len(j))
            d0, cbest, vb, hb = Df[ar, j], cf[ar, j], Vf[ar, j], Hf[ar, j]
        if cbest is None:                                                      # no refinement stage: coarse answer
            g_i = np.argmax(J, axis=1)
            ar = np.arange(len(g_i))
            cbest = np.einsum('tcd,td->tc', Qi[g_i], m[ar, g_i])
            vb, hb = V[g_i], H[g_i]
        pe = 0.5 * np.arctan2(2 * cbest[:, 1], cbest[:, 0] - cbest[:, 2])
        Eh = np.cos(pe)[:, None] * vb + np.sin(pe)[:, None] * hb
        out["pol"][ts] = np.degrees(np.arccos(np.clip(np.abs(np.sum(Eh * E[ts], axis=1)), 0, 1)))
        out["ang"][ts] = np.degrees(np.arccos(np.clip(np.sum(d0 * d[ts], axis=1), -1, 1)))
        eaz = np.degrees(np.arctan2(d0[:, 1], d0[:, 0]))
        eel = np.degrees(np.arcsin(np.clip(d0[:, 2], -1, 1)))
        out["daz"][ts] = np.abs((eaz - az[ts] + 180) % 360 - 180)
        out["del_"][ts] = np.abs(eel - el[ts])
        out["est_el"][ts] = eel
    return out


def summ(o, key="ang"):
    e = o[key]
    return float(np.median(e)), float(np.percentile(e, 95)), float(np.mean(o["ang"] > 10) * 100)


def sphere_sources(T, seed):
    rng = np.random.default_rng(seed)
    az = rng.uniform(0, 360, T)
    el = np.degrees(np.arcsin(rng.uniform(-1, 1, T)))
    return az, el, rng.uniform(0, 180, T)


# ---- orientation paths -------------------------------------------------------------------------------------
def tri(x, lo, hi, start):
    span = hi - lo
    u = np.mod(x + (start - lo), 2 * span)
    return lo + np.where(u <= span, u, 2 * span - u)


def sample_wp(wp, K, L):
    """Waypoints (M, 3) deg joined by straight moves, every axis at up to the same rate: a segment takes
    max |delta| 'degrees of time'. Returns the K orientations at times 0 .. L (holding the last waypoint)."""
    wp = np.asarray(wp, float)
    seg = np.max(np.abs(np.diff(wp, axis=0)), axis=1) if len(wp) > 1 else np.zeros(0)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    cum = cum + np.arange(len(cum)) * 1e-9
    t = np.linspace(0, L, K)
    return np.column_stack([np.interp(t, cum, wp[:, i]) for i in range(3)])


PANEL_ROWS = {1: [0.0], 2: [-35.0, 35.0], 3: [-50.0, 0.0, 50.0], 4: [-60.0, -20.0, 20.0, 60.0]}
PLATE_ROWS = {1: [35.0], 2: [20.0, 50.0], 3: [15.0, 40.0, 65.0], 4: [10.0, 30.0, 50.0, 70.0]}


def path(kind, name, L, K):
    """(K, 3) degrees. kind 'plate' | 'panel'. See the module docstring for the axes."""
    t = np.linspace(0, L, K)
    z = np.zeros(K)
    panel = kind == "panel"
    if name == "fixed":
        return np.tile([0.0, 90.0 if panel else 0.0, 0.0], (K, 1))
    if name == "az1":
        return np.column_stack([t, z + (90 if panel else 0), z])
    if name == "roll1":                       # panel: about boresight; plate: tilt only (about a horizontal axis)
        return np.column_stack([z, z + 90, t]) if panel else np.column_stack([z, tri(t, 0, 90, 0), z])
    if name == "raster":
        n = int(max(1, min(4, round(L / 360))))
        rows = (PANEL_ROWS if panel else PLATE_ROWS)[n]
        beta = [90 - e for e in rows] if panel else rows
        A = max((L - (abs(beta[-1] - beta[0]))) / n, 0.0)
        wp, a = [], 0.0
        for i, bb in enumerate(beta):
            wp.append([a, bb, 0.0])
            a += A if i % 2 == 0 else -A
            wp.append([a, bb, 0.0])
        return sample_wp(wp, K, L)
    if name == "spiral":
        if panel:
            e0 = min(70.0, L / 2)
            return np.column_stack([t, 90 - (-e0 + 2 * e0 * t / max(L, 1e-9)), z])
        return np.column_stack([t, min(90.0, L) * t / max(L, 1e-9), z])
    if name in ("tumble2", "tumble3"):
        b = 90 - tri(0.618 * t, -75, 75, 0) if panel else tri(0.618 * t, 0, 90, 0)
        return np.column_stack([t, b, 0.382 * t if name == "tumble3" else z])
    raise ValueError(name)


N_WP = 6


def opt_bounds(kind, L):
    lo = np.tile([-min(L, 360.0), 0.0, -min(L, 180.0)], N_WP)
    hi = np.tile([min(L, 360.0), 180.0 if kind == "panel" else 90.0, min(L, 180.0)], N_WP)
    return lo, hi


def opt_path(x, K, L):
    return sample_wp(np.asarray(x).reshape(N_WP, 3), K, L)


def make_path_score(arr, f, L, *, K=16, trials=64, seed=0):
    def score(X, gen):
        az, el, psi = sphere_sources(trials, [seed, gen, 1])
        out = np.empty(len(X))
        for i, x in enumerate(X):
            o = sphere_mc(arr, f, opt_path(x, K, L), az, el, psi, seed=[seed, gen, 2], spacing=5.0,
                          stages=((2.5, 0.625),), chunk=trials)
            out[i] = np.mean(np.minimum(o["ang"], 30.0)) + 0.2 * np.mean(np.minimum(o["pol"], 45.0))
        return out
    return score


# =========================================================================================================
# sections
# =========================================================================================================
def sec0(res, quick):
    log("## 0. Cross-checks")
    T = 300 if quick else 1500
    rows = []
    lam = C0 / 1e9
    cases = [("4 whips, r 0.25 wl", EL.ring(4, 0.25 * lam), 1e9, 15.0, 12345, "0.93 / 2.62 / 0.0 %"),
             ("4 whips, r 0.354 wl", EL.ring(4, 0.354 * lam), 1e9, 15.0, 12345, "0.75 / 170.7 / 10.6 %"),
             ("8 whips, r 0.5 wl", EL.ring(8, 0.5 * lam), 1e9, 15.0, 12345, "0.33 / 0.93 / 0.0 %"),
             ("8 whips, r 225 mm, 2.3 GHz, 0-30 deg up", EL.ring(8, 0.225), 2.3e9, 30.0, 4242, "(30 % gross worst in band)")]
    for name, p, f, elm, sd, pub in cases:
        ref = EL.check(p, [f], el_max=elm, trials=T, seed=sd)[0]
        mine = az_mc(p, f, [0.0], el_max=elm, trials=T, seed=sd)
        rows.append([name, pub, fmt3(ref), fmt3(mine)])
        res.setdefault("crosscheck", {})[name] = dict(evolve_layout=list(map(float, ref)), this=mine)
    table("az-only model against evolve_layout.check() (median / 95th / gross)",
          ["case", "evolve_study.md", "evolve_layout.check now", "rotation_study.az_mc, K = 1"], rows)
    # sphere model: effective snapshot against 64-sample covariance; panel on boresight against REDESIGN section 9
    T2 = 100 if quick else 300
    rows = []
    panel, plate = make_panel(), make_plate()
    rng = np.random.default_rng(5)
    az, el, psi = rng.uniform(-45, 45, T2), rng.uniform(-30, 30, T2), np.zeros(T2)
    for f in (2.44e9, 3.5e9, 5.8e9):
        a = sphere_mc(panel, f, panel["home"][None], az, el, psi, seed=1)
        c = sphere_mc(panel, f, panel["home"][None], az, el, psi, seed=1, covariance=True, chunk=50)
        rows.append([f"panel {f / 1e9:.2f} GHz, source within 45 deg az / 30 deg el of boresight",
                     f"{np.median(a['ang']):.2f} / {np.percentile(a['ang'], 95):.2f}",
                     f"{np.median(c['ang']):.2f} / {np.percentile(c['ang'], 95):.2f}"])
    az, el = rng.uniform(0, 360, T2), rng.uniform(20, 70, T2)
    a = sphere_mc(plate, 1.575e9, plate["home"][None], az, el, psi, seed=1)
    c = sphere_mc(plate, 1.575e9, plate["home"][None], az, el, psi, seed=1, covariance=True, chunk=50)
    rows.append(["plate 1.575 GHz, sources 20-70 deg up", f"{np.median(a['ang']):.2f} / {np.percentile(a['ang'], 95):.2f}",
                 f"{np.median(c['ang']):.2f} / {np.percentile(c['ang'], 95):.2f}"])
    table("sphere model, one orientation: great-circle error median / 95th (deg); REDESIGN section 9 quotes "
          "0.12/0.49, 0.12/0.49, 0.10/0.46 for the panel", ["case", "one snapshot at SNR + 18 dB", "64 samples, covariance"], rows)
    res["crosscheck"]["sphere"] = rows


def scen_a():
    lam = C0 / 1e9
    c4 = (np.arange(4) - 1.5) * 0.164
    gx, gy = np.meshgrid(c4, c4, indexing="ij")
    return [
        ("ring of 8, 225 mm, 2.3 GHz, sources 0-30 deg up", EL.ring(8, 0.225), 2.3e9, 30.0, 8),
        ("ring of 8, 225 mm, 1.8 GHz, sources 0-30 deg up", EL.ring(8, 0.225), 1.8e9, 30.0, 8),
        ("4 whips, r 0.354 wl, sources 0-15 deg up", EL.ring(4, 0.354 * lam), 1e9, 15.0, 4),
        ("4 whips, r 0.45 wl", EL.ring(4, 0.45 * lam), 1e9, 15.0, 4),
        ("4 whips, r 0.60 wl", EL.ring(4, 0.60 * lam), 1e9, 15.0, 4),
        ("4 whips, r 1.0 wl", EL.ring(4, 1.0 * lam), 1e9, 15.0, 4),
        ("4x4 grid, 164 mm pitch at 2.3 GHz (1.26 wl)", np.column_stack([gx.ravel(), gy.ravel()]), 2.3e9, 15.0, 4),
        ("ring of 8, 225 mm at 5.8 GHz (spacing 3.3 wl)", EL.ring(8, 0.225), 5.8e9, 15.0, 8),
    ]


def sector_steps(K, sym):
    return np.arange(K) * (360.0 / sym) / K


def sec_a(res, quick):
    log()
    log("## a. Fragile / ambiguous layouts turned about the vertical axis, K orientations combined")
    T = 300 if quick else 1500
    Ks = (1, 2, 4, 8)
    out = {}
    rows_e, rows_t, rows_n = [], [], []
    for name, p, f, elm, sym in scen_a():
        re, rt, rn = [name], [name], [name]
        for K in Ks:
            st = sector_steps(K, sym)
            a = az_mc(p, f, st, el_max=elm, trials=T, k_each=64, seed=101)
            b = az_mc(p, f, st, el_max=elm, trials=T, k_each=64 // K, seed=101)
            c = az_mc(p, f, np.arange(K) * 360.0 / K, el_max=elm, trials=T, k_each=64 // K, seed=101)
            out.setdefault(name, {})[K] = dict(each64=a, total64=b, naive_total64=c)
            re.append(fmt3(a)); rt.append(fmt3(b)); rn.append(fmt3(c))
        rows_e.append(re); rows_t.append(rt); rows_n.append(rn)
    hdr = ["layout"] + [f"K = {K}" for K in Ks]
    table("a1. Equal time: 64 samples in total, split over K orientations that fill one symmetry sector "
          "(step = 360 / (symmetry x K): 45/K deg for a ring of 8, 90/K deg for 4 whips or a square grid). "
          "median / 95th percentile (deg) / gross", hdr, rows_t)
    table("a2. K times the time: 64 samples at every orientation, same steps", hdr, rows_e)
    table("a3. Equal time, naive steps of 360 / K deg (a symmetric layout lands on itself)", hdr, rows_n)
    # is rotation even needed for the elevated-source failure? the same single capture with a 2-D search
    rows = []
    for name, p, f, elm, sym in scen_a()[:2]:
        a = az_mc(p, f, [0.0], el_max=elm, trials=T, seed=101)
        b = az_mc(p, f, [0.0], el_max=elm, trials=T, seed=101, est2d=True)
        c = az_mc(p, f, sector_steps(4, sym), el_max=elm, trials=T, k_each=16, seed=101, est2d=True)
        rows.append([name, fmt3(a), fmt3(b), fmt3(c)])
        out[name]["est2d"] = dict(az_only=a, azel=b, azel_K4=c)
    table("a4. The elevated-source failure of the 225 mm ring is an estimator fault: the same data searched over "
          "azimuth and elevation (0-60 deg)", ["layout", "az-only, K = 1", "az-el search, K = 1", "az-el search, K = 4"], rows)
    res["a"] = out


def sec_b(res, quick):
    log()
    log("## b. Orientation error and cable-flex phase drift")
    T = 300 if quick else 1500
    sc = scen_a()
    pick = [sc[0], sc[3], sc[6]]
    conds = [("perfect", {}), ("angle 0.1 deg", dict(orient_err=0.1)), ("angle 0.5 deg", dict(orient_err=0.5)),
             ("angle 1 deg", dict(orient_err=1.0)), ("flex 1 deg", dict(flex=1.0)), ("flex 3 deg", dict(flex=3.0)),
             ("flex 10 deg", dict(flex=10.0)), ("angle 0.5 deg + flex 3 deg", dict(orient_err=0.5, flex=3.0))]
    out = {}
    for name, p, f, elm, sym in pick:
        rows = []
        for cn, kw in conds:
            r = [cn]
            for K in (1, 2, 4, 8):
                a = az_mc(p, f, sector_steps(K, sym), el_max=elm, trials=T, k_each=64 // K, seed=202, **kw)
                out.setdefault(name, {}).setdefault(cn, {})[K] = a
                r.append(fmt3(a))
            rows.append(r)
        table(f"b. {name}: 64 samples in total; median / 95th (deg) / gross. The first orientation is the surveyed "
              "home (exact, no drift); errors apply to the moves", ["condition", "K = 1", "K = 2", "K = 4", "K = 8"], rows)
    # the same with a well-behaved layout, where there is nothing to fix: does motion hurt?
    p = np.vstack([EL.ring(8, 0.225), EL.ring(8, 0.112, np.pi / 8)])
    rows = []
    for cn, kw in conds:
        r = [cn]
        for K in (1, 2, 4, 8):
            a = az_mc(p, 1.575e9, sector_steps(K, 8), el_max=15.0, trials=T, k_each=64 // K, seed=203, **kw)
            out.setdefault("two rings of 8", {}).setdefault(cn, {})[K] = a
            r.append(fmt3(a))
        rows.append(r)
    table("b. two rings of 8 (225 / 112 mm), 1.575 GHz, sources 0-15 deg up: a layout with nothing to fix",
          ["condition", "K = 1", "K = 2", "K = 4", "K = 8"], rows)
    res["b"] = out


GAIN0 = {2.44e9: 6.4, 3.5e9: 10.7, 5.8e9: 15.5}      # REDESIGN section 4, realised gain on boresight (dBi)


def sec_c(res, quick):
    log()
    log("## c. The panel on one axis (azimuth)")
    T = 100 if quick else 300
    panel = make_panel()
    out = {}
    offs = (0, 15, 30, 45, 60, 75, 90, 120, 180)
    rows = []
    for f in (2.44e9, 3.5e9, 5.8e9):
        r1, r2 = [f"{f / 1e9:.2f} GHz: az error median / 95th (deg) / gross"], [f"{f / 1e9:.2f} GHz: gain (dBi)"]
        for o in offs:
            rng = np.random.default_rng([3, o])
            az, el = np.full(T, float(o)) + rng.uniform(-2, 2, T), rng.uniform(0, 15, T)
            s = summ(sphere_mc(panel, f, panel["home"][None], az, el, np.zeros(T), seed=[31, o]), "daz")
            c = np.cos(np.radians(o))
            gain = GAIN0[f] + 15 * np.log10(c) if c > 1e-3 else None
            out.setdefault(f"{f / 1e9:.2f}", {})[o] = dict(daz=s, gain_dbi=gain)
            r1.append(fmt3(s))
            r2.append(f"{gain:.1f}" if gain is not None and gain > GAIN0[f] - 20 else f"<= {GAIN0[f] - 20:.1f}")
        rows += [r1, r2]
    table("c1. Fixed panel: azimuth DF and gain against the source's azimuth off boresight (source 0-15 deg up, "
          "V-pol, one capture; gain = REDESIGN boresight gain + 15 log10 cos, a model, with a -20 dB floor behind)",
          ["", *[f"{o} deg" for o in offs]], rows)
    rows = []
    for f in (2.44e9, 5.8e9):
        for K in (1, 2, 3, 4, 6, 8):
            rng = np.random.default_rng([4, K])
            az, el = rng.uniform(0, 360, T * 2), rng.uniform(0, 15, T * 2)
            ang = np.column_stack([np.arange(K) * 360.0 / K, np.full(K, 90.0), np.zeros(K)])
            o = sphere_mc(panel, f, ang, az, el, np.zeros(2 * T), seed=[32, K])
            ok = float(np.mean(o["daz"] < 1.0) * 100)
            worst = GAIN0[f] + 15 * np.log10(np.cos(np.radians(min(180.0 / K, 89.9)))) if K > 2 else None
            out.setdefault(f"cover_{f / 1e9:.2f}", {})[K] = dict(daz=summ(o, "daz"), within1deg_pct=ok, worst_gain=worst)
            rows.append([f"{f / 1e9:.2f} GHz", K, f"{360 / K:.0f} deg", fmt3(summ(o, "daz")), f"{ok:.0f} %",
                         f"{worst:.1f} dBi" if worst is not None else "none (gaps)"])
    table("c2. Panel stepped to K azimuths (one capture each), sources anywhere in azimuth, 0-15 deg up",
          ["band", "K", "step", "az error median / 95th / gross", "azimuths found within 1 deg", "worst-azimuth gain"], rows)
    rows = []
    for rate in (6.0, 30.0, 90.0):
        for K in (3, 4, 6, 8):
            r = [f"{rate:.0f} deg/s", K]
            for dwell in (0.05, 0.5, 5.0):
                r.append(f"{K * dwell + (K - 1) * (360.0 / K) / rate:.1f} s")
            r.append(f"{360.0 / rate:.0f} s")
            rows.append(r)
    table("c3. Time for one 360 deg look (K stops, moving between them at the given rate; acceleration ignored); "
          "last column: one continuous turn capturing on the fly", ["slew rate", "K stops", "dwell 0.05 s", "0.5 s", "5 s", "continuous"], rows)
    res["c"] = out


def sec_d(res, quick):
    log()
    log("## d. Two axes: elevation near the horizon")
    T = 100 if quick else 300
    plate, panel = make_plate(), make_panel()
    out = {}
    cfgs = [("flat, fixed", [(0, 0, 0)] * 4), ("flat, 4 azimuths (one axis)", [(k * 95.6, 0, 0) for k in range(4)]),
            ("tilted 30 deg, fixed", [(0, 30, 0)] * 4), ("tilted 30 deg, 4 azimuths", [(k * 90, 30, 0) for k in range(4)]),
            ("tilted 60 deg, 4 azimuths", [(k * 90, 60, 0) for k in range(4)]),
            ("on edge (90 deg), 4 azimuths", [(k * 90, 90, 0) for k in range(4)]),
            ("tilts 0/30/60/90 at 4 azimuths", [(k * 90, 30 * k, 0) for k in range(4)])]
    els = (0, 2, 5, 10, 20, 30, 45, 70)
    for f in (1.575e9, 0.433e9):
        rows = []
        for cn, ang in cfgs:
            r = [cn]
            for e in els:
                rng = np.random.default_rng([5, e])
                az = rng.uniform(0, 360, T)
                o = sphere_mc(plate, f, np.array(ang, float), az, np.full(T, float(e)), np.zeros(T), seed=[41, e])
                fold = np.abs(np.abs(o["est_el"]) - e)
                wrong = float(np.mean(o["est_el"] < -max(e, 1) / 2) * 100)
                out.setdefault(f"plate_{f / 1e9:.3f}", {}).setdefault(cn, {})[e] = dict(
                    el=float(np.median(o["del_"])), folded=float(np.median(fold)), az=float(np.median(o["daz"])),
                    below_pct=wrong, gross=float(np.mean(o["ang"] > 10) * 100))
                r.append(f"{np.median(fold):.2f} ({wrong:.0f} %) / {np.median(o['daz']):.2f}")
            rows.append(r)
        table(f"d1. Whip plate (two rings of 8), {f / 1e9:.3f} GHz, V-pol source at the given elevation, 4 captures: "
              "median elevation error with the sign folded (share of answers below the horizon) / median azimuth "
              "error, deg", ["plate", *[f"el {e}" for e in els]], rows)
    cfgs = [("boresight on the horizon, fixed", [(0, 90, 0)] * 4), ("boresight 45 deg up, fixed", [(0, 45, 0)] * 4),
            ("elevation axis: 0 / 30 / 60 / 90 deg", [(0, 90 - 30 * k, 0) for k in range(4)])]
    els = (0, 15, 30, 45, 60, 75, 90)
    for f in (2.44e9, 5.8e9):
        rows = []
        for cn, ang in cfgs:
            r = [cn]
            for e in els:
                rng = np.random.default_rng([6, e])
                az = rng.uniform(-20, 20, T)
                o = sphere_mc(panel, f, np.array(ang, float), az, np.full(T, float(e)), np.zeros(T), seed=[42, e])
                out.setdefault(f"panel_{f / 1e9:.2f}", {}).setdefault(cn, {})[e] = dict(
                    el=float(np.median(o["del_"])), ang=summ(o))
                r.append(f"{np.median(o['del_']):.2f} / {np.mean(o['ang'] > 10) * 100:.0f} %")
            rows.append(r)
        table(f"d2. Panel, {f / 1e9:.2f} GHz, source within 20 deg of the boresight azimuth, 4 captures: median "
              "elevation error (deg) / gross", ["panel", *[f"el {e}" for e in els]], rows)
    res["d"] = out


def sec_e(res, quick):
    log()
    log("## e. Panel polarisation")
    T = 100 if quick else 300
    panel = make_panel()
    out = {}
    rows = []
    psis = (0, 30, 45, 60, 75, 85, 90)
    rng = np.random.default_rng(7)
    az, el = rng.uniform(-20, 20, T), rng.uniform(0, 15, T)
    for f in (2.44e9, 5.8e9):
        r = [f"{f / 1e9:.2f} GHz: great-circle error median / 95th (deg) / gross"]
        for ps in psis:
            s = summ(sphere_mc(panel, f, panel["home"][None], az, el, np.full(T, float(ps)), seed=[51, ps]))
            out.setdefault(f"loss_{f / 1e9:.2f}", {})[ps] = s
            r.append(fmt3(s))
        rows.append(r)
    rows.insert(0, ["polarisation loss (dB) = 20 log10 cos(angle)"] +
                [f"{-20 * np.log10(max(np.cos(np.radians(p)), 1e-9)):.1f}" if p < 90 else "total" for p in psis])
    table("e1. Fixed V-pol panel, linearly polarised source turned from vertical (ideal flare: no cross-polar "
          "response; one capture, source near boresight)", ["", *[f"{p} deg" for p in psis]], rows)
    rows = []
    for f in (5.8e9,):
        for snr in (20.0, 0.0):
            for K in (1, 2, 3, 4, 8):
                for pfl in ((0.0, 1.0, 3.0) if K in (4, 8) and snr == 20.0 else (0.0,)):
                    ang = np.column_stack([np.zeros(K), np.full(K, 90.0), np.arange(K) * 180.0 / K])
                    psi = np.random.default_rng(8).uniform(0, 180, T)
                    o = sphere_mc(panel, f, ang, az, el, psi, snr_db=snr, pfluct_db=pfl, seed=[52, K])
                    out.setdefault("roll", {})[f"{snr:.0f}dB_K{K}_pf{pfl}"] = dict(pol=summ(o, "pol"), ang=summ(o))
                    rows.append([f"{snr:.0f} dB", K, f"{180 / K:.0f} deg" if K > 1 else "-", f"{pfl:.0f} dB",
                                 f"{np.median(o['pol']):.1f} / {np.percentile(o['pol'], 95):.1f}", fmt3(summ(o))])
    table("e2. Panel rolled about boresight to K angles over 180 deg, 5.8 GHz, source of random linear "
          "polarisation near boresight: polarisation-angle error and direction error",
          ["SNR", "K", "roll step", "source power change between captures (1 sigma)",
           "polarisation angle error median / 95th (deg)", "direction median / 95th / gross"], rows)
    res["e"] = out


def sec_f(res, quick):
    log()
    log("## f. Arbitrary orientation paths, sources anywhere on the sphere, random linear polarisation")
    T = 100 if quick else 300
    K = 24
    Ls = (90, 180, 360, 720, 1440) if not quick else (180, 720)
    arrs = [("plate", make_plate(), 1.575e9), ("panel", make_panel(), 5.8e9), ("panel", make_panel(), 2.44e9)]
    names = ["fixed", "az1", "roll1", "raster", "spiral", "tumble2", "tumble3"]
    label = dict(fixed="fixed", az1="one axis: azimuth", roll1="one axis: {}", raster="two axes: az-el raster",
                 spiral="two axes: spiral", tumble2="two axes: tumble", tumble3="three axes: tumble (gyroscopic)",
                 evolved="three axes: evolved path", random="three axes: best random path")
    out = {}
    az, el, psi = sphere_sources(T, 99)
    opt_L = (180, 720) if not quick else (180,)
    for kind, arr, f in arrs:
        key = f"{kind}_{f / 1e9:.3f}"
        lab = dict(label)
        lab["roll1"] = "one axis: roll about boresight" if kind == "panel" else "one axis: tilt 0-90 deg"
        # (iii) paths designed with the repo's optimiser
        designs = {}
        if (kind, f) != ("panel", 2.44e9):
            for L in opt_L:
                t0 = time.time()
                score = make_path_score(arr, f, L, seed=7)
                lo, hi = opt_bounds(kind, L)
                init = path(kind, "tumble3", L, N_WP).reshape(-1)
                gens, pop = (4, 8) if quick else (30, 16)
                ga = evolve(score, lo, hi, pop=pop, generations=gens, elite=2, sigma=0.15, crossover=0.5, seed=L, init=init)
                rs = random_search(score, lo, hi, batches=ga.evaluations // pop, batch=pop, seed=L)
                for nm, r in (("evolved", ga), ("random", rs)):
                    X, F = reselect(score, r.finalists[:6], draws=3 if quick else 8)
                    designs[(nm, L)] = X[0]
                    out.setdefault(key, {}).setdefault("designs", {})[f"{nm}_L{L}"] = dict(
                        waypoints_deg=np.round(X[0].reshape(N_WP, 3), 1).tolist(), search_score=float(F[0]))
                Xs, Fs = reselect(score, [init], draws=3 if quick else 8)
                out[key]["designs"][f"tumble3_L{L}"] = dict(search_score=float(Fs[0]))
                log(f"  {key} L = {L}: search scores (lower is better) evolved {out[key]['designs'][f'evolved_L{L}']['search_score']:.2f}, "
                    f"random {out[key]['designs'][f'random_L{L}']['search_score']:.2f}, tumble3 {Fs[0]:.2f}  "
                    f"({ga.evaluations} evaluations each, {time.time() - t0:.0f} s)")
        rows = []
        for nm in names + ["evolved", "random"]:
            r = [lab[nm]]
            for L in Ls:
                if nm in ("evolved", "random"):
                    if (nm, L) not in designs:
                        r.append("")
                        continue
                    ang = opt_path(designs[(nm, L)], K, L)
                else:
                    ang = path(kind, nm, L, K)
                o = sphere_mc(arr, f, ang, az, el, psi, seed=[61, L])
                s = summ(o)
                out.setdefault(key, {}).setdefault(nm, {})[L] = dict(ang=s, pol=summ(o, "pol"), el=float(np.median(o["del_"])))
                r.append(f"{s[0]:.2f} / {s[1]:.1f} / {s[2]:.0f} % / {np.median(o['pol']):.1f}")
            rows.append(r)
        table(f"f1. {kind}, {f / 1e9:.3f} GHz, perfect orientation knowledge, {K} captures spread over the path: "
              "great-circle error median / 95th (deg) / gross / median polarisation-angle error (deg); "
              "45 deg = polarisation not recovered", ["path \\ length L (time at 30 deg/s)",
                                                       *[f"{L} deg ({L / 30:.0f} s)" for L in Ls]], rows)
        # imperfect orientation knowledge at one length
        L = 720
        conds = [("perfect", {}), ("jitter 0.2 deg", dict(jitter=0.2)), ("jitter 1 deg", dict(jitter=1.0)),
                 ("axis zero offset 0.5 deg", dict(bias=0.5)), ("flex 3 deg", dict(flex=3.0)),
                 ("jitter 1 + offset 0.5 + flex 3", dict(jitter=1.0, bias=0.5, flex=3.0))]
        rows = []
        for nm in names + ["evolved"]:
            if nm == "evolved" and (nm, L) not in designs:
                continue
            ang = opt_path(designs[(nm, L)], K, L) if nm == "evolved" else path(kind, nm, L, K)
            r = [lab[nm]]
            for cn, kw in conds:
                o = sphere_mc(arr, f, ang, az, el, psi, seed=[62, L], **kw)
                s = summ(o)
                out[key].setdefault("imperfect", {}).setdefault(nm, {})[cn] = dict(ang=s, pol=summ(o, "pol"))
                r.append(f"{s[0]:.2f} / {s[1]:.1f} / {s[2]:.0f} % / {np.median(o['pol']):.1f}")
            rows.append(r)
        table(f"f2. {kind}, {f / 1e9:.3f} GHz, L = {L} deg: imperfect orientation knowledge (per-axis angle errors "
              "unknown to the estimator; flex = per-element phase drift with travel). Same columns as f1",
              ["path", *[c for c, _ in conds]], rows)
        # sources restricted to where a fixed array already works: does motion still help there?
        rng = np.random.default_rng(77)
        if kind == "plate":
            az2, el2 = rng.uniform(0, 360, T), rng.uniform(15, 60, T)
            where = "sources 15-60 deg above the plate, any azimuth"
        else:
            az2, el2 = rng.uniform(-45, 45, T), rng.uniform(-30, 30, T)
            where = "sources within 45 deg az / 30 deg el of the home boresight"
        rows = []
        for nm in names:
            r = [lab[nm]]
            for L in (180, 720):
                o = sphere_mc(arr, f, path(kind, nm, L, K), az2, el2, psi, seed=[63, L])
                s = summ(o)
                out[key].setdefault("home_sector", {}).setdefault(nm, {})[L] = dict(ang=s, pol=summ(o, "pol"))
                r.append(f"{s[0]:.2f} / {s[1]:.1f} / {s[2]:.0f} % / {np.median(o['pol']):.1f}")
            rows.append(r)
        if kind == "plate":      # free-space whips cannot tell up from down: give the estimator the hemisphere
            for nm in ("fixed", "az1"):
                r = [lab[nm] + ", estimator told the source is above the horizon"]
                for L in (180, 720):
                    o = sphere_mc(arr, f, path(kind, nm, L, K), az2, el2, psi, seed=[63, L], upper=True)
                    s = summ(o)
                    out[key]["home_sector"].setdefault(nm + "_upper", {})[L] = dict(ang=s, pol=summ(o, "pol"))
                    r.append(f"{s[0]:.2f} / {s[1]:.1f} / {s[2]:.0f} % / {np.median(o['pol']):.1f}")
                rows.append(r)
        table(f"f3. {kind}, {f / 1e9:.3f} GHz, {where} (where the fixed array already works). Same columns",
              ["path", "L = 180 deg", "L = 720 deg"], rows)
    res["f"] = out


def plots(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if "a" in res:
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        for name, v in res["a"].items():
            Ks = [k for k in v if isinstance(k, int) or str(k).isdigit()]
            ax[0].plot([int(k) for k in Ks], [v[k]["total64"][2] for k in Ks], "o-", label=name)
            ax[1].plot([int(k) for k in Ks], [v[k]["total64"][0] for k in Ks], "o-", label=name)
        ax[0].set_ylabel("gross errors (%)"); ax[1].set_ylabel("median error (deg)"); ax[1].set_yscale("log")
        for a in ax:
            a.set_xlabel("orientations K (64 samples in total)"); a.set_xscale("log", base=2); a.grid(alpha=0.3)
        ax[0].legend(fontsize=6)
        fig.tight_layout(); fig.savefig(OUT / "a_gross_vs_K.png", dpi=120); plt.close(fig)
    if "f" in res:
        keys = [k for k in res["f"]]
        fig, axes = plt.subplots(2, len(keys), figsize=(5 * len(keys), 7), squeeze=False)
        for j, key in enumerate(keys):
            for nm, v in res["f"][key].items():
                if nm in ("designs", "imperfect", "home_sector"):
                    continue
                Ls = sorted(int(x) for x in v)
                g = lambda L: v[L] if L in v else v[str(L)]
                axes[0, j].plot(Ls, [g(L)["ang"][2] for L in Ls], "o-", label=nm)
                axes[1, j].plot(Ls, [g(L)["pol"][0] for L in Ls], "o-", label=nm)
            axes[0, j].set_title(f"{key}: gross errors (%), sources on the whole sphere", fontsize=9)
            axes[1, j].set_title(f"{key}: median polarisation-angle error (deg)", fontsize=9)
            for a in axes[:, j]:
                a.set_xscale("log", base=2); a.set_xlabel("path length (deg of axis travel)"); a.grid(alpha=0.3)
            axes[0, j].legend(fontsize=6)
        fig.tight_layout(); fig.savefig(OUT / "f_paths.png", dpi=120); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--only", default="0,a,b,c,d,e,f")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    res: dict = {}
    secs = dict(a=sec_a, b=sec_b, c=sec_c, d=sec_d, e=sec_e, f=sec_f)
    secs["0"] = sec0
    t00 = time.time()
    for s in args.only.split(","):
        t0 = time.time()
        secs[s](res, args.quick)
        log(f"[section {s}: {time.time() - t0:.0f} s]")
    log(f"[total {time.time() - t00:.0f} s]")
    sfx = ("_quick" if args.quick else "") + ("" if args.only == "0,a,b,c,d,e,f" else "_" + args.only.replace(",", ""))
    (OUT / f"results{sfx}.json").write_text(json.dumps(res, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
    (OUT / f"tables{sfx}.md").write_text("\n".join(LINES) + "\n")
    try:
        plots(res)
    except Exception as e:  # plotting must not lose the tables
        log(f"plot failed: {e}")


if __name__ == "__main__":
    main()
