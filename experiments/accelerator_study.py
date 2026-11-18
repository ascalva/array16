"""Charged particles in fields this project could make: what contains a particle's path, what accelerates it.

The owner's idea (2026-11-16): a field that acts as a small particle accelerator, with a knotted field containing
the particle path, and the particles' own fields refining the geometry. This script puts numbers on it.

  scales     what the array's own fields do to an electron (nothing useful in free space), and what it would take
  validate   the integrator against theory: gyration, Mathieu stability, secular frequency, damped trapping
  ecr        electron cyclotron resonance: a magnet plus RF, the one resonant way to pump electrons at milliwatts
  dust       a trefoil-knot quadrupole guide for charged dust in air (60 Hz, kV): the table-top version
  rfguide    the same knot at 1 GHz for electrons in vacuum (after Hoffrogge et al. 2011)
  solenoid   a solenoid wound on a trefoil knot: electrons follow the knotted magnetic axis; a plain ring loses them
  ring       four RF gaps on that knot driven by one chirp with fixed phase offsets: a toy synchrotron
  selffield  where a particle's own field helps and where it hurts (closed-form numbers)

    python experiments/accelerator_study.py            # everything, about 4 minutes; results in experiments/accelerator/

Models: point charges, relativistic Boris integrator with exact gyro-rotation, thin-wire electrodes (line charges),
filament coils (Biot-Savart), no image charges, no gas collisions, no space charge in the tracking. Idealised: upper
bounds on what hardware would do.
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

E, ME, C0, EPS0, MU0, KB = 1.602176634e-19, 9.1093837015e-31, 299792458.0, 8.8541878128e-12, 4e-7 * math.pi, 1.380649e-23
QM_E = -E / ME
OUT = Path(__file__).resolve().parent / "accelerator"
RES: dict = {}


def log(*a):
    print(*a, flush=True)


def v_of(w_ev, m=ME):
    g = 1 + w_ev * E / (m * C0 ** 2)
    return C0 * math.sqrt(1 - 1 / g ** 2)


# ---- integrator ---------------------------------------------------------------------------------------------
def cross(a, b):
    return np.stack([a[..., 1] * b[..., 2] - a[..., 2] * b[..., 1],
                     a[..., 2] * b[..., 0] - a[..., 0] * b[..., 2],
                     a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]], axis=-1)


def push(x, u, qm, dt, ef, bf, drag=0.0, acc=None):
    """One relativistic Boris step for N particles (u = gamma v). The magnetic rotation uses the tangent of the
    exact half angle, so the gyration frequency has no time-step error. drag: Stokes rate (1/s); acc: (3,) m/s^2."""
    k = (qm * dt / 2)[:, None]
    um = u + k * ef
    if acc is not None:
        um = um + acc * (dt / 2)
    gam = np.sqrt(1 + np.sum(um * um, 1) / C0 ** 2)
    bn = np.sqrt(np.sum(bf * bf, 1))
    tv = bf * (np.tan(qm * dt / 2 * bn / gam) / np.where(bn > 0, bn, 1.0))[:, None]
    up = um + cross(um, tv)
    u = um + cross(up, tv * (2 / (1 + np.sum(tv * tv, 1)))[:, None]) + k * ef
    if acc is not None:
        u = u + acc * (dt / 2)
    if drag:
        u = u * math.exp(-drag * dt)
    gam = np.sqrt(1 + np.sum(u * u, 1) / C0 ** 2)
    return x + u / gam[:, None] * dt, u


# ---- curves, frames, electrodes, coils ----------------------------------------------------------------------
def resample(c, n):
    d = np.linalg.norm(np.roll(c, -1, 0) - c, axis=1)
    s = np.concatenate([[0], np.cumsum(d)])
    cc = np.vstack([c, c[:1]])
    si = np.arange(n) * s[-1] / n
    return np.stack([np.interp(si, s, cc[:, k]) for k in range(3)], 1), float(s[-1])


def torus_knot(R, r, n, p=2, q=3):
    """The (2,3) torus knot is the trefoil. Returns n points equally spaced in arc length, and the length."""
    t = np.linspace(0, 2 * np.pi, 40 * n, endpoint=False)
    return resample(np.stack([(R + r * np.cos(q * t)) * np.cos(p * t), (R + r * np.cos(q * t)) * np.sin(p * t),
                              r * np.sin(q * t)], 1), n)


def circle(R, n):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.stack([R * np.cos(t), R * np.sin(t), 0 * t], 1), 2 * math.pi * R


def frame(c, close_mod):
    """Unit tangent and a rotation-minimising normal frame (A, B) along a closed curve, with a uniform extra twist so
    that after one lap the frame returns rotated by a multiple k of close_mod (pi: a quadrupole's like-signed rods
    meet; 2 pi: a winding closes)."""
    n = len(c)
    T = np.roll(c, -1, 0) - np.roll(c, 1, 0)
    T /= np.linalg.norm(T, axis=1)[:, None]
    n1 = np.empty_like(c)
    v = np.cross(T[0], [0.0, 0.0, 1.0])
    if np.linalg.norm(v) < 1e-6:
        v = np.cross(T[0], [1.0, 0.0, 0.0])
    n1[0] = v / np.linalg.norm(v)
    for i in range(1, n):
        v = n1[i - 1] - np.dot(n1[i - 1], T[i]) * T[i]
        n1[i] = v / np.linalg.norm(v)
    v = n1[-1] - np.dot(n1[-1], T[0]) * T[0]
    v /= np.linalg.norm(v)
    alpha = math.atan2(np.dot(v, np.cross(T[0], n1[0])), np.dot(v, n1[0]))
    k = round(alpha / close_mod)
    tw = -(alpha - k * close_mod) * np.arange(n) / n
    n2 = np.cross(T, n1)
    ct, st = np.cos(tw)[:, None], np.sin(tw)[:, None]
    return T, ct * n1 + st * n2, -st * n1 + ct * n2, k


def curvature_radius(c, L):
    ds = L / len(c)
    k = np.linalg.norm(np.roll(c, -1, 0) - 2 * c + np.roll(c, 1, 0), axis=1) / ds ** 2
    return float(1 / k.max())


def quad_wires(c, r0):
    """Four thin rods at r0 round the curve, alternating sign. Returns segment midpoints (4, n, 3), lengths (4, n), signs."""
    T, A, B, k = frame(c, math.pi)
    W = np.array([c + r0 * (math.cos(j * math.pi / 2) * A + math.sin(j * math.pi / 2) * B) for j in range(4)])
    nxt = np.empty_like(W)
    for j in range(4):
        nxt[j, :-1] = W[j, 1:]
        nxt[j, -1] = W[(j + 2 * k) % 4, 0]            # after a lap a rod meets the like-signed rod opposite
    return T, (W + nxt) / 2, np.linalg.norm(nxt - W, axis=2), np.array([1.0, -1.0, 1.0, -1.0]), (A, B)


def e_wires(x, I, mid, dl, sign, lam, offs):
    """Field (N, 3) of the four line charges (sign * lam C/m) from the segments at stations I + offs."""
    idx = (I[:, None] + offs[None, :]) % mid.shape[1]
    out = np.zeros_like(x)
    for j in range(4):
        r = x[:, None, :] - mid[j][idx]
        d2 = np.sum(r * r, axis=2)
        out += sign[j] * np.sum(r * (dl[j][idx] / (d2 * np.sqrt(d2)))[:, :, None], axis=1)
    return out * (np.asarray(lam) / (4 * math.pi * EPS0)).reshape(-1, 1)


def helix(c, a_c, turns, per_turn=36):
    """A single-layer winding of `turns` turns, radius a_c, around the closed curve c. Returns midpoints, dl vectors."""
    T, A, B, _ = frame(c, 2 * math.pi)
    n, m = len(c), turns * per_turn
    f = np.arange(m) / m * n
    i0 = np.floor(f).astype(int)
    w = (f - i0)[:, None]
    i1 = (i0 + 1) % n
    ax, aa, bb = c[i0] * (1 - w) + c[i1] * w, A[i0] * (1 - w) + A[i1] * w, B[i0] * (1 - w) + B[i1] * w
    ph = 2 * math.pi * turns * np.arange(m) / m
    P = ax + a_c * (np.cos(ph)[:, None] * aa + np.sin(ph)[:, None] * bb)
    nxt = np.roll(P, -1, 0)
    return (P + nxt) / 2, nxt - P


def b_coil(x, mid, dlv):
    """Biot-Savart field (N, 3) per ampere of the filament segments."""
    r = x[:, None, :] - mid[None, :, :]
    d2 = np.sum(r * r, axis=2)
    k3 = 1 / (d2 * np.sqrt(d2))
    return MU0 / (4 * math.pi) * np.stack([
        np.sum((dlv[None, :, 1] * r[:, :, 2] - dlv[None, :, 2] * r[:, :, 1]) * k3, 1),
        np.sum((dlv[None, :, 2] * r[:, :, 0] - dlv[None, :, 0] * r[:, :, 2]) * k3, 1),
        np.sum((dlv[None, :, 0] * r[:, :, 1] - dlv[None, :, 1] * r[:, :, 0]) * k3, 1)], 1)


def track_station(x, I, c, adv, iters=1):
    n = len(c)
    for _ in range(iters):
        d0 = ((x - c[I]) ** 2).sum(1)
        ip, im = (I + 1) % n, (I - 1) % n
        fw = ((x - c[ip]) ** 2).sum(1) < d0
        bw = (~fw) & (((x - c[im]) ** 2).sum(1) < d0)
        I = np.where(fw, ip, np.where(bw, im, I))
        adv = adv + fw - bw
    return I, adv


def offset(x, I, c, T):
    d = x - c[I]
    d = d - (d * T[I]).sum(1)[:, None] * T[I]
    return d, np.sqrt((d * d).sum(1))


# ---- 1. scales ----------------------------------------------------------------------------------------------
def scales():
    log("\n== scales ==")
    rows = {}
    for f in (0.4e9, 2.45e9, 5.8e9):
        w = 2 * math.pi * f
        coef = E / (4 * ME * w * w)                      # ponderomotive energy in eV per (V/m)^2
        rows[f"{f / 1e9:g} GHz"] = dict(
            up_eV_at_17Vm=coef * 17 ** 2, up_eV_at_247Vm=coef * 247 ** 2,
            E_for_25meV_kVm=math.sqrt(0.025 / coef) / 1e3, E_for_1eV_kVm=math.sqrt(1 / coef) / 1e3,
            cyclotron_field_mT=w * ME / E * 1e3)
        r = rows[f"{f / 1e9:g} GHz"]
        log(f"  {f / 1e9:4.2f} GHz: ponderomotive depth {r['up_eV_at_17Vm']:.1e} eV at 17 V/m (array at 0.3 m), "
            f"{r['up_eV_at_247Vm']:.1e} eV at 247 V/m (all power in one spot); needs {r['E_for_25meV_kVm']:.0f} kV/m for "
            f"25 meV, {r['E_for_1eV_kVm']:.0f} kV/m for 1 eV; cyclotron resonance at {r['cyclotron_field_mT']:.1f} mT")
    # a resonant gap: V = sqrt(2 P R) with an assumed shunt resistance (estimate)
    gap = {f"{p * 1e3:g} mW": [math.sqrt(2 * p * r) for r in (0.25e6, 1e6)] for p in (3.2e-3, 1.0)}
    log("  gap voltage for an assumed shunt resistance of 0.25-1 Mohm (estimate):",
        {k: [round(x) for x in v] for k, v in gap.items()}, "V")
    # magnetic rigidity: the field-times-radius that bends an electron
    rig = {f"{w} eV": ME * v_of(w) / math.sqrt(1 - (v_of(w) / C0) ** 2) / E for w in (100, 1e3, 1e4, 1e5)}
    log("  magnetic rigidity B*rho (T m):", {k: f"{v:.2e}" for k, v in rig.items()})
    # mean free path of electrons in air (estimate: cross-section 1e-19 m^2)
    mfp = {p: KB * 293 / (p * 100 * 1e-19) for p in (1013, 1e-3, 1e-5, 1e-7)}      # p in mbar
    log("  electron mean free path (estimate, 1e-19 m^2):", {f"{p:g} mbar": f"{v:.2g} m" for p, v in mfp.items()})
    RES["scales"] = dict(ponderomotive=rows, gap_voltage=gap, rigidity_Tm=rig, mean_free_path_m={str(k): v for k, v in mfp.items()})


# ---- 2. validation ------------------------------------------------------------------------------------------
def quad_field(g):
    """Ideal quadrupole: E = g cos(wt) (-x, +y, 0), one g per particle; returns a function of (x, c)."""
    return lambda x, c: np.stack([-g * x[:, 0], g * x[:, 1], 0 * g], 1) * c


def validate():
    log("\n== validate ==")
    ok = {}
    # (a) gyration: 100 eV electron in 10 mT
    v, b = v_of(100), 0.01
    T = 2 * math.pi / (E / ME * b) / math.sqrt(1 - (v / C0) ** 2)      # relativistic period
    x, u = np.zeros((1, 3)), np.array([[v / math.sqrt(1 - (v / C0) ** 2), 0, 0]])
    dt, ys = T / 40, []
    for _ in range(40 * 50):
        x, u = push(x, u, np.array([QM_E]), dt, np.zeros((1, 3)), np.array([[0, 0, b]]))
        ys.append(x[0, 1])
    diam, want = max(ys) - min(ys), 2 * ME * v / math.sqrt(1 - (v / C0) ** 2) / (E * b)
    ang = math.atan2(u[0, 1], u[0, 0])
    ok["gyration"] = abs(diam / want - 1) < 0.01 and abs(ang) < 1e-6
    log(f"  gyration: diameter {diam * 1e3:.3f} mm (theory {want * 1e3:.3f}), velocity angle after 50 turns {ang:.1e} rad: "
        f"{'pass' if ok['gyration'] else 'FAIL'}")
    # (b) Mathieu: stable below q = 0.908, unstable above; secular frequency q Omega / (2 sqrt 2)
    om, qs = 2 * math.pi, np.array([0.2, 0.85, 1.0])
    fld = quad_field(qs * om ** 2 / 2)                    # charge/mass 1, r0 1: q = 2 g / Omega^2
    x, u = np.tile([1e-3, 1e-3, 0.0], (3, 1)), np.zeros((3, 3))
    per, nper, hist = 100, 400, []
    for it in range(per * nper):
        t = it / per
        x, u = push(x, u, np.ones(3), 1 / per, fld(x, math.cos(om * t)), np.zeros((3, 3)))
        if it % per == 0:
            hist.append(x[:, 0].copy())
    h = np.array(hist)
    grow = np.abs(h[-100:]).max(0) / 1e-3
    up = [i + h[i, 0] / (h[i, 0] - h[i + 1, 0]) for i in range(len(h) - 1) if h[i, 0] < 0 <= h[i + 1, 0]]
    beta = 2 * (len(up) - 1) / (up[-1] - up[0])           # secular frequency = beta Omega / 2
    ok["mathieu"] = grow[1] < 50 and grow[2] > 1e3 and abs(beta / (0.2 / math.sqrt(2)) - 1) < 0.02
    log(f"  Mathieu: growth at q = 0.85 {grow[1]:.1f}x (stable), at q = 1.0 {grow[2]:.1e}x (unstable); beta at q = 0.2 "
        f"{beta:.4f} (first order {0.2 / math.sqrt(2):.4f}): {'pass' if ok['mathieu'] else 'FAIL'}")
    # (c) dust in air: stability with Stokes drag, and gravity sag against the damped pseudo-potential
    m, qm, gam, f = 1.17e-12, 0.0176, 1450.0, 60.0        # 10 um borosilicate sphere, Vinitsky et al. 2014
    om, r0, a = 2 * math.pi * f, 6e-3, 0.75e-3
    geff = lambda V: 2 * V / math.log(r0 / a) / r0 ** 2    # field gradient of four thin rods, rod-to-rod amplitude V
    Vs = np.array([250, 500, 1000, 2000, 4000, 8000, 16000, 32000.0])
    fld = quad_field(geff(Vs))
    N = len(Vs)
    x, u = np.tile([2e-4, 2e-4, 0.0], (N, 1)), np.zeros((N, 3))
    per, nper = 200, 300
    dt, sag, amp = 1 / (f * per), np.zeros(N), np.zeros(N)
    for it in range(per * nper):
        x, u = push(x, u, np.full(N, qm), dt, fld(x, math.cos(om * it * dt)), np.zeros((N, 3)), drag=gam,
                    acc=np.array([0.0, -9.81, 0.0]))
        if it >= per * (nper - 50):
            sag += x[:, 1] / (per * 50)
            amp = np.maximum(amp, np.abs(x).max(1))
    qM = 2 * qm * geff(Vs) / om ** 2
    pred = -9.81 / ((qm * geff(Vs)) ** 2 / (2 * (om ** 2 + gam ** 2)))
    stable = amp < 0.8 * r0
    log("  dust (10 um glass, 60 Hz, r0 6 mm, air): V, Mathieu q, held?, sag simulated / damped pseudo-potential (mm)")
    for V, q_, s_, sg, pr in zip(Vs, qM, stable, sag, pred):
        log(f"     {V:6.0f} V  q {q_:6.1f}  {'held' if s_ else 'LOST'}  {sg * 1e3:8.3f} / {pr * 1e3:8.3f}")
    i = int(np.argmin(np.abs(Vs - 1000)))
    ok["dust"] = bool(stable[i]) and abs(sag[i] / pred[i] - 1) < 0.5
    RES["validate"] = dict(passed=ok, dust=dict(V=Vs.tolist(), q=qM.tolist(), held=stable.tolist(),
                                                 sag_mm=(sag * 1e3).tolist(), sag_pred_mm=(pred * 1e3).tolist()))
    assert all(ok.values()), ok
    return dict(m=m, qm=qm, gam=gam, f=f, r0=r0, a=a)


# ---- 3. electron cyclotron resonance -------------------------------------------------------------------------
def ecr_rwa(e0, delta, f=2.45e9, tmax=None, n=40000):
    """Rotating-wave equations for an electron in a uniform field B = B_res (1 + delta) and a linearly polarised RF
    field e0 (V/m, peak). a = p_perp / (m c) in the frame rotating at the RF frequency. Returns t, W (eV)."""
    w = 2 * math.pi * f
    kap = E * e0 / (2 * ME * C0)
    tstar = (12 * math.pi * ME ** 2 * C0 ** 2 / (w * E ** 2 * e0 ** 2)) ** (1 / 3)
    tmax = tmax or 6 * tstar
    dt, a, ts, ws = tmax / n, 0j, [], []
    rhs = lambda z: kap - 1j * (w * (1 + delta) / math.sqrt(1 + abs(z) ** 2) - w) * z
    for i in range(n):
        k1 = rhs(a); k2 = rhs(a + dt / 2 * k1); k3 = rhs(a + dt / 2 * k2); k4 = rhs(a + dt * k3)
        a += dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if i % 40 == 0:
            ts.append((i + 1) * dt); ws.append((math.sqrt(1 + abs(a) ** 2) - 1) * ME * C0 ** 2 / E)
    return np.array(ts), np.array(ws), tstar


def ecr():
    log("\n== ecr ==")
    f = 2.45e9
    w = 2 * math.pi * f
    b0 = w * ME / E
    out = {}
    for e0 in (17.0, 250.0, 1e3, 1e4):
        row = {}
        for d in (-1e-3, -1e-4, 0.0, 1e-4, 1e-3, 1e-2):
            t, wv, ts = ecr_rwa(e0, d)
            row[f"{d:g}"] = float(wv.max())
        t, wv, ts = ecr_rwa(e0, 0.0)
        out[f"{e0:g}"] = dict(wmax_eV=row, t_peak_us=float(t[np.argmax(wv)] * 1e6), tstar_us=ts * 1e6)
        log(f"  E0 {e0:7.0f} V/m: peak energy {row['0']:8.1f} eV after {out[f'{e0:g}']['t_peak_us']:.3f} us (uniform field); "
            f"field off by -1e-3/-1e-4/+1e-4/+1e-3/+1e-2: " + " / ".join(f"{row[k]:.1f}" for k in ("-0.001", "-0.0001", "0.0001", "0.001", "0.01")))
    # full-orbit cross-check at 30 kV/m
    e0 = 3e4
    t, wv, ts = ecr_rwa(e0, 0.0)
    per = 100
    dt = 1 / (f * per)
    x, u, wmax = np.zeros((1, 3)), np.zeros((1, 3)), 0.0
    for it in range(int(6 * ts / dt)):
        ef = np.array([[e0 * math.cos(w * (it + 0.5) * dt), 0, 0]])
        x, u = push(x, u, np.array([QM_E]), dt, ef, np.array([[0, 0, b0]]))
        wmax = max(wmax, (math.sqrt(1 + float(np.sum(u * u)) / C0 ** 2) - 1) * ME * C0 ** 2 / E)
    agree = abs(wmax / wv.max() - 1)
    log(f"  cross-check at 30 kV/m: full orbit {wmax:.0f} eV, rotating-wave {wv.max():.0f} eV ({agree * 100:.1f} % apart)")
    assert agree < 0.05
    # passing through a resonance zone in a field gradient: one kick, random phase (closed form)
    kick = {}
    for e0 in (17.0, 1e3, 1e4):
        vpar, lb = v_of(5.0), 0.05                     # 5 eV along the field, gradient length 50 mm
        teff = math.sqrt(2 * math.pi / (w / lb * vpar))
        kick[f"{e0:g}"] = 0.5 * ME * (E / ME * e0 / 2 * teff) ** 2 / E
    log("  one pass through the resonance layer of a 50 mm-gradient magnet, 5 eV electron: energy step (eV)",
        {k: f"{v:.2e}" for k, v in kick.items()})
    RES["ecr"] = dict(b_res_mT=b0 * 1e3, uniform=out, boris_check=dict(e0=e0, full_orbit_eV=wmax, rwa_eV=float(wv.max())),
                      gradient_pass_eV=kick)
    return {e0: ecr_rwa(e0, 0.0)[:2] for e0 in (250.0, 1e3, 1e4)}


# ---- 4. knotted quadrupole guide: dust in air ------------------------------------------------------------------
def dust(p):
    log("\n== dust: trefoil quadrupole guide, 10 um glass spheres in air ==")
    n = 601
    c, L = torus_knot(0.060, 0.016, n)
    T, mid, dl, sign, _ = quad_wires(c, p["r0"])
    rc = curvature_radius(c, L)
    V, vph, lam_tw = 1000.0, 0.05, L / round(L / 0.03)
    lam = math.pi * EPS0 * V / math.log(p["r0"] / p["a"])
    om = 2 * math.pi * p["f"]
    slope = float(np.abs(T[:, 2]).max())
    need = (p["gam"] * vph + 9.81 * slope) / p["qm"]
    ea = 1.6 * need
    N = 3
    qm = p["qm"] * np.array([0.5, 1.0, 2.0])            # charge spread
    x, u = np.tile(c[0], (N, 1)) + 3e-4, np.zeros((N, 3))
    I, adv = np.zeros(N, int), np.zeros(N)
    offs = np.arange(n) - n // 2
    per = 200
    dt, ds = 1 / (p["f"] * per), L / n
    kw, ww = 2 * math.pi / lam_tw, 2 * math.pi * vph / lam_tw
    alive, maxoff, traj, zsum, cnt = np.ones(N, bool), np.zeros(N), [], np.zeros(N), 0
    tlap = np.full(N, np.nan)
    t0, steps = time.time(), int(1.4 * L / vph / dt)
    for it in range(steps):
        t = it * dt
        ef = e_wires(x, I, mid, dl, sign, lam, offs) * math.cos(om * t)
        s = I * ds + ((x - c[I]) * T[I]).sum(1)
        ef += T[I] * (ea * np.sin(kw * s - ww * t))[:, None]
        xn, un = push(x, u, qm, dt, ef, np.zeros((N, 3)), drag=p["gam"], acc=np.array([0.0, 0.0, -9.81]))
        x, u = np.where(alive[:, None], xn, x), np.where(alive[:, None], un, u)
        I, adv = track_station(x, I, c, adv)
        if it % 20 == 0:
            d, off = offset(x, I, c, T)
            maxoff = np.where(alive, np.maximum(maxoff, off), maxoff)
            alive &= off < 0.8 * p["r0"]
            zsum += d[:, 2]; cnt += 1
            tlap = np.where(np.isnan(tlap) & (adv >= n), t, tlap)
            if it % 400 == 0:
                traj.append(x[1].copy())
            if np.all(~alive | (adv >= n)):
                break
    log(f"  knot: length {L * 1e3:.0f} mm, tightest bend {rc * 1e3:.0f} mm, steepest slope {slope:.2f}; rods +/-{V / 2:.0f} V at "
        f"{p['f']:.0f} Hz; conveyor wave {ea / 1e3:.1f} kV/m ({ea * lam_tw / (2 * math.pi):.0f} V on {lam_tw * 1e3:.0f} mm), {vph * 100:.0f} cm/s")
    for i in range(N):
        log(f"     charge x{qm[i] / p['qm']:.1f}: {'held' if alive[i] else 'LOST'}, laps {adv[i] / n:.2f}, lap time {tlap[i]:.1f} s, "
            f"largest offset {maxoff[i] * 1e3:.2f} mm of {p['r0'] * 1e3:.0f}, mean height offset {zsum[i] / cnt * 1e3:+.3f} mm")
    log(f"  ({time.time() - t0:.0f} s)")
    # Coulomb chain: spacing below which a line of grains buckles into a zigzag (closed form)
    Q = p["qm"] * p["m"]
    kt = p["m"] * (p["qm"] * 2 * V / math.log(p["r0"] / p["a"]) / p["r0"] ** 2) ** 2 / (2 * (om ** 2 + p["gam"] ** 2))
    dc = (4.207 * Q * Q / (4 * math.pi * EPS0) / kt) ** (1 / 3)
    log(f"  grains closer than {dc * 1e3:.2f} mm along the guide buckle from a line into a zigzag (their own fields): "
        f"more than {L / dc:.0f} grains on this knot")
    RES["dust"] = dict(length_mm=L * 1e3, bend_mm=rc * 1e3, V=V, f=p["f"], conveyor_kVm=ea / 1e3, conveyor_V=ea * lam_tw / (2 * math.pi),
                       held=alive.tolist(), laps=(adv / n).tolist(), lap_s=tlap.tolist(), max_offset_mm=(maxoff * 1e3).tolist(),
                       zigzag_spacing_mm=dc * 1e3)
    return c, np.array(traj)


# ---- 5. the same knot at 1 GHz for electrons -------------------------------------------------------------------
def rfguide():
    log("\n== rfguide: trefoil quadrupole guide at 1 GHz, electrons in vacuum ==")
    n, r0, f = 2001, 0.5e-3, 1e9
    c, L = torus_knot(0.020, 0.006, n)
    T, mid, dl, sign, (A, _) = quad_wires(c, r0)
    rc = curvature_radius(c, L)
    om = 2 * math.pi * f
    ws = np.array([0.5, 1, 2, 4, 8, 0.5, 1, 2, 4, 8, 16, 32.0])
    qM = np.array([0.3] * 5 + [0.6] * 7)
    lam = math.pi * EPS0 * (qM * om ** 2 / (2 * E / ME)) * r0 ** 2 / 2      # line charge for that Mathieu q
    vrod = lam / (math.pi * EPS0) * math.log(r0 / 25e-6)                    # rod-to-rod volts for 50 um wire
    N = len(ws)
    v = np.array([v_of(w) for w in ws])
    x, u = np.tile(c[0], (N, 1)) + 3e-5 * A[0], v[:, None] * T[0][None, :]
    I, adv, alive, maxoff = np.zeros(N, int), np.zeros(N), np.ones(N, bool), np.zeros(N)
    offs = np.arange(-60, 61)
    dt, laps_want, traj = 1 / (f * 30), 5, []
    t0 = time.time()
    for it in range(int(laps_want * L / v.min() / dt)):
        ef = e_wires(x, I, mid, dl, sign, lam, offs) * math.cos(om * it * dt)
        xn, un = push(x, u, np.full(N, QM_E), dt, ef, np.zeros((N, 3)))
        x, u = np.where(alive[:, None], xn, x), np.where(alive[:, None], un, u)
        I, adv = track_station(x, I, c, adv)
        if it % 10 == 0:
            _, off = offset(x, I, c, T)
            maxoff = np.where(alive, np.maximum(maxoff, off), maxoff)
            alive &= off < 0.8 * r0
            if it % 200 == 0 and adv[1] < n:
                traj.append(x[1].copy())
            if not alive.any():
                break
    wsec = qM / (2 * math.sqrt(2)) * om
    depth = 0.5 * ME * wsec ** 2 * r0 ** 2 / E
    wlim = 0.25 * ME * rc * wsec ** 2 * r0 / E          # centrifugal shift of r0 / 2
    log(f"  knot: length {L * 1e3:.0f} mm, tightest bend {rc * 1e3:.1f} mm, rods {r0 * 1e3:.1f} mm from the axis")
    for q in (0.3, 0.6):
        m_ = qM == q
        log(f"   q {q}: rod-to-rod {vrod[m_][0]:.0f} V, secular {wsec[m_][0] / 2 / math.pi / 1e6:.0f} MHz, depth {depth[m_][0]:.2f} eV, "
            f"centrifugal estimate of the largest energy {wlim[m_][0]:.1f} eV")
        for w_, al, la, mo in zip(ws[m_], alive[m_], adv[m_] / n, maxoff[m_]):
            log(f"      {w_:4.1f} eV: {'held' if al else 'LOST'} after {la:5.1f} laps, largest offset {mo * 1e6:5.0f} um")
    log(f"  ({time.time() - t0:.0f} s)")
    RES["rfguide"] = dict(length_mm=L * 1e3, bend_mm=rc * 1e3, r0_mm=r0 * 1e3, energies_eV=ws.tolist(), q=qM.tolist(),
                          rod_volts=vrod.tolist(), held=alive.tolist(), laps=(adv / n).tolist(), max_offset_um=(maxoff * 1e6).tolist(),
                          depth_eV=depth.tolist(), centrifugal_limit_eV=wlim.tolist())
    return c, np.array(traj)


# ---- 6. solenoid wound on a knot: the magnetic guide ------------------------------------------------------------
def solenoid():
    log("\n== solenoid: a winding on a trefoil knot, electrons in vacuum ==")
    n, a_c, turns = 600, 0.011, 180
    out, keep = {}, {}
    ck, Lk = torus_knot(0.060, 0.028, n)
    shapes = {"trefoil knot": (ck, Lk), "plain ring, same length": circle(Lk / (2 * math.pi), n)}
    for name, (c, L) in shapes.items():
        T = frame(c, 2 * math.pi)[0]
        mid, dlv = helix(c, a_c, turns)
        bax = np.vstack([b_coil(c[i:i + 40], mid, dlv) for i in range(0, n, 40)])
        bpar = (bax * T).sum(1)
        if bpar.mean() < 0:
            dlv, bax, bpar = -dlv, -bax, -bpar
        per_amp = float(np.linalg.norm(bax, axis=1).mean())
        log(f"  {name}: length {L * 1e3:.0f} mm, tightest bend {curvature_radius(c, L) * 1e3:.0f} mm, winding radius {a_c * 1e3:.0f} mm, "
            f"{turns} turns; axis field {per_amp * 1e3:.4f} mT/A (ideal {MU0 * turns / L * 1e3:.4f}), ripple "
            f"{np.linalg.norm(bax, axis=1).min() / per_amp:.2f}-{np.linalg.norm(bax, axis=1).max() / per_amp:.2f}, along the axis {(bpar / np.linalg.norm(bax, axis=1)).min():.3f}")
        ws = np.array([100.0, 300.0, 1000.0, 3000.0])
        rows = {}
        for bt in (5e-3, 10e-3, 20e-3):
            cur = bt / per_amp
            N = len(ws)
            v = np.array([v_of(w) for w in ws])
            A0 = frame(c, 2 * math.pi)[1][0]
            x = np.tile(c[0], (N, 1)) + 2e-3 * A0
            b0 = b_coil(x, mid, dlv)
            u = (v / np.sqrt(1 - (v / C0) ** 2))[:, None] * b0 / np.linalg.norm(b0, axis=1)[:, None]
            I, adv, alive, maxoff = np.zeros(N, int), np.zeros(N), np.ones(N, bool), np.zeros(N)
            dt = 2 * math.pi / (E / ME * bt) / 20
            laps_want, traj = 20, []
            for it in range(int(laps_want * L / v.min() / dt)):
                bf = b_coil(x, mid, dlv) * cur
                xn, un = push(x, u, np.full(N, QM_E), dt, np.zeros((N, 3)), bf)
                x, u = np.where(alive[:, None], xn, x), np.where(alive[:, None], un, u)
                I, adv = track_station(x, I, c, adv, iters=4)
                _, off = offset(x, I, c, T)
                maxoff = np.where(alive, np.maximum(maxoff, off), maxoff)
                alive &= off < 0.9 * a_c
                if bt == 10e-3 and it % 3 == 0 and alive[1]:
                    traj.append(x[1].copy())
                if not alive.any():
                    break
            rows[f"{bt * 1e3:g} mT"] = dict(held=alive.tolist(), laps=(adv / n).tolist(), max_offset_mm=(maxoff * 1e3).tolist(), amps=cur)
            log(f"     {bt * 1e3:4.0f} mT ({cur:5.1f} A in {turns} turns): " + "; ".join(
                f"{w_:.0f} eV {'held' if al else 'LOST'} {la:.1f} laps, {mo * 1e3:.1f} mm" for w_, al, la, mo in zip(ws, alive, adv / n, maxoff)))
            if bt == 10e-3:
                keep[name] = (c, np.array(traj))
        out[name] = dict(length_mm=L * 1e3, mT_per_A=per_amp * 1e3, cases=rows)
    # copper for a continuous 10 mT winding (estimate): 6 mm build, 60 % fill
    J = (10e-3 / MU0) / (0.006 * 0.6)
    vol = 2 * math.pi * (a_c + 0.003) * 0.006 * 0.6 * Lk
    out["copper_10mT"] = dict(A_per_mm2=J / 1e6, watts=J ** 2 * 1.72e-8 * vol, kg=vol * 8960)
    log(f"  continuous 10 mT winding (estimate): {J / 1e6:.1f} A/mm^2, {out['copper_10mT']['watts']:.0f} W, {out['copper_10mT']['kg']:.1f} kg of copper")
    RES["solenoid"] = out
    return keep, Lk


# ---- 7. four RF gaps on the knot: a toy synchrotron driven by one chirp ---------------------------------------
def ring(L):
    log("\n== ring: four gaps, one chirp, fixed phase offsets ==")
    w0, vg, phis, ngap, band = 50.0, 0.5, math.radians(30), 4, 54e6
    mc2 = ME * C0 ** 2 / E
    vel = lambda w: C0 * np.sqrt(1 - 1 / (1 + w / mc2) ** 2)
    h = round(60e6 / (vel(w0) / L))
    f0 = h * vel(w0) / L
    # synchronous particle: kicked by vg sin(phis) at each gap; phase table Phi(t_j) = j h / ngap cycles
    ts, wsyn, t, w = [-(L / ngap) / vel(w0), 0.0], [w0, w0 + vg * math.sin(phis)], 0.0, w0 + vg * math.sin(phis)
    while h * vel(w) / L < f0 + band:
        t += (L / ngap) / vel(w)
        w += vg * math.sin(phis)
        ts.append(t); wsyn.append(w)
    ts, wsyn = np.array(ts), np.array(wsyn)
    phi_tab = (np.arange(len(ts)) - 1) * h / ngap
    N = 4000
    rng = np.random.default_rng(1)
    tp = rng.uniform(-0.5, 0.5, N) / f0                    # a continuous beam: every RF phase
    wp = w0 + rng.uniform(-1, 1, N)
    sample = []
    for j in range(len(ts) - 1):
        ph = 2 * math.pi * (np.interp(tp, ts, phi_tab) - j * h / ngap) + phis      # phase seen at this gap
        wp = np.maximum(wp + vg * np.sin(ph), 0.5)
        tp = tp + (L / ngap) / vel(wp)
        if j % ngap == 0:
            sample.append(wp[:40].copy())
    cap = np.abs(wp - wsyn[-1]) < 0.05 * wsyn[-1]
    out = dict(harmonic=int(h), f_start_MHz=f0 / 1e6, f_end_MHz=h * vel(wsyn[-1]) / L / 1e6, laps=(len(ts) - 2) / ngap,
               time_us=float(ts[-1] * 1e6), w_start_eV=w0, w_end_eV=float(wsyn[-1]), captured=float(cap.mean()),
               spread_eV=float(wp[cap].std()) if cap.any() else None, gap_volts=vg, phase_step_deg=360 * (h / ngap % 1))
    log(f"  injected {w0:.0f} eV continuous beam, {ngap} gaps of {vg} V at harmonic {h}: chirp {out['f_start_MHz']:.1f} -> "
        f"{out['f_end_MHz']:.1f} MHz in {out['time_us']:.0f} us ({out['laps']:.0f} laps); synchronous energy {out['w_end_eV']:.0f} eV; "
        f"captured {out['captured'] * 100:.0f} % with {out['spread_eV']:.1f} eV rms spread")
    # what the beam itself puts on a pickup
    for ne in (1e5, 1e6, 1e7):
        i_avg = ne * E * vel(wsyn[-1]) / L
        log(f"     {ne:.0e} electrons circulating: {i_avg * 1e6:.2f} uA; about {20 * math.log10(0.3 * 2 * i_avg * 50 / math.sqrt(2) / 0.2236):.0f} dBm "
            f"on a 50 ohm pickup coupling 30 % (estimate)")
    RES["ring"] = out
    return np.array(sample), wsyn[1::ngap]


# ---- 8. self-fields (closed form) -----------------------------------------------------------------------------
def selffield():
    log("\n== selffield ==")
    out = {}
    for b in (5e-3, 10e-3, 20e-3):
        nb = EPS0 * b * b / (2 * ME)                      # Brillouin limit of a magnetised electron column
        out[f"{b * 1e3:g} mT"] = dict(brillouin_m3=nb, max_mA_100eV_3mm=nb * E * v_of(100) * math.pi * 0.003 ** 2 * 1e3)
    i, a, v = 1e-6, 0.003, v_of(100)
    er = i / (2 * math.pi * EPS0 * a * v)
    out["beam_1uA_100eV"] = dict(self_E_Vm=er, guide_vxB_Vm_at_10mT=v * 0.01, beta2=(v / C0) ** 2)
    n_line, kT = 1e15, 2 * E                              # Bennett pinch: mu0 I^2 = 8 pi N k (Te + Ti)
    out["bennett_A"] = math.sqrt(8 * math.pi * n_line * 2 * kT / MU0)
    log(f"  Brillouin limit: " + ", ".join(f"{k}: {v_['brillouin_m3']:.1e} /m^3 ({v_['max_mA_100eV_3mm']:.0f} mA in a 3 mm 100 eV beam)"
                                             for k, v_ in out.items() if "mT" in k))
    log(f"  1 uA, 100 eV, 3 mm beam: its own field {er:.2f} V/m against {v * 0.01:.1e} V/m of magnetic guiding; magnetic self-attraction is "
        f"{(v / C0) ** 2:.1e} of the electric self-repulsion")
    log(f"  a neutral plasma column of 1e15 particles/m at 2 eV pinches itself at {out['bennett_A']:.0f} A")
    RES["selffield"] = out


# ---- figures -----------------------------------------------------------------------------------------------------
def figures(dust_r, rf_r, sol_r, ecr_r, ring_r):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 4, figsize=(15, 4))
    for a_, (c, tr), ttl in ((ax[0], dust_r, "dust grain, 60 Hz knot guide in air: one lap"),
                             (ax[1], rf_r, "1 eV electron, 1 GHz knot guide: one lap"),
                             (ax[2], sol_r["trefoil knot"], "300 eV electron, 10 mT knotted solenoid: 20 laps"),
                             (ax[3], sol_r["plain ring, same length"], "same electron, plain ring solenoid: lost")):
        a_.plot(c[:, 0] * 1e3, c[:, 1] * 1e3, color="0.75", lw=4)
        if len(tr):
            a_.plot(tr[:, 0] * 1e3, tr[:, 1] * 1e3, color="C3", lw=0.7)
        a_.set_aspect("equal"); a_.set_title(ttl, fontsize=8); a_.tick_params(labelsize=7); a_.set_xlabel("mm", fontsize=7)
    fig.tight_layout(); fig.savefig(OUT / "guides.png", dpi=130); plt.close(fig)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    for e0, (t, w) in ecr_r.items():
        ax[0].loglog(t * 1e6, np.maximum(w, 1e-3), label=f"{e0:g} V/m")
    ax[0].set_xlabel("microseconds"); ax[0].set_ylabel("eV"); ax[0].set_title("Cyclotron resonance at 2.45 GHz, 87.5 mT, uniform field", fontsize=9)
    ax[0].legend(fontsize=7); ax[0].grid(alpha=0.3, which="both")
    samp, wsyn = ring_r
    ax[1].plot(samp, color="C0", lw=0.4, alpha=0.6)
    ax[1].plot(wsyn[:len(samp)], "k--", lw=1.2, label="synchronous energy")
    ax[1].set_xlabel("lap"); ax[1].set_ylabel("eV"); ax[1].set_title("Four 0.5 V gaps on the knot, one chirp: 40 of 4000 electrons", fontsize=9)
    ax[1].legend(fontsize=7); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(OUT / "resonance_and_ring.png", dpi=130); plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    t0 = time.time()
    scales()
    p = validate()
    ecr_r = ecr()
    dust_r = dust(p)
    rf_r = rfguide()
    sol_r, Lk = solenoid()
    ring_r = ring(Lk)
    selffield()
    (OUT / "results.json").write_text(json.dumps(RES, indent=1, default=float))
    figures(dust_r, rf_r, sol_r, ecr_r, ring_r)
    log(f"\ndone in {time.time() - t0:.0f} s; results in {OUT}")


if __name__ == "__main__":
    sys.exit(main())
