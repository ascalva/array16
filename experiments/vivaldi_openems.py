#!/usr/bin/env python
"""openEMS model and optimisation of the rev B panel element: a connected antipodal Vivaldi on 0.8 mm FR-4.

    .env/bin/python experiments/vivaldi_openems.py launch     # start `all` detached, Mac kept awake
    .env/bin/python experiments/vivaldi_openems.py status     # progress, best design so far, ETA
    .env/bin/python experiments/vivaldi_openems.py all        # validate -> search -> check -> report (resumes)
    .env/bin/python experiments/vivaldi_openems.py validate | search | check | report

Everything is resumable. Each finished openEMS run is one JSON line in experiments/vivaldi/evals.jsonl, keyed by
a hash of (model, mesh and solver settings, rounded parameters, eps_r). A design is looked up there before it is
simulated, and the optimiser is seeded, so running the same command again after a crash replays the cache in
seconds and carries on where it stopped. `report` uses only the files on disk.

Coordinates (mm): x along a card (the E-plane; the E-field across a flare gap is along x), y normal to the cards
(the H-plane; cards repeat every 25 mm in y), z up from the tile ground (boresight). Time dependence exp(+jwt).

Element (one 25 mm cell, repeated six times along a 150 mm card):
  front copper (y = +0.4):  50 ohm microstrip up from the card foot, narrowing to the throat width, then fin A,
                            which flares towards +x;
  back copper (y = -0.4):   a ground strip along the card foot (soldered to the tile ground), a ground taper
                            under the microstrip (microstrip -> parallel strip, no balun), then fin B towards -x;
  cell edge:                fin A of one cell overlaps fin B of the next by 1 mm and a row of vias (modelled as
                            a solid wall through the laminate) joins them above height z_j. Below z_j the two
                            outlines part and leave an open window down to the ground strip. The joined fins
                            make the array "connected"; the window is the Vivaldi cavity (a shunt inductance).
  port:                     a 50 ohm lumped port through the laminate between strip and card ground, 1 mm above
                            the tile. The right-angle joint to the tile is NOT modelled.

Models
  row   one card, six cells, PMC walls at y = +/-12.5 mm (an infinite stack of cards, all phased alike: H-plane
        broadside; note the image cards are mirror images, so the fins of neighbouring cards face each other),
        PEC ground at z = 0, PML at +/-x and above. Cell 2 (0-based) is excited, the others are loads; cells 0
        and 5 are the 49.9 ohm dummies. One run gives column c of S and, by reciprocity, row c.
  tile  the finite 6 x 6 tile: six translated cards, a 160 mm square ground with air below, PML on all six sides,
        element (card 2, cell 2) excited, the outer ring in 49.9 ohm, NF2FF at 2.4, 3.5 and 5.8 GHz.

Active reflection of element c with the active set A (the inner four cells of a card, or the inner 4 x 4) driven
with a_j = exp(-j k r_j . u), u = (sin(th) cos(ph), sin(th) sin(ph)), which steers the beam to (th, ph):
    Gamma_c(th, ph) = sum_{j in A} S_cj exp(-j k (r_j - r_c) . u)
By the mirror symmetry of a card, -th for cell 2 is +th for cell 3, so both centre cells are covered by +/-th.

Score (lower is better): worst |Gamma_c| in dB over 2.3-6 GHz and th = 0, +/-30, +/-45 deg in the E-plane, plus
0.3 x the mean of |Gamma_c| in dB over the same set, plus 10 dB per 0.1 mm of manufacturing-rule violation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT = os.path.join(HERE, "vivaldi")
EVALS = os.path.join(OUT, "evals.jsonl")
SCRATCH = os.environ.get(
    "VIVALDI_SCRATCH",
    "/private/tmp/claude-501/-Users-ascalva--claude-spares-spare-iFcbz1/8b6bf486-3081-448c-8f46-b00ba3ebceb0/scratchpad/vivaldi_runs")
sys.path.insert(0, REPO)

# ----------------------------------------------------------------------------- fixed quantities
C0 = 299792458.0
EPS0 = 8.854187817e-12
P = 25.0            # mm, cell pitch
NCELL = 6           # cells per card, cards per tile
CARD = P * NCELL    # 150 mm
T = 0.8             # mm, laminate
TAND = 0.02
F_KAPPA = 4.0e9     # the loss is a constant conductivity giving tan_d = 0.02 at this frequency
W_MS = 1.5          # mm, 50 ohm microstrip on 0.8 mm FR-4
OV = 0.5            # mm, each fin reaches this far past the cell edge (1 mm overlap for the via lands)
GND_STRIP = 1.0     # mm, height of the ground strip along the card foot
MIN_FEATURE = 0.15  # mm
FREQ = np.linspace(2.3e9, 6.0e9, 38)
F_PAT = [2.4e9, 3.5e9, 5.8e9]
SCAN = [0.0, 30.0, -30.0, 45.0, -45.0]
EXC_CELL = 2
EXC_CARD = 2
GROUND_HALF = 80.0  # mm, half side of the tile ground in the tile model

NAMES = ["H", "R", "w_t", "z_t", "a_w", "c_j", "W_g"]
UNITS = ["mm", "1/mm", "mm", "mm", "mm", "-", "mm"]
WHAT = ["fin (card) height", "exponential opening rate of the flare", "throat (parallel-strip) width",
        "height of the throat above the tile = length of the feed taper", "aperture width at the top",
        "start of the joined section, as a fraction of the flare length (z_j = z_t + c_j (H - z_t))",
        "width of the ground taper at the card foot"]
LO = np.array([35.0, 0.02, 1.0, 3.0, 14.0, 0.05, 3.0])
HI = np.array([45.0, 0.30, 3.0, 15.0, 24.0, 0.90, 24.5])
STEP = np.array([0.1, 0.005, 0.05, 0.1, 0.1, 0.01, 0.1])   # designs are snapped to this lattice
X0 = np.array([40.0, 0.10, 1.5, 8.0, 22.0, 0.40, 12.0])    # starting design, see vivaldi_study.md

# ----------------------------------------------------------------------------- solver and run settings
ENGINE = os.environ.get("VIVALDI_ENGINE", "gpu")       # 15x faster here than the CPU engines; cross-checked in validate
THREADS = int(os.environ.get("VIVALDI_THREADS", "4"))      # per row-model run
JOBS = int(os.environ.get("VIVALDI_JOBS", "1"))            # row-model runs side by side (one GPU)
TILE_THREADS = int(os.environ.get("VIVALDI_TILE_THREADS", "8"))


def mesh_cfg(nh=35, nsub=3, gap=30.0, end=1e-4, engine=None):
    """nh: mesh cells per half pitch along the card (35 -> 0.357 mm; z uses the same size); nsub: cells through the
    laminate; gap: air between the structure and the start of each PML (mm); end: energy end criterion."""
    return dict(nh=int(nh), nsub=int(nsub), gap=float(gap), end=float(end), engine=engine or ENGINE,
                dair=2.5, dy=1.5, ratio=1.3, pml=8, f0=4.15e9, fc=2.4e9)


M_BASE = mesh_cfg()
M_COARSE = mesh_cfg(nh=25, nsub=2)
M_FINE1 = mesh_cfg(nh=50, nsub=4)
M_FINE2 = mesh_cfg(nh=70, nsub=4)
M_GAP = mesh_cfg(gap=50.0)
M_END = mesh_cfg(end=1e-5)
M_TILE = mesh_cfg(nh=40, nsub=3)
ALT_ENGINE = "multithreaded" if ENGINE != "multithreaded" else "basic"

CLIMB_STEPS, CLIMB_LAM, CLIMB_SIGMA, CLIMB_SEED = 25, 8, 0.08, 1
EVO_POP, EVO_GENS, EVO_ELITE, EVO_SIGMA, EVO_SEED = 12, 10, 2, 0.08, 2
N_FINAL = 4
EPS_NOM, EPS_LO, EPS_HI = 4.4, 4.2, 4.6
PLAN_SEARCH = 1 + CLIMB_STEPS * CLIMB_LAM + EVO_POP + EVO_GENS * (EVO_POP - EVO_ELITE) + 2 * (N_FINAL + 1)

# validation gates (decided before looking at results)
TOL_MESH_DB = 1.0      # |score(fine1) - score(fine2)| must be below this ("stopped moving")
TOL_MESH_G = 0.08      # and the largest change of |Gamma_act| between fine1 and fine2
TOL_BND_DB = 0.5       # score change when the PML moves from 30 to 50 mm away
TOL_BND_G = 0.05
TOL_PASSIVE = 1.005    # column power sum and |S| limit
TOL_ENGINE = 0.02      # largest |S| difference between two FDTD engines on the same mesh


# ----------------------------------------------------------------------------- small utilities
def log(*a):
    print(time.strftime("%Y-%m-%d %H:%M:%S"), *a, flush=True)


def atomic_write(path, text, mode="w"):
    tmp = "%s.tmp%d" % (path, os.getpid())
    with open(tmp, mode) as f:
        f.write(text)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def write_json(path, obj):
    atomic_write(path, json.dumps(obj, indent=1, default=_js))


def read_json(path):
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def _js(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    raise TypeError(type(o))


def snap(x):
    x = np.clip(np.asarray(x, float), LO, HI)
    return np.clip(np.round(x / STEP) * STEP, LO, HI).round(6)


def db(v):
    return 20 * np.log10(np.maximum(np.abs(v), 1e-12))


# ----------------------------------------------------------------------------- geometry
def spans(p, z):
    """Copper of one cell at height z, relative to the cell centre: (front span or None, back span)."""
    H, R, wt, zt, aw, cj, Wg = p
    zj = zt + cj * (H - zt)
    if z < zt:
        hw = W_MS / 2 + (wt / 2 - W_MS / 2) * z / zt
        hb = P / 2 if z < GND_STRIP else Wg / 2 + (wt / 2 - Wg / 2) * (z - GND_STRIP) / (zt - GND_STRIP)
        return (-hw, hw), (-hb, hb)
    s = math.expm1(R * (z - zt)) / math.expm1(R * (H - zt))
    xin = -wt / 2 + (aw / 2 + wt / 2) * s
    xout = P / 2 + OV if z >= zj else wt / 2 + (P / 2 + OV - wt / 2) * (z - zt) / (zj - zt)
    return (xin, xout), (-xout, -xin)


def penalty(p):
    """Millimetres by which the outline breaks the 0.15 mm copper/gap rule or the via-land rule (0 if legal).
    The parameter bounds were chosen so that this is normally zero; it guards the corners of the box."""
    H, R, wt, zt, aw, cj, Wg = p
    v = 0.0
    for z in np.linspace(0.0, H, 400, endpoint=False):
        f, b = spans(p, z)
        v = max(v, MIN_FEATURE - (f[1] - f[0]), MIN_FEATURE - (b[1] - b[0]))
        # same-layer gap to the neighbour's 1 mm via tab (front: tab at -P/2 .. -P/2+OV of this cell)
        v = max(v, MIN_FEATURE - (f[0] - (-P / 2 + OV)))
        if GND_STRIP <= z < zt and b[1] < P / 2:           # gap between neighbouring ground tapers
            v = max(v, MIN_FEATURE - 2 * (P / 2 - b[1]))
    v = max(v, 1.0 - (1 - cj) * (H - zt))                  # at least 1 mm of via row
    v = max(v, 0.6 - (P / 2 + OV - aw / 2))                # via land: 0.3 mm drill plus 0.15 mm ring each side
    v = max(v, H - 45.0)
    return max(v, 0.0)


def raster(p, zl, res):
    """Staircase the outline on the mesh: per z-row, copper from the snapped inner to the snapped outer edge.
    Returns (front boxes, back boxes, via walls) in card coordinates; boxes are (x1, x2, z1, z2)."""
    front, back, vias = [], [], []
    H, R, wt, zt, aw, cj, Wg = p
    zj = zt + cj * (H - zt)
    kj = int(np.argmin(np.abs(zl - zj)))
    sn = lambda v: round(v / res) * res
    for c in range(NCELL):
        xc = (c - (NCELL - 1) / 2) * P
        for k in range(len(zl) - 1):
            f, b = spans(p, 0.5 * (zl[k] + zl[k + 1]))
            for lst, s, skip in ((front, f, k == 0), (back, b, False)):   # the strip starts one row above the tile
                if skip:
                    continue
                a, e = sn(s[0]), sn(s[1])
                if e - a < res / 2:
                    e = a + res
                a, e = max(xc + a, -CARD / 2), min(xc + e, CARD / 2)
                if e > a:
                    if lst and lst[-1][3] == zl[k] and abs(lst[-1][0] - a) < 1e-9 and abs(lst[-1][1] - e) < 1e-9:
                        lst[-1] = (a, e, lst[-1][2], zl[k + 1])
                    else:
                        lst.append((a, e, zl[k], zl[k + 1]))
        if c < NCELL - 1:
            vias.append((xc + P / 2, zl[kj], zl[-1]))
    return front, back, vias


def grade(d0, dmax, length, ratio):
    """Offsets 0 < o_1 < ... = length with steps growing from about d0 by `ratio` up to dmax."""
    steps, d = [], d0
    while sum(steps) < length - 1e-9:
        d = min(d * ratio, dmax)
        steps.append(d)
    steps = np.array(steps) * (length / sum(steps))
    return np.cumsum(steps)


def mesh_lines(model, p, m):
    H = p[0]
    res = P / 2 / m["nh"]
    far = m["gap"] + m["pml"] * m["dair"]
    out = grade(res, m["dair"], far, m["ratio"])
    xin = np.linspace(-CARD / 2, CARD / 2, 2 * NCELL * m["nh"] + 1)
    x = np.concatenate([-CARD / 2 - out[::-1], xin, CARD / 2 + out])
    zin = np.linspace(0.0, H, int(round(H / res)) + 1)
    z = np.concatenate([zin, H + out])
    ysub = np.linspace(-T / 2, T / 2, m["nsub"] + 1)
    side = grade(T / m["nsub"], m["dy"], P / 2 - T / 2, m["ratio"])
    ycell = np.concatenate([-T / 2 - side[::-1], ysub, T / 2 + side])        # -12.5 .. 12.5
    if model == "row":
        y = ycell
    else:
        ys = [ycell[:-1] + (c - (NCELL - 1) / 2) * P for c in range(NCELL)]
        yout = grade(m["dy"], m["dair"], far, m["ratio"])
        y = np.concatenate([-CARD / 2 - yout[::-1]] + ys + [[CARD / 2], CARD / 2 + yout])
        z = np.concatenate([-out[::-1], z])
    return x, y, z, zin, res


def geom_hash(model, p, eps, m):
    x, y, z, zin, res = mesh_lines(model, p, m)
    fr, bk, vi = raster(p, zin, res)
    blob = json.dumps([model, m, round(float(eps), 4), np.round(zin, 5).tolist(),
                       np.round(np.array(fr), 5).tolist(), np.round(np.array(bk), 5).tolist(),
                       np.round(np.array(vi), 5).tolist()], sort_keys=True)
    return hashlib.sha1(blob.encode()).hexdigest()[:16]


def eval_key(model, p, eps, m):
    blob = json.dumps([model, m, [round(float(v), 6) for v in snap(p)], round(float(eps), 4)], sort_keys=True)
    return hashlib.sha1(blob.encode()).hexdigest()[:16]


def port_positions(model):
    """(x, y) of every port in mm, in port order (card by card), and the index of the excited one."""
    xs = (np.arange(NCELL) - (NCELL - 1) / 2) * P
    if model == "row":
        return np.array([[x, 0.0] for x in xs]), EXC_CELL
    pos = np.array([[x, (c - (NCELL - 1) / 2) * P] for c in range(NCELL) for x in xs])
    return pos, EXC_CARD * NCELL + EXC_CELL


def active_set(model):
    if model == "row":
        return [1, 2, 3, 4]
    return [c * NCELL + k for c in range(1, NCELL - 1) for k in range(1, NCELL - 1)]


# ----------------------------------------------------------------------------- the openEMS run (worker)
def run_openems(spec):
    """Build and run one model in spec['dir']; return the result dict. Runs in a worker process."""
    from CSXCAD import ContinuousStructure
    from openEMS import openEMS

    model, p, eps, m, wd = spec["model"], np.array(spec["x"], float), float(spec["eps"]), spec["mesh"], spec["dir"]
    x, y, z, zin, res = mesh_lines(model, p, m)
    front, back, vias = raster(p, zin, res)
    H = p[0]
    FDTD = openEMS(NrTS=80000, EndCriteria=m["end"])
    FDTD.SetGaussExcite(m["f0"], m["fc"])
    pml = "PML_%d" % m["pml"]
    FDTD.SetBoundaryCond([pml, pml, "PMC", "PMC", "PEC", pml] if model == "row" else [pml] * 6)
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    g = CSX.GetGrid()
    g.SetDeltaUnit(1e-3)
    g.AddLine("x", x)
    g.AddLine("y", y)
    g.AddLine("z", z)
    fr4 = CSX.AddMaterial("fr4", epsilon=eps, kappa=2 * np.pi * F_KAPPA * EPS0 * eps * TAND)
    cu = CSX.AddMetal("cu")
    cards = [0.0] if model == "row" else [(c - (NCELL - 1) / 2) * P for c in range(NCELL)]
    zp = float(zin[2])
    ports = []
    pos, exc = port_positions(model)
    act = set(active_set(model))
    for ic, yc in enumerate(cards):
        fr4.AddBox([-CARD / 2, yc - T / 2, 0], [CARD / 2, yc + T / 2, H], priority=1)
        for (a, e, z1, z2) in front:
            cu.AddBox([a, yc + T / 2, z1], [e, yc + T / 2, z2], priority=10)
        for (a, e, z1, z2) in back:
            cu.AddBox([a, yc - T / 2, z1], [e, yc - T / 2, z2], priority=10)
        for (xb, z1, z2) in vias:
            cu.AddBox([xb, yc - T / 2, z1], [xb, yc + T / 2, z2], priority=10)
        for k in range(NCELL):
            n = ic * NCELL + k
            xc = (k - (NCELL - 1) / 2) * P
            R = 50.0 if n in act else 49.9
            ports.append(FDTD.AddLumpedPort(n + 1, R, [xc, yc - T / 2, zp], [xc, yc + T / 2, zp], "y",
                                            excite=1.0 if n == exc else 0, priority=5))
    nf2ff = None
    if model == "tile":
        gh = float(x[np.argmin(np.abs(x - GROUND_HALF))])
        gy = float(y[np.argmin(np.abs(y - GROUND_HALF))])
        cu.AddBox([-gh, -gy, 0], [gh, gy, 0], priority=10)
        if spec.get("far", True):
            nf2ff = FDTD.CreateNF2FFBox(frequency=F_PAT, opt_resolution=[3.0, 3.0, 3.0])
    t0 = time.time()
    kw = dict(verbose=1, numThreads=int(spec["threads"]))
    if m["engine"] != "fastest":
        kw["engine"] = m["engine"]
    FDTD.Run(wd, cleanup=True, **kw)
    wall = time.time() - t0
    for pt in ports:
        pt.CalcPort(wd, FREQ)
    inc = ports[exc].uf_inc
    S = np.array([pt.uf_ref / inc for pt in ports])
    leak = max(float(np.max(np.abs(pt.uf_inc / inc))) for i, pt in enumerate(ports) if i != exc)
    res_d = dict(S_re=np.round(S.real, 6).tolist(), S_im=np.round(S.imag, 6).tolist(), wall=round(wall, 1),
                 cells=int((len(x) - 1) * (len(y) - 1) * (len(z) - 1)), passive_inc=leak,
                 nboxes=len(front) + len(back))
    if nf2ff is not None:
        res_d["far"] = farfield(nf2ff, ports[exc], wd)
    return res_d


def farfield(nf2ff, port, wd):
    """Embedded element pattern: realised gain (co and cross, Ludwig 3, co = x) on three cuts, and the power
    budget, at F_PAT. Realised gain = 4 pi r^2 S(th, ph) / P_incident: mismatch, coupling into the other ports
    and dielectric loss are all counted against the element."""
    th = np.arange(0.0, 180.1, 2.0)
    ph = np.arange(0.0, 360.0, 5.0)
    r = nf2ff.CalcNF2FF(wd, F_PAT, th, ph, center=[0, 0, 0], read_cached=False)
    fpat = np.array(F_PAT)
    port.CalcPort(wd, fpat)
    TH, PH = np.meshgrid(np.deg2rad(th), np.deg2rad(ph), indexing="ij")
    out = dict(theta=th.tolist(), cuts=[0, 45, 90], f=F_PAT, rows=[])
    for i in range(len(F_PAT)):
        Et, Ep = np.asarray(r.E_theta[i]), np.asarray(r.E_phi[i])
        co = Et * np.cos(PH) - Ep * np.sin(PH)
        cx = Et * np.sin(PH) + Ep * np.cos(PH)
        pw = np.abs(Et) ** 2 + np.abs(Ep) ** 2
        prad, pinc, pacc = float(r.Prad[i]), float(port.P_inc[i]), float(port.P_acc[i])
        k = float(r.Dmax[i]) / pw.max() * prad / pinc          # |E|^2 -> realised gain
        # cross-check of Dmax: integrate |E|^2 over the sphere
        w = np.sin(TH) * np.deg2rad(2.0) * np.deg2rad(5.0)
        d_int = 4 * np.pi * pw.max() / float((pw * w).sum())
        imax = np.unravel_index(np.argmax(pw), pw.shape)
        row = dict(f=F_PAT[i], Prad=prad, Pinc=pinc, Pacc=pacc, Dmax_dBi=10 * np.log10(float(r.Dmax[i])),
                   Dmax_int_dBi=10 * np.log10(d_int), peak_theta=float(th[imax[0]]), peak_phi=float(ph[imax[1]]),
                   peak_gain_dBi=10 * np.log10(k * pw.max()), cuts={})
        for c in out["cuts"]:
            j1, j2 = int(c // 5), int((c + 180) // 5)
            n = int(np.searchsorted(th, 90.0)) + 1           # theta 0..90
            tt = np.concatenate([-th[1:n][::-1], th[:n]])
            gco = np.concatenate([np.abs(co[1:n, j2][::-1]) ** 2, np.abs(co[:n, j1]) ** 2]) * k
            gcx = np.concatenate([np.abs(cx[1:n, j2][::-1]) ** 2, np.abs(cx[:n, j1]) ** 2]) * k
            row["cuts"][str(c)] = dict(theta=tt.tolist(), co_dBi=np.round(10 * np.log10(gco + 1e-30), 2).tolist(),
                                       cx_dBi=np.round(10 * np.log10(gcx + 1e-30), 2).tolist())
        out["rows"].append(row)
    return out


def worker_main(specfile):
    spec = read_json(specfile)
    parent = os.getppid()

    def watchdog():                      # die with the parent, so a killed run leaves no orphan solver
        while True:
            time.sleep(2)
            if os.getppid() != parent:
                os._exit(3)
    threading.Thread(target=watchdog, daemon=True).start()
    os.environ["MPLBACKEND"] = "Agg"
    res = run_openems(spec)
    write_json(os.path.join(spec["dir"], "result.json"), res)


# ----------------------------------------------------------------------------- cache and evaluation
class Cache:
    def __init__(self):
        os.makedirs(OUT, exist_ok=True)
        self.by_key, self.by_geom, self.order = {}, {}, []
        if os.path.exists(EVALS):
            good, bad = [], 0
            with open(EVALS) as f:
                for line in f:
                    try:
                        rec = json.loads(line)
                        rec["key"]
                    except (ValueError, KeyError, TypeError):
                        bad += 1
                        continue
                    good.append(line if line.endswith("\n") else line + "\n")
                    self._index(rec)
            if bad:                                           # a torn line from a kill: rewrite the file without it
                atomic_write(EVALS, "".join(good))
                log("cache: dropped %d corrupt line(s) from evals.jsonl" % bad)

    def _index(self, rec):
        if rec["key"] not in self.by_key:
            self.order.append(rec["key"])
        self.by_key[rec["key"]] = rec
        if rec.get("ok"):
            self.by_geom.setdefault(rec["geom"], rec)

    def put(self, rec):
        line = json.dumps(rec, default=_js, separators=(",", ":")) + "\n"
        fd = os.open(EVALS, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o644)
        try:
            os.write(fd, line.encode())                       # one write call: the line lands whole or not at all
            os.fsync(fd)
        finally:
            os.close(fd)
        self._index(rec)


_children = set()


def _simulate(spec, timeout):
    wd = spec["dir"]
    shutil.rmtree(wd, ignore_errors=True)
    os.makedirs(wd)
    sf = os.path.join(wd, "spec.json")
    write_json(sf, spec)
    err = ""
    for attempt in (1, 2):
        with open(os.path.join(wd, "solver.log"), "w") as lf:
            pr = subprocess.Popen([sys.executable, "-u", os.path.abspath(__file__), "_worker", sf],
                                  stdout=lf, stderr=subprocess.STDOUT, env=dict(os.environ, MPLBACKEND="Agg"))
            _children.add(pr)
            try:
                rc = pr.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                pr.kill()
                pr.wait()
                rc = -9
            _children.discard(pr)
        res = read_json(os.path.join(wd, "result.json"))
        if rc == 0 and res is not None and np.all(np.isfinite(res["S_re"])) and np.all(np.isfinite(res["S_im"])):
            try:
                with open(os.path.join(wd, "solver.log")) as lf:
                    txt = lf.read()
                ts = re.findall(r"Timestep:\s+(\d+).*?\(\s*(-?[\d.]+)dB\)", txt)
                if ts:
                    res["nts"], res["end_dB"] = int(ts[-1][0]), float(ts[-1][1])
            except OSError:
                pass
            shutil.rmtree(wd, ignore_errors=True)
            return res, ""
        try:
            with open(os.path.join(wd, "solver.log")) as lf:
                err = "rc=%s: %s" % (rc, lf.read()[-400:])
        except OSError:
            err = "rc=%s" % rc
        log("  run failed (attempt %d) in %s: %s" % (attempt, wd, err.replace("\n", " | ")[-300:]))
    return None, err                                          # the scratch folder is kept for inspection


def evaluate(cache, items, tag, jobs=None, threads=None, timeout=3600):
    """items: list of (model, x, eps, mesh). Returns their records, from the cache where possible."""
    out, todo = [None] * len(items), {}
    for i, (model, x, eps, m) in enumerate(items):
        x = snap(x)
        key = eval_key(model, x, eps, m)
        rec = cache.by_key.get(key)
        if rec is None:
            gh = geom_hash(model, x, eps, m)
            twin = cache.by_geom.get(gh)
            if twin is not None:                              # a different parameter vector, the same staircase
                rec = dict(twin, key=key, x=x.tolist(), tag=tag, t=time.time(), wall=0.0, alias_of=twin["key"])
                cache.put(rec)
                log("  [%s] %s alias of %s (identical mesh geometry)" % (tag, key, twin["key"]))
            elif key in todo:
                todo[key][1].append(i)
                continue
            else:
                todo[key] = (dict(model=model, x=x.tolist(), eps=float(eps), mesh=m, geom=gh, key=key), [i])
                continue
        else:
            log("  [%s] %s cached  score %.2f" % (tag, key, rec.get("score", float("nan"))))
        out[i] = rec
    if todo:
        nj = jobs or JOBS

        def one(job):
            spec = dict(job, dir=os.path.join(SCRATCH, "%s_%s" % (job["model"], job["key"])),
                        threads=threads or THREADS)
            res, err = _simulate(spec, timeout)
            rec = dict(key=job["key"], geom=job["geom"], model=job["model"], mesh=job["mesh"], x=job["x"],
                       eps=job["eps"], tag=tag, t=time.time(), ok=res is not None)
            if res is None:
                rec.update(err=err[-300:], score=1e3)
            else:
                rec.update(res)
                rec.update(metrics(rec))
            return rec
        with ThreadPoolExecutor(max_workers=nj) as ex:
            futs = {ex.submit(one, job): idx for job, idx in todo.values()}
            for fu in as_completed(futs):
                rec = fu.result()
                cache.put(rec)
                log("  [%s] %s ran %.0f s  score %.2f  worst %.2f dB  x=%s" % (
                    tag, rec["key"], rec.get("wall", 0), rec["score"], rec.get("worst", float("nan")),
                    " ".join("%g" % v for v in rec["x"])))
                for i in futs[fu]:
                    out[i] = rec
    return out


def _terminate(signum, frame):
    for pr in list(_children):
        try:
            pr.kill()
        except OSError:
            pass
    log("stopped by signal %d" % signum)
    os._exit(1)


# ----------------------------------------------------------------------------- post-processing
def S_of(rec):
    return np.array(rec["S_re"]) + 1j * np.array(rec["S_im"])


def gamma_act(rec, theta_deg, phi_deg=0.0, active=None):
    """Active reflection of the excited element against frequency for a beam steered to (theta, phi)."""
    S = S_of(rec)
    pos, c = port_positions(rec["model"])
    act = active_set(rec["model"]) if active is None else active
    th, ph = np.deg2rad(theta_deg), np.deg2rad(phi_deg)
    u = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph)])
    k = 2 * np.pi * FREQ / C0
    d = (pos[act] - pos[c]) @ u * 1e-3
    return np.sum(S[act] * np.exp(-1j * k[None, :] * d[:, None]), axis=0)


def scan_table(rec, planes=((0.0, "E"),)):
    """|Gamma_act| in dB, shape (planes, scans, freq)."""
    return np.array([[db(gamma_act(rec, t, ph)) for t in SCAN] for ph, _ in planes])


def metrics(rec, planes=((0.0, "E"),)):
    g = scan_table(rec, planes)
    S = S_of(rec)
    pen = penalty(rec["x"])
    worst, mean = float(g.max()), float(g.mean())
    iw = np.unravel_index(np.argmax(g), g.shape)
    return dict(worst=worst, mean=mean, pen_mm=pen, score=worst + 0.3 * mean + 100.0 * pen,
                worst_f=float(FREQ[iw[2]]), worst_scan=float(SCAN[iw[1]]),
                smax=float(np.abs(S).max()), colsum=float((np.abs(S) ** 2).sum(axis=0).max()))


def score_planes(rec, both):
    pl = ((0.0, "E"), (90.0, "H")) if both else ((0.0, "E"),)
    g = scan_table(rec, pl)
    return float(g.max() + 0.3 * g.mean() + 100.0 * penalty(rec["x"])), float(g.max()), float(g.mean())


# ----------------------------------------------------------------------------- stages
def stage_validate(cache):
    log("validate: starting design on the row model")
    spec = lambda m: ("row", X0, EPS_NOM, m)
    # coarse first (also the timing reference), then the expensive ones one at a time with more threads
    base = evaluate(cache, [spec(M_BASE)], "validate")[0]
    others = evaluate(cache, [spec(M_GAP), spec(M_END), spec(mesh_cfg(engine=ALT_ENGINE))], "validate")
    f1 = evaluate(cache, [spec(M_FINE1)], "validate", jobs=1, threads=TILE_THREADS)[0]
    f2 = evaluate(cache, [spec(M_FINE2)], "validate", jobs=1, threads=TILE_THREADS, timeout=4 * 3600)[0]
    gapr, endr, engr = others
    crs = evaluate(cache, [spec(M_COARSE)], "validate")[0]
    runs = dict(coarse=crs, base=base, fine1=f1, fine2=f2, gap50=gapr, end1e5=endr, engine=engr)
    if not all(r.get("ok") for r in runs.values()):
        v = dict(complete=True, trustworthy=False, reason="a validation run failed", keys={k: r["key"] for k, r in runs.items()})
        write_json(os.path.join(OUT, "validate.json"), v)
        return v

    def dg(a, b):
        return float(np.max(np.abs(10 ** (scan_table(a) / 20) - 10 ** (scan_table(b) / 20))))
    num = dict(
        score={k: r["score"] for k, r in runs.items()}, worst={k: r["worst"] for k, r in runs.items()},
        mean={k: r["mean"] for k, r in runs.items()},
        wall={k: r.get("wall") for k, r in runs.items()}, cells={k: r.get("cells") for k, r in runs.items()},
        nts={k: r.get("nts") for k, r in runs.items()},
        d_score_base_fine1=base["score"] - f1["score"], d_score_fine1_fine2=f1["score"] - f2["score"],
        d_score_base_fine2=base["score"] - f2["score"],
        dG_base_fine1=dg(base, f1), dG_fine1_fine2=dg(f1, f2), dG_base_fine2=dg(base, f2),
        d_score_gap=base["score"] - gapr["score"], dG_gap=dg(base, gapr),
        d_score_end=base["score"] - endr["score"], dG_end=dg(base, endr),
        dS_engine=float(np.max(np.abs(S_of(base) - S_of(engr)))), alt_engine=ALT_ENGINE, engine=ENGINE,
        smax=max(r["smax"] for r in runs.values()), colsum=max(r["colsum"] for r in runs.values()),
        passive_inc=max(r.get("passive_inc", 0) for r in runs.values()))
    S = S_of(f1)
    adj = db(S[[1, 3]])
    num["lit"] = dict(coupling_adjacent_dB=[float(adj.min()), float(adj.max())],
                      coupling_second_dB=[float(db(S[[0, 4]]).min()), float(db(S[[0, 4]]).max())],
                      vswr_broadside=_vswr_stats(gamma_act(f1, 0.0)), vswr_45=_vswr_stats(gamma_act(f1, 45.0)),
                      s11_dB=[float(db(S[2]).min()), float(db(S[2]).max())],
                      absorbed_or_radiated=[float(1 - (np.abs(S) ** 2).sum(axis=0).max()),
                                            float(1 - (np.abs(S) ** 2).sum(axis=0).min())])
    gates = dict(
        mesh=abs(num["d_score_fine1_fine2"]) <= TOL_MESH_DB and num["dG_fine1_fine2"] <= TOL_MESH_G,
        boundary=abs(num["d_score_gap"]) <= TOL_BND_DB and num["dG_gap"] <= TOL_BND_G,
        passive=num["smax"] <= TOL_PASSIVE and num["colsum"] <= TOL_PASSIVE,
        engine=num["dS_engine"] <= TOL_ENGINE,
        end=abs(num["d_score_end"]) <= TOL_BND_DB and num["dG_end"] <= TOL_BND_G)
    v = dict(complete=True, trustworthy=all(gates.values()), gates=gates, numbers=num,
             keys={k: r["key"] for k, r in runs.items()},
             tolerances=dict(mesh_dB=TOL_MESH_DB, mesh_G=TOL_MESH_G, bnd_dB=TOL_BND_DB, bnd_G=TOL_BND_G,
                             passive=TOL_PASSIVE, engine=TOL_ENGINE))
    write_json(os.path.join(OUT, "validate.json"), v)
    log("validate: gates %s -> trustworthy=%s" % (gates, v["trustworthy"]))
    return v


def _vswr_stats(g):
    a = np.minimum(np.abs(g), 0.999)
    v = (1 + a) / (1 - a)
    return dict(max=float(v.max()), median=float(np.median(v)), frac_below_2=float((v <= 2).mean()),
                frac_below_3=float((v <= 3).mean()), at_2p3=float(v[0]), at_6=float(v[-1]))


def stage_search(cache):
    from array16.evolve import climb, evolve
    seen = {}

    def score(X, gen):
        recs = evaluate(cache, [("row", x, EPS_NOM, M_BASE) for x in X], "search")
        for r in recs:
            seen[r["key"]] = r
        return np.array([r["score"] for r in recs])

    def progress(phase, hist):
        write_json(os.path.join(OUT, "search_progress.json"), dict(phase=phase, t=time.time(), history=hist))

    log("search: climb from the starting design (%d steps x %d)" % (CLIMB_STEPS, CLIMB_LAM))
    hist1 = []
    r1 = climb(score, X0, LO, HI, steps=CLIMB_STEPS, lam=CLIMB_LAM, sigma=CLIMB_SIGMA, seed=CLIMB_SEED, repair=snap,
               callback=lambda s, x, b: (hist1.append(dict(step=s, best=b, x=list(x))), progress("climb", hist1),
                                         log("search: climb step %d best %.2f" % (s, b))))
    log("search: evolve, population %d, %d generations, seeded with the climb result and the start" % (EVO_POP, EVO_GENS))
    hist2 = []
    r2 = evolve(score, LO, HI, pop=EVO_POP, generations=EVO_GENS, elite=EVO_ELITE, sigma=EVO_SIGMA, seed=EVO_SEED,
                init=np.array([r1.x, X0]), repair=snap,
                callback=lambda g, x, b: (hist2.append(dict(gen=g, best=b, x=list(x))), progress("evolve", hist2),
                                          log("search: evolve generation %d best %.2f" % (g, b))))
    # finalists: the best distinct designs seen, re-scored at both ends of the eps_r range
    ranked = sorted((r for r in seen.values() if r.get("ok")), key=lambda r: r["score"])
    fin = []
    for r in ranked:
        if all(np.max(np.abs((np.array(r["x"]) - np.array(q["x"])) / (HI - LO))) > 0.03 for q in fin):
            fin.append(r)
        if len(fin) == N_FINAL:
            break
    start = seen[eval_key("row", X0, EPS_NOM, M_BASE)]
    cands = fin + ([start] if all(q["key"] != start["key"] for q in fin) else [])
    log("search: re-scoring %d finalists and the start at eps_r %.1f and %.1f" % (len(fin), EPS_LO, EPS_HI))
    items = [("row", r["x"], e, M_BASE) for r in cands for e in (EPS_LO, EPS_HI)]
    rr = evaluate(cache, items, "finalist")
    table = []
    for i, r in enumerate(cands):
        s = [rr[2 * i]["score"], r["score"], rr[2 * i + 1]["score"]]
        table.append(dict(x=r["x"], key=r["key"], score_eps=dict(zip(["4.2", "4.4", "4.6"], s)), robust=max(s),
                          is_start=r["key"] == start["key"]))
    table.sort(key=lambda d: d["robust"])
    win = next(d for d in table if True)
    res = dict(complete=True, x0=X0.tolist(), winner=win["x"], winner_is_start=win["is_start"], finalists=table,
               climb=dict(x=r1.x.tolist(), score=r1.score, history=r1.history),
               evolve=dict(x=r2.x.tolist(), score=r2.score, history=r2.history),
               settings=dict(climb=[CLIMB_STEPS, CLIMB_LAM, CLIMB_SIGMA, CLIMB_SEED],
                             evolve=[EVO_POP, EVO_GENS, EVO_ELITE, EVO_SIGMA, EVO_SEED]))
    write_json(os.path.join(OUT, "search.json"), res)
    log("search: winner %s robust score %.2f (start %.2f)" % (
        win["x"], win["robust"], next(d["robust"] for d in table if d["is_start"])))
    return res


def stage_check(cache):
    s = read_json(os.path.join(OUT, "search.json"))
    xw = np.array(s["winner"])
    out = dict(complete=False, x0=X0.tolist(), winner=xw.tolist(), keys={})

    def save():
        write_json(os.path.join(OUT, "check.json"), out)
    log("check: winner on the finer row mesh")
    for name, x in (("start", X0), ("winner", xw)):
        r = evaluate(cache, [("row", x, EPS_NOM, M_FINE1)], "check", jobs=1, threads=TILE_THREADS)[0]
        out["keys"]["row_fine_%s" % name] = r["key"]
        save()
    for eps in (EPS_LO, EPS_HI, EPS_NOM):
        for name, x in (("start", X0), ("winner", xw)):
            log("check: 6 x 6 tile, %s, eps_r %.1f" % (name, eps))
            r = evaluate(cache, [("tile", x, eps, M_TILE)], "check", jobs=1, threads=TILE_THREADS, timeout=6 * 3600)[0]
            out["keys"]["tile_%s_%.1f" % (name, eps)] = r["key"]
            save()
    out["complete"] = True
    save()
    return out


# ----------------------------------------------------------------------------- report
def _fmt_x(x):
    return ", ".join("%s %g" % (n, v) for n, v in zip(NAMES, x))


def outline(p, n=300):
    """Outline polygons of one cell for drawing: (front xs, zs), (back xs, zs)."""
    zs = np.linspace(0.0, p[0], n)
    f = np.array([spans(p, min(z, p[0] - 1e-6))[0] for z in zs])
    b = np.array([spans(p, min(z, p[0] - 1e-6))[1] for z in zs])
    return zs, f, b


def stage_report(cache=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cache = cache or Cache()
    K = cache.by_key
    val = read_json(os.path.join(OUT, "validate.json"))
    sea = read_json(os.path.join(OUT, "search.json"))
    chk = read_json(os.path.join(OUT, "check.json"))
    fg = FREQ / 1e9
    L = ["# Connected antipodal Vivaldi, openEMS: machine-written summary", "",
         "Written by `experiments/vivaldi_openems.py report` on %s from the files in this folder. "
         "%d evaluations in `evals.jsonl`." % (time.strftime("%Y-%m-%d %H:%M"), len(K)), "",
         "Stages complete: validate %s, search %s, check %s." % tuple(
             "yes" if d and d.get("complete") else "NO (pending)" for d in (val, sea, chk)), "",
         "Score = worst |Gamma_act| (dB) over 2.3-6 GHz and E-plane scan 0, +/-30, +/-45 deg + 0.3 x mean (dB) "
         "+ manufacturing penalty. Lower is better. Reference: 50 ohm.", ""]

    def savefig(fig, name):
        tmp = os.path.join(OUT, "." + name + ".tmp.png")
        fig.savefig(tmp, dpi=130)
        plt.close(fig)
        os.replace(tmp, os.path.join(OUT, name))

    # ---- parameters
    L += ["## Parameters", "", "| Name | Meaning | Unit | Low | High | Start |" + (" Winner |" if sea else ""),
          "| --- | --- | --- | --- | --- | --- |" + (" --- |" if sea else "")]
    for i, n in enumerate(NAMES):
        L.append("| %s | %s | %s | %g | %g | %g |" % (n, WHAT[i], UNITS[i], LO[i], HI[i], X0[i])
                 + (" %g |" % sea["winner"][i] if sea else ""))
    L.append("")

    # ---- validation
    if val and val.get("numbers"):
        n, kk = val["numbers"], val["keys"]
        L += ["## Validation (starting design, row model, eps_r %.1f)" % EPS_NOM, "",
              "| Run | Mesh along card and z (mm) | Cells through laminate | Air to PML (mm) | Cells | Time steps | Wall (s) | Worst (dB) | Mean (dB) | Score |",
              "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
        for name in ("coarse", "base", "fine1", "fine2", "gap50", "end1e5", "engine"):
            r = K.get(kk[name])
            if r:
                m = r["mesh"]
                L.append("| %s | %.3f | %d | %g | %s | %s | %s | %.2f | %.2f | %.2f |" % (
                    name + (" (%s)" % m["engine"] if name == "engine" else "") + (" (end 1e-5)" if name == "end1e5" else ""),
                    P / 2 / m["nh"], m["nsub"], m["gap"], r.get("cells"), r.get("nts"), r.get("wall"),
                    r["worst"], r["mean"], r["score"]))
        g = val["gates"]
        t = val["tolerances"]
        yn = lambda b: "pass" if b else "**FAIL**"
        L += ["", "| Check | Number | Limit | Result |", "| --- | --- | --- | --- |",
              "| Mesh: score, fine1 to fine2 | %+.2f dB | %.1f dB | %s |" % (n["d_score_fine1_fine2"], t["mesh_dB"], yn(g["mesh"])),
              "| Mesh: largest change of abs(Gamma_act), fine1 to fine2 | %.3f | %.2f | %s |" % (n["dG_fine1_fine2"], t["mesh_G"], yn(g["mesh"])),
              "| Mesh: score, base to fine2 (the error of the search mesh) | %+.2f dB | reported | - |" % n["d_score_base_fine2"],
              "| Mesh: largest change of abs(Gamma_act), base to fine2 | %.3f | reported | - |" % n["dG_base_fine2"],
              "| Boundary: score, PML 30 to 50 mm away | %+.2f dB | %.1f dB | %s |" % (n["d_score_gap"], t["bnd_dB"], yn(g["boundary"])),
              "| Boundary: largest change of abs(Gamma_act) | %.3f | %.2f | %s |" % (n["dG_gap"], t["bnd_G"], yn(g["boundary"])),
              "| End criterion 1e-4 to 1e-5: score | %+.2f dB | %.1f dB | %s |" % (n["d_score_end"], t["bnd_dB"], yn(g["end"])),
              "| End criterion: largest change of abs(Gamma_act) | %.3f | %.2f | %s |" % (n["dG_end"], t["bnd_G"], yn(g["end"])),
              "| Passivity: largest abs(S) | %.4f | %.3f | %s |" % (n["smax"], t["passive"], yn(g["passive"])),
              "| Passivity: largest column power sum | %.4f | %.3f | %s |" % (n["colsum"], t["passive"], yn(g["passive"])),
              "| Engine %s against %s: largest abs(dS) | %.4f | %.2f | %s |" % (n["engine"], n["alt_engine"], n["dS_engine"], t["engine"], yn(g["engine"])),
              "| Incident wave at the passive ports (should be 0) | %.1e | - | - |" % n["passive_inc"],
              "", "**Model trustworthy by these gates: %s.**" % ("yes" if val["trustworthy"] else "NO"), ""]
        lit = n["lit"]
        L += ["Behaviour of the starting design on the fine1 mesh, for comparison with published arrays "
              "(sources and discussion in `../vivaldi_study.md`):", "",
              "| Quantity | Model |", "| --- | --- |",
              "| Coupling to the adjacent cells, 2.3-6 GHz | %.1f to %.1f dB |" % tuple(lit["coupling_adjacent_dB"]),
              "| Coupling to the cells two away | %.1f to %.1f dB |" % tuple(lit["coupling_second_dB"]),
              "| Passive S11 of the centre cell | %.1f to %.1f dB |" % tuple(lit["s11_dB"]),
              "| Active VSWR at broadside: max, median, share of band below 2, below 3 | %.2f, %.2f, %.0f %%, %.0f %% |" % (
                  lit["vswr_broadside"]["max"], lit["vswr_broadside"]["median"], 100 * lit["vswr_broadside"]["frac_below_2"], 100 * lit["vswr_broadside"]["frac_below_3"]),
              "| Active VSWR at 45 deg E-plane: max, median, share below 2, below 3 | %.2f, %.2f, %.0f %%, %.0f %% |" % (
                  lit["vswr_45"]["max"], lit["vswr_45"]["median"], 100 * lit["vswr_45"]["frac_below_2"], 100 * lit["vswr_45"]["frac_below_3"]),
              "| Active VSWR at broadside at 2.3 and 6 GHz | %.2f, %.2f |" % (lit["vswr_broadside"]["at_2p3"], lit["vswr_broadside"]["at_6"]),
              "| Share of the incident power radiated or lost in FR-4 (1 - column sum) | %.2f to %.2f |" % tuple(lit["absorbed_or_radiated"]), ""]
        fig, ax = plt.subplots(1, 2, figsize=(11, 4), tight_layout=True)
        for name, st in (("coarse", "-."), ("base", "-"), ("fine1", "-"), ("fine2", "-"), ("gap50", "--"), ("end1e5", ":"), ("engine", ":")):
            r = K.get(kk[name])
            if r:
                ax[0].plot(fg, db(gamma_act(r, 0)), st, label=name)
                ax[1].plot(fg, db(gamma_act(r, 45)), st, label=name)
        for a, t_ in zip(ax, ("broadside", "45 deg, E-plane")):
            a.set(xlabel="GHz", ylabel="|Gamma_act| (dB)", title="Start design, row model: " + t_)
            a.grid(alpha=.3)
            a.legend(fontsize=7)
        savefig(fig, "validation.png")
        L += ["![validation](validation.png)", ""]
    elif val:
        L += ["## Validation", "", "Validation did not finish: %s" % val.get("reason"), ""]

    # ---- search
    srch = [K[k] for k in cache.order if K[k].get("tag") == "search" and K[k].get("ok")]
    r0 = K.get(eval_key("row", X0, EPS_NOM, M_BASE))
    if srch:
        best = np.minimum.accumulate([r["score"] for r in srch])
        L += ["## Search (row model, base mesh, eps_r %.1f)" % EPS_NOM, "",
              "%d distinct designs evaluated so far (planned: about %d including finalists). Best score %.2f against %.2f for the start." % (
                  len(srch), PLAN_SEARCH, best[-1], r0["score"] if r0 else float("nan")), ""]
        fig, ax = plt.subplots(figsize=(7, 4), tight_layout=True)
        ax.plot(np.arange(1, len(srch) + 1), [r["score"] for r in srch], ".", ms=3, color="0.6", label="each evaluation")
        ax.plot(np.arange(1, len(srch) + 1), best, "-", label="best so far")
        if r0:
            ax.axhline(r0["score"], ls="--", color="k", lw=.8, label="start")
        ax.set(xlabel="new evaluations (cache order)", ylabel="score (dB)", title="Convergence")
        ax.set_ylim(best[-1] - 1, (r0["score"] if r0 else best[0]) + 4)
        ax.grid(alpha=.3)
        ax.legend()
        savefig(fig, "convergence.png")
        L += ["![convergence](convergence.png)", "", "| Evaluations | Best score (dB) |", "| --- | --- |"]
        for n_ in sorted(set([1, 10, 25, 50, 75, 100, 150, 200, 250, 300, len(srch)])):
            if n_ <= len(srch):
                L.append("| %d | %.2f |" % (n_, best[n_ - 1]))
        L.append("")
    if sea:
        L += ["### Finalists, re-scored at both ends of the eps_r range (row model, base mesh)", "",
              "| Design | " + " | ".join(NAMES) + " | Score eps 4.2 | 4.4 | 4.6 | Robust (worst of the three) |",
              "| --- |" + " --- |" * (len(NAMES) + 4)]
        for i, d in enumerate(sea["finalists"]):
            L.append("| %s | %s | %.2f | %.2f | %.2f | %.2f |" % (
                "start" if d["is_start"] else ("winner" if d["x"] == sea["winner"] else "finalist %d" % (i + 1)),
                " | ".join("%g" % v for v in d["x"]), d["score_eps"]["4.2"], d["score_eps"]["4.4"], d["score_eps"]["4.6"], d["robust"]))
        L += ["", "Winner = best robust score. " + ("**The start design was not beaten.**" if sea["winner_is_start"] else ""), "",
              "| Parameter | Start | Winner | Change |", "| --- | --- | --- | --- |"]
        for i, n_ in enumerate(NAMES):
            L.append("| %s (%s) | %g | %g | %+g |" % (n_, UNITS[i], X0[i], sea["winner"][i], round(sea["winner"][i] - X0[i], 4)))
        L.append("")
        rw = K.get(eval_key("row", sea["winner"], EPS_NOM, M_BASE))
        if r0 and rw:
            fig, ax = plt.subplots(1, 3, figsize=(14, 4), tight_layout=True, sharey=True)
            for a, th in zip(ax, (0, 30, 45)):
                for r, lab, col in ((r0, "start", "C0"), (rw, "winner", "C3")):
                    a.plot(fg, db(gamma_act(r, th)), color=col, label="%s +%d" % (lab, th))
                    if th:
                        a.plot(fg, db(gamma_act(r, -th)), "--", color=col, label="%s -%d" % (lab, th))
                a.axhline(-10, color="k", lw=.5)
                a.set(xlabel="GHz", title="row model, E-plane scan %d deg" % th)
                a.grid(alpha=.3)
                a.legend(fontsize=7)
            ax[0].set_ylabel("|Gamma_act| (dB)")
            savefig(fig, "active_reflection_row.png")
            L += ["### Active reflection per frequency, row model (dB)", "", "![row](active_reflection_row.png)", "",
                  "| GHz | start 0 | start +30 | start -30 | start +45 | start -45 | winner 0 | winner +30 | winner -30 | winner +45 | winner -45 |",
                  "| --- |" + " --- |" * 10]
            g0, gw = scan_table(r0)[0], scan_table(rw)[0]
            for i in range(0, len(FREQ), 2):
                L.append("| %.1f | " % fg[i] + " | ".join("%.1f" % v for v in list(g0[:, i]) + list(gw[:, i])) + " |")
            L.append("")
        # outlines
        fig, ax = plt.subplots(1, 2, figsize=(10, 5), tight_layout=True)
        for a, (p, name) in zip(ax, ((X0, "start"), (np.array(sea["winner"]), "winner"))):
            zs, f, b = outline(p)
            for c in (-1, 0, 1):
                a.fill_betweenx(zs, f[:, 0] + c * P, f[:, 1] + c * P, color="C1", alpha=.6, lw=0, label="front copper" if c == 0 else None)
                a.fill_betweenx(zs, b[:, 0] + c * P, b[:, 1] + c * P, color="C0", alpha=.45, lw=0, label="back copper" if c == 0 else None)
            zj = p[3] + p[5] * (p[0] - p[3])
            for xb in (-P / 2, P / 2):
                a.plot([xb, xb], [zj, p[0]], "k-", lw=2, label="via row" if xb < 0 else None)
            a.set(xlim=(-1.5 * P, 1.5 * P), ylim=(0, 46), xlabel="x (mm)", ylabel="z (mm)", title="%s: %s" % (name, _fmt_x(p)))
            a.title.set_fontsize(7)
            a.set_aspect("equal")
            a.legend(fontsize=7, loc="lower right")
        savefig(fig, "outline.png")
        L += ["![outline](outline.png)", ""]

    # ---- check
    if chk:
        kk = chk["keys"]
        L += ["## Check%s" % ("" if chk.get("complete") else " (INCOMPLETE: %d of 8 runs done)" % len(kk)), "",
              "Tile model: finite 6 x 6, mesh %.3f mm, element (card %d, cell %d) excited, inner 4 x 4 driven for "
              "the active reflection, scan 0, +/-30, +/-45 deg." % (P / 2 / M_TILE["nh"], EXC_CARD, EXC_CELL), "",
              "| Model | Design | eps_r | Worst E-plane (dB) | Mean E (dB) | Score E | Worst E+H (dB) | Score E+H | Worst at |",
              "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
        rows = {}
        for lab, mk in [("row, base mesh", lambda d, e: eval_key("row", chk["x0"] if d == "start" else chk["winner"], e, M_BASE)),
                        ("row, fine1 mesh", lambda d, e: kk.get("row_fine_%s" % d) if e == EPS_NOM else None),
                        ("tile 6 x 6", lambda d, e: kk.get("tile_%s_%.1f" % (d, e)))]:
            for e in (EPS_LO, EPS_NOM, EPS_HI):
                for d in ("start", "winner"):
                    r = K.get(mk(d, e) or "")
                    if not r or not r.get("ok"):
                        continue
                    sE, wE, mE = score_planes(r, False)
                    if r["model"] == "tile":
                        sB, wB, _ = score_planes(r, True)
                        extra = "%.2f | %.2f" % (wB, sB)
                    else:
                        extra = "- | -"
                    mt = metrics(r)
                    rows[(lab, d, e)] = (sE, wE, sB if r["model"] == "tile" else None)
                    L.append("| %s | %s | %.1f | %.2f | %.2f | %.2f | %s | %.1f GHz, %+.0f deg |" % (
                        lab, d, e, wE, mE, sE, extra, mt["worst_f"] / 1e9, mt["worst_scan"]))
        L.append("")
        # verdict
        gains = {}
        for lab in ("row, base mesh", "row, fine1 mesh", "tile 6 x 6"):
            for idx, nm in ((0, "E"), (2, "E+H")):
                pairs = [(rows[(lab, "start", e)][idx], rows[(lab, "winner", e)][idx]) for e in (EPS_LO, EPS_NOM, EPS_HI)
                         if (lab, "start", e) in rows and (lab, "winner", e) in rows and rows[(lab, "start", e)][idx] is not None]
                if pairs:
                    gains[(lab, nm)] = (max(p_[0] for p_ in pairs), max(p_[1] for p_ in pairs), len(pairs))
        L += ["| Model | Planes | Start, worst score over eps_r | Winner, worst score over eps_r | Improvement (dB) | eps_r values |",
              "| --- | --- | --- | --- | --- | --- |"]
        for (lab, nm), (a, b, n_) in gains.items():
            L.append("| %s | %s | %.2f | %.2f | %+.2f | %d |" % (lab, nm, a, b, a - b, n_))
        L.append("")
        mesh_err = abs(val["numbers"]["d_score_base_fine2"]) if val and val.get("numbers") else float("nan")
        tg = gains.get(("tile 6 x 6", "E"))
        tb = gains.get(("tile 6 x 6", "E+H"))
        fgain = gains.get(("row, fine1 mesh", "E"))
        if tg and chk.get("complete"):
            imp, imp2 = tg[0] - tg[1], tb[0] - tb[1]
            impf = fgain[0] - fgain[1] if fgain else float("nan")
            if imp > mesh_err and imp2 > 0 and impf > 0:
                verdict = ("**The winner holds on the check model**: %.2f dB better than the start on the tile (E-plane score, worst over eps_r), "
                           "%.2f dB with the H-plane included, %.2f dB on the finer row mesh; the search-mesh error is %.2f dB." % (imp, imp2, impf, mesh_err))
            elif imp > 0 and imp2 > 0 and impf > 0:
                verdict = ("**The winner is better on every model, but by less than the mesh error on the tile**: %.2f dB on the tile (E-plane), %.2f dB with "
                           "the H-plane, %.2f dB on the finer row mesh, against a search-mesh error of %.2f dB. Treat it as not proven." % (imp, imp2, impf, mesh_err))
            else:
                verdict = ("**The winner only wins on the search model: it is not a winner.** Tile improvement %.2f dB (E-plane), %.2f dB (E and H), "
                           "finer row mesh %.2f dB; search-mesh error %.2f dB. Keep the starting design, or re-run the search with the tile in the loop." % (imp, imp2, impf, mesh_err))
            L += [verdict, ""]
        else:
            L += ["Verdict pending: the check stage has not finished.", ""]
        # tile figures and tables
        tiles = [(d, e, K.get(kk.get("tile_%s_%.1f" % (d, e), ""))) for e in (EPS_LO, EPS_NOM, EPS_HI) for d in ("start", "winner")]
        tiles = [t_ for t_ in tiles if t_[2] and t_[2].get("ok")]
        if tiles:
            fig, ax = plt.subplots(2, 3, figsize=(14, 7.5), tight_layout=True, sharey=True)
            for row_, (ph, pn) in enumerate(((0.0, "E"), (90.0, "H"))):
                for a, th in zip(ax[row_], (0, 30, 45)):
                    for d, e, r in tiles:
                        col = "C0" if d == "start" else "C3"
                        ls = {EPS_LO: "--", EPS_NOM: "-", EPS_HI: ":"}[e]
                        g = np.maximum(db(gamma_act(r, th, ph)), db(gamma_act(r, -th, ph)))
                        a.plot(fg, g, ls, color=col, label="%s eps %.1f" % (d, e))
                    a.axhline(-10, color="k", lw=.5)
                    a.set(xlabel="GHz", title="tile, %s-plane scan %d deg (worse of +/-)" % (pn, th))
                    a.grid(alpha=.3)
                ax[row_][0].set_ylabel("|Gamma_act| (dB)")
            ax[0][0].legend(fontsize=7)
            savefig(fig, "active_reflection_tile.png")
            ths = np.arange(0, 61, 5)
            fig, ax = plt.subplots(1, 3, figsize=(14, 4), tight_layout=True, sharey=True)
            for a, f_ in zip(ax, F_PAT):
                i = int(np.argmin(np.abs(FREQ - f_)))
                for d, e, r in tiles:
                    if e != tiles[0][1]:
                        continue
                    for ph, ls in ((0.0, "-"), (90.0, "--")):
                        g = [max(db(gamma_act(r, t_, ph))[i], db(gamma_act(r, -t_, ph))[i]) for t_ in ths]
                        a.plot(ths, g, ls, color="C0" if d == "start" else "C3", label="%s %s-plane" % (d, "E" if ph == 0 else "H"))
                a.set(xlabel="scan angle (deg)", title="tile, %.1f GHz, eps_r %.1f" % (FREQ[i] / 1e9, tiles[0][1]))
                a.grid(alpha=.3)
            ax[0].set_ylabel("|Gamma_act| (dB)")
            ax[0].legend(fontsize=7)
            savefig(fig, "active_vs_scan_tile.png")
            L += ["![tile](active_reflection_tile.png)", "", "![scan](active_vs_scan_tile.png)", "",
                  "### Active reflection on the tile at the pattern frequencies (dB, worse of +/- scan)", "",
                  "| Design | eps_r | GHz | E 0 | E 30 | E 45 | H 30 | H 45 | passive S11 |", "| --- |" + " --- |" * 8]
            for d, e, r in tiles:
                S = S_of(r)
                c = port_positions("tile")[1]
                for f_ in F_PAT:
                    i = int(np.argmin(np.abs(FREQ - f_)))
                    v = [db(gamma_act(r, 0))[i]] + [max(db(gamma_act(r, t_, ph))[i], db(gamma_act(r, -t_, ph))[i])
                                                    for ph in (0.0, 90.0) for t_ in (30, 45)]
                    L.append("| %s | %.1f | %.1f | " % (d, e, FREQ[i] / 1e9) + " | ".join("%.1f" % q for q in v) + " | %.1f |" % db(S[c])[i])
            L.append("")
            far = [t_ for t_ in tiles if t_[2].get("far")]
            if far:
                L += ["### Embedded element pattern (centre element excited, all others terminated)", "",
                      "Realised gain = 4 pi r^2 S / P_incident: mismatch, power coupled into the other 35 ports and FR-4 loss all count against it. "
                      "Co-polar = x (Ludwig 3). Cross-polar levels are relative to the co-polar gain at the same angle; "
                      "'D' is the 45 deg cut. Efficiency = P_radiated / P_accepted.", "",
                      "| Design | eps_r | GHz | Realised gain, broadside (dBi) | Unit-cell limit 4 pi A / lambda^2 (dBi) | Peak (dBi) at theta, phi | Directivity (dBi), solver / integrated | Reflected | Coupled to other ports | FR-4 loss | Radiated / incident | Efficiency (rad/acc) | X-pol broadside (dB) | X-pol worst within 45 deg: E, H, D (dB) |",
                      "| --- |" + " --- |" * 13]
                for d, e, r in far:
                    S = S_of(r)
                    c = port_positions("tile")[1]
                    for row_ in r["far"]["rows"]:
                        i = int(np.argmin(np.abs(FREQ - row_["f"])))
                        refl = abs(S[c, i]) ** 2
                        coup = float((np.abs(S[:, i]) ** 2).sum() - refl)
                        radf = row_["Prad"] / row_["Pinc"]
                        xp = []
                        for cname in ("0", "90", "45"):
                            cu_ = row_["cuts"][cname]
                            th_ = np.array(cu_["theta"])
                            sel = np.abs(th_) <= 45
                            xp.append(float(np.max(np.array(cu_["cx_dBi"])[sel] - np.array(cu_["co_dBi"])[sel])))
                        c0 = row_["cuts"]["0"]
                        i0 = int(np.argmin(np.abs(np.array(c0["theta"]))))
                        L.append("| %s | %.1f | %.1f | %.1f | %.1f | %.1f at %.0f, %.0f | %.1f / %.1f | %.0f %% | %.0f %% | %.0f %% | %.0f %% | %.0f %% (%.1f dB) | %.1f | %.1f, %.1f, %.1f |" % (
                            d, e, row_["f"] / 1e9, c0["co_dBi"][i0], 10 * np.log10(4 * np.pi * (P * 1e-3 * row_["f"] / C0) ** 2), row_["peak_gain_dBi"], row_["peak_theta"], row_["peak_phi"],
                            row_["Dmax_dBi"], row_["Dmax_int_dBi"], 100 * refl, 100 * coup, 100 * max(1 - refl - coup - radf, 0), 100 * radf,
                            100 * row_["Prad"] / row_["Pacc"], 10 * np.log10(row_["Prad"] / row_["Pacc"]),
                            c0["cx_dBi"][i0] - c0["co_dBi"][i0], xp[0], xp[1], xp[2]))
                L.append("")
                e_show = far[0][1]
                fig, ax = plt.subplots(3, 3, figsize=(14, 10), tight_layout=True, sharey=True)
                for i, f_ in enumerate(F_PAT):
                    for j, cname in enumerate(("0", "90", "45")):
                        a = ax[i][j]
                        for d, e, r in far:
                            if e != e_show:
                                continue
                            cu_ = r["far"]["rows"][i]["cuts"][cname]
                            col = "C0" if d == "start" else "C3"
                            a.plot(cu_["theta"], cu_["co_dBi"], "-", color=col, label="%s co" % d)
                            a.plot(cu_["theta"], cu_["cx_dBi"], "--", color=col, label="%s cross" % d)
                        a.set(title="%.1f GHz, %s, eps_r %.1f" % (f_ / 1e9, {"0": "E-plane", "90": "H-plane", "45": "D-plane"}[cname], e_show),
                              xlabel="theta (deg)", ylim=(-30, 10))
                        a.grid(alpha=.3)
                    ax[i][0].set_ylabel("realised gain (dBi)")
                ax[0][0].legend(fontsize=7)
                savefig(fig, "embedded_pattern.png")
                L += ["![pattern](embedded_pattern.png)", ""]
    # ---- compute
    tot = sum(r.get("wall", 0) or 0 for r in K.values())
    L += ["## Compute", "", "| Tag | Runs | Solver wall time (h) |", "| --- | --- | --- |"]
    for tag in ("validate", "search", "finalist", "check"):
        rs = [r for r in K.values() if r.get("tag") == tag and not r.get("alias_of")]
        L.append("| %s | %d | %.2f |" % (tag, len(rs), sum(r.get("wall", 0) or 0 for r in rs) / 3600))
    L += ["| all | %d | %.2f |" % (len(K), tot / 3600), "",
          "Failed runs: %d. Staircase aliases (different parameters, identical mesh geometry, not re-run): %d." % (
              sum(1 for r in K.values() if not r.get("ok")), sum(1 for r in K.values() if r.get("alias_of"))), ""]
    atomic_write(os.path.join(OUT, "summary.md"), "\n".join(L) + "\n")
    log("report: wrote %s" % os.path.join(OUT, "summary.md"))


# ----------------------------------------------------------------------------- status, all, launch
def pid_alive(pid):
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, ValueError, TypeError):
        return False


def stage_status():
    cache = Cache()
    K = cache.by_key
    val, sea, chk = (read_json(os.path.join(OUT, n)) for n in ("validate.json", "search.json", "check.json"))
    done = lambda d: bool(d and d.get("complete"))
    stage = "validate" if not done(val) else "search" if not done(sea) else "check" if not done(chk) else "finished"
    if done(val) and not val.get("trustworthy"):
        stage = "STOPPED: validation failed (see summary.md)"
    try:
        pid = open(os.path.join(OUT, "run.pid")).read().strip()
    except OSError:
        pid = None
    print("run process: %s" % ("pid %s, running" % pid if pid_alive(pid) else "not running (last pid %s)" % pid))
    print("stage: %s" % stage)
    by = {}
    for r in K.values():
        by.setdefault(r.get("tag"), []).append(r)
    print("evaluations in cache: %d  (%s)  failed: %d" % (
        len(K), ", ".join("%s %d" % (k, len(v)) for k, v in by.items()), sum(1 for r in K.values() if not r.get("ok"))))
    srch = [r for r in by.get("search", []) if r.get("ok")]
    r0 = K.get(eval_key("row", X0, EPS_NOM, M_BASE))
    if r0:
        print("start design:  score %.2f  worst %.2f dB  mean %.2f dB" % (r0["score"], r0["worst"], r0["mean"]))
    if srch:
        b = min(srch, key=lambda r: r["score"])
        print("best so far:   score %.2f  worst %.2f dB  mean %.2f dB\n   %s" % (b["score"], b["worst"], b["mean"], _fmt_x(b["x"])))
    real = sorted((r for r in K.values() if not r.get("alias_of") and r.get("wall")), key=lambda r: r["t"])
    recent = [r for r in real if r["t"] > time.time() - 3600 and r["model"] == "row" and r["mesh"]["nh"] == M_BASE["nh"]]
    if len(recent) >= 2:
        rate = (len(recent) - 1) / max(recent[-1]["t"] - recent[0]["t"], 1.0)
        print("rate, last hour: %.1f row evaluations per minute (mean solver wall %.0f s)" % (
            rate * 60, np.mean([r["wall"] for r in recent])))
    else:
        rate = None
    if stage == "search":
        left = max(PLAN_SEARCH - len(srch) - len(by.get("finalist", [])), 0)
        tile_w = [r["wall"] for r in real if r["model"] == "tile"]
        chk_s = 6 * (np.mean(tile_w) if tile_w else 1500.0) + 600
        if rate:
            print("ETA (rough): search %.1f h for at most %d more evaluations, then check about %.1f h" % (
                left / rate / 3600, left, chk_s / 3600))
    elif stage == "check":
        n = len((chk or {}).get("keys", {}))
        tile_w = [r["wall"] for r in real if r["model"] == "tile"]
        print("check: %d of 8 runs done; ETA (rough) %.1f h" % (n, (8 - n) * (np.mean(tile_w) if tile_w else 1500.0) / 3600))
    try:
        with open(os.path.join(OUT, "run.log")) as f:
            tail = f.readlines()[-3:]
        print("last log lines:\n   " + "   ".join(tail), end="")
    except OSError:
        pass


def stage_all():
    os.makedirs(OUT, exist_ok=True)
    pidf = os.path.join(OUT, "run.pid")
    try:
        old = open(pidf).read().strip()
    except OSError:
        old = None
    if old and old != str(os.getpid()) and pid_alive(old):
        log("another run is alive (pid %s); not starting a second one" % old)
        return 2
    atomic_write(pidf, "%d\n" % os.getpid())
    log("=== all: pid %d, engine %s, %d jobs x %d threads ===" % (os.getpid(), ENGINE, JOBS, THREADS))
    cache = Cache()
    log("cache: %d evaluations on disk" % len(cache.by_key))
    val = read_json(os.path.join(OUT, "validate.json"))
    if not (val and val.get("complete")):
        val = stage_validate(cache)
    else:
        log("validate: already complete, skipped")
    if not val.get("trustworthy"):
        log("STOP: the model failed validation; no search was run. See summary.md.")
        stage_report(cache)
        return 1
    sea = read_json(os.path.join(OUT, "search.json"))
    if not (sea and sea.get("complete")):
        stage_search(cache)
        stage_report(cache)
    else:
        log("search: already complete, skipped")
    chk = read_json(os.path.join(OUT, "check.json"))
    if not (chk and chk.get("complete")):
        stage_check(cache)
    else:
        log("check: already complete, skipped")
    stage_report(cache)
    log("=== all: finished ===")
    return 0


def launch():
    """Start `all` in its own session (it survives the terminal or agent that started it), under caffeinate."""
    os.makedirs(OUT, exist_ok=True)
    if os.fork():
        time.sleep(1.5)
        print("started; log: %s\npid file: %s" % (os.path.join(OUT, "run.log"), os.path.join(OUT, "run.pid")))
        return
    os.setsid()
    if os.fork():
        os._exit(0)
    os.chdir(REPO)
    fd = os.open(os.path.join(OUT, "run.log"), os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o644)
    nul = os.open(os.devnull, os.O_RDONLY)
    os.dup2(nul, 0)
    os.dup2(fd, 1)
    os.dup2(fd, 2)
    os.execvp("caffeinate", ["caffeinate", "-i", "-s", sys.executable, "-u", os.path.abspath(__file__), "all"])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cmd", choices=["validate", "search", "check", "report", "status", "all", "launch", "_worker"])
    ap.add_argument("arg", nargs="?")
    a = ap.parse_args()
    if a.cmd == "_worker":
        return worker_main(a.arg)
    if a.cmd == "status":
        return stage_status()
    if a.cmd == "report":
        return stage_report()
    if a.cmd == "launch":
        return launch()
    signal.signal(signal.SIGTERM, _terminate)
    signal.signal(signal.SIGINT, _terminate)
    os.makedirs(SCRATCH, exist_ok=True)
    if a.cmd == "all":
        return stage_all()
    cache = Cache()
    if a.cmd == "validate":
        stage_validate(cache)
    elif a.cmd == "search":
        stage_search(cache)
    elif a.cmd == "check":
        stage_check(cache)
    stage_report(cache)


if __name__ == "__main__":
    sys.exit(main() or 0)
