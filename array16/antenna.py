"""Antenna models: patch design, and the 16 patches as one network coupled through space.

Each patch is the transmission-line model of Balanis (Antenna Theory, 4th ed., 14.2.1): a
length L of wide microstrip between two radiating edges. Each edge is a slot, a uniform
magnetic line current of length W on an infinite ground plane. Every slot couples to every
other slot of every patch through space, with mutual admittances from the reaction integral
between them, so the array is one network. Solving it gives the S-matrix at the 16 probes
and, for a wave into each probe, the voltage on all 32 edges, hence the embedded element
patterns the simulator can use in place of its ideal element.

Included: radiation and coupling through space; dielectric and conductor loss in the patch;
surface-wave loss as extra conductance on each patch's space-wave conductance (Jackson and
Alexopoulos, IEEE TAP 1991); the probe's inductance; eps_r and h in the resonance.

Accuracy, against a moment-method solution with the grounded-slab Green's function (the model
review of 2026-11-15): on 1.28 mm FR-4 at the 26 mm pitch the nearest E-plane coupling is about
right (surface waves replace the grazing space wave the slot model uses); H-plane coupling is
1.5-2 dB too weak and diagonal coupling 2-3 dB too weak, because the non-radiating edges carry no
magnetic current here; on lower-eps_r substrates E-plane coupling is overstated by 2-5 dB. The
model does not hold when facing edges are closer than their fringing fields (pitch - L below two
open-end extensions, as for air patches at 26 mm): there a moment method gives E-plane coupling
about 6 dB stronger and the patterns differ; PatchArray.coupling_reliable says which case applies.
Bandwidth is about 20 % narrower than the Jackson-Alexopoulos estimate (conservative).
Not included: the finite ground plane, radiation from the probe, the cross-polarised mode.
All patches are polarised along x (E-plane along the rows of the array).
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache

import numpy as np

from . import board as B

ETA0 = 376.730313668
MU0 = 1.25663706212e-6
EULER = 0.5772156649


# ---- substrates and closed-form microstrip results ---------------------------------------------
@dataclass(frozen=True)
class Substrate:
    name: str
    eps_r: float
    tan_d: float
    h: float                     # m, from the patch to its ground plane
    sigma: float = 5.8e7         # conductor conductivity, S/m
    rough: float = 1.0           # conductor-loss multiplier for surface roughness


def eps_eff(eps_r, h, W):
    """Effective permittivity of a microstrip of width W >= h (Hammerstad; Balanis 14-1)."""
    return (eps_r + 1) / 2 + (eps_r - 1) / 2 / np.sqrt(1 + 12 * h / W)


def open_end(eps_r, h, W):
    """How far the fringing field extends each open end (Hammerstad; Balanis 14-2)."""
    ee = eps_eff(eps_r, h, W)
    return 0.412 * h * (ee + 0.3) * (W / h + 0.264) / ((ee - 0.258) * (W / h + 0.8))


def z_line(eps_r, h, W):
    """Characteristic impedance of a microstrip of width W >= h (Pozar 3.196)."""
    return ETA0 / np.sqrt(eps_eff(eps_r, h, W)) / (W / h + 1.393 + 0.667 * np.log(W / h + 1.444))


def hed_efficiency(eps_r, h, f):
    """Share of a horizontal electric dipole's power that leaves as space wave rather than
    surface wave (Jackson and Alexopoulos 1991). 1 for air."""
    c1 = 1 - 1 / eps_r + 0.4 / eps_r ** 2
    return 1 / (1 + 0.75 * np.pi * (2 * np.pi * f / B.C0 * h) / c1 * (1 - 1 / eps_r) ** 3)


def jackson_alexopoulos(sub: Substrate, W, L, f0=B.F0) -> dict:
    """Q factors, radiation efficiency and VSWR-2 bandwidth of a resonant patch from the
    cavity-model formulas of Jackson and Alexopoulos (IEEE TAP 39(3), 1991): an independent
    estimate to set beside the transmission-line network."""
    er, h, k0 = sub.eps_r, sub.h, 2 * np.pi * f0 / B.C0
    c1 = 1 - 1 / er + 0.4 / er ** 2
    p = 1 - 0.16605 / 20 * (k0 * W) ** 2 + 0.00761 * 3 / 560 * (k0 * W) ** 4 - 0.0914153 / 10 * (k0 * L) ** 2
    q_sp = 3 / 16 * er / (p * c1) * L / W / (h * f0 / B.C0)
    e = hed_efficiency(er, h, f0)
    loss = (1 - e) / (q_sp * e) + sub.tan_d + sub.rough / (h * np.sqrt(np.pi * f0 * MU0 * sub.sigma))
    q = 1 / (1 / q_sp + loss)
    return dict(q_space=q_sp, q_total=q, efficiency=q / q_sp, bw_vswr2=f0 / (np.sqrt(2) * q))


def textbook_patch(eps_r, h, f0=B.F0, w_max=None):
    """(W, L) from the closed-form design procedure (Balanis 14.2.1)."""
    W = B.C0 / (2 * f0) * np.sqrt(2 / (eps_r + 1))
    if w_max is not None:
        W = min(W, w_max)
    return W, B.C0 / (2 * f0 * np.sqrt(eps_eff(eps_r, h, W))) - 2 * open_end(eps_r, h, W)


# ---- slots ----------------------------------------------------------------------------------------
@lru_cache(maxsize=None)
def _legendre(n):
    return np.polynomial.legendre.leggauss(n)


def _gauss(W, n):
    t, w = _legendre(n)
    return t * W / 2, w * W / 2


def slot_conductance(k, W, n=64):
    """Radiation conductance of a slot of length W with uniform voltage, into one half space."""
    y, w = _gauss(W, n)
    kern = k * np.sinc(k * np.abs(y[:, None] - y[None, :]) / np.pi)       # sin(kR) / R
    return (k * (w @ kern @ w) - 2 + 2 * np.sinc(k * W / np.pi)) / (2 * np.pi * ETA0)


def slot_mutual(k, W, dx, dy, n=64):
    """Mutual admittance of two parallel slots of length W along y, radiating into the same half
    space; the second is displaced by (dx, dy) and must not overlap the first.

    Reaction of the field of one uniform magnetic line current (doubled by its image in the
    ground plane) on the other, including the magnetic charges at the slot ends.
    """
    y, w = _gauss(W, n)
    G = lambda r: np.exp(-1j * k * r) / r
    A = w @ G(np.hypot(dx, dy + y[:, None] - y[None, :])) @ w
    e = W / 2
    X = sum(s * G(np.hypot(dx, dy + a - b)) for s, a, b in ((1, e, e), (-1, e, -e), (-1, -e, e), (1, -e, -e)))
    return (1j * k / ETA0 * A + X / (1j * k * ETA0)) / (2 * np.pi)


# ---- one patch -------------------------------------------------------------------------------------
def probe_reactance(sub: Substrate, probe_r, f):
    k0 = 2 * np.pi * f / B.C0
    return ETA0 * k0 * sub.h / (2 * np.pi) * (np.log(2 / (k0 * np.sqrt(sub.eps_r) * probe_r)) - EULER)


@dataclass(frozen=True)
class Patch:
    """A probe-fed rectangular patch, resonant along x and so polarised along x."""

    sub: Substrate
    W: float                     # m, along y
    L: float                     # m, along x
    feed: float                  # m, probe position along x from the patch centre
    probe_r: float = 0.2e-3      # m
    x_comp: float = 0.0          # ohm: reactance a series capacitor cancels at F0 (capacitive or L-probe feed)

    def line(self, f):
        """Characteristic admittance and propagation constant of the patch as a wide microstrip."""
        s, k0 = self.sub, 2 * np.pi * f / B.C0
        ee, z0 = eps_eff(s.eps_r, s.h, self.W), z_line(s.eps_r, s.h, self.W)
        a_d = 0.0 if s.eps_r == 1 else k0 * s.eps_r * (ee - 1) * s.tan_d / (2 * np.sqrt(ee) * (s.eps_r - 1))
        a_c = s.rough * np.sqrt(np.pi * f * MU0 / s.sigma) / (z0 * self.W)          # Pozar 3.199
        return 1 / z0, a_d + a_c + 1j * k0 * np.sqrt(ee)

    def probe_x(self, f):
        """Reactance at the probe: the probe crossing the substrate (Garg et al., Microstrip Antenna
        Design Handbook) less the fixed series capacitor, if any."""
        return probe_reactance(self.sub, self.probe_r, f) - self.x_comp * B.F0 / f

    def edge(self, f):
        """(admittance of one radiating edge on its own, its radiation conductance)."""
        s, k0 = self.sub, 2 * np.pi * f / B.C0
        g = slot_conductance(k0, self.W)
        y0, gam = self.line(f)
        b = y0 * np.tan(gam.imag * open_end(s.eps_r, s.h, self.W))      # fringing, as an open-end extension
        return g / hed_efficiency(s.eps_r, s.h, f) + 1j * b, g


# ---- the array as one network -----------------------------------------------------------------------
class PatchArray:
    """Patches at the given positions (default: the board's 16), coupled through space.

    solve(f) returns the S-matrix at the probes and, for a unit incident wave (peak, so 0.5 W)
    into each probe in turn with the others matched, the voltage on every radiating edge.
    """

    def __init__(self, patch: Patch, x=B.X, y=B.Y, z_ref=50.0):
        self.p, self.z_ref = patch, z_ref
        self.x, self.y = np.atleast_1d(np.asarray(x, float)), np.atleast_1d(np.asarray(y, float))
        self.n = len(self.x)
        self.xs = np.concatenate([self.x - patch.L / 2, self.x + patch.L / 2])     # the -x edges, then the +x edges
        self.ys = np.concatenate([self.y, self.y])
        # An edge's magnetic current is -2 n x E: the two edges face opposite ways, so for the same
        # voltage their currents point opposite ways (TM10's opposite edge voltages then radiate in phase).
        self.sign = np.concatenate([-np.ones(self.n), np.ones(self.n)])
        # Facing radiating edges of neighbouring patches closer than their fringing fields, or than
        # about 0.2 wavelength: there the slot model's coupling is 5-7 dB too weak and its patterns are
        # off (moment-method check, model review 2026-11-15).
        gap = np.abs(self.x[:, None] - self.x[None, :]) - patch.L
        same_row = (np.abs(self.y[:, None] - self.y[None, :]) < patch.W) & ~np.eye(self.n, dtype=bool)
        limit = max(2 * open_end(patch.sub.eps_r, patch.sub.h, patch.W), 0.2 * B.C0 / B.F0)
        self.coupling_reliable = bool(not np.any(same_row & (gap < limit)))
        self._solved: dict[float, dict] = {}

    def solve(self, f: float) -> dict:
        f = float(f)
        if f in self._solved:
            return self._solved[f]
        p, n, m = self.p, self.n, 2 * self.n
        k0 = 2 * np.pi * f / B.C0
        y_edge, g_edge = p.edge(f)
        e_sw = hed_efficiency(p.sub.eps_r, p.sub.h, f)
        ys, gr, memo = np.empty((m, m), complex), np.empty((m, m)), {}
        for i in range(m):
            ys[i, i], gr[i, i] = y_edge, g_edge
            for j in range(i + 1, m):
                key = (round(abs(self.xs[j] - self.xs[i]) * 1e8), round(abs(self.ys[j] - self.ys[i]) * 1e8))
                if key not in memo:
                    memo[key] = slot_mutual(k0, p.W, key[0] * 1e-8, key[1] * 1e-8)
                y = self.sign[i] * self.sign[j] * memo[key]
                gr[i, j] = gr[j, i] = y.real
                if j == i + n:                    # a patch's own two edges: its surface wave scales with
                    y = y.real / e_sw + 1j * y.imag      # its whole space-wave conductance, G1 + G12
                ys[i, j] = ys[j, i] = y
        y0, gam = p.line(f)
        l1, l2 = p.L / 2 + p.feed, p.L / 2 - p.feed                  # -x edge to probe, probe to +x edge
        Y = np.zeros((m + n, m + n), complex)                       # nodes: edges, then probes
        Y[:m, :m] = ys
        a, b, c = np.arange(n), np.arange(n) + n, np.arange(n) + m
        Y[a, a] += y0 / np.tanh(gam * l1)
        Y[b, b] += y0 / np.tanh(gam * l2)
        Y[c, c] += y0 / np.tanh(gam * l1) + y0 / np.tanh(gam * l2)
        Y[a, c] = Y[c, a] = -y0 / np.sinh(gam * l1)
        Y[b, c] = Y[c, b] = -y0 / np.sinh(gam * l2)
        T = -np.linalg.solve(Y[:m, :m], Y[:m, m:])                  # edge voltages per volt on the probes
        z_patch = np.linalg.inv(Y[m:, m:] + Y[m:, :m] @ T)
        z = z_patch + 1j * p.probe_x(f) * np.eye(n)
        eye = np.eye(n)
        S = (z - self.z_ref * eye) @ np.linalg.inv(z + self.z_ref * eye)
        v = T @ z_patch @ (eye - S) / np.sqrt(self.z_ref)            # edge voltages for a unit wave into each port
        out = self._solved[f] = dict(S=S, v=v, z=z, g_rad=gr)
        return out

    def eep(self, u, v, f=B.F0) -> np.ndarray:
        """Embedded element patterns, co-polar (x) component: complex amplitude per element whose
        squared magnitude is the realised gain toward (u, v). Shape (*u.shape, n)."""
        u, v = np.asarray(u, float), np.asarray(v, float)
        k = 2 * np.pi * f / B.C0
        V = self.solve(f)["v"]
        s2 = u * u + v * v
        c = np.sqrt(np.clip(1 - s2, 0, None))
        pol = np.where(s2 > 0, (u * u + c * v * v) / np.where(s2 > 0, s2, 1), 1.0)        # Ludwig-3 co-polar
        ly = self._ly(u, v, k) @ V
        return -1j * k / np.sqrt(4 * np.pi * ETA0) * ly * pol[..., None]

    def _ly(self, u, v, k):
        """Radiation vector (y component) of each edge for a unit edge voltage, shape (*u.shape, 2n)."""
        ph = self.sign * np.exp(1j * k * (u[..., None] * self.xs + v[..., None] * self.ys))
        return 2 * self.p.W * np.sinc(k * self.p.W * v / (2 * np.pi))[..., None] * ph

    def radiated_power(self, f, weights, n=90) -> float:
        """Power radiated into the upper half space for incident waves `weights` (peak), by integrating the far field."""
        k = 2 * np.pi * f / B.C0
        th, wt = np.polynomial.legendre.leggauss(n)
        th, wt = (th + 1) * np.pi / 4, wt * np.pi / 4
        ph = np.linspace(0, 2 * np.pi, 2 * n, endpoint=False)
        T, P = np.meshgrid(th, ph, indexing="ij")
        u, v = np.sin(T) * np.cos(P), np.sin(T) * np.sin(P)
        ly = self._ly(u, v, k) @ (self.solve(f)["v"] @ np.asarray(weights, complex))
        U = k * k * np.abs(ly) ** 2 * (np.cos(P) ** 2 + np.cos(T) ** 2 * np.sin(P) ** 2) / (32 * np.pi ** 2 * ETA0)
        return float(np.sum(U * np.sin(T) * wt[:, None]) * (2 * np.pi / (2 * n)))


# ---- design and characterisation of one element --------------------------------------------------
def design(sub: Substrate, f0=B.F0, w_max=B.PITCH - 4e-3, probe_r=0.2e-3, probe_comp=False) -> Patch:
    """Width from the textbook formula (at most w_max, so neighbours keep a gap), then length
    and probe position for the best match of the isolated patch at f0. The probe stays at least
    0.1 mm inside the patch edge; if no position inside matches, the result shows it in S11.
    probe_comp: a series capacitor, chosen here and then fixed, cancels the probe's reactance at F0."""
    W, L0 = textbook_patch(sub.eps_r, sub.h, f0, w_max)
    x_comp = probe_reactance(sub, probe_r, B.F0) if probe_comp else 0.0
    x_max = lambda L: 0.5 - (probe_r + 0.1e-3) / L                # probe centre as a fraction of L from the centre
    s11 = lambda L, x: abs(PatchArray(Patch(sub, W, L, -x * L, probe_r, x_comp), [0.0], [0.0]).solve(f0)["S"][0, 0])
    Ls, xs = L0 * np.linspace(0.85, 1.08, 24), np.linspace(0.02, x_max(L0 * 0.85), 24)
    for _ in range(4):
        grid = np.array([[s11(L, min(x, x_max(L))) for x in xs] for L in Ls])
        i, j = np.unravel_index(np.argmin(grid), grid.shape)
        dL, dx = Ls[1] - Ls[0], xs[1] - xs[0]
        Ls = np.linspace(Ls[i] - 2 * dL, Ls[i] + 2 * dL, 13)
        xs = np.clip(np.linspace(xs[j] - 2 * dx, xs[j] + 2 * dx, 13), 0.005, x_max(Ls[0]))
    L, x = Ls[6], min(xs[6], x_max(Ls[6]))
    return Patch(sub, W, L, -x * L, probe_r, x_comp)


def characterise(patch: Patch, f0=B.F0) -> dict:
    """Resonance, bandwidth, efficiency and gain of the isolated patch, and its sensitivity to eps_r."""
    def sweep(p):
        iso = PatchArray(p, [0.0], [0.0])
        f = np.linspace(0.85 * f0, 1.15 * f0, 871)
        s = np.array([abs(iso.solve(x)["S"][0, 0]) for x in f])
        return f, s
    f, s = sweep(patch)
    k = int(np.argmin(np.abs(f - f0)))
    lo, hi = k, k
    while lo > 0 and s[lo - 1] < 10 ** (-10 / 20):
        lo -= 1
    while hi < len(f) - 1 and s[hi + 1] < 10 ** (-10 / 20):
        hi += 1
    iso = PatchArray(patch, [0.0], [0.0])
    r = iso.solve(f0)
    p_in = (1 - abs(r["S"][0, 0]) ** 2) / 2
    p_rad = iso.radiated_power(f0, [1.0])
    g = float(abs(iso.eep(np.array(0.0), np.array(0.0), f0)[0]) ** 2)
    def f_shift(**change):
        fd, sd = sweep(replace(patch, sub=replace(patch.sub, **change)))
        return fd[np.argmin(sd)] - f[np.argmin(s)]
    eps_shift = (f_shift(eps_r=patch.sub.eps_r + 0.1) - f_shift(eps_r=max(1.0, patch.sub.eps_r - 0.1))) / (
        0.2 if patch.sub.eps_r > 1 else 0.1) * 0.1
    h_shift = (f_shift(h=patch.sub.h * 1.1) - f_shift(h=patch.sub.h * 0.9)) / 2
    ja = jackson_alexopoulos(patch.sub, patch.W, patch.L, f0)
    return dict(
        cavity_model_efficiency=float(ja["efficiency"]),
        cavity_model_bw_vswr2_MHz=float(ja["bw_vswr2"] / 1e6),
        W_mm=patch.W * 1e3, L_mm=patch.L * 1e3, feed_mm=abs(patch.feed) * 1e3,
        f_res_GHz=float(f[np.argmin(s)] / 1e9),
        s11_f0_dB=float(20 * np.log10(s[k])),
        bw_10dB_MHz=float((f[hi] - f[lo]) / 1e6) if s[k] < 10 ** (-10 / 20) else 0.0,
        covers_ism=bool(s[k] < 10 ** (-10 / 20) and f[lo] <= 5.725e9 and f[hi] >= 5.875e9),
        efficiency=float(p_rad / p_in),
        surface_wave_share=float(1 - hed_efficiency(patch.sub.eps_r, patch.sub.h, f0)),
        gain_dbi=float(10 * np.log10(g)),
        directivity_dbi=float(10 * np.log10(g * 0.5 / p_rad)),          # realised gain x available power / radiated power
        probe_x_ohm=float(patch.probe_x(f0)),
        f_shift_MHz_per_0p1_eps=float(eps_shift / 1e6),          # resonance moves this much for eps_r + 0.1
        f_shift_MHz_per_10pct_h=float(h_shift / 1e6),            # ... and for the height 10 % larger
    )


# ---- a coupling model with chosen coefficients -----------------------------------------------------
class CoupledIdeal:
    """The simulator's ideal element (cos^q, the default realised array gain) with chosen complex
    coupling: element n receives sum_m C[n, m] times what element m alone would. Give C whole
    (for example fit_coupling() of a PatchArray, scaled), or E-plane (along x), H-plane (along y)
    and diagonal nearest-neighbour coefficients. For sensitivity studies."""

    def __init__(self, c_e=0j, c_h=0j, c_d=0j, q=B.ELEM_Q, gain_dbi=B.ARRAY_GAIN_DBI, C=None):
        dx = np.round((B.X[:, None] - B.X[None, :]) / B.PITCH)
        dy = np.round((B.Y[:, None] - B.Y[None, :]) / B.PITCH)
        self.C = np.asarray(C, complex) if C is not None else (
            np.eye(B.N_ELEM) + c_e * ((abs(dx) == 1) & (dy == 0)) + c_h * ((dx == 0) & (abs(dy) == 1))
            + c_d * ((abs(dx) == 1) & (abs(dy) == 1)))
        self.q, self.g0 = q, 10 ** (gain_dbi / 10) / B.N_ELEM

    def eep(self, u, v, f=B.F0) -> np.ndarray:
        u, v = np.asarray(u, float), np.asarray(v, float)
        c = np.sqrt(np.clip(1 - u * u - v * v, 0, None))
        base = np.sqrt(self.g0 * c ** self.q)[..., None] * np.exp(2j * np.pi * f / B.C0 * (u[..., None] * B.X + v[..., None] * B.Y))
        return base @ self.C.T


def fit_coupling(arr: PatchArray, f=B.F0, theta_max_deg=70.0, n=24) -> tuple[np.ndarray, float]:
    """The receive-coupling matrix C of a patch array: its embedded patterns fitted, over directions
    out to theta_max, as sum_m C[n, m] times the isolated patch's pattern moved to element m.
    Returns (C, relative rms residual of the fit)."""
    t, p = np.meshgrid(np.radians(np.linspace(0, theta_max_deg, n)), np.radians(np.linspace(0, 360, 2 * n, endpoint=False)))
    u, v = (np.sin(t) * np.cos(p)).ravel(), (np.sin(t) * np.sin(p)).ravel()
    k = 2 * np.pi * f / B.C0
    iso = PatchArray(arr.p, [0.0], [0.0]).eep(u, v, f)[:, :1] * np.exp(1j * k * (u[:, None] * arr.x + v[:, None] * arr.y))
    e = arr.eep(u, v, f)
    C = np.linalg.lstsq(iso, e, rcond=None)[0].T
    return C, float(np.linalg.norm(iso @ C.T - e) / np.linalg.norm(e))
