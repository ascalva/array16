"""Simulator: a virtual board, scene and radios behind the backend interface.

The host library talks to SimBackend exactly as it will talk to the hardware:
controller frames in, I/Q samples out. Everything physical lives here:
element patterns, switch and line losses, LNA and radio noise, unit-to-unit
errors, unknown oscillator phases, stream start offsets and USB latency.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import board as B
from .backend import Capture
from .board import ControlWord
from .fwsim import ControllerModel


def db2amp(db):
    return 10 ** (np.asarray(db) / 20)


def dbm2w(dbm):
    return 10 ** ((np.asarray(dbm) - 30) / 10)


def w2dbm(w):
    return 10 * np.log10(w) + 30


class Noise:
    """Thermal noise in units of sqrt(W) for one sample rate (noise bandwidth = fs)."""

    def __init__(self, rng, bandwidth):
        self.rng, self.sigma = rng, np.sqrt(B.K_BOLTZ * B.T0 * bandwidth)

    def cn(self, shape, power=1.0):
        g = self.rng.standard_normal(shape) + 1j * self.rng.standard_normal(shape)
        return self.sigma * np.sqrt(np.asarray(power) / 2) * g

    def atten(self, x, loss_db):
        """Matched attenuator at 290 K: attenuates and tops the noise back up."""
        a2 = 10 ** (-np.asarray(loss_db) / 10)
        return x * np.sqrt(a2) + self.cn(x.shape, 1 - a2)

    def amp(self, x, gain_db, nf_db):
        return (x + self.cn(x.shape, 10 ** (nf_db / 10) - 1)) * db2amp(gain_db)


class BoardModel:
    """RF behaviour of one physical board, including its unit-to-unit errors."""

    def __init__(self, seed=1, amp_sigma_db=0.5, phase_sigma_deg=10.0, stage_sigma_deg=3.0,
                 bypass_sigma_deg=5.0, port_sigma_deg=40.0, antenna=None):
        # antenna: anything with eep(u, v, f) -> (16,) complex embedded element patterns (antenna.PatchArray,
        # antenna.CoupledIdeal). None: identical ideal elements, cos(theta)**ELEM_Q, no coupling.
        self.antenna = antenna
        r = np.random.default_rng(seed)
        ph = lambda sig, n: np.exp(1j * np.radians(r.normal(0, sig, n)))
        common = db2amp(r.normal(0, amp_sigma_db, 16)) * ph(phase_sigma_deg, 16)
        self.rx_err = common * ph(stage_sigma_deg, 16)            # path through the LNA
        self.tx_err = common * ph(bypass_sigma_deg, 16)           # path through the bypass line
        self.stage_err = np.stack([ph(stage_sigma_deg, 16), ph(stage_sigma_deg, 16)])  # long path of [180, 90]
        self.port_err = ph(port_sigma_deg, 4)                     # quadrant receive path and its cable
        self.tx_port_err = ph(port_sigma_deg, 4)
        self.g0 = 10 ** (B.ARRAY_GAIN_DBI / 10) / B.N_ELEM        # embedded element gain at boresight
        rr = np.hypot(B.X, B.Y)                                   # centre probe: nearer elements couple harder,
        self.probe = (db2amp(-28 - 20 * np.log10(rr / rr.min()))  # the lower half sees the opposite sign
                      * np.exp(-2j * np.pi * B.F0 / B.C0 * 1.6 * rr) * np.sign(B.Y))

    # ---- pieces --------------------------------------------------------------
    def elem_amp(self, u, v):
        cos_t = np.sqrt(max(0.0, 1 - u * u - v * v))
        return np.sqrt(self.g0 * cos_t ** B.ELEM_Q)

    def elem_resp(self, u, v, f):
        """Each antenna port's response to a plane wave from (u, v): sqrt(realised gain) with phase."""
        if self.antenna is None:
            return self.elem_amp(u, v) * B.steering(u, v, f)
        return self.antenna.eep(np.asarray(u, float), np.asarray(v, float), f)

    def ps_loss(self, w: ControlWord):
        return (4 * B.IL_SW + B.LINE_DB_PER_CM * (B.LEN_PS_BASE_CM + B.LEN_180_CM * w.p180 + B.LEN_90_CM * w.p90))

    def ps_phase(self, w: ControlWord, f):
        lag = (np.pi * w.p180 + np.pi / 2 * w.p90) * f / B.F0     # delay lines: lag scales with frequency
        return np.exp(-1j * lag) * np.where(w.p180, self.stage_err[0], 1) * np.where(w.p90, self.stage_err[1], 1)

    # ---- receive -------------------------------------------------------------
    def receive(self, w: ControlWord, x_el, f, nz: Noise):
        """Antenna ports (16, n) -> quadrant combiner outputs (4, n)."""
        x = x_el + nz.cn(x_el.shape)                              # sky/room noise at 290 K
        x = nz.atten(x, B.FEED_LOSS + B.IL_SW)
        x = np.where(w.rx[:, None], nz.amp(x, B.LNA_GAIN, B.LNA_NF), x)
        x = nz.atten(x, B.IL_SW)
        x = nz.atten(x, self.ps_loss(w)[:, None]) * (self.ps_phase(w, f) * self.rx_err)[:, None]
        quad = np.stack([x[B.QUAD_ELEMS[q]].sum(0) / 2 for q in range(4)])
        return nz.atten(quad, 2 * B.WILK_EXCESS + B.LINE_DB_PER_CM * B.LEN_QUAD_TREE_CM)

    def to_sockets(self, w: ControlWord, quad, nz: Noise):
        """Quadrant outputs -> (four quadrant RX sockets (4, n), HackRF socket (1, n))."""
        iso = db2amp(-B.ISO_SW)
        thru, leak = nz.atten(quad, B.IL_SW), quad * iso
        to_quad, to_shared = (leak, thru) if w.hackrf_mode else (thru, leak)
        lna = (B.LNA_GAIN, B.LNA_NF) if w.lna2 else (B.LNA_OFF_GAIN, -B.LNA_OFF_GAIN)
        y = nz.atten(to_quad, B.IL_SW) * np.where(w.qrx, 1, iso)[:, None]
        y = nz.amp(y, *lna) * self.port_err[:, None]
        h = nz.atten(to_shared, B.LINE_DB_PER_CM * B.LEN_DIAG_CM).sum(0, keepdims=True) / 2
        h = nz.atten(h, 2 * B.WILK_EXCESS + B.IL_SW) * (1 if w.hrx else iso)
        h = nz.atten(nz.amp(h, *lna), B.IL_SW + B.LINE_DB_PER_CM * B.LEN_HACKRF_CM)
        return y, h

    # ---- transmit ------------------------------------------------------------
    def transmit(self, w: ControlWord, f, drive):
        """Wave amplitude (sqrt W) reaching each antenna port.

        drive: complex amplitude (sqrt W) at each of the four quadrant TX sockets
        (four-radio mode) or at the HackRF socket (HackRF mode).
        """
        iso = db2amp(-B.ISO_SW)
        el = (db2amp(-(2 * B.IL_SW + B.FEED_LOSS + self.ps_loss(w))) * self.ps_phase(w, f) * self.tx_err
              * np.where(w.rx, db2amp(B.LNA_OFF_GAIN), 1))        # an element left in receive blocks transmit
        tree = db2amp(-(2 * B.WILK_EXCESS + B.LINE_DB_PER_CM * B.LEN_QUAD_TREE_CM)) / 2
        if w.hackrf_mode:
            loss = B.LINE_DB_PER_CM * (B.LEN_HACKRF_CM + B.LEN_DIAG_CM) + 3 * B.IL_SW + 2 * B.WILK_EXCESS
            q_in = np.full(4, drive * db2amp(-loss) / 2 * (iso if w.hrx else 1), dtype=complex)
        else:
            q_in = np.asarray(drive) * db2amp(-2 * B.IL_SW) * self.tx_port_err * np.where(w.qrx, iso, 1)
        return q_in[B.QUAD] * tree * el


@dataclass
class Source:
    """An emitter. power_dbm is what an isotropic antenna at the array would receive.
    range_m: distance from the array centre; None means far away (a plane wave). A finite range
    adds the wavefront's curvature across the aperture, as a calibration source on the bench does."""

    theta_deg: float
    phi_deg: float
    power_dbm: float = -70.0
    kind: str = "tone"            # tone | qpsk
    f_offset: float = 20e3        # Hz from the tuned frequency
    symbol_rate: float = 200e3
    seed: int = 7
    range_m: float | None = None

    def curvature(self, f) -> np.ndarray | float:
        """Per-element phase factor beyond the plane wave, for a source at range_m."""
        if self.range_m is None:
            return 1.0
        u, v = B.direction(self.theta_deg, self.phi_deg)
        w = np.sqrt(max(0.0, 1 - u * u - v * v))
        src = self.range_m * np.array([u, v, w])
        rn = np.sqrt((src[0] - B.X) ** 2 + (src[1] - B.Y) ** 2 + src[2] ** 2)
        return np.exp(-2j * np.pi * f / B.C0 * (rn - self.range_m + B.X * u + B.Y * v))

    def waveform(self, t):
        carrier = np.exp(2j * np.pi * (self.f_offset * t + self.seed * 0.137))
        if self.kind == "tone":
            return carrier
        table = np.random.default_rng(self.seed).integers(0, 4, 65536)
        sym = np.exp(1j * (np.pi / 4 + np.pi / 2 * table[np.floor(t * self.symbol_rate).astype(np.int64) % 65536]))
        return carrier * sym


@dataclass
class Reflector:
    """A point target lit by the auxiliary patch; motion(t) is its displacement along range, m."""

    theta_deg: float
    phi_deg: float
    range_m: float
    rcs_m2: float = 1.0
    motion: callable = lambda t: 0.0 * t


class _Radio:
    def __init__(self, fs, bits, nf_db, boards):
        self.fs, self.bits, self.nf_db, self.boards = fs, bits, nf_db, np.asarray(boards)
        self.channels = len(boards)
        self.lo = np.ones(self.channels, complex)

    def retune(self, rng):
        """Every radio board locks with its own arbitrary oscillator phase, and its stream restarts."""
        self.lo = np.exp(1j * rng.uniform(0, 2 * np.pi, self.boards.max() + 1))[self.boards]
        self.tx_lo = np.exp(1j * rng.uniform(0, 2 * np.pi, self.boards.max() + 1))[self.boards]
        self.restart(rng)

    def restart(self, rng):
        """Each board's sample stream starts at its own instant; the offset then stays put."""
        self.skew = rng.integers(0, 40, self.boards.max() + 1)[self.boards] * (self.boards > 0)

    def digitise(self, x, skew, n, nz):
        x = nz.amp(x, 0.0, self.nf_db) * self.lo[:, None]
        y = np.stack([x[c, skew[c]:skew[c] + n] for c in range(self.channels)])
        full = 6 * np.sqrt(np.mean(np.abs(y) ** 2, axis=1).max()) + 1e-30
        step, top = full / 2 ** (self.bits - 1), 2 ** (self.bits - 1)
        q = lambda a: np.clip(np.round(a / step), -top, top - 1) * step
        return ((q(y.real) + 1j * q(y.imag)) / full).astype(np.complex64), float(full)


class SimBackend:
    """Stands in for the board, the radios and the room."""

    def __init__(self, radio="quad", seed=0, board: BoardModel | None = None, fs=2e6):
        self.rng = np.random.default_rng(seed)
        self.board = board or BoardModel(seed + 1)
        self.ctrl = ControllerModel()
        self.kind = radio
        self.radio = (_Radio(fs, 12, B.RADIO_NF["quad"], [0, 0, 1, 1]) if radio == "quad"
                      else _Radio(fs, 8, B.RADIO_NF["hackrf"], [0]))
        self.fs, self.f, self.t = fs, B.F0, 0.0
        self.sources: list[Source] = []
        self.reflectors: list[Reflector] = []
        self.cal_tx_dbm: float | None = None       # level of the probe signal at the CAL socket, None = off
        self.cal_offset = 100e3
        self.aux_tx_dbm: float | None = None       # tone into the auxiliary patch (radar illuminator)
        self.aux_gain_dbi, self.aux_leak_db = 5.0, -38.0
        self.tx_weights = np.ones(4, complex)
        self.radio.retune(self.rng)

    # ---- the three calls the host library uses --------------------------------
    def controller(self, frame: bytes) -> bytes:
        self.t += self.rng.uniform(0.2e-3, 2e-3)   # USB latency
        self.ctrl.now = self.t
        return self.ctrl.handle(frame)

    def tune(self, f: float):
        self.f = f
        self.radio.retune(self.rng)
        self.t += 5e-3

    def capture(self, n: int) -> Capture:
        skew = self.radio.skew
        m = n + int(skew.max())
        self.t += self.rng.uniform(0.1e-3, 1e-3)
        t = self.t + np.arange(m) / self.fs
        nz = Noise(self.rng, self.fs)
        x = np.zeros((B.N_ELEM, m), complex)
        for s in self.sources:
            u, v = B.direction(s.theta_deg, s.phi_deg)
            x += (np.sqrt(dbm2w(s.power_dbm)) * self.board.elem_resp(u, v, self.f) * s.curvature(self.f))[:, None] * s.waveform(t)
        if self.aux_tx_dbm is not None:
            lam, eirp = B.C0 / self.f, dbm2w(self.aux_tx_dbm + self.aux_gain_dbi)
            tone = np.exp(2j * np.pi * self.cal_offset * t)
            x += (np.sqrt(dbm2w(self.aux_tx_dbm)) * db2amp(self.aux_leak_db) * np.exp(1j * 0.7 * np.arange(B.N_ELEM)))[:, None] * tone
            for r in self.reflectors:
                u, v = B.direction(r.theta_deg, r.phi_deg)
                p_iso = eirp * lam ** 2 * r.rcs_m2 / ((4 * np.pi) ** 3 * r.range_m ** 4)
                echo = tone * np.exp(-4j * np.pi / lam * (r.range_m + r.motion(t)))
                x += (np.sqrt(p_iso) * self.board.elem_resp(u, v, self.f))[:, None] * echo
        widx, words = self.ctrl.word_index(t), self.ctrl.words()
        out = np.zeros((self.radio.channels, m), complex)
        for i in np.unique(widx):
            sel, w = widx == i, words[i]
            xi = x[:, sel]
            if self.cal_tx_dbm is not None:             # probe code; the CAL switch isolates it when off
                level = self.cal_tx_dbm - B.IL_SW - (0 if w.cal else B.ISO_SW)
                xi = xi + (np.sqrt(dbm2w(level)) * self.board.probe)[:, None] * self.cal_wave(t[sel])
            y, h = self.board.to_sockets(w, self.board.receive(w, xi, self.f, nz), nz)
            out[:, sel] = y if self.radio.channels == 4 else h
        edges = np.flatnonzero(np.diff(widx)) + 1              # switching transients: two bad samples
        for e in edges:
            out[:, e:e + 2] *= 0.3
        samples, scale = self.radio.digitise(out, skew, n, nz)
        self.t = float(t[-1])
        return Capture(samples, self.fs, scale)

    CAL_CHIP_RATE, CAL_CHIPS = 500e3, 256

    def _cal_chips(self):
        return np.exp(1j * (np.pi / 4 + np.pi / 2 * np.random.default_rng(99).integers(0, 4, self.CAL_CHIPS)))

    def cal_signal(self, on: bool):
        """Start or stop the pseudo-noise code on the transmit channel that feeds the centre probe.

        Returns one period of the code at the sample rate, the replica the host correlates against.
        """
        self.cal_tx_dbm = -12.0 if on else None
        return np.repeat(self._cal_chips(), int(round(self.fs / self.CAL_CHIP_RATE)))

    def cal_wave(self, t):
        return self._cal_chips()[np.floor(t * self.CAL_CHIP_RATE).astype(np.int64) % self.CAL_CHIPS]

    # ---- virtual test equipment -----------------------------------------------
    def set_tx(self, weights):
        """Baseband weight (amplitude up to 1, and phase) on each of the four transmit channels."""
        self.tx_weights = np.asarray(weights, dtype=complex)

    def eirp_dbm(self, theta_deg, phi_deg, drive_dbm=0.0) -> float:
        """Far-field power meter: EIRP toward a direction for the word now latched."""
        u, v = B.direction(theta_deg, phi_deg)
        amp = np.sqrt(dbm2w(drive_dbm))
        drive = amp if self.kind == "hackrf" else amp * self.tx_weights * self.radio.tx_lo
        a = self.board.transmit(self.ctrl.word, self.f, drive)
        field = np.sum(a * self.board.elem_resp(u, v, self.f))
        return float(w2dbm(abs(field) ** 2 + 1e-30))
