"""High-level driver. It only ever talks to a backend (backend.Backend) through:

    backend.controller(frame) -> reply     beam-controller protocol frames
    backend.tune(freq_hz)                  retune the radio(s)
    backend.capture(n) -> Capture          n I/Q samples per channel
    backend.cal_signal(on) -> replica      four-channel: probe code on the CAL IN transmit channel
    backend.set_tx(weights)                four-channel: baseband transmit weights
    backend.kind, backend.fs               "quad" or "hackrf"; sample rate

so the same code runs against the simulator (sim.SimBackend) and the hardware.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import board as B
from . import dsp
from . import protocol as p
from .board import ControlWord
from .control import phase_states, snapshot_words, state_lag


@dataclass
class Calibration:
    elem_ref: np.ndarray = field(default_factory=lambda: np.ones(B.N_ELEM, complex))   # raw snapshot of a boresight source
    probe_ref: np.ndarray | None = None       # channel phasors against the probe code at that time
    chan: np.ndarray = field(default_factory=lambda: np.ones(4, complex))              # channel change since then

    def quad_phase(self) -> np.ndarray:
        """One phase per quadrant channel, as seen when elem_ref was taken."""
        return np.array([np.exp(1j * np.angle(self.elem_ref[e].sum())) for e in B.QUAD_ELEMS])

    def elem_offsets(self, hackrf: bool) -> np.ndarray:
        """Phase of each element path at state 0, for the analogue beam."""
        ref = self.elem_ref if hackrf else self.elem_ref * np.conj(self.quad_phase())[B.QUAD]
        return np.angle(ref)

    def chan_weights(self) -> np.ndarray:
        """Weights that bring the four quadrant channels into phase."""
        c = self.quad_phase() * self.chan
        return np.conj(c / np.abs(c))


class Array:
    def __init__(self, backend, slot_us: int = 50):
        self.b, self.slot_us = backend, slot_us
        self.hackrf = backend.kind == "hackrf"
        self.slot = int(round(slot_us * 1e-6 * backend.fs))
        self.cal = Calibration()
        self.f, self.skew = B.F0, 0
        self.tx_ref = None                    # (channel weights, quadrant phases) from calibrate_tx()
        if not self.hackrf:
            self.replica = backend.cal_signal(False)  # one period of the probe code, for correlation

    # ---- controller link -------------------------------------------------------
    def _cmd(self, cmd: int, payload: bytes = b"") -> bytes:
        rcmd, status, data = p.parse_reply(self.b.controller(p.request(cmd, payload)))
        if rcmd != cmd or status != p.OK:
            raise RuntimeError("controller refused command 0x%02x (status %d)" % (cmd, status))
        return data

    def ping(self) -> int:
        return self._cmd(p.PING)[0]

    def status(self) -> p.Status:
        """Controller state, reference clock and the current the USB-C source offers."""
        return p.Status.from_payload(self._cmd(p.STATUS))

    def set_word(self, w: ControlWord):
        self._cmd(p.SET_WORD, w.to_bytes())

    def run(self, words: list[ControlWord]):
        self._cmd(p.LOAD_SEQ, p.Sequence(words, self.slot_us).payload())
        self._cmd(p.RUN)

    def stop(self):
        self._cmd(p.STOP)

    def tune(self, f: float):
        """Retune. In four-radio mode the channels are re-aligned straight away."""
        self.f = f
        self.b.tune(f)
        if not self.hackrf and self.cal.probe_ref is not None:
            self.align()

    # ---- captures ----------------------------------------------------------------
    def _capture(self, words: list[ControlWord]) -> np.ndarray:
        """Play a sequence and capture three repetitions of it, in sqrt(W) at the radio input."""
        self.run(words)
        cap = self.b.capture((3 * (len(words) + 1) + 2) * self.slot + 64)
        self.stop()
        return cap.samples.astype(complex) * cap.scale

    def _slots(self, words: list[ControlWord]):
        """Capture a sequence and return the aligned slots of two repetitions."""
        y = self._capture(words)
        if not self.hackrf:
            y = dsp.shift_boards(y, self.skew)
        period = len(words) + 1
        return dsp.cut_slots(y, dsp.find_marker(y, self.slot, period), self.slot, period)

    def _measure_channels(self) -> np.ndarray:
        """Send the probe code and time and phase the four radio channels against it."""
        n_cal = 4
        self.b.cal_signal(True)
        try:
            y = self._capture([ControlWord.receive(0, cal=True)] * n_cal)
        finally:
            self.b.cal_signal(False)
        m = dsp.find_marker(y[:2], self.slot, n_cal + 1)
        self.skew, chan = dsp.cal_slot(y, m + (n_cal + 2) * self.slot, n_cal * self.slot, self.replica)
        return chan

    def raw_snapshot(self, base_state=0, blind: bool = False, probe: bool = False) -> np.ndarray:
        """Uncalibrated complex value per element."""
        words, schedule = snapshot_words(base_state, self.hackrf, blind)
        for w in words:
            w.cal = probe
        seg, idx = self._slots(words)
        if blind:
            return dsp.decode_blind(seg, schedule, self.f, base_state)
        tone = dsp.estimate_tone(seg, idx, self.b.fs)
        return dsp.decode_coherent(dsp.slot_means(seg, idx, self.b.fs, tone).mean(0), self.hackrf, self.f, base_state)

    # ---- calibration ---------------------------------------------------------------
    def calibrate(self, range_m: float | None = None):
        """With a steady source straight ahead: store the element reference.

        In four-radio mode the channel phases are measured first and stored
        with it, so that align() can undo any later retune. range_m: the source's
        distance from the array centre, to remove its wavefront curvature (at 2 m and
        5.8 GHz that is 0.6 / 2.9 / 5.3 degrees on inner / edge / corner elements). The
        error is symmetric about the centre, so it barely biases single-source direction
        finding, but it defocuses beams and moves nulls.
        """
        self.cal = Calibration()
        if not self.hackrf:
            self.cal.probe_ref = self._measure_channels()
        self.cal.elem_ref = self.raw_snapshot()
        if range_m is not None:
            extra = np.sqrt(range_m ** 2 + B.X ** 2 + B.Y ** 2) - range_m          # path beyond the centre's
            self.cal.elem_ref = self.cal.elem_ref * np.exp(2j * np.pi * self.f / B.C0 * extra)

    def align(self) -> np.ndarray:
        """After a retune or a stream restart: re-time and re-phase the four radio channels."""
        c = self._measure_channels() / self.cal.probe_ref
        self.cal.chan = c / np.abs(c)
        return self.cal.chan

    def probe_snapshot(self) -> np.ndarray:
        """Raw per-element response to the centre probe: a health check against drift or a dead part."""
        return self.raw_snapshot(blind=not self.hackrf, probe=True)

    # ---- measurements ----------------------------------------------------------------
    def snapshot(self, blind: bool = False) -> np.ndarray:
        """Calibrated complex amplitude of the incoming field at each element."""
        z = self.raw_snapshot(blind=blind) / self.cal.elem_ref
        return z if self.hackrf else z / self.cal.chan[B.QUAD]

    def locate(self, blind: bool = False) -> tuple[float, float]:
        """(theta_deg, phi_deg) of the strongest source."""
        return B.angles(*dsp.locate(self.snapshot(blind), self.f))

    def _quad_phase(self, state, u, v) -> np.ndarray:
        """Phase each quadrant adds toward (u, v) for these phase states, calibration included."""
        ph = self.cal.elem_offsets(False) - state_lag(state, self.f) + 2 * np.pi * self.f / B.C0 * (B.X * u + B.Y * v)
        r = np.array([np.exp(1j * ph[e]).sum() for e in B.QUAD_ELEMS])
        return r / np.abs(r)

    def point(self, theta_deg: float, phi_deg: float, tx: bool = False) -> ControlWord:
        """Steer the analogue beam; returns the word that was latched.

        For a four-channel transmit beam the baseband phase of each channel is
        set as well, which needs calibrate_tx() to have run since the last retune.
        """
        u, v = B.direction(theta_deg, phi_deg)
        s = phase_states(u, v, self.f, self.cal.elem_offsets(self.hackrf))
        w = ControlWord.transmit(s, self.hackrf) if tx else ControlWord.receive(s, self.hackrf)
        self.set_word(w)
        if tx and not self.hackrf and self.tx_ref is not None:
            weights, ref = self.tx_ref
            self.b.set_tx(weights * ref * np.conj(self._quad_phase(s, u, v)))
        return w

    def calibrate_tx(self, power_meter):
        """Bring the four transmit channels into phase, using a power meter straight ahead.

        power_meter() returns the received power for the weights now set. Each
        channel is turned through three phases against channel 0 and the best
        phase read off the three powers. Repeat after every retune.
        """
        s0 = phase_states(0, 0, self.f, self.cal.elem_offsets(False))
        self.set_word(ControlWord.transmit(s0))
        weights, turns = np.ones(4, complex), np.exp(2j * np.pi * np.arange(3) / 3)
        for q in (1, 2, 3):
            readings = []
            for t in turns:
                trial = np.zeros(4, complex)
                trial[0], trial[q] = 1, t
                self.b.set_tx(trial)
                readings.append(power_meter())
            weights[q] = np.exp(1j * np.angle(np.sum(np.array(readings) * turns)))
        self.tx_ref = (weights, self._quad_phase(s0, 0, 0))
        self.b.set_tx(weights)

    def beam(self, word: ControlWord, n: int = 2000, toward: tuple[float, float] | None = None) -> np.ndarray:
        """Samples of the formed beam for one word: (n,) after combining the channels.

        The word is played as a one-slot sequence so that the marker lets the
        two radio boards be lined up. Four-radio mode: with toward=(theta_deg, phi_deg)
        each quadrant channel is phased toward that direction, which also corrects the
        2-bit rounding left in each quadrant (as point(tx=True) does); without it the
        channels are only aligned for broadside.
        """
        slot_us, self.slot_us = self.slot_us, int(np.ceil((n + 2 * dsp.GUARD) / self.b.fs * 1e6))
        slot, self.slot = self.slot, int(round(self.slot_us * 1e-6 * self.b.fs))
        try:
            seg, _ = self._slots([word])
        finally:
            self.slot_us, self.slot = slot_us, slot
        y = seg[0, 0]
        if self.hackrf:
            return y[0]
        if toward is None:
            return self.cal.chan_weights() @ y / 2
        return np.conj(self.quadrant_response(word.state, *B.direction(*toward))) @ y / 2

    def quadrant_response(self, state, u: float, v: float) -> np.ndarray:
        """Expected phase of each quadrant channel for a plane wave from (u, v) with these phase states.

        It folds in the calibration, the channel phases and the 2-bit rounding of each state.
        """
        r = self._quad_phase(state, u, v) * np.conj(self.cal.chan_weights())
        return r / np.abs(r)

    def monopulse(self, theta_deg: float, phi_deg: float) -> tuple[float, float]:
        """Refine a direction from one burst, using the four quadrants (four-radio mode)."""
        u0, v0 = B.direction(theta_deg, phi_deg)
        s = phase_states(u0, v0, self.f, self.cal.elem_offsets(False))
        seg, _ = self._slots([ControlWord.receive(s)])
        y = seg[0, 0] * np.conj(self.quadrant_response(s, u0, v0))[:, None]
        k = 2 * np.pi * self.f / B.C0 * 2 * B.PITCH
        du = np.angle(np.vdot(y[0] + y[2], y[1] + y[3])) / k      # right half against left half
        dv = np.angle(np.vdot(y[2] + y[3], y[0] + y[1])) / k      # top half against bottom half
        return B.angles(u0 + du, v0 + dv)
