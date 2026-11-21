"""Virtual acceptance tests: the host library against the simulated board.

Each test mirrors a bench check from the design doc's build order.
"""
import struct

import numpy as np
import pytest

from array16 import Array, board as B, dsp, protocol as p
from array16.board import ControlWord
from array16.fwsim import ControllerModel
from array16.control import phase_states
from array16.sim import BoardModel, Reflector, SimBackend, Source

IDEAL = dict(amp_sigma_db=0, phase_sigma_deg=0, stage_sigma_deg=0, bypass_sigma_deg=0, port_sigma_deg=0)


def bench(kind, seed=1, calibrate=True, ideal=False):
    """A board on the virtual bench, calibrated against a source straight ahead."""
    b = SimBackend(kind, seed=seed, board=BoardModel(seed + 1, **IDEAL) if ideal else None)
    a = Array(b)
    if calibrate:
        b.sources = [Source(0, 0, -75)]
        a.calibrate()
        b.sources = []
    return b, a


# ---- control word and protocol ---------------------------------------------------
def test_control_word_round_trip_and_layout():
    rng = np.random.default_rng(0)
    for _ in range(50):
        v = int(rng.integers(0, 2 ** 56, dtype=np.uint64))
        assert ControlWord.from_bytes(ControlWord.from_int(v).to_bytes()).to_int() == v
    assert ControlWord().to_int() == 0                              # power-up state: everything off
    w = ControlWord.receive(3, hackrf_mode=True)
    assert w.to_int() & 0xFFFF == 0xFFFF and (w.to_int() >> B.BIT_MODE) & 1 == 1
    assert list(w.state) == [3] * 16 and len(w.to_bytes()) == 7


def test_controller_protocol_and_timeline():
    c = ControllerModel()
    assert p.parse_reply(c.handle(b"\xa5\x01\x00\x00\x00\x00"))[1] == p.ERR_CRC
    assert p.parse_reply(c.handle(p.request(0x77)))[1] == p.ERR_CMD
    assert p.parse_reply(c.handle(p.request(p.RUN)))[1] == p.ERR_ARG   # nothing loaded yet
    words = [ControlWord.receive(k) for k in range(3)]
    assert p.parse_reply(c.handle(p.request(p.LOAD_SEQ, p.Sequence(words, 50).payload())))[1] == p.OK
    c.now = 1.0
    c.handle(p.request(p.RUN))
    idx = c.word_index(1.0 + np.array([10, 60, 110, 160, 210]) * 1e-6)
    assert list(idx) == [0, 1, 2, 3, 0]                             # marker, three words, marker again
    assert c.words()[0].to_int() == 0 and c.words()[2].to_int() == words[1].to_int()


def test_finite_sequence_ends_like_stop():
    c = ControllerModel()
    static = ControlWord.receive(1)
    c.handle(p.request(p.SET_WORD, static.to_bytes()))
    words = [ControlWord.receive(k) for k in (2, 3)]
    c.handle(p.request(p.LOAD_SEQ, p.Sequence(words, 50, repeat=2).payload()))
    c.now = 1.0
    c.handle(p.request(p.RUN))
    latched = lambda t: [c.words()[i].to_int() for i in c.word_index(np.asarray(t))]
    slots = 1.0 + (np.arange(8) * 50 + 25) * 1e-6
    w0, w1, s = words[0].to_int(), words[1].to_int(), static.to_int()
    assert latched(slots) == [0, w0, w1, 0, w0, w1, s, s]           # two periods, then the static word again
    status = lambda: p.Status.from_payload(p.parse_reply(c.handle(p.request(p.STATUS)))[2])
    assert status().running
    c.now = 1.0 + 6 * 50e-6
    st = status()
    assert not st.running and st.count == 2 and st.slot_us == 50
    assert latched([c.now + 1e-3]) == [s]


def test_status_reports_reference_and_usb_power():
    b = SimBackend("quad")
    a = Array(b)
    b.ctrl.usb_ma, b.ctrl.ref_ok = 1500, False
    a._cmd(p.CLOCK_SRC, b"\x01")
    st = a.status()
    assert (st.version, st.running, st.clock_src, st.ref_ok, st.usb_ma) == (p.FW_VERSION, False, 1, False, 1500)


class Version1Controller:
    """A controller still on protocol version 1: PING says 1, and STATUS has 9 bytes, its count u8."""

    kind, fs = "hackrf", 2e6

    def controller(self, frame: bytes) -> bytes:
        cmd = frame[1]
        if cmd == p.PING:
            return p.reply(cmd, p.OK, bytes([1]))
        if cmd == p.STATUS:
            return p.reply(cmd, p.OK, struct.pack("<BBBHBBH", 1, 0, 3, 50, 0, 1, 1500))
        return p.reply(cmd, p.ERR_ARG)


def test_host_refuses_a_controller_of_another_protocol_version():
    assert Array(SimBackend("hackrf")).check_version() == p.FW_VERSION == 2
    a = Array(Version1Controller())
    assert a.ping() == 1
    with pytest.raises(RuntimeError, match="protocol version 1, this host version 2"):
        a.check_version()
    with pytest.raises(ValueError, match="9-byte STATUS from protocol version 1: this host reads version 2's 10 bytes"):
        a.status()
    with pytest.raises(ValueError, match="0-byte STATUS from a controller"):
        p.Status.from_payload(b"")


# ---- step 2 of the build order: the board with only the HackRF -----------------------
@pytest.mark.parametrize("direction", [(10, 0), (25, 30), (40, 135), (35, -100)])
def test_hackrf_mode_finds_direction(direction):
    b, a = bench("hackrf")
    b.sources = [Source(*direction, power_dbm=-75)]
    assert dsp.angle_between(a.locate(), direction) < 1.5


def test_hackrf_mode_beam_points_where_commanded():
    b, a = bench("hackrf")
    b.sources = [Source(30, 45, -75)]
    power = lambda th, ph: np.mean(np.abs(a.beam(a.point(th, ph))) ** 2)
    on, off = power(30, 45), power(30, 225)
    assert 10 * np.log10(on / off) > 8                              # beam on the source against 60 degrees away
    cut = [power(th, 45) for th in range(0, 61, 5)]
    assert abs(5 * int(np.argmax(cut)) - 30) <= 5


def test_calibration_source_distance_is_corrected():
    """A boresight source 1.5 m away bends the wavefront by up to 7 degrees across the aperture
    (corners); calibrate(range_m) gives back the reference a far source would. The error is
    symmetric about the centre, so it hardly biases single-source DF, but it distorts beams."""
    refs = {}
    for label, rng, kw in (("far", None, {}), ("near", 1.5, {}), ("corrected", 1.5, {"range_m": 1.5})):
        b, a = bench("hackrf", seed=12, calibrate=False)            # same seed: the same noise in all three
        b.sources = [Source(0, 0, -75, range_m=rng)]
        a.calibrate(**kw)
        refs[label] = a.cal.elem_ref / a.cal.elem_ref[5]
    worst = lambda k: np.degrees(np.max(np.abs(np.angle(refs[k] / refs["far"]))))
    assert 5 < worst("near") < 9
    assert worst("corrected") < 1.5                                 # what is left is 8-bit quantisation of a different signal


def test_four_radio_receive_beam_phased_toward_its_direction():
    gains = []
    for seed in (1, 2, 3, 8):
        b, a = bench("quad", seed=seed)
        d = (20 + 2 * seed, (35 * seed) % 360)
        b.sources = [Source(*d, -70)]
        w = a.point(*d)
        plain = np.mean(np.abs(a.beam(w)) ** 2)
        toward = np.mean(np.abs(a.beam(w, toward=d)) ** 2)
        gains.append(10 * np.log10(toward / plain))
    assert min(gains) > -0.05 and np.mean(gains) > 0.05              # corrects each quadrant's 2-bit residue


# ---- step 3: four radio channels -------------------------------------------------------
def test_four_radio_snapshot_survives_a_retune():
    b, a = bench("quad")
    b.sources = [Source(20, 60, -75)]
    assert dsp.angle_between(a.locate(), (20, 60)) < 1.5
    a.tune(B.F0)                                                    # both radio boards relock; tune() re-aligns
    assert np.max(np.abs(np.angle(a.cal.chan / a.cal.chan[0]))) > np.radians(20)   # the retune really moved a channel
    assert dsp.angle_between(a.locate(), (20, 60)) < 1.5


def test_alignment_holds_with_a_strong_signal_on_air():
    b, a = bench("quad", seed=5)
    b.sources = [Source(0, 0, -60)]                                 # present while the channels are re-aligned
    a.tune(B.F0)
    z = a.snapshot()
    assert np.degrees(np.std(np.angle(z * np.conj(z.mean())))) < 3


def test_lost_stream_alignment_is_caught_by_realigning():
    b, a = bench("quad", seed=6)
    b.sources = [Source(25, 10, -70, kind="qpsk")]
    b.radio.restart(b.rng)                                          # a dropped buffer: the boards slip apart
    a.align()
    assert dsp.angle_between(a.locate(blind=True), (25, 10)) < 2


def test_blind_snapshot_handles_a_modulated_signal():
    b, a = bench("quad", seed=2)
    b.sources = [Source(30, -45, -70, kind="qpsk")]
    assert dsp.angle_between(a.locate(blind=True), (30, -45)) < 2


def test_single_burst_direction_from_four_quadrants():
    b, a = bench("quad", seed=4)
    b.sources = [Source(24, 20, -70, kind="qpsk")]
    assert dsp.angle_between(a.monopulse(20, 0), (24, 20)) < 2.5


# ---- off the design frequency (as if the RF parts were wideband) ----------------------------
def test_delay_line_states_scale_with_frequency():
    from array16.control import state_lag, toggle_factor
    assert np.allclose(toggle_factor(0, B.F0), -1) and np.allclose(toggle_factor(2, B.F0), -1)
    assert np.allclose(toggle_factor(0, B.F0 / 2), -1j)                # the 180-degree line is 90 degrees at half frequency
    assert np.allclose(state_lag(np.arange(4), 2.9e9), np.arange(4) * np.pi / 4)


@pytest.mark.parametrize("kind", ["quad", "hackrf"])
def test_snapshots_decode_with_the_actual_toggle_off_frequency(kind):
    b, a = bench(kind, seed=9, calibrate=False)
    a.tune(3.0e9)
    b.sources = [Source(0, 0, -75)]
    a.calibrate()                                                   # calibration is per frequency
    b.sources = [Source(25, 70, -75)]
    assert dsp.angle_between(a.locate(), (25, 70)) < 2.5


def test_beam_steering_uses_the_delay_lines_actual_phases():
    b, a = bench("hackrf", seed=7, calibrate=False)
    a.tune(4.2e9)
    b.sources = [Source(0, 0, -75)]
    a.calibrate()
    b.sources = []
    a.point(20, 0, tx=True)
    cut = [b.eirp_dbm(th, 0) for th in range(0, 61, 2)]
    assert abs(2 * int(np.argmax(cut)) - 20) <= 6


# ---- link budget against the design doc ---------------------------------------------------
def test_noise_figure_of_the_four_radio_chain():
    b, a = bench("quad", calibrate=False, ideal=True)
    p_dbm = -110.0
    b.sources = [Source(0, 0, p_dbm, f_offset=20e3)]
    # the two radio boards have unknown phases here, so combine the channels on the tone
    seg, _ = a._slots([ControlWord.receive(0)] * 8)
    ch = seg[0].transpose(1, 0, 2).reshape(4, -1)
    tone = np.exp(-2j * np.pi * 20e3 / b.fs * np.arange(seg.shape[-1]))
    amp = np.array([[np.mean(s[c] * tone) for c in range(4)] for s in seg[0]])
    sig = np.mean(np.abs(amp) ** 2, axis=0)
    noise = np.mean(np.abs(ch) ** 2, axis=1) - sig
    snr = 4 * np.mean(sig / noise)                                  # four channels combined in phase
    kTB = 10 * np.log10(B.K_BOLTZ * B.T0 * b.fs) + 30
    nf = p_dbm + B.ARRAY_GAIN_DBI - kTB - 10 * np.log10(snr)
    assert 2.0 < nf < 3.0, nf


def test_transmit_eirp_matches_the_budget():
    b, a = bench("quad", calibrate=False, ideal=True)
    a.calibrate_tx(lambda: 10 ** (b.eirp_dbm(0, 0) / 10))           # the two radio boards' transmit phases
    a.point(0, 0, tx=True)
    assert abs(b.eirp_dbm(0, 0) - 14.7) < 0.3
    b, a = bench("hackrf", calibrate=False, ideal=True)
    a.set_word(ControlWord.transmit(0, hackrf_mode=True))
    assert abs(b.eirp_dbm(0, 0) - 4.4) < 0.3


def test_transmit_beam_steers_in_hackrf_mode():
    b, a = bench("hackrf", seed=7)
    a.point(30, 0, tx=True)
    cut = [b.eirp_dbm(th, 0) for th in range(0, 61, 2)]
    assert abs(2 * int(np.argmax(cut)) - 30) <= 6
    assert max(cut) > 4.4 - 3.5                                     # scan loss, 2-bit loss and unit errors


def test_four_channel_transmit_needs_and_gets_its_own_alignment():
    b, a = bench("quad", seed=7)
    a.point(0, 0, tx=True)
    before = b.eirp_dbm(0, 0)                                       # four transmit channels at unknown phases
    a.calibrate_tx(lambda: 10 ** (b.eirp_dbm(0, 0) / 10))
    a.point(0, 0, tx=True)
    assert b.eirp_dbm(0, 0) > 14.7 - 1.5 and b.eirp_dbm(0, 0) > before + 0.5
    a.point(30, 0, tx=True)
    cut = [b.eirp_dbm(th, 0) for th in range(0, 61, 2)]
    assert abs(2 * int(np.argmax(cut)) - 30) <= 6 and max(cut) > 14.7 - 3.5


# ---- radar practice: motion from the phase of an echo -----------------------------------------
def test_millimetre_motion_shows_in_the_echo_phase():
    b, a = bench("quad", seed=8)
    lam, f_tone = B.C0 / B.F0, b.cal_offset
    u, v = B.direction(15, 0)
    s = phase_states(u, v, a.f, a.cal.elem_offsets(False))
    b.aux_tx_dbm = 0.0                                              # steady tone from the auxiliary patch

    def phasors():                                                  # tone amplitude on each quadrant channel
        seg, idx = a._slots([ControlWord.receive(s)])
        return dsp.slot_means(seg, idx, b.fs, f_tone)[0, 0]

    leak = phasors()                                                # empty room: direct leakage only
    target = a.quadrant_response(s, u, v)
    b.reflectors = [Reflector(15, 0, 2.0, rcs_m2=0.5, motion=lambda t: 2e-3 * np.sin(2 * np.pi * 0.5 * t))]
    basis, times, ratio = np.stack([leak, target], axis=1), [], []
    for _ in range(40):                                             # two seconds, one capture every 50 ms
        b.t += 0.05
        coef = np.linalg.lstsq(basis, phasors(), rcond=None)[0]     # split each capture into leakage and echo
        times.append(b.t)
        ratio.append(coef[1] / coef[0])                             # the leakage is the phase reference
    ph = np.unwrap(np.angle(ratio))
    moved = -(ph - ph.mean()) * lam / (4 * np.pi)                   # metres along range
    truth = 2e-3 * np.sin(2 * np.pi * 0.5 * np.array(times))
    assert np.sqrt(np.mean((moved - (truth - truth.mean())) ** 2)) < 0.3e-3
