"""The firmware's C core against its specification, fwsim.ControllerModel.

The core (firmware/core/ctrl.c) is built for this machine and driven through ctypes. A
fake board stands in for the RP2040 layer: it records what the core tells the hardware.
Every frame goes to both the core and the model; replies must match byte for byte, and
the hardware must have been told to latch or play what the model says is latched.
"""
import ctypes
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from array16 import protocol as p
from array16.board import ControlWord
from array16.fwsim import ControllerModel

CORE = Path(__file__).resolve().parents[1] / "firmware" / "core"
U8P = ctypes.POINTER(ctypes.c_uint8)
LATCH = ctypes.CFUNCTYPE(None, ctypes.c_void_p, U8P)
RUN = ctypes.CFUNCTYPE(None, ctypes.c_void_p, U8P, ctypes.c_uint8, ctypes.c_uint16, ctypes.c_uint16, U8P)
PLAYING = ctypes.CFUNCTYPE(ctypes.c_bool, ctypes.c_void_p)
CLOCK = ctypes.CFUNCTYPE(ctypes.c_bool, ctypes.c_void_p, ctypes.c_uint8)
USB_MA = ctypes.CFUNCTYPE(ctypes.c_uint16, ctypes.c_void_p)


class HW(ctypes.Structure):
    _fields_ = [("latch", LATCH), ("run", RUN), ("playing", PLAYING), ("clock", CLOCK), ("usb_ma", USB_MA),
                ("ctx", ctypes.c_void_p)]


@pytest.fixture(scope="module")
def lib(tmp_path_factory):
    cc = shutil.which("cc") or shutil.which("gcc") or shutil.which("clang")
    if cc is None:
        pytest.skip("no C compiler")
    so = tmp_path_factory.mktemp("core") / "libctrl.so"
    subprocess.run([cc, "-std=c11", "-O2", "-Wall", "-Wextra", "-Wconversion", "-Werror", "-shared", "-fPIC",
                    "-o", str(so), str(CORE / "ctrl.c")], check=True)
    L = ctypes.CDLL(str(so))
    L.ctrl_size.restype = ctypes.c_size_t
    L.ctrl_init.argtypes = [ctypes.c_void_p, ctypes.POINTER(HW)]
    L.ctrl_handle.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_size_t, ctypes.c_void_p]
    L.ctrl_handle.restype = ctypes.c_size_t
    L.ctrl_feed.argtypes = [ctypes.c_void_p, ctypes.c_uint8, ctypes.c_uint32, ctypes.c_void_p]
    L.ctrl_feed.restype = ctypes.c_size_t
    L.ctrl_pack.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_uint32)]
    return L


class FakeBoard:
    """The RP2040 layer as the core sees it. Time, the reference clock and USB power come from the model."""

    def __init__(self, lib, model: ControllerModel):
        self.m, self.lib = model, lib
        self.mode, self.static, self.run_args = "static", bytes(7), None
        self.hw = HW(LATCH(self._latch), RUN(self._run), PLAYING(lambda ctx: self.playing()),
                     CLOCK(lambda ctx, src: self.m.ref_ok), USB_MA(lambda ctx: self.m.usb_ma), None)
        self.ctrl = ctypes.create_string_buffer(lib.ctrl_size())
        self.reply = ctypes.create_string_buffer(p.MAX_PAYLOAD)
        lib.ctrl_init(self.ctrl, ctypes.byref(self.hw))

    def _latch(self, ctx, word):
        self.mode, self.static = "static", bytes(word[:7])

    def _run(self, ctx, words, n, slot_us, repeat, after):
        self.mode = "run"
        self.run_args = dict(table=[bytes(7)] + [bytes(words[7 * i:7 * i + 7]) for i in range(n)],
                             slot_us=slot_us, repeat=repeat, after=bytes(after[:7]), t0=self.m.now)

    def playing(self) -> bool:
        if self.mode != "run":
            return False
        a = self.run_args
        period = len(a["table"]) * a["slot_us"] * 1e-6
        return not (a["repeat"] and self.m.now - a["t0"] >= a["repeat"] * period - 1e-12)

    def handle(self, frame: bytes) -> bytes:
        n = self.lib.ctrl_handle(self.ctrl, frame, len(frame), self.reply)
        return self.reply.raw[:n]

    def feed(self, data: bytes, now_ms: int) -> bytes:
        out = b""
        for b in data:
            n = self.lib.ctrl_feed(self.ctrl, b, now_ms, self.reply)
            out += self.reply.raw[:n]
        return out

    def check_against(self, m: ControllerModel):
        """The hardware holds what the model says is latched."""
        if m.running:
            a = self.run_args
            assert self.mode == "run" and self.playing()
            assert a["table"] == [w.to_bytes() for w in m.words()[:-1]]
            assert (a["after"], a["slot_us"], a["repeat"], a["t0"]) == (m.word.to_bytes(), m.seq.slot_us, m.seq.repeat, m.t_start)
        else:
            held = self.static if self.mode == "static" else self.run_args["after"]
            assert self.mode == "static" or not self.playing()
            assert held == m.words()[0].to_bytes()


def rand_word(rng) -> ControlWord:
    return ControlWord.from_int(int(rng.integers(0, 2 ** 56, dtype=np.uint64)))


def rand_frame(rng) -> bytes:
    k = int(rng.integers(0, 13))
    if k == 0:
        return p.request(p.PING, rng.bytes(int(rng.integers(0, 3))))
    if k == 1:
        return p.request(p.SET_WORD, rand_word(rng).to_bytes())
    if k == 2:
        words = [rand_word(rng) for _ in range(int(rng.integers(1, 6)))]
        return p.request(p.LOAD_SEQ, p.Sequence(words, int(rng.choice([10, 37, 50, 100])), int(rng.choice([0, 0, 1, 2, 3]))).payload())
    if k in (3, 4):
        return p.request(p.RUN)
    if k == 5:
        return p.request(p.STOP)
    if k == 6:
        return p.request(p.STATUS)
    if k == 7:
        return p.request(p.CLOCK_SRC, bytes([int(rng.integers(0, 3))]) * int(rng.choice([1, 1, 1, 0, 2])))
    if k == 8:                                                     # any command code, junk payload
        return p.request(int(rng.integers(0, 256)), rng.bytes(int(rng.integers(0, 9))))
    if k == 9:                                                     # malformed LOAD_SEQ
        count = int(rng.choice([0, 1, 2, 64, 65]))
        head = int(rng.choice([5, 9, 10])).to_bytes(2, "little") + bytes(2) + bytes([count])
        body = rng.bytes(7 * count + int(rng.choice([0, 0, -1, 1])))
        return p.request(p.LOAD_SEQ, (head + body)[:int(rng.integers(0, 6))] if rng.random() < 0.2 else head + body)
    if k == 10:
        return p.request(p.SET_WORD, rng.bytes(int(rng.choice([0, 6, 8]))))
    f = bytearray(p.request(p.STATUS if k == 11 else p.PING))     # damaged on the wire
    i = int(rng.integers(0, len(f)))
    f[i] ^= 1 << int(rng.integers(0, 8))
    return bytes(f) if rng.random() < 0.7 else bytes(f[:int(rng.integers(0, len(f)))])


@pytest.mark.parametrize("seed", range(4))
def test_core_answers_every_frame_like_the_model(lib, seed):
    rng = np.random.default_rng(seed)
    m = ControllerModel()
    m.ref_ok, m.usb_ma = bool(seed % 2), [0, 500, 1500, 3000][seed]
    board = FakeBoard(lib, m)
    for _ in range(1500):
        m.now += float(rng.uniform(0, 300e-6))                     # finite sequences end between some frames
        frame = rand_frame(rng)
        assert board.handle(frame) == m.handle(frame), frame.hex()
        board.check_against(m)


def test_stream_reassembles_split_and_noisy_frames(lib):
    rng = np.random.default_rng(10)
    m, ref = ControllerModel(), ControllerModel()
    board = FakeBoard(lib, m)
    frames = []
    for _ in range(300):
        f = rand_frame(rng)
        if p.parse_request(f) is not None and len(f) - 6 <= p.MAX_PAYLOAD:   # whole, undamaged, not oversized
            frames.append(f)
    stream, want = b"", b""
    for f in frames:
        noise = bytes(b for b in rng.bytes(int(rng.integers(0, 6))) if b != p.SOF_REQ)
        if rng.random() < 0.2:
            noise += bytes([p.SOF_REQ, 0x10, 0xFF, 0xFF])         # a start byte followed by an impossible length
        stream += noise + f
        want += ref.handle(f)
    got, now, i = b"", 0, 0
    while i < len(stream):                                         # arbitrary USB packet boundaries, short gaps
        n = int(rng.integers(1, 70))
        got += board.feed(stream[i:i + n], now)
        i, now = i + n, now + int(rng.integers(0, 20))
    assert got == want


def test_stream_drops_a_stale_partial_frame(lib):
    board = FakeBoard(lib, ControllerModel())
    ping, status = p.request(p.PING), p.request(p.STATUS)
    assert board.feed(ping[:3], 0) == b""
    assert board.feed(status, 51) == ControllerModel().handle(status)          # the half PING was dropped
    assert board.feed(ping[:3], 100) + board.feed(ping[3:], 140) == ControllerModel().handle(ping)
    bad = bytearray(ping)
    bad[-1] ^= 0xFF
    assert p.parse_reply(board.feed(bytes(bad), 200))[:2] == (0, p.ERR_CRC)
    assert board.feed(ping, 201) == ControllerModel().handle(ping)
    # Longer than any request: dropped without a reply (fwsim, given it whole, would say ERR_ARG)
    huge = p.request(p.LOAD_SEQ, bytes([50, 0, 0, 0, 65]) + bytes(7 * 65))
    assert len(huge) - 6 > p.MAX_PAYLOAD and board.feed(huge + ping, 300) == ControllerModel().handle(ping)


def test_pack_puts_every_bit_on_its_register_output(lib):
    rng = np.random.default_rng(3)
    for _ in range(200):
        v = int(rng.integers(0, 2 ** 56, dtype=np.uint64))
        fifo = (ctypes.c_uint32 * 2)()
        lib.ctrl_pack(v.to_bytes(7, "little"), fifo)
        assert fifo[1] & 0xFF == 0
        # The PIO shifter sends 32 bits of the first FIFO word, then 24 of the second, MSB first,
        # on SR_DATA, one SR_CLK edge each. On every edge each 74HC595 moves QA..QH up by one and
        # takes the previous register's QH (via QH') into QA; register 0 takes SR_DATA.
        bits = [(fifo[0] >> (31 - i)) & 1 for i in range(32)] + [(fifo[1] >> (31 - i)) & 1 for i in range(24)]
        chain = [0] * 56                                           # chain[8 * k + j] = output Q(j) of register k
        for b in bits:
            chain = [b] + chain[:-1]
        assert chain == [(v >> i) & 1 for i in range(56)]          # board.py: bit i on Q(i % 8) of register i // 8
