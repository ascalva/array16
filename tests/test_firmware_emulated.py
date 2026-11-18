"""The compiled RP2040 firmware, run in an emulator, against its specification fwsim.ControllerModel.

firmware/emu/harness.mjs boots the firmware image in rp2040js (Wokwi's RP2040 emulator), plays the
scripts below over emulated USB CDC and records the 74HC595 pins with emulated timestamps. Every
frame comes from protocol.py and also goes through the model, at the moment the firmware received
it: replies must match byte for byte, every word latched into the register chain must be the word
the model holds at that moment, and the slot edges of a running sequence must fall on the model's
slot boundaries, counted from the firmware's marker latch.

Known, documented differences, which the checks allow for:
  - RUN: the marker latches 11-12 us after the firmware starts the sequence (12 edges of the 1 MHz
    slot clock; firmware/README.md), not at once. TRIG rises at RUN and falls within about 1 us of
    the marker.
  - SET_WORD, STOP and LOAD_SEQ while running latch about 7 us after the frame is handled (one shift).
  - CLOCK_SRC while a sequence runs may move later slot edges by up to one slot-clock period.
  - If the reference timing the slots goes away, no slot edges come until the firmware notices
    (a few ms at most); then the sequence goes on from the crystal where it stopped.
  - The model gets each frame when the device has it; the firmware may handle it up to 1 ms later
    (stdio_usb puts off USB work while it is busy sending a reply). The checks refuse a script in
    which that could change an answer, such as a STATUS at the end of a finite sequence.

The emulator's slot clock is exact, so slot edges must land within SLOT_TOL_NS of the model's.
What the emulator does not show is in firmware/emu/README.md.

Skips unless node, rp2040js (~/opt/src/rp2040js, or $RP2040JS, with npm install done) and a firmware
image (firmware/build/array16_controller.uf2, or $ARRAY16_FW_IMAGE) are all present.
"""
from __future__ import annotations

import copy
import json
import math
import os
import shutil
import struct
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest

from array16 import protocol as p
from array16.board import ControlWord
from array16.fwsim import ControllerModel

ROOT = Path(__file__).resolve().parents[1]
HARNESS = ROOT / "firmware" / "emu" / "harness.mjs"
RP2040JS = Path(os.environ.get("RP2040JS") or Path.home() / "opt" / "src" / "rp2040js")
IMAGE = Path(os.environ.get("ARRAY16_FW_IMAGE") or ROOT / "firmware" / "build" / "array16_controller.uf2")
NODE = shutil.which("node")


def _missing() -> str | None:
    if NODE is None:
        return "node not found"
    if not (RP2040JS / "src" / "rp2040.ts").exists() or not (RP2040JS / "node_modules" / "tsx").exists():
        return f"rp2040js not installed in {RP2040JS} (git clone https://github.com/wokwi/rp2040js; npm install)"
    if not IMAGE.is_file():
        return f"no firmware image at {IMAGE} (build firmware/, or set ARRAY16_FW_IMAGE)"
    return None


pytestmark = pytest.mark.skipif(_missing() is not None, reason=str(_missing()))

SLOT_TOL_NS = 50                   # the emulated slot clock is exact; the PIO reacts within a few 8 ns cycles
MARKER_AFTER_TRIG_NS = (11_000, 12_200)
CLOCK_SWITCH_SHIFT_NS = 1_000      # CLOCK_SRC while running: one slot-clock period at most
REF_LOSS_GAP_NS = 5_000_000        # reference gone: slots stop until the firmware moves to the crystal
ENV = dict(ref40_hz=40e6, ref10_hz=10e6, cc1_mv=1680, cc2_mv=0)
# USB Type-C sink: voltage across the 5.1k CC pull-down for each advertised current (vRd ranges)
CC_RANGES = ((0, 150, 0), (250, 610, 500), (700, 1160, 1500), (1310, 2040, 3000))
CC_MV = {0: 0, 500: 420, 1500: 930, 3000: 1680}


def advertised_ma(env: dict) -> int:
    v = max(env["cc1_mv"], env["cc2_mv"])
    for lo, hi, ma in CC_RANGES:
        if lo <= v <= hi:
            return ma
    raise ValueError(f"{v} mV lies between Type-C ranges: pick an unambiguous voltage")


def ref_ok(env: dict, src: int) -> bool:
    """The selected reference is present and within 1 % (firmware/README.md)."""
    hz, nominal = (env["ref40_hz"], 40e6) if src == 0 else (env["ref10_hz"], 10e6)
    return abs(hz - nominal) < 0.01 * nominal


# ---- scripts ---------------------------------------------------------------------------------------
class Script:
    """Steps for the harness, plus the request frames the firmware must act on (in order)."""

    def __init__(self, name: str, **env):
        self.name, self.env = name, {**ENV, **env}
        self.steps: list[dict] = []
        self.frames: list[tuple[int, bytes]] = []   # (step that completes the frame, frame)

    def send(self, data: bytes, frames: list[bytes] | None = None, packets: list[int] | None = None,
             gap_us: int | None = None) -> "Script":
        self.frames += [(len(self.steps), f) for f in ([data] if frames is None else frames)]
        step = {"send": data.hex(), "wait_replies": len(self.frames), "timeout_us": 50_000}
        if packets is not None:
            step["packets"] = packets
        if gap_us is not None:
            step["gap_us"] = gap_us
        self.steps.append(step)
        return self

    def req(self, cmd: int, payload: bytes = b"", **kw) -> "Script":
        return self.send(p.request(cmd, payload), **kw)

    def wait(self, us: int) -> "Script":
        self.steps.append({"wait_us": us})
        return self

    def set_env(self, **env) -> "Script":
        self.steps.append({"env": env})
        return self

    def to_json(self) -> dict:
        return {"name": self.name, "env": self.env, "steps": self.steps}


def words(rng, n: int) -> list[ControlWord]:
    return [ControlWord.from_int(int(rng.integers(0, 2 ** 56, dtype=np.uint64))) for _ in range(n)]


def seq(ws, slot_us: int, repeat: int) -> bytes:
    return p.Sequence(ws, slot_us, repeat).payload()


def scripts() -> list[Script]:
    rng = np.random.default_rng(2026)
    out = []

    s = Script("boot")
    s.wait(5000).req(p.STATUS).req(p.PING)
    out.append(s)

    s = Script("set_word")
    a, b, c = words(rng, 3)
    s.req(p.SET_WORD, a.to_bytes()).req(p.STATUS).req(p.SET_WORD, b.to_bytes())
    s.req(p.SET_WORD, c.to_bytes()[:6]).req(p.SET_WORD, c.to_bytes() + b"\0")       # wrong length: ERR_ARG, no latch
    s.req(p.SET_WORD, ControlWord.from_int(2 ** 56 - 1).to_bytes()).req(p.SET_WORD, bytes(7))
    s.req(p.SET_WORD, c.to_bytes()).req(p.PING, b"\x01\x02").req(p.STOP)            # STOP when idle: the static word again
    out.append(s)

    s = Script("status_power", cc1_mv=0)
    for ma in (0, 500, 1500, 3000):
        s.set_env(cc1_mv=CC_MV[ma], cc2_mv=0).req(p.STATUS)
    for ma in (500, 1500, 3000):                                                    # cable the other way round
        s.set_env(cc1_mv=0, cc2_mv=CC_MV[ma]).req(p.STATUS)
    out.append(s)

    s = Script("seq_continuous")
    static, w5 = words(rng, 1)[0], words(rng, 5)
    s.req(p.SET_WORD, static.to_bytes()).req(p.LOAD_SEQ, seq(w5, 50, 0)).req(p.STATUS).req(p.RUN)
    s.wait(1600).req(p.STATUS).wait(200).req(p.STOP).req(p.STATUS).wait(300)
    out.append(s)

    s = Script("seq_finite")
    static, w3, w1 = words(rng, 1)[0], words(rng, 3), words(rng, 1)
    s.req(p.SET_WORD, static.to_bytes()).req(p.LOAD_SEQ, seq(w3, 20, 2)).req(p.RUN).wait(1000).req(p.STATUS)
    # STATUS straight after RUN: the firmware may take up to 1 ms to see a frame that arrives while it is
    # still sending a reply (stdio_usb defers tud_task to a 1 ms timer then), so the sequence runs 3 ms.
    s.req(p.LOAD_SEQ, seq(w1, 500, 3)).req(p.RUN).req(p.STATUS).wait(4000).req(p.STATUS)
    s.req(p.LOAD_SEQ, seq(w3, 10, 1)).req(p.RUN).wait(500).req(p.STATUS).wait(200)
    out.append(s)

    s = Script("seq_max")
    static, w64, v64 = words(rng, 1)[0], words(rng, 64), words(rng, 64)
    s.req(p.SET_WORD, static.to_bytes()).req(p.LOAD_SEQ, seq(w64, 10, 1)).req(p.RUN).wait(1000).req(p.STATUS)
    s.req(p.LOAD_SEQ, seq(v64, 10, 0)).req(p.RUN)
    for _ in range(4):                                                               # USB traffic while 10 us slots play
        s.req(p.STATUS).req(p.PING)
    s.wait(500).req(p.STOP)
    out.append(s)

    s = Script("run_stop_while_running")
    a, b = words(rng, 2)
    w4, w2 = words(rng, 4), words(rng, 2)
    s.req(p.SET_WORD, a.to_bytes()).req(p.LOAD_SEQ, seq(w4, 30, 0)).req(p.RUN).wait(400)
    s.req(p.RUN).wait(400)                                       # restart: the marker, then the first word again
    s.req(p.STOP).req(p.STATUS).req(p.RUN).wait(200)
    s.req(p.LOAD_SEQ, seq(w2, 40, 0))                            # while running: stops, the static word returns
    s.req(p.STATUS).req(p.RUN).wait(300)
    s.req(p.SET_WORD, b.to_bytes()).req(p.STATUS).req(p.STOP).wait(300)
    out.append(s)

    s = Script("clock_src")
    w3 = words(rng, 3)
    s.req(p.STATUS).req(p.CLOCK_SRC, b"\x01").req(p.STATUS)
    s.req(p.LOAD_SEQ, seq(w3, 25, 0)).req(p.RUN).wait(500)
    s.req(p.CLOCK_SRC, b"\x00").wait(500).req(p.STATUS).req(p.STOP)          # switched while running
    s.set_env(ref10_hz=0).req(p.CLOCK_SRC, b"\x01").req(p.STATUS)              # HackRF clock missing: crystal
    s.req(p.RUN).wait(300).req(p.STOP)                                         # sequences still play
    s.set_env(ref40_hz=40.2e6).req(p.CLOCK_SRC, b"\x00").req(p.STATUS)         # 0.5 % off: accepted
    s.set_env(ref40_hz=40.6e6).req(p.CLOCK_SRC, b"\x00").req(p.STATUS)         # 1.5 % off: refused
    s.set_env(ref40_hz=40e6).req(p.CLOCK_SRC, b"\x00").req(p.STATUS)
    s.req(p.CLOCK_SRC, b"\x02").req(p.CLOCK_SRC).req(p.CLOCK_SRC, b"\x00\x00").req(p.STATUS)
    out.append(s)

    s = Script("boot_without_reference", ref40_hz=0)
    w2 = words(rng, 2)
    s.req(p.STATUS).req(p.LOAD_SEQ, seq(w2, 20, 0)).req(p.RUN).wait(300).req(p.STOP)
    out.append(s)

    s = Script("reference_lost")
    w3 = words(rng, 3)
    s.req(p.CLOCK_SRC, b"\x01").req(p.LOAD_SEQ, seq(w3, 20, 0)).req(p.RUN).wait(300)
    s.set_env(ref10_hz=0).wait(6000)                   # the HackRF clock goes away while the sequence plays
    s.req(p.STATUS).req(p.STOP)                        # ref_ok = 0, and the slots run on from the crystal
    s.req(p.LOAD_SEQ, seq(w3, 20, 2)).req(p.RUN).wait(500).req(p.STATUS)
    s.set_env(ref10_hz=10e6).req(p.CLOCK_SRC, b"\x01").req(p.STATUS)   # back: CLOCK_SRC locks on again
    out.append(s)

    # The RP2040's frequency counter may wait for ever for the first edge of a clock that is not
    # running (FC0_STATUS.WAITING). The firmware must not hang on it at boot or at CLOCK_SRC.
    s = Script("no_reference_counter_waits", ref40_hz=0, ref10_hz=0, fc_dead_clock="wait")
    w2 = words(rng, 2)
    s.req(p.STATUS).req(p.LOAD_SEQ, seq(w2, 20, 0)).req(p.RUN).wait(300).req(p.STOP)
    s.req(p.CLOCK_SRC, b"\x01").req(p.STATUS)
    s.set_env(ref40_hz=40e6).req(p.CLOCK_SRC, b"\x00").req(p.STATUS).req(p.RUN).wait(300)
    s.set_env(ref40_hz=0).wait(6000).req(p.STATUS).req(p.STOP)
    out.append(s)

    s = Script("crc_error")
    w = words(rng, 1)[0]
    bad = bytearray(p.request(p.PING))
    bad[-1] ^= 0x40                                              # CRC damaged
    s.send(bytes(bad))
    bad = bytearray(p.request(p.SET_WORD, w.to_bytes()))
    bad[6] ^= 0x01                                               # payload damaged: no latch
    s.send(bytes(bad)).req(p.STATUS).req(p.SET_WORD, w.to_bytes())
    out.append(s)

    s = Script("oversized_length")
    ping, status = p.request(p.PING), p.request(p.STATUS)
    noise = bytes(x for x in rng.bytes(40) if x != p.SOF_REQ)
    huge = p.request(p.LOAD_SEQ, bytes([50, 0, 0, 0, 65]) + bytes(7 * 65))          # 460-byte payload
    over = p.request(p.LOAD_SEQ, bytes([50, 0, 0, 0, 64]) + bytes(7 * 64 + 1))      # one byte over the limit
    for f in (huge, over):
        assert len(f) - 6 > p.MAX_PAYLOAD and p.SOF_REQ not in f[1:]                 # its remains are plain noise
    s.send(bytes([p.SOF_REQ, p.SET_WORD, 0xFF, 0xFF]) + ping, frames=[ping])
    s.send(noise + status, frames=[status])
    s.send(huge + ping, frames=[ping])
    s.send(over, frames=[]).wait(2000).req(p.PING)
    s.req(p.LOAD_SEQ, seq(words(rng, 64), 50, 0))                                    # the longest legal frame
    out.append(s)

    s = Script("split_frame")
    w5, w = words(rng, 5), words(rng, 1)[0]
    f = p.request(p.LOAD_SEQ, seq(w5, 40, 0))
    s.send(f, packets=[1, 2, 1, 10, len(f) - 14], gap_us=2000)
    g = p.request(p.SET_WORD, w.to_bytes())
    s.send(g, packets=[4, 9], gap_us=45_000)                     # 45 ms between packets: still one frame
    s.send(p.request(p.PING), packets=[3, 3])                    # two packets back to back
    s.req(p.RUN).wait(500).req(p.STOP)                           # the reassembled sequence plays
    out.append(s)

    s = Script("partial_frame_timeout")
    ping, status = p.request(p.PING), p.request(p.STATUS)
    w = words(rng, 1)[0]
    sw = p.request(p.SET_WORD, w.to_bytes())
    s.send(ping[:3], frames=[]).wait(60_000).send(status)       # stale half PING dropped: only STATUS answered
    s.send(ping[:3], frames=[]).wait(30_000).send(ping[3:], frames=[ping])
    s.send(sw[:8], frames=[]).wait(70_000).req(p.PING)          # stale SET_WORD: never latched
    out.append(s)

    s = Script("invalid_requests")
    w2 = words(rng, 2)
    body = b"".join(x.to_bytes() for x in w2)
    s.req(p.RUN)                                                 # nothing loaded
    s.req(0x55, b"\x00").req(0x00)                               # unknown commands
    s.req(p.LOAD_SEQ, struct.pack("<HHB", 9, 0, 2) + body)       # slot below 10 us
    s.req(p.LOAD_SEQ, struct.pack("<HHB", 50, 0, 0))             # no words
    s.req(p.LOAD_SEQ, struct.pack("<HHB", 50, 0, 3) + body)      # count and length disagree
    s.req(p.LOAD_SEQ, b"\x32\x00").req(p.STATUS).req(p.STOP)
    out.append(s)
    return out


# ---- running the harness -------------------------------------------------------------------------
@dataclass
class Trace:
    name: str
    events: list[dict]

    def all(self, kind: str) -> list[dict]:
        return [e for e in self.events if e["ev"] == kind]

    def rx_time(self, step: int) -> int:
        """When the device had the last packet of a step."""
        return [e for e in self.all("out") if e["step"] == step][-1]["t_ns"]

    def edges(self, pin: str) -> list[tuple[int, int | None]]:
        return [(e["t_ns"], e["level"]) for e in self.all("pin") if e["pin"] == pin]


def _run(batch: list[Script], tmp: Path) -> list[dict]:
    path = tmp / f"{batch[0].name}.json"
    path.write_text(json.dumps({"scenarios": [s.to_json() for s in batch]}))
    proc = subprocess.run([NODE, str(HARNESS), "--image", str(IMAGE), "--rp2040js", str(RP2040JS), "--script", str(path)],
                          capture_output=True, text=True, timeout=900)
    events = [json.loads(line) for line in proc.stdout.splitlines() if line.startswith("{")]
    errors = [e for e in events if e["ev"] in ("error", "timeout")]
    assert proc.returncode == 0 and not errors, f"harness failed: {errors} {proc.stderr[-3000:]}"
    return events


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    ss = scripts()
    tmp = tmp_path_factory.mktemp("emu")
    jobs = min(len(ss), max(1, (os.cpu_count() or 2) // 2), 6)
    batches = [ss[i::jobs] for i in range(jobs)]
    with ThreadPoolExecutor(jobs) as pool:
        events = [e for batch in pool.map(lambda b: _run(b, tmp), batches) for e in batch]
    return {s.name: (s, Trace(s.name, [e for e in events if e.get("scenario") == s.name])) for s in ss}


NAMES = [s.name for s in scripts()]


# ---- the model as the oracle -----------------------------------------------------------------------
@dataclass
class Expect:
    lo: float
    hi: float
    word: int
    what: str
    optional: bool = False
    nominal: float | None = None   # slot edge time, for sequence words


def model_running(m: ControllerModel, t: float) -> bool:
    c = copy.deepcopy(m)
    c.now = t
    return p.Status.from_payload(p.parse_reply(c.handle(p.request(p.STATUS)))[2]).running


class Segment:
    """A running sequence: the model's word for each slot, anchored on the firmware's marker latch."""

    def __init__(self, m: ControllerModel, t0_ns: float, name: str):
        self.m = copy.deepcopy(m)
        self.t0, self.slot, self.k, self.done, self.name = t0_ns, m.seq.slot_us * 1000, 0, False, name
        period = (len(m.seq.words) + 1) * self.slot
        self.end = t0_ns + m.seq.repeat * period if m.seq.repeat else math.inf

    def edges(self, until_ns: float, optional: bool = False) -> list[Expect]:
        out = []
        while not self.done and self.t0 + self.k * self.slot < until_ns:
            t = self.t0 + self.k * self.slot
            ws = self.m.words()
            i = int(self.m.word_index(np.array([self.m.t_start + (self.k + 0.5) * self.slot * 1e-9]))[0])
            out.append(Expect(t - SLOT_TOL_NS, t + SLOT_TOL_NS, ws[i].to_int(), f"{self.name}: slot {self.k}", optional, t))
            self.done = bool(self.m.seq.repeat) and i == len(ws) - 1    # the static word after the last period
            self.k += 1
        return out


@dataclass
class Checked:
    replies: list[bytes]
    want_replies: list[bytes]
    slot_errors_ns: list[float]
    clock_shifts_ns: list[float]
    ref_loss_gaps_ns: list[float]
    immediate_ns: list[float]
    marker_ns: list[float]


def check_against_model(script: Script, tr: Trace) -> Checked:
    m = ControllerModel()
    env, env_at = dict(script.env), {}
    for i, st in enumerate(script.steps):
        env.update(st.get("env", {}))
        env_at[i] = dict(env)
    latches = tr.all("latch")[1:]                       # the first is the boot latch, checked on its own
    replies = [bytes.fromhex(e["hex"]) for e in tr.all("reply")]
    reply_t = [e["t_first_ns"] for e in tr.all("reply")]   # the firmware had handled the frame by then
    rises = [t for t, v in tr.edges("TRIG") if v == 1]
    gpout = [e for e in tr.all("gpout") if e["n"] == 0]
    envs = [(e["t_ns"], e["env"]) for e in tr.all("env")]
    cur_env = dict(script.env)
    t_end = tr.all("end")[0]["t_ns"]
    exp: list[Expect] = []
    seg: Segment | None = None
    want, shifts, marker_ns, runs, gaps = [], [], [], 0, []

    for j, (step, frame) in enumerate(script.frames):
        t_rx = tr.rx_time(step)
        t_rep = reply_t[j] if j < len(reply_t) else math.inf
        e = env_at[step]
        while envs and envs[0][0] < t_rx:
            t_env, new_env = envs.pop(0)
            if seg is not None and not seg.done and ref_ok(cur_env, m.clock_src) and not ref_ok(new_env, m.clock_src):
                # The reference timing the slots has gone: no slot edges until the firmware notices
                # and moves the slot clock to the crystal; then the sequence carries on where it was.
                exp += seg.edges(t_env)
                gap = [x["t_ns"] for x in latches if x["t_ns"] > t_env][0] - (seg.t0 + seg.k * seg.slot)
                assert 0 <= gap <= REF_LOSS_GAP_NS, f"{script.name}: slots stopped for {gap} ns"
                seg.t0 += gap
                gaps.append(gap)
            cur_env = new_env
        m.usb_ma, m.ref_ok = advertised_ma(e), ref_ok(e, m.clock_src)
        req = p.parse_request(frame)
        cmd = req[0] if req else None
        if seg is not None:
            assert not (t_rx - 2_000 <= seg.end <= t_rep + 2_000), f"{script.name}: frame {j} lands on the end of a sequence"
            exp += seg.edges(t_rx)
        if cmd == p.RUN and m.seq is not None:
            trig = [t for t in rises if t_rx <= t <= t_rep]
            assert len(trig) == 1, f"{script.name}: frame {j} (RUN): TRIG rose {len(trig)} times"
            lo, hi = trig[0] + MARKER_AFTER_TRIG_NS[0], trig[0] + MARKER_AFTER_TRIG_NS[1]
            marker = [x for x in latches if lo <= x["t_ns"] <= hi]
            assert marker, f"{script.name}: frame {j} (RUN): no marker latch 11-12 us after TRIG rose"
            if seg is not None:
                exp += seg.edges(trig[0], optional=True)
            m.now = marker[0]["t_ns"] * 1e-9
            want.append(m.handle(frame))
            seg = Segment(m, marker[0]["t_ns"], f"{script.name} RUN {runs}")
            marker_ns.append(marker[0]["t_ns"] - trig[0])
            runs += 1
            continue
        m.now = t_rx * 1e-9
        was_running = model_running(m, m.now)
        want.append(m.handle(frame))
        ok = p.parse_reply(want[-1])[1] == p.OK
        if ok and (cmd in (p.SET_WORD, p.STOP) or (cmd == p.LOAD_SEQ and was_running)):
            if seg is not None:
                exp += seg.edges(t_rep, optional=True)
                seg = None
            exp.append(Expect(t_rx, t_rep + 20_000, m.words()[0].to_int(), f"{script.name}: frame {j} ({cmd:#04x}) latches"))
        elif ok and cmd == p.CLOCK_SRC and seg is not None and not seg.done:
            # The slot clock changes source when the firmware reprograms GPOUT0; later edges may move.
            t_sw = [g["t_ns"] for g in gpout if t_rx <= g["t_ns"] <= t_rep and g["hz"] > 0][-1]
            exp += seg.edges(t_sw)
            nominal = seg.t0 + seg.k * seg.slot
            after = [x["t_ns"] for x in latches if x["t_ns"] > t_sw]
            shift = after[0] - nominal
            assert abs(shift) <= CLOCK_SWITCH_SHIFT_NS + SLOT_TOL_NS, f"{script.name}: slot edges moved {shift} ns at CLOCK_SRC"
            seg.t0 += shift
            shifts.append(shift)
    if seg is not None:
        exp += seg.edges(t_end - 2_000)

    # Every expected latch happened, with the model's word; nothing else was latched. A stopping
    # frame's latch is the last one in its window: slot edges the sequence still played before the
    # firmware stopped it come first. So those are matched first, then the slot edges.
    left = list(latches)
    slot_err, immediate = [], []
    order = sorted(exp, key=lambda x: (x.nominal is not None, x.lo))
    for x in order:
        hit = [y for y in left if x.lo <= y["t_ns"] <= x.hi]
        if not hit:
            assert x.optional, f"missing latch: {x.what}, expected {x.word:014x} in [{x.lo:.0f}, {x.hi:.0f}] ns"
            continue
        y = hit[-1] if x.nominal is None else min(hit, key=lambda y: abs(y["t_ns"] - x.nominal))
        left.remove(y)
        assert int(y["word"], 16) == x.word, f"{x.what}: latched {y['word']}, model says {x.word:014x}"
        assert y["oe_n"] == 0, f"{x.what}: outputs disabled"
        assert y["clocks"] >= 56, f"{x.what}: only {y['clocks']} bits shifted"
        if x.nominal is not None:
            slot_err.append(y["t_ns"] - x.nominal)
            if x.what.split(": slot ")[-1] != "0":
                assert y["clocks"] == 56, f"{x.what}: {y['clocks']} SR_CLK edges"
        else:
            immediate.append(y["t_ns"] - x.lo)
    assert not left, f"{script.name}: latches the model does not explain: {left[:5]}"

    # TRIG: one pulse per RUN, falling as the marker latches.
    edges = tr.edges("TRIG")
    pulses = [(t0, t1) for (t0, v0), (t1, v1) in zip(edges, edges[1:]) if v0 == 1 and v1 == 0]
    assert len(rises) == runs == len(pulses), f"{script.name}: {len(rises)} TRIG pulses for {runs} RUNs"
    for (t0, t1), dm in zip(pulses, marker_ns):
        # busy_wait_us_32(12) on the 1 us timer: 11 to 12 us. The marker comes from the slot clock.
        assert 11_000 <= t1 - t0 <= 12_500 and abs((t1 - t0) - dm) <= 1_500, (t0, t1, dm)
    return Checked(replies, want, slot_err, shifts, gaps, immediate, marker_ns)


def check_boot(tr: Trace):
    """SR_OE_N held high (the board pulls it up), zeros shifted in and latched, then SR_OE_N low for good."""
    first = tr.all("latch")[0]
    assert (first["word"], first["clocks"], first["oe_n"]) == ("0" * 14, 56, 1)
    oe = tr.edges("SR_OE_N")
    driven = [(t, v) for t, v in oe if v is not None]
    assert [v for _, v in driven] == [1, 0], f"SR_OE_N went {oe}"
    assert driven[0][0] < first["t_ns"] < driven[1][0] and oe[-1] == driven[-1]
    assert [v for _, v in tr.edges("TRIG") if v is not None][0] == 0


# ---- tests -------------------------------------------------------------------------------------------
@pytest.mark.parametrize("name", NAMES)
def test_firmware_behaves_like_the_model(runs, name):
    script, tr = runs[name]
    check_boot(tr)
    r = check_against_model(script, tr)
    assert r.replies == r.want_replies, "\n".join(
        f"{i}: firmware {a.hex() if a else None} model {b.hex() if b else None}"
        for i, (a, b) in enumerate(zip(r.replies + [None] * 9, r.want_replies + [None] * 9)) if a != b)
    if r.slot_errors_ns:
        assert max(abs(e) for e in r.slot_errors_ns) <= SLOT_TOL_NS


def test_frames_are_protocol_frames(runs):
    """The harness sent exactly the bytes protocol.py built, split where the script said."""
    for script, tr in runs.values():
        for i, st in enumerate(script.steps):
            if "send" in st:
                got = "".join(e["hex"] for e in tr.all("out") if e["step"] == i)
                assert got == st["send"], (script.name, i)


def slot_clocks(tr: Trace) -> list[tuple[str, float, float]]:
    """GPOUT0 settings that lasted (clock_gpio_init writes CTRL, then DIV: 48 ns of a mixed setting)."""
    g = [e for e in tr.all("gpout") if e["n"] == 0]
    return [(a["src"], a["div"], a["hz"]) for a, b in zip(g, g[1:] + [None])
            if a["hz"] > 0 and (b is None or b["t_ns"] - a["t_ns"] > 1_000)]


def test_slot_clock_comes_from_the_selected_reference(runs):
    """GPOUT0 (GPIO21) runs at 1 MHz: 40 MHz / 40, HackRF 10 MHz / 10, or the 12 MHz crystal / 12."""
    ref40, ref10, xosc = ("gpin0", 40, 1e6), ("gpin1", 10, 1e6), ("xosc", 12, 1e6)
    # 40.2 MHz is accepted (slots 0.5 % short); at 40.6 MHz it is refused and the crystal takes over.
    assert slot_clocks(runs["clock_src"][1]) == [ref40, ref10, ref40, xosc, ("gpin0", 40, 1.005e6),
                                                 ("gpin0", 40, 1.015e6), xosc, ref40]
    assert slot_clocks(runs["boot_without_reference"][1]) == [xosc]
    assert slot_clocks(runs["reference_lost"][1]) == [ref40, ref10, xosc, ref10]
    assert slot_clocks(runs["no_reference_counter_waits"][1]) == [xosc, ref40, xosc]


def test_register_signal_timing(runs):
    """What the 74HC595s see from the PIO program at 125 MHz: 64 ns setup and pulse widths."""
    for _, tr in runs.values():
        t = tr.all("timing")[0]
        assert t["clk_sys_hz"] == 125e6
        assert t["data_clk_same_instant"] == 0 and t["clk_latch_same_instant"] == 0
        for k in ("min_data_setup_ns", "min_clk_high_ns", "min_clk_low_ns", "min_latch_high_ns", "min_clk_to_latch_ns"):
            assert t[k] >= 60, (tr.name, k, t[k])


def test_emulator_gaps_are_only_the_known_ones(runs):
    """The firmware touches nothing rp2040js lacks except flash/XIP set-up and the clock-divider restart."""
    known = ("IO_QSPI_BASE", "SSI:", "address: 1400", "clkDivRestart")
    for _, tr in runs.values():
        other = [w["what"] for w in tr.all("warnings")[0]["items"] if not any(k in w["what"] for k in known)]
        assert not other, (tr.name, other)
