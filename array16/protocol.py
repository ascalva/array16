"""Host <-> beam controller protocol (USB CDC serial, little endian).

Request : 0xA5 | cmd | len_lo len_hi | payload | crc16
Reply   : 0x5A | cmd | status | len_lo len_hi | payload | crc16
crc16 is CRC-16/CCITT-FALSE over everything before it.

The same frames drive the firmware model in fwsim.py and the real RP2040.
On the serial stream the firmware hunts for 0xA5, drops a frame whose length
field exceeds MAX_PAYLOAD (and resumes hunting one byte later), drops a partial
frame after FRAME_TIMEOUT_MS without a byte, and answers a CRC failure with
command 0 and ERR_CRC, as fwsim does.

FW_VERSION 2 (2026-11-21) widened the word count of LOAD_SEQ and of STATUS to u16, for
up to MAX_SEQ = 256 words; STATUS grew from 9 bytes to 10. PING is the same in both versions
and tells the host which one it talks to. A version 2 controller answers any version 1
LOAD_SEQ with ERR_ARG: a payload of 5 + 7m bytes never has the 6 + 7n bytes its count asks
for. A version 1 controller answers a version 2 LOAD_SEQ of 1 to 63 words with ERR_ARG the
same way, but one of 64 words or more is over version 1's 453-byte limit: it drops that as
noise, without a reply, and may then take a 0xA5 among the words for the start of a frame,
which ends in a stray ERR_CRC reply or swallows the next request. So a host reads PING
before it sends LOAD_SEQ (api.Array.check_version).
"""
from __future__ import annotations

import struct
from dataclasses import dataclass

from .board import WORD_BYTES, ControlWord

SOF_REQ, SOF_REP = 0xA5, 0x5A
PING, SET_WORD, LOAD_SEQ, RUN, STOP, STATUS, CLOCK_SRC = 0x01, 0x10, 0x11, 0x12, 0x13, 0x20, 0x30
OK, ERR_CRC, ERR_CMD, ERR_ARG = 0, 1, 2, 3
MAX_SEQ = 256                                # words after the marker
SEQ_HEAD = 6                                 # LOAD_SEQ: slot_us u16, repeat u16, count u16, then the words
MAX_PAYLOAD = SEQ_HEAD + MAX_SEQ * WORD_BYTES  # the largest request is a full LOAD_SEQ, 1798 bytes
FRAME_TIMEOUT_MS = 50
FW_VERSION = 2


def crc16(data: bytes) -> int:
    crc = 0xFFFF
    for b in data:
        crc ^= b << 8
        for _ in range(8):
            crc = ((crc << 1) ^ 0x1021 if crc & 0x8000 else crc << 1) & 0xFFFF
    return crc


def request(cmd: int, payload: bytes = b"") -> bytes:
    body = bytes([SOF_REQ, cmd]) + struct.pack("<H", len(payload)) + payload
    return body + struct.pack("<H", crc16(body))


def reply(cmd: int, status: int = OK, payload: bytes = b"") -> bytes:
    body = bytes([SOF_REP, cmd, status]) + struct.pack("<H", len(payload)) + payload
    return body + struct.pack("<H", crc16(body))


def parse_request(frame: bytes) -> tuple[int, bytes] | None:
    if len(frame) < 6 or frame[0] != SOF_REQ:
        return None
    (n,) = struct.unpack_from("<H", frame, 2)
    if len(frame) != 6 + n or crc16(frame[:-2]) != struct.unpack_from("<H", frame, 4 + n)[0]:
        return None
    return frame[1], frame[4:4 + n]


def parse_reply(frame: bytes) -> tuple[int, int, bytes]:
    if len(frame) < 7 or frame[0] != SOF_REP:
        raise ValueError("not a reply frame")
    (n,) = struct.unpack_from("<H", frame, 3)
    if len(frame) != 7 + n or crc16(frame[:-2]) != struct.unpack_from("<H", frame, 5 + n)[0]:
        raise ValueError("bad reply length or CRC")
    return frame[1], frame[2], frame[5:5 + n]


@dataclass
class Status:
    """STATUS reply payload, 10 bytes; count is u16 from version 2 on.

    clock_src: reference the host selected (0 = on-board 40 MHz, 1 = HackRF 10 MHz).
    ref_ok: that reference is present and times the slots. If 0 the controller has fallen
    back to its own 12 MHz crystal: sequences still play, but slot edges drift against the radios.
    usb_ma: current the USB-C source advertises on CC (500, 1500 or 3000; 0 = no CC seen).
    Both radio boards need 1500.
    """

    version: int
    running: bool
    count: int
    slot_us: int
    clock_src: int = 0
    ref_ok: bool = True
    usb_ma: int = 0

    FORMAT = "<BBHHBBH"

    def payload(self) -> bytes:
        return struct.pack(self.FORMAT, self.version, int(self.running), self.count, self.slot_us,
                           self.clock_src, int(self.ref_ok), self.usb_ma)

    @classmethod
    def from_payload(cls, b: bytes) -> "Status":
        size = struct.calcsize(cls.FORMAT)
        if len(b) != size:
            sender = "protocol version %d" % b[0] if b else "a controller"
            raise ValueError("a %d-byte STATUS from %s: this host reads version %d's %d bytes"
                             % (len(b), sender, FW_VERSION, size))
        v, run, n, slot, src, ok, ma = struct.unpack(cls.FORMAT, b)
        return cls(v, bool(run), n, slot, src, bool(ok), ma)


@dataclass
class Sequence:
    """A list of control words stepped at a fixed slot time.

    The controller plays one blank marker slot, then the words, and repeats.
    repeat = 0 means until STOP.
    """

    words: list[ControlWord]
    slot_us: int = 50
    repeat: int = 0

    def payload(self) -> bytes:
        if not 1 <= len(self.words) <= MAX_SEQ:
            raise ValueError("sequence length must be 1..%d" % MAX_SEQ)
        head = struct.pack("<HHH", self.slot_us, self.repeat, len(self.words))
        return head + b"".join(w.to_bytes() for w in self.words)

    @classmethod
    def from_payload(cls, p: bytes) -> "Sequence":
        slot_us, repeat, n = struct.unpack_from("<HHH", p, 0)
        if len(p) != SEQ_HEAD + n * WORD_BYTES or not 1 <= n <= MAX_SEQ or slot_us < 10:
            raise ValueError("bad sequence payload")
        body = p[SEQ_HEAD:]
        words = [ControlWord.from_bytes(body[i * WORD_BYTES:(i + 1) * WORD_BYTES]) for i in range(n)]
        return cls(words, slot_us, repeat)
