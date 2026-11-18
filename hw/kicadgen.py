"""A small KiCad 7 schematic writer.

Sheets are described as parts with a pin -> net mapping. Every pin gets a short
wire stub and a label, so connectivity is by net name and no wire routing is
needed. The writer also computes the netlist it intends, which check.py compares
with the netlist KiCad itself exports from the generated files.
"""
from __future__ import annotations

import os
import re
import uuid
from pathlib import Path

# KiCad's symbol libraries: $KICAD_SYMBOL_DIR, else KiCad's own variables, else the usual install places.
LIB_ENV = ("KICAD_SYMBOL_DIR", "KICAD9_SYMBOL_DIR", "KICAD8_SYMBOL_DIR", "KICAD7_SYMBOL_DIR")
LIB_DIRS = ("/usr/share/kicad/symbols", "/usr/local/share/kicad/symbols",
            "/Applications/KiCad/KiCad.app/Contents/SharedSupport/symbols")
PROJECT = "array16"
POWER_NETS = {"GND", "+5V", "+3V3", "+3V3_RF", "+1V1"}       # global labels; everything else is local to its sheet
NC = "NC"
GRID = 2.54


def uid(*key) -> str:
    return str(uuid.uuid5(uuid.NAMESPACE_URL, "array16/" + "/".join(map(str, key))))


# ---- s-expressions ------------------------------------------------------------------
def parse(text: str):
    tokens = re.findall(r'\(|\)|"(?:\\.|[^"\\])*"|[^\s()]+', text)
    stack, cur = [], []
    for t in tokens:
        if t == "(":
            stack.append(cur)
            cur = []
        elif t == ")":
            done, cur = cur, stack.pop()
            cur.append(done)
        else:
            cur.append(t)
    return cur[0]


def dump(x) -> str:
    return "(" + " ".join(dump(i) if isinstance(i, list) else i for i in x) + ")" if isinstance(x, list) else x


def q(s) -> str:
    return '"' + str(s).replace('"', "'") + '"'


def find(x, name):
    return [i for i in x if isinstance(i, list) and i and i[0] == name]


# ---- symbols --------------------------------------------------------------------------
class Symbol:
    """A library symbol: its s-expression and the position of every pin."""

    def __init__(self, lib_id: str, sx):
        self.lib_id, self.name = lib_id, lib_id.split(":")[1]
        sx[1] = q(lib_id)
        self.sx, self.pins = sx, []
        for unit in find(sx, "symbol"):
            for pn in find(unit, "pin"):
                at = find(pn, "at")[0]
                self.pins.append(dict(name=find(pn, "name")[0][1].strip('"'), number=find(pn, "number")[0][1].strip('"'),
                                      x=float(at[1]), y=float(at[2]), a=int(float(at[3])) if len(at) > 3 else 0))
        xs, ys = [p["x"] for p in self.pins], [p["y"] for p in self.pins]
        self.box = (min(xs), min(ys), max(xs), max(ys))


_libs: dict[str, list] = {}


def lib_dir() -> Path:
    for var in LIB_ENV:
        if os.environ.get(var):
            return Path(os.environ[var])
    for d in map(Path, LIB_DIRS):
        if d.is_dir():
            return d
    raise FileNotFoundError("KiCad symbol libraries not found; install KiCad 7+ or set KICAD_SYMBOL_DIR")


def library(lib_id: str) -> Symbol:
    """Load a symbol from the installed KiCad libraries, flattening `extends`."""
    lib, name = lib_id.split(":")
    if lib not in _libs:
        _libs[lib] = parse((lib_dir() / f"{lib}.kicad_sym").read_text())
    syms = {s[1].strip('"'): s for s in find(_libs[lib], "symbol")}
    sx = syms[name]
    ext = find(sx, "extends")
    if ext:
        parent = parse(dump(syms[ext[0][1].strip('"')]))
        pname = parent[1].strip('"')
        props = {p[1]: p for p in find(sx, "property")}
        body = [props.get(i[1], i) if isinstance(i, list) and i[0] == "property" else i for i in parent]
        for u in find(body, "symbol"):
            u[1] = q(u[1].strip('"').replace(pname, name, 1))
        sx = body
    return Symbol(lib_id, parse(dump(sx)))


def custom(name: str, ref: str, pins: list[tuple], width: float = 15.24, footprint: str = "") -> Symbol:
    """A box symbol. pins: (number, name, side) with side L, R, T or B; order sets position."""
    sides = {s: [p for p in pins if p[2] == s] for s in "LRTB"}
    rows = max(len(sides["L"]), len(sides["R"]), 1)
    h = max(rows + 1, 3) * GRID
    w = max(width, (max(len(sides["T"]), len(sides["B"])) + 1) * GRID)
    w, h = round(w / (2 * GRID)) * 2 * GRID, round(h / (2 * GRID)) * 2 * GRID
    font = "(effects (font (size 1.27 1.27)))"
    out = [f'(symbol {q(name)} (in_bom yes) (on_board yes)',
           f'(property "Reference" {q(ref)} (at {-w / 2} {h / 2 + 1.27} 0) (effects (font (size 1.27 1.27)) (justify left)))',
           f'(property "Value" {q(name)} (at {-w / 2} {-h / 2 - 1.27} 0) (effects (font (size 1.27 1.27)) (justify left)))',
           f'(property "Footprint" {q(footprint)} (at 0 0 0) (effects (font (size 1.27 1.27)) hide))',
           f'(property "Datasheet" "~" (at 0 0 0) (effects (font (size 1.27 1.27)) hide))',
           f'(symbol {q(name + "_0_1")} (rectangle (start {-w / 2} {h / 2}) (end {w / 2} {-h / 2}) (stroke (width 0.254) (type default)) (fill (type background))))',
           f'(symbol {q(name + "_1_1")}']
    for side, (dx, dy, ang) in {"L": (-1, 0, 0), "R": (1, 0, 180), "T": (0, 1, 270), "B": (0, -1, 90)}.items():
        for i, (num, pname, _) in enumerate(sides[side]):
            n = len(sides[side])
            if side in "LR":
                x, y = dx * (w / 2 + GRID), (n - 1) / 2 * GRID - i * GRID
                y = round(y / GRID) * GRID if n % 2 else round((y - GRID / 2) / GRID) * GRID + GRID
            else:
                x, y = (i - (n - 1) / 2) * GRID, dy * (h / 2 + GRID)
                x = round(x / GRID) * GRID if n % 2 else round((x - GRID / 2) / GRID) * GRID + GRID
            out.append(f'(pin passive line (at {x} {y} {ang}) (length {GRID}) (name {q(pname)} {font}) (number {q(num)} {font}))')
    out.append("))")
    return Symbol(f"{PROJECT}:{name}", parse(" ".join(out)))


# ---- sheets ---------------------------------------------------------------------------------
class Part:
    def __init__(self, sym: Symbol, prefix: str, value: str, nets: dict, footprint: str = "", note: str = ""):
        self.sym, self.prefix, self.value, self.footprint, self.note = sym, prefix, value, footprint, note
        self.pin_net = {}
        for p in sym.pins:
            key = p["name"] if p["name"] in nets else p["number"]
            if key not in nets:
                raise KeyError(f"{sym.lib_id}: pin {p['number']} '{p['name']}' has no net")
            self.pin_net[p["number"]] = nets[key]
        unused = set(nets) - {p["name"] for p in sym.pins} - {p["number"] for p in sym.pins}
        if unused:
            raise KeyError(f"{sym.lib_id}: no such pin(s) {sorted(unused)}")


class Sheet:
    def __init__(self, name: str, title: str, ports: list[str], paper: str = "A3"):
        self.name, self.title, self.ports, self.paper = name, title, ports, paper
        self.parts: list[Part] = []
        self.children: list[tuple[str, Sheet, dict]] = []      # (instance name, sheet, port -> local net)
        self.notes: list[str] = []

    def add(self, sym, prefix, value, nets, footprint="", note=""):
        self.parts.append(Part(sym, prefix, value, nets, footprint, note))

    def child(self, inst: str, sheet: "Sheet", nets: dict):
        missing = set(sheet.ports) - set(nets)
        if missing:
            raise KeyError(f"{inst}: unconnected sheet ports {sorted(missing)}")
        self.children.append((inst, sheet, nets))


class Design:
    """Walks the hierarchy, assigns references, writes the files and the intended netlist."""

    def __init__(self, root: Sheet, out_dir: Path):
        self.root, self.out = root, Path(out_dir)
        self.instances: dict[str, list[tuple[str, str]]] = {}    # sheet name -> [(uuid path, name path)]
        self.refs: dict[tuple[str, int, str], str] = {}          # (sheet name, part index, uuid path) -> reference
        self.counter: dict[str, int] = {}
        self.sheets: dict[str, Sheet] = {}
        self.parent = {}                                         # union-find over (name path, net)
        self.nodes = {}                                          # reference.pin -> key
        self._walk(root, "/" + uid("root"), "")

    def _find(self, k):
        self.parent.setdefault(k, k)
        while self.parent[k] != k:
            self.parent[k] = self.parent[self.parent[k]]
            k = self.parent[k]
        return k

    def _key(self, npath, net):
        return ("", net) if net in POWER_NETS else (npath, net)

    def _walk(self, sh: Sheet, upath: str, npath: str):
        self.sheets[sh.name] = sh
        self.instances.setdefault(sh.name, []).append((upath, npath))
        for i, part in enumerate(sh.parts):
            n = self.counter[part.prefix] = self.counter.get(part.prefix, 0) + 1
            ref = self.refs[(sh.name, i, upath)] = f"{part.prefix}{n}"
            for num, net in part.pin_net.items():
                if net != NC:
                    self.nodes[f"{ref}.{num}"] = self._find(self._key(npath, net))
        for inst, child, nets in sh.children:
            cpath = npath + "/" + inst
            for port, net in nets.items():
                a, b = self._find(self._key(npath, net)), self._find(self._key(cpath, port))
                self.parent[a] = b
            self._walk(child, upath + "/" + uid(sh.name, inst), cpath)

    def netlist(self) -> set[frozenset]:
        groups: dict = {}
        for node, key in self.nodes.items():
            groups.setdefault(self._find(key), set()).add(node)
        return {frozenset(g) for g in groups.values() if len(g) > 1}

    # ---- writing ---------------------------------------------------------------------------
    def write(self):
        self.out.mkdir(parents=True, exist_ok=True)
        pages = {}
        for name in self.sheets:
            for upath, _ in self.instances[name]:
                pages[upath] = len(pages) + 1
        for name, sh in self.sheets.items():
            fname = f"{PROJECT}.kicad_sch" if sh is self.root else f"{name}.kicad_sch"
            (self.out / fname).write_text(self._sheet_text(sh, pages))
        (self.out / f"{PROJECT}.kicad_pro").write_text('{\n  "meta": {"filename": "%s.kicad_pro", "version": 1}\n}\n' % PROJECT)

    def _label(self, kind, net, x, y, ang):
        just = {0: "left", 180: "right", 90: "left", 270: "right"}[ang]
        eff = f'(effects (font (size 1.27 1.27)) (justify {just}))'
        if kind == "global":
            return (f'(global_label {q(net)} (shape passive) (at {x:.2f} {y:.2f} {ang}) {eff} (uuid {uid("lbl", net, x, y)}) '
                    f'(property "Intersheetrefs" "${{INTERSHEET_REFS}}" (at 0 0 0) (effects (font (size 1.27 1.27)) hide)))')
        if kind == "hier":
            return f'(hierarchical_label {q(net)} (shape passive) (at {x:.2f} {y:.2f} {ang}) {eff} (uuid {uid("lbl", net, x, y)}))'
        return f'(label {q(net)} (at {x:.2f} {y:.2f} {ang}) (effects (font (size 1.27 1.27)) (justify {just} bottom)) (uuid {uid("lbl", net, x, y)}))'

    def _stub(self, sh, net, x, y, dx, dy, key):
        """Wire stub from a pin end, with the right kind of label at its far end."""
        if net == NC:
            return f'(no_connect (at {x:.2f} {y:.2f}) (uuid {uid("nc", key)}))'
        x2, y2 = x + dx * GRID, y + dy * GRID
        ang = {(1, 0): 0, (-1, 0): 180, (0, -1): 90, (0, 1): 270}[(dx, dy)]
        kind = "global" if net in POWER_NETS else "hier" if net in sh.ports else "local"
        wire = f'(wire (pts (xy {x:.2f} {y:.2f}) (xy {x2:.2f} {y2:.2f})) (stroke (width 0) (type default)) (uuid {uid("w", key)}))'
        return wire + "\n  " + self._label(kind, net, x2, y2, ang)

    def _sheet_text(self, sh: Sheet, pages) -> str:
        libs, body = {}, []
        w_page, h_page = {"A4": (297, 210), "A3": (420, 297), "A2": (594, 420), "A1": (841, 594)}[sh.paper]
        # child sheets as boxes in columns on the left, ports on the right edge of each box
        x, y = 40.0, 45.0
        for inst, child, nets in sh.children:
            h = (len(child.ports) + 1) * GRID
            if y + h > h_page - 25 and y > 45:
                x, y = x + 100, 45.0
            sx, sy, w = round(x / GRID) * GRID, round(y / GRID) * GRID, 25.4
            su = uid(sh.name, inst)
            pins = []
            for k, port in enumerate(child.ports):
                py = sy + (k + 1) * GRID
                pins.append(f'    (pin {q(port)} passive (at {sx + w:.2f} {py:.2f} 0) (effects (font (size 1.27 1.27)) (justify right)) (uuid {uid("spin", sh.name, inst, port)}))')
                body.append(self._stub(sh, nets[port], sx + w, py, 1, 0, (sh.name, inst, port)))
            paths = " ".join(f'(path {q(up)} (page {q(pages[up + "/" + su])}))' for up, _ in self.instances[sh.name])
            body.append(
                f'(sheet (at {sx:.2f} {sy:.2f}) (size {w} {h:.2f}) (fields_autoplaced) (stroke (width 0.1524) (type solid)) (fill (color 0 0 0 0.0)) (uuid {su})\n'
                f'    (property "Sheetname" {q(inst)} (at {sx:.2f} {sy - 0.7:.2f} 0) (effects (font (size 1.27 1.27)) (justify left bottom)))\n'
                f'    (property "Sheetfile" {q(child.name + ".kicad_sch")} (at {sx:.2f} {sy + h + 0.6:.2f} 0) (effects (font (size 1.27 1.27)) (justify left top)))\n'
                + "\n".join(pins) + f'\n    (instances (project {q(PROJECT)} {paths})))')
            y += h + 14
        x_left = x + 100 if sh.children else 40.0
        x, y, row_h = x_left, 45.0, 0.0
        for i, part in enumerate(sh.parts):
            s = part.sym
            libs[s.lib_id] = s
            x0, y0, x1, y1 = s.box
            side = 9 if all(p["a"] % 180 == 90 for p in s.pins) else 34   # room for stubs and labels at the sides
            left, right = side - x0, side + x1
            if x + left + right > w_page - 25 and x > x_left:
                x, y, row_h = x_left, y + row_h, 0.0
            px, py = round((x + left) / GRID) * GRID, round((y + y1 + 12) / GRID) * GRID
            row_h = max(row_h, (y1 - y0) + 30)
            x += left + right
            inst = self.instances[sh.name]
            ref0 = self.refs[(sh.name, i, inst[0][0])]
            paths = " ".join(f'(path {q(up)} (reference {q(self.refs[(sh.name, i, up)])}) (unit 1))' for up, _ in inst)
            hide = "(effects (font (size 1.27 1.27)) hide)"
            show = "(effects (font (size 1.27 1.27)) (justify left))"
            body.append(
                f'(symbol (lib_id {q(s.lib_id)}) (at {px:.2f} {py:.2f} 0) (unit 1) (in_bom yes) (on_board yes) (dnp no) (uuid {uid("sym", sh.name, i)})\n'
                f'    (property "Reference" {q(ref0)} (at {px + x0:.2f} {py - y1 - 7.62:.2f} 0) {show})\n'
                f'    (property "Value" {q(part.value)} (at {px + x0:.2f} {py - y1 - 5.08:.2f} 0) {show})\n'
                f'    (property "Footprint" {q(part.footprint)} (at {px:.2f} {py:.2f} 0) {hide})\n'
                f'    (property "Datasheet" "~" (at {px:.2f} {py:.2f} 0) {hide})\n'
                + "".join(f'    (pin {q(p["number"])} (uuid {uid("pin", sh.name, i, p["number"])}))\n' for p in s.pins)
                + f'    (instances (project {q(PROJECT)} {paths})))')
            for p in s.pins:
                d = {0: (-1, 0), 180: (1, 0), 90: (0, 1), 270: (0, -1)}[p["a"] % 360]
                body.append(self._stub(sh, part.pin_net[p["number"]], px + p["x"], py - p["y"], d[0], d[1], (sh.name, i, p["number"])))
        if sh.parts and y + row_h > h_page - 35:
            raise ValueError(f"sheet {sh.name}: parts run off the {sh.paper} page; use a larger paper size")
        for k, note in enumerate(sh.notes):
            body.append(f'(text {q(note)} (at 30 {20 + 4 * k:.2f} 0) (effects (font (size 1.6 1.6)) (justify left bottom)) (uuid {uid("note", sh.name, k)}))')
        root = sh is self.root
        return ("(kicad_sch (version 20230121) (generator array16_kicadgen)\n"
                f"  (uuid {uid('root') if root else uid('file', sh.name)})\n  (paper {q(sh.paper)})\n"
                f'  (title_block (title {q(sh.title)}) (date "2026-11-14") (rev "0.1") (company "array16: 4x4 phased-array SDR board"))\n'
                "  (lib_symbols\n    " + "\n    ".join(dump(s.sx) for s in libs.values()) + "\n  )\n  "
                + "\n  ".join(body) + "\n"
                + ('  (sheet_instances (path "/" (page "1")))\n' if root else "") + ")\n")
