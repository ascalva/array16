"""Round-trip check: KiCad's own netlist of the generated schematics must equal the intended one.

    python3 hw/check.py          (needs kicad-cli 7 or newer)
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import kicadgen as K            # noqa: E402
from design import build        # noqa: E402


def kicad_netlist(sch: Path) -> set:
    net = sch.with_suffix(".net")
    subprocess.run(["kicad-cli", "sch", "export", "netlist", "--format", "kicadsexpr", "-o", str(net), str(sch)], check=True, capture_output=True)
    tree = K.parse(net.read_text())
    groups = set()
    for n in K.find(K.find(tree, "nets")[0], "net"):
        nodes = frozenset(K.find(x, "ref")[0][1].strip('"') + "." + K.find(x, "pin")[0][1].strip('"') for x in K.find(n, "node"))
        if len(nodes) > 1:
            groups.add(nodes)
    net.unlink()
    return groups


def main():
    design, summary = build()
    want, got = design.netlist(), kicad_netlist(HERE / "kicad" / "array16.kicad_sch")
    missing, extra = want - got, got - want
    print(f"parts {summary['total']}, intended nets {len(want)}, KiCad nets {len(got)}, differing {len(missing) + len(extra)}")
    for label, nets in (("intended but not in KiCad", missing), ("in KiCad but not intended", extra)):
        for n in sorted(nets, key=lambda s: sorted(s))[:5]:
            print(" ", label + ":", sorted(n)[:8])
    return 0 if not missing and not extra else 1


if __name__ == "__main__":
    sys.exit(main())
