"""The laws ledger and what it cites, checked against each other in both directions (docs/LAWS.md)."""

from __future__ import annotations

import ast
import re
from dataclasses import dataclass
from functools import cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "docs" / "LAWS.md"
FORMAL = ROOT / "formal" / "Array16"
LABELS = {"Proposition", "Claim", "Hypothesis"}
LAW_ID = r"[A-Z]{3}-H?\d+"
SPAN = re.compile(r"`([^`]+)`")
DECLARATION = re.compile(
    r"/--(?P<doc>.*?)-/\s*(?:noncomputable\s+)?(?:theorem|lemma|def)\s+(?P<name>[\w.']+)", re.DOTALL
)
LAW_TEST = re.compile(r"^test_([a-z]{3})_(h?\d+)_")


@dataclass(frozen=True)
class Row:
    law: str
    statement: str
    structure: str
    anchor: str
    label: str
    falsifier: str
    formal: str

    def cites(self, cell: str) -> list[str]:
        """The `file::symbol` references in one of this row's cells."""
        return [span for span in SPAN.findall(cell) if "::" in span]


@cache
def rows() -> tuple[Row, ...]:
    """The ledger's rows: the table lines whose first cell is a law id. A cell may hold an escaped bar."""
    found = []
    for line in LEDGER.read_text().splitlines():
        cells = [cell.strip() for cell in re.split(r"(?<!\\)\|", line.strip())[1:-1]]
        if len(cells) == 7 and re.fullmatch(LAW_ID, cells[0]):
            found.append(Row(*cells))
    return tuple(found)


def defines(path: Path, dotted: str) -> bool:
    """Whether the Python module at `path` defines the dotted name: a function, a class, or a member of a class."""
    scope: list[ast.stmt] = ast.parse(path.read_text()).body
    for part in dotted.split("."):
        inner = next((n.body for n in scope if isinstance(n, (ast.FunctionDef, ast.ClassDef)) and n.name == part), None)
        if inner is None:
            return False
        scope = inner
    return True


@cache
def declarations() -> dict[str, str | None]:
    """`file::name` -> the law id in bold at the head of its docstring, for each documented Lean declaration."""
    marked: dict[str, str | None] = {}
    for source in sorted(FORMAL.rglob("*.lean")):
        for found in DECLARATION.finditer(source.read_text()):
            head = re.match(rf"\*\*({LAW_ID})\*\*", found["doc"].strip())
            marked[f"{source.relative_to(FORMAL).as_posix()}::{found['name']}"] = head[1] if head else None
    return marked


def test_the_ledger_has_rows_with_distinct_ids_and_known_labels() -> None:
    laws = [row.law for row in rows()]
    assert len(laws) >= 20
    assert len(set(laws)) == len(laws)
    assert {row.label for row in rows()} <= LABELS
    assert all(("-H" in row.law) == (row.label == "Hypothesis") for row in rows())


def test_every_anchor_exists() -> None:
    for row in rows():
        assert row.cites(row.anchor), f"{row.law} names no anchor"
        for cited in row.cites(row.anchor):
            path, symbol = cited.split("::")
            assert defines(ROOT / path, symbol), f"{row.law}: no {symbol} in {path}"


def test_every_falsifier_exists() -> None:
    for row in rows():
        assert row.cites(row.falsifier), f"{row.law} names no falsifier"
        for cited in row.cites(row.falsifier):
            path, test = cited.split("::")
            assert test.startswith("test_"), f"{row.law}: {test} is not a test"
            assert defines(ROOT / path, test), f"{row.law}: no {test} in {path}"


def test_every_cited_declaration_exists_and_carries_its_law() -> None:
    for row in rows():
        for cited in row.cites(row.formal):
            assert cited in declarations(), f"{row.law}: no declaration {cited}"
            assert declarations()[cited] == row.law, f"{cited} is marked {declarations()[cited]}, cited by {row.law}"


def test_every_marked_declaration_is_cited_by_its_law() -> None:
    cited = {(row.law, name) for row in rows() for name in row.cites(row.formal)}
    for name, law in declarations().items():
        if law is not None:
            assert (law, name) in cited, f"{name} is marked {law}, and no row of that law cites it"


def test_a_proposition_is_proved_and_nothing_else_claims_to_be() -> None:
    for row in rows():
        if row.label == "Proposition":
            assert row.cites(row.formal), f"{row.law} is a Proposition with no proof"
        else:
            assert row.formal == "—", f"{row.law} is a {row.label} and cites a proof"


def test_every_test_named_for_a_law_is_one_of_its_falsifiers() -> None:
    falsifiers = {row.law: row.cites(row.falsifier) for row in rows()}
    for source in sorted((ROOT / "tests").glob("test_*.py")):
        for node in ast.parse(source.read_text()).body:
            named = LAW_TEST.match(node.name) if isinstance(node, ast.FunctionDef) else None
            if named and isinstance(node, ast.FunctionDef):
                law = f"{named[1].upper()}-{named[2].upper()}"
                assert f"tests/{source.name}::{node.name}" in falsifiers.get(law, []), f"{node.name}: not in {law}"


def test_the_checks_have_teeth() -> None:
    # A reference to something that is not there is seen to be missing, by the same functions the checks use.
    optics = ROOT / "array16" / "particles" / "optics.py"
    assert defines(optics, "GuideOptics.closed_orbit")
    assert not defines(optics, "GuideOptics.no_such_method")
    assert not defines(optics, "closed_orbit")  # a method is not a module function
    assert declarations()["Guide.lean::lap_closedOrbit"] == "GDE-1"
    assert declarations()["Guide.lean::lap_sub_closedOrbit"] is None  # documented, and marked for no law
    assert "Guide.lean::no_such_theorem" not in declarations()
