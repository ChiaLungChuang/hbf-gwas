"""Committed records of earlier criteria versions stay exactly as recorded.

CI gates only the current criteria version. Older runs are records: their
results and gate outputs must never be edited or re-judged. These hashes were
taken from the files as committed on 2026-11-10.
"""

import hashlib
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

RECORDS = {
    "runs/2026-11-10-sim-stratification/results.json": (
        "53308e1e606a8f76aa4259711e2119e1f907d3f1d542a7a05e50569183d0ee39"
    ),
    "runs/2026-11-10-sim-stratification/gate.json": (
        "959c29f0ef316b03f360da1ae25f79ef92c5f4ba345ffe05e5167de5ecfe0ea7"
    ),
    "runs/2026-11-10-sim-v2/results.json": (
        "a3e68a0750498d1f5535a59efc3e050b67775d8752127a90fe86d655131b6da2"
    ),
    "runs/2026-11-10-sim-v2/gate.json": (
        "1e85d9e0e1a22bf84793e51f7714c900f4ac30eb3e6e49d71694f27a3d4a5685"
    ),
}


@pytest.mark.parametrize("path,digest", sorted(RECORDS.items()))
def test_record_unchanged(path, digest):
    assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == digest


def test_sim_v1_passed_and_sim_v2_failed_as_recorded():
    v1 = json.loads((ROOT / "runs/2026-11-10-sim-stratification/gate.json").read_text())
    v2 = json.loads((ROOT / "runs/2026-11-10-sim-v2/gate.json").read_text())
    assert v1["passed"] is True
    assert v2["passed"] is False
    failed = [c["name"] for c in v2["checks"] if c["status"] == "FAIL"]
    assert failed == ["locus HBS1L-MYB recovery"]
