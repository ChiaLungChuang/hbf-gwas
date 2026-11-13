"""Pinned dependency versions (requirements/constraints.txt, D-058).

Every dependency, direct and transitive, has one exact version per supported
Python version. CI and scripts/rebuild.sh install with these constraints and
set HBFGWAS_REQUIRE_LOCK=1, which makes the last test compare what is
installed with the pins.
"""

from __future__ import annotations

import os
import re
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
LOCK = ROOT / "requirements" / "constraints.txt"
PYTHONS = ("3.11", "3.12")
LINE = re.compile(
    r'^([a-z0-9][a-z0-9._-]*)==([0-9][^;\s]*)(?:; python_version (<|>=) "([0-9.]+)")?$'
)


def _pins() -> list[tuple[str, str, str | None, str | None]]:
    out = []
    for raw in LOCK.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        m = LINE.match(line)
        assert m, f"not an exact pin: {line!r}"
        out.append(m.groups())
    return out


def _applies(op, ver, py: str) -> bool:
    if op is None:
        return True
    a = tuple(int(x) for x in py.split("."))
    b = tuple(int(x) for x in ver.split("."))
    return a < b if op == "<" else a >= b


def _resolved(py: str) -> dict[str, str]:
    got: dict[str, str] = {}
    for name, ver, op, mv in _pins():
        if _applies(op, mv, py):
            assert name not in got, f"{name} pinned twice for Python {py}"
            got[name] = ver
    return got


def _direct() -> set[str]:
    import tomllib

    proj = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    reqs = list(proj["dependencies"]) + list(proj["optional-dependencies"]["dev"])
    return {re.split(r"[<>=!~ ;\[]", r, maxsplit=1)[0].lower() for r in reqs}


def test_every_line_is_an_exact_pin():
    assert len(_pins()) >= 20


@pytest.mark.parametrize("py", PYTHONS)
def test_one_version_per_package_for_each_supported_python(py):
    a, b = _resolved("3.11"), _resolved("3.12")
    assert set(a) == set(b)  # the same packages for both
    assert _resolved(py)


@pytest.mark.parametrize("py", PYTHONS)
def test_every_direct_dependency_is_pinned(py):
    missing = _direct() - set(_resolved(py))
    assert not missing, missing


def test_pyproject_pins_ruff_like_the_lock():
    assert _resolved("3.12")["ruff"] == "0.17.4"
    assert "ruff==0.17.4" in (ROOT / "pyproject.toml").read_text()


@pytest.mark.skipif(
    os.environ.get("HBFGWAS_REQUIRE_LOCK") != "1",
    reason="set HBFGWAS_REQUIRE_LOCK=1 (CI, rebuild.sh) to compare installed versions",
)
def test_installed_versions_equal_the_pins():
    py = f"{sys.version_info.major}.{sys.version_info.minor}"
    wrong = {}
    for name, want in _resolved(py).items():
        try:
            have = version(name)
        except PackageNotFoundError:
            have = None
        if have != want:
            wrong[name] = (have, want)
    assert not wrong, wrong
