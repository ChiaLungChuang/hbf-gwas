#!/usr/bin/env python3
"""Fail if forbidden wording appears in tracked files or commit messages.

docs/BRIEF.md bans three phrasings because 1000 Genomes has no HbF phenotype
and every phenotype in this project is simulated. docs/BRIEF.md itself quotes
the banned phrases, so it is the only file exempt from the scan (plus this
script and its test, which must spell out the patterns).

Usage:
    python scripts/check_language.py            # scan tracked files
    python scripts/check_language.py --commits  # also scan every commit message
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

PATTERNS = [
    re.compile(r"\bre-?discover\w*", re.IGNORECASE),
    re.compile(r"\breplicat\w*", re.IGNORECASE),
    re.compile(r"\bvalidat\w*\s+against\s+biolog\w*", re.IGNORECASE),
]
EXEMPT = {
    "docs/BRIEF.md",
    "scripts/check_language.py",
}


def find_violations(text: str) -> list[str]:
    """Return every banned phrase found in ``text``."""
    hits: list[str] = []
    for pattern in PATTERNS:
        hits.extend(m.group(0) for m in pattern.finditer(text))
    return hits


def _git(*args: str) -> str:
    return subprocess.run(["git", *args], check=True, capture_output=True, text=True).stdout


def scan_files() -> list[str]:
    problems: list[str] = []
    for rel in _git("ls-files").splitlines():
        if rel in EXEMPT or not Path(rel).is_file():
            continue
        try:
            lines = Path(rel).read_text(encoding="utf-8").splitlines()
        except UnicodeDecodeError:
            continue  # binary file (e.g. a committed PNG)
        for lineno, line in enumerate(lines, start=1):
            for hit in find_violations(line):
                problems.append(f"{rel}:{lineno}: {hit!r}")
    return problems


def scan_commits() -> list[str]:
    problems: list[str] = []
    log = _git("log", "--all", "--format=%H%x00%B%x01")
    for record in log.split("\x01"):
        if "\x00" not in record:
            continue
        sha, message = record.strip().split("\x00", 1)
        for hit in find_violations(message):
            problems.append(f"commit {sha[:10]}: {hit!r}")
    return problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--commits", action="store_true", help="also scan commit messages")
    args = parser.parse_args(argv)

    problems = scan_files()
    if args.commits:
        problems += scan_commits()
    if problems:
        print("Forbidden wording found (see docs/BRIEF.md):")
        print("\n".join(f"  {p}" for p in problems))
        return 1
    print("Wording check passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
