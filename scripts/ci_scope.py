#!/usr/bin/env python3
"""Decide what a CI run checks (used by .github/workflows/ci.yml).

Full checks: lint, tests on Python 3.11 and 3.12, and the acceptance gate.
They run on pushes to main, on manual runs, and on pull requests into main
that are not drafts. Stacked pull requests are opened as drafts; drafts, and
pull requests into any other branch, get lint plus tests on Python 3.12
only. GitHub cannot merge a draft, and marking it "Ready for review" starts a
full run. So every pull request gets a full run before it can be merged.

Tests and the gate are skipped when every changed file is documentation
that no test, script or gate reads. Lint, which includes the
forbidden-wording check, always runs. When the changed files cannot be
worked out, the tests run.

Which files count as changed:

* a pull request that can be merged (not a draft): every file it changes
  against its base (``git diff base...head``), however many pushes it took;
* a draft's new push (``synchronize``): only that push's files. This is safe
  because the draft cannot be merged without the full run above;
* a push to main: that push's files.

Prints ``code=``, ``full=`` and ``pythons=`` lines for ``$GITHUB_OUTPUT``,
and a human-readable note on stderr.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import PurePosixPath

PYTHONS_FULL = ["3.11", "3.12"]
PYTHONS_REDUCED = ["3.12"]

# Documentation that a test, script or the CI gate reads. A change to one of
# these runs the tests. tests/test_ci_scope.py fails if a test reads another
# .md file that is missing from this list.
READ_BY_CHECKS = frozenset(
    {
        "docs/locus-signoff.md",  # tests/test_signoff.py: must match the generator
        "docs/figures/sim_stratification_summary.json",  # CI gate: must reproduce
    }
)

NO_COMMIT = "0" * 40


def is_doc_only(path: str) -> bool:
    """True for a documentation file that no check reads."""
    if path in READ_BY_CHECKS:
        return False
    p = PurePosixPath(path)
    if path == "LICENSE":
        return True
    if path.startswith("docs/figures/"):
        return False  # figures and their summaries stay with the code
    return p.suffix == ".md"


def needs_tests(paths: list[str] | None) -> bool:
    """Tests run unless the change is known and every file is doc-only."""
    if not paths:
        return True
    return not all(is_doc_only(p) for p in paths)


def is_full(event: str, base_ref: str, draft: bool) -> bool:
    """Full checks everywhere except drafts and PRs into a branch other than main."""
    if event != "pull_request":
        return True
    return base_ref == "main" and not draft


def _git_diff(*args: str) -> list[str] | None:
    try:
        out = subprocess.run(
            ["git", "diff", "--name-only", *args],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    return [line for line in out.splitlines() if line.strip()]


def changed_files(env: dict[str, str]) -> list[str] | None:
    """Changed files for this event, or None when they cannot be worked out."""
    event = env.get("EVENT", "")
    if event == "pull_request":
        head = env.get("HEAD_SHA", "")
        before = env.get("BEFORE", "")
        draft = env.get("DRAFT", "false") == "true"
        if draft and env.get("ACTION") == "synchronize" and before and before != NO_COMMIT:
            return _git_diff(before, head)
        base = env.get("BASE_SHA", "")
        if not base or not head:
            return None
        return _git_diff(f"{base}...{head}")
    if event == "push":
        before, after = env.get("BEFORE", ""), env.get("AFTER", "")
        if not before or before == NO_COMMIT or not after:
            return None
        return _git_diff(before, after)
    return None


def decide(env: dict[str, str]) -> dict[str, object]:
    event = env.get("EVENT", "")
    full = is_full(event, env.get("BASE_REF", ""), env.get("DRAFT", "false") == "true")
    paths = changed_files(env)
    code = needs_tests(paths)
    return {
        "code": code,
        "full": full,
        "pythons": PYTHONS_FULL if full else PYTHONS_REDUCED,
        "paths": paths,
    }


def main() -> int:
    d = decide(dict(os.environ))
    print(f"code={'true' if d['code'] else 'false'}")
    print(f"full={'true' if d['full'] else 'false'}")
    print(f"pythons={json.dumps(d['pythons'])}")
    paths = d["paths"]
    shown = "unknown (tests run)" if paths is None else f"{len(paths)} file(s)"
    what = "full checks" if d["full"] else "draft checks (lint + tests on 3.12)"
    tests = "tests run" if d["code"] else "documentation only: tests and gate skipped"
    print(f"ci_scope: {what}; changed: {shown}; {tests}", file=sys.stderr)
    for p in (paths or [])[:50]:
        print(f"  {p}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
