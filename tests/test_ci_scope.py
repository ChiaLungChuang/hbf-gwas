"""CI scope rules (scripts/ci_scope.py) and the workflow files that use them.

The rules cut CI minutes without weakening what is checked before a merge to
main: a pull request that can be merged always gets the full checks.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import ci_scope  # noqa: E402

WORKFLOWS = ROOT / ".github" / "workflows"
PLINK2_ACTION = ROOT / ".github" / "actions" / "install-plink2" / "action.yml"


def _workflow(name: str) -> dict:
    wf = yaml.safe_load((WORKFLOWS / name).read_text())
    wf["on"] = wf.get("on", wf.pop(True, None))  # PyYAML reads the key `on` as True
    return wf


# --- which runs get the full checks -------------------------------------------


@pytest.mark.parametrize(
    ("event", "base", "draft", "full"),
    [
        ("push", "", False, True),  # push to main
        ("workflow_dispatch", "", False, True),
        ("pull_request", "main", False, True),  # can be merged: full checks
        ("pull_request", "main", True, False),  # draft (stacked): lint + 3.12
        ("pull_request", "signoff-3", False, False),  # into another branch
    ],
)
def test_full_checks_only_where_a_merge_to_main_can_follow(event, base, draft, full):
    assert ci_scope.is_full(event, base, draft) is full


def test_reduced_checks_still_run_tests_on_one_python():
    assert ci_scope.PYTHONS_REDUCED == ["3.12"]
    assert ci_scope.PYTHONS_FULL == ["3.11", "3.12"]
    assert set(ci_scope.PYTHONS_REDUCED) <= set(ci_scope.PYTHONS_FULL)


# --- documentation-only changes ------------------------------------------------


@pytest.mark.parametrize(
    ("paths", "tests"),
    [
        (["README.md"], False),
        (["CHANGELOG.md", "docs/decisions.md", "docs/power.md"], False),
        (["runs/2026-11-11/README.md"], False),
        (["LICENSE"], False),
        (["README.md", "src/hbfgwas/assoc.py"], True),
        (["docs/locus-signoff.md"], True),  # test_signoff compares it with the generator
        (["docs/figures/sim_stratification_summary.json"], True),  # the gate reads it
        (["docs/figures/sim_stratification_qq.png"], True),
        (["runs/2026-11-11/prereg.json"], True),
        (["config/analysis.yaml"], True),
        ([".github/workflows/ci.yml"], True),
        ([], True),  # nothing known: run the tests
        (None, True),  # changes could not be worked out: run the tests
    ],
)
def test_tests_are_skipped_only_for_docs_no_check_reads(paths, tests):
    assert ci_scope.needs_tests(paths) is tests


def test_every_markdown_file_a_test_reads_is_listed_as_read_by_checks():
    # A quoted *.md name in a test means the test may read that file; a change
    # to it must then run the tests.
    listed = {Path(p).name for p in ci_scope.READ_BY_CHECKS}
    found = set()
    for f in (ROOT / "tests").glob("*.py"):
        if f.name == Path(__file__).name:
            continue  # this file names doc files only as examples
        found |= set(re.findall(r"""["']([\w./-]+\.md)["']""", f.read_text()))
    assert {Path(p).name for p in found} <= listed, found


def test_every_docs_file_the_workflows_read_is_listed_as_read_by_checks():
    for f in [*WORKFLOWS.glob("*.yml"), PLINK2_ACTION]:
        for p in re.findall(r"docs/[\w./-]+", f.read_text()):
            if p.endswith(".md") and p != "docs/ci.md":  # docs/ci.md: a comment
                raise AssertionError(f"{f.name} reads {p}")
            if p.startswith("docs/figures/"):
                assert p in ci_scope.READ_BY_CHECKS, (f.name, p)


# --- changed files, from a real git history ------------------------------------


def _git(cwd: Path, *args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=cwd, check=True, capture_output=True, text=True
    ).stdout.strip()


@pytest.fixture
def repo(tmp_path, monkeypatch):
    _git(tmp_path, "init", "-q", "-b", "main")
    _git(tmp_path, "config", "user.email", "ci@example.invalid")
    _git(tmp_path, "config", "user.name", "ci")
    shas = {}

    def commit(name: str, files: dict[str, str]) -> str:
        for rel, text in files.items():
            (tmp_path / rel).parent.mkdir(parents=True, exist_ok=True)
            (tmp_path / rel).write_text(text)
        _git(tmp_path, "add", "-A")
        _git(tmp_path, "commit", "-q", "-m", name)
        shas[name] = _git(tmp_path, "rev-parse", "HEAD")
        return shas[name]

    commit("base", {"src/a.py": "x = 1\n", "README.md": "a\n"})
    _git(tmp_path, "checkout", "-q", "-b", "feature")
    commit("code", {"src/a.py": "x = 2\n"})
    commit("docs", {"README.md": "b\n"})
    monkeypatch.chdir(tmp_path)
    return shas


def test_mergeable_pr_counts_every_file_against_its_base(repo):
    env = {
        "EVENT": "pull_request",
        "ACTION": "synchronize",
        "BASE_REF": "main",
        "DRAFT": "false",
        "BASE_SHA": repo["base"],
        "HEAD_SHA": repo["docs"],
        "BEFORE": repo["code"],
    }
    d = ci_scope.decide(env)
    assert sorted(d["paths"]) == ["README.md", "src/a.py"]
    assert d["code"] is True and d["full"] is True


def test_draft_push_counts_only_that_push(repo):
    env = {
        "EVENT": "pull_request",
        "ACTION": "synchronize",
        "BASE_REF": "main",
        "DRAFT": "true",
        "BASE_SHA": repo["base"],
        "HEAD_SHA": repo["docs"],
        "BEFORE": repo["code"],
    }
    d = ci_scope.decide(env)
    assert d["paths"] == ["README.md"]
    assert d["code"] is False and d["full"] is False


def test_draft_marked_ready_gets_the_full_checks_on_every_file(repo):
    env = {
        "EVENT": "pull_request",
        "ACTION": "ready_for_review",
        "BASE_REF": "main",
        "DRAFT": "false",
        "BASE_SHA": repo["base"],
        "HEAD_SHA": repo["docs"],
        "BEFORE": "",
    }
    d = ci_scope.decide(env)
    assert d["full"] is True and d["code"] is True
    assert d["pythons"] == ["3.11", "3.12"]


def test_push_to_main_counts_that_push(repo):
    d = ci_scope.decide({"EVENT": "push", "BEFORE": repo["code"], "AFTER": repo["docs"]})
    assert d["paths"] == ["README.md"] and d["code"] is False and d["full"] is True


def test_unknown_history_runs_the_tests(repo):
    unknown = "f" * 40
    d = ci_scope.decide({"EVENT": "push", "BEFORE": unknown, "AFTER": repo["docs"]})
    assert d["paths"] is None and d["code"] is True
    d = ci_scope.decide({"EVENT": "push", "BEFORE": "0" * 40, "AFTER": repo["docs"]})
    assert d["code"] is True
    assert ci_scope.decide({"EVENT": "workflow_dispatch"})["code"] is True


def test_output_lines_for_github(repo, capsys, monkeypatch):
    monkeypatch.setattr(
        ci_scope.os,
        "environ",
        {"EVENT": "push", "BEFORE": repo["code"], "AFTER": repo["docs"]},
    )
    assert ci_scope.main() == 0
    out = capsys.readouterr().out.splitlines()
    assert out == ["code=false", "full=true", 'pythons=["3.11", "3.12"]']
    assert json.loads(out[2].split("=", 1)[1]) == ["3.11", "3.12"]


# --- the workflow files --------------------------------------------------------


def test_ci_cancels_superseded_runs_and_starts_a_full_run_on_ready_for_review():
    wf = _workflow("ci.yml")
    assert wf["concurrency"]["cancel-in-progress"] is True
    assert "pull_request.number" in wf["concurrency"]["group"]
    assert "ready_for_review" in wf["on"]["pull_request"]["types"]
    assert wf["on"]["push"]["branches"] == ["main"]


def test_ci_jobs_follow_the_scope_and_have_timeouts():
    jobs = _workflow("ci.yml")["jobs"]
    assert set(jobs) == {"lint", "test", "gate"}
    assert all("timeout-minutes" in j for j in jobs.values())
    assert "if" not in jobs["lint"]  # lint (with the wording check) always runs
    assert "needs.lint.outputs.code == 'true'" in jobs["test"]["if"]
    assert "needs.lint.outputs.pythons" in jobs["test"]["strategy"]["matrix"]["python-version"]
    assert "needs.lint.outputs.full == 'true'" in jobs["gate"]["if"]
    steps = [s.get("run", "") for s in jobs["lint"]["steps"]]
    assert any("scripts/ci_scope.py" in s for s in steps)
    assert any("check_language.py --commits" in s for s in steps)


def test_ci_caches_dependencies_and_the_plink2_download():
    jobs = _workflow("ci.yml")["jobs"]
    for name in ("lint", "test", "gate"):
        setup = [s for s in jobs[name]["steps"] if "setup-python" in s.get("uses", "")]
        assert setup and setup[0]["with"]["cache"] == "pip", name
    uses = [s.get("uses", "") for s in jobs["test"]["steps"]]
    assert "./.github/actions/install-plink2" in uses
    action = yaml.safe_load(PLINK2_ACTION.read_text())
    cache = [s for s in action["runs"]["steps"] if "actions/cache" in s.get("uses", "")][0]
    sha = re.search(r"PLINK2_SHA256: (\w+)", PLINK2_ACTION.read_text()).group(1)
    assert cache["with"]["key"].endswith(sha)  # a new pin never reuses an old zip
    assert "sha256sum -c" in PLINK2_ACTION.read_text()  # checked even when cached


def test_rebuild_runs_only_on_manual_trigger_and_weekly():
    wf = _workflow("rebuild.yml")
    assert set(wf["on"]) == {"workflow_dispatch", "schedule"}
    assert len(wf["on"]["schedule"]) == 1
    steps = wf["jobs"]["rebuild"]["steps"]
    plink = [s for s in steps if s.get("uses") == "./.github/actions/install-plink2"][0]
    assert plink["with"]["use-cache"] == "false"  # a clean rebuild downloads fresh
    assert any(s.get("run") == "scripts/rebuild.sh" for s in steps)


def test_rebuild_script_starts_from_a_clean_clone_without_the_pip_cache():
    text = (ROOT / "scripts" / "rebuild.sh").read_text()
    assert "git clone" in text and "--no-cache-dir" in text
    assert "hbfgwas.gate" in text and "pytest" in text
    assert (ROOT / "scripts" / "rebuild.sh").stat().st_mode & 0o111
