"""The forbidden-wording check must catch the phrases banned by docs/BRIEF.md.

The phrases are assembled from fragments so that this file does not itself
trip the repository-wide scan.
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import check_language  # noqa: E402

BANNED_EXAMPLES = [
    "the pipeline " + "re" + "discovers BCL11A",
    "this " + "valid" + "ates against " + "biology",
    "results " + "repli" + "cate prior work",
    "an independent " + "repli" + "cation cohort",
]
ALLOWED_EXAMPLES = [
    "the pipeline recovers planted effects",
    "validate the config file",
    "simulation repetitions",
    "biology is not claimed",
]


def test_banned_phrases_are_flagged():
    for text in BANNED_EXAMPLES:
        assert check_language.find_violations(text), text


def test_allowed_phrases_pass():
    for text in ALLOWED_EXAMPLES:
        assert not check_language.find_violations(text), text


def test_repository_files_and_commit_messages_are_clean():
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "check_language.py"), "--commits"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
