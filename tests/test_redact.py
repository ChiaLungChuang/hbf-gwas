"""Run records are redacted at write time: no hostname, no absolute local paths.

The end-to-end check on the synthetic fileset (every file the scripts write
into RUN_DIR) is in tests/test_scripts.py.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

from hbfgwas import experiment, redact

ROOT = Path(__file__).resolve().parents[1]
REPO = "/Users/alice/projects/hbf-gwas"
KW = {
    "root": REPO,
    "hostnames": ["lab-mac-01"],
    "home": "/Users/alice",
    "user": "alice",
    "tmp": ["/var/folders/xy/T"],
}

PLINK_LOG = f"""PLINK v2.0.0-b.1.4 64-bit (7 Nov 2026)
Options in effect:
  --out {REPO}/data/1kg/qc.1_geno
  --pfile {REPO}/data/1kg/afr

Hostname: lab-mac-01
Working directory: {REPO}
Start time: Wed Nov 11 17:50:01 2026
"""


def test_plink2_log_loses_hostname_and_absolute_paths():
    out = redact.redact_text(PLINK_LOG, **KW)
    assert "Hostname: [redacted]" in out
    assert "Working directory: .\n" in out
    assert "--out data/1kg/qc.1_geno" in out and "--pfile data/1kg/afr" in out
    assert "lab-mac-01" not in out and "/Users/" not in out and "alice" not in out
    # everything else is as produced
    assert "PLINK v2.0.0-b.1.4 64-bit (7 Nov 2026)" in out
    assert "Start time: Wed Nov 11 17:50:01 2026" in out


@pytest.mark.parametrize(
    ("text", "want"),
    [
        ("ran on lab-mac-01 today", "ran on [host] today"),
        (f"{REPO}/runs/2026-11-12/x.json", "runs/2026-11-12/x.json"),
        (f"cd {REPO} && ls", "cd . && ls"),
        (f"'{REPO}/config/analysis.yaml'", "'config/analysis.yaml'"),
        ("/Users/alice/.local/bin/plink2", "~/.local/bin/plink2"),
        ("/Users/alice/projects/hbf-gwas-old/x", "~/projects/hbf-gwas-old/x"),  # not the repo
        ("/Users/bob/data", "/Users/[user]/data"),
        ("/home/carol/x", "/home/[user]/x"),
        ("/tmp/pytest-of-alice/p0/run", "/tmp/pytest-of-[user]/p0/run"),
        ("~/alice-notes/x [TMP]/pytest-of-alice/y", "~/[user]-notes/x [TMP]/pytest-of-[user]/y"),
        (
            "file:///var/folders/xy/T/pytest-of-alice/serve/a.zst",
            "file://[TMP]/pytest-of-[user]/serve/a.zst",
        ),
        ("alice wrote 2 files", "alice wrote 2 files"),  # a name outside a path stays
        ("lambda 1.3422; 652 samples; NA20314", "lambda 1.3422; 652 samples; NA20314"),
        ("https://www.cog-genomics.org/plink/2.0/", "https://www.cog-genomics.org/plink/2.0/"),
    ],
)
def test_redaction_rules(text, want):
    assert redact.redact_text(text, **KW) == want


def test_data_and_run_dirs_outside_the_repo_are_named():
    aliases = {"DATA_DIR": "/scratch/x/data", "RUN_DIR": f"{REPO}/runs/2026-11-12"}
    out = redact.redact_text("/scratch/x/data/afr.pgen and /scratch/x/data", aliases=aliases, **KW)
    assert out == "[DATA_DIR]/afr.pgen and [DATA_DIR]"
    # inside the repository a named directory is simply repo-relative
    assert redact.redact_text(f"{REPO}/runs/2026-11-12/a", aliases=aliases, **KW) == (
        "runs/2026-11-12/a"
    )


def test_resolved_paths_are_redacted_too(tmp_path):
    real = tmp_path / "real"
    real.mkdir()
    link = tmp_path / "link"
    link.symlink_to(real)
    kw = {"root": link, "hostnames": [], "home": "/nonexistent-home", "user": None, "tmp": []}
    assert redact.redact_text(f"{real}/data/x", **kw) == "data/x"
    assert redact.redact_text(f"{link}/data/x", **kw) == "data/x"


def test_short_names_are_never_replaced():
    kw = {"root": REPO, "hostnames": ["box"], "home": "/Users/al", "user": "al", "tmp": []}
    assert redact.redact_text("box of al", **kw) == "box of al"


def test_redaction_is_idempotent_and_keeps_json_valid():
    doc = {"config_path": f"{REPO}/config/analysis.yaml", "lambda": 1.3422, "host": "lab-mac-01"}
    once = redact.redact_text(json.dumps(doc, indent=2), **KW)
    assert redact.redact_text(once, **KW) == once
    parsed = json.loads(once)
    assert parsed == {"config_path": "config/analysis.yaml", "lambda": 1.3422, "host": "[host]"}


def test_this_machine_is_redacted_with_the_defaults(tmp_path):
    text = f"{ROOT}/data/1kg/afr {Path.home()}/x " + " ".join(redact.local_hostnames())
    out = redact.redact_text(text)
    assert str(ROOT) not in out and str(Path.home()) not in out
    for h in redact.local_hostnames():
        assert h not in out
    assert out.startswith("data/1kg/afr ~/x")


def test_cli_writes_a_redacted_copy_and_the_note(tmp_path):
    src = tmp_path / "in.log"
    src.write_text(f"Hostname: x-host-1\n{tmp_path}/data/a\n")
    dest = tmp_path / "run" / "logs" / "in.log"
    cmd = [sys.executable, "-m", "hbfgwas.redact", "--repo-root", str(ROOT)]
    subprocess.run(
        [*cmd, "--alias", f"DATA_DIR={tmp_path / 'data'}", str(src), str(dest)], check=True
    )
    assert dest.read_text() == "Hostname: [redacted]\n[DATA_DIR]/a\n"
    subprocess.run([*cmd, "--note", str(tmp_path / "run")], check=True)
    note = (tmp_path / "run" / redact.NOTE_FILE).read_text()
    assert "redacted" in note and "when they were\nwritten" in note


def test_capture_shows_and_saves_redacted_output(tmp_path, capsys):
    with redact.capture(tmp_path, "step", **KW):
        print(f"wrote {REPO}/runs/x.json on lab-mac-01")
        print("warning", file=sys.stderr)
    shown = capsys.readouterr()
    assert shown.out == f"wrote {REPO}/runs/x.json on lab-mac-01\n"  # terminal unchanged
    assert shown.err == "warning\n"
    assert (tmp_path / "logs" / "step.stdout.txt").read_text() == "wrote runs/x.json on [host]\n"
    assert (tmp_path / "logs" / "step.stderr.txt").read_text() == "warning\n"
    assert (tmp_path / redact.NOTE_FILE).exists()


def test_simulated_results_record_a_repo_relative_config_path(tmp_path):
    prov = experiment.provenance()
    assert Path(prov["config_path"]).is_absolute()  # the in-memory value is unchanged
    path = redact.write_record(tmp_path / "results.json", json.dumps({"provenance": prov}))
    assert json.loads(path.read_text())["provenance"]["config_path"] == "config/analysis.yaml"
    src = (ROOT / "src" / "hbfgwas" / "experiment.py").read_text()
    assert 'redact.write_record(paths["results"]' in src and "redact.write_note(out_dir)" in src


def test_every_script_captures_its_output_and_redacts_its_copies():
    for f in sorted([*(ROOT / "scripts").glob("0*.sh"), ROOT / "scripts" / "step_e.sh"]):
        text = f.read_text()
        assert '\nsetup_dirs\ncapture_logs "$@"\n' in text, f.name
        assert not re.search(r'^\s*cp .*"\$RUN_DIR', text, re.M), f.name
    common = (ROOT / "scripts" / "_common.sh").read_text()
    assert 'save_log() { if [ -f "$1.log" ]; then redact_to' in common
    for f in ("03b_pc_plot.py", "05_report.py"):
        assert "redact.capture(" in (ROOT / "scripts" / f).read_text(), f
