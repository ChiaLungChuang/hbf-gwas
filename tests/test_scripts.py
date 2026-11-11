"""Tests for scripts/00-05.

* Shell syntax and the guard rails (download checks, 04 refusing to run on
  unverified loci or without the prereg-v1 tag) run anywhere.
* The end-to-end run of 00 -> 05 on the synthetic fileset (tests/synthetic_1kg.py)
  needs plink2 and zstd on PATH. It is skipped without them, unless
  REQUIRE_PLINK2=1 (set in CI), where a missing tool is a failure.
"""

from __future__ import annotations

import copy
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

from hbfgwas.config import load_config

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(Path(__file__).resolve().parent))
import synthetic_1kg as synth  # noqa: E402

HAVE_TOOLS = bool(shutil.which("plink2") and shutil.which("zstd"))
REQUIRE = os.environ.get("REQUIRE_PLINK2") == "1"


def _env(tmp: Path, **extra) -> dict:
    env = dict(os.environ)
    env.update(
        PYTHON=sys.executable,
        DATA_DIR=str(tmp / "data"),
        RUN_DIR=str(tmp / "run"),
        CONFIG=str(ROOT / "config" / "analysis.yaml"),
    )
    env.update({k: str(v) for k, v in extra.items()})
    return env


def _run(script: str, env: dict) -> subprocess.CompletedProcess:
    cmd = (
        [str(SCRIPTS / script)]
        if script.endswith(".sh")
        else [sys.executable, str(SCRIPTS / script)]
    )
    return subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True)


# --- anywhere --------------------------------------------------------------------


@pytest.mark.parametrize("script", sorted(p.name for p in SCRIPTS.glob("0*.sh")))
def test_shell_syntax(script):
    subprocess.run(["bash", "-n", str(SCRIPTS / script)], check=True)


def _serve(tmp: Path, pgen_bytes: bytes) -> dict:
    d = tmp / "serve"
    d.mkdir()
    (d / "pgen").write_bytes(pgen_bytes)
    (d / "pvar").write_bytes(b"\x28\xb5\x2f\xfd" + b"\0" * 2000)
    (d / "psam").write_text("#IID\tSuperPop\n" + "x\tAFR\n" * 50)
    (d / "rel").write_text("#IID\n" + "a\n" * 50)
    return {f"{k.upper()}_URL": (d / k).as_uri() for k in ("pgen", "pvar", "psam", "rel")}


def test_download_rejects_a_small_file(tmp_path):
    urls = _serve(tmp_path, b"<html>gone</html>")
    r = _run("00_download_1kg.sh", _env(tmp_path, PLINK2="true", MIN_FREE_GB=0, **urls))
    assert r.returncode != 0
    assert "bytes" in r.stderr and "re-check the URL" in r.stderr


def test_download_rejects_an_html_page_even_when_large(tmp_path):
    urls = _serve(tmp_path, b"<!DOCTYPE html><html>" + b"x" * 5000)
    r = _run(
        "00_download_1kg.sh",
        _env(tmp_path, PLINK2="true", MIN_FREE_GB=0, MIN_BYTES_LARGE=10, **urls),
    )
    assert r.returncode != 0
    assert "HTML" in r.stderr


def test_download_rejects_a_non_zstd_file(tmp_path):
    urls = _serve(tmp_path, b"PK\x03\x04" + b"\1" * 5000)
    r = _run(
        "00_download_1kg.sh",
        _env(tmp_path, PLINK2="true", MIN_FREE_GB=0, MIN_BYTES_LARGE=10, **urls),
    )
    assert r.returncode != 0
    assert "zstd" in r.stderr


def test_download_checks_free_disk(tmp_path):
    urls = _serve(tmp_path, b"\x28\xb5\x2f\xfd" + b"\0" * 5000)
    r = _run("00_download_1kg.sh", _env(tmp_path, PLINK2="true", MIN_FREE_GB=10**9, **urls))
    assert r.returncode != 0
    assert "free" in r.stderr


def test_assoc_refuses_unverified_loci(tmp_path):
    r = _run("04_assoc.sh", _env(tmp_path, PLINK2="true"))
    assert r.returncode != 0
    assert "not signed off" in r.stderr


def test_assoc_refuses_without_prereg_tag(tmp_path):
    loci = tmp_path / "loci.tsv"
    subprocess.run([sys.executable, "-m", "hbfgwas.loci", "export", "--out", str(loci)], check=True)
    empty = tmp_path / "empty.git"  # a repository with no tags, whatever this checkout has
    subprocess.run(["git", "init", "-q", "--bare", str(empty)], check=True)
    r = _run("04_assoc.sh", _env(tmp_path, PLINK2="true", HBFGWAS_LOCI_TSV=loci, GIT_DIR=empty))
    assert r.returncode != 0
    assert "prereg-v1" in r.stderr


# --- end to end on the synthetic fileset (plink2) --------------------------------


@pytest.fixture(scope="module")
def e2e(tmp_path_factory):
    if not HAVE_TOOLS:
        if REQUIRE:
            pytest.fail("REQUIRE_PLINK2=1 but plink2 or zstd is not on PATH")
        pytest.skip("plink2 and zstd are needed for the end-to-end script test")
    tmp = tmp_path_factory.mktemp("e2e")
    built = synth.build(tmp / "synthetic")
    cfg = copy.deepcopy(load_config())
    cfg["qc"]["expected_afr_samples"] = {"min": 130, "max": 140}
    cfg["concordance"]["n_random_variants"] = 300
    cfg_path = tmp / "synthetic_config.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    loci = tmp / "loci.tsv"
    subprocess.run([sys.executable, "-m", "hbfgwas.loci", "export", "--out", str(loci)], check=True)
    s = built["serve"]
    env = _env(
        tmp,
        CONFIG=cfg_path,
        PGEN_URL=(s / "all_phase3.pgen.zst").as_uri(),
        PVAR_URL=(s / "all_phase3.pvar.zst").as_uri(),
        PSAM_URL=(s / "all_phase3.psam").as_uri(),
        REL_URL=(s / "deg2_phase3.king.cutoff.out.id").as_uri(),
        MIN_FREE_GB=0,
        MIN_BYTES_LARGE=100,
        MIN_BYTES_SMALL=10,
        HBFGWAS_LOCI_TSV=loci,
        HBFGWAS_SKIP_PREREG_CHECK=1,
    )
    out = {}
    for script in (
        "00_download_1kg.sh",
        "01_subset_afr.sh",
        "02_qc.sh",
        "03_pca.sh",
        "04_assoc.sh",
        "05_report.py",
    ):
        r = _run(script, env)
        out[script] = r
        if r.returncode != 0:
            pytest.fail(
                f"{script} failed ({r.returncode})\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
            )
    return {
        "tmp": tmp,
        "data": tmp / "data",
        "run": tmp / "run",
        "out": out,
        "built": built,
        "serve": s,
        "config": cfg_path,
    }


def _ids(pvar: Path) -> set[str]:
    return {
        line.split("\t")[2] for line in pvar.read_text().splitlines() if not line.startswith("#")
    }


def _n(path: Path) -> int:
    return sum(1 for line in path.read_text().splitlines() if line and not line.startswith("#"))


def test_e2e_download_manifest(e2e):
    m = pd.read_csv(e2e["run"] / "download_manifest.tsv", sep="\t")
    assert list(m.columns) == ["file", "url", "downloaded_utc", "bytes", "sha256"]
    assert len(m) == 4
    for _, r in m.iterrows():
        src = e2e["serve"] / r["file"]
        assert r["sha256"] == hashlib.sha256(src.read_bytes()).hexdigest()
        assert r["bytes"] == src.stat().st_size
    assert (e2e["data"] / "all_phase3.pgen").read_bytes() == e2e["built"]["pgen_plain"].read_bytes()
    assert (e2e["run"] / "plink2_version.txt").read_text().strip()


def test_e2e_afr_subset(e2e):
    assert _n(e2e["data"] / "afr.psam") == synth.EXPECTED_AFR_AFTER_REL
    ids = _ids(e2e["data"] / "afr.pvar")
    for vid in synth.REMOVED_BY_01:
        assert vid not in ids, vid
    assert not any(v.startswith(("rsX", "rsINDEL", "rsMULTI", "rsSPLIT")) for v in ids)
    assert sum(bool(re.fullmatch(r"1:\d+:C:T", v)) for v in ids) == 10  # '.' IDs filled
    counts = pd.read_csv(e2e["run"] / "01_subset_counts.tsv", sep="\t")
    assert counts["step"].tolist()[-1] == "afr_dedup"


def test_e2e_qc_table_and_missingness_message(e2e):
    t = pd.read_csv(e2e["run"] / "qc_filters.tsv", sep="\t")
    assert t["filter"].str.split().str[0].tolist() == ["--geno", "--mind", "--maf", "--hwe"]
    assert t.loc[0, "variants_removed"] == 0 and t.loc[1, "samples_removed"] == 0
    assert t.loc[2, "variants_removed"] >= len(synth.REMOVED_BY_MAF)
    assert t.loc[3, "variants_removed"] >= len(synth.REMOVED_BY_HWE)
    assert "fully called" in e2e["out"]["02_qc.sh"].stderr
    ids = _ids(e2e["data"] / "afr_qc.pvar")
    for vid in synth.REMOVED_BY_MAF + synth.REMOVED_BY_HWE:
        assert vid not in ids, vid
    for vid in synth.KEPT:
        assert vid in ids, vid


def test_e2e_sentinel_report(e2e):
    s = pd.read_csv(e2e["run"] / "sentinels.tsv", sep="\t").set_index("gene")
    assert s.at["HBB cluster", "found_by"] == "position"
    assert s.at["BCL11A", "found_by"] == "rsid"
    assert s["survived_qc"].all()
    assert s["effect_allele_freq"].between(0, 1).all()


def test_e2e_pca(e2e):
    assert _n(e2e["run"] / "pca_eigenval.txt") == 10
    eig = pd.read_csv(e2e["data"] / "pca.eigenvec", sep="\t")
    assert [f"PC{i}" for i in range(1, 11)] == [c for c in eig.columns if c.startswith("PC")]
    assert (e2e["run"] / "het_summary.tsv").exists()
    counts = pd.read_csv(e2e["run"] / "03_pca_counts.tsv", sep="\t")
    assert counts["step"].tolist() == [
        "input",
        "outside_long_range_ld_and_maf",
        "ld_pruned",
        "king_unrelated",
    ]
    # some synthetic SNPs sit inside the long-range LD regions and must be gone
    assert counts.loc[1, "variants"] < counts.loc[0, "variants"]


def test_e2e_association_outputs(e2e):
    d = e2e["data"]
    for f in (
        "glm_nopc.HBF_SIM.glm.linear",
        "glm_pc.HBF_SIM.glm.linear",
        "pheno.tsv",
        "concord.raw",
    ):
        assert (d / f).exists(), f
    meta = json.loads((e2e["run"] / "phenotype_meta.json").read_text())
    assert [m["gene"] for m in meta["planted"]] == ["BCL11A", "HBS1L-MYB", "HBB cluster"]
    assert "simulated" in meta["phenotype"]


def test_e2e_report_and_plink2_concordance(e2e):
    res = json.loads((e2e["run"] / "results.json").read_text())
    print(json.dumps({k: res[k] for k in ("unadjusted", "adjusted", "concordance")}, indent=1))
    print(pd.read_csv(e2e["run"] / "qc_filters.tsv", sep="\t").to_string(index=False))
    assert res["concordance"]["passed"], res["concordance"]
    for side in ("unadjusted", "adjusted"):
        assert res["concordance"][side]["n_compared"] >= 300
        assert res[side]["lambda_gc"] > 0
    assert len(res["loci"]) == 3
    assert all("NOT pre-registered" in x["predicted_power_source"] for x in res["loci"])
    assert (e2e["run"] / "qq_real.png").stat().st_size > 10_000
    assert (e2e["run"] / "pc1_pc2_by_population.png").stat().st_size > 10_000


def test_e2e_gate_refuses_real_data_before_preregistration(e2e):
    r = subprocess.run(
        [
            sys.executable,
            "-m",
            "hbfgwas.gate",
            "--config",
            str(e2e["config"]),
            "--experiment",
            "real",
            "--results",
            str(e2e["run"] / "results.json"),
        ],
        capture_output=True,
        text=True,
    )
    assert r.returncode == 2


def test_e2e_every_plink2_log_is_kept(e2e):
    logs = {p.name for p in (e2e["run"] / "logs").glob("*.log")}
    for expected in ("afr.log", "afr_qc.log", "pca.3.log", "glm_pc.log", "glm_nopc.log"):
        assert expected in logs, sorted(logs)
