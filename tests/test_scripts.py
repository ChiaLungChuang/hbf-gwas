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


def _fake_plink2(tmp: Path, version: str = "PLINK v2.0.0-b.1.4 fake stub (7 Nov 2026)") -> Path:
    """A stand-in plink2 that only answers --version (for tests that stop earlier)."""
    p = tmp / "fake_plink2"
    p.write_text(f"#!/bin/sh\nif [ \"$1\" = --version ]; then echo '{version}'; fi\nexit 0\n")
    p.chmod(0o755)
    return p


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
    r = _run(
        "00_download_1kg.sh", _env(tmp_path, PLINK2=_fake_plink2(tmp_path), MIN_FREE_GB=0, **urls)
    )
    assert r.returncode != 0
    assert "bytes" in r.stderr and "re-check the URL" in r.stderr


def test_download_rejects_an_html_page_even_when_large(tmp_path):
    urls = _serve(tmp_path, b"<!DOCTYPE html><html>" + b"x" * 5000)
    r = _run(
        "00_download_1kg.sh",
        _env(tmp_path, PLINK2=_fake_plink2(tmp_path), MIN_FREE_GB=0, MIN_BYTES_LARGE=10, **urls),
    )
    assert r.returncode != 0
    assert "HTML" in r.stderr


def test_download_rejects_a_non_zstd_file(tmp_path):
    urls = _serve(tmp_path, b"PK\x03\x04" + b"\1" * 5000)
    r = _run(
        "00_download_1kg.sh",
        _env(tmp_path, PLINK2=_fake_plink2(tmp_path), MIN_FREE_GB=0, MIN_BYTES_LARGE=10, **urls),
    )
    assert r.returncode != 0
    assert "zstd" in r.stderr


def test_download_checks_free_disk(tmp_path):
    urls = _serve(tmp_path, b"\x28\xb5\x2f\xfd" + b"\0" * 5000)
    r = _run(
        "00_download_1kg.sh",
        _env(tmp_path, PLINK2=_fake_plink2(tmp_path), MIN_FREE_GB=10**9, **urls),
    )
    assert r.returncode != 0
    assert "free" in r.stderr


def _unsigned_loci(tmp_path: Path) -> Path:
    out = tmp_path / "unsigned_loci.tsv"
    subprocess.run([sys.executable, "-m", "hbfgwas.loci", "export", "--out", str(out)], check=True)
    df = pd.read_csv(out, sep="\t")
    df.loc[0, "analysis_ready"] = False  # one locus not signed off
    df.to_csv(out, sep="\t", index=False)
    return out


def test_assoc_refuses_unverified_loci(tmp_path):
    r = _run(
        "04_assoc.sh",
        _env(tmp_path, PLINK2=_fake_plink2(tmp_path), HBFGWAS_LOCI_TSV=_unsigned_loci(tmp_path)),
    )
    assert r.returncode != 0
    assert "not signed off" in r.stderr


def _signed_loci(tmp_path: Path) -> Path:
    out = tmp_path / "signed_loci.tsv"
    subprocess.run([sys.executable, "-m", "hbfgwas.loci", "export", "--out", str(out)], check=True)
    return out


def test_assoc_refuses_an_empty_step_e_seed(tmp_path):
    empty = tmp_path / "seed.yaml"
    empty.write_text("seed: null\n")
    r = _run(
        "04_assoc.sh",
        _env(
            tmp_path,
            PLINK2=_fake_plink2(tmp_path),
            HBFGWAS_LOCI_TSV=_signed_loci(tmp_path),
            HBFGWAS_SKIP_PREREG_CHECK=1,
            STEP_E_SEED_FILE=empty,
        ),
    )
    assert r.returncode != 0
    assert "step e seed" in r.stderr


def test_assoc_refuses_without_prereg_tag(tmp_path):
    loci = tmp_path / "loci.tsv"
    subprocess.run([sys.executable, "-m", "hbfgwas.loci", "export", "--out", str(loci)], check=True)
    empty = tmp_path / "empty.git"  # a repository with no tags, whatever this checkout has
    subprocess.run(["git", "init", "-q", "--bare", str(empty)], check=True)
    r = _run(
        "04_assoc.sh",
        _env(tmp_path, PLINK2=_fake_plink2(tmp_path), HBFGWAS_LOCI_TSV=loci, GIT_DIR=empty),
    )
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
    cfg["prediction"].update(n_null_variants=500, power_shift_draws=200, lambda_shift_draws=50)
    # test-only declared shift and sensitivity exclusion on the synthetic samples
    cfg["real_data_simulation"]["population_shift"] = {
        "ACB": -0.3,
        "ASW": -0.5,
        "ESN": 0.2,
        "GWD": 0.6,
        "LWK": -0.1,
        "MSL": 0.4,
        "YRI": 0.1,
    }
    cfg["sensitivity"]["exclude_samples"] = ["SYN_ASW_05"]
    # the real-data criteria (predictions for the 652 real genotypes) do not apply
    # to the synthetic fileset
    cfg["acceptance"]["real"] = None
    seed_file = tmp / "step_e_seed.yaml"
    seed_file.write_text("seed: 777\n")
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
        STEP_E_SEED_FILE=seed_file,
    )
    out = {}
    for script in (
        "00_download_1kg.sh",
        "01_subset_afr.sh",
        "02_qc.sh",
        "03_pca.sh",
        "03b_pc_plot.py",
        "03d_sensitivity_pca.sh",
        "03c_prereg_predictions.sh",
    ):
        r = _run(script, env)
        out[script] = r
        if r.returncode != 0:
            pytest.fail(
                f"{script} failed ({r.returncode})\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
            )
    # Test-only real-data criteria for the synthetic fileset, written from its own
    # step d predictions (the real ones are pre-registered in the step d PR). The
    # lambda and type I bounds are deliberately wide: this checks the step e code
    # path, not the synthetic numbers. The power calibration is the real rule.
    pred = json.loads((tmp / "run" / "prereg_predictions.json").read_text())
    cfg["acceptance"]["real"] = {
        "criteria_version": "synthetic-test",
        "unadjusted_lambda_gc": {"at_least": 0.0, "at_most": 100.0},
        "adjusted_lambda_gc": {"at_least": 0.0, "at_most": 100.0},
        "adjusted_type1_error": {"alpha": 0.05, "at_least": 0.0, "at_most": 1.0},
        "power_calibration": {"max_abs_diff": 0.05, "n_redraws": 2000},
        "predicted_power": {
            g: round(v["pc_adjusted"]["power"], 4) for g, v in pred["conditional_power"].items()
        },
        "required_loci": [],
    }
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    r = _run("step_e.sh", env)
    out["step_e.sh"] = r
    if r.returncode not in (0, 3):  # 3: a gated criterion failed, which the tests inspect
        pytest.fail(f"step_e.sh failed ({r.returncode})\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}")
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
    for vid in synth.KEPT + synth.KEPT_BY_KEEP_FEWHET:
        assert vid in ids, vid
    assert t.loc[3, "filter"].endswith("keep-fewhet")


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
    # review item E: the synthetic samples are unrelated, so KING removes nobody and says so
    assert "--king-cutoff 0.0884 removed 0 samples" in e2e["out"]["03_pca.sh"].stderr


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
    sources = {x["predicted_power_source"] for x in res["loci"]}
    assert sources == {"pre-registered (acceptance.real)"}
    assert (e2e["run"] / "qq_real.png").stat().st_size > 10_000


def test_e2e_pc_plot_at_step_c(e2e):
    # review item F: its own script, run right after 03_pca.sh
    assert (e2e["run"] / "pc1_pc2_by_population.png").stat().st_size > 10_000
    counts = json.loads((e2e["run"] / "pc_plot_counts.json").read_text())
    assert set(counts) == set(synth.AFR)  # the seven AFR populations, not SuperPop
    assert sum(counts.values()) == synth.EXPECTED_AFR_AFTER_REL


def test_e2e_gate_refuses_real_data_before_preregistration(e2e, tmp_path):
    cfg = yaml.safe_load(e2e["config"].read_text())
    cfg["acceptance"]["real"] = None  # as on main before step d
    unregistered = tmp_path / "unregistered.yaml"
    unregistered.write_text(yaml.safe_dump(cfg, sort_keys=False))
    r = subprocess.run(
        [
            sys.executable,
            "-m",
            "hbfgwas.gate",
            "--config",
            str(unregistered),
            "--experiment",
            "real",
            "--results",
            str(e2e["run"] / "results.json"),
        ],
        capture_output=True,
        text=True,
    )
    assert r.returncode == 2


def test_e2e_step_e_power_redraws(e2e):
    rd = json.loads((e2e["run"] / "power_redraws.json").read_text())
    assert rd["n_redraws"] == 2000 and rd["seed_child"] == "power_reps"
    assert rd["step_e_seed"] == 777
    declared = yaml.safe_load(e2e["config"].read_text())["acceptance"]["real"]["predicted_power"]
    assert set(rd["loci"]) == set(declared)
    for gene, r in rd["loci"].items():
        print(gene, r["empirical_power"], r["empirical_power_ci95"], declared[gene])
        # the run recomputes the prediction from the same genotypes: it equals the declared value
        assert round(r["predicted_power_recomputed"], 4) == declared[gene]
        assert abs(r["empirical_power"] - declared[gene]) <= 0.05


def test_e2e_step_e_sensitivity_scan(e2e):
    for f in ("glm_sens_nopc.HBF_SIM.glm.linear", "glm_sens_pc.HBF_SIM.glm.linear"):
        assert (e2e["data"] / f).exists(), f
    res = json.loads((e2e["run"] / "results.json").read_text())
    sens = res["sensitivity"]
    assert sens["status"] == "run" and sens["gated"] is False
    assert sens["excluded_samples"] == ["SYN_ASW_05"]
    assert sens["n_samples_analysed"] == res["n_samples_analysed"] - 1
    assert set(sens["loci"]) == {x["locus"] for x in res["loci"]}
    for side in ("unadjusted", "adjusted"):
        assert sens[side]["lambda_gc"] > 0


def test_e2e_step_e_gate_judges_every_criterion(e2e):
    gate = json.loads((e2e["run"] / "gate.json").read_text())
    names = [c["name"] for c in gate["checks"]]
    for want in ("unadjusted lambda_GC", "PC-adjusted lambda_GC", "PC-adjusted type I error"):
        assert want in names
    calib = [c for c in gate["checks"] if c["name"].endswith("power calibration")]
    assert len(calib) == 3
    assert all(c["status"] == "PASS" for c in calib), calib
    assert e2e["out"]["step_e.sh"].returncode == (0 if gate["passed"] else 3)


def _step_e_config(tmp_path: Path, real=True, shift=True) -> Path:
    cfg = load_config()
    cfg["acceptance"]["real"] = (
        {
            "unadjusted_lambda_gc": {"at_least": 0.0, "at_most": 100.0},
            "adjusted_lambda_gc": {"at_least": 0.0, "at_most": 100.0},
            "adjusted_type1_error": {"alpha": 0.05, "at_least": 0.0, "at_most": 1.0},
            "power_calibration": {"max_abs_diff": 0.05, "n_redraws": 10},
        }
        if real
        else None
    )
    cfg["real_data_simulation"]["population_shift"] = (
        {p: 0.0 for p in ("ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI")} if shift else None
    )
    path = tmp_path / f"step_e_config_{real}_{shift}.yaml"  # one file per variant
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return path


def _seed(tmp_path: Path, value="1234") -> Path:
    p = tmp_path / f"seed_{value}.yaml"  # one file per value
    p.write_text(f"seed: {value}\n")
    return p


def _step_e_env(tmp_path: Path, **extra) -> dict:
    base = {
        "PLINK2": _fake_plink2(tmp_path),
        "CONFIG": _step_e_config(tmp_path),
        "HBFGWAS_LOCI_TSV": _signed_loci(tmp_path),
        "HBFGWAS_SKIP_PREREG_CHECK": 1,
        "STEP_E_SEED_FILE": _seed(tmp_path),
    }
    base.update(extra)
    return _env(tmp_path, **base)


def test_step_e_refuses_without_the_prereg_tag(tmp_path):
    empty = tmp_path / "empty.git"  # a repository with no tags, whatever this checkout has
    subprocess.run(["git", "init", "-q", "--bare", str(empty)], check=True)
    env = _step_e_env(tmp_path, GIT_DIR=empty)
    env.pop("HBFGWAS_SKIP_PREREG_CHECK")
    r = _run("step_e.sh", env)
    assert r.returncode == 1 and "prereg-v1" in r.stderr


def test_step_e_refuses_an_unsigned_locus_table(tmp_path):
    r = _run("step_e.sh", _step_e_env(tmp_path, HBFGWAS_LOCI_TSV=_unsigned_loci(tmp_path)))
    assert r.returncode == 1 and "not signed off" in r.stderr


def test_step_e_refuses_an_empty_seed(tmp_path):
    r = _run("step_e.sh", _step_e_env(tmp_path, STEP_E_SEED_FILE=_seed(tmp_path, "null")))
    assert r.returncode == 1 and "step e seed" in r.stderr


def test_step_e_refuses_without_real_criteria(tmp_path):
    r = _run("step_e.sh", _step_e_env(tmp_path, CONFIG=_step_e_config(tmp_path, real=False)))
    assert r.returncode == 1 and "acceptance.real is empty" in r.stderr


def test_step_e_refuses_an_undeclared_shift(tmp_path):
    r = _run("step_e.sh", _step_e_env(tmp_path, CONFIG=_step_e_config(tmp_path, shift=False)))
    assert r.returncode == 1 and "population_shift" in r.stderr


def test_step_e_runs_once(tmp_path):
    (tmp_path / "run").mkdir()
    (tmp_path / "run" / "results.json").write_text("{}")
    r = _run("step_e.sh", _step_e_env(tmp_path))
    assert r.returncode == 1 and "runs once" in r.stderr


def test_step_e_refusals_happen_before_any_phenotype(tmp_path):
    r = _run("step_e.sh", _step_e_env(tmp_path, STEP_E_SEED_FILE=_seed(tmp_path, "null")))
    assert r.returncode == 1
    assert not (tmp_path / "data" / "pheno.tsv").exists()
    assert not (tmp_path / "run" / "phenotype_meta.json").exists()


def test_e2e_every_plink2_log_is_kept(e2e):
    logs = {p.name for p in (e2e["run"] / "logs").glob("*.log")}
    for expected in ("afr.log", "afr_qc.log", "pca.3.log", "glm_pc.log", "glm_nopc.log"):
        assert expected in logs, sorted(logs)


def test_e2e_run_records_are_redacted_at_write_time(e2e):
    from hbfgwas import redact

    run, tmp = e2e["run"], e2e["tmp"]
    forbidden = {str(tmp), os.path.realpath(tmp), str(ROOT), str(Path.home())}
    forbidden |= set(redact.local_hostnames())
    user = redact.local_user()
    checked = 0
    for f in sorted(run.rglob("*")):
        if not f.is_file() or f.suffix == ".png":
            continue
        text = f.read_text(errors="replace")
        for s in forbidden:
            assert s not in text, (f.relative_to(run), s)
        assert not re.search(r"/(Users|home)/(?!\[user\])", text), f.relative_to(run)
        if user:
            for token in re.findall(r"/[^\s\"'`]+", text):
                assert not re.search(rf"(?<![A-Za-z0-9]){re.escape(user)}(?![A-Za-z0-9])", token)
        checked += 1
    assert checked > 30
    assert (run / redact.NOTE_FILE).exists()
    logs = run / "logs"
    for script in e2e["out"]:
        name = script.rsplit(".", 1)[0]
        assert (logs / f"{name}.stdout.txt").exists(), name
        assert (logs / f"{name}.stderr.txt").exists(), name
    plink_logs = "".join(p.read_text() for p in logs.glob("*.log"))
    assert "Hostname: [redacted]" in plink_logs
    assert "Working directory: .\n" in plink_logs  # the scripts ran from the repository
    assert "[DATA_DIR]/afr" in plink_logs  # DATA_DIR is outside the repository here


def test_e2e_captured_logs_are_what_the_script_printed_redacted(e2e):
    from hbfgwas import redact

    aliases = {"DATA_DIR": e2e["data"], "RUN_DIR": e2e["run"]}
    r = e2e["out"]["02_qc.sh"]
    logs = e2e["run"] / "logs"
    assert (logs / "02_qc.stderr.txt").read_text() == redact.redact_text(r.stderr, aliases=aliases)
    assert (logs / "02_qc.stdout.txt").read_text() == redact.redact_text(r.stdout, aliases=aliases)
    assert "missingness" in (logs / "02_qc.stderr.txt").read_text()


def test_e2e_prereg_predictions_before_any_phenotype(e2e):
    # review item H: predicted lambda and conditional power exist before 04 runs
    pred = json.loads((e2e["run"] / "prereg_predictions.json").read_text())
    assert pred["computed_before_any_phenotype"] is True
    assert set(pred["conditional_power"]) == {"BCL11A", "HBS1L-MYB", "HBB cluster"}
    for side in ("unadjusted", "pc_adjusted"):
        lam = pred["predicted_lambda_gc"][side]
        assert 0.5 < lam["q05"] <= lam["median"] <= lam["q95"] < 10
    assert pred["n_null_variants"] > 400
    assert pred["not_testable"] == []
    mtime_pred = (e2e["run"] / "prereg_predictions.json").stat().st_mtime
    mtime_pheno = (e2e["data"] / "pheno.tsv").stat().st_mtime
    assert mtime_pred < mtime_pheno  # written before the phenotype existed


def test_scripts_refuse_an_unpinned_plink2(tmp_path):
    # review item G: config/analysis.yaml pins plink2 v2.0.0-b.1.4 (7 Nov 2026)
    other = _fake_plink2(tmp_path, "PLINK v2.0.0-a.7.15 64-bit (7 Nov 2026)")
    r = _run("02_qc.sh", _env(tmp_path, PLINK2=other))
    assert r.returncode != 0
    assert "pins v2.0.0-b.1.4 (7 Nov 2026)" in r.stderr


def test_e2e_sensitivity_pca_and_predictions(e2e):
    # NA20314-style sensitivity analysis, on a synthetic sample (not gated)
    assert _n(e2e["run"] / "sens_pca_eigenval.txt") == 10
    sens = json.loads((e2e["run"] / "prereg_predictions_sensitivity.json").read_text())
    prim = json.loads((e2e["run"] / "prereg_predictions.json").read_text())
    assert sens["excluded_samples"] == ["SYN_ASW_05"]
    assert sens["n_samples"] == prim["n_samples"] - 1
    assert prim["shift_mode"] == "fixed" and sens["shift_mode"] == "fixed"
    assert prim["population_shift"]["GWD"] == 0.6
