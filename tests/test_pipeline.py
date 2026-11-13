"""Unit tests for the helpers the real-data scripts call (no plink2 needed)."""

import json
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from hbfgwas import pipeline
from hbfgwas.config import load_config

CFG = load_config()

# A test-only declared shift for the seven AFR populations (the real values are
# pre-registered in config/analysis.yaml in step d).
TEST_SHIFT = {"ACB": -0.3, "ASW": -0.5, "ESN": 0.2, "GWD": 0.6, "LWK": -0.1, "MSL": 0.4, "YRI": 0.1}


def _cfg_with_shift():
    import copy

    cfg = copy.deepcopy(CFG)
    cfg["real_data_simulation"]["population_shift"] = TEST_SHIFT
    return cfg


@pytest.fixture
def seed_file(tmp_path):
    p = tmp_path / "step_e_seed.yaml"
    p.write_text("seed: 4242\n")
    return p


PVAR = """##fileformat=PVARv1.0
#CHROM\tPOS\tID\tREF\tALT
2\t60718000\trsA\tC\tT
2\t60718043\trs1427407\tT\tG
6\t135419018\trs9399137\tT\tC
11\t5276169\trs41322148\tG\tA
11\t5300000\trsB\tA\tG
"""


@pytest.fixture
def pvar(tmp_path):
    p = tmp_path / "x.pvar"
    p.write_text(PVAR)
    return p


@pytest.fixture
def loci(tmp_path):
    out = tmp_path / "loci.tsv"
    subprocess.run([sys.executable, "-m", "hbfgwas.loci", "export", "--out", str(out)], check=True)
    return pd.read_csv(out, sep="\t", dtype={"chromosome": str, "rsid": str})


def test_config_get_and_format():
    assert pipeline.format_scalar(pipeline.config_get(CFG, "qc.maf")) == "0.01"
    assert pipeline.format_scalar(pipeline.config_get(CFG, "qc.hwe_p")) == "1e-06"
    assert pipeline.format_scalar(pipeline.config_get(CFG, "pca.n_pcs")) == "10"
    assert pipeline.format_scalar(pipeline.config_get(CFG, "acceptance.real")) == "null"
    with pytest.raises(KeyError):
        pipeline.config_get(CFG, "qc.nope")


def test_config_get_cli():
    out = subprocess.run(
        [sys.executable, "-m", "hbfgwas.pipeline", "config-get", "qc.expected_afr_samples.min"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    assert out == "640"


def test_find_sentinels_by_rsid_and_by_position(pvar, loci):
    smap = pipeline.find_sentinels(pvar, loci).set_index("gene")
    assert smap.at["BCL11A", "found_by"] == "rsid"
    assert smap.at["HBS1L-MYB", "found_by"] == "rsid"
    # rs7482144 is stored under a synonym: found by GRCh37 position
    assert smap.at["HBB cluster", "found_by"] == "position"
    assert smap.at["HBB cluster", "id"] == "rs41322148"
    assert smap["position_matches_table"].all()
    assert smap["alleles_match_table"].all()


def test_find_sentinels_reports_missing_and_mismatch(tmp_path, loci):
    p = tmp_path / "y.pvar"
    p.write_text("#CHROM\tPOS\tID\tREF\tALT\n2\t99\trs1427407\tT\tG\n")
    smap = pipeline.find_sentinels(p, loci).set_index("gene")
    assert smap.at["BCL11A", "position_matches_table"] == False  # noqa: E712
    assert smap.at["HBS1L-MYB", "found_by"] == "not found"


def test_sentinel_report(pvar, loci):
    smap = pipeline.find_sentinels(pvar, loci)
    afreq = pd.DataFrame(
        {
            "ID": ["rs1427407", "rs9399137", "rs41322148"],
            "ALT_FREQS": ["0.76", "0.04", "0.005"],
            "OBS_CT": ["1300"] * 3,
        }
    )
    steps = {
        "geno": {"rs1427407", "rs9399137", "rs41322148"},
        "maf": {"rs1427407", "rs9399137"},
        "hwe": {"rs1427407", "rs9399137"},
    }
    rep = pipeline.sentinel_report(smap, afreq, steps).set_index("gene")
    # effect allele T is REF for rs1427407 -> 1 - ALT freq
    assert rep.at["BCL11A", "effect_allele_freq"] == pytest.approx(0.24)
    assert rep.at["HBS1L-MYB", "effect_allele_freq"] == pytest.approx(0.04)
    assert rep.at["HBB cluster", "removed_by"] == "maf"
    assert rep.at["BCL11A", "survived_qc"]
    assert not rep.at["HBB cluster", "survived_qc"]
    # review item C: failing QC means "not testable in this callset", no proxy
    assert rep.at["HBB cluster", "status"] == "not testable in this callset"
    assert rep.at["BCL11A", "status"] == "testable"


def test_absent_sentinel_is_not_testable(tmp_path, loci):
    p = tmp_path / "no_hbb.pvar"
    p.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\n2\t60718043\trs1427407\tT\tG\n"
        "6\t135419018\trs9399137\tT\tC\n11\t5276170\trsNEAR\tC\tT\n"
    )
    smap = pipeline.find_sentinels(p, loci)
    rep = pipeline.sentinel_report(smap, None, {"geno": {"rs1427407", "rs9399137", "rsNEAR"}})
    rep = rep.set_index("gene")
    # a variant 1 bp away is not used as a proxy
    assert rep.at["HBB cluster", "found_by"] == "not found"
    assert rep.at["HBB cluster", "status"] == pipeline.NOT_TESTABLE
    assert rep.at["HBB cluster", "removed_by"] == "not in data"


def test_sample_ids_is_reproducible_and_includes(pvar):
    a = pipeline.sample_ids(pvar, 2, seed=1, include={"rs9399137"})
    b = pipeline.sample_ids(pvar, 2, seed=1, include={"rs9399137"})
    assert a == b
    assert "rs9399137" in a
    assert len(set(a)) in (2, 3)


def _raw(tmp_path, counted=("G", "C", "A"), n=60, seed=0):
    rng = np.random.default_rng(seed)
    ids = ["rs1427407", "rs9399137", "rs41322148"]
    g = rng.integers(0, 3, size=(n, 3))
    iids = [f"S{i}" for i in range(n)]
    df = pd.DataFrame({"FID": iids, "IID": iids, "PAT": 0, "MAT": 0, "SEX": 1, "PHENOTYPE": -9})
    for j, (vid, a) in enumerate(zip(ids, counted, strict=True)):
        df[f"{vid}_{a}"] = g[:, j]
    path = tmp_path / "planted.raw"
    df.to_csv(path, sep="\t", index=False)
    pops = ["YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"]
    psam = pd.DataFrame(
        {"IID": iids, "SuperPop": "AFR", "Population": [pops[i % 7] for i in range(n)]}
    )
    return path, psam, g


def test_simulate_phenotype_flips_to_the_effect_allele(tmp_path, pvar, loci, seed_file):
    path, psam, g = _raw(tmp_path)  # export counted G (ALT) for rs1427407; effect allele is T
    smap = pipeline.find_sentinels(pvar, loci)
    pheno, meta = pipeline.simulate_phenotype_from_raw(
        _cfg_with_shift(), pipeline.read_raw(path), smap, psam, seed_file
    )
    planted = {m["gene"]: m for m in meta["planted"]}
    assert planted["BCL11A"]["flipped_to_effect_allele"] is True
    assert planted["HBS1L-MYB"]["flipped_to_effect_allele"] is False
    assert planted["BCL11A"]["effect_allele_freq"] == pytest.approx((2 - g[:, 0]).mean() / 2)
    assert list(pheno.columns) == ["FID", "IID", "HBF_SIM"] or list(pheno.columns) == [
        "IID",
        "HBF_SIM",
    ]
    assert meta["step_e_seed"] == 4242
    assert meta["population_shift"] == TEST_SHIFT  # exactly the declared values
    again, _ = pipeline.simulate_phenotype_from_raw(
        _cfg_with_shift(), pipeline.read_raw(path), smap, psam, seed_file
    )
    np.testing.assert_array_equal(pheno["HBF_SIM"], again["HBF_SIM"])


def test_simulate_phenotype_rejects_an_allele_mismatch(tmp_path, pvar, loci, seed_file):
    path, psam, _ = _raw(tmp_path)
    smap = pipeline.find_sentinels(pvar, loci)
    smap.loc[smap["gene"] == "HBS1L-MYB", ["ref", "alt"]] = ["A", "G"]  # strand flip
    with pytest.raises(ValueError, match="neither REF"):
        pipeline.simulate_phenotype_from_raw(
            _cfg_with_shift(), pipeline.read_raw(path), smap, psam, seed_file
        )


def test_write_pheno_header(tmp_path):
    out = tmp_path / "p.tsv"
    pipeline.write_pheno(pd.DataFrame({"IID": ["a", "b"], "HBF_SIM": [0.1, -0.2]}), out)
    assert out.read_text().splitlines()[0] == "#IID\tHBF_SIM"


def test_het_summary():
    het = pd.DataFrame(
        {"IID": [f"s{i}" for i in range(50)] + ["odd"], "F": [0.0] * 25 + [0.01] * 25 + [0.9]}
    )
    s = pipeline.het_summary(het)
    assert s.at[0, "n_beyond_3sd"] == 1
    assert s.at[0, "ids_beyond_3sd"] == "odd"


def test_null_mask_excludes_window():
    glm = pd.DataFrame(
        {"CHROM": ["2", "2", "2", "6"], "POS": [60_000_000, 60_718_043, 62_000_000, 60_718_043]}
    )
    planted = pd.DataFrame({"chrom": ["2"], "pos": [60_718_043]})
    assert pipeline.null_mask(glm, planted, 1_000_000).tolist() == [False, False, True, True]
    assert pipeline.null_mask(glm, planted, 100).tolist() == [True, False, True, True]


def test_read_raw_counted_alleles(tmp_path):
    path, _, g = _raw(tmp_path)
    raw = pipeline.read_raw(path)
    assert raw.ids == ["rs1427407", "rs9399137", "rs41322148"]
    assert raw.counted == ["G", "C", "A"]
    np.testing.assert_array_equal(raw.genotypes, g)
    json.dumps(raw.ids)


def test_sample_ids_can_avoid_the_sentinels(pvar):
    ids = pipeline.sample_ids(pvar, 10, seed=3, exclude_windows=[("2", 60718043)], window_bp=1000)
    assert "rs1427407" not in ids and "rsA" not in ids  # rsA is 43 bp away
    assert set(ids) == {"rs9399137", "rs41322148", "rsB"}


def test_real_phenotype_refuses_an_empty_seed_or_undeclared_shift(tmp_path, pvar, loci):
    path, psam, _ = _raw(tmp_path)
    smap = pipeline.find_sentinels(pvar, loci)
    empty = tmp_path / "empty.yaml"
    empty.write_text("seed: null\n")
    with pytest.raises(ValueError, match="seed is empty"):
        pipeline.simulate_phenotype_from_raw(
            _cfg_with_shift(), pipeline.read_raw(path), smap, psam, empty
        )
    seeded = tmp_path / "s.yaml"
    seeded.write_text("seed: 1\n")
    with pytest.raises(ValueError, match="not pre-registered"):
        pipeline.simulate_phenotype_from_raw(CFG, pipeline.read_raw(path), smap, psam, seeded)


def test_step_e_seed_children_are_stable(tmp_path):
    p = tmp_path / "s.yaml"
    p.write_text("seed: 99\n")
    a = [pipeline.seed_child_int(p, n) for n in pipeline.STEP_E_SPAWN]
    assert a == [pipeline.seed_child_int(p, n) for n in pipeline.STEP_E_SPAWN]
    assert len(set(a)) == 3


# --- step e: power redraws at the sentinels -----------------------------------------


def _eigenvec(iids, seed=1):
    rng = np.random.default_rng(seed)
    pcs = rng.standard_normal((len(iids), 10))
    df = pd.DataFrame(pcs, columns=[f"PC{i}" for i in range(1, 11)])
    df.insert(0, "IID", iids)
    return df


def _cfg_for_redraws(n_redraws=1000):
    cfg = _cfg_with_shift()
    cfg["acceptance"]["real"] = {
        "power_calibration": {"max_abs_diff": 0.05, "n_redraws": n_redraws}
    }
    return cfg


def test_power_redraws_match_the_conditional_prediction(tmp_path, pvar, loci, seed_file):
    path, psam, _ = _raw(tmp_path, n=150, seed=3)
    raw = pipeline.read_raw(path)
    smap = pipeline.find_sentinels(pvar, loci)
    res = pipeline.power_redraws(
        _cfg_for_redraws(1000), raw, smap, psam, _eigenvec(raw.iid), seed_file
    )
    assert res["n_redraws"] == 1000 and res["seed_child"] == "power_reps"
    assert res["step_e_seed"] == 4242
    assert set(res["loci"]) == {"BCL11A", "HBS1L-MYB", "HBB cluster"}
    for gene, r in res["loci"].items():
        lo, hi = r["empirical_power_ci95"]
        assert lo <= r["empirical_power"] <= hi
        # empirical power over redraws agrees with the conditional prediction
        # (SE <= 0.016 at 1000 redraws)
        assert abs(r["empirical_power"] - r["predicted_power_recomputed"]) < 0.07, gene
        assert r["hits_pc_adjusted"] == round(r["empirical_power"] * 1000)
    again = pipeline.power_redraws(
        _cfg_for_redraws(1000), raw, smap, psam, _eigenvec(raw.iid), seed_file
    )
    assert again == res  # the step e seed fixes every redraw


def test_power_redraws_use_their_own_seed_child(tmp_path, pvar, loci, seed_file):
    path, psam, _ = _raw(tmp_path, n=150, seed=3)
    raw = pipeline.read_raw(path)
    smap = pipeline.find_sentinels(pvar, loci)
    other = tmp_path / "other_seed.yaml"
    other.write_text("seed: 4243\n")
    a = pipeline.power_redraws(
        _cfg_for_redraws(300), raw, smap, psam, _eigenvec(raw.iid), seed_file
    )
    b = pipeline.power_redraws(_cfg_for_redraws(300), raw, smap, psam, _eigenvec(raw.iid), other)
    assert a["loci"] != b["loci"]
    assert pipeline.seed_child_int(seed_file, "power_reps") != pipeline.seed_child_int(
        seed_file, "phenotype"
    )


def test_power_redraws_refuse_without_a_count_or_a_declared_shift(tmp_path, pvar, loci, seed_file):
    path, psam, _ = _raw(tmp_path, n=60)
    raw = pipeline.read_raw(path)
    smap = pipeline.find_sentinels(pvar, loci)
    no_count = _cfg_with_shift()
    no_count["acceptance"]["real"] = None
    with pytest.raises(ValueError, match="n_redraws"):
        pipeline.power_redraws(no_count, raw, smap, psam, _eigenvec(raw.iid), seed_file)
    import copy

    no_shift = copy.deepcopy(_cfg_for_redraws(10))
    no_shift["real_data_simulation"]["population_shift"] = None
    with pytest.raises(ValueError, match="population_shift"):
        pipeline.power_redraws(no_shift, raw, smap, psam, _eigenvec(raw.iid), seed_file)


# --- step e integrity check ------------------------------------------------------------


def _declared_cfg():
    import copy

    cfg = copy.deepcopy(CFG)
    cfg["acceptance"]["real"] = {
        "predicted_power": {"BCL11A": 0.9896, "HBS1L-MYB": 0.1694},
        "predicted_lambda_gc": {
            "unadjusted": {"median": 1.2001, "q01": 1.1248, "q99": 1.2874},
            "pc_adjusted": {"median": 1.0094},
        },
    }
    return cfg


def _recomputed(**over):
    pred = {
        "conditional_power": {
            "BCL11A": {"pc_adjusted": {"power": 0.989576}},
            "HBS1L-MYB": {"pc_adjusted": {"power": 0.169391}},
        },
        "predicted_lambda_gc": {
            "unadjusted": {"median": 1.1999},
            "pc_adjusted": {"median": 1.0045},
        },
        "null_distribution": {
            "unadjusted": {"lambda_gc": {"median": 1.20012, "q01": 1.124814, "q99": 1.287446}},
            "pc_adjusted": {"lambda_gc": {"median": 1.00941}},
        },
    }
    pred.update(over)
    return pred


def test_integrity_check_passes_when_the_recomputed_values_match():
    res = pipeline.check_declared(_declared_cfg(), _recomputed())
    assert res["passed"], res
    assert len(res["checks"]) == 6
    assert "not an acceptance criterion" in res["what"]
    # the newer null_distribution is used, not the older summary (1.1999)
    row = [c for c in res["checks"] if c["name"] == "predicted lambda_GC unadjusted median"][0]
    assert row["recomputed"] == 1.20012


def test_integrity_check_fails_on_any_difference_or_missing_value():
    pred = _recomputed()
    pred["conditional_power"]["HBS1L-MYB"]["pc_adjusted"]["power"] = 0.1699
    res = pipeline.check_declared(_declared_cfg(), pred)
    assert not res["passed"]
    assert [c["name"] for c in res["checks"] if not c["ok"]] == [
        "conditional power HBS1L-MYB (PC-adjusted)"
    ]
    pred = _recomputed(null_distribution=None)  # no q01/q99 to compare
    res = pipeline.check_declared(_declared_cfg(), pred)
    bad = sorted(c["name"] for c in res["checks"] if not c["ok"])
    assert "predicted lambda_GC unadjusted q01" in bad and not res["passed"]


def test_integrity_check_tolerance_is_the_4_decimal_rounding():
    assert pipeline.INTEGRITY_TOL == 1e-4
    pred = _recomputed()
    pred["conditional_power"]["BCL11A"]["pc_adjusted"]["power"] = 0.9896 + 1.5e-4
    assert not pipeline.check_declared(_declared_cfg(), pred)["passed"]
