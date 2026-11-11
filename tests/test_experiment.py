"""Tests for the stratification experiment runner, on a SMALL configuration.

These use their own seeds and sizes; the pre-registered experiment
(config/analysis.yaml, simulated_experiment) is run by
scripts/run_stratification.py, not by the test suite.
"""

import copy
import json

import matplotlib

matplotlib.use("Agg")

import pytest

from hbfgwas import experiment, gate
from hbfgwas.config import load_config


@pytest.fixture(scope="module")
def small_config():
    cfg = copy.deepcopy(load_config())
    exp = cfg["simulated_experiment"]
    exp["n_per_subpop"] = [60] * 7
    exp["n_null_variants"] = 3000
    exp["power_reps"] = 40
    exp["prediction_shift_draws"] = 50
    exp["predicted"]["conditional_power_pc_adjusted"] = {}
    exp["seed_key"] = "test_seed"
    cfg["seeds"] = {"test_seed": 11, "real_phenotype": 14, "real_concordance_sample": 15}
    return cfg


@pytest.fixture(scope="module")
def small_results(small_config):
    return experiment.run_stratification(small_config)


def test_result_schema(small_results):
    r = small_results
    for side in ("unadjusted", "adjusted"):
        assert {"lambda_gc", "type1_error", "type1_alpha"} <= set(r[side])
    assert r["adjusted"]["n_covariates"] == 10
    assert len(r["loci"]) == 3
    for loc in r["loci"]:
        assert {
            "locus",
            "predicted_power",
            "p_adjusted",
            "p_unadjusted",
            "recovered",
            "empirical_power",
            "empirical_power_ci95",
        } <= set(loc)
    assert r["design"]["n"] == 420
    assert r["design"]["n_null_variants"] == 3000


def test_structure_is_detected_and_removed(small_results):
    assert small_results["unadjusted"]["lambda_gc"] > small_results["adjusted"]["lambda_gc"]
    assert small_results["fst_estimate"] > 0


def test_gate_accepts_the_result_format(small_config, small_results):
    checks = gate.evaluate(small_config, small_results, experiment="simulated")
    assert {c.name for c in checks} >= {"unadjusted lambda_GC", "PC-adjusted lambda_GC"}


def test_deterministic(small_config, small_results):
    again = experiment.run_stratification(small_config)
    assert again["adjusted"]["lambda_gc"] == small_results["adjusted"]["lambda_gc"]
    assert [x["p_adjusted"] for x in again["loci"]] == [
        x["p_adjusted"] for x in small_results["loci"]
    ]


def test_predictions_are_recomputed_from_the_genotypes(small_config, small_results):
    # sizes changed, so predicted power must come from the code, not the YAML numbers
    for loc in small_results["loci"]:
        assert 0 < loc["predicted_power"] <= 1
    assert small_results["predicted"]["unadjusted_lambda_gc"] > 1
    assert small_results["predictions_match_declared"] is False  # nothing declared here
    cond = experiment.conditional_predictions(small_config)["pc_adjusted"]
    assert [cond[x["locus"]]["power"] for x in small_results["loci"]] == [
        x["predicted_power"] for x in small_results["loci"]
    ]


def test_planted_effects_are_the_per_allele_betas(small_config, small_results):
    betas = [p["beta_sd"] for p in small_config["simulated_experiment"]["planted"]]
    assert [x["beta_planted"] for x in small_results["loci"]] == betas


def test_writes_outputs(tmp_path, small_config, small_results):
    paths = experiment.write_outputs(small_results, small_config, tmp_path)
    assert paths["figure"].stat().st_size > 10_000
    summary = json.loads(paths["summary"].read_text())
    assert summary["unadjusted"]["lambda_gc"] == small_results["unadjusted"]["lambda_gc"]
    full = json.loads(paths["results"].read_text())
    assert "provenance" in full


def test_summary_comparison(small_results):
    s = experiment.summarize(small_results)
    assert experiment.compare_summaries(s, s) == []
    other = copy.deepcopy(s)
    other["adjusted"]["lambda_gc"] += 0.01
    assert experiment.compare_summaries(s, other)


def test_unregistered_seed_runs_nothing(small_config):
    cfg = copy.deepcopy(small_config)
    cfg["seeds"]["test_seed"] = None
    with pytest.raises(experiment.NotPreregisteredError):
        experiment.run_stratification(cfg)


def test_seeds_are_spawned_from_one_seed(small_config, small_results):
    seeds = experiment.experiment_seeds(small_config)
    assert seeds["seed"] == 11
    assert small_results["design"]["seeds"]["spawned_in_order"] == [
        "genotypes",
        "phenotype",
        "power_reps",
        "prediction",
    ]
    other = copy.deepcopy(small_config)
    other["seeds"]["test_seed"] = 12
    a = experiment.simulate_genotypes_and_pcs(small_config)["plant"].genotypes
    b = experiment.simulate_genotypes_and_pcs(other)["plant"].genotypes
    assert (a != b).any()
