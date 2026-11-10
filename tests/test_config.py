"""The pre-registration file must say exactly what the brief says.

These literal values come from docs/BRIEF.md. If one of these tests fails,
an acceptance criterion or analysis setting was edited: that is not allowed
after a result has been seen.
"""

import numpy as np
import pytest

from hbfgwas import power, simulate
from hbfgwas.config import load_config

CFG = load_config()


def test_qc_thresholds_from_the_brief():
    qc = CFG["qc"]
    assert qc["geno"] == 0.02
    assert qc["mind"] == 0.02
    assert qc["maf"] == 0.01
    assert qc["hwe_p"] == 1e-6
    assert qc["hwe_mode"] == "midp"
    assert qc["expected_afr_samples"] == {"min": 640, "max": 665}


def test_pca_settings_from_the_brief():
    assert CFG["pca"]["n_pcs"] == 10
    assert CFG["pca"]["method"] == "exact"
    assert CFG["pca"]["king_cutoff"] == 0.0884


def test_significance_threshold():
    assert CFG["association"]["genome_wide_alpha"] == 5e-8
    assert CFG["association"]["type1_alpha"] == 0.05


def test_simulated_acceptance_criteria_from_the_brief():
    sim = CFG["acceptance"]["simulated"]
    assert sim["unadjusted_lambda_gc"] == {"greater_than": 1.2}
    assert sim["adjusted_lambda_gc"] == {"at_least": 0.97, "at_most": 1.05}
    assert sim["adjusted_type1_error"] == {"alpha": 0.05, "at_least": 0.04, "at_most": 0.06}


def test_locus_recovery_criterion_does_not_demand_all_loci():
    rec = CFG["acceptance"]["locus_recovery"]
    assert rec["required_if_predicted_power_at_least"] == 0.8
    assert rec["not_required_label"] == "underpowered, not required"
    assert not any("all" in str(k).lower() for k in rec)


def test_real_data_criteria_not_written_before_phase2():
    assert CFG["acceptance"]["real"] is None


def test_seeds_are_distinct_integers():
    seeds = CFG["seeds"]
    assert all(isinstance(v, int) for v in seeds.values())
    assert len(set(seeds.values())) == len(seeds)


def test_declared_predictions_match_the_code():
    exp = CFG["simulated_experiment"]
    n = sum(exp["n_per_subpop"])
    pred = exp["predicted"]
    lam = simulate.expected_null_chi2(
        n=n,
        fst=exp["fst"],
        shift_var_fraction=exp["shift_var_fraction"],
        n_pops=len(exp["n_per_subpop"]),
    )
    assert pred["unadjusted_lambda_gc"] == pytest.approx(lam, abs=5e-4)
    alpha = CFG["association"]["genome_wide_alpha"]
    for planted in exp["planted"]:
        r2 = power.r2_from_beta(planted["beta_sd"], planted["ancestral_freq"])
        assert planted["r2"] == pytest.approx(r2, abs=5e-6)
        p = power.power(n=n, alpha=alpha, r2=r2, covariate_r2=exp["shift_var_fraction"])
        assert pred["power_pc_adjusted"][planted["locus"]] == pytest.approx(p, abs=5e-4)


def test_planted_loci_match_the_locus_table():
    from hbfgwas.loci import LOCI

    names = [p["locus"] for p in CFG["simulated_experiment"]["planted"]]
    assert names == [loc.gene.value for loc in LOCI]
    assert np.all([0 < p["ancestral_freq"] < 1 for p in CFG["simulated_experiment"]["planted"]])
