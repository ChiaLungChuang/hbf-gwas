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
    assert qc["hwe_keep_fewhet"] is True  # review item D (user, D-029)
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
    # The brief's rule was 0.8 ("underpowered, not required"). Amended by the user
    # after the failed sim-v2 result (issue #14, D-024): REQUIRED only at >= 0.99,
    # every other locus reported, not gated.
    rec = CFG["acceptance"]["locus_recovery"]
    assert rec["required_if_predicted_power_at_least"] == 0.99
    assert rec["not_required_label"] == "reported, not gated"
    assert not any("all" in str(k).lower() for k in rec)


def test_real_data_criteria_follow_the_v3_rule():
    # Phase 2 step d (proposed; final when merged and tagged prereg-v1), amended
    # by the user before any real-data phenotype existed (D-048 to D-050)
    real = CFG["acceptance"]["real"]
    sim = CFG["acceptance"]["simulated"]
    for key in ("adjusted_lambda_gc", "adjusted_type1_error", "power_calibration"):
        assert {k: v for k, v in real[key].items() if k != "n_redraws"} == sim[key], key
    assert real["power_calibration"]["n_redraws"] == 2000
    assert real["predicted_power"] == {"BCL11A": 0.9896, "HBS1L-MYB": 0.1694}
    threshold = CFG["acceptance"]["locus_recovery"]["required_if_predicted_power_at_least"]
    assert real["required_loci"] == [
        g for g, p in real["predicted_power"].items() if p >= threshold
    ]  # flag 2: no change, BCL11A (0.9896) is reported, not gated
    assert real["required_loci"] == []
    assert real["not_testable"] == ["HBB cluster"]


def test_unadjusted_lambda_criterion_is_the_pre_stated_1_99_interval():
    # Flag 1 (D-048): "> 1.2" replaced by the 1st-99th percentile prediction
    # interval of the observed value, rounded outward to 4 decimals.
    import json
    import math
    from pathlib import Path

    real = CFG["acceptance"]["real"]
    run = Path(__file__).resolve().parents[1] / "runs/2026-11-12-stepd-amend"
    pred = json.loads((run / "prereg_predictions.json").read_text())
    q = pred["null_distribution"]["unadjusted"]["lambda_gc"]
    crit = real["unadjusted_lambda_gc"]
    assert set(crit) == {"at_least", "at_most"}
    assert crit["at_least"] == math.floor(q["q01"] * 1e4) / 1e4
    assert crit["at_most"] == math.ceil(q["q99"] * 1e4) / 1e4
    assert "greater_than" not in crit


def test_every_gated_criterion_has_a_pass_probability_and_low_ones_are_flagged():
    import json
    from pathlib import Path

    real = CFG["acceptance"]["real"]
    pp = real["pass_probability"]
    run = Path(__file__).resolve().parents[1] / "runs/2026-11-12-stepd-amend"
    rec = json.loads((run / "prereg_predictions.json").read_text())["pass_probability"]
    per = rec["per_criterion"]
    assert pp["unadjusted_lambda_gc"] == round(per["unadjusted_lambda_gc_in_q01_q99"], 4)
    assert pp["adjusted_lambda_gc"] == round(per["adjusted_lambda_gc"], 4)
    assert pp["adjusted_type1_error"] == round(per["adjusted_type1_error"], 4)
    for gene in real["predicted_power"]:
        got = rec["power_calibration"][gene]["if_prediction_exact"]
        assert pp[f"power_calibration_{gene}"] == round(got, 4)
    assert pp["all_gated_criteria"] == round(rec["all_gated_criteria"], 4)
    low = sorted(
        k for k, v in pp.items() if k not in ("all_gated_criteria", "below_0_99") and v < 0.99
    )
    assert sorted(pp["below_0_99"]) == low == ["unadjusted_lambda_gc"]


def test_real_predictions_match_the_committed_step_d_records():
    import json
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    real = CFG["acceptance"]["real"]
    for run in ("runs/2026-11-12-stepd", "runs/2026-11-12-stepd-amend"):
        pred = json.loads((root / run / "prereg_predictions.json").read_text())
        for gene, p in real["predicted_power"].items():
            got = pred["conditional_power"][gene]["pc_adjusted"]["power"]
            assert p == pytest.approx(got, abs=5e-5), (run, gene)
        assert pred["population_shift"] == CFG["real_data_simulation"]["population_shift"]
    amend = json.loads((root / "runs/2026-11-12-stepd-amend/prereg_predictions.json").read_text())
    for side in ("unadjusted", "pc_adjusted"):
        q = amend["null_distribution"][side]["lambda_gc"]
        for k in ("median", "q01", "q99"):
            assert real["predicted_lambda_gc"][side][k] == pytest.approx(q[k], abs=5e-5), side


def test_population_shift_is_centred_with_variance_02():
    import numpy as np

    shift = CFG["real_data_simulation"]["population_shift"]
    counts = {"ACB": 95, "ASW": 55, "ESN": 99, "GWD": 113, "LWK": 97, "MSL": 85, "YRI": 108}
    v = np.repeat([shift[k] for k in counts], list(counts.values()))
    assert abs(v.mean()) < 1e-6 and v.var() == pytest.approx(0.2, abs=1e-6)


def test_step_e_seed_slot_is_empty_until_after_the_tag():
    from pathlib import Path

    import yaml

    seed = yaml.safe_load(
        (Path(__file__).resolve().parents[1] / "config/step_e_seed.yaml").read_text()
    )
    assert seed == {"seed": None}


def test_seeds_are_distinct_integers():
    # The current experiment's seed may be null only until it is pre-registered.
    key = CFG["simulated_experiment"]["seed_key"]
    seeds = {k: v for k, v in CFG["seeds"].items() if not (k == key and v is None)}
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
        assert pred["unconditional_power_pc_adjusted"][planted["locus"]] == pytest.approx(
            p, abs=5e-4
        )


def test_predictions_cannot_be_declared_before_the_seed():
    exp = CFG["simulated_experiment"]
    if CFG["seeds"][exp["seed_key"]] is None:
        assert exp["predicted"]["conditional_power_pc_adjusted"] is None


def test_declared_conditional_predictions_match_the_code():
    from hbfgwas import experiment

    declared = CFG["simulated_experiment"]["predicted"]["conditional_power_pc_adjusted"]
    if declared is None:
        pytest.skip("sim-v3 predictions are declared together with the sim-v3 seed")
    cond = experiment.conditional_predictions(CFG)["pc_adjusted"]
    assert set(declared) == set(cond)
    for locus, value in declared.items():
        assert value == pytest.approx(cond[locus]["power"], abs=5e-4), locus


def test_power_calibration_is_gated_at_005():
    # issue #10 proposed it; issue #14 (sim-v3) makes it a gated criterion
    assert CFG["acceptance"]["simulated"]["power_calibration"] == {"max_abs_diff": 0.05}


def test_criteria_version_is_sim_v3():
    assert CFG["prereg"]["version"] == "sim-v3"
    assert CFG["simulated_experiment"]["version"] == "sim-v3"
    assert CFG["simulated_experiment"]["seed_key"] == "sim_v3"


def test_planted_loci_match_the_locus_table():
    from hbfgwas.loci import LOCI

    names = [p["locus"] for p in CFG["simulated_experiment"]["planted"]]
    assert names == [loc.gene.value for loc in LOCI]
    assert np.all([0 < p["ancestral_freq"] < 1 for p in CFG["simulated_experiment"]["planted"]])


def test_planted_effects_agree_between_sections():
    effects = {e["locus"]: e for e in CFG["planted_effects"]}
    for planted in CFG["simulated_experiment"]["planted"]:
        assert planted["beta_sd"] == effects[planted["locus"]]["beta_sd"]


def test_planted_effect_alleles_match_the_locus_table():
    from hbfgwas.loci import LOCI

    table = {loc.gene.value: loc.effect_allele.value for loc in LOCI}
    for e in CFG["planted_effects"]:
        assert e["effect_allele"] == table[e["locus"]]


def test_sentinel_rule_has_no_proxies():
    # review item C (user, D-028)
    from hbfgwas.pipeline import NOT_TESTABLE

    assert CFG["sentinels"] == {"not_testable_label": NOT_TESTABLE, "proxies": "none"}


def test_prediction_settings():
    # review item H: Phase 2 step d predictions, computed before any phenotype
    assert CFG["prediction"] == {
        "n_null_variants": 20000,
        "power_shift_draws": 2000,
        "lambda_shift_draws": 200,
        # step d amendment (D-049), before any real-data phenotype
        "null_draws": 2000,
        "n_effective_variants": 215522,
        "power_model_error": 0.013,
    }
    assert isinstance(CFG["seeds"]["real_prediction"], int)


def test_plink2_pin_matches_ci():
    # review item G: CI downloads exactly the Linux build recorded in the config
    from pathlib import Path

    # CI and the weekly rebuild both install plink2 through this one action
    action = Path(__file__).resolve().parents[1] / ".github/actions/install-plink2/action.yml"
    ci = action.read_text()
    pin = CFG["plink2"]
    assert pin["version"] == "v2.0.0-b.1.4" and pin["build_date"] == "7 Nov 2026"
    linux = pin["builds"]["linux_x86_64"]
    assert f"PLINK2_URL: {linux['url']}" in ci
    assert f"PLINK2_SHA256: {linux['zip_sha256']}" in ci
    assert "20261107" in pin["builds"]["mac_arm64"]["url"]


def test_bcl11a_beta_is_the_derived_sd_value_not_wonkams_trait_units():
    # user note (2026-11-12): Wonkam 2025's beta = 0.26 is in trait units; keep 0.456435 SD
    import math

    b = {e["locus"]: e["beta_sd"] for e in CFG["planted_effects"]}["BCL11A"]
    assert b == 0.456435
    assert b == pytest.approx(math.sqrt(0.076 / (2 * 0.24 * 0.76)), abs=1e-6)
