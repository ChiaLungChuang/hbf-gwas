"""Tests for the acceptance gate (logic only, on hand-made results)."""

import copy
import json
import subprocess
import sys

import pytest

from hbfgwas import gate

CONFIG = {
    "association": {"genome_wide_alpha": 5e-8, "type1_alpha": 0.05},
    "acceptance": {
        "simulated": {
            "unadjusted_lambda_gc": {"greater_than": 1.2},
            "adjusted_lambda_gc": {"at_least": 0.97, "at_most": 1.05},
            "adjusted_type1_error": {"alpha": 0.05, "at_least": 0.04, "at_most": 0.06},
        },
        "locus_recovery": {
            "required_if_predicted_power_at_least": 0.8,
            "not_required_label": "underpowered, not required",
        },
        "real": None,
    },
}

GOOD = {
    "experiment": "simulated",
    "unadjusted": {"lambda_gc": 1.4, "type1_error": 0.12, "type1_alpha": 0.05},
    "adjusted": {"lambda_gc": 1.01, "type1_error": 0.051, "type1_alpha": 0.05},
    "loci": [
        {"locus": "A", "predicted_power": 0.95, "p_adjusted": 1e-12, "recovered": True},
        {"locus": "B", "predicted_power": 0.40, "p_adjusted": 1e-3, "recovered": False},
        {"locus": "C", "predicted_power": 0.10, "p_adjusted": 1e-9, "recovered": True},
    ],
}


def _status(checks, name):
    return next(c.status for c in checks if c.name == name)


def test_all_pass():
    checks = gate.evaluate(CONFIG, GOOD, experiment="simulated")
    assert not gate.any_failed(checks)
    assert _status(checks, "locus A recovery") == "PASS"
    assert _status(checks, "locus B recovery") == "NOT REQUIRED"
    assert _status(checks, "locus C recovery") == "NOT REQUIRED"


@pytest.mark.parametrize(
    "path,value,expected",
    [
        (("unadjusted", "lambda_gc"), 1.2, "FAIL"),  # strictly greater than 1.2
        (("unadjusted", "lambda_gc"), 1.2001, "PASS"),
        (("adjusted", "lambda_gc"), 0.9699, "FAIL"),
        (("adjusted", "lambda_gc"), 0.97, "PASS"),  # closed interval
        (("adjusted", "lambda_gc"), 1.05, "PASS"),
        (("adjusted", "lambda_gc"), 1.0501, "FAIL"),
        (("adjusted", "type1_error"), 0.0399, "FAIL"),
        (("adjusted", "type1_error"), 0.04, "PASS"),
        (("adjusted", "type1_error"), 0.06, "PASS"),
        (("adjusted", "type1_error"), 0.0601, "FAIL"),
    ],
)
def test_boundaries(path, value, expected):
    res = copy.deepcopy(GOOD)
    res[path[0]][path[1]] = value
    names = {
        ("unadjusted", "lambda_gc"): "unadjusted lambda_GC",
        ("adjusted", "lambda_gc"): "PC-adjusted lambda_GC",
        ("adjusted", "type1_error"): "PC-adjusted type I error",
    }
    checks = gate.evaluate(CONFIG, res, experiment="simulated")
    assert _status(checks, names[path]) == expected


def test_required_locus_not_recovered_fails():
    res = copy.deepcopy(GOOD)
    res["loci"][0].update(p_adjusted=1e-6, recovered=False)
    checks = gate.evaluate(CONFIG, res, experiment="simulated")
    assert _status(checks, "locus A recovery") == "FAIL"
    assert gate.any_failed(checks)


def test_power_exactly_at_threshold_is_required():
    res = copy.deepcopy(GOOD)
    res["loci"][1].update(predicted_power=0.8)
    checks = gate.evaluate(CONFIG, res, experiment="simulated")
    assert _status(checks, "locus B recovery") == "FAIL"


def test_not_required_detail_states_power():
    checks = gate.evaluate(CONFIG, GOOD, experiment="simulated")
    detail = next(c.detail for c in checks if c.name == "locus B recovery")
    assert "underpowered, not required" in detail
    assert "0.40" in detail


def test_inconsistent_recovered_flag_is_an_error():
    res = copy.deepcopy(GOOD)
    res["loci"][0]["recovered"] = False  # p says recovered
    with pytest.raises(ValueError):
        gate.evaluate(CONFIG, res, experiment="simulated")


def test_alpha_mismatch_is_an_error():
    res = copy.deepcopy(GOOD)
    res["adjusted"]["type1_alpha"] = 0.01
    with pytest.raises(ValueError):
        gate.evaluate(CONFIG, res, experiment="simulated")


def test_real_without_preregistered_criteria_refuses():
    with pytest.raises(gate.NotPreregisteredError):
        gate.evaluate(CONFIG, GOOD, experiment="real")


def test_cli_exit_codes(tmp_path):
    cfg = tmp_path / "c.yaml"
    import yaml

    cfg.write_text(yaml.safe_dump(CONFIG))
    good = tmp_path / "good.json"
    good.write_text(json.dumps(GOOD))
    bad_res = copy.deepcopy(GOOD)
    bad_res["adjusted"]["lambda_gc"] = 1.2
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(bad_res))

    def run(results, experiment="simulated"):
        return subprocess.run(
            [
                sys.executable,
                "-m",
                "hbfgwas.gate",
                "--config",
                str(cfg),
                "--results",
                str(results),
                "--experiment",
                experiment,
            ],
            capture_output=True,
            text=True,
        )

    ok = run(good)
    assert ok.returncode == 0, ok.stdout + ok.stderr
    assert "PASS" in ok.stdout
    fail = run(bad)
    assert fail.returncode == 1
    assert "FAIL" in fail.stdout
    assert run(good, "real").returncode == 2


def _with_calibration(tol=0.05):
    cfg = copy.deepcopy(CONFIG)
    cfg["acceptance"]["simulated"]["power_calibration"] = {
        "max_abs_diff": tol,
        "status": "proposed",
    }
    res = copy.deepcopy(GOOD)
    for loc, emp in zip(res["loci"], (0.96, 0.44, 0.12), strict=True):
        loc["empirical_power"] = emp
        loc["empirical_power_ci95"] = [emp - 0.02, emp + 0.02]
    return cfg, res


def test_power_calibration_passes_within_tolerance():
    cfg, res = _with_calibration()
    checks = gate.evaluate(cfg, res, experiment="simulated")
    assert _status(checks, "locus B power calibration") == "PASS"  # 0.44 vs 0.40
    assert not gate.any_failed(checks)


def test_power_calibration_fails_outside_tolerance():
    cfg, res = _with_calibration()
    res["loci"][1]["empirical_power"] = 0.4501  # 0.0501 above 0.40
    checks = gate.evaluate(cfg, res, experiment="simulated")
    assert _status(checks, "locus B power calibration") == "FAIL"
    res["loci"][1]["empirical_power"] = 0.45  # exactly at the tolerance passes
    checks = gate.evaluate(cfg, res, experiment="simulated")
    assert _status(checks, "locus B power calibration") == "PASS"


def test_declared_predictions_must_be_reproduced():
    res = copy.deepcopy(GOOD)
    res["predictions_match_declared"] = False
    checks = gate.evaluate(CONFIG, res, experiment="simulated")
    assert _status(checks, "declared predictions reproduced") == "FAIL"
