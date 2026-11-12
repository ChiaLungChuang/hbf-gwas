"""Acceptance gate: compare a results JSON with the pre-registered criteria.

The criteria live in ``config/analysis.yaml`` under ``acceptance`` and were
declared before any experiment ran. This module never changes them; it only
reads them and reports PASS / FAIL for each gated criterion and REPORTED for
results that are shown but not gated.

Exit codes of ``python -m hbfgwas.gate``: 0 all criteria met, 1 at least one
criterion failed, 2 the requested experiment has no pre-registered criteria.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass

from hbfgwas import redact
from hbfgwas.config import load_config


class NotPreregisteredError(RuntimeError):
    """The experiment has no acceptance criteria in the config."""


@dataclass
class Check:
    name: str
    status: str  # PASS, FAIL, REPORTED (shown, not gated) or NOT TESTABLE
    observed: object
    criterion: str
    detail: str = ""


def _interval(name: str, value: float, rule: dict) -> Check:
    ok = True
    parts = []
    if "greater_than" in rule:
        ok &= value > rule["greater_than"]
        parts.append(f"> {rule['greater_than']}")
    if "at_least" in rule:
        ok &= value >= rule["at_least"]
        parts.append(f">= {rule['at_least']}")
    if "at_most" in rule:
        ok &= value <= rule["at_most"]
        parts.append(f"<= {rule['at_most']}")
    if not parts:
        raise ValueError(f"criterion for {name} has no bound")
    return Check(name, "PASS" if ok else "FAIL", value, " and ".join(parts))


def _locus_checks(config: dict, loci: list[dict]) -> list[Check]:
    rec = config["acceptance"]["locus_recovery"]
    threshold = rec["required_if_predicted_power_at_least"]
    label = rec["not_required_label"]
    alpha = config["association"]["genome_wide_alpha"]
    checks = []
    for loc in loci:
        recovered = loc["p_adjusted"] < alpha
        if bool(loc["recovered"]) != recovered:
            raise ValueError(f"locus {loc['locus']}: 'recovered' disagrees with p < {alpha}")
        power = loc["predicted_power"]
        name = f"locus {loc['locus']} recovery"
        observed = f"p = {loc['p_adjusted']:.3g}"
        crit = f"p < {alpha} if predicted power >= {threshold}"
        hit = "hit" if recovered else "miss"
        if power >= threshold:
            status = "PASS" if recovered else "FAIL"
            detail = f"REQUIRED (predicted power {power:.4f}); single scan: {hit}"
        else:
            status = "REPORTED"
            detail = f"{label} (predicted power {power:.4f}); single scan: {hit}"
        checks.append(Check(name, status, observed, crit, detail))
    return checks


def evaluate(config: dict, results: dict, experiment: str = "simulated") -> list[Check]:
    """Evaluate every pre-registered criterion for ``experiment``."""
    criteria = config["acceptance"].get(experiment)
    if not criteria:
        raise NotPreregisteredError(
            f"no acceptance criteria for '{experiment}' in the config; "
            "real-data criteria are written in Phase 2 step d and tagged prereg-v1"
        )
    t1 = criteria["adjusted_type1_error"]
    if results["adjusted"]["type1_alpha"] != t1["alpha"]:
        raise ValueError("results used a different alpha for type I error than the config")
    checks = [
        _interval(
            "unadjusted lambda_GC",
            results["unadjusted"]["lambda_gc"],
            criteria["unadjusted_lambda_gc"],
        ),
        _interval(
            "PC-adjusted lambda_GC",
            results["adjusted"]["lambda_gc"],
            criteria["adjusted_lambda_gc"],
        ),
        _interval(
            "PC-adjusted type I error",
            results["adjusted"]["type1_error"],
            {k: v for k, v in t1.items() if k != "alpha"},
        ),
    ]
    checks[-1].criterion = f"at alpha {t1['alpha']}: " + checks[-1].criterion
    checks += _locus_checks(config, results.get("loci", []))
    label = (config.get("sentinels") or {}).get(
        "not_testable_label", "not testable in this callset"
    )
    for name in results.get("not_testable", []):
        checks.append(
            Check(
                f"locus {name} recovery",
                "NOT TESTABLE",
                "-",
                "sentinel absent from the callset or failed QC; no proxy",
                label,
            )
        )
    calib = criteria.get("power_calibration")
    if calib:
        checks += _calibration_checks(calib, results.get("loci", []))
    if "predictions_match_declared" in results:
        ok = bool(results["predictions_match_declared"])
        checks.append(
            Check(
                "declared predictions reproduced",
                "PASS" if ok else "FAIL",
                ok,
                "predicted power recomputed by the run equals the value declared in the config",
            )
        )
    return checks


def _calibration_checks(calib: dict, loci: list[dict]) -> list[Check]:
    tol = calib["max_abs_diff"]
    label = f" ({calib['status']})" if calib.get("status") else ""
    checks = []
    for loc in loci:
        e, p = loc["empirical_power"], loc["predicted_power"]
        lo, hi = loc.get("empirical_power_ci95", [float("nan")] * 2)
        ok = abs(e - p) <= tol
        checks.append(
            Check(
                f"locus {loc['locus']} power calibration",
                "PASS" if ok else "FAIL",
                f"{e:.3f} vs {p:.3f}",
                f"|empirical - predicted| <= {tol}{label}",
                f"empirical 95% CI [{lo:.3f}, {hi:.3f}]; difference {e - p:+.3f}",
            )
        )
    return checks


def any_failed(checks: list[Check]) -> bool:
    return any(c.status == "FAIL" for c in checks)


def format_table(checks: list[Check]) -> str:
    rows = [("criterion", "status", "observed", "rule", "detail")]
    for c in checks:
        obs = f"{c.observed:.4f}" if isinstance(c.observed, float) else str(c.observed)
        rows.append((c.name, c.status, obs, c.criterion, c.detail))
    widths = [max(len(r[i]) for r in rows) for i in range(5)]
    return "\n".join("  ".join(v.ljust(w) for v, w in zip(r, widths, strict=True)) for r in rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m hbfgwas.gate", description=__doc__)
    parser.add_argument("--config", default=None, help="default: config/analysis.yaml")
    parser.add_argument("--results", required=True, help="results JSON from an experiment")
    parser.add_argument("--experiment", default="simulated", choices=["simulated", "real"])
    parser.add_argument("--out", default=None, help="write the checks as JSON here")
    args = parser.parse_args(argv)

    config = load_config(args.config)
    with open(args.results) as fh:
        results = json.load(fh)
    try:
        checks = evaluate(config, results, args.experiment)
    except NotPreregisteredError as exc:
        print(f"GATE NOT RUN: {exc}")
        return 2
    print(format_table(checks))
    failed = any_failed(checks)
    print("\nGATE " + ("FAILED" if failed else "PASSED"))
    if args.out:
        doc = {"passed": not failed, "checks": [asdict(c) for c in checks]}
        redact.write_record(args.out, json.dumps(doc, indent=2))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
