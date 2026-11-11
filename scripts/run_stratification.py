#!/usr/bin/env python3
"""Run the pre-registered stratification experiment on simulated genotypes.

Every setting comes from config/analysis.yaml (simulated_experiment, seeds,
association). The phenotype is simulated. Writes, under --out-dir:
results.json (full, with provenance), sim_stratification_summary.json and
sim_stratification_qq.png; copies the last two to --figures-dir.

    python scripts/run_stratification.py                      # the recorded run
    python scripts/run_stratification.py --out-dir /tmp/sim --figures-dir "" \
        --compare-to docs/figures/sim_stratification_summary.json   # CI check

Exit codes: 0 ok; 3 the summary differs from --compare-to beyond tolerance.
The acceptance gate is separate: python -m hbfgwas.gate --results <results.json>
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from hbfgwas import experiment  # noqa: E402
from hbfgwas.config import DEFAULT_CONFIG, load_config  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    today = time.strftime("%Y-%m-%d")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--out-dir", default=f"runs/{today}-sim-stratification")
    parser.add_argument("--figures-dir", default="docs/figures", help='"" to skip copying')
    parser.add_argument("--compare-to", default=None, help="summary JSON to reproduce")
    args = parser.parse_args(argv)

    config = load_config(args.config)
    results = experiment.run_stratification(config)
    paths = experiment.write_outputs(results, config, Path(args.out_dir), Path(args.config))

    s = experiment.summarize(results)
    print(
        f"n = {s['design']['n']}, null variants = {s['design']['n_null_variants']}, "
        f"Fst = {s['design']['fst']} (estimated {s['fst_estimate']:.5f})"
    )
    print(f"predicted unadjusted lambda_GC = {s['predicted_unadjusted_lambda_gc']:.4f}")
    for side in ("unadjusted", "adjusted"):
        r = s[side]
        print(
            f"{side:>10}: lambda_GC = {r['lambda_gc']:.4f}, "
            f"type I error at {r['type1_alpha']} = {r['type1_error']:.4f}"
        )
    for loc in s["loci"]:
        lo, hi = loc["empirical_power_ci95"]
        print(
            f"{loc['locus']:>12}: predicted power {loc['predicted_power']:.4f}, "
            f"empirical {loc['empirical_power']:.3f} [{lo:.3f}, {hi:.3f}], "
            f"primary p_adjusted = {loc['p_adjusted']:.3g} "
            f"({'recovered' if loc['recovered'] else 'not recovered'})"
        )
    print(f"predictions match declared config values: {results['predictions_match_declared']}")
    print(f"runtime {results['runtime_s']:.1f} s; wrote {paths['results']}")

    if args.figures_dir:
        fig_dir = Path(args.figures_dir)
        fig_dir.mkdir(parents=True, exist_ok=True)
        for key in ("summary", "figure"):
            shutil.copy2(paths[key], fig_dir / paths[key].name)
        print(f"copied summary and figure to {fig_dir}/")

    if args.compare_to:
        reference = json.loads(Path(args.compare_to).read_text())
        diffs = experiment.compare_summaries(reference, s)
        if diffs:
            print("SUMMARY DIFFERS from " + args.compare_to + ":\n  " + "\n  ".join(diffs))
            return 3
        print(f"summary reproduces {args.compare_to}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
