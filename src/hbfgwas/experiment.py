"""Experiment 1: population stratification on simulated genotypes.

Reads every setting from the pre-registration (``config/analysis.yaml``):
Balding-Nichols genotypes for seven subpopulations, a phenotype with planted
effects at three variants plus a subpopulation shift, association without and
with 10 PCs, and empirical power over repeated phenotype draws. The phenotype
is simulated; nothing here is a measured trait.

``run_stratification`` returns a results dict in the format the gate reads;
``write_outputs`` writes the full results JSON, a compact summary JSON, and the
two-panel QQ figure.
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy import stats

import hbfgwas
from hbfgwas import assoc, power, simulate
from hbfgwas.config import DEFAULT_CONFIG


def _binomial_ci95(k: int, n: int) -> list[float]:
    """Clopper-Pearson 95% interval."""
    lo = 0.0 if k == 0 else float(stats.beta.ppf(0.025, k, n - k + 1))
    hi = 1.0 if k == n else float(stats.beta.ppf(0.975, k + 1, n - k))
    return [lo, hi]


def predictions(config: dict) -> dict:
    """Analytic predictions recomputed from the config's design."""
    exp = config["simulated_experiment"]
    n = int(sum(exp["n_per_subpop"]))
    alpha = config["association"]["genome_wide_alpha"]
    h = exp["shift_var_fraction"]
    return {
        "unadjusted_lambda_gc": simulate.expected_null_chi2(
            n=n, fst=exp["fst"], shift_var_fraction=h, n_pops=len(exp["n_per_subpop"])
        ),
        "power_pc_adjusted": {
            p["locus"]: float(power.power(n=n, alpha=alpha, r2=p["r2"], covariate_r2=h))
            for p in exp["planted"]
        },
    }


def run_stratification(config: dict) -> dict:
    t0 = time.time()
    exp = config["simulated_experiment"]
    seeds = config["seeds"]
    alpha_gw = config["association"]["genome_wide_alpha"]
    alpha_t1 = config["association"]["type1_alpha"]
    n_pcs = exp["n_pcs"]
    planted = exp["planted"]

    rng_g = np.random.default_rng(seeds["sim_genotypes"])
    null = simulate.simulate_genotypes(
        exp["n_per_subpop"],
        fst=exp["fst"],
        rng=rng_g,
        n_variants=exp["n_null_variants"],
        freq_range=tuple(exp["ancestral_freq_range"]),
    )
    plant = simulate.simulate_genotypes(
        exp["n_per_subpop"],
        fst=exp["fst"],
        rng=rng_g,
        ancestral_freq=np.array([p["ancestral_freq"] for p in planted]),
    )
    pop = null.pop
    n = pop.size
    monomorphic_null = int(np.sum(np.ptp(null.genotypes, axis=0) == 0))

    pcs = assoc.pca(null.genotypes, n_pcs=n_pcs)

    r2 = [p["r2"] for p in planted]
    h = exp["shift_var_fraction"]
    pheno = simulate.simulate_phenotype(
        plant.genotypes,
        pop,
        rng=np.random.default_rng(seeds["sim_phenotype"]),
        r2=r2,
        shift_var_fraction=h,
    )

    raw = assoc.linear_assoc(null.genotypes, pheno.y)
    adj = assoc.linear_assoc(null.genotypes, pheno.y, pcs.pcs)
    raw_pl = assoc.linear_assoc(plant.genotypes, pheno.y)
    adj_pl = assoc.linear_assoc(plant.genotypes, pheno.y, pcs.pcs)

    reps = simulate.simulate_phenotype(
        plant.genotypes,
        pop,
        rng=np.random.default_rng(seeds["sim_power_reps"]),
        r2=r2,
        shift_var_fraction=h,
        n_reps=exp["power_reps"],
    )
    rep_adj = assoc.linear_assoc(plant.genotypes, reps.y, pcs.pcs)
    rep_raw = assoc.linear_assoc(plant.genotypes, reps.y)
    hits = (rep_adj.p < alpha_gw).sum(axis=1)
    hits_raw = (rep_raw.p < alpha_gw).sum(axis=1)

    pred = predictions(config)
    declared = exp.get("predicted", {})
    loci = []
    for j, p in enumerate(planted):
        k = int(hits[j])
        loci.append(
            {
                "locus": p["locus"],
                "r2_planted": float(pheno.r2_realized[j]),
                "beta_planted": float(pheno.betas[j]),
                "sample_freq": float(plant.genotypes[:, j].mean() / 2),
                "predicted_power": pred["power_pc_adjusted"][p["locus"]],
                "p_unadjusted": float(raw_pl.p[j]),
                "p_adjusted": float(adj_pl.p[j]),
                "beta_adjusted": float(adj_pl.beta[j]),
                "recovered": bool(adj_pl.p[j] < alpha_gw),
                "empirical_power": k / exp["power_reps"],
                "empirical_power_ci95": _binomial_ci95(k, exp["power_reps"]),
                "empirical_power_unadjusted": int(hits_raw[j]) / exp["power_reps"],
            }
        )

    declared_power = declared.get("power_pc_adjusted", {})
    match = bool(declared_power) and all(
        abs(declared_power[loc] - pred["power_pc_adjusted"][loc]) < 5e-4 for loc in declared_power
    )
    return {
        "experiment": "simulated",
        "description": "stratification on simulated genotypes; phenotype simulated",
        "design": {
            "n": int(n),
            "n_per_subpop": list(exp["n_per_subpop"]),
            "n_null_variants": int(exp["n_null_variants"]),
            "monomorphic_null_variants": monomorphic_null,
            "fst": exp["fst"],
            "shift_var_fraction": h,
            "n_pcs": n_pcs,
            "power_reps": exp["power_reps"],
            "seeds": {k: seeds[k] for k in ("sim_genotypes", "sim_phenotype", "sim_power_reps")},
        },
        "fst_estimate": simulate.estimate_fst(null.genotypes, pop),
        "pca_eigenvalues": pcs.eigenvalues.tolist(),
        "phenotype_variance": {
            "genetic": float(pheno.genetic.var()),
            "shift": float(pheno.shift.var()),
            "noise": float(pheno.noise.var()),
            "total": float(pheno.y.var()),
        },
        "unadjusted": {
            "lambda_gc": assoc.lambda_gc(raw.p),
            "type1_error": assoc.type1_error(raw.p, alpha_t1),
            "type1_alpha": alpha_t1,
            "n_covariates": 0,
        },
        "adjusted": {
            "lambda_gc": assoc.lambda_gc(adj.p),
            "type1_error": assoc.type1_error(adj.p, alpha_t1),
            "type1_alpha": alpha_t1,
            "n_covariates": n_pcs,
        },
        "loci": loci,
        "predicted": pred,
        "predictions_match_declared": match,
        "_pvalues": {"unadjusted": raw.p, "adjusted": adj.p},  # for the figure; not serialised
        "runtime_s": time.time() - t0,
    }


def summarize(results: dict) -> dict:
    """Compact, figure-adjacent summary (what docs/figures/ keeps)."""
    return {
        "experiment": results["experiment"],
        "design": results["design"],
        "unadjusted": results["unadjusted"],
        "adjusted": results["adjusted"],
        "fst_estimate": results["fst_estimate"],
        "loci": [
            {
                k: loc[k]
                for k in (
                    "locus",
                    "predicted_power",
                    "p_adjusted",
                    "recovered",
                    "empirical_power",
                    "empirical_power_ci95",
                )
            }
            for loc in results["loci"]
        ],
        "predicted_unadjusted_lambda_gc": results["predicted"]["unadjusted_lambda_gc"],
    }


def compare_summaries(a: dict, b: dict, tol: float = 1e-4) -> list[str]:
    """Differences between two summaries beyond ``tol`` (log10 scale for p-values)."""
    diffs = []
    for side in ("unadjusted", "adjusted"):
        for key in ("lambda_gc", "type1_error"):
            if abs(a[side][key] - b[side][key]) > tol:
                diffs.append(f"{side}.{key}: {a[side][key]} vs {b[side][key]}")
    for la, lb in zip(a["loci"], b["loci"], strict=True):
        if la["locus"] != lb["locus"]:
            diffs.append(f"locus order: {la['locus']} vs {lb['locus']}")
            continue
        if abs(np.log10(la["p_adjusted"]) - np.log10(lb["p_adjusted"])) > 1e-3:
            diffs.append(f"{la['locus']}.p_adjusted: {la['p_adjusted']} vs {lb['p_adjusted']}")
        if la["recovered"] != lb["recovered"]:
            diffs.append(f"{la['locus']}.recovered differs")
        if abs(la["empirical_power"] - lb["empirical_power"]) > 2e-3:
            diffs.append(
                f"{la['locus']}.empirical_power: {la['empirical_power']} vs {lb['empirical_power']}"
            )
    return diffs


def _git_commit() -> str | None:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            cwd=Path(__file__).resolve().parent,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def provenance(config_path: Path | None = None) -> dict:
    config_path = Path(config_path) if config_path else DEFAULT_CONFIG
    return {
        "hbfgwas_version": hbfgwas.__version__,
        "git_commit": _git_commit(),
        "config_path": str(config_path),
        "config_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "platform": platform.platform(),
        "run_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }


def make_figure(results: dict, path: Path) -> Path:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), sharex=True, sharey=True)
    for ax, side, label, color in (
        (axes[0], "unadjusted", "No ancestry adjustment", "#c0392b"),
        (axes[1], "adjusted", f"Adjusted for {results['design']['n_pcs']} PCs", "#2471a3"),
    ):
        assoc.qq_plot(ax, results["_pvalues"][side], title=label, color=color)
        t1 = results[side]["type1_error"]
        ax.text(
            0.04,
            0.96,
            f"type I error at 0.05: {t1:.4f}",
            transform=ax.transAxes,
            va="top",
            fontsize=9,
        )
    d = results["design"]
    fig.suptitle(
        f"Null variants, simulated genotypes (n = {d['n']}, {d['n_null_variants']:,} variants, "
        f"{len(d['n_per_subpop'])} subpopulations, Fst = {d['fst']}); phenotype simulated",
        fontsize=10,
    )
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _jsonable(results: dict) -> dict:
    return {k: v for k, v in results.items() if not k.startswith("_")}


def write_outputs(
    results: dict, config: dict, out_dir: Path, config_path: Path | None = None
) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    full = _jsonable(results) | {"provenance": provenance(config_path)}
    paths = {
        "results": out_dir / "results.json",
        "summary": out_dir / "sim_stratification_summary.json",
        "figure": out_dir / "sim_stratification_qq.png",
    }
    paths["results"].write_text(json.dumps(full, indent=2) + "\n")
    paths["summary"].write_text(json.dumps(summarize(results), indent=2) + "\n")
    make_figure(results, paths["figure"])
    return paths
