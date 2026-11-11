"""Experiment 1: population stratification on simulated genotypes.

Reads every setting from the pre-registration (``config/analysis.yaml``):
Balding-Nichols genotypes for seven subpopulations, a phenotype with planted
effects at three variants plus a subpopulation shift, association without and
with 10 PCs, and empirical power over repeated phenotype draws. The phenotype
is simulated; nothing here is a measured trait.

sim-v2 (after the Phase 1 review, issue #10): per-allele effects are planted
directly (beta mode, as on real genotypes), predicted power is conditional on
the realised genotypes and PCs (``hbfgwas.predict``) and computed before any
phenotype is drawn, and empirical power comes from ``power_reps`` phenotype
redraws tested at the planted variants only.

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
from hbfgwas import assoc, power, predict, simulate
from hbfgwas.config import DEFAULT_CONFIG


def _binomial_ci95(k: int, n: int) -> list[float]:
    """Clopper-Pearson 95% interval."""
    lo = 0.0 if k == 0 else float(stats.beta.ppf(0.025, k, n - k + 1))
    hi = 1.0 if k == n else float(stats.beta.ppf(0.975, k + 1, n - k))
    return [lo, hi]


def unconditional_predictions(config: dict) -> dict:
    """sim-v1 analytic predictions (noncentral chi-square, partial R2); reference only."""
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


def simulate_genotypes_and_pcs(config: dict) -> dict:
    """The experiment's genotypes (seed sim_genotypes) and PCs; no phenotype involved."""
    exp = config["simulated_experiment"]
    rng_g = np.random.default_rng(config["seeds"]["sim_genotypes"])
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
        ancestral_freq=np.array([p["ancestral_freq"] for p in exp["planted"]]),
    )
    pcs = assoc.pca(null.genotypes, n_pcs=exp["n_pcs"])
    return {"null": null, "plant": plant, "pcs": pcs, "pop": null.pop}


def conditional_predictions(config: dict, sim: dict | None = None) -> dict:
    """Power at each planted variant, conditional on the realised genotypes and PCs.

    Computed before any phenotype is drawn (genotypes, PCs, betas, and the
    shift's law only). The PC-adjusted values decide which loci are REQUIRED
    and are the predictions the power-calibration criterion is judged against.
    """
    exp = config["simulated_experiment"]
    sim = sim or simulate_genotypes_and_pcs(config)
    g = sim["plant"].genotypes.astype(float)
    betas = np.array([p["beta_sd"] for p in exp["planted"]])
    alpha = config["association"]["genome_wide_alpha"]
    h = exp["shift_var_fraction"]
    out = {}
    for side, cov in (("pc_adjusted", sim["pcs"].pcs), ("unadjusted", None)):
        res = predict.conditional_power(
            g,
            betas,
            alpha,
            covariates=cov,
            pop=sim["pop"],
            shift_var_fraction=h,
            n_shift_draws=exp["prediction_shift_draws"],
            rng=np.random.default_rng(config["seeds"]["sim_prediction"]),
        )
        out[side] = {p["locus"]: r for p, r in zip(exp["planted"], res, strict=True)}
    return out


def run_stratification(config: dict) -> dict:
    t0 = time.time()
    exp = config["simulated_experiment"]
    seeds = config["seeds"]
    alpha_gw = config["association"]["genome_wide_alpha"]
    alpha_t1 = config["association"]["type1_alpha"]
    n_pcs = exp["n_pcs"]
    planted = exp["planted"]
    tol = (config["acceptance"]["simulated"].get("power_calibration") or {}).get("max_abs_diff")

    sim = simulate_genotypes_and_pcs(config)
    null, plant, pcs, pop = sim["null"], sim["plant"], sim["pcs"], sim["pop"]
    n = pop.size
    monomorphic_null = int(np.sum(np.ptp(null.genotypes, axis=0) == 0))
    betas = [p["beta_sd"] for p in planted]
    h = exp["shift_var_fraction"]

    # Predictions first: genotypes and PCs only, no phenotype yet.
    cond = conditional_predictions(config, sim)

    # One genome-wide scan with the pre-registered phenotype seed.
    pheno = simulate.simulate_phenotype(
        plant.genotypes,
        pop,
        rng=np.random.default_rng(seeds["sim_phenotype"]),
        beta=betas,
        shift_var_fraction=h,
    )
    raw = assoc.linear_assoc(null.genotypes, pheno.y)
    adj = assoc.linear_assoc(null.genotypes, pheno.y, pcs.pcs)
    raw_pl = assoc.linear_assoc(plant.genotypes, pheno.y)
    adj_pl = assoc.linear_assoc(plant.genotypes, pheno.y, pcs.pcs)

    # Empirical power: phenotype redraws, tested at the planted variants only.
    reps = simulate.simulate_phenotype(
        plant.genotypes,
        pop,
        rng=np.random.default_rng(seeds["sim_power_reps_v2"]),
        beta=betas,
        shift_var_fraction=h,
        n_reps=exp["power_reps"],
    )
    rep_adj = assoc.linear_assoc(plant.genotypes, reps.y, pcs.pcs)
    rep_raw = assoc.linear_assoc(plant.genotypes, reps.y)
    hits = (rep_adj.p < alpha_gw).sum(axis=1)
    hits_raw = (rep_raw.p < alpha_gw).sum(axis=1)

    uncond = unconditional_predictions(config)
    declared = exp.get("predicted", {}).get("conditional_power_pc_adjusted", {})
    loci = []
    for j, p in enumerate(planted):
        name = p["locus"]
        k = int(hits[j])
        emp = k / exp["power_reps"]
        pa, pu = cond["pc_adjusted"][name], cond["unadjusted"][name]
        loci.append(
            {
                "locus": name,
                "beta_planted": float(pheno.betas[j]),
                "r2_realized": float(pheno.r2_realized[j]),
                "sample_freq": float(plant.genotypes[:, j].mean() / 2),
                "predicted_power": pa["power"],
                "predicted_power_mc_se": pa["power_mc_se"],
                "predicted_detail": {
                    k2: pa[k2]
                    for k2 in (
                        "expected_beta",
                        "leakage",
                        "genotype_variance_kept",
                        "mean_noncentrality",
                        "mean_effective_critical_t",
                    )
                },
                "predicted_power_unadjusted": pu["power"],
                "predicted_power_unconditional_v1": uncond["power_pc_adjusted"][name],
                "p_unadjusted": float(raw_pl.p[j]),
                "p_adjusted": float(adj_pl.p[j]),
                "beta_adjusted": float(adj_pl.beta[j]),
                "recovered": bool(adj_pl.p[j] < alpha_gw),
                "empirical_power": emp,
                "empirical_power_ci95": _binomial_ci95(k, exp["power_reps"]),
                "empirical_minus_predicted": emp - pa["power"],
                "power_calibration_tolerance": tol,
                "empirical_power_unadjusted": int(hits_raw[j]) / exp["power_reps"],
                "empirical_power_unadjusted_ci95": _binomial_ci95(
                    int(hits_raw[j]), exp["power_reps"]
                ),
            }
        )

    match = bool(declared) and all(
        abs(declared[name] - cond["pc_adjusted"][name]["power"]) < 5e-4 for name in declared
    )
    return {
        "experiment": "simulated",
        "prereg_version": config.get("prereg", {}).get("version"),
        "description": "stratification on simulated genotypes; phenotype simulated",
        "design": {
            "n": int(n),
            "n_per_subpop": list(exp["n_per_subpop"]),
            "n_null_variants": int(exp["n_null_variants"]),
            "monomorphic_null_variants": monomorphic_null,
            "fst": exp["fst"],
            "shift_var_fraction": h,
            "n_pcs": n_pcs,
            "planted_effect_mode": "beta (per-allele beta_sd)",
            "power_reps": exp["power_reps"],
            "prediction_shift_draws": exp["prediction_shift_draws"],
            "seeds": {
                k: seeds[k]
                for k in ("sim_genotypes", "sim_phenotype", "sim_power_reps_v2", "sim_prediction")
            },
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
        "predicted": {"unadjusted_lambda_gc": uncond["unadjusted_lambda_gc"]},
        "predictions_match_declared": match,
        "_pvalues": {"unadjusted": raw.p, "adjusted": adj.p},  # for the figure; not serialised
        "runtime_s": time.time() - t0,
    }


def summarize(results: dict) -> dict:
    """Compact, figure-adjacent summary (what docs/figures/ keeps)."""
    keys = (
        "locus",
        "predicted_power",
        "empirical_power",
        "empirical_power_ci95",
        "empirical_minus_predicted",
        "p_adjusted",
        "recovered",
    )
    return {
        "experiment": results["experiment"],
        "prereg_version": results.get("prereg_version"),
        "design": results["design"],
        "unadjusted": results["unadjusted"],
        "adjusted": results["adjusted"],
        "fst_estimate": results["fst_estimate"],
        "loci": [{k: loc[k] for k in keys} for loc in results["loci"]],
        "predicted_unadjusted_lambda_gc": results["predicted"]["unadjusted_lambda_gc"],
        "predictions_match_declared": results["predictions_match_declared"],
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
        for key in ("empirical_power", "predicted_power"):
            if abs(la[key] - lb[key]) > 2e-3:
                diffs.append(f"{la['locus']}.{key}: {la[key]} vs {lb[key]}")
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
