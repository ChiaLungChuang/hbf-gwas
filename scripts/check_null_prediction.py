#!/usr/bin/env python3
"""Check the step d prediction of lambda_GC and type I error on SIMULATED genotypes.

No real genotypes are read, and no phenotype touches real data. The check:

1. Simulate Balding-Nichols genotypes for seven subpopulations with the real
   sample counts (652 in all), plus two planted variants at the real sentinel
   frequencies.
2. Take the step d phenotype model: planted effects (BCL11A 0.456435, HBS1L-MYB
   0.571), the declared population shift, and noise.
3. Predict the distribution of the observed statistics with
   ``hbfgwas.predict.predicted_null_distribution``, from a 20,000-variant
   sample and 10 PCs, exactly as step d does.
4. Draw many phenotypes, run the genome-wide scan each time (without and with
   the PCs), and count how often each observed statistic falls inside each
   predicted interval.

If the prediction is right, the observed unadjusted lambda_GC falls inside the
predicted 1st-99th percentile interval in 98% of draws, inside the 5th-95th in
90%, and below the median in 50%.

    python scripts/check_null_prediction.py --out runs/<date>/null_prediction_check_simulated.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time

import numpy as np
from scipy import stats

from hbfgwas import assoc, predict, redact, simulate
from hbfgwas.config import load_config

COUNTS = {"ACB": 95, "ASW": 55, "ESN": 99, "GWD": 113, "LWK": 97, "MSL": 85, "YRI": 108}
SEED = 20261113  # this check's own seed; it touches no real data


def observed_stats(g, y, covariates, alpha):
    """lambda_GC and type I error of the scan of each column of y over the columns of g."""
    n = g.shape[0]
    q = assoc._covariate_basis(n, covariates)
    g_r = g - q @ (q.T @ g)
    g_r /= np.linalg.norm(g_r, axis=0)
    y_r = y - q @ (q.T @ y)
    y_r /= np.linalg.norm(y_r, axis=0)
    df = n - q.shape[1] - 1
    t_crit = stats.t.ppf(1 - alpha / 2, df)
    lam, t1 = [], []
    for a in range(0, y.shape[1], 50):
        r = g_r.T @ y_r[:, a : a + 50]
        t2 = r**2 * df / (1 - r**2)
        lam.append(np.median(t2, axis=0) / predict.MEDIAN_CHI2_1DF)
        t1.append(np.mean(t2 > t_crit**2, axis=0))
    return np.concatenate(lam), np.concatenate(t1)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-variants", type=int, default=60_000)
    ap.add_argument("--n-phenotypes", type=int, default=1000)
    ap.add_argument("--n-draws", type=int, default=2000)
    ap.add_argument("--fst", type=float, default=0.01)
    a = ap.parse_args(argv)
    cfg = load_config()
    shift = cfg["real_data_simulation"]["population_shift"]
    rng = np.random.default_rng(SEED)
    t0 = time.time()

    labels = sorted(COUNTS)
    sizes = [COUNTS[k] for k in labels]
    sim = simulate.simulate_genotypes(sizes, a.fst, rng, n_variants=a.n_variants)
    g = sim.genotypes.astype(float)
    pop = np.array(labels)[sim.pop]
    plant = simulate.simulate_genotypes(
        sizes, a.fst, rng, ancestral_freq=np.array([0.2393, 0.0391])
    ).genotypes.astype(float)
    betas = np.array([0.456435, 0.571])
    svec = np.array([shift[p] for p in pop])
    genetic = plant @ betas
    noise_var = predict.noise_variance(plant, betas, shift_var_fraction=float(svec.var()))
    fixed = genetic + svec

    pcs = assoc.pca(g[:, rng.choice(g.shape[1], 20_000, replace=False)], 10).pcs
    indicators = (pop[:, None] == np.array(labels)[None, :]).astype(float)
    sample = rng.choice(g.shape[1], 20_000, replace=False)
    pred = predict.predicted_null_distribution(
        g[:, sample],
        fixed,
        noise_var,
        np.hstack([indicators, pcs]),
        pcs,
        n_draws=a.n_draws,
        rng=np.random.default_rng(SEED + 1),
        n_effective=g.shape[1],
    )

    # many phenotypes with the same fixed part, new noise each time
    y = fixed[:, None] + rng.standard_normal((g.shape[0], a.n_phenotypes)) * np.sqrt(noise_var)
    results = {
        "what": "check of the step d null-statistic prediction on SIMULATED genotypes only",
        "seed": SEED,
        "n_samples": int(g.shape[0]),
        "population_counts": COUNTS,
        "fst": a.fst,
        "n_variants_observed": int(g.shape[1]),
        "n_variants_prediction": 20_000,
        "n_phenotypes": a.n_phenotypes,
        "n_prediction_draws": a.n_draws,
        "noise_var": noise_var,
        "sides": {},
    }
    bounds = {
        "pc_adjusted": {"lambda": (0.97, 1.05), "type1": (0.04, 0.06)},
        "unadjusted": {"lambda": None, "type1": None},
    }
    for side, cov in (("unadjusted", None), ("pc_adjusted", pcs)):
        lam_o, t1_o = observed_stats(g, y, cov, 0.05)
        p = pred[side]
        row = {"predicted": {"lambda_gc": p["lambda_gc"], "type1_error": p["type1_error"]}}
        for name, obs, summ in (
            ("lambda_gc", lam_o, p["lambda_gc"]),
            ("type1_error", t1_o, p["type1_error"]),
        ):
            row[name] = {
                "observed_quantiles": predict._quantile_summary(obs),
                "inside_q01_q99": float(np.mean((obs >= summ["q01"]) & (obs <= summ["q99"]))),
                "inside_q05_q95": float(np.mean((obs >= summ["q05"]) & (obs <= summ["q95"]))),
                "below_median": float(np.mean(obs < summ["median"])),
            }
        b = bounds[side]
        if b["lambda"]:
            lo, hi = b["lambda"]
            row["lambda_gc"]["observed_pass_rate"] = float(np.mean((lam_o >= lo) & (lam_o <= hi)))
            row["lambda_gc"]["predicted_pass_probability"] = predict.prob_in_interval(
                p["lambda_draws"], p["lambda_sd"], lo, hi
            )
            lo, hi = b["type1"]
            row["type1_error"]["observed_pass_rate"] = float(np.mean((t1_o >= lo) & (t1_o <= hi)))
            row["type1_error"]["predicted_pass_probability"] = predict.prob_in_interval(
                p["type1_draws"], p["type1_sd"], lo, hi
            )
        results["sides"][side] = row
    results["seconds"] = round(time.time() - t0, 1)
    redact.write_record(a.out, json.dumps(results, indent=2) + "\n")
    for side, row in results["sides"].items():
        lg = row["lambda_gc"]
        print(
            f"{side}: predicted lambda q01-q99 [{row['predicted']['lambda_gc']['q01']:.4f}, "
            f"{row['predicted']['lambda_gc']['q99']:.4f}], "
            f"observed inside {lg['inside_q01_q99']:.3f}"
            f" (q05-q95: {lg['inside_q05_q95']:.3f}; below median {lg['below_median']:.3f})"
        )
    print(f"wrote {a.out} ({results['seconds']} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
