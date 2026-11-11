"""Tests for predictions made from genotypes before any phenotype is drawn.

conditional_power must agree with Monte Carlo power of the actual test
(hbfgwas.assoc.linear_assoc on phenotypes from hbfgwas.simulate) for a FIXED
genotype matrix, with and without covariates. These tests use their own seeds.
"""

import numpy as np
import pytest
from scipy import stats

from hbfgwas import assoc, predict
from hbfgwas import simulate as sim


@pytest.fixture(scope="module")
def scenario():
    rng = np.random.default_rng(404)
    sizes = [150, 150, 150]
    null = sim.simulate_genotypes(sizes, fst=0.05, rng=rng, n_variants=3000)
    planted = sim.simulate_genotypes(
        sizes, fst=0.05, rng=rng, ancestral_freq=np.array([0.3, 0.1, 0.2])
    )
    pcs = assoc.pca(null.genotypes, n_pcs=5).pcs
    return dict(
        g=planted.genotypes.astype(float),
        pop=null.pop,
        pcs=pcs,
        betas=np.array([0.25, 0.30, 0.20]),
        h=0.15,
        alpha=1e-4,
    )


def _empirical(s, covariates, reps=4000, seed=405):
    ph = sim.simulate_phenotype(
        s["g"],
        s["pop"],
        rng=np.random.default_rng(seed),
        beta=s["betas"],
        shift_var_fraction=s["h"],
        n_reps=reps,
    )
    res = assoc.linear_assoc(s["g"], ph.y, covariates)
    return (res.p < s["alpha"]).mean(axis=1), reps


@pytest.mark.parametrize("adjusted", [True, False])
def test_conditional_power_matches_monte_carlo(scenario, adjusted):
    s = scenario
    cov = s["pcs"] if adjusted else None
    pred = predict.conditional_power(
        s["g"],
        s["betas"],
        s["alpha"],
        covariates=cov,
        pop=s["pop"],
        shift_var_fraction=s["h"],
        n_shift_draws=2000,
        rng=np.random.default_rng(406),
    )
    emp, reps = _empirical(s, cov)
    for j, p in enumerate(pred):
        se_emp = np.sqrt(max(p["power"] * (1 - p["power"]), 1e-4) / reps)
        se = np.hypot(se_emp, p["power_mc_se"])
        assert abs(emp[j] - p["power"]) < 4 * se, (j, emp[j], p["power"], se)


def test_leakage_from_a_correlated_planted_locus():
    rng = np.random.default_rng(1)
    g1 = rng.binomial(2, 0.3, 500).astype(float)
    g = np.column_stack([g1, g1])  # identical columns: each picks up the other's effect
    pred = predict.conditional_power(g, np.array([0.2, 0.1]), 0.05)
    assert pred[0]["expected_beta"] == pytest.approx(0.3)
    assert pred[0]["leakage"] == pytest.approx(0.1)


def test_single_variant_no_shift_matches_noncentral_t_by_hand():
    rng = np.random.default_rng(2)
    n = 400
    g = rng.binomial(2, 0.25, n).astype(float)[:, None]
    beta, alpha = 0.15, 1e-3
    pred = predict.conditional_power(g, np.array([beta]), alpha)[0]
    gc = g[:, 0] - g[:, 0].mean()
    noise = 1 - (beta * gc).var()
    df = n - 2
    delta = beta * np.sqrt(gc @ gc / noise)
    c = stats.t.isf(alpha / 2, df)
    expect = stats.nct.sf(c, df, delta) + stats.nct.cdf(-c, df, delta)
    assert pred["power"] == pytest.approx(expect, rel=1e-9)
    assert pred["df"] == df


def test_default_noise_variance_matches_simulate_phenotype(scenario):
    s = scenario
    ph = sim.simulate_phenotype(
        s["g"], s["pop"], rng=np.random.default_rng(3), beta=s["betas"], shift_var_fraction=s["h"]
    )
    assert predict.noise_variance(s["g"], s["betas"], s["h"]) == pytest.approx(
        1 - ph.genetic.var() - s["h"]
    )


def test_shift_vectors_have_the_simulation_law(scenario):
    s = scenario
    sh = predict.shift_vectors(s["pop"], s["h"], 50, np.random.default_rng(4))
    assert sh.shape == (len(s["pop"]), 50)
    np.testing.assert_allclose(sh.var(axis=0), s["h"], rtol=1e-10)
    for k in np.unique(s["pop"]):
        assert np.ptp(sh[s["pop"] == k], axis=0).max() < 1e-12


# --- predicted lambda_GC (review item H) ------------------------------------------


@pytest.fixture(scope="module")
def lam_scenario():
    rng = np.random.default_rng(808)
    sizes = [150, 150, 150]
    g = sim.simulate_genotypes(sizes, fst=0.05, rng=rng, n_variants=4000)
    pcs = assoc.pca(g.genotypes[:, :2000], n_pcs=5).pcs
    return dict(g=g.genotypes[:, 2000:].astype(float), pop=g.pop, pcs=pcs, h=0.15)


def test_mixture_median_of_central_chi2_is_the_chi2_median():
    m = predict.mixture_median_chi2(np.zeros(100))
    assert m == pytest.approx(stats.chi2.ppf(0.5, 1), rel=1e-6)


def test_mixture_median_matches_simulation():
    rng = np.random.default_rng(9)
    d2 = rng.uniform(0, 4, 500)
    draws = (rng.standard_normal((400, 500)) + np.sqrt(d2)) ** 2
    assert predict.mixture_median_chi2(d2) == pytest.approx(np.median(draws), rel=0.01)


@pytest.mark.parametrize("adjusted", [False, True])
def test_predicted_lambda_matches_observed_lambda_over_phenotype_draws(lam_scenario, adjusted):
    """Prediction (genotypes only) vs lambda_GC observed on 40 simulated phenotypes."""
    s = lam_scenario
    cov = s["pcs"] if adjusted else None
    pred = predict.predicted_lambda(
        s["g"], s["pop"], s["h"], covariates=cov, n_shift_draws=400, rng=np.random.default_rng(12)
    )
    ph = sim.simulate_phenotype(
        s["g"][:, :1],
        s["pop"],
        rng=np.random.default_rng(13),
        r2=[0.0],
        shift_var_fraction=s["h"],
        n_reps=40,
    )
    res = assoc.linear_assoc(s["g"], ph.y, cov)
    observed = np.array([assoc.lambda_gc(res.p[:, r]) for r in range(40)])
    se = observed.std(ddof=1) / np.sqrt(40)
    assert abs(observed.mean() - pred["mean"]) < 3 * se + 0.01, (observed.mean(), pred)
    # the interval plus the sampling error of a median over 2000 variants
    samp = predict.lambda_sampling_sd(pred["median"], pred["n_variants"])
    lo, hi = pred["q05"] - 1.645 * samp, pred["q95"] + 1.645 * samp
    inside = np.mean((observed >= lo) & (observed <= hi))
    assert inside >= 0.75, (inside, lo, hi, observed.min(), observed.max())


def test_predicted_lambda_weak_structure_matches_observed():
    """Weak structure (the experiment's own design: 7 x 93, Fst 0.01, h 0.2).

    Sensitive to double counting the within-subpopulation sampling part of
    the genotype means, which shifts lambda by about h.
    """
    rng = np.random.default_rng(909)
    g = sim.simulate_genotypes([93] * 7, fst=0.01, rng=rng, n_variants=6000)
    null = g.genotypes[:, :4000].astype(float)
    pred = predict.predicted_lambda(
        null, g.pop, 0.2, n_shift_draws=300, rng=np.random.default_rng(910)
    )
    ph = sim.simulate_phenotype(
        null[:, :1],
        g.pop,
        rng=np.random.default_rng(911),
        r2=[0.0],
        shift_var_fraction=0.2,
        n_reps=40,
    )
    res = assoc.linear_assoc(null, ph.y)
    observed = np.array([assoc.lambda_gc(res.p[:, r]) for r in range(40)])
    se = observed.std(ddof=1) / np.sqrt(40)
    assert abs(observed.mean() - pred["mean"]) < 3 * se + 0.01, (observed.mean(), pred["mean"], se)


def test_lambda_sampling_sd():
    # 100,000 independent null variants: SD ~ 0.0074 (the sim-v1 design note)
    assert predict.lambda_sampling_sd(1.0, 100_000) == pytest.approx(0.00738, abs=2e-4)


def test_predicted_lambda_summary(lam_scenario):
    s = lam_scenario
    out = predict.predicted_lambda(
        s["g"], s["pop"], s["h"], n_shift_draws=50, rng=np.random.default_rng(11)
    )
    adj = predict.predicted_lambda(
        s["g"],
        s["pop"],
        s["h"],
        covariates=s["pcs"],
        n_shift_draws=50,
        rng=np.random.default_rng(11),
    )
    assert out["q05"] <= out["median"] <= out["q95"]
    assert out["median"] > 1.5
    assert adj["median"] == pytest.approx(1.0, abs=0.06)
    assert out["n_variants"] == 2000 and out["n_shift_draws"] == 50
