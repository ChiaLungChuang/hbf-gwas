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
