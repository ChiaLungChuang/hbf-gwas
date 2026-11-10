"""Tests for analytic power of a single-variant additive test.

For one degree of freedom the noncentral chi-square is the square of a shifted
standard normal, chi2_1(ncp) = (Z + sqrt(ncp))^2, so power has the closed form
    Phi(sqrt(ncp) - sqrt(c)) + Phi(-sqrt(ncp) - sqrt(c)),
with c the central chi-square critical value. The hand checks below use that
closed form (a different code path from scipy's ncx2).
"""

import numpy as np
import pytest
from scipy import stats

from hbfgwas import power as pw


def closed_form_power(ncp: float, alpha: float) -> float:
    c = stats.chi2.isf(alpha, 1)
    return stats.norm.cdf(np.sqrt(ncp) - np.sqrt(c)) + stats.norm.cdf(-np.sqrt(ncp) - np.sqrt(c))


def test_critical_value_genome_wide():
    # chi2_1 critical value at 5e-8 equals the squared two-sided z of 5.45131...
    assert pw.critical_value(5e-8) == pytest.approx(stats.norm.isf(2.5e-8) ** 2, rel=1e-10)
    assert pw.critical_value(5e-8) == pytest.approx(29.7168, abs=1e-3)


def test_ncp_formula():
    assert pw.ncp(n=1000, r2=0.05) == pytest.approx(1000 * 0.05 / 0.95)
    # covariates that explain 20% of phenotypic variance raise the partial R^2
    assert pw.ncp(n=1000, r2=0.04, covariate_r2=0.2) == pytest.approx(1000 * 0.05 / 0.95)


def test_textbook_80_percent_at_alpha_005():
    # ncp = (z_{0.975} + z_{0.80})^2 = 7.8489 gives 80% power at alpha = 0.05
    target_ncp = (stats.norm.isf(0.025) + stats.norm.isf(0.20)) ** 2
    n = 1000
    r2 = target_ncp / (n + target_ncp)  # inverts ncp = n r2 / (1 - r2)
    assert pw.power(n=n, alpha=0.05, r2=r2) == pytest.approx(0.80, abs=1e-4)


def test_brief_case_n650_genome_wide_r2_045():
    p = pw.power(n=650, alpha=5e-8, r2=0.045)
    assert 0.45 < p < 0.55  # "about 50% power near R2 = 4.5%"
    assert p == pytest.approx(closed_form_power(650 * 0.045 / 0.955, 5e-8), rel=1e-6)


def test_r2_for_half_power_at_n650_is_near_4_5_percent():
    r2 = pw.required_r2(n=650, alpha=5e-8, target_power=0.5)
    assert 0.04 < r2 < 0.05
    assert pw.power(n=650, alpha=5e-8, r2=r2) == pytest.approx(0.5, abs=1e-6)


def test_zero_effect_gives_alpha():
    for alpha in (0.05, 1e-3, 5e-8):
        assert pw.power(n=650, alpha=alpha, r2=0.0) == pytest.approx(alpha, rel=1e-9)


def test_beta_and_maf_give_r2():
    assert pw.r2_from_beta(beta=0.5, maf=0.2) == pytest.approx(2 * 0.2 * 0.8 * 0.25)
    assert pw.r2_from_beta(beta=0.5, maf=0.2, phenotype_var=2.0) == pytest.approx(0.04)
    assert pw.power(n=650, alpha=5e-8, beta=0.5, maf=0.2) == pytest.approx(
        pw.power(n=650, alpha=5e-8, r2=0.08)
    )


def test_power_increases_with_n_and_r2():
    ns = np.array([200, 400, 650, 1000, 5000])
    p_n = pw.power(n=ns, alpha=5e-8, r2=0.03)
    assert np.all(np.diff(p_n) > 0)
    r2s = np.array([0.0, 0.01, 0.03, 0.06, 0.1])
    p_r = pw.power(n=650, alpha=5e-8, r2=r2s)
    assert np.all(np.diff(p_r) > 0)


def test_closed_form_agreement_over_a_grid():
    for n in (100, 650, 2000):
        for r2 in (0.001, 0.02, 0.08):
            for alpha in (0.05, 5e-8):
                ncp = n * r2 / (1 - r2)
                assert pw.power(n=n, alpha=alpha, r2=r2) == pytest.approx(
                    closed_form_power(ncp, alpha), rel=1e-6, abs=1e-12
                )


def test_required_n_inverts_power():
    n = pw.required_n(r2=0.02, alpha=5e-8, target_power=0.8)
    assert pw.power(n=n, alpha=5e-8, r2=0.02) == pytest.approx(0.8, abs=1e-6)


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(n=650, alpha=5e-8),  # no effect given
        dict(n=650, alpha=5e-8, r2=0.1, beta=0.3, maf=0.2),  # both given
        dict(n=650, alpha=5e-8, beta=0.3),  # maf missing
        dict(n=650, alpha=5e-8, r2=1.0),  # r2 out of range
        dict(n=650, alpha=0.0, r2=0.1),
        dict(n=0, alpha=0.05, r2=0.1),
        dict(n=650, alpha=0.05, beta=0.3, maf=0.0),
    ],
)
def test_bad_inputs_raise(kwargs):
    with pytest.raises(ValueError):
        pw.power(**kwargs)


def test_analytic_power_matches_simulated_regression():
    """Empirical power of OLS on simulated data agrees with the analytic value."""
    rng = np.random.default_rng(20261110)
    n, r2, alpha, reps, maf = 400, 0.03, 1e-3, 2000, 0.3
    g = rng.binomial(2, maf, size=(reps, n)).astype(float)
    g -= g.mean(axis=1, keepdims=True)
    g /= g.std(axis=1, keepdims=True)
    y = np.sqrt(r2) * g + np.sqrt(1 - r2) * rng.standard_normal((reps, n))
    y -= y.mean(axis=1, keepdims=True)
    r = (g * y).sum(axis=1) / np.sqrt((g * g).sum(axis=1) * (y * y).sum(axis=1))
    t = r * np.sqrt((n - 2) / (1 - r**2))
    p = 2 * stats.t.sf(np.abs(t), n - 2)
    empirical = np.mean(p < alpha)
    analytic = pw.power(n=n, alpha=alpha, r2=r2)
    se = np.sqrt(analytic * (1 - analytic) / reps)
    assert abs(empirical - analytic) < 4 * se, (empirical, analytic)
