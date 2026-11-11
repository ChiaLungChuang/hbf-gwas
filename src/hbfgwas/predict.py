"""Predictions made from the genotypes alone, before any phenotype is drawn.

``conditional_power`` (issue #10): power of the additive t-test at each
planted variant, CONDITIONAL on the realised genotype matrix and covariates.
For variant j, with Q an orthonormal basis of [1, covariates] and residuals
taken after projecting on Q:

* the planted variant keeps only its residual genotype g_r (what the PCs leave);
* the expected estimate is E[beta_hat_j] = g_r . (G beta)_r / (g_r . g_r), i.e.
  its own beta plus "leakage" from the other planted loci it is correlated with;
* beta_hat's sampling spread comes from the noise only, but the residual
  variance the t-test ESTIMATES also contains whatever of the other loci (and
  of the subpopulation shift) the model leaves behind; that fixed part raises
  the effective critical value instead of adding spread;
* the statistic is then a noncentral t with df = n - rank([1, C]) - 1, judged
  against the t critical value (not the normal one).

The subpopulation shift is random (its subpopulation means are redrawn with
each phenotype), so power is averaged over draws from the same law that
``simulate.simulate_phenotype`` uses. Noise is handled exactly; only the
shift is averaged by Monte Carlo. Details and derivation: docs/power.md.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from hbfgwas.assoc import _covariate_basis
from hbfgwas.simulate import scale_shift_draws


def noise_variance(genotypes, betas, shift_var_fraction: float = 0.0) -> float:
    """Noise variance simulate_phenotype uses in beta mode: 1 - Var(G beta) - h."""
    g = np.asarray(genotypes, dtype=float)
    gen = g @ np.asarray(betas, dtype=float)
    return float(1.0 - gen.var() - shift_var_fraction)


def shift_vectors(pop, shift_var_fraction: float, n_draws: int, rng: np.random.Generator):
    """(n, n_draws) subpopulation shifts with the law of simulate_phenotype."""
    labels, idx = np.unique(np.asarray(pop), return_inverse=True)
    draws = rng.standard_normal((n_draws, labels.size))
    return scale_shift_draws(draws, idx, shift_var_fraction)[0]


def two_sided_power(crit, df: int, delta) -> np.ndarray:
    """P(|T| > crit) for T ~ noncentral t(df, delta), robust where scipy returns NaN.

    scipy's noncentral t can return NaN far in the tails (e.g. |delta| ~ 38);
    there the tail probability is 0 or 1 to double precision.
    """
    delta = np.asarray(delta, dtype=float)
    upper = stats.nct.sf(crit, df, delta)
    lower = stats.nct.sf(crit, df, -delta)
    upper = np.where(np.isnan(upper), (delta > crit).astype(float), upper)
    lower = np.where(np.isnan(lower), (-delta > crit).astype(float), lower)
    return upper + lower


def conditional_power(
    genotypes,
    betas,
    alpha: float,
    covariates=None,
    pop=None,
    shift_var_fraction: float = 0.0,
    noise_var: float | None = None,
    n_shift_draws: int = 2000,
    rng: np.random.Generator | None = None,
) -> list[dict]:
    """Power of the 1-df additive t-test at each planted variant, given the genotypes.

    ``genotypes`` (n, k) are counts of each planted variant's effect allele,
    ``betas`` (k,) their per-allele effects in phenotype SD units.
    """
    g = np.asarray(genotypes, dtype=float)
    if g.ndim == 1:
        g = g[:, None]
    betas = np.asarray(betas, dtype=float)
    n, k = g.shape
    if betas.shape != (k,):
        raise ValueError("one beta per planted variant")
    h = float(shift_var_fraction)
    if noise_var is None:
        noise_var = noise_variance(g, betas, h)
    if noise_var <= 0:
        raise ValueError("noise variance must be positive")

    q = _covariate_basis(n, covariates)
    df = n - q.shape[1] - 1
    crit = float(stats.t.isf(alpha / 2, df))
    g_r = g - q @ (q.T @ g)
    gen = g @ betas
    gen_r = gen - q @ (q.T @ gen)
    if h > 0:
        if pop is None:
            raise ValueError("a shift needs subpopulation labels")
        rng = rng if rng is not None else np.random.default_rng(0)
        s = shift_vectors(pop, h, n_shift_draws, rng)
        s_r = s - q @ (q.T @ s)
    else:
        s_r = np.zeros((n, 1))
    f = gen_r[:, None] + s_r  # fixed part of the residualised phenotype, per shift draw
    f_ss = np.sum(f**2, axis=0)

    out = []
    for j in range(k):
        gj = g_r[:, j]
        ss = float(gj @ gj)
        raw_ss = float(np.sum((g[:, j] - g[:, j].mean()) ** 2))
        if ss <= 1e-10 * max(raw_ss, 1.0):
            out.append({"power": float("nan"), "df": df, "note": "no genotype variance left"})
            continue
        m = (gj @ f) / ss  # E[beta_hat | shift draw]
        rss_fixed = np.maximum(f_ss - m**2 * ss, 0.0)
        # beta_hat ~ N(m, noise_var / ss) exactly; the residual sum of squares is
        # noise_var * noncentral chi2(df, lam) with lam = rss_fixed / noise_var,
        # independent of beta_hat. The fixed part (other loci, unabsorbed shift)
        # inflates the ESTIMATED standard error but not beta_hat's spread, so it
        # raises the effective critical value: |t| > crit  <=>
        # |(Z + delta) / sqrt(V / df)| > crit * sqrt(1 + lam / df), with V matched
        # to a central chi2(df) in mean.
        delta = m / np.sqrt(noise_var / ss)
        crit_eff = crit * np.sqrt(1.0 + rss_fixed / (noise_var * df))
        pw = two_sided_power(crit_eff, df, delta)
        m_gen = float(gj @ gen_r / ss)
        out.append(
            {
                "power": float(np.mean(pw)),
                # Monte Carlo error from averaging over shift draws (0 without a shift)
                "power_mc_se": float(np.std(pw) / np.sqrt(pw.size)) if pw.size > 1 else 0.0,
                "beta": float(betas[j]),
                "expected_beta": m_gen,
                "leakage": m_gen - float(betas[j]),
                "genotype_variance_kept": ss / raw_ss,
                "mean_noncentrality": float(np.mean(delta)),
                "mean_effective_critical_t": float(np.mean(crit_eff)),
                "noise_var": float(noise_var),
                "df": int(df),
                "critical_t": crit,
                "n_shift_draws": int(s_r.shape[1]) if h > 0 else 0,
            }
        )
    return out


# --- predicted lambda_GC (review item H) --------------------------------------------

MEDIAN_CHI2_1DF = float(stats.chi2.ppf(0.5, 1))


def mixture_median_chi2(delta2) -> float:
    """Median of an equal-weight mixture of noncentral chi-square(1, delta2_j).

    Uses chi2_1(d) = (Z + sqrt(d))^2, so CDF(x) = Phi(sqrt(x) - sqrt(d)) - Phi(-sqrt(x) - sqrt(d)).
    """
    d = np.sqrt(np.maximum(np.asarray(delta2, dtype=float), 0.0))

    def cdf(x):
        r = np.sqrt(x)
        return float(np.mean(stats.norm.cdf(r - d) - stats.norm.cdf(-r - d)))

    lo, hi = 0.0, max(1.0, float(np.max(d)) + 10.0) ** 2
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if cdf(mid) < 0.5:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-10 * max(1.0, hi):
            break
    return 0.5 * (lo + hi)


def _residualise(x: np.ndarray, q: np.ndarray) -> np.ndarray:
    return x - q @ (q.T @ x)


def _lambdas(genotypes, shifts, covariates, phenotype_var):
    """lambda_GC for each column of ``shifts`` (effective subpopulation shifts).

    For a null variant, t_j ~ N(delta_j, tau^2) with
    delta_j = sqrt(n) cov(g_r, s_r) / (sd(g_r) sd(y_r)) computed from the
    REALISED genotypes (so it already contains the within-subpopulation
    sampling deviations of the genotype means), and tau^2 = (phenotype_var -
    Var(s)) / Var(y_r): only the non-shift part of the phenotype is noise.
    Hence t_j^2 / tau^2 ~ chi2_1(delta_j^2 / tau^2), and lambda is the median
    of that mixture times tau^2, over 0.4549.
    """
    g = np.asarray(genotypes, dtype=float)
    n = g.shape[0]
    keep = g.var(axis=0) > 0
    g = g[:, keep]
    q = _covariate_basis(n, covariates)
    g_r = _residualise(g, q)
    sd_g = np.sqrt(np.sum(g_r**2, axis=0) / n)
    s = shifts - shifts.mean(axis=0)
    s_r = _residualise(s, q)
    nonshift = phenotype_var - s.var(axis=0)
    if np.any(nonshift <= 0):
        raise ValueError("shift variance exceeds the phenotype variance")
    var_y_r = nonshift + s_r.var(axis=0)
    tau2 = nonshift / var_y_r
    cov = (g_r.T @ s_r) / n  # (M, B)
    delta = np.sqrt(n) * cov / sd_g[:, None] / np.sqrt(var_y_r)[None, :]
    lam = [tau2[b] * mixture_median_chi2(delta[:, b] ** 2 / tau2[b]) for b in range(delta.shape[1])]
    return np.array(lam) / MEDIAN_CHI2_1DF


def lambda_sampling_sd(lam: float, n_variants: int) -> float:
    """SD of lambda_GC estimated as a median over ``n_variants`` independent statistics.

    Var(sample median) ~ 1 / (4 M f(m)^2) with f the chi2_1 density at its
    median m; relative to lambda this is 2.333 / sqrt(M). LD makes the
    effective M smaller than the variant count.
    """
    m = MEDIAN_CHI2_1DF
    f = float(stats.chi2.pdf(m, 1))
    return float(lam * 0.5 / (np.sqrt(n_variants) * f * m))


def predicted_lambda(
    genotypes,
    pop,
    shift_var_fraction: float,
    covariates=None,
    n_shift_draws: int = 200,
    rng: np.random.Generator | None = None,
    phenotype_var: float = 1.0,
    noise_var: float | None = None,
) -> dict:
    """Predicted lambda_GC from genotypes alone, before any phenotype is drawn.

    ``genotypes`` (n, M) are null variants (e.g. a random subset away from the
    planted sentinels). Each null t statistic is ~ N(delta_j, tau^2) with
    delta_j = sqrt(n) cov(g_r, s_r) / (sd(g_r) sd(y_r)) from the realised
    genotypes and tau^2 the non-shift share of the residual phenotype
    variance; lambda is the median of the implied chi-square mixture over the
    variants divided by 0.4549 (see ``_lambdas``).

    Two random parts act on ALL null variants at once and are drawn
    ``n_shift_draws`` times: the subpopulation shift (the simulation's own law),
    and the subpopulation means of the noise, N(0, noise_var / n_k), which act
    like a second, small shift. Their spread is why the summary is a median
    with 5th / 95th percentiles rather than one number.
    """
    rng = rng if rng is not None else np.random.default_rng(0)
    pop = np.asarray(pop)
    if noise_var is None:
        noise_var = phenotype_var - shift_var_fraction
    shifts = shift_vectors(pop, shift_var_fraction, n_shift_draws, rng)
    labels, idx = np.unique(pop, return_inverse=True)
    nk = np.bincount(idx).astype(float)
    noise_means = rng.standard_normal((n_shift_draws, labels.size)) * np.sqrt(noise_var / nk)
    effective = shifts + noise_means[:, idx].T
    lam = _lambdas(genotypes, effective, covariates, phenotype_var)
    n_var = int(np.sum(np.asarray(genotypes, dtype=float).var(axis=0) > 0))
    med = float(np.median(lam))
    return {
        "median": med,
        "q05": float(np.quantile(lam, 0.05)),
        "q95": float(np.quantile(lam, 0.95)),
        "mean": float(lam.mean()),
        "n_variants": n_var,
        "n_shift_draws": int(n_shift_draws),
        "note": (
            "q05-q95 reflect the random shift and noise means only; estimating lambda_GC "
            "as a median over M independent variants adds an SD of about "
            "2.333 * lambda / sqrt(M) (lambda_sampling_sd)"
        ),
    }
