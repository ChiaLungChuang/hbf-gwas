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
    shift=None,
) -> list[dict]:
    """Power of the 1-df additive t-test at each planted variant, given the genotypes.

    ``genotypes`` (n, k) are counts of each planted variant's effect allele,
    ``betas`` (k,) their per-allele effects in phenotype SD units. ``shift``
    (n,) is an exact, fixed shift (pre-registered on real genotypes); then
    nothing is averaged and ``pop``/``shift_var_fraction`` are not used.
    """
    g = np.asarray(genotypes, dtype=float)
    if g.ndim == 1:
        g = g[:, None]
    betas = np.asarray(betas, dtype=float)
    n, k = g.shape
    if betas.shape != (k,):
        raise ValueError("one beta per planted variant")
    fixed = None if shift is None else np.asarray(shift, dtype=float)
    h = float(fixed.var()) if fixed is not None else float(shift_var_fraction)
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
    if fixed is not None:
        if fixed.shape != (n,):
            raise ValueError("shift needs one value per individual")
        s_r = (fixed - q @ (q.T @ fixed))[:, None]
    elif h > 0:
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
                "n_shift_draws": 0 if (fixed is not None or h == 0) else int(s_r.shape[1]),
                "shift_mode": "fixed" if fixed is not None else ("random" if h > 0 else "none"),
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
    shift=None,
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
    fixed = None if shift is None else np.asarray(shift, dtype=float)
    if fixed is not None:
        shift_var_fraction = float(fixed.var())
    if noise_var is None:
        noise_var = phenotype_var - shift_var_fraction
    if fixed is not None:  # exact, pre-registered shift: only the noise means vary
        shifts = np.repeat(fixed[:, None], n_shift_draws, axis=1)
    else:
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
        "shift_mode": "fixed" if fixed is not None else "random",
        "note": (
            "q05-q95 reflect the random shift and noise means only; estimating lambda_GC "
            "as a median over M independent variants adds an SD of about "
            "2.333 * lambda / sqrt(M) (lambda_sampling_sd)"
        ),
    }


# --- predictive distribution of the null statistics (step d amendment) ------------


def _mixture_median_block(delta2: np.ndarray, iters: int = 48) -> np.ndarray:
    """Median of the equal-weight noncentral chi2(1, delta2[:, b]) mixture, per column."""
    from scipy.special import ndtr

    d = np.sqrt(np.maximum(np.asarray(delta2, dtype=float), 0.0))
    lo = np.zeros(d.shape[1])
    hi = (np.max(d, axis=0) + 10.0) ** 2
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        r = np.sqrt(mid)[None, :]
        cdf = np.mean(ndtr(r - d) - ndtr(-r - d), axis=0)
        below = cdf < 0.5
        lo = np.where(below, mid, lo)
        hi = np.where(below, hi, mid)
    return 0.5 * (lo + hi)


def _null_stats(g_r, sd_g, s_r, nonshift, z_crit, idx=None):
    """lambda_GC and expected type I error per column of s_r (residualised fixed parts)."""
    from scipy.special import ndtr

    n = g_r.shape[0]
    if idx is not None:
        g_r, sd_g = g_r[:, idx], sd_g[idx]
    var_y_r = nonshift + s_r.var(axis=0)
    tau2 = nonshift / var_y_r
    delta = np.sqrt(n) * ((g_r.T @ s_r) / n) / sd_g[:, None] / np.sqrt(var_y_r)[None, :]
    tau = np.sqrt(tau2)[None, :]
    lam = tau2 * _mixture_median_block(delta**2 / tau2[None, :]) / MEDIAN_CHI2_1DF
    t1 = np.mean(ndtr((-z_crit - delta) / tau) + 1.0 - ndtr((z_crit - delta) / tau), axis=0)
    return lam, t1


def predicted_null_distribution(
    genotypes,
    fixed,
    noise_var: float,
    structure,
    pcs,
    n_draws: int = 2000,
    rng: np.random.Generator | None = None,
    n_effective: int | None = None,
    type1_alpha: float = 0.05,
    n_bootstrap: int = 30,
    block: int = 50,
) -> dict:
    """Predictive distribution of lambda_GC and type I error, before any phenotype.

    The step e phenotype is ``y = fixed + e``. ``fixed`` (n,) is set once the
    genotypes and the pre-registered shift are set: planted effects plus the
    population shift. ``e ~ N(0, noise_var I)`` is the only random part. Its
    projection on ``structure`` (n, k), the population indicators and the PCs,
    is drawn exactly in each of ``n_draws`` draws and added to ``fixed``. The
    remainder acts on each null variant as independent noise, so its t statistic
    is ~ N(delta_jb, tau_b^2) (see ``_lambdas``). Per draw and per side
    (without and with the PCs as covariates):

    * lambda_b, the median of that chi-square mixture / 0.4549;
    * type1_b, the mean over variants of P(|t| > z_(1 - alpha/2)).

    The observed statistics differ from these by estimation noise, added once
    per draw as a normal deviate:

    * from taking a median (or fraction) over a finite number of variants:
      ``lambda_sampling_sd(lambda, n_effective)`` and
      sqrt(t1 (1 - t1) / n_effective);
    * from predicting with a sample of null variants: the bootstrap SD over
      variants, at the median draw.

    Returns per-draw arrays (``*_draws`` without estimation noise, ``*_obs``
    with it) and quantile summaries of the ``*_obs`` distributions.
    """
    rng = rng if rng is not None else np.random.default_rng(0)
    g = np.asarray(genotypes, dtype=float)
    g = g[:, g.var(axis=0) > 0]
    n, m = g.shape
    f = np.asarray(fixed, dtype=float)
    phenotype_var = float(f.var() + noise_var)
    basis = _covariate_basis(n, structure)  # [1, structure], orthonormal
    z = rng.standard_normal((basis.shape[1], n_draws))
    eff = f[:, None] + np.sqrt(noise_var) * (basis @ z)
    s = eff - eff.mean(axis=0)
    nonshift = phenotype_var - s.var(axis=0)
    if np.any(nonshift <= 0):
        raise ValueError("the fixed part and the drawn noise exceed the phenotype variance")
    z_crit = float(stats.norm.ppf(1 - type1_alpha / 2))
    n_eff = int(n_effective) if n_effective else m
    out: dict = {
        "n_draws": int(n_draws),
        "n_null_variants": int(m),
        "n_effective": n_eff,
        "noise_var": float(noise_var),
        "fixed_var": float(f.var()),
        "structure_dim": int(basis.shape[1]),
        "type1_alpha": type1_alpha,
    }
    eps = rng.standard_normal((2, n_draws))  # estimation noise, shared by both sides
    boot_idx = [rng.integers(0, m, m) for _ in range(n_bootstrap)]
    for side, cov in (("unadjusted", None), ("pc_adjusted", pcs)):
        q = _covariate_basis(n, cov)
        g_r = _residualise(g, q)
        sd_g = np.sqrt(np.sum(g_r**2, axis=0) / n)
        s_r = _residualise(s, q)
        lam = np.empty(n_draws)
        t1 = np.empty(n_draws)
        for a in range(0, n_draws, block):
            sl = slice(a, min(a + block, n_draws))
            lam[sl], t1[sl] = _null_stats(g_r, sd_g, s_r[:, sl], nonshift[sl], z_crit)
        b_mid = int(np.argsort(lam)[n_draws // 2])
        boots = [
            _null_stats(g_r, sd_g, s_r[:, [b_mid]], nonshift[[b_mid]], z_crit, idx)
            for idx in boot_idx
        ]
        boot_lam_sd = float(np.std([b[0][0] for b in boots], ddof=1))
        boot_t1_sd = float(np.std([b[1][0] for b in boots], ddof=1))
        lam_sd = np.sqrt(
            np.array([lambda_sampling_sd(x, n_eff) for x in lam]) ** 2 + boot_lam_sd**2
        )
        t1_sd = np.sqrt(t1 * (1 - t1) / n_eff + boot_t1_sd**2)
        lam_obs = lam + lam_sd * eps[0]
        t1_obs = t1 + t1_sd * eps[1]
        out[side] = {
            "lambda_draws": lam,
            "lambda_sd": lam_sd,
            "lambda_obs": lam_obs,
            "type1_draws": t1,
            "type1_sd": t1_sd,
            "type1_obs": t1_obs,
            "bootstrap_sd": {"lambda": boot_lam_sd, "type1": boot_t1_sd},
            "lambda_gc": _quantile_summary(lam_obs),
            "type1_error": _quantile_summary(t1_obs),
        }
    return out


def _quantile_summary(x) -> dict:
    x = np.asarray(x, dtype=float)
    qs = {"q01": 0.01, "q05": 0.05, "median": 0.5, "q95": 0.95, "q99": 0.99}
    out = {k: float(np.quantile(x, v)) for k, v in qs.items()}
    out["mean"] = float(x.mean())
    return out


def prob_in_interval(draws, sd, lo: float | None, hi: float | None) -> float:
    """P(observed in [lo, hi]) when observed ~ N(draw, sd^2), averaged over draws."""
    d, s = np.asarray(draws, dtype=float), np.asarray(sd, dtype=float)
    up = stats.norm.cdf((hi - d) / s) if hi is not None else 1.0
    dn = stats.norm.cdf((lo - d) / s) if lo is not None else 0.0
    return float(np.mean(up - dn))


def prob_calibration_passes(p_true: float, predicted: float, n: int, tol: float) -> float:
    """P(|X / n - predicted| <= tol) for X ~ Binomial(n, p_true): the power calibration."""
    lo = int(np.ceil((predicted - tol) * n - 1e-9))
    hi = int(np.floor((predicted + tol) * n + 1e-9))
    lo, hi = max(lo, 0), min(hi, n)
    if hi < lo:
        return 0.0
    return float(stats.binom.cdf(hi, n, p_true) - stats.binom.cdf(lo - 1, n, p_true))
