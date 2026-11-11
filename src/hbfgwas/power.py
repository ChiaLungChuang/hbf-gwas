"""Analytic power for a single-variant additive association test.

The 1-df test statistic for a variant explaining a fraction ``r2`` of the
residual phenotypic variance is approximately noncentral chi-square with

    ncp = n * r2_partial / (1 - r2_partial),   r2_partial = r2 / (1 - covariate_r2)

where ``covariate_r2`` is the fraction of phenotypic variance explained by the
covariates (e.g. the part of a non-genetic shift that ancestry PCs absorb).
``ncp = n f^2`` with Cohen's ``f^2 = R^2 / (1 - R^2)``.

Power is ``P(chi2_1(ncp) > c)`` with ``c`` the central chi-square critical
value at ``alpha``. For a per-allele effect ``beta`` (in phenotype units) at
minor allele frequency ``maf`` under Hardy-Weinberg proportions,
``r2 = 2 maf (1 - maf) beta^2 / phenotype_var``.

All functions broadcast over numpy arrays.
"""

from __future__ import annotations

import numpy as np
from scipy import optimize, stats


def critical_value(alpha: float, df: int = 1) -> float:
    """Central chi-square critical value for a test of size ``alpha``."""
    if not 0 < alpha < 1:
        raise ValueError("alpha must be in (0, 1)")
    return float(stats.chi2.isf(alpha, df))


def r2_from_beta(beta, maf, phenotype_var: float = 1.0):
    """Variance explained by an additive variant: 2 p (1 - p) beta^2 / Var(y)."""
    beta = np.asarray(beta, dtype=float)
    maf = np.asarray(maf, dtype=float)
    if np.any((maf <= 0) | (maf >= 1)):
        raise ValueError("maf must be in (0, 1)")
    if phenotype_var <= 0:
        raise ValueError("phenotype_var must be positive")
    out = 2.0 * maf * (1.0 - maf) * beta**2 / phenotype_var
    return out.item() if out.ndim == 0 else out


def _check_r2(r2) -> np.ndarray:
    r2 = np.asarray(r2, dtype=float)
    if np.any((r2 < 0) | (r2 >= 1)):
        raise ValueError("r2 must be in [0, 1)")
    return r2


def ncp(n, r2, covariate_r2: float = 0.0):
    """Noncentrality parameter n * r2_partial / (1 - r2_partial)."""
    n = np.asarray(n, dtype=float)
    if np.any(n <= 0):
        raise ValueError("n must be positive")
    r2 = _check_r2(r2)
    if not 0 <= covariate_r2 < 1:
        raise ValueError("covariate_r2 must be in [0, 1)")
    r2_partial = _check_r2(r2 / (1.0 - covariate_r2))
    out = n * r2_partial / (1.0 - r2_partial)
    return out.item() if out.ndim == 0 else out


def power(
    n,
    alpha: float,
    r2=None,
    beta=None,
    maf=None,
    phenotype_var: float = 1.0,
    covariate_r2: float = 0.0,
):
    """Power of the 1-df additive test.

    Give either ``r2`` (variance explained) or both ``beta`` and ``maf``.
    """
    if r2 is None:
        if beta is None or maf is None:
            raise ValueError("give r2, or both beta and maf")
        r2 = r2_from_beta(beta, maf, phenotype_var)
    elif beta is not None or maf is not None:
        raise ValueError("give r2 or (beta, maf), not both")
    c = critical_value(alpha)
    lam = np.asarray(ncp(n, r2, covariate_r2), dtype=float)
    out = np.where(lam > 0, stats.ncx2.sf(c, 1, np.maximum(lam, 1e-300)), stats.chi2.sf(c, 1))
    return out.item() if out.ndim == 0 else out


def required_r2(n: float, alpha: float, target_power: float = 0.8, covariate_r2: float = 0.0):
    """Smallest variance explained that reaches ``target_power``."""
    if not alpha < target_power < 1:
        raise ValueError("target_power must be in (alpha, 1)")
    hi = 1.0 - 1e-12
    return float(
        optimize.brentq(
            lambda r2: power(n, alpha, r2=r2, covariate_r2=covariate_r2) - target_power,
            0.0,
            hi * (1.0 - covariate_r2),
            xtol=1e-14,
        )
    )


def required_n(r2: float, alpha: float, target_power: float = 0.8, covariate_r2: float = 0.0):
    """Sample size (not rounded) that reaches ``target_power``."""
    if not alpha < target_power < 1:
        raise ValueError("target_power must be in (alpha, 1)")
    if r2 <= 0:
        raise ValueError("r2 must be positive")
    return float(
        optimize.brentq(
            lambda n: power(n, alpha, r2=r2, covariate_r2=covariate_r2) - target_power,
            1e-6,
            1e9,
            xtol=1e-8,
        )
    )
