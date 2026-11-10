"""Reference implementation of additive linear association, PCA and summaries.

This module is the tested reference, not the genome-wide engine: plink2
``--glm`` does the genome-wide scan, and :func:`compare_to_plink2` checks that
plink2 and this code agree on the same variants within a tolerance.

Model, per variant ``j``: ``y = a + b_j g_j + C c + e`` with ``C`` the
covariates. Fitting uses the Frisch-Waugh-Lovell theorem: residualise ``y``
and every ``g_j`` on ``[1, C]`` once, then each ``b_j`` is a 1-D regression of
residuals, so all variants (and many phenotypes) are fitted with a few matrix
products. Residual degrees of freedom are ``n - rank([1, C]) - 1`` and p-values
use the t distribution, as for ordinary least squares.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

MEDIAN_CHI2_1DF = float(stats.chi2.ppf(0.5, 1))  # 0.454936...


@dataclass
class AssocResult:
    beta: np.ndarray  # (m,) or (m, R)
    se: np.ndarray
    t: np.ndarray
    p: np.ndarray
    df: int
    n: int

    def to_frame(self, ids=None) -> pd.DataFrame:
        if self.beta.ndim != 1:
            raise ValueError("to_frame needs a single phenotype")
        m = self.beta.size
        return pd.DataFrame(
            {
                "ID": list(ids) if ids is not None else np.arange(m),
                "BETA": self.beta,
                "SE": self.se,
                "T_STAT": self.t,
                "P": self.p,
                "OBS_CT": self.n,
                "DF": self.df,
            }
        )


def _covariate_basis(n: int, covariates) -> np.ndarray:
    """Orthonormal basis of the column space of [1, covariates]."""
    cols = [np.ones((n, 1))]
    if covariates is not None:
        c = np.asarray(covariates, dtype=float)
        if c.ndim == 1:
            c = c[:, None]
        if c.shape[0] != n:
            raise ValueError("covariates must have one row per individual")
        if np.isnan(c).any():
            raise ValueError("covariates contain missing values")
        cols.append(c)
    design = np.hstack(cols)
    u, s, _ = np.linalg.svd(design, full_matrices=False)
    keep = s > s.max() * max(design.shape) * np.finfo(float).eps
    return u[:, keep]


def linear_assoc(genotypes, y, covariates=None, chunk: int = 10_000) -> AssocResult:
    """Additive linear regression of ``y`` on each genotype column, with covariates.

    ``genotypes`` is (n, m) allele counts (no missing values); ``y`` is (n,) or
    (n, R) for R phenotypes at once. Variants with no variance left after
    covariate adjustment get NaN.
    """
    g_all = np.asarray(genotypes)
    if g_all.ndim != 2:
        raise ValueError("genotypes must be (n, m)")
    yy = np.asarray(y, dtype=float)
    single = yy.ndim == 1
    if single:
        yy = yy[:, None]
    n, m = g_all.shape
    if yy.shape[0] != n:
        raise ValueError("y must have one row per individual")
    if np.isnan(yy).any():
        raise ValueError("y contains missing values")

    q = _covariate_basis(n, covariates)
    df = n - q.shape[1] - 1
    if df < 1:
        raise ValueError("not enough individuals for this many covariates")
    y_r = yy - q @ (q.T @ yy)
    y_ss = np.sum(y_r**2, axis=0)  # (R,)

    shape = (m, yy.shape[1])
    beta, se, tstat, pval = (np.full(shape, np.nan) for _ in range(4))
    for c0 in range(0, m, chunk):
        c1 = min(c0 + chunk, m)
        g = g_all[:, c0:c1].astype(float)
        if np.isnan(g).any():
            raise ValueError("genotypes contain missing values")
        g_r = g - q @ (q.T @ g)
        g_ss = np.sum(g_r**2, axis=0)
        raw_ss = np.sum((g - g.mean(axis=0)) ** 2, axis=0)
        ok = g_ss > 1e-10 * np.maximum(raw_ss, 1.0)
        if not ok.any():
            continue
        gy = g_r[:, ok].T @ y_r  # (k, R)
        b = gy / g_ss[ok, None]
        rss = y_ss[None, :] - b * gy
        s = np.sqrt(np.maximum(rss, 0.0) / df / g_ss[ok, None])
        t = b / s
        idx = np.arange(c0, c1)[ok]
        beta[idx], se[idx], tstat[idx] = b, s, t
        pval[idx] = 2.0 * stats.t.sf(np.abs(t), df)
    if single:
        beta, se, tstat, pval = beta[:, 0], se[:, 0], tstat[:, 0], pval[:, 0]
    return AssocResult(beta=beta, se=se, t=tstat, p=pval, df=int(df), n=int(n))


@dataclass
class PCAResult:
    pcs: np.ndarray  # (n, k) unit-norm eigenvectors of the GRM
    eigenvalues: np.ndarray  # (k,) descending
    variance_explained: np.ndarray  # eigenvalues / trace
    n_variants_used: int


def pca(genotypes, n_pcs: int = 10, chunk: int = 10_000) -> PCAResult:
    """PCA of the standardised genotype matrix via the n x n GRM.

    Each variant is standardised as ``(g - 2p) / sqrt(2p(1-p))`` with ``p`` the
    sample allele frequency; monomorphic variants are dropped. The GRM
    ``Z Z^T / m`` is accumulated in chunks and fully eigendecomposed (exact,
    not randomised).
    """
    g_all = np.asarray(genotypes)
    n, m = g_all.shape
    if not 1 <= n_pcs < n:
        raise ValueError("n_pcs must be in [1, n)")
    grm = np.zeros((n, n))
    used = 0
    for c0 in range(0, m, chunk):
        g = g_all[:, c0 : c0 + chunk].astype(float)
        p = g.mean(axis=0) / 2.0
        keep = (p > 0) & (p < 1)
        if not keep.any():
            continue
        z = (g[:, keep] - 2.0 * p[keep]) / np.sqrt(2.0 * p[keep] * (1.0 - p[keep]))
        grm += z @ z.T
        used += int(keep.sum())
    if used == 0:
        raise ValueError("no polymorphic variants")
    grm /= used
    w, v = np.linalg.eigh(grm)
    order = np.argsort(w)[::-1][:n_pcs]
    return PCAResult(
        pcs=v[:, order],
        eigenvalues=w[order],
        variance_explained=w[order] / np.trace(grm),
        n_variants_used=used,
    )


def lambda_gc(p) -> float:
    """Genomic control inflation: median 1-df chi-square over its null median."""
    p = np.asarray(p, dtype=float).ravel()
    p = p[np.isfinite(p)]
    if p.size == 0:
        raise ValueError("no finite p-values")
    return float(np.median(stats.chi2.isf(p, 1)) / MEDIAN_CHI2_1DF)


def type1_error(p_null, alpha: float = 0.05) -> float:
    """Fraction of null variants with p < alpha."""
    p = np.asarray(p_null, dtype=float).ravel()
    p = p[np.isfinite(p)]
    return float(np.mean(p < alpha))


def empirical_power(p_planted, threshold: float = 5e-8) -> np.ndarray:
    """Per-variant fraction of repetitions with p < threshold; input (variants, reps)."""
    p = np.atleast_2d(np.asarray(p_planted, dtype=float))
    return np.mean(p < threshold, axis=1)


def qq_plot(ax, p, title: str | None = None, ci: float = 0.95, color: str = "C0", label=None):
    """QQ plot of -log10 p against the uniform expectation, with a pointwise band."""
    p = np.asarray(p, dtype=float).ravel()
    p = np.sort(p[np.isfinite(p) & (p > 0)])
    n = p.size
    i = np.arange(1, n + 1)
    expected = -np.log10(i / (n + 1))
    observed = -np.log10(p)
    lo = -np.log10(stats.beta.ppf((1 + ci) / 2, i, n - i + 1))
    hi = -np.log10(stats.beta.ppf((1 - ci) / 2, i, n - i + 1))
    ax.fill_between(expected, lo, hi, color="0.85", lw=0, label=f"{int(ci * 100)}% band")
    top = max(expected.max(), observed.max())
    ax.plot([0, top], [0, top], color="0.4", lw=0.8, ls="--")
    ax.scatter(expected, observed, s=4, color=color, label=label, rasterized=True)
    ax.set_xlabel(r"expected $-\log_{10} p$")
    ax.set_ylabel(r"observed $-\log_{10} p$")
    lam = lambda_gc(p)
    head = f"{title}\n" if title else ""
    ax.set_title(f"{head}" + r"$\lambda_{GC}$" + f" = {lam:.3f}")
    return ax


# --- plink2 concordance ----------------------------------------------------


def read_plink2_glm(path: str | Path) -> pd.DataFrame:
    """Read a plink2 ``--glm`` linear output file, keeping only ADD rows."""
    df = pd.read_csv(path, sep="\t", dtype={"ID": str, "A1": str})
    df = df.rename(columns={"#CHROM": "CHROM"})
    if "TEST" in df.columns:
        df = df[df["TEST"] == "ADD"].copy()
    for col in ("BETA", "SE", "T_STAT", "P"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    missing = {"ID", "A1", "BETA", "P"} - set(df.columns)
    if missing:
        raise ValueError(f"plink2 glm file lacks columns: {sorted(missing)}")
    return df.reset_index(drop=True)


@dataclass
class ConcordanceResult:
    passed: bool
    n_compared: int
    n_flipped: int
    missing: list = field(default_factory=list)
    mismatches: pd.DataFrame = field(default_factory=pd.DataFrame)
    max_abs_beta_diff: float = float("nan")
    max_rel_beta_diff: float = float("nan")
    max_abs_log10p_diff: float = float("nan")
    tolerances: dict = field(default_factory=dict)

    def summary(self) -> str:
        state = "PASS" if self.passed else "FAIL"
        return (
            f"{state}: {self.n_compared} variants compared ({self.n_flipped} allele flips), "
            f"{len(self.missing)} missing, {len(self.mismatches)} beyond tolerance; "
            f"max |dBETA| = {self.max_abs_beta_diff:.3g}, "
            f"max rel dBETA = {self.max_rel_beta_diff:.3g}, "
            f"max |d(-log10 p)| = {self.max_abs_log10p_diff:.3g}"
        )

    def to_dict(self) -> dict:
        return {
            "passed": bool(self.passed),
            "n_compared": int(self.n_compared),
            "n_flipped": int(self.n_flipped),
            "missing": list(self.missing),
            "n_beyond_tolerance": int(len(self.mismatches)),
            "beyond_tolerance_ids": self.mismatches.get("ID", pd.Series(dtype=str)).tolist(),
            "max_abs_beta_diff": float(self.max_abs_beta_diff),
            "max_rel_beta_diff": float(self.max_rel_beta_diff),
            "max_abs_log10p_diff": float(self.max_abs_log10p_diff),
            "tolerances": dict(self.tolerances),
        }


def compare_to_plink2(
    python_results: pd.DataFrame,
    plink_glm,
    beta_rtol: float = 1e-3,
    beta_atol: float = 1e-8,
    log10p_rtol: float = 1e-3,
    log10p_atol: float = 1e-2,
) -> ConcordanceResult:
    """Compare this module's BETA and P with a plink2 ``--glm`` output.

    ``python_results`` needs ``ID``, ``BETA``, ``P`` and, to detect allele
    flips, ``A1`` (the allele whose count was regressed). ``plink_glm`` is a
    path to a plink2 ``.glm.linear`` file or an equivalent DataFrame. Where
    the counted alleles differ the Python beta is sign-flipped before
    comparison. A variant passes when
    ``|dBETA| <= beta_atol + beta_rtol |BETA_plink2|`` and
    ``|d(-log10 P)| <= log10p_atol + log10p_rtol |-log10 P_plink2|``.
    Every Python variant must be present in the plink2 results.
    """
    glm = plink_glm if isinstance(plink_glm, pd.DataFrame) else read_plink2_glm(plink_glm)
    if glm["ID"].duplicated().any():
        raise ValueError("duplicate variant IDs in the plink2 results")
    tol = dict(
        beta_rtol=beta_rtol, beta_atol=beta_atol, log10p_rtol=log10p_rtol, log10p_atol=log10p_atol
    )
    py = python_results.copy()
    cols = ["ID", "BETA", "P"] + (["A1"] if "A1" in glm.columns else [])
    merged = py.merge(glm[cols], on="ID", how="left", suffixes=("_py", "_pl"))
    missing = merged.loc[merged["P_pl"].isna() | merged["BETA_pl"].isna(), "ID"].tolist()
    merged = merged[~merged["ID"].isin(missing)].copy()

    flip = np.zeros(len(merged), dtype=bool)
    if "A1_py" in merged.columns and "A1_pl" in merged.columns:
        flip = (merged["A1_py"].astype(str) != merged["A1_pl"].astype(str)).to_numpy()
    b_py = np.where(flip, -1.0, 1.0) * merged["BETA_py"].to_numpy(dtype=float)
    b_pl = merged["BETA_pl"].to_numpy(dtype=float)
    l_py = -np.log10(merged["P_py"].to_numpy(dtype=float))
    l_pl = -np.log10(merged["P_pl"].to_numpy(dtype=float))

    d_beta = np.abs(b_py - b_pl)
    d_logp = np.abs(l_py - l_pl)
    bad = (d_beta > beta_atol + beta_rtol * np.abs(b_pl)) | (
        d_logp > log10p_atol + log10p_rtol * np.abs(l_pl)
    )
    rel = d_beta / np.maximum(np.abs(b_pl), 1e-300)
    mism = merged.loc[bad, ["ID"]].assign(
        BETA_python=b_py[bad],
        BETA_plink2=b_pl[bad],
        log10p_python=l_py[bad],
        log10p_plink2=l_pl[bad],
    )
    n = len(merged)
    return ConcordanceResult(
        passed=bool(n > 0 and not missing and not bad.any()),
        n_compared=n,
        n_flipped=int(flip.sum()),
        missing=missing,
        mismatches=mism.reset_index(drop=True),
        max_abs_beta_diff=float(d_beta.max()) if n else float("nan"),
        max_rel_beta_diff=float(rel.max()) if n else float("nan"),
        max_abs_log10p_diff=float(d_logp.max()) if n else float("nan"),
        tolerances=tol,
    )
