"""Tests for the reference association implementation."""

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from hbfgwas import assoc
from hbfgwas import simulate as sim


def _rng(s=0):
    return np.random.default_rng(101 + s)


def _toy(n=300, m=40, c=3, seed=0):
    rng = _rng(seed)
    g = rng.binomial(2, rng.uniform(0.1, 0.9, m), size=(n, m)).astype(float)
    cov = rng.standard_normal((n, c))
    y = 0.3 * g[:, 0] + cov @ rng.standard_normal(c) + rng.standard_normal(n)
    return g, y, cov


def _loop_ols(g, y, cov):
    """Per-variant OLS with numpy lstsq: the slow, obviously-correct reference."""
    n = len(y)
    out = []
    for j in range(g.shape[1]):
        x = (
            np.column_stack([np.ones(n), g[:, j], cov])
            if cov is not None
            else np.column_stack([np.ones(n), g[:, j]])
        )
        coef, *_ = np.linalg.lstsq(x, y, rcond=None)
        resid = y - x @ coef
        df = n - x.shape[1]
        sigma2 = resid @ resid / df
        se = np.sqrt(sigma2 * np.linalg.inv(x.T @ x)[1, 1])
        t = coef[1] / se
        out.append((coef[1], se, t, 2 * stats.t.sf(abs(t), df)))
    return np.array(out)


def test_no_covariates_matches_linregress():
    g, y, _ = _toy()
    res = assoc.linear_assoc(g, y)
    for j in range(g.shape[1]):
        lr = stats.linregress(g[:, j], y)
        assert res.beta[j] == pytest.approx(lr.slope, rel=1e-10)
        assert res.se[j] == pytest.approx(lr.stderr, rel=1e-8)
        assert res.p[j] == pytest.approx(lr.pvalue, rel=1e-6)
    assert res.df == len(y) - 2


def test_with_covariates_matches_per_variant_lstsq():
    g, y, cov = _toy()
    res = assoc.linear_assoc(g, y, cov)
    ref = _loop_ols(g, y, cov)
    np.testing.assert_allclose(res.beta, ref[:, 0], rtol=1e-9)
    np.testing.assert_allclose(res.se, ref[:, 1], rtol=1e-8)
    np.testing.assert_allclose(res.t, ref[:, 2], rtol=1e-8)
    np.testing.assert_allclose(res.p, ref[:, 3], rtol=1e-6)
    assert res.df == len(y) - 2 - cov.shape[1]


def test_chunking_does_not_change_results():
    g, y, cov = _toy(m=55)
    a = assoc.linear_assoc(g, y, cov, chunk=7)
    b = assoc.linear_assoc(g, y, cov, chunk=10_000)
    np.testing.assert_allclose(a.beta, b.beta, rtol=1e-12)
    np.testing.assert_allclose(a.p, b.p, rtol=1e-10)


def test_many_phenotypes_at_once_equals_one_at_a_time():
    g, _, cov = _toy()
    ys = _rng(5).standard_normal((g.shape[0], 4))
    res = assoc.linear_assoc(g, ys, cov)
    assert res.beta.shape == (g.shape[1], 4)
    for r in range(4):
        one = assoc.linear_assoc(g, ys[:, r], cov)
        np.testing.assert_allclose(res.beta[:, r], one.beta, rtol=1e-10)
        np.testing.assert_allclose(res.p[:, r], one.p, rtol=1e-8)


def test_monomorphic_variant_gives_nan():
    g, y, cov = _toy()
    g[:, 3] = 1.0
    res = assoc.linear_assoc(g, y, cov)
    assert np.isnan(res.beta[3]) and np.isnan(res.p[3])
    assert np.isfinite(res.beta[2])


def test_to_frame():
    g, y, cov = _toy(m=5)
    df = assoc.linear_assoc(g, y, cov).to_frame(ids=[f"v{i}" for i in range(5)])
    assert list(df.columns[:5]) == ["ID", "BETA", "SE", "T_STAT", "P"]
    assert df["ID"].tolist() == ["v0", "v1", "v2", "v3", "v4"]


def test_bad_inputs_raise():
    g, y, cov = _toy()
    with pytest.raises(ValueError):
        assoc.linear_assoc(g[:-1], y)
    g2 = g.copy()
    g2[0, 0] = np.nan
    with pytest.raises(ValueError):
        assoc.linear_assoc(g2, y)


# --- PCA -----------------------------------------------------------------


def test_pca_matches_direct_eigendecomposition():
    s = sim.simulate_genotypes([60] * 4, n_variants=1500, fst=0.03, rng=_rng(6))
    res = assoc.pca(s.genotypes, n_pcs=5, chunk=400)
    x = s.genotypes.astype(float)
    p = x.mean(axis=0) / 2
    z = (x - 2 * p) / np.sqrt(2 * p * (1 - p))
    w, v = np.linalg.eigh(z @ z.T / z.shape[1])
    np.testing.assert_allclose(res.eigenvalues, w[::-1][:5], rtol=1e-8)
    for k in range(5):  # eigenvectors match up to sign
        assert abs(abs(res.pcs[:, k] @ v[:, -1 - k]) - 1) < 1e-8
    np.testing.assert_allclose(res.pcs.T @ res.pcs, np.eye(5), atol=1e-10)
    assert res.n_variants_used == 1500


def test_pca_drops_monomorphic_variants():
    s = sim.simulate_genotypes([50] * 3, n_variants=500, fst=0.03, rng=_rng(7))
    g = s.genotypes.copy()
    g[:, :10] = 0
    assert assoc.pca(g, n_pcs=3).n_variants_used == 490


def test_pc_adjustment_removes_stratification():
    n_per_pop, fst, h = [120] * 5, 0.05, 0.2
    s = sim.simulate_genotypes(n_per_pop, n_variants=6000, fst=fst, rng=_rng(8))
    ph = sim.simulate_phenotype(
        s.genotypes[:, :1], s.pop, r2=[0.0], shift_var_fraction=h, rng=_rng(9)
    )
    null = s.genotypes[:, 1:]
    raw = assoc.linear_assoc(null, ph.y)
    pcs = assoc.pca(s.genotypes, n_pcs=10).pcs
    adj = assoc.linear_assoc(null, ph.y, pcs)
    assert assoc.lambda_gc(raw.p) > 2.0
    assert 0.9 < assoc.lambda_gc(adj.p) < 1.1


# --- summaries -----------------------------------------------------------


def test_lambda_gc_is_one_for_exact_uniform_quantiles():
    p = (np.arange(1, 100_001) - 0.5) / 100_000
    assert assoc.lambda_gc(p) == pytest.approx(1.0, abs=1e-4)


def test_lambda_gc_recovers_a_scaled_chi_square():
    chi2 = stats.chi2.isf((np.arange(1, 20_001) - 0.5) / 20_000, 1) * 1.3
    assert assoc.lambda_gc(stats.chi2.sf(chi2, 1)) == pytest.approx(1.3, rel=1e-3)


def test_type1_error_and_power():
    p = np.array([0.001, 0.02, 0.04, 0.2, 0.5, 0.06, 0.9, 0.03, 0.7, 0.051])
    assert assoc.type1_error(p, alpha=0.05) == pytest.approx(0.4)
    reps = np.array([[1e-9, 1e-3, 1e-10], [0.5, 1e-8, 1e-12]])  # variants x reps
    np.testing.assert_allclose(assoc.empirical_power(reps, threshold=5e-8), [2 / 3, 2 / 3])


def test_qq_plot_draws(tmp_path):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    p = _rng(10).uniform(size=5000)
    assoc.qq_plot(ax, p, title="null")
    out = tmp_path / "qq.png"
    fig.savefig(out)
    plt.close(fig)
    assert out.stat().st_size > 1000
    assert "lambda" in ax.get_title().lower() or any(
        "lambda" in t.get_text().lower() for t in ax.texts
    )


# --- plink2 concordance --------------------------------------------------

GLM_HEADER = (
    "#CHROM\tPOS\tID\tREF\tALT\tPROVISIONAL_REF?\tA1\tOMITTED\tA1_FREQ\tTEST\tOBS_CT"
    "\tBETA\tSE\tT_STAT\tP\tERRCODE\n"
)


def _write_glm(path, rows):
    with open(path, "w") as fh:
        fh.write(GLM_HEADER)
        for r in rows:
            fh.write(
                f"1\t{r['pos']}\t{r['id']}\tA\tG\tN\t{r['a1']}\t{r['omit']}\t0.3\tADD\t300"
                f"\t{r['beta']:.6g}\t{r['se']:.6g}\t{r['t']:.6g}\t{r['p']:.6g}\t.\n"
            )


def _python_results():
    g, y, cov = _toy(m=6)
    res = assoc.linear_assoc(g, y, cov)
    ids = [f"rs{i}" for i in range(6)]
    df = res.to_frame(ids=ids)
    df["A1"] = "G"
    return df


def _glm_rows(df, flip=()):
    rows = []
    for i, r in df.iterrows():
        sign = -1 if r["ID"] in flip else 1
        rows.append(
            dict(
                pos=1000 + i,
                id=r["ID"],
                a1="A" if r["ID"] in flip else "G",
                omit="G" if r["ID"] in flip else "A",
                beta=sign * r["BETA"],
                se=r["SE"],
                t=sign * r["T_STAT"],
                p=r["P"],
            )
        )
    return rows


def test_read_plink2_glm(tmp_path):
    df = _python_results()
    path = tmp_path / "x.PHENO1.glm.linear"
    _write_glm(path, _glm_rows(df))
    with open(path, "a") as fh:  # a covariate row that must be filtered out
        fh.write("1\t9999\trs0\tA\tG\tN\tG\tA\t0.3\tPC1\t300\t0.1\t0.1\t1\t0.3\t.\n")
    glm = assoc.read_plink2_glm(path)
    assert len(glm) == 6
    assert set(glm.columns) >= {"ID", "A1", "BETA", "SE", "P"}


def test_concordance_passes_on_matching_results(tmp_path):
    df = _python_results()
    path = tmp_path / "x.glm.linear"
    _write_glm(path, _glm_rows(df))
    result = assoc.compare_to_plink2(df, path)
    assert result.passed, result.summary()
    assert result.n_compared == 6


def test_concordance_handles_allele_flips(tmp_path):
    df = _python_results()
    path = tmp_path / "x.glm.linear"
    _write_glm(path, _glm_rows(df, flip={"rs2", "rs4"}))
    result = assoc.compare_to_plink2(df, path)
    assert result.passed, result.summary()
    assert result.n_flipped == 2


def test_concordance_fails_on_a_wrong_beta(tmp_path):
    df = _python_results()
    rows = _glm_rows(df)
    rows[3]["beta"] *= 1.05
    path = tmp_path / "x.glm.linear"
    _write_glm(path, rows)
    result = assoc.compare_to_plink2(df, path)
    assert not result.passed
    assert result.mismatches["ID"].tolist() == ["rs3"]


def test_concordance_fails_on_a_missing_variant(tmp_path):
    df = _python_results()
    path = tmp_path / "x.glm.linear"
    _write_glm(path, _glm_rows(df)[:-1])
    result = assoc.compare_to_plink2(df, path)
    assert not result.passed
    assert result.missing == ["rs5"]


def test_concordance_accepts_a_dataframe():
    df = _python_results()
    glm = pd.DataFrame(_glm_rows(df)).rename(
        columns={"id": "ID", "a1": "A1", "beta": "BETA", "se": "SE", "p": "P"}
    )
    assert assoc.compare_to_plink2(df, glm).passed
