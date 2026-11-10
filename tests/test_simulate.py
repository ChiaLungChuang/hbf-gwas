"""Tests for genotype and phenotype simulation.

The calibration tests fail if the simulation stops producing population
structure: estimated Fst must track the configured Fst, the top principal
components must carry the subpopulation signal, and a subpopulation-correlated
shift must inflate null test statistics (and must not when Fst = 0).
These tests use their own seeds, not the pre-registered experiment seeds.
"""

import numpy as np
import pytest

from hbfgwas import simulate as sim

SEED = 7


def _rng(offset=0):
    return np.random.default_rng(SEED + offset)


def _null_chi2(genotypes, y):
    """Marginal 1-df chi-square (n r^2) for each column, computed directly."""
    g = genotypes.astype(float)
    g = g - g.mean(axis=0)
    yc = y - y.mean()
    r = (g * yc[:, None]).sum(axis=0) / np.sqrt((g * g).sum(axis=0) * (yc @ yc))
    return len(y) * r**2


# --- genotypes ---------------------------------------------------------------


def test_pop_labels():
    pop = sim.pop_labels([3, 2, 4])
    assert pop.tolist() == [0, 0, 0, 1, 1, 2, 2, 2, 2]


def test_genotype_shapes_and_values():
    g = sim.simulate_genotypes([50, 60, 70], n_variants=500, fst=0.05, rng=_rng())
    assert g.genotypes.shape == (180, 500)
    assert g.genotypes.dtype == np.int8
    assert set(np.unique(g.genotypes)) <= {0, 1, 2}
    assert g.pop.shape == (180,)
    assert g.subpop_freq.shape == (3, 500)
    assert g.ancestral_freq.shape == (500,)


def test_same_seed_same_genotypes():
    a = sim.simulate_genotypes([40, 40], n_variants=300, fst=0.05, rng=_rng())
    b = sim.simulate_genotypes([40, 40], n_variants=300, fst=0.05, rng=_rng())
    np.testing.assert_array_equal(a.genotypes, b.genotypes)


def test_given_ancestral_freq_is_used():
    p = np.array([0.1, 0.5, 0.9])
    g = sim.simulate_genotypes([2000, 2000], ancestral_freq=p, fst=0.0, rng=_rng())
    np.testing.assert_allclose(g.genotypes.mean(axis=0) / 2, p, atol=0.02)


def test_fst_zero_means_identical_subpop_freqs():
    g = sim.simulate_genotypes([30, 30, 30], n_variants=200, fst=0.0, rng=_rng())
    np.testing.assert_array_equal(g.subpop_freq[0], g.subpop_freq[1])


def test_balding_nichols_moments():
    # E[p_k] = p and Var[p_k] = F p (1 - p)
    p = np.full(200_000, 0.3)
    f = sim.balding_nichols_freqs(p, fst=0.02, n_pops=1, rng=_rng())[0]
    assert f.mean() == pytest.approx(0.3, abs=2e-3)
    assert f.var() == pytest.approx(0.02 * 0.3 * 0.7, rel=0.03)


@pytest.mark.parametrize("fst", [0.01, 0.05])
def test_calibration_fst_estimate_tracks_configured_fst(fst):
    g = sim.simulate_genotypes([150] * 6, n_variants=8000, fst=fst, rng=_rng(1))
    est = sim.estimate_fst(g.genotypes, g.pop)
    assert est == pytest.approx(fst, rel=0.15), est


def test_calibration_no_structure_when_fst_zero():
    g = sim.simulate_genotypes([150] * 6, n_variants=8000, fst=0.0, rng=_rng(2))
    assert abs(sim.estimate_fst(g.genotypes, g.pop)) < 0.001


def test_calibration_top_pcs_carry_the_subpopulations():
    k = 5
    g = sim.simulate_genotypes([80] * k, n_variants=4000, fst=0.03, rng=_rng(3))
    x = g.genotypes.astype(float)
    p = x.mean(axis=0) / 2
    z = (x - 2 * p) / np.sqrt(2 * p * (1 - p))
    eig = np.sort(np.linalg.eigvalsh(z @ z.T / z.shape[1]))[::-1]
    # k subpopulations give k - 1 structured axes, then a clear drop
    assert eig[k - 2] / eig[k - 1] > 1.5
    # and the first PC separates subpopulations: between-group variance dominates
    u = np.linalg.eigh(z @ z.T)[1][:, -1]
    groups = [u[g.pop == j] for j in range(k)]
    between = np.var([grp.mean() for grp in groups])
    within = np.mean([grp.var() for grp in groups])
    assert between > 5 * within


# --- phenotype ---------------------------------------------------------------


def _planted(n_per_pop=(200, 200, 200), fst=0.02, offset=0):
    return sim.simulate_genotypes(
        list(n_per_pop), ancestral_freq=np.array([0.3, 0.1, 0.5]), fst=fst, rng=_rng(offset)
    )


def test_phenotype_components_sum_and_have_target_variances():
    g = _planted()
    ph = sim.simulate_phenotype(
        g.genotypes, g.pop, r2=[0.05, 0.02, 0.01], shift_var_fraction=0.2, rng=_rng(10)
    )
    np.testing.assert_allclose(ph.y, ph.genetic + ph.shift + ph.noise)
    # planted variance is exact by construction (beta from the sample genotype variance)
    per_variant = ph.betas[None, :] * g.genotypes
    np.testing.assert_allclose(per_variant.var(axis=0), [0.05, 0.02, 0.01], rtol=1e-10)
    # the shift variance is exact by construction; noise is random
    assert ph.shift.var() == pytest.approx(0.2, rel=1e-10)
    assert ph.noise.var() == pytest.approx(1 - 0.08 - 0.2, rel=0.15)
    assert ph.y.var() == pytest.approx(1.0, rel=0.15)


def test_shift_is_constant_within_subpopulation():
    g = _planted()
    ph = sim.simulate_phenotype(
        g.genotypes, g.pop, r2=[0.05, 0.0, 0.0], shift_var_fraction=0.3, rng=_rng(11)
    )
    for k in np.unique(g.pop):
        assert np.ptp(ph.shift[g.pop == k]) == pytest.approx(0.0, abs=1e-12)
    assert len(ph.pop_shift) == 3


def test_noise_is_not_subpopulation_correlated():
    g = _planted(n_per_pop=(2000, 2000, 2000))
    ph = sim.simulate_phenotype(
        g.genotypes, g.pop, r2=[0.0, 0.0, 0.0], shift_var_fraction=0.0, rng=_rng(12)
    )
    means = [ph.noise[g.pop == k].mean() for k in range(3)]
    assert np.ptp(means) < 0.1


def test_beta_mode_uses_given_effects():
    g = _planted()
    ph = sim.simulate_phenotype(g.genotypes, g.pop, beta=[0.3, 0.0, -0.1], rng=_rng(13))
    np.testing.assert_allclose(ph.betas, [0.3, 0.0, -0.1])
    np.testing.assert_allclose(
        ph.genetic,
        g.genotypes @ np.array([0.3, 0.0, -0.1])
        - (g.genotypes @ np.array([0.3, 0.0, -0.1])).mean(),
    )


def test_accepts_a_real_genotype_matrix_with_string_populations():
    # stands in for genotypes read from a plink2 --export A file
    rng = _rng(14)
    real = rng.integers(0, 3, size=(120, 2)).astype(float)
    pops = np.array(["YRI", "LWK", "GWD", "MSL", "ESN", "ASW"] * 20)
    ph = sim.simulate_phenotype(real, pops, r2=[0.04, 0.02], shift_var_fraction=0.1, rng=rng)
    assert ph.y.shape == (120,)
    assert set(ph.pop_shift) == {"YRI", "LWK", "GWD", "MSL", "ESN", "ASW"}
    np.testing.assert_allclose(ph.y, ph.genetic + ph.shift + ph.noise)


def test_many_repetitions_at_once():
    g = _planted()
    ph = sim.simulate_phenotype(
        g.genotypes, g.pop, r2=[0.05, 0.02, 0.01], shift_var_fraction=0.2, rng=_rng(15), n_reps=50
    )
    assert ph.y.shape == (600, 50)
    assert ph.shift.shape == (600, 50)
    # betas are fixed by the genotypes; shift and noise are redrawn per repetition
    assert not np.allclose(ph.noise[:, 0], ph.noise[:, 1])
    np.testing.assert_allclose(ph.shift.var(axis=0), 0.2, rtol=1e-10)


def test_same_seed_same_phenotype():
    g = _planted()
    kw = dict(r2=[0.05, 0.02, 0.01], shift_var_fraction=0.2)
    a = sim.simulate_phenotype(g.genotypes, g.pop, rng=_rng(16), **kw)
    b = sim.simulate_phenotype(g.genotypes, g.pop, rng=_rng(16), **kw)
    np.testing.assert_array_equal(a.y, b.y)


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(r2=[0.5, 0.3, 0.1], shift_var_fraction=0.2),  # fractions sum >= 1
        dict(r2=[0.05, 0.02], shift_var_fraction=0.0),  # wrong length
        dict(r2=[0.05, 0.02, 0.01], beta=[0.1, 0.1, 0.1]),  # both
        dict(),  # neither
        dict(r2=[0.05, 0.02, 0.01], shift_var_fraction=-0.1),
    ],
)
def test_bad_phenotype_inputs_raise(kwargs):
    g = _planted()
    with pytest.raises(ValueError):
        sim.simulate_phenotype(g.genotypes, g.pop, rng=_rng(17), **kwargs)


def test_monomorphic_planted_variant_raises_in_r2_mode():
    geno = np.zeros((10, 1))
    with pytest.raises(ValueError):
        sim.simulate_phenotype(geno, np.zeros(10, int), r2=[0.05], rng=_rng(18))


def test_missing_genotypes_raise():
    geno = np.ones((10, 1))
    geno[0, 0] = np.nan
    with pytest.raises(ValueError):
        sim.simulate_phenotype(geno, np.zeros(10, int), beta=[0.1], rng=_rng(19))


# --- calibration: structure plus shift inflates null statistics ---------------


def test_calibration_shift_with_structure_inflates_null_tests():
    n_per_pop, fst, h = [200] * 5, 0.05, 0.2
    g = sim.simulate_genotypes(n_per_pop, n_variants=3000, fst=fst, rng=_rng(20))
    ph = sim.simulate_phenotype(
        g.genotypes[:, :1], g.pop, r2=[0.0], shift_var_fraction=h, rng=_rng(21)
    )
    chi2 = _null_chi2(g.genotypes[:, 1:], ph.y)
    predicted = sim.expected_null_chi2(
        n=sum(n_per_pop), fst=fst, shift_var_fraction=h, n_pops=len(n_per_pop)
    )
    assert predicted > 3
    assert chi2.mean() == pytest.approx(predicted, rel=0.2)


def test_calibration_no_inflation_without_structure():
    g = sim.simulate_genotypes([200] * 5, n_variants=3000, fst=0.0, rng=_rng(22))
    ph = sim.simulate_phenotype(
        g.genotypes[:, :1], g.pop, r2=[0.0], shift_var_fraction=0.2, rng=_rng(23)
    )
    chi2 = _null_chi2(g.genotypes[:, 1:], ph.y)
    assert chi2.mean() == pytest.approx(1.0, abs=0.08)
