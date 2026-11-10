"""Genotype simulation with population structure, and phenotype simulation.

Genotypes follow the Balding-Nichols model: each variant has an ancestral
allele frequency ``p``; subpopulation ``k`` has frequency
``p_k ~ Beta(p (1 - F) / F, (1 - p) (1 - F) / F)`` (mean ``p``, variance
``F p (1 - p)``), and individuals draw ``Binomial(2, p_k)`` allele counts, so
Hardy-Weinberg proportions hold inside each subpopulation.

The phenotype has three separable components on a scale whose total variance
is 1 (in expectation):

* ``genetic``: planted additive effects of the given variants,
* ``shift``: a non-genetic mean shift per subpopulation (the confounder), and
* ``noise``: independent normal noise.

``simulate_phenotype`` accepts any genotype matrix, simulated here or real
(e.g. read from a plink2 ``--export A`` file), with any subpopulation labels.
The phenotype is always simulated; it is never a measured trait.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

_CHUNK = 20_000  # variants generated per block; fixed so the RNG stream is reproducible


@dataclass
class SimulatedGenotypes:
    genotypes: np.ndarray  # (n, m) int8 counts of the counted allele
    pop: np.ndarray  # (n,) subpopulation index
    ancestral_freq: np.ndarray  # (m,)
    subpop_freq: np.ndarray  # (K, m)
    fst: float


@dataclass
class Phenotype:
    y: np.ndarray  # (n,) or (n, reps)
    genetic: np.ndarray
    shift: np.ndarray
    noise: np.ndarray
    betas: np.ndarray  # per-allele effect of each planted variant, phenotype units
    r2_realized: np.ndarray  # variance of beta_j * g_j in this sample
    pop_shift: dict = field(default_factory=dict)  # label -> shift (array if reps)
    shift_var_fraction: float = 0.0


def pop_labels(n_per_pop: Sequence[int]) -> np.ndarray:
    """Subpopulation index for each individual: 0 repeated n_0 times, then 1, ..."""
    return np.repeat(np.arange(len(n_per_pop)), n_per_pop)


def balding_nichols_freqs(
    ancestral_freq: np.ndarray, fst: float, n_pops: int, rng: np.random.Generator
) -> np.ndarray:
    """Draw subpopulation allele frequencies, shape (n_pops, m)."""
    p = np.asarray(ancestral_freq, dtype=float)
    if np.any((p <= 0) | (p >= 1)):
        raise ValueError("ancestral frequencies must be in (0, 1)")
    if not 0 <= fst < 1:
        raise ValueError("fst must be in [0, 1)")
    if fst == 0:
        return np.tile(p, (n_pops, 1))
    scale = (1 - fst) / fst
    return rng.beta(p * scale, (1 - p) * scale, size=(n_pops, p.size))


def simulate_genotypes(
    n_per_pop: Sequence[int],
    fst: float,
    rng: np.random.Generator,
    n_variants: int | None = None,
    ancestral_freq: np.ndarray | None = None,
    freq_range: tuple[float, float] = (0.05, 0.95),
) -> SimulatedGenotypes:
    """Balding-Nichols genotypes for subpopulations of the given sizes.

    Ancestral frequencies are ``ancestral_freq`` if given, otherwise
    ``Uniform(freq_range)`` draws for ``n_variants`` variants.
    """
    if ancestral_freq is None:
        if n_variants is None:
            raise ValueError("give n_variants or ancestral_freq")
        lo, hi = freq_range
        ancestral_freq = rng.uniform(lo, hi, size=n_variants)
    p = np.asarray(ancestral_freq, dtype=float)
    sizes = [int(s) for s in n_per_pop]
    subpop_freq = balding_nichols_freqs(p, fst, len(sizes), rng)
    geno = np.empty((sum(sizes), p.size), dtype=np.int8)
    starts = np.concatenate([[0], np.cumsum(sizes)])
    for c0 in range(0, p.size, _CHUNK):
        c1 = min(c0 + _CHUNK, p.size)
        for k, n_k in enumerate(sizes):
            geno[starts[k] : starts[k + 1], c0:c1] = rng.binomial(
                2, subpop_freq[k, c0:c1], size=(n_k, c1 - c0)
            )
    return SimulatedGenotypes(geno, pop_labels(sizes), p, subpop_freq, float(fst))


def estimate_fst(genotypes: np.ndarray, pop: np.ndarray) -> float:
    """Moment estimate of Balding-Nichols F (ratio of averages over variants).

    Between-subpopulation variance of sample frequencies, minus its expected
    binomial sampling part, divided by an estimate of ``p (1 - p)``. Assumes
    Hardy-Weinberg proportions within subpopulations.
    """
    pop = np.asarray(pop)
    labels = np.unique(pop)
    k = labels.size
    if k < 2:
        raise ValueError("need at least two subpopulations")
    g = np.asarray(genotypes, dtype=float)
    freqs = np.vstack([g[pop == lab].mean(axis=0) / 2 for lab in labels])  # (K, m)
    n_alleles = np.array([2 * np.sum(pop == lab) for lab in labels], dtype=float)[:, None]
    s2 = freqs.var(axis=0, ddof=1)
    sampling = np.mean(freqs * (1 - freqs) / (n_alleles - 1), axis=0)
    pbar = freqs.mean(axis=0)
    num = s2 - sampling
    den = pbar * (1 - pbar) + s2 / k
    return float(num.sum() / den.sum())


def expected_null_chi2(
    n: int,
    fst: float,
    shift_var_fraction: float,
    n_pops: int,
) -> float:
    """Expected 1-df chi-square of a null variant without ancestry adjustment.

    For equal-sized subpopulations, a subpopulation shift explaining a fraction
    ``h`` of phenotypic variance, and Balding-Nichols differentiation ``F``:

        E[chi2] ~= 1 + 2 n F h / (K (1 + F (1 - 2/K)))

    Derivation: the sample covariance of genotype and shift is
    ``2 sum_k w_k (p_k - p)(s_k - sbar)``, whose variance over variants is
    ``4 F p (1 - p) sum_k w_k^2 (s_k - sbar)^2 = 4 F p (1 - p) h / K`` for
    ``w_k = 1/K``; the genotype variance is ``2 p (1 - p) (1 + F (1 - 2/K))``.
    Because chi2 is roughly ``(Z + delta)^2`` with ``delta`` independent of
    ``Z``, the same number approximates lambda_GC.
    """
    k = n_pops
    return 1.0 + 2.0 * n * fst * shift_var_fraction / (k * (1.0 + fst * (1.0 - 2.0 / k)))


def scale_shift_draws(draws: np.ndarray, idx: np.ndarray, shift_var_fraction: float):
    """Turn N(0, 1) subpopulation draws (reps, K) into individual-level shifts.

    Each repetition is centred over individuals and scaled so its variance is
    exactly ``shift_var_fraction``. Returns (shift (n, reps), pop_means (K, reps)).
    Shared by ``simulate_phenotype`` and ``hbfgwas.predict`` so both use one law.
    """
    reps = draws.shape[0]
    shift = np.zeros((idx.size, reps))
    pop_means = np.zeros((draws.shape[1], reps))
    for r in range(reps):
        s = draws[r][idx]
        centre = s.mean()
        scale = np.sqrt(shift_var_fraction / (s - centre).var())
        shift[:, r] = (s - centre) * scale
        pop_means[:, r] = (draws[r] - centre) * scale
    return shift, pop_means


def _validate_genotypes(genotypes) -> np.ndarray:
    g = np.asarray(genotypes, dtype=float)
    if g.ndim == 1:
        g = g[:, None]
    if g.ndim != 2:
        raise ValueError("genotypes must be (n, k)")
    if np.isnan(g).any():
        raise ValueError("genotypes contain missing values; impute or drop them first")
    return g


def simulate_phenotype(
    genotypes,
    pop,
    rng: np.random.Generator,
    r2: Sequence[float] | None = None,
    beta: Sequence[float] | None = None,
    shift_var_fraction: float = 0.0,
    n_reps: int | None = None,
) -> Phenotype:
    """Simulate y = genetic + shift + noise.

    Parameters
    ----------
    genotypes : (n, k) allele counts of the planted variants (simulated or real)
    pop : (n,) subpopulation labels of any type
    r2 : target variance explained per planted variant; beta_j is set to
        ``sqrt(r2_j / Var(g_j))`` using the sample genotype variance, so the
        realized variance of each planted term equals ``r2_j`` exactly
    beta : per-allele effects in phenotype units (alternative to ``r2``)
    shift_var_fraction : variance of the subpopulation shift; subpopulation
        means are drawn N(0, 1), centred, and scaled to this variance exactly
    n_reps : if given, draw this many phenotypes at once (shift and noise are
        redrawn; betas are fixed by the genotypes); outputs are (n, n_reps)
    """
    g = _validate_genotypes(genotypes)
    n, k = g.shape
    pop = np.asarray(pop)
    if pop.shape != (n,):
        raise ValueError("pop must have one label per individual")
    if (r2 is None) == (beta is None):
        raise ValueError("give exactly one of r2 or beta")
    if not 0 <= shift_var_fraction < 1:
        raise ValueError("shift_var_fraction must be in [0, 1)")

    gvar = g.var(axis=0)
    if r2 is not None:
        r2_arr = np.asarray(r2, dtype=float)
        if r2_arr.shape != (k,):
            raise ValueError("r2 needs one value per planted variant")
        if np.any(r2_arr < 0):
            raise ValueError("r2 must be non-negative")
        if np.any((gvar == 0) & (r2_arr > 0)):
            raise ValueError("a planted variant is monomorphic in this sample")
        ratio = np.divide(r2_arr, gvar, out=np.zeros_like(r2_arr), where=gvar > 0)
        betas = np.sqrt(ratio)
    else:
        betas = np.asarray(beta, dtype=float)
        if betas.shape != (k,):
            raise ValueError("beta needs one value per planted variant")

    raw = g @ betas
    genetic = raw - raw.mean()
    genetic_var = genetic.var()
    noise_var = 1.0 - genetic_var - shift_var_fraction
    if noise_var <= 0:
        raise ValueError("planted and shift variance fractions must sum to less than 1")

    labels, idx = np.unique(pop, return_inverse=True)
    if shift_var_fraction > 0 and labels.size < 2:
        raise ValueError("a subpopulation shift needs at least two subpopulations")

    reps = 1 if n_reps is None else int(n_reps)
    shift = np.zeros((n, reps))
    pop_means = np.zeros((labels.size, reps))
    if shift_var_fraction > 0:
        draws = rng.standard_normal((reps, labels.size))
        shift, pop_means = scale_shift_draws(draws, idx, shift_var_fraction)
    noise = rng.standard_normal((n, reps)) * np.sqrt(noise_var)

    r2_realized = (betas[None, :] * g).var(axis=0)
    if n_reps is None:
        shift, noise = shift[:, 0], noise[:, 0]
        y = genetic + shift + noise
        pop_shift = {lab: float(pop_means[i, 0]) for i, lab in enumerate(labels.tolist())}
        genetic_out = genetic
    else:
        genetic_out = np.broadcast_to(genetic[:, None], (n, reps))
        y = genetic_out + shift + noise
        pop_shift = {lab: pop_means[i].copy() for i, lab in enumerate(labels.tolist())}
    return Phenotype(
        y=y,
        genetic=genetic_out,
        shift=shift,
        noise=noise,
        betas=betas,
        r2_realized=r2_realized,
        pop_shift=pop_shift,
        shift_var_fraction=float(shift_var_fraction),
    )
