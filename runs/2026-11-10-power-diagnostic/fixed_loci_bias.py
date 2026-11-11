"""Post hoc, on its own scenario (seeds 404-406, 1000+, 2000+; not the
pre-registered seeds): the second cause found while building the conditional
power formula.

With genotypes held fixed, the other planted loci contribute a FIXED part to
the residual. It inflates the t-test's estimated standard error but adds no
sampling spread to beta_hat. Treating it as random noise (scaling the
noncentrality, "old") over-predicts power for weak loci; raising the critical
value instead ("new", hbfgwas.predict.conditional_power) does not.
40,000 phenotype redraws per scenario.

    python runs/2026-11-10-power-diagnostic/fixed_loci_bias.py
"""

import numpy as np

from hbfgwas import assoc, predict
from hbfgwas import simulate as sim

rng = np.random.default_rng(404)
sizes = [150] * 3
null = sim.simulate_genotypes(sizes, fst=0.05, rng=rng, n_variants=3000)
pl = sim.simulate_genotypes(sizes, fst=0.05, rng=rng, ancestral_freq=np.array([0.3, 0.1, 0.2]))
pcs = assoc.pca(null.genotypes, 5).pcs
G = pl.genotypes.astype(float)
BETAS = np.array([0.25, 0.30, 0.20])
ALPHA, R = 1e-4, 40_000


def old_formula(g, b, cov, h):
    """Fixed loci folded into the noise: noncentrality scaled down, critical value unchanged."""
    q = assoc._covariate_basis(g.shape[0], cov)
    df = g.shape[0] - q.shape[1] - 1
    from scipy import stats

    crit = stats.t.isf(ALPHA / 2, df)
    gr = g - q @ (q.T @ g)
    gen = g @ b
    genr = gen - q @ (q.T @ gen)
    nv = predict.noise_variance(g, b, h)
    out = []
    for j in range(g.shape[1]):
        ss = gr[:, j] @ gr[:, j]
        m = gr[:, j] @ genr / ss
        lam = (genr @ genr - m * m * ss) / nv
        d = m / np.sqrt(nv / ss) / np.sqrt(1 + lam / df)
        out.append(float(predict.two_sided_power(crit, df, d)))
    return out


def run(label, g, b, h, cov, seed0):
    new = predict.conditional_power(
        g,
        b,
        ALPHA,
        covariates=cov,
        pop=null.pop,
        shift_var_fraction=h,
        n_shift_draws=20000,
        rng=np.random.default_rng(7),
    )
    old = old_formula(g, b, cov, h) if h == 0 else [float("nan")] * g.shape[1]
    hits, done = np.zeros(g.shape[1]), 0
    for blk in range(R // 4000):
        ph = sim.simulate_phenotype(
            g,
            null.pop,
            rng=np.random.default_rng(seed0 + blk),
            beta=b,
            shift_var_fraction=h,
            n_reps=4000,
        )
        hits += (assoc.linear_assoc(g, ph.y, cov).p < ALPHA).sum(1)
        done += 4000
    for j in range(g.shape[1]):
        e = hits[j] / done
        p = new[j]["power"]
        se = np.hypot(np.sqrt(p * (1 - p) / done), new[j]["power_mc_se"])
        print(
            f"{label:30} locus {j}: empirical {e:.4f} | new {p:.4f} ({(e - p) / se:+.2f} SE)"
            f" | old {old[j]:.4f}"
        )


run("no shift, 5 PCs", G, BETAS, 0.0, pcs, 2000)
run("no shift, no covariates", G, BETAS, 0.0, None, 2000)
run("no shift, locus 1 alone", G[:, [1]], BETAS[[1]], 0.0, None, 2000)
run("shift h=0.15, 5 PCs", G, BETAS, 0.15, pcs, 1000)
run("shift h=0.15, no covariates", G, BETAS, 0.15, None, 1000)
