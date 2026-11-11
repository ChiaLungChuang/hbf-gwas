"""EXPLORATORY diagnostic, written AFTER the pre-registered run; not part of it.

The secondary outcome (empirical power over 1000 phenotype redraws) disagreed
with the analytic prediction for HBS1L-MYB and HBB cluster. The redraws keep
ONE genotype realization fixed, so they estimate power conditional on those
genotypes, while the analytic formula averages over genotype realizations.
This script rebuilds the same genotypes from the pre-registered seeds and
computes the conditional expected z of each planted variant after PC
adjustment, given the fixed effects of the other planted variants.

    python runs/2026-11-10-sim-stratification/diagnostic_conditional_power.py
"""

import numpy as np
from scipy import stats

from hbfgwas import assoc, power, simulate
from hbfgwas.config import load_config

cfg = load_config()
exp, s = cfg["simulated_experiment"], cfg["seeds"]
rng = np.random.default_rng(s["sim_genotypes"])
null = simulate.simulate_genotypes(
    exp["n_per_subpop"],
    fst=exp["fst"],
    rng=rng,
    n_variants=exp["n_null_variants"],
    freq_range=tuple(exp["ancestral_freq_range"]),
)
pl = simulate.simulate_genotypes(
    exp["n_per_subpop"],
    fst=exp["fst"],
    rng=rng,
    ancestral_freq=np.array([p["ancestral_freq"] for p in exp["planted"]]),
)
pcs = assoc.pca(null.genotypes, exp["n_pcs"]).pcs
G = pl.genotypes.astype(float)
n, h = G.shape[0], exp["shift_var_fraction"]
r2 = np.array([p["r2"] for p in exp["planted"]])
beta = np.sqrt(r2 / G.var(axis=0))
print("sample correlation between planted genotypes:")
print(np.round(np.corrcoef(G.T), 3))
Q, _ = np.linalg.qr(np.column_stack([np.ones(n), pcs]))
Gr = G - Q @ (Q.T @ G)
SS, SStot = (Gr**2).sum(0), ((G - G.mean(0)) ** 2).sum(0)
gen = G @ beta
genr = gen - Q @ (Q.T @ gen)
c = stats.norm.isf(2.5e-8)
for j, p in enumerate(exp["planted"]):
    mean_b = (Gr[:, j] @ genr) / SS[j]
    mu = mean_b / np.sqrt((1 - h - r2[j]) / SS[j])
    cond = stats.norm.sf(c - mu) + stats.norm.cdf(-c - mu)
    print(
        f"{p['locus']:>12}: genotype variance kept after PCs {SS[j] / SStot[j]:.3f}; "
        f"E[z | genotypes] {mu:.3f} vs unconditional sqrt(ncp) "
        f"{np.sqrt(power.ncp(n, r2[j], h)):.3f}; conditional power {cond:.3f}; "
        f"unconditional {power.power(n, 5e-8, r2=r2[j], covariate_r2=h):.3f}"
    )
