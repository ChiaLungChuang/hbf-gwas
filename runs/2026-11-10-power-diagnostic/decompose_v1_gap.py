"""Post hoc (after the Phase 1 result): why the Phase 1 conditional diagnostic
over-predicted power for HBS1L-MYB and HBB cluster.

Uses the Phase 1 genotypes (seed sim_genotypes) and the Phase 1 realised
betas, and draws NO new phenotypes on them: every number below is analytic and
is compared with the empirical power already recorded in
runs/2026-11-10-sim-stratification/results.json (1000 redraws).

    python runs/2026-11-10-power-diagnostic/decompose_v1_gap.py
"""

import json
from pathlib import Path

import numpy as np
from scipy import stats

from hbfgwas import assoc, predict, simulate
from hbfgwas.config import load_config

V1 = Path("runs/2026-11-10-sim-stratification/results.json")
cfg = load_config()
exp, seeds = cfg["simulated_experiment"], cfg["seeds"]
alpha = cfg["association"]["genome_wide_alpha"]
v1 = {loc["locus"]: loc for loc in json.loads(V1.read_text())["loci"]}

rng = np.random.default_rng(seeds["sim_genotypes"])
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
names = [p["locus"] for p in exp["planted"]]
betas = np.array([v1[k]["beta_planted"] for k in names])  # Phase 1 used r2 mode
r2 = np.array([p["r2"] for p in exp["planted"]])

Q, _ = np.linalg.qr(np.column_stack([np.ones(n), pcs]))
Gr = G - Q @ (Q.T @ G)
SS = (Gr**2).sum(0)
gen = G @ betas
genr = gen - Q @ (Q.T @ gen)
df = n - Q.shape[1] - 1
z_crit = stats.norm.isf(alpha / 2)
t_crit = stats.t.isf(alpha / 2, df)

full = predict.conditional_power(
    G,
    betas,
    alpha,
    covariates=pcs,
    pop=null.pop,
    shift_var_fraction=h,
    n_shift_draws=20000,
    rng=np.random.default_rng(1),
)
print(f"n = {n}, df = {df}, normal critical value {z_crit:.4f}, t critical value {t_crit:.4f}\n")
rows = []
for j, name in enumerate(names):
    m = Gr[:, j] @ genr / SS[j]
    mu_old = m / np.sqrt((1 - h - r2[j]) / SS[j])  # Phase 1 diagnostic: fixed loci treated as noise
    a = stats.norm.sf(z_crit - mu_old) + stats.norm.cdf(-z_crit - mu_old)
    b = stats.norm.sf(t_crit - mu_old) + stats.norm.cdf(-t_crit - mu_old)
    emp = v1[name]["empirical_power"]
    lo, hi = v1[name]["empirical_power_ci95"]
    se = np.sqrt(emp * (1 - emp) / 1000)
    c = full[j]["power"]
    rows.append((name, a, b, c, emp, lo, hi, se))
    print(
        f"{name:>12}: (a) Phase 1 diagnostic {a:.3f} | (b) + t critical value {b:.3f} | "
        f"(c) + fixed loci raise the critical value, not the spread {c:.3f} | "
        f"recorded empirical {emp:.3f} [{lo:.3f}, {hi:.3f}]; (c) - empirical = "
        f"{(c - emp):+.3f} ({(emp - c) / se:+.2f} SE)"
    )
