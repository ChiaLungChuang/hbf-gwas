"""Post hoc audit for issue #14: why HBS1L-MYB's predicted power is 0.809 in sim-v2.

Analytic only, on the genotype draw shared by sim-v1 and sim-v2 (seed
sim_genotypes); no phenotype is drawn. Compares with the committed records in
runs/2026-11-10-sim-stratification and runs/2026-11-10-sim-v2.

    python runs/2026-11-11-hbs1l-myb-audit/audit.py
"""

import json
from pathlib import Path

import numpy as np

from hbfgwas import assoc, power, predict, simulate
from hbfgwas.config import load_config

cfg = load_config()
exp, seeds = cfg["simulated_experiment"], cfg["seeds"]
v1 = json.loads(Path("runs/2026-11-10-sim-stratification/results.json").read_text())
v2 = json.loads(Path("runs/2026-11-10-sim-v2/results.json").read_text())
loci_v1 = {x["locus"]: x for x in v1["loci"]}
loci_v2 = {x["locus"]: x for x in v2["loci"]}
names = [p["locus"] for p in exp["planted"]]
J = names.index("HBS1L-MYB")

rng = np.random.default_rng(seeds["sim_genotypes"])  # the sim-v1 / sim-v2 genotype seed
null = simulate.simulate_genotypes(
    exp["n_per_subpop"],
    fst=exp["fst"],
    rng=rng,
    n_variants=exp["n_null_variants"],
    freq_range=tuple(exp["ancestral_freq_range"]),
)
planted = simulate.simulate_genotypes(
    exp["n_per_subpop"],
    fst=exp["fst"],
    rng=rng,
    ancestral_freq=np.array([p["ancestral_freq"] for p in exp["planted"]]),
)
pcs = assoc.pca(null.genotypes, exp["n_pcs"]).pcs
G = planted.genotypes.astype(float)

print("sample allele frequencies (v1 record, v2 record, recomputed):")
print(
    " ",
    [round(loci_v1[n]["sample_freq"], 4) for n in names],
    [round(loci_v2[n]["sample_freq"], 4) for n in names],
    list(np.round(G.mean(0) / 2, 4)),
)
var_g = G.var(0)
p_anc = exp["planted"][J]["ancestral_freq"]
print(
    f"HBS1L-MYB genotype variance in the sample {var_g[J]:.5f}; "
    f"2p(1-p) at the ancestral frequency {p_anc}: {2 * p_anc * (1 - p_anc):.5f}"
)

b_v1 = np.array([loci_v1[n]["beta_planted"] for n in names])
b_v2 = np.array([p["beta_sd"] for p in exp["planted"]])
print(
    f"HBS1L-MYB beta v1 {b_v1[J]:.4f} (r2 {b_v1[J] ** 2 * var_g[J]:.5f}); "
    f"v2 {b_v2[J]:.4f} (r2 {b_v2[J] ** 2 * var_g[J]:.5f})"
)
print("phenotype variance v1", v1["phenotype_variance"], "v2", v2["phenotype_variance"])


def conditional(betas):
    return predict.conditional_power(
        G,
        betas,
        5e-8,
        covariates=pcs,
        pop=null.pop,
        shift_var_fraction=0.2,
        n_shift_draws=2000,
        rng=np.random.default_rng(seeds["sim_prediction"]),
    )[J]


b_mixed = b_v1.copy()
b_mixed[J] = b_v2[J]
print(
    f"\nv1 declared (unconditional chi-square, r2 0.026, partial R2): "
    f"{power.power(650, 5e-8, r2=0.026, covariate_r2=0.2):.4f}"
)
for label, betas in (
    ("v1 betas, corrected conditional", b_v1),
    ("v1 betas, HBS1L-MYB beta = 0.571 only", b_mixed),
    ("v2 betas (all beta_sd)", b_v2),
):
    r = conditional(betas)
    print(
        f"{label:40}: power {r['power']:.4f}; E[beta_hat] {r['expected_beta']:.4f} "
        f"= beta {r['beta']:.4f} + leakage {r['leakage']:+.4f}; "
        f"noncentrality {r['mean_noncentrality']:.3f} vs effective critical "
        f"{r['mean_effective_critical_t']:.3f}; noise {r['noise_var']:.4f}"
    )
print(
    f"\nv2 provenance: commit {v2['provenance']['git_commit'][:7]}, "
    f"run_utc {v2['provenance']['run_utc']}"
)
