# 2026-11-10: why predicted power disagreed in Phase 1 (post hoc)

**Post hoc.** Everything here was done after the Phase 1 result, at the
user's request (issue #10). No phenotype was drawn on the Phase 1
genotypes: `decompose_v1_gap.py` is analytic and compares with the empirical
power already recorded in `runs/2026-11-10-sim-stratification/results.json`.
`fixed_loci_bias.py` uses its own scenario and seeds.

## Question

Phase 1 recorded empirical power (1000 redraws, genotypes fixed) of 0.375
for HBS1L-MYB and 0.419 for HBB cluster. The Phase 1 diagnostic's
"conditional" predictions, 0.406 and 0.461, were still 2-3 SE too high.

## Causes

1. **Normal instead of t critical value** (`decompose_v1_gap.txt`). The test is
   a t-test with df = 650 - 11 - 1 = 638; at alpha = 5e-8 its two-sided
   critical value is 5.5176, not the normal 5.4513. Near the threshold this
   alone lowers power by 0.026 at both loci (0.406 → 0.380, 0.461 → 0.435).
   The noncentral chi-square of `hbfgwas.power` has the same approximation
   (see `docs/power.md`).
2. **Fixed loci treated as random noise** (`fixed_loci_bias.txt`). With the
   genotypes held fixed, the other planted loci are a fixed part of the
   residual. They inflate the standard error the t-test *estimates* but add
   no spread to beta_hat. Treating them as noise (scaling the noncentrality
   down, as the Phase 1 diagnostic and a first version of the new formula
   did) over-predicts weak loci: by 0.005-0.007 at 40,000 redraws in the test
   scenario (4-5 SE). The correct treatment raises the effective critical
   value by sqrt(1 + lambda/df) and keeps the noncentrality unscaled.

## Result with both fixed (`hbfgwas.predict.conditional_power`)

| locus | Phase 1 diagnostic | + t critical value | full conditional model | recorded empirical [95% CI] | gap |
|---|---|---|---|---|---|
| BCL11A | 0.999 | 0.998 | 0.999 | 0.998 [0.993, 1.000] | 0.6 SE |
| HBS1L-MYB | 0.406 | 0.380 | 0.389 | 0.375 [0.345, 0.406] | 0.9 SE |
| HBB cluster | 0.461 | 0.435 | 0.436 | 0.419 [0.388, 0.450] | 1.1 SE |

The remaining gaps are within Monte Carlo error. On the separate scenario,
the corrected formula is within 2 SE of 40,000-redraw empirical power in
all 13 comparisons (with and without PCs, with and without a shift).

The proof of the redesigned check is the pre-registered run
`runs/2026-11-10-sim-v2/`, on fresh phenotype redraws.
