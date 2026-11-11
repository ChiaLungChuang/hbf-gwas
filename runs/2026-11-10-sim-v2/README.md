# Run 2026-11-10: sim-v2 (stratification + redesigned power check) — GATE FAILED

The single pre-registered run of sim-v2 (`config/analysis.yaml`, prereg
`sim-v2`, committed in `df3b117` before this run). Genotypes and phenotype
are simulated. Code at commit `df3b117` (see `results.json` → `provenance`).

    python scripts/run_stratification.py --out-dir runs/2026-11-10-sim-v2
    python -m hbfgwas.gate --results runs/2026-11-10-sim-v2/results.json

## Gate (`gate_log.txt`)

| criterion | observed | rule | status |
|---|---|---|---|
| unadjusted lambda_GC | 1.3825 (predicted 1.369) | > 1.2 | PASS |
| PC-adjusted lambda_GC | 1.0030 | [0.97, 1.05] | PASS |
| PC-adjusted type I error at 0.05 | 0.0509 | [0.04, 0.06] | PASS |
| BCL11A recovery (predicted power 0.9996, REQUIRED) | p = 6.81e-22 | p < 5e-8 | PASS |
| **HBS1L-MYB recovery (predicted power 0.8090, REQUIRED)** | **p = 1.24e-06** | p < 5e-8 | **FAIL** |
| HBB cluster recovery (predicted power 0.4621) | p = 6.81e-07 | not required | underpowered, not required |
| BCL11A power calibration | 0.999 vs 0.9996 | \|diff\| <= 0.05 (proposed) | PASS |
| HBS1L-MYB power calibration | 0.806 [0.788, 0.823] vs 0.8090 | \|diff\| <= 0.05 (proposed) | PASS |
| HBB cluster power calibration | 0.451 [0.429, 0.473] vs 0.4621 | \|diff\| <= 0.05 (proposed) | PASS |
| declared predictions reproduced | yes | | PASS |

**GATE FAILED** on one criterion: HBS1L-MYB was REQUIRED (conditional
predicted power 0.809 >= 0.8) and the single genome-wide scan did not reach
genome-wide significance at it. Per the brief, nothing was changed after this
result; it is reported to the user.

## What the numbers say

* The redesigned power check (issue #10) agrees: empirical power over 2000
  redraws is within 0.011 of the conditional prediction at every locus and
  every prediction lies inside its 95% interval.
* The failure is the single-run recovery of a locus whose predicted power is
  0.81: a REQUIRED locus at that power is expected to be missed about 19% of
  the time. HBS1L-MYB became REQUIRED only in sim-v2, because sim-v2 plants
  the per-allele effect directly (D-022); under sim-v1's fixed-r2 scaling it
  would have been ~0.39 and not required.
