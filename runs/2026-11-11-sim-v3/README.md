# Run 2026-11-11: sim-v3 — GATE PASSED

The single pre-registered run of sim-v3 (criteria amended after the failed
sim-v2, D-024). Genotypes and phenotype are simulated.

Order of events (all on GitHub before the run):
1. `#17` (criteria) reached `main` via `#18`; `main`'s tree verified identical
   to the `#17` tip.
2. Seed commit `c2ca787` (`seeds.sim_v3: 20261111`) pushed and confirmed on
   GitHub.
3. Predictions commit `88d4530` (computed from the genotypes and PCs only)
   pushed and confirmed on GitHub.
4. This run, at commit `88d4530` (see `results.json` → `provenance`):

       python scripts/run_stratification.py --out-dir runs/2026-11-11-sim-v3
       python -m hbfgwas.gate --results runs/2026-11-11-sim-v3/results.json

## Gate (`gate_log.txt`)

| criterion | observed | rule | status |
|---|---|---|---|
| unadjusted lambda_GC | 1.3422 (predicted 1.369) | > 1.2 | PASS |
| PC-adjusted lambda_GC | 0.9980 | [0.97, 1.05] | PASS |
| PC-adjusted type I error at 0.05 | 0.0505 | [0.04, 0.06] | PASS |
| BCL11A single scan (predicted power 0.9958, REQUIRED) | p = 3.41e-09, hit | p < 5e-8 | PASS |
| HBS1L-MYB single scan (predicted power 0.0311) | p = 7.51e-06, miss | reported, not gated | REPORTED |
| HBB cluster single scan (predicted power 0.2160) | p = 2.94e-06, miss | reported, not gated | REPORTED |
| BCL11A power calibration | 0.995 [0.990, 0.997] vs 0.9958 | \|diff\| <= 0.05 | PASS |
| HBS1L-MYB power calibration | 0.030 [0.023, 0.038] vs 0.0311 | \|diff\| <= 0.05 | PASS |
| HBB cluster power calibration | 0.203 [0.186, 0.222] vs 0.2160 | \|diff\| <= 0.05 | PASS |
| declared predictions reproduced | yes | | PASS |

**GATE PASSED.** Estimated Fst 0.0100 (configured 0.01). In this genotype
draw HBS1L-MYB has sample frequency 0.0377 and negative leakage (-0.092)
from the other planted loci, so its predicted power is low.

Observation, not a criterion: across sim-v2 and sim-v3 all six empirical
powers sit at or below their predictions (-0.001 to -0.013), each inside its
95% interval and inside the 0.05 tolerance.
