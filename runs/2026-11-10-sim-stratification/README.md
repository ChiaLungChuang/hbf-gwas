# Run 2026-11-10: stratification experiment on simulated genotypes

The single pre-registered run of experiment 1 (`config/analysis.yaml`,
`simulated_experiment`). Genotypes **and** phenotype are simulated.

- Code: commit `fbc75d9` (see `results.json` → `provenance`), run with
  `python scripts/run_stratification.py --out-dir runs/2026-11-10-sim-stratification`
  on macOS arm64, Python 3.12, then
  `python -m hbfgwas.gate --results runs/2026-11-10-sim-stratification/results.json`.
- The configuration was committed (`f63e8bd`) before the runner existed.

## Results (`run_log.txt`, `gate_log.txt`)

| | predicted in advance | observed | criterion | status |
|---|---|---|---|---|
| unadjusted lambda_GC | 1.369 | 1.3806 | > 1.2 | PASS |
| PC-adjusted lambda_GC | | 0.9995 | [0.97, 1.05] | PASS |
| PC-adjusted type I error at 0.05 | | 0.0507 | [0.04, 0.06] | PASS |
| unadjusted type I error at 0.05 | | 0.0966 | (none) | |
| BCL11A recovery (power 0.997, REQUIRED) | | p = 5.69e-21 | p < 5e-8 | PASS |
| HBB cluster (power 0.531) | | p = 1.02e-06 | not required | underpowered, not required |
| HBS1L-MYB (power 0.218) | | p = 1.55e-04 | not required | underpowered, not required |

Estimated Fst 0.00997 (configured 0.01). **Gate: PASSED.**

## Secondary outcome (reported, not gating): empirical vs predicted power

| locus | predicted (analytic) | empirical over 1000 redraws [95% CI] | agrees? |
|---|---|---|---|
| BCL11A | 0.997 | 0.998 [0.993, 1.000] | yes |
| HBS1L-MYB | 0.218 | 0.375 [0.345, 0.406] | **no** (higher) |
| HBB cluster | 0.531 | 0.419 [0.388, 0.450] | **no** (lower) |

**Why (exploratory, investigated after the run, nothing changed):** the 1000
redraws re-simulate the shift and noise but keep one genotype realization
fixed, so they measure power *conditional on those genotypes*; the analytic
formula averages over genotype realizations. In this realization the
BCL11A and HBS1L-MYB genotypes have a chance sample correlation of 0.081
(about 2 SE from 0), which adds about 13% to HBS1L-MYB's expected effect
estimate, and the PCs absorb 6.2% of the HBB genotype variance (its frequency
varies across subpopulations). Conditional power computed from the fixed
genotypes is 0.406 and 0.461, close to the empirical 0.375 and 0.419
(`diagnostic_conditional_power.py`, output in `diagnostic_conditional_power.txt`).
This is a weakness in how Claude designed the secondary outcome, not a
failure of a pre-registered criterion. Any re-design (e.g. redrawing
genotypes in every repetition) needs a new pre-registration approved by the
user; see the open decision issue.
