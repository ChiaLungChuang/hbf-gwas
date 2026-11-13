# Run 2026-11-12: step e benchmark on SIMULATED genotypes of the real size

`scripts/benchmark_step_e.py`, seed 20261114. **No real genotypes were read
and no phenotype touched real data.**

**The simulated fileset.** plink2 `--dummy` wrote random genotypes for **652
samples x 14,698,808 variants**, the real counts after QC: a 1.8 GB fileset
in a scratch directory, deleted afterwards. Each expensive part of
`scripts/step_e.sh` then ran on it, in its own process. The prediction code
was the complete stack's, including the step d amendment's 2000-draw
predictive distribution.

**The machine.** macOS arm64, 8 cores, 16 GB RAM, plink2
`PLINK v2.0.0-b.1.4 M1 (7 Nov 2026)` with its default threads and memory.

| Step e part | Wall time | Peak memory | Written |
|---|---|---|---|
| integrity check: export 20,000 null variants | 2.1 s | 1.0 GB | 25 MB |
| integrity check: recompute the declared predictions | 83.6 s | 1.1 GB | none |
| 04: `--glm`, no covariates | 6.3 s | 0.9 GB | 1.24 GB |
| 04: `--glm`, 10 PCs | 10.0 s | 1.0 GB | 1.24 GB |
| 04: concordance export | 0.7 s | 1.0 GB | 2.5 MB |
| 04b: 2000 redraws at 2 sentinels | 1.8 s | 0.2 GB | none |
| 04c: the two scans without NA20314 (**not run**; assumed equal to 04's) | 16.3 s | about 1.0 GB | 2.61 GB |
| 05: read both genome-wide results, lambda_GC, type I error, QQ figure | **321.2 s** | **2.5 GB** | 0.1 MB |
| **step e, expected** | **442 s, about 7.4 min** | **2.5 GB** | **5.25 GB** |

**Notes:**
- **05 dominates.** It reads two 1.24 GB tables with pandas and draws the QQ
  figure. The sensitivity scan's two tables are read on top of that, so add
  roughly another 2.5 minutes, with peak memory about the same, since they
  are read after the primary ones.
- **This benchmark leaves out one step of the integrity check.** On real data
  it also reads the 14.7M-line `.pvar` in Python to sample the null variants
  (`pipeline sample-ids`); the benchmark used plink2's `--thin-count`
  instead. On the real genotypes, the primary 03c run of the step d amendment
  took **117 s** from start to the primary predictions (its log,
  `runs/2026-11-12-stepd-amend/logs/03c_prereg_predictions.stderr.txt`),
  compared with 86 s here.
- **Realistic expectation:** about **8 to 10 minutes**, a peak of about
  **2.5 GB** of memory, and about **5.3 GB** of new files in `data/1kg/`
  (the four genome-wide result tables). `runs/<date>-stepe/` itself stays
  small (JSON, two figures, logs).
- plink2 `--dummy` genotypes have no LD and no structure. That does not
  change the cost of `--glm` or of reading its output, which depend on the
  number of samples and variants.
