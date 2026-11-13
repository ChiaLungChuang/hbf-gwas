# Run 2026-11-12: step d amendment, predictions only (no phenotype)

Made at the user's request (2026-11-12), **before any real-data phenotype
existed**, to amend the step d pre-registration (D-048 to D-050):

* flag 1: "unadjusted lambda_GC > 1.2" is replaced by "inside the pre-stated
  1st-99th percentile prediction interval";
* for every gated criterion, the probability that a correct pipeline passes.

The original step d records in `runs/2026-11-12-stepd/` are not edited.

**No phenotype was formed and no association was run on the real
genotypes.** The run reads the 652 real genotypes (sentinels and 20,000
random null variants more than 1 Mb from them), their PCs and their
population labels. It draws only the noise's projection on 17 directions
(the population indicators and the 10 PCs), never a phenotype. Same code
path as step d: `scripts/03c_prereg_predictions.sh`, plink2
`PLINK v2.0.0-b.1.4 M1 (7 Nov 2026)`. Logs are redacted when written
(`REDACTION.md`).

## What changed in the method (D-049), and why

The original step d prediction (`predicted_lambda`) gave a median and a
5th-95th percentile over 200 draws, and its own note said it left out the
estimation noise of the observed median. A 1st-99th percentile interval for
the **observed** value needs that, and more draws in the tails.
`predict.predicted_null_distribution` makes four changes:

1. **The fixed part is the whole fixed part.** It is planted effects plus the
   declared shift (variance 0.3044), not the shift alone. The noise variance
   is 0.6999, the one `simulate_phenotype` uses.
2. **The structured part of the noise is drawn exactly.** That is the
   noise's projection on the population indicators **and the 10 PCs** (17
   dimensions), rather than on the population means only. 2000 draws, from
   seed `real_prediction`'s child 2. Children 0 and 1 are unchanged, so the
   conditional power reproduces exactly: 0.9896 and 0.1694.
3. **Estimation noise is added.** It comes from taking a median, or for the
   type I error a fraction, over a finite number of variants. The effective
   number is the 215,522 LD-pruned variants of step c. A variant bootstrap
   adds the error of predicting from a 20,000-variant sample.
4. Type I error is predicted the same way, as the mean over variants of
   P(|t| > 1.96).

## Results (`prereg_predictions.json`)

| Observed statistic | median | 1st | 5th | 95th | 99th percentile |
|---|---|---|---|---|---|
| unadjusted lambda_GC | 1.2001 | **1.1248** | 1.1459 | 1.2630 | **1.2874** |
| PC-adjusted lambda_GC | 1.0094 | 0.9966 | 1.0001 | 1.0186 | 1.0226 |
| unadjusted type I error (alpha 0.05) | 0.0736 | 0.0644 | | | 0.0848 |
| PC-adjusted type I error | 0.0511 | 0.0498 | | | 0.0523 |

Estimation SD at the median: unadjusted lambda 0.0081 (variant bootstrap
0.0054 of it), PC-adjusted 0.0054 (bootstrap 0.0020). The unadjusted median
is essentially the original 1.1999; the interval is wider because it now
covers the observed value.

**Probability that a correct pipeline passes each gated criterion:**

| Gated criterion | P(pass) | How |
|---|---|---|
| unadjusted lambda_GC in [1.1248, 1.2875] | **0.98, below 0.99** | by construction of a 1st-99th interval: the fraction of the 2000 predictive draws inside it |
| PC-adjusted lambda_GC in [0.97, 1.05] | 1 − 2.6e-12 | mean over draws of P(draw + estimation noise in the range) |
| PC-adjusted type I error in [0.04, 0.06] | 1.0 (double precision) | the same |
| power calibration BCL11A, \|emp − 0.9896\| ≤ 0.05 over 2000 redraws | 1 − <1e-15 | Binomial(2000, 0.9896); still 1.0 if the true power is off by ±0.013 |
| power calibration HBS1L-MYB, \|emp − 0.1694\| ≤ 0.05 | 1 − 4.8e-9 | Binomial(2000, 0.1694); 0.99998 if the true power is off by ±0.013 |
| **all gated criteria together** | **0.98** | the three null criteria jointly over the draws, times the two calibrations |

The calibration is the only place where the model's own error enters
explicitly. 0.013 is the largest gap between empirical and predicted power
seen in sim-v2 and sim-v3.

**Sensitivity without NA20314** (`prereg_predictions_sensitivity.json`,
never gated): conditional power 0.9897 and 0.1645. Unadjusted lambda_GC
1.1989 [1.1240, 1.2919], PC-adjusted 1.0091 [0.9962, 1.0225]. Same pass
probabilities.

## Check of the prediction on simulated genotypes (`null_prediction_check_simulated.json`)

`scripts/check_null_prediction.py`, seed 20261113, **simulated genotypes
only**. The setting:
- Balding-Nichols, Fst 0.01, the real population counts, 60,000 variants;
- planted variants at the real sentinel frequencies, and the declared shift;
- the prediction made from a 20,000-variant sample, exactly as above;
- then 1000 phenotypes scanned genome-wide.

| | observed inside 1st-99th (target 0.98) | inside 5th-95th (target 0.90) | below predicted median (target 0.50) |
|---|---|---|---|
| unadjusted lambda_GC | 0.993 | 0.947 | 0.601 |
| PC-adjusted lambda_GC | 0.983 | 0.895 | 0.488 |
| unadjusted type I error | 0.992 | 0.943 | 0.629 |
| PC-adjusted type I error | 0.989 | 0.899 | 0.565 |

PC-adjusted lambda_GC in [0.97, 1.05]: observed 0.999, predicted 0.9995.

**Reading:** the PC-adjusted predictions are well calibrated. The unadjusted
ones are slightly too wide and centred about 0.011 too high (predicted
median 1.3977, observed 1.3867). In this setting the 1st-99th coverage was
0.993, so the 0.98 above is if anything conservative.

**What this check cannot test:**
- real linkage disequilibrium (the simulated variants are independent), so
  the effective number of variants is an assumption;
- real admixture beyond what the 10 PCs capture.

Those are the parts of the real prediction that rest on judgment.

## Order of events

1. The user's request: flag 1 replaced, flag 2 unchanged, pass probabilities.
2. Code and tests: `predict.predicted_null_distribution`, the pass
   probabilities in `pipeline.prereg_predictions`, and the settings in
   `config/analysis.yaml` `prediction` (`null_draws`, `n_effective_variants`,
   `power_model_error`).
3. The simulated check (above).
4. This run on the real genotypes.
5. The criterion written from its 1st-99th interval, rounded outward to 4
   decimals: [1.1248, 1.2875].

No phenotype exists for the real genotypes at any point. The step e seed is
still empty.
