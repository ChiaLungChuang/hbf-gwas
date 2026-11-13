# Run 2026-11-12: Phase 2 step d predictions (no phenotype)

Everything here was computed from the 652 real genotypes, their PCs and their
population labels. **No phenotype was formed and no association was run**:
the conditional power is analytic in the noise, the fixed shift has no
randomness, and the lambda_GC prediction averages only over draws of the
noise's population means (Monte Carlo integration). Code: branch
`stepd-prereg`; plink2 `PLINK v2.0.0-b.1.4 M1 (7 Nov 2026)`.

Order of events (all pushed to GitHub before the next step):
1. `67a8df2`: code for a fixed population shift, an empty step e seed slot,
   and the sensitivity PCA.
2. `f0c9a5f`: seed `real_shift_draw: 20261112` declared, pushed, and
   confirmed on GitHub before use.
3. `ed3d3ce`: the seven shift values computed from that seed and the
   population counts only, then committed.
4. This run: `03d_sensitivity_pca.sh`, then `03c_prereg_predictions.sh`.

## The fixed non-genetic shift (SD units of the simulated phenotype)

| population | n | raw N(0,1) draw | shift |
|---|---|---|---|
| ACB | 95 | 1.370097 | +0.646472 |
| ASW | 55 | -1.322758 | -0.391917 |
| ESN | 99 | -0.133478 | +0.066680 |
| GWD | 113 | -1.321810 | -0.391551 |
| LWK | 97 | -1.757012 | -0.559369 |
| MSL | 85 | 1.336254 | +0.633422 |
| YRI | 108 | -0.349562 | -0.016644 |

Centred and scaled over the 652 samples: mean 0, variance 0.200000.

## Predictions (`prereg_predictions.json`)

| | BCL11A rs1427407 | HBS1L-MYB rs9399137 |
|---|---|---|
| effect-allele frequency, 652 samples | T 0.2393 | C 0.0391 |
| planted beta (SD units) | 0.456435 | 0.571 |
| expected beta-hat = beta + leakage (PC-adjusted) | 0.4606 (+0.0041) | 0.5863 (+0.0153) |
| genotype variance kept after 10 PCs | 0.949 | 0.968 |
| noncentrality / effective critical t (df 640) | 8.160 / 5.821 | 5.015 / 5.988 |
| **conditional power, PC-adjusted (the criterion's prediction)** | **0.9896** | **0.1694** |
| conditional power, unadjusted | 0.9931 | 0.0262 |
| analytic power, `hbfgwas.power` (reference) | 0.9975 | 0.1822 |
| REQUIRED in the single scan (>= 0.99)? | no | no |

HBB cluster: not testable in this callset (rs7482144 absent).

Predicted lambda_GC over 20,000 random null variants more than 1 Mb from
the sentinels (200 draws of the noise's population means), as median
[5th, 95th percentile]: **unadjusted 1.1999 [1.1449, 1.2703]**; **PC-adjusted
1.0045 [1.0028, 1.0066]**. Noise variance 0.6999 = 1 - genetic 0.1001 -
shift 0.2.

## Sensitivity analysis without NA20314 (`prereg_predictions_sensitivity.json`)

PCs recomputed on 651 samples, same 215,522 pruned variants
(`sens_pca_eigenval.txt`: 4.03, 3.00, 2.64, ... vs 4.19, 3.21, 2.67 with
NA20314). Shift variance 0.200071 with the same per-population values.
Conditional power: BCL11A 0.9897, HBS1L-MYB 0.1645. Predicted lambda_GC:
unadjusted 1.1984 [1.1444, 1.2690], PC-adjusted 1.0045 [1.0028, 1.0063].
