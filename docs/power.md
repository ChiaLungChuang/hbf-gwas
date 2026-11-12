# Power: what is planted, and how power is predicted

This page answers review question B (2026-11-11). Every number on it can be
reproduced with `hbfgwas.power`, `hbfgwas.predict` and the committed run
records. **The phenotype is always simulated.**

## 1. The scale: what "beta in SD units" is relative to

The simulated phenotype has three parts, and its total variance is 1 (in
expectation; `hbfgwas.simulate.simulate_phenotype`):

    y = genetic + shift + noise
    Var(shift) = 0.2  (exactly, each draw)
    Var(noise) = 1 - Var(genetic) - 0.2

`beta_sd` is the per-allele effect of the **HbF-increasing allele**, in units
of the standard deviation of this **total** simulated phenotype (so the
shift is included in the SD). Since sim-v2 it is planted exactly as given:
each copy of the effect allele adds `beta_sd` to y. This applies to both
simulated genotypes and `scripts/04_assoc.sh` on real genotypes.

Where the numbers come from (`src/hbfgwas/loci.py`, every value still
`verified: false`):

| locus | beta_sd | origin |
|---|---|---|
| BCL11A | 0.456435 | sqrt(PVE / (2 p (1 - p))) with PVE 0.076 and MAF 0.24 (Wonkam et al. 2025 global meta-analysis, cube-root HbF) |
| HBS1L-MYB | 0.571 | Lettre et al. 2008 Table 1: "standard deviation units for the minor allele" in CSSCD |
| HBB cluster | 0.407 | Lettre et al. 2008 Table 1, same units |

**Assumption.** Lettre et al.'s SD is the SD of their transformed,
covariate-adjusted HbF in sickle cell patients. Treating it as an SD of our
total simulated phenotype is a modelling choice, not an equivalence.

sim-v1 was different. It fixed the variance explained,
r2 = 2 p (1 - p) beta_sd^2 at the **ancestral** frequency, and rescaled beta to
whatever gave that r2 in the sample. HBS1L-MYB therefore received 0.469 per
allele, not 0.571, because its sample frequency (0.0623) was above the
ancestral 0.0416. sim-v2 switched to planting `beta_sd` directly (D-022: Claude's
decision, approved by the user after review, D-027). The `r2` values still in
the config are informational: the variance explained at the ancestral
frequency.

## 2. Allele frequency of each simulated locus, and why

Each planted variant is drawn by the Balding-Nichols model. It starts from an
**ancestral frequency** p; each of the 7 subpopulations then gets p_k ~ Beta
(mean p, variance Fst p (1 - p), Fst = 0.01), and individuals draw
Binomial(2, p_k). The sample frequency therefore differs from p in every
genotype draw.

| locus | ancestral frequency (effect allele) | why | sample frequency, sim-v1/v2 draw | sample frequency, sim-v3 draw |
|---|---|---|---|---|
| BCL11A | 0.237519 (T) | 1000 Genomes phase 3 AFR, via Ensembl GRCh37 | 0.2462 | 0.2508 |
| HBS1L-MYB | 0.041604 (C) | 1000 Genomes phase 3 AFR, via Ensembl GRCh37 | 0.0623 | 0.0377 |
| HBB cluster | 0.1238 (A) | gnomAD exomes AFR via Ensembl GRCh37, because Ensembl has no 1000 Genomes phase 3 value. **Suspect:** the user noted far too few homozygotes for that frequency, consistent with HBG1/HBG2 paralogy. | 0.1254 | 0.1177 |

The predictions the criteria use (section 4) are computed from the
**realised** genotypes, so they use the sample frequency, not p. On real
genotypes the frequency is the callset's own. HBB is "not testable in this
callset" if rs7482144 is absent or fails QC.

## 3. sim-v1's analytic formula (`hbfgwas.power`), and why it gave 0.997

For a 1-df additive test with variance explained r2, the statistic is
approximately noncentral chi-square with

    ncp = n r2_p / (1 - r2_p),   r2_p = r2 / (1 - covariate_r2)

and power = P(chi2_1(ncp) > 29.7168) at alpha = 5e-8, i.e. the normal
critical value 5.4513, squared.

**Yes: in sim-v1 the PC-removed structure variance was excluded from the
denominator.** The prediction assumed the 10 PCs absorb the whole shift, so
`covariate_r2 = 0.2` and the variant's share of the *remaining* variance is
r2 / 0.8 (D-010). Your hand calculation with plain n r2 / (1 - r2) is the
same formula without that step:

| locus | r2 = 2p(1-p)beta^2 | plain ncp | plain power | partial r2 = r2 / 0.8 | ncp | sim-v1 declared power |
|---|---|---|---|---|---|---|
| BCL11A | 0.07546 | 53.05 | **0.9666** | 0.09432 | 67.70 | **0.9973** |
| HBS1L-MYB | 0.02600 | 17.35 | 0.0993 | 0.03250 | 21.84 | 0.2181 |
| HBB cluster | 0.03594 | 24.23 | 0.2984 | 0.04492 | 30.57 | 0.5310 |

Worked example, BCL11A: r2 = 2 x 0.237519 x 0.762481 x 0.456435^2 = 0.07546.
Plain: ncp = 650 x 0.07546 / 0.92454 = 53.05, so power ~ Phi(sqrt(53.05) -
5.4513) = Phi(1.833) = 0.9666. With the shift removed: 650 x 0.09432 /
0.90568 = 67.70, so power ~ Phi(8.228 - 5.4513) = 0.9973.

Two further approximations in this formula, both found after sim-v1:

* It uses the **normal** critical value. The actual test is a t-test with
  df = 650 - 11 - 1 = 638, whose critical value is 5.5176. Near the threshold
  this overstates power by about 0.026.
* It is **unconditional**: it averages over genotype draws, while every run
  uses one fixed draw, with its own sample frequencies and chance
  correlations between planted variants.

`hbfgwas.power` is kept because the brief asks for the noncentral
chi-square. The criteria no longer use it.

## 4. The prediction the criteria use (`hbfgwas.predict.conditional_power`)

This is the power of the actual t-test at each planted variant, **conditional
on the realised genotype matrix and covariates**, computed before any
phenotype is drawn. With Q an orthonormal basis of [1, PCs] and residuals
taken after projecting on Q:

1. g_r = residual genotype of the planted variant (what the PCs leave);
   SS = g_r . g_r.
2. Expected estimate: E[beta_hat] = g_r . (G beta)_r / SS. This is its own beta plus
   **leakage** from the other planted variants it is correlated with.
3. beta_hat ~ N(E[beta_hat], sigma_e^2 / SS) exactly, where sigma_e^2 is the
   noise variance, 1 - Var(G beta) - 0.2.
4. The residual sum of squares is sigma_e^2 times a noncentral chi-square
   with df = n - rank([1, PCs]) - 1 and noncentrality lambda = RSS_fixed /
   sigma_e^2. RSS_fixed is the part of the other loci and of the shift that
   the model leaves behind. It is fixed for fixed genotypes, so it raises the
   effective critical value to t_crit sqrt(1 + lambda / df) rather than adding
   spread.
5. Power = P(|T| > effective critical value), with T noncentral t(df,
   E[beta_hat] / sqrt(sigma_e^2 / SS)), averaged over 2000 draws of the
   subpopulation shift from the simulation's own law.

Checked against Monte Carlo (40,000 redraws; within 2 SE in all 13 scenario
comparisons, `runs/2026-11-10-power-diagnostic/`) and against the 2000-redraw
empirical power in every run since:

| run | locus | sample freq | planted beta | predicted (conditional, PC-adjusted) | empirical [95% CI] |
|---|---|---|---|---|---|
| sim-v2 | BCL11A | 0.2462 | 0.4564 | 0.9996 | 0.999 [0.996, 1.000] |
| sim-v2 | HBS1L-MYB | 0.0623 | 0.571 | 0.8090 | 0.806 [0.788, 0.823] |
| sim-v2 | HBB cluster | 0.1254 | 0.407 | 0.4621 | 0.451 [0.429, 0.473] |
| sim-v3 | BCL11A | 0.2508 | 0.4564 | 0.9958 | 0.995 [0.990, 0.997] |
| sim-v3 | HBS1L-MYB | 0.0377 | 0.571 | 0.0311 | 0.030 [0.023, 0.038] |
| sim-v3 | HBB cluster | 0.1177 | 0.407 | 0.2160 | 0.203 [0.186, 0.222] |

All six empirical values are at or slightly below prediction (-0.001 to
-0.013); each is inside its interval and inside the 0.05 tolerance.

## 5. How the criteria use it (sim-v3 rule, D-024)

* Gated: |empirical - predicted| <= 0.05 per locus over 2000 redraws.
* Gated: a locus with predicted power >= 0.99 must reach p < 5e-8 in the
  single genome-wide scan.
* Reported, not gated: every other locus's single-scan hit or miss, next to
  its predicted power.

The Phase 2 step d real-data pre-registration uses the same rule, with
predictions computed from the real genotypes and PCs.

## 6. Predicted lambda_GC from genotypes (review item H)

`hbfgwas.predict.predicted_lambda` predicts the unadjusted (and PC-adjusted)
lambda_GC from the genotypes alone, before any phenotype is drawn. For each
draw of the random parts that act on all null variants at once (the
subpopulation shift, and the subpopulation means of the noise), each null
variant's t statistic is N(delta_j, tau^2). Here delta_j = sqrt(n) cov(g_r,
s_r) / (sd(g_r) sd(y_r)), computed from the **realised** genotypes, and
tau^2 is the non-shift share of the residual phenotype variance. lambda is
the median of the implied chi-square mixture over the null variants,
divided by 0.4549. The summary is the median and the 5th-95th percentile
over draws. Estimating lambda as a median over M independent variants adds
an SD of about 2.333 lambda / sqrt(M).

Two corrections were needed while building it, both found by tests and
checks:

* The noise's own subpopulation means act like a second, small shift on all
  variants at once, so they are drawn as well.
* The within-subpopulation sampling part of the genotype means is already in
  delta_j, so the noise term has variance 1 - h, not 1. With this,
  E[chi2] ~ (1 - h) + 2nFh/K + h(1 - 1/K). For the simulated design that is
  1.342, which shows that `simulate.expected_null_chi2` (1.369, declared in
  sim-v1) over-counts by about h/K. That declared value stays as it was; it
  was informational and never a criterion.

Check on the sim-v3 genotypes (20,000 random null variants, 200 draws; no
phenotype): predicted unadjusted lambda_GC 1.372 [1.292, 1.471], with PCs
1.000 [1.000, 1.000]. sim-v3 observed 1.342 and 0.998.

On real data, `scripts/03c_prereg_predictions.sh` computes the same
quantities from the real genotypes and PCs, together with the conditional
power of each testable sentinel, for the Phase 2 step d pre-registration.


## 7. The step d amendment: the predictive distribution of the observed statistics

Added on 2026-11-12, before any real-data phenotype existed, when the user
replaced "unadjusted lambda_GC > 1.2" with an interval and asked for each
gated criterion's pass probability (D-048, D-049).
`hbfgwas.predict.predicted_null_distribution` predicts the **observed**
lambda_GC and type I error, not only their central value.

* **The fixed and the random parts.** The step e phenotype is
  y = f + e. f, planted effects plus the declared shift, is fixed once the
  genotypes are. e ~ N(0, sigma_e^2) is the only random part.
* **The structured noise is drawn exactly.** e is split into its projection
  on the population indicators and the 10 PCs, drawn in each of 2000 draws,
  and the rest. In each draw every null variant's t statistic is
  N(delta_j, tau^2), with delta_j computed from the realised genotypes and
  f + the drawn projection, as in section 6. lambda is the mixture median /
  0.4549; the type I error is the mean of P(|t| > 1.96).
* **Estimation noise is added.** The observed median over M_eff
  effectively independent variants has SD 2.333 lambda / sqrt(M_eff), with
  M_eff = 215,522 (the LD-pruned variants of step c). The observed type I
  fraction has SD sqrt(t1 (1 - t1) / M_eff). Using a 20,000-variant sample
  adds a bootstrap SD.
* **Pass probabilities:**
  - an interval criterion on lambda or type I: the average over draws of
    P(draw + estimation noise inside it);
  - the power calibration: Binomial(2000, predicted power), and the worst
    case with the true power off by up to 0.013.

The numbers and a check on simulated genotypes are in
`runs/2026-11-12-stepd-amend/README.md`. In that check the PC-adjusted
predictions were well calibrated. The unadjusted ones covered the observed
value slightly more often than nominal: 0.993 inside the 1st-99th
percentile.
