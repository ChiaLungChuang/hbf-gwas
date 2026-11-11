# Changelog

All notable changes to this project are recorded here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Decisions and their
reasons are in [docs/decisions.md](docs/decisions.md).

## [Unreleased]

### Added
- `docs/BRIEF.md` (the brief and addendum, verbatim) and `CLAUDE.md`.
- Item 1 scaffold: `hbfgwas` package (src layout, Python 3.11+), pytest,
  ruff, GitHub Actions CI (lint, shellcheck, forbidden-wording check on files
  and commit messages, tests on Python 3.11 and 3.12).
- Item 2 `hbfgwas.loci`: the three loci (BCL11A rs1427407, HBS1L-MYB
  rs9399137, HBB cluster rs7482144) with GRCh37 positions, alleles,
  1000 Genomes AFR frequencies and literature effect sizes, each with its
  source. Every value is `verified: false`; the 1000 Genomes AFR frequency of
  rs7482144 is null (not found). `python -m hbfgwas.loci check` and
  `pytest -m realdata` fail until the user signs the table off.
- Item 3 `hbfgwas.power`: analytic power of the 1-df additive test via the
  noncentral chi-square (ncp = n R2/(1-R2), optional covariate R2), R2 from
  beta + MAF, and inverses `required_r2` / `required_n`.
- Item 4 `hbfgwas.simulate`: Balding-Nichols genotypes for several
  subpopulations with configurable Fst, an Fst moment estimator, the analytic
  expected null chi-square under stratification, and a three-component
  phenotype (planted effects, subpopulation shift, noise) that accepts real
  genotype matrices and many repetitions at once. Calibration tests fail if
  the simulation stops producing structure.
- Item 5 `hbfgwas.assoc`: vectorised additive linear regression with
  covariates (Frisch-Waugh-Lovell, t-distribution p-values, many phenotypes
  at once), exact PCA of the standardised genotype matrix, lambda_GC, type I
  error, empirical power, QQ plot with a pointwise band, a plink2 `--glm`
  reader, and `compare_to_plink2` (allele-flip aware, tolerance on BETA and
  -log10 P).
- Item 6 `config/analysis.yaml` (pre-registration, declared before any
  experiment ran) and `hbfgwas.gate` (`python -m hbfgwas.gate`): QC, PCA,
  significance, concordance tolerances, seeds, planted effects, simulated
  experiment design with analytic predictions, and the brief's acceptance
  criteria. `tests/test_config.py` pins every value the brief specifies.
- Item 7 `hbfgwas.experiment` and `scripts/run_stratification.py`: the
  pre-registered stratification experiment on simulated genotypes (lambda_GC
  and type I error without and with 10 PCs, planted-locus recovery, empirical
  power over 1000 phenotype redraws), two-panel QQ figure and summary JSON in
  `docs/figures/`, full results under `runs/`. CI re-runs it, applies the
  gate, and checks it reproduces the committed summary.
- Item 8 real-data scripts `scripts/00_download_1kg.sh` to `05_report.py`,
  `scripts/_common.sh`, `hbfgwas.pipeline` (sentinel lookup by rsID or
  GRCh37 position, phenotype simulation on exported real genotypes, report
  with Python-vs-plink2 concordance), `resources/long_range_ld_grch37.tsv`
  (all entries unverified), and `tests/synthetic_1kg.py`, a synthetic
  plink2 fileset on which CI runs 00 -> 05 end to end with a pinned official
  plink2 build. Not run on real data.
- Item 9 README (claim and non-claim first; status table with separate
  simulated-genotype and real-genotype columns, the real column "not yet
  run") and `docs/running-real-data.md` (Phase 2 walkthrough).

### Changed after the Phase 1 review
- A (issue #10): `hbfgwas.predict.conditional_power` (power conditional on
  the realised genotypes and covariates, t critical value, leakage from other
  planted loci); experiment sim-v2 with 2000 redraws at the sentinels, a
  proposed per-locus power-calibration criterion (0.05) and per-allele
  effects planted directly; post hoc diagnostic of the Phase 1 gap in
  `runs/2026-11-10-power-diagnostic/`.
- Criteria amendment sim-v3 (issue #14), made after the failed sim-v2
  result: power calibration gated at 0.05; a locus is required in the single
  scan only at conditional predicted power >= 0.99; other loci reported. One
  pre-registered seed for sim-v3 (not yet set). CI skips the gate until it is
  set. `tests/test_records.py` pins the sim-v1 and sim-v2 records. Post hoc
  audit for #14 in `runs/2026-11-11-hbs1l-myb-audit/`. README history of
  sim-v1 (passed), sim-v2 (failed) and sim-v3.
- D-027: D-022 recorded as Claude's decision, approved by the user on
  2026-11-11 after review; CLAUDE.md now requires asking before any change to
  what is simulated or how a criterion is judged.
- sim-v3: seed 20261111 pre-registered (c2ca787) and predictions declared
  (88d4530), each pushed before the single run; gate PASSED
  (`runs/2026-11-11-sim-v3/`). `docs/figures/` now shows sim-v3.
- B: `docs/power.md`: what beta_sd is relative to, the allele frequency of each
  simulated locus and why, sim-v1's analytic formula (with and without the
  partial-R2 step), and the conditional formula the criteria use.
- C: a sentinel absent from the callset or failing QC is "not testable in this
  callset" (no proxy), in the sentinel report, the real-data report and the
  gate (status NOT TESTABLE); the gnomAD AFR frequency of rs7482144 is flagged
  suspect in the locus table and the config.
- D: HWE QC uses `--hwe 1e-6 midp keep-fewhet` (only excess heterozygosity
  filtered); the README explains why. The synthetic fileset gains
  heterozygote-deficit records that must be kept.
- E: `03_pca.sh` states explicitly when the KING step removes no samples
  (and lists any it does remove), like the missingness filters in 02.
- F: `scripts/03b_pc_plot.py` draws PC1 vs PC2 by the seven AFR populations at
  step c, right after the PCA; `05_report.py` no longer draws it.
- G: `01_subset_afr.sh` reads the compressed variant file with the documented
  `--pfile <prefix> vzs` modifier instead of naming `--pgen/--pvar/--psam`.
  CI's plink2 pin will be changed to the build the user installs locally once
  the version string is known.
- H: `hbfgwas.predict.predicted_lambda` (lambda_GC from genotypes before any
  phenotype) and `scripts/03c_prereg_predictions.sh` (step d: predicted
  lambda_GC and conditional power from the real genotypes and PCs). Checked on
  simulated phenotypes and on the sim-v3 genotypes; documented in
  `docs/power.md`, including that the sim-v1 analytic lambda formula
  over-counts by about h/K.
- G (follow-up): plink2 pinned to Beta 1.4 (7 Nov 2026) in `config/analysis.yaml`
  for the user's Mac (arm64) and CI (Linux x86_64); CI checks the version
  string and every real-data script refuses another version.
