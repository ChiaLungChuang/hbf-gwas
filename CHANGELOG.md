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
