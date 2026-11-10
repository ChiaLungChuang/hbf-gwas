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
