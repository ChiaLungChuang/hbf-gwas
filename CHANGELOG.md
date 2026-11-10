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
- Item 3 `hbfgwas.power`: analytic power of the 1-df additive test via the
  noncentral chi-square (ncp = n R2/(1-R2), optional covariate R2), R2 from
  beta + MAF, and inverses `required_r2` / `required_n`.
