# CLAUDE.md

## Read the brief first, every time

At the start of **every session**, and again at the start of **every phase**
(Phase 1, Phase 2, and each Phase 2 step a-f), re-read `docs/BRIEF.md` in full
before doing anything else. It holds the original brief followed by an
ADDENDUM; the addendum overrides the brief wherever they conflict.

After the brief, read `docs/decisions.md` and `CHANGELOG.md` to see what has
been decided and where the work stands, and check the open GitHub issues.

## Hard rules (summary only; docs/BRIEF.md is authoritative)

- The phenotype is always simulated. 1000 Genomes has no HbF phenotype.
  Never use the forbidden wording listed in docs/BRIEF.md in code, docs, or
  commit messages.
- Never commit genotype data (`data/` is git-ignored). Do commit everything
  that documents a run, under `runs/<date>/`.
- Never report a test as passing unless it was run and seen to pass.
- Never tune a simulation to hit a number. Fix seeds; report what comes out.
- Never edit an acceptance criterion in `config/analysis.yaml` after seeing a
  result. If a criterion fails, stop and report it to the user.
- Ask the user BEFORE any change to what is simulated (effect sizes,
  frequencies, variance components, the phenotype model, seeds' roles) or to
  how a criterion is judged (which prediction decides REQUIRED, thresholds,
  gated vs reported). Declaring it in the config before a run is not enough
  (docs/decisions.md D-027).
- Locus table values stay `verified: false` until the user signs them off by
  hand.
- All work goes through a branch and a PR whose description says what was
  done and what the output was. Push after every working step. One GitHub
  issue per phase and per open decision.
- Record each decision in `docs/decisions.md` (who made it, and why) and each
  change in `CHANGELOG.md`.
- Phase gates: stop after Phase 1 and report; do not start Phase 2 until the
  user replies. Inside Phase 2, stop and wait at every point the addendum
  names.

## Commands

```bash
/opt/homebrew/bin/python3.12 -m venv .venv && .venv/bin/pip install -e ".[dev]"
.venv/bin/python -m pytest              # unit tests (plink2 tests skip without plink2)
.venv/bin/python -m pytest -m realdata  # real-data readiness; fails until loci are signed off
.venv/bin/ruff check . && .venv/bin/ruff format --check .
.venv/bin/python scripts/check_language.py --commits
```
