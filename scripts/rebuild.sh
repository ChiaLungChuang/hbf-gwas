#!/usr/bin/env bash
# Rebuild everything that can be rebuilt without real data, from a clean clone.
#
#   scripts/rebuild.sh            # rebuilds the committed HEAD of this checkout
#   REQUIRE_PLINK2=1 scripts/rebuild.sh   # fail instead of skipping plink2 tests
#
# It clones the committed HEAD into a new temporary directory (uncommitted
# changes, data/ and any existing .venv are not used). Then it creates a new
# virtualenv without the pip cache, installs the package, and runs:
# - ruff and shellcheck (if installed);
# - the forbidden-wording check on files and commit messages;
# - every test;
# - the pre-registered simulated experiment with its acceptance gate;
# - a check that the result reproduces docs/figures/sim_stratification_summary.json.
# Nothing is written to this checkout. No phenotype is drawn on real genotypes.
# PYTHON selects the interpreter (default python3; it must be 3.11 or newer).
set -euo pipefail

src=$(git -C "$(dirname "$0")/.." rev-parse --show-toplevel)
head=$(git -C "$src" rev-parse HEAD)
work=$(mktemp -d "${TMPDIR:-/tmp}/hbfgwas-rebuild.XXXXXX")
trap 'rm -rf "$work"' EXIT

echo "rebuild: commit $head in a clean clone"
git clone --quiet --no-local "$src" "$work/repo"
cd "$work/repo"
git checkout --quiet "$head"

"${PYTHON:-python3}" -m venv .venv
.venv/bin/python -m pip install --quiet --no-cache-dir --upgrade pip
.venv/bin/python -m pip install --quiet --no-cache-dir -e ".[dev]"

echo "rebuild: lint"
.venv/bin/ruff check .
.venv/bin/ruff format --check .
if command -v shellcheck >/dev/null; then
  # shellcheck disable=SC2046
  shellcheck $(git ls-files '*.sh')
else
  echo "rebuild: shellcheck not installed, skipped (CI runs it)"
fi
.venv/bin/python scripts/check_language.py --commits

echo "rebuild: tests"
.venv/bin/python -m pytest -q

echo "rebuild: pre-registered simulated experiment and gate"
key=$(.venv/bin/python -m hbfgwas.pipeline config-get simulated_experiment.seed_key)
if [ "$(.venv/bin/python -m hbfgwas.pipeline config-get "seeds.$key")" = "null" ]; then
  echo "rebuild: seeds.$key is not pre-registered yet, so there is nothing to gate"
  echo "rebuild: OK ($head)"
  exit 0
fi
.venv/bin/python scripts/run_stratification.py --out-dir "$work/sim" --figures-dir ""
.venv/bin/python -m hbfgwas.gate --results "$work/sim/results.json"
SIM_DIR="$work/sim" .venv/bin/python - <<'PY'
import json
import os
import sys

from hbfgwas import experiment

ref = json.load(open("docs/figures/sim_stratification_summary.json"))
new = json.load(open(os.path.join(os.environ["SIM_DIR"], "sim_stratification_summary.json")))
diffs = experiment.compare_summaries(ref, new)
print("\n".join(diffs) or "summary reproduces docs/figures/sim_stratification_summary.json")
sys.exit(1 if diffs else 0)
PY

echo "rebuild: OK ($head)"
