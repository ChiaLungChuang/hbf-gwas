#!/usr/bin/env bash
# Phase 2 step e in one command: the pre-registered real-data experiment, run once.
# The phenotype is SIMULATED by hbfgwas; 1000 Genomes has no HbF phenotype.
#
#   scripts/step_e.sh
#
# Refuses to start, before any phenotype is simulated, unless:
#   * the tag prereg-v1 exists here and on GitHub (same commit), and
#     config/analysis.yaml is unchanged since it;
#   * acceptance.real (the real-data criteria) is filled in;
#   * the locus table passes the analysis check (signed off by the user);
#   * the step e seed (config/step_e_seed.yaml) is filled;
#   * the population shift is declared;
#   * the PCA (03) and, if samples are excluded for the sensitivity analysis,
#     the sensitivity PCA (03d) have been run;
#   * RUN_DIR has no results.json yet (step e runs once);
#   * integrity check: each locus's conditional power and the predicted lambda_GC,
#     recomputed from the genotypes (03c), equal the values declared in
#     acceptance.real (within 1e-4). Not an acceptance criterion.
# Then it runs, in order:
#   04_assoc.sh            simulated phenotype; plink2 --glm without and with the PCs
#   04b_power_redraws.sh   redraws at the sentinels (acceptance.real.power_calibration)
#   04c_sensitivity_scan.sh  the same phenotype without the excluded samples (never gated)
#   05_report.py           results.json, figures, Python-vs-plink2 concordance
#   python -m hbfgwas.gate --experiment real   PASS / FAIL per pre-registered criterion
# Exit status: 0 every gated criterion passed; 3 at least one gated criterion
# failed (change nothing and report it); 4 the gate passed but the
# Python-vs-plink2 concordance check failed; 1 refused or an error.
# Overrides, for the synthetic test fileset only (as in 04_assoc.sh):
#   HBFGWAS_LOCI_TSV=<file>, HBFGWAS_SKIP_PREREG_CHECK=1
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=scripts/_common.sh
. "$here/_common.sh"
need "$PLINK2"
setup_dirs
capture_logs "$@"
check_plink_version

# --- preconditions: nothing below runs unless all of these hold ------------------
[ ! -e "$RUN_DIR/results.json" ] ||
	die "$RUN_DIR/results.json already exists: step e runs once"
[ "$(cfg acceptance.real)" != "null" ] || die "acceptance.real is empty in $CONFIG"
require_step_e_ready
if [ "${HBFGWAS_SKIP_PREREG_CHECK:-0}" != 1 ]; then
	local_tag=$(git rev-parse 'prereg-v1^{commit}')
	remote_tag=$(git ls-remote origin 'refs/tags/prereg-v1^{}' 'refs/tags/prereg-v1' |
		awk '{print $1}' | tail -n 1) || true
	[ -n "$remote_tag" ] || die "tag prereg-v1 is not on GitHub (origin): push it before step e"
	[ "$remote_tag" = "$local_tag" ] ||
		die "tag prereg-v1 points to $local_tag here but $remote_tag on GitHub"
fi
[ -s "$DATA_DIR/pca.eigenvec" ] || die "no $DATA_DIR/pca.eigenvec: run 03_pca.sh first"
if [ "$(cfg sensitivity.exclude_samples)" != "[]" ] && [ "$(cfg sensitivity.exclude_samples)" != "null" ]; then
	[ -s "$DATA_DIR/pca_sens.eigenvec" ] ||
		die "no $DATA_DIR/pca_sens.eigenvec: run 03d_sensitivity_pca.sh first"
fi
export PYTHON PLINK2 CONFIG DATA_DIR RUN_DIR

# --- integrity check, still before any phenotype ---------------------------------
# Recompute each locus's conditional power and the predicted lambda_GC from the
# genotypes (03c, primary analysis only) and stop if any differs from the value
# declared in acceptance.real. Not an acceptance criterion (D-057).
log "integrity check: recomputing the declared predictions from the genotypes (no phenotype)"
RUN_DIR="$RUN_DIR/integrity" PREREG_SKIP_SENSITIVITY=1 bash "$here/03c_prereg_predictions.sh"
"$PYTHON" -m hbfgwas.pipeline check-declared --config "$CONFIG" \
	--predictions "$RUN_DIR/integrity/prereg_predictions.json" \
	--out "$RUN_DIR/integrity_check.json" ||
	die "integrity check failed: the recomputed predictions differ from the declared values ($RUN_DIR/integrity_check.json); no phenotype was drawn"
log "integrity check passed; running 04, 04b, 04c, 05 and the gate into $RUN_DIR"

bash "$here/04_assoc.sh"
bash "$here/04b_power_redraws.sh"
if [ -s "$DATA_DIR/pca_sens.eigenvec" ]; then
	bash "$here/04c_sensitivity_scan.sh"
fi
report_status=0
"$PYTHON" "$here/05_report.py" || report_status=$?
[ "$report_status" = 0 ] || [ "$report_status" = 4 ] ||
	die "05_report.py failed (exit $report_status)"
gate_status=0
"$PYTHON" -m hbfgwas.gate --config "$CONFIG" --experiment real \
	--results "$RUN_DIR/results.json" --out "$RUN_DIR/gate.json" || gate_status=$?
case "$gate_status" in
0) ;;
1)
	log "gate: at least one pre-registered criterion failed; change nothing and report it"
	exit 3
	;;
*) die "the gate could not evaluate the results (exit $gate_status)" ;;
esac
[ "$report_status" = 0 ] || die "gate passed, but the Python-vs-plink2 concordance check failed"
log "step e done: every gated criterion passed ($RUN_DIR/gate.json)"
