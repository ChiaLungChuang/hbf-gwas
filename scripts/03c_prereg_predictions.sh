#!/usr/bin/env bash
# 03c (Phase 2 step d): predictions from the real genotypes and PCs, computed
# BEFORE any phenotype is drawn, for the real-data pre-registration PR:
#   * conditional power per testable sentinel (PC-adjusted and unadjusted), and
#     which loci the sim-v3 rule makes REQUIRED (predicted power >= 0.99);
#   * predicted lambda_GC, unadjusted and PC-adjusted, over random null variants
#     more than association.null_exclusion_window_bp from every sentinel
#     (median and 5th-95th percentile over draws of the population shift).
# A sentinel absent from the QC'd callset is "not testable in this callset".
#
#   HBFGWAS_LOCI_TSV=<file>   use this locus table instead of exporting loci.py
#                             (synthetic test fileset); it must still be analysis_ready
#   PREREG_SKIP_SENSITIVITY=1 skip the sensitivity predictions (step e's integrity check)
#
# In:  DATA_DIR/afr_qc.*, DATA_DIR/unrelated.id, DATA_DIR/pca.eigenvec
# Out: RUN_DIR/prereg_predictions.json (+ _sensitivity.json), prereg_sentinels.tsv, logs/
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=scripts/_common.sh
. "$here/_common.sh"
need "$PLINK2"
setup_dirs
capture_logs "$@"
check_plink_version

in="$DATA_DIR/afr_qc"
keep="$DATA_DIR/unrelated.id"
loci="${HBFGWAS_LOCI_TSV:-}"
if [ -n "$loci" ]; then
	log "WARNING: HBFGWAS_LOCI_TSV is set; using $loci instead of exporting src/hbfgwas/loci.py"
	"$PYTHON" -m hbfgwas.pipeline loci-ready --loci "$loci" ||
		die "the locus table $loci is not signed off (analysis_ready is false for a locus)"
else
	loci="$DATA_DIR/loci.tsv"
	"$PYTHON" -m hbfgwas.loci export --require-analysis-ready --out "$loci" ||
		die "the locus table is not signed off by the user (python -m hbfgwas.loci check --analysis)"
fi

"$PYTHON" -m hbfgwas.pipeline sentinels --pvar "$in.pvar" --loci "$loci" \
	--out-ids "$DATA_DIR/prereg_ids.txt" --out-map "$DATA_DIR/prereg_map.tsv"
redact_to "$DATA_DIR/prereg_map.tsv" "$RUN_DIR/prereg_sentinels.tsv"
[ -s "$DATA_DIR/prereg_ids.txt" ] || die "no sentinel is testable in this callset"
plink --pfile "$in" --keep "$keep" --extract "$DATA_DIR/prereg_ids.txt" --export A \
	--out "$DATA_DIR/prereg_sentinels"
save_log "$DATA_DIR/prereg_sentinels"

"$PYTHON" -m hbfgwas.pipeline sample-ids --pvar "$in.pvar" \
	--n "$(cfg prediction.n_null_variants)" --seed "$(cfg seeds.real_prediction)" \
	--exclude-around "$DATA_DIR/prereg_map.tsv" \
	--window-bp "$(cfg association.null_exclusion_window_bp)" \
	--out "$DATA_DIR/prereg_null_ids.txt"
plink --pfile "$in" --keep "$keep" --extract "$DATA_DIR/prereg_null_ids.txt" --export A \
	--out "$DATA_DIR/prereg_null"
save_log "$DATA_DIR/prereg_null"

"$PYTHON" -m hbfgwas.pipeline predictions --config "$CONFIG" \
	--raw-sentinels "$DATA_DIR/prereg_sentinels.raw" --map "$DATA_DIR/prereg_map.tsv" \
	--raw-null "$DATA_DIR/prereg_null.raw" --psam "$in.psam" --eigenvec "$DATA_DIR/pca.eigenvec" \
	--out "$RUN_DIR/prereg_predictions.json"
log "predictions (no phenotype drawn): $RUN_DIR/prereg_predictions.json"

# NA20314 sensitivity analysis (not gated): the same predictions without the
# excluded samples, with the PCs that 03d_sensitivity_pca.sh recomputed.
excl=$("$PYTHON" -c 'import json,sys; print(",".join(json.loads(sys.argv[1]) or []))' \
	"$(cfg sensitivity.exclude_samples)")
if [ "${PREREG_SKIP_SENSITIVITY:-0}" = 1 ]; then
	excl=""  # step e's integrity check recomputes the primary predictions only
fi
if [ -n "$excl" ]; then
	[ -s "$DATA_DIR/pca_sens.eigenvec" ] || die "run scripts/03d_sensitivity_pca.sh first"
	"$PYTHON" -m hbfgwas.pipeline predictions --config "$CONFIG" \
		--raw-sentinels "$DATA_DIR/prereg_sentinels.raw" --map "$DATA_DIR/prereg_map.tsv" \
		--raw-null "$DATA_DIR/prereg_null.raw" --psam "$in.psam" \
		--eigenvec "$DATA_DIR/pca_sens.eigenvec" --exclude "$excl" \
		--out "$RUN_DIR/prereg_predictions_sensitivity.json"
	log "sensitivity predictions without $excl: $RUN_DIR/prereg_predictions_sensitivity.json"
fi
