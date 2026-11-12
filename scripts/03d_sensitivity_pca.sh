#!/usr/bin/env bash
# 03d (Phase 2 step d): PCs for the sensitivity analysis. Repeats the exact
# PCA of 03_pca.sh on the same LD-pruned variants and KING-unrelated samples,
# minus config sensitivity.exclude_samples (NA20314). No phenotype involved.
#
# In:  DATA_DIR/afr_qc.*, DATA_DIR/unrelated.id, DATA_DIR/pca.1_prune.prune.in
# Out: DATA_DIR/pca_sens.eigenvec; RUN_DIR/sens_pca_eigenval.txt, sens_pca_counts.tsv, logs/
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=scripts/_common.sh
. "$here/_common.sh"
need "$PLINK2"
setup_dirs
check_plink_version

in="$DATA_DIR/afr_qc"
n_pcs=$(cfg pca.n_pcs)
"$PYTHON" -c 'import json,sys; ids=json.loads(sys.argv[1]) or []; print("#IID"); print("\n".join(ids))' \
	"$(cfg sensitivity.exclude_samples)" >"$DATA_DIR/sens_exclude.id"
n_excl=$(count_records "$DATA_DIR/sens_exclude.id")
[ "$n_excl" -gt 0 ] || die "config sensitivity.exclude_samples is empty"

p="$DATA_DIR/pca_sens"
plink --pfile "$in" --keep "$DATA_DIR/unrelated.id" --remove "$DATA_DIR/sens_exclude.id" \
	--extract "$DATA_DIR/pca.1_prune.prune.in" --pca "$n_pcs" --out "$p.3"
save_log "$p.3"
cp "$p.3.eigenvec" "$DATA_DIR/pca_sens.eigenvec"
cp "$p.3.eigenval" "$RUN_DIR/sens_pca_eigenval.txt"
n_left=$(count_records "$DATA_DIR/pca_sens.eigenvec")
printf 'step\tsamples\tvariants\tnote\n' >"$RUN_DIR/sens_pca_counts.tsv"
printf 'sensitivity_pca\t%s\t%s\tunrelated minus %s; same LD-pruned variants as 03\n' \
	"$n_left" "$(count_records "$DATA_DIR/pca.1_prune.prune.in")" \
	"$(awk '!/^#/' "$DATA_DIR/sens_exclude.id" | paste -sd, -)" >>"$RUN_DIR/sens_pca_counts.tsv"
log "sensitivity PCA: $n_left samples, $n_pcs PCs; eigenvalues in $RUN_DIR/sens_pca_eigenval.txt"
