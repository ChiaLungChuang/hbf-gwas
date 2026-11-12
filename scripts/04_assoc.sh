#!/usr/bin/env bash
# 04: simulate an HbF-like phenotype on the real genotypes, then plink2 --glm
# without and with the PCs. The phenotype is SIMULATED by hbfgwas (planted
# effects at the sentinel variants + a shift by 1000 Genomes population +
# noise); 1000 Genomes has no HbF phenotype.
#
# Refuses to run unless:
#   * the locus table passes the analysis check (src/hbfgwas/loci.py, signed off by the user),
#   * the step e seed (config/step_e_seed.yaml) is filled and the population shift is declared, and
#   * the tag prereg-v1 exists, config/analysis.yaml is unchanged since it, and
#     acceptance.real is filled in (Phase 2 step d).
# Overrides, for the synthetic test fileset only:
#   HBFGWAS_LOCI_TSV=<file>       use this locus table (it must still be analysis_ready)
#   HBFGWAS_SKIP_PREREG_CHECK=1   skip the prereg-v1 checks
#
# In:  DATA_DIR/afr_qc.*, DATA_DIR/pca.eigenvec, DATA_DIR/unrelated.id
# Out: DATA_DIR/pheno.tsv, glm_nopc.HBF_SIM.glm.linear, glm_pc.HBF_SIM.glm.linear,
#      concord.raw (subset for the Python-vs-plink2 check);
#      RUN_DIR/phenotype_meta.json, planted.tsv, logs/
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
n_pcs=$(cfg pca.n_pcs)

# --- preconditions ---------------------------------------------------------
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
if [ "${HBFGWAS_SKIP_PREREG_CHECK:-0}" = 1 ]; then
	log "WARNING: HBFGWAS_SKIP_PREREG_CHECK=1; not checking the prereg-v1 tag"
else
	git rev-parse -q --verify refs/tags/prereg-v1 >/dev/null ||
		die "tag prereg-v1 does not exist: real-data criteria must be pre-registered first"
	git diff --quiet prereg-v1 -- "$CONFIG" ||
		die "$CONFIG differs from the prereg-v1 tag"
	[ "$(cfg acceptance.real)" != "null" ] || die "acceptance.real is empty in $CONFIG"
fi

seed_file="${STEP_E_SEED_FILE:-$(cfg real_data_simulation.seed_file)}"
"$PYTHON" -m hbfgwas.pipeline seed-child --seed-file "$seed_file" --name phenotype >/dev/null ||
	die "the step e seed in $seed_file is empty: it is filled only after the prereg-v1 tag"
[ "$(cfg real_data_simulation.population_shift)" != "null" ] ||
	die "real_data_simulation.population_shift is not pre-registered in $CONFIG"

# --- planted variants: sentinels present after QC ----------------------------
# A sentinel absent from the callset (by rsID and by position) or removed by QC
# is "not testable in this callset": it is not planted and no proxy is chosen
# (review item C, docs/decisions.md D-028).
"$PYTHON" -m hbfgwas.pipeline sentinels --pvar "$in.pvar" --loci "$loci" \
	--out-ids "$DATA_DIR/planted_ids.txt" --out-map "$DATA_DIR/planted_map.tsv"
redact_to "$DATA_DIR/planted_map.tsv" "$RUN_DIR/planted.tsv"
[ -s "$DATA_DIR/planted_ids.txt" ] || die "no sentinel variant survived QC; nothing to plant"
plink --pfile "$in" --keep "$keep" --extract "$DATA_DIR/planted_ids.txt" --export A --out "$DATA_DIR/planted"
save_log "$DATA_DIR/planted"

"$PYTHON" -m hbfgwas.pipeline simulate-phenotype --config "$CONFIG" --seed-file "$seed_file" \
	--raw "$DATA_DIR/planted.raw" --map "$DATA_DIR/planted_map.tsv" --psam "$in.psam" \
	--out "$DATA_DIR/pheno.tsv" --meta "$RUN_DIR/phenotype_meta.json"

# --- genome-wide association with plink2 -------------------------------------
plink --pfile "$in" --keep "$keep" --pheno "$DATA_DIR/pheno.tsv" --pheno-name HBF_SIM \
	--glm allow-no-covars --out "$DATA_DIR/glm_nopc"
save_log "$DATA_DIR/glm_nopc"
plink --pfile "$in" --keep "$keep" --pheno "$DATA_DIR/pheno.tsv" --pheno-name HBF_SIM \
	--covar "$DATA_DIR/pca.eigenvec" --covar-name "PC1-PC$n_pcs" --covar-variance-standardize \
	--glm hide-covar --out "$DATA_DIR/glm_pc"
save_log "$DATA_DIR/glm_pc"

# --- subset for the Python reference check (05_report.py) --------------------
"$PYTHON" -m hbfgwas.pipeline sample-ids --pvar "$in.pvar" \
	--n "$(cfg concordance.n_random_variants)" \
	--seed "$("$PYTHON" -m hbfgwas.pipeline seed-child --seed-file "$seed_file" --name concordance)" \
	--include "$DATA_DIR/planted_ids.txt" --out "$DATA_DIR/concord_ids.txt"
plink --pfile "$in" --keep "$keep" --extract "$DATA_DIR/concord_ids.txt" --export A --out "$DATA_DIR/concord"
save_log "$DATA_DIR/concord"
log "association done: $DATA_DIR/glm_nopc.HBF_SIM.glm.linear, $DATA_DIR/glm_pc.HBF_SIM.glm.linear"
