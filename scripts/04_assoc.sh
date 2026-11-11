#!/usr/bin/env bash
# 04: simulate an HbF-like phenotype on the real genotypes, then plink2 --glm
# without and with the PCs. The phenotype is SIMULATED by hbfgwas (planted
# effects at the sentinel variants + a shift by 1000 Genomes population +
# noise); 1000 Genomes has no HbF phenotype.
#
# Refuses to run unless:
#   * every value in the locus table is verified by the user (src/hbfgwas/loci.py), and
#   * the tag prereg-v1 exists, config/analysis.yaml is unchanged since it, and
#     acceptance.real is filled in (Phase 2 step d).
# Overrides, for the synthetic test fileset only:
#   HBFGWAS_LOCI_TSV=<file>       use this locus table, skip the sign-off check
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

in="$DATA_DIR/afr_qc"
keep="$DATA_DIR/unrelated.id"
n_pcs=$(cfg pca.n_pcs)

# --- preconditions ---------------------------------------------------------
loci="${HBFGWAS_LOCI_TSV:-}"
if [ -n "$loci" ]; then
	log "WARNING: HBFGWAS_LOCI_TSV is set; using $loci without the sign-off check"
else
	loci="$DATA_DIR/loci.tsv"
	"$PYTHON" -m hbfgwas.loci export --require-verified --out "$loci" ||
		die "the locus table is not signed off by the user (python -m hbfgwas.loci check)"
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

# --- planted variants: sentinels present after QC ----------------------------
"$PYTHON" -m hbfgwas.pipeline sentinels --pvar "$in.pvar" --loci "$loci" \
	--out-ids "$DATA_DIR/planted_ids.txt" --out-map "$DATA_DIR/planted_map.tsv"
cp "$DATA_DIR/planted_map.tsv" "$RUN_DIR/planted.tsv"
[ -s "$DATA_DIR/planted_ids.txt" ] || die "no sentinel variant survived QC; nothing to plant"
plink --pfile "$in" --keep "$keep" --extract "$DATA_DIR/planted_ids.txt" --export A --out "$DATA_DIR/planted"
save_log "$DATA_DIR/planted"

"$PYTHON" -m hbfgwas.pipeline simulate-phenotype --config "$CONFIG" \
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
	--n "$(cfg concordance.n_random_variants)" --seed "$(cfg seeds.real_concordance_sample)" \
	--include "$DATA_DIR/planted_ids.txt" --out "$DATA_DIR/concord_ids.txt"
plink --pfile "$in" --keep "$keep" --extract "$DATA_DIR/concord_ids.txt" --export A --out "$DATA_DIR/concord"
save_log "$DATA_DIR/concord"
log "association done: $DATA_DIR/glm_nopc.HBF_SIM.glm.linear, $DATA_DIR/glm_pc.HBF_SIM.glm.linear"
