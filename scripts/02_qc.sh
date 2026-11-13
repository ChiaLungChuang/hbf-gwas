#!/usr/bin/env bash
# 02: sample and variant QC, one filter at a time, in this order:
#   --geno, --mind, --maf, --hwe (midp keep-fewhet), thresholds from config/analysis.yaml (qc).
# Each filter is computed on the output of the previous one and its counts in
# and out are written to RUN_DIR/qc_filters.tsv. 1000 Genomes phase 3 is fully
# called (no missing genotype calls), so --geno and --mind are expected to
# remove nothing; the script says so explicitly.
#
# Also reports each sentinel variant (src/hbfgwas/loci.py, looked up by rsID and
# by GRCh37 position): its allele frequency in this AFR sample and the filter,
# if any, that removed it -> RUN_DIR/sentinels.tsv.
#
#   HBFGWAS_LOCI_TSV  use this locus table instead of exporting src/hbfgwas/loci.py
#
# In:  DATA_DIR/afr.{pgen,pvar,psam}
# Out: DATA_DIR/afr_qc.{pgen,pvar,psam}; RUN_DIR/qc_filters.tsv, sentinels.tsv, logs/
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=scripts/_common.sh
. "$here/_common.sh"
need "$PLINK2"
setup_dirs
capture_logs "$@"
check_plink_version

in="$DATA_DIR/afr"
geno=$(cfg qc.geno)
mind=$(cfg qc.mind)
maf=$(cfg qc.maf)
hwe=$(cfg qc.hwe_p)
hwe_mode=$(cfg qc.hwe_mode)
hwe_args=("$hwe" "$hwe_mode")
if [ "$(cfg qc.hwe_keep_fewhet)" = "true" ]; then hwe_args+=(keep-fewhet); fi
q="$DATA_DIR/qc"
table="$RUN_DIR/qc_filters.tsv"
printf 'step\tfilter\tsamples_in\tsamples_out\tsamples_removed\tvariants_in\tvariants_out\tvariants_removed\n' >"$table"
row() { printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$(($3 - $4))" "$5" "$6" "$(($5 - $6))" >>"$table"; }

s0=$(count_records "$in.psam")
v0=$(count_records "$in.pvar")
log "input: $s0 samples, $v0 variants"

plink --pfile "$in" --geno "$geno" --write-snplist --out "$q.1_geno"
save_log "$q.1_geno"
v1=$(count_records "$q.1_geno.snplist")
row 1 "--geno $geno" "$s0" "$s0" "$v0" "$v1"

plink --pfile "$in" --extract "$q.1_geno.snplist" --mind "$mind" --write-samples --out "$q.2_mind"
save_log "$q.2_mind"
s2=$(count_records "$q.2_mind.id")
row 2 "--mind $mind" "$s0" "$s2" "$v1" "$v1"

if [ "$((v0 - v1))" -eq 0 ] && [ "$((s0 - s2))" -eq 0 ]; then
	log "NOTE: --geno removed 0 variants and --mind removed 0 samples. This is expected:"
	log "      1000 Genomes phase 3 genotypes are fully called (no missing genotype calls),"
	log "      so the missingness filters cannot remove anything from this dataset."
else
	log "WARNING: missingness filters removed $((v0 - v1)) variants and $((s0 - s2)) samples."
	log "         Fully called phase 3 data should lose none: check the input files."
fi

plink --pfile "$in" --keep "$q.2_mind.id" --extract "$q.1_geno.snplist" --maf "$maf" \
	--write-snplist --out "$q.3_maf"
save_log "$q.3_maf"
v3=$(count_records "$q.3_maf.snplist")
row 3 "--maf $maf" "$s2" "$s2" "$v1" "$v3"

# keep-fewhet: only excess heterozygosity is filtered (Wahlund effect in pooled AFR)
plink --pfile "$in" --keep "$q.2_mind.id" --extract "$q.3_maf.snplist" --hwe "${hwe_args[@]}" \
	--write-snplist --out "$q.4_hwe"
save_log "$q.4_hwe"
v4=$(count_records "$q.4_hwe.snplist")
row 4 "--hwe ${hwe_args[*]}" "$s2" "$s2" "$v3" "$v4"

plink --pfile "$in" --keep "$q.2_mind.id" --extract "$q.4_hwe.snplist" --make-pgen --out "$DATA_DIR/afr_qc"
save_log "$DATA_DIR/afr_qc"
log "after QC: $(count_records "$DATA_DIR/afr_qc.psam") samples, $(count_records "$DATA_DIR/afr_qc.pvar") variants"

# Sentinels: frequency in this AFR sample (all samples passing --mind) and QC survival.
loci="${HBFGWAS_LOCI_TSV:-}"
if [ -z "$loci" ]; then
	loci="$DATA_DIR/loci.tsv"
	"$PYTHON" -m hbfgwas.loci export --out "$loci"
fi
"$PYTHON" -m hbfgwas.pipeline sentinels --pvar "$in.pvar" --loci "$loci" \
	--out-ids "$DATA_DIR/sentinel_ids.txt" --out-map "$DATA_DIR/sentinel_map.tsv"
if [ -s "$DATA_DIR/sentinel_ids.txt" ]; then
	plink --pfile "$in" --keep "$q.2_mind.id" --extract "$DATA_DIR/sentinel_ids.txt" --freq \
		--out "$DATA_DIR/sentinels"
	save_log "$DATA_DIR/sentinels"
	freq_arg=(--afreq "$DATA_DIR/sentinels.afreq")
else
	freq_arg=()
fi
"$PYTHON" -m hbfgwas.pipeline sentinel-report --map "$DATA_DIR/sentinel_map.tsv" \
	${freq_arg[@]+"${freq_arg[@]}"} \
	--step "geno=$q.1_geno.snplist" --step "maf=$q.3_maf.snplist" --step "hwe=$q.4_hwe.snplist" \
	--out "$RUN_DIR/sentinels.tsv"
log "QC table: $table; sentinel report: $RUN_DIR/sentinels.tsv"
