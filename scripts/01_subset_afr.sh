#!/usr/bin/env bash
# 01: keep 1000 Genomes SuperPop AFR, remove the related-sample list, keep
# autosomal biallelic ACGT SNPs, and remove duplicate variants.
#
# Duplicates: every record at a chromosome:position that carries more than one
# biallelic SNP record (e.g. a split multiallelic site), and every record whose
# ID occurs more than once. Missing IDs are first set to chrom:pos:ref:alt.
# Fails unless the remaining sample count lies in qc.expected_afr_samples
# (640-665 in config/analysis.yaml).
#
# In:  DATA_DIR/all_phase3.{pgen,pvar.zst,psam}, DATA_DIR/deg2_phase3.king.cutoff.out.id
# Out: DATA_DIR/afr.{pgen,pvar,psam}; RUN_DIR/01_subset_counts.tsv; RUN_DIR/logs/
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=scripts/_common.sh
. "$here/_common.sh"
need "$PLINK2"
setup_dirs

in="$DATA_DIR/all_phase3"
rel="$DATA_DIR/deg2_phase3.king.cutoff.out.id"
min=$(cfg qc.expected_afr_samples.min)
max=$(cfg qc.expected_afr_samples.max)
counts="$RUN_DIR/01_subset_counts.tsv"
printf 'step\tsamples\tvariants\tnote\n' >"$counts"

# 1. AFR sample list from the SuperPop column (found by header name).
awk -F'\t' '
	NR == 1 {
		for (i = 1; i <= NF; i++) { h = $i; sub(/^#/, "", h); col[h] = i }
		if (!("IID" in col) || !("SuperPop" in col)) {
			print "psam lacks an IID or SuperPop column" > "/dev/stderr"; exit 2
		}
		print "#IID"; next
	}
	$col["SuperPop"] == "AFR" { print $col["IID"] }
' "$in.psam" >"$DATA_DIR/afr_keep.txt"
n_all=$(count_records "$in.psam")
n_afr=$(count_records "$DATA_DIR/afr_keep.txt")
n_rel_afr=$(awk 'NR == FNR { if (!/^#/) keep[$1] = 1; next } !/^#/ && ($NF in keep) { n++ } END { print n + 0 }' \
	"$DATA_DIR/afr_keep.txt" "$rel")
printf 'all_phase3\t%s\tNA\tall samples in the .psam\n' "$n_all" >>"$counts"
printf 'superpop_afr\t%s\tNA\tSuperPop == AFR\n' "$n_afr" >>"$counts"
log "AFR samples: $n_afr of $n_all; $n_rel_afr of them are on the related-sample list"

# 2. AFR minus related samples; autosomes; biallelic ACGT SNPs; fill missing IDs.
# --pfile does not find a .pvar.zst on its own, so name the three files.
# shellcheck disable=SC2016 # $r and $a are plink2 template fields, not shell variables
plink --pgen "$in.pgen" --pvar "$in.pvar.zst" --psam "$in.psam" \
	--keep "$DATA_DIR/afr_keep.txt" --remove "$rel" \
	--autosome --snps-only just-acgt --max-alleles 2 \
	--set-missing-var-ids '@:#:$r:$a' \
	--make-pgen --out "$DATA_DIR/afr_step1"
save_log "$DATA_DIR/afr_step1"
v_loaded=$(grep -Eo '[0-9]+ variants loaded' "$DATA_DIR/afr_step1.log" | head -n 1 | awk '{print $1}' || true)
printf 'loaded\tNA\t%s\tvariants in the input .pvar (from the plink2 log)\n' "${v_loaded:-NA}" >>"$counts"
printf 'afr_unrelated_autosomal_biallelic_acgt_snps\t%s\t%s\tminus related list; --autosome --snps-only just-acgt --max-alleles 2\n' \
	"$(count_records "$DATA_DIR/afr_step1.psam")" "$(count_records "$DATA_DIR/afr_step1.pvar")" >>"$counts"

# 3. Duplicates. Same-position records are adjacent in a sorted .pvar.
pvar1="$DATA_DIR/afr_step1.pvar"
awk -F'\t' '
	/^#/ { next }
	{
		key = $1 ":" $2
		if ($1 == chrom && $2 + 0 < pos + 0) { print "unsorted .pvar at " key > "/dev/stderr"; exit 3 }
		if (key == prev) { buf = buf "\n" $3; run++ }
		else { if (run > 1) print buf; buf = $3; run = 1 }
		prev = key; chrom = $1; pos = $2
	}
	END { if (run > 1) print buf }
' "$pvar1" >"$DATA_DIR/dup_position_ids.txt"
awk -F'\t' '!/^#/ { print $3 }' "$pvar1" | LC_ALL=C sort | uniq -d >"$DATA_DIR/dup_ids.txt"
LC_ALL=C sort -u "$DATA_DIR/dup_position_ids.txt" "$DATA_DIR/dup_ids.txt" >"$DATA_DIR/dup_exclude.txt"
n_pos=$(awk 'NF' "$DATA_DIR/dup_position_ids.txt" | wc -l | tr -d ' ')
n_ids=$(awk 'NF' "$DATA_DIR/dup_ids.txt" | wc -l | tr -d ' ')
log "duplicates: $n_pos records at shared positions; $n_ids IDs used more than once"

exclude=()
if [ -s "$DATA_DIR/dup_exclude.txt" ]; then exclude=(--exclude "$DATA_DIR/dup_exclude.txt"); fi
plink --pfile "$DATA_DIR/afr_step1" ${exclude[@]+"${exclude[@]}"} --make-pgen --out "$DATA_DIR/afr"
save_log "$DATA_DIR/afr"
n_final=$(count_records "$DATA_DIR/afr.psam")
v_final=$(count_records "$DATA_DIR/afr.pvar")
printf 'afr_dedup\t%s\t%s\tremoved %s records at shared positions and all records of %s duplicated IDs\n' \
	"$n_final" "$v_final" "$n_pos" "$n_ids" >>"$counts"

if [ "${KEEP_INTERMEDIATE:-0}" != 1 ]; then rm -f "$DATA_DIR"/afr_step1.{pgen,pvar,psam}; fi

if [ "$n_final" -lt "$min" ] || [ "$n_final" -gt "$max" ]; then
	die "$n_final AFR samples remain; the pre-registered range is $min-$max"
fi
log "AFR subset: $n_final samples (expected $min-$max), $v_final variants; counts in $counts"
