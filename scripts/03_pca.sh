#!/usr/bin/env bash
# 03: ancestry PCs. Exclude long-range LD regions (resources/, GRCh37), keep
# MAF >= pca.maf_for_pca, LD-prune, remove second-degree relatives with
# KING-robust (cutoff pca.king_cutoff), exact PCA (plink2 --pca N, not approx),
# and a --het report (no filter). Settings from config/analysis.yaml (pca).
#
# In:  DATA_DIR/afr_qc.{pgen,pvar,psam}
# Out: DATA_DIR/pca.eigenvec, DATA_DIR/unrelated.id (samples kept by KING);
#      RUN_DIR/pca_eigenval.txt, 03_pca_counts.tsv, het_summary.tsv, logs/
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=scripts/_common.sh
. "$here/_common.sh"
need "$PLINK2"
setup_dirs

in="$DATA_DIR/afr_qc"
n_pcs=$(cfg pca.n_pcs)
regions=$(cfg pca.long_range_ld_regions)
pca_maf=$(cfg pca.maf_for_pca)
win=$(cfg pca.ld_prune.window_kb)
step=$(cfg pca.ld_prune.step)
r2=$(cfg pca.ld_prune.r2)
king=$(cfg pca.king_cutoff)
p="$DATA_DIR/pca"
counts="$RUN_DIR/03_pca_counts.tsv"
printf 'step\tsamples\tvariants\tnote\n' >"$counts"
printf 'input\t%s\t%s\tafter QC\n' "$(count_records "$in.psam")" "$(count_records "$in.pvar")" >>"$counts"

# Region file -> plink2 bed1 (chrom, start, end, name; 1-based inclusive).
awk -F'\t' 'NR > 1 && NF >= 4 { print $1 "\t" $2 "\t" $3 "\t" $4 }' "$regions" >"$p.lrld.bed1"
log "long-range LD regions: $(wc -l <"$p.lrld.bed1" | tr -d ' ') (all flagged unverified in $regions)"

plink --pfile "$in" --exclude bed1 "$p.lrld.bed1" --maf "$pca_maf" --write-snplist --out "$p.0_candidates"
save_log "$p.0_candidates"
printf 'outside_long_range_ld_and_maf\t%s\t%s\t--exclude bed1 regions; --maf %s\n' \
	"$(count_records "$in.psam")" "$(count_records "$p.0_candidates.snplist")" "$pca_maf" >>"$counts"

plink --pfile "$in" --extract "$p.0_candidates.snplist" --indep-pairwise "${win}kb" "$step" "$r2" --out "$p.1_prune"
save_log "$p.1_prune"
printf 'ld_pruned\t%s\t%s\t--indep-pairwise %skb %s %s\n' \
	"$(count_records "$in.psam")" "$(count_records "$p.1_prune.prune.in")" "$win" "$step" "$r2" >>"$counts"

plink --pfile "$in" --extract "$p.1_prune.prune.in" --king-cutoff "$king" --out "$p.2_king"
save_log "$p.2_king"
cp "$p.2_king.king.cutoff.in.id" "$DATA_DIR/unrelated.id"
# plink2 may not write the .out.id file when nobody is removed
if [ ! -f "$p.2_king.king.cutoff.out.id" ]; then printf '#IID\n' >"$p.2_king.king.cutoff.out.id"; fi
n_unrel=$(count_records "$DATA_DIR/unrelated.id")
n_out=$(count_records "$p.2_king.king.cutoff.out.id")
printf 'king_unrelated\t%s\t%s\t--king-cutoff %s removed %s samples\n' \
	"$n_unrel" "$(count_records "$p.1_prune.prune.in")" "$king" "$n_out" >>"$counts"
cp "$p.2_king.king.cutoff.out.id" "$RUN_DIR/king_removed.id"

plink --pfile "$in" --keep "$DATA_DIR/unrelated.id" --extract "$p.1_prune.prune.in" \
	--pca "$n_pcs" --out "$p.3"
save_log "$p.3"
cp "$p.3.eigenvec" "$DATA_DIR/pca.eigenvec"
cp "$p.3.eigenval" "$RUN_DIR/pca_eigenval.txt"
[ "$(count_records "$p.3.eigenval")" -eq "$n_pcs" ] || die "expected $n_pcs eigenvalues"

plink --pfile "$in" --keep "$DATA_DIR/unrelated.id" --extract "$p.1_prune.prune.in" --het --out "$p.4"
save_log "$p.4"
"$PYTHON" -m hbfgwas.pipeline het-summary --het "$p.4.het" --out "$RUN_DIR/het_summary.tsv"
log "PCA: $n_unrel unrelated samples, $n_pcs PCs; eigenvalues in $RUN_DIR/pca_eigenval.txt"
