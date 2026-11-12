#!/usr/bin/env bash
# 04c: sensitivity scan without the samples in sensitivity.exclude_samples
# (NA20314; D-040). Reported, never gated.
# The phenotype is SIMULATED by hbfgwas; 1000 Genomes has no HbF phenotype.
#
# Uses the step e phenotype as drawn by 04_assoc.sh (not redrawn), drops the
# excluded samples, and runs plink2 --glm without covariates and with the PCs
# recomputed without them (03d_sensitivity_pca.sh).
#
# In:  DATA_DIR/afr_qc.*, unrelated.id, pheno.tsv (04), pca_sens.eigenvec and
#      sens_exclude.id (03d)
# Out: DATA_DIR/glm_sens_nopc.HBF_SIM.glm.linear, glm_sens_pc.HBF_SIM.glm.linear;
#      RUN_DIR/logs/
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=scripts/_common.sh
. "$here/_common.sh"
need "$PLINK2"
setup_dirs
capture_logs "$@"
check_plink_version

in="$DATA_DIR/afr_qc"
n_pcs=$(cfg pca.n_pcs)
[ -s "$DATA_DIR/pheno.tsv" ] || die "no $DATA_DIR/pheno.tsv: run 04_assoc.sh first"
[ -s "$DATA_DIR/pca_sens.eigenvec" ] || die "no $DATA_DIR/pca_sens.eigenvec: run 03d_sensitivity_pca.sh first"
[ -s "$DATA_DIR/sens_exclude.id" ] || die "no $DATA_DIR/sens_exclude.id: run 03d_sensitivity_pca.sh first"

plink --pfile "$in" --keep "$DATA_DIR/unrelated.id" --remove "$DATA_DIR/sens_exclude.id" \
	--pheno "$DATA_DIR/pheno.tsv" --pheno-name HBF_SIM \
	--glm allow-no-covars --out "$DATA_DIR/glm_sens_nopc"
save_log "$DATA_DIR/glm_sens_nopc"
plink --pfile "$in" --keep "$DATA_DIR/unrelated.id" --remove "$DATA_DIR/sens_exclude.id" \
	--pheno "$DATA_DIR/pheno.tsv" --pheno-name HBF_SIM \
	--covar "$DATA_DIR/pca_sens.eigenvec" --covar-name "PC1-PC$n_pcs" --covar-variance-standardize \
	--glm hide-covar --out "$DATA_DIR/glm_sens_pc"
save_log "$DATA_DIR/glm_sens_pc"
log "sensitivity scan done (without $(awk '!/^#/' "$DATA_DIR/sens_exclude.id" | paste -sd, -))"
