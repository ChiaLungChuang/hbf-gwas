#!/usr/bin/env bash
# 04b: empirical power at the planted sentinels over phenotype redraws (step e).
# The phenotype is SIMULATED by hbfgwas; 1000 Genomes has no HbF phenotype.
#
# Each redraw is the step e phenotype model with new noise: the same planted
# effects on the same genotypes and the same fixed population shift; the noise
# comes from the step e seed's child "power_reps". Every redraw is tested at
# the sentinels only, with the PCs (the gated power calibration) and without
# them (reported). The number of redraws is
# acceptance.real.power_calibration.n_redraws (2000 in the proposal).
#
# In:  DATA_DIR/planted.raw and planted_map.tsv (from 04_assoc.sh),
#      DATA_DIR/afr_qc.psam, DATA_DIR/pca.eigenvec
# Out: RUN_DIR/power_redraws.json
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=scripts/_common.sh
. "$here/_common.sh"
setup_dirs
capture_logs "$@"

seed_file="${STEP_E_SEED_FILE:-$(cfg real_data_simulation.seed_file)}"
[ -s "$DATA_DIR/planted.raw" ] || die "no $DATA_DIR/planted.raw: run 04_assoc.sh first"
"$PYTHON" -m hbfgwas.pipeline seed-child --seed-file "$seed_file" --name power_reps >/dev/null ||
	die "the step e seed in $seed_file is empty: it is filled only after the prereg-v1 tag"

"$PYTHON" -m hbfgwas.pipeline power-redraws --config "$CONFIG" --seed-file "$seed_file" \
	--raw "$DATA_DIR/planted.raw" --map "$DATA_DIR/planted_map.tsv" \
	--psam "$DATA_DIR/afr_qc.psam" --eigenvec "$DATA_DIR/pca.eigenvec" \
	--out "$RUN_DIR/power_redraws.json"
