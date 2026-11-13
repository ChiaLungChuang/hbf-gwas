# shellcheck shell=bash
# Shared helpers for scripts/0[0-4]_*.sh. Sourced, not executed.
# Written for bash 3.2 (macOS /bin/bash) as well as bash 5.
#
# Environment (all optional):
#   PYTHON      python with hbfgwas installed      (default: python3)
#   PLINK2      plink2 executable                  (default: plink2)
#   CONFIG      pre-registration file              (default: config/analysis.yaml)
#   DATA_DIR    genotype files; never committed    (default: data/1kg)
#   RUN_DIR     run documentation; committed       (default: runs/<UTC date>)
#   PLINK_THREADS, PLINK_MEMORY_MB                 passed to plink2 if set

: "${PYTHON:=python3}"
: "${PLINK2:=plink2}"
: "${CONFIG:=config/analysis.yaml}"
: "${DATA_DIR:=data/1kg}"
: "${RUN_DIR:=runs/$(date -u +%Y-%m-%d)}"

log() { printf '[%s] %s\n' "$(date -u +%H:%M:%S)" "$*" >&2; }
die() {
	printf 'ERROR: %s\n' "$*" >&2
	exit 1
}
need() { command -v "$1" >/dev/null 2>&1 || die "'$1' not found on PATH"; }

# Make DATA_DIR and RUN_DIR absolute so later cd's cannot break them.
setup_dirs() {
	mkdir -p "$DATA_DIR" "$RUN_DIR/logs"
	DATA_DIR=$(cd "$DATA_DIR" && pwd)
	RUN_DIR=$(cd "$RUN_DIR" && pwd)
}

# Value from the pre-registration, e.g. cfg qc.maf
cfg() { "$PYTHON" -m hbfgwas.pipeline config-get --config "$CONFIG" "$1"; }

sha256_of() {
	if command -v sha256sum >/dev/null 2>&1; then
		sha256sum "$1" | awk '{print $1}'
	else
		shasum -a 256 "$1" | awk '{print $1}'
	fi
}

file_size() { wc -c <"$1" | tr -d ' '; }

# Lines not starting with '#': samples in a .psam/.id file, variants in a .pvar,
# entries in a .snplist.
count_records() { awk '!/^#/ {n++} END {print n + 0}' "$1"; }

plink() {
	local extra=()
	if [ -n "${PLINK_THREADS:-}" ]; then extra+=(--threads "$PLINK_THREADS"); fi
	if [ -n "${PLINK_MEMORY_MB:-}" ]; then extra+=(--memory "$PLINK_MEMORY_MB"); fi
	"$PLINK2" ${extra[@]+"${extra[@]}"} "$@"
}

record_plink_version() { "$PLINK2" --version >"$RUN_DIR/plink2_version.txt"; }

# Refuse to run unless plink2 is the build pinned in the config (review item G):
# the version string must contain plink2.version and "(plink2.build_date)".
check_plink_version() {
	local want_v want_d got
	want_v=$(cfg plink2.version)
	want_d=$(cfg plink2.build_date)
	got=$("$PLINK2" --version 2>/dev/null | head -n 1 || true)
	case "$got" in
	*"$want_v"*"($want_d)"*) log "plink2: $got (pinned $want_v, $want_d)" ;;
	*)
		if [ "${PLINK2_ALLOW_UNPINNED:-0}" = 1 ]; then
			log "WARNING: plink2 is '$got', not the pinned $want_v ($want_d); PLINK2_ALLOW_UNPINNED=1"
		else
			die "plink2 is '$got'; config/analysis.yaml pins $want_v ($want_d)"
		fi
		;;
	esac
}

# Preconditions for step e (04_assoc.sh and step_e.sh), checked before any
# phenotype is simulated. Sets $loci (the locus table used) and $seed_file.
#   * the locus table passes the analysis check (signed off by the user);
#   * the tag prereg-v1 exists, $CONFIG is unchanged since it, and
#     acceptance.real is filled in (unless HBFGWAS_SKIP_PREREG_CHECK=1, which is
#     for the synthetic test fileset only);
#   * the step e seed is filled and the population shift is declared.
require_step_e_ready() {
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
}

# Copy the plink2 log for output prefix $1 into RUN_DIR/logs/.
save_log() { if [ -f "$1.log" ]; then cp "$1.log" "$RUN_DIR/logs/"; fi; }
