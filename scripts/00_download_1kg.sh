#!/usr/bin/env bash
# 00: download the 1000 Genomes phase 3 (GRCh37) PLINK2 fileset and record a manifest.
#
# The Dropbox links rotate. Re-check every URL against
# https://www.cog-genomics.org/plink/2.0/resources before running, and override
# any that changed with the environment variables below.
#
#   PGEN_URL PVAR_URL PSAM_URL REL_URL   download URLs (defaults below)
#   MIN_FREE_GB      required free disk in DATA_DIR       (default 20)
#   MIN_BYTES_LARGE  minimum size of the .pgen.zst/.pvar.zst (default 1048576 = 1 MB)
#   MIN_BYTES_SMALL  minimum size of the .psam and the related-sample list (default 50)
#
# The brief asks for a 1 MB minimum on every file; the .psam (2504 samples) and
# the related-sample list (14 samples) are far smaller than 1 MB by nature, so
# they get a small minimum plus content checks instead (docs/decisions.md D-015).
# Any file that looks like an HTML page (a dead link) is rejected.
#
# Writes DATA_DIR/all_phase3.{pgen.zst,pvar.zst,psam}, the related-sample list,
# then DATA_DIR/all_phase3.pgen via plink2 --zst-decompress, and
# RUN_DIR/download_manifest.tsv (file, URL, UTC date, bytes, sha256) and
# RUN_DIR/plink2_version.txt.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=scripts/_common.sh
. "$here/_common.sh"

: "${PGEN_URL:=https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1}"
: "${PVAR_URL:=https://www.dropbox.com/s/c95n8quqwqww4s0/all_phase3_noannot.pvar.zst?dl=1}"
: "${PSAM_URL:=https://www.dropbox.com/scl/fi/haqvrumpuzfutklstazwk/phase3_corrected.psam?rlkey=0yyifzj2fb863ddbmsv4jkeq6&dl=1}"
: "${REL_URL:=https://www.dropbox.com/s/zj8d14vv9mp6x3c/deg2_phase3.king.cutoff.out.id?dl=1}"
: "${MIN_FREE_GB:=20}"
: "${MIN_BYTES_LARGE:=1048576}"
: "${MIN_BYTES_SMALL:=50}"

need curl
need "$PLINK2"
setup_dirs

free_kb=$(df -Pk "$DATA_DIR" | awk 'NR == 2 {print $4}')
need_kb=$((MIN_FREE_GB * 1024 * 1024))
if [ "$free_kb" -lt "$need_kb" ]; then
	die "only $((free_kb / 1024 / 1024)) GB free in $DATA_DIR; need >= $MIN_FREE_GB GB"
fi
log "free disk in $DATA_DIR: $((free_kb / 1024 / 1024)) GB (need >= $MIN_FREE_GB GB)"

manifest="$RUN_DIR/download_manifest.tsv"
printf 'file\turl\tdownloaded_utc\tbytes\tsha256\n' >"$manifest"

# fetch URL DEST MIN_BYTES KIND   (KIND: zst | psam | idlist)
fetch() {
	local url=$1 dest=$2 min=$3 kind=$4
	local out="$DATA_DIR/$dest"
	log "downloading $dest"
	curl -L --fail --retry 3 --retry-delay 5 --silent --show-error -o "$out.part" "$url"
	mv "$out.part" "$out"
	local bytes
	bytes=$(file_size "$out")
	if [ "$bytes" -lt "$min" ]; then
		die "$dest is only $bytes bytes (minimum $min). A dead or rotated link usually returns a small HTML page: re-check the URL at https://www.cog-genomics.org/plink/2.0/resources"
	fi
	if head -c 1024 "$out" | LC_ALL=C grep -qi -e '<!doctype html' -e '<html'; then
		die "$dest looks like an HTML page, not data: re-check the URL"
	fi
	case $kind in
	zst)
		if [ "$(head -c 4 "$out" | od -An -tx1 | tr -d ' \n')" != "28b52ffd" ]; then
			die "$dest is not zstd-compressed (bad magic bytes)"
		fi
		;;
	psam | idlist)
		if ! head -n 1 "$out" | grep -q '^#'; then
			die "$dest has no '#' header line; expected a plink2 sample file"
		fi
		;;
	esac
	printf '%s\t%s\t%s\t%s\t%s\n' "$dest" "$url" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
		"$bytes" "$(sha256_of "$out")" >>"$manifest"
	log "  $dest: $bytes bytes"
}

fetch "$PGEN_URL" all_phase3.pgen.zst "$MIN_BYTES_LARGE" zst
fetch "$PVAR_URL" all_phase3.pvar.zst "$MIN_BYTES_LARGE" zst
fetch "$PSAM_URL" all_phase3.psam "$MIN_BYTES_SMALL" psam
fetch "$REL_URL" deg2_phase3.king.cutoff.out.id "$MIN_BYTES_SMALL" idlist
cp "$manifest" "$DATA_DIR/download_manifest.tsv"

record_plink_version
log "plink2: $(cat "$RUN_DIR/plink2_version.txt")"
log "decompressing all_phase3.pgen.zst"
(cd "$DATA_DIR" && "$PLINK2" --zst-decompress all_phase3.pgen.zst all_phase3.pgen)
[ -s "$DATA_DIR/all_phase3.pgen" ] || die "decompression produced no all_phase3.pgen"
log "done: manifest at $manifest"
