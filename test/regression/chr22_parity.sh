#!/usr/bin/env bash
# test/regression/chr22_parity.sh
#
# Regression: bwa and bwa-mem3 produce identical SAM output on a ~50k PE
# holodeck simulation against hg38 chr22.
#
# Was: the "Index chr22 with bwa and align" + "Align chr22 with bwa-mem3"
# + "Compare chr22 bwa vs bwa-mem3 (parity)" steps inline in ci.yml.
#
# Inputs:
#   BWA_MEM3       — path to the bwa-mem3 binary under test
#   CHR22_FA       — path to chr22.fa (pre-indexed with bwa-mem3 by caller)
#   BWA_CHR22_FA   — path to a chr22.fa copy pre-indexed with bwa
#                    (sidecars: .amb .ann .bwt .pac .sa). Lives next to
#                    CHR22_FA so the bwa-mem3 sidecars don't collide.
#   CHR22_SIM_DIR  — directory containing holodeck-simulated
#                    reads.r1.fastq.gz / reads.r2.fastq.gz
#
# Optional:
#   BWA_CHR22_NORMALIZED_SAM — a precomputed upstream-bwa result in exactly
#                    the normalized, sorted form compared below (as written
#                    by CHR22_UPSTREAM_OUT mode). Upstream bwa's output never
#                    changes with bwa-mem3 commits, so CI computes it once and
#                    caches it (.github/actions/chr22-fixture) instead of
#                    re-running bwa on every row. It is used only if it passes
#                    a structural sanity check (one primary record per read);
#                    otherwise, or when unset, bwa is run here as before.
#   CHR22_UPSTREAM_OUT — producer mode: run only the upstream-bwa side, write
#                    its normalized form to this path and exit. BWA_MEM3 and
#                    CHR22_FA are not needed in this mode.
set -euo pipefail

: "${BWA_CHR22_FA:?BWA_CHR22_FA must be set}"
: "${CHR22_SIM_DIR:?CHR22_SIM_DIR must be set}"

r1="$CHR22_SIM_DIR/reads.r1.fastq.gz"
r2="$CHR22_SIM_DIR/reads.r2.fastq.gz"

# Drop the tags that legitimately differ and the @PG/@HD lines, then sort.
# LC_ALL=C makes the sort order locale-independent, so a cached upstream file
# compares correctly whatever locale the consuming step runs under.
normalize() {
    grep -v '^@PG' "$1" | grep -v '^@HD' | sed 's/\tMQ:i:[0-9]*//' | sed 's/\tHN:i:[0-9]*//' | sed 's/\tXS:i:[-0-9]*//' | LC_ALL=C sort
}

# Upstream bwa, normalized into $1. The bwa command line (including -t, which
# sets the batch size and so the per-batch insert-size estimate) is part of
# this script, so the cache key hashes this file.
run_upstream() {
    bwa mem -t 4 "$BWA_CHR22_FA" "$r1" "$r2" \
        > "$CHR22_SIM_DIR/bwa.sam" 2> "$CHR22_SIM_DIR/bwa.log"
    normalize "$CHR22_SIM_DIR/bwa.sam" > "$1"
}

# A precomputed upstream file is trusted only if it holds exactly one primary
# record (no 0x100 secondary / 0x800 supplementary bit) per input read: two
# per pair in r1. This catches an empty, truncated, or wrong-fixture file; it
# deliberately does not detect a record whose content differs, which instead
# surfaces as a parity FAIL.
upstream_is_valid() {
    local f=$1 pairs primaries
    [[ -s $f ]] || return 1
    grep -q '^@SQ' "$f" || return 1
    pairs=$(($(gzip -dc "$r1" | wc -l) / 4))
    primaries=$(grep -v '^@' "$f" | awk -F'\t' 'NF >= 11 && int($2 / 256) % 2 == 0 && int($2 / 2048) % 2 == 0' | wc -l)
    ((pairs > 0 && primaries == 2 * pairs))
}

if [[ -n ${CHR22_UPSTREAM_OUT:-} ]]; then
    run_upstream "$CHR22_UPSTREAM_OUT"
    if ! upstream_is_valid "$CHR22_UPSTREAM_OUT"; then
        echo "FAIL: upstream bwa output failed its sanity check" >&2
        exit 1
    fi
    echo "upstream bwa records: $(grep -cv '^@' "$CHR22_UPSTREAM_OUT")"
    exit 0
fi

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${CHR22_FA:?CHR22_FA must be set}"

upstream="$CHR22_SIM_DIR/bwa.normalized.sam"
if [[ -n ${BWA_CHR22_NORMALIZED_SAM:-} ]] && upstream_is_valid "$BWA_CHR22_NORMALIZED_SAM"; then
    echo "using precomputed upstream bwa result: $BWA_CHR22_NORMALIZED_SAM"
    upstream=$BWA_CHR22_NORMALIZED_SAM
else
    if [[ -n ${BWA_CHR22_NORMALIZED_SAM:-} ]]; then
        echo "note: precomputed upstream bwa result missing or invalid; running bwa" >&2
    fi
    run_upstream "$upstream"
fi

"$BWA_MEM3" mem -t 4 "$CHR22_FA" "$r1" "$r2" \
    > "$CHR22_SIM_DIR/bwamem3.sam" 2> "$CHR22_SIM_DIR/bwamem3.log"
normalize "$CHR22_SIM_DIR/bwamem3.sam" > "$CHR22_SIM_DIR/bwamem3.normalized.sam"

echo "bwa records:      $(grep -cv '^@' "$upstream")"
echo "bwa-mem3 records: $(grep -cv '^@' "$CHR22_SIM_DIR/bwamem3.normalized.sam")"
if diff "$upstream" "$CHR22_SIM_DIR/bwamem3.normalized.sam" > /dev/null 2>&1; then
    echo "PASS: bwa == bwa-mem3 on chr22 (~50k PE holodeck reads)"
else
    echo "FAIL: chr22 parity diff; first 40 lines:"
    diff "$upstream" "$CHR22_SIM_DIR/bwamem3.normalized.sam" | head -40
    exit 1
fi
