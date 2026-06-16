#!/usr/bin/env bash
# test/regression/longread_pe_parity.sh
#
# Regression: the recovered long-read 8-bit banded-SW path must produce
# byte-identical SAM to the 16-bit reference on LONG paired-end reads.
#
# chr22_parity.sh covers parity vs upstream bwa at 151 bp. This test targets the
# recovered >=128 bp 8-bit extension path specifically: it simulates 400 bp PE
# reads (whose seed-extension segments land in the 128-1088 bp range routed to
# the 8-bit kernel) and diffs `bwa-mem3 mem` with the 8-bit path ENABLED (default)
# against the same binary with BWAMEM3_DISABLE_BSW8=1 (every pair forced to the
# 16-bit, then scalar, path -- i.e. the pre-recovery behavior). The two must be
# byte-identical; any divergence is a recovered-8-bit-path regression.
#
# Self-contained A/B (one binary, only the routing gate flips), so it does not
# depend on the upstream bwa version's long-read behavior.
#
# Inputs (env vars):
#   BWA_MEM3       - path to the bwa-mem3 binary under test
#   CHR22_FA       - chr22.fa, pre-indexed with bwa-mem3 by the caller
#   CHR22_SIM_DIR  - directory for fixture-private intermediates
#   PIXI_MANIFEST  - path to test/holodeck/pixi.toml (holodeck + samtools)
set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${CHR22_FA:?CHR22_FA must be set}"
: "${CHR22_SIM_DIR:?CHR22_SIM_DIR must be set}"
: "${PIXI_MANIFEST:?PIXI_MANIFEST must be set}"

mkdir -p "$CHR22_SIM_DIR"

# 400 bp PE reads: long enough that seed-extension segments exceed the int8
# position ceiling (>127) and route to the recovered diagonal-offset 8-bit path,
# but well under MAX_SEQ_LEN8 (1088). ~0.3x of 50.8 Mb chr22 at 400 bp PE is a
# few tens of thousands of pairs -- enough coverage of the path for a regression
# gate without dominating CI wall time.
pfx="$CHR22_SIM_DIR/longread_pe"
pixi run --frozen --manifest-path "$PIXI_MANIFEST" -- \
  holodeck simulate \
    -r "$CHR22_FA" \
    -o "$pfx" \
    --coverage 0.3 \
    --read-length 400 \
    --fragment-mean 900 --fragment-stddev 100 \
    --seed 42 \
    --simple-names \
    --min-error-rate 0.002 --max-error-rate 0.002

r1="$pfx.r1.fastq.gz"
r2="$pfx.r2.fastq.gz"
[ -s "$r1" ] && [ -s "$r2" ] || { echo "FAIL: holodeck produced no PE fastqs ($r1, $r2)" >&2; exit 1; }

# Normalize: drop @PG (carries the command line) and @HD; sort so thread
# interleaving does not affect the comparison.
normalize() { grep -v '^@PG' "$1" | grep -v '^@HD' | sort; }

sam8="$CHR22_SIM_DIR/longread_pe.8bit.sam"
sam16="$CHR22_SIM_DIR/longread_pe.16bit.sam"

"$BWA_MEM3" mem -t 4 "$CHR22_FA" "$r1" "$r2" > "$sam8.raw" 2>"$CHR22_SIM_DIR/longread_pe.8bit.log"
BWAMEM3_DISABLE_BSW8=1 "$BWA_MEM3" mem -t 4 "$CHR22_FA" "$r1" "$r2" > "$sam16.raw" 2>"$CHR22_SIM_DIR/longread_pe.16bit.log"

normalize "$sam8.raw"  > "$sam8"
normalize "$sam16.raw" > "$sam16"

n8=$(grep -cv '^@' "$sam8" || true)
echo "long-read PE records: $n8 (read-length 400, chr22 ~0.3x)"

# Guard against a trivially-passing empty run.
if [ "$n8" -lt 1000 ]; then
    echo "FAIL: long-read PE SAM has $n8 records, expected tens of thousands (holodeck/align regression?)" >&2
    exit 1
fi

if diff -q "$sam8" "$sam16" > /dev/null 2>&1; then
    echo "PASS: recovered 8-bit path == 16-bit reference on 400 bp PE reads ($n8 records)"
else
    echo "FAIL: 8-bit vs 16-bit SAM diverged on long PE reads; first 40 diff lines:" >&2
    diff "$sam8" "$sam16" | head -40 >&2
    exit 1
fi
