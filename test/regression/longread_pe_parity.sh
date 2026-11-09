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

# Checked explicitly rather than with ${VAR:?}, which exits with a bare shell
# error (and before any ERR trap runs), so the failure would carry no FAIL: line.
for v in BWA_MEM3 CHR22_FA CHR22_SIM_DIR PIXI_MANIFEST; do
    if [[ -z ${!v:-} ]]; then
        echo "FAIL: $v must be set" >&2
        exit 1
    fi
done

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
[ -s "$r1" ] && [ -s "$r2" ] || {
    echo "FAIL: holodeck produced no PE fastqs ($r1, $r2)" >&2
    exit 1
}

# Normalize: drop @PG (carries the command line) and @HD; sort so thread
# interleaving does not affect the comparison.
normalize() { grep -v '^@PG' "$1" | grep -v '^@HD' | sort; }

sam8="$CHR22_SIM_DIR/longread_pe.8bit.sam"
sam16="$CHR22_SIM_DIR/longread_pe.16bit.sam"
pairs8="$CHR22_SIM_DIR/longread_pe.8bit.pairs.tsv"
pairs16="$CHR22_SIM_DIR/longread_pe.16bit.pairs.tsv"

# Pin both legs, even when the caller disabled the 8-bit tier. Remove stale dumps
# so a binary that no longer implements the hook cannot pass on an earlier run.
rm -f "$pairs8" "$pairs16"
BWAMEM3_DISABLE_BSW8=0 BWAMEM3_DUMP_PAIRS="$pairs8" "$BWA_MEM3" mem -t 4 "$CHR22_FA" "$r1" "$r2" > "$sam8.raw" 2> "$CHR22_SIM_DIR/longread_pe.8bit.log"
BWAMEM3_DISABLE_BSW8=1 BWAMEM3_DUMP_PAIRS="$pairs16" "$BWA_MEM3" mem -t 4 "$CHR22_FA" "$r1" "$r2" > "$sam16.raw" 2> "$CHR22_SIM_DIR/longread_pe.16bit.log"

for pairs in "$pairs8" "$pairs16"; do
    if [[ ! -s $pairs ]]; then
        echo "FAIL: extension routing dump is missing or empty: $pairs" >&2
        exit 1
    fi
    if ! mawk -F '\t' '
        NR == 1 { if ($0 != "len1\tlen2\th0\tw\ttier") exit 1; next }
        {
            if (NF != 5) exit 1
            for (i = 1; i <= 5; i++) if ($i !~ /^-?[0-9]+$/) exit 1
            if ($1 < 0 || $2 < 0 || $4 < 0 || ($5 != 1 && $5 != 8 && $5 != 16)) exit 1
            n++
        }
        END { if (n == 0) exit 1 }
    ' "$pairs"; then
        echo "FAIL: extension routing dump has invalid rows or no data: $pairs" >&2
        exit 1
    fi
done
# The dump records the actual tier decision; len2 is the query extension length.
nlong8=$(mawk 'NR > 1 && $2 >= 128 && $5 == 8 { n++ } END { print n+0 }' "$pairs8")
nref8=$(mawk 'NR > 1 && $5 == 8 { n++ } END { print n+0 }' "$pairs16")
if [[ $nlong8 -eq 0 || $nref8 -ne 0 ]]; then
    echo "FAIL: expected long extensions in the 8-bit tier only in the enabled run (long enabled=$nlong8, all disabled=$nref8)" >&2
    exit 1
fi
echo "long-read PE routing: $nlong8 query extensions >=128 bp in the 8-bit tier; disabled run has none"

normalize "$sam8.raw" > "$sam8"
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
