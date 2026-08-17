#!/usr/bin/env bash
# test/regression/chain_flt_empty_compat.sh
#
# Regression: when `mem_chain_flt`'s weight filter drops EVERY chain for a
# read, the two upstreams disagree about the read, and each `--compat` target
# must reproduce its own (fg-labs/bwa-mem3#310).
#
# bwa returns zero survivors -- its tail loops are bounded by `n_chn == 0`, so
# the read goes out unmapped. bwa-mem2's seqid-range machinery instead
# synthesizes a count of 1 for the emptied array, the unconditional
# `kept[0] = 3` becomes load-bearing, and the chain the filter just rejected is
# extended into an alignment. bwa-mem3 inherited bwa-mem2's behavior at the
# fork point and keeps it on the default path, because that path is the
# drop-in; `--compat=bwa-mem` selects bwa's.
#
# Reachable only when `min_chain_weight > 0`, which is never the default: `-W`,
# or the `-x pacbio`/`pbref`/`ont2d` presets. This script uses `-W` because it
# makes the trigger exact rather than probabilistic -- with `-W` above the read
# length NO chain can clear the threshold, so every read takes the path. (On
# real long reads the presets essentially never trigger it: measured over 500
# HiFi reads at `-x pacbio`, zero reads had every chain dropped, which is why
# this is pinned with `-W` and not with a preset.)
#
# The three assertions:
#   1. default        -- every read mapped (bwa-mem2 behavior, drop-in intact)
#   2. --compat=bwa-mem2 -- every read mapped (same; the target IS bwa-mem2)
#   3. --compat=bwa-mem  -- every read UNMAPPED (bwa's behavior)
#
# Guard against a vacuous pass: the same reads at `-W 0` must map under all
# three, proving the fixture aligns at all and that assertion 3 is the weight
# filter firing rather than reads that simply do not align.
#
# Fixture (deterministic, no PRNG): 150 bp SE reads sliced from the committed
# phix.fa, same construction as compat_byte_identical.sh. `-W 200` exceeds the
# read length, so every chain's weight is below the threshold.
#
# Inputs (env vars):
#   BWA_MEM3           — path to the bwa-mem3 binary under test
#   CHAIN_FLT_PHIX_FA  — path to test/fixtures/phix.fa (the reference source)
#   CHAIN_FLT_WORK_DIR — directory for fixture-private intermediates

set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${CHAIN_FLT_PHIX_FA:?CHAIN_FLT_PHIX_FA must be set}"
: "${CHAIN_FLT_WORK_DIR:?CHAIN_FLT_WORK_DIR must be set}"

mkdir -p "$CHAIN_FLT_WORK_DIR"

ref="$CHAIN_FLT_WORK_DIR/phix.fa"
cp "$CHAIN_FLT_PHIX_FA" "$ref"

# --- Build SE FASTQ by slicing phix (deterministic, no PRNG). ---
seq="$(grep -v '^>' "$ref" | tr -d '\n' | tr 'acgt' 'ACGT')"
r1="$CHAIN_FLT_WORK_DIR/r1.fq"
: > "$r1"
L=150
n_reads=0
for off in 200 700 1200 1700 2200 2700 3200 3700 4200; do
    s="${seq:$off:$L}"
    [ "${#s}" -eq "$L" ] || continue
    q="$(printf 'I%.0s' $(seq 1 "$L"))"
    printf '@read%d\n%s\n+\n%s\n' "$n_reads" "$s" "$q" >> "$r1"
    n_reads=$((n_reads + 1))
done

if [ "$n_reads" -lt 5 ]; then
    echo "FAIL: fixture built only $n_reads reads (phix slicing regression?)" >&2
    exit 1
fi
echo "fixture: $n_reads SE reads sliced from phix"

"$BWA_MEM3" index "$ref" > "$CHAIN_FLT_WORK_DIR/index.log" 2>&1

# mapped_count <label> <-W value> [extra flags...]
mapped_count() {
    local label=$1 wt=$2
    shift 2
    local sam="$CHAIN_FLT_WORK_DIR/$label.sam"
    "$BWA_MEM3" mem "$@" -W "$wt" "$ref" "$r1" > "$sam" 2> /dev/null
    # FLAG bit 0x4 (unmapped) via arithmetic: `and()` is a gawk builtin and is
    # absent from BSD awk (macOS) and mawk, both of which run this suite.
    awk '!/^@/ && int($2 / 4) % 2 == 0' "$sam" | wc -l | tr -d ' '
}

# --- Guard: at -W 0 the weight filter cannot drop anything, so all three
# --- configurations must map every read. Without this, assertion 3 below
# --- would also pass on a fixture that simply fails to align.
for cfg in "default:" "compat-mem2:--compat=bwa-mem2" "compat-mem:--compat=bwa-mem"; do
    label=${cfg%%:*}
    flag=${cfg#*:}
    if [ -n "$flag" ]; then
        got=$(mapped_count "w0-$label" 0 "$flag")
    else got=$(mapped_count "w0-$label" 0); fi
    if [ "$got" -ne "$n_reads" ]; then
        echo "FAIL: guard -- at -W 0, $label mapped $got/$n_reads reads (expected all)" >&2
        exit 1
    fi
done
echo "guard: all $n_reads reads map at -W 0 under all three configurations"

# --- The load-bearing assertions, at -W 200 (> read length). ---
def_mapped=$(mapped_count "w200-default" 200)
mem2_mapped=$(mapped_count "w200-compat-mem2" 200 --compat=bwa-mem2)
mem_mapped=$(mapped_count "w200-compat-mem" 200 --compat=bwa-mem)

fail=0
if [ "$def_mapped" -ne "$n_reads" ]; then
    echo "FAIL: default at -W 200 mapped $def_mapped/$n_reads (expected all -- bwa-mem2 drop-in)" >&2
    fail=1
fi
if [ "$mem2_mapped" -ne "$n_reads" ]; then
    echo "FAIL: --compat=bwa-mem2 at -W 200 mapped $mem2_mapped/$n_reads (expected all)" >&2
    fail=1
fi
if [ "$mem_mapped" -ne 0 ]; then
    echo "FAIL: --compat=bwa-mem at -W 200 mapped $mem_mapped/$n_reads (expected 0 -- bwa leaves them unmapped)" >&2
    fail=1
fi
[ "$fail" -eq 0 ] || exit 1

echo "PASS: default=$def_mapped/$n_reads mapped, --compat=bwa-mem2=$mem2_mapped/$n_reads, --compat=bwa-mem=0/$n_reads (#310)"
