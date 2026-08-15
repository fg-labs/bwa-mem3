#!/usr/bin/env bash
# test/regression/hic_alias_identical.sh
#
# Regression: `bwa-mem3 mem --hic` must be exactly `-5SP` — same flag bits,
# same output — and must keep being so as the three letters evolve.
#
# `--hic` exists because `-5SP` reads as three unrelated single letters, and
# the one that dominates cost on Hi-C (`-S`) is the one a reader is least
# likely to infer: mate rescue runs BEFORE the pairing bail-out in
# mem_sam_pe, so `-P` alone skips pairing while leaving the full rescue
# Smith-Waterman in place. An alias only helps if it cannot drift from what
# it aliases, which is what this script pins.
#
# The load-bearing assertion is byte-identity between a `--hic` run and a
# `-5SP` run of the same inputs. If `--hic` ever sets a different bit set
# (a letter added to one and not the other, a flag renumbered), the streams
# diverge and this fails. `@PG` is excluded from both sides: its `CL:`
# records the argv, which necessarily differs between the two invocations.
#
# Two guards keep the identity check from passing vacuously:
#   1. `--hic` output must actually DIFFER from a plain `mem` run. If the
#      flags stopped taking effect entirely, `--hic` == `-5SP` would still
#      hold trivially while the option did nothing.
#   2. `--hic` must clear the 0x2 proper-pair bit that a plain run sets on
#      this fixture — the observable signature of -P (no pairing), and the
#      cheapest positive proof the flags reached the aligner.
#
# Fixture (deterministic, no PRNG): PE reads sliced directly from the
# committed phix.fa, same construction as compat_byte_identical.sh, so the
# FASTQ is byte-identical across awk implementations and every pair aligns
# concordantly (=> a plain run flags them proper, giving guard 2 something
# to observe).
#
# Inputs (env vars):
#   BWA_MEM3        — path to the bwa-mem3 binary under test
#   HIC_PHIX_FA     — path to test/fixtures/phix.fa (the reference source)
#   HIC_WORK_DIR    — directory for fixture-private intermediates

set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${HIC_PHIX_FA:?HIC_PHIX_FA must be set}"
: "${HIC_WORK_DIR:?HIC_WORK_DIR must be set}"

mkdir -p "$HIC_WORK_DIR"

ref="$HIC_WORK_DIR/phix.fa"
cp "$HIC_PHIX_FA" "$ref"

# --- Build the PE FASTQs by slicing phix (deterministic, no PRNG). ---
seq="$(grep -v '^>' "$ref" | tr -d '\n' | tr 'acgt' 'ACGT')"
r1="$HIC_WORK_DIR/r1.fq"
r2="$HIC_WORK_DIR/r2.fq"
: > "$r1"
: > "$r2"
L=120 # read length
G=300 # inner gap between the two windows (fragment ~ G+L)
# Enough pairs to clear mem_pestat's MIN_DIR_CNT (10 unique pairs per
# orientation) with margin. Below that threshold no orientation gets an
# insert-size distribution, so a plain run flags nothing proper and -5SP
# becomes a no-op -- both guards below would then fire on a perfectly
# correct --hic. phix is ~5.4 kb, so a 200 bp stride fits ~24 pairs.
n_pairs=0
for off in $(seq 100 200 4700); do
    s1="${seq:$off:$L}"
    s2fwd="${seq:$((off + G)):$L}"
    s2="$(printf '%s' "$s2fwd" | rev | tr 'ACGT' 'TGCA')"
    [ "${#s1}" -eq "$L" ] && [ "${#s2}" -eq "$L" ] || continue
    q="$(printf 'I%.0s' $(seq 1 "$L"))"
    printf '@pair%d/1\n%s\n+\n%s\n' "$n_pairs" "$s1" "$q" >> "$r1"
    printf '@pair%d/2\n%s\n+\n%s\n' "$n_pairs" "$s2" "$q" >> "$r2"
    n_pairs=$((n_pairs + 1))
done

if [ "$n_pairs" -lt 3 ]; then
    echo "FAIL: fixture built only $n_pairs pairs (phix slicing regression?)" >&2
    exit 1
fi
echo "fixture: $n_pairs PE pairs sliced from phix"

# --- Index the reference. ---
"$BWA_MEM3" index "$ref" > "$HIC_WORK_DIR/index.log" 2>&1

# --- Three runs: plain, -5SP, --hic. ---
plain_sam="$HIC_WORK_DIR/plain.sam"
long_sam="$HIC_WORK_DIR/5SP.sam"
hic_sam="$HIC_WORK_DIR/hic.sam"
"$BWA_MEM3" mem "$ref" "$r1" "$r2" > "$plain_sam" 2> /dev/null
"$BWA_MEM3" mem -5SP "$ref" "$r1" "$r2" > "$long_sam" 2> /dev/null
"$BWA_MEM3" mem --hic "$ref" "$r1" "$r2" > "$hic_sam" 2> /dev/null

# Guard 1: --hic must actually change something versus a plain run, else the
# identity check below would hold trivially for a no-op option.
if diff -q <(grep -v '^@PG' "$plain_sam") <(grep -v '^@PG' "$hic_sam") > /dev/null 2>&1; then
    echo "FAIL: --hic output is identical to a plain run — the flags are not taking effect" >&2
    exit 1
fi

# Guard 2: a plain run flags these concordant pairs proper (0x2); --hic must
# not, since -P skips pairing. This is the observable signature of the alias
# having reached the aligner, not just the parser.
# int($2/2)%2 rather than and($2,2): `and()` is a gawk builtin and is absent
# from BSD awk (macOS) and mawk, both of which run this suite.
plain_proper=$(awk '!/^@/ && int($2 / 2) % 2 == 1' "$plain_sam" | wc -l | tr -d ' ')
hic_proper=$(awk '!/^@/ && int($2 / 2) % 2 == 1' "$hic_sam" | wc -l | tr -d ' ')
if [ "$plain_proper" -lt 1 ]; then
    echo "FAIL: plain run flagged no proper pairs — fixture too weak for guard 2" >&2
    exit 1
fi
if [ "$hic_proper" -ne 0 ]; then
    echo "FAIL: --hic emitted $hic_proper proper-pair (0x2) records; -P should leave none" >&2
    exit 1
fi
echo "guards: plain=$plain_proper proper-pair records, --hic=0"

# The load-bearing check: --hic is byte-identical to -5SP. @PG is excluded
# from BOTH sides because its CL: records the argv, which differs by
# construction between the two invocations.
if ! diff -u <(grep -v '^@PG' "$long_sam") <(grep -v '^@PG' "$hic_sam") > "$HIC_WORK_DIR/hic.diff"; then
    echo "FAIL: --hic is not byte-identical to -5SP:" >&2
    head -40 "$HIC_WORK_DIR/hic.diff" >&2
    exit 1
fi

echo "PASS: --hic is byte-identical to -5SP ($n_pairs pairs, phix)"
