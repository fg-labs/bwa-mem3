#!/usr/bin/env bash
# test/regression/compat_byte_identical.sh
#
# Regression: `bwa-mem3 mem --compat` must suppress exactly the two
# bwa-mem3-only record additions that keep the drop-in profile from being
# a byte-for-byte match to bwa-mem2 — the `MQ:i` and `HN:i` tags — and
# must change NOTHING else.
#
# `@PG` is deliberately NOT suppressed: bwa-mem2 emits its own, so
# dropping ours would turn a changed line into a missing one, and `CL:`
# embeds the invocation either way. It is excluded from the comparison
# below (the default and --compat runs have different argv) and asserted
# to still be present.
#
# The load-bearing assertion is byte-identity: a default run with its
# MQ:i / HN:i tags stripped must be identical, byte-for-byte, to a
# `--compat` run of the same inputs. If --compat ever suppresses too much
# (drops another tag, drops @HD/@SQ) or too little (leaks MQ/HN), or
# perturbs an alignment, the streams diverge and this fails. Both the
# SAM-text and --bam paths are checked, since the tags are emitted from
# two independent code paths (bwamem.cpp and bam_writer.cpp) gated on the
# same MEM_F_COMPAT flag.
#
# Fixture (deterministic, no PRNG): PE reads are sliced directly from the
# committed phix.fa so the FASTQ is byte-identical across awk
# implementations. read1 is a forward window; read2 is the reverse
# complement of a downstream window, so every pair aligns concordantly
# and carries a mate (=> MQ:i is emitted) and a primary hit count
# (=> HN:i is emitted). Several pairs guarantee non-zero tag counts.
#
# Inputs (env vars):
#   BWA_MEM3         — path to the bwa-mem3 binary under test
#   COMPAT_PHIX_FA   — path to test/fixtures/phix.fa (the reference source)
#   COMPAT_WORK_DIR  — directory for fixture-private intermediates

set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${COMPAT_PHIX_FA:?COMPAT_PHIX_FA must be set}"
: "${COMPAT_WORK_DIR:?COMPAT_WORK_DIR must be set}"

mkdir -p "$COMPAT_WORK_DIR"

ref="$COMPAT_WORK_DIR/phix.fa"
cp "$COMPAT_PHIX_FA" "$ref"

# --- Build the PE FASTQs by slicing phix (deterministic, no PRNG). ---
# Flatten the reference to a single sequence, then emit N pairs at fixed
# offsets: read1 = ref[off .. off+L], read2 = revcomp(ref[off+G .. off+G+L]).
seq="$(grep -v '^>' "$ref" | tr -d '\n' | tr 'acgt' 'ACGT')"
r1="$COMPAT_WORK_DIR/r1.fq"
r2="$COMPAT_WORK_DIR/r2.fq"
: > "$r1"; : > "$r2"
L=120          # read length
G=300          # inner gap between the two windows (fragment ~ G+L)
n_pairs=0
for off in 200 900 1600 2300 3000 3700; do
    s1="${seq:$off:$L}"
    s2fwd="${seq:$((off+G)):$L}"
    # reverse complement of s2fwd
    s2="$(printf '%s' "$s2fwd" | rev | tr 'ACGT' 'TGCA')"
    [ "${#s1}" -eq "$L" ] && [ "${#s2}" -eq "$L" ] || continue
    q="$(printf 'I%.0s' $(seq 1 "$L"))"
    printf '@pair%d/1\n%s\n+\n%s\n' "$n_pairs" "$s1" "$q" >> "$r1"
    printf '@pair%d/2\n%s\n+\n%s\n' "$n_pairs" "$s2" "$q" >> "$r2"
    n_pairs=$((n_pairs+1))
done

if [ "$n_pairs" -lt 3 ]; then
    echo "FAIL: fixture built only $n_pairs pairs (phix slicing regression?)" >&2
    exit 1
fi
echo "fixture: $n_pairs PE pairs sliced from phix"

# --- Index the reference. ---
"$BWA_MEM3" index "$ref" > "$COMPAT_WORK_DIR/index.log" 2>&1

# --- SAM-text path: default vs --compat. ---
def_sam="$COMPAT_WORK_DIR/default.sam"
cmp_sam="$COMPAT_WORK_DIR/compat.sam"
"$BWA_MEM3" mem          "$ref" "$r1" "$r2" > "$def_sam" 2>/dev/null
"$BWA_MEM3" mem --compat "$ref" "$r1" "$r2" > "$cmp_sam" 2>/dev/null

# The default run must actually contain the tags we claim to suppress,
# else the byte-identity check below would pass vacuously.
def_mq=$(grep -c 'MQ:i:'  "$def_sam" || true)
def_hn=$(grep -c 'HN:i:'  "$def_sam" || true)
if [ "$def_mq" -lt 1 ] || [ "$def_hn" -lt 1 ]; then
    echo "FAIL: default SAM lacks MQ/HN to suppress (MQ=$def_mq HN=$def_hn) — fixture too weak" >&2
    exit 1
fi

# --compat must emit neither tag.
cmp_mq=$(grep -c 'MQ:i:'  "$cmp_sam" || true)
cmp_hn=$(grep -c 'HN:i:'  "$cmp_sam" || true)
if [ "$cmp_mq" -ne 0 ] || [ "$cmp_hn" -ne 0 ]; then
    echo "FAIL: --compat SAM still emits MQ=$cmp_mq HN=$cmp_hn (expected 0/0)" >&2
    exit 1
fi

# --compat must preserve the header, @PG included -- suppressing @PG would
# only trade a changed line for a missing one (bwa-mem2 emits its own).
if ! grep -q '^@HD' "$cmp_sam" || ! grep -q '^@SQ' "$cmp_sam"; then
    echo "FAIL: --compat SAM dropped @HD/@SQ header lines" >&2
    exit 1
fi
if ! grep -q '^@PG.*ID:bwa-mem3' "$cmp_sam"; then
    echo "FAIL: --compat SAM dropped the bwa-mem3 @PG line (it must be preserved)" >&2
    exit 1
fi

# The load-bearing check: default (minus MQ/HN) is byte-identical to --compat.
# @PG is excluded from BOTH sides -- its CL: records the actual argv, which
# necessarily differs between the default and --compat invocations.
def_stripped="$COMPAT_WORK_DIR/default.stripped.sam"
cmp_stripped="$COMPAT_WORK_DIR/compat.stripped.sam"
grep -v '^@PG' "$def_sam" \
  | sed -E 's/\tMQ:i:[0-9-]+//; s/\tHN:i:[0-9-]+//' > "$def_stripped"
grep -v '^@PG' "$cmp_sam" > "$cmp_stripped"
if ! diff -q "$def_stripped" "$cmp_stripped" > /dev/null; then
    echo "FAIL: --compat SAM differs from default beyond MQ/HN:" >&2
    diff "$def_stripped" "$cmp_stripped" | head -20 >&2
    exit 1
fi
echo "PASS: SAM-text --compat is byte-identical to default minus MQ/HN (@PG excluded)"

# --- BAM path: same assertion, if samtools is available. ---
if command -v samtools > /dev/null 2>&1; then
    def_bam="$COMPAT_WORK_DIR/default.bam"
    cmp_bam="$COMPAT_WORK_DIR/compat.bam"
    "$BWA_MEM3" mem --bam          "$ref" "$r1" "$r2" > "$def_bam" 2>/dev/null
    "$BWA_MEM3" mem --bam --compat "$ref" "$r1" "$r2" > "$cmp_bam" 2>/dev/null

    # --no-PG so samtools does not inject its own @PG line into the header dump.
    if [ "$(samtools view -H --no-PG "$cmp_bam" | grep -c '^@PG.*ID:bwa-mem3')" -ne 1 ]; then
        echo "FAIL: --compat BAM lacks the bwa-mem3 @PG line (it must be preserved)" >&2
        exit 1
    fi
    bam_def_rec="$COMPAT_WORK_DIR/default.bam.rec"
    bam_cmp_rec="$COMPAT_WORK_DIR/compat.bam.rec"
    samtools view --no-PG "$def_bam" \
      | sed -E 's/\tMQ:i:[0-9-]+//; s/\tHN:i:[0-9-]+//' > "$bam_def_rec"
    samtools view --no-PG "$cmp_bam" > "$bam_cmp_rec"
    if ! diff -q "$bam_def_rec" "$bam_cmp_rec" > /dev/null; then
        echo "FAIL: --compat BAM records differ from default beyond MQ/HN:" >&2
        diff "$bam_def_rec" "$bam_cmp_rec" | head -20 >&2
        exit 1
    fi
    echo "PASS: BAM --compat records are byte-identical to default minus MQ/HN"
else
    echo "SKIP: samtools not on PATH — BAM-path assertion skipped"
fi

# --- --compat and --fast are mutually exclusive: must be a hard error. ---
# (A diff-clean-looking stream over --fast's changed alignments would defeat
# the parity-validation purpose of --compat, so combining them is rejected.)
if "$BWA_MEM3" mem --compat --fast "$ref" "$r1" "$r2" > /dev/null 2>"$COMPAT_WORK_DIR/mutex.log"; then
    echo "FAIL: --compat --fast exited 0 (expected a hard error)" >&2
    exit 1
fi
if ! grep -q 'mutually exclusive' "$COMPAT_WORK_DIR/mutex.log"; then
    echo "FAIL: --compat --fast failed without the expected 'mutually exclusive' message:" >&2
    cat "$COMPAT_WORK_DIR/mutex.log" >&2
    exit 1
fi
echo "PASS: --compat --fast is rejected as a hard error"

# --- --compat is non-meth only: --compat --meth must be a hard error. ---
# bwa-mem2 has no methylation mode, so byte-identity is undefined under --meth.
# The guard fires during option validation, before any index load, so the
# non-meth phix reference here is sufficient to exercise it.
if "$BWA_MEM3" mem --compat --meth "$ref" "$r1" "$r2" > /dev/null 2>"$COMPAT_WORK_DIR/meth.log"; then
    echo "FAIL: --compat --meth exited 0 (expected a hard error)" >&2
    exit 1
fi
if ! grep -q 'not supported with --meth' "$COMPAT_WORK_DIR/meth.log"; then
    echo "FAIL: --compat --meth failed without the expected 'not supported with --meth' message:" >&2
    cat "$COMPAT_WORK_DIR/meth.log" >&2
    exit 1
fi
echo "PASS: --compat --meth is rejected as a hard error"

echo "PASS: compat byte-identical regression"
