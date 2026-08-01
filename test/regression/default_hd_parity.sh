#!/usr/bin/env bash
# test/regression/default_hd_parity.sh
#
# Regression: the DEFAULT `@HD` record must be byte-identical on every output
# path — SAM text, `--bam`, and `--meth`.
#
# It was not. Each of the three writers hardcoded its own literal, and the two
# BAM ones had drifted, so the same run emitted a different `@HD` depending on
# `--bam` (fg-labs/bwa-mem3#288):
#
#   src/bwa.cpp        @HD VN:1.5 SO:unsorted GO:query
#   src/bam_writer.cpp @HD VN:1.6 SO:unsorted
#   src/meth_bam.cpp   @HD VN:1.6 SO:unsorted
#
# All three now emit BWAMEM3_DEFAULT_HD_LINE (src/compat_target.h), which is
# byte-identical to upstream bwa's literal at bwa.c:426. This test pins that:
# a fourth writer, or a "fix" to any one of the three, has to keep them equal.
#
# The assertion is equality ACROSS PATHS *and* against the expected string.
# Equality alone would pass if all three drifted together; the literal alone
# would not catch one path diverging on a value the test did not think to
# check. Both are cheap, so assert both.
#
# `--compat` is out of scope here: it suppresses the default `@HD` entirely,
# which compat_byte_identical.sh covers. This is about DEFAULT output.
#
# Inputs (env vars):
#   BWA_MEM3        — path to the bwa-mem3 binary under test
#   HD_PARITY_PHIX_FA — path to test/fixtures/phix.fa
#   HD_PARITY_WORK_DIR — directory for fixture-private intermediates

set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${HD_PARITY_PHIX_FA:?HD_PARITY_PHIX_FA must be set}"
: "${HD_PARITY_WORK_DIR:?HD_PARITY_WORK_DIR must be set}"

# The one true default, duplicated here ON PURPOSE: a test that read the
# constant from the source it is testing could not catch the constant changing.
EXPECT=$'@HD\tVN:1.5\tSO:unsorted\tGO:query'

mkdir -p "$HD_PARITY_WORK_DIR"
cd "$HD_PARITY_WORK_DIR"
# Remove only this script's own artifacts (never `rm -rf` a caller-supplied path).
rm -f phix.fa phix.fa.* phix.dict phix.meth.* r1.fq r2.fq ./*.sam ./*.bam ./*.hdr ./*.log

cp "$HD_PARITY_PHIX_FA" phix.fa
ref=phix.fa

# The reference must have NO sidecar: a <prefix>.hdr / <baseprefix>.dict @HD
# takes precedence over the default and would make every assertion vacuous.
rm -f "$ref.hdr" phix.dict
"$BWA_MEM3" index "$ref" > index.log 2>&1

# PE reads sliced from phix — deterministic, no PRNG.
python3 - << 'PY'
seq = "".join(l.strip() for l in open("phix.fa") if not l.startswith(">")).upper()
comp = str.maketrans("ACGT", "TGCA")
with open("r1.fq", "w") as a, open("r2.fq", "w") as b:
    for i, off in enumerate([200, 900, 1600]):
        x = seq[off:off+120]
        y = seq[off+300:off+420][::-1].translate(comp)
        a.write("@p%d/1\n%s\n+\n%s\n" % (i, x, "I"*120))
        b.write("@p%d/2\n%s\n+\n%s\n" % (i, y, "I"*120))
PY

check() { # <label> <actual-@HD>
    if [ -z "$2" ]; then
        echo "FAIL: $1 emitted no default @HD" >&2
        exit 1
    fi
    if [ "$2" != "$EXPECT" ]; then
        echo "FAIL: $1 default @HD differs from the expected default" >&2
        printf '  expected: %s\n  actual:   %s\n' \
            "$(printf '%s' "$EXPECT" | sed 's/\t/\\t/g')" \
            "$(printf '%s' "$2" | sed 's/\t/\\t/g')" >&2
        exit 1
    fi
    printf '  %-12s %s\n' "$1" "$(printf '%s' "$2" | sed 's/\t/\\t/g')"
}

echo "default @HD by output path:"

# --- SAM text ---
"$BWA_MEM3" mem "$ref" r1.fq r2.fq > sam.sam 2> sam.err
sam_hd=$(grep -m1 '^@HD' sam.sam || true)
check "SAM text" "$sam_hd"

if ! command -v samtools > /dev/null 2>&1; then
    echo "SKIP: samtools not on PATH — --bam and --meth paths not checked"
    echo "PASS: default @HD parity (SAM text only)"
    exit 0
fi

# --- --bam ---
"$BWA_MEM3" mem --bam "$ref" r1.fq r2.fq > bam.bam 2> bam.err
bam_hd=$(samtools view -H --no-PG bam.bam | grep -m1 '^@HD' || true)
check "--bam" "$bam_hd"

# --- --meth (needs the D3 seed index) ---
if [[ ! -s "$ref.meth.bwt.2bit.64" || ! -s "$ref.meth.amb" ||
    ! -s "$ref.meth.ann" || ! -s "$ref.meth.pac" ]]; then
    "$BWA_MEM3" index --meth "$ref" > index-meth.log 2>&1 \
        || {
            echo "FAIL: bwa-mem3 index --meth failed" >&2
            cat index-meth.log >&2
            exit 1
        }
fi
"$BWA_MEM3" mem --meth "$ref" r1.fq > meth.bam 2> meth.err \
    || {
        echo "FAIL: bwa-mem3 mem --meth exited non-zero" >&2
        cat meth.err >&2
        exit 1
    }
meth_hd=$(samtools view -H --no-PG meth.bam | grep -m1 '^@HD' || true)
check "--meth" "$meth_hd"

# Equality across paths is implied by all three matching EXPECT, but assert it
# directly so a future change that relaxes EXPECT still cannot let them drift.
if [ "$sam_hd" != "$bam_hd" ] || [ "$sam_hd" != "$meth_hd" ]; then
    echo "FAIL: default @HD differs across output paths" >&2
    exit 1
fi

echo "PASS: default @HD is byte-identical across SAM text, --bam and --meth"
