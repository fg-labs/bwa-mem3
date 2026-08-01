#!/usr/bin/env bash
# test/regression/header_parity.sh
#
# Regression for the two header behaviors that compat_byte_identical.sh
# cannot reach with its phix fixture, because both need a reference that
# has ALT contigs and a .hdr/.dict sidecar:
#
#   1. AH:* on ALT contigs is emitted by the GENERATED @SQ block on BOTH
#      output paths. Upstream bwa (bwa.c:432) and bwa-mem2 (bwa.cpp:538)
#      both emit it; bwa-mem3's SAM text path did, but the BAM writer was
#      ported as an SN+LN-only loop and dropped it for EVERY ALT-aware
#      reference, sidecar or not (fg-labs/bwa-mem3#281).
#
#   2. `--compat=bwa-mem2` ignores the sidecar entirely, so @SQ regenerates
#      as bare SN/LN (+AH:*). The <prefix>.hdr / <baseprefix>.dict sidecar
#      is bwa-mem3-only -- a port of lh3/bwa#348, which lh3 closed unmerged
#      -- so neither upstream has anything to load, and its M5/AS/UR/SP
#      would break byte-identity.
#
# The sidecar's own @SQ is AUTHORITATIVE and is never rewritten: lh3/bwa#348
# was designed for the block to be produced complete by an external dict
# tool, and `samtools dict --alt` (samtools/samtools#1676) exists to do it.
# So a sidecar WITHOUT AH legitimately yields output without AH -- asserted
# below, along with the warning that says so, and with the `samtools dict
# --alt` path that fixes it.
#
# Fixture is generated here, not committed (repo convention): a 3-contig FASTA
# whose ALT contig name EXTENDS a primary one (chr1 / chr1_alt), a matching
# .alt, and PE reads sliced from the primary contig. Sequence comes from
# python3's Mersenne Twister at a fixed seed, so the fixture is identical on
# every run and platform.
#
# Inputs (env vars):
#   BWA_MEM3         — path to the bwa-mem3 binary under test
#   HEADER_PARITY_WORK_DIR — directory for fixture-private intermediates.
#                     This script DELETES the fixture files it owns inside
#                     it, so give it a directory of its own.

set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${HEADER_PARITY_WORK_DIR:?HEADER_PARITY_WORK_DIR must be set}"

mkdir -p "$HEADER_PARITY_WORK_DIR"
cd "$HEADER_PARITY_WORK_DIR"
# Remove only what a previous run of THIS script created. A blanket
# `rm -rf "$HEADER_PARITY_WORK_DIR"` would obey a mistyped or over-broad value
# of a caller-supplied variable; the `:?` guard above catches unset and empty,
# but not `/` or `$HOME`.
rm -f ref.fasta ref.fasta.* ref.dict index.log probe.sam \
    ./*.sq ./*.hdr ./*.out ./*.err r1.fq r2.fq

ALT_CONTIG=chr1_alt

# --- Generate the ALT-aware fixture (deterministic; fixed PRNG seed). ---
python3 - << 'PY'
import random
random.seed(282)                      # fixed: fixture must be reproducible
def rnd(n): return "".join(random.choice("ACGT") for _ in range(n))
# chr1_alt deliberately EXTENDS chr1. The sidecar then carries both SN:chr1
# and SN:chr1_alt, so the AH scan in bwa.cpp must match SN as a whole token to
# report the right contig -- a prefix match would conflate them. Cases 2 and 3
# below exercise exactly that through the missing-AH warning.
contigs = [("chr1", 20000), ("chr2", 15000), ("chr1_alt", 5000)]
seqs = {n: rnd(l) for n, l in contigs}
with open("ref.fasta", "w") as f:
    for n, _ in contigs:
        f.write(">%s\n" % n)
        s = seqs[n]
        for i in range(0, len(s), 60):
            f.write(s[i:i+60] + "\n")
# <prefix>.alt is a SAM file; bwa reads only column 1 (QNAME) to set is_alt.
with open("ref.fasta.alt", "w") as f:
    f.write("chr1_alt\t0\tchr1\t1\t60\t5000M\t*\t0\t0\t*\t*\n")
comp = str.maketrans("ACGT", "TGCA")
with open("r1.fq", "w") as a, open("r2.fq", "w") as b:
    for i, off in enumerate(range(500, 4500, 400)):
        x = seqs["chr1"][off:off+120]
        y = seqs["chr1"][off+300:off+420][::-1].translate(comp)
        a.write("@p%d/1\n%s\n+\n%s\n" % (i, x, "I"*120))
        b.write("@p%d/2\n%s\n+\n%s\n" % (i, y, "I"*120))
PY

[ -s ref.fasta ] && [ -s ref.fasta.alt ] && [ -s r1.fq ] && [ -s r2.fq ] \
    || {
        echo "FAIL: fixture generation produced empty files" >&2
        exit 1
    }

"$BWA_MEM3" index ref.fasta > index.log 2>&1

# The index must actually mark the ALT contig, else every AH assertion below
# passes vacuously. `mem` prints no is_alt summary, so infer it from the one
# path that has always emitted AH: SAM text with no sidecar.
rm -f ref.dict ref.fasta.hdr
"$BWA_MEM3" mem ref.fasta r1.fq r2.fq 2> /dev/null > probe.sam
grep -q "^@SQ.*SN:${ALT_CONTIG}.*AH:\*" probe.sam \
    || {
        echo "FAIL: .alt not honored — no AH:* on $ALT_CONTIG even on the SAM/no-sidecar path" >&2
        grep '^@SQ' probe.sam >&2
        exit 1
    }
echo "fixture: 3 contigs, $ALT_CONTIG marked ALT via ref.fasta.alt"

have_samtools=0
command -v samtools > /dev/null 2>&1 && have_samtools=1
[ "$have_samtools" -eq 1 ] || echo "NOTE: samtools not on PATH — BAM and sidecar cases will be skipped"

# Run one invocation, capturing its stderr as "$tag.err", its raw output as
# "$tag.out", the full header as "$tag.hdr" and the @SQ block as "$tag.sq".
# Extracting the header once here means the assertions below are plain greps
# over a file and none of them needs to know which output path ran.
# ('@' cannot start a SAM record -- it is excluded from QNAME -- so grep '^@'
# over SAM text picks up header lines only.)
sq_of() { # <tag> [bwa-mem3 args...]
    local tag="$1"
    shift
    local a is_bam=0
    for a in "$@"; do [ "$a" = "--bam" ] && is_bam=1; done
    "$BWA_MEM3" mem "$@" ref.fasta r1.fq r2.fq 2> "$tag.err" > "$tag.out"
    if [ "$is_bam" -eq 1 ]; then
        samtools view -H --no-PG "$tag.out" > "$tag.hdr"
    else
        grep '^@' "$tag.out" > "$tag.hdr"
    fi
    grep '^@SQ' "$tag.hdr" > "$tag.sq" \
        || {
            echo "FAIL: $tag emitted no @SQ header block" >&2
            exit 1
        }
}

# Assert the run emitted no @HD at all.
assert_no_HD() { # <tag> <human description>
    if grep -q '^@HD' "$1.hdr"; then
        echo "FAIL: $2 — @HD emitted under compat (bwa-mem2 emits none)" >&2
        exit 1
    fi
}

assert_alt_has_AH() { # <tag> <human description>
    grep -q "SN:${ALT_CONTIG}.*AH:\*" "$1.sq" \
        || {
            echo "FAIL: $2 — AH:* missing on $ALT_CONTIG:" >&2
            cat "$1.sq" >&2
            exit 1
        }
}
assert_no_primary_AH() { # <tag> <human description>
    # Spec: AH "must not be present on sequences in the primary assembly".
    # Match SN as a whole tab-delimited token. Not load-bearing for the current
    # fixture (chr1_alt is the longer name, so a substring filter happens to
    # behave), but it keeps this assertion correct if the ALT contig is ever
    # renamed to a PREFIX of another contig, where a loose filter would drop
    # both lines and the check would silently stop firing.
    if grep -v "SN:${ALT_CONTIG}$(printf '\t')" "$1.sq" | grep -q 'AH:'; then
        echo "FAIL: $2 — AH on a primary-assembly contig:" >&2
        cat "$1.sq" >&2
        exit 1
    fi
}

# --- 1. No sidecar: AH:* on BOTH paths, default AND compat. ---------------
# The BAM/no-sidecar cell is fg-labs/bwa-mem3#281's broader half: it dropped
# AH for every ALT-aware reference, independent of --compat.
rm -f ref.dict ref.fasta.hdr
sq_of nosc_sam_def
assert_alt_has_AH nosc_sam_def "SAM text, no sidecar, default"
assert_no_primary_AH nosc_sam_def "SAM text, no sidecar, default"
sq_of nosc_sam_cmp --compat=bwa-mem2
assert_alt_has_AH nosc_sam_cmp "SAM text, no sidecar, --compat=bwa-mem2"
if [ "$have_samtools" -eq 1 ]; then
    sq_of nosc_bam_def --bam
    assert_alt_has_AH nosc_bam_def "BAM, no sidecar, default (issue #281)"
    assert_no_primary_AH nosc_bam_def "BAM, no sidecar, default"
    sq_of nosc_bam_cmp --bam --compat=bwa-mem2
    assert_alt_has_AH nosc_bam_cmp "BAM, no sidecar, --compat=bwa-mem2"
fi
echo "PASS: AH:* emitted on generated @SQ (SAM+BAM, default and --compat)"

# The remaining cases need a sidecar, which we build with `samtools dict`.
if [ "$have_samtools" -eq 0 ]; then
    echo "SKIP: sidecar cases need samtools dict"
    echo "PASS: compat header parity regression (partial — no samtools)"
    exit 0
fi

# --- 2. Sidecar WITHOUT AH (what `samtools dict` emits by default). -------
samtools dict -a testasm -s testus -o ref.dict ref.fasta
if grep -q 'AH:' ref.dict; then
    echo "FAIL: fixture sidecar unexpectedly carries AH — case 2 is vacuous" >&2
    exit 1
fi

# Default: the sidecar is authoritative, so its tags appear and AH does not.
# This is DESIGNED behavior (lh3/bwa#348), not a bug -- do not "fix" it.
assert_sidecar_verbatim() { # <tag>
    grep -q 'M5:' "$1.sq" \
        || {
            echo "FAIL: $1 — sidecar identity tags missing; sidecar not honored" >&2
            exit 1
        }
    if grep -q 'AH:' "$1.sq"; then
        echo "FAIL: $1 — AH present, but the sidecar is authoritative and has none" >&2
        exit 1
    fi
}
sq_of sc_sam_def
assert_sidecar_verbatim sc_sam_def
sq_of sc_bam_def --bam
assert_sidecar_verbatim sc_bam_def
# ...and the run must SAY so, rather than losing ALT status silently.
grep -q 'sidecar supplies @SQ without an AH tag' sc_sam_def.err \
    || {
        echo "FAIL: no warning about the sidecar's missing AH:" >&2
        cat sc_sam_def.err >&2
        exit 1
    }
grep -q 'samtools dict --alt' sc_sam_def.err \
    || {
        echo "FAIL: missing-AH warning does not name the remedy" >&2
        cat sc_sam_def.err >&2
        exit 1
    }
echo "PASS: sidecar @SQ is authoritative (AH absent), and the gap is warned about"

# compat: sidecar skipped entirely -> bare SN/LN + AH:*, no identity tags.
assert_compat_ignores_sidecar() { # <tag>
    assert_alt_has_AH "$1" "$1"
    if grep -qE "\t(M5|AS|UR|SP):" "$1.sq"; then
        echo "FAIL: $1 — sidecar identity tags leaked into compat output:" >&2
        head -3 "$1.sq" >&2
        exit 1
    fi
    # The sidecar's own @HD must not sneak through either.
    assert_no_HD "$1" "$1"
    # No warning: the sidecar was never read.
    if grep -q 'sidecar supplies @SQ without an AH tag' "$1.err"; then
        echo "FAIL: $1 — warned about a sidecar that compat does not read" >&2
        exit 1
    fi
}
sq_of sc_sam_cmp --compat=bwa-mem2
assert_compat_ignores_sidecar sc_sam_cmp
sq_of sc_bam_cmp --bam --compat=bwa-mem2
assert_compat_ignores_sidecar sc_bam_cmp
echo "PASS: --compat=bwa-mem2 ignores the sidecar (bare SN/LN + AH:*, no @HD)"

# --- 3. Sidecar WITH AH (`samtools dict --alt`, the supported remedy). ----
samtools dict -a testasm -s testus --alt ref.fasta.alt -o ref.dict ref.fasta
grep -q "SN:${ALT_CONTIG}.*AH:" ref.dict \
    || {
        echo "FAIL: 'samtools dict --alt' produced no AH — fixture/samtools mismatch" >&2
        exit 1
    }
assert_sidecar_AH_passthrough() { # <tag>
    assert_alt_has_AH "$1" "$1 (sidecar carries AH)"
    assert_no_primary_AH "$1" "$1"
    grep -q 'M5:' "$1.sq" \
        || {
            echo "FAIL: $1 — sidecar identity tags missing" >&2
            exit 1
        }
    # No warning: the sidecar supplies AH, so there is no gap to report.
    if grep -q 'sidecar supplies @SQ without an AH tag' "$1.err"; then
        echo "FAIL: $1 — warned despite the sidecar carrying AH" >&2
        exit 1
    fi
}
sq_of ah_sam
assert_sidecar_AH_passthrough ah_sam
sq_of ah_bam --bam
assert_sidecar_AH_passthrough ah_bam
echo "PASS: 'samtools dict --alt' sidecar carries AH:* through both paths, no warning"

echo "PASS: compat header parity regression"
