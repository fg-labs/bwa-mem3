#!/usr/bin/env bash
# test/regression/meth_collapsed_scoring.sh
#
# Regression: --meth-scoring collapsed reproduces bwameth's collapsed-space
# behavior end-to-end by ALSO freeing the conversion MIRROR cell, where genomic
# keeps it a real mismatch.
#
# For an OT (forward) read, the conversion cell ref-C x read-T is freed in BOTH
# modes; the MIRROR cell ref-T x read-C is a real mismatch under `genomic`
# (variant-aware) but freed under `collapsed` (C/T interchangeable). So a read
# that differs from the reference only by a single ref-T -> read-C substitution
# scores (a+b) higher under collapsed than under genomic. NM is derived from the
# same matrix, so the variant is a mismatch under genomic (NM=1) and hidden under
# collapsed (NM=0) -- the observable cost of bwameth-compatible placement. This
# exercises bandedSWA's general (>=2 freed cell) matrix path that collapsed uses.
#
# The mirror cell alone cannot separate `neutral` from `genomic` -- both leave it
# at -b -- so a second read exercises the CONVERSION cell itself, which genomic
# scores +a and neutral scores 0. There NM must be 0 in both (neither penalizes
# it, and NM is "mismatch iff the matrix penalizes"), while AS differs by exactly
# one match score. Without that case a regression collapsing neutral into genomic
# would pass this fixture.

set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"
command -v samtools > /dev/null 2>&1 || {
    echo "SKIP: samtools not on PATH (--meth emits BAM)"
    exit 0
}

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
fail() {
    echo "FAIL: $*" >&2
    exit 1
}

# Deterministic 1500 bp reference (same generator as meth_output_integrity).
REF=TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTAGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGAATCGTCCCTTCTAACGTAGTTCGCATAGTTCCCGTCCGTAGCCGGACTATTCGAACACCCAGTATTCGATTAACTCGGGCTTGACGTATTAGAGGCGTTAGTGTGCCAGGTAAGATACGCCAACGGAATTAACCTCTGTGACACTCCGCGGAGCCTTCGGACATATAAGTGATCGGGTCTACGTTTGTTAGACTTGAGACGTCTGTTAAGAGTTGGGTCTAATAAATCGCCTACACGTGGAGTCTAACGGGGAAGCGTCGAATCCTGATACATCATATAATGGAGCGTGTTATGAAAAAAGAGCATTCCATTGTACGAGCCGTGCCAGAAACGGCTTGACTACGTGAGCGTAGTGTTAGATAAACAGGAAACACTGACGCGGTTAGAAGGCGGATTGCCGGTAGGTTTTGGAAACATAAATACACACGGTATCATGTTGGGTCACGATTCCTATCACCGCACAGGGCCAACCATAGAAGAACTGAAAGAACTAATCTGGCGGCGGGCTCGGTGCTTATATTTTCCACCCAACATCGTGCACATTAGGCTCACCGCGCCCTACGGGCGAAGGGTGCGTACGGTGTTTATAAGGCGTGACGGCCCCAAGTAGAGGGTAATTCTGTGAAAGAATCTCAGGACGGTGGCATGAATTCAATTCCTTTTAAACCTATCGTTCCGACCTTATGCAATCCTTCAATGAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTACTCTTCACTGCCCAAGCGAGGCTAACCTTTTGAAATGTCACAGTCGAAGCATATCTCCCGTACATCTTTTTCGGAGATCGCAGCTCGCGGAGCTATAAGCGACTTAAGCCCTTGTGTCGGTGATCCCAAGGGTCTGACTCCTGTACCAGGGTTACTGTTTCGCTTTACGGAGTAGCCTGTGAGGTGAACTGAAAGGAGCATATTTGAGATCTAAGATAGGGTCCTCCTCTGCGTCTACGTTCTCTCCGTTACGTACGGCTTCGCACCGGAGTGCATCTTGGCCCCGAAACGCACTGTGTGTGCTGATACAGCGTCCCTGGCCGGCCATGGGTTCAGAACTCCCGGGAACGCTTTTCAACTTAGAGGAACCCCGTCATGGAAGTAGATCGCGTCGAATGAGGGAGTTAGTCCTCGTTCCAGCTGGTAATTGTTTTACCGCTTGGGACCACTATAGGCCGCGGGTAGAAGTTGCTGGGTGTTGATTCCAACCCTCGAACCACGATACGACCTGCCATTTATGGCACAGTAAGGTTCAAACAGCATAATGAATACAGTTATAGTAACTTCCTCACGTACGATTAGGACGCAGCCTTG
printf '>chrA\n%s\n' "$REF" > ref.fa
"$BWA_MEM3" index --meth ref.fa > /dev/null 2>&1 || fail "index --meth nonzero exit"

# 60 bp forward window; flip one internal ref-T to C (no other change, so no
# C->T conversions — the single ref-T x read-C diff is the only off-diagonal).
START=300
LEN=60
SUB=${REF:START:LEN}
READ=""
FLIPPED=0
for ((i = 0; i < LEN; i++)); do
    b=${SUB:i:1}
    if [ "$FLIPPED" -eq 0 ] && [ "$b" = "T" ] && [ "$i" -ge 18 ] && [ "$i" -le 42 ]; then
        READ+="C"
        FLIPPED=1
    else
        READ+="$b"
    fi
done
[ "$FLIPPED" -eq 1 ] || fail "no internal ref-T found to flip in window"
Q=$(printf 'I%.0s' $(seq 1 $LEN))
printf '@m\n%s\n+\n%s\n' "$READ" "$Q" > r.fq

# Same window, but flip one internal ref-C to T: a pure OT conversion, the only
# off-diagonal column, and the cell whose SCORE differs between genomic (+a) and
# neutral (0).
CONV=""
CONVERTED=0
for ((i = 0; i < LEN; i++)); do
    b=${SUB:i:1}
    if [ "$CONVERTED" -eq 0 ] && [ "$b" = "C" ] && [ "$i" -ge 18 ] && [ "$i" -le 42 ]; then
        CONV+="T"
        CONVERTED=1
    else
        CONV+="$b"
    fi
done
[ "$CONVERTED" -eq 1 ] || fail "no internal ref-C found to convert in window"
printf '@c\n%s\n+\n%s\n' "$CONV" "$Q" > c.fq

tag() { mawk -v k="$2" '{for(i=12;i<=NF;i++) if(substr($i,1,5)==k":i:"){print substr($i,6); exit}}' <<< "$1"; }
as_of() { # $1 = scoring mode -> echoes "AS NM" of the primary alignment
    "$BWA_MEM3" mem --meth --meth-scoring "$1" -t 1 ref.fa r.fq > o.bam 2> /dev/null || fail "$1 mem --meth nonzero exit"
    samtools quickcheck o.bam || fail "$1 invalid BAM"
    local line
    line=$(samtools view -F 0x104 o.bam | mawk 'NR==1')
    [ -n "$line" ] || fail "$1: no primary alignment"
    echo "$(tag "$line" AS) $(tag "$line" NM)"
}

read -r AS_G NM_G < <(as_of genomic)
read -r AS_C NM_C < <(as_of collapsed)

echo "[meth_collapsed] genomic: AS=$AS_G NM=$NM_G   collapsed: AS=$AS_C NM=$NM_C"

# a=1, b=4 (collapsed default sets b=2, but we want a fixed delta: pin -B 4 in
# BOTH so the only variable is the matrix). Re-run with explicit -B 4.
as_of_b4() { # $1 = scoring mode, $2 = fastq -> echoes "AS NM" of the primary alignment
    "$BWA_MEM3" mem --meth --meth-scoring "$1" -B 4 -t 1 ref.fa "$2" > o.bam 2> /dev/null || fail "$1 -B4 nonzero exit on $2"
    local line
    line=$(samtools view -F 0x104 o.bam | mawk 'NR==1')
    [ -n "$line" ] || fail "$1 -B4: no primary alignment for $2"
    echo "$(tag "$line" AS) $(tag "$line" NM)"
}
read -r AS_Gb NM_Gb <<< "$(as_of_b4 genomic r.fq)"
read -r AS_Cb NM_Cb <<< "$(as_of_b4 collapsed r.fq)"

# The mirror cell is a mismatch under genomic, a match under collapsed:
#   collapsed AS - genomic AS == a + b == 1 + 4 == 5.
[ $((AS_Cb - AS_Gb)) -eq 5 ] || fail "collapsed should free the ref-T x read-C mirror: AS_collapsed=$AS_Cb AS_genomic=$AS_Gb (want diff 5)"
# NM follows the same matrix, so the mirror cell is visible in NM too: a real
# ref-T x read-C variant is a mismatch under genomic and hidden under collapsed.
# This is the observable cost of bwameth-compatible placement (issue #327).
[ "$NM_Gb" = "1" ] || fail "genomic should keep the ref-T x read-C variant as a mismatch: NM=$NM_Gb, want 1"
[ "$NM_Cb" = "0" ] || fail "collapsed should hide the ref-T x read-C variant (C/T interchangeable): NM=$NM_Cb, want 0"

# `neutral` (the --meth=taps default) is the third mode and the only one whose
# conversion cell is 0 rather than +a, so it exercises a distinct edge of the
# "mismatch iff the matrix penalizes it" predicate. Like `genomic` it leaves the
# mirror cell at -b, so the real ref-T x read-C variant must stay visible.
read -r AS_Nb NM_Nb <<< "$(as_of_b4 neutral r.fq)"
[ "$NM_Nb" = "1" ] || fail "neutral should keep the ref-T x read-C variant as a mismatch: NM=$NM_Nb, want 1"
# This read has no conversion column, so the one cell where neutral and genomic
# differ is never exercised: their scores must agree exactly. (The conversion
# read below is what separates them.)
[ "$AS_Nb" = "$AS_Gb" ] || fail "neutral and genomic must score the mirror cell identically: AS_neutral=$AS_Nb AS_genomic=$AS_Gb"

# ---------------------------------------------------------------------------
# Conversion cell (ref-C x read-T), the one cell genomic and neutral score
# DIFFERENTLY: +a under genomic, 0 under neutral. Neither value is negative, so
# "mismatch iff the matrix penalizes it" makes it a match for NM in both -- the
# NM=0 the whole issue is about. The scores must still diverge by exactly one
# match score, which is what pins neutral apart from genomic.
read -r AS_Gc NM_Gc <<< "$(as_of_b4 genomic c.fq)"
read -r AS_Nc NM_Nc <<< "$(as_of_b4 neutral c.fq)"

echo "[meth_collapsed] conversion read -> genomic: AS=$AS_Gc NM=$NM_Gc   neutral: AS=$AS_Nc NM=$NM_Nc"

[ "$NM_Gc" = "0" ] || fail "genomic should not count the ref-C x read-T conversion: NM=$NM_Gc, want 0"
[ "$NM_Nc" = "0" ] || fail "neutral should not count the ref-C x read-T conversion (cell scores 0, not < 0): NM=$NM_Nc, want 0"
# a=1: genomic credits the conversion +1, neutral credits it 0, and the read has
# no other off-diagonal column, so the whole AS difference is that one cell.
[ $((AS_Gc - AS_Nc)) -eq 1 ] || fail "neutral must score the conversion 0 where genomic scores +a=1: AS_genomic=$AS_Gc AS_neutral=$AS_Nc (want diff 1)"

# ---------------------------------------------------------------------------
# Degenerate -B under --meth. "Mismatch iff the matrix penalizes it" only says
# something while SOME cell is penalized: bwa_fill_scmat stores -b, so -B 0
# makes every substitution free and -B <0 makes it a positive reward. Either
# way NM collapses to 0 and MD to a bare match run for REAL variants, not just
# conversions -- output that looks clean because scoring is degenerate. Both
# must be refused, so the bound is `<= 0` rather than `== 0`.
#
# Assert the guard's own message, not merely a non-zero exit: `mem` exits 1 for
# plenty of unrelated reasons (missing index, bad path), so an exit-code-only
# check would pass with the guard never firing.
reject_meth() { # $1 = human label, $2... = args before the positional operands
    local label="$1"
    shift
    local rc=0
    "$BWA_MEM3" mem --meth "$@" -t 1 ref.fa r.fq > /dev/null 2> badB.err || rc=$?
    [ "$rc" -eq 1 ] || fail "$label must be rejected with exit 1, got exit $rc"
    grep -q 'requires a positive mismatch penalty' badB.err \
        || fail "$label must be rejected by the mismatch-penalty guard, but stderr was: $(cat badB.err)"
}
# Explicit -B: bwa_fill_scmat stores -b, so 0 frees every substitution and a
# negative -B rewards it.
for badB in 0 -1; do
    reject_meth "--meth -B $badB" -B "$badB"
done
# ...and the -A route, which is why the message reports the EFFECTIVE penalty
# rather than echoing "-B": update_a() scales b by -A when -B was not given, so
# -A 0 drives b to 0 without the user ever passing -B.
reject_meth "--meth -A 0 (b scaled to 0 by update_a)" -A 0

# ...and the rejection is scoped to --meth: the non-meth path keeps the literal
# base comparison, where -B 0 only affects placement, never NM/MD truthfulness.
# A separate non-doubled index is required -- `index --meth` builds the doubled
# per-strand reference, which plain `mem` must not be pointed at.
cp ref.fa ref_nonmeth.fa
"$BWA_MEM3" index ref_nonmeth.fa > /dev/null 2>&1 || fail "index (non-meth) nonzero exit"
"$BWA_MEM3" mem -B 0 -t 1 ref_nonmeth.fa r.fq > nm.sam 2> /dev/null \
    || fail "non-meth -B 0 must still be accepted (the --meth guard must not leak into the default path)"
NM_LIT=$(mawk '!/^@/{for(i=12;i<=NF;i++) if(substr($i,1,5)=="NM:i:"){print substr($i,6); exit}}' nm.sam)
# The read carries one real ref-T x read-C substitution; literal NM must see it
# even at -B 0, which is exactly why -B 0 is harmless off the --meth path.
[ "$NM_LIT" = "1" ] || fail "non-meth -B 0 should keep literal NM: NM=$NM_LIT, want 1"

echo "PASS: meth_collapsed_scoring (collapsed-space --meth scoring, and literal NM preserved off the --meth path at -B 0)"
