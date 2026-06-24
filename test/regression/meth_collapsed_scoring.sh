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
# scores (a+b) higher under collapsed than under genomic, while the literal
# mismatch is preserved in NM either way. This exercises bandedSWA's general
# (>=2 freed cell) matrix path that collapsed mode uses.

set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"
command -v samtools >/dev/null 2>&1 || { echo "SKIP: samtools not on PATH (--meth emits BAM)"; exit 0; }

WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT; cd "$WORK"
fail() { echo "FAIL: $*" >&2; exit 1; }

# Deterministic 1500 bp reference (same generator as meth_output_integrity).
REF=TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTAGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGAATCGTCCCTTCTAACGTAGTTCGCATAGTTCCCGTCCGTAGCCGGACTATTCGAACACCCAGTATTCGATTAACTCGGGCTTGACGTATTAGAGGCGTTAGTGTGCCAGGTAAGATACGCCAACGGAATTAACCTCTGTGACACTCCGCGGAGCCTTCGGACATATAAGTGATCGGGTCTACGTTTGTTAGACTTGAGACGTCTGTTAAGAGTTGGGTCTAATAAATCGCCTACACGTGGAGTCTAACGGGGAAGCGTCGAATCCTGATACATCATATAATGGAGCGTGTTATGAAAAAAGAGCATTCCATTGTACGAGCCGTGCCAGAAACGGCTTGACTACGTGAGCGTAGTGTTAGATAAACAGGAAACACTGACGCGGTTAGAAGGCGGATTGCCGGTAGGTTTTGGAAACATAAATACACACGGTATCATGTTGGGTCACGATTCCTATCACCGCACAGGGCCAACCATAGAAGAACTGAAAGAACTAATCTGGCGGCGGGCTCGGTGCTTATATTTTCCACCCAACATCGTGCACATTAGGCTCACCGCGCCCTACGGGCGAAGGGTGCGTACGGTGTTTATAAGGCGTGACGGCCCCAAGTAGAGGGTAATTCTGTGAAAGAATCTCAGGACGGTGGCATGAATTCAATTCCTTTTAAACCTATCGTTCCGACCTTATGCAATCCTTCAATGAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTACTCTTCACTGCCCAAGCGAGGCTAACCTTTTGAAATGTCACAGTCGAAGCATATCTCCCGTACATCTTTTTCGGAGATCGCAGCTCGCGGAGCTATAAGCGACTTAAGCCCTTGTGTCGGTGATCCCAAGGGTCTGACTCCTGTACCAGGGTTACTGTTTCGCTTTACGGAGTAGCCTGTGAGGTGAACTGAAAGGAGCATATTTGAGATCTAAGATAGGGTCCTCCTCTGCGTCTACGTTCTCTCCGTTACGTACGGCTTCGCACCGGAGTGCATCTTGGCCCCGAAACGCACTGTGTGTGCTGATACAGCGTCCCTGGCCGGCCATGGGTTCAGAACTCCCGGGAACGCTTTTCAACTTAGAGGAACCCCGTCATGGAAGTAGATCGCGTCGAATGAGGGAGTTAGTCCTCGTTCCAGCTGGTAATTGTTTTACCGCTTGGGACCACTATAGGCCGCGGGTAGAAGTTGCTGGGTGTTGATTCCAACCCTCGAACCACGATACGACCTGCCATTTATGGCACAGTAAGGTTCAAACAGCATAATGAATACAGTTATAGTAACTTCCTCACGTACGATTAGGACGCAGCCTTG
printf '>chrA\n%s\n' "$REF" > ref.fa
"$BWA_MEM3" index --meth ref.fa >/dev/null 2>&1 || fail "index --meth nonzero exit"

# 60 bp forward window; flip one internal ref-T to C (no other change, so no
# C->T conversions — the single ref-T x read-C diff is the only off-diagonal).
START=300; LEN=60
SUB=${REF:START:LEN}
READ=""; FLIPPED=0
for ((i=0; i<LEN; i++)); do
    b=${SUB:i:1}
    if [ "$FLIPPED" -eq 0 ] && [ "$b" = "T" ] && [ "$i" -ge 18 ] && [ "$i" -le 42 ]; then
        READ+="C"; FLIPPED=1
    else
        READ+="$b"
    fi
done
[ "$FLIPPED" -eq 1 ] || fail "no internal ref-T found to flip in window"
Q=$(printf 'I%.0s' $(seq 1 $LEN))
printf '@m\n%s\n+\n%s\n' "$READ" "$Q" > r.fq

tag() { mawk -v k="$2" '{for(i=12;i<=NF;i++) if(substr($i,1,5)==k":i:"){print substr($i,6); exit}}' <<<"$1"; }
as_of() { # $1 = scoring mode -> echoes "AS NM" of the primary alignment
    "$BWA_MEM3" mem --meth --meth-scoring "$1" -t 1 ref.fa r.fq > o.bam 2>/dev/null || fail "$1 mem --meth nonzero exit"
    samtools quickcheck o.bam || fail "$1 invalid BAM"
    local line; line=$(samtools view -F 0x104 o.bam | mawk 'NR==1')
    [ -n "$line" ] || fail "$1: no primary alignment"
    echo "$(tag "$line" AS) $(tag "$line" NM)"
}

read -r AS_G NM_G < <(as_of genomic)
read -r AS_C NM_C < <(as_of collapsed)

echo "[meth_collapsed] genomic: AS=$AS_G NM=$NM_G   collapsed: AS=$AS_C NM=$NM_C"

# a=1, b=4 (collapsed default sets b=2, but we want a fixed delta: pin -B 4 in
# BOTH so the only variable is the matrix). Re-run with explicit -B 4.
as_of_b4() {
    "$BWA_MEM3" mem --meth --meth-scoring "$1" -B 4 -t 1 ref.fa r.fq > o.bam 2>/dev/null || fail "$1 -B4 nonzero exit"
    local line; line=$(samtools view -F 0x104 o.bam | mawk 'NR==1')
    echo "$(tag "$line" AS) $(tag "$line" NM)"
}
read -r AS_Gb NM_Gb <<<"$(as_of_b4 genomic)"
read -r AS_Cb NM_Cb <<<"$(as_of_b4 collapsed)"

# The mirror cell is a mismatch under genomic, a match under collapsed:
#   collapsed AS - genomic AS == a + b == 1 + 4 == 5.
[ $((AS_Cb - AS_Gb)) -eq 5 ] || fail "collapsed should free the ref-T x read-C mirror: AS_collapsed=$AS_Cb AS_genomic=$AS_Gb (want diff 5)"
# Literal mismatch preserved in NM regardless of scoring mode.
[ "$NM_Gb" = "$NM_Cb" ] || fail "NM should be identical across modes (literal mismatch preserved): genomic=$NM_Gb collapsed=$NM_Cb"
[ "$NM_Cb" = "1" ] || fail "expected exactly one literal mismatch (the flipped base): NM=$NM_Cb"

echo "OK"
