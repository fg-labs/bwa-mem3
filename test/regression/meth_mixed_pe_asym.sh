#!/usr/bin/env bash
# test/regression/meth_mixed_pe_asym.sh
#
# Regression (D3 A1/PR-4 + PR-6): a MIXED-hypothesis paired-end batch — R1 on the
# OT strand (C->T conversions) interleaved with R2 on the OB strand (G->A
# conversions) — is scored with the correct per-hypothesis asymmetric matrix on
# the live `mem --meth` path. Before A1 a mixed batch fell back to the SYMMETRIC
# matrix for the batched-SIMD placement scorer, so bisulfite conversions were
# penalized as mismatches (deflated AS/MAPQ). After A1's per-hypothesis partition
# pass each mate is scored against its OT/OB object, so conversions are FREE in
# scoring (AS == perfect) while remaining literal mismatches in NM/MD (for
# downstream Revelio-style masking).
#
# The decisive assertion per mate: AS == read length (60) AND NM == #conversions
# (> 0). A symmetric-scored conversion would cost (a + b) per event, dropping AS
# well below 60.
#
# Inputs:
#   BWA_MEM3 — path to the bwa-mem3 binary under test
set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"

command -v samtools >/dev/null 2>&1 || { echo "SKIP: samtools not on PATH (--meth emits BAM)"; exit 0; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

fail() { echo "FAIL: $*" >&2; exit 1; }

# Deterministic 1500 bp reference (PRNG seed 4242; see test generation notes).
REF=TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTAGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGAATCGTCCCTTCTAACGTAGTTCGCATAGTTCCCGTCCGTAGCCGGACTATTCGAACACCCAGTATTCGATTAACTCGGGCTTGACGTATTAGAGGCGTTAGTGTGCCAGGTAAGATACGCCAACGGAATTAACCTCTGTGACACTCCGCGGAGCCTTCGGACATATAAGTGATCGGGTCTACGTTTGTTAGACTTGAGACGTCTGTTAAGAGTTGGGTCTAATAAATCGCCTACACGTGGAGTCTAACGGGGAAGCGTCGAATCCTGATACATCATATAATGGAGCGTGTTATGAAAAAAGAGCATTCCATTGTACGAGCCGTGCCAGAAACGGCTTGACTACGTGAGCGTAGTGTTAGATAAACAGGAAACACTGACGCGGTTAGAAGGCGGATTGCCGGTAGGTTTTGGAAACATAAATACACACGGTATCATGTTGGGTCACGATTCCTATCACCGCACAGGGCCAACCATAGAAGAACTGAAAGAACTAATCTGGCGGCGGGCTCGGTGCTTATATTTTCCACCCAACATCGTGCACATTAGGCTCACCGCGCCCTACGGGCGAAGGGTGCGTACGGTGTTTATAAGGCGTGACGGCCCCAAGTAGAGGGTAATTCTGTGAAAGAATCTCAGGACGGTGGCATGAATTCAATTCCTTTTAAACCTATCGTTCCGACCTTATGCAATCCTTCAATGAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTACTCTTCACTGCCCAAGCGAGGCTAACCTTTTGAAATGTCACAGTCGAAGCATATCTCCCGTACATCTTTTTCGGAGATCGCAGCTCGCGGAGCTATAAGCGACTTAAGCCCTTGTGTCGGTGATCCCAAGGGTCTGACTCCTGTACCAGGGTTACTGTTTCGCTTTACGGAGTAGCCTGTGAGGTGAACTGAAAGGAGCATATTTGAGATCTAAGATAGGGTCCTCCTCTGCGTCTACGTTCTCTCCGTTACGTACGGCTTCGCACCGGAGTGCATCTTGGCCCCGAAACGCACTGTGTGTGCTGATACAGCGTCCCTGGCCGGCCATGGGTTCAGAACTCCCGGGAACGCTTTTCAACTTAGAGGAACCCCGTCATGGAAGTAGATCGCGTCGAATGAGGGAGTTAGTCCTCGTTCCAGCTGGTAATTGTTTTACCGCTTGGGACCACTATAGGCCGCGGGTAGAAGTTGCTGGGTGTTGATTCCAACCCTCGAACCACGATACGACCTGCCATTTATGGCACAGTAAGGTTCAAACAGCATAATGAATACAGTTATAGTAACTTCCTCACGTACGATTAGGACGCAGCCTTG

# R1: OT (C->T) read over chrA[100:160] (1-based pos 101), forward.  NM == 10 conversions.
R1=ATTCTGGCGATTTAGGTGACTCACAGAATCGTTCTTTCTAACGTAGTTTGTATAGTTCTC
# R2: OB (G->A) read, reverse mate over chrA[300:360] (1-based pos 301).  NM == 5 conversions.
R2=TAAGCGATTTATTAAACCCAACTCTTAACAAACGTCTCAAATCTAACAAACGTAAACCCG

printf '>chrA\n%s\n' "$REF" > ref.fa
Q=$(printf 'I%.0s' $(seq 1 60))
printf '@p\n%s\n+\n%s\n' "$R1" "$Q" > r1.fq
printf '@p\n%s\n+\n%s\n' "$R2" "$Q" > r2.fq

"$BWA_MEM3" index --meth ref.fa >/dev/null 2>&1 || fail "index --meth nonzero exit"
"$BWA_MEM3" mem --meth --meth-scoring genomic -t 1 ref.fa r1.fq r2.fq > pe.bam 2>/dev/null || fail "mem --meth nonzero exit"
samtools quickcheck pe.bam || fail "mem --meth produced an invalid BAM"

# Extract one decoded line per mate (READ1 = flag&64, READ2 = flag&128).
get() { # $1 = mate-flag-bit (64=READ1, 128=READ2)  -> echo "RNAME POS MAPQ rev CIGAR AS NM XR"
    # mawk has no bitwise and(); extract a flag bit via integer arithmetic.
    samtools view pe.bam | mawk -v bit="$1" '
        (int($2/bit) % 2) == 1 {
            rev = (int($2/16) % 2) == 1 ? 1 : 0; as=""; nm=""; xr="";
            for (i=12;i<=NF;i++){ if($i~/^AS:i:/){as=substr($i,6)} if($i~/^NM:i:/){nm=substr($i,6)} if($i~/^XR:Z:/){xr=substr($i,6)} }
            print $3, $4, $5, rev, $6, as, nm, xr; exit
        }'
}

check_mate() { # $1=label $2=mateflag $3=want_pos $4=want_rev $5=want_as $6=want_nm $7=want_xr
    read -r rn pos mapq rev cig as nm xr < <(get "$2")
    [ "$rn" = "chrA" ]      || fail "$1: RNAME $rn, want chrA"
    [ "$pos" = "$3" ]       || fail "$1: POS $pos, want $3"
    [ "$rev" = "$4" ]       || fail "$1: strand rev=$rev, want $4"
    [ "$cig" = "60M" ]      || fail "$1: CIGAR $cig, want 60M"
    [ "$mapq" = "60" ]      || fail "$1: MAPQ $mapq, want 60 (asym scoring must not deflate MAPQ)"
    [ "$xr" = "$7" ]        || fail "$1: XR $xr, want $7"
    # The A1 contract: conversions are FREE in scoring (AS == full length) yet
    # appear as literal mismatches in NM (> 0). A symmetric fallback would give AS < $5.
    [ "$as" = "$5" ]        || fail "$1: AS $as, want $5 (conversions must be scored free under the per-hypothesis matrix)"
    [ "$nm" = "$6" ]        || fail "$1: NM $nm, want $6 (conversions must remain literal mismatches for downstream masking)"
    [ "$nm" -gt 0 ]         || fail "$1: NM must be > 0 (test must actually exercise conversions)"
}

# R1 OT forward @101, AS 60 (10 C->T free), NM 10, XR CT.
check_mate "R1/OT" 64  101 0 60 10 CT
# R2 OB reverse @301, AS 60 (5 G->A free), NM 5, XR GA.
check_mate "R2/OB" 128 301 1 60 5  GA

echo "PASS: meth_mixed_pe_asym (mixed OT+OB PE batch scored per-hypothesis; conversions free in AS, literal in NM)"
