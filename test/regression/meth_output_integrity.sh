#!/usr/bin/env bash
# test/regression/meth_output_integrity.sh
#
# Regression (D3 B6): the --meth BAM is original-alphabet and Bismark-compatible.
# Five invariants downstream methylation + variant callers depend on:
#
#  1. Real-SNP vs bisulfite-conversion. Under --meth, NM/MD are derived from the
#     per-hypothesis asymmetric matrix rather than from literal base inequality,
#     so a column is a mismatch iff the matrix penalizes it. An OT read with 10
#     C->T conversions plus one real A->G SNP must score AS = 60 - (a+b) (only the
#     SNP penalized, conversions free), NM = 1, and an MD string with exactly one
#     ref-A mismatch (the real variant) and NO ref-C entries -- the conversions
#     are matches for NM/MD exactly as they already are for the DP (issue #327).
#
#     (The complementary case -- that `collapsed` additionally hides a real T->C
#     variant while `genomic` keeps it -- is covered by meth_collapsed_scoring.sh,
#     which owns that fixture.)
#
#  2. Bismark four-strand XR/XG. XR is the per-read conversion (R1=CT, R2=GA); XG
#     is the genome strand shared by both mates of a fragment. A top-strand
#     fragment -> OT/CTOT (both XG:CT); a bottom-strand fragment -> OB/CTOB (both
#     XG:GA), tracking R1's mapping strand, exactly as bwameth/Bismark emit.
#
#  3. SEQ<->CIGAR orientation. A reverse-mapped mate stores SEQ as the reverse
#     complement of the input read with a consistent CIGAR length (no S1
#     inconsistent-BAM); samtools quickcheck passes.
#
#  4. MQ:i / HN:i tag-set parity with the non-meth writer. src/meth_bam.cpp is a
#     second BAM emitter alongside src/bam_writer.cpp, and a --meth BAM must not
#     be a subset of a non-meth one: a record whose mate is mapped carries MQ:i
#     equal to that mate's MAPQ, and a mapped record carries HN:i equal to its
#     XA hit count. Checked on a pair with asymmetric MAPQ and asymmetric hit
#     counts, so both tags are pinned to a value and not merely to presence.
#
#  5. SA:Z on a split alignment. The two records a chimeric read splits into must
#     cross-reference each other's placement, with all six SA fields intact --
#     the meth writer builds SA:Z itself and nothing else in the suite splits.
#
# Inputs:
#   BWA_MEM3 — path to the bwa-mem3 binary under test
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

# Deterministic 1500 bp reference (PRNG seed 4242).
REF=TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTAGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGAATCGTCCCTTCTAACGTAGTTCGCATAGTTCCCGTCCGTAGCCGGACTATTCGAACACCCAGTATTCGATTAACTCGGGCTTGACGTATTAGAGGCGTTAGTGTGCCAGGTAAGATACGCCAACGGAATTAACCTCTGTGACACTCCGCGGAGCCTTCGGACATATAAGTGATCGGGTCTACGTTTGTTAGACTTGAGACGTCTGTTAAGAGTTGGGTCTAATAAATCGCCTACACGTGGAGTCTAACGGGGAAGCGTCGAATCCTGATACATCATATAATGGAGCGTGTTATGAAAAAAGAGCATTCCATTGTACGAGCCGTGCCAGAAACGGCTTGACTACGTGAGCGTAGTGTTAGATAAACAGGAAACACTGACGCGGTTAGAAGGCGGATTGCCGGTAGGTTTTGGAAACATAAATACACACGGTATCATGTTGGGTCACGATTCCTATCACCGCACAGGGCCAACCATAGAAGAACTGAAAGAACTAATCTGGCGGCGGGCTCGGTGCTTATATTTTCCACCCAACATCGTGCACATTAGGCTCACCGCGCCCTACGGGCGAAGGGTGCGTACGGTGTTTATAAGGCGTGACGGCCCCAAGTAGAGGGTAATTCTGTGAAAGAATCTCAGGACGGTGGCATGAATTCAATTCCTTTTAAACCTATCGTTCCGACCTTATGCAATCCTTCAATGAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTACTCTTCACTGCCCAAGCGAGGCTAACCTTTTGAAATGTCACAGTCGAAGCATATCTCCCGTACATCTTTTTCGGAGATCGCAGCTCGCGGAGCTATAAGCGACTTAAGCCCTTGTGTCGGTGATCCCAAGGGTCTGACTCCTGTACCAGGGTTACTGTTTCGCTTTACGGAGTAGCCTGTGAGGTGAACTGAAAGGAGCATATTTGAGATCTAAGATAGGGTCCTCCTCTGCGTCTACGTTCTCTCCGTTACGTACGGCTTCGCACCGGAGTGCATCTTGGCCCCGAAACGCACTGTGTGTGCTGATACAGCGTCCCTGGCCGGCCATGGGTTCAGAACTCCCGGGAACGCTTTTCAACTTAGAGGAACCCCGTCATGGAAGTAGATCGCGTCGAATGAGGGAGTTAGTCCTCGTTCCAGCTGGTAATTGTTTTACCGCTTGGGACCACTATAGGCCGCGGGTAGAAGTTGCTGGGTGTTGATTCCAACCCTCGAACCACGATACGACCTGCCATTTATGGCACAGTAAGGTTCAAACAGCATAATGAATACAGTTATAGTAACTTCCTCACGTACGATTAGGACGCAGCCTTG

printf '>chrA\n%s\n' "$REF" > ref.fa
Q=$(printf 'I%.0s' $(seq 1 60))
emit() { printf '@%s\n%s\n+\n%s\n' "$1" "$2" "$Q" > "$3"; }
"$BWA_MEM3" index --meth ref.fa > /dev/null 2>&1 || fail "index --meth nonzero exit"

tag() { mawk -v k="$2" '{for(i=12;i<=NF;i++) if(substr($i,1,5)==k":Z:"||substr($i,1,5)==k":i:") {print substr($i,6); exit}}' <<< "$1"; }

# --- 1. real SNP vs conversion ---------------------------------------------
# OT fwd @101: 10 C->T conversions (free) + 1 real A->G SNP at read pos 22 (penalized).
SNP_READ=ATTCTGGCGATTTAGGTGACTCGCAGAATCGTTCTTTCTAACGTAGTTTGTATAGTTCTC
emit snp "$SNP_READ" snp.fq
"$BWA_MEM3" mem --meth --meth-scoring genomic -t 1 ref.fa snp.fq > snp.bam 2> /dev/null || fail "snp mem --meth nonzero exit"
samtools quickcheck snp.bam || fail "snp invalid BAM"
line=$(samtools view snp.bam | mawk 'NR==1')
as=$(tag "$line" AS)
nm=$(tag "$line" NM)
md=$(tag "$line" MD)
# --meth keeps the bwa default mismatch penalty b=4 (not bwameth's b=2), so the single real SNP costs -4
# (conversions free): 59 match/freed columns (+59) - 1 SNP (-4) = 55.
[ "$as" = "55" ] || fail "SNP/conv: AS $as, want 55 (only the real SNP penalized at b=4; conversions free)"
[ "$nm" = "1" ] || fail "SNP/conv: NM $nm, want 1 (the real SNP only; conversions are matches for NM/MD)"
# MD reference-base mismatch letters: exactly one ref-A (the SNP), no ref-C (the
# conversions must not appear at all).
nC=$(printf '%s' "$md" | tr -cd 'C' | wc -c | tr -d ' ')
nA=$(printf '%s' "$md" | tr -cd 'A' | wc -c | tr -d ' ')
nG=$(printf '%s' "$md" | tr -cd 'G' | wc -c | tr -d ' ')
[ "$nC" = "0" ] || fail "SNP/conv: MD has $nC ref-C mismatches, want 0 (conversions must be hidden); MD=$md"
[ "$nA" = "1" ] || fail "SNP/conv: MD has $nA ref-A mismatches, want 1 (the real SNP); MD=$md"
[ "$nG" = "0" ] || fail "SNP/conv: MD has $nG ref-G mismatches, want 0; MD=$md"

# --- 2. Bismark four-strand XR/XG ------------------------------------------
# Forward (top-strand) fragment: R1 OT fwd, R2 CTOT rev -> both XG:CT.
FWD_R1=ATTCTGGCGATTTAGGTGACTCACAGAATCGTTCTTTCTAACGTAGTTTGTATAGTTCTC
FWD_R2=TAAGCGATTTATTAAACCCAACTCTTAACAAACGTCTCAAATCTAACAAACGTAAACCCG
# Reverse (bottom-strand) fragment: R1 OB rev, R2 CTOB fwd -> both XG:GA.
REV_R1=AGATTCTTTCACAGAATTACTCTCTATTTGGGGCTGTCACGCTTTATAAATATCGTACGC
REV_R2=CTAACGCGATTAAAAGACAGATTGCCAGTAAGTTTTAGAAACATAAATACACACAGTATC

xrxg() { # $1 bam  $2 mateflag -> "XR XG"
    samtools view "$1" | mawk -v bit="$2" '(int($2/bit)%2)==1{xr="";xg="";for(i=12;i<=NF;i++){if($i~/^XR:Z:/)xr=substr($i,6);if($i~/^XG:Z:/)xg=substr($i,6)}print xr,xg;exit}'
}
check_strand() { # $1 bam $2 mateflag $3 label $4 wantXR $5 wantXG
    read -r xr xg < <(xrxg "$1" "$2")
    [ "$xr" = "$4" ] || fail "$3: XR $xr, want $4"
    [ "$xg" = "$5" ] || fail "$3: XG $xg, want $5 (Bismark genome strand)"
}
emit f "$FWD_R1" f1.fq
emit f "$FWD_R2" f2.fq
"$BWA_MEM3" mem --meth --meth-scoring genomic -t 1 ref.fa f1.fq f2.fq > fwd.bam 2> /dev/null || fail "fwd nonzero exit"
check_strand fwd.bam 64 "OT   R1" CT CT
check_strand fwd.bam 128 "CTOT R2" GA CT
emit r "$REV_R1" r1.fq
emit r "$REV_R2" r2.fq
"$BWA_MEM3" mem --meth --meth-scoring genomic -t 1 ref.fa r1.fq r2.fq > rev.bam 2> /dev/null || fail "rev nonzero exit"
check_strand rev.bam 64 "OB   R1" CT GA
check_strand rev.bam 128 "CTOB R2" GA GA

# --- 3. SEQ<->CIGAR orientation on a reverse-mapped mate -------------------
SEQ_R1=ATCCCGGCGACTTAGGCGACCCACAGAATCGTCCCTTCTAACGTAGTTCGCATAGTTCCC
SEQ_R2=TAGGCGATTTATTAGACCCAACTCTTAACAGACGTCTCAAGTCTAACAAACGTAGACCCG
SEQ_R2_RC=CGGGTCTACGTTTGTTAGACTTGAGACGTCTGTTAAGAGTTGGGTCTAATAAATCGCCTA
emit pp "$SEQ_R1" s1.fq
emit pp "$SEQ_R2" s2.fq
"$BWA_MEM3" mem --meth --meth-scoring genomic -t 1 ref.fa s1.fq s2.fq > seq.bam 2> /dev/null || fail "seq nonzero exit"
samtools quickcheck seq.bam || fail "seq invalid BAM (S1 inconsistent-BAM trap)"
read -r rev seq cig < <(samtools view seq.bam | mawk '(int($2/128)%2)==1{print (int($2/16)%2),$10,$6;exit}')
[ "$rev" = "1" ] || fail "SEQ orient: R2 expected reverse-mapped (rev=$rev)"
[ "$cig" = "60M" ] || fail "SEQ orient: CIGAR $cig, want 60M"
[ "${#seq}" = "60" ] || fail "SEQ orient: SEQ length ${#seq}, want 60 (CIGAR consistency)"
[ "$seq" = "$SEQ_R2_RC" ] || fail "SEQ orient: reverse-mate SEQ is not revcomp(input read)"

# --- 4. MQ:i / HN:i parity with the non-meth writer ------------------------
# The fixture is deliberately ASYMMETRIC on both tags: chrB repeats chrA:241-420,
# which contains R2's locus, so R2 has two equally good hits (MAPQ 0, one XA hit)
# while R1 stays unique (MAPQ 60, no XA hits). That pins both tags to a VALUE
# rather than to presence — a writer emitting the record's own MAPQ as MQ:i would
# survive a symmetric MAPQ-60 pair, and an HN:i wired to a constant 0 would
# survive a fixture with no repeats. Neither survives this one.
#
# The MAPQ preconditions are asserted first so that a future scoring change which
# flattens the asymmetry is reported as a stale fixture rather than passing
# vacuously.
printf '>chrA\n%s\n>chrB\n%s\n' "$REF" "${REF:240:180}" > dup.fa
"$BWA_MEM3" index --meth dup.fa > /dev/null 2>&1 || fail "dup index --meth nonzero exit"
"$BWA_MEM3" mem --meth --meth-scoring genomic -t 1 dup.fa f1.fq f2.fq > dup.bam 2> /dev/null \
    || fail "dup mem --meth nonzero exit"
samtools quickcheck dup.bam || fail "dup invalid BAM"

mq_hn() { # $1 bam  $2 mate flag bit -> "MAPQ MQ HN"
    samtools view "$1" | mawk -v bit="$2" '
        (int($2/bit)%2)==1 {
            mq = ""; hn = ""
            for (i = 12; i <= NF; i++) {
                if ($i ~ /^MQ:i:/) mq = substr($i, 6)
                if ($i ~ /^HN:i:/) hn = substr($i, 6)
            }
            print $5, mq, hn
            exit
        }'
}
read -r mapq1 mq1 hn1 < <(mq_hn dup.bam 64)
read -r mapq2 mq2 hn2 < <(mq_hn dup.bam 128)
[ "$mapq1" = "60" ] || fail "MQ/HN: stale fixture — R1 MAPQ $mapq1, want 60 (unique locus)"
[ "$mapq2" = "0" ] || fail "MQ/HN: stale fixture — R2 MAPQ $mapq2, want 0 (locus duplicated on chrB)"
# Non-empty first, so an absent tag is reported as absent rather than as a value
# mismatch against the empty string.
[ -n "$mq1" ] || fail "MQ/HN: R1 has no MQ:i (the non-meth BAM writer emits it)"
[ -n "$mq2" ] || fail "MQ/HN: R2 has no MQ:i (the non-meth BAM writer emits it)"
[ "$mq1" = "$mapq2" ] || fail "MQ/HN: R1 MQ:i $mq1, want R2's MAPQ $mapq2 (not R1's own $mapq1)"
[ "$mq2" = "$mapq1" ] || fail "MQ/HN: R2 MQ:i $mq2, want R1's MAPQ $mapq1 (not R2's own $mapq2)"
# HN counts the hits XA enumerates: 0 for the unique mate, 1 for the mate whose
# locus is duplicated on chrB.
[ -n "$hn1" ] || fail "MQ/HN: R1 has no HN:i (the non-meth BAM writer emits it)"
[ -n "$hn2" ] || fail "MQ/HN: R2 has no HN:i (the non-meth BAM writer emits it)"
[ "$hn1" = "0" ] || fail "MQ/HN: R1 HN:i $hn1, want 0 (unique mapper, no XA hits)"
[ "$hn2" = "1" ] || fail "MQ/HN: R2 HN:i $hn2, want 1 (the chrA copy of the duplicated locus)"

# --- 5. SA:Z on a split (chimeric) alignment -------------------------------
# The meth writer builds SA:Z itself, and nothing else in the --meth suite
# produces a split alignment, so the whole builder is otherwise untested. A read
# whose halves come from chrA:101 and chrA:1001 splits into two records that must
# cross-reference each other's position. Asserting the cross-reference and the
# field count — rather than a literal CIGAR — keeps this robust to scoring
# changes that move the split point, while still catching a truncated SA:Z.
CHIM="${REF:100:75}${REF:1000:75}"
printf '@chim\n%s\n+\n%s\n' "$CHIM" "$(printf 'I%.0s' $(seq 1 150))" > chim.fq
"$BWA_MEM3" mem --meth --meth-scoring genomic -t 1 ref.fa chim.fq > chim.bam 2> /dev/null \
    || fail "chim mem --meth nonzero exit"
samtools quickcheck chim.bam || fail "chim invalid BAM"
n_chim=$(samtools view -c chim.bam)
[ "$n_chim" = "2" ] || fail "SA: stale fixture — chimeric read gave $n_chim records, want 2 (a split alignment)"

sa_of() { # $1 record number -> "RNAME POS SA:Z-value"
    samtools view chim.bam | mawk -v n="$1" 'NR==n {
        sa = ""
        for (i = 12; i <= NF; i++) if ($i ~ /^SA:Z:/) sa = substr($i, 6)
        print $3, $4, sa
    }'
}
read -r rname1 pos1 sa1 < <(sa_of 1)
read -r rname2 pos2 sa2 < <(sa_of 2)
[ -n "$sa1" ] || fail "SA: split record 1 has no SA:Z"
[ -n "$sa2" ] || fail "SA: split record 2 has no SA:Z"
# SA:Z is `rname,pos,strand,cigar,mapq,NM;` — a truncated build loses the
# trailing ';' or a field, which is exactly the silent-truncation failure mode.
for sa in "$sa1" "$sa2"; do
    case "$sa" in
        *\;) ;;
        *) fail "SA: '$sa' does not end in ';' (truncated?)" ;;
    esac
    nf=$(printf '%s' "${sa%;}" | mawk -F, '{print NF}')
    [ "$nf" = "6" ] || fail "SA: '$sa' has $nf comma-separated fields, want 6 (truncated?)"
done
# Each half points at the other's placement.
[ "$(printf '%s' "$sa1" | cut -d, -f1)" = "$rname2" ] \
    || fail "SA: record 1's SA:Z names $(printf '%s' "$sa1" | cut -d, -f1), want record 2's RNAME $rname2"
[ "$(printf '%s' "$sa1" | cut -d, -f2)" = "$pos2" ] \
    || fail "SA: record 1's SA:Z points at $(printf '%s' "$sa1" | cut -d, -f2), want record 2's POS $pos2"
[ "$(printf '%s' "$sa2" | cut -d, -f1)" = "$rname1" ] \
    || fail "SA: record 2's SA:Z names $(printf '%s' "$sa2" | cut -d, -f1), want record 1's RNAME $rname1"
[ "$(printf '%s' "$sa2" | cut -d, -f2)" = "$pos1" ] \
    || fail "SA: record 2's SA:Z points at $(printf '%s' "$sa2" | cut -d, -f2), want record 1's POS $pos1"

echo "PASS: meth_output_integrity (real-SNP vs conversion in AS/NM/MD; Bismark four-strand XR/XG; reverse SEQ orientation; MQ:i/HN:i tag parity; SA:Z cross-reference on a split alignment)"
