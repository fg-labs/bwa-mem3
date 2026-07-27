#!/usr/bin/env bash
# test/regression/meth_output_integrity.sh
#
# Regression (D3 B6): the --meth BAM is original-alphabet and Bismark-compatible.
# Three invariants downstream methylation + variant callers depend on:
#
#  1. Real-SNP vs bisulfite-conversion. The asymmetric matrix frees conversions
#     in SCORING but NM/MD count every literal mismatch. An OT read with 10 C->T
#     conversions plus one real A->G SNP must score AS = 60 - (a+b) (only the SNP
#     penalized, conversions free), NM = 11, and an MD string with ten ref-C
#     mismatches (the conversions, the Revelio masking substrate) and exactly one
#     ref-A mismatch (the real variant) -- so a downstream caller sees the SNP and
#     a Revelio double-mask hides the conversions.
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
# Inputs:
#   BWA_MEM3 — path to the bwa-mem3 binary under test
set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"
command -v samtools >/dev/null 2>&1 || { echo "SKIP: samtools not on PATH (--meth emits BAM)"; exit 0; }

WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT; cd "$WORK"
fail() { echo "FAIL: $*" >&2; exit 1; }

# Deterministic 1500 bp reference (PRNG seed 4242).
REF=TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTAGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGAATCGTCCCTTCTAACGTAGTTCGCATAGTTCCCGTCCGTAGCCGGACTATTCGAACACCCAGTATTCGATTAACTCGGGCTTGACGTATTAGAGGCGTTAGTGTGCCAGGTAAGATACGCCAACGGAATTAACCTCTGTGACACTCCGCGGAGCCTTCGGACATATAAGTGATCGGGTCTACGTTTGTTAGACTTGAGACGTCTGTTAAGAGTTGGGTCTAATAAATCGCCTACACGTGGAGTCTAACGGGGAAGCGTCGAATCCTGATACATCATATAATGGAGCGTGTTATGAAAAAAGAGCATTCCATTGTACGAGCCGTGCCAGAAACGGCTTGACTACGTGAGCGTAGTGTTAGATAAACAGGAAACACTGACGCGGTTAGAAGGCGGATTGCCGGTAGGTTTTGGAAACATAAATACACACGGTATCATGTTGGGTCACGATTCCTATCACCGCACAGGGCCAACCATAGAAGAACTGAAAGAACTAATCTGGCGGCGGGCTCGGTGCTTATATTTTCCACCCAACATCGTGCACATTAGGCTCACCGCGCCCTACGGGCGAAGGGTGCGTACGGTGTTTATAAGGCGTGACGGCCCCAAGTAGAGGGTAATTCTGTGAAAGAATCTCAGGACGGTGGCATGAATTCAATTCCTTTTAAACCTATCGTTCCGACCTTATGCAATCCTTCAATGAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTACTCTTCACTGCCCAAGCGAGGCTAACCTTTTGAAATGTCACAGTCGAAGCATATCTCCCGTACATCTTTTTCGGAGATCGCAGCTCGCGGAGCTATAAGCGACTTAAGCCCTTGTGTCGGTGATCCCAAGGGTCTGACTCCTGTACCAGGGTTACTGTTTCGCTTTACGGAGTAGCCTGTGAGGTGAACTGAAAGGAGCATATTTGAGATCTAAGATAGGGTCCTCCTCTGCGTCTACGTTCTCTCCGTTACGTACGGCTTCGCACCGGAGTGCATCTTGGCCCCGAAACGCACTGTGTGTGCTGATACAGCGTCCCTGGCCGGCCATGGGTTCAGAACTCCCGGGAACGCTTTTCAACTTAGAGGAACCCCGTCATGGAAGTAGATCGCGTCGAATGAGGGAGTTAGTCCTCGTTCCAGCTGGTAATTGTTTTACCGCTTGGGACCACTATAGGCCGCGGGTAGAAGTTGCTGGGTGTTGATTCCAACCCTCGAACCACGATACGACCTGCCATTTATGGCACAGTAAGGTTCAAACAGCATAATGAATACAGTTATAGTAACTTCCTCACGTACGATTAGGACGCAGCCTTG

printf '>chrA\n%s\n' "$REF" > ref.fa
Q=$(printf 'I%.0s' $(seq 1 60))
emit() { printf '@%s\n%s\n+\n%s\n' "$1" "$2" "$Q" > "$3"; }
"$BWA_MEM3" index --meth ref.fa >/dev/null 2>&1 || fail "index --meth nonzero exit"

tag() { mawk -v k="$2" '{for(i=12;i<=NF;i++) if(substr($i,1,5)==k":Z:"||substr($i,1,5)==k":i:") {print substr($i,6); exit}}' <<<"$1"; }

# --- 1. real SNP vs conversion ---------------------------------------------
# OT fwd @101: 10 C->T conversions (free) + 1 real A->G SNP at read pos 22 (penalized).
SNP_READ=ATTCTGGCGATTTAGGTGACTCGCAGAATCGTTCTTTCTAACGTAGTTTGTATAGTTCTC
emit snp "$SNP_READ" snp.fq
"$BWA_MEM3" mem --meth --meth-scoring genomic -t 1 ref.fa snp.fq > snp.bam 2>/dev/null || fail "snp mem --meth nonzero exit"
samtools quickcheck snp.bam || fail "snp invalid BAM"
line=$(samtools view snp.bam | mawk 'NR==1')
as=$(tag "$line" AS); nm=$(tag "$line" NM); md=$(tag "$line" MD)
# --meth keeps the bwa default mismatch penalty b=4 (not bwameth's b=2), so the single real SNP costs -4
# (conversions free): 59 match/freed columns (+59) - 1 SNP (-4) = 55.
[ "$as" = "55" ] || fail "SNP/conv: AS $as, want 55 (only the real SNP penalized at b=4; conversions free)"
[ "$nm" = "11" ] || fail "SNP/conv: NM $nm, want 11 (10 conversions + 1 SNP)"
# MD reference-base mismatch letters: ten ref-C (conversions) + exactly one ref-A (SNP).
nC=$(printf '%s' "$md" | tr -cd 'C' | wc -c | tr -d ' ')
nA=$(printf '%s' "$md" | tr -cd 'A' | wc -c | tr -d ' ')
nG=$(printf '%s' "$md" | tr -cd 'G' | wc -c | tr -d ' ')
[ "$nC" = "10" ] || fail "SNP/conv: MD has $nC ref-C mismatches, want 10 (conversions); MD=$md"
[ "$nA" = "1" ]  || fail "SNP/conv: MD has $nA ref-A mismatches, want 1 (the real SNP); MD=$md"
[ "$nG" = "0" ]  || fail "SNP/conv: MD has $nG ref-G mismatches, want 0; MD=$md"

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
emit f "$FWD_R1" f1.fq; emit f "$FWD_R2" f2.fq
"$BWA_MEM3" mem --meth --meth-scoring genomic -t 1 ref.fa f1.fq f2.fq > fwd.bam 2>/dev/null || fail "fwd nonzero exit"
check_strand fwd.bam 64  "OT   R1" CT CT
check_strand fwd.bam 128 "CTOT R2" GA CT
emit r "$REV_R1" r1.fq; emit r "$REV_R2" r2.fq
"$BWA_MEM3" mem --meth --meth-scoring genomic -t 1 ref.fa r1.fq r2.fq > rev.bam 2>/dev/null || fail "rev nonzero exit"
check_strand rev.bam 64  "OB   R1" CT GA
check_strand rev.bam 128 "CTOB R2" GA GA

# --- 3. SEQ<->CIGAR orientation on a reverse-mapped mate -------------------
SEQ_R1=ATCCCGGCGACTTAGGCGACCCACAGAATCGTCCCTTCTAACGTAGTTCGCATAGTTCCC
SEQ_R2=TAGGCGATTTATTAGACCCAACTCTTAACAGACGTCTCAAGTCTAACAAACGTAGACCCG
SEQ_R2_RC=CGGGTCTACGTTTGTTAGACTTGAGACGTCTGTTAAGAGTTGGGTCTAATAAATCGCCTA
emit pp "$SEQ_R1" s1.fq; emit pp "$SEQ_R2" s2.fq
"$BWA_MEM3" mem --meth --meth-scoring genomic -t 1 ref.fa s1.fq s2.fq > seq.bam 2>/dev/null || fail "seq nonzero exit"
samtools quickcheck seq.bam || fail "seq invalid BAM (S1 inconsistent-BAM trap)"
read -r rev seq cig < <(samtools view seq.bam | mawk '(int($2/128)%2)==1{print (int($2/16)%2),$10,$6;exit}')
[ "$rev" = "1" ]            || fail "SEQ orient: R2 expected reverse-mapped (rev=$rev)"
[ "$cig" = "60M" ]         || fail "SEQ orient: CIGAR $cig, want 60M"
[ "${#seq}" = "60" ]       || fail "SEQ orient: SEQ length ${#seq}, want 60 (CIGAR consistency)"
[ "$seq" = "$SEQ_R2_RC" ]  || fail "SEQ orient: reverse-mate SEQ is not revcomp(input read)"

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
"$BWA_MEM3" index --meth dup.fa >/dev/null 2>&1 || fail "dup index --meth nonzero exit"
"$BWA_MEM3" mem --meth --meth-scoring genomic -t 1 dup.fa f1.fq f2.fq > dup.bam 2>/dev/null \
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
[ "$mapq2" = "0" ]  || fail "MQ/HN: stale fixture — R2 MAPQ $mapq2, want 0 (locus duplicated on chrB)"
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

echo "PASS: meth_output_integrity (real-SNP vs conversion in AS/NM/MD; Bismark four-strand XR/XG; reverse SEQ orientation; MQ:i/HN:i tag parity)"
