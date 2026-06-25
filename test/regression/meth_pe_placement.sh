#!/usr/bin/env bash
# test/regression/meth_pe_placement.sh
#
# Regression (D3 B4): the seed->original remap places paired-end --meth reads at
# the correct ORIGINAL locus + strand + bisulfite hypothesis for ALL FOUR
# combinations of {OT, OB} x {forward, reverse}. The reverse-strand re-encode in
# meth_seed_to_orig ((orig_bns->l_pac<<1)-1-orig_fwd) is the single most
# placement-critical step in --meth; this exercises it on both strands for both
# hypotheses via the live `mem --meth` path.
#
#   Forward fragment:  R1 OT forward @101,  R2 OB reverse @301
#   Reverse fragment:  R1 OT reverse @701,  R2 OB forward @501
#
# Directional contract: R1 -> OT (XR:CT), R2 -> OB (XR:GA), regardless of which
# genomic strand the mate maps to. Conversions are scored free (AS == 60) but
# remain literal mismatches (NM == #conversions), confirming placement and
# original-alphabet scoring are jointly correct on every strand/hypothesis.
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

# Forward fragment: R1 OT fwd @101 (10 C->T), R2 OB rev @301 (5 G->A).
FWD_R1=ATTCTGGCGATTTAGGTGACTCACAGAATCGTTCTTTCTAACGTAGTTTGTATAGTTCTC
FWD_R2=TAAGCGATTTATTAAACCCAACTCTTAACAAACGTCTCAAATCTAACAAACGTAAACCCG
# Reverse fragment: R1 OT rev @701 (6 C->T), R2 OB fwd @501 (9 G->A).
REV_R1=AGATTCTTTCACAGAATTACTCTCTATTTGGGGCTGTCACGCTTTATAAATATCGTACGC
REV_R2=CTAACGCGATTAAAAGACAGATTGCCAGTAAGTTTTAGAAACATAAATACACACAGTATC

printf '>chrA\n%s\n' "$REF" > ref.fa
Q=$(printf 'I%.0s' $(seq 1 60))
emit() { printf '@%s\n%s\n+\n%s\n' "$1" "$2" "$Q" > "$3"; }

"$BWA_MEM3" index --meth ref.fa >/dev/null 2>&1 || fail "index --meth nonzero exit"

# Decode one mate (mate-flag bit 64=READ1, 128=READ2) from $1.bam.
get() { samtools view "$1" | mawk -v bit="$2" '
    (int($2/bit) % 2) == 1 {
        rev = (int($2/16) % 2); as=""; nm=""; xr="";
        for (i=12;i<=NF;i++){ if($i~/^AS:i:/)as=substr($i,6); if($i~/^NM:i:/)nm=substr($i,6); if($i~/^XR:Z:/)xr=substr($i,6) }
        print $3, $4, $5, rev, $6, as, nm, xr; exit }'
}
check() { # $1 bam  $2 mateflag  $3 label  $4 pos  $5 rev  $6 xr  $7 as  $8 nm
    read -r rn pos mapq rev cig as nm xr < <(get "$1" "$2")
    [ "$rn" = "chrA" ]  || fail "$3: RNAME $rn, want chrA"
    [ "$pos" = "$4" ]   || fail "$3: POS $pos, want $4 (remap placement)"
    [ "$rev" = "$5" ]   || fail "$3: strand rev=$rev, want $5 (reverse-strand re-encode)"
    [ "$cig" = "60M" ]  || fail "$3: CIGAR $cig, want 60M"
    [ "$mapq" = "60" ]  || fail "$3: MAPQ $mapq, want 60"
    [ "$xr" = "$6" ]    || fail "$3: XR $xr, want $6 (hypothesis label)"
    [ "$as" = "$7" ]    || fail "$3: AS $as, want $7 (conversions scored free)"
    [ "$nm" = "$8" ]    || fail "$3: NM $nm, want $8 (conversions remain literal mismatches)"
}

# --- forward fragment ---
emit f "$FWD_R1" f1.fq; emit f "$FWD_R2" f2.fq
"$BWA_MEM3" mem --meth --meth-scoring genomic -t 1 ref.fa f1.fq f2.fq > fwd.bam 2>/dev/null || fail "fwd mem --meth nonzero exit"
samtools quickcheck fwd.bam || fail "fwd produced an invalid BAM"
check fwd.bam 64  "FWD R1 OT-fwd" 101 0 CT 60 10
check fwd.bam 128 "FWD R2 OB-rev" 301 1 GA 60 5

# --- reverse fragment ---
emit r "$REV_R1" r1.fq; emit r "$REV_R2" r2.fq
"$BWA_MEM3" mem --meth --meth-scoring genomic -t 1 ref.fa r1.fq r2.fq > rev.bam 2>/dev/null || fail "rev mem --meth nonzero exit"
samtools quickcheck rev.bam || fail "rev produced an invalid BAM"
check rev.bam 64  "REV R1 OT-rev" 701 1 CT 60 6
check rev.bam 128 "REV R2 OB-fwd" 501 0 GA 60 9

echo "PASS: meth_pe_placement (OT/OB x forward/reverse remap placement + free-conversion scoring)"
