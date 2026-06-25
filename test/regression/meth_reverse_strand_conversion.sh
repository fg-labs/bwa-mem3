#!/usr/bin/env bash
# test/regression/meth_reverse_strand_conversion.sh
#
# Regression: the per-strand asymmetric scoring matrix MUST be strand-adjusted
# for reverse-strand extension. bwa-mem extends a reverse-strand seed against the
# reverse-COMPLEMENTED reference window, so a conversion the read's hypothesis
# matrix frees in the forward frame (OT frees ref-C x read-T) presents in the RC
# extension frame as the OTHER cell (ref-G x read-A). Selecting the matrix from
# the raw seed hypothesis with no strand term therefore feeds reverse-strand reads
# the WRONG matrix: their real bisulfite conversions score as mismatches, the
# extension collapses, and the read soft-clips or drops out entirely. The fix
# carries meth_strand_hyp = meth_hypothesis XOR is_rev and selects the matrix from
# it at every extension site.
#
# This guards a case the OT/OB-x-strand placement regression (meth_pe_placement.sh)
# does NOT catch — its reverse-strand reads carry conversions sparse enough that the
# unflipped matrix still scores them free. Here a reverse OT read (R1, XR:CT) carries
# 25 dense bottom-strand conversions: under the FIX it maps full-length (75M) with
# every conversion free; under the unflipped matrix it FAILS TO MAP at all.
#
# Fixture: a real EM-seq read pair from the chr22 holodeck simulation at
# chr22:20564334 (R1 reverse) / chr22:20563995 (R2 forward), with a 601 bp chr22
# window (chr22:20563900-20564500) as the reference, reduced to a self-contained
# inline fixture. The forward mate R2 (OB, XR:GA) maps cleanly under BOTH the buggy
# and fixed binaries — a built-in control showing the defect is reverse-strand only.
#
# Verified empirically: PASS on the fixed binary (d71365b), FAIL on the pre-fix
# binary (498521b, R1 unmapped: CIGAR=*, MAPQ=0).
#
# Inputs:
#   BWA_MEM3 — path to the bwa-mem3 binary under test
set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"
command -v samtools >/dev/null 2>&1 || { echo "SKIP: samtools not on PATH (--meth emits BAM)"; exit 0; }

WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT; cd "$WORK"
fail() { echo "FAIL: $*" >&2; exit 1; }

# 601 bp reference window (chr22:20563900-20564500, GRCh38).
REF=TCAGAAAATCTTGTAATTGACTGTTCATATAAGCAGTTTTCATAATAGCCACAAAGTAGAAGCAGCGCAAATGTTCATCACAGGAAATCGGTGTAGTAGAATATTATTTGGCGATAAAAATAACACAGAAAAACTGATAACATGCCACGACATGGGCGAACTCGTAAATATTATGCTAAAGAAGCTAGTCACAAAAGAGCACATGTTCGATGCCATTCGTACAAAGAGTCCGGAAGTGGCCAATCCATAGAGACAGAGAGGAGATTACTGATTGCCAGAGTCTAGGCTTCTTATTTATCCAAAAGACTTAGTTGTCCCTTTTCTTTTGTCTTTGGTTATTATAGAGTAACTCATGATAGGAAATCCCAAAATCAACACAAATGCTACTTCGTATTCTATCTTTCTGTCTGTGGTAAATGGAACGTTCAGATTCCAGCGGCAGCCGTGGCAGTGGGGCTTTTGCTGGCTGTTTTGTCCCTTGCTGTGCAGCCCTGCAGCGTTTCTGGGAATCTGCCCTGTGGACTGACTGGCGACTCTGGTCTTTTCTCAGCCCAGCTGCAGCTCCAGCAGGTGGCGCTGCAGCAGCAGCAGCAACAGCAGC

# R1: reverse OT read (XR:CT) with 25 dense bottom-strand conversions — the affected read.
R1=TTTTTAGAAACGTTGTAGGGTTGTATAGTAAGGGATACAATAGTTAGTAAAAGTTTTATTGTTACGGTTGTCGTT
# R2: forward OB read (XR:GA) — proper mate + strand-unaffected control.
R2=ATAAAATATTATTTAACAATAAAAATAACACAAAAAAACTAATAACATACCACAACATAAACGAACTCGTAAATA

printf '>chrT\n%s\n' "$REF" > ref.fa
Q=$(printf 'I%.0s' $(seq 1 75))
printf '@p\n%s\n+\n%s\n' "$R1" "$Q" > r1.fq
printf '@p\n%s\n+\n%s\n' "$R2" "$Q" > r2.fq

"$BWA_MEM3" index --meth ref.fa >/dev/null 2>&1 || fail "index --meth nonzero exit"
"$BWA_MEM3" mem --meth --meth-scoring genomic -t 1 ref.fa r1.fq r2.fq > pe.bam 2>/dev/null || fail "mem --meth nonzero exit"
samtools quickcheck pe.bam || fail "produced an invalid BAM"

# Decode one mate (mate-flag bit 64=READ1, 128=READ2) into: rname pos mapq rev cigar as nm xr.
get() { samtools view "$1" | mawk -v bit="$2" '
    (int($2/bit) % 2) == 1 {
        unm = (int($2/4) % 2); rev = (int($2/16) % 2); as=""; nm=""; xr="";
        for (i=12;i<=NF;i++){ if($i~/^AS:i:/)as=substr($i,6); if($i~/^NM:i:/)nm=substr($i,6); if($i~/^XR:Z:/)xr=substr($i,6) }
        print $3, $4, $5, rev, $6, as, nm, xr, unm; exit }'
}

# --- R1: reverse-strand OT read maps full-length with conversions scored free ---
# This is the assertion the bug fails: under the unflipped matrix R1 is unmapped.
read -r rn pos mapq rev cig as nm xr unm < <(get pe.bam 64)
[ "$unm" = "0" ]  || fail "R1 reverse-OT is UNMAPPED (flag 0x4) — the strand-matrix bug: conversions scored as mismatches collapsed the extension"
[ "$rn" = "chrT" ] || fail "R1: RNAME $rn, want chrT"
[ "$pos" = "435" ] || fail "R1: POS $pos, want 435 (reverse-strand placement)"
[ "$rev" = "1" ]   || fail "R1: strand rev=$rev, want 1 (reverse)"
[ "$cig" = "75M" ] || fail "R1: CIGAR $cig, want 75M (no soft-clip — conversions freed by strand-adjusted matrix)"
[ "$mapq" = "60" ] || fail "R1: MAPQ $mapq, want 60"
[ "$xr" = "CT" ]   || fail "R1: XR $xr, want CT (OT hypothesis)"
[ "$nm" = "25" ]   || fail "R1: NM $nm, want 25 (24 conversions + 1 real mismatch; all literal in NM)"
# Conversions are scored free; the 1 real mismatch costs the bwa default penalty b=4 (--meth keeps b=4):
# 74 match/freed columns (+74) - 1 real mismatch (-4) = 70.
[ "$as" = "70" ]   || fail "R1: AS $as, want 70 (conversions free; 1 real mismatch at b=4)"

# --- R2: forward-strand OB mate — strand-unaffected control (passes on buggy + fixed) ---
read -r rn2 pos2 mapq2 rev2 cig2 as2 nm2 xr2 unm2 < <(get pe.bam 128)
[ "$unm2" = "0" ]  || fail "R2 forward-OB should map (control)"
[ "$rn2" = "chrT" ] || fail "R2: RNAME $rn2, want chrT"
[ "$pos2" = "96" ]  || fail "R2: POS $pos2, want 96"
[ "$rev2" = "0" ]   || fail "R2: strand rev=$rev2, want 0 (forward)"
[ "$cig2" = "75M" ] || fail "R2: CIGAR $cig2, want 75M"
[ "$xr2" = "GA" ]   || fail "R2: XR $xr2, want GA (OB hypothesis)"

echo "PASS: meth_reverse_strand_conversion (reverse-strand asymmetric matrix flip; R1 reverse-OT maps full-length, conversions free)"
