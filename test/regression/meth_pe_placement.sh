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
# genomic strand the mate maps to. Conversions are scored free (AS == 60), and
# because NM is derived from the same matrix they are hidden from NM too
# (NM == 0; issue #327) -- with XM's lowercase call count proving the read really
# carries #conversions > 0, so NM == 0 cannot pass vacuously. Together these
# confirm placement and original-alphabet scoring are jointly correct on every
# strand/hypothesis. MD stays literal on all four records: each must carry one
# MD:Z, and samtools calmd must leave it unchanged while reporting the
# conversion-free NM as different.
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

# Forward fragment: R1 OT fwd @101 (10 C->T), R2 OB rev @301 (5 G->A).
FWD_R1=ATTCTGGCGATTTAGGTGACTCACAGAATCGTTCTTTCTAACGTAGTTTGTATAGTTCTC
FWD_R2=TAAGCGATTTATTAAACCCAACTCTTAACAAACGTCTCAAATCTAACAAACGTAAACCCG
# Reverse fragment: R1 OT rev @701 (6 C->T), R2 OB fwd @501 (9 G->A).
REV_R1=AGATTCTTTCACAGAATTACTCTCTATTTGGGGCTGTCACGCTTTATAAATATCGTACGC
REV_R2=CTAACGCGATTAAAAGACAGATTGCCAGTAAGTTTTAGAAACATAAATACACACAGTATC

printf '>chrA\n%s\n' "$REF" > ref.fa
Q=$(printf 'I%.0s' $(seq 1 60))
emit() { printf '@%s\n%s\n+\n%s\n' "$1" "$2" "$Q" > "$3"; }

"$BWA_MEM3" index --meth ref.fa > /dev/null 2>&1 || fail "index --meth nonzero exit"

# Decode one mate (mate-flag bit 64=READ1, 128=READ2) from $1.bam.
#
# mawk keeps the first matching record and prints it from END rather than
# `exit`ing on the match. Quitting early closed the pipe under samtools, and the
# old `read ... < <(get ...)` discarded samtools' status on top of that -- so a
# truncated or unparseable BAM whose first record happened to satisfy every
# assertion below would pass. Draining the stream lets pipefail see a real
# samtools failure, and check() tests get()'s status before parsing fields.
get() {
    samtools view "$1" | mawk -v bit="$2" '
    out == "" && (int($2/bit) % 2) == 1 {
        rev = (int($2/16) % 2); as=""; nm=""; xr=""; nconv=0;
        for (i=12;i<=NF;i++){ if($i~/^AS:i:/)as=substr($i,6); if($i~/^NM:i:/)nm=substr($i,6); if($i~/^XR:Z:/)xr=substr($i,6); if($i~/^XM:Z:/){xm=substr($i,6); nconv=gsub(/[zxhu]/,"",xm)} }
        out = $3 " " $4 " " $5 " " rev " " $6 " " as " " nm " " xr " " nconv }
    END { if (out != "") print out }'
}
check() { # $1 bam  $2 mateflag  $3 label  $4 pos  $5 rev  $6 xr  $7 as  $8 nconv
    local fields
    fields="$(get "$1" "$2")" || fail "$3: reading $1 failed (samtools/mawk exited non-zero)"
    [ -n "$fields" ] || fail "$3: no record with mate flag $2 in $1"
    read -r rn pos mapq rev cig as nm xr nconv <<< "$fields"
    [ "$rn" = "chrA" ] || fail "$3: RNAME $rn, want chrA"
    [ "$pos" = "$4" ] || fail "$3: POS $pos, want $4 (remap placement)"
    [ "$rev" = "$5" ] || fail "$3: strand rev=$rev, want $5 (reverse-strand re-encode)"
    [ "$cig" = "60M" ] || fail "$3: CIGAR $cig, want 60M"
    [ "$mapq" = "60" ] || fail "$3: MAPQ $mapq, want 60"
    [ "$xr" = "$6" ] || fail "$3: XR $xr, want $6 (hypothesis label)"
    [ "$as" = "$7" ] || fail "$3: AS $as, want $7 (conversions scored free)"
    [ "$nconv" = "$8" ] || fail "$3: XM shows $nconv converted bases, want $8 (fixture must actually exercise conversions)"
    [ "$nconv" -gt 0 ] || fail "$3: converted-base count must be > 0 (test must actually exercise conversions)"
    [ "$nm" = "0" ] || fail "$3: NM $nm, want 0 ($nconv conversions are matrix-freed, so they are not edits for NM)"
}

# MD is literal under --meth: samtools calmd, recomputing MD from ref.fa, must
# leave every record's MD unchanged. calmd adds a missing MD without reporting
# "different MD", so every record must first carry exactly one MD:Z, or the MD
# check would pass on a record that has none. NM excludes conversions, so calmd
# must instead report a different NM on both mates; that also proves calmd
# compared the records rather than passing vacuously.
md_matches_calmd() { # $1 = BAM
    local counts n_rec n_md1
    counts="$(samtools view "$1" | mawk '
        { n = 0; for (i = 12; i <= NF; i++) if ($i ~ /^MD:Z:/) n++; rec++; if (n == 1) ok++ }
        END { print rec + 0, ok + 0 }')" || fail "$1: reading records for the MD:Z count failed"
    read -r n_rec n_md1 <<< "$counts"
    [ "$n_rec" = "2" ] || fail "$1: want 2 records, got $n_rec"
    [ "$n_md1" = "$n_rec" ] || fail "$1: only $n_md1 of $n_rec records carry exactly one MD:Z tag; calmd would add a missing MD silently"
    samtools calmd "$1" ref.fa > /dev/null 2> "$1.calmd.err" || fail "$1: samtools calmd nonzero exit"
    ! grep -q 'different MD' "$1.calmd.err" || fail "$1: MD differs from samtools calmd: $(cat "$1.calmd.err")"
    local n_nm
    n_nm=$(grep -c 'different NM' "$1.calmd.err" || true)
    [ "$n_nm" = "2" ] || fail "$1: want samtools calmd to report NM on both mates, got $n_nm: $(cat "$1.calmd.err")"
}

# --- forward fragment ---
emit f "$FWD_R1" f1.fq
emit f "$FWD_R2" f2.fq
"$BWA_MEM3" mem --meth --meth-scoring genomic -t 1 ref.fa f1.fq f2.fq > fwd.bam 2> /dev/null || fail "fwd mem --meth nonzero exit"
samtools quickcheck fwd.bam || fail "fwd produced an invalid BAM"
check fwd.bam 64 "FWD R1 OT-fwd" 101 0 CT 60 10
check fwd.bam 128 "FWD R2 OB-rev" 301 1 GA 60 5
md_matches_calmd fwd.bam

# --- reverse fragment ---
emit r "$REV_R1" r1.fq
emit r "$REV_R2" r2.fq
"$BWA_MEM3" mem --meth --meth-scoring genomic -t 1 ref.fa r1.fq r2.fq > rev.bam 2> /dev/null || fail "rev mem --meth nonzero exit"
samtools quickcheck rev.bam || fail "rev produced an invalid BAM"
check rev.bam 64 "REV R1 OT-rev" 701 1 CT 60 6
check rev.bam 128 "REV R2 OB-fwd" 501 0 GA 60 9
md_matches_calmd rev.bam

echo "PASS: meth_pe_placement (OT/OB x forward/reverse remap placement + free-conversion scoring + literal MD)"
