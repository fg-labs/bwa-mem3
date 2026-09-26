#!/usr/bin/env bash
# test/regression/meth_nul_base.sh
#
# Regression: under --meth, a control byte (0, TAB, ...) in a read's sequence
# line must not corrupt the record. The FASTQ reader keeps the line's full
# parsed length, but --meth saved the original bases with strdup(), which stops
# at the first 0: extension, mate rescue and the SEQ writer then read l_seq
# bytes of a shorter copy -- a heap over-read that also wrote the bytes past it
# into SEQ (a 100 bp read came out as 20 real bases plus 80 N) and aligned
# against them. The same bytes went into the read's YS:Z/YC:Z carrier comment,
# which is read as a tab-separated C string, so a 0 byte silently dropped XR:Z
# and every -C tag, and a TAB split the carrier and could forge a tag. --meth
# now writes each control byte as N before copying anything.
#
# The fixture is built programmatically: a constructed bisulfite reference and
# three PE pairs aligned with -C, every read carrying a BC:Z FASTQ comment.
#   p1: R1 (forward, OT) has a 0 byte at base 20, where the reference has A.
#       It keeps its methylated CpGs as C, so its original bases differ from
#       the converted ones and SEQ restoration is tested, not just length.
#   p2: R2 (maps to the reverse strand, OB) has a 0 byte at base 31 and a TAB
#       at base 61.
#   p3: R1 (forward, OT) has a TAB at base 86 followed by the tag-shaped bytes
#       ZZ:Z:ACGT, so a TAB left in the carrier would split off a well-formed
#       ZZ:Z tag and emit it (the writer silently drops a split-off token that
#       is not tag-shaped, so p2's TAB alone cannot show a missed TAB).
# For each affected read: it is mapped, on its expected strand; SEQ equals the
# input with control bytes (and p3's non-base bytes) read as N; the alignment
# scores that base as N too (p1's NM/MD record a mismatch against the
# reference A, where a raw 0 byte would be scored as base A and match); the
# record carries XR:Z (CT for R1, GA for R2) and BC:Z:ACGT; the internal
# YS:Z/YC:Z carriers do not leak; no forged ZZ:Z tag appears; and every aux
# field is a well-formed tag.
#
# Inputs:
#   BWA_MEM3 — path to the bwa-mem3 binary under test
set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"
command -v samtools > /dev/null 2>&1 || {
    echo "SKIP: samtools not on PATH (needed to filter records by FLAG)"
    exit 0
}

case "$BWA_MEM3" in
    /*) BIN="$BWA_MEM3" ;;
    *) BIN="$PWD/$BWA_MEM3" ;;
esac
[ -x "$BIN" ] || {
    echo "FAIL: BWA_MEM3 ($BIN) is not executable" >&2
    exit 1
}

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
fail() {
    echo "FAIL: $*" >&2
    exit 1
}

# Deterministic 1500 bp reference (same PRNG seed 4242 as the sibling meth tests).
REF=TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTAGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGAATCGTCCCTTCTAACGTAGTTCGCATAGTTCCCGTCCGTAGCCGGACTATTCGAACACCCAGTATTCGATTAACTCGGGCTTGACGTATTAGAGGCGTTAGTGTGCCAGGTAAGATACGCCAACGGAATTAACCTCTGTGACACTCCGCGGAGCCTTCGGACATATAAGTGATCGGGTCTACGTTTGTTAGACTTGAGACGTCTGTTAAGAGTTGGGTCTAATAAATCGCCTACACGTGGAGTCTAACGGGGAAGCGTCGAATCCTGATACATCATATAATGGAGCGTGTTATGAAAAAAGAGCATTCCATTGTACGAGCCGTGCCAGAAACGGCTTGACTACGTGAGCGTAGTGTTAGATAAACAGGAAACACTGACGCGGTTAGAAGGCGGATTGCCGGTAGGTTTTGGAAACATAAATACACACGGTATCATGTTGGGTCACGATTCCTATCACCGCACAGGGCCAACCATAGAAGAACTGAAAGAACTAATCTGGCGGCGGGCTCGGTGCTTATATTTTCCACCCAACATCGTGCACATTAGGCTCACCGCGCCCTACGGGCGAAGGGTGCGTACGGTGTTTATAAGGCGTGACGGCCCCAAGTAGAGGGTAATTCTGTGAAAGAATCTCAGGACGGTGGCATGAATTCAATTCCTTTTAAACCTATCGTTCCGACCTTATGCAATCCTTCAATGAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTACTCTTCACTGCCCAAGCGAGGCTAACCTTTTGAAATGTCACAGTCGAAGCATATCTCCCGTACATCTTTTTCGGAGATCGCAGCTCGCGGAGCTATAAGCGACTTAAGCCCTTGTGTCGGTGATCCCAAGGGTCTGACTCCTGTACCAGGGTTACTGTTTCGCTTTACGGAGTAGCCTGTGAGGTGAACTGAAAGGAGCATATTTGAGATCTAAGATAGGGTCCTCCTCTGCGTCTACGTTCTCTCCGTTACGTACGGCTTCGCACCGGAGTGCATCTTGGCCCCGAAACGCACTGTGTGTGCTGATACAGCGTCCCTGGCCGGCCATGGGTTCAGAACTCCCGGGAACGCTTTTCAACTTAGAGGAACCCCGTCATGGAAGTAGATCGCGTCGAATGAGGGAGTTAGTCCTCGTTCCAGCTGGTAATTGTTTTACCGCTTGGGACCACTATAGGCCGCGGGTAGAAGTTGCTGGGTGTTGATTCCAACCCTCGAACCACGATACGACCTGCCATTTATGGCACAGTAAGGTTCAAACAGCATAATGAATACAGTTATAGTAACTTCCTCACGTACGATTAGGACGCAGCCTTG

printf '>chrA\n%s\n' "$REF" > ref.fa
"$BIN" index --meth ref.fa > /dev/null 2>&1 || fail "index --meth nonzero exit"

# Three FR pairs from the reference, each R1 = ref[a:a+100] and
# R2 = revcomp(ref[b:b+100]), bisulfite-converted (R1 C->T except CpG, R2 G->A)
# so they map under --meth. Written with python so the control bytes land in
# the files; expected.tsv holds each affected read's SEQ in read order with
# those bytes as N, plus the strand it must map to (16 = reverse).
python3 - "$REF" << 'PY' || fail "fixture generation failed"
import sys
ref = sys.argv[1]
comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
def rc(s):
    return ''.join(comp[b] for b in reversed(s))
def convert_ot(s):
    # C->T, except a C followed by G (a methylated CpG) stays C.
    return ''.join('T' if b == 'C' and (i + 1 == len(s) or s[i + 1] != 'G') else b
                   for i, b in enumerate(s))
def put(s, i, byte):
    return s[:i] + byte + s[i + 1:]
q = 'I' * 100
assert ref[219] == 'A'
p1_r1 = put(convert_ot(ref[200:300]), 19, '\x00')
p1_r2 = rc(ref[450:550]).replace('G', 'A')
p2_r1 = ref[700:800].replace('C', 'T')
p2_r2 = put(put(rc(ref[950:1050]).replace('G', 'A'), 30, '\x00'), 60, '\t')
forged = '\tZZ:Z:ACGT'
p3_r1 = convert_ot(ref[1150:1250])
p3_r1 = p3_r1[:85] + forged + p3_r1[85 + len(forged):]
p3_r2 = rc(ref[1380:1480]).replace('G', 'A')
assert len(p3_r1) == 100
assert 'C' in p1_r1  # restoration must be observable
with open('r1.fq', 'wb') as f1, open('r2.fq', 'wb') as f2:
    for name, r1, r2 in (('p1', p1_r1, p1_r2), ('p2', p2_r1, p2_r2), ('p3', p3_r1, p3_r2)):
        f1.write(('@%s/1\tBC:Z:ACGT\n%s\n+\n%s\n' % (name, r1, q)).encode('latin-1'))
        f2.write(('@%s/2\tBC:Z:ACGT\n%s\n+\n%s\n' % (name, r2, q)).encode('latin-1'))
def as_n(s):
    return ''.join(b if b in 'ACGTN' else 'N' for b in s)
with open('expected.tsv', 'w') as f:
    f.write('p1\t64\t0\tCT\t%s\n' % as_n(p1_r1))
    f.write('p2\t128\t16\tGA\t%s\n' % as_n(p2_r2))
    f.write('p3\t64\t0\tCT\t%s\n' % as_n(p3_r1))
PY

"$BIN" mem --meth -C ref.fa r1.fq r2.fq > out.sam 2> mem.log || {
    cat mem.log >&2
    fail "mem --meth nonzero exit"
}
samtools view -F 0x900 out.sam > primary.sam || fail "samtools view failed on the output"

# check_read NAME READ_FLAG STRAND XR EXPECTED_SEQ: the primary record of that
# read, which must be mapped with FLAG & 0x10 equal to STRAND.
check_read() {
    local name=$1 read_flag=$2 strand=$3 xr=$4 want=$5
    local line="" qname flag rest
    while IFS=$'\t' read -r qname flag rest; do
        if [ "$qname" = "$name" ] && ((flag & read_flag)); then
            [ -z "$line" ] || fail "$name: more than one primary record for read flag $read_flag"
            line="$qname"$'\t'"$flag"$'\t'"$rest"
        fi
    done < primary.sam
    [ -n "$line" ] || fail "$name: no primary record for read flag $read_flag"

    local seq
    flag=$(printf '%s\n' "$line" | cut -f2)
    (((flag & 4) == 0)) || fail "$name: read flag $read_flag is unmapped (FLAG $flag)"
    (((flag & 16) == strand)) || fail "$name: read flag $read_flag mapped to the wrong strand (FLAG $flag, want FLAG & 0x10 = $strand)"
    seq=$(printf '%s\n' "$line" | cut -f10)
    # SEQ is stored reverse-complemented for a reverse-strand hit; compare in
    # read order.
    if ((flag & 16)); then
        seq=$(printf '%s\n' "$seq" | rev | tr ACGTN TGCAN)
    fi
    [ "$seq" = "$want" ] || fail "$name: SEQ differs from the input (control bytes as N): got $seq, want $want"

    local tags field
    tags=$'\t'$(printf '%s\n' "$line" | cut -f12-)$'\t'
    [[ $tags == *$'\t'"XR:Z:$xr"$'\t'* ]] || fail "$name: missing XR:Z:$xr; tags:$tags"
    [[ $tags == *$'\tBC:Z:ACGT\t'* ]] || fail "$name: lost the -C tag BC:Z:ACGT; tags:$tags"
    [[ $tags != *$'\tYS:Z:'* && $tags != *$'\tYC:Z:'* ]] || fail "$name: internal YS:Z/YC:Z carrier leaked; tags:$tags"
    [[ $tags != *$'\tZZ:'* ]] || fail "$name: a TAB in the carrier forged a ZZ tag; tags:$tags"
    while IFS= read -r field; do
        [[ $field =~ ^[A-Za-z][A-Za-z0-9]:[AifZHB]: ]] || fail "$name: malformed aux field '$field'"
    done < <(printf '%s\n' "$line" | cut -f12- | tr '\t' '\n')
}

while IFS=$'\t' read -r name read_flag strand xr want; do
    check_read "$name" "$read_flag" "$strand" "$xr" "$want"
done < expected.tsv

# p1 R1 is otherwise an exact match, so the N at base 20 is its only edit.
while IFS=$'\t' read -r qname flag rest; do
    if [ "$qname" != p1 ] || ! ((flag & 64)); then continue; fi
    tags=$'\t'$(printf '%s\n' "$rest" | cut -f10-)$'\t'
    [[ $tags == *$'\tNM:i:1\t'* && $tags == *$'\tMD:Z:19A80\t'* ]] \
        || fail "p1: the control byte was not scored as N (want NM:i:1 MD:Z:19A80); tags:$tags"
done < primary.sam
echo "PASS: meth_nul_base"
