#!/usr/bin/env bash
# Synthetic meth option/error checks and single-end success; paired-end stays disabled.
# Usage: bash test/tools/aln_meth_cli_test.sh BWA_MEM3_BINARY
set -euo pipefail
[ $# -eq 1 ] || {
    echo "usage: $0 BWA_MEM3_BINARY" >&2
    exit 2
}
bin=$1
[[ $bin = /* ]] || bin=$PWD/$bin
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
printf '>chr1\nACGTTGCATGTCGCATGATGCATGAGAGCTACGTTGCATGTCGCATGATGCATGAGAGCT\n' > "$tmp/ref.fa"
printf '@r1\nACGTTGCATGTCGCATGATG\n+\nIIIIIIIIIIIIIIIIIIII\n' > "$tmp/r.fq"
"$bin" aln --help > "$tmp/help"
[ "$(grep -c -- '--set-as-failed f|r' "$tmp/help")" = 1 ]
echo 'PASS shared-help-once'
check_error() {
    local name=$1 pattern=$2 rc=0
    shift 2
    "$bin" "$@" > "$tmp/$name.out" 2> "$tmp/$name.err" || rc=$?
    [ "$rc" = 1 ] || {
        echo "FAIL $name: exit $rc" >&2
        exit 1
    }
    grep -F -- "$pattern" "$tmp/$name.err" > /dev/null || {
        echo "FAIL $name: missing diagnostic $pattern" >&2
        sed -n '1,20p' "$tmp/$name.err" >&2
        exit 1
    }
    [ ! -s "$tmp/$name.out" ] || {
        echo "FAIL $name: output before refusal" >&2
        exit 1
    }
    [ ! -e "$tmp/out.sam" ] || {
        echo "FAIL $name: created output file" >&2
        exit 1
    }
    echo "PASS $name"
}
check_error stray 'use --meth=taps' aln --meth taps "$tmp/ref.fa" "$tmp/r.fq"
check_error stray-tags '--meth-tags' aln --meth --meth-tags XR XG "$tmp/ref.fa" "$tmp/r.fq"
check_error invalid-failed "requires 'f' or 'r'" aln --meth --set-as-failed x "$tmp/ref.fa" "$tmp/r.fq"
check_error mem-invalid-failed "requires 'f' or 'r'" mem --meth --set-as-failed x "$tmp/ref.fa" "$tmp/r.fq"
check_error failed-without-meth '--set-as-failed requires --meth' aln --set-as-failed f "$tmp/ref.fa" "$tmp/r.fq"
for bam in '' --bam; do
    check_error "sq${bam}" 'ERROR: --meth does not accept @SQ' aln --meth ${bam:+$bam} -f "$tmp/out.sam" \
        -r '@RG\tID:x\n@SQ\tSN:chr1\tLN:56' "$tmp/ref.fa" "$tmp/r.fq"
done
check_error missing-original 'original reference index file' aln --meth -f "$tmp/out.sam" "$tmp/ref.fa" "$tmp/r.fq"
grep -F 'bwa-mem3 index --meth' "$tmp/missing-original.err" > /dev/null
# A bogus converted index must not be opened when the original reference is missing.
printf 'bad index\n' > "$tmp/ref.fa.meth.bwt.2bit.64"
check_error original-before-converted 'original reference index file' aln --meth "$tmp/ref.fa.meth" "$tmp/r.fq"
rm "$tmp/ref.fa.meth.bwt.2bit.64"
"$bin" index "$tmp/ref.fa" > "$tmp/index.out" 2> "$tmp/index.err"
check_error missing-converted 'converted index' aln --meth -f "$tmp/out.sam" "$tmp/ref.fa" "$tmp/r.fq"
grep -F 'bwa-mem3 index --meth' "$tmp/missing-converted.err" > /dev/null
# shellcheck disable=SC2016 # Match literal diagnostic backticks, not command substitution.
if grep -F 'run `bwa-mem3 index`' "$tmp/missing-converted.err"; then exit 1; fi
"$bin" index --meth "$tmp/ref.fa" > "$tmp/index-meth.out" 2> "$tmp/index-meth.err"
for kind in original converted; do
    prefix=$tmp/ref.fa
    if [ "$kind" = converted ]; then prefix=$prefix.meth; fi
    for suffix in ann amb pac; do
        mv "$prefix.$suffix" "$prefix.$suffix.saved"
        check_error "$kind-$suffix" "$kind" aln --meth "$tmp/ref.fa" "$tmp/r.fq"
        grep -F 'bwa-mem3 index --meth' "$tmp/$kind-$suffix.err" > /dev/null
        mv "$prefix.$suffix.saved" "$prefix.$suffix"
    done
done
for chem in emseq taps 5base; do
    "$bin" aln "--meth=$chem" "$tmp/ref.fa" "$tmp/r.fq" > "$tmp/se-$chem.sam" 2> "$tmp/se-$chem.err"
    grep -F 'XR:Z:CT' "$tmp/se-$chem.sam" > /dev/null
    echo "PASS se-$chem"
done
"$bin" aln --meth "$tmp/ref.fa.meth" "$tmp/r.fq" > "$tmp/suffix.sam" 2> "$tmp/suffix.err"
grep -v '^@' "$tmp/se-emseq.sam" > "$tmp/base.records"
grep -v '^@' "$tmp/suffix.sam" > "$tmp/suffix.records"
cmp "$tmp/base.records" "$tmp/suffix.records"
echo 'PASS suffix-equivalence'
for value in f r; do
    "$bin" aln --meth --set-as-failed "$value" "$tmp/ref.fa" "$tmp/r.fq" > "$tmp/failed-$value.sam" 2> "$tmp/failed-$value.err"
    grep -F 'XR:Z:CT' "$tmp/failed-$value.sam" > /dev/null
    echo "PASS failed-$value"
done
check_error paired-refusal 'paired-end --meth is not yet supported' aln --meth "$tmp/ref.fa" "$tmp/r.fq" "$tmp/r.fq"
check_error cpg-refusal '--meth-context cpg is not yet supported' aln --meth=taps --meth-context cpg "$tmp/ref.fa" "$tmp/r.fq"
cp "$tmp/r.fq" "$tmp/taps"
(cd "$tmp" && "$bin" aln --meth "$tmp/ref.fa" taps > "$tmp/path.sam" 2> "$tmp/path.err")
grep -v '^@' "$tmp/path.sam" > "$tmp/path.records"
cmp "$tmp/base.records" "$tmp/path.records"
echo 'PASS real-path'
# Reject a converted index whose annotation still parses but has the wrong contig name.
cp "$tmp/ref.fa.meth.ann" "$tmp/converted.ann.saved"
sed 's/rchr1/rwrong/' "$tmp/converted.ann.saved" > "$tmp/ref.fa.meth.ann"
check_error incompatible-index 'converted index does not match' aln --meth "$tmp/ref.fa" "$tmp/r.fq"
cp "$tmp/converted.ann.saved" "$tmp/ref.fa.meth.ann"
