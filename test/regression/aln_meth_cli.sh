#!/usr/bin/env bash
# test/regression/aln_meth_cli.sh
#
# `bwa-mem3 aln --meth` from the command line, on a synthetic reference: option and index errors
# (each refused before any output), prefix resolution, output aliasing of index files and header
# sidecars, single-end records for each chemistry and --set-as-failed, -v, paired-end records and
# their name and FASTQ errors, --meth-context cpg and its option checks, and, on Linux with nothing
# else staged, a converted index attached from shared memory with its disk files gone.
#
# Env:
#   BWA_MEM3  the binary (default: the repository's ./bwa-mem3)
set -euo pipefail
# A bare assertion that fails still names itself.
trap 'echo "FAIL: line $LINENO: $BASH_COMMAND" >&2' ERR
bin="${BWA_MEM3:-$(dirname "$0")/../../bwa-mem3}"
[[ $bin = /* ]] || bin=$PWD/$bin
[[ -x $bin ]] || {
    echo "FAIL: BWA_MEM3 ('$bin') is not an executable binary" >&2
    exit 1
}
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
printf '>chr1\nACGTTGCATGTCGCATGATGCATGAGAGCTACGTTGCATGTCGCATGATGCATGAGAGCT\n' > "$tmp/ref.fa"
printf '@r1\nACGTTGCATGTCGCATGATG\n+\nIIIIIIIIIIIIIIIIIIII\n' > "$tmp/r.fq"
"$bin" aln --help > "$tmp/help"
[ "$(grep -c -- '--set-as-failed f|r' "$tmp/help")" = 1 ] || {
    echo "FAIL: shared-help-once: expected one --set-as-failed entry" >&2
    exit 1
}
echo 'PASS: shared-help-once'
# The one alignment record of a run must be r1 mapped forward at either copy (chr1:1 or chr1:31,
# a two-copy repeat) as 20M, with the expected flag.
check_record() {
    local name=$1 sam=$2 flag=$3
    awk -v flag="$flag" '!/^@/ { n++; if ($1 != "r1" || $2 != flag || $3 != "chr1" || ($4 != 1 && $4 != 31) || $6 != "20M" || $0 !~ /\tX0:i:2\t/ || $0 !~ /\tXR:Z:CT/) bad = 1 }
        END { exit !(n == 1 && !bad) }' "$sam" || {
        echo "FAIL: $name: expected one r1 record with flag $flag at chr1:1 or chr1:31 (20M, X0:i:2, XR:Z:CT)" >&2
        grep -v '^@' "$sam" >&2 || true
        exit 1
    }
}
check_error() {
    local name=$1 pattern=$2 rc=0
    shift 2
    "$bin" "$@" > "$tmp/$name.out" 2> "$tmp/$name.err" || rc=$?
    [ "$rc" = 1 ] || {
        echo "FAIL: $name: exit $rc" >&2
        exit 1
    }
    grep -F -- "$pattern" "$tmp/$name.err" > /dev/null || {
        echo "FAIL: $name: missing diagnostic $pattern" >&2
        sed -n '1,20p' "$tmp/$name.err" >&2
        exit 1
    }
    [ ! -s "$tmp/$name.out" ] || {
        echo "FAIL: $name: output before refusal" >&2
        exit 1
    }
    [ ! -e "$tmp/out.sam" ] || {
        echo "FAIL: $name: created output file" >&2
        exit 1
    }
    echo "PASS: $name"
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
grep -F 'bwa-mem3 index --meth' "$tmp/missing-original.err" > /dev/null || {
    echo "FAIL: missing-original: missing converted-index command" >&2
    exit 1
}
# A bogus converted index must not be opened when the original reference is missing.
printf 'bad index\n' > "$tmp/ref.fa.meth.bwt.2bit.64"
check_error original-before-converted 'original reference index file' aln --meth "$tmp/ref.fa.meth" "$tmp/r.fq"
rm "$tmp/ref.fa.meth.bwt.2bit.64"
"$bin" index "$tmp/ref.fa" > "$tmp/index.out" 2> "$tmp/index.err"
check_error missing-converted 'converted index' aln --meth -f "$tmp/out.sam" "$tmp/ref.fa" "$tmp/r.fq"
grep -F 'bwa-mem3 index --meth' "$tmp/missing-converted.err" > /dev/null || {
    echo "FAIL: missing-converted: missing converted-index command" >&2
    exit 1
}
# shellcheck disable=SC2016 # Match literal diagnostic backticks, not command substitution.
if grep -F 'run `bwa-mem3 index`' "$tmp/missing-converted.err"; then
    echo 'FAIL: missing-converted: recommends a plain index' >&2
    exit 1
fi
"$bin" index --meth "$tmp/ref.fa" > "$tmp/index-meth.out" 2> "$tmp/index-meth.err"
for kind in original converted; do
    prefix=$tmp/ref.fa
    if [ "$kind" = converted ]; then prefix=$prefix.meth; fi
    for suffix in ann amb pac; do
        mv "$prefix.$suffix" "$prefix.$suffix.saved"
        check_error "$kind-$suffix" "$kind" aln --meth "$tmp/ref.fa" "$tmp/r.fq"
        grep -F 'bwa-mem3 index --meth' "$tmp/$kind-$suffix.err" > /dev/null || {
            echo "FAIL: $kind-$suffix: missing converted-index command" >&2
            exit 1
        }
        mv "$prefix.$suffix.saved" "$prefix.$suffix"
    done
done
for chem in emseq taps 5base; do
    "$bin" aln "--meth=$chem" "$tmp/ref.fa" "$tmp/r.fq" > "$tmp/se-$chem.sam" 2> "$tmp/se-$chem.err"
    check_record "se-$chem" "$tmp/se-$chem.sam" 0
    echo "PASS: se-$chem"
done
"$bin" aln --meth "$tmp/ref.fa.meth" "$tmp/r.fq" > "$tmp/suffix.sam" 2> "$tmp/suffix.err"
grep -v '^@' "$tmp/se-emseq.sam" > "$tmp/base.records" || {
    echo "FAIL: suffix-equivalence: no baseline alignment records" >&2
    exit 1
}
grep -v '^@' "$tmp/suffix.sam" > "$tmp/suffix.records" || {
    echo "FAIL: suffix-equivalence: no converted-prefix alignment records" >&2
    exit 1
}
cmp "$tmp/base.records" "$tmp/suffix.records" || {
    echo "FAIL: suffix-equivalence: alignment records differ" >&2
    exit 1
}
echo 'PASS: suffix-equivalence'
for level in 0 4; do
    "$bin" aln --meth -v "$level" "$tmp/ref.fa" "$tmp/r.fq" > "$tmp/verbose.sam" 2> "$tmp/verbose.err"
    grep -v '^@' "$tmp/verbose.sam" > "$tmp/verbose.records" || {
        echo "FAIL: verbosity-$level: no alignment records" >&2
        exit 1
    }
    cmp "$tmp/base.records" "$tmp/verbose.records" || {
        echo "FAIL: verbosity-$level: alignment records differ" >&2
        exit 1
    }
    if [ "$level" = 4 ]; then
        grep -E '[1-9][0-9]* SA lookups.*[1-9][0-9]* unique windows scored' "$tmp/verbose.err" > /dev/null || {
            echo "FAIL: verbosity-$level: missing lookup/window diagnostic" >&2
            exit 1
        }
    elif grep -q 'SA lookups' "$tmp/verbose.err"; then
        echo "FAIL: verbosity-$level: SA lookup report printed" >&2
        exit 1
    fi
    echo "PASS: verbosity-$level"
done
for value in -1 5 abc; do
    check_error "verbosity-invalid-$value" '-v' aln -v "$value" "$tmp/ref.fa" "$tmp/r.fq"
done
for value in f r; do
    "$bin" aln --meth --set-as-failed "$value" "$tmp/ref.fa" "$tmp/r.fq" > "$tmp/failed-$value.sam" 2> "$tmp/failed-$value.err"
    # r1 is a forward read-1 placement: --set-as-failed f marks it QC-fail (0x200), r does not.
    if [ "$value" = f ]; then want=512; else want=0; fi
    check_record "failed-$value" "$tmp/failed-$value.sam" "$want"
    echo "PASS: failed-$value"
done
"$bin" aln --meth "$tmp/ref.fa" "$tmp/r.fq" "$tmp/r.fq" > "$tmp/paired.sam" 2> "$tmp/paired.err"
# Both mates are written: one r1 record with the read-1 flag (XR:Z:CT) and one with the read-2
# flag (XR:Z:GA).
awk '!/^@/ { n++
        if ($1 != "r1") bad = 1
        else if (int($2 / 64) % 4 == 1 && $0 ~ /\tXR:Z:CT(\t|$)/) r1++
        else if (int($2 / 64) % 4 == 2 && $0 ~ /\tXR:Z:GA(\t|$)/) r2++
        else bad = 1 }
    END { exit !(n == 2 && r1 == 1 && r2 == 1 && !bad) }' "$tmp/paired.sam" || {
    echo 'FAIL: paired-meth: expected one r1 read-1 record (XR:Z:CT) and one r1 read-2 record (XR:Z:GA)' >&2
    grep -v '^@' "$tmp/paired.sam" >&2 || true
    exit 1
}
echo 'PASS: paired-meth'
printf '@wrong\nACGTTGCATGTCGCATGATG\n+\nIIIIIIIIIIIIIIIIIIII\n' > "$tmp/wrong.fq"
rc=0
"$bin" aln --meth "$tmp/ref.fa" "$tmp/r.fq" "$tmp/wrong.fq" > "$tmp/wrong.sam" 2> "$tmp/wrong.err" || rc=$?
[ "$rc" = 1 ] || {
    echo "FAIL: paired-name-error: exit $rc (expected 1)" >&2
    exit 1
}
grep -F 'paired reads have different names' "$tmp/wrong.err" > /dev/null || {
    echo "FAIL: paired-name-error: missing name-mismatch diagnostic" >&2
    exit 1
}
echo 'PASS: paired-name-error'
printf '@r1\nACGT\n+\nII\n' > "$tmp/malformed.fq"
rc=0
"$bin" aln --meth "$tmp/ref.fa" "$tmp/r.fq" "$tmp/malformed.fq" > "$tmp/malformed.sam" 2> "$tmp/malformed.err" || rc=$?
[ "$rc" = 1 ] || {
    echo "FAIL: paired-malformed: exit $rc (expected 1)" >&2
    exit 1
}
grep -F 'malformed' "$tmp/malformed.err" > /dev/null || {
    echo "FAIL: paired-malformed: missing malformed-input diagnostic" >&2
    exit 1
}
echo 'PASS: paired-malformed'
"$bin" aln --meth=taps --meth-context cpg "$tmp/ref.fa" "$tmp/r.fq" > "$tmp/cpg.sam" 2> "$tmp/cpg.err"
grep -F 'XG:Z:' "$tmp/cpg.sam" > /dev/null || {
    echo "FAIL: cpg-success: missing XG tag" >&2
    exit 1
}
echo 'PASS: cpg-success'
check_error cpg-emseq 'requires --meth=taps' aln --meth=emseq --meth-context cpg "$tmp/ref.fa" "$tmp/r.fq"
check_error cpg-collapsed 'incompatible with --meth-scoring collapsed' aln --meth=taps --meth-context cpg --meth-scoring collapsed "$tmp/ref.fa" "$tmp/r.fq"
cp "$tmp/r.fq" "$tmp/taps"
(cd "$tmp" && "$bin" aln --meth "$tmp/ref.fa" taps > "$tmp/path.sam" 2> "$tmp/path.err")
grep -v '^@' "$tmp/path.sam" > "$tmp/path.records" || {
    echo "FAIL: real-path: no alignment records" >&2
    exit 1
}
cmp "$tmp/base.records" "$tmp/path.records" || {
    echo "FAIL: real-path: alignment records differ" >&2
    exit 1
}
echo 'PASS: real-path'
# Reject a converted index whose annotation still parses but has the wrong contig name.
cp "$tmp/ref.fa.meth.ann" "$tmp/converted.ann.saved"
sed 's/rchr1/rwrong/' "$tmp/converted.ann.saved" > "$tmp/ref.fa.meth.ann"
check_error incompatible-index 'converted index does not match' aln --meth "$tmp/ref.fa" "$tmp/r.fq"
cp "$tmp/converted.ann.saved" "$tmp/ref.fa.meth.ann"
# -f naming an original-reference index file is refused before the index is read.
cp "$tmp/ref.fa.pac" "$tmp/pac.saved"
check_error out-original-index 'is the index file' aln --meth -f "$tmp/ref.fa.pac" "$tmp/ref.fa" "$tmp/r.fq"
cmp "$tmp/ref.fa.pac" "$tmp/pac.saved" || {
    echo "FAIL: out-original-index: original index changed" >&2
    exit 1
}
# A reference whose own name ends in .meth is found under its own prefix.
cp "$tmp/ref.fa" "$tmp/x.meth"
"$bin" index --meth "$tmp/x.meth" > "$tmp/index-x.out" 2> "$tmp/index-x.err"
"$bin" aln --meth "$tmp/x.meth" "$tmp/r.fq" > "$tmp/x.sam" 2> "$tmp/x.err"
check_record meth-named-reference "$tmp/x.sam" 0
echo 'PASS: meth-named-reference'
# Once a reference x is indexed as well, x.meth names either reference: refuse to pick one.
cp "$tmp/ref.fa" "$tmp/x"
"$bin" index "$tmp/x" > "$tmp/index-xx.out" 2> "$tmp/index-xx.err"
check_error ambiguous-prefix 'rename one of them' aln --meth "$tmp/x.meth" "$tmp/r.fq"
"$bin" aln --meth "$tmp/x.meth.meth" "$tmp/r.fq" > "$tmp/x2.sam" 2> "$tmp/x2.err"
check_record meth-named-converted "$tmp/x2.sam" 0
echo 'PASS: meth-named-converted'
# An index prefix spelled like an option value is not a stray value when its FASTA is gone.
mkdir "$tmp/named"
cp "$tmp/ref.fa" "$tmp/named/taps"
"$bin" index --meth "$tmp/named/taps" > "$tmp/index-taps.out" 2> "$tmp/index-taps.err"
rm "$tmp/named/taps"
(cd "$tmp/named" && "$bin" aln --meth taps "$tmp/r.fq" > "$tmp/named.sam" 2> "$tmp/named.err")
check_record index-named-taps "$tmp/named.sam" 0
echo 'PASS: index-named-taps'
# -f naming the original reference's header sidecar (.hdr, or the .dict it falls back to), or a
# link to it, is refused before the sidecar is read.
printf '@CO\tsidecar\n' > "$tmp/ref.fa.hdr"
cp "$tmp/ref.fa.hdr" "$tmp/hdr.saved"
check_error out-hdr 'is the header sidecar' aln --meth -f "$tmp/ref.fa.hdr" "$tmp/ref.fa" "$tmp/r.fq"
cmp "$tmp/ref.fa.hdr" "$tmp/hdr.saved" || {
    echo "FAIL: out-hdr: header sidecar changed" >&2
    exit 1
}
rm "$tmp/ref.fa.hdr"
printf '@CO\tdict\n' > "$tmp/ref.dict"
cp "$tmp/ref.dict" "$tmp/dict.saved"
ln -s ref.dict "$tmp/dict-link.sam"
check_error out-dict-link 'is the header sidecar' aln --meth -f "$tmp/dict-link.sam" "$tmp/ref.fa" "$tmp/r.fq"
cmp "$tmp/ref.dict" "$tmp/dict.saved" || {
    echo "FAIL: out-dict-link: header sidecar changed" >&2
    exit 1
}
rm "$tmp/ref.dict" "$tmp/dict-link.sam"
# A converted index staged with `bwa-mem3 shm` is attached even with its disk files gone (Linux
# only: macOS POSIX shm cannot hold an index). `shm -d` drops every staged index, so this runs only
# when the registry is readable and empty, and drops at the end only while its own index is the
# sole entry; a failed drop, or an entry left behind by it, fails the test.
# shellcheck source=test/framework/shm_registry_guard.sh
source "$(dirname "$0")/../framework/shm_registry_guard.sh"
if [ "$(uname -s)" != Linux ]; then
    echo 'SKIP: shm-staged-converted (macOS POSIX shm cannot hold an index)'
elif ! shm_registry_holds_only "$bin"; then
    echo 'SKIP: shm-staged-converted (another index is already staged, or shm -l failed; run bwa-mem3 shm -l to inspect)'
else
    shm_drop_own() {
        if ! shm_registry_holds_only "$bin" ref.fa.meth s.meth.meth; then
            echo 'FAIL: shm-staged-converted: the registry lists another index or shm -l failed; leaving shared memory as is' >&2
            return 1
        fi
        "$bin" shm -d > "$tmp/shm-d.out" 2> "$tmp/shm-d.err" || {
            echo 'FAIL: shm-staged-converted: shm -d failed' >&2
            cat "$tmp/shm-d.err" >&2
            return 1
        }
        shm_registry_holds_only "$bin" || {
            echo 'FAIL: shm-staged-converted: an index is still staged after shm -d' >&2
            return 1
        }
    }
    trap 'shm_drop_own || true; rm -rf "$tmp"' EXIT
    "$bin" shm --meth "$tmp/ref.fa" > "$tmp/shm.out" 2> "$tmp/shm.err"
    mkdir "$tmp/conv.saved"
    mv "$tmp"/ref.fa.meth.* "$tmp/conv.saved/"
    "$bin" aln --meth "$tmp/ref.fa" "$tmp/r.fq" > "$tmp/shm.sam" 2> "$tmp/shm.aln.err"
    grep -F 'attached from shm' "$tmp/shm.aln.err" > /dev/null || {
        echo "FAIL: shm-staged-converted: missing shared-memory attachment diagnostic" >&2
        exit 1
    }
    grep -v '^@' "$tmp/shm.sam" > "$tmp/shm.records" || {
        echo "FAIL: shm-staged-converted: no alignment records" >&2
        exit 1
    }
    cmp "$tmp/base.records" "$tmp/shm.records" || {
        echo "FAIL: shm-staged-converted: alignment records differ" >&2
        exit 1
    }
    mv "$tmp"/conv.saved/* "$tmp/"
    echo 'PASS: shm-staged-converted'
    # A reference whose own name ends in .meth, with its converted index staged and off disk, is
    # still the original prefix: the staged <prefix>.meth counts like one on disk.
    cp "$tmp/ref.fa" "$tmp/s.meth"
    "$bin" index --meth "$tmp/s.meth" > "$tmp/index-s.out" 2> "$tmp/index-s.err"
    # shm --meth keeps a prefix that already ends in .meth, so name the converted index itself.
    "$bin" shm --meth "$tmp/s.meth.meth" > "$tmp/shm-s.out" 2> "$tmp/shm-s.err"
    mkdir "$tmp/s.saved"
    mv "$tmp"/s.meth.meth.* "$tmp/s.saved/"
    "$bin" aln --meth "$tmp/s.meth" "$tmp/r.fq" > "$tmp/shm-s.sam" 2> "$tmp/shm-s.aln.err"
    grep -F 'attached from shm' "$tmp/shm-s.aln.err" > /dev/null || {
        echo "FAIL: shm-staged-meth-named-original: missing shared-memory attachment diagnostic" >&2
        exit 1
    }
    grep -v '^@' "$tmp/shm-s.sam" > "$tmp/shm-s.records" || {
        echo "FAIL: shm-staged-meth-named-original: no alignment records" >&2
        exit 1
    }
    cmp "$tmp/base.records" "$tmp/shm-s.records" || {
        echo "FAIL: shm-staged-meth-named-original: alignment records differ" >&2
        exit 1
    }
    trap 'rm -rf "$tmp"' EXIT
    shm_drop_own
    echo 'PASS: shm-staged-meth-named-original'
fi
