#!/usr/bin/env bash
# test/regression/meth_tags.sh
#
# Regression (issue #331): --meth-tags selects which Bismark tags reach the BAM.
#
#  1. The default is the full set -- a bare `--meth` emits XR:Z, XG:Z and XM:Z,
#     so existing methylation pipelines are unaffected by the flag's addition.
#  2. Every spec form resolves to the right set: `all`, `none`, an inclusion
#     list (`XR,XG`), and `^`-prefixed exclusions (`^XM`, the motivating case).
#  3. Deselecting a tag removes ONLY that tag. The alignment itself is
#     untouched: fields 1-11 plus NM/MD/AS are byte-identical between a default
#     run and a `--meth-tags none` run, so the flag is an emission filter and
#     not an alignment-affecting knob.
#  4. Malformed specs are rejected with a non-zero exit rather than silently
#     falling back to a default -- a typo'd `--meth-tags ^XN` must not quietly
#     emit everything.
#
# Inputs:
#   BWA_MEM3 — path to the bwa-mem3 binary under test
set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"
command -v samtools >/dev/null 2>&1 || { echo "SKIP: samtools not on PATH (--meth emits BAM)"; exit 0; }

WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT; cd "$WORK"
fail() { echo "FAIL: $*" >&2; exit 1; }

# Deterministic 1500 bp reference (same generator as meth_output_integrity).
REF=TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTAGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGAATCGTCCCTTCTAACGTAGTTCGCATAGTTCCCGTCCGTAGCCGGACTATTCGAACACCCAGTATTCGATTAACTCGGGCTTGACGTATTAGAGGCGTTAGTGTGCCAGGTAAGATACGCCAACGGAATTAACCTCTGTGACACTCCGCGGAGCCTTCGGACATATAAGTGATCGGGTCTACGTTTGTTAGACTTGAGACGTCTGTTAAGAGTTGGGTCTAATAAATCGCCTACACGTGGAGTCTAACGGGGAAGCGTCGAATCCTGATACATCATATAATGGAGCGTGTTATGAAAAAAGAGCATTCCATTGTACGAGCCGTGCCAGAAACGGCTTGACTACGTGAGCGTAGTGTTAGATAAACAGGAAACACTGACGCGGTTAGAAGGCGGATTGCCGGTAGGTTTTGGAAACATAAATACACACGGTATCATGTTGGGTCACGATTCCTATCACCGCACAGGGCCAACCATAGAAGAACTGAAAGAACTAATCTGGCGGCGGGCTCGGTGCTTATATTTTCCACCCAACATCGTGCACATTAGGCTCACCGCGCCCTACGGGCGAAGGGTGCGTACGGTGTTTATAAGGCGTGACGGCCCCAAGTAGAGGGTAATTCTGTGAAAGAATCTCAGGACGGTGGCATGAATTCAATTCCTTTTAAACCTATCGTTCCGACCTTATGCAATCCTTCAATGAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTACTCTTCACTGCCCAAGCGAGGCTAACCTTTTGAAATGTCACAGTCGAAGCATATCTCCCGTACATCTTTTTCGGAGATCGCAGCTCGCGGAGCTATAAGCGACTTAAGCCCTTGTGTCGGTGATCCCAAGGGTCTGACTCCTGTACCAGGGTTACTGTTTCGCTTTACGGAGTAGCCTGTGAGGTGAACTGAAAGGAGCATATTTGAGATCTAAGATAGGGTCCTCCTCTGCGTCTACGTTCTCTCCGTTACGTACGGCTTCGCACCGGAGTGCATCTTGGCCCCGAAACGCACTGTGTGTGCTGATACAGCGTCCCTGGCCGGCCATGGGTTCAGAACTCCCGGGAACGCTTTTCAACTTAGAGGAACCCCGTCATGGAAGTAGATCGCGTCGAATGAGGGAGTTAGTCCTCGTTCCAGCTGGTAATTGTTTTACCGCTTGGGACCACTATAGGCCGCGGGTAGAAGTTGCTGGGTGTTGATTCCAACCCTCGAACCACGATACGACCTGCCATTTATGGCACAGTAAGGTTCAAACAGCATAATGAATACAGTTATAGTAACTTCCTCACGTACGATTAGGACGCAGCCTTG

printf '>chrA\n%s\n' "$REF" > ref.fa
"$BWA_MEM3" index --meth ref.fa >/dev/null 2>&1 || fail "index --meth nonzero exit"

# OT read over chrA[100:160] with 10 C->T conversions -- carries all three tags.
READ=ATTCTGGCGATTTAGGTGACTCACAGAATCGTTCTTTCTAACGTAGTTTGTATAGTTCTC
Q=$(printf 'I%.0s' $(seq 1 60))
printf '@t\n%s\n+\n%s\n' "$READ" "$Q" > r.fq

# Echo the meth tags present on the first record, space-separated and sorted,
# so the assertion does not depend on aux emission order.
tags_of() { # $1 = bam -> e.g. "XG XM XR"
    samtools view "$1" | mawk 'NR==1{for(i=12;i<=NF;i++) if($i~/^X[RGM]:Z:/) print substr($i,1,2)}' \
        | sort | tr '\n' ' ' | sed 's/ $//'
}
run() { # $1 = out.bam, rest = extra args
    local out="$1"; shift
    "$BWA_MEM3" mem --meth "$@" -t 1 ref.fa r.fq > "$out" 2>/dev/null \
        || fail "mem --meth $* nonzero exit"
    samtools quickcheck "$out" || fail "mem --meth $* produced an invalid BAM"
}

# --- 1. default is the full Bismark set ------------------------------------
run default.bam
got=$(tags_of default.bam)
[ "$got" = "XG XM XR" ] || fail "default: tags '$got', want 'XG XM XR' (a bare --meth must emit all three)"

# --- 2. every spec form resolves correctly ---------------------------------
check_spec() { # $1 = spec, $2 = expected sorted tag list
    run "spec.bam" --meth-tags "$1"
    local g; g=$(tags_of spec.bam)
    [ "$g" = "$2" ] || fail "--meth-tags $1: tags '$g', want '$2'"
}
check_spec all       "XG XM XR"
check_spec none      ""
check_spec "^XM"     "XG XR"          # the motivating case from issue #331
check_spec "XR,XG"   "XG XR"          # equivalent inclusion spelling
check_spec "XM"      "XM"
check_spec "^XR,^XG" "XM"
check_spec "xr,xg"   "XG XR"          # case-insensitive

# --- 3. tag selection does not perturb the alignment ------------------------
# Fields 1-11 plus the alignment-bearing tags must match the default run; only
# the meth tags may differ.
run none.bam --meth-tags none
core() { samtools view "$1" | mawk '{o=$1;for(i=2;i<=11;i++)o=o"\t"$i;
    for(i=12;i<=NF;i++) if($i~/^(NM:i:|MD:Z:|AS:i:)/) o=o"\t"$i; print o}'; }
core default.bam > core.default
core none.bam    > core.none
cmp -s core.default core.none \
    || fail "--meth-tags none changed the alignment: fields 1-11/NM/MD/AS differ from the default run"
[ -s core.default ] || fail "core comparison was vacuous (no records extracted)"

# --- 4. malformed specs are rejected, not silently defaulted ----------------
for bad in "XR,^XM" "XZ" "XR," ",XR" "XR,,XG" "" "^" "XRX"; do
    if "$BWA_MEM3" mem --meth --meth-tags "$bad" -t 1 ref.fa r.fq >/dev/null 2>/dev/null; then
        fail "--meth-tags '$bad' exited 0; a malformed spec must be rejected"
    fi
done

echo "PASS: meth_tags (default emits XR/XG/XM; all|none|list|^exclusion resolve correctly; selection does not perturb the alignment; malformed specs rejected)"
