#!/usr/bin/env bash
# test/regression/meth_flag_footguns.sh
#
# Regression: a stray word from a mis-spelled --meth-family flag must be
# diagnosed, not silently consumed as the reference.
#
# Several spellings orphan their value on the command line:
#   --meth taps                -- OPTIONAL argument; getopt_long only binds it with '='
#   --meth-seed-prune baseline -- OPTIONAL argument; same, silently selects spec30
#   --meth-tags XR XG          -- REQUIRED argument; 'XR' binds, 'XG' does not
#
# The orphan lands in the <idxbase> slot. With a single-end read file the
# invocation still has three positionals, which is a well-formed paired-end
# call whose reference happens to be named "taps" or "XG" -- so the arity check
# passes and the user gets a missing-index error naming a token they never
# meant as a path (previously: "--meth seed index 'XG.meth.*' not found. Run:
# bwa-mem3 index --meth XG").
#
# Asserted here:
#  1. Both spellings are rejected with a message that names the offending token
#     AND the flag it belongs to, so the fix is readable off the error.
#  2. The check is gated on the token not existing: a reference genuinely named
#     `taps` (a value in the vocabulary) still aligns normally. A guard that
#     blocks a legitimate filename would be worse than the footgun.
#  3. `-XM` is accepted as a synonym for `^XM`, so the exclusion form needs no
#     shell quoting (a bare ^XM is a negated glob under zsh EXTENDED_GLOB).
#
# Inputs:
#   BWA_MEM3 — path to the bwa-mem3 binary under test
set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"
command -v samtools > /dev/null 2>&1 || {
    echo "SKIP: samtools not on PATH (--meth emits BAM)"
    exit 0
}
command -v mawk > /dev/null 2>&1 || {
    echo "SKIP: mawk not on PATH (required to inspect BAM tags)"
    exit 0
}

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
fail() {
    echo "FAIL: $*" >&2
    exit 1
}

REF=TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTAGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGAATCGTCCCTTCTAACGTAGTTCGCATAGTTCCCGTCCGTAGCCGGACTATTCGAACACCCAGTATTCGATTAACTCGGGCTTGACGTATTAGAGGCGTTAGTGTGCCAGGTAAGATACGCCAACGGAATTAACCTCTGTGACACTCCGCGGAGCCTTCGGACATATAAGTGATCGGGTCTACGTTTGTTAGACTTGAGACGTCTGTTAAGAGTTGGGTCTAATAAATCGCCTACACGTGGAGTCTAACGGGGAAGCGTCGAATCCTGATACATCATATAATGGAGCGTGTTATGAAAAAAGAGCATTCCATTGTACGAGCCGTGCCAGAAACGGCTTGACTACGTGAGCGTAGTGTTAGATAAACAGGAAACACTGACGCGGTTAGAAGGCGGATTGCCGGTAGGTTTTGGAAACATAAATACACACGGTATCATGTTGGGTCACGATTCCTATCACCGCACAGGGCCAACCATAGAAGAACTGAAAGAACTAATCTGGCGGCGGGCTCGGTGCTTATATTTTCCACCCAACATCGTGCACATTAGGCTCACCGCGCCCTACGGGCGAAGGGTGCGTACGGTGTTTATAAGGCGTGACGGCCCCAAGTAGAGGGTAATTCTGTGAAAGAATCTCAGGACGGTGGCATGAATTCAATTCCTTTTAAACCTATCGTTCCGACCTTATGCAATCCTTCAATGAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTACTCTTCACTGCCCAAGCGAGGCTAACCTTTTGAAATGTCACAGTCGAAGCATATCTCCCGTACATCTTTTTCGGAGATCGCAGCTCGCGGAGCTATAAGCGACTTAAGCCCTTGTGTCGGTGATCCCAAGGGTCTGACTCCTGTACCAGGGTTACTGTTTCGCTTTACGGAGTAGCCTGTGAGGTGAACTGAAAGGAGCATATTTGAGATCTAAGATAGGGTCCTCCTCTGCGTCTACGTTCTCTCCGTTACGTACGGCTTCGCACCGGAGTGCATCTTGGCCCCGAAACGCACTGTGTGTGCTGATACAGCGTCCCTGGCCGGCCATGGGTTCAGAACTCCCGGGAACGCTTTTCAACTTAGAGGAACCCCGTCATGGAAGTAGATCGCGTCGAATGAGGGAGTTAGTCCTCGTTCCAGCTGGTAATTGTTTTACCGCTTGGGACCACTATAGGCCGCGGGTAGAAGTTGCTGGGTGTTGATTCCAACCCTCGAACCACGATACGACCTGCCATTTATGGCACAGTAAGGTTCAAACAGCATAATGAATACAGTTATAGTAACTTCCTCACGTACGATTAGGACGCAGCCTTG
printf '>chrA\n%s\n' "$REF" > ref.fa
READ=ATTCTGGCGATTTAGGTGACTCACAGAATCGTTCTTTCTAACGTAGTTTGTATAGTTCTC
printf '@t\n%s\n+\n%s\n' "$READ" "$(printf 'I%.0s' $(seq 1 60))" > r.fq
"$BWA_MEM3" index --meth ref.fa > /dev/null 2>&1 || fail "index --meth nonzero exit"

# --- 1. both orphan spellings are diagnosed --------------------------------
# Single-end read file, so the stray token makes exactly three positionals and
# the arity check cannot catch it -- this is the case that used to slip through.
expect_diagnosed() { # $1 = label, $2 = stray token, $3 = flag it belongs to, rest = argv
    local label="$1" token="$2" flag="$3"
    shift 3
    if "$BWA_MEM3" "$@" > /dev/null 2> err.txt; then
        fail "$label: exited 0; the stray '$token' must be rejected"
    fi
    grep -q "'$token' was taken as a positional argument" err.txt \
        || fail "$label: error does not name the stray token '$token'; got: $(head -2 err.txt | tr '\n' ' ')"
    # Match the whole clause, not just the flag: a bare "--meth" would also be
    # satisfied by "--meth-tags", so a value attributed to the wrong member of
    # the family would slip through.
    grep -Fq -- "the value for $flag and no such file exists" err.txt \
        || fail "$label: error does not attribute '$token' to $flag; got: $(head -2 err.txt | tr '\n' ' ')"
    # The old, misleading advice must be gone: it told the user to build an
    # index for the stray token.
    if grep -q "index --meth $token" err.txt; then
        fail "$label: error still advises building an index for the stray token '$token'"
    fi
}
expect_diagnosed "--meth-tags space-separated" XG --meth-tags mem --meth --meth-tags XR XG ref.fa r.fq
expect_diagnosed "--meth separated argument" taps --meth mem --meth taps ref.fa r.fq
expect_diagnosed "--meth-scoring orphan" genomic --meth-scoring mem --meth --meth-scoring collapsed genomic ref.fa r.fq
# --meth-seed-prune is also OPTIONAL-argument, so `--meth-seed-prune baseline`
# (space, no '=') silently selects spec30 and orphans `baseline` into <idxbase>.
# Diagnose it by name rather than letting it slide into a missing-index error.
expect_diagnosed "--meth-seed-prune orphan" baseline --meth-seed-prune mem --meth --meth-seed-prune baseline ref.fa r.fq

# --- 2. no false positive on a real file with a vocabulary name -------------
# A reference genuinely named `taps` must still work: the guard fires only when
# the token names nothing on disk.
cp ref.fa taps
"$BWA_MEM3" index --meth taps > /dev/null 2>&1 || fail "index --meth taps nonzero exit"
"$BWA_MEM3" mem --meth taps r.fq > real.bam 2> /dev/null \
    || fail "a reference genuinely named 'taps' was blocked by the stray-value guard"
[ "$(samtools view -c real.bam)" = "1" ] || fail "reference named 'taps': expected 1 record"
# Same file, this time with the chemistry actually selected.
"$BWA_MEM3" mem --meth=taps taps r.fq > real2.bam 2> /dev/null \
    || fail "--meth=taps with a reference named 'taps' was rejected"
[ "$(samtools view -c real2.bam)" = "1" ] || fail "--meth=taps + ref 'taps': expected 1 record"

# --- 3. '-XM' is a shell-safe synonym for '^XM' -----------------------------
tags_of() { samtools view "$1" | mawk 'NR==1{for(i=12;i<=NF;i++) if($i~/^X[RGM]:Z:/) print substr($i,1,2)}' \
    | sort | tr '\n' ' ' | sed 's/ $//'; }
"$BWA_MEM3" mem --meth --meth-tags -XM -t 1 ref.fa r.fq > dash.bam 2> /dev/null \
    || fail "--meth-tags -XM nonzero exit"
"$BWA_MEM3" mem --meth --meth-tags '^XM' -t 1 ref.fa r.fq > caret.bam 2> /dev/null \
    || fail "--meth-tags ^XM nonzero exit"
d=$(tags_of dash.bam)
c=$(tags_of caret.bam)
[ "$d" = "XG XR" ] || fail "--meth-tags -XM: tags '$d', want 'XG XR'"
[ "$d" = "$c" ] || fail "-XM and ^XM disagree: '$d' vs '$c'"

echo "PASS: meth_flag_footguns (orphaned --meth/--meth-tags/--meth-scoring/--meth-seed-prune values diagnosed by name; real files with vocabulary names unaffected; -XM == ^XM)"
