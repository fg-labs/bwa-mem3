#!/usr/bin/env bash
# test/index_alt_sidecar_warn_test.sh
#
# Asserts that `bwa-mem3 index` reports the ALT/AH gap: the reference has ALT
# contigs (via <prefix>.alt) but the header sidecar (<baseprefix>.dict /
# <prefix>.hdr) supplies @SQ records with no AH tag.
#
# Why this matters: the sidecar's @SQ is authoritative by design (a port of
# lh3/bwa#348), so bwa-mem3 emits it verbatim and does NOT inject AH:*. A
# Picard/GATK-style .dict carries no AH -- Picard has no notion of a .alt file
# -- so the standard Broad reference layout silently loses ALT status from every
# output header. `mem` reports this too; indexing is the earlier place to hear
# it, when the remedy (`samtools dict --alt`) still costs nothing.
#
# Scope: this covers the `index` call site -- that it is reached, and that it
# stays quiet in each of the three ways it should (no sidecar, no ALT contig,
# sidecar already carrying AH). The remaining gating inside
# bwa_warn_sidecar_missing_AH (--compat, the --bam path) is covered at `mem`
# level by test/regression/header_parity.sh, which uses a real
# `samtools dict --alt` sidecar rather than a hand-written one. This script
# deliberately needs no samtools, so it runs on every CI leg.
#
# Usage: test/index_alt_sidecar_warn_test.sh <bwa-mem3-binary>

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 <bwa-mem3-binary>" >&2
    exit 2
fi

bin="$1"
[[ -x "$bin" ]] || { echo "FAIL: bwa-mem3 binary not executable at $bin" >&2; exit 1; }

fail () {  # $1 = stderr file to dump, $2... = message
    local err="$1"; shift
    echo "FAIL [index alt sidecar warn]: $*" >&2
    echo "--- index stderr:" >&2
    cat "$err" >&2
    exit 1
}

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

ALT_SN="chrX_KI000001v1_alt"
# The warning is emitted by bwa_warn_sidecar_missing_AH (src/bwa.cpp); match its
# text, not a bare "ALT", so the assertion cannot pass on unrelated output.
WARN_RE="sidecar supplies @SQ without an AH tag"

# One line per sequence: `fold` drops the final newline when the length is not a
# multiple of the width, gluing the next '>' onto the sequence -- silently, and
# the index then builds with a single contig.
emit_contig () {  # $1 = name, $2 = length
    printf '>%s\n' "$1"
    head -c "$2" /dev/zero | tr '\0' 'A'
    printf '\n'
}

make_ref () {  # $1 = destination fasta
    {
        emit_contig chr1      2000
        emit_contig chr2      1500
        emit_contig "$ALT_SN" 1000
    } > "$1"
}

# Picard/GATK-style sequence dictionary (SN/LN/M5/AS/UR/SP) -- what
# CreateSequenceDictionary emits and what the Broad bundle ships as
# <basename>.dict. $2 appends tags to the ALT record: '' for Picard (no AH).
write_dict () {  # $1 = destination .dict, $2 = extra tags on the ALT record
    {
        printf '@HD\tVN:1.5\tSO:unsorted\n'
        printf '@SQ\tSN:chr1\tLN:2000\tM5:aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\tAS:test\tUR:file:/refs/ref.fa\tSP:synthetic\n'
        printf '@SQ\tSN:chr2\tLN:1500\tM5:bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb\tAS:test\tUR:file:/refs/ref.fa\tSP:synthetic\n'
        printf '@SQ\tSN:%s\tLN:1000\tM5:cccccccccccccccccccccccccccccccc\tAS:test\tUR:file:/refs/ref.fa\tSP:synthetic%s\n' \
            "$ALT_SN" "${2-}"
    } > "$1"
}

write_alt () {  # $1 = destination .alt
    printf '%s\t0\tchr1\t1\t60\t1000M\t*\t0\t0\t*\t*\n' "$ALT_SN" > "$1"
}

# --------------------------------------------------------------------------
# Case A: .alt present, sidecar present without AH -> MUST warn
# --------------------------------------------------------------------------
a="$tmp/a"; mkdir -p "$a"
make_ref   "$a/ref.fa"
write_alt  "$a/ref.fa.alt"   # .alt is <prefix>.alt
write_dict "$a/ref.dict"     # .dict is <baseprefix>.dict
"$bin" index "$a/ref.fa" > /dev/null 2> "$a/err.txt" \
    || fail "$a/err.txt" "case A: index exited non-zero"

grep -q "$WARN_RE" "$a/err.txt" \
    || fail "$a/err.txt" "case A: expected the missing-AH warning at index time"
grep -q "samtools dict" "$a/err.txt" \
    || fail "$a/err.txt" "case A: warning must name the \`samtools dict --alt\` remedy"
grep -q "$ALT_SN" "$a/err.txt" \
    || fail "$a/err.txt" "case A: warning should name the offending contig ($ALT_SN)"

# --------------------------------------------------------------------------
# Case B: .alt present, NO sidecar -> must NOT warn. bwa-mem3 generates its own
# @SQ with AH:*, so nothing is lost. Exercises the no-sidecar early return.
# --------------------------------------------------------------------------
b="$tmp/b"; mkdir -p "$b"
make_ref  "$b/ref.fa"
write_alt "$b/ref.fa.alt"
"$bin" index "$b/ref.fa" > /dev/null 2> "$b/err.txt" \
    || fail "$b/err.txt" "case B: index exited non-zero"

! grep -q "$WARN_RE" "$b/err.txt" \
    || fail "$b/err.txt" "case B: no sidecar present; must not warn"

# --------------------------------------------------------------------------
# Case C: sidecar without AH but NO .alt -> must NOT warn. No contig is ALT, so
# a missing AH is correct rather than a loss. This is the common case for most
# references, where a spurious warning would be pure noise.
# --------------------------------------------------------------------------
c="$tmp/c"; mkdir -p "$c"
make_ref   "$c/ref.fa"
write_dict "$c/ref.dict"
"$bin" index "$c/ref.fa" > /dev/null 2> "$c/err.txt" \
    || fail "$c/err.txt" "case C: index exited non-zero"

! grep -q "$WARN_RE" "$c/err.txt" \
    || fail "$c/err.txt" "case C: reference has no ALT contigs; must not warn"

# --------------------------------------------------------------------------
# Case D: .alt present AND the sidecar's ALT @SQ already carries AH:* -> must
# NOT warn. This is the state case A tells the user to reach, so it is the one
# case where a false positive would train them to ignore the warning. It is
# also the only negative case that exercises the AH-detection itself rather
# than an early return: cases B and C bail before any @SQ is scanned.
# --------------------------------------------------------------------------
d="$tmp/d"; mkdir -p "$d"
make_ref   "$d/ref.fa"
write_alt  "$d/ref.fa.alt"
# $'...' rather than '\tAH:*': write_dict passes the extra tags through a
# printf %s, which does NOT expand escapes, so a plain '\tAH:*' would leave a
# literal backslash-t in the record. The AH token would not be tab-delimited,
# the sidecar would still be AH-less as far as the scanner is concerned, and the
# case would "fail" on a fixture bug rather than on the code.
write_dict "$d/ref.dict" $'\tAH:*'   # what `samtools dict --alt` produces
# Verify the fixture is what the case claims: AH:* must be a whole tab-delimited
# field on the ALT record, which is the only form bwa_warn_sidecar_missing_AH
# recognizes.
awk -F'\t' -v sn="SN:$ALT_SN" '
    $1 == "@SQ" {
        for (i = 2; i <= NF; ++i) if ($i == sn) alt = 1
        if (alt) { for (i = 2; i <= NF; ++i) if ($i == "AH:*") ok = 1; alt = 0 }
    }
    END { exit !ok }' "$d/ref.dict" \
    || fail "$d/ref.dict" "case D: fixture bug -- AH:* is not a tab-delimited field on the ALT @SQ"
"$bin" index "$d/ref.fa" > /dev/null 2> "$d/err.txt" \
    || fail "$d/err.txt" "case D: index exited non-zero"

! grep -q "$WARN_RE" "$d/err.txt" \
    || fail "$d/err.txt" "case D: sidecar already carries AH:* on the ALT contig; must not warn"

echo "PASS: index reports an AH-less sidecar only when the reference has ALT contigs and the sidecar omits AH"
