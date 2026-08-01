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
# sidecar already carrying AH). Every case that reads a sidecar runs against
# both names bwa_load_hdr_from_index accepts (<prefix>.hdr and
# <baseprefix>.dict), plus one case pinning that .hdr wins when both exist.
# The remaining gating inside bwa_warn_sidecar_missing_AH (--compat, the --bam
# path) is covered at `mem` level by test/regression/header_parity.sh, which
# uses a real `samtools dict --alt` sidecar rather than a hand-written one.
# This script deliberately needs no samtools, so it runs on every CI leg.
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
write_dict () {  # $1 = destination sidecar, $2 = extra tags on the ALT record
    {
        printf '@HD\tVN:1.5\tSO:unsorted\n'
        printf '@SQ\tSN:chr1\tLN:2000\tM5:aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\tAS:test\tUR:file:/refs/ref.fa\tSP:synthetic\n'
        printf '@SQ\tSN:chr2\tLN:1500\tM5:bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb\tAS:test\tUR:file:/refs/ref.fa\tSP:synthetic\n'
        printf '@SQ\tSN:%s\tLN:1000\tM5:cccccccccccccccccccccccccccccccc\tAS:test\tUR:file:/refs/ref.fa\tSP:synthetic%s\n' \
            "$ALT_SN" "${2-}"
    } > "$1"
}

# bwa_load_hdr_from_index (src/bwa.cpp) resolves <prefix>.hdr FIRST and only
# falls back to <baseprefix>.dict. Both names must therefore be exercised: a
# regression in either lookup leaves the other passing. For `index ref.fa` the
# prefix is `ref.fa`, so the two names are `ref.fa.hdr` and `ref.dict`.
SIDECAR_KINDS=(dict hdr)

# Assigns the path to the global SIDECAR rather than printing it: a `$(...)`
# substitution runs in a subshell, where the unknown-kind `exit` below would
# kill only that subshell and leave the script running with an empty path.
set_sidecar () {  # $1 = case dir, $2 = kind (dict|hdr)
    case "$2" in
        dict) SIDECAR="$1/ref.dict"   ;;   # <baseprefix>.dict, the fallback
        hdr)  SIDECAR="$1/ref.fa.hdr" ;;   # <prefix>.hdr, tried first
        *)    echo "internal error: unknown sidecar kind '$2'" >&2; exit 2 ;;
    esac
}

write_alt () {  # $1 = destination .alt
    printf '%s\t0\tchr1\t1\t60\t1000M\t*\t0\t0\t*\t*\n' "$ALT_SN" > "$1"
}

# --------------------------------------------------------------------------
# Case A: .alt present, sidecar present without AH -> MUST warn
# --------------------------------------------------------------------------
case_A () {  # $1 = sidecar kind
    local kind="$1"
    local a="$tmp/a-$kind"; mkdir -p "$a"
    make_ref   "$a/ref.fa"
    write_alt  "$a/ref.fa.alt"                 # .alt is <prefix>.alt
    set_sidecar "$a" "$kind"
    write_dict "$SIDECAR"                      # no AH on the ALT record
    "$bin" index "$a/ref.fa" > /dev/null 2> "$a/err.txt" \
        || fail "$a/err.txt" "case A[$kind]: index exited non-zero"

    grep -q "$WARN_RE" "$a/err.txt" \
        || fail "$a/err.txt" "case A[$kind]: expected the missing-AH warning at index time"
    grep -q "samtools dict" "$a/err.txt" \
        || fail "$a/err.txt" "case A[$kind]: warning must name the \`samtools dict --alt\` remedy"
    grep -q "$ALT_SN" "$a/err.txt" \
        || fail "$a/err.txt" "case A[$kind]: warning should name the offending contig ($ALT_SN)"
}

for kind in "${SIDECAR_KINDS[@]}"; do case_A "$kind"; done

# --------------------------------------------------------------------------
# Case B: .alt present, NO sidecar -> must NOT warn. bwa-mem3 generates its own
# @SQ with AH:*, so nothing is lost. Exercises the no-sidecar early return.
#
# Not parameterized over SIDECAR_KINDS: the case is defined by the absence of
# every sidecar name, so there is only one variant of it to run.
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
case_C () {  # $1 = sidecar kind
    local kind="$1"
    local c="$tmp/c-$kind"; mkdir -p "$c"
    make_ref   "$c/ref.fa"
    set_sidecar "$c" "$kind"
    write_dict "$SIDECAR"
    "$bin" index "$c/ref.fa" > /dev/null 2> "$c/err.txt" \
        || fail "$c/err.txt" "case C[$kind]: index exited non-zero"

    ! grep -q "$WARN_RE" "$c/err.txt" \
        || fail "$c/err.txt" "case C[$kind]: reference has no ALT contigs; must not warn"
}

for kind in "${SIDECAR_KINDS[@]}"; do case_C "$kind"; done

# --------------------------------------------------------------------------
# Case D: .alt present AND the sidecar's ALT @SQ already carries AH:* -> must
# NOT warn. This is the state case A tells the user to reach, so it is the one
# case where a false positive would train them to ignore the warning. It is
# also the only negative case that exercises the AH-detection itself rather
# than an early return: cases B and C bail before any @SQ is scanned.
# --------------------------------------------------------------------------
case_D () {  # $1 = sidecar kind
    local kind="$1"
    local d="$tmp/d-$kind"; mkdir -p "$d"
    set_sidecar "$d" "$kind"
    make_ref   "$d/ref.fa"
    write_alt  "$d/ref.fa.alt"
    # $'...' rather than '\tAH:*': write_dict passes the extra tags through a
    # printf %s, which does NOT expand escapes, so a plain '\tAH:*' would leave a
    # literal backslash-t in the record. The AH token would not be tab-delimited,
    # the sidecar would still be AH-less as far as the scanner is concerned, and the
    # case would "fail" on a fixture bug rather than on the code.
    write_dict "$SIDECAR" $'\tAH:*'   # what `samtools dict --alt` produces
    # Verify the fixture is what the case claims: AH:* must be a whole tab-delimited
    # field on the ALT record, which is the only form bwa_warn_sidecar_missing_AH
    # recognizes.
    awk -F'\t' -v sn="SN:$ALT_SN" '
        $1 == "@SQ" {
            for (i = 2; i <= NF; ++i) if ($i == sn) alt = 1
            if (alt) { for (i = 2; i <= NF; ++i) if ($i == "AH:*") ok = 1; alt = 0 }
        }
        END { exit !ok }' "$SIDECAR" \
        || fail "$SIDECAR" "case D[$kind]: fixture bug -- AH:* is not a tab-delimited field on the ALT @SQ"
    "$bin" index "$d/ref.fa" > /dev/null 2> "$d/err.txt" \
        || fail "$d/err.txt" "case D[$kind]: index exited non-zero"

    ! grep -q "$WARN_RE" "$d/err.txt" \
        || fail "$d/err.txt" "case D[$kind]: sidecar already carries AH:* on the ALT contig; must not warn"
}

for kind in "${SIDECAR_KINDS[@]}"; do case_D "$kind"; done

# --------------------------------------------------------------------------
# Case E: .alt present, an AH-less <prefix>.hdr AND an AH-carrying
# <baseprefix>.dict -> MUST warn. Pins the resolution ORDER rather than either
# lookup in isolation: the cases above pass whichever name wins, so only a
# disagreeing pair proves .hdr is consulted first. If the precedence ever
# flipped to .dict, the AH-carrying .dict would silence the warning here.
# --------------------------------------------------------------------------
e="$tmp/e"; mkdir -p "$e"
make_ref   "$e/ref.fa"
write_alt  "$e/ref.fa.alt"
write_dict "$e/ref.fa.hdr"            # no AH -- the winner, so this must warn
write_dict "$e/ref.dict" $'\tAH:*'    # has AH -- must be ignored
"$bin" index "$e/ref.fa" > /dev/null 2> "$e/err.txt" \
    || fail "$e/err.txt" "case E: index exited non-zero"

grep -q "$WARN_RE" "$e/err.txt" \
    || fail "$e/err.txt" "case E: <prefix>.hdr outranks <baseprefix>.dict; the AH-less .hdr must warn"

# --------------------------------------------------------------------------
# Case F: an AH-less sidecar whose ALT contig has a LONG name (>= 256 bytes)
# -> MUST still warn. bwa_warn_sidecar_missing_AH has to copy each sidecar SN
# out of the header text to get a NUL-terminated hash key; when that copy went
# into a fixed 256-byte stack buffer, any longer SN was skipped outright. That
# silently disabled this very check for exactly the contigs it is about --
# contig names are strdup()ed from the FASTA name with no length bound
# (bntseq.cpp), so a long-named ALT contig is representable in bns.
#
# `local ALT_SN` rather than a subshell: make_ref / write_alt / write_dict read
# $ALT_SN at call time, so bash's dynamic scoping makes them all use the long
# name here, while `fail`'s exit still terminates the script (it would only
# leave a subshell).
# --------------------------------------------------------------------------
case_F () {
    local ALT_SN
    ALT_SN="chrX_$(printf 'L%.0s' $(seq 1 300))_alt"   # ~306 bytes
    [[ ${#ALT_SN} -ge 256 ]] \
        || { echo "FAIL [index alt sidecar warn]: case F fixture bug -- ALT_SN is only ${#ALT_SN} bytes, need >= 256" >&2; exit 1; }
    local f="$tmp/f"; mkdir -p "$f"
    make_ref   "$f/ref.fa"
    write_alt  "$f/ref.fa.alt"
    write_dict "$f/ref.dict"          # no AH on the long-named ALT record
    "$bin" index "$f/ref.fa" > /dev/null 2> "$f/err.txt" \
        || fail "$f/err.txt" "case F: index exited non-zero"

    grep -q "$WARN_RE" "$f/err.txt" \
        || fail "$f/err.txt" "case F: a ${#ALT_SN}-byte ALT contig name must not be skipped by the ALT/AH check"
    grep -q "$ALT_SN" "$f/err.txt" \
        || fail "$f/err.txt" "case F: warning should name the offending contig"
}

case_F

echo "PASS: index reports an AH-less sidecar only when the reference has ALT contigs and the sidecar omits AH (both .dict and .hdr, .hdr first; long contig names included)"
