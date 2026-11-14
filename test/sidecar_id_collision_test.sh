#!/usr/bin/env bash
# test/sidecar_id_collision_test.sh
#
# Asserts that an @RG or @PG supplied by the user (-R / -H) REPLACES a
# <prefix>.hdr / <baseprefix>.dict sidecar record of the same type and ID,
# rather than being emitted next to it.
#
# Bug: with a sidecar carrying `@RG ID:g SM:old` and `-R '@RG\tID:g\tSM:new'`,
# `mem` wrote both `@RG ID:g` lines. The SAM spec requires @RG (and @PG) IDs to
# be unique, and the documented header precedence is user > index sidecar >
# default. The sidecar's records with other IDs must still be forwarded, in
# their original order.
#
# Checked on every writer: SAM text, --bam, --meth (SAM text) and
# --meth --bam. The BAM outputs are read back with samtools; without it, those
# cases are skipped with a warning (CI has samtools, so they stay enforced).
# The SAM writer is also checked with the sidecar as <baseprefix>.dict.
#
# The test copies phix.fa to a temp dir and indexes it there, so the shared
# fixtures (and their sidecars) are untouched.
#
# Usage: test/sidecar_id_collision_test.sh <bwa-mem3-binary> <fixtures-dir>

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 <bwa-mem3-binary> <fixtures-dir>" >&2
    exit 2
fi

bin="$1"
fixtures="$2"
src_ref="$fixtures/phix.fa"
reads="$fixtures/reads.fa"

[[ -x "$bin" ]] || {
    echo "FAIL: bwa-mem3 binary not executable at $bin" >&2
    exit 1
}
[[ -s "$src_ref" ]] || {
    echo "FAIL: phix.fa missing at $src_ref" >&2
    exit 1
}
[[ -s "$reads" ]] || {
    echo "FAIL: reads.fa missing at $reads" >&2
    exit 1
}

fail() {
    echo "FAIL [sidecar id collision]: $*" >&2
    exit 1
}

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

ref="$tmp/phix.fa"
cp "$src_ref" "$ref"
# `index --meth` also builds the plain index the non---meth cases align against.
"$bin" index --meth "$ref" > /dev/null 2>&1 || fail "bwa-mem3 index --meth on phix.fa failed"

# The sidecar: a colliding @RG and @PG between two records that must survive.
printf '%s\n' \
    $'@RG\tID:keep1\tSM:k1' \
    $'@RG\tID:g\tSM:old' \
    $'@PG\tID:p\tPN:old' \
    $'@RG\tID:keep2\tSM:k2' \
    > "$ref.hdr"

RG=$'@RG\tID:g\tSM:new'
USER_PG=$'@PG\tID:p\tPN:new'

# check_header <label> <header-file>
check_header() {
    local label="$1" hdr="$2" n
    n="$(grep -cE $'^@RG\t(.*\t)?ID:g(\t|$)' "$hdr" || true)"
    [[ "$n" -eq 1 ]] || {
        sed $'s/\t/<TAB>/g' "$hdr" >&2
        fail "[$label] expected exactly 1 '@RG ID:g', got $n"
    }
    grep -qE $'^@RG\tID:g\tSM:new$' "$hdr" || fail "[$label] the user's @RG ID:g (SM:new) is missing"
    ! grep -q 'SM:old' "$hdr" || fail "[$label] the sidecar's colliding @RG (SM:old) was emitted"

    n="$(grep -cE $'^@PG\t(.*\t)?ID:p(\t|$)' "$hdr" || true)"
    [[ "$n" -eq 1 ]] || {
        sed $'s/\t/<TAB>/g' "$hdr" >&2
        fail "[$label] expected exactly 1 '@PG ID:p', got $n"
    }
    grep -qE $'^@PG\tID:p\tPN:new$' "$hdr" || fail "[$label] the user's @PG ID:p (PN:new) is missing"
    ! grep -q 'PN:old' "$hdr" || fail "[$label] the sidecar's colliding @PG (PN:old) was emitted"

    # The non-colliding sidecar @RGs survive, in order, before the user's.
    local order
    order="$(grep -oE $'^@RG\tID:[^\t]*' "$hdr" | cut -f2 | tr '\n' ' ' || true)"
    [[ "$order" == "ID:keep1 ID:keep2 ID:g " ]] \
        || fail "[$label] @RG order is '$order', expected 'ID:keep1 ID:keep2 ID:g '"
    echo "OK:   [$label] user @RG/@PG replaced the sidecar's; other sidecar @RGs kept in order"
}

# run_case <label> <is-bam> <mem args...>
run_case() {
    local label="$1" is_bam="$2"
    shift 2
    "$bin" mem "$@" -R "$RG" -H "$USER_PG" "$ref" "$reads" > "$tmp/out" 2> "$tmp/err" || {
        cat "$tmp/err" >&2
        fail "[$label] mem exited non-zero"
    }
    if [[ "$is_bam" -eq 1 ]]; then
        samtools view -H "$tmp/out" > "$tmp/hdr" 2> /dev/null \
            || fail "[$label] samtools could not read the BAM header"
    else
        grep '^@' "$tmp/out" > "$tmp/hdr" || true
    fi
    check_header "$label" "$tmp/hdr"
}

run_case "SAM" 0
run_case "--meth" 0 --meth
if command -v samtools > /dev/null 2>&1; then
    run_case "--bam" 1 --bam
    run_case "--meth --bam" 1 --meth --bam
    bam_note=""
else
    echo "SKIP: samtools not found; --bam and --meth --bam not checked" >&2
    bam_note=" (SKIPPED: --bam, --meth --bam; samtools not found)"
fi

# The same sidecar as <baseprefix>.dict, the loader's fallback when no .hdr
# exists, on the SAM writer.
mv "$ref.hdr" "$tmp/phix.dict"
run_case "SAM, .dict sidecar" 0

echo "PASS: user @RG/@PG replace same-ID sidecar records on each writer checked$bam_note"
