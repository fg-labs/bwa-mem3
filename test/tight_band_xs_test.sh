#!/usr/bin/env bash
# test/tight_band_xs_test.sh
#
# Regression for the retry ladder's tight_band early-accept. A right extension
# with a tight_band proof of 91 used to stop at w = 100 while the full-width
# ladder (max_off 84 >= 3w/4) goes on to w = 200. Score and extent agree, but
# the region records a->w = 100 instead of 200, which narrows the contained-seed
# purge window: a seed the full ladder purges is extended instead, and its
# score-91 region became the read's XS. The exact ladders now make their own
# stop decision, so this read reports XS:i:37, as bwa does.
#
# Asserts, on one synthetic 306 bp read (fixture from #529), the record bwa
# 0.7.19 emits for it:
#   * FLAG 16, POS 1001, MAPQ 60, CIGAR 143M18I145M, AS:i:193;
#   * the secondary score is XS:i:37.
# Modes: the default certified ladder, --compat=bwa-mem2, the full-width
# ladder (--no-band-cert, --no-adaptive-band: ACCEPT_PAIR's non-cert branch),
# and --adaptive-band, whose 8-bit tier (where this read's short extensions
# run) must stay the exact ladder so the flag remains a no-op on short reads.
#
# Usage: test/tight_band_xs_test.sh <bwa-mem3-binary> <fixtures-dir>

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 <bwa-mem3-binary> <fixtures-dir>" >&2
    exit 2
fi

bin="$1"
fixtures="$2"
ref_src="$fixtures/tight_band_xs_ref.fa"
reads="$fixtures/tight_band_xs_read.fq"

[[ -x "$bin" ]] || {
    echo "FAIL: bwa-mem3 binary not executable at $bin" >&2
    exit 1
}
[[ -s "$ref_src" && -s "$reads" ]] || {
    echo "FAIL: tight_band_xs fixtures missing under $fixtures" >&2
    exit 1
}

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

# Index a copy so the fixtures directory stays clean.
cp "$ref_src" "$tmp/ref.fa"
"$bin" index "$tmp/ref.fa" > "$tmp/index.log" 2>&1 || {
    echo "FAIL: bwa-mem3 index on tight_band_xs_ref.fa failed" >&2
    cat "$tmp/index.log" >&2
    exit 1
}

# $1: label; remaining args: extra `mem` options.
check() {
    local label="$1"
    shift
    local sam="$tmp/$label.sam"
    "$bin" mem "$@" "$tmp/ref.fa" "$reads" > "$sam" 2> "$tmp/$label.log" || {
        echo "FAIL: [$label] bwa-mem3 mem exited non-zero" >&2
        cat "$tmp/$label.log" >&2
        exit 1
    }
    local rec
    rec="$(grep -v '^@' "$sam" || true)"
    if [[ -z "$rec" || "$(printf '%s\n' "$rec" | wc -l | tr -d ' ')" != 1 ]]; then
        echo "FAIL: [$label] expected exactly one SAM record" >&2
        exit 1
    fi
    local core tags
    core="$(printf '%s\n' "$rec" | cut -f2,4,5,6)"
    tags="$(printf '%s\n' "$rec" | cut -f12- | tr '\t' '\n')"
    [[ "$core" == $'16\t1001\t60\t143M18I145M' ]] || {
        echo "FAIL: [$label] FLAG/POS/MAPQ/CIGAR $(tr '\t' ' ' <<< "$core"), expected 16 1001 60 143M18I145M" >&2
        exit 1
    }
    grep -qx 'AS:i:193' <<< "$tags" || {
        echo "FAIL: [$label] primary score changed: $(tr '\n' ' ' <<< "$tags")" >&2
        exit 1
    }
    grep -qx 'XS:i:37' <<< "$tags" || {
        echo "FAIL: [$label] expected XS:i:37: $(tr '\n' ' ' <<< "$tags")" >&2
        exit 1
    }
    echo "PASS: [$label] 16 1001 60 143M18I145M AS:i:193 XS:i:37"
}

check native
check compat --compat=bwa-mem2
check no_band_cert --no-band-cert
check no_adaptive_band --no-adaptive-band
check adaptive_band --adaptive-band
