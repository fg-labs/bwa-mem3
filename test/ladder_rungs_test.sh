#!/usr/bin/env bash
# test/ladder_rungs_test.sh
#
# Regression for the extension retry ladder's rung count. Upstream bwa and
# bwa-mem2 try w = 100 and then 200 and stop (MAX_BAND_TRY 2). bwa-mem3 carried
# a four-rung ladder [100, 200, 400, 800] from #58, so an extension whose
# stop test failed at 200 (score changed, max_off >= 150) ran on to 400 and
# could commit a different alignment and a wider a->w.
#
# Three reads, built from phix slices by fixtures/ladder_proofs/build_fixture.awk:
#
#   ladder     900 bp: a 210 bp exact seed and then three segments that each net
#              +19 over an insertion, at diagonal offsets 80, 150 and 250. At
#              w = 200 the alignment reaches offset 150 with a changed score, so
#              the stop test fails there; only a third rung reaches the segment
#              at offset 250. The four-rung ladder reported
#              209M80I141M70I130M100I170M with AS:i:267.
#   ladder_rc  the same read reverse-complemented: the seed is at the 3' end of
#              the query, so the staircase is a left extension.
#   adaptive   554 bp: gains at offsets 16, 32 and 64, just past 3/4 of the
#              --adaptive-band narrowing rungs 20, 40 and 80. That ladder keeps
#              four rungs (ADAPTIVE_BAND_TRY); with two it stops at 40 and
#              clips the last segment.
#
# Each read must give bwa 0.7.19's record (FLAG, POS, MAPQ, CIGAR, NM, MD, AS,
# XS) under the default certified ladder, --compat=bwa-mem2, --no-band-cert and
# --no-adaptive-band (the full-width ladder via either flag), and the adaptive
# read also under --adaptive-band. The 900 bp reads are not checked under
# --adaptive-band, which is not byte-identical by design on them.
#
# Usage: test/ladder_rungs_test.sh <bwa-mem3-binary> <fixtures-dir>

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 <bwa-mem3-binary> <fixtures-dir>" >&2
    exit 2
fi

bin="$1"
fixtures="$2"
phix="$fixtures/phix.fa"
awk_src="$fixtures/ladder_proofs/build_fixture.awk"

[[ -x "$bin" ]] || {
    echo "FAIL: bwa-mem3 binary not executable at $bin" >&2
    exit 1
}
[[ -s "$phix" && -s "$awk_src" ]] || {
    echo "FAIL: phix.fa or ladder_proofs/build_fixture.awk missing under $fixtures" >&2
    exit 1
}

AWK="${AWK:-$(command -v mawk || command -v gawk || command -v awk)}"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

"$AWK" -v MODE=ref -f "$awk_src" "$phix" > "$tmp/ref.fa"
for c in ladder ladder_rc adaptive; do
    "$AWK" -v MODE=read -v CASE="$c" -f "$awk_src" "$phix" > "$tmp/$c.fq"
done
"$bin" index "$tmp/ref.fa" > "$tmp/index.log" 2>&1 || {
    echo "FAIL: bwa-mem3 index failed" >&2
    cat "$tmp/index.log" >&2
    exit 1
}

# $1: read case; $2: label; $3: expected FLAG/POS/MAPQ/CIGAR (tab-separated);
# $4: expected tags (space-separated, each must be present); remaining args:
# extra `mem` options.
check() {
    local case="$1" label="$2" want_core="$3" want_tags="$4"
    shift 4
    local sam="$tmp/$case.$label.sam"
    "$bin" mem "$@" "$tmp/ref.fa" "$tmp/$case.fq" > "$sam" 2> "$tmp/$case.$label.log" || {
        echo "FAIL: [$case $label] bwa-mem3 mem exited non-zero" >&2
        cat "$tmp/$case.$label.log" >&2
        exit 1
    }
    local rec
    rec="$(grep -v '^@' "$sam" || true)"
    if [[ -z "$rec" || "$(printf '%s\n' "$rec" | wc -l | tr -d ' ')" != 1 ]]; then
        echo "FAIL: [$case $label] expected exactly one SAM record" >&2
        exit 1
    fi
    local core tags t
    core="$(printf '%s\n' "$rec" | cut -f2,4,5,6)"
    tags="$(printf '%s\n' "$rec" | cut -f12- | tr '\t' '\n')"
    [[ "$core" == "$want_core" ]] || {
        echo "FAIL: [$case $label] FLAG/POS/MAPQ/CIGAR $(tr '\t' ' ' <<< "$core"), expected $(tr '\t' ' ' <<< "$want_core")" >&2
        exit 1
    }
    for t in $want_tags; do
        grep -qxF "$t" <<< "$tags" || {
            echo "FAIL: [$case $label] missing $t: $(tr '\n' ' ' <<< "$tags")" >&2
            exit 1
        }
    done
    echo "PASS: [$case $label] $(tr '\t' ' ' <<< "$core") $want_tags"
}

stair_md='MD:Z:227C17G17T17T17T17T17G31A17T17A17T17C17A17'
stair_tags="NM:i:163 $stair_md AS:i:248 XS:i:0"
adapt_md='MD:Z:117G17G17C17G17A17T29G17A17G17G17A17C29A17C17T17G17C17A17G17C6'
adapt_tags="NM:i:84 $adapt_md AS:i:308 XS:i:0"
for mode in "native" "compat --compat=bwa-mem2" "no_band_cert --no-band-cert" \
    "no_adaptive_band --no-adaptive-band"; do
    read -r label opts <<< "$mode"
    # shellcheck disable=SC2086 # opts is one flag or empty
    check ladder "$label" $'0\t1001\t60\t209M80I141M70I125M275S' "$stair_tags" $opts
    # shellcheck disable=SC2086
    check ladder_rc "$label" $'16\t1001\t60\t209M80I141M70I125M275S' "$stair_tags" $opts
    # shellcheck disable=SC2086
    check adaptive "$label" $'0\t1001\t60\t100M16I120M16I120M32I150M' "$adapt_tags" $opts
done
check adaptive adaptive_band $'0\t1001\t60\t100M16I120M16I120M32I150M' "$adapt_tags" --adaptive-band
