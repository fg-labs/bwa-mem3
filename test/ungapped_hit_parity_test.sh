#!/usr/bin/env bash
# test/ungapped_hit_parity_test.sh
#
# Regression for the ungapped extension fast path (src/ungapped_ext.h). When it
# reports a HIT the banded-SW ladder is skipped and the region is committed from
# the diagonal walk. Contract: a HIT commits exactly the fields the extension
# kernel would, so each read below gives bwa 0.7.19's record. Three reads
# (150 bp, phix slices with substitutions in their first bases; see
# fixtures/ladder_proofs/build_fixture.awk) each isolate one channel, at the
# non-default parameter that exposes it:
#
#   L0tie   -L 0   the diagonal returns to its record without exceeding it; the
#                  kernels keep the earlier row (qle 0), so bwa clips: 5S145M.
#                  A walk that keeps the later tied row reports 150M NM:1.
#   O20tie  -O 20  the same tie followed by two more mismatches, so the clip
#                  decision is branch A even at the default -L 5: 7S143M
#                  (2S148M from a later-row walk).
#   d3zdrop -d 3   one mismatch drops the kernels' score by 4 > zdrop, so they
#                  stop: 11S139M AS:139 (a walk without the z-drop guard goes on
#                  to 150M AS:145).
#
# Asserts bwa 0.7.19's record for each (FLAG/POS/MAPQ/CIGAR, NM, MD, AS, XS),
# taken from `bwa mem <opt>` on the same reference and read, in the default
# mode, --compat=bwa-mem2, --no-band-cert and --adaptive-band.
#
# Usage: test/ungapped_hit_parity_test.sh <bwa-mem3-binary> <fixtures-dir>

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
"$bin" index "$tmp/ref.fa" > "$tmp/index.log" 2>&1 || {
    echo "FAIL: bwa-mem3 index failed" >&2
    cat "$tmp/index.log" >&2
    exit 1
}

# $1: case; $2: expected FLAG/POS/MAPQ/CIGAR (tab-separated); $3: expected tags
# (space-separated, each must be present); $4: label; remaining args: `mem`
# options.
check() {
    local case_="$1" expect_core="$2" expect_tags="$3" label="$4"
    shift 4
    local reads="$tmp/$case_.fq" sam="$tmp/$case_.$label.sam"
    [[ -s "$reads" ]] || "$AWK" -v MODE=read -v CASE="$case_" -f "$awk_src" "$phix" > "$reads"
    "$bin" mem "$@" "$tmp/ref.fa" "$reads" > "$sam" 2> "$sam.log" || {
        echo "FAIL: [$case_ $label] bwa-mem3 mem exited non-zero" >&2
        cat "$sam.log" >&2
        exit 1
    }
    local rec
    rec="$(grep -v '^@' "$sam" || true)"
    if [[ -z "$rec" || "$(printf '%s\n' "$rec" | wc -l | tr -d ' ')" != 1 ]]; then
        echo "FAIL: [$case_ $label] expected exactly one SAM record" >&2
        exit 1
    fi
    local core tags t
    core="$(printf '%s\n' "$rec" | cut -f2,4,5,6)"
    tags="$(printf '%s\n' "$rec" | cut -f12- | tr '\t' '\n')"
    [[ "$core" == "$expect_core" ]] || {
        echo "FAIL: [$case_ $label] FLAG/POS/MAPQ/CIGAR $(tr '\t' ' ' <<< "$core"), expected $(tr '\t' ' ' <<< "$expect_core")" >&2
        exit 1
    }
    for t in $expect_tags; do
        grep -qxF "$t" <<< "$tags" || {
            echo "FAIL: [$case_ $label] missing $t: $(tr '\n' ' ' <<< "$tags")" >&2
            exit 1
        }
    done
    echo "PASS: [$case_ $label] $(tr '\t' ' ' <<< "$core") $expect_tags"
}

for mode in "native" "compat --compat=bwa-mem2" "no_band_cert --no-band-cert" "adaptive_band --adaptive-band"; do
    # shellcheck disable=SC2086  # $mode splits into label + option on purpose
    set -- $mode
    label="$1"
    shift
    check L0tie $'0\t1006\t60\t5S145M' "NM:i:0 MD:Z:145 AS:i:145 XS:i:0" "$label" -L 0 "$@"
    check O20tie $'0\t1008\t60\t7S143M' "NM:i:0 MD:Z:143 AS:i:143 XS:i:0" "$label" -O 20 "$@"
    check d3zdrop $'0\t1012\t60\t11S139M' "NM:i:0 MD:Z:139 AS:i:139 XS:i:0" "$label" -d 3 "$@"
done
echo "PASS: ungapped_hit_parity_test (3 reads x 4 modes)"
