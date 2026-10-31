#!/usr/bin/env bash
# test/regression/mem_synth_ref_check.sh
#
# The committed synthetic reference test/fixtures/mem_synth.fa is what its
# generator (test/tools/mem_synth_ref.c) writes, with the planted features and
# the 19-mer uniqueness that the mem CLI scripts' expected values rely on;
# `mem_synth_ref --check` asserts all three. A checker that stopped failing
# would pass the committed file whatever it held, so first confirm that
# --check rejects three corrupted copies, each failing a different check:
#
#   copy_broken   one base of the ctg2 copy of block R changed: the bytes
#                 differ and ctg2:1001-1300 is no longer a copy of R.
#   dup_19mer     ctg1:101-119 written over ctg1:4001-4019, outside every
#                 feature: that 19-mer now occurs twice.
#   short_contig  the last line of ctg4 dropped: ctg4 is shorter than 600 bp.
#
# then that it accepts the committed file.
#
# Env: MEM_SYNTH_REF -- the built generator (test/mem_synth_ref, from
# `make test-binaries`); TMPDIR (optional) for the scratch directory.

set -Eeuo pipefail
# -E: a command failing inside a function (mutate's awk) also reaches this trap,
# so an unexpected failure is labelled FAIL rather than ending the script silently.
# shellcheck disable=SC2154  # rc is assigned by `rc=$?` at the start of the trap body
trap 'rc=$?; echo "FAIL: mem_synth_ref_check: unexpected command failure at line $LINENO (exit $rc)" >&2; exit "$rc"' ERR

: "${MEM_SYNTH_REF:?MEM_SYNTH_REF must be set to the built test/mem_synth_ref}"
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fixture="$repo/test/fixtures/mem_synth.fa"
work="$(mktemp -d "${TMPDIR:-/tmp}/mem_synth_ref_check.XXXXXX")"
trap 'rm -rf "$work"' EXIT

# mutate <mode> <out.fa>: rewrites the fixture's contigs, 60 bases per line.
mutate() {
    awk -v mode="$1" '
        /^>/ { name[++n] = substr($0, 2); next }
        { seq[n] = seq[n] $0 }
        END {
            next_base["A"] = "C"; next_base["C"] = "G"; next_base["G"] = "T"; next_base["T"] = "A"
            if (mode == "copy_broken") {
                b = substr(seq[2], 1150, 1)
                seq[2] = substr(seq[2], 1, 1149) next_base[b] substr(seq[2], 1151)
            } else if (mode == "dup_19mer") {
                seq[1] = substr(seq[1], 1, 4000) substr(seq[1], 101, 19) substr(seq[1], 4020)
            } else if (mode == "short_contig") {
                seq[4] = substr(seq[4], 1, length(seq[4]) - 60)
            }
            for (i = 1; i <= n; i++) {
                print ">" name[i]
                for (j = 1; j <= length(seq[i]); j += 60) print substr(seq[i], j, 60)
            }
        }' "$fixture" > "$2"
}

fails=0
# expect_reject <mode> <pattern>: --check must exit 1 and print <pattern>.
expect_reject() {
    local fa="$work/$1.fa" rc=0
    mutate "$1" "$fa"
    "$MEM_SYNTH_REF" --check "$fa" > "$work/$1.out" 2>&1 || rc=$?
    if [[ $rc -eq 1 ]] && grep -q -- "$2" "$work/$1.out"; then
        echo "ok: --check rejects $1"
    else
        echo "FAIL: --check on $1 exited $rc without '$2':" >&2
        head -n 5 "$work/$1.out" >&2
        fails=$((fails + 1))
    fi
}

expect_reject copy_broken "ctg2:1001-1300 is not a copy of R"
expect_reject dup_19mer "19-mer at ctg1:4001 occurs 2 times forward"
expect_reject short_contig "ctg4 is 540 bp, expected 600"

if ((fails > 0)); then
    echo "FAIL: mem_synth_ref_check: --check accepted $fails corrupted cop(y/ies); its verdict on the fixture means nothing"
    exit 1
fi
if ! "$MEM_SYNTH_REF" --check "$fixture"; then
    echo "FAIL: mem_synth_ref_check: $fixture does not pass --check"
    exit 1
fi
echo "PASS: mem_synth_ref_check: --check rejects 3 corrupted copies and accepts the committed fixture"
