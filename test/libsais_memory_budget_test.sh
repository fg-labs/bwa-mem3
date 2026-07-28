#!/usr/bin/env bash
# Memory-budget-honored test. Builds the 1 Mbp synthetic fixture under
# several --max-memory settings, measures actual peak RSS via /usr/bin/time,
# and fails if peak exceeds the budget by more than 10% slack. The preflight
# should gate any case where the libsais estimate would exceed the budget.
#
# Parses /usr/bin/time -l output on macOS and /usr/bin/time -v on Linux.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
BWAMEM3="$ROOT/bwa-mem3"
FA_SRC="$HERE/fixtures/synthetic_1mb.fa"
[[ -s "$FA_SRC" ]] || { echo "FAIL: $FA_SRC missing"; exit 1; }

UNAME_S="$(uname -s)"
if [[ "$UNAME_S" == "Darwin" ]]; then
    TIME_CMD=(/usr/bin/time -l)
else
    TIME_CMD=(/usr/bin/time -v)
fi

parse_peak() {
    # $1: timing output path. echoes peak RSS in bytes.
    local out="$1"
    if [[ "$UNAME_S" == "Darwin" ]]; then
        grep -E '^ *[0-9]+ +maximum resident set size' "$out" \
            | awk '{print $1}'
    else
        grep -E 'Maximum resident set size' "$out" \
            | awk -F': +' '{print $2 * 1024}'
    fi
}

# Budgets chosen (128 / 512 / 2048 MiB) exercise the preflight + actual
# peak-RSS contract above the fixed thread-pool / libomp / mimalloc-arena
# infrastructure cost (~15 MiB on this host that's independent of the
# budget). Tighter budgets below ~64 MiB would fail the +10% slack check
# even on a trivially small input, so we lean on the preflight test
# (libsais_index_diff) to cover the "does the budget get respected"
# contract at larger scales.
FAIL=0
for budget_mib in 128 512 2048; do
    budget_bytes=$(( budget_mib * 1024 * 1024 ))
    slack_bytes=$(( budget_bytes * 11 / 10 ))
    TD="$(mktemp -d)"
    # Cleanup unconditionally on exit (parse-failure / set -e / signal)
    # so $TD doesn't leak between iterations or on the early-exit paths.
    trap 'rm -rf "$TD"' EXIT
    cp "$FA_SRC" "$TD/t.fa"
    TIMING="$TD/time.out"
    "${TIME_CMD[@]}" "$BWAMEM3" index --max-memory "${budget_mib}M" "$TD/t.fa" >"$TD/stdout" 2>"$TIMING" || {
        echo "FAIL: build failed at --max-memory ${budget_mib}M"
        cat "$TIMING" | tail -20
        exit 1
    }
    peak="$(parse_peak "$TIMING")"
    [[ -n "$peak" ]] || { echo "FAIL: could not parse peak RSS"; cat "$TIMING"; exit 1; }
    if [[ "$peak" -gt "$slack_bytes" ]]; then
        echo "FAIL: --max-memory ${budget_mib}M -> peak $(( peak / 1024 / 1024 )) MiB (budget ${budget_mib}M + 10% = $(( slack_bytes / 1024 / 1024 )) MiB)"
        FAIL=1
    else
        printf "OK:   --max-memory %dM -> peak %d MiB (under budget+10%%)\n" \
            "$budget_mib" "$(( peak / 1024 / 1024 ))"
    fi
    rm -rf "$TD"
    trap - EXIT
done

# `--meth` builds two indexes sequentially in one process: the original over
# N = 2*l_pac, then a seed over a per-strand-converted text twice as long
# (4*l_pac), so the seed costs ~2x the original and decides whether the
# invocation can run. Pick a budget between the two estimates: the original
# would fit, the seed cannot. The preflight must refuse UP FRONT and leave no
# index behind -- previously the seed's own check fired only once that build
# started, i.e. after the original index had been built and written (an hour
# into an hg38 run, leaving a half-populated index directory).
TD="$(mktemp -d)"
trap 'rm -rf "$TD"' EXIT
cp "$FA_SRC" "$TD/m.fa"
# The fixture is ~1 Mbp -> original est ~11.4 MiB, seed est ~22.9 MiB.
if "$BWAMEM3" index --meth --max-memory 16M "$TD/m.fa" >"$TD/out" 2>"$TD/err"; then
    echo "FAIL: --meth --max-memory 16M should have been refused (seed needs ~23 MiB)"
    FAIL=1
else
    rc=$?
    if [[ "$rc" -ne 3 ]]; then
        echo "FAIL: --meth budget refusal exited $rc, expected 3"
        cat "$TD/err"
        FAIL=1
    elif ! grep -q 'seed index dominates' "$TD/err"; then
        echo "FAIL: --meth refusal did not explain that the seed dominates"
        cat "$TD/err"
        FAIL=1
    else
        # The point of the up-front check: nothing was built before refusing.
        leaked="$(find "$TD" -name 'm.fa.*' -not -name '*.fai' | head -5)"
        if [[ -n "$leaked" ]]; then
            echo "FAIL: --meth refusal left index artifacts behind:"
            echo "$leaked"
            FAIL=1
        else
            echo "OK:   --meth --max-memory 16M refused up front (exit 3, nothing written)"
        fi
    fi
fi
rm -rf "$TD"
trap - EXIT

if [[ $FAIL -ne 0 ]]; then
    echo "libsais_memory_budget_test: FAILED"
    exit 1
fi
echo "OK: libsais_memory_budget (all budgets respected)"
