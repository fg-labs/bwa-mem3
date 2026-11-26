#!/usr/bin/env bash
# test/regression/version_banner.sh
#
# Asserts that 'bwa-mem3 version' prints the SIMD floor and runtime banner
# lines added by the multi-arch-deployment plan.
#
# Inputs:
#   BWA_MEM3        — path to the bwa-mem3 binary under test
#   EXPECTED_FLOOR  — (optional) expected build floor (e.g. "avx2"); if set,
#                     the floor tier is checked to match exactly, and for avx2
#                     and avx512bw so is the whole floor line (and, for avx2,
#                     every example of it quoted in docs/src)

set -euo pipefail

# Checked explicitly rather than with ${VAR:?}, which exits with a bare shell
# error (and before any ERR trap runs), so the failure would carry no FAIL: line.
if [[ -z ${BWA_MEM3:-} ]]; then
    echo "FAIL: BWA_MEM3 must be set" >&2
    exit 1
fi

OUT="$(mktemp -t bwamem3-version-XXXXXX)"
trap 'rm -f "$OUT"' EXIT

# --- Scenario 1: default invocation prints floor + runtime, no warning ---
"$BWA_MEM3" version > "$OUT" 2>&1
rc=$?

if [[ $rc -ne 0 ]]; then
    echo "FAIL: bwa-mem3 version exited $rc (expected 0)" >&2
    cat "$OUT" >&2
    exit 1
fi

if ! grep -q '^SIMD floor: ' "$OUT"; then
    echo "FAIL: 'SIMD floor:' line not found in version output" >&2
    cat "$OUT" >&2
    exit 1
fi

if ! grep -q '^SIMD runtime: ' "$OUT"; then
    echo "FAIL: 'SIMD runtime:' line not found in version output" >&2
    cat "$OUT" >&2
    exit 1
fi

floor_line="$(grep '^SIMD floor: ' "$OUT" | head -1)"

# The parenthetical names the oldest CPU the floor runs on, never a psABI level:
# the floor check verifies only the tier's own subsets, not a full x86-64-v2/-v3/-v4
# set, so such a label would overstate what the binary was checked for.
if [[ $floor_line == *x86-64-v* ]]; then
    echo "FAIL: floor line names a psABI level: $floor_line" >&2
    exit 1
fi

if [[ -n "${EXPECTED_FLOOR:-}" ]]; then
    actual="$(sed -E 's/^SIMD floor: ([a-z0-9]+).*/\1/' <<< "$floor_line")"
    if [[ "$actual" != "$EXPECTED_FLOOR" ]]; then
        echo "FAIL: floor mismatch — expected '$EXPECTED_FLOOR', got '$actual'" >&2
        cat "$OUT" >&2
        exit 1
    fi

    # The whole line for the x86 floors a build can have (arch=avx2, arch=avx512bw).
    case "$EXPECTED_FLOOR" in
        avx2) floor_label="Haswell 2013+" ;;
        avx512bw) floor_label="Skylake-X 2017+" ;;
        *) floor_label="" ;;
    esac
    if [[ -n $floor_label ]]; then
        want="SIMD floor: $EXPECTED_FLOOR ($floor_label); kernels: sse41 sse42 avx avx2 avx512bw"
        if [[ "$floor_line" != "$want" ]]; then
            echo "FAIL: floor line mismatch — expected '$want', got '$floor_line'" >&2
            exit 1
        fi

        # The docs quote the floor line of a default (avx2) build as example output;
        # each quote must match what the binary prints (runs of spaces aside).
        docs_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/docs/src"
        if [[ $EXPECTED_FLOOR == avx2 && -d $docs_dir ]]; then
            n_quotes=0
            while IFS= read -r quote; do
                n_quotes=$((n_quotes + 1))
                quote_file="${quote%%:*}"
                quote_line="$(tr -s ' ' <<< "${quote#*:}")"
                if [[ "$quote_line" != "$want" ]]; then
                    echo "FAIL: $quote_file quotes '$quote_line'; the binary prints '$want'" >&2
                    exit 1
                fi
            done < <(grep -rE --include='*.md' '^SIMD floor: +avx2 ' "$docs_dir" || true)
            if ((n_quotes == 0)); then
                echo "FAIL: no documented 'SIMD floor: avx2' example found under $docs_dir" >&2
                exit 1
            fi
        fi
    fi
fi

# --- Scenario 2: BWAMEM3_FORCE_TIER downgrades the runtime line ---
#
# Proves that the runtime line reads g_tier (FORCE_TIER-influenced) while
# the warning reads g_host_capability (raw host capability). The forced
# tier MUST appear in the runtime line; the [W::bwa-mem3] warning MUST
# NOT appear, because g_host_capability is unchanged.
#
# Only meaningful on x86 hosts where FORCE_TIER is honoured. On arm64,
# cross-family FORCE_TIER requests are rejected by the existing
# simd_dispatch.cpp logic, so skip the test.
if uname -m | grep -qE '^(x86_64|amd64)$'; then
    BWAMEM3_FORCE_TIER=sse41 "$BWA_MEM3" version > "$OUT" 2>&1
    rc=$?
    if [[ $rc -ne 0 ]]; then
        echo "FAIL: bwa-mem3 version with FORCE_TIER=sse41 exited $rc (expected 0)" >&2
        cat "$OUT" >&2
        exit 1
    fi
    if ! grep -q 'SIMD runtime: sse41 (BWAMEM3_FORCE_TIER=sse41)' "$OUT"; then
        echo "FAIL: runtime line did not show forced tier with parenthetical" >&2
        cat "$OUT" >&2
        exit 1
    fi
    if grep -q '\[W::bwa-mem3\]' "$OUT"; then
        echo "FAIL: warning line appeared under FORCE_TIER (g_host_capability should be unchanged)" >&2
        cat "$OUT" >&2
        exit 1
    fi
    echo "PASS: version banner has SIMD floor + runtime, FORCE_TIER tested on x86"
else
    echo "PASS: version banner has SIMD floor + runtime (FORCE_TIER scenario skipped on $(uname -m))"
fi
