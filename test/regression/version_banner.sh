#!/usr/bin/env bash
# test/regression/version_banner.sh
#
# Asserts that 'bwa-mem3 version' prints the SIMD floor and runtime banner
# lines added by the multi-arch-deployment plan.
#
# Inputs:
#   BWA_MEM3        — path to the bwa-mem3 binary under test
#   EXPECTED_FLOOR  — (optional) expected build floor (e.g. "avx2"); if set,
#                     the floor line is checked to match exactly

set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"

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

if [[ -n "${EXPECTED_FLOOR:-}" ]]; then
    actual="$(grep '^SIMD floor: ' "$OUT" | sed -E 's/^SIMD floor: ([a-z0-9]+).*/\1/' | head -1)"
    if [[ "$actual" != "$EXPECTED_FLOOR" ]]; then
        echo "FAIL: floor mismatch — expected '$EXPECTED_FLOOR', got '$actual'" >&2
        cat "$OUT" >&2
        exit 1
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
