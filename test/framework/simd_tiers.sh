#!/usr/bin/env bash
# test/framework/simd_tiers.sh
#
# Source-only helper for regression scripts that run one bwa-mem3 binary under
# each tier with a batched kswv kernel that it can use on this host
# (BWAMEM3_FORCE_TIER), so that every such script sweeps the same tiers and
# parses the dispatcher the same way. Used by all_tiers_parity.sh,
# band_cert_identity.sh and, through test/framework/mem_cli.sh,
# mem_cli_basics.sh. (dash_read_as_n.sh and rescue_prune_identity.sh sweep a
# different set, the sub-AVX2 tiers included, and parse the dispatcher
# themselves.) Keeps to bash 3.2 (macOS /bin/bash).
#
#     source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../framework" && pwd)/simd_tiers.sh"
#     simd_sweep_tiers "$BWA_MEM3" || exit $?
#     for t in "${SIMD_TIERS[@]}"; do ...; done

# The tiers a sweep covers, low to high.
#
# sse41/sse42/avx are DELIBERATELY ABSENT. There is no batched kswv kernel below
# AVX2: src/kswv.cpp's SSE-only getScores8/getScores16 are stubs that
# exit(EXIT_FAILURE). On a multi-tier build the main TUs are compiled at
# BASELINE_ARCH (avx2 by default), and the batched mate-rescue path is the only
# one there (the scalar pre-AVX2 fallback was removed), so mem_sam_pe_batch_run
# calls those stubs the moment a pair needs mate rescue.
# Forcing any of those three tiers on paired-end input therefore kills the
# aligner: exit 1, one [E::getScores8] per thread, header-only SAM. Measured
# 2026-07-26 on a c7i.4xlarge, all three tiers, hg38 + 32k pairs.
#
# This is not a user-facing gap. Pre-AVX2 x86 (sse41/sse42/avx) is no longer a
# supported tier: the batched mate-rescue kernels require AVX2+, the scalar
# (pre-AVX2) mate-rescue fallback was removed, and the Makefile now refuses
# BASELINE_ARCH<avx2 at configure time (so no crashing binary can be built). A
# default build also refuses to start below its floor (the simd_dispatch.cpp
# host precheck, covered by host_floor_enforce.sh). Restore a tier below avx2
# here only alongside a runtime scalar fallback for tiers with no batched
# kernel.
SIMD_SWEEP_TIERS=(avx2 avx512bw neon)

# Every tier name the dispatcher can report, batched or not, low to high. Kept
# separate from SIMD_SWEEP_TIERS so that narrowing the sweep does not turn a
# correctly-detected low-tier host into a bogus "detection broke" failure: an
# sse41-baseline build reports sse41, which is a perfectly good detection and
# simply has no tiers to sweep.
SIMD_KNOWN_TIERS=(sse41 sse42 avx avx2 avx512bw neon)

# simd_sweep_tiers <bwa-mem3>
#   Sets SIMD_HOST_TIER to the tier the binary dispatches to on this host
#   (with any BWAMEM3_FORCE_TIER in the caller's environment ignored), and
#   SIMD_TIERS to the sweepable tiers at or below it, lowest first: avx2 and,
#   where the host has it, avx512bw on x86; neon on arm64. SIMD_TIERS is empty
#   on a host below the avx2 floor; the caller decides whether that is a SKIP
#   or a failure. Returns 2, with the dispatcher's output on stderr, when the
#   host tier cannot be read.
simd_sweep_tiers() {
    local raw
    # bwa-mem3 with no args prints usage and exits non-zero; the dispatcher's
    # debug line goes to stderr.
    raw="$(env -u BWAMEM3_FORCE_TIER BWAMEM3_DEBUG_SIMD=1 "$1" 2>&1 || true)"
    SIMD_HOST_TIER="$(printf '%s\n' "$raw" | sed -n 's/.*SIMD tier: \([a-z0-9]*\).*/\1/p' | head -n 1)"
    SIMD_TIERS=()

    # Rank against SIMD_KNOWN_TIERS rather than walking SIMD_SWEEP_TIERS until
    # the host tier matches: an sse41 host never appears in the sweep list, so
    # such a walk would never stop and would "test" tiers above the host, which
    # BWAMEM3_FORCE_TIER (downgrade-only) refuses.
    local host_rank=-1 i t
    for i in "${!SIMD_KNOWN_TIERS[@]}"; do
        if [[ ${SIMD_KNOWN_TIERS[$i]} == "$SIMD_HOST_TIER" ]]; then
            host_rank=$i
            break
        fi
    done
    if [[ $host_rank -lt 0 ]]; then
        echo "FAIL: could not detect host SIMD tier from BWAMEM3_DEBUG_SIMD output" >&2
        echo "  parsed tier '$SIMD_HOST_TIER' (expected one of: ${SIMD_KNOWN_TIERS[*]})" >&2
        echo "  raw output was:" >&2
        printf '%s\n' "$raw" >&2
        return 2
    fi

    if [[ $SIMD_HOST_TIER == neon ]]; then
        # arm64 has exactly one tier; the dispatcher always picks NEON.
        SIMD_TIERS=(neon)
        return 0
    fi
    for t in "${SIMD_SWEEP_TIERS[@]}"; do
        [[ $t == neon ]] && continue
        for i in "${!SIMD_KNOWN_TIERS[@]}"; do
            if [[ ${SIMD_KNOWN_TIERS[$i]} == "$t" ]]; then
                if [[ $i -le $host_rank ]]; then SIMD_TIERS+=("$t"); fi
                break
            fi
        done
    done
}
