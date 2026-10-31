#!/usr/bin/env bash
# test/framework/simd_tiers.sh
#
# Source-only helper for regression scripts that read which SIMD tier a
# bwa-mem3 binary dispatches to, or run it under each tier with a batched kswv
# kernel that it can use on this host (BWAMEM3_FORCE_TIER), so that every such
# script parses the dispatcher, ranks tiers and sweeps them the same way. Used
# by all_tiers_parity.sh, band_cert_identity.sh, dash_read_as_n.sh,
# meth_context_identity.sh, rescue_prune_identity.sh and, through
# test/framework/mem_cli.sh, mem_cli_basics.sh. Keeps to bash 3.2 (macOS /bin/bash).
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
# simply has no tiers to sweep. scalar is BWAMEM3_TIER_NONE (no SIMD at all;
# src/simd_dispatch.cpp's tier names).
SIMD_KNOWN_TIERS=(scalar sse41 sse42 avx avx2 avx512bw neon)

# simd_tier_rank <tier>
#   Prints the tier's index in SIMD_KNOWN_TIERS (low to high), or -1 for a name
#   the dispatcher never reports.
simd_tier_rank() {
    local i
    for i in "${!SIMD_KNOWN_TIERS[@]}"; do
        if [[ ${SIMD_KNOWN_TIERS[$i]} == "$1" ]]; then
            echo "$i"
            return 0
        fi
    done
    echo -1
}

# simd_parse_tier <file>
#   Prints the tier named by the first BWAMEM3_DEBUG_SIMD "SIMD tier: <tier>"
#   line in <file>, or nothing. A run's own stderr, kept with
#   BWAMEM3_DEBUG_SIMD=1, shows the tier it ran at after any force.
simd_parse_tier() {
    sed -n 's/.*SIMD tier: \([a-z0-9]*\).*/\1/p' "$1" | head -n 1
}

# simd_detect_tier <command...>
#   Runs <command...> (a bwa-mem3 binary, optionally behind `env ...`) with
#   BWAMEM3_DEBUG_SIMD=1 and sets SIMD_HOST_TIER to the tier its dispatcher
#   reports. A BWAMEM3_FORCE_TIER in the environment counts, so SIMD_HOST_TIER
#   is then the forced tier, not the host's own (simd_sweep_tiers unsets it to
#   read the host's). Returns 2, with the dispatcher's output on stderr, when
#   that is not a tier in SIMD_KNOWN_TIERS.
simd_detect_tier() {
    local raw
    # bwa-mem3 with no args prints usage and exits non-zero; the dispatcher's
    # debug line goes to stderr.
    raw="$(BWAMEM3_DEBUG_SIMD=1 "$@" 2>&1 || true)"
    SIMD_HOST_TIER="$(simd_parse_tier <(printf '%s\n' "$raw"))"
    if [[ $(simd_tier_rank "$SIMD_HOST_TIER") -lt 0 ]]; then
        echo "FAIL: could not detect host SIMD tier from BWAMEM3_DEBUG_SIMD output" >&2
        echo "  parsed tier '$SIMD_HOST_TIER' (expected one of: ${SIMD_KNOWN_TIERS[*]})" >&2
        echo "  raw output was:" >&2
        printf '%s\n' "$raw" >&2
        return 2
    fi
}

# simd_sweep_tiers <bwa-mem3>
#   Sets SIMD_HOST_TIER to the tier the binary dispatches to on this host
#   (with any BWAMEM3_FORCE_TIER in the caller's environment ignored), and
#   SIMD_TIERS to the sweepable tiers at or below it, lowest first: avx2 and,
#   where the host has it, avx512bw on x86; neon on arm64. SIMD_TIERS is empty
#   on a host below the avx2 floor; the caller decides whether that is a SKIP
#   or a failure. Returns 2, with the dispatcher's output on stderr, when the
#   host tier cannot be read.
simd_sweep_tiers() {
    SIMD_TIERS=()
    simd_detect_tier env -u BWAMEM3_FORCE_TIER "$1" || return 2

    if [[ $SIMD_HOST_TIER == neon ]]; then
        # arm64 has exactly one tier; the dispatcher always picks NEON.
        SIMD_TIERS=(neon)
        return 0
    fi
    # Rank against SIMD_KNOWN_TIERS rather than walking SIMD_SWEEP_TIERS until
    # the host tier matches: an sse41 host never appears in the sweep list, so
    # such a walk would never stop and would "test" tiers above the host, which
    # BWAMEM3_FORCE_TIER (downgrade-only) refuses.
    # A sweep entry missing from SIMD_KNOWN_TIERS (rank -1) is skipped, not
    # swept as if it were below every tier.
    local host_rank rank t
    host_rank="$(simd_tier_rank "$SIMD_HOST_TIER")"
    for t in "${SIMD_SWEEP_TIERS[@]}"; do
        [[ $t == neon ]] && continue
        rank="$(simd_tier_rank "$t")"
        if [[ $rank -ge 0 && $rank -le $host_rank ]]; then SIMD_TIERS+=("$t"); fi
    done
}
