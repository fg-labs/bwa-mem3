#!/usr/bin/env bash
# test/regression/all_tiers_parity.sh
#
# Regression: validates that one binary, on one host, produces byte-identical
# SAM whichever of its usable kernels BWAMEM3_FORCE_TIER selects.
#
# SCOPE: avx2 and avx512bw only (plus neon on arm64, where there is one tier).
# sse41/sse42/avx are excluded because there is no batched kswv kernel below
# AVX2 and forcing them kills the aligner — see the ALL_TIERS comment for the
# full reasoning and the measurement behind it.
#
# This script is therefore NOT the cross-build parity check, and should not be
# mistaken for one. "Does a low-baseline build agree with an avx2 build?" is
# answered by CI building at sse41/avx2/arm64 baselines and running
# chr22_parity.sh on each, which asserts byte-identity against bwa itself — a
# stronger oracle than self-comparison, and transitive across the three builds.
#
# Real coverage needs an avx512bw host: on an avx2-only x86 runner, or on arm64,
# there is exactly one sweepable tier and the script reports SKIPPED rather than
# claiming a pass it did not earn.
#
# Inputs:
#   BWA_MEM3       — path to the multi-tier bwa-mem3 binary
#   PARITY_FA      — pre-indexed reference fasta
#   PARITY_R1      — paired-end FASTQ R1
#   PARITY_R2      — paired-end FASTQ R2

set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${PARITY_FA:?PARITY_FA must be set}"
: "${PARITY_R1:?PARITY_R1 must be set}"
: "${PARITY_R2:?PARITY_R2 must be set}"

OUT_DIR="$(mktemp -d -t bwamem3-parity-XXXXXX)"
trap 'rm -rf "$OUT_DIR"' EXIT

# Detect host max tier. The dispatcher's debug message goes to stderr.
# bwa-mem3 with no args prints usage and exits non-zero; capture the output
# under `|| true` so set -euo pipefail doesn't abort the whole script.
HOST_TIER_RAW="$(BWAMEM3_DEBUG_SIMD=1 "$BWA_MEM3" 2>&1 || true)"
HOST_TIER="$(printf '%s\n' "$HOST_TIER_RAW" | sed -n 's/.*SIMD tier: \([a-z0-9]*\).*/\1/p' | head -1)"

# Tiers in low-to-high order. Skip any tier above the host (BWAMEM3_FORCE_TIER
# would refuse the request anyway).
#
# sse41/sse42/avx are DELIBERATELY ABSENT. There is no batched kswv kernel below
# AVX2: src/kswv.cpp's SSE-only getScores8/getScores16 are stubs that
# exit(EXIT_FAILURE). On a multi-tier build the main TUs are compiled at
# BASELINE_ARCH (avx2 by default), so BWAMEM_BATCHED_MATESW is 1 and
# mem_sam_pe_batch_run calls those stubs the moment a pair needs mate rescue.
# Forcing any of those three tiers on paired-end input therefore kills the
# aligner: exit 1, one [E::getScores8] per thread, header-only SAM. Measured
# 2026-07-26 on a c7i.4xlarge, all three tiers, hg38 + 32k pairs.
#
# This is not a user-facing gap. Those tiers are unreachable in production --
# a default build refuses to start below its floor (the simd_dispatch.cpp host
# precheck, covered by host_floor_enforce.sh), and a genuine pre-AVX2 user
# builds BASELINE_ARCH=sse41, which compiles BWAMEM_BATCHED_MATESW to 0 and
# takes the correct scalar mem_sam_pe path. That low-baseline build is covered:
# CI builds it and runs chr22_parity.sh on it, asserting byte-identity against
# bwa itself -- as it does for the avx2 and arm64 builds, so the three builds
# agree with each other transitively. That is a stronger oracle than this
# script's self-comparison, and it is where cross-build parity actually lives.
#
# What remains here is the reachable-tier question: one binary, on one host,
# must produce the same SAM whichever of its usable kernels it dispatches.
# Restore a tier below avx2 here only alongside a runtime scalar fallback for
# tiers with no batched kernel.
ALL_TIERS=(avx2 avx512bw neon)

# Every tier name the dispatcher can report, batched or not. Kept separate from
# ALL_TIERS so that narrowing the sweep does not turn a correctly-detected
# low-tier host into a bogus "detection broke" failure: an sse41-baseline build
# reports HOST_TIER=sse41, which is a perfectly good detection and simply has no
# tiers to sweep.
KNOWN_TIERS=(sse41 sse42 avx avx2 avx512bw neon)

# Fail fast if HOST_TIER didn't parse, or names a tier we don't know about —
# silent fallthrough below would then test every x86 tier on a non-AVX-512 host
# and produce a confusing SIGILL rather than a clear "detection broke" message.
detected_ok=0
for t in "${KNOWN_TIERS[@]}"; do
    if [[ "$t" == "$HOST_TIER" ]]; then
        detected_ok=1
        break
    fi
done
if [[ "$detected_ok" -ne 1 ]]; then
    echo "FAIL: could not detect host SIMD tier from BWAMEM3_DEBUG_SIMD output" >&2
    echo "  parsed HOST_TIER='$HOST_TIER' (expected one of: ${KNOWN_TIERS[*]})" >&2
    echo "  raw output was:" >&2
    printf '%s\n' "$HOST_TIER_RAW" >&2
    exit 2
fi
echo "Host tier: $HOST_TIER"

# Select the sweepable tiers: those that both have a batched kswv kernel
# (ALL_TIERS) and are at or below the host (BWAMEM3_FORCE_TIER is downgrade-only
# and refuses up-tier requests).
#
# Rank against KNOWN_TIERS rather than walking ALL_TIERS until HOST_TIER matches.
# The old walk relied on the host tier appearing in the list it was iterating;
# once sub-AVX2 tiers were dropped from that list, an sse41 host would never
# match, the loop would never break, and it would end up "testing" avx2 and
# avx512bw -- tiers ABOVE the host, which the dispatcher refuses.
host_rank=-1
for i in "${!KNOWN_TIERS[@]}"; do
    if [[ "${KNOWN_TIERS[$i]}" == "$HOST_TIER" ]]; then
        host_rank=$i
        break
    fi
done

if [[ "$HOST_TIER" == "neon" ]]; then
    # arm64 has exactly one tier; the dispatcher always picks NEON.
    TIERS=("neon")
else
    TIERS=()
    for t in "${ALL_TIERS[@]}"; do
        if [[ "$t" == "neon" ]]; then continue; fi
        for i in "${!KNOWN_TIERS[@]}"; do
            if [[ "${KNOWN_TIERS[$i]}" == "$t" ]]; then
                if [[ "$i" -le "$host_rank" ]]; then TIERS+=("$t"); fi
                break
            fi
        done
    done
fi

# A host below the batched floor (an sse41/sse42/avx baseline build) has nothing
# to sweep. Say so plainly instead of reporting an empty pass.
if [[ "${#TIERS[@]}" -eq 0 ]]; then
    echo "ALL TIERS PARITY: SKIPPED -- host tier '$HOST_TIER' is below the batched-kswv"
    echo "  floor (avx2), so this build has no sweepable tiers. Its correctness is covered"
    echo "  by chr22_parity.sh, which runs on this build in CI and asserts byte-identity"
    echo "  against bwa. This is NOT a pass."
    exit 0
fi

echo "Testing tiers: ${TIERS[*]}"

# Pick a reference SAM (the highest tier — usually the host's actual tier).
# Use ${#TIERS[@]}-1 instead of negative indexing so we work on macOS bash 3.2.
REF_TIER="${TIERS[$((${#TIERS[@]} - 1))]}"

# Scoring variants to diff under. Each entry is a label plus extra `mem`
# args. The default run drives the 8-bit kswv mate-rescue kernel; the
# "match4" run raises the match score so mate-rescue pairs clear the
# byte-mode gate at src/bwamem_pair.cpp:217 (KSW_XBYTE set iff
# l_ms*opt->a < 250) and route to getScores16 — without this no end-to-end
# test exercises the 16-bit kswv kernels across tiers.
#
# The gate scales with l_ms (the rescued read length), so the match score
# needed to reach 16-bit depends on the fixture: -A 4 clears it for any
# read >= 63bp, covering all realistic short-read fixtures. If a future
# fixture uses reads shorter than that, raise -A here or the match4 variant
# silently degenerates into a second 8-bit run.
VARIANT_LABELS=(default match4)
VARIANT_ARGS=("" "-A 4")

EXIT=0
COMPARISONS=0
for vi in "${!VARIANT_LABELS[@]}"; do
    vlabel="${VARIANT_LABELS[$vi]}"
    # shellcheck disable=SC2206  # intentional word-splitting of the arg string
    vargs=(${VARIANT_ARGS[$vi]})
    echo "=== Scoring variant: $vlabel (mem args: ${VARIANT_ARGS[$vi]:-none}) ==="
    REF_SAM="$OUT_DIR/$vlabel.$REF_TIER.sam"

    for t in "${TIERS[@]}"; do
        echo ">>> Generating SAM for tier=$t"
        # Check the exit status explicitly. Relying on `set -e` here aborted the
        # whole script with no message at all -- the aligner's stderr goes to a
        # per-tier log, so the operator saw "Generating SAM for tier=X" and a
        # bare exit code. Say what happened and show the log.
        set +e
        BWAMEM3_FORCE_TIER="$t" "$BWA_MEM3" mem -t 1 ${vargs[@]+"${vargs[@]}"} \
            "$PARITY_FA" "$PARITY_R1" "$PARITY_R2" \
            > "$OUT_DIR/$vlabel.$t.sam" 2> "$OUT_DIR/$vlabel.$t.log"
        rc=$?
        set -e
        if [[ $rc -ne 0 ]]; then
            echo "FAIL: [$vlabel] bwa-mem3 exited $rc at tier=$t" >&2
            sed 's/^/      /' "$OUT_DIR/$vlabel.$t.log" | tail -10 >&2
            exit 1
        fi
        echo "    Produced $(wc -c < "$OUT_DIR/$vlabel.$t.sam") bytes"
    done

    for t in "${TIERS[@]}"; do
        if [[ "$t" == "$REF_TIER" ]]; then
            continue
        fi
        if ! diff -q "$REF_SAM" "$OUT_DIR/$vlabel.$t.sam" > /dev/null; then
            echo "FAIL: [$vlabel] tier $t differs from $REF_TIER"
            # `|| true`: we are already inside the "files differ" branch, so
            # this diff exits 1 by construction. Under `set -o pipefail` that
            # status would end the script right here -- losing EXIT=1 and every
            # remaining tier comparison -- turning a real parity failure into a
            # silent partial run.
            diff "$REF_SAM" "$OUT_DIR/$vlabel.$t.sam" | head -20 || true
            EXIT=1
        else
            echo "OK: [$vlabel] tier $t matches $REF_TIER"
            COMPARISONS=$((COMPARISONS + 1))
        fi
    done
done

# NON-VACUITY. With a single testable tier the loop above `continue`s on the
# reference tier and compares nothing -- and this script used to print
# "ALL TIERS PARITY: PASS" on the back of zero comparisons. That is the normal
# case on an avx2-only x86 runner and on every arm64 host (one NEON tier), so
# the reassuring line would have been the usual outcome rather than the
# exception. Report the count, and never say PASS for zero.
if [[ "$EXIT" -ne 0 ]]; then
    exit "$EXIT"
fi
if [[ "$COMPARISONS" -eq 0 ]]; then
    echo "ALL TIERS PARITY: SKIPPED -- only one testable tier on this host (${TIERS[*]}),"
    echo "  so nothing was compared. This is expected on arm64 and on x86 hosts without"
    echo "  AVX-512BW; it is NOT a pass. Run on an avx512bw host for real coverage."
    exit 0
fi
echo "ALL TIERS PARITY: PASS ($COMPARISONS comparison(s) across ${#TIERS[@]} tiers: ${TIERS[*]})"
exit 0
