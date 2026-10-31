#!/usr/bin/env bash
# test/regression/all_tiers_parity.sh
#
# Regression: validates that one binary, on one host, produces byte-identical
# SAM whichever of its usable kernels BWAMEM3_FORCE_TIER selects.
#
# SCOPE: avx2 and avx512bw only (plus neon on arm64, where there is one tier).
# sse41/sse42/avx are excluded because there is no batched kswv kernel below
# AVX2 and forcing them kills the aligner — see SIMD_SWEEP_TIERS in
# test/framework/simd_tiers.sh for the full reasoning and the measurement
# behind it.
#
# This script is therefore NOT the cross-build parity check, and should not be
# mistaken for one. "Does a low-baseline build agree with an avx2 build?" is
# answered by CI building at sse41/avx2/arm64 baselines and running
# chr22_parity.sh on each, which asserts byte-identity against bwa itself — a
# stronger oracle than self-comparison, and transitive across the three builds.
#
# Real coverage needs an avx512bw host: on an avx2-only x86 runner, or on arm64,
# there is exactly one sweepable tier and the script reports SKIP: rather than
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

# The sweepable tiers (avx2/avx512bw on x86, neon on arm64) at or below the
# host's. The helper's comment says why sse41/sse42/avx are excluded. The
# supported tiers (avx2/avx512bw and arm64) are each run through
# chr22_parity.sh, asserting byte-identity against bwa itself, so the builds
# agree transitively -- a stronger oracle than this script's self-comparison,
# and where cross-build parity actually lives. Caveat on scope: chr22_parity.sh
# strips MQ:i, HN:i, and XS:i before its diff against bwa (see
# proto-neon-kswv.yml), so that transitive cross-build claim does NOT cover
# those three emitted tags. This script compares the FULL, unstripped SAM
# record, so it is what establishes equality of MQ:i/HN:i/XS:i across a single
# build's reachable tiers.
#
# What remains here is the reachable-tier question: one binary, on one host,
# must produce the same SAM whichever of its usable kernels it dispatches.
# shellcheck source=test/framework/simd_tiers.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../framework" && pwd)/simd_tiers.sh"
simd_sweep_tiers "$BWA_MEM3" || exit $?
HOST_TIER="$SIMD_HOST_TIER"
# ${a[@]+...}: an empty array is "unbound" under set -u in bash 3.2.
TIERS=(${SIMD_TIERS[@]+"${SIMD_TIERS[@]}"})
echo "Host tier: $HOST_TIER"

# A host below the batched floor (an sse41/sse42/avx baseline build) has nothing
# to sweep. Say so plainly instead of reporting an empty pass.
if [[ "${#TIERS[@]}" -eq 0 ]]; then
    echo "SKIP: all_tiers_parity -- host tier '$HOST_TIER' is below the batched-kswv"
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
        rc=0
        BWAMEM3_FORCE_TIER="$t" "$BWA_MEM3" mem -t 1 ${vargs[@]+"${vargs[@]}"} \
            "$PARITY_FA" "$PARITY_R1" "$PARITY_R2" \
            > "$OUT_DIR/$vlabel.$t.sam" 2> "$OUT_DIR/$vlabel.$t.log" || rc=$?
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
    echo "SKIP: all_tiers_parity -- only one testable tier on this host (${TIERS[*]}),"
    echo "  so nothing was compared. This is expected on arm64 and on x86 hosts without"
    echo "  AVX-512BW; it is NOT a pass. Run on an avx512bw host for real coverage."
    exit 0
fi
echo "PASS: all_tiers_parity ($COMPARISONS comparison(s) across ${#TIERS[@]} tiers: ${TIERS[*]})"
exit 0
