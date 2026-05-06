#!/usr/bin/env bash
# test/regression/all_tiers_parity.sh
#
# Regression: validates that running bwa-mem3 with each value of
# BWAMEM3_FORCE_TIER produces byte-identical SAM output.
#
# Must be run on a host that can host every tier — e.g. an AVX-512BW host,
# which can downgrade to avx2/avx/sse42/sse41 via the env var. Lower-tier
# hosts skip the tiers above their host floor (BWAMEM3_FORCE_TIER refuses
# up-tier requests, so we'd just be retesting the host tier).
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

# x86 tiers in low-to-high order. Skip any tier above the host (BWAMEM3_FORCE_TIER
# would refuse the request anyway).
ALL_TIERS=(sse41 sse42 avx avx2 avx512bw neon)

# Fail fast if HOST_TIER didn't parse, or names a tier we don't know about —
# silent fallthrough below would then test every x86 tier on a non-AVX-512 host
# and produce a confusing SIGILL rather than a clear "detection broke" message.
detected_ok=0
for t in "${ALL_TIERS[@]}"; do
    if [[ "$t" == "$HOST_TIER" ]]; then
        detected_ok=1
        break
    fi
done
if [[ "$detected_ok" -ne 1 ]]; then
    echo "FAIL: could not detect host SIMD tier from BWAMEM3_DEBUG_SIMD output" >&2
    echo "  parsed HOST_TIER='$HOST_TIER' (expected one of: ${ALL_TIERS[*]})" >&2
    echo "  raw output was:" >&2
    printf '%s\n' "$HOST_TIER_RAW" >&2
    exit 2
fi
echo "Host tier: $HOST_TIER"

# On arm64 there's only one tier; the dispatcher always picks NEON.
if [[ "$HOST_TIER" == "neon" ]]; then
    TIERS=("neon")
else
    TIERS=()
    for t in "${ALL_TIERS[@]}"; do
        if [[ "$t" == "neon" ]]; then continue; fi
        TIERS+=("$t")
        if [[ "$t" == "$HOST_TIER" ]]; then
            break
        fi
    done
fi

echo "Testing tiers: ${TIERS[*]}"

# Pick a reference SAM (the highest tier — usually the host's actual tier).
# Use ${#TIERS[@]}-1 instead of negative indexing so we work on macOS bash 3.2.
REF_TIER="${TIERS[$((${#TIERS[@]}-1))]}"
REF_SAM="$OUT_DIR/$REF_TIER.sam"

for t in "${TIERS[@]}"; do
    echo ">>> Generating SAM for tier=$t"
    BWAMEM3_FORCE_TIER="$t" "$BWA_MEM3" mem -t 1 \
        "$PARITY_FA" "$PARITY_R1" "$PARITY_R2" \
        > "$OUT_DIR/$t.sam" 2>"$OUT_DIR/$t.log"
    echo "    Produced $(wc -c < "$OUT_DIR/$t.sam") bytes"
done

EXIT=0
for t in "${TIERS[@]}"; do
    if [[ "$t" == "$REF_TIER" ]]; then
        continue
    fi
    if ! diff -q "$REF_SAM" "$OUT_DIR/$t.sam" >/dev/null; then
        echo "FAIL: tier $t differs from $REF_TIER"
        diff "$REF_SAM" "$OUT_DIR/$t.sam" | head -20
        EXIT=1
    else
        echo "OK: tier $t matches $REF_TIER"
    fi
done

if [[ "$EXIT" -eq 0 ]]; then
    echo "ALL TIERS PARITY: PASS"
fi
exit "$EXIT"
