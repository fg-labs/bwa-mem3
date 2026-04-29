#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

BIN=${BIN:-./bwa-mem3}
INNER=${INNER:-./shm_pack_round_trip_test}
PREFIX=test/fixtures/phix.fa

if [[ ! -x "$BIN" ]]; then
    echo "FAIL: $BIN not built. Run 'make -j4' first." >&2
    exit 2
fi
if [[ ! -x "$INNER" ]]; then
    echo "FAIL: $INNER not built. Run 'make -j4 shm_pack_round_trip_test' first." >&2
    exit 2
fi

# Build the phix index if any of the index files are missing.
need_index=0
for ext in .0123 .amb .ann .bwt.2bit.64 .pac; do
    if [[ ! -s "${PREFIX}${ext}" ]]; then
        need_index=1; break
    fi
done
if [[ "$need_index" -eq 1 ]]; then
    echo "[setup] Building phix index..."
    "$BIN" index "$PREFIX" >/dev/null 2>&1
fi

echo "[run] $INNER $PREFIX"
"$INNER" "$PREFIX"
