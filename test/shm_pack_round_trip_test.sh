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

# Build the phix index if any of the index files are missing. No `.0123`: it is
# not built by default (mem pac-fetches from `.pac`) and never staged in shm.
need_index=0
for ext in .amb .ann .bwt.2bit.64 .pac; do
    if [[ ! -s "${PREFIX}${ext}" ]]; then
        need_index=1
        break
    fi
done
if [[ "$need_index" -eq 1 ]]; then
    echo "[setup] Building phix index..."
    "$BIN" index "$PREFIX" > /dev/null 2>&1
fi

# Build a denser oracle (index -u 2) so the inner test can check the `shm -u`
# densify path against an independent from-scratch SA construction.
ORACLE_DIR="$(mktemp -d)"
trap 'rm -rf "$ORACLE_DIR"' EXIT
cp "$PREFIX" "$ORACLE_DIR/phix.fa"
"$BIN" index -u 2 "$ORACLE_DIR/phix.fa" > /dev/null 2>&1

echo "[run] $INNER $PREFIX $ORACLE_DIR/phix.fa"
"$INNER" "$PREFIX" "$ORACLE_DIR/phix.fa"
