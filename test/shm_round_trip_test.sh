#!/usr/bin/env bash
# End-to-end byte-identical regression: bwa-mem3 mem with a staged shm
# segment vs. a disk-load run on the same phix index. Pass criterion is
# (a) staged run shows the attach log line, (b) disk run does not, and
# (c) `diff` between the two SAM outputs produces no output.

set -euo pipefail
cd "$(dirname "$0")/.."

BIN=${BIN:-./bwa-mem3}
PREFIX=test/fixtures/phix.fa
READS=test/fixtures/reads.fa

if [[ ! -x "$BIN" ]]; then
    echo "FAIL: $BIN not built. Run 'make -j4' first." >&2
    exit 2
fi

# Clean up any stale registry from a prior run (and on exit).
"$BIN" shm -d >/dev/null 2>&1 || true
trap '"$BIN" shm -d >/dev/null 2>&1 || true' EXIT

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

if [[ ! -f "$READS" ]]; then
    echo "FAIL: reads fixture missing: $READS" >&2
    exit 2
fi

DISK_SAM=$(mktemp -t shm_disk.XXXXXX).sam
SHM_SAM=$(mktemp -t shm_attached.XXXXXX).sam
DISK_ERR=$(mktemp -t shm_disk.XXXXXX).err
SHM_ERR=$(mktemp -t shm_attached.XXXXXX).err
trap 'rm -f "$DISK_SAM" "$SHM_SAM" "$DISK_ERR" "$SHM_ERR"; "$BIN" shm -d >/dev/null 2>&1 || true' EXIT

echo "[run] disk baseline"
"$BIN" mem "$PREFIX" "$READS" > "$DISK_SAM" 2> "$DISK_ERR"

echo "[run] stage and re-align"
"$BIN" shm "$PREFIX"
"$BIN" mem "$PREFIX" "$READS" > "$SHM_SAM" 2> "$SHM_ERR"

# (a) Attach log line present in shm run.
if ! grep -q "attached from shm" "$SHM_ERR"; then
    echo "FAIL: shm run did not show 'attached from shm' in stderr" >&2
    echo "----- shm stderr (last 40 lines) -----" >&2
    tail -40 "$SHM_ERR" >&2
    exit 1
fi

# (b) Attach log line absent in disk run.
if grep -q "attached from shm" "$DISK_ERR"; then
    echo "FAIL: disk run unexpectedly showed 'attached from shm' in stderr" >&2
    echo "----- disk stderr (last 40 lines) -----" >&2
    tail -40 "$DISK_ERR" >&2
    exit 1
fi

# (c) Byte-identical SAM.
if ! diff -q "$DISK_SAM" "$SHM_SAM" >/dev/null; then
    echo "FAIL: SAM differs between disk and shm runs" >&2
    diff -u "$DISK_SAM" "$SHM_SAM" | head -40 >&2
    exit 1
fi

# (d) After dropping the segment, mem reverts to the disk path.
"$BIN" shm -d
DROP_ERR=$(mktemp -t shm_drop.XXXXXX).err
"$BIN" mem "$PREFIX" "$READS" >/dev/null 2>"$DROP_ERR"
if grep -q "attached from shm" "$DROP_ERR"; then
    echo "FAIL: mem still shows 'attached from shm' after drop" >&2
    rm -f "$DROP_ERR"
    exit 1
fi
rm -f "$DROP_ERR"

echo "OK"
