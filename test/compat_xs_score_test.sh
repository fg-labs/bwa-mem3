#!/usr/bin/env bash
# Synthetic compatibility routing and secondary-score check.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BWAMEM3="${BWAMEM3:-$ROOT/bwa-mem3}"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/bwamem3-compat-xs.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT

cp "$ROOT/test/fixtures/compat_xs_ref.fa" "$TMP/ref.fa"
"$BWAMEM3" index "$TMP/ref.fa" > "$TMP/index.log" 2>&1

"$BWAMEM3" mem -v 4 --compat=bwa-mem2 \
    "$TMP/ref.fa" "$ROOT/test/fixtures/compat_xs_read.fq" \
    > "$TMP/compat.sam" 2> "$TMP/compat.log"
"$BWAMEM3" mem -v 4 \
    "$TMP/ref.fa" "$ROOT/test/fixtures/compat_xs_read.fq" \
    > "$TMP/native.sam" 2> "$TMP/native.log"

python3 - "$TMP/compat.sam" "$TMP/compat.log" "$TMP/native.log" << 'PY'
import pathlib
import sys

sam, compat_log, native_log = map(lambda p: pathlib.Path(p).read_text(), sys.argv[1:])
records = [line.split('\t') for line in sam.splitlines() if not line.startswith('@')]
assert len(records) == 1, f"expected one SAM record, got {len(records)}"
tags = records[0][11:]
assert 'AS:i:75' in tags, f"unexpected primary score: {tags}"
assert 'XS:i:70' in tags, f"unexpected secondary score: {tags}"
assert 'Ungapped fast-path' not in compat_log, 'compatibility path took ungapped shortcut'
assert 'Ungapped fast-path' in native_log, 'fixture did not exercise native shortcut'
PY
