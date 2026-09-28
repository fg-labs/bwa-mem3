#!/usr/bin/env bash
# Synthetic XA cutoff regression at secondary scores 32 and 31 of primary 40.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BWAMEM3="${BWAMEM3:-$ROOT/bwa-mem3}"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/bwamem3-compat-xa.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT

cp "$ROOT/test/fixtures/compat_xa_ref.fa" "$TMP/ref.fa"
"$BWAMEM3" index "$TMP/ref.fa" > "$TMP/index.log" 2>&1
"$BWAMEM3" mem --compat=bwa-mem2 \
    "$TMP/ref.fa" "$ROOT/test/fixtures/compat_xa_reads.fq" \
    > "$TMP/out.sam" 2> "$TMP/mem.log"

python3 - "$TMP/out.sam" << 'PY'
import pathlib
import sys

records = {}
for line in pathlib.Path(sys.argv[1]).read_text().splitlines():
    if not line.startswith('@'):
        fields = line.split('\t')
        records[fields[0]] = fields[11:]

assert set(records) == {'boundary32', 'boundary31'}, records
assert 'AS:i:40' in records['boundary32'], records['boundary32']
assert 'XS:i:32' in records['boundary32'], records['boundary32']
assert 'XA:Z:boundary32_alt,+1,20M2D20M,2;' in records['boundary32'], records['boundary32']
assert 'AS:i:40' in records['boundary31'], records['boundary31']
assert 'XS:i:31' in records['boundary31'], records['boundary31']
assert not any(tag.startswith('XA:Z:') for tag in records['boundary31']), records['boundary31']
PY
