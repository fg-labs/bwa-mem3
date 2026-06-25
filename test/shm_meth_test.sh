#!/usr/bin/env bash
# Regression: `bwa-mem3 shm --meth` resolves the `.meth` seed-index prefix the
# same way `bwa-mem3 mem --meth` does, so an end user passes the same plain
# FASTA path to all three commands (`index --meth`, `shm --meth`, `mem --meth`)
# without juggling the `.meth` suffix by hand.
#
# D3 (dual index): `index --meth` builds BOTH the original-alphabet index at the
# bare prefix (`ref.fa.{ann,amb,pac,bwt.2bit.64}`) AND the converted seed FM-index at
# `ref.fa.meth.*`. `shm --meth <plain>` stages the seed segment under the
# `ref.fa.meth` basename; `mem --meth <plain>` attaches it from shm. Because the
# bare index now exists too, plain `shm <plain>` also succeeds (stages the
# original index) — the dual index removes the old meth-only failure mode.

set -euo pipefail
cd "$(dirname "$0")/.."

BIN=${BIN:-./bwa-mem3}

if [[ ! -x "$BIN" ]]; then
    echo "FAIL: $BIN not built. Run 'make -j4' first." >&2
    exit 2
fi

# Isolated dir so the staged segments are unambiguous.
WORKDIR="$(mktemp -d -t shm_meth.XXXXXX)"
cp test/fixtures/phix.fa "$WORKDIR/ref.fa"
PLAIN_PREFIX="$WORKDIR/ref.fa"
METH_PREFIX="$WORKDIR/ref.fa.meth"

# Drop any stale registry from prior runs and on every exit.
"$BIN" shm -d >/dev/null 2>&1 || true
ERR="$(mktemp)"
trap '"$BIN" shm -d >/dev/null 2>&1 || true; rm -rf "$WORKDIR" "$ERR"' EXIT

echo "[setup] bwa-mem3 index --meth $PLAIN_PREFIX"
"$BIN" index --meth "$PLAIN_PREFIX" >/dev/null 2>&1

# D3 dual-index invariant: `index --meth` writes the bare original index AND the
# .meth seed index. Neither carries the unpacked `.0123` by default — `mem`
# pac-fetches the original reference from `.pac`, and `mem --meth` never extends
# against the seed at all (saves ~6.4 GB original + ~13 GB seed on hg38).
for ext in .ann .amb .pac .bwt.2bit.64; do
    if [[ ! -e "${PLAIN_PREFIX}${ext}" ]]; then
        echo "FAIL: dual index missing bare original index ${PLAIN_PREFIX}${ext}" >&2
        exit 1
    fi
done
for ext in .ann .amb .pac .bwt.2bit.64; do
    if [[ ! -e "${METH_PREFIX}${ext}" ]]; then
        echo "FAIL: dual index missing seed index ${METH_PREFIX}${ext}" >&2
        exit 1
    fi
done
if [[ -e "${PLAIN_PREFIX}.0123" ]]; then
    echo "FAIL: original index ${PLAIN_PREFIX}.0123 was built; it must not be by default (mem pac-fetches from .pac)" >&2
    exit 1
fi
if [[ -e "${METH_PREFIX}.0123" ]]; then
    echo "FAIL: seed index ${METH_PREFIX}.0123 was built; it must not be (never read in --meth)" >&2
    exit 1
fi

# --- A: `shm --meth <plain>` stages the seed under the .meth basename ----
echo "[A] bwa-mem3 shm --meth $PLAIN_PREFIX"
"$BIN" shm --meth "$PLAIN_PREFIX" >/dev/null 2>&1 \
    || { echo "FAIL A: shm --meth <plain> exited non-zero" >&2; exit 1; }

LIST="$("$BIN" shm -l 2>&1)"
COUNT_A="$(echo "$LIST" | awk '{print $1}' | grep -xc 'ref.fa.meth' || true)"
if [[ "$COUNT_A" -ne 1 ]]; then
    echo "FAIL A: expected exactly one 'ref.fa.meth' entry, got $COUNT_A" >&2
    echo "----- shm -l -----" >&2; echo "$LIST" >&2
    exit 1
fi
"$BIN" shm -d >/dev/null 2>&1

# --- B: mem --meth attaches from the shm segment staged by shm --meth ---
echo "[B] shm --meth + mem --meth ⇒ 'attached from shm'"
"$BIN" shm --meth "$PLAIN_PREFIX" >/dev/null 2>&1
: > "$ERR"
"$BIN" mem --meth "$PLAIN_PREFIX" test/fixtures/reads.fa >/dev/null 2>"$ERR" \
    || { echo "FAIL B: mem --meth exited non-zero" >&2; cat "$ERR" >&2; exit 1; }
if ! grep -q "attached from shm" "$ERR"; then
    echo "FAIL B: mem --meth did not attach from shm (registry mismatch?)" >&2
    tail -40 "$ERR" >&2
    exit 1
fi
"$BIN" shm -d >/dev/null 2>&1

# --- C: `shm --meth <.meth-suffixed-prefix>` is idempotent (no double-append) ---
echo "[C] bwa-mem3 shm --meth <prefix>.meth (no double-append)"
"$BIN" shm --meth "$METH_PREFIX" >/dev/null 2>&1 \
    || { echo "FAIL C: shm --meth on already-.meth prefix exited non-zero" >&2; exit 1; }
LIST="$("$BIN" shm -l 2>&1)"
COUNT_C="$(echo "$LIST" | awk '{print $1}' | grep -xc 'ref.fa.meth' || true)"
if [[ "$COUNT_C" -ne 1 ]]; then
    echo "FAIL C: expected exactly one 'ref.fa.meth' entry, got $COUNT_C" >&2
    echo "----- shm -l -----" >&2; echo "$LIST" >&2
    exit 1
fi
if echo "$LIST" | awk '{print $1}' | grep -qx 'ref.fa.meth.meth'; then
    echo "FAIL C: double-appended 'ref.fa.meth.meth' entry present" >&2
    exit 1
fi
"$BIN" shm -d >/dev/null 2>&1

# --- D: plain `shm <plain>` (no --meth) now SUCCEEDS on the dual index ----
# The dual index writes the bare original index, so the old "meth-only directory"
# failure mode is gone: `shm <plain>` stages the original index.
echo "[D] bwa-mem3 shm <plain>  (dual index ⇒ plain staging works)"
"$BIN" shm "$PLAIN_PREFIX" >/dev/null 2>&1 \
    || { echo "FAIL D: shm <plain> on a dual index unexpectedly failed" >&2; exit 1; }
LIST="$("$BIN" shm -l 2>&1)"
COUNT_D="$(echo "$LIST" | awk '{print $1}' | grep -xc 'ref.fa' || true)"
if [[ "$COUNT_D" -ne 1 ]]; then
    echo "FAIL D: expected exactly one 'ref.fa' entry, got $COUNT_D" >&2
    echo "----- shm -l -----" >&2; echo "$LIST" >&2
    exit 1
fi
"$BIN" shm -d >/dev/null 2>&1

# --- E: extra positional arguments after idxbase are rejected -----------
echo "[E] bwa-mem3 shm --meth <prefix> <typo>  (extra positional rejected)"
: > "$ERR"
if "$BIN" shm --meth "$PLAIN_PREFIX" stray-arg >/dev/null 2>"$ERR"; then
    echo "FAIL E: shm --meth <prefix> <stray-arg> unexpectedly succeeded" >&2
    cat "$ERR" >&2
    exit 1
fi
if ! grep -qiE "positional|too many|unexpected" "$ERR"; then
    echo "FAIL E: error message did not flag the extra positional arg" >&2
    echo "----- stderr -----" >&2; cat "$ERR" >&2
    exit 1
fi

echo "OK"
