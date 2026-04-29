#!/usr/bin/env bash
# Regression: `bwa-mem3 shm --meth` resolves the c2t-suffixed prefix the
# same way `bwa-mem3 mem --meth` does, so an end user passes the same plain
# FASTA path to all three commands (`index --meth`, `shm --meth`,
# `mem --meth`) without juggling the `.bwameth.c2t` suffix by hand.
#
# Pre-fix, `bwa-mem3 shm <plain>` opens `bns_restore(<plain>)` and fails on
# `<plain>.ann` when only the meth artifacts are on disk; staging instead
# required typing `bwa-mem3 shm <plain>.bwameth.c2t`. This test pins the
# new flag-driven behavior in place.

set -euo pipefail
cd "$(dirname "$0")/.."

BIN=${BIN:-./bwa-mem3}

if [[ ! -x "$BIN" ]]; then
    echo "FAIL: $BIN not built. Run 'make -j4' first." >&2
    exit 2
fi

# Isolated dir so the meth-only artifacts are the *only* index in scope.
# test/fixtures/phix.fa already has a plain index from other tests, which
# would mask the negative case below.
WORKDIR="$(mktemp -d -t shm_meth.XXXXXX)"
cp test/fixtures/phix.fa "$WORKDIR/ref.fa"
PLAIN_PREFIX="$WORKDIR/ref.fa"
C2T_PREFIX="$WORKDIR/ref.fa.bwameth.c2t"

# Drop any stale registry from prior runs and on every exit. Allocate the
# stderr-capture file once and truncate between cases — avoids mktemp churn
# and ensures the EXIT trap (which expands $ERR at fire time) cleans up the
# single live path rather than leaking earlier ones.
"$BIN" shm -d >/dev/null 2>&1 || true
ERR="$(mktemp)"
trap '"$BIN" shm -d >/dev/null 2>&1 || true; rm -rf "$WORKDIR" "$ERR"' EXIT

# Build only the meth index. The plain ref.fa.{0123,ann,...} are NOT built,
# so this is the exact scenario where the pre-fix `bwa-mem3 shm ref.fa` UX
# would have failed with `fail to open '...ref.fa.ann'`.
echo "[setup] bwa-mem3 index --meth $PLAIN_PREFIX"
"$BIN" index --meth "$PLAIN_PREFIX" >/dev/null 2>&1

# Sanity: the plain index files must be absent (so test cases C and D below
# actually exercise the wrong-flag path).
for ext in .0123 .ann .amb .pac .bwt.2bit.64; do
    if [[ -e "${PLAIN_PREFIX}${ext}" ]]; then
        echo "FAIL: precondition violated — ${PLAIN_PREFIX}${ext} exists" >&2
        exit 1
    fi
done

# --- A: `shm --meth <plain>` stages under the c2t basename --------------
echo "[A] bwa-mem3 shm --meth $PLAIN_PREFIX"
"$BIN" shm --meth "$PLAIN_PREFIX" >/dev/null 2>&1 \
    || { echo "FAIL A: shm --meth <plain> exited non-zero" >&2; exit 1; }

LIST="$("$BIN" shm -l 2>&1)"
COUNT_A="$(echo "$LIST" | awk '{print $1}' | grep -xc 'ref.fa.bwameth.c2t' || true)"
if [[ "$COUNT_A" -ne 1 ]]; then
    echo "FAIL A: expected exactly one 'ref.fa.bwameth.c2t' entry, got $COUNT_A" >&2
    echo "----- shm -l -----" >&2
    echo "$LIST" >&2
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

# --- C: `shm --meth <c2t-suffixed-prefix>` is idempotent ----------------
echo "[C] bwa-mem3 shm --meth <prefix>.bwameth.c2t (no double-append)"
"$BIN" shm --meth "$C2T_PREFIX" >/dev/null 2>&1 \
    || { echo "FAIL C: shm --meth on already-c2t prefix exited non-zero" >&2; exit 1; }
LIST="$("$BIN" shm -l 2>&1)"
COUNT_C="$(echo "$LIST" | awk '{print $1}' | grep -xc 'ref.fa.bwameth.c2t' || true)"
if [[ "$COUNT_C" -ne 1 ]]; then
    echo "FAIL C: expected exactly one 'ref.fa.bwameth.c2t' entry, got $COUNT_C" >&2
    echo "----- shm -l -----" >&2
    echo "$LIST" >&2
    exit 1
fi
if echo "$LIST" | awk '{print $1}' | grep -qx 'ref.fa.bwameth.c2t.bwameth.c2t'; then
    echo "FAIL C: double-appended 'ref.fa.bwameth.c2t.bwameth.c2t' entry present" >&2
    exit 1
fi
"$BIN" shm -d >/dev/null 2>&1

# --- D: `shm <plain>` (no --meth) errors with a helpful hint ------------
echo "[D] bwa-mem3 shm <plain>  (only meth artifacts on disk → hint about --meth)"
: > "$ERR"
if "$BIN" shm "$PLAIN_PREFIX" >/dev/null 2>"$ERR"; then
    echo "FAIL D: shm <plain> on a meth-only directory unexpectedly succeeded" >&2
    cat "$ERR" >&2
    exit 1
fi
if ! grep -q -- "--meth" "$ERR"; then
    echo "FAIL D: error message did not mention '--meth' as a hint" >&2
    echo "----- stderr -----" >&2
    cat "$ERR" >&2
    exit 1
fi

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
    echo "----- stderr -----" >&2
    cat "$ERR" >&2
    exit 1
fi

echo "OK"
