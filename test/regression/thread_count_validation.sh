#!/usr/bin/env bash
# test/regression/thread_count_validation.sh
#
# Regression: `-t` must reject malformed values instead of silently
# substituting a working thread count, and must clamp (not reject) an
# in-range-but-too-large value.
#
# `-t` was parsed with atoi(), which maps every unparseable string to 0 and
# ignores trailing garbage: `-t nope` silently became 1 thread, and
# `-t 300oops` silently became 300 threads, with no diagnostic either way.
# Separately, per-thread profiling arrays (tprof[][tid]) are sized to
# MAX_THREADS (256) and indexed by tid up to n_threads-1, so a valid-looking
# `-t` above 256 previously wrote out of bounds; that is now clamped to 256
# with a warning rather than rejected, since it is a valid (if oversized)
# request rather than a malformed one.
#
# This test pins the contract:
#   * `-t` with a non-numeric or partially-numeric value -> hard error
#   * `-t` with a value <= 0 -> floors to 1 (unchanged, not an error)
#   * `-t` with a value > 256 -> clamped to 256, with a warning, run succeeds
#   * `-t` with a value in 1..256 -> used as given
#
# Inputs:
#   BWA_MEM3 — path to bwa-mem3 binary
#   FIXTURES — directory containing phix.fa and reads.fa (default: test/fixtures)
set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
FIXTURES="${FIXTURES:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../fixtures" && pwd)}"

src_ref="$FIXTURES/phix.fa"
reads="$FIXTURES/reads.fa"
[[ -x "$BWA_MEM3" ]] || {
    echo "FAIL: binary not executable: $BWA_MEM3" >&2
    exit 1
}
[[ -s "$src_ref" ]] || {
    echo "FAIL: phix.fa missing: $src_ref" >&2
    exit 1
}
[[ -s "$reads" ]] || {
    echo "FAIL: reads.fa missing: $reads" >&2
    exit 1
}

# Index a private copy so the test never writes into the fixtures tree.
mdir="$(mktemp -d)"
err="$mdir/err.log"
trap 'rm -rf "$mdir"' EXIT
ref="$mdir/phix.fa"
cp "$src_ref" "$ref"
"$BWA_MEM3" index "$ref" > /dev/null 2>&1 || {
    echo "FAIL: index phix.fa" >&2
    exit 1
}

fails=0
t_err='ERROR: -t requires an integer thread count'

# $1 = value for -t, $2 = expected diagnostic substring.
reject() {
    local val="$1" want="$2" rc=0
    "$BWA_MEM3" mem -t "$val" "$ref" "$reads" > /dev/null 2> "$err" || rc=$?
    if [[ "$rc" -eq 0 ]]; then
        echo "FAIL: '-t $val' was accepted; expected a non-zero exit"
        fails=$((fails + 1))
    elif ! grep -q "$want" "$err"; then
        echo "FAIL: '-t $val' exited $rc but without the expected diagnostic ('$want')"
        cat "$err" >&2
        fails=$((fails + 1))
    else
        echo "  ok: '-t $val' rejected"
    fi
}

# $1 = value for -t that must be ACCEPTED (the run itself must succeed).
# $2 (optional) = substring that must appear on stderr (e.g. the clamp warning).
accept() {
    local val="$1" want="${2:-}"
    if ! "$BWA_MEM3" mem -t "$val" "$ref" "$reads" > /dev/null 2> "$err"; then
        echo "FAIL: '-t $val' was rejected; expected it to run"
        cat "$err" >&2
        fails=$((fails + 1))
        return
    fi
    if [[ -n "$want" ]] && ! grep -q "$want" "$err"; then
        echo "FAIL: '-t $val' ran but stderr lacked the expected substring ('$want')"
        cat "$err" >&2
        fails=$((fails + 1))
    else
        echo "  ok: '-t $val' accepted"
    fi
}

# --- malformed values: not a complete integer ---------------------------------
reject "nope" "$t_err"                 # not a number at all -> would have meant 1 thread
reject "300oops" "$t_err"              # trailing junk: partial parse must not pass
reject "" "$t_err"                     # empty
reject "99999999999999999999" "$t_err" # overflows a 64-bit parse (ERANGE)

# --- valid but <= 0: floors to 1, unchanged, not an error ----------------------
accept "0"
accept "-5"

# --- valid and in range ---------------------------------------------------------
accept "1"
accept "4"
accept "256" # MAX_THREADS, the high edge -- must be accepted, not clamped

# --- valid but above MAX_THREADS: clamped with a warning, run still succeeds ---
accept "257" "clamping to 256"
accept "999999" "clamping to 256"

if [[ "$fails" -ne 0 ]]; then
    echo "FAIL: -t option validation ($fails failure(s))"
    exit 1
fi
echo "PASS: -t rejects malformed values and clamps out-of-range values"
