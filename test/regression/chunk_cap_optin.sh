#!/usr/bin/env bash
# test/regression/chunk_cap_optin.sh
#
# Regression: the batch-size cap must stay OFF by default.
#
# bwa (fastmap.c) and bwa-mem2 v2.2.1 (fastmap.cpp) both compute the default
# batch as `opt->chunk_size * opt->n_threads` with no upper bound. bwa-mem3
# briefly capped that default at 256M bases, which silently re-partitioned the
# input for every `-t >= 26` (10M * 26 > 256M). The partition is not cosmetic:
# mem_pestat() infers the insert-size distribution from whatever reads land in a
# batch, and those bounds feed pairing, mate rescue and MAPQ — so a capped run
# is not byte-identical to bwa/bwa-mem2 at the same -t. The regression hid
# because the benchmark aligns at -t 16, where the cap never engaged.
#
# This test pins the contract:
#   * default              -> never capped, at any -t (identical batching to bwa-mem2)
#   * --chunk-cap N        -> caps at N (opt-in)
#   * --fast               -> implies the cap (--fast already is not byte-identical)
#   * BWA_MEM3_CHUNK_CAP   -> wins over both --chunk-cap and --fast's implied cap,
#                             including setting 0 to switch an explicit cap back off;
#                             an empty value is ignored
#   * -K                   -> always wins, is never capped, and does not consult
#                             the env var at all (it takes a different branch)
#
# Inputs:
#   BWA_MEM3 — path to bwa-mem3 binary
#   FIXTURES — directory containing phix.fa and reads.fa (default: test/fixtures)
set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
FIXTURES="${FIXTURES:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../fixtures" && pwd)}"

src_ref="$FIXTURES/phix.fa"
reads="$FIXTURES/reads.fa"
[[ -x "$BWA_MEM3" ]] || { echo "FAIL: binary not executable: $BWA_MEM3" >&2; exit 1; }
[[ -s "$src_ref" ]]  || { echo "FAIL: phix.fa missing: $src_ref" >&2; exit 1; }
[[ -s "$reads" ]]    || { echo "FAIL: reads.fa missing: $reads" >&2; exit 1; }

# Index a private copy so the test never writes into the fixtures tree.
mdir="$(mktemp -d)"; err="$mdir/err.log"
trap 'rm -rf "$mdir"' EXIT
ref="$mdir/phix.fa"
cp "$src_ref" "$ref"
"$BWA_MEM3" index "$ref" >/dev/null 2>&1 || { echo "FAIL: index phix.fa" >&2; exit 1; }

fails=0

# The effective batch size is reported by the reader as `read_chunk: <bases>`.
# Note the args must be passed as separate words: under zsh an unquoted "$1"
# holding "-t 32 -K 100" is ONE argv entry and atoi() silently eats the rest.
batch_size() {
    "$BWA_MEM3" mem "$@" "$ref" "$reads" >/dev/null 2>"$err" \
        || { echo "FAIL: mem exited nonzero (args: $*)" >&2; cat "$err" >&2; exit 1; }
    grep -m1 -oE 'read_chunk: [0-9]+' "$err" | grep -oE '[0-9]+$'
}

expect() {  # $1 = expected batch size, rest = mem args
    local want="$1"; shift
    local got; got="$(batch_size "$@")"
    if [[ "$got" != "$want" ]]; then
        echo "FAIL: 'mem $*' batch=$got, expected $want"
        fails=$((fails + 1))
    else
        echo "  ok: 'mem $*' batch=$got"
    fi
}

# As expect(), but with BWA_MEM3_CHUNK_CAP set for just this invocation.
#
# Written as an explicit helper rather than as a `BWA_MEM3_CHUNK_CAP=N expect ...`
# prefix for two reasons. The prefix form works on a bash function (the value
# reaches the binary and does not leak into later calls -- verified) but POSIX
# leaves that unspecified for functions, so it is not something to rely on. And
# the failure message from expect() prints only the mem args, so a prefixed
# assignment would be invisible in exactly the output someone reads to debug it.
#
# $1 = env value (may be empty), $2 = expected batch size, rest = mem args
expect_env() {
    local envval="$1" want="$2"; shift 2
    local got
    got="$(BWA_MEM3_CHUNK_CAP="$envval" batch_size "$@")"
    if [[ "$got" != "$want" ]]; then
        echo "FAIL: 'BWA_MEM3_CHUNK_CAP=$envval mem $*' batch=$got, expected $want"
        fails=$((fails + 1))
    else
        echo "  ok: 'BWA_MEM3_CHUNK_CAP=$envval mem $*' batch=$got"
    fi
}

# --- the actual regression: default must never cap, at any -t ---------------
expect 320000000 -t 32
expect 640000000 -t 64
expect 160000000 -t 16

# --- opt-in cap -------------------------------------------------------------
expect 256000000 -t 32 --chunk-cap 256000000
expect 160000000 -t 16 --chunk-cap 256000000   # below the cap: untouched
expect 320000000 -t 32 --chunk-cap 0           # 0 disables

# --- --fast implies the cap; explicit flags still win -----------------------
# Both argument orders, and both a non-zero cap and an explicit 0. --fast's cap is
# applied after the option loop and only when --chunk-cap was absent
# (`if (!chunk_cap_set)`), so precedence is order-independent by construction --
# these cases pin that, so a parser change that started consuming --fast after the
# flag (and clobbering it) cannot silently re-enable the implied cap.
expect 256000000 -t 32 --fast
expect  64000000 -t 32 --fast --chunk-cap 64000000
expect  64000000 -t 32 --chunk-cap 64000000 --fast
expect 320000000 -t 32 --fast --chunk-cap 0     # explicit 0 beats --fast's cap
expect 320000000 -t 32 --chunk-cap 0 --fast

# --- BWA_MEM3_CHUNK_CAP wins over everything except -K ----------------------
# The env override exists for cap sweeps, so it has to beat both the CLI flag and
# --fast's implied cap, and be able to switch a cap back OFF.
expect_env 128000000 128000000 -t 32                          # enables a cap the default lacks
expect_env 128000000 128000000 -t 32 --chunk-cap 256000000     # beats an explicit cap
expect_env 128000000 128000000 -t 32 --fast                    # beats --fast's implied cap
expect_env 0         320000000 -t 32 --chunk-cap 256000000     # 0 switches an explicit cap off
expect_env 0         320000000 -t 32 --fast                    # ... and --fast's implied one
expect_env 512000000 320000000 -t 32                           # cap above scaled: untouched

# An EMPTY value is ignored rather than parsed as 0 -- `cap_env && *cap_env` in
# main_mem. Without this case, dropping the `*cap_env` check would go unnoticed:
# atoll("") is 0, which would silently disable an explicit cap.
expect_env "" 256000000 -t 32 --chunk-cap 256000000

# A MALFORMED value is reported and ignored, not parsed as 0. Parsing it as 0
# would silently switch an explicit cap off -- a typo in a sweep variable
# quietly changing how the input is partitioned is the failure mode this whole
# option exists to prevent.
expect_env abc 320000000 -t 32                          # falls back to no cap
expect_env xyz 256000000 -t 32 --chunk-cap 256000000    # falls back to the CLI cap
BWA_MEM3_CHUNK_CAP=abc batch_size -t 32 >/dev/null
grep -q "BWA_MEM3_CHUNK_CAP='abc' is not a non-negative integer" "$err" \
    || { echo "FAIL: a malformed BWA_MEM3_CHUNK_CAP must be reported, not silently ignored"
         fails=$((fails + 1)); }

# --- -K always wins and is never capped ------------------------------------
expect 100000000 -t 32 -K 100000000
expect 100000000 -t 32 --fast -K 100000000
# -K takes the fixed_chunk_size branch, which never reads BWA_MEM3_CHUNK_CAP, so
# even an env cap far below -K must not clip it.
expect_env 1000000 100000000 -t 32 -K 100000000

# --- the warning fires only when the cap actually engages -------------------
# Three cases, because the warning is gated on `scaled > cap`, not on "a cap is
# configured": an implementation that warns whenever --chunk-cap is present
# would pass the first and third checks and still be wrong.
batch_size -t 32 --chunk-cap 256000000 >/dev/null   # 320M > 256M: clipped
grep -q 'chunk cap engaged' "$err" || { echo "FAIL: expected 'chunk cap engaged' warning"; fails=$((fails + 1)); }
batch_size -t 16 --chunk-cap 256000000 >/dev/null   # 160M < 256M: configured but inert
if grep -q 'chunk cap engaged' "$err"; then
    echo "FAIL: an unengaged chunk cap must not report a warning"
    fails=$((fails + 1))
fi
batch_size -t 32 >/dev/null                          # no cap at all
if grep -q 'chunk cap engaged' "$err"; then
    echo "FAIL: default run must not report a chunk cap"
    fails=$((fails + 1))
fi

# --- an unparseable --chunk-cap is rejected, not silently read as 0 ---------
# atoll() maps every unparseable value to 0, and 0 means "no cap" -- so without
# validation `--chunk-cap 100M` would look accepted and quietly leave batching
# uncapped. Uses `=` form so a leading '-' is not taken for another option.
reject() {  # $1 = the --chunk-cap value that must be refused
    # `|| rc=$?` rather than toggling errexit off around the call: a `set +e`
    # window silently stops checking every command inside it, including any
    # later added between the call and `rc=$?`.
    local rc=0
    "$BWA_MEM3" mem "--chunk-cap=$1" "$ref" "$reads" >/dev/null 2>"$err" || rc=$?
    if [[ "$rc" -eq 0 ]]; then
        echo "FAIL: '--chunk-cap=$1' was accepted; expected a non-zero exit"
        fails=$((fails + 1))
    elif ! grep -q 'ERROR: --chunk-cap requires a non-negative integer' "$err"; then
        echo "FAIL: '--chunk-cap=$1' exited $rc but without the expected diagnostic"
        fails=$((fails + 1))
    else
        echo "  ok: '--chunk-cap=$1' rejected"
    fi
}
reject abc      # not a number at all
reject 100M     # trailing suffix: the option takes plain bases, not 100M
reject -5       # negative

if [[ "$fails" -ne 0 ]]; then
    echo "FAIL: chunk-cap opt-in regression ($fails failure(s))"
    exit 1
fi
echo "PASS: chunk-cap is opt-in; default batching matches bwa/bwa-mem2 at every -t"
