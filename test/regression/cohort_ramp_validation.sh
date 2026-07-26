#!/usr/bin/env bash
# test/regression/cohort_ramp_validation.sh
#
# Regression: --cohort-ramp-first / --cohort-ramp-ratio and their environment
# spellings must reject a malformed value instead of silently parsing a prefix.
#
# --cohort-ramp-first takes plain bases, and the C conversion the obvious
# implementation reaches for accepts a prefix: atoll("16M") is 16. That is
# positive, so it passes a `> 0` check and configures a SIXTEEN BYTE first
# slice where the user asked for 16 Mbases -- a millionfold error that produces
# no diagnostic and no crash, only a run that reads the input in dribbles. The
# ratio has the same shape via atof(): "1.5x" becomes 1.5 and "abc" becomes 0.0.
#
# The two spellings are validated the same way but fail differently, on purpose:
#
#   flag  -> hard error, non-zero exit. Nothing is loaded yet, so refusing costs
#            the user nothing and a typo in a script must not run for an hour.
#   env   -> WARNING and fall back to the flag/default value. These are sweep
#            knobs read after the index is loaded; aborting there would throw
#            away minutes of work over a variable the run does not require.
#            Same choice BWA_MEM3_CHUNK_CAP and BWA_MEM3_COHORT_SLICES make.
#
# A silent fallback is what this test rules out in the env case: without the
# WARNING, `BWA_MEM3_COHORT_RAMP_FIRST=16M` and an unset variable are
# indistinguishable from the outside.
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

run_mem() {   # rest = mem args; env comes from the caller. Returns mem's exit code.
    # `|| rc=$?` rather than toggling errexit, matching every other exit-code-
    # capturing helper below. Disabling `set -e` around the call would also
    # suppress failures from any command later added inside the window.
    local rc=0
    "$BWA_MEM3" mem "$@" "$ref" "$reads" >/dev/null 2>"$err" || rc=$?
    return "$rc"
}

# --- the flags reject a malformed value outright ----------------------------
# `=` form so a leading '-' is not taken for another option.
reject_flag() {   # $1 = option name, $2 = value that must be refused
    local opt="$1" val="$2" rc=0
    run_mem "--$opt=$val" || rc=$?
    if [[ "$rc" -eq 0 ]]; then
        echo "FAIL: '--$opt=$val' was accepted; expected a non-zero exit"
        fails=$((fails + 1))
    elif ! grep -q "ERROR: --$opt requires" "$err"; then
        echo "FAIL: '--$opt=$val' exited $rc but without the expected diagnostic"
        fails=$((fails + 1))
    else
        echo "  ok: '--$opt=$val' rejected"
    fi
}

reject_flag cohort-ramp-first 16M     # the millionfold case: bases take no suffix
reject_flag cohort-ramp-first abc     # not a number at all
reject_flag cohort-ramp-first -5      # negative
reject_flag cohort-ramp-ratio 1.5x    # trailing junk after a valid prefix
reject_flag cohort-ramp-ratio abc     # not a number at all

# --- the flags accept a well-formed value -----------------------------------
accept_flag() {   # $1 = option name, $2 = value that must be accepted
    local opt="$1" val="$2" rc=0
    run_mem "--$opt=$val" || rc=$?
    if [[ "$rc" -ne 0 ]]; then
        echo "FAIL: '--$opt=$val' exited $rc; expected success"
        cat "$err" >&2
        fails=$((fails + 1))
    else
        echo "  ok: '--$opt=$val' accepted"
    fi
}

accept_flag cohort-ramp-first 4000000
accept_flag cohort-ramp-first 0        # 0 selects the fractional ramp shape
accept_flag cohort-ramp-ratio 2.5

# --- the env spellings warn and fall back, rather than parsing a prefix -----
env_warns() {   # $1 = variable, $2 = malformed value
    local var="$1" val="$2" rc=0
    env "$var=$val" bash -c '"$0" mem "$1" "$2" >/dev/null 2>"$3"' \
        "$BWA_MEM3" "$ref" "$reads" "$err" || rc=$?
    if [[ "$rc" -ne 0 ]]; then
        echo "FAIL: '$var=$val' exited $rc; a bad sweep variable must not abort the run"
        fails=$((fails + 1))
    elif ! grep -q "WARNING: $var='$val'" "$err"; then
        echo "FAIL: '$var=$val' was accepted silently; a prefix-parsed sweep"
        echo "      variable is indistinguishable from an unset one"
        fails=$((fails + 1))
    else
        echo "  ok: '$var=$val' reported and ignored"
    fi
}

env_quiet() {   # $1 = variable, $2 = well-formed value
    local var="$1" val="$2" rc=0
    env "$var=$val" bash -c '"$0" mem "$1" "$2" >/dev/null 2>"$3"' \
        "$BWA_MEM3" "$ref" "$reads" "$err" || rc=$?
    if [[ "$rc" -ne 0 ]]; then
        echo "FAIL: '$var=$val' exited $rc; expected success"
        cat "$err" >&2
        fails=$((fails + 1))
    elif grep -q "WARNING: $var=" "$err"; then
        echo "FAIL: '$var=$val' is well-formed but was reported as malformed"
        fails=$((fails + 1))
    else
        echo "  ok: '$var=$val' accepted"
    fi
}

env_warns BWA_MEM3_COHORT_RAMP_FIRST 16M
env_warns BWA_MEM3_COHORT_RAMP_FIRST abc
env_warns BWA_MEM3_COHORT_RAMP_FIRST -5
env_warns BWA_MEM3_COHORT_RAMP_RATIO 1.5x
env_warns BWA_MEM3_COHORT_RAMP_RATIO abc
env_quiet BWA_MEM3_COHORT_RAMP_FIRST 4000000
env_quiet BWA_MEM3_COHORT_RAMP_FIRST 0
env_quiet BWA_MEM3_COHORT_RAMP_RATIO 2.5
# An EMPTY value leaves the flag/default alone rather than parsing as 0/0.0.
env_quiet BWA_MEM3_COHORT_RAMP_FIRST ""
env_quiet BWA_MEM3_COHORT_RAMP_RATIO ""

if [[ "$fails" -ne 0 ]]; then
    echo "FAIL: cohort ramp validation ($fails failure(s))"
    exit 1
fi
echo "PASS: cohort ramp values are validated identically by flag and environment"
