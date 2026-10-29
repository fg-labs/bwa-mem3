#!/usr/bin/env bash
# test/regression/mem_cli_selftest.sh
#
# Self-test for the shared shell helpers the mem CLI scripts are built on
# (test/framework/mem_cli.sh and test/framework/simd_tiers.sh). A helper that
# stopped failing would make every script on it pass without checking
# anything, so each one is driven here with inputs whose answer is known:
#
#   simd_sweep_tiers   stub binaries reporting each dispatcher tier, and one
#                      reporting garbage: avx512bw -> avx2 avx512bw, avx2 ->
#                      avx2, sse41 -> none, neon -> neon, garbage -> status 2;
#                      a BWAMEM3_FORCE_TIER in the caller's environment is
#                      ignored.
#   mem_cli_init       MEM_CLI_TIERS that is whitespace only, names an
#                      unknown tier or names a tier twice ends the script
#                      with FAIL; an empty one counts as unset; a host below
#                      the floor reports SKIP and exits 0.
#   mem_cli_verdict    zero checks is a failure, a failed check is a failure.
#   mem_cli_expect     a mismatch is counted as a failure.
#   mem_cli_same_file  differing, missing and header-only files all fail;
#                      equal files that differ only in @PG pass.
#   mem_cli_primary    a read with two primary records yields nothing.
#   mem_cli_tag        a present tag's value, nothing for an absent one.
#   mem_cli_seq        an out-of-range slice, or an unknown contig even for
#                      an empty slice, returns 1 with a FAIL line.
#
# Needs no bwa-mem3 build. Env: TMPDIR (optional) for the scratch directory.

set -euo pipefail

framework="$(cd "$(dirname "${BASH_SOURCE[0]}")/../framework" && pwd)"
work="$(mktemp -d "${TMPDIR:-/tmp}/mem_cli_selftest.XXXXXX")"
trap 'rm -rf "$work"' EXIT

fails=0
check() { # check <label> <expected> <actual>
    if [[ $2 == "$3" ]]; then
        echo "ok: $1"
    else
        echo "FAIL: $1: expected '$2', got '$3'" >&2
        fails=$((fails + 1))
    fi
}

# A stub bwa-mem3 whose dispatcher reports <tier>, and which records the
# BWAMEM3_FORCE_TIER it saw.
stub() {
    local path="$work/stub-$1"
    # The stub's own ${BWAMEM3_FORCE_TIER:-} is expanded when the stub runs.
    # shellcheck disable=SC2016
    printf '#!/bin/sh\necho "[M::main] SIMD tier: %s (build baseline: x)" >&2\necho "${BWAMEM3_FORCE_TIER:-}" > "%s.forced"\nexit 1\n' \
        "$1" "$path" > "$path"
    chmod +x "$path"
    printf '%s\n' "$path"
}

# --- simd_sweep_tiers ------------------------------------------------------
for case in "avx512bw:avx2 avx512bw" "avx2:avx2" "sse41:" "neon:neon"; do
    host="${case%%:*}"
    want="${case#*:}"
    bin="$(stub "$host")"
    got="$(
        # shellcheck source=test/framework/simd_tiers.sh
        source "$framework/simd_tiers.sh"
        BWAMEM3_FORCE_TIER=avx2 simd_sweep_tiers "$bin"
        echo "$SIMD_HOST_TIER:${SIMD_TIERS[*]-}"
    )"
    check "simd_sweep_tiers on a $host host" "$host:$want" "$got"
    check "simd_sweep_tiers ignores the caller's BWAMEM3_FORCE_TIER ($host)" "" "$(cat "$bin.forced")"
done
bin="$(stub garbage)"
rc=0
(
    # shellcheck source=test/framework/simd_tiers.sh
    source "$framework/simd_tiers.sh"
    simd_sweep_tiers "$bin"
) 2> /dev/null || rc=$?
check "simd_sweep_tiers on an unparseable tier returns 2" 2 "$rc"

# --- mem_cli_init and mem_cli_verdict (each in its own script) -------------
# run_script <name> <env...> -- <body>: runs a script that sources mem_cli.sh;
# prints "<status>|<last stdout+stderr line>".
run_script() {
    local name="$1" body
    shift
    local envs=()
    while [[ $1 != -- ]]; do
        envs+=("$1")
        shift
    done
    body="$2"
    printf '#!/usr/bin/env bash\nset -euo pipefail\nsource "%s/mem_cli.sh"\n%s\n' "$framework" "$body" > "$work/$name.sh"
    local rc=0
    env ${envs[@]+"${envs[@]}"} bash "$work/$name.sh" > "$work/$name.out" 2>&1 || rc=$?
    printf '%s|%s\n' "$rc" "$(tail -n 1 "$work/$name.out")"
}
neon="$(stub neon)"
sse41="$(stub sse41)"

out="$(run_script blank BWA_MEM3="$neon" "MEM_CLI_TIERS= " -- 'mem_cli_init blank')"
check "mem_cli_init: blank MEM_CLI_TIERS fails" "1|FAIL: blank: MEM_CLI_TIERS=' ' names no tier" "$out"
# The generated script expands MEM_CLI_TIER_LIST itself.
# shellcheck disable=SC2016
out="$(run_script empty BWA_MEM3="$neon" "MEM_CLI_TIERS=" -- 'mem_cli_init empty; echo "tiers: ${MEM_CLI_TIER_LIST[*]}"')"
check "mem_cli_init: an empty MEM_CLI_TIERS counts as unset (detected tiers)" "0|tiers: neon" "$out"
out="$(run_script unknown BWA_MEM3="$neon" MEM_CLI_TIERS=sse41 -- 'mem_cli_init unknown')"
check "mem_cli_init: an unknown MEM_CLI_TIERS tier fails" \
    "1|FAIL: unknown: MEM_CLI_TIERS names 'sse41', not one of: avx2 avx512bw neon" "$out"
out="$(run_script repeated BWA_MEM3="$neon" "MEM_CLI_TIERS=neon neon" -- 'mem_cli_init repeated')"
check "mem_cli_init: a repeated MEM_CLI_TIERS tier fails" \
    "1|FAIL: repeated: MEM_CLI_TIERS names 'neon' more than once" "$out"
out="$(run_script floor BWA_MEM3="$sse41" -- 'mem_cli_init floor; echo REACHED')"
check "mem_cli_init: a host below the floor is a SKIP, exit 0" 0 "${out%%|*}"
check "mem_cli_init: the below-floor SKIP says so" yes \
    "$(grep -q '^SKIP: floor: host tier .sse41. is below the avx2 floor' "$work/floor.out" && echo yes || echo no)"
# Expanded by the generated script, not here.
# shellcheck disable=SC2016
verdict='if mem_cli_verdict; then echo "PASS: $MEM_CLI_VERDICT"; else echo "FAIL: $MEM_CLI_VERDICT"; exit 1; fi'
out="$(run_script none BWA_MEM3="$neon" -- "mem_cli_init none; $verdict")"
check "mem_cli_verdict: zero checks is a FAIL" "1|FAIL: none: no check ran (tiers: neon )" "$out"
out="$(run_script bad BWA_MEM3="$neon" -- "mem_cli_init bad; mem_cli_expect x a b; $verdict")"
check "mem_cli_verdict: a failed check is a FAIL" "1|FAIL: bad: 1 of 1 checks failed (tiers: neon )" "$out"
out="$(run_script good BWA_MEM3="$neon" -- "mem_cli_init good; mem_cli_expect x a a; mem_cli_skip y z; $verdict")"
check "mem_cli_verdict: passing checks are a PASS with skips counted" \
    "0|PASS: good: 1 checks on tiers: neon (1 skipped)" "$out"

# --- the checks, in this shell ---------------------------------------------
# shellcheck source=test/framework/mem_cli.sh
source "$framework/mem_cli.sh"
MEM_CLI_NAME=selftest
MEM_CLI_WORK="$work"
MEM_CLI_CHECKS=0
MEM_CLI_FAILS=0

mem_cli_expect "a mismatch" a b 2> /dev/null
check "mem_cli_expect counts a mismatch" "1/1" "$MEM_CLI_FAILS/$MEM_CLI_CHECKS"

tab=$'\t'
printf '@HD\tVN:1.6\n@PG\tID:a\nr1\t0\tctg1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tAS:i:4\n' > "$work/a.sam"
printf '@HD\tVN:1.6\n@PG\tID:b\nr1\t0\tctg1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tAS:i:4\n' > "$work/b.sam"
printf '@HD\tVN:1.6\nr1\t0\tctg1\t6\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tAS:i:4\n' > "$work/c.sam"
printf '@HD\tVN:1.6\n@PG\tID:a\n' > "$work/h.sam"
for case in "a b 0 equal but for @PG" "a c 1 differing" "h h 1 header-only" "nope1 nope2 1 missing"; do
    read -r x y want label <<< "$case"
    before=$MEM_CLI_FAILS
    mem_cli_same_file "$label" "$work/$x.sam" "$work/$y.sam" > /dev/null 2>&1
    check "mem_cli_same_file: $label" "$want" "$((MEM_CLI_FAILS - before))"
done

printf 'r1\t0\tctg1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\nr1\t16\tctg2\t9\t60\t4M\t*\t0\t0\tACGT\tIIII\nr1\t256\tctg3\t1\t0\t4M\t*\t0\t0\tACGT\tIIII\n' > "$work/two.sam"
check "mem_cli_primary: two primaries yield nothing" "" "$(mem_cli_primary "$work/two.sam" r1 2> /dev/null)"
check "mem_cli_primary: the one primary of a read" "r1${tab}0${tab}ctg1" \
    "$(mem_cli_primary "$work/a.sam" r1 | cut -f 1-3)"
rec="$(mem_cli_primary "$work/a.sam" r1)"
check "mem_cli_tag: a present tag" 4 "$(mem_cli_tag "$rec" AS)"
check "mem_cli_tag: an absent tag" "" "$(mem_cli_tag "$rec" XS)"

printf '>c\nACGTACGT\n' > "$work/ref.fa"
check "mem_cli_seq: an in-range slice" GTAC "$(mem_cli_seq "$work/ref.fa" c 3 4)"
rc=0
mem_cli_seq "$work/ref.fa" c 6 4 > /dev/null 2> "$work/seq.err" || rc=$?
check "mem_cli_seq: an out-of-range slice returns 1" 1 "$rc"
rc=0
mem_cli_seq "$work/ref.fa" nosuch 1 0 > /dev/null 2>&1 || rc=$?
check "mem_cli_seq: an unknown contig returns 1, even for an empty slice" 1 "$rc"
check "mem_cli_seq: and says so" "FAIL: selftest: mem_cli_seq c:6+4 is outside the reference $work/ref.fa" \
    "$(cat "$work/seq.err")"

if ((fails > 0)); then
    echo "FAIL: mem_cli_selftest: $fails check(s) failed"
    exit 1
fi
echo "PASS: mem_cli_selftest"
