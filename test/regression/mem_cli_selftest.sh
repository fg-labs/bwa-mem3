#!/usr/bin/env bash
# test/regression/mem_cli_selftest.sh
#
# Self-test for the shared shell helpers the mem CLI scripts are built on
# (test/framework/mem_cli.sh, test/framework/simd_tiers.sh and
# test/framework/dna.sh). A helper that
# stopped failing would make every script on it pass without checking
# anything, so each one is driven here with inputs whose answer is known:
#
#   simd_sweep_tiers   stub binaries reporting each dispatcher tier, and one
#                      reporting garbage: avx512bw -> avx2 avx512bw, avx2 ->
#                      avx2, sse41 -> none, scalar -> none, neon -> neon,
#                      garbage -> status 2; a BWAMEM3_FORCE_TIER in the
#                      caller's environment is ignored, and a sweep entry
#                      with no rank is skipped.
#   simd_detect_tier   a forced tier is reported as forced (avx2 on an
#                      avx512bw stub).
#   simd_tier_rank     scalar < sse41 < ... < avx512bw < neon; -1 for an
#                      unknown name.
#   simd_parse_tier    the tier of the first "SIMD tier:" line, or nothing.
#   dna_revcomp        reverses and complements, keeping case, N and '-'.
#   mem_cli_init       MEM_CLI_TIERS that is whitespace only, names an
#                      unknown tier or names a tier twice ends the script
#                      with FAIL; an empty one counts as unset; a host below
#                      the floor reports SKIP and exits 0.
#   mem_cli_verdict    zero checks is a failure, a failed check is a failure,
#                      and so is a mem_cli_seq slice that failed inside a
#                      command substitution, where its status is lost.
#   mem_cli_expect     a mismatch is counted as a failure.
#   mem_cli_same_file  differing, missing and header-only files all fail;
#                      equal files that differ only in @PG pass.
#   mem_cli_primary    a read with two primary records yields nothing.
#   mem_cli_tag        a present tag's value, nothing for an absent one.
#   mem_cli_ref_name   names the script, the stem and the pid.
#   mem_cli_tags       the optional fields of a record, TAB-joined.
#   mem_cli_expect_fail a nonzero exit with the expected stderr passes; exit
#                      0, or a nonzero exit with other stderr, fails.
#   mem_cli_count      counts the records an awk condition selects, never a
#                      header line, and 0 when none match.
#   mem_cli_mem_status an error exit is recorded, not fatal; a signal death
#                      and a successful run in the wrong tier end the
#                      script with FAIL.
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
bin="$(stub scalar)"
got="$(
    # shellcheck source=test/framework/simd_tiers.sh
    source "$framework/simd_tiers.sh"
    simd_sweep_tiers "$bin" && echo "$SIMD_HOST_TIER:${SIMD_TIERS[*]-}"
)"
check "simd_sweep_tiers on a scalar host: detected, nothing to sweep" "scalar:" "$got"
bin="$(stub avx2)"
got="$(
    # shellcheck source=test/framework/simd_tiers.sh
    source "$framework/simd_tiers.sh"
    SIMD_SWEEP_TIERS=(avx2 bogus neon)
    simd_sweep_tiers "$bin" && echo "${SIMD_TIERS[*]-}"
)"
check "simd_sweep_tiers skips a sweep entry it cannot rank" avx2 "$got"

# --- simd_detect_tier and simd_tier_rank ------------------------------------
# A stub whose dispatcher reports the forced tier, or avx512bw unforced.
bin="$work/stub-forceable"
# The stub's own ${BWAMEM3_FORCE_TIER:-} is expanded when the stub runs.
# shellcheck disable=SC2016
printf '#!/bin/sh\necho "[M::main] SIMD tier: ${BWAMEM3_FORCE_TIER:-avx512bw} (build baseline: x)" >&2\nexit 1\n' > "$bin"
chmod +x "$bin"
got="$(
    # shellcheck source=test/framework/simd_tiers.sh
    source "$framework/simd_tiers.sh"
    BWAMEM3_FORCE_TIER=avx2 simd_detect_tier "$bin" && echo "$SIMD_HOST_TIER"
)"
check "simd_detect_tier honours the caller's BWAMEM3_FORCE_TIER" avx2 "$got"
got="$(
    # shellcheck source=test/framework/simd_tiers.sh
    source "$framework/simd_tiers.sh"
    BWAMEM3_FORCE_TIER=avx2 simd_sweep_tiers "$bin" && echo "$SIMD_HOST_TIER"
)"
check "simd_sweep_tiers reads the unforced tier from the same binary" avx512bw "$got"
got="$(
    # shellcheck source=test/framework/simd_tiers.sh
    source "$framework/simd_tiers.sh"
    for t in scalar sse41 sse42 avx avx2 avx512bw neon bogus; do printf '%s ' "$(simd_tier_rank "$t")"; done
)"
check "simd_tier_rank orders the tiers and gives -1 for an unknown name" "0 1 2 3 4 5 6 -1 " "$got"
printf 'x\n[M::main] SIMD tier: avx2 (build baseline: x)\n[M::main] SIMD tier: neon\n' > "$work/tiers.log"
got="$(
    # shellcheck source=test/framework/simd_tiers.sh
    source "$framework/simd_tiers.sh"
    printf '%s|%s' "$(simd_parse_tier "$work/tiers.log")" "$(simd_parse_tier /dev/null)"
)"
check "simd_parse_tier: the first SIMD tier line's tier, else nothing" "avx2|" "$got"

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
# The generated script expands MEM_CLI_WORK itself.
# shellcheck disable=SC2016
out="$(run_script lostseq BWA_MEM3="$neon" -- 'mem_cli_init lostseq; printf ">c\nACGT\n" > "$MEM_CLI_WORK/r.fa"; x="$(echo "$(mem_cli_seq "$MEM_CLI_WORK/r.fa" c 3 9)")"; mem_cli_expect x a a; '"$verdict")"
check "mem_cli_verdict: a slice that failed inside \$(...) is a FAIL" \
    "1|FAIL: lostseq: 1 mem_cli_seq slice(s) fell outside the reference (tiers: neon )" "$out"
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
check "mem_cli_tags: every optional field" "NM:i:0${tab}AS:i:4" "$(mem_cli_tags "$rec")"

check "mem_cli_ref_name: script, stem and pid" "selftest.u3.$$.fa" "$(mem_cli_ref_name u3)"

# Expanded by awk, not here.
# shellcheck disable=SC2016
check "mem_cli_count: records matching a condition" 2 "$(mem_cli_count "$work/two.sam" '$1 == "r1" && $5 == 60')"
# shellcheck disable=SC2016
check "mem_cli_count: header lines never count" 0 "$(mem_cli_count "$work/a.sam" '$1 ~ /^@/')"
# shellcheck disable=SC2016
check "mem_cli_count: no match is 0" 0 "$(mem_cli_count "$work/two.sam" '$3 == "ctg9"')"

# A stub `mem` that exits with <status>, after printing the banner of <tier>.
status_stub() {
    local path="$work/status-stub-$1-$2"
    printf '#!/bin/sh
echo "Executing in %s mode" >&2
%s
' "$2" "$3" > "$path"
    chmod +x "$path"
    printf '%s\n' "$path"
}
# run_status <stub-args...>: "<script status>|<MEM_CLI_STATUS>|<first FAIL: line>"
# (bash itself reports a killed job on stderr first)
# for mem_cli_mem_status at tier avx2 in a subshell.
run_status() {
    local bin rc=0
    bin="$(status_stub "$@")"
    (
        BWA_MEM3="$bin"
        mem_cli_mem_status avx2 "$work/status.out"
        echo "$MEM_CLI_STATUS" > "$work/status.code"
    ) 2> "$work/status.err" || rc=$?
    printf '%s|%s|%s\n' "$rc" "$(cat "$work/status.code" 2> /dev/null)" "$(grep -m 1 '^FAIL:' "$work/status.err" || true)"
    rm -f "$work/status.code"
}
check "mem_cli_mem_status: an error exit is recorded" "0|1|" "$(run_status exit1 avx2 'exit 1')"
check "mem_cli_mem_status: success in the forced tier" "0|0|" "$(run_status ok avx2 'exit 0')"
check "mem_cli_mem_status: a signal death ends the script" \
    "1||FAIL: selftest: bwa-mem3 mem  (tier avx2) was killed by signal 9:" \
    "$(run_status kill avx2 'kill -9 $$')"
# run_fail <stub-args...> <pattern>: the failures mem_cli_expect_fail counted.
run_fail() {
    local bin pattern="${*: -1}" before=$MEM_CLI_FAILS
    bin="$(status_stub "${@:1:3}")"
    BWA_MEM3="$bin" mem_cli_expect_fail "x" "$pattern" avx2 "$work/fail.out" > /dev/null 2>&1
    echo $((MEM_CLI_FAILS - before))
}
check "mem_cli_expect_fail: nonzero with the expected stderr passes" 0 \
    "$(run_fail fmsg avx2 'echo "[E::x] bad thing" >&2; exit 1' 'bad thing')"
check "mem_cli_expect_fail: nonzero with other stderr fails" 1 \
    "$(run_fail fother avx2 'echo "[E::x] other" >&2; exit 1' 'bad thing')"
check "mem_cli_expect_fail: exit 0 fails" 1 "$(run_fail fzero avx2 'exit 0' 'bad thing')"
check "mem_cli_mem_status: success in another tier ends the script" \
    "1||FAIL: selftest: BWAMEM3_FORCE_TIER=avx2 was not honoured:" "$(run_status wrong neon 'exit 0')"

check "dna_revcomp: reverses, complements, keeps case, N and '-'" "-nacgtNACGT" "$(dna_revcomp "ACGTNacgtn-")"

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
