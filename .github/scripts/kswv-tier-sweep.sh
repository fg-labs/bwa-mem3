#!/usr/bin/env bash
# .github/scripts/kswv-tier-sweep.sh
#
# Run the unit/kswv doctest suite once per SIMD tier this host can actually
# host, and record which tiers were exercised.
#
# Why this exists
# ---------------
# The kswv job used to run the suite exactly once, at whatever tier the runner
# happened to expose: `grep -q avx512bw /proc/cpuinfo` picked MAKE_ARCH and
# silently fell back to avx2. AVX-512BW is one of the two tiers that carried
# the PR #290 query-padding bug, so a runner-fleet refresh could drop precisely
# the coverage that matters most, and no artifact would say so. Worse, a single
# run only ever exercises ONE tier, while the kernels differ per tier by
# construction.
#
# How tier selection works
# ------------------------
# BWAMEM3_FORCE_TIER can only downgrade. An up-tier (or cross-family, or
# unknown) request is ignored with a warning and leaves the tier unchanged --
# see the runtime-line logic in src/simd_dispatch.cpp. Forcing blindly would
# therefore run the host tier N times and cheerfully report N tiers.
#
# So the gate is bwa-mem3's own banner: it annotates an ignored request as
# ", ignored" on the `SIMD runtime:` line. That makes the dispatcher the single
# authority on which tier a run actually used -- no second copy of the CPU
# detection logic to drift out of lockstep with it. (Same idiom as
# test/regression/all_tiers_parity.sh, which sweeps tiers for SAM parity.)
#
# Inputs (environment):
#   BWA_MEM3      path to the built bwa-mem3 binary, used as the tier oracle
#   TEST_BIN      path to bwa_mem3_tests_unit
#   OUT_DIR       directory to write per-tier JUnit XML into
#   ARCH_LABEL    "x86" | "arm" -- selects the candidate tier list, names files
# Optional:
#   REQUIRE_TIERS space-separated tiers that MUST have run; exit 1 otherwise.
#                 Unset by default: a runner without AVX-512 cannot be blamed
#                 for lacking it, and failing every such run would train people
#                 to ignore this job. Set it on a runner class you control to
#                 turn "we lost the avx512bw leg" into a hard error.
#   SUMMARY_FILE  file to append a markdown summary to (default $GITHUB_STEP_SUMMARY)
#   SUITE_ARGS    doctest --test-suite filter selecting which unit suites to sweep
#                 per tier (default: "--test-suite=unit/kswv"). doctest ORs a single
#                 --test-suite's comma-separated values (repeating the flag does NOT
#                 OR -- the last one wins), so sweep several kernel suites at once with
#                 "--test-suite=unit/kswv,unit/bandedswa". This is how a vector kernel
#                 other than kswv (e.g. the banded-SW getScores16 wrappers) gets
#                 exercised on the avx512bw tier, which no unit-test matrix row builds.
#                 Each comma-separated suite is run and validated SEPARATELY: every
#                 named suite must match >= 1 test case per tier or the tier is
#                 reported EMPTY (exit 1). A union run cannot enforce this -- doctest
#                 does not tag its JUnit testcases with the owning suite, so one
#                 misspelled suite would ride on a sibling's count and be dropped
#                 silently. Any non---test-suite tokens in SUITE_ARGS are preserved
#                 and passed through to every per-suite run.

set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${TEST_BIN:?TEST_BIN must be set}"
: "${OUT_DIR:?OUT_DIR must be set}"
: "${ARCH_LABEL:?ARCH_LABEL must be set}"

SUMMARY_FILE="${SUMMARY_FILE:-${GITHUB_STEP_SUMMARY:-/dev/null}}"
REQUIRE_TIERS="${REQUIRE_TIERS:-}"
SUITE_ARGS="${SUITE_ARGS:---test-suite=unit/kswv}"

mkdir -p "$OUT_DIR"

# Candidates, highest first so the most capable (and most bug-prone) kernel is
# reported at the top of the summary. arm64 has exactly one tier; the x86 list
# covers the three that have distinct kswv kernels plus sse41, whose kernels
# are exit() stubs and which therefore checks the skip path rather than the
# kernels (the unit tests gate on bwamem3_simd_tier() and skip cleanly there).
case "$ARCH_LABEL" in
    arm) CANDIDATES="neon" ;;
    x86) CANDIDATES="avx512bw avx2 sse41" ;;
    *)
        echo "::error::unknown ARCH_LABEL '$ARCH_LABEL' (expected x86 or arm)"
        exit 1
        ;;
esac

# Ask the dispatcher what a given force actually selects. Prints the effective
# tier name and returns 0; prints nothing and returns 0 when the request was
# refused (host cannot reach that tier); returns 2 when the oracle itself is
# unusable.
#
# The 2 matters. This runs inside a command substitution, i.e. a subshell, so
# an `exit` here would only leave the subshell and the loop would carry on --
# reporting every tier as "not available" and finishing with "no SIMD tier
# could be exercised", which reads as "old runner" and buries the real fault.
# Returning a distinct status lets the caller tell "can't" from "broken".
effective_tier() {
    local want="$1" line version_output
    # Take the binary's own exit status before looking at what it printed.
    # Piping straight into grep discards that status, so a `version` that
    # printed the line and then died would be read as a real tier answer
    # instead of the broken oracle it is.
    version_output="$(BWAMEM3_FORCE_TIER="$want" "$BWA_MEM3" version 2> /dev/null)" \
        || return 2
    line="$(printf '%s\n' "$version_output" | grep '^SIMD runtime:' || true)"
    [ -n "$line" ] || return 2
    case "$line" in
        *", ignored"*) return 0 ;; # refused: host cannot reach this tier
    esac
    printf '%s\n' "$line" | sed -n 's/^SIMD runtime: \([a-z0-9]*\).*/\1/p'
}

# doctest with --reporters=junit --out=FILE writes the report to FILE and
# prints nothing to stdout, so the run itself leaves no readable trace in the
# CI log. Recover the counts from the report so both the log and the summary
# say what actually happened per tier.
#
# Emits four fields: assertions failures errors testcases. The <testsuite>
# `tests` attribute is doctest's *assertion* count (p.numAsserts), NOT the
# number of test cases -- so it reads 0 for an empty selection but is also a
# poor liveness signal in general. The count of executed test cases is the
# number of <testcase> child elements, which doctest emits one-per-case at
# test_case_start (before the body runs), so a case that detects its tier and
# skips still shows up. That is the value the empty-selection gate keys on.
counts_from_junit() {
    local f="$1" line cases
    if [ ! -s "$f" ]; then
        echo "? ? ? ?"
        return
    fi
    # grep returns nonzero when it finds nothing (no <testsuite> line, zero
    # <testcase> children); under errexit that would abort the sweep, so absorb
    # it -- an empty match is a real state the callers below handle, not a fault.
    line="$(grep -m1 -oE '<testsuite [^>]*' "$f" || true)"
    cases="$(grep -c '<testcase ' "$f" || true)"
    printf '%s %s %s %s\n' \
        "$(printf '%s' "$line" | sed -n 's/.*tests="\([0-9]*\)".*/\1/p')" \
        "$(printf '%s' "$line" | sed -n 's/.*failures="\([0-9]*\)".*/\1/p')" \
        "$(printf '%s' "$line" | sed -n 's/.*errors="\([0-9]*\)".*/\1/p')" \
        "$cases"
}

# Suite names requested by SUITE_ARGS: split the LAST --test-suite=<csv> token on
# commas, one per line. doctest itself is last-flag-wins for a repeated
# --test-suite (only the final flag selects), so honor that here rather than
# unioning every token -- otherwise this per-suite runner would exercise suites a
# direct doctest invocation with the same args would ignore. Empty when SUITE_ARGS
# carries no --test-suite= token (a caller passing some other doctest filter --
# handled by the :ALL: fallback).
suites_from_args() {
    local tok suites=""
    # shellcheck disable=SC2086  # SUITE_ARGS is a deliberately word-split arg list
    for tok in $SUITE_ARGS; do
        case "$tok" in
            --test-suite=*) suites="${tok#--test-suite=}" ;; # last flag wins
        esac
    done
    [ -z "$suites" ] || printf '%s\n' "$suites" | tr ',' '\n'
}

# SUITE_ARGS with every --test-suite= token removed, so each per-suite run can
# re-add exactly one --test-suite filter while preserving any other doctest args.
other_args_from_suite_args() {
    local tok out=""
    # shellcheck disable=SC2086  # SUITE_ARGS is a deliberately word-split arg list
    for tok in $SUITE_ARGS; do
        case "$tok" in
            --test-suite=*) ;; # dropped; re-added one suite at a time
            *) out="$out $tok" ;;
        esac
    done
    printf '%s' "${out# }"
}

OTHER_ARGS="$(other_args_from_suite_args)"

ran=""
skipped=""
rows=""
overall=0

for tier in $CANDIDATES; do
    # effective_tier returns 2 (oracle broken) as well as 0 (tier printed, or
    # empty when refused). Capture that status through `if` so errexit does not
    # abort before the rc==2 "broken oracle" branch below can report it.
    if got="$(effective_tier "$tier")"; then
        rc=0
    else
        rc=$?
    fi
    if [ "$rc" -eq 2 ]; then
        echo "::error::no 'SIMD runtime:' banner from $BWA_MEM3 -- tier oracle is broken, not a low-tier runner"
        exit 1
    fi
    if [ -z "$got" ]; then
        skipped="$skipped $tier"
        echo "--- tier $tier: not available on this host, skipping"
        continue
    fi
    if [ "$got" != "$tier" ]; then
        # Banner accepted the request but reported a different tier. Never
        # expected; bail rather than attribute one tier's results to another.
        echo "::error::asked for tier '$tier' but dispatcher reports '$got'"
        exit 1
    fi

    echo "--- tier $tier: running $SUITE_ARGS"

    # Validate each requested suite INDEPENDENTLY. A multi-suite filter such as
    # --test-suite=unit/kswv,unit/bandedswa selects the UNION, and doctest does
    # not tag its JUnit <testcase> elements with the suite that owns them -- so a
    # single union run whose aggregate testcase count is positive cannot tell
    # that one named suite (renamed or misspelled) matched nothing while a sibling
    # carried the count. That would silently drop the coverage the missing suite
    # was added to guarantee. Run each suite on its own and require every one to
    # match >= 1 test case before recording a pass.
    suites=()
    while IFS= read -r s; do [ -n "$s" ] && suites+=("$s"); done < <(suites_from_args)
    if [ "${#suites[@]}" -eq 0 ]; then
        # SUITE_ARGS carried no --test-suite= token; run it as-is and gate on the
        # aggregate non-empty count (there are no individual suites to split).
        suites=(":ALL:")
    fi

    tier_status=pass # pass | FAILED | EMPTY
    empty_suites=""
    agg_tests=0
    agg_fail=0
    agg_err=0
    agg_cases=0
    for s in "${suites[@]}"; do
        if [ "$s" = ":ALL:" ]; then
            slabel="all"
            sfilter=()
            report="$OUT_DIR/kswv-unit-${ARCH_LABEL}-${tier}.xml"
        else
            slabel="$s"
            sfilter=(--test-suite="$s")
            report="$OUT_DIR/kswv-unit-${ARCH_LABEL}-${tier}-$(printf '%s' "$s" | tr '/,' '__').xml"
        fi
        # A failing suite exits nonzero; capture that through `if` so errexit
        # does not abort before the FAILED aggregation below records it.
        # shellcheck disable=SC2086  # OTHER_ARGS is a deliberately word-split arg list
        if BWAMEM3_FORCE_TIER="$tier" "$TEST_BIN" $OTHER_ARGS "${sfilter[@]}" \
            --reporters=junit --out="$report"; then
            s_status=0
        else
            s_status=$?
        fi
        read -r s_tests s_fail s_err s_cases << EOF
$(counts_from_junit "$report")
EOF
        # counts_from_junit emits "?" placeholders when the report is missing;
        # normalize to 0 so the aggregate arithmetic below never chokes.
        [ "$s_tests" -ge 0 ] 2> /dev/null || s_tests=0
        [ "$s_fail" -ge 0 ] 2> /dev/null || s_fail=0
        [ "$s_err" -ge 0 ] 2> /dev/null || s_err=0
        [ "$s_cases" -ge 0 ] 2> /dev/null || s_cases=0
        agg_tests=$((agg_tests + s_tests))
        agg_fail=$((agg_fail + s_fail))
        agg_err=$((agg_err + s_err))
        agg_cases=$((agg_cases + s_cases))
        if [ "$s_status" -ne 0 ]; then
            tier_status=FAILED
            echo "::error::suite '$slabel' failed at tier $tier (exit $s_status, ${s_fail} failures, ${s_err} errors)"
        elif [ "$s_cases" -le 0 ]; then
            # doctest exits 0 for a filter that matches no test cases, emitting
            # tests="0" with no <testcase> children. Refuse to record that as a
            # pass -- it means a suite was renamed out from under us or typo'd.
            [ "$tier_status" = pass ] && tier_status=EMPTY
            empty_suites="$empty_suites $slabel"
            echo "::error::suite '$slabel' selected 0 test cases at tier $tier -- empty selection, not a pass"
        else
            echo "    tier $tier / suite $slabel: $s_cases test cases, $s_tests assertions, $s_fail failures, $s_err errors"
        fi
    done

    if [ "$tier_status" = FAILED ]; then
        overall=1
        echo "::error::unit suite ($SUITE_ARGS) failed at tier $tier"
        rows="$rows| \`$tier\` | $agg_tests | $agg_fail | $agg_err | **FAILED** |"$'\n'
    elif [ "$tier_status" = EMPTY ]; then
        overall=1
        echo "::error::tier $tier: empty suite selection(s) --$empty_suites -- not a pass"
        rows="$rows| \`$tier\` | $agg_tests | $agg_fail | $agg_err | **EMPTY** |"$'\n'
    else
        echo "    tier $tier: ${#suites[@]} suite(s), $agg_cases test cases, $agg_tests assertions, $agg_fail failures, $agg_err errors -- pass"
        rows="$rows| \`$tier\` | $agg_tests | $agg_fail | $agg_err | pass |"$'\n'
    fi
    ran="$ran $tier"
done

# The point of the whole exercise: make the tier set a visible, durable fact
# rather than something you have to reconstruct from a ccache key.
{
    echo "### vector-kernel unit suites ($SUITE_ARGS) — SIMD tier sweep (${ARCH_LABEL})"
    echo
    echo "| tier | tests | failures | errors | status |"
    echo "|---|---|---|---|---|"
    printf '%s' "$rows"
    for t in $skipped; do echo "| \`$t\` | — | — | — | not available on this runner |"; done
    echo
} >> "$SUMMARY_FILE"

if [ -z "$ran" ]; then
    echo "::error::no SIMD tier could be exercised on this runner"
    exit 1
fi

case " $ran " in
    *" avx512bw "*) ;;
    *)
        if [ "$ARCH_LABEL" = "x86" ]; then
            echo "::warning::AVX-512BW was NOT exercised on this runner; the avx512bw vector kernels are untested by this run"
            echo "> :warning: \`avx512bw\` not exercised — that kernel is untested by this run." >> "$SUMMARY_FILE"
        fi
        ;;
esac

for req in $REQUIRE_TIERS; do
    case " $ran " in
        *" $req "*) ;;
        *)
            echo "::error::required tier '$req' was not exercised"
            overall=1
            ;;
    esac
done

echo "tiers exercised:${ran:- none}"
echo "tiers skipped:${skipped:- none}"
exit "$overall"
