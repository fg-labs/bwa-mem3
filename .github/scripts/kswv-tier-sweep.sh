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

set -uo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${TEST_BIN:?TEST_BIN must be set}"
: "${OUT_DIR:?OUT_DIR must be set}"
: "${ARCH_LABEL:?ARCH_LABEL must be set}"

SUMMARY_FILE="${SUMMARY_FILE:-${GITHUB_STEP_SUMMARY:-/dev/null}}"
REQUIRE_TIERS="${REQUIRE_TIERS:-}"

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
counts_from_junit() {
    local f="$1" line
    if [ ! -s "$f" ]; then
        echo "? ? ?"
        return
    fi
    line="$(grep -m1 -oE '<testsuite [^>]*' "$f")"
    printf '%s %s %s\n' \
        "$(printf '%s' "$line" | sed -n 's/.*tests="\([0-9]*\)".*/\1/p')" \
        "$(printf '%s' "$line" | sed -n 's/.*failures="\([0-9]*\)".*/\1/p')" \
        "$(printf '%s' "$line" | sed -n 's/.*errors="\([0-9]*\)".*/\1/p')"
}

ran=""
skipped=""
rows=""
overall=0

for tier in $CANDIDATES; do
    got="$(effective_tier "$tier")"
    rc=$?
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

    echo "--- tier $tier: running unit/kswv"
    report="$OUT_DIR/kswv-unit-${ARCH_LABEL}-${tier}.xml"
    BWAMEM3_FORCE_TIER="$tier" "$TEST_BIN" --test-suite="unit/kswv" \
        --reporters=junit --out="$report"
    status=$?
    read -r n_tests n_fail n_err << EOF
$(counts_from_junit "$report")
EOF
    if [ "$status" -ne 0 ]; then
        overall=$status
        echo "::error::unit/kswv failed at tier $tier (exit $status, ${n_fail} failures, ${n_err} errors)"
        rows="$rows| \`$tier\` | $n_tests | $n_fail | $n_err | **FAILED** |"$'\n'
    else
        echo "    tier $tier: $n_tests tests, $n_fail failures, $n_err errors -- pass"
        rows="$rows| \`$tier\` | $n_tests | $n_fail | $n_err | pass |"$'\n'
    fi
    ran="$ran $tier"
done

# The point of the whole exercise: make the tier set a visible, durable fact
# rather than something you have to reconstruct from a ccache key.
{
    echo "### kswv unit suite — SIMD tier sweep (${ARCH_LABEL})"
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
            echo "::warning::AVX-512BW was NOT exercised on this runner; the avx512bw kswv kernel is untested by this run"
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
