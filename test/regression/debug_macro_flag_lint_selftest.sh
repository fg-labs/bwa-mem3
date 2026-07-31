#!/usr/bin/env bash
# test/regression/debug_macro_flag_lint_selftest.sh
#
# Regression: debug_macro_flag_lint.sh still detects a list that has drifted.
#
# The lint prints PASS when the two lists agree -- and also when its extraction
# has quietly stopped extracting, because an empty flag set trivially satisfies
# every check over it. Those two states are indistinguishable from CI, which is
# the same "green while doing no work" failure the lint exists to catch. The
# only way to separate them is to hand it lists that have drifted and confirm
# it says so.
#
# Fixtures are generated per case: a minimal .github/workflows/ci.yml carrying
# an EXTRA_CXXFLAGS line, plus a src/ holding preprocessor conditionals.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."

lint="$PWD/test/regression/debug_macro_flag_lint.sh"
fixture_root="$(mktemp -d)"
trap 'rm -rf "$fixture_root"' EXIT

failures=0
case_seq=0

# Build a fixture repo: $1 = the EXTRA_CXXFLAGS line body, $2 = src/ content,
# $3 = optional text appended after the assignment's closing quote.
#
# Each call gets its own directory rather than reusing one path. Several cases
# below build a fixture and then append to its ci.yml; sharing a path would make
# those mutations survive only as long as no other case is built in between, and
# a clobbered fixture still reports `ok` against whatever replaced it.
#
# The flags body goes through printf, not echo: some of it carries a literal
# backslash-newline to exercise a wrapped assignment, and echo may interpret it.
build_fixture() {
    local flags="$1" src_body="$2" line_suffix="${3-}" dir

    case_seq=$((case_seq + 1))
    dir="$fixture_root/case-$case_seq"
    mkdir -p "$dir/.github/workflows" "$dir/src"
    {
        printf 'jobs:\n  build:\n    steps:\n      - run: |\n'
        printf '          make EXTRA_CXXFLAGS="%s"%s\n' "$flags" "$line_suffix"
    } > "$dir/.github/workflows/ci.yml"
    printf '%s\n' "$src_body" > "$dir/src/fixture.cpp"
    printf '%s' "$dir"
}

# $1 = description, $2 = PASS|FAIL, $3 = flags line, $4 = src content,
# $5 = optional substring the output must contain.
#
# The exit status alone is not discriminating: the lint exits 1 from several
# distinct checks, and it also exits 1 when it aborts outright (an unbound
# variable under `set -u`, say). Without $5 a case that means to exercise one
# branch is satisfied by any other -- including a crash -- so each drift case
# below pins the message that names its own branch.
check_case() {
    local description="$1" expected="$2" expect_text="${5-}" dir
    dir="$(build_fixture "$3" "$4")"

    local output status verdict
    output="$(bash "$lint" "$dir" 2>&1)" && status=0 || status=$?
    case "$status" in
        0) verdict=PASS ;;
        1) verdict=FAIL ;;
        *) verdict="ERROR(exit $status)" ;;
    esac

    if [[ "$verdict" == "$expected" ]] \
        && [[ -z $expect_text || "$output" == *"$expect_text"* ]]; then
        echo "  ok   $description -> $verdict"
    else
        echo "  FAIL $description: expected $expected${expect_text:+ matching \"$expect_text\"}, got $verdict" >&2
        printf '%s\n' "$output" | sed 's/^/       /' >&2
        failures=$((failures + 1))
    fi
}

# $1 = description, $2 = fixture dir, $3 = expected exit status,
# $4 = optional substring the output must contain.
#
# Sibling of check_case for the cases that have to mutate a fixture after
# build_fixture returns, which check_case cannot express. Same reason for
# pinning the message: "nonzero" alone is satisfied by a syntax error or by an
# unrelated early guard firing first, neither of which is the branch under test.
check_dir() {
    local description="$1" dir="$2" expect_status="$3" expect_text="${4-}"

    local output status
    output="$(bash "$lint" "$dir" 2>&1)" && status=0 || status=$?

    if (( status == expect_status )) \
        && [[ -z $expect_text || "$output" == *"$expect_text"* ]]; then
        echo "  ok   $description -> exit $status"
    else
        echo "  FAIL $description: expected exit $expect_status${expect_text:+ matching \"$expect_text\"}, got exit $status" >&2
        printf '%s\n' "$output" | sed 's/^/       /' >&2
        failures=$((failures + 1))
    fi
}

echo "== lists that agree =="
check_case "flag names a real macro" PASS \
    '-DBWA_MEM3_DEBUG_ALPHA' \
    '#ifdef BWA_MEM3_DEBUG_ALPHA
int alpha;
#endif'
check_case "several flags, all matched" PASS \
    '-DBWA_MEM3_DEBUG_ALPHA -DBWA_MEM3_DEBUG_BETA' \
    '#ifdef BWA_MEM3_DEBUG_ALPHA
#endif
#if defined(BWA_MEM3_DEBUG_BETA)
#endif'
check_case "non-family macro needs no flag" PASS \
    '-DBWA_MEM3_DEBUG_ALPHA' \
    '#ifdef BWA_MEM3_DEBUG_ALPHA
#endif
#ifdef STAGE_PROF
#endif'

echo "== forward drift: a flag that enables nothing =="
check_case "flag names no macro" FAIL \
    '-DBWA_MEM3_DEBUG_ALPHA -DBWA_MEM3_DEBUG_TYPO' \
    '#ifdef BWA_MEM3_DEBUG_ALPHA
#endif' \
    'passes -DBWA_MEM3_DEBUG_TYPO'
# The whole point: a mention in prose is not a compiled gate. The unrelated
# conditional keeps the src/ macro set non-empty, so this reaches the
# forward-drift branch instead of stopping at the "nothing was checked" guard.
check_case "macro only named in a comment" FAIL \
    '-DBWA_MEM3_DEBUG_ALPHA' \
    '// BWA_MEM3_DEBUG_ALPHA would guard this, if anything did
int alpha;
#ifdef STAGE_PROF
#endif' \
    'passes -DBWA_MEM3_DEBUG_ALPHA'
# The same rule on the directive line itself, which the case above does not
# reach. A trailing comment there is not a gate either, so it must not be
# demanded of the flag list: reporting that src/ gates code on it would be false,
# and adding the -D to silence that would enable a macro nothing tests.
check_case "macro named in a comment on the #ifdef line" PASS \
    '-DBWA_MEM3_DEBUG_ALPHA' \
    '#ifdef BWA_MEM3_DEBUG_ALPHA  // see also BWA_MEM3_DEBUG_BETA, not yet used
#endif
#if defined(BWA_MEM3_DEBUG_ALPHA) /* not BWA_MEM3_DEBUG_GAMMA */
#endif
#ifdef BWA_MEM3_DEBUG_ALPHA /* a * inside must not defeat the strip: BWA_MEM3_DEBUG_DELTA */
#endif' \
    'PASS: debug_macro_flag_lint'
# The other half of the same rule, and the half that matters more: stripping a
# comment must never take real code with it. A closed comment mid-condition is
# left in place precisely so the macro after it survives -- getting this wrong
# drops a gate silently, which is the failure mode this lint exists to catch.
check_case "macro after a mid-condition comment still counts" PASS \
    '-DBWA_MEM3_DEBUG_ALPHA -DBWA_MEM3_DEBUG_BETA' \
    '#if defined(BWA_MEM3_DEBUG_ALPHA) /* 2*3 */ || defined(BWA_MEM3_DEBUG_BETA)
#endif' \
    'PASS: debug_macro_flag_lint'

# The exempt macro is read out of the lint's own pending_flags array rather than
# hardcoded: that array is documented as something to empty once the pending
# check lands, and a hardcoded name would turn that deletion into a spurious
# failure here. No entries left means there is nothing to exercise.
#
# No `| head -1` on the end: under `set -o pipefail` head exits after the first
# line, the upstream stage dies of SIGPIPE, and the pipeline reports 141 -- so
# the assignment fails and `set -e` kills this script. That is the same trap the
# lint documents above its membership tests. Take the first line in the shell
# instead, and `|| true` for the no-entries case, where grep exits 1.
# Strip the indent with sed rather than `tr -d '[:space:]'`: tr would delete the
# newlines too, collapsing a multi-entry array into one concatenated token that
# the first-line split below cannot separate.
pending_flags_in_lint="$(sed -n '/^pending_flags=(/,/^)/p' "$lint" \
    | sed -n 's/^[[:space:]]*\(BWA_MEM3_DEBUG_[A-Za-z0-9_]*\).*/\1/p')"
pending_flag="${pending_flags_in_lint%%$'\n'*}"
if [[ -n $pending_flag ]]; then
    check_case "pending flag is exempt" PASS \
        "-DBWA_MEM3_DEBUG_ALPHA -D$pending_flag" \
        '#ifdef BWA_MEM3_DEBUG_ALPHA
#endif' \
        "note: $pending_flag names no macro"
else
    echo "  skip pending flag is exempt (lint has no pending_flags entries)"
fi

# The deletion the array documents, exercised now rather than discovered later:
# an empty pending_flags must still reach the FAIL branch. Bash before 4.4
# aborts on the plain expansion of an empty array under `set -u`, and that abort
# also exits 1 -- so this case only means something because it pins the message
# rather than the status.
#
# Run through the shebang instead of `bash "$lint"`. The lint declares
# `#!/bin/bash`, which is 3.2 on macOS, while `bash` on PATH is whatever the
# developer installed -- and under a 4.4+ bash this case cannot observe the
# abort it exists to catch.
empty_pending_lint="$fixture_root/lint_no_pending.sh"
sed '/^pending_flags=(/,/^)/{ /^pending_flags=(/!{ /^)/!d; }; }' \
    "$lint" > "$empty_pending_lint"
chmod +x "$empty_pending_lint"
remaining="$(sed -n '/^pending_flags=(/,/^)/p' "$empty_pending_lint" \
    | grep -cE '^[[:space:]]*BWA_MEM3_DEBUG_' || true)"
empty_dir="$(build_fixture '-DBWA_MEM3_DEBUG_ALPHA -DBWA_MEM3_DEBUG_TYPO' \
    '#ifdef BWA_MEM3_DEBUG_ALPHA
#endif')"
empty_out="$("$empty_pending_lint" "$empty_dir" 2>&1)" \
    && empty_status=0 || empty_status=$?
if (( remaining != 0 )); then
    echo "  FAIL could not empty pending_flags; the empty-array case checked nothing" >&2
    failures=$((failures + 1))
elif (( empty_status == 1 )) \
    && [[ "$empty_out" == *'passes -DBWA_MEM3_DEBUG_TYPO'* ]]; then
    echo "  ok   empty pending_flags -> FAIL (no unbound-variable abort)"
else
    echo "  FAIL empty pending_flags: expected the forward-drift message, got exit $empty_status" >&2
    printf '%s\n' "$empty_out" | sed 's/^/       /' >&2
    failures=$((failures + 1))
fi

echo "== reverse drift: a macro no CI row compiles =="
check_case "src macro missing from flag list" FAIL \
    '-DBWA_MEM3_DEBUG_ALPHA' \
    '#ifdef BWA_MEM3_DEBUG_ALPHA
#endif
#ifdef BWA_MEM3_DEBUG_ORPHAN
#endif' \
    'gates code on BWA_MEM3_DEBUG_ORPHAN'

echo "== a backslash continuation is one logical line =="
# Both extractions are line-oriented, and both of the constructs they read can
# be split across physical lines. Reading only the first one drops whatever the
# continuation carries, and drops it silently -- the same "green while doing no
# work" failure the lint exists to catch, reopened inside the lint itself.
#
# src/ already splits preprocessor conditions this way today, so the reverse
# direction is the live risk: a continued condition's second macro would never
# be required to have a -D.
check_case "macro on a continued #if line" FAIL \
    '-DBWA_MEM3_DEBUG_ALPHA' \
    '#if defined(BWA_MEM3_DEBUG_ALPHA) || \
    defined(BWA_MEM3_DEBUG_CONTINUED)
#endif' \
    'gates code on BWA_MEM3_DEBUG_CONTINUED'
# Forward direction, same shape. The flag list is one long line today; wrapping
# it is the natural next edit, and a continuation keeps the single-line anchor's
# count at 1 -- so the trailing flags would go unchecked while every guard in the
# script still reported success.
check_case "flag on a continued EXTRA_CXXFLAGS line" FAIL \
    '-DBWA_MEM3_DEBUG_ALPHA \
            -DBWA_MEM3_DEBUG_WRAPPED' \
    '#ifdef BWA_MEM3_DEBUG_ALPHA
#endif' \
    'passes -DBWA_MEM3_DEBUG_WRAPPED'

echo "== extraction that stopped extracting must not report PASS =="
check_case "no flag list in the workflow" FAIL \
    '-DSTAGE_PROF' \
    '#ifdef BWA_MEM3_DEBUG_ALPHA
#endif' \
    'expected exactly one -DBWA_MEM3_DEBUG flag list'
check_case "src with no preprocessor conditionals" FAIL \
    '-DBWA_MEM3_DEBUG_ALPHA' \
    'int main(void) { return 0; }' \
    'no preprocessor conditionals found under'
echo "== the anchor is scoped to the flag line, not the whole workflow =="
# Two independent mechanisms keep prose out of the flag set -- requiring
# EXTRA_CXXFLAGS on the line, and dropping comment lines outright -- so there is
# a case per mechanism. One fixture carrying both would let either mechanism
# cover for the other's absence.
#
# Mechanism 1: a comment that names only a macro is filtered by the missing
# EXTRA_CXXFLAGS. The macro is one that appears nowhere else, not the one on the
# flag line: repeating the flag-line macro would leave this blind to an extractor
# that scans the whole workflow, since it would `sort -u` both mentions back into
# the same one-flag set and still print PASS. A distinct name lands in the flag
# set only if the comment was read, and then names no macro in src/, so the
# forward-drift branch fails the case.
prefix_comment_dir="$(build_fixture '-DBWA_MEM3_DEBUG_ALPHA' '#ifdef BWA_MEM3_DEBUG_ALPHA
#endif')"
echo '      # note: -DBWA_MEM3_DEBUG_COMMENT_ONLY is deliberately not enabled' \
    >> "$prefix_comment_dir/.github/workflows/ci.yml"
check_dir "comment naming only a macro" "$prefix_comment_dir" 0 \
    'PASS: debug_macro_flag_lint'

# Mechanism 2: a comment that quotes the whole assignment carries EXTRA_CXXFLAGS
# too, so requiring both terms does not filter it -- it would count as a second
# flag list and fail claiming the step was restructured. Only the comment filter
# catches this one.
full_comment_dir="$(build_fixture '-DBWA_MEM3_DEBUG_ALPHA' '#ifdef BWA_MEM3_DEBUG_ALPHA
#endif')"
echo '      # was: EXTRA_CXXFLAGS="-DBWA_MEM3_DEBUG_ALPHA" before the rename' \
    >> "$full_comment_dir/.github/workflows/ci.yml"
check_dir "comment quoting the whole assignment" "$full_comment_dir" 0 \
    'PASS: debug_macro_flag_lint'

# A trailing comment on the real flag line is part of that line, not a comment
# line, so the filter must leave it alone. Filtering on "contains #" rather than
# "starts with #" would take the flag list with it and report an absent list.
trailing_comment_dir="$(build_fixture '-DBWA_MEM3_DEBUG_ALPHA' '#ifdef BWA_MEM3_DEBUG_ALPHA
#endif' '  # keep this list in sync with src/')"
check_dir "trailing comment on the flag line survives" "$trailing_comment_dir" 0 \
    'PASS: debug_macro_flag_lint'

# A `#` comment ends at its own newline in YAML and in shell alike, so a trailing
# backslash inside one continues nothing. Dropping comments only after joining
# would let this comment absorb the flag list below it and discard both, leaving
# the list looking absent.
comment_backslash_dir="$fixture_root/comment-backslash"
mkdir -p "$comment_backslash_dir/.github/workflows" "$comment_backslash_dir/src"
{
    printf 'jobs:\n  build:\n    steps:\n      - run: |\n'
    printf '      # a note that happens to end in a backslash \\\n'
    printf '          make EXTRA_CXXFLAGS="-DBWA_MEM3_DEBUG_ALPHA"\n'
} > "$comment_backslash_dir/.github/workflows/ci.yml"
printf '#ifdef BWA_MEM3_DEBUG_ALPHA\n#endif\n' \
    > "$comment_backslash_dir/src/fixture.cpp"
check_dir "comment ending in a backslash keeps the next line" \
    "$comment_backslash_dir" 0 'PASS: debug_macro_flag_lint'

# Two flag lines means the step was restructured and the single-line anchor no
# longer describes it; guessing which one is authoritative would be the bug.
two_line_dir="$(build_fixture '-DBWA_MEM3_DEBUG_ALPHA' '#ifdef BWA_MEM3_DEBUG_ALPHA
#endif')"
echo '          make EXTRA_CXXFLAGS="-DBWA_MEM3_DEBUG_ALPHA"' \
    >> "$two_line_dir/.github/workflows/ci.yml"
check_dir "two flag lines" "$two_line_dir" 1 \
    'expected exactly one -DBWA_MEM3_DEBUG flag list'

missing_dir="$fixture_root/missing"
mkdir -p "$missing_dir/src"
check_dir "missing workflow" "$missing_dir" 1 \
    '.github/workflows/ci.yml missing under'

# A path of the right name but the wrong type exists, so the existence check
# above passes it through. Both of these already exited nonzero before the type
# checks were added -- but blaming a restructured step and an empty preprocessor
# set respectively, which is the wrong bug. These cases pin the diagnostic, not
# just the status.
wrong_workflow_dir="$fixture_root/wrong-type-workflow"
mkdir -p "$wrong_workflow_dir/.github/workflows/ci.yml" "$wrong_workflow_dir/src"
check_dir "workflow path is a directory" "$wrong_workflow_dir" 1 \
    '.github/workflows/ci.yml is not a regular file'

wrong_src_dir="$fixture_root/wrong-type-src"
mkdir -p "$wrong_src_dir/.github/workflows"
printf 'jobs:\n          make EXTRA_CXXFLAGS="-DBWA_MEM3_DEBUG_ALPHA"\n' \
    > "$wrong_src_dir/.github/workflows/ci.yml"
printf '#ifdef BWA_MEM3_DEBUG_ALPHA\n#endif\n' > "$wrong_src_dir/src"
check_dir "src path is a regular file" "$wrong_src_dir" 1 \
    'src is not a directory'

# The expected text is derived from $lint rather than spelled out: the lint
# prints its own basename, so a hardcoded name would make this case fail for a
# renamed copy of the script instead of for the behaviour under test.
usage_out="$(bash "$lint" one two 2>&1)" && usage_status=0 || usage_status=$?
if (( usage_status == 2 )) && [[ "$usage_out" == *"usage: ${lint##*/}"* ]]; then
    echo "  ok   too many arguments -> usage error (exit 2)"
else
    echo "  FAIL too many arguments: expected the usage message and exit 2, got exit $usage_status" >&2
    printf '%s\n' "$usage_out" | sed 's/^/       /' >&2
    failures=$((failures + 1))
fi

if (( failures > 0 )); then
    echo "FAIL: debug_macro_flag_lint_selftest ($failures case(s))" >&2
    exit 1
fi

echo "PASS: debug_macro_flag_lint_selftest"
