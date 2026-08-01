#!/usr/bin/env bash
# test/regression/regression_coverage_lint_selftest.sh
#
# Regression: regression_coverage_lint.sh still notices a script nothing runs.
#
# The lint prints PASS when every script is covered, and would print the same
# PASS if its resolution silently stopped resolving -- an empty script list, or
# a covered_text that matches everything, both satisfy the check vacuously.
# From CI those are indistinguishable, which is the failure this whole family
# of lints exists to catch. So the cases below hand it trees that have real
# gaps and confirm it says so.
#
# Two of them guard bugs this lint actually had:
#
#   - the `make test-binaries` substring trap. Resolving invoked targets by
#     substring pulls the `test:` recipe in on any `make test-binaries` line,
#     which marks as covered precisely the scripts only `make test` runs --
#     the ones the lint was written to find.
#
#   - the SIGPIPE case. `printf "$big" | grep -q` under `set -o pipefail`
#     returns 141 when grep short-circuits on an early match, so a covered
#     script reads as uncovered once the text outgrows the pipe buffer. It
#     bit for real at ~87 KB of workflow text.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."

lint="$PWD/test/regression/regression_coverage_lint.sh"
fixture_root="$(mktemp -d)"
trap 'rm -rf "$fixture_root"' EXIT

failures=0

# $1 = workflow body, $2 = Makefile body, $3.. = regression script basenames
# Set FIXTURE_EXEMPTIONS to the body of the fixture's coverage_exemptions.txt;
# left empty, no file is written (the lint must cope with that too).
FIXTURE_EXEMPTIONS=""
build_fixture() {
    local workflow="$1" makefile="$2"
    shift 2
    local dir="$fixture_root/case" name

    rm -rf "$dir"
    mkdir -p "$dir/.github/workflows" "$dir/test/regression"
    printf '%s\n' "$workflow" > "$dir/.github/workflows/ci.yml"
    printf '%s\n' "$makefile" > "$dir/Makefile"
    for name in "$@"; do
        printf '#!/bin/bash\necho %s\n' "$name" > "$dir/test/regression/$name.sh"
    done
    if [[ -n $FIXTURE_EXEMPTIONS ]]; then
        printf '%s\n' "$FIXTURE_EXEMPTIONS" > "$dir/test/regression/coverage_exemptions.txt"
    fi
    printf '%s' "$dir"
}

run_case() {
    local description="$1" expected="$2" dir="$3"
    local output status verdict
    output="$(bash "$lint" "$dir" 2>&1)" && status=0 || status=$?
    case "$status" in
        0) verdict=PASS ;;
        1) verdict=FAIL ;;
        *) verdict="ERROR(exit $status)" ;;
    esac
    if [[ "$verdict" == "$expected" ]]; then
        echo "  ok   $description -> $verdict"
    else
        echo "  FAIL $description: expected $expected, got $verdict" >&2
        printf '%s\n' "$output" | sed 's/^/       /' >&2
        failures=$((failures + 1))
    fi
}

echo "== covered scripts =="
run_case "named directly in a workflow" PASS \
    "$(build_fixture 'jobs:
  build:
    steps:
      - run: bash test/regression/alpha.sh' \
        'test:
	echo nothing' alpha)"

run_case "in a Makefile target a workflow invokes" PASS \
    "$(build_fixture 'jobs:
  build:
    steps:
      - run: make check' \
        'check:
	./test/regression/alpha.sh' alpha)"

# One rule line may name several targets, and the recipe belongs to each. Read
# with a leading-name anchor, only `fmt` resolves and `lint` names nothing, so a
# workflow that invokes `make lint` reports alpha as run by no one.
run_case "every target on a multi-target rule line resolves" PASS \
    "$(build_fixture 'jobs:
  build:
    steps:
      - run: make lint' \
        'fmt lint:
	./test/regression/alpha.sh' alpha)"

echo "== uncovered scripts =="

# A target name is a literal, not a pattern. This Makefile really does define
# `.PHONY`, `.SUFFIXES`, `.c.o` and `.cpp.o`, so a `.` left unescaped in the
# invocation ERE lets a target match a `make` line that names something else --
# here `.check` matching `make xcheck`, which drags in a recipe no workflow runs
# and marks its scripts covered.
run_case "a dotted target does not wildcard onto another make invocation" FAIL \
    "$(build_fixture 'jobs:
  build:
    steps:
      - run: make xcheck' \
        'xcheck:
	echo building
.check:
	./test/regression/alpha.sh' alpha)"

# The same literalness on the recipe side: `.check` must attach to its own
# recipe, not to every rule line the `.` happens to wildcard over. beta is run
# by `xcheck`, which no workflow invokes, so it has to read as uncovered.
run_case "a dotted target attaches only to its own recipe" FAIL \
    "$(build_fixture 'jobs:
  build:
    steps:
      - run: make .check' \
        'xcheck:
	./test/regression/beta.sh
.check:
	./test/regression/alpha.sh' alpha beta)"
# `make` has to be the whole tool name, not a tail of one. Unanchored it matches
# inside `cmake`, `gmake` and `nmake`, so a line that never invokes this
# Makefile at all -- `cmake --build check` here -- resolves `check` and pulls in
# a recipe nothing runs, marking its scripts covered.
run_case "a make-suffixed tool name does not invoke a target" FAIL \
    "$(build_fixture 'jobs:
  build:
    steps:
      - run: cmake --build check' \
        'check:
	./test/regression/alpha.sh' alpha)"

# The anchor must not go the other way and lose a real invocation: make is
# routinely spelled with a path, and `/` is not part of the tool name.
run_case "make invoked by path still resolves its target" PASS \
    "$(build_fixture 'jobs:
  build:
    steps:
      - run: /usr/bin/make check' \
        'check:
	./test/regression/alpha.sh' alpha)"

run_case "referenced nowhere at all" FAIL \
    "$(build_fixture 'jobs:
  build:
    steps:
      - run: echo hi' \
        'test:
	echo nothing' orphan)"

run_case "only in a target no workflow invokes" FAIL \
    "$(build_fixture 'jobs:
  build:
    steps:
      - run: make test-binaries' \
        'test-binaries:
	echo building
test:
	./test/regression/alpha.sh' alpha)"

# The substring trap: `make test-binaries` must not count as invoking `test`.
run_case "make test-binaries does not pull in the test recipe" FAIL \
    "$(build_fixture 'jobs:
  build:
    steps:
      - run: make -j4 test-binaries' \
        'test-binaries:
	echo building
test:
	./test/regression/alpha.sh' alpha)"

# The same trap one level down, in the script names themselves. Every lint here
# ships next to an `_selftest` sibling whose name contains it, so a substring
# match reads `alpha` as covered on the strength of `alpha_selftest.sh` alone --
# and dropping alpha's own reference would then go unnoticed.
run_case "sibling _selftest reference does not cover the shorter name" FAIL \
    "$(build_fixture 'jobs:
  build:
    steps:
      - run: bash test/regression/alpha_selftest.sh' \
        'test:
	echo nothing' alpha alpha_selftest)"

run_case "both siblings referenced" PASS \
    "$(build_fixture 'jobs:
  build:
    steps:
      - run: bash test/regression/alpha_selftest.sh
      - run: bash test/regression/alpha.sh' \
        'test:
	echo nothing' alpha alpha_selftest)"

# Tightening the match must not go so far that ordinary prose stops counting: a
# reference is often the last thing in a sentence, and the Makefile already
# carries one written that way.
run_case "reference followed by punctuation still counts" PASS \
    "$(build_fixture 'jobs:
  build:
    steps:
      - run: echo see test/regression/alpha.sh.' \
        'test:
	echo nothing' alpha)"

echo "== exemptions =="
# Driven by the fixture's own coverage_exemptions.txt, so these cases keep
# testing the mechanism no matter what the repository's real list holds --
# including when it is empty.
FIXTURE_EXEMPTIONS='# reason goes here

orphan'
run_case "exempt script may be uncovered" PASS \
    "$(build_fixture 'jobs:
  build:
    steps:
      - run: echo hi' \
        'test:
	echo nothing' orphan)"

# An exemption that is no longer true has to be removed, or the list decays
# into a place where covered scripts go to be forgotten.
run_case "exempt script that IS covered is rejected" FAIL \
    "$(build_fixture 'jobs:
  build:
    steps:
      - run: bash test/regression/orphan.sh' \
        'test:
	echo nothing' orphan)"

# An exemption whose script no longer exists is the same decay in the other
# direction: it names nothing, so it can never be contradicted, and it outlives
# the rename or deletion that made it meaningless.
FIXTURE_EXEMPTIONS='# script was renamed away
ghost'
run_case "exemption naming no script is rejected" FAIL \
    "$(build_fixture 'jobs:
  build:
    steps:
      - run: bash test/regression/alpha.sh' \
        'test:
	echo nothing' alpha)"

# Comments and blank lines are not script names.
FIXTURE_EXEMPTIONS='# orphan
'
run_case "commented-out entry does not exempt" FAIL \
    "$(build_fixture 'jobs:
  build:
    steps:
      - run: echo hi' \
        'test:
	echo nothing' orphan)"
FIXTURE_EXEMPTIONS=""

echo "== resolution that stopped resolving must not report PASS =="
empty_wf="$(build_fixture 'jobs: {}' 'test:
	echo nothing' alpha)"
rm -f "$empty_wf/.github/workflows/ci.yml"
run_case "no workflow files" FAIL "$empty_wf"

no_scripts="$(build_fixture 'jobs:
  build:
    steps:
      - run: echo hi' 'test:
	echo nothing')"
run_case "no regression scripts" FAIL "$no_scripts"

no_makefile="$(build_fixture 'jobs:
  build:
    steps:
      - run: bash test/regression/alpha.sh' 'test:
	echo nothing' alpha)"
rm -f "$no_makefile/Makefile"
run_case "no Makefile" FAIL "$no_makefile"

# Regression for the SIGPIPE bug: the name appears at the very top, then far
# more text than a pipe buffer holds. A piped `grep -q` short-circuits, printf
# takes SIGPIPE, and pipefail turns the match into a miss.
big_dir="$(build_fixture 'jobs:
  build:
    steps:
      - run: bash test/regression/alpha.sh' \
    'test:
	echo nothing' alpha)"
{
    printf '# filler to push the workflow text past the pipe buffer\n'
    for _ in $(seq 1 4000); do
        printf '# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n'
    done
} >> "$big_dir/.github/workflows/ci.yml"
workflow_bytes="$(wc -c < "$big_dir/.github/workflows/ci.yml" | tr -d ' ')"
run_case "early match in >64KB of workflow text ($workflow_bytes bytes)" PASS "$big_dir"

bash "$lint" one two > /dev/null 2>&1 && usage_status=0 || usage_status=$?
if ((usage_status == 2)); then
    echo "  ok   too many arguments -> usage error (exit 2)"
else
    echo "  FAIL too many arguments: expected exit 2, got $usage_status" >&2
    failures=$((failures + 1))
fi

# Checked on the message, not just the status: an unguarded `cd` into a missing
# directory also dies with exit 1 under `set -e`, so a status-only assertion
# passes either way. The point of the case is that the failure arrives in this
# family's own format rather than as bash's raw `cd:` diagnostic.
absent_dir="$fixture_root/absent"
absent_out="$(bash "$lint" "$absent_dir" 2>&1)" && absent_status=0 || absent_status=$?
if ((absent_status == 1)) && [[ $absent_out == FAIL:* ]]; then
    echo "  ok   argument that is not a directory -> FAIL"
else
    echo "  FAIL argument that is not a directory: expected exit 1 and a FAIL: message," >&2
    echo "       got exit $absent_status: $absent_out" >&2
    failures=$((failures + 1))
fi

if ((failures > 0)); then
    echo "FAIL: regression_coverage_lint_selftest ($failures case(s))" >&2
    exit 1
fi

echo "PASS: regression_coverage_lint_selftest"
