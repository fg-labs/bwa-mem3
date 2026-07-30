#!/usr/bin/env bash
# test/regression/ndebug_gate_lint_selftest.sh
#
# Regression: ndebug_gate_lint.sh still recognises an NDEBUG gate.
#
# The lint it tests has exactly one observable behaviour on a healthy tree --
# it prints PASS -- and that is also what it prints when its matcher has
# stopped matching anything at all. Those two states are indistinguishable from
# CI, which is the same "green while doing no work" failure the lint itself
# exists to catch upstream. The only way to tell them apart is to hand it
# directives it must flag and check that it does.
#
# That is not hypothetical. The matcher spells its word boundaries as explicit
# character classes because `\b` is a GNU extension, not POSIX ERE: an
# implementation that treats `\b` as a literal `b` makes the whole pattern
# unmatchable, so the lint reports PASS on a tree full of gates -- silently,
# and only on the machines whose tools differ from CI's. A single wrong
# character in the pattern reproduces that by hand.
#
# Fixtures are generated here rather than committed: each case is one line, and
# a file of one-line fixtures away from the expectations that name them is
# harder to read than the table below.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."

lint="test/regression/ndebug_gate_lint.sh"
fixture_root="$(mktemp -d)"
trap 'rm -rf "$fixture_root"' EXIT

failures=0

# Run the lint over a fixture directory holding a single source line, and
# check the verdict. `expected` is FLAG (the lint must reject it, exit 1) or
# PASS (the lint must accept it, exit 0).
check_line() {
    local description="$1" expected="$2" content="$3"
    local dir="$fixture_root/case"

    rm -rf "$dir"
    mkdir -p "$dir"
    printf '%s\n' "$content" > "$dir/fixture.cpp"

    local output status verdict
    output="$(bash "$lint" "$dir" 2>&1)" && status=0 || status=$?
    case "$status" in
        0) verdict=PASS ;;
        1) verdict=FLAG ;;
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

echo "== gates the lint must reject =="
check_line "plain #ifndef"              FLAG '#ifndef NDEBUG'
check_line "plain #ifdef"               FLAG '#ifdef NDEBUG'
check_line "#if !defined()"             FLAG '#if !defined(NDEBUG)'
check_line "space between # and elif"   FLAG '#  elif defined NDEBUG'
check_line "leading tab"                FLAG $'\t#ifndef NDEBUG'
check_line "NDEBUG last on the line"    FLAG '#if defined(FOO) && !defined NDEBUG'
# The compiler splices a backslash-continued line before it sees a directive,
# so this is one `#if defined(NDEBUG)`. A physical-line matcher misses it.
check_line "backslash continuation"     FLAG '#if \
defined(NDEBUG)'
check_line "continuation mid-expression" FLAG '#if defined(FOO) || \
    defined(NDEBUG)'

echo "== non-gates the lint must accept =="
# bntseq.cpp documents the historical bug in exactly this shape; a matcher that
# flags prose makes the lint unusable in the file that motivated it.
check_line "NDEBUG in block comment"    PASS ' * this was once gated on #ifndef NDEBUG and never compiled out'
check_line "NDEBUG in line comment"     PASS '// see the #ifndef NDEBUG note above'
check_line "unrelated directive"        PASS '#define NDEBUG_FOO 1'
check_line "NDEBUG as a name prefix"    PASS '#ifndef NDEBUG_TRACE'
check_line "NDEBUG as a name suffix"    PASS '#ifdef MY_NDEBUG'
check_line "#endif trailing comment"    PASS '#endif // NDEBUG'
check_line "ordinary code"              PASS 'int main(void) { return 0; }'

echo "== a lint that scanned nothing must not report PASS =="
empty_dir="$fixture_root/empty"
mkdir -p "$empty_dir"
if bash "$lint" "$empty_dir" >/dev/null 2>&1; then
    echo "  FAIL directory with no sources reported PASS" >&2
    failures=$((failures + 1))
else
    echo "  ok   directory with no sources -> FAIL"
fi

if bash "$lint" "$fixture_root/does-not-exist" >/dev/null 2>&1; then
    echo "  FAIL missing directory reported PASS" >&2
    failures=$((failures + 1))
else
    echo "  ok   missing directory -> FAIL"
fi

# Exit 2 rather than 1, so a caller that mis-invokes the lint cannot read the
# result as "gates found" -- or, worse, as a clean run.
bash "$lint" one two >/dev/null 2>&1 && usage_status=0 || usage_status=$?
if (( usage_status == 2 )); then
    echo "  ok   too many arguments -> usage error (exit 2)"
else
    echo "  FAIL too many arguments: expected exit 2, got $usage_status" >&2
    failures=$((failures + 1))
fi

# Build artifacts share src/ with the sources in a working tree. Feeding a .o
# to a text scanner is an error, not a finding, so they must be skipped rather
# than aborting the run.
artifact_dir="$fixture_root/artifacts"
mkdir -p "$artifact_dir"
printf '%s\n' 'int main(void) { return 0; }' > "$artifact_dir/fixture.cpp"
head -c 4096 /dev/urandom > "$artifact_dir/fixture.o"
if bash "$lint" "$artifact_dir" >/dev/null 2>&1; then
    echo "  ok   binary artifact alongside sources -> PASS"
else
    echo "  FAIL binary artifact alongside sources aborted the lint" >&2
    failures=$((failures + 1))
fi

if (( failures > 0 )); then
    echo "FAIL: ndebug_gate_lint_selftest ($failures case(s))" >&2
    exit 1
fi

echo "PASS: ndebug_gate_lint_selftest"
