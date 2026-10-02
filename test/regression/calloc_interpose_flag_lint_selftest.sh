#!/usr/bin/env bash
# test/regression/calloc_interpose_flag_lint_selftest.sh
#
# Regression: calloc_interpose_flag_lint.sh still rejects a calloc-defining test
# whose object is built without $(CALLOC_INTERPOSE_CXXFLAGS).
#
# On a healthy tree the lint prints PASS, which is also what it would print if
# its matcher stopped recognising calloc definitions or recipes. The only way to
# tell those apart is to hand it fixture trees it must reject and check that it
# does.
# The fixtures are Makefile text, so `$(...)` in single quotes is meant literally.
# shellcheck disable=SC2016
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."

lint="$PWD/test/regression/calloc_interpose_flag_lint.sh"
fixture_root="$(mktemp -d)"
trap 'rm -rf "$fixture_root"' EXIT

failures=0

# Build a fixture tree with one calloc-defining test/t.cpp and the given
# Makefile body, run the lint over it, and check the verdict. `expected` is
# FLAG (the lint must reject it, exit 1) or PASS (it must accept it, exit 0).
check_case() {
    local description="$1" expected="$2" definition="$3" makefile_body="$4"
    local dir="$fixture_root/case"

    rm -rf "$dir"
    mkdir -p "$dir/test"
    printf '%s\n' "$definition" '{ return 0; }' > "$dir/test/t.cpp"
    printf 'CALLOC_INTERPOSE_CXXFLAGS = -fno-builtin-malloc\n\n%b' "$makefile_body" > "$dir/Makefile"

    local rc=0
    "$lint" "$dir" > /dev/null 2>&1 || rc=$?
    local got
    case "$rc" in
        0) got=PASS ;;
        1) got=FLAG ;;
        *) got="exit $rc" ;;
    esac
    if [[ $got != "$expected" ]]; then
        echo "FAIL: $description: expected $expected, got $got" >&2
        failures=$((failures + 1))
    fi
}

ext='extern "C" void *calloc(size_t nmemb, size_t size)'
good_rule='test/t.o: test/t.cpp\n\t$(CXX) -c $(CXXFLAGS) $(CALLOC_INTERPOSE_CXXFLAGS) $< -o $@\n'

check_case "explicit rule passes the flag" PASS "$ext" "$good_rule"
check_case "plain C definition, rule passes the flag" PASS 'void *calloc(size_t n, size_t s)' "$good_rule"
check_case "no explicit rule (suffix rule)" FLAG "$ext" 'other.o: other.cpp\n\t$(CXX) -c $< -o $@\n'
check_case "rule without the flag" FLAG "$ext" 'test/t.o: test/t.cpp\n\t$(CXX) -c $(CXXFLAGS) $< -o $@\n'
check_case "flag only in a later, unrelated recipe" FLAG "$ext" \
    'test/t.o: test/t.cpp\n\t$(CXX) -c $(CXXFLAGS) $< -o $@\n\nx.o: x.cpp\n\t$(CXX) $(CALLOC_INTERPOSE_CXXFLAGS) $< -o $@\n'
# A bare declaration is not a definition, so there is nothing to check -- and
# scanning nothing must fail rather than pass.
check_case "declaration only (nothing to check)" FLAG 'extern "C" void *calloc(size_t, size_t);' "$good_rule"

if ((failures > 0)); then
    exit 1
fi
echo "PASS: calloc_interpose_flag_lint.sh accepts and rejects every fixture as expected"
