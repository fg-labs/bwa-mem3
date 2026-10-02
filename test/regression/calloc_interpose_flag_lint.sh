#!/bin/bash
# Every test source that defines its own calloc must be compiled with
# $(CALLOC_INTERPOSE_CXXFLAGS).
#
# A test that interposes calloc (to make calloc(0, ...) return NULL, say)
# implements the ordinary case as malloc + memset. GCC at -O2 and above folds
# that pair back into a call to calloc -- into the replacement itself -- so the
# test spins forever instead of failing. -fno-builtin-malloc stops the fold; the
# Makefile carries it as CALLOC_INTERPOSE_CXXFLAGS (see the comment there).
#
# Nothing else catches a missing flag. The CI gcc lane uses Ubuntu's g++, whose
# default _FORTIFY_SOURCE=3 turns memset into __memset_chk and so dodges the fold
# by accident, and clang never folds inside a function named calloc. The hang
# only shows up on an upstream GCC, after the fact, as a `make test` that never
# returns. So this checks the wiring statically instead: for each test/ source
# that defines calloc, its object must have an explicit Makefile rule (the
# .cpp.o suffix rule does not pass the flag) whose recipe names the variable.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

if (($# > 1)); then
    echo "usage: ${BASH_SOURCE[0]##*/} [repo-root-to-check]" >&2
    exit 2
fi

# Optional argument so calloc_interpose_flag_lint_selftest.sh can aim the same
# matcher at a fixture tree.
if (($# == 1)); then
    if [[ ! -d $1 ]]; then
        echo "FAIL: '$1' is not a directory -- nothing was checked" >&2
        exit 1
    fi
    repo_root="$(cd "$1" && pwd)"
fi
cd "$repo_root"

if [[ ! -f Makefile || ! -d test ]]; then
    echo "FAIL: no Makefile or test/ under $repo_root -- nothing was checked" >&2
    exit 1
fi

if ! grep -Eq '^CALLOC_INTERPOSE_CXXFLAGS[[:space:]]*:?=.*-fno-builtin-malloc' Makefile; then
    echo "FAIL: Makefile does not define CALLOC_INTERPOSE_CXXFLAGS with -fno-builtin-malloc" >&2
    exit 1
fi

# A calloc *definition*: `[extern "C"] void *calloc(` on a line that is not a
# bare declaration ending in `;`.
definers=()
while IFS= read -r file; do
    definers+=("$file")
done < <(grep -rlE --include='*.c' --include='*.cc' --include='*.cpp' --include='*.cxx' \
    '^[[:space:]]*(extern[[:space:]]+"C"[[:space:]]+)?void[[:space:]]*\*[[:space:]]*calloc[[:space:]]*\([^;]*$' \
    test | sort)

# Scanning nothing must not read as a PASS.
if ((${#definers[@]} == 0)); then
    echo "FAIL: no test source defines calloc -- the matcher found nothing to check" >&2
    exit 1
fi

failures=0
for src in "${definers[@]}"; do
    obj="${src%.*}.o"
    # Recipe lines of the explicit `<obj>:` rule: the tab-indented lines that
    # follow it, up to the first line that is not one.
    recipe="$(mawk -v obj="$obj" '
        index($0, obj ":") == 1 { in_rule = 1; found = 1; next }
        in_rule && /^\t/        { print; next }
        in_rule                 { in_rule = 0 }
        END                     { if (!found) exit 3 }
    ' Makefile)" || {
        echo "FAIL: $src defines calloc but $obj has no explicit Makefile rule;" \
            "the .cpp.o suffix rule does not pass \$(CALLOC_INTERPOSE_CXXFLAGS)" >&2
        failures=$((failures + 1))
        continue
    }
    # shellcheck disable=SC2016 # a literal Makefile variable reference, not shell
    if ! grep -qF '$(CALLOC_INTERPOSE_CXXFLAGS)' <<< "$recipe"; then
        echo "FAIL: $src defines calloc but the recipe for $obj does not pass" \
            "\$(CALLOC_INTERPOSE_CXXFLAGS); GCC -O2+ folds its malloc + memset" \
            "into a call to itself" >&2
        failures=$((failures + 1))
    fi
done

if ((failures > 0)); then
    exit 1
fi
echo "PASS: all ${#definers[@]} calloc-defining test sources build with \$(CALLOC_INTERPOSE_CXXFLAGS)"
