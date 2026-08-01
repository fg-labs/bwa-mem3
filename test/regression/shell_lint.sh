#!/usr/bin/env bash
# test/regression/shell_lint.sh
#
# Lint and format-check every tracked shell script: shellcheck for
# correctness, shfmt for layout.
#
# Shell is a load-bearing part of the test surface here -- test/regression/
# alone is ~35 scripts, and they are what stands between a kernel change and a
# silent correctness regression. Until this lint existed nothing checked them,
# so the usual shell footguns (unquoted expansions, `local x=$(...)` swallowing
# an exit status, variables that are assigned and never read) could land
# unnoticed in the very code whose job is to notice things.
#
# Two tools, because they catch disjoint classes of problem:
#
#   - shellcheck, for correctness. Run at DEFAULT severity, which on this tree
#     is a quiet ~15-finding baseline that was cleared before the gate went in.
#     `--enable=all` was measured and rejected: 3214 findings, 2875 of them
#     SC2250 ("use ${braces}"), a whole-tree rewrite for no correctness gain.
#   - shfmt, for layout, so review comments are about behaviour rather than
#     indentation. See SHFMT_FLAGS below for why each flag is there.
#
# (Both list entries start with `- ` because a comment line whose first word is
# "shellcheck" is parsed as a directive, not prose -- SC1072/SC1073. This gate
# caught that in its own header the first time it ran.)
#
# Both tools are PINNED (see the versions below and the CI job that installs
# them). This is not incidental: a new upstream release that adds a check would
# otherwise turn CI red on a pull request that never touched shell. The gate
# should only get stricter when someone deliberately bumps the pin.
#
# Behaviour when the tools are missing depends on where this runs:
#
#   locally   SKIP with a visible notice and exit 0, so `make test` does not
#             break for a contributor who has not installed them.
#   in CI     SHELL_LINT_REQUIRED=1 makes a missing tool a hard error. Without
#             that, a broken install step would turn the gate into a permanent
#             silent pass -- the same "green while doing no work" failure that
#             ndebug_gate_lint_selftest.sh exists to rule out next door.
#
# Usage:
#   shell_lint.sh [--fix] [directory-to-scan]
#
#   --fix   rewrite files with shfmt instead of only reporting, then report
#           whatever findings remain (those need a human).
#
# Environment:
#   SHELL_LINT_REQUIRED=1   missing tools are a hard error rather than a SKIP.

set -euo pipefail

# Pinned versions. Keep in lockstep with the installer in
# .github/workflows/ci.yml; a mismatch between them is what makes CI and a
# local run disagree about the same tree.
readonly SHELLCHECK_VERSION="0.11.0"
readonly SHFMT_VERSION="v3.13.1"

# Chosen by measuring churn against this tree rather than by preference:
#
#   -i 4  matches the existing indent ladder (1439 lines at indent 4, then
#         8/12/16/20). -i 2 was measured at roughly twice the churn.
#   -ci   indent switch cases.
#   -bn   binary operators lead the continuation line, matching the existing
#         `\` + `| mawk` style used throughout test/regression/.
#   -sr   space after redirect operators; the tree already leans this way.
#
# `-s` (--simplify) is deliberately absent. It rewrites code rather than
# layout, and a formatter that can change behaviour is not one you want
# wired to a gate that people will run with --fix.
readonly SHFMT_FLAGS=(-i 4 -ci -bn -sr)

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

fix=0
scan_root=""
while (($# > 0)); do
    case "$1" in
        --fix) fix=1 ;;
        -h | --help)
            echo "usage: ${BASH_SOURCE[0]##*/} [--fix] [directory-to-scan]"
            exit 0
            ;;
        -*)
            echo "FAIL: unknown option '$1'" >&2
            echo "usage: ${BASH_SOURCE[0]##*/} [--fix] [directory-to-scan]" >&2
            exit 2
            ;;
        *)
            # Exit 2 rather than 1 so a caller that mis-invokes the lint cannot
            # read the result as "findings" -- or, worse, as a clean run.
            if [[ -n $scan_root ]]; then
                echo "FAIL: at most one directory may be given" >&2
                echo "usage: ${BASH_SOURCE[0]##*/} [--fix] [directory-to-scan]" >&2
                exit 2
            fi
            scan_root="$1"
            ;;
    esac
    shift
done

# ---------------------------------------------------------------------------
# Tool discovery
# ---------------------------------------------------------------------------

# Report every missing tool at once. Telling someone to install shellcheck,
# watching them do it, and only then mentioning shfmt is a poor trade for the
# three lines it saves here.
missing=()
command -v shellcheck > /dev/null 2>&1 || missing+=("shellcheck (pinned: $SHELLCHECK_VERSION)")
command -v shfmt > /dev/null 2>&1 || missing+=("shfmt (pinned: $SHFMT_VERSION)")

if ((${#missing[@]} > 0)); then
    if [[ -n ${SHELL_LINT_REQUIRED:-} ]]; then
        echo "FAIL: SHELL_LINT_REQUIRED is set but these tools are missing:" >&2
        printf '  - %s\n' "${missing[@]}" >&2
        echo >&2
        echo "This is set by CI, where a missing tool means the install step" >&2
        echo "broke. Skipping here would make the gate a silent pass." >&2
        exit 1
    fi
    echo "SKIP: shell_lint -- not installed:"
    printf '  - %s\n' "${missing[@]}"
    echo "  (install both, or set SHELL_LINT_REQUIRED=1 to make this a failure)"
    exit 0
fi

# A version skew changes what the tools report, so a local run can disagree
# with CI over an unmodified tree. Warn rather than fail: being unable to lint
# at all is worse than linting with a near-enough version.
#
# Compared with any leading `v` stripped from BOTH sides. The upstream release
# binary CI installs reports "v3.13.1", but the Homebrew build reports
# "3.13.1" -- and `brew install shfmt` is what the contributing guide tells
# people to run, so without this the documented install path warns about CI
# skew on an exactly-correct version. A note that fires when nothing is wrong
# is one people learn to ignore.
actual_shellcheck="$(shellcheck --version | awk '/^version:/ {print $2}')"
actual_shfmt="$(shfmt --version)"
if [[ ${actual_shellcheck#v} != "${SHELLCHECK_VERSION#v}" ]]; then
    echo "note: shellcheck $actual_shellcheck, pinned $SHELLCHECK_VERSION -- results may differ from CI" >&2
fi
if [[ ${actual_shfmt#v} != "${SHFMT_VERSION#v}" ]]; then
    echo "note: shfmt $actual_shfmt, pinned $SHFMT_VERSION -- results may differ from CI" >&2
fi

# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

# Default is every TRACKED *.sh, so a new script is covered the moment it is
# added -- there is no opt-in list to forget to update. `git ls-files` rather
# than `find` because a working tree here also holds untracked scratch scripts
# and sibling worktrees, and linting those would make the result depend on
# whatever the developer happens to have lying around.
#
# The optional directory argument exists so shell_lint_selftest.sh can aim the
# same checks at a fixture tree of known-good and known-bad scripts, which is
# the only way to tell "the tree is clean" apart from "the lint stopped
# checking". Fixtures are not tracked, hence find(1) in that mode.
scripts=()
if [[ -z $scan_root ]]; then
    while IFS= read -r -d '' file; do
        scripts+=("$file")
    done < <(git ls-files -z -- '*.sh')
    scan_label="tracked *.sh"
else
    if [[ ! -d $scan_root ]]; then
        echo "FAIL: no such directory: $scan_root -- nothing was linted" >&2
        exit 1
    fi
    while IFS= read -r -d '' file; do
        scripts+=("$file")
    done < <(find "$scan_root" -type f -name '*.sh' -print0)
    scan_label="$scan_root"
fi

# Every way this lint can end up checking nothing has to be a failure, not a
# PASS. A lint that silently checks nothing reads as green while doing no work,
# which is precisely the state it is supposed to make impossible.
if ((${#scripts[@]} == 0)); then
    echo "FAIL: no shell scripts found under $scan_label -- nothing was linted" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# --fix
# ---------------------------------------------------------------------------

if ((fix)); then
    # Same exit-status ambiguity as the check path below: shfmt -w returns 1 on
    # a file it cannot parse, which `set -e` would turn into a silent abort
    # part-way through rewriting the tree.
    if ! shfmt -w "${SHFMT_FLAGS[@]}" "${scripts[@]}"; then
        echo "FAIL: shfmt could not rewrite every file; see the errors above" >&2
        exit 1
    fi
    echo "shfmt: rewrote ${#scripts[@]} file(s) under $scan_label"
    if shellcheck -s bash -x -f gcc "${scripts[@]}"; then
        echo "PASS: shell_lint --fix (formatting applied, no shellcheck findings)"
        exit 0
    fi
    echo >&2
    echo "Formatting is fixed, but the shellcheck findings above need a human:" >&2
    echo "either correct the code, or add a scoped '# shellcheck disable=SCxxxx'" >&2
    echo "with a comment saying why the diagnostic does not apply." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Check
# ---------------------------------------------------------------------------

status=0

# -x follows `source`d files so a helper's definitions are in scope for the
# script that sources it; without it every such variable reads as undefined.
if ! shellcheck -s bash -x -f gcc "${scripts[@]}"; then
    echo >&2
    echo "shellcheck found the issue(s) above in $scan_label." >&2
    echo "Fix them, or add a scoped '# shellcheck disable=SCxxxx' with a" >&2
    echo "comment saying why the diagnostic does not apply here." >&2
    status=1
fi

# shfmt exits 1 for "files differ", for "could not parse a file", and for
# "no such file" alike, so the exit status alone cannot tell a finding from a
# failure -- and `unformatted=$(shfmt -l ...)` under `set -e` would abort the
# script with no message at all. Split the streams instead: -l names differing
# files on stdout and reports real errors on stderr. A file shfmt cannot parse
# has to fail loudly rather than be filed as a formatting nit.
shfmt_stderr="$(mktemp)"
trap 'rm -f "$shfmt_stderr"' EXIT
unformatted="$(shfmt -l "${SHFMT_FLAGS[@]}" "${scripts[@]}" 2> "$shfmt_stderr")" || true
if [[ -s $shfmt_stderr ]]; then
    echo "FAIL: shfmt did not complete, so the format check is not trustworthy:" >&2
    cat "$shfmt_stderr" >&2
    exit 1
fi

if [[ -n $unformatted ]]; then
    if ((status != 0)); then echo >&2; fi
    echo "FAIL: not shfmt-formatted:" >&2
    while IFS= read -r file; do
        printf '  %s\n' "$file" >&2
    done <<< "$unformatted"
    echo >&2
    echo "Run:  ./test/regression/shell_lint.sh --fix" >&2
    status=1
fi

if ((status != 0)); then
    exit 1
fi

echo "PASS: shell_lint (${#scripts[@]} script(s) in $scan_label: shellcheck clean, shfmt-formatted)"
