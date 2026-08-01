#!/usr/bin/env bash
# test/regression/shell_lint_selftest.sh
#
# Regression: shell_lint.sh still rejects the things it is supposed to reject.
#
# The lint it tests has exactly one observable behaviour on a healthy tree --
# it prints PASS -- and that is also what it would print if it had stopped
# checking anything: if the file discovery returned nothing, if a tool went
# missing and the SKIP path swallowed it, if a flag were dropped. From CI those
# states are indistinguishable, and they are the same "green while doing no
# work" failure that ndebug_gate_lint_selftest.sh exists to rule out for the
# lint next door. The only way to tell them apart is to hand the lint scripts
# it must flag and check that it does.
#
# The SKIP path makes this sharper than it looks. shell_lint.sh exits 0 when
# either tool is absent, on purpose, so that `make test` works on a box without
# them -- which means a mis-wired CI job gets a passing gate forever.
# The SHELL_LINT_REQUIRED check below is what makes that impossible, so it is
# tested here rather than assumed.
#
# Fixtures are generated rather than committed: each case is a few lines, and a
# tracked fixture would itself be linted by the very gate under test.
#
# File-scoped rather than per-call: the fixtures below are shell source held in
# single-quoted strings, so almost every one of them contains a `$` that must
# NOT expand here -- expanding it in this script is precisely the bug that
# would stop the fixture from reproducing the finding it is named after.
# shellcheck disable=SC2016
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."

lint="test/regression/shell_lint.sh"
fixture_root="$(mktemp -d)"
trap 'rm -rf "$fixture_root"' EXIT

failures=0

# The selftest asserts on rejection, which requires the tools to be present.
# If they are not, skip in the same shape shell_lint.sh does -- asserting that
# a lint rejects things is meaningless when the lint is itself skipping.
if ! command -v shellcheck > /dev/null 2>&1 || ! command -v shfmt > /dev/null 2>&1; then
    if [[ -n ${SHELL_LINT_REQUIRED:-} ]]; then
        echo "FAIL: SHELL_LINT_REQUIRED is set but shellcheck/shfmt are missing" >&2
        exit 1
    fi
    echo "SKIP: shell_lint_selftest -- shellcheck and/or shfmt not installed"
    exit 0
fi

# Run the lint over a fixture directory holding a single script and check the
# verdict. `expected` is FLAG (the lint must reject it, exit 1) or PASS.
#
# A rejection must also SAY something. The exit status alone is too weak an
# assertion: `unformatted=$(shfmt -l ...)` under `set -e` once aborted this
# lint with exit 1 and no output whatsoever, which every exit-code-only case
# here scored as a correct FLAG. A gate that fails a pull request without
# naming the file or the fix is barely better than no gate, so the third
# argument is a pattern the output must contain.
check_script() {
    local description="$1" expected="$2" content="$3" want_output="${4:-}"
    local dir="$fixture_root/case"

    rm -rf "$dir"
    mkdir -p "$dir"
    printf '%s\n' "$content" > "$dir/fixture.sh"

    local output status verdict
    output="$(bash "$lint" "$dir" 2>&1)" && status=0 || status=$?
    case "$status" in
        0) verdict=PASS ;;
        1) verdict=FLAG ;;
        *) verdict="ERROR(exit $status)" ;;
    esac

    if [[ "$verdict" != "$expected" ]]; then
        echo "  FAIL $description: expected $expected, got $verdict" >&2
        printf '%s\n' "$output" | sed 's/^/       /' >&2
        failures=$((failures + 1))
        return
    fi

    if [[ -n $want_output && $output != *"$want_output"* ]]; then
        echo "  FAIL $description: $verdict as expected, but the output never" >&2
        echo "       mentioned '$want_output' -- a silent verdict is not a verdict" >&2
        printf '%s\n' "$output" | sed 's/^/       /' >&2
        failures=$((failures + 1))
        return
    fi

    echo "  ok   $description -> $verdict"
}

echo "== shellcheck findings the lint must reject =="
# SC2086: the classic. An unquoted expansion that word-splits on a path with a
# space is the single most common way a test script silently does the wrong
# thing rather than failing.
check_script "unquoted expansion (SC2086)" FLAG '#!/usr/bin/env bash
f="a b"
cat $f' 'SC2086'
# SC2164: `cd` that can fail, after which every subsequent path is resolved
# against the wrong directory -- a test that then "passes" is meaningless.
check_script "unchecked cd (SC2164)" FLAG '#!/usr/bin/env bash
cd /nonexistent
echo hi' 'SC2164'
# SC2155: the finding this tree actually had, in index_alt_sidecar_warn_test.sh.
check_script "local masks exit status (SC2155)" FLAG '#!/usr/bin/env bash
f() {
    local x="$(false)"
    echo "$x"
}
f' 'SC2155'

echo "== formatting the lint must reject =="
# The output must name the offending file and the --fix command, not merely
# exit non-zero: see the note on check_script above.
check_script "wrong indent width" FLAG '#!/usr/bin/env bash
if true; then
  echo hi
fi' 'fixture.sh'
check_script "wrong indent names the fix" FLAG '#!/usr/bin/env bash
if true; then
  echo hi
fi' '--fix'
check_script "missing space after redirect" FLAG '#!/usr/bin/env bash
echo hi >"/dev/null"' 'fixture.sh'

echo "== clean scripts the lint must accept =="
check_script "clean and formatted" PASS '#!/usr/bin/env bash
set -euo pipefail
f="a b"
cat -- "$f"' 'PASS: shell_lint'
check_script "scoped disable is honoured" PASS '#!/usr/bin/env bash
set -euo pipefail
# shellcheck disable=SC2016  # single quotes are deliberate
echo '"'"'$notavariable'"'"'' 'PASS: shell_lint'

echo "== the version note fires on real skew, not on spelling =="
# The upstream release binary reports "v3.13.1"; the Homebrew build reports
# "3.13.1". They are the same release, and `brew install shfmt` is what the
# contributing guide tells people to run -- so the unprefixed spelling is the
# common case, and warning on it would train contributors to ignore a note
# that is supposed to mean something.
#
# The pin is read out of the lint rather than repeated here, so bumping the
# version in one place does not quietly stop testing this.
pinned_shfmt="$(mawk -F'"' '/^readonly SHFMT_VERSION=/ {print $2}' "$lint")"
real_shfmt="$(command -v shfmt)"
version_bin="$fixture_root/verbin"
version_dir="$fixture_root/verclean"
mkdir -p "$version_bin" "$version_dir"
printf '%s\n' '#!/usr/bin/env bash' 'set -euo pipefail' 'echo hi' > "$version_dir/fixture.sh"

# $1 = description, $2 = what the stub reports for --version, $3 = want-note (yes|no)
check_version_note() {
    local description="$1" reported="$2" want_note="$3" output
    # The real shfmt path reaches the stub through the environment rather than
    # being interpolated into it: a path containing whitespace would otherwise
    # word-split at exec, so the stub would run the wrong command and every
    # version case would fail on the clean fixture instead of on the note.
    printf '%s\n' \
        '#!/usr/bin/env bash' \
        "if [[ \${1:-} == --version ]]; then echo '$reported'; exit 0; fi" \
        'exec "$REAL_SHFMT" "$@"' > "$version_bin/shfmt"
    chmod +x "$version_bin/shfmt"

    output="$(REAL_SHFMT="$real_shfmt" PATH="$version_bin:$PATH" bash "$lint" "$version_dir" 2>&1)" || {
        echo "  FAIL $description: the lint did not pass on a clean fixture" >&2
        printf '%s\n' "$output" | sed 's/^/       /' >&2
        failures=$((failures + 1))
        return
    }

    if [[ $output == *"note: shfmt"* && $want_note == no ]]; then
        echo "  FAIL $description: warned about CI skew on a matching version" >&2
        printf '%s\n' "$output" | sed 's/^/       /' >&2
        failures=$((failures + 1))
    elif [[ $output != *"note: shfmt"* && $want_note == yes ]]; then
        echo "  FAIL $description: a genuinely different version went unreported" >&2
        failures=$((failures + 1))
    else
        echo "  ok   $description"
    fi
}

check_version_note "pinned spelling ($pinned_shfmt)" "$pinned_shfmt" no
check_version_note "brew spelling (${pinned_shfmt#v})" "${pinned_shfmt#v}" no
check_version_note "a genuinely older version still warns" "v0.0.1" yes
rm -rf "$version_bin"

echo "== a file shfmt cannot parse must fail loudly, not as a nit =="
# shfmt returns 1 for "differs", "cannot parse" and "no such file" alike, so a
# gate that keys on the exit status alone would file a syntax error as a
# formatting nit and tell the author to run --fix, which cannot help.
check_script "unparseable script" FLAG '#!/usr/bin/env bash
if true; then
echo unterminated' 'fixture.sh'

# As check_script(), but for a run that is expected to fail without any
# fixture script to point at. Asserts on the message and not only on the exit
# status, for the reason given above check_script: an abort with exit 1 and no
# output is what these two cases exist to catch, and it scores as a correct
# rejection under an exit-code-only assertion.
check_rejects() { # $1 = description, $2 = directory argument, $3 = expected message
    local description="$1" dir="$2" want_output="$3"
    local output status

    output="$(bash "$lint" "$dir" 2>&1)" && status=0 || status=$?
    if ((status == 0)); then
        echo "  FAIL $description: reported PASS" >&2
        failures=$((failures + 1))
    elif [[ $output != *"$want_output"* ]]; then
        echo "  FAIL $description: exited $status, but the output never" >&2
        echo "       mentioned '$want_output' -- a silent verdict is not a verdict" >&2
        printf '%s\n' "$output" | sed 's/^/       /' >&2
        failures=$((failures + 1))
    else
        echo "  ok   $description -> FAIL"
    fi
}

echo "== a lint that checked nothing must not report PASS =="
empty_dir="$fixture_root/empty"
mkdir -p "$empty_dir"
check_rejects "directory with no scripts" "$empty_dir" "no shell scripts found"
check_rejects "missing directory" "$fixture_root/does-not-exist" "no such directory"

echo "== misuse must not read as a clean run =="
# Exit 2 rather than 1, so a caller that mis-invokes the lint cannot mistake
# the result for "findings" or for success.
bash "$lint" one two > /dev/null 2>&1 && usage_status=0 || usage_status=$?
if ((usage_status == 2)); then
    echo "  ok   two directories -> usage error (exit 2)"
else
    echo "  FAIL two directories: expected exit 2, got $usage_status" >&2
    failures=$((failures + 1))
fi

bash "$lint" --bogus > /dev/null 2>&1 && usage_status=0 || usage_status=$?
if ((usage_status == 2)); then
    echo "  ok   unknown option -> usage error (exit 2)"
else
    echo "  FAIL unknown option: expected exit 2, got $usage_status" >&2
    failures=$((failures + 1))
fi

echo "== SHELL_LINT_REQUIRED turns a missing tool into a failure =="
# The SKIP path is what keeps `make test` working without the tools installed,
# and it is also what would turn a mis-wired CI job into a permanent silent
# pass. Simulate absence with a PATH holding neither tool, and check both
# halves of the contract.
#
# Both cases set SHELL_LINT_REQUIRED explicitly -- `env -u` for the unset half
# -- rather than inheriting whatever the caller has. The CI job sets it at the
# job level, so it is exported into this script too: reading the ambient value
# made the "unset" case assert the opposite of what it says, and it passed
# locally and failed in CI for that reason alone. A test of an environment
# contract has to own that variable.
tool_free_bin="$fixture_root/nobin"
mkdir -p "$tool_free_bin"
# dirname/pwd are here because shell_lint.sh resolves its own repo root with
# them before it checks for anything else; without them the script fails for a
# reason that has nothing to do with the tools being absent.
for cmd in bash env git awk sed cat find mktemp rm mkdir printf dirname pwd; do
    target="$(command -v "$cmd" 2> /dev/null)" && ln -sf "$target" "$tool_free_bin/$cmd"
done
clean_dir="$fixture_root/clean"
mkdir -p "$clean_dir"
printf '%s\n' '#!/usr/bin/env bash' 'echo hi' > "$clean_dir/fixture.sh"

if PATH="$tool_free_bin" env -u SHELL_LINT_REQUIRED bash "$lint" "$clean_dir" > /dev/null 2>&1; then
    echo "  ok   tools absent, unset       -> SKIP (exit 0)"
else
    echo "  FAIL tools absent should SKIP so 'make test' still works" >&2
    failures=$((failures + 1))
fi

if PATH="$tool_free_bin" SHELL_LINT_REQUIRED=1 bash "$lint" "$clean_dir" > /dev/null 2>&1; then
    echo "  FAIL tools absent with SHELL_LINT_REQUIRED=1 reported success" >&2
    failures=$((failures + 1))
else
    echo "  ok   tools absent, REQUIRED=1  -> FAIL (exit 1)"
fi

echo "== --fix repairs formatting =="
fix_dir="$fixture_root/fix"
mkdir -p "$fix_dir"
printf '%s\n' '#!/usr/bin/env bash' 'if true; then' '  echo hi' 'fi' > "$fix_dir/fixture.sh"
if bash "$lint" "$fix_dir" > /dev/null 2>&1; then
    echo "  FAIL unformatted fixture was accepted before --fix" >&2
    failures=$((failures + 1))
elif ! bash "$lint" --fix "$fix_dir" > /dev/null 2>&1; then
    echo "  FAIL --fix did not resolve the formatting" >&2
    failures=$((failures + 1))
elif ! bash "$lint" "$fix_dir" > /dev/null 2>&1; then
    echo "  FAIL tree still unformatted after --fix" >&2
    failures=$((failures + 1))
else
    echo "  ok   unformatted -> --fix -> PASS"
fi

# --fix formats, but shellcheck findings need a human; it must not exit 0 and
# leave the impression that it fixed everything.
fixable_dir="$fixture_root/fixable"
mkdir -p "$fixable_dir"
printf '%s\n' '#!/usr/bin/env bash' 'f="a b"' 'cat $f' > "$fixable_dir/fixture.sh"
if bash "$lint" --fix "$fixable_dir" > /dev/null 2>&1; then
    echo "  FAIL --fix reported success despite a shellcheck finding" >&2
    failures=$((failures + 1))
else
    echo "  ok   --fix with a shellcheck finding -> FAIL"
fi

if ((failures > 0)); then
    echo "FAIL: shell_lint_selftest ($failures case(s))" >&2
    exit 1
fi

echo "PASS: shell_lint_selftest"
