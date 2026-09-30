#!/usr/bin/env bash
# test/regression/rescue_docs_lint_selftest.sh
#
# Regression: rescue_docs_lint.sh still detects each drift it exists to catch.
#
# The lint prints PASS when the pages, the knob list and the source agree -- and
# also when an extraction has quietly stopped extracting, because an empty set
# satisfies every check over it. The only way to tell the two apart is to hand
# the lint trees that have drifted and confirm it says so, each case pinning the
# message of the branch it means to exercise (the exit status alone is shared by
# every branch, and by a crash).
#
# Each case builds a minimal fixture tree -- a knob list, a source that reads it,
# the user-facing table, a SUMMARY, one page, a workflow, a unit test and a
# regression script -- then applies one mutation to it.
# The fixtures are Markdown and C, whose backticks and quotes are literal text:
# shellcheck disable=SC2016  # single-quoted backticks here are Markdown code spans
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."

command -v python3 > /dev/null 2>&1 || {
    echo "SKIP: python3 not on PATH (the self-test edits its fixtures with it)"
    exit 0
}

lint="$PWD/test/regression/rescue_docs_lint.sh"
gen="$PWD/scripts/rescue_knobs.sh"
fixture_root="$(mktemp -d)"
trap 'rm -rf "$fixture_root"' EXIT

failures=0

# The page every fixture starts from. Its claim names all three kinds of gate.
readonly PAGE='# Rescue

The filter lives in `src/rescue_x.cpp` and its entry point is `rescue_x_run()`;
`BWA3_RESCUE_ALPHA=0` turns it off.

Output is identical with it on or off. Gates: `rescue_x_identity.sh`,
`rescue x: the filter matches the reference`, `Run the rescue identity`.

```
Code blocks are not prose: identical, `no_such_symbol`.
```'

# Build a fixture tree and print its path. Called in a command substitution, so
# the tree gets a fresh mktemp directory: a counter would not survive the subshell,
# and a reused path would let one case's mutation leak into the next.
build_fixture() {
    local dir
    dir="$(mktemp -d "$fixture_root/case-XXXXXX")"
    mkdir -p "$dir/src" "$dir/docs/src/developer-guide" "$dir/docs/src/whats-different" \
        "$dir/docs/_generated/rescue" "$dir/.github/workflows" "$dir/test/unit" "$dir/test/regression"
    cat > "$dir/src/rescue_env.h" << 'EOF'
/* Knobs. Design notes: docs/src/developer-guide/rescue.md.
 *
 * rescue-knobs:begin
 * BWA3_RESCUE_ALPHA           1     the alpha filter, a meaning that runs
 *                                   onto a second line
 * BWA3_RESCUE_BETA            7     the beta threshold
 * rescue-knobs:end */
EOF
    cat > "$dir/src/rescue_x.cpp" << 'EOF'
/* See docs/src/developer-guide/rescue.md. */
#include "rescue_env.h"
static bool rescue_x_on() { return rescue_env_on("BWA3_RESCUE_ALPHA"); }
static int rescue_x_beta() { return rescue_env_int("BWA3_RESCUE_BETA", 7); }
int rescue_x_run() { return rescue_x_on() ? rescue_x_beta() : 0; }
EOF
    cat > "$dir/docs/src/whats-different/performance.md" << 'EOF'
| Variable | Effect |
|---|---|
| `BWA3_RESCUE_ALPHA=0` | off |
| `BWA3_RESCUE_BETA=<n>` | threshold |
EOF
    printf '# Summary\n\n- [Rescue](developer-guide/rescue.md)\n' > "$dir/docs/src/SUMMARY.md"
    printf '%s\n' "$PAGE" > "$dir/docs/src/developer-guide/rescue.md"
    cat > "$dir/.github/workflows/ci.yml" << 'EOF'
jobs:
  build:
    steps:
      - name: Run the rescue identity
        run: bash test/regression/rescue_x_identity.sh
EOF
    printf 'TEST_CASE("rescue x: the filter matches the reference" * doctest::test_suite("unit/pair")) {}\n' \
        > "$dir/test/unit/test_rescue_x.cpp"
    printf '#!/usr/bin/env bash\necho PASS\n' > "$dir/test/regression/rescue_x_identity.sh"
    bash "$gen" "$dir" > "$dir/docs/_generated/rescue/knobs.md"
    printf '%s' "$dir"
}

# $1 = description, $2 = PASS|FAIL, $3 = fixture dir, $4 = substring the output must contain.
check() {
    local description="$1" expected="$2" dir="$3" expect_text="${4-}" output status verdict
    output="$(bash "$lint" "$dir" 2>&1)" && status=0 || status=$?
    case "$status" in
        0) verdict=PASS ;;
        1) verdict=FAIL ;;
        *) verdict="ERROR(exit $status)" ;;
    esac
    if [[ $verdict == "$expected" ]] && [[ -z $expect_text || $output == *"$expect_text"* ]]; then
        echo "  ok   $description -> $verdict"
    else
        echo "  FAIL $description: expected $expected${expect_text:+ matching \"$expect_text\"}, got $verdict" >&2
        printf '%s\n' "$output" | sed 's/^/       /' >&2
        failures=$((failures + 1))
    fi
}

# Replace the first occurrence of $2 by $3 in file $1 (literal strings).
edit() {
    python3 - "$1" "$2" "$3" << 'PY'
import sys
path, old, new = sys.argv[1:4]
text = open(path).read()
if old not in text:
    sys.exit("edit: %r not found in %s" % (old, path))
open(path, 'w').write(text.replace(old, new, 1))
PY
}

d="$(build_fixture)"
check "an agreeing tree" PASS "$d" "PASS: rescue_docs_lint (1 page(s), 2 knobs"

d="$(build_fixture)"
check "a missing directory" FAIL "$d/nope" "is not a directory"

d="$(build_fixture)"
printf '/* no list */\n' > "$d/src/rescue_env.h"
check "no knob list" FAIL "$d" "could not render the knob list"

# 1. the generated table
d="$(build_fixture)"
edit "$d/src/rescue_env.h" "the beta threshold" "the beta threshold, reworded"
check "a stale generated table" FAIL "$d" "knobs.md is stale"

d="$(build_fixture)"
rm "$d/docs/_generated/rescue/knobs.md"
check "a missing generated table" FAIL "$d" "knobs.md is missing"

# 2. knobs read in src/
d="$(build_fixture)"
printf 'static int g() { return rescue_env_int("BWA3_RESCUE_GAMMA", 3); }\n' > "$d/src/rescue_y.cpp"
printf '/* docs/src/developer-guide/rescue.md */\n' >> "$d/src/rescue_y.cpp"
check "src/ reads an unlisted knob" FAIL "$d" "reads BWA3_RESCUE_GAMMA, which is not in the knob list"

d="$(build_fixture)"
edit "$d/src/rescue_x.cpp" 'rescue_env_int("BWA3_RESCUE_BETA", 7)' '7'
bash "$gen" "$d" > "$d/docs/_generated/rescue/knobs.md"
check "a listed knob nothing reads" FAIL "$d" "names BWA3_RESCUE_BETA, which nothing in src/ reads"

d="$(build_fixture)"
edit "$d/src/rescue_x.cpp" '"BWA3_RESCUE_BETA", 7' '"BWA3_RESCUE_BETA", 9'
check "an integer default that disagrees" FAIL "$d" "says default 7, rescue_env_int reads it with 9"

d="$(build_fixture)"
edit "$d/src/rescue_x.cpp" '"BWA3_RESCUE_BETA", 7' '"BWA3_RESCUE_BETA", -1'
check "a computed default listed as a number" FAIL "$d" "default is computed (-1)"

d="$(build_fixture)"
edit "$d/src/rescue_env.h" "BWA3_RESCUE_ALPHA           1" "BWA3_RESCUE_ALPHA           0"
bash "$gen" "$d" > "$d/docs/_generated/rescue/knobs.md"
check "a toggle listed as off" FAIL "$d" "is read with rescue_env_on (default on)"

d="$(build_fixture)"
edit "$d/src/rescue_x.cpp" 'rescue_env_on("BWA3_RESCUE_ALPHA")' 'rescue_env_opt_in("BWA3_RESCUE_ALPHA")'
check "an opt-in listed as on" FAIL "$d" "is read with rescue_env_opt_in (default off)"

d="$(build_fixture)"
edit "$d/src/rescue_x.cpp" 'rescue_env_on("BWA3_RESCUE_ALPHA")' 'getenv("BWA3_RESCUE_ALPHA") != 0'
check "a bare getenv" FAIL "$d" "read with a bare getenv"

d="$(build_fixture)"
edit "$d/src/rescue_x.cpp" 'rescue_env_on("BWA3_RESCUE_ALPHA")' 'getenv (
    "BWA3_RESCUE_ALPHA") != 0'
check "a bare getenv spaced and split across lines" FAIL "$d" "read with a bare getenv"

d="$(build_fixture)"
edit "$d/src/rescue_x.cpp" 'rescue_env_int("BWA3_RESCUE_BETA", 7)' 'rescue_env_int ( "BWA3_RESCUE_BETA" , 7 )'
check "a spaced reader call is a read" PASS "$d"

d="$(build_fixture)"
printf 'static int g() {\n    return rescue_env_int(\n        "BWA3_RESCUE_GAMMA", 3);\n}\n' > "$d/src/rescue_y.cpp"
printf '/* docs/src/developer-guide/rescue.md */\n' >> "$d/src/rescue_y.cpp"
check "src/ reads an unlisted knob in a call split across lines" FAIL "$d" "reads BWA3_RESCUE_GAMMA, which is not in the knob list"

d="$(build_fixture)"
edit "$d/src/rescue_x.cpp" 'rescue_env_int("BWA3_RESCUE_BETA", 7)' '7 /* rescue_env_int("BWA3_RESCUE_BETA", 7) */'
check "a block-commented reader call is not a read" FAIL "$d" "names BWA3_RESCUE_BETA, which nothing in src/ reads"

d="$(build_fixture)"
edit "$d/src/rescue_x.cpp" 'return rescue_env_int("BWA3_RESCUE_BETA", 7); }' 'return 7; } // rescue_env_int("BWA3_RESCUE_BETA", 7)'
check "a line-commented reader call is not a read" FAIL "$d" "names BWA3_RESCUE_BETA, which nothing in src/ reads"

# 3. the user-facing table
d="$(build_fixture)"
edit "$d/docs/src/whats-different/performance.md" '| `BWA3_RESCUE_BETA=<n>` | threshold |' ''
check "a knob missing from the user-facing table" FAIL "$d" "has no row for BWA3_RESCUE_BETA"

d="$(build_fixture)"
printf '| `BWA3_RESCUE_OLD=0` | gone |\n' >> "$d/docs/src/whats-different/performance.md"
check "a documented knob that is gone" FAIL "$d" "documents BWA3_RESCUE_OLD, which is not in the knob list"

# 4. pages and their names
d="$(build_fixture)"
printf '# Summary\n' > "$d/docs/src/SUMMARY.md"
check "a page missing from the SUMMARY" FAIL "$d" "is not linked from docs/src/SUMMARY.md"

d="$(build_fixture)"
edit "$d/docs/src/developer-guide/rescue.md" '`rescue_x_run()`' '`rescue_x_renamed()`'
check "a renamed function" FAIL "$d" "'rescue_x_renamed' is not a word anywhere"

d="$(build_fixture)"
edit "$d/docs/src/developer-guide/rescue.md" '`rescue_x_run()`' '`ns::rescue_x_gone`'
check "a scoped name with a missing part" FAIL "$d" "'rescue_x_gone' is not a word anywhere"

d="$(build_fixture)"
edit "$d/docs/src/developer-guide/rescue.md" '`src/rescue_x.cpp`' '`src/rescue_z.cpp:12`'
check "a path that does not exist" FAIL "$d" "and src/rescue_z.cpp does not exist"

d="$(build_fixture)"
edit "$d/docs/src/developer-guide/rescue.md" '`src/rescue_x.cpp`' '`src/rescue_x.{h,cpp}`'
check "a path family" FAIL "$d" "names the path family"

d="$(build_fixture)"
edit "$d/docs/src/developer-guide/rescue.md" '`src/rescue_x.cpp`' '`rescue_z.cpp`'
check "a bare file name that does not exist" FAIL "$d" "no file of that name exists"

d="$(build_fixture)"
edit "$d/docs/src/developer-guide/rescue.md" '`BWA3_RESCUE_ALPHA=0`' '`BWA3_RESCUE_DELTA=0`'
check "a page naming an unlisted knob" FAIL "$d" "names BWA3_RESCUE_DELTA, which is not in the knob list"

# 5. links between pages and source
d="$(build_fixture)"
edit "$d/src/rescue_x.cpp" '/* See docs/src/developer-guide/rescue.md. */' '/* no page */'
check "a rescue source naming no page" FAIL "$d" "src/rescue_x.cpp names no rescue developer-guide page"

d="$(build_fixture)"
edit "$d/src/rescue_x.cpp" 'docs/src/developer-guide/rescue.md' 'docs/src/developer-guide/rescue-gone.md'
edit "$d/src/rescue_env.h" 'docs/src/developer-guide/rescue.md' 'the pages'
check "a source naming a missing page" FAIL "$d" "names docs/src/developer-guide/rescue-gone.md, which does not exist"

d="$(build_fixture)"
cp "$d/docs/src/developer-guide/rescue.md" "$d/docs/src/developer-guide/rescue-extra.md"
printf -- '- [Extra](developer-guide/rescue-extra.md)\n' >> "$d/docs/src/SUMMARY.md"
check "a page no source names" FAIL "$d" "no file under src/ names docs/src/developer-guide/rescue-extra.md"

# 6. claims and gates
d="$(build_fixture)"
printf '\nThe result is unchanged whatever the knob says.\n' >> "$d/docs/src/developer-guide/rescue.md"
check "a claim with no gate" FAIL "$d" "makes an exactness claim without naming its gate"

d="$(build_fixture)"
printf '\nBoth kernels give the same outputs.\n' >> "$d/docs/src/developer-guide/rescue.md"
check "a same-outputs claim with no gate" FAIL "$d" "makes an exactness claim without naming its gate"

d="$(build_fixture)"
printf '\n| a claim | byte-for-byte |\n' >> "$d/docs/src/developer-guide/rescue.md"
check "a table-row claim with no gate" FAIL "$d" "makes an exactness claim without naming its gate"

d="$(build_fixture)"
printf '\nEvery shortcut is exact.\n' >> "$d/docs/src/developer-guide/rescue.md"
check "an is-exact claim with no gate" FAIL "$d" "makes an exactness claim without naming its gate"

d="$(build_fixture)"
printf '\nThe two forms are equivalent.\n' >> "$d/docs/src/developer-guide/rescue.md"
check "an equivalence claim with no gate" FAIL "$d" "makes an exactness claim without naming its gate"

d="$(build_fixture)"
printf '\nThe filter counts exact 5-mer hits and the exactness section explains why.\n' >> "$d/docs/src/developer-guide/rescue.md"
check "exact as a plain adjective is no claim" PASS "$d"

d="$(build_fixture)"
printf '\nIt is identical. Gate: `no such test`.\n' >> "$d/docs/src/developer-guide/rescue.md"
check "a gate that names nothing" FAIL "$d" "names the gate 'no such test'"

d="$(build_fixture)"
edit "$d/.github/workflows/ci.yml" 'run: bash test/regression/rescue_x_identity.sh' 'run: true'
check "a regression gate no workflow runs" FAIL "$d" "names the gate 'rescue_x_identity.sh'"

d="$(build_fixture)"
edit "$d/.github/workflows/ci.yml" 'run: bash test/regression/rescue_x_identity.sh' 'run: true # was test/regression/rescue_x_identity.sh'
check "a regression gate only a workflow comment names" FAIL "$d" "names the gate 'rescue_x_identity.sh'"

d="$(build_fixture)"
edit "$d/test/unit/test_rescue_x.cpp" 'TEST_CASE(' '// TEST_CASE('
check "a unit gate whose TEST_CASE is commented out" FAIL "$d" "names the gate 'rescue x: the filter matches the reference'"

d="$(build_fixture)"
edit "$d/test/unit/test_rescue_x.cpp" 'TEST_CASE(' '/*
TEST_CASE('
printf '*/\n' >> "$d/test/unit/test_rescue_x.cpp"
check "a unit gate whose TEST_CASE is block-commented" FAIL "$d" "names the gate 'rescue x: the filter matches the reference'"

d="$(build_fixture)"
edit "$d/test/unit/test_rescue_x.cpp" 'TEST_CASE(' 'TEST_CASE (
    '
check "a unit gate whose TEST_CASE is split across lines" PASS "$d"

d="$(build_fixture)"
edit "$d/docs/src/developer-guide/rescue.md" '`rescue_x_run()`' '`rescue_x_comment_only`'
printf '/* rescue_x_comment_only */\n' >> "$d/src/rescue_x.cpp"
check "a name only a C comment spells" FAIL "$d" "'rescue_x_comment_only' is not a word anywhere"

d="$(build_fixture)"
edit "$d/docs/src/developer-guide/rescue.md" '`rescue_x_run()`' '`rescue_x_script_only`'
printf '# rescue_x_script_only\n' >> "$d/test/regression/rescue_x_identity.sh"
check "a name only a script comment spells" FAIL "$d" "'rescue_x_script_only' is not a word anywhere"

d="$(build_fixture)"
printf '/* "BWA3_RESCUE_GHOST" is only mentioned here */\n' >> "$d/src/rescue_x.cpp"
check "a knob named only in a comment is not a read" PASS "$d"

d="$(build_fixture)"
printf '\nNothing here claims anything; it only mentions `Gate:` in passing.\n' >> "$d/docs/src/developer-guide/rescue.md"
check "a block with no claim" PASS "$d"

if ((failures > 0)); then
    echo "FAIL: rescue_docs_lint_selftest ($failures case(s) did not behave)" >&2
    exit 1
fi
echo "PASS: rescue_docs_lint_selftest (every case behaved)"
