#!/usr/bin/env bash
# test/regression/readme_contract_lint_selftest.sh
#
# Regression: readme_contract_lint.sh still recognises a stale README.
#
# The lint it tests has one observable behaviour on a healthy tree -- it prints
# PASS -- and that is also what it prints when its own matching has stopped
# matching anything. Those two states are indistinguishable from CI, which is
# the same "green while doing no work" failure the lint exists to catch in the
# README. The only way to tell them apart is to hand it trees it must reject
# and check that it does.
#
# That is not hypothetical for this lint in particular. Deciding whether a
# script "reads no environment" is a heuristic over shell text, and it was
# wrong twice while being written: it first read `MAKE="${MAKE:-make}"` as an
# assignment rather than an environment input, and then read the same
# expansion inside its own header comment as live code. Either mistake moves a
# script between the two sets silently.
#
# Fixtures are generated here rather than committed: each case is a two-file
# tree, and a directory of those away from the expectations that name them is
# harder to read than the cases below.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."

lint="test/regression/readme_contract_lint.sh"
fixture_root="$(mktemp -d)"
trap 'rm -rf "$fixture_root"' EXIT

failures=0

readonly BEGIN_MARKER='<!-- source-only lints: begin -->'
readonly END_MARKER='<!-- source-only lints: end -->'

# Build a fixture tree: a README, plus one script per "name:body" pair given
# after it. Returns the tree's path on stdout.
make_tree() {
    local name="$1" readme_body="$2"
    shift 2

    local dir="$fixture_root/$name"
    rm -rf "$dir"
    mkdir -p "$dir/test/regression"
    printf '%s\n' "$readme_body" > "$dir/test/regression/README.md"

    local spec
    for spec in "$@"; do
        printf '%s\n' "${spec#*:}" > "$dir/test/regression/${spec%%:*}"
    done

    printf '%s' "$dir"
}

# Run the lint over a fixture tree and check the verdict. `expected` is FLAG
# (the lint must reject it, exit 1) or PASS (must accept it, exit 0).
check_tree() {
    local description="$1" expected="$2" dir="$3"

    local output status
    output="$(bash "$lint" "$dir" 2>&1)" && status=0 || status=$?

    local verdict
    case "$status" in
        0) verdict=PASS ;;
        1) verdict=FLAG ;;
        *) verdict="EXIT$status" ;;
    esac

    if [[ $verdict == "$expected" ]]; then
        echo "  ok   $description -> $expected"
    else
        echo "  FAIL $description -> got $verdict, wanted $expected" >&2
        printf '       %s\n' "$output" >&2
        failures=1
    fi
}

# A README whose source-only block names exactly the scripts passed in.
readme_naming() {
    local rows="" name
    for name in "$@"; do
        rows+="- \`$name\` takes an optional positional directory"$'\n'
    done
    printf '%s\n\n%s\n%s%s\n' \
        "Every script but the source-only lints reads its inputs from the environment." \
        "$BEGIN_MARKER" "$rows" "$END_MARKER"
}

# The fixture bodies below are shell source held as data, so every `$` in them
# must survive into the fixture file unexpanded. Single quotes are the point,
# not an oversight -- and the lint under test has to read them as data too,
# which is what several of the cases check.
# shellcheck disable=SC2016

# A script that takes the named environment variable, and one that takes none.
readonly ENV_SCRIPT='set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"
echo "PASS: fixture"'
# shellcheck disable=SC2016
readonly SOURCE_ONLY_SCRIPT='set -euo pipefail
scan_root="${1:-src}"
echo "PASS: fixture $scan_root"'
# shellcheck disable=SC2016
readonly DEFAULTED_ENV_SCRIPT='set -euo pipefail
MAKE="${MAKE:-make}"
"$MAKE" --version'

# A script that holds shell text as data -- the shape of the self-tests next to
# this one, which build fixture scripts from single-quoted constants spanning
# several lines. Everything between the inner quotes is a string; nothing in it
# is an environment read. Written with the `'"'"'` idiom so the inner quotes
# survive into the fixture file.
# shellcheck disable=SC2016
readonly QUOTED_DATA_SCRIPT='set -euo pipefail
payload='"'"'set -euo pipefail
MAKE="${MAKE:-make}"
: "${BWA_MEM3:?BWA_MEM3 must be set}"'"'"'
printf '"'"'%s\n'"'"' "$payload" > /dev/null
echo "PASS: fixture"'

echo "== stale READMEs the lint must reject =="

check_tree "README names a script that was deleted" FLAG \
    "$(make_tree deleted "$(readme_naming lint.sh)

See also \`retired.sh\` for the older check." \
        "lint.sh:$SOURCE_ONLY_SCRIPT")"

check_tree "source-only script missing from the block" FLAG \
    "$(make_tree unlisted "$(readme_naming lint.sh)

Also here: \`other_lint.sh\`." \
        "lint.sh:$SOURCE_ONLY_SCRIPT" \
        "other_lint.sh:$SOURCE_ONLY_SCRIPT")"

# The same omission, with the two names swapped so one is a suffix of the other.
# A block searched as raw text answers "is `lint.sh` in here?" with the
# `other_lint.sh` row, so this tree -- whose block names only `other_lint.sh`
# while `lint.sh` is source-only too -- reports PASS. The case above does not
# catch it: there the listed name is the shorter one, which no row contains.
check_tree "source-only script whose name is a suffix of a listed one" FLAG \
    "$(make_tree suffix_unlisted "$(readme_naming other_lint.sh)" \
        "lint.sh:$SOURCE_ONLY_SCRIPT" \
        "other_lint.sh:$SOURCE_ONLY_SCRIPT")"

check_tree "env-reading script listed in the block" FLAG \
    "$(make_tree miscategorised "$(readme_naming lint.sh parity.sh)" \
        "lint.sh:$SOURCE_ONLY_SCRIPT" \
        "parity.sh:$ENV_SCRIPT")"

# The bug that shipped twice while this lint was being written: a variable both
# assigned and defaulted from the environment in one expansion is an input.
check_tree "script whose only input is a defaulted env var" FLAG \
    "$(make_tree defaulted "$(readme_naming lint.sh build.sh)" \
        "lint.sh:$SOURCE_ONLY_SCRIPT" \
        "build.sh:$DEFAULTED_ENV_SCRIPT")"

echo "== healthy trees the lint must accept =="

check_tree "block names exactly the source-only scripts" PASS \
    "$(make_tree healthy "$(readme_naming lint.sh other_lint.sh)

Also here: \`parity.sh\`, which needs a binary." \
        "lint.sh:$SOURCE_ONLY_SCRIPT" \
        "other_lint.sh:$SOURCE_ONLY_SCRIPT" \
        "parity.sh:$ENV_SCRIPT")"

# The loud half of the same collision. `lint.sh` reads the environment and is
# correctly absent from the block, but a raw-text search finds it in the
# `other_lint.sh` row and reports it as wrongly listed. That direction fails a
# healthy tree rather than passing a stale one, so it is asserted as PASS here.
check_tree "env-reading script whose name is a suffix of a listed one" PASS \
    "$(make_tree suffix_env "$(readme_naming other_lint.sh)

Also here: \`lint.sh\`, which needs a binary." \
        "lint.sh:$ENV_SCRIPT" \
        "other_lint.sh:$SOURCE_ONLY_SCRIPT")"

# `test/regression/*.sh` is a glob, not a reference to a script named `*.sh`.
check_tree "glob in prose is not a dangling reference" PASS \
    "$(make_tree glob "$(readme_naming lint.sh)

\`ls test/regression/*.sh\` is the authoritative list." \
        "lint.sh:$SOURCE_ONLY_SCRIPT")"

# An env var named only in a header comment is documentation, not an input.
check_tree "env var named only in a comment is not an input" PASS \
    "$(make_tree commented "$(readme_naming lint.sh)" \
        "lint.sh:# Inputs: none. Contrast \${BWA_MEM3:-} in the parity scripts.
$SOURCE_ONLY_SCRIPT")"

# Shell text held as data. This is the shape of the self-tests next to this
# one, which build fixture scripts from single-quoted constants that run over
# several lines: `${MAKE:-make}` inside them is a string, not an expansion, and
# only tracking the quote state across lines can tell the difference. Dropping
# that tracking leaves every other case here passing, so without this one the
# self-test does not cover it.
check_tree "env var inside a multi-line single-quoted string is not an input" PASS \
    "$(make_tree quoted_data "$(readme_naming lint.sh)" \
        "lint.sh:$QUOTED_DATA_SCRIPT")"

# An apostrophe inside a double-quoted message is a literal character, not the
# start of a quoted span. Reading it as one shifts the quote state for
# everything after, which turns the single-quoted payload on the next line into
# live code and invents an input the script does not take.
check_tree "apostrophe in a double-quoted string does not shift quote state" PASS \
    "$(make_tree apostrophe "$(readme_naming lint.sh)" \
        "lint.sh:set -euo pipefail
label=\"the lint's own message\"
payload='\${BWA_MEM3:?BWA_MEM3 must be set}'
echo \"\$label \$payload\" > /dev/null
echo \"PASS: fixture\"")"

echo "== a lint that checked nothing must not report PASS =="

no_readme_dir="$(make_tree no-readme "$(readme_naming lint.sh)" \
    "lint.sh:$SOURCE_ONLY_SCRIPT")"
rm -f "$no_readme_dir/test/regression/README.md"
check_tree "no README" FLAG "$no_readme_dir"

check_tree "no scripts at all" FLAG \
    "$(make_tree no-scripts "$(readme_naming lint.sh)")"

check_tree "no source-only script in the tree" FLAG \
    "$(make_tree none-source-only "$(readme_naming parity.sh)" \
        "parity.sh:$ENV_SCRIPT")"

check_tree "README with no block markers" FLAG \
    "$(make_tree no-markers "Every script reads its inputs from the environment." \
        "lint.sh:$SOURCE_ONLY_SCRIPT")"

check_tree "README with an empty block" FLAG \
    "$(make_tree empty-block "$(printf '%s\n%s\n' "$BEGIN_MARKER" "$END_MARKER")" \
        "lint.sh:$SOURCE_ONLY_SCRIPT")"

# Reported as an ordering problem rather than as an empty block: the end marker
# above its begin marker also leaves nothing between the two, so the emptiness
# test alone would name the wrong edit.
check_tree "README with its block markers the wrong way round" FLAG \
    "$(make_tree swapped-markers "$(printf '%s\n%s\n%s\n' \
        "Prose naming lint.sh." "$END_MARKER" "$BEGIN_MARKER")" \
        "lint.sh:$SOURCE_ONLY_SCRIPT")"

check_tree "README naming no script at all" FLAG \
    "$(make_tree no-names "$(printf '%s\n%s\nnothing here\n%s\n' \
        "Prose with no script names." "$BEGIN_MARKER" "$END_MARKER")" \
        "lint.sh:$SOURCE_ONLY_SCRIPT")"

# Checked on the message, not just the status: an unguarded `cd` into a missing
# directory also dies with exit 1 under `set -e`, so a status-only assertion
# passes either way. The point of the case is that the failure arrives in this
# family's own format rather than as bash's raw `cd:` diagnostic.
missing_dir="$fixture_root/does-not-exist"
missing_out="$(bash "$lint" "$missing_dir" 2>&1)" && missing_status=0 || missing_status=$?
if ((missing_status == 1)) && [[ $missing_out == FAIL:* ]]; then
    echo "  ok   missing directory -> FLAG"
else
    echo "  FAIL missing directory -> expected exit 1 and a FAIL: message," >&2
    echo "       got exit $missing_status: $missing_out" >&2
    failures=1
fi

bash "$lint" one two > /dev/null 2>&1 && usage_status=0 || usage_status=$?
if ((usage_status == 2)); then
    echo "  ok   too many arguments -> usage error (exit 2)"
else
    echo "  FAIL too many arguments -> exit $usage_status, wanted 2" >&2
    failures=1
fi

if ((failures != 0)); then
    echo "FAIL: readme_contract_lint_selftest" >&2
    exit 1
fi

echo "PASS: readme_contract_lint_selftest"
