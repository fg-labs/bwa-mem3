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

# The step name every fixture's workflow defines, and that readme_naming puts
# in the table's origin column. Check 3 compares the two, so the default tree
# has to satisfy it or every case below would fail for the wrong reason.
readonly DEFAULT_STEP='Run the regression scripts'
readonly FIXTURE_JOB_NAME='The build job'

# Extra step names for the fixture workflow, one per line, written into the
# workflow verbatim so a case can give a step the YAML quoting a real one needs.
# Set around the case that wants it; reset after, the way FIXTURE_EXEMPTIONS is
# used in the sibling coverage-lint self-test.
FIXTURE_STEPS=""

# Build a fixture tree: a README, a workflow, plus one script per "name:body"
# pair given after it. Returns the tree's path on stdout.
make_tree() {
    local name="$1" readme_body="$2"
    shift 2

    local dir="$fixture_root/$name"
    rm -rf "$dir"
    mkdir -p "$dir/test/regression" "$dir/.github/workflows"
    printf '%s\n' "$readme_body" > "$dir/test/regression/README.md"

    # Every tree carries a workflow, because check 3 treats a missing one as
    # "nothing was checked" rather than as a pass -- so a tree without it would
    # make every case here fail on the workflow instead of on what it tests.
    # The job carries a name of its own, so a case can aim the origin column at
    # it and confirm a job name does not satisfy a column that names steps.
    {
        printf 'jobs:\n  build:\n    name: %s\n    steps:\n' "$FIXTURE_JOB_NAME"
        printf '      - name: %s\n' "$DEFAULT_STEP"
        local step
        while IFS= read -r step; do
            [[ -n $step ]] && printf '      - name: %s\n' "$step"
        done <<< "$FIXTURE_STEPS"
    } > "$dir/.github/workflows/ci.yml"

    local spec
    for spec in "$@"; do
        printf '%s\n' "${spec#*:}" > "$dir/test/regression/${spec%%:*}"
    done

    printf '%s' "$dir"
}

# Run the lint over a fixture tree and check the verdict. `expected` is FLAG
# (the lint must reject it, exit 1) or PASS (must accept it, exit 0).
#
# The optional fourth argument is a phrase the report has to carry. A verdict
# alone cannot tell two rejection branches apart, so a case that exists to move
# a tree from one branch to another reads the same before and after the change
# it pins -- green while checking nothing, which is the failure this whole file
# is here to prevent. Such a case names its branch.
check_tree() {
    local description="$1" expected="$2" dir="$3" must_say="${4-}"

    local output status
    output="$(bash "$lint" "$dir" 2>&1)" && status=0 || status=$?

    local verdict
    case "$status" in
        0) verdict=PASS ;;
        1) verdict=FLAG ;;
        *) verdict="EXIT$status" ;;
    esac

    if [[ $verdict != "$expected" ]]; then
        echo "  FAIL $description -> got $verdict, wanted $expected" >&2
        printf '       %s\n' "$output" >&2
        failures=1
        return
    fi

    if [[ -n $must_say && $output != *"$must_say"* ]]; then
        echo "  FAIL $description -> $expected, but the report never says" >&2
        echo "       '$must_say'" >&2
        printf '       %s\n' "$output" >&2
        failures=1
        return
    fi

    echo "  ok   $description -> $expected"
}

# A README whose source-only block names exactly the scripts passed in, and
# whose table gives each one an origin step the fixture workflow defines.
#
# `ORIGIN_OVERRIDE` replaces the origin cell of every row, so a case can aim a
# stale or unquoted value at check 3 without rebuilding the README by hand.
ORIGIN_OVERRIDE=""
readme_naming() {
    local rows="" table="" name origin
    # Spelled out rather than as `${ORIGIN_OVERRIDE:-...}`: check 2 of the lint
    # under test reads a name used in a `:-` expansion as an environment input
    # by design, which would move this self-test out of the source-only set it
    # belongs in -- and the lint would then correctly fail on its own README.
    if [[ -n $ORIGIN_OVERRIDE ]]; then
        origin="$ORIGIN_OVERRIDE"
    else
        origin="\"$DEFAULT_STEP\""
    fi
    for name in "$@"; do
        rows+="- \`$name\` takes an optional positional directory"$'\n'
        table+="| \`$name\` | what it checks | $origin |"$'\n'
    done
    printf '%s\n\n| Script | What it checks | Origin in ci.yml |\n|---|---|---|\n%s\n%s\n%s%s\n' \
        "Every script but the source-only lints reads its inputs from the environment." \
        "$table" \
        "$BEGIN_MARKER" "$rows" "$END_MARKER"
}

# The fixture bodies below are shell source held as data, so every `$` in them
# must survive into the fixture file unexpanded. Single quotes are the point,
# not an oversight -- and the lint under test has to read them as data too,
# which is what several of the cases check.
# shellcheck disable=SC2016

# A script that takes the named environment variable, and one that takes none.
#
# Every fixture below carries both a PASS: and a FAIL: literal, because check 4
# requires them of any script in the tree. The cases here are aimed at checks 1
# and 2, so the markers are scaffolding rather than the thing under test -- a
# fixture missing one would FLAG for a reason the case never meant to exercise.
# The dedicated check-4 cases further down drop a marker on purpose.
readonly ENV_SCRIPT='set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"
[[ -n ${BWA_MEM3:-} ]] || { echo "FAIL: fixture" >&2; exit 1; }
echo "PASS: fixture"'
# shellcheck disable=SC2016
readonly SOURCE_ONLY_SCRIPT='set -euo pipefail
scan_root="${1:-src}"
[[ -n $scan_root ]] || { echo "FAIL: fixture" >&2; exit 1; }
echo "PASS: fixture $scan_root"'
# shellcheck disable=SC2016
readonly DEFAULTED_ENV_SCRIPT='set -euo pipefail
MAKE="${MAKE:-make}"
"$MAKE" --version || { echo "FAIL: fixture" >&2; exit 1; }
echo "PASS: fixture"'

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
printf '"'"'%s\n'"'"' "$payload" > /dev/null || { echo "FAIL: fixture" >&2; exit 1; }
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
echo \"\$label \$payload\" > /dev/null || { echo \"FAIL: fixture\" >&2; exit 1; }
echo \"PASS: fixture\"")"

echo "== the ci.yml origin column ==" # check 3

# The drift this check was added for: two rows named steps that no longer
# existed, because a PR merged the two steps they described into one and only
# the workflow was updated. A reader following the column looks for a step name
# that is not there.
ORIGIN_OVERRIDE='"A step no workflow defines"'
check_tree "origin column names a step ci.yml does not define" FLAG \
    "$(make_tree stale-origin "$(readme_naming lint.sh)" \
        "lint.sh:$SOURCE_ONLY_SCRIPT")"

# An unquoted cell is rejected rather than skipped: a row that opts out of the
# check silently is how the column rots again.
ORIGIN_OVERRIDE='Run the regression scripts'
check_tree "origin cell that is not quoted is rejected" FLAG \
    "$(make_tree unquoted-origin "$(readme_naming lint.sh)" \
        "lint.sh:$SOURCE_ONLY_SCRIPT")"

# An empty pair of quotes names no step, so it belongs with the unquoted cells
# rather than downstream, where the name strips to an empty grep pattern and is
# reported as a blank line under a heading about missing steps. Both branches
# exit 1, so the branch is named: on the verdict alone this case reads green
# against a lint that still routes the cell the wrong way.
ORIGIN_OVERRIDE='""'
check_tree "origin cell holding an empty pair of quotes is rejected" FLAG \
    "$(make_tree empty-origin "$(readme_naming lint.sh)" \
        "lint.sh:$SOURCE_ONLY_SCRIPT")" \
    "does not hold a"
ORIGIN_OVERRIDE=""

# The column names a *step*, not a job. A job name would let a row point at
# something that runs no script, so it must not satisfy the check.
ORIGIN_OVERRIDE="\"$FIXTURE_JOB_NAME\""
check_tree "a job name does not satisfy the origin column" FLAG \
    "$(make_tree job-name-origin "$(readme_naming lint.sh)" \
        "lint.sh:$SOURCE_ONLY_SCRIPT")"
ORIGIN_OVERRIDE=""

# The README has more than one three-column table, and the other one's third
# column holds a script name rather than a step name. Locating the table by
# position instead of by its header reads those rows as origins and reports
# every one of them as a step that does not exist.
check_tree "a second table is not read as origins" PASS \
    "$(make_tree second-table "$(readme_naming lint.sh)

| Lint | Positional argument | Self-test |
|------|---------------------|-----------|
| \`lint.sh\` | directory to scan | \`lint.sh\` |" \
        "lint.sh:$SOURCE_ONLY_SCRIPT")"

# A data cell may contain three hyphens. Recognising the separator row by
# "contains ---" rather than by its shape drops such a row from the check and
# from the row count -- a row opting out of check 3 in silence, which is what an
# unquoted cell is rejected for two cases above. Only the second row here is
# stale, so the case can only reach FLAG if that row was read.
emdash_readme="$(
    cat << EOF
Every script but the source-only lints reads its inputs from the environment.

| Script | What it checks | Origin in ci.yml |
|---|---|---|
| \`lint.sh\` | a plain cell | "$DEFAULT_STEP" |
| \`other.sh\` | a cell holding --- three hyphens | "A step no workflow defines" |

$BEGIN_MARKER
- \`lint.sh\` takes an optional positional directory
- \`other.sh\` takes an optional positional directory
$END_MARKER
EOF
)"
check_tree "a data cell holding --- is still read as a row" FLAG \
    "$(make_tree emdash-cell "$emdash_readme" \
        "lint.sh:$SOURCE_ONLY_SCRIPT" "other.sh:$SOURCE_ONLY_SCRIPT")"

# YAML's quoting is not part of the step name. A name containing `: ` cannot be
# written unquoted at all, so a check that compares the quoted spelling reports
# a step the workflow plainly defines as missing -- and tells the reader to fix
# a column that is already correct.
FIXTURE_STEPS='"Regression: option validation"'
ORIGIN_OVERRIDE='"Regression: option validation"'
check_tree "a double-quoted step name satisfies the column" PASS \
    "$(make_tree double-quoted-step "$(readme_naming lint.sh)" \
        "lint.sh:$SOURCE_ONLY_SCRIPT")"

FIXTURE_STEPS="'Regression: single quoted'"
ORIGIN_OVERRIDE='"Regression: single quoted"'
check_tree "a single-quoted step name satisfies the column" PASS \
    "$(make_tree single-quoted-step "$(readme_naming lint.sh)" \
        "lint.sh:$SOURCE_ONLY_SCRIPT")"
FIXTURE_STEPS=""
ORIGIN_OVERRIDE=""

# The header string is not unique in the file: the real README quotes it inside
# the row that describes this very check. Locating the table with a plain match
# takes the first prose mention instead and reports the table under it as empty.
check_tree "the column header named in prose above the table" PASS \
    "$(make_tree prose-header "$(printf '%s\n\n%s\n' \
        "The Origin in ci.yml column names the step that runs each script." \
        "$(readme_naming lint.sh)")" \
        "lint.sh:$SOURCE_ONLY_SCRIPT")"

echo "== the PASS:/FAIL: marker contract ==" # check 4

# Every case below names this branch rather than resting on the verdict. Each
# one hands the lint a tree that is stale in exactly one way, and checks 1-3
# reject a tree for three other reasons with the same exit 1 -- so a case that
# tripped one of those by accident reads green while never reaching check 4.
readonly MARKER_REPORT='cannot emit the markers'

# The three shapes that were actually in the tree when check 4 was written: no
# PASS: at all (meth_collapsed_scoring.sh ended on a bare `echo "OK"`), a PASS
# without the colon (all_tiers_parity.sh printed "ALL TIERS PARITY: PASS ("),
# and neither marker (meth_oracle.sh `exec`d its child and said nothing).
# shellcheck disable=SC2016
readonly NO_PASS_SCRIPT='set -euo pipefail
scan_root="${1:-src}"
[[ -n $scan_root ]] || { echo "FAIL: fixture" >&2; exit 1; }
echo "OK"'
# shellcheck disable=SC2016
readonly UNCOLONED_PASS_SCRIPT='set -euo pipefail
scan_root="${1:-src}"
[[ -n $scan_root ]] || { echo "FAIL: fixture" >&2; exit 1; }
echo "FIXTURE PARITY: PASS (1 comparison)"'
readonly NO_MARKER_SCRIPT='set -euo pipefail
exec bash some/other/harness'
# The asymmetric case. Every shape above is missing PASS: or missing both, so a
# check that dropped the FAIL: half entirely would still turn all of them FLAG
# and this file would report green while covering half the contract.
# shellcheck disable=SC2016
readonly NO_FAIL_SCRIPT='set -euo pipefail
scan_root="${1:-src}"
[[ -n $scan_root ]] || exit 1
echo "PASS: fixture"'
# Both markers present, and neither reachable: a script that documents the
# contract in its header and prints nothing. That is what a raw substring match
# accepts, and what the header of the lint itself would have satisfied.
readonly COMMENT_ONLY_SCRIPT='set -euo pipefail
# Prints PASS: on success and FAIL: on failure.
exec bash some/other/harness'
# The same shape one column over: the markers trail live code instead of
# starting the line, so a `^[[:space:]]*#` filter leaves them in place while
# the shell reads them as comment text and the script still prints nothing.
# Distinguishing this from the quoted `"... # PASS:"` on the line below it is
# what the lint needs quote state for, and the `PASS` case pins that both ways:
# the marker in the comment must not count, and the one in the string must.
readonly TRAILING_COMMENT_SCRIPT='set -euo pipefail
echo "ran # PASS: fixture"
exec bash some/other/harness # reports FAIL: on a bad exit'
# The accepted half of the same rule, in the other quote. Dropping comments has
# to leave both kinds of quoted span intact, and single quotes are the ones a
# strip written for `$NAME` expansion would throw away wholesale -- check 2
# drops them by design, and the two checks share a walk. shell_lint_selftest.sh
# is this shape today, so getting it wrong would fail a healthy tree.
# shellcheck disable=SC2016
readonly SINGLE_QUOTED_MARKER_SCRIPT='set -euo pipefail
scan_root="${1:-src}"
[[ -n $scan_root ]] || { printf '"'"'FAIL: fixture\n'"'"' >&2; exit 1; }
printf '"'"'PASS: fixture\n'"'"''

check_tree "script with no PASS: marker" FLAG \
    "$(make_tree no-pass-marker "$(readme_naming lint.sh)" \
        "lint.sh:$NO_PASS_SCRIPT")" \
    "$MARKER_REPORT"

check_tree "script with no FAIL: marker" FLAG \
    "$(make_tree no-fail-marker "$(readme_naming lint.sh)" \
        "lint.sh:$NO_FAIL_SCRIPT")" \
    "$MARKER_REPORT"

check_tree "script whose PASS lacks the colon" FLAG \
    "$(make_tree uncoloned-pass "$(readme_naming lint.sh)" \
        "lint.sh:$UNCOLONED_PASS_SCRIPT")" \
    "$MARKER_REPORT"

# `other.sh` carries the bad shape while `lint.sh` stays healthy, and both are
# named in the block because neither reads the environment. Putting the shape on
# a second script is what keeps check 2 from flagging the tree first, which
# would let the case pass without ever reaching check 4.
check_tree "markers only in a comment" FLAG \
    "$(make_tree comment-only-markers "$(readme_naming lint.sh other.sh)" \
        "lint.sh:$SOURCE_ONLY_SCRIPT" "other.sh:$COMMENT_ONLY_SCRIPT")" \
    "$MARKER_REPORT"

# Split across two scripts for the same reason. The report is checked down to
# the marker this tree is missing, because the script's other marker is inside
# a string on the line above: a lint that dropped the whole line -- comment,
# quoted text and all -- would flag the tree for both markers and read green
# here while rejecting every script that prints a marker mid-line.
check_tree "marker only in a trailing comment" FLAG \
    "$(make_tree trailing-comment-markers "$(readme_naming lint.sh other.sh)" \
        "lint.sh:$SOURCE_ONLY_SCRIPT" "other.sh:$TRAILING_COMMENT_SCRIPT")" \
    "other.sh (FAIL:)"

# Both markers missing at once -- the wrapper shape. Split across two scripts
# for the same reason as the case above.
check_tree "wrapper emitting neither marker" FLAG \
    "$(make_tree no-markers "$(readme_naming lint.sh other.sh)" \
        "lint.sh:$SOURCE_ONLY_SCRIPT" "other.sh:$NO_MARKER_SCRIPT")" \
    "$MARKER_REPORT"

check_tree "markers only inside single quotes" PASS \
    "$(make_tree single-quoted-markers "$(readme_naming lint.sh)" \
        "lint.sh:$SINGLE_QUOTED_MARKER_SCRIPT")"

echo "== a lint that checked nothing must not report PASS =="

no_workflows_dir="$(make_tree no-workflows "$(readme_naming lint.sh)" \
    "lint.sh:$SOURCE_ONLY_SCRIPT")"
rm -rf "$no_workflows_dir/.github"
check_tree "no workflows directory" FLAG "$no_workflows_dir"

empty_workflow_dir="$(make_tree empty-workflow "$(readme_naming lint.sh)" \
    "lint.sh:$SOURCE_ONLY_SCRIPT")"
printf 'jobs: {}\n' > "$empty_workflow_dir/.github/workflows/ci.yml"
check_tree "workflow defining no steps" FLAG "$empty_workflow_dir"

# A README with no origin table at all leaves check 3 with nothing to compare,
# which is the same green-while-doing-no-work state the rest of this block
# guards against.
#
# Named rather than left to the status, because the header search is a grep
# pipeline: with no table to match, an unguarded grep exits 1 and `pipefail`
# carries that out of the command substitution, killing the script on the
# assignment with exit 1 -- the same status this branch reports, and no message
# at all. Status alone reads FLAG either way.
check_tree "README with no origin table" FLAG \
    "$(make_tree no-origin-table "$(printf '%s\n\n%s\n- %s takes a directory\n%s\n' \
        "Prose naming lint.sh but carrying no table." "$BEGIN_MARKER" "lint.sh" "$END_MARKER")" \
        "lint.sh:$SOURCE_ONLY_SCRIPT")" \
    "no table column headed 'Origin in ci.yml'"

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
