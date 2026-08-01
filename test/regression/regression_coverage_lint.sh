#!/bin/bash
# Every script in test/regression/ must be wired into CI, not just into
# `make test`. What that resolves to in practice -- a workflow naming the
# script, or invoking a Makefile target whose recipe runs it -- is spelled out
# below; this lint checks the wiring, not that the run happened.
#
# The regression scripts are enumerated in two independent, hand-maintained
# places: the `test:` target in the Makefile, and explicit steps and loops in
# .github/workflows/. A script listed in one and not the other still looks
# wired up -- it has a home, it passes when run by hand -- while in fact no
# automated run ever executes it. Nothing fails; the coverage just is not
# there.
#
# CI does not invoke `make test`. It builds with `make test-binaries` and then
# runs the regression scripts from its own list, so the `test:` recipe is a
# developer convenience rather than the thing CI executes. That makes the
# Makefile list and the workflow list genuinely independent, and drift between
# them invisible. Both scripts this lint was written for were found that way:
# one reachable only from `test:`, one only from `test-injection`, a target no
# workflow names.
#
# A script counts as covered when a workflow names it, or when it is in the
# recipe of a Makefile target a workflow actually invokes. Comments do not
# count, on either side: a step's comment usually names the script the step
# runs, so counting it would let the coverage outlive the invocation that
# produced it. Anything else has to be listed in
# test/regression/coverage_exemptions.txt with the reason -- so a script nothing
# runs is a reviewable line there rather than silence.
#
# Sibling lints: ndebug_gate_lint.sh, debug_macro_flag_lint.sh. Same failure
# in each case -- a check that reads as green while doing no work.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

if (($# > 1)); then
    echo "usage: ${BASH_SOURCE[0]##*/} [repo-root-to-check]" >&2
    exit 2
fi

# Optional argument so regression_coverage_lint_selftest.sh can aim the same
# resolution at a fixture tree; a lint that has stopped resolving anything
# reports the same PASS as one whose lists agree.
#
# Checked rather than left to `cd` to fail under `set -e`: bash's own `cd:` line
# carries no FAIL: prefix, so the one way of pointing this lint at nothing that
# is reachable from the command line would be the one failure a caller
# scanning for the family's format does not see.
if (($# == 1)); then
    if [[ ! -d $1 ]]; then
        echo "FAIL: '$1' is not a directory -- nothing was checked" >&2
        exit 1
    fi
    cd "$1"
fi

workflow_dir=".github/workflows"
script_dir="test/regression"
makefile="Makefile"

for required in "$workflow_dir" "$script_dir" "$makefile"; do
    if [[ ! -e $required ]]; then
        echo "FAIL: $required missing under $PWD -- nothing was checked" >&2
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# Scripts deliberately not run by CI.
#
# Kept in a data file rather than an array in here, so adding an exemption is a
# reviewable one-line diff next to its reason, and so this list can be exercised
# by the self-test against fixture trees that carry their own.
# ---------------------------------------------------------------------------
exemptions_file="$script_dir/coverage_exemptions.txt"
exempt_scripts=()
if [[ -f $exemptions_file ]]; then
    while IFS= read -r entry; do
        entry="${entry%%#*}"
        entry="${entry//[[:space:]]/}"
        if [[ -n $entry ]]; then
            exempt_scripts+=("$entry")
        fi
    done < "$exemptions_file"
fi

# Read loops rather than `mapfile`, which is bash 4+: macOS still ships bash
# 3.2 as /bin/bash, and the sibling lints run there. A lint that cannot run on
# a contributor's shell is the same missing coverage by another route.
workflow_files=()
while IFS= read -r file; do
    [[ -n $file ]] && workflow_files+=("$file")
done < <(find "$workflow_dir" -type f \( -name '*.yml' -o -name '*.yaml' \) | sort)
if ((${#workflow_files[@]} == 0)); then
    echo "FAIL: no workflow files under $workflow_dir/ -- nothing was checked" >&2
    exit 1
fi

scripts=()
while IFS= read -r file; do
    [[ -n $file ]] && scripts+=("$file")
done < <(find "$script_dir" -type f -name '*.sh' | sort)
if ((${#scripts[@]} == 0)); then
    echo "FAIL: no scripts under $script_dir/ -- nothing was checked" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Comments are stripped before anything is matched against this text.
#
# A script counts as covered when a workflow names it. Without this, a workflow
# *comment* naming a script satisfies that -- and comments here routinely name
# the very script the step below runs, because the comment explains why it runs.
# So the coverage survives deleting the thing that produces it: remove the
# invocation, leave the comment that describes it, and the lint stays green
# while nothing runs the script. That is the same read-as-green-while-doing-no-
# work failure the lint exists to catch, one level up.
#
# `#` opens a comment only where a word starts and only outside quotes, which
# is the rule in both YAML and the shell inside a `run:` block. Both halves
# matter, and they pull in opposite directions:
#
#   - A word starts at the start of a line, after whitespace, and after a
#     control operator -- `;`, `&`, `|`, `(`, `)` -- which ends the command
#     before it. Treating only whitespace as the boundary reads `echo hi;# bash
#     foo.sh` as a token, so the commented-out invocation counts as coverage:
#     this lint's own bug, one level down.
#   - Inside '...' or "..." a `#` is a literal, whatever precedes it. Truncating
#     at the `#` in `echo "counted # of reads" && bash foo.sh` drops the real
#     invocation after it and calls a script CI does run uncovered.
#
# A `#` inside a token (`sha256#...`, a URL fragment) is left alone by the same
# word-start rule, and a backslash escapes the character after it.
#
# A quote opens a quoted run only when it does not follow a word character and
# its partner is on the same line, and state resets at each newline rather than
# carrying across the file. All three rules exist for the same reason: an
# apostrophe in prose -- `- name: Don't rebuild` -- is far more common here than
# a string spanning lines, and treating one as an opening quote would make every
# `#` after it on that line read as a literal. That is the false-coverage
# direction, which is the one that fails silently.
#
# Pairing alone does not get there, because the partner it finds may be another
# prose apostrophe on the far side of the comment: in `- name: Don't rebuild  #
# runs alpha.sh when it's stale` the two apostrophes pair across the `#`, the
# comment reads as quoted text, and the script it names counts as covered. So a
# quote must also start a word. An apostrophe inside one is prose (`Don't`,
# `it's`); a real opening quote follows whitespace or punctuation (`bash -c
# '...'`, `--flag='a # b'`, `$'...'`), never a letter or digit.
#
# The Makefile recipes gathered below are stripped the same way, and so is the
# text the target-invocation search reads: a comment mentioning
# `make test-injection` must not count as invoking that target either.
# ---------------------------------------------------------------------------
strip_comments() {
    awk '
    {
        n = length($0)
        in_single = 0
        in_double = 0
        kept = ""
        for (i = 1; i <= n; i++) {
            c = substr($0, i, 1)
            if (in_single) {
                if (c == "\047") in_single = 0
                kept = kept c
                continue
            }
            if (c == "\\" && i < n) {
                kept = kept c substr($0, i + 1, 1)
                i++
                continue
            }
            if (in_double) {
                if (c == "\"") in_double = 0
                kept = kept c
                continue
            }
            if (c == "\047" || c == "\"") {
                prev = (i == 1) ? "" : substr($0, i - 1, 1)
                if (prev !~ /[A-Za-z0-9]/ && index(substr($0, i + 1), c)) {
                    if (c == "\047") in_single = 1; else in_double = 1
                }
                kept = kept c
                continue
            }
            if (c == "#") {
                prev = (i == 1) ? "" : substr($0, i - 1, 1)
                if (prev == "" || prev ~ /[ \t;&|()]/) break
            }
            kept = kept c
        }
        print kept
    }'
}

workflow_text="$(cat "${workflow_files[@]}" | strip_comments)"

# ---------------------------------------------------------------------------
# Recipes of Makefile targets that a workflow actually invokes.
#
# Matched as `make [flags/vars...] <target>` so a target is only counted when a
# workflow really calls it. Substring matching would be worse than useless
# here: `test` would match every `make test-binaries` line, pull the `test:`
# recipe in, and mark as covered exactly the scripts this lint exists to find.
#
# Both halves therefore treat the target name as a literal. In the ERE it is
# escaped; the recipe scan compares whole names instead of building a regex out
# of one. Interpolated raw, this Makefile's own `.PHONY`, `.SUFFIXES`, `.c.o`
# and `.cpp.o` each turn their `.` into a wildcard, and a target that matches
# the wrong line attaches the wrong recipe -- over-counting coverage in the
# same silent direction as the substring trap above.
# ---------------------------------------------------------------------------
invoked_recipes=""
while IFS= read -r target; do
    [[ -z $target ]] && continue

    # Target names are drawn from [A-Za-z0-9_.-] by the enumeration below, and
    # `.` is the only ERE metacharacter in that set.
    target_re="${target//./\\.}"
    # `make` must be the whole tool name. Unanchored it also matches the tail of
    # `cmake`, `gmake` and `nmake`, so `cmake --build <target>` -- a line that
    # never runs this Makefile -- would resolve that target and count its recipe
    # as covered. Spelled out rather than as `\b`/`[[:<:]]`, neither of which is
    # portable across GNU and BSD grep. `/` is deliberately not excluded:
    # `/usr/bin/make check` is a real invocation.
    if ! grep -qE \
        "(^|[^A-Za-z0-9_])make([[:space:]]+[^[:space:]]+)*[[:space:]]+${target_re}([[:space:]]|\$)" \
        <<< "$workflow_text"; then
        continue
    fi
    recipe="$(awk -v t="$target" '
        # A rule line opens a recipe; a later unindented, non-comment line
        # closes it. One rule line may name several targets ("fmt lint: deps")
        # and the recipe belongs to each, so the whole list left of ":" is
        # searched -- anchoring on a single leading name misses every target
        # but the first, and reports its scripts as run by nothing.
        /^[^\t#]/ && index($0, ":") > 0 {
            in_recipe = 0
            n = split(substr($0, 1, index($0, ":") - 1), names, /[ \t]+/)
            for (i = 1; i <= n; i++) {
                if (names[i] == t) in_recipe = 1
            }
            next
        }
        in_recipe && /^\t/ { print; next }
        in_recipe && !/^\t/ && !/^#/ && NF { in_recipe = 0 }
    ' "$makefile" | strip_comments)"
    invoked_recipes+="$recipe"$'\n'
done < <(awk '
    # Rule lines only. Plenty of other lines carry a ":": assignments
    # (`.DEFAULT_GOAL := ...`), and directives whose arguments hide one inside a
    # substitution reference (`-include $(OBJS:.o=.d)`). Both would otherwise
    # enumerate as targets, with whatever follows read as their recipe.
    #
    # So the whole list left of ":" must consist of nothing but target names --
    # one token that is not a name disqualifies the line, rather than being
    # dropped while its neighbours are kept.
    /^[a-zA-Z0-9_.-]/ && index($0, ":") > 0 {
        head = substr($0, 1, index($0, ":") - 1)
        if (substr($0, index($0, ":") + 1, 1) == "=") next
        if (substr($0, index($0, ":") + 1, 2) == ":=") next
        n = split(head, names, /[ \t]+/)
        for (i = 1; i <= n; i++) {
            if (names[i] !~ /^[a-zA-Z0-9_.-]+$/) next
        }
        for (i = 1; i <= n; i++) print names[i]
    }
' "$makefile" | sort -u)

covered_text="$workflow_text"$'\n'"$invoked_recipes"

# ---------------------------------------------------------------------------
# Coverage is matched on whole reference tokens, not substrings.
#
# Every lint here ships beside an `_selftest` sibling whose name contains it, so
# a substring match reads `regression_coverage_lint` as covered on the strength
# of `regression_coverage_lint_selftest.sh` alone -- and dropping the shorter
# script's own reference would then leave it uncovered and still green. Same
# class of bug as the `make test-binaries` trap above, one level down.
#
# Splitting the text on everything that cannot appear inside a script reference
# leaves whole tokens (`alpha.sh`, `alpha_selftest.sh`) to compare exactly,
# which needs no regex escaping of the names. `.` and `-` stay inside tokens so
# a name carrying either is still matched as a unit rather than in pieces --
# `alpha-beta` is as much its own script as `alpha_selftest` is.
#
# Keeping `.` and `-` means a reference ending a sentence would otherwise
# tokenize as `alpha.sh.` and match nothing, so trailing runs of them are
# trimmed. The Makefile already writes one reference that way.
# ---------------------------------------------------------------------------
covered_tokens="$(LC_ALL=C tr -c 'A-Za-z0-9_.-' '\n' <<< "$covered_text" \
    | sed 's/[.-]*$//' | sort -u)"

failures=0
exempt_used=()

for path in "${scripts[@]}"; do
    name="$(basename "$path" .sh)"

    # `${a[@]+"${a[@]}"}` rather than `"${a[@]}"`: under `set -u`, bash 3.2
    # treats an empty array's expansion as an unbound variable, and an empty
    # exemption list is the normal case.
    is_exempt=0
    for exempt in ${exempt_scripts[@]+"${exempt_scripts[@]}"}; do
        if [[ $name == "$exempt" ]]; then
            is_exempt=1
            break
        fi
    done

    # Herestring, not `printf ... | grep`: with `set -o pipefail`, `grep -q`
    # exits on its first match, printf dies of SIGPIPE, and the pipeline
    # reports 141 -- turning a found script into "not covered" for whichever
    # names happen to appear early enough in the text to short-circuit grep.
    #
    # Both spellings count: references carry the extension (`alpha.sh`), while
    # a bare name is how a script gets mentioned without one.
    if grep -qxF -- "$name" <<< "$covered_tokens" \
        || grep -qxF -- "$name.sh" <<< "$covered_tokens"; then
        if ((is_exempt)); then
            echo "FAIL: $path is listed as exempt but IS run by CI -- drop the" >&2
            echo "      exemption so the list keeps meaning something." >&2
            failures=$((failures + 1))
        fi
        continue
    fi

    if ((is_exempt)); then
        exempt_used+=("$name")
        echo "note: $name is not run by CI (listed as exempt)"
        continue
    fi

    echo "FAIL: $path is run by no workflow, and by no Makefile target that a" >&2
    echo "      workflow invokes. It passes only when someone runs it by hand." >&2
    echo "      Add it to a workflow, or list it in $exemptions_file with a" >&2
    echo "      reason." >&2
    failures=$((failures + 1))
done

# The loop above only reads the exemption list from the script side, so it
# catches an entry that has become untrue but never one that has stopped naming
# anything. An exemption whose script was renamed or deleted can no longer be
# contradicted by any check, and so outlives the change that made it
# meaningless -- silence again, in the one file whose whole point is that
# leaving a script uncovered costs a reviewable line.
for exempt in ${exempt_scripts[@]+"${exempt_scripts[@]}"}; do
    found=0
    for path in "${scripts[@]}"; do
        if [[ "$(basename "$path" .sh)" == "$exempt" ]]; then
            found=1
            break
        fi
    done
    if ((!found)); then
        echo "FAIL: $exemptions_file exempts '$exempt', which is not a script" >&2
        echo "      under $script_dir/. Drop the entry, or fix the name." >&2
        failures=$((failures + 1))
    fi
done

if ((failures > 0)); then
    echo "FAIL: regression_coverage_lint ($failures problem(s))" >&2
    exit 1
fi

echo "PASS: regression_coverage_lint (${#scripts[@]} scripts, ${#exempt_used[@]} exempt)"
