#!/bin/bash
# Keep the opt-in debug/instrumentation macro list in ci.yml and the macros in
# src/ pointing at each other.
#
# ci.yml has one step that rebuilds with every opt-in debug macro enabled and
# reruns the chr22 parity check, so the guarded code is compiled and exercised
# rather than left to rot (see the "Opt-in debug/instrumentation macro build"
# step). That step is only worth anything if its `-D` flags name macros that
# actually exist. Two ways it silently stops being worth anything:
#
#   - A flag names no macro. The build compiles it into nothing, the parity
#     check passes, and the check the flag was supposed to enable never runs.
#     This is not hypothetical: the flags for three checks were added ahead of
#     the PRs introducing them, so for a while the step really was enabling
#     nothing for those three. Two have since landed and matched; one has not.
#     A rename of one character during review would have been invisible.
#
#   - A macro exists but no flag names it. A new BWA_MEM3_DEBUG_* check added
#     to src/ without a corresponding -D is never compiled by any CI row --
#     the exact gap the step was added to close, reopened one macro at a time.
#
# So this lint checks both directions. It is the same failure this repository
# keeps rediscovering: a check that reads as green while doing no work. See
# test/regression/ndebug_gate_lint.sh for the sibling case.
#
# Scope note: the reverse direction is restricted to the BWA_MEM3_DEBUG_ prefix
# because that is the one unambiguous naming convention for this family. Build
# mode switches (STAGE_PROF, DISABLE_OUTPUT, COVERAGE) have their own CI rows
# and are deliberately not in scope here.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

if (( $# > 1 )); then
    echo "usage: ${BASH_SOURCE[0]##*/} [repo-root-to-check]" >&2
    exit 2
fi

# Optional argument so debug_macro_flag_lint_selftest.sh can aim the same
# checks at a fixture tree; there is no other way to tell "the lists agree"
# apart from "the extraction stopped extracting".
if (( $# == 1 )); then
    cd "$1"
fi

workflow=".github/workflows/ci.yml"
src_dir="src"

for required in "$workflow" "$src_dir"; do
    if [[ ! -e $required ]]; then
        echo "FAIL: $required missing under $PWD -- nothing was checked" >&2
        exit 1
    fi
done

# Existing is not the same as usable, and the type errors do not announce
# themselves: `grep` over a directory named ci.yml reads nothing, so the flag
# list looks absent and the failure below blames a restructured step; `find`
# over a regular file named src finds no sources, so the failure blames an empty
# preprocessor set. Both exit nonzero either way, but both name the wrong bug.
if [[ ! -f $workflow ]]; then
    echo "FAIL: $workflow is not a regular file -- nothing was checked" >&2
    exit 1
fi
if [[ ! -d $src_dir ]]; then
    echo "FAIL: $src_dir is not a directory -- nothing was checked" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Deliberately-ahead-of-the-code flags.
#
# A flag may name a macro that does not exist yet when it is added ahead of the
# PR introducing the check, which avoids a follow-up CI edit. That is fine, but
# it has to be stated here with a reason, so the exemption is a line a reviewer
# sees rather than silence. Delete the entry when the check lands.
# ---------------------------------------------------------------------------
pending_flags=(
    # Enabled ahead of the ungapped-CIGAR fast-path change that introduces it.
    BWA_MEM3_DEBUG_UNGAPPED_XCHECK
)

# ---------------------------------------------------------------------------
# Fold `\`-continued lines into one logical line.
#
# Every match below is line-oriented, and both of the constructs they match --
# a preprocessor conditional and a shell assignment -- can be split across
# physical lines. src/ already splits conditions that way today. A grep that
# reads only the first physical line silently drops what the rest carries,
# which is the same set-shrinking failure the `\b` note further down guards
# against: the dropped macro is never required to have a -D, so the reverse
# check stops covering it and says nothing.
#
# A file whose last line ends in `\` would join across a file boundary below,
# but that is undefined behaviour in C anyway, and the only consequence here is
# one spurious identifier in the set -- which fails loudly rather than quietly.
# ---------------------------------------------------------------------------
join_continuations() {
    sed -e ':a' -e '/\\$/{$!{N;s/\\\n/ /;ba' -e '}' -e '}'
}

# ---------------------------------------------------------------------------
# The -D flags the opt-in macro build passes.
#
# Anchored on the single EXTRA_CXXFLAGS line carrying the family prefix rather
# than on a step name or line number, either of which drifts silently. Zero or
# several matches means the step was restructured and this extraction no longer
# describes it -- a hard error, not a quiet empty list.
#
# Comment lines are dropped first, then both terms are required. Requiring both
# terms alone is not enough: a comment that quotes the assignment carries both,
# and would count as a second "flag list" -- a spurious restructured-step failure
# pointing at the wrong cause. Dropping comments is right under either reading of
# such a line, since `#` opens a comment in YAML and in the shell of a `run: |`
# block alike. A trailing `# ...` on a real flag line is untouched; only a line
# whose first character is `#` is one.
#
# Comments are dropped before joining, not after. A `#` comment runs to the end
# of its physical line in both languages -- a trailing `\` inside one continues
# nothing -- so joining first would let a comment absorb the real flag list and
# then discard both, reporting the list as absent.
#
# Numbered before either step: the number a wrapped list reports has to be the
# line the reader can go and look at, and numbering the joined stream would count
# logical lines the file does not have. The comment filter therefore has to skip
# that `N:` prefix to find the line's first real character.
# ---------------------------------------------------------------------------
flag_lines="$(grep -n '' "$workflow" \
    | grep -vE '^[0-9]+:[[:space:]]*#' \
    | join_continuations \
    | grep 'EXTRA_CXXFLAGS' | grep 'DBWA_MEM3_DEBUG' || true)"
flag_line_count="$(printf '%s' "$flag_lines" | grep -c . || true)"

if [[ $flag_line_count -ne 1 ]]; then
    echo "FAIL: expected exactly one -DBWA_MEM3_DEBUG flag list in $workflow," >&2
    echo "found $flag_line_count. The opt-in macro build step was restructured;" >&2
    echo "update the extraction in ${BASH_SOURCE[0]##*/} to match." >&2
    [[ -n $flag_lines ]] && printf '%s\n' "$flag_lines" >&2
    exit 1
fi

ci_flags="$(printf '%s' "$flag_lines" \
    | grep -oE '\-D[A-Za-z_][A-Za-z0-9_]*' \
    | sed 's/^-D//' | sort -u)"

if [[ -z $ci_flags ]]; then
    echo "FAIL: no -D flags parsed out of $workflow -- nothing was checked" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Macro names that appear in a real preprocessor conditional under src/.
#
# Identifiers are pulled with an explicit character class, never `\b` -- that
# is a GNU extension an implementation may read as a literal `b`, which would
# empty this set and make every check below vacuously pass.
#
# C comments are stripped off the directive line first, the same way YAML
# comments are stripped off the flag line above. A trailing `// see also
# BWA_MEM3_DEBUG_X` is routine, and without this the name lands in the set as
# though a conditional tested it -- so the reverse check demands a -D for a macro
# that gates nothing, reporting that src/ gates code on it when src/ does not.
# The obvious way to silence that is to add the flag, which then enables a macro
# no conditional tests: the forward-direction hole this lint exists to close,
# reopened by its own false alarm.
#
# Stripped after the directive-line match, not before, so a commented-out
# directive stays excluded rather than becoming eligible once its `//` is gone.
#
# Only comments that run to the end of the line are stripped: a `//` tail, a
# `/* ... */` that ends the line, and an unterminated `/*` tail. A closed
# `/* ... */` with code after it is deliberately left alone, because no portable
# BRE distinguishes it from the trailing case without risking the code that
# follows -- and getting that wrong drops a real macro, which is silent, while
# leaving a comment in only ever keeps an extra identifier, which is loud. Bias
# the residual case toward the failure a reader will see.
#
# No `\|` alternation, for the reason the `\b` note gives: it is a GNU extension
# that BSD sed reads as a literal `|`, and here that silently dropped the second
# macro of a two-macro condition. The `[[:space:]]*$` anchor and the negated
# `/\*\//` address below are portable BRE and behave identically on both.
# ---------------------------------------------------------------------------
src_macros="$(find "$src_dir" -type f \( \
        -name '*.c' -o -name '*.cc' -o -name '*.cpp' -o -name '*.cxx' -o \
        -name '*.h' -o -name '*.hh' -o -name '*.hpp' -o -name '*.hxx' \) \
        -exec cat {} + 2>/dev/null \
    | join_continuations \
    | grep -E '^[[:space:]]*#[[:space:]]*(if|ifdef|ifndef|elif)' \
    | sed -e 's://.*::' \
          -e 's:/\*.*\*/[[:space:]]*$: :' \
          -e '/\*\//!s:/\*.*$: :' \
    | grep -oE '[A-Za-z_][A-Za-z0-9_]*' | sort -u || true)"

if [[ -z $src_macros ]]; then
    echo "FAIL: no preprocessor conditionals found under $src_dir/ -- nothing was checked" >&2
    exit 1
fi

# Membership tests use a herestring rather than `printf ... | grep -q`. Under
# `set -o pipefail`, `grep -q` exits on its first match, printf is killed by
# SIGPIPE, and the pipeline reports 141 -- so a match reads as a miss. It only
# bites once the text outgrows the pipe buffer, which is exactly the kind of
# bug that lies dormant until the list it scans gets longer.
failures=0

# Forward: every flag the CI step passes must name a macro src/ actually tests.
while IFS= read -r flag; do
    [[ -z $flag ]] && continue
    if grep -qx -- "$flag" <<< "$src_macros"; then
        continue
    fi
    exempt=0
    # `${a[@]+"${a[@]}"}` rather than `"${a[@]}"`: the pending_flags block above
    # tells the maintainer to delete the entry when the check lands, and bash
    # before 4.4 (3.2 ships on macOS, and this script runs under `#!/bin/bash`)
    # treats the plain expansion of an empty array under `set -u` as an unbound
    # variable and aborts.
    for pending in ${pending_flags[@]+"${pending_flags[@]}"}; do
        [[ $flag == "$pending" ]] && exempt=1 && break
    done
    if (( exempt )); then
        echo "note: $flag names no macro in $src_dir/ yet (listed as pending)"
    else
        echo "FAIL: $workflow passes -D$flag, but no preprocessor conditional in" >&2
        echo "      $src_dir/ tests it. The opt-in build compiles it into nothing," >&2
        echo "      so whatever check it was meant to enable never runs. Fix the" >&2
        echo "      name, or add it to pending_flags with a reason if the check" >&2
        echo "      has not landed yet." >&2
        failures=$((failures + 1))
    fi
done <<< "$ci_flags"

# Reverse: every BWA_MEM3_DEBUG_* macro in src/ must be enabled by that step.
while IFS= read -r macro; do
    [[ -z $macro ]] && continue
    case "$macro" in BWA_MEM3_DEBUG_*) ;; *) continue ;; esac
    if grep -qx -- "$macro" <<< "$ci_flags"; then
        continue
    fi
    echo "FAIL: $src_dir/ gates code on $macro, but the opt-in macro build in" >&2
    echo "      $workflow does not pass -D$macro. No CI row compiles that code," >&2
    echo "      so it can rot unnoticed. Add it to the EXTRA_CXXFLAGS list." >&2
    failures=$((failures + 1))
done <<< "$src_macros"

if (( failures > 0 )); then
    echo "FAIL: debug_macro_flag_lint ($failures problem(s))" >&2
    exit 1
fi

echo "PASS: debug_macro_flag_lint (ci.yml -D list and src/ macros agree)"
