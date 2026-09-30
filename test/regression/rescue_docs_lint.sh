#!/usr/bin/env bash
# test/regression/rescue_docs_lint.sh
#
# Keep the mate-rescue developer-guide pages (docs/src/developer-guide/rescue*.md),
# the BWA3_RESCUE_* knob list and the source pointing at each other.
#
# The rescue shortcuts are exact by argument, and the arguments live in prose:
# the pages, the knob list in src/rescue_env.h, and the user-facing knob table in
# docs/src/whats-different/performance.md. Prose does not fail a build when the
# code under it moves, so it rots quietly -- a renamed function the page still
# names, a knob whose default changed in one place, an exactness claim whose test
# was deleted. This lint makes each of those a failure:
#
#   1. The knob list parses (scripts/rescue_knobs.sh), and the committed
#      docs/_generated/rescue/knobs.md is what the generator prints now.
#   2. Every BWA3_RESCUE_* name src/ reads is in the list, and every listed knob
#      is read somewhere. Every read goes through a rescue_env.h reader, and the
#      listed default is the reader's: the literal of rescue_env_int (a negative
#      literal means a computed default, listed as a word such as `auto`), 1 for
#      rescue_env_on, 0 for rescue_env_opt_in.
#   3. The knob table in performance.md lists exactly the knobs in the list.
#   4. Every page is in docs/src/SUMMARY.md. Every backticked name in a page
#      exists: a path is a file or directory of the repository (a trailing
#      `:<line>` is dropped; a brace or star family is refused, since it cannot
#      be checked), a bare file name is a file under src/, test/ or scripts/, a
#      BWA3_RESCUE_* name is in the knob list, and an identifier (`a_b`, `A::b`,
#      `f()`) is a whole word somewhere under src/ or test/. Anything else in
#      backticks (an expression, an option, a value) is not checked.
#   5. Pages and source link both ways: every src/rescue_* file names at least
#      one page, every page is named by at least one source file, and every page
#      path a source file names exists.
#   6. Claims name their gates. A paragraph, list item or table row of a page
#      that says something is identical, unchanged or byte-for-byte, or gives the
#      same score, output, result or decision, must say
#      `Gate:` (or `Gates:`) followed by one or more backticked gates, and every
#      gate named anywhere must exist: `<name>.sh` is a script in test/regression/
#      that a workflow runs; anything else is the exact name of a unit TEST_CASE
#      under test/unit/ or of a workflow step.
#
# Source-only: no build, no binary. The optional argument aims the same checks
# at another tree, which is how rescue_docs_lint_selftest.sh shows each check
# still fails on the drift it exists to catch.

set -euo pipefail

if (($# > 1)); then
    echo "usage: ${BASH_SOURCE[0]##*/} [repo-root-to-check]" >&2
    exit 2
fi
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
root="$repo_root"
if (($# == 1)); then
    if [[ ! -d $1 ]]; then
        echo "FAIL: '$1' is not a directory -- nothing was checked" >&2
        exit 1
    fi
    root="$(cd "$1" && pwd)"
fi
cd "$root"

env_h="src/rescue_env.h"
generated="docs/_generated/rescue/knobs.md"
perf_md="docs/src/whats-different/performance.md"
summary="docs/src/SUMMARY.md"
guide_dir="docs/src/developer-guide"
for required in "$env_h" "$perf_md" "$summary"; do
    if [[ ! -f $required ]]; then
        echo "FAIL: $required is not a regular file under $PWD -- nothing was checked" >&2
        exit 1
    fi
done

failures=0
fail() {
    echo "FAIL: $*" >&2
    failures=$((failures + 1))
}

# --- 1. The knob list and its generated table. --------------------------------
# The generator is the one parser of the list; the lint reads the table it prints.
table="$(bash "$repo_root/scripts/rescue_knobs.sh" "$root")" || {
    echo "FAIL: scripts/rescue_knobs.sh could not render the knob list in $env_h" >&2
    exit 1
}
knob_rows="$(grep -E '^\| `BWA3_RESCUE_' <<< "$table" || true)"
if [[ -z $knob_rows ]]; then
    echo "FAIL: the knob list in $env_h renders no knob -- nothing was checked" >&2
    exit 1
fi
listed="$(cut -d'`' -f2 <<< "$knob_rows" | sort)"
listed_default() { grep -F "| \`$1\` |" <<< "$knob_rows" | cut -d'`' -f4; }

if [[ ! -f $generated ]]; then
    fail "$generated is missing; run \`make docs-rescue-knobs\` and commit it"
elif [[ "$(cat "$generated")" != "$table" ]]; then
    fail "$generated is stale; run \`make docs-rescue-knobs\` and commit it"
    diff <(printf '%s\n' "$table") "$generated" | head -20 >&2 || true
fi

# --- 2. Knobs read in src/, their readers and defaults. -----------------------
src_files=()
while IFS= read -r f; do
    [[ -n $f ]] && src_files+=("$f")
done < <(find src -maxdepth 1 -type f \( -name '*.c' -o -name '*.cpp' -o -name '*.h' \) | sort)
if ((${#src_files[@]} == 0)); then
    echo "FAIL: no sources under $PWD/src -- nothing was checked" >&2
    exit 1
fi

# C/C++ with its comments removed (line and block, a block may span lines) and its
# string and character literals kept, so a reader call, a TEST_CASE or a name inside
# a comment does not count as live code.
strip_c_comments() {
    awk '
    FNR == 1 { in_block = 0 }
    {
        n = length($0); out = ""; in_str = 0; q = ""
        for (i = 1; i <= n; i++) {
            c = substr($0, i, 1)
            if (in_block) {
                if (c == "*" && substr($0, i + 1, 1) == "/") { in_block = 0; i++; out = out " " }
                continue
            }
            if (in_str) {
                out = out c
                if (c == "\\" && i < n) { out = out substr($0, i + 1, 1); i++; continue }
                if (c == q) in_str = 0
                continue
            }
            if (c == "\"" || c == "\047") { in_str = 1; q = c; out = out c; continue }
            if (c == "/" && substr($0, i + 1, 1) == "/") break
            if (c == "/" && substr($0, i + 1, 1) == "*") { in_block = 1; i++; continue }
            out = out c
        }
        print out
    }' "$@"
}
# The sources' code on one line, so a call split across lines (or spaced before its
# parenthesis) reads like any other.
src_code="$(strip_c_comments "${src_files[@]}" | tr '\n\t' '  ')"
# The (knob, rest) calls of reader $1 in src_code, one "name<TAB>arguments" per line.
reader_calls() {
    { grep -oE "$1"'[[:space:]]*\([[:space:]]*"BWA3_RESCUE_[A-Z0-9_]+"[[:space:]]*(,[^)]*)?\)' <<< "$src_code" || true; } \
        | sed -E 's/^[^"]*"([A-Z0-9_]+)"[[:space:]]*,?[[:space:]]*([^)]*[^)[:space:]])?[[:space:]]*\)$/\1\t\2/'
}

# A knob counts as read only where a reader is called on it (or a bare getenv,
# refused below): a name quoted in a comment is not a read.
read_in_src="$(reader_calls '(rescue_env_(int|on|opt_in)|getenv)' | cut -f1 | sort -u)"
while IFS= read -r k; do
    [[ -z $k ]] && continue
    grep -qxF "$k" <<< "$listed" || fail "src/ reads $k, which is not in the knob list in $env_h"
done <<< "$read_in_src"
while IFS= read -r k; do
    [[ -z $k ]] && continue
    grep -qxF "$k" <<< "$read_in_src" || fail "the knob list names $k, which nothing in src/ reads"
done <<< "$listed"

bare="$(for f in "${src_files[@]}"; do
    strip_c_comments "$f" | tr '\n\t' '  ' | grep -qE 'getenv[[:space:]]*\([[:space:]]*"BWA3_RESCUE_' && echo "$f"
done | head -5 || true)"
[[ -z $bare ]] || fail "BWA3_RESCUE_* read with a bare getenv (use a rescue_env.h reader):"$'\n'"$bare"

while IFS=$'\t' read -r k lit; do
    [[ -z $k ]] && continue
    d="$(listed_default "$k")"
    if [[ $lit =~ ^[0-9]+$ ]]; then
        [[ $d == "$lit" ]] || fail "$k: the knob list says default $d, rescue_env_int reads it with $lit"
    elif [[ $lit =~ ^-[0-9]+$ ]]; then
        [[ ! $d =~ ^[0-9]+$ ]] \
            || fail "$k: rescue_env_int's default is computed ($lit), the knob list says the number $d"
    else
        fail "$k: rescue_env_int's default '$lit' is not an integer literal the lint can check"
    fi
done < <(reader_calls rescue_env_int | sort -u)
while IFS= read -r k; do
    [[ -z $k ]] && continue
    [[ "$(listed_default "$k")" == 1 ]] || fail "$k is read with rescue_env_on (default on); the knob list says $(listed_default "$k")"
done < <(reader_calls rescue_env_on | cut -f1 | sort -u)
while IFS= read -r k; do
    [[ -z $k ]] && continue
    [[ "$(listed_default "$k")" == 0 ]] || fail "$k is read with rescue_env_opt_in (default off); the knob list says $(listed_default "$k")"
done < <(reader_calls rescue_env_opt_in | cut -f1 | sort -u)

# --- 3. The user-facing knob table. --------------------------------------------
documented="$({ grep -oE '^\| `BWA3_RESCUE_[A-Z0-9_]+' "$perf_md" || true; } | cut -d'`' -f2 | sort -u)"
while IFS= read -r k; do
    [[ -z $k ]] && continue
    grep -qxF "$k" <<< "$documented" || fail "$perf_md's knob table has no row for $k"
done <<< "$listed"
while IFS= read -r k; do
    [[ -z $k ]] && continue
    grep -qxF "$k" <<< "$listed" || fail "$perf_md documents $k, which is not in the knob list"
done <<< "$documented"

# --- 4. Pages: in the SUMMARY, and every backticked name exists. ---------------
pages=()
while IFS= read -r f; do
    [[ -n $f ]] && pages+=("$f")
done < <(find "$guide_dir" -maxdepth 1 -type f -name 'rescue*.md' 2> /dev/null | sort)
if ((${#pages[@]} == 0)); then
    echo "FAIL: no $guide_dir/rescue*.md pages under $PWD -- nothing was checked" >&2
    exit 1
fi

# The text of a page outside fenced code blocks and {{#include}} lines.
prose() { awk '/^[[:space:]]*```/ { fence = !fence; next } !fence && !/^\{\{#include/' "$1"; }

# Code only: a comment naming a symbol -- in a regression script, in this lint, or in
# the C/C++ itself -- must not keep a page's reference to it alive after the symbol
# is gone.
tree_code="$(find src test -type f \( -name '*.c' -o -name '*.cpp' -o -name '*.h' \) -print0 2> /dev/null \
    | xargs -0 cat /dev/null | strip_c_comments)"
# A name also exists as the stem of a C/C++ file under src/ or test/ (a harness named
# after its source, such as a test binary).
word_in_tree() {
    grep -qw -- "$1" <<< "$tree_code" \
        || [[ -n "$(find src test -type f \( -name "$1.c" -o -name "$1.cpp" -o -name "$1.h" \) -print -quit 2> /dev/null)" ]]
}

for page in "${pages[@]}"; do
    grep -qF "(developer-guide/${page##*/})" "$summary" || fail "$page is not linked from $summary"
    # shellcheck disable=SC2016  # the backticks are Markdown code spans, not command substitution
    while IFS= read -r tok; do
        [[ -z $tok ]] && continue
        if [[ $tok == BWA3_RESCUE_* ]]; then
            k="${tok%%=*}"
            grep -qxF "$k" <<< "$listed" || fail "$page names $k, which is not in the knob list"
        elif [[ $tok =~ ^[A-Za-z0-9_.-]+(/[A-Za-z0-9_.-]*)+(:[0-9]+)?$ ]]; then
            path="$tok"
            [[ $tok =~ :[0-9]+$ ]] && path="${tok%:*}"
            [[ -e $path ]] || fail "$page names $tok, and $path does not exist"
        elif [[ $tok == */* && $tok =~ [*{] ]]; then
            fail "$page names the path family $tok; spell each path out so the lint can check it"
        elif [[ $tok =~ ^[A-Za-z0-9_-]+\.(sh|c|cpp|h)$ ]]; then
            [[ -n "$(find src test scripts -name "$tok" -print -quit 2> /dev/null)" ]] \
                || fail "$page names $tok, and no file of that name exists under src/, test/ or scripts/"
        elif [[ $tok =~ ^[A-Za-z_][A-Za-z0-9_]*(::[A-Za-z_][A-Za-z0-9_]*)*(\(\))?$ ]]; then
            ident="${tok%()}"
            IFS=':' read -r -a parts <<< "${ident//::/:}"
            for p in "${parts[@]}"; do
                [[ -z $p ]] && continue
                word_in_tree "$p" || fail "$page names \`$tok\`, and '$p' is not a word anywhere under src/ or test/"
            done
        fi
    done < <(prose "$page" | grep -oE '`[^`]+`' | tr -d '`' | sort -u)
done

# --- 5. Pages and source link both ways. ---------------------------------------
links="$({ grep -ohE 'docs/src/developer-guide/rescue[a-z0-9-]*\.md' "${src_files[@]}" || true; } | sort -u)"
while IFS= read -r l; do
    [[ -z $l ]] && continue
    [[ -f $l ]] || fail "src/ names $l, which does not exist"
done <<< "$links"
for page in "${pages[@]}"; do
    grep -qxF "$page" <<< "$links" || fail "no file under src/ names $page (put the path in the header of the code it documents)"
done
for f in "${src_files[@]}"; do
    [[ ${f##*/} == rescue_* ]] || continue
    grep -qE 'docs/src/developer-guide/rescue[a-z0-9-]*\.md' "$f" || fail "$f names no rescue developer-guide page"
done

# --- 6. Claims name their gates, and the gates exist. --------------------------
workflows=()
while IFS= read -r f; do
    [[ -n $f ]] && workflows+=("$f")
done < <(find .github/workflows -type f \( -name '*.yml' -o -name '*.yaml' \) 2> /dev/null | sort)
# Workflow text with comments removed, so a comment naming a script or step does
# not count as running it. The same quote-aware stripper as
# regression_coverage_lint.sh's strip_comments.
strip_comments() {
    awk '
    {
        n = length($0); in_single = 0; in_double = 0; kept = ""
        for (i = 1; i <= n; i++) {
            c = substr($0, i, 1)
            if (in_single) { if (c == "\047") in_single = 0; kept = kept c; continue }
            if (c == "\\" && i < n) { kept = kept c substr($0, i + 1, 1); i++; continue }
            if (in_double) { if (c == "\"") in_double = 0; kept = kept c; continue }
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
workflow_text="$(cat "${workflows[@]}" /dev/null | strip_comments)"
step_names="$({ grep -E '^[[:space:]]*-[[:space:]]+name:[[:space:]]' <<< "$workflow_text" || true; } \
    | sed -E -e 's/^[[:space:]]*-[[:space:]]+name:[[:space:]]*//' -e 's/[[:space:]]+$//' \
        -e 's/^"(.*)"$/\1/' -e "s/^'(.*)'\$/\1/" | sort -u)"
# A TEST_CASE inside a comment (line or block) is not a live test.
unit_tests=()
while IFS= read -r f; do
    [[ -n $f ]] && unit_tests+=("$f")
done < <(find test/unit -maxdepth 1 -type f -name '*.cpp' 2> /dev/null | sort)
test_cases="$({ strip_c_comments "${unit_tests[@]}" /dev/null | tr '\n\t' '  ' \
    | grep -oE 'TEST_CASE[[:space:]]*\([[:space:]]*"[^"]+"' || true; } \
    | sed -E 's/^TEST_CASE[[:space:]]*\([[:space:]]*"(.*)"$/\1/' | sort -u)"

gate_ok() {
    local g="$1"
    if [[ $g == *.sh ]]; then
        [[ -f test/regression/$g ]] || return 1
        grep -qF "test/regression/$g" <<< "$workflow_text"
        return
    fi
    grep -qxF -- "$g" <<< "$test_cases" && return 0
    grep -qxF -- "$g" <<< "$step_names"
}

# What counts as an exactness claim: the words the pages use for "the output does
# not change".
claim_words='identical|unchanged|byte-for-byte|byte for byte|(is|are) exact([^a-z]|$)|exactly|equivalent|the same (score|output|result|decision)s?'

# One block per paragraph, list item or table row, as "<first line number><TAB><text>".
blocks() {
    awk '
        function flush() { if (buf != "") print start "\t" buf; buf = "" }
        /^[[:space:]]*```/ { flush(); fence = !fence; next }
        fence { next }
        /^[[:space:]]*$/ || /^#/ || /^\{\{#include/ { flush(); next }
        /^[[:space:]]*([-*]|[0-9]+\.)[[:space:]]/ || /^\|/ { flush() }
        { if (buf == "") { start = NR; buf = $0 } else buf = buf " " $0 }
        END { flush() }
    ' "$1"
}

for page in "${pages[@]}"; do
    while IFS=$'\t' read -r line text; do
        [[ -z $line ]] && continue
        gates=""
        if [[ $text =~ Gates?:(.*)$ ]]; then
            # shellcheck disable=SC2016  # the backticks are Markdown code spans, not command substitution
            gates="$(grep -oE '`[^`]+`' <<< "${BASH_REMATCH[1]}" | tr -d '`' || true)"
        fi
        if grep -qiE "$claim_words" <<< "$text" && [[ -z $gates ]]; then
            fail "$page:$line makes an exactness claim without naming its gate (Gate: \`<test>\`)"
        fi
        while IFS= read -r g; do
            [[ -z $g ]] && continue
            gate_ok "$g" || fail "$page:$line names the gate '$g', which is no CI-run regression script, unit TEST_CASE or workflow step"
        done <<< "$gates"
    done < <(blocks "$page")
done

if ((failures > 0)); then
    echo "FAIL: rescue_docs_lint found $failures problem(s)" >&2
    exit 1
fi
echo "PASS: rescue_docs_lint (${#pages[@]} page(s), $(wc -l <<< "$listed" | tr -d ' ') knobs, generated table current, symbols, links and gates resolve)"
