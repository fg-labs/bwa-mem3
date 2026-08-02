#!/usr/bin/env bash
# test/regression/readme_contract_lint.sh
#
# Check that test/regression/README.md still describes the scripts next to it.
#
# The README is prose plus a table, and the two rot differently. A new script
# reliably gets a table row, because adding one is part of writing the script.
# The prose around the table -- which jobs run what, and which scripts take
# environment variables versus a positional directory -- reliably does not,
# because nothing points at it. That has happened twice:
#
#   - the README claimed every script "reads its inputs from environment
#     variables" for as long as ndebug_gate_lint.sh, which reads none, had
#     existed;
#   - when the debug_macro_flag_lint pair landed it got its two table rows,
#     while the intro kept naming "the two ndebug_gate_lint* scripts" and the
#     single job that ran them.
#
# Both were caught by a human reading the file for an unrelated reason. This
# lint makes the two claims that actually drifted checkable:
#
#   1. Every test/regression/*.sh the README names still exists. The reverse --
#      requiring every script to appear -- is deliberately NOT checked: the
#      README says in as many words that its table is "a reading guide, not an
#      inventory", and a lint that forced it to become one would be arguing
#      with the file it is checking. A row for a script that was deleted is a
#      different thing: it is not a judgement call about coverage, it is a
#      broken reference. The developer guide carried one for the better part of
#      a year, telling readers to run a phix_parity.sh that no longer existed.
#
#   2. The set of scripts that read no environment at all is exactly the set
#      named in the README's source-only-lint block. Add a lint and forget the
#      prose, or give an existing lint an env var without moving it out, and
#      this fails.
#
# Check 2 needs to know where the claim lives, so the README marks it:
#
#     <!-- source-only lints: begin -->
#     ... naming each script that takes no environment ...
#     <!-- source-only lints: end -->
#
# The markers are the contract. Prose inside them can be rewritten freely --
# only the script names it mentions are checked.
#
# "Reads no environment" is decided by reading the script: a variable that is
# referenced but never assigned anywhere in the file is an input from the
# environment. Shell-provided names (BASH_SOURCE, PWD, IFS, ...) are excluded
# by the allowlist below, since no caller sets those to configure a run.
#
# Input:
#   $1  optional repository root to check; defaults to this repository. The
#       argument exists so readme_contract_lint_selftest.sh can aim the same
#       checks at a fixture tree of known-good and known-bad input, which is
#       the only way to tell "the README is accurate" apart from "the checks
#       stopped checking".

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

if (($# > 1)); then
    echo "usage: ${BASH_SOURCE[0]##*/} [repo-root-to-check]" >&2
    exit 2
fi

# Checked rather than left to `cd` to fail under `set -e`: bash's own `cd:` line
# carries no FAIL: prefix, so the one way of pointing this lint at nothing that
# is reachable from the command line would be the one failure a caller scanning
# for the family's format does not see.
if (($# == 1)); then
    if [[ ! -d $1 ]]; then
        echo "FAIL: '$1' is not a directory -- nothing was checked" >&2
        exit 1
    fi
    cd "$1"
fi

regression_dir="test/regression"
readme="$regression_dir/README.md"

# Every way this lint can end up checking nothing has to be a failure, not a
# PASS. A lint that silently checks nothing reads as green while doing no work,
# which is the same class of bug as the drift it exists to catch.
if [[ ! -d $regression_dir ]]; then
    echo "FAIL: no $regression_dir/ under $PWD -- nothing was checked" >&2
    exit 1
fi
if [[ ! -f $readme ]]; then
    echo "FAIL: $readme missing under $PWD -- nothing was checked" >&2
    exit 1
fi

scripts=()
while IFS= read -r path; do
    scripts+=("${path##*/}")
done < <(find "$regression_dir" -maxdepth 1 -type f -name '*.sh' | sort)

if ((${#scripts[@]} == 0)); then
    echo "FAIL: no *.sh under $regression_dir/ -- nothing was checked" >&2
    exit 1
fi

failures=0

# What counts as a reference to a script in this README. Both checks read names
# out of the file with it -- check 1 from the whole README, check 2 from inside
# the source-only block -- so it is written once: two copies would let check 2
# go on matching an older shape after check 1 was widened, which is this lint's
# own failure mode turned inward.
script_name_pattern='[A-Za-z_][A-Za-z0-9_.-]*\.sh'

# --- Check 1: every script the README names still exists. ------------------

# Names are taken from the whole file rather than just the table, so a script
# retired out of a prose sentence is caught the same way as one retired out of
# a row. A bare `*.sh` cannot come out of the pattern above, which has to start
# on a letter or underscore, so a glob like `test/regression/*.sh` yields no
# name at all; the skip below states that intent for anyone who widens it.
mentioned=()
while IFS= read -r name; do
    [[ $name == '*.sh' ]] && continue
    mentioned+=("$name")
done < <(grep -oE "$script_name_pattern" "$readme" | sort -u)

if ((${#mentioned[@]} == 0)); then
    echo "FAIL: $readme names no *.sh at all -- check 1 detected nothing" >&2
    exit 1
fi

dangling=()
for name in "${mentioned[@]}"; do
    if [[ ! -f "$regression_dir/$name" ]]; then
        dangling+=("$name")
    fi
done

if ((${#dangling[@]} > 0)); then
    echo "FAIL: $readme names scripts that do not exist under $regression_dir/:" >&2
    printf '  %s\n' "${dangling[@]}" >&2
    echo >&2
    echo "A reader following the README runs a script that is not there. Drop the" >&2
    echo "reference, or point it at whatever replaced the script." >&2
    failures=1
fi

# --- Check 2: the source-only-lint block lists exactly the env-free scripts. -

begin_marker='<!-- source-only lints: begin -->'
end_marker='<!-- source-only lints: end -->'

for marker in "$begin_marker" "$end_marker"; do
    count="$(grep -cF -- "$marker" "$readme" || true)"
    if [[ $count != 1 ]]; then
        echo "FAIL: $readme has $count copies of '$marker' -- expected exactly 1" >&2
        echo "      (check 2 cannot locate the source-only-lint block)" >&2
        exit 1
    fi
done

begin_line="$(grep -nF -- "$begin_marker" "$readme" | cut -d: -f1)"
end_line="$(grep -nF -- "$end_marker" "$readme" | cut -d: -f1)"

# Ordering first, because the emptiness test below cannot tell the two apart:
# an end marker above its begin marker also leaves nothing between them, and
# reporting that as an empty block sends the reader to look at contents that are
# right there in the file.
if ((end_line < begin_line)); then
    echo "FAIL: the source-only-lint markers in $readme are the wrong way round --" >&2
    echo "      end on line $end_line, begin on line $begin_line" >&2
    exit 1
fi

if ((end_line <= begin_line + 1)); then
    echo "FAIL: the source-only-lint block in $readme is empty -- nothing to check" >&2
    exit 1
fi

block="$(sed -n "$((begin_line + 1)),$((end_line - 1))p" "$readme")"

# Names the shell itself provides. A script reading only these takes no
# configuration from its caller, so it still counts as source-only.
shell_provided_names='^(BASH|BASHPID|BASH_ARGC|BASH_ARGV|BASH_COMMAND|BASH_LINENO|BASH_REMATCH|BASH_SOURCE|BASH_SUBSHELL|BASH_VERSINFO|BASH_VERSION|COLUMNS|EUID|FUNCNAME|GROUPS|HISTFILE|HOME|HOSTNAME|HOSTTYPE|IFS|LINENO|LINES|MACHTYPE|OLDPWD|OPTARG|OPTIND|OSTYPE|PATH|PIPESTATUS|PPID|PS1|PS2|PS4|PWD|RANDOM|REPLY|SECONDS|SHELL|SHLVL|UID)$'

# Emit only the parts of a script where a `$NAME` would actually expand:
# comments, single-quoted spans, and backslash-escaped dollars are dropped.
#
# This is the difference between reading a script and grepping it. Two of these
# scripts hold shell text as data -- fixture bodies for their self-tests --
# inside single quotes that may run over several lines, so a line-at-a-time
# filter cannot see where the quote opened. `'${MAKE:-make}'` is a string; the
# same characters unquoted are an environment read. Tracking the quote state
# across the whole file is what tells them apart.
#
# Not handled: quoted heredocs (`<<'EOF'`), whose body is also literal. Their
# contents are still scanned, which can only over-detect inputs -- a script
# would be kept out of the source-only block rather than slipped into it, so
# the failure is loud.
strip_non_expanding() {
    # Both quote states are tracked, not just single. A single quote inside a
    # double-quoted string is a literal character and opens nothing -- the
    # `'"'"'` idiom these scripts use to embed a quote is exactly that shape,
    # and reading its quote as an opener inverts the state for the whole rest
    # of the file.
    awk '
    BEGIN { single = 0; double = 0 }
    {
        line = $0
        out = ""
        i = 1
        n = length(line)
        while (i <= n) {
            c = substr(line, i, 1)
            if (single) {
                if (c == "\047") single = 0
                i++
                continue
            }
            if (double) {
                # Inside double quotes a backslash still escapes, and `$` still
                # expands, so the contents are kept.
                if (c == "\\") { out = out " "; i += 2; continue }
                if (c == "\"") double = 0
                else out = out c
                i++
                continue
            }
            if (c == "\047") { single = 1; i++; continue }
            if (c == "\"") { double = 1; i++; continue }
            # A backslash quotes the next character, so `\$FOO` is a literal
            # dollar sign. Drop both, leaving a space so adjacent words do not
            # fuse into a new token.
            if (c == "\\") { out = out " "; i += 2; continue }
            # A `#` starting a word begins a comment that runs to end of line.
            if (c == "#" && (i == 1 || substr(line, i - 1, 1) ~ /[ \t;&|()]/)) break
            out = out c
            i++
        }
        print out
    }' "$1"
}

# A variable that is referenced but never assigned anywhere in the file is an
# input from the environment. Assignment covers the plain, prefix, `local`,
# `export`, `declare`/`readonly`/`typeset`, `for NAME in`, and `read -r NAME`
# forms -- anything the script sets for itself is not something a caller passes
# in.
#
# One shape defeats that rule and has to be recognised on its own:
#
#     MAKE="${MAKE:-make}"
#
# is an assignment, but the value assigned is whatever the caller exported. The
# `:-` `-` `:=` `:?` `?` `:+` `+` expansions exist precisely to consume an
# outside value, so a name used in one is an environment input no matter what
# else the script does with it. Without this, make_default_goal.sh and
# make_header_deps.sh -- whose only inputs are MAKE and MAKE_ARGS, both
# defaulted this way -- would be misread as source-only.
env_inputs_of() {
    local file="$1"

    local code
    code="$(strip_non_expanding "$file")"

    local referenced assigned defaulted
    referenced="$(grep -oE '\$\{[A-Za-z_][A-Za-z0-9_]*|\$[A-Za-z_][A-Za-z0-9_]*' <<< "$code" \
        | sed -E 's/^\$\{?//' | sort -u || true)"
    assigned="$({
        grep -oE '(^|[[:space:]]|[;&|(])[A-Za-z_][A-Za-z0-9_]*=' <<< "$code" \
            | sed -E 's/^[^A-Za-z_]*//; s/=$//'
        grep -oE '\bfor[[:space:]]+[A-Za-z_][A-Za-z0-9_]*' <<< "$code" \
            | sed -E 's/^for[[:space:]]+//'
        grep -oE '\bread[[:space:]]+(-[A-Za-z]+[[:space:]]+)*[A-Za-z_][A-Za-z0-9_]*' <<< "$code" \
            | sed -E 's/.*[[:space:]]//'
    } | sort -u || true)"
    defaulted="$(grep -oE '\$\{[A-Za-z_][A-Za-z0-9_]*:?[-=?+]' <<< "$code" \
        | sed -E 's/^\$\{//; s/:?[-=?+]$//' | sort -u || true)"

    {
        comm -23 <(printf '%s\n' "$referenced") <(printf '%s\n' "$assigned")
        printf '%s\n' "$defaulted"
    } | sort -u | grep -E '^[A-Z][A-Z0-9_]*$' | grep -vE "$shell_provided_names" || true
}

source_only=()
env_reading=()
for script in "${scripts[@]}"; do
    if [[ -z "$(env_inputs_of "$regression_dir/$script")" ]]; then
        source_only+=("$script")
    else
        env_reading+=("$script")
    fi
done

# Not a possible state today, and a bug in env_inputs_of if it ever is: it
# would make check 2 vacuous while still printing PASS.
if ((${#source_only[@]} == 0)); then
    echo "FAIL: no script under $regression_dir/ reads zero environment variables" >&2
    echo "      -- check 2 detected nothing, which means it is no longer checking" >&2
    exit 1
fi

# Membership is decided against the block's script names, extracted with the
# shared pattern above, rather than against its raw text. A substring search
# reads one name as another whenever one is a suffix of the other: with
# `other_lint.sh` in the block, `grep -qF lint.sh` matches its row, so an
# omitted source-only `lint.sh` reports as present and an env-reading `lint.sh`
# reports as wrongly listed. Both directions are wrong, and the first is silent.
block_names="$(grep -oE "$script_name_pattern" <<< "$block" | sort -u || true)"

block_names_contains() {
    grep -qxF -- "$1" <<< "$block_names"
}

missing_from_block=()
for script in "${source_only[@]}"; do
    if ! block_names_contains "$script"; then
        missing_from_block+=("$script")
    fi
done

wrongly_in_block=()
for script in "${env_reading[@]}"; do
    if block_names_contains "$script"; then
        wrongly_in_block+=("$script")
    fi
done

if ((${#missing_from_block[@]} > 0)); then
    echo "FAIL: these scripts read no environment but the source-only-lint block" >&2
    echo "      in $readme does not name them:" >&2
    printf '  %s\n' "${missing_from_block[@]}" >&2
    echo >&2
    echo "The block is what tells a reader which scripts need no setup. Name them" >&2
    echo "there, with the positional argument each one takes." >&2
    failures=1
fi

if ((${#wrongly_in_block[@]} > 0)); then
    echo "FAIL: the source-only-lint block in $readme names these scripts, but they" >&2
    echo "      do read environment variables:" >&2
    for script in "${wrongly_in_block[@]}"; do
        echo "  $script ($(env_inputs_of "$regression_dir/$script" | tr '\n' ' '))" >&2
    done
    echo >&2
    echo "Either move them out of the block or drop the environment input." >&2
    failures=1
fi

if ((failures != 0)); then
    exit 1
fi

echo "PASS: readme_contract_lint (${#mentioned[@]} scripts named and present; ${#source_only[@]} source-only)"
