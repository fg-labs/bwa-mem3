#!/bin/bash
# Reject NDEBUG preprocessor gates in production sources under src/ -- any
# `#if`, `#ifdef`, `#ifndef` or `#elif` directive that tests NDEBUG.
#
# Why this is a hard error rather than a style preference: no target in this
# build system has ever defined NDEBUG -- not the Makefile, not any workflow in
# .github/, including the release build. So a gate written as
#
#     #ifndef NDEBUG
#         ... expensive debug-only work ...
#     #endif
#
# never compiles out of anything this repository produces. The work runs in
# every shipped binary while the surrounding comment claims it does not. That
# has bitten this codebase repeatedly:
#
#   - the pac-fetch poison in bntseq.cpp, which memset whole buffers on every
#     fetch of a production run before it was converted to BWA_MEM3_DEBUG_POISON;
#   - a chaining cross-check that ran the O(n^2) reference pass -- the exact work
#     its own optimization existed to skip -- on every read.
#
# The idiom that actually works here is a dedicated opt-in macro, off by
# default and named in the comment so it can be turned on deliberately:
#
#     #ifdef BWA_MEM3_DEBUG_MY_CHECK
#         ... expensive debug-only work ...
#     #endif
#
# built with `make EXTRA_CXXFLAGS=-DBWA_MEM3_DEBUG_MY_CHECK`. Existing examples:
# BSW8_ASSERT_ENVELOPE, BWA_MEM3_DEBUG_POISON, BWAMEM3_UGP_PROFILE.
#
# Scope is src/ only. test/ is exempt: test/unit/test_xassert_macro.cpp gates on
# NDEBUG deliberately, because asserting on assert-vs-xassert semantics is the
# whole point of that test.
#
# Note this lint says nothing about plain `assert()` calls, which are live in
# every build here and are fine. `xassert()` (utils.h) remains the right guard
# for a check that must survive a hypothetical -DNDEBUG build.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

if (( $# > 1 )); then
    echo "usage: ${BASH_SOURCE[0]##*/} [directory-to-scan]" >&2
    exit 2
fi

# The directory to scan. Defaults to this repository's src/; the optional
# argument exists so ndebug_gate_lint_selftest.sh can aim the same matcher at a
# fixture tree of known-good and known-bad directives, which is the only way to
# tell "no gates in src/" apart from "the matcher stopped matching".
scan_root="${1:-src}"

# Every way this lint can end up scanning nothing has to be a failure, not a
# PASS. A lint that silently checks nothing is the same class of bug as the
# NDEBUG gate it exists to catch: it reads as green while doing no work.
if [[ ! -d $scan_root ]]; then
    echo "FAIL: no $scan_root/ directory under $repo_root -- nothing was linted" >&2
    exit 1
fi

# Select by extension rather than scanning every regular file: in a working
# tree src/ also holds untracked build artifacts (.o and friends), and feeding
# those to a text scanner is a hard error, not a lint finding. The list is
# deliberately wider than the extensions the tree uses today so adding a .cc or
# a .inc does not quietly fall outside the lint.
sources=()
while IFS= read -r -d '' file; do
    sources+=("$file")
done < <(find "$scan_root" -type f \( \
    -name '*.c' -o -name '*.cc' -o -name '*.cpp' -o -name '*.cxx' -o \
    -name '*.h' -o -name '*.hh' -o -name '*.hpp' -o -name '*.hxx' -o \
    -name '*.inc' -o -name '*.ipp' -o -name '*.tcc' \) -print0)

if (( ${#sources[@]} == 0 )); then
    echo "FAIL: no C/C++ sources under $scan_root/ -- nothing was linted" >&2
    exit 1
fi

# Flag only real preprocessor conditionals: a line whose first non-blank
# character is `#`, an if/ifdef/ifndef/elif directive, mentioning NDEBUG as a
# whole word. Prose inside block comments (` * ... #ifndef NDEBUG ...`) does not
# start with `#` and so does not match -- bntseq.cpp documents the historical
# bug in exactly that shape and must keep passing.
#
# Two details this gets right that a plain `grep -E` does not:
#
#   - Word boundaries are spelled out as explicit character classes. `\b` is a
#     GNU extension, not POSIX ERE; an implementation without it treats `\b` as
#     a literal `b`, so `(if|ifdef|ifndef|elif)\b` matches nothing at all and
#     the lint passes everything -- silently, on exactly the machines whose
#     grep differs from CI's.
#
#   - A trailing backslash splices the next physical line onto this one before
#     the compiler sees a directive, so `#if \` / `defined(NDEBUG)` is one
#     `#if defined(NDEBUG)`. grep matches physical lines and would miss it;
#     this joins continuations first and reports the line the directive starts
#     on.
# shellcheck disable=SC2016  # single quotes are deliberate: $0 is awk's, not the shell's
scan_file='
{
    line = $0
    start = FNR
    while (line ~ /\\[ \t]*$/) {
        sub(/\\[ \t]*$/, "", line)
        spliced = (getline continuation)
        # getline returns 0 at end of file and -1 on a read error. Treating the
        # error as end of file would abandon the rest of the file and still let
        # the lint print PASS.
        if (spliced < 0) {
            print "awk: read error while splicing a continued line" > "/dev/stderr"
            exit 2
        }
        if (spliced == 0) break
        line = line continuation
    }
    if (line ~ /^[ \t]*#[ \t]*(ifndef|ifdef|elif|if)([^A-Za-z0-9_]|$)/ &&
        line ~ /(^|[^A-Za-z0-9_])NDEBUG([^A-Za-z0-9_]|$)/) {
        printf "%s:%d:%s\n", FILENAME, start, line
    }
}'

# One awk per file so a continuation at end-of-file cannot splice across into
# the next one, and so a read error on any single file fails the whole lint
# instead of being mistaken for "that file had no matches".
hits=""
for file in "${sources[@]}"; do
    if ! file_hits="$(awk "$scan_file" "$file")"; then
        echo "FAIL: unable to scan $file -- lint did not complete" >&2
        exit 1
    fi
    if [[ -n $file_hits ]]; then
        hits+="$file_hits"$'\n'
    fi
done

if [[ -n $hits ]]; then
    echo "FAIL: NDEBUG preprocessor gate(s) found in $scan_root/:" >&2
    echo >&2
    printf '%s' "$hits" >&2
    echo >&2
    echo "This build never defines NDEBUG, so the guarded code runs in every" >&2
    echo "binary produced here -- including the release build. Replace the gate" >&2
    echo "with a dedicated opt-in macro (see BWA_MEM3_DEBUG_POISON in" >&2
    echo "src/bntseq.cpp) and document how to enable it." >&2
    exit 1
fi

echo "PASS: ndebug_gate_lint (no NDEBUG preprocessor gates in $scan_root/)"
