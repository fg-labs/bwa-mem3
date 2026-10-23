#!/usr/bin/env bash
# test/regression/make_c_opt_level.sh
#
# Regression: the FASTQ reader's C TUs (src/fast_reader*.c, src/fr_fastq.c)
# compile at the C++ TUs' -O3 with debug info.
#
# Was: the .c.o rule passed only $(CFLAGS), which the Makefile never gave an
# optimization level, so a default build compiled the reader at the compiler
# default -O0 with no -g while every C++ TU got -g -O3. Nothing failed; the
# shipped reader was just unoptimized. This pins the compile line instead:
#
#   1. default build: the last -O on the cc line is -O3, and -g is present;
#   2. an -O level in CFLAGS (environment) does not lower it;
#   3. ASAN=1: the last -O is the sanitizer's -O1;
#   4. EXTRA_CFLAGS=-O0 (the documented override): the last -O is -O0.
#
# Checked with `make -n -B` on one reader object, so it needs no built tree.
#
# Inputs:
#   MAKE      — make binary (default: make)
#   MAKE_ARGS — extra make arguments, e.g. "arch=avx2 CXX=clang++ CC=clang"
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."

MAKE="${MAKE:-make}"
read -r -a make_args <<< "${MAKE_ARGS:-}"

# Dry-run src/fr_fastq.o and set LINE to its compile line and OPT to the last
# -O token on it. Run in the current shell (not $(...)) so a FAIL can print it.
compile_line() { # <env assignments...> -- <make args...>
    local envs=()
    while [[ $1 != -- ]]; do
        envs+=("$1")
        shift
    done
    shift
    local out
    out=$(env -u CFLAGS -u EXTRA_CFLAGS "${envs[@]}" "$MAKE" -n -B "${make_args[@]}" "$@" src/fr_fastq.o 2>&1) || true
    LINE=$(grep -E -- '-c .*src/fr_fastq\.c' <<< "$out" | tail -n 1) || true
    if [[ -z $LINE ]]; then
        echo "FAIL: no compile line for src/fr_fastq.c in \`make -n -B ${MAKE_ARGS:-} $* src/fr_fastq.o\`"
        tail -n 20 <<< "$out"
        exit 1
    fi
    OPT=$(grep -oE -- '(^| )-O[0-3sgz]?( |$)' <<< "$LINE" | tr -d ' ' | tail -n 1) || true
}

check() { # <case> <expected last -O> <env assignments...> -- <make args...>
    local name=$1 want=$2
    shift 2
    compile_line "$@"
    if [[ $OPT != "$want" ]]; then
        echo "FAIL: $name: last optimization flag on the reader's compile line is '${OPT:-none}', expected $want"
        echo "  $LINE"
        exit 1
    fi
}

check "default build" -O3 --
if ! grep -qE -- '(^| )-g( |$)' <<< "$LINE"; then
    echo "FAIL: default build: no -g on the reader's compile line"
    echo "  $LINE"
    exit 1
fi
check "environment CFLAGS=-O0" -O3 CFLAGS=-O0 --
check "ASAN=1" -O1 -- ASAN=1
check "EXTRA_CFLAGS=-O0" -O0 -- EXTRA_CFLAGS=-O0

echo "PASS: make_c_opt_level (reader C TUs: -g -O3 by default, environment CFLAGS -O0 overridden, ASAN=1 -O1, EXTRA_CFLAGS=-O0 honoured)"
