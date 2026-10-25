#!/usr/bin/env bash
# test/regression/make_c_opt_level.sh
#
# Regression: the optimization level on the aligner's compile and link lines.
#
# Was (1): the .c.o rule passed only $(CFLAGS), which the Makefile never gave an
# optimization level, so a default build compiled the FASTQ reader's C TUs
# (src/fast_reader*.c, src/fr_fastq.c) at the compiler default -O0 with no -g
# while every C++ TU got -g -O3. Nothing failed; the shipped reader was just
# unoptimized.
#
# Was (2): the aligner's compile rules and the bwa-mem3 link put $(CPPFLAGS)
# after the Makefile's -O level, and conda-forge's compiler activation exports
# CPPFLAGS="-DNDEBUG -D_FORTIFY_SOURCE=2 -O2", so a conda build compiled the
# aligner, SIMD kernels included, at -O2 (and an LTO link generated code at
# -O2). The same -DNDEBUG compiled every plain assert() out of that build,
# though no target here defines NDEBUG. The Makefile now strips -O levels and
# -DNDEBUG from CPPFLAGS.
#
# This pins the lines instead:
#   1. default build: the last -O on the reader's cc line is -O3, and -g is present;
#   2. an -O level in CFLAGS (environment) does not lower it;
#   3. ASAN=1: the last -O is the sanitizer's -O1;
#   4. EXTRA_CFLAGS=-O0 (the documented override): the last -O is -O0;
#   5. an -O level in CPPFLAGS (environment) does not lower the reader, a C++
#      TU, or a per-tier kernel TU below -O3, does not replace ASAN=1's -O1 or
#      COVERAGE=1's COVERAGE_OPT (default -O0), and does not lower the bwa-mem3
#      link line under LTO; the rest of CPPFLAGS (-D_FORTIFY_SOURCE=2) and the
#      Makefile's own defines (-DV17=1) are still on the line. A command-line
#      CPPFLAGS (the Makefile's `override` branch) also loses its -O and keeps
#      -D_FORTIFY_SOURCE=2;
#   6. a -DNDEBUG (or -DNDEBUG=<v>) in CPPFLAGS, environment or command line, is
#      stripped from the reader, C++, kernel and link lines, while the explicit
#      opt-out (EXTRA_CXXFLAGS=-DNDEBUG / EXTRA_CFLAGS=-DNDEBUG) still reaches
#      the line, and pgo-generate / pgo-use (EXTRA_CXXFLAGS) and lto-build
#      (EXTRA_CFLAGS) forward it to their sub-makes.
#
# Checked with `make -n -B` on single objects and on the bwa-mem3 link (with
# the prebuilt libraries held old), so it needs no built tree. The LTO link
# case sets BWA_LTO_MODE=full on the command line, which skips the toolchain
# probe, so it dry-runs on any host.
#
# Inputs:
#   MAKE      — make binary (default: make)
#   MAKE_ARGS — extra make arguments, e.g. "arch=avx2 CXX=clang++ CC=clang"
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."

MAKE="${MAKE:-make}"
read -r -a make_args <<< "${MAKE_ARGS:-}"

# The make target to dry-run and the text that identifies its line in the
# output: the source file for a compile, the output name for the link. The
# cases below reassign these to move from the reader to other objects.
OBJ=src/fr_fastq.o
SRC=src/fr_fastq.c
KIND=compile

# Dry-run $OBJ and set LINE to the $SRC compile (or link) line and OPT to the
# last -O token on it. Run in the current shell (not $(...)) so a FAIL can
# print it.
compile_line() { # <env assignments...> -- <make args...>
    local envs=()
    while [[ $1 != -- ]]; do
        envs+=("$1")
        shift
    done
    shift
    # For the link, hold the prebuilt libraries old (-o): under -n -B their
    # recipes would otherwise run, and recipe lines that call $(MAKE) execute
    # even under -n (the vendored cmake / configure sub-builds). The print-%
    # calls hold the included ext/htslib/htslib_static.mk old as well: on a
    # clean tree make would otherwise remake it, and so build htslib, before
    # printing the variable.
    local hold=() v w words
    if [[ $KIND == link ]]; then
        for v in BWA_LIB HTS_LIB LIBSAIS_OBJS MIMALLOC_LIB ZLIBNG_LIB; do
            words=$("$MAKE" -s -o ext/htslib/htslib_static.mk "${make_args[@]}" "print-$v" 2> /dev/null) || true
            # These three always name something to hold; empty means the
            # print-% helper or the variable changed, and the dry run would then
            # rebuild (and run the sub-builds of) the libraries.
            if [[ -z $words && $v =~ ^(BWA_LIB|HTS_LIB|LIBSAIS_OBJS)$ ]]; then
                echo "FAIL: \`make -s -o ext/htslib/htslib_static.mk print-$v\` printed nothing; cannot hold the prebuilt libraries for the link dry run"
                exit 1
            fi
            for w in $words; do
                hold+=(-o "$w")
            done
        done
    fi
    local out
    out=$(env -u CFLAGS -u CXXFLAGS -u CPPFLAGS -u EXTRA_CFLAGS -u EXTRA_CXXFLAGS -u COVERAGE_OPT -u LDFLAGS "${envs[@]}" "$MAKE" -n -B "${hold[@]}" "${make_args[@]}" "$@" "$OBJ" 2>&1) || true
    if [[ $KIND == link ]]; then
        LINE=$(grep -E -- " -o $SRC( |\$)" <<< "$out" | tail -n 1) || true
    else
        LINE=$(grep -F -- "$SRC" <<< "$out" | grep -E -- ' -c ' | tail -n 1) || true
    fi
    if [[ -z $LINE ]]; then
        echo "FAIL: no $KIND line for $SRC in \`make -n -B ${MAKE_ARGS:-} $* $OBJ\`"
        tail -n 20 <<< "$out"
        exit 1
    fi
    OPT=$(grep -oE -- '(^| )-O([0-3sgz]|fast)?( |$)' <<< "$LINE" | tr -d ' ' | tail -n 1) || true
}

check() { # <case> <expected last -O> <env assignments...> -- <make args...>
    local name=$1 want=$2
    shift 2
    compile_line "$@"
    if [[ $OPT != "$want" ]]; then
        echo "FAIL: $name: last optimization flag on the $KIND line for $SRC is '${OPT:-none}', expected $want"
        echo "  $LINE"
        exit 1
    fi
}

# The CPPFLAGS filter must remove only -O levels: the environment's other
# flags and the Makefile's own defines stay on the line.
check_keeps() { # <case>
    local f
    for f in -D_FORTIFY_SOURCE=2 -DV17=1; do
        if ! grep -qE -- "(^| )$f( |\$)" <<< "$LINE"; then
            echo "FAIL: $1: $f missing from the $KIND line for $SRC"
            echo "  $LINE"
            exit 1
        fi
    done
}

# The CPPFLAGS filter must also remove -DNDEBUG (conda's CPPFLAGS defines it):
# plain assert() stays live in every build of this repository without an
# explicit -DNDEBUG opt-out (EXTRA_CXXFLAGS=-DNDEBUG / EXTRA_CFLAGS=-DNDEBUG).
# Only the single-token -DNDEBUG / -DNDEBUG=<v> is stripped; a two-token
# `-D NDEBUG` in CPPFLAGS is not matched and still reaches the compiler.
check_no_ndebug() { # <case>
    if grep -qE -- '(^| )-DNDEBUG(=[^ ]*)?( |$)' <<< "$LINE"; then
        echo "FAIL: $1: -DNDEBUG from CPPFLAGS still on the $KIND line for $SRC"
        echo "  $LINE"
        exit 1
    fi
}

check_has_ndebug() { # <case>
    if ! grep -qE -- '(^| )-DNDEBUG( |$)' <<< "$LINE"; then
        echo "FAIL: $1: explicit -DNDEBUG opt-out missing from the $KIND line for $SRC"
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

CONDA_CPPFLAGS="-DNDEBUG -D_FORTIFY_SOURCE=2 -O2"
check "environment CPPFLAGS with -O2" -O3 "CPPFLAGS=$CONDA_CPPFLAGS" --
check_keeps "environment CPPFLAGS with -O2"
check_no_ndebug "environment CPPFLAGS with -O2"
check "EXTRA_CFLAGS=-DNDEBUG with environment CPPFLAGS" -O3 "CPPFLAGS=$CONDA_CPPFLAGS" -- EXTRA_CFLAGS=-DNDEBUG
check_has_ndebug "EXTRA_CFLAGS=-DNDEBUG with environment CPPFLAGS"
check "environment CPPFLAGS with -Ofast" -O3 "CPPFLAGS=-D_FORTIFY_SOURCE=2 -Ofast" --
check_keeps "environment CPPFLAGS with -Ofast"
check "environment CPPFLAGS with -O2, ASAN=1" -O1 "CPPFLAGS=$CONDA_CPPFLAGS" -- ASAN=1

OBJ=src/bwamem.o SRC=src/bwamem.cpp
check "C++ TU, environment CPPFLAGS with -O2" -O3 "CPPFLAGS=$CONDA_CPPFLAGS" --
check_keeps "C++ TU, environment CPPFLAGS with -O2"
check_no_ndebug "C++ TU, environment CPPFLAGS with -O2"
check "C++ TU, environment CPPFLAGS with -DNDEBUG=1" -O3 "CPPFLAGS=-DNDEBUG=1 -D_FORTIFY_SOURCE=2" --
check_no_ndebug "C++ TU, environment CPPFLAGS with -DNDEBUG=1"
check "C++ TU, EXTRA_CXXFLAGS=-DNDEBUG with environment CPPFLAGS" -O3 "CPPFLAGS=$CONDA_CPPFLAGS" -- EXTRA_CXXFLAGS=-DNDEBUG
check_has_ndebug "C++ TU, EXTRA_CXXFLAGS=-DNDEBUG with environment CPPFLAGS"
check "C++ TU, environment CPPFLAGS with -O2, COVERAGE=1" -O0 "CPPFLAGS=$CONDA_CPPFLAGS" -- COVERAGE=1
# A command-line CPPFLAGS takes the Makefile's `override` branch. It replaces
# the Makefile's own defines by design (-DV17=1 is not expected here), but its
# -O must still be dropped and the rest of it kept.
check "C++ TU, command-line CPPFLAGS with -O2" -O3 -- "CPPFLAGS=$CONDA_CPPFLAGS"
check_no_ndebug "C++ TU, command-line CPPFLAGS with -O2"
if ! grep -qE -- '(^| )-D_FORTIFY_SOURCE=2( |$)' <<< "$LINE"; then
    echo "FAIL: C++ TU, command-line CPPFLAGS with -O2: -D_FORTIFY_SOURCE=2 missing from the $KIND line for $SRC"
    echo "  $LINE"
    exit 1
fi

# The per-tier rule exists on every host (the x86 build is what uses it).
OBJ=src/bandedSWA.avx2.o SRC=src/bandedSWA.cpp
check "kernel TU (avx2 tier), environment CPPFLAGS with -O2" -O3 "CPPFLAGS=$CONDA_CPPFLAGS" --
check_keeps "kernel TU (avx2 tier), environment CPPFLAGS with -O2"
check_no_ndebug "kernel TU (avx2 tier), environment CPPFLAGS with -O2"

OBJ=bwa-mem3 SRC=bwa-mem3 KIND=link
check "bwa-mem3 link, BWA_LTO_MODE=full, environment CPPFLAGS with -O2" -O3 "CPPFLAGS=$CONDA_CPPFLAGS" -- BWA_LTO_MODE=full
if ! grep -qE -- '(^| )-flto=full( |$)' <<< "$LINE"; then
    echo "FAIL: bwa-mem3 link, BWA_LTO_MODE=full: no -flto=full on the link line"
    echo "  $LINE"
    exit 1
fi
check_no_ndebug "bwa-mem3 link, BWA_LTO_MODE=full, environment CPPFLAGS with -O2"

# The PGO and lto-build targets re-invoke make with their own EXTRA_CXXFLAGS on
# the sub-make's command line, which replaces the caller's value unless the
# recipe forwards it; lto-build also clears the command-line variables, so it
# must forward EXTRA_CFLAGS itself. Dry-run each target with MAKE replaced by
# a probe that prints its arguments (a $(MAKE) line runs even under -n, so
# this shows the sub-make's command line without building), check the caller's -DNDEBUG opt-out is forwarded, then
# dry-run src/bwamem.o with pgo-use's forwarded EXTRA_CXXFLAGS and check it
# reaches the compile line.
# A variable the caller set on the command line but the recipe does not pass
# reaches the sub-make through MAKEFLAGS instead, so MAKE is a shell that also
# prints the MAKEFLAGS it was started with into SUBFLAGS (in its own
# environment, so a recipe that clears MAKEFLAGS for the sub-make shows that),
# and the CPPFLAGS in that environment into SUBCPPFLAGS (unset prints <unset>).
submake_args() { # <target> <make args...>: set SUB to the echoed sub-make line, SUBFLAGS / SUBCPPFLAGS to its MAKEFLAGS / CPPFLAGS
    local target=$1 out
    shift
    # $$MAKEFLAGS, $${CPPFLAGS...} and $$@ are for make and the probe shell to
    # expand, not this one.
    # shellcheck disable=SC2016
    out=$(env -u CFLAGS -u CXXFLAGS -u CPPFLAGS -u EXTRA_CFLAGS -u EXTRA_CXXFLAGS "$MAKE" -n "${make_args[@]}" "$target" 'MAKE=sh -c '\''echo "SUBFLAGS=$$MAKEFLAGS"; echo "SUBCPPFLAGS=$${CPPFLAGS-<unset>}"; echo "$$@"'\'' submake' "$@" 2>&1) || true
    SUB=$(grep -E -- '(^| )all$' <<< "$out" | grep -v -e '^echo ' -e 'SUBFLAGS=' | tail -n 1) || true
    SUBFLAGS=$(grep -- '^SUBFLAGS=' <<< "$out" | tail -n 1) || true
    SUBCPPFLAGS=$(grep -- '^SUBCPPFLAGS=' <<< "$out" | tail -n 1) || true
    if [[ -z $SUB ]]; then
        echo "FAIL: $target: no sub-make line in \`make -n $target $*\` with the MAKE probe"
        tail -n 20 <<< "$out"
        exit 1
    fi
}
for target in pgo-generate pgo-use; do
    submake_args "$target" EXTRA_CXXFLAGS=-DNDEBUG EXTRA_CFLAGS=-DNDEBUG
    if ! grep -qE -- '(^| )EXTRA_CXXFLAGS=-DNDEBUG -fprofile-' <<< "$SUB"; then
        echo "FAIL: $target: the caller's EXTRA_CXXFLAGS=-DNDEBUG is not forwarded to its sub-make"
        echo "  $SUB"
        exit 1
    fi
    # The recipe passes no EXTRA_CFLAGS, so the caller's arrives via MAKEFLAGS;
    # an EXTRA_CFLAGS on the sub-make's own command line would replace it.
    if grep -qE -- '(^| )EXTRA_CFLAGS=' <<< "$SUB" \
        || ! grep -qE -- '(^| )EXTRA_CFLAGS=-DNDEBUG( |$)' <<< "$SUBFLAGS"; then
        echo "FAIL: $target: the caller's EXTRA_CFLAGS=-DNDEBUG does not reach its sub-make"
        echo "  $SUB"
        echo "  $SUBFLAGS"
        exit 1
    fi
done
PGO_USE_EXTRA=$(grep -oE -- 'EXTRA_CXXFLAGS=.* CXX=' <<< "$SUB" | sed -e 's/^EXTRA_CXXFLAGS=//' -e 's/ CXX=$//')
submake_args lto-build EXTRA_CFLAGS=-DNDEBUG
if ! grep -qE -- '(^| )EXTRA_CFLAGS=-DNDEBUG( |$)' <<< "$SUB"; then
    echo "FAIL: lto-build: the caller's EXTRA_CFLAGS=-DNDEBUG is not forwarded to its sub-make"
    echo "  $SUB"
    exit 1
fi
# lto-build clears MAKEFLAGS / MAKEOVERRIDES, so a command-line CPPFLAGS reaches
# its sub-make only through the environment: the filtered value must be
# exported (kept: -D_FORTIFY_SOURCE=2; stripped: -O2, -DNDEBUG).
submake_args lto-build "CPPFLAGS=$CONDA_CPPFLAGS"
if ! grep -qE -- '(^|[= ])-D_FORTIFY_SOURCE=2( |$)' <<< "$SUBCPPFLAGS" \
    || grep -qE -- '(^|[= ])(-O[^ ]*|-DNDEBUG(=[^ ]*)?)( |$)' <<< "$SUBCPPFLAGS"; then
    echo "FAIL: lto-build: a command-line CPPFLAGS does not reach its sub-make filtered"
    echo "  $SUBCPPFLAGS"
    exit 1
fi
OBJ=src/bwamem.o SRC=src/bwamem.cpp KIND=compile
check "C++ TU, pgo-use EXTRA_CXXFLAGS=-DNDEBUG, environment CPPFLAGS" -O3 "CPPFLAGS=$CONDA_CPPFLAGS" -- "EXTRA_CXXFLAGS=$PGO_USE_EXTRA"
check_has_ndebug "C++ TU, pgo-use EXTRA_CXXFLAGS=-DNDEBUG, environment CPPFLAGS"
if ! grep -qE -- '(^| )-fprofile-use=' <<< "$LINE"; then
    echo "FAIL: C++ TU, pgo-use: -fprofile-use missing from the compile line for $SRC"
    echo "  $LINE"
    exit 1
fi

echo "PASS: make_c_opt_level (reader C TUs: -g -O3 by default, environment CFLAGS -O0 overridden, ASAN=1 -O1, EXTRA_CFLAGS=-O0 honoured; an -O in CPPFLAGS ignored on reader, C++ and kernel TUs and the LTO link, with ASAN/COVERAGE levels and other flags kept; -DNDEBUG in CPPFLAGS stripped, EXTRA_CXXFLAGS/EXTRA_CFLAGS=-DNDEBUG honoured and forwarded through pgo-generate, pgo-use and lto-build; a command-line CPPFLAGS reaches lto-build's sub-make filtered)"
