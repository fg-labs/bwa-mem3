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
# Dry-run checks use `make -n -B` on objects and links, holding libraries old.
# Native checks also compile/run the differential reader, profiling off/on:
# selected compression libraries and their generated headers must already exist.
# The LTO link check forces BWA_LTO_MODE=full, so it dry-runs on every host.
#
# Inputs:
#   MAKE      — make binary (default: make)
#   MAKE_ARGS — extra make arguments, e.g. "arch=avx2 CXX=clang++ CC=clang".
#               Native probe preserves compiler/library configuration and the
#               include paths (-I/-isystem/-idirafter/-iquote) of CPPFLAGS from
#               the environment or MAKE_ARGS; overrides instrumentation, native
#               tuning, and sentinel CFLAGS/CPPFLAGS.
#   CC, CXX   — compiler defaults, overridden by MAKE_ARGS as with make itself.
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
    if [[ $KIND == link || $KIND == reader_diff ]]; then
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
    if [[ $KIND == link || $KIND == reader_diff ]]; then
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

# A distro compiler can default to x86-64-v3. Enabling SSE4.1 does not
# disable that baseline's AVX2: verify the actual predefined ISA macros for
# every tier under an intentionally stronger inherited baseline.
read -r -a tier_cxx <<< "$("$MAKE" -s -o ext/htslib/htslib_static.mk "${make_args[@]}" print-CXX)"
tier_target=()
tier_names=(sse41 sse42 avx avx2 avx512bw)
if [[ $(uname -m) == arm64 || $(uname -m) == aarch64 ]]; then
    # Clang can preprocess for x86 without an SDK; native ARM GCC cannot.
    tier_version=$("${tier_cxx[@]}" --version)
    if [[ $tier_version == *clang* ]]; then
        tier_target=(--target=x86_64-linux-gnu)
    else
        echo "SKIP: x86 kernel ISA macros require Clang when testing on ARM"
        tier_names=()
    fi
fi
for tier in "${tier_names[@]}"; do
    read -r -a tier_flags <<< "$("$MAKE" -s -o ext/htslib/htslib_static.mk "${make_args[@]}" "print-KERNEL_FLAGS_$tier")"
    tier_macros=$("${tier_cxx[@]}" "${tier_target[@]}" -march=x86-64-v4 \
        "${tier_flags[@]}" -dM -E -x c++ /dev/null)
    for feature in SSE4_1 SSE4_2 AVX AVX2 AVX512F AVX512BW; do
        want=0
        case "$tier:$feature" in
            *:SSE4_1 | sse42:SSE4_2 | avx:SSE4_2 | avx2:SSE4_2 | avx512bw:SSE4_2 | \
                avx:AVX | avx2:AVX | avx512bw:AVX | avx2:AVX2 | avx512bw:AVX2 | \
                avx512bw:AVX512F | avx512bw:AVX512BW) want=1 ;;
        esac
        got=0
        if grep -q "^#define __${feature}__ 1$" <<< "$tier_macros"; then got=1; fi
        if [[ $got != "$want" ]]; then
            echo "FAIL: $tier kernel: __${feature}__=$got, expected $want after resetting a stronger baseline"
            exit 1
        fi
    done
done
if ((${#tier_names[@]} > 0)); then
    echo "PASS: kernel tier ISA macros do not inherit a stronger compiler baseline"
fi

# The differential reader compiles C directly rather than using .c.o. Check
# that it sees the same effective preprocessor and ARM tuning flags as production.
OBJ=fr_fastq_diff_test SRC=fr_fastq_diff_test KIND=reader_diff
for tuning in ARM_CPU=neoverse-n1 ARM_TUNE=neoverse-v2; do
    check "differential reader, $tuning" -O3 'CPPFLAGS=-DREADER_FLAGS_PROBE=1' -- \
        UNAME_M=aarch64 UNAME_S=Linux ARM_TUNE= "$tuning" LTO=off
    # $(CFLAGS) precedes $(SRC_CFLAGS) on the reader line: an environment -O2
    # must not replace the reader's -O3.
    check "differential reader, $tuning, environment CFLAGS=-O2" -O3 \
        CFLAGS=-O2 'CPPFLAGS=-DREADER_FLAGS_PROBE=1' -- \
        UNAME_M=aarch64 UNAME_S=Linux ARM_TUNE= "$tuning" LTO=off
    case "$tuning" in
        ARM_CPU=*) tune_flag=-mcpu=${tuning#*=} ;;
        ARM_TUNE=*) tune_flag=-mtune=${tuning#*=} ;;
    esac
    for flag in -DREADER_FLAGS_PROBE=1 "$tune_flag"; do
        if [[ " $LINE " != *" $flag "* ]]; then
            echo "FAIL: differential reader: $flag missing from C invocation"
            echo "  $LINE"
            exit 1
        fi
    done
    echo "ok: differential reader flags: $LINE"
done

# Also execute the C recipe, not just its dry-run rendering. Keep its binary
# and flags stamp in a tiny workspace; source and dependency symlinks are read-only.
root=$PWD
probe_dir=$(mktemp -d -t bwamem3-reader-flags-XXXXXX)
trap 'rm -rf "$probe_dir"' EXIT
for path in src test ext scripts version.txt; do ln -s "$root/$path" "$probe_dir/$path"; done
cat > "$probe_dir/cc" << 'EOF'
#!/usr/bin/env bash
set -euo pipefail
printf 'reader C invocation:'
printf ' %s' "$@"
printf '\n'
read -r -a compiler <<< "$READER_REAL_CC"
case ${1:-} in --version | -dumpversion) exec "${compiler[@]}" "$@" ;; esac
for flag in -DREADER_FLAGS_PROBE=1 -DREADER_CFLAGS_PROBE=1 "$READER_NATIVE_FLAG"; do
    found=0
    for arg in "$@"; do [[ $arg == "$flag" ]] && found=1; done
    if ((found == 0)); then
        echo "FAIL: actual differential C invocation lacks $flag" >&2
        exit 1
    fi
done
exec "${compiler[@]}" "$@"
EOF
chmod +x "$probe_dir/cc"
export READER_REAL_CC
if ! compiler_output=$("$MAKE" -s -o ext/htslib/htslib_static.mk "${make_args[@]}" LTO=off print-CC); then
    echo 'FAIL: cannot resolve C compiler for selected MAKE_ARGS' >&2
    exit 1
fi
READER_REAL_CC=$(tail -n 1 <<< "$compiler_output")
if [[ -z $READER_REAL_CC ]]; then
    echo 'FAIL: selected make configuration has no C compiler' >&2
    exit 1
fi
case $(uname -m) in
    arm64 | aarch64)
        native_settings=(ARM_CPU=native)
        export READER_NATIVE_FLAG=-mcpu=native
        ;;
    *)
        native_settings=(arch=native)
        export READER_NATIVE_FLAG=-march=native
        ;;
esac
# Replace the caller's CPPFLAGS (environment or MAKE_ARGS) with the probe
# sentinels so profiling and flag checks stay isolated, but keep its include
# paths: a system zlib-ng header may be found only through them.
native_args=()
caller_cppflags=${CPPFLAGS:-}
for arg in "${make_args[@]}"; do
    case $arg in
        CPPFLAGS=*) caller_cppflags+=" ${arg#CPPFLAGS=}" ;;
        CFLAGS=*) ;;
        *) native_args+=("$arg") ;;
    esac
done
read -r -a caller_cppflags_words <<< "$caller_cppflags"
include_cppflags=()
for ((i = 0; i < ${#caller_cppflags_words[@]}; i++)); do
    word=${caller_cppflags_words[i]}
    case $word in
        -I | -isystem | -idirafter | -iquote)
            include_cppflags+=("$word")
            ((i + 1 < ${#caller_cppflags_words[@]})) && include_cppflags+=("${caller_cppflags_words[++i]}")
            ;;
        -I* | -isystem* | -idirafter* | -iquote*) include_cppflags+=("$word") ;;
    esac
done
(
    cd "$probe_dir"
    for profiling in 0 1; do
        rm -f fr_fastq_diff_test
        profile_cppflags=
        [[ $profiling == 1 ]] && profile_cppflags=-DSTAGE_PROF=1
        CPPFLAGS="${include_cppflags[*]} -DREADER_FLAGS_PROBE=1 $profile_cppflags" "$MAKE" -s -f "$root/Makefile" "${native_args[@]}" \
            -o ext/htslib/htslib_static.mk -o ext/zlib-ng/build/libz-ng.a \
            fr_fastq_diff_test "CC=$probe_dir/cc" \
            CFLAGS=-DREADER_CFLAGS_PROBE=1 \
            VERSION_STRING=reader-flags LTO=off ARM_TUNE= "${native_settings[@]}" "STAGE_PROF=$profiling"
        ./fr_fastq_diff_test
        echo "ok: executed reader recipe STAGE_PROF=$profiling with selected MAKE_ARGS compiler/library mode"
    done
)
echo "PASS: make_c_opt_level (production and differential reader flags, executed native C recipe, optimization levels, CPPFLAGS filtering, recursive forwarding and kernel ISA isolation)"
