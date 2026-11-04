#!/usr/bin/env bash
# test/regression/make_lto_mode.sh
#
# Regression: how the Makefile resolves LTO= into the mode the build uses
# (BWA_LTO_MODE), read back with `make -s print-BWA_LTO_MODE`. That goal runs the
# same parse-time logic as a real build, including its small capability probes.
#
# The default (LTO=auto) must never fail a build: it picks full LTO only when
# the toolchain can do it and otherwise builds without LTO, saying so on stderr.
# An explicit LTO=full|thin must fail loudly instead. Checked:
#   - LTO= validation: an unknown value and a multi-word value fail; an invalid
#     value inherited from the environment is ignored with a warning.
#   - LTO=off, ASAN=1, and a gcc toolchain resolve to off; LTO=full with gcc fails.
#   - On Linux: a clang that cannot link with lld resolves auto to off with a
#     NOTE on stderr and nothing extra on stdout, and fails LTO=full; a clang
#     CXX with a non-clang CC resolves auto to off with a NOTE.
#   - On macOS: auto resolves to off and LTO=full fails.
#   - BWA_LTO_MODE inherited from the environment does not override LTO=; on
#     the command line (what lto-build passes its inner make) it does.
#   - With EXPECT_LTO set (CI's clang-19 + lld-19 rows set it to `full`), the
#     real toolchain must resolve auto to that mode and LTO=thin to thin; and
#     a clang whose lld link only succeeds with a directory that ld.bfd (not
#     lld) searches must still resolve auto (and LTO=full) to full, with that
#     directory at the end of LTO_LDFLAGS as -Wl,-L=<dir> (BWA_LTO_LIBDIRS), and
#     with no directory added when the first probe succeeds.
#   - On every host, native archive consumers exercise retries for ordinary
#     and sysroot-only directories, and rejection of whitespace/shell-special
#     paths without exporting LTO flags.
#   - A recursive make inherits BWA_LTO_MODE and BWA_LTO_LIBDIRS; BWA_LTO_MODE
#     on the command line drops an inherited BWA_LTO_LIBDIRS.
#
# Stub compilers stand in for gcc and for a clang without lld, so only the
# EXPECT_LTO and clang-CXX/gcc-CC checks need a real clang. Archive checks use
# the real native compilers through wrappers that strip the ELF-only LTO flags.
#
# Inputs:
#   CXX, CC     — the real toolchain (default clang++ / clang)
#   EXPECT_LTO  — optional: the mode the real toolchain must resolve auto to
#   MAKE        — make binary (default: make)
#   AR, NM      — archive tools (defaults: ar / nm; must match the LTO compiler)

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."

MAKE="${MAKE:-make}"
CXX="${CXX:-clang++}"
CC="${CC:-clang}"
OS="$(uname -s)"

stub_dir="$(mktemp -d -t bwamem3-lto-XXXXXX)"
trap 'rm -rf "$stub_dir"' EXIT

# gcc stubs: a GCC --version banner, a floor-passing -dumpversion, and every
# compile fails, so they can only ever resolve to off.
for tool in g++ gcc; do
    cat > "$stub_dir/$tool" << 'EOF'
#!/bin/sh
case "$1" in
    --version)
        echo "gcc (GCC) 15.1.0"
        echo "Copyright (C) 2025 Free Software Foundation, Inc."
        exit 0
        ;;
    -dumpversion)
        echo "15"
        exit 0
        ;;
esac
exit 1
EOF
    chmod +x "$stub_dir/$tool"
done

fail=0
out_file="$stub_dir/out"
err_file="$stub_dir/err"

# ARM_TUNE= (empty) on every call: on aarch64 Linux the Makefile otherwise
# probes -mtune with these stub compilers and prints its own NOTE to stdout,
# which is not what this script checks.
#
# resolve <make args...>: run `make -s print-BWA_LTO_MODE`; sets rc, and leaves
# stdout and stderr in $out_file / $err_file. The value is stdout's last line: on
# a tree whose htslib is not configured yet, the first call also runs htslib's
# configure, which prints to stdout ahead of it.
resolve() {
    rc=0
    "$MAKE" -s print-BWA_LTO_MODE ARM_TUNE= "$@" > "$out_file" 2> "$err_file" || rc=$?
}

# expect_mode <mode> <description> <make args...>
expect_mode() {
    local want=$1 what=$2
    shift 2
    resolve "$@"
    local got
    got="$(tail -n 1 "$out_file")"
    if ((rc != 0)) || [[ $got != "$want" ]]; then
        echo "FAIL: $what: expected '$want' on stdout and exit 0, got '$got' (exit $rc)" >&2
        cat "$err_file" >&2
        fail=1
    else
        echo "ok: $what -> $got"
    fi
}

# expect_error <stderr substring> <description> <make args...>
expect_error() {
    local want=$1 what=$2
    shift 2
    resolve "$@"
    if ((rc == 0)) || ! grep -qF -- "$want" "$err_file"; then
        echo "FAIL: $what: expected a failure mentioning '$want', got exit $rc" >&2
        cat "$out_file" "$err_file" >&2
        fail=1
    else
        echo "ok: $what fails ($want)"
    fi
}

# expect_note <description>: after expect_mode, the NOTE must be on stderr only.
expect_note() {
    if ! grep -qF "NOTE: building without LTO" "$err_file"; then
        echo "FAIL: $1: no 'building without LTO' NOTE on stderr" >&2
        fail=1
    elif grep -qF "NOTE" "$out_file"; then
        echo "FAIL: $1: the NOTE leaked onto stdout, where print-<VAR> output is parsed" >&2
        fail=1
    else
        echo "ok: $1 prints its NOTE on stderr only"
    fi
}

# resolve_vars <description> <extra PATH dir or ""> <make args...>: one parse
# printing BWA_LTO_MODE, BWA_LTO_LIBDIRS and LTO_LDFLAGS; sets rc, got_mode,
# got_libdirs and got_ldflags (the last three lines of stdout).
resolve_vars() {
    local what=$1 extra_path=$2
    shift 2
    local path=$PATH
    [[ -n $extra_path ]] && path="$extra_path:$PATH"
    rc=0
    PATH=$path "$MAKE" -s print-BWA_LTO_MODE print-BWA_LTO_LIBDIRS print-LTO_LDFLAGS ARM_TUNE= "$@" \
        > "$out_file" 2> "$err_file" || rc=$?
    got_mode="$(tail -n 3 "$out_file" | sed -n 1p)"
    got_libdirs="$(tail -n 3 "$out_file" | sed -n 2p)"
    got_ldflags="$(tail -n 3 "$out_file" | sed -n 3p)"
    : "$what"
}

bfd_retry_checked=0

gcc_args=(CXX="$stub_dir/g++" CC="$stub_dir/gcc" ALLOW_UNSUPPORTED_COMPILER=1)

expect_error "LTO must be one of" "LTO=bogus" "${gcc_args[@]}" LTO=bogus
expect_error "LTO must be one of" "LTO='auto off'" "${gcc_args[@]}" "LTO=auto off"
rc=0
LTO=1 "$MAKE" -s print-BWA_LTO_MODE ARM_TUNE= "${gcc_args[@]}" > "$out_file" 2> "$err_file" || rc=$?
if ((rc != 0)) || ! grep -qF "from the environment is not auto" "$err_file" || [[ "$(tail -n 1 "$out_file")" != off ]]; then
    echo "FAIL: an invalid LTO in the environment should warn and resolve, not fail" >&2
    cat "$out_file" "$err_file" >&2
    fail=1
else
    echo "ok: LTO=1 in the environment warns and is ignored"
fi
expect_mode off "gcc toolchain, LTO=auto" "${gcc_args[@]}"
# macOS is refused before the compiler is looked at.
if [[ $OS == Darwin ]]; then gcc_full_err="not supported on macOS"; else gcc_full_err="needs clang for both CXX and CC"; fi
expect_error "$gcc_full_err" "gcc toolchain, LTO=full" "${gcc_args[@]}" LTO=full
expect_mode off "LTO=off" CXX="$CXX" CC="$CC" LTO=off
expect_mode off "ASAN=1" CXX="$CXX" CC="$CC" ASAN=1
# lto-build's inner make: BWA_LTO_MODE=off on the command line skips resolution,
# so even a forced LTO=full is not applied (or checked).
expect_mode off "BWA_LTO_MODE=off on the command line wins over LTO=full" \
    CXX="$CXX" CC="$CC" LTO=full BWA_LTO_MODE=off
rc=0
BWA_LTO_MODE=full "$MAKE" -s print-BWA_LTO_MODE ARM_TUNE= "${gcc_args[@]}" > "$out_file" 2> "$err_file" || rc=$?
if ((rc != 0)) || [[ "$(tail -n 1 "$out_file")" != off ]]; then
    echo "FAIL: top-level resolution with inherited BWA_LTO_MODE=full: expected off and exit 0, got '$(tail -n 1 "$out_file")' (exit $rc)" >&2
    cat "$err_file" >&2
    fail=1
else
    echo "ok: BWA_LTO_MODE in the environment is ignored at the top level"
fi
# A recursive make (MAKELEVEL > 0) takes the mode and the ld.bfd directories
# from its parent's environment without probing again.
rc=0
MAKELEVEL=1 BWA_LTO_MODE=full BWA_LTO_LIBDIRS=-Wl,-L=/inherited "$MAKE" -s print-LTO_LDFLAGS ARM_TUNE= "${gcc_args[@]}" \
    > "$out_file" 2> "$err_file" || rc=$?
if ((rc != 0)) || [[ "$(tail -n 1 "$out_file")" != "-flto=full -ffat-lto-objects -fuse-ld=lld -Wl,-L=/inherited" ]]; then
    echo "FAIL: sub-make inheritance: expected full LTO flags with the inherited directory and exit 0, got '$(tail -n 1 "$out_file")' (exit $rc)" >&2
    cat "$err_file" >&2
    fail=1
else
    echo "ok: a sub-make inherits the mode and the ld.bfd directories"
fi
# BWA_LTO_MODE on the command line (lto-build's inner make) drops an inherited
# BWA_LTO_LIBDIRS.
rc=0
BWA_LTO_LIBDIRS=-Wl,-L=/stale "$MAKE" -s print-BWA_LTO_LIBDIRS ARM_TUNE= "${gcc_args[@]}" BWA_LTO_MODE=off \
    > "$out_file" 2> "$err_file" || rc=$?
if ((rc != 0)); then
    echo "FAIL: make failed in the BWA_LTO_MODE command-line override check (exit $rc)" >&2
    cat "$err_file" >&2
    fail=1
elif [[ -n "$(tail -n 1 "$out_file")" ]]; then
    echo "FAIL: BWA_LTO_MODE on the command line kept an inherited BWA_LTO_LIBDIRS ('$(tail -n 1 "$out_file")')" >&2
    fail=1
else
    echo "ok: BWA_LTO_MODE on the command line drops an inherited BWA_LTO_LIBDIRS"
fi

if [[ $OS == Darwin ]]; then
    expect_mode off "macOS, LTO=auto" CXX="$CXX" CC="$CC"
    expect_error "not supported on macOS" "macOS, LTO=full" CXX="$CXX" CC="$CC" LTO=full
elif "$CXX" --version 2> /dev/null | head -1 | grep -q clang; then
    # A clang whose -fuse-ld=lld link fails: wrap the real one.
    for pair in "clang++:$CXX" "clang:$CC"; do
        name=${pair%%:*}
        real=${pair#*:}
        cat > "$stub_dir/nolld-$name" << EOF
#!/bin/sh
for a in "\$@"; do
    case "\$a" in -fuse-ld=lld) echo "clang: error: invalid linker name in argument '-fuse-ld=lld'" >&2; exit 1 ;; esac
done
exec "$real" "\$@"
EOF
        chmod +x "$stub_dir/nolld-$name"
    done
    expect_mode off "clang without lld, LTO=auto" CXX="$stub_dir/nolld-clang++" CC="$stub_dir/nolld-clang"
    expect_note "clang without lld, LTO=auto"
    expect_error "needs clang with -ffat-lto-objects and ld.lld" "clang without lld, LTO=full" \
        CXX="$stub_dir/nolld-clang++" CC="$stub_dir/nolld-clang" LTO=full
    expect_mode off "clang CXX with a gcc CC, LTO=auto" CXX="$CXX" CC="$stub_dir/gcc"
    expect_note "clang CXX with a gcc CC, LTO=auto"
    if [[ -n ${EXPECT_LTO:-} ]]; then
        expect_mode "$EXPECT_LTO" "real toolchain, LTO=auto" CXX="$CXX" CC="$CC"
        expect_mode thin "real toolchain, LTO=thin" CXX="$CXX" CC="$CC" LTO=thin
        # Test the installed archive tools with actual optimized fat objects,
        # not the native-only wrappers below. Both languages must be extracted
        # and run correctly by a native consumer AND an LTO consumer.
        printf 'int archive_c(const int *p, int n){int s=0; for(int i=0;i<n;++i)s+=p[i]; return s;}\n' > "$stub_dir/member.c"
        printf 'int archive_cpp(const int *p, int n){int s=0; for(int i=0;i<n;++i)s+=p[i]*2; return s;}\n' > "$stub_dir/member.cpp"
        printf 'extern "C" int archive_c(const int*,int);\nint archive_cpp(const int*,int);\nint main(){int p[]={1,2,3};return archive_c(p,3)+archive_cpp(p,3)!=18;}\n' > "$stub_dir/consumer.cpp"
        "$CXX" -O3 -c "$stub_dir/consumer.cpp" -o "$stub_dir/consumer.o"
        for mode in full thin; do
            "$CC" -g -O3 -fopenmp -flto="$mode" -ffat-lto-objects -c "$stub_dir/member.c" -o "$stub_dir/member-c.o"
            "$CXX" -g -O3 -fopenmp -flto="$mode" -ffat-lto-objects -c "$stub_dir/member.cpp" -o "$stub_dir/member-cpp.o"
            "${AR:-ar}" rcs "$stub_dir/members.a" "$stub_dir/member-c.o" "$stub_dir/member-cpp.o"
            "${AR:-ar}" --version
            "${NM:-nm}" --print-armap "$stub_dir/members.a"
            "$CXX" "$stub_dir/consumer.o" "$stub_dir/members.a" -o "$stub_dir/native-consumer"
            "$stub_dir/native-consumer"
            "$CXX" -flto="$mode" -ffat-lto-objects -fuse-ld=lld "$stub_dir/consumer.o" "$stub_dir/members.a" -o "$stub_dir/lto-consumer"
            "$stub_dir/lto-consumer"
            echo "ok: real $mode fat archive: native and LTO consumers link and run"
        done
        # The first probe succeeded, so no ld.bfd directories were added.
        resolve_vars "real toolchain: no ld.bfd directories added" "" CXX="$CXX" CC="$CC"
        if ((rc == 0)) && [[ -z $got_libdirs ]]; then
            echo "ok: real toolchain adds no ld.bfd directories"
        else
            echo "FAIL: real toolchain: expected an empty BWA_LTO_LIBDIRS, got '$got_libdirs' (exit $rc)" >&2
            fail=1
        fi
        # A library only ld.bfd finds (libdeflate under /usr/local/lib64 on
        # Amazon Linux 2023): a fake ld.bfd whose only existing SEARCH_DIR is
        # $bfd_dir, and a clang whose lld link fails unless -Wl,-L=$bfd_dir is
        # given. The probe must retry with ld.bfd's directories, keep LTO on,
        # and put them on the LTO link as -Wl,-L (searched after every other
        # -L, so they cannot shadow the vendored libraries), keeping the `=`.
        bfd_dir="$stub_dir/bfd-only-lib"
        bfd_flag="-Wl,-L=$bfd_dir"
        mkdir -p "$bfd_dir" "$stub_dir/bin"
        cat > "$stub_dir/bin/ld.bfd" << EOF
#!/bin/sh
echo 'SEARCH_DIR("=$bfd_dir"); SEARCH_DIR("=$stub_dir/no-such-dir"); SEARCH_DIR("$stub_dir/has space");'
EOF
        chmod +x "$stub_dir/bin/ld.bfd"
        mkdir -p "$stub_dir/has space"
        for pair in "clang++:$CXX" "clang:$CC"; do
            name=${pair%%:*}
            real=${pair#*:}
            cat > "$stub_dir/bfdonly-$name" << EOF
#!/bin/sh
lld=0
found=0
for a in "\$@"; do
    case "\$a" in -fuse-ld=lld) lld=1 ;; "$bfd_flag") found=1 ;; esac
done
if [ \$lld = 1 ] && [ \$found = 0 ]; then
    echo "ld.lld: error: unable to find library -lstub" >&2
    exit 1
fi
exec "$real" "\$@"
EOF
            chmod +x "$stub_dir/bfdonly-$name"
        done
        bfd_args=(CXX="$stub_dir/bfdonly-clang++" CC="$stub_dir/bfdonly-clang")
        resolve_vars "a library only ld.bfd finds" "$stub_dir/bin" "${bfd_args[@]}"
        if ((rc == 0)) && [[ $got_mode == full && $got_libdirs == "$bfd_flag" && $got_ldflags == *" $bfd_flag" ]] \
            && ! grep -qF "NOTE" "$err_file"; then
            echo "ok: a library only ld.bfd finds -> full, with $bfd_flag at the end of LTO_LDFLAGS"
            bfd_retry_checked=1
        else
            echo "FAIL: a library only ld.bfd finds: expected 'full', BWA_LTO_LIBDIRS '$bfd_flag' and LTO_LDFLAGS ending in it, with no NOTE; got '$got_mode' / '$got_libdirs' / '$got_ldflags' (exit $rc)" >&2
            cat "$err_file" >&2
            fail=1
        fi
        resolve_vars "a library only ld.bfd finds, LTO=full" "$stub_dir/bin" "${bfd_args[@]}" LTO=full
        if ((rc == 0)) && [[ $got_mode == full && $got_libdirs == "$bfd_flag" ]]; then
            echo "ok: a library only ld.bfd finds, LTO=full -> full"
        else
            echo "FAIL: a library only ld.bfd finds, LTO=full: got '$got_mode' / '$got_libdirs' (exit $rc)" >&2
            cat "$err_file" >&2
            fail=1
        fi
    else
        echo "skip: EXPECT_LTO unset, not asserting the real toolchain's resolved mode"
    fi
else
    echo "skip: CXX=$CXX is not clang; the clang-only cases need a real clang"
fi

# Exercise the archive gate on every host with real native objects. Only the
# ELF-only LTO flags are stripped: these cases test probe control flow, not
# whether this host supports fat objects (EXPECT_LTO covers that on Linux).
real_cxx="$(command -v "$CXX")"
real_ar="$(command -v ar)"
cat > "$stub_dir/native-clang" << 'EOF'
#!/usr/bin/env bash
set -euo pipefail
case "${1:-}" in
    --version) echo 'clang version 21.0.0'; exit 0 ;;
    -dumpversion) echo 21; exit 0 ;;
esac
if [[ ${1:-} == -fuse-ld=lld && ${2:-} == -### ]]; then
    printf ' "ld.lld" "--sysroot=%s"\n' "${RETRY_SYSROOT:-}"
    exit 0
fi
if [[ ${CHECK_PRODUCTION_FLAGS:-0} == 1 ]]; then
    for arg in "$@"; do
        case "$arg" in
            */a.c) language=c ;;
            */b.cpp) language=cpp ;;
            *) continue ;;
        esac
        flags=" $* "
        for flag in -g -O3 "-DPROBE_${language}"; do
            if [[ $flags != *" $flag "* ]]; then
                echo "probe $language compile missing production flag $flag" >&2
                exit 1
            fi
        done
    done
fi
args=()
archive=0
lto=0
retry_found=0
for arg in "$@"; do
    case "$arg" in
        -flto=*) lto=1 ;;
        -ffat-lto-objects | -fuse-ld=lld) ;;
        -Wl,-L*)
            [[ $arg == "${RETRY_FLAG:-}" ]] && retry_found=1
            ;;
        *.a) archive=1; args+=("$arg") ;;
        *) args+=("$arg") ;;
    esac
done
if ((archive && lto)) && [[ -n ${RETRY_FLAG:-} ]] && ((retry_found == 0)); then
    exit 1
fi
if ((archive)) && [[ ${REJECT_ARCHIVE_LINK:-} == native && $lto == 0 || ${REJECT_ARCHIVE_LINK:-} == lto && $lto == 1 ]]; then
    exit 1
fi
for arg in "${args[@]}"; do
    case "$arg" in *.c) exec "$REAL_CC" "${args[@]}" ;; esac
done
exec "$REAL_CXX" "${args[@]}"
EOF
cat > "$stub_dir/archive" << 'EOF'
#!/usr/bin/env bash
set -euo pipefail
# Simulate an archiver that succeeds but omits one language's symbols.
args=()
for arg in "$@"; do
    case "${OMIT_MEMBER:-}:$arg" in
        c:*/a.o | cpp:*/b.o) ;;
        *) args+=("$arg") ;;
    esac
done
exec "$REAL_AR" "${args[@]}"
EOF
chmod +x "$stub_dir/native-clang" "$stub_dir/archive"
REAL_CC="$(command -v "$CC")"
export REAL_CXX="$real_cxx" REAL_CC REAL_AR="$real_ar"
probe_args=(CXX="$stub_dir/native-clang" CC="$stub_dir/native-clang"
    UNAME_S=Linux UNAME_M=x86_64 arch=avx2 BWA_LTO_PROBE_LIBS=
    ARCH_FLAGS= LIBSAIS_OPENMP_CFLAGS=
    AR="$stub_dir/archive")
expect_mode full "archive and both consumers work" "${probe_args[@]}"
CHECK_PRODUCTION_FLAGS=1 expect_mode full "probe compiles with production flags" "${probe_args[@]}" \
    'EXTRA_CFLAGS=-DPROBE_c' 'EXTRA_CXXFLAGS=-DPROBE_cpp' ARCH_FLAGS=
for member in c cpp; do
    OMIT_MEMBER="$member" expect_mode off "archive omits $member entry" "${probe_args[@]}"
    expect_note "archive omits $member entry"
    OMIT_MEMBER="$member" expect_error "archive" "forced LTO with missing $member entry" "${probe_args[@]}" LTO=full
done
for consumer in native lto; do
    REJECT_ARCHIVE_LINK="$consumer" expect_mode off "archive rejects $consumer consumer" "${probe_args[@]}"
    expect_note "archive rejects $consumer consumer"
    REJECT_ARCHIVE_LINK="$consumer" expect_error "archive" "forced LTO with rejected $consumer consumer" "${probe_args[@]}" LTO=thin
done

# Retry path checks run on every host using the native archive consumer above.
# ld.bfd output is data: unsafe entries must not become shell recipe text.
mkdir -p "$stub_dir/retry-bin"
cat > "$stub_dir/retry-bin/ld.bfd" << 'EOF'
#!/bin/sh
printf 'SEARCH_DIR("%s");\n' "$RETRY_DIR"
EOF
chmod +x "$stub_dir/retry-bin/ld.bfd"
for retry_dir in "$stub_dir/plain+lib" "=$stub_dir/equal-lib"; do
    mkdir -p "${retry_dir#=}"
    RETRY_DIR="$retry_dir" RETRY_FLAG="-Wl,-L$retry_dir" resolve_vars "safe retry directory" "$stub_dir/retry-bin" "${probe_args[@]}"
    if ((rc != 0)) || [[ $got_mode != full || $got_libdirs != "-Wl,-L$retry_dir" ]]; then
        echo "FAIL: safe retry directory '$retry_dir': got '$got_mode' / '$got_libdirs' (exit $rc)" >&2
        fail=1
    else
        echo 'ok: safe retry directory retains LTO'
    fi
done
sysroot="$stub_dir/sysroot"
mkdir -p "$sysroot/bwa-retry-only/lib"
RETRY_SYSROOT="$sysroot" RETRY_DIR='=/bwa-retry-only/lib' RETRY_FLAG='-Wl,-L=/bwa-retry-only/lib' \
    resolve_vars 'sysroot-only retry directory' "$stub_dir/retry-bin" "${probe_args[@]}"
if ((rc != 0)) || [[ $got_mode != full || $got_libdirs != '-Wl,-L=/bwa-retry-only/lib' ]]; then
    echo "FAIL: sysroot-only retry directory: got '$got_mode' / '$got_libdirs' (exit $rc)" >&2
    fail=1
else
    echo 'ok: sysroot-only retry directory retains LTO'
fi
for suffix in "has space" "O'Reilly" "dollar\$lib" 'semi;lib' 'comma,lib' 'back`tick' 'hash#lib'; do
    retry_dir="$stub_dir/$suffix"
    mkdir -p "$retry_dir"
    RETRY_DIR="$retry_dir" RETRY_FLAG="-Wl,-L$retry_dir" resolve_vars "unsafe retry directory" "$stub_dir/retry-bin" "${probe_args[@]}"
    if ((rc != 0)) || [[ $got_mode != off || -n $got_libdirs || -n $got_ldflags ]]; then
        echo "FAIL: unsafe retry directory '$retry_dir' did not fall back to off (exit $rc)" >&2
        fail=1
    else
        echo 'ok: unsafe retry directory falls back without LTO flags'
    fi
done

# The actual library recipe must honor AR, not just the capability probe.
cat > "$stub_dir/reject-ar" << 'EOF'
#!/bin/sh
echo 'selected archiver invoked' >&2
exit 19
EOF
chmod +x "$stub_dir/reject-ar"
printf 'int archive_override_entry(void){return 7;}\n' > "$stub_dir/override.cpp"
"$real_cxx" -c "$stub_dir/override.cpp" -o "$stub_dir/override.o"
rc=0
"$MAKE" -s "$stub_dir/override.a" LTO=off CXX="$CXX" CC="$CC" \
    BWA_LIB="$stub_dir/override.a" OBJS="$stub_dir/override.o" KERNEL_TIER_OBJS= \
    AR="$stub_dir/reject-ar" > "$out_file" 2> "$err_file" || rc=$?
if ((rc == 0)) || ! grep -qF 'selected archiver invoked' "$err_file"; then
    echo 'FAIL: libbwa archive recipe ignored AR override' >&2
    fail=1
else
    echo 'ok: libbwa archive recipe honors AR override'
fi

if ((fail != 0)); then
    echo "FAIL: make_lto_mode"
    exit 1
fi
if ((bfd_retry_checked)); then retry_note="ld.bfd search-dir retry, "; else retry_note=""; fi
echo "PASS: make_lto_mode (LTO= validation, auto fallbacks, ${retry_note}forced-mode errors, BWA_LTO_MODE/BWA_LTO_LIBDIRS inheritance)"
