#!/usr/bin/env bash
# test/regression/make_lto_mode.sh
#
# Regression: how the Makefile resolves LTO= into the mode the build uses
# (BWA_LTO_MODE), read back with `make -s print-BWA_LTO_MODE`. That goal runs the
# same parse-time logic as a real build but compiles nothing.
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
#     real toolchain must resolve auto to that mode and LTO=thin to thin.
#
# Stub compilers stand in for gcc and for a clang without lld, so only the
# EXPECT_LTO and clang-CXX/gcc-CC checks need a real clang.
#
# Inputs:
#   CXX, CC     — the real toolchain (default clang++ / clang)
#   EXPECT_LTO  — optional: the mode the real toolchain must resolve auto to
#   MAKE        — make binary (default: make)

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

gcc_args=(CXX="$stub_dir/g++" CC="$stub_dir/gcc" ALLOW_UNSUPPORTED_COMPILER=1)

expect_error "LTO must be one of" "LTO=bogus" "${gcc_args[@]}" LTO=bogus
expect_error "LTO must be one of" "LTO='auto off'" "${gcc_args[@]}" "LTO=auto off"
LTO=1 "$MAKE" -s print-BWA_LTO_MODE ARM_TUNE= "${gcc_args[@]}" > "$out_file" 2> "$err_file" || true
if ! grep -qF "from the environment is not auto" "$err_file" || [[ "$(tail -n 1 "$out_file")" != off ]]; then
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
BWA_LTO_MODE=full "$MAKE" -s print-BWA_LTO_MODE ARM_TUNE= "${gcc_args[@]}" > "$out_file" 2> "$err_file" || true
if [[ "$(tail -n 1 "$out_file")" != off ]]; then
    echo "FAIL: BWA_LTO_MODE=full in the environment overrode the top-level resolution (got '$(tail -n 1 "$out_file")')" >&2
    fail=1
else
    echo "ok: BWA_LTO_MODE in the environment is ignored at the top level"
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
    else
        echo "skip: EXPECT_LTO unset, not asserting the real toolchain's resolved mode"
    fi
else
    echo "skip: CXX=$CXX is not clang; the clang-only cases need a real clang"
fi

if ((fail != 0)); then
    echo "FAIL: make_lto_mode"
    exit 1
fi
echo "PASS: make_lto_mode (LTO= validation, auto fallbacks, forced-mode errors, BWA_LTO_MODE inheritance)"
