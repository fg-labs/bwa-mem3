#!/bin/sh
# Compiler shim that routes one compiler through ccache. The setup-ccache action
# symlinks each compiler name the jobs use (clang++-19, clang-19, g++, gcc, cc,
# ...) to this script in a directory it puts first on PATH, so CXX/CC keep their
# plain values (`clang++-19`, not `ccache clang++-19`) and nothing that inspects
# them -- the Makefile's CC pairing, compiler floor, LTO probe and .build-flags
# stamp, htslib's configure, the vendored cmake builds -- sees a difference.
#
# The real compiler is the first one of the same name on PATH after this
# directory. It is passed to ccache by absolute path, so ccache never has to
# search for it itself.
#
# -march=native / -mcpu=native / -mtune=native: ccache keys on the command line,
# not on the CPU the compiler resolves "native" to, and hosted runners vary
# (with and without AVX-512, different CPU models). An object built native on
# one runner and restored on another could carry instructions the second host
# lacks. For those compiles only, a description of what "native" means on this
# host -- the compiler's own expansion of the flag plus the CPU's feature list
# -- is hashed in through CCACHE_EXTRAFILES, so native objects are shared only
# between hosts that resolve "native" identically. Every other compile stays
# host-independent and shared across runners.
#
# POSIX sh and no subshells on the common path: it runs once per compile.

name=${0##*/}
self_dir=${0%/*}

real=
saved_ifs=$IFS
IFS=:
set -f
for dir in $PATH; do
    [ "$dir" = "$self_dir" ] && continue
    [ -n "$dir" ] || dir=.
    if [ -f "$dir/$name" ] && [ -x "$dir/$name" ]; then
        real=$dir/$name
        break
    fi
done
set +f
IFS=$saved_ifs

if [ -z "$real" ]; then
    echo "compiler-shim: no '$name' on PATH outside $self_dir" >&2
    exit 127
fi

for arg in "$@"; do
    case $arg in
        -march=native | -mcpu=native | -mtune=native)
            isa_file=$self_dir/native-isa.$name.${arg%%=*}
            if [ ! -s "$isa_file" ]; then
                tmp=$isa_file.$$
                {
                    echo "$name $arg"
                    # -E: preprocess only, so the -### plan names no temp files.
                    # From / so the plan's compilation-dir entries do not
                    # depend on which tree compiled the first native object.
                    (cd / && "$real" "$arg" -### -E -x c /dev/null -o /dev/null 2>&1)
                    if [ -r /proc/cpuinfo ]; then
                        grep -m1 -E '^(flags|Features)[[:space:]]*:' /proc/cpuinfo
                    else
                        sysctl -n machdep.cpu.brand_string hw.cpufamily 2> /dev/null
                    fi
                } > "$tmp" 2>&1
                mv -f "$tmp" "$isa_file"
            fi
            CCACHE_EXTRAFILES=$isa_file
            export CCACHE_EXTRAFILES
            break
            ;;
    esac
done

exec ccache "$real" "$@"
