#!/usr/bin/env bash
# test/regression/make_test_relink.sh
#
# Regression: every standalone binary test/Makefile builds is relinked when
# ../libbwa.a changes.
#
# Was: the standalone rules depended only on their own object, so `make <test>`
# after a src/ change said "up to date" and the test ran against the OLD
# library. A parity test then passes or fails against code that is no longer
# in the tree, and nothing says so.
#
# The check is a dry run: with the tree up to date, `make -n <bin>` must print
# no link for the binary; after ../libbwa.a is touched, it must print one for
# every binary. Nothing is compiled or linked. The library's mtime is put back
# afterwards so the next real build sees no spurious change.
#
# Inputs:
#   MAKE — make binary (default: make). No MAKE_ARGS: test/Makefile is driven
#          directly and only under -n, so the flags the tree was built with do
#          not matter here.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."

MAKE="${MAKE:-make}"
LIB=libbwa.a

if [[ ! -f "$LIB" ]]; then
    echo "FAIL: $LIB is not built — run this against a built tree"
    exit 1
fi

# Whole-second mtime, for the same reason as make_header_deps.sh: GNU Make 3.81
# (Apple's stock make) compares whole seconds, so the touched library has to
# land a full second past the binaries or make sees nothing to do.
if stat -c %Y . > /dev/null 2>&1; then
    mtime_sec() { stat -c %Y "$1"; } # GNU coreutils
else
    mtime_sec() { stat -f %m "$1"; } # BSD / macOS
fi

# Enumerated from the Makefile, not from a copy of its list, so a binary added
# to EXE without the ../libbwa.a prerequisite fails here.
read -r -a all_bins <<< "$("$MAKE" -s -C test print-EXE)"
if ((${#all_bins[@]} == 0)); then
    echo "FAIL: test/Makefile EXE is empty — is print-% missing?"
    exit 1
fi

# Only binaries that exist can be up to date; an unbuilt one would show a link
# line either way and prove nothing.
bins=()
for bin in "${all_bins[@]}"; do
    [[ -x "test/$bin" ]] && bins+=("$bin")
done
if ((${#bins[@]} == 0)); then
    echo "FAIL: none of the standalone binaries (${all_bins[*]}) is built — run 'make -C test' first"
    exit 1
fi

# would_link <bin>: true when a dry run of test/Makefile would link <bin>.
would_link() {
    "$MAKE" -n -C test "$1" 2> /dev/null | grep -qE -- "(^| )-o $1( |\$)"
}

fail=0

# --- Precondition: the tree is up to date, so a link line below can only come
# --- from the library prerequisite and not from a stale object.
for bin in "${bins[@]}"; do
    if would_link "$bin"; then
        echo "FAIL: test/$bin is already out of date before the library is touched — build test/ first"
        fail=1
    fi
done
((fail == 0)) || exit 1

MTIME_REF="$(mktemp)"
touch -r "$LIB" "$MTIME_REF"
trap 'touch -r "$MTIME_REF" "$LIB" 2> /dev/null || true; rm -f "$MTIME_REF"' EXIT

newest=0
for bin in "${bins[@]}"; do
    m="$(mtime_sec "test/$bin")"
    ((m > newest)) && newest=$m
done
for _ in 1 2 3; do
    touch "$LIB"
    (($(mtime_sec "$LIB") > newest)) && break
    sleep 1
done
if (($(mtime_sec "$LIB") <= newest)); then
    echo "FAIL: could not make $LIB newer than the test binaries (timestamp granularity?)"
    exit 1
fi

# --- The check: every standalone binary relinks after the library changed.
relinked=0
for bin in "${bins[@]}"; do
    if would_link "$bin"; then
        relinked=$((relinked + 1))
    else
        echo "FAIL: test/$bin would NOT relink after $LIB changed — its rule lacks the ../libbwa.a prerequisite"
        fail=1
    fi
done

if ((fail != 0)); then
    echo "FAIL: test/Makefile standalone binaries do not all track $LIB"
    exit 1
fi
echo "PASS: all $relinked built standalone test binaries relink after $LIB changes (${#all_bins[@]} listed in EXE)"
