#!/usr/bin/env bash
# test/regression/make_clean_no_compiler_floor.sh
#
# Regression: the compiler-version floor gates *builds*, not clean-only goals.
#
# The floor (Makefile "Compiler version floor" block) evaluates at parse time,
# so it runs for every goal make is asked to make -- including `clean`, which
# invokes no compiler and only `rm`s artifacts. Left ungated it FAILs
# `make clean` on GNU make's default CXX (`g++`, i.e. gcc-13 on Ubuntu 24.04),
# which is below the floor. That broke both contributor cleanups and any CI
# step that runs `make clean` before a build with a pinned CXX -- concretely the
# kswv ASan steps in .github/workflows/proto-neon-kswv.yml, which run
# `make clean` bare and then `make ... CXX="ccache clang++-19"`: the bare clean
# tripped the floor and aborted the step under `bash -e` before the pinned build
# ever ran.
#
# This asserts the goal-aware gate both ways, so neither half can silently
# regress: a below-floor compiler must NOT block the clean-only goals, and must
# STILL block a real build (otherwise the exemption would be a hole that lets an
# accidentally-slow build through -- the very thing the floor exists to stop).
#
# Checked with `make -n` (dry run) so nothing is compiled or removed: the floor
# is a parse-time $(error), so it fires -- or is skipped -- before any recipe
# would run. A stub `g++` reporting version 13 stands in for a below-floor
# toolchain, so the test needs no particular compiler installed.
#
# Inputs:
#   MAKE — make binary (default: make)

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."

MAKE="${MAKE:-make}"

# A stub that answers the two probes the floor uses: a gcc `--version` banner
# (the "Free Software Foundation" line is what the Makefile keys off to call it
# gcc) and `-dumpversion` = 13, one below the gcc floor of 15.
stub_dir="$(mktemp -d -t bwamem3-floor-XXXXXX)"
trap 'rm -rf "$stub_dir"' EXIT

stub_cxx="$stub_dir/g++"
cat > "$stub_cxx" << 'EOF'
#!/bin/sh
case "$1" in
    --version)
        echo "g++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0"
        echo "Copyright (C) 2023 Free Software Foundation, Inc."
        ;;
    -dumpversion) echo "13" ;;
esac
exit 0
EOF
chmod +x "$stub_cxx"

FLOOR_MARKER="unsupported compiler"
fail=0

# Run `make -n <goals...> CXX=<stub>` and report whether the floor tripped.
# stdout+stderr are captured together so the marker is caught wherever make
# writes it; the recipe text `make -n` prints is irrelevant here.
floor_tripped() {
    local out rc=0
    out="$("$MAKE" -n CXX="$stub_cxx" "$@" 2>&1)" || rc=$?
    # Herestring, not `printf | grep`: under `set -o pipefail` a `grep -q` that
    # matches early closes the pipe, printf dies of SIGPIPE, and the pipeline
    # reports 141 -- flipping a tripped floor to "not tripped". Same trap the
    # sibling lints avoid the same way.
    if grep -qF "$FLOOR_MARKER" <<< "$out"; then
        return 0 # floor tripped
    fi
    # A build goal that failed for some other reason must not read as "not
    # tripped" -- that would let a broken invocation pass the enforce checks.
    if ((rc != 0)); then
        printf '%s\n' "$out" >&2
        echo "FAIL: 'make -n CXX=<stub> $*' exited $rc without the floor message -- unexpected build error" >&2
        fail=1
    fi
    return 1 # floor did not trip
}

# Clean-only goals must be exempt: no floor, and a clean dry run succeeds.
for goal in clean pgo-clean profile-clean lto-clean; do
    if floor_tripped "$goal"; then
        echo "FAIL: 'make $goal' tripped the compiler floor on a below-floor CXX (clean invokes no compiler)" >&2
        fail=1
    else
        echo "exempt: 'make $goal' is not blocked by the floor"
    fi
done

# A real build must STILL be gated: the bare default goal, an explicit build
# target, and a mixed goal list that contains a build target. Each aborts at the
# floor's parse-time $(error) -- which sits near the top of the Makefile, before
# the `-include`d generated makefiles -- so no recipe (and no submodule build)
# ever runs here.
for goal in "" bwa-mem3 "clean bwa-mem3"; do
    # goal is intentionally word-split: it carries zero or more goals.
    # shellcheck disable=SC2086
    if floor_tripped $goal; then
        echo "enforced: 'make ${goal:-<default goal>}' is blocked by the floor"
    else
        echo "FAIL: 'make ${goal:-<default goal>}' was NOT blocked by the floor on a below-floor CXX (a build must be gated)" >&2
        fail=1
    fi
done

if ((fail != 0)); then
    echo "FAIL: the compiler floor does not correctly distinguish clean-only goals from builds" >&2
    exit 1
fi
echo "PASS: clean-only goals bypass the compiler floor; builds remain gated"
