#!/usr/bin/env bash
# test/regression/make_default_goal.sh
#
# Regression: a bare `make` builds the aligner, whatever else the root Makefile
# grows above `all:`.
#
# Was: the default goal was never stated, so make used its fallback — the first
# target of the first rule in the file. `myall` gets a rule only in the two
# no-`arch=` branches; once `arch=` is set on x86 (`make arch=avx2 CXX=g++
# USE_MIMALLOC=1`, exactly what CI's build step and `single:`'s own recursion
# run) no myall rule exists and the goal fell through to whichever rule came
# next. `all:` held that slot by ordering alone.
#
# Adding the objects-depend-on-$(FLAGS_STAMP) rule above it took the slot away:
# the default goal silently became `src/fastmap.o`, the first entry of $(OBJS).
# `make arch=avx2 …` then compiled one object and exited 0 — a *successful*
# build that produced no binary, so every x86 row failed one step later at
# `./bwa-mem3 version` with "No such file or directory" while arm64 (which does
# reach a myall rule, and whose rule sits above the new one) stayed green.
#
# Checked by name rather than by building: the failure is a mis-selected goal,
# and the goal make would pick is the whole finding. `make -n` cannot see it —
# on the built tree this runs against it reports nothing to do either way.
#
# Inputs:
#   MAKE      — make binary (default: make)
#   MAKE_ARGS — extra make arguments, e.g. "CXX=g++ USE_MIMALLOC=1". Must not
#               contain `arch=`; this test sets that itself, per case.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."

MAKE="${MAKE:-make}"
MAKE_ARGS="${MAKE_ARGS:-}"

if [[ "$MAKE_ARGS" == *arch=* ]]; then
    echo "FAIL: MAKE_ARGS must not set arch= — this test varies it per case"
    exit 1
fi

# `myall` is the no-arch goal on both platforms: it recurses into `single` on x86
# and `arm64` on Apple Silicon. Both end at a linked binary, which `all` also
# does; what must never appear here is an object file.
#
# The arch= cases use tiers real invocations pass, since those are the ones whose
# goal went wrong. Only ARCH_FLAGS varies between them, so the platform's own
# tier is enough to cover the branch.
case "$(uname -m)" in
    arm64 | aarch64) arch_cases=(arm64 native) ;;
    *)               arch_cases=(sse41 avx2 avx512bw native) ;;
esac

fail=0

default_goal_for() {
    # MAKE_ARGS is intentionally word-split: it carries zero or more arguments.
    # shellcheck disable=SC2086
    "$MAKE" $MAKE_ARGS "$@" -s print-.DEFAULT_GOAL
}

check() {
    local want="$1" got label
    shift
    label="${*:-<no arch>}"
    if ! got="$(default_goal_for "$@")"; then
        echo "FAIL: 'make $label print-.DEFAULT_GOAL' failed — a build error, not a goal finding"
        fail=1
        return
    fi
    if [[ "$got" != "$want" ]]; then
        echo "FAIL: default goal for 'make $label' is '$got', want '$want'"
        # The symptom is worth naming: an object here builds and exits 0.
        [[ "$got" == *.o ]] && echo "      (an object goal makes a bare build succeed while producing no binary)"
        fail=1
        return
    fi
    echo "goal: 'make $label' -> $got"
}

check myall
for arch in "${arch_cases[@]}"; do
    check all "arch=$arch"
done

if (( fail != 0 )); then
    echo "FAIL: root Makefile default goal does not build the aligner"
    exit 1
fi
echo "PASS: default goal builds the aligner for a bare make and every arch= tier"
