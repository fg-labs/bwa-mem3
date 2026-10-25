#!/usr/bin/env bash
# .github/scripts/ci-depth.sh
#
# Decide how deep this CI run checks, and export the answer as CI_DEPTH.
#
#   CI_DEPTH=full   the merge queue (merge_group), any event other than
#                   pull_request, and a pull request whose diff touches a path
#                   below (or whose diff cannot be computed)
#   CI_DEPTH=light  every other pull request
#
# Why
# ---
# Every change lands through the merge queue, whose merge_group run tests the
# exact commit that lands, so that run keeps every check at full depth. Pull
# request runs repeat on every push to a PR branch and are where most of the
# compute goes, so there the deep safety nets for the SIMD kernels (the full
# rescue kernel sweeps, the AVX-512BW kswv unit tree, the MSan padded-lane
# gate) run at full depth only on the pull requests that can break them: those
# that touch the kernels, their harnesses and tests, the build, the compiler
# install or this workflow. A light PR run still runs every other check, and the
# kernel sweeps at their reduced job counts; nothing that lands is less tested.
#
# Fail safe: anything unexpected (not a two-parent merge commit, the base not
# fetchable, git failing) gives full, never light. The script always exits 0.
#
# How the diff is computed
# ------------------------
# On pull_request, actions/checkout checks out refs/pull/<N>/merge at depth 1:
# a merge commit GitHub made of the PR head onto the base branch, whose first
# parent is that base. The shallow clone hides the parents from rev-parse, but
# the commit object still names them, so the first parent is read from it,
# fetched at depth 1 (a public repository needs no credential), and diffed
# against HEAD. That is exactly the PR's change as it would merge, whatever the
# base branch has done since the PR branched. --no-renames lists both sides of
# a rename, so moving a kernel file out of a matched path still counts.
#
# Inputs (environment):
#   EVENT_NAME     github.event_name
# Outputs:
#   CI_DEPTH=<full|light> appended to $GITHUB_ENV, and depth=<full|light> to
#   $GITHUB_OUTPUT, each when that variable is set.

set -euo pipefail
# Fail safe under errexit too: a command that fails outside the explicit
# fallbacks below gives full depth rather than a failed job.
trap 'emit full "unexpected failure at line $LINENO"' ERR

# Paths whose change makes a PR run at full depth. Shell case patterns, where
# `*` also matches `/`. The kernels (every TU the Makefile compiles at each ISA
# tier, its KERNEL_SRCS: bandedSWA, kswv, ksw, sam_encode) and every header
# they include (the shared macro.h and utils.h too), the rescue planner, the
# harnesses and unit tests of all of them, the build files (flags, targets),
# the compiler install (codegen), this workflow and this script.
full_paths=(
    'src/rescue_*'
    'src/ksw*'
    'src/bandedSWA*'
    'src/sam_encode*'
    'src/bsw_*'
    'src/kernel_dispatch.h'
    'src/macro.h'
    'src/utils.h'
    'src/simd_*'
    'src/neon_*'
    'src/x86_soa_pack.h'
    'src/matesw_u8.h'
    'test/rescue_*'
    'test/kernel_padded_lane_uninit_test.cpp'
    'test/framework/kswv_runner.*'
    'test/framework/rescue_prune_test.h'
    'test/unit/test_ksw*'
    'test/unit/test_bandedswa_*'
    'test/unit/test_sam_encode*'
    'test/unit/test_rescue_*'
    'test/unit/test_simd_*'
    'Makefile'
    'test/Makefile'
    '.github/actions/install-clang-19/*'
    '.github/workflows/ci.yml'
    '.github/scripts/ci-depth.sh'
)

# A failed write leaves CI_DEPTH / depth unset, which every consumer reads as
# full, so it costs time, never depth; it is reported as a warning rather than
# failing the job.
emit() { # <depth> <reason>
    echo "CI depth: $1 ($2)"
    if [ -n "${GITHUB_ENV:-}" ] && ! echo "CI_DEPTH=$1" >> "$GITHUB_ENV"; then
        echo "::warning::could not write CI_DEPTH to GITHUB_ENV; the steps will run at full depth"
    fi
    if [ -n "${GITHUB_OUTPUT:-}" ] && ! echo "depth=$1" >> "$GITHUB_OUTPUT"; then
        echo "::warning::could not write the depth output; the jobs that read it will run at full depth"
    fi
    exit 0
}

event=${EVENT_NAME:-}
case "$event" in
    pull_request) ;;
    merge_group) emit full "merge queue: every check at full depth" ;;
    *) emit full "event '${event:-unset}' is not pull_request" ;;
esac

# `git cat-file -p` prints the commit's own parent lines even in a shallow
# clone, where rev-parse HEAD^1 fails.
parents=$(git cat-file -p HEAD 2> /dev/null | sed -n 's/^parent //p') \
    || emit full "could not read HEAD"
n=$(printf '%s\n' "$parents" | grep -c . || true)
[ "$n" = 2 ] || emit full "HEAD has $n parents, not a pull request merge commit"
base=$(printf '%s\n' "$parents" | head -n 1)

if ! git cat-file -e "$base^{commit}" 2> /dev/null; then
    git fetch --quiet --no-tags --no-recurse-submodules --depth=1 origin "$base" \
        || emit full "could not fetch the base commit $base"
fi
changed=$(git diff --name-only --no-renames "$base" HEAD) \
    || emit full "could not diff against the base commit $base"

hits=()
while IFS= read -r f; do
    [ -n "$f" ] || continue
    for p in "${full_paths[@]}"; do
        # shellcheck disable=SC2254 # $p is a pattern on purpose
        case "$f" in $p)
            hits+=("$f")
            break
            ;;
        esac
    done
done <<< "$changed"

nfiles=$(printf '%s\n' "$changed" | grep -c . || true)
if [ "${#hits[@]}" -gt 0 ]; then
    printf '  %s\n' "${hits[@]}"
    emit full "pull request touches ${#hits[@]} kernel, harness or build path(s), listed above"
fi
emit light "pull request changes $nfiles file(s), none on a kernel, harness or build path"
