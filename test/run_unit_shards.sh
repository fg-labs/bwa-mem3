#!/usr/bin/env bash
# test/run_unit_shards.sh
#
# Run a doctest binary (in CI: test/bwa_mem3_tests_unit) as concurrent shards,
# so its wall time drops by up to the core count.
#
# Why: the unit binary runs its ~500 cases one after another on one core, and
# no one case dominates, so nothing short of running cases concurrently
# shortens it. That is ~16 s on every matrix row and over two minutes in the
# -O0 coverage build, while the other cores sit idle. The CPU spent is the
# same either way; the wall time (which is what a runner is billed for) is not.
#
# How: doctest numbers the cases that pass the filters 1..N in a fixed order
# (by file, then line) and runs only the cases whose number falls in
# --first..--last (inclusive, 1-based). Every shard is the same binary with the
# same arguments, so every shard numbers the cases identically; this script
# asks for N once (--count with the same filters) and hands each shard one
# contiguous slice of 1..N. The slices partition 1..N, so every case runs
# exactly once.
#
# There are more shards than CPUs (four per job slot by default), handed out to
# a fixed number of concurrent jobs as each finishes. The cost is not spread
# evenly over the file order: the banded-SW and kswv files hold most of it, so
# one slice per CPU leaves the whole run waiting on the one or two slices that
# landed on them. Smaller slices let the other jobs keep working meanwhile; a
# process start costs a few milliseconds, so they are cheap.
#
# Each shard is a separate process, so the unit tests must not share fixed
# file paths (the ones that write files use mkstemp(3) or a pid-keyed name),
# and a case cannot rely on an earlier case having run in the same process.
# Under gcov (the coverage job), concurrent processes writing the same .gcda
# files is safe: libgcov locks each file and merges its counters into it.
#
# Usage:
#   test/run_unit_shards.sh BINARY [doctest-args...]
#
#   doctest-args (e.g. -ts=unit/kswv) are passed to the count and to every
#   shard. The options that would break the partition or redirect a shard's
#   output (--first/--last/--count/--list-*/--order-by/--rand-seed/--out) are
#   rejected.
#
# Inputs (env):
#   UNIT_JOBS          shards run concurrently (default: the CPU count).
#   UNIT_SHARDS        number of slices (default: 4 x UNIT_JOBS). Both are
#                      capped at the number of selected test cases.
#   UNIT_JUNIT_PREFIX  if set, each shard runs with the JUnit reporter and
#                      writes "${UNIT_JUNIT_PREFIX}-shard<k>.xml". Unset, each
#                      shard uses doctest's console reporter.
#
# Output: one PASS:/FAIL: line per shard; a failing shard's full output
# (doctest's failure report on the console reporter, its JUnit XML otherwise)
# is printed after all shards finish. Exits non-zero when any shard exits
# non-zero (a failed assertion or a crash) or did not report an exit status,
# when a shard's JUnit file is missing, or when the filters select no test
# cases at all.
set -euo pipefail

fail() {
    echo "FAIL: $*" >&2
    exit 1
}

# Worker mode: one shard, invoked by the xargs pool below as
#   run_unit_shards.sh <k> <first> <last> [doctest-args...]
# with the binary and the scratch directory in the environment. It records the
# shard's exit status in a file and itself exits 0, so xargs keeps going and the
# parent sees every shard's status, not just the first failure.
if [[ -n ${RUN_UNIT_SHARDS_WORKER_DIR:-} ]]; then
    k="$1" first="$2" last="$3"
    shift 3
    shard_args=(--first="$first" --last="$last")
    if [[ -n ${UNIT_JUNIT_PREFIX:-} ]]; then
        shard_args+=(--reporters=junit --out="${UNIT_JUNIT_PREFIX}-shard${k}.xml")
    fi
    status=0
    "$RUN_UNIT_SHARDS_BIN" "${shard_args[@]}" "$@" \
        > "$RUN_UNIT_SHARDS_WORKER_DIR/shard${k}.log" 2>&1 || status=$?
    echo "$status" > "$RUN_UNIT_SHARDS_WORKER_DIR/shard${k}.status"
    exit 0
fi

(($# >= 1)) || fail "usage: ${BASH_SOURCE[0]##*/} BINARY [doctest-args...]"
bin="$1"
shift
[[ -x $bin ]] || fail "test binary '$bin' is not an executable file"
[[ $bin == */* ]] || bin="./$bin"

for arg in "$@"; do
    # doctest accepts -x, --x and the --dt-x prefixed spelling of each option.
    opt="${arg#-}"
    opt="${opt#-}"
    opt="${opt#dt-}"
    opt="${opt%%=*}"
    case "$opt" in
        first | f | last | l | count | c | list-test-cases | ltc | list-test-suites | lts | \
            order-by | ob | rand-seed | rs | out | o)
            fail "'$arg' conflicts with sharding; the shard runner sets the range, order and outputs itself"
            ;;
    esac
done

positive_int() {
    [[ $2 =~ ^[1-9][0-9]*$ ]] || fail "$1='$2' is not a positive integer"
}
if [[ -n ${UNIT_JOBS:-} ]]; then
    positive_int UNIT_JOBS "$UNIT_JOBS"
    jobs="$UNIT_JOBS"
else
    jobs="$(nproc 2> /dev/null || sysctl -n hw.ncpu 2> /dev/null || echo 1)"
fi
if [[ -n ${UNIT_SHARDS:-} ]]; then
    positive_int UNIT_SHARDS "$UNIT_SHARDS"
    shards="$UNIT_SHARDS"
else
    shards=$((4 * jobs))
fi

# The count runs with exactly the shards' filters, so N is the same N every
# shard numbers against. A count that cannot be parsed is a failure, not 0.
count_out="$("$bin" --count --no-colors "$@" 2> /dev/null)" \
    || fail "'$bin --count $*' exited non-zero"
total="$(printf '%s\n' "$count_out" \
    | sed -n 's/^\[doctest\] unskipped test cases passing the current filters: \([0-9][0-9]*\)$/\1/p')"
[[ -n $total ]] || fail "could not read the test-case count from '$bin --count $*'"
((total > 0)) || fail "the filters ($*) select no test cases -- nothing would run"
((shards <= total)) || shards="$total"
((jobs <= shards)) || jobs="$shards"

workdir="$(mktemp -d "${TMPDIR:-/tmp}/run_unit_shards.XXXXXX")"
trap 'rm -rf "$workdir"' EXIT

echo "Running $total test cases from $bin as $shards shard(s), $jobs at a time"
start=$SECONDS
firsts=()
lasts=()
for ((k = 0; k < shards; k++)); do
    firsts+=($((k * total / shards + 1)))
    lasts+=($(((k + 1) * total / shards)))
    if [[ -n ${UNIT_JUNIT_PREFIX:-} ]]; then
        rm -f "${UNIT_JUNIT_PREFIX}-shard${k}.xml"
    fi
done

# xargs -P is the job pool: it starts the next shard as soon as one exits. Each
# shard's arguments are one NUL-separated record of a fixed size, so arguments
# containing spaces or glob characters reach the binary unchanged.
xargs_status=0
for ((k = 0; k < shards; k++)); do
    printf '%s\0' "$k" "${firsts[k]}" "${lasts[k]}" "$@"
done | RUN_UNIT_SHARDS_WORKER_DIR="$workdir" RUN_UNIT_SHARDS_BIN="$bin" \
    xargs -0 -n $((3 + $#)) -P "$jobs" bash "${BASH_SOURCE[0]}" || xargs_status=$?
elapsed=$((SECONDS - start))
((xargs_status == 0)) || echo "FAIL: the shard pool (xargs) exited $xargs_status" >&2

failed=0
statuses=()
for ((k = 0; k < shards; k++)); do
    label="shard $k/$shards (test cases ${firsts[k]}-${lasts[k]})"
    junit="${UNIT_JUNIT_PREFIX:+${UNIT_JUNIT_PREFIX}-shard${k}.xml}"
    status="$(cat "$workdir/shard${k}.status" 2> /dev/null || echo missing)"
    statuses+=("$status")
    if [[ $status == missing ]]; then
        echo "FAIL: $label never reported an exit status" >&2
    elif [[ $status != 0 ]]; then
        echo "FAIL: $label exited $status" >&2
    elif [[ -n $junit && ! -s $junit ]]; then
        echo "FAIL: $label exited 0 but wrote no JUnit file at $junit" >&2
        statuses[k]=no-junit
    else
        echo "PASS: $label"
        continue
    fi
    failed=$((failed + 1))
done

for ((k = 0; k < shards; k++)); do
    [[ ${statuses[k]} != 0 ]] || continue
    echo "---- output of shard $k (test cases ${firsts[k]}-${lasts[k]}) ----" >&2
    cat "$workdir/shard${k}.log" >&2 2> /dev/null || true
    junit="${UNIT_JUNIT_PREFIX:+${UNIT_JUNIT_PREFIX}-shard${k}.xml}"
    if [[ -n $junit && -s $junit ]]; then
        cat "$junit" >&2
    fi
    echo "---- end of shard $k; rerun it alone with: $bin --first=${firsts[k]} --last=${lasts[k]} $* ----" >&2
done

if ((failed > 0 || xargs_status != 0)); then
    fail "$failed of $shards shard(s) failed (${elapsed}s)"
fi
echo "PASS: all $total test cases ran across $shards shard(s) (${elapsed}s)"
