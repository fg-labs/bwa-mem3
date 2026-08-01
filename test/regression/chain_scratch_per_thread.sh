#!/usr/bin/env bash
# test/regression/chain_scratch_per_thread.sh
#
# Regression: the chaining scratch buffers (worker_t::chain_scratch and
# worker_t::seed_scratch) must be sized PER THREAD, i.e. by
# nthreads * BATCH_SIZE, and must NOT scale with the batch size (-K) or with
# the read count of a chunk.
#
# Why this needs its own test. These buffers were originally sized by
# `nreads = chunk_bases / NREADS_ESTIMATE_AVG_BASES + 10` and indexed by
# seq_id, left over from when worker_bwt and worker_aln were two separate
# kt_for passes and the chains plus their seeds had to survive the barrier
# between them. Fusing them into worker_bwt_aln made that unnecessary, but the
# allocation was not shrunk, so seed_scratch reached 13.1 GB at -t 64 and
# 39.3 GB at -t 192 -- growth linear in -t that made a large -t untenable.
#
# The crucial point for testing: re-introducing the nreads sizing is
# COMPLETELY INVISIBLE to output. The alignments are byte-identical either
# way, so every parity/determinism test in this suite would still pass while
# the memory regression silently returned. Nothing else guards this, which is
# why the assertion is on the allocation report rather than on the SAM.
#
# Method: worker_alloc reports the scratch allocation separately from the
# per-read `regs` array:
#   "   per-thread chaining scratch: <N> MB (<T> threads)"
# Because that figure depends only on nthreads and BATCH_SIZE, three things must
# hold, and the test checks all three:
#
#   1. It is EXACTLY equal across a 16x change in -K. Asserting equality rather
#      than a ratio matters: the combined "Memory pre-allocation for Chaining"
#      total also contains `regs`, which legitimately tracks -K, so any ratio
#      test against the total needs a fudge factor whose safe range shifts with
#      the chosen -K values. Exact equality on the isolated figure has none of
#      that fragility.
#   2. It is EXACTLY linear in -t. Without this, a sizing that mixed in a
#      constant, or scaled with something other than the thread count, would
#      still pass (1).
#   3. The scratch is never resized per chunk. This one cannot be checked from
#      the figures at all: worker_alloc prints them before the first chunk is
#      parsed, so a reallocation in kt_pipeline step 1 -- which is exactly where
#      the old code sized chain_ar/seedBuf from ret->n_seqs -- would leave both
#      figures untouched and this test green. It is instead enforced in the code
#      (a hard size check in that step, surviving NDEBUG) and observed here via
#      the exit status plus the absence of its error message.
#
# Inputs (env vars):
#   BWA_MEM3  -- path to the bwa-mem3 binary under test
#   CHR22_FA  -- path to chr22.fa, pre-indexed with bwa-mem3 by the caller
#   TMPDIR    -- optional; scratch location (must be disk-backed, not tmpfs)
set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${CHR22_FA:?CHR22_FA must be set}"

WORK=$(mktemp -d "${TMPDIR:-/tmp}/chain_scratch.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

# A tiny read set: the assertion is about ALLOCATION, not alignment, so the
# reads only need to be real enough to get through worker_alloc. Generated
# here rather than committed, per the project's no-test-data-files rule.
python3 - "$WORK" << 'PY'
import gzip, random, sys
out = sys.argv[1]
random.seed(20260725)          # fixed seed: fixture must be reproducible
bases = "ACGT"
with gzip.open(f"{out}/r1.fq.gz", "wt") as f1, gzip.open(f"{out}/r2.fq.gz", "wt") as f2:
    for i in range(2000):
        s1 = "".join(random.choice(bases) for _ in range(100))
        s2 = "".join(random.choice(bases) for _ in range(100))
        f1.write(f"@read{i}/1\n{s1}\n+\n{'I'*100}\n")
        f2.write(f"@read{i}/2\n{s2}\n+\n{'I'*100}\n")
PY

# Report the isolated per-thread scratch allocation, in MB, for a given -t / -K.
scratch_mb_for() {
    local t="$1" k="$2" err="$WORK/t$1k$2.err"
    local rc=0 mb

    # The exit status is asserted, not assumed, for two reasons. worker_alloc
    # reports the scratch line early, so a run that printed it and then died
    # would otherwise yield a perfectly comparable number from a failed run. And
    # the per-chunk resize guard in kt_pipeline step 1 reports a violation by
    # exiting non-zero -- see the note below on why that guard exists.
    "$BWA_MEM3" mem -t "$t" -K "$k" "$CHR22_FA" \
        "$WORK/r1.fq.gz" "$WORK/r2.fq.gz" > /dev/null 2> "$err" || rc=$?
    if [ "$rc" -ne 0 ]; then
        echo "FAIL: 'bwa-mem3 mem -t $t -K $k' exited $rc" >&2
        tail -20 "$err" >&2
        exit 1
    fi

    # grep exits 1 when the line is absent, which under `set -e` (plus pipefail)
    # would abort the script on the assignment itself and never reach the
    # diagnostic below -- so the failure is caught explicitly instead.
    mb="$(grep -m1 "per-thread chaining scratch" "$err" | sed 's/.*scratch: *//; s/ MB.*//')" || mb=""
    if [ -z "$mb" ]; then
        echo "FAIL: could not find the per-thread scratch line in $err" >&2
        echo "      (worker_alloc must report it; see src/fastmap.cpp)" >&2
        tail -20 "$err" >&2
        exit 1
    fi
    echo "$mb"
}

# --- 1. Invariant across -K -------------------------------------------------
SMALL=$(scratch_mb_for 4 4000000)  #   4 Mbase
LARGE=$(scratch_mb_for 4 64000000) #  64 Mbase, 16x

echo "per-thread chaining scratch:  -K 4M => ${SMALL} MB    -K 64M => ${LARGE} MB"

if [ "$SMALL" != "$LARGE" ]; then
    echo "FAIL: chaining scratch changed with -K (${SMALL} MB -> ${LARGE} MB)." >&2
    echo "      It must depend only on nthreads * BATCH_SIZE, never on the read" >&2
    echo "      count or batch size. See worker_alloc() in src/fastmap.cpp." >&2
    exit 1
fi

# --- 2. Exactly linear in -t ------------------------------------------------
# -K invariance alone does not pin the formula: a sizing that mixed in any
# constant, or that scaled with something other than the thread count, could
# still be -K-invariant. Doubling -t must double the figure exactly.
T8=$(scratch_mb_for 8 4000000)

echo "per-thread chaining scratch:  -t 4 => ${SMALL} MB    -t 8 => ${T8} MB"

awk -v a="$SMALL" -v b="$T8" 'BEGIN {
    # Both come from the same %0.4lf format, so an exact ratio is safe to
    # demand; the tolerance only absorbs that 4-decimal rounding.
    exit !(a > 0 && b > 0 && (b / a) > 1.9995 && (b / a) < 2.0005)
}' || {
    echo "FAIL: chaining scratch is not linear in -t (${SMALL} MB at -t 4," >&2
    echo "      ${T8} MB at -t 8; expected exactly 2x). The sizing must be" >&2
    echo "      nthreads * BATCH_SIZE. See worker_alloc() in src/fastmap.cpp." >&2
    exit 1
}

# --- 3. The per-chunk resize guard never fired ------------------------------
# kt_pipeline step 1 exits non-zero if the scratch was resized per chunk, which
# the exit-status check in scratch_mb_for already turns into a failure. Assert
# the message is absent as well, so this stays a failure even if that guard is
# ever softened to a warning: a reallocation at THAT site is invisible to the
# figures above, because worker_alloc reports them before the first chunk is
# read.
if grep -l "per-thread seed scratch was resized" "$WORK"/*.err > /dev/null 2>&1; then
    echo "FAIL: the per-chunk scratch resize guard fired:" >&2
    grep -h "per-thread seed scratch was resized" "$WORK"/*.err >&2
    exit 1
fi

echo "PASS: chaining scratch is per-thread ($SMALL MB at -t 4), invariant across a 16x -K change, exactly linear in -t, and never resized per chunk."
