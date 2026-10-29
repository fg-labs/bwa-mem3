#!/usr/bin/env bash
# test/proc_freq_calibration_test.sh
#
# `proc_freq` (the __rdtsc() tick rate that every timing report divides by)
# is resolved at the first timing report by proc_freq_resolve()
# (src/profiling.cpp): read from CNTFRQ_EL0 on aarch64, measured over the run
# elsewhere. It used to be measured with a literal `sleep(1)` at startup,
# which charged a full second of wall time to EVERY invocation -- including
# `version`, `shm -l`, and usage errors, none of which ever read proc_freq.
# This test guards the replacement end to end; the branches of
# proc_freq_resolve() are pinned one by one in test/unit/test_proc_freq.cpp.
#
#   1. Startup is cheap: N back-to-back `version` calls finish in well
#      under N seconds. A reintroduced sleep(1) fails here immediately.
#   2. A real `mem` run reports a processor frequency inside a generous
#      band, which catches only a grossly mis-scaled rate. A modest
#      mis-scaling may pass: check 4 catches only a rate more than ~1.5x too
#      low, and a modestly high rate passes this script.
#   3. At -v 4, that `mem` run says where the rate came from: CNTFRQ_EL0 on
#      aarch64 (or the run window, where that register reads 0), the run
#      window elsewhere. A regression to sleeping for the
#      rate at report time (e.g. main() no longer opening the window) fails
#      here; checks 2 and 4 would still pass, as the sleep is accurate.
#   4. `index` resolves the rate before its "Total time taken" line divides
#      by it (an unresolved proc_freq of 0 prints `inf`), and the time it
#      prints does not exceed the wall time this script measures around it,
#      so a rate too low by more than the margin given in check 4 fails here
#      on every architecture.
#
# Usage: test/proc_freq_calibration_test.sh <bwa-mem3-binary> <fixtures-dir>

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 <bwa-mem3-binary> <fixtures-dir>" >&2
    exit 2
fi

abspath() { (cd "$(dirname "$1")" && printf '%s/%s\n' "$(pwd)" "$(basename "$1")"); }
bin="$(abspath "$1")"
fixtures="$(abspath "$2")"
ref="$fixtures/phix.fa"
reads="$fixtures/reads.fa"

[[ -x "$bin" ]] || {
    echo "FAIL: bwa-mem3 binary not executable at $bin" >&2
    exit 1
}
[[ -s "$ref" ]] || {
    echo "FAIL: phix.fa missing at $ref" >&2
    exit 1
}
[[ -s "$reads" ]] || {
    echo "FAIL: reads.fa missing at $reads" >&2
    exit 1
}

# Build the phiX FMI index if not already present.
if [[ ! -s "$ref.bwt.2bit.64" || ! -s "$ref.amb" ||
    ! -s "$ref.ann" || ! -s "$ref.pac" ]]; then
    "$bin" index "$ref" > /dev/null 2>&1 \
        || {
            echo "FAIL: bwa-mem3 index on phix.fa failed" >&2
            exit 1
        }
fi

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

# --- 1. Startup does not sleep ----------------------------------------------
# `version` reads no index and prints a handful of lines, so its entire
# runtime is process startup. Timed with bash's integer-second SECONDS
# rather than a high-resolution clock so the test needs no non-POSIX date
# or extra interpreter. 8 invocations against a 3 s budget sits between the
# two costs it has to separate: ~25x above the measured cost (~0.11 s for
# all 8, i.e. ~14 ms each) and ~2.7x below the old sleep(1) cost (8 s). The
# large upper margin is what keeps it from flaking on a loaded CI runner.
runs=8
budget=3
SECONDS=0
for _ in $(seq "$runs"); do
    "$bin" version > "$tmpdir/version.out" 2> "$tmpdir/version.err" \
        || {
            echo "FAIL: bwa-mem3 version exited non-zero" >&2
            cat "$tmpdir/version.err" >&2
            exit 1
        }
done
elapsed=$SECONDS
if ((elapsed > budget)); then
    echo "FAIL: $runs 'bwa-mem3 version' calls took ${elapsed}s (budget ${budget}s)." >&2
    echo "      Startup is sleeping again -- see proc_freq_mark_start() / proc_freq_resolve() in src/profiling.cpp." >&2
    exit 1
fi

# --- 2. The resolved rate is a plausible frequency --------------------------
# display_stats() (src/profiling.cpp) calls proc_freq_resolve(), which never
# leaves the rate at 0, and prints it as "Processor is running @<freq> MHz".
# Require a value inside a deliberately generous band that no real tick
# source (a ~24 MHz arm64 CNTVCT_EL0 through a multi-GHz x86 TSC) falls
# outside of. The band only catches gross errors; check 4 and the unit test
# pin the scale.
"$bin" mem "$ref" "$reads" > "$tmpdir/mem.sam" 2> "$tmpdir/mem.err" \
    || {
        echo "FAIL: bwa-mem3 mem on phix.fa exited non-zero" >&2
        cat "$tmpdir/mem.err" >&2
        exit 1
    }

# Capture the field verbatim ([^ ]+, not just digits) so whatever %lf
# actually printed -- an absurd magnitude or a non-numeric nan/inf -- lands
# in the failure message instead of the match silently failing and reporting
# a missing line.
freq_mhz="$(sed -nE 's/^Processor is running @([^ ]+) MHz$/\1/p' "$tmpdir/mem.err" | head -n1)"
[[ -n "$freq_mhz" ]] \
    || {
        echo "FAIL: 'Processor is running @... MHz' line missing from mem stderr" >&2
        cat "$tmpdir/mem.err" >&2
        exit 1
    }

# The regex rejects nan/inf before awk sees them: macOS awk evaluates a
# comparison with nan as true, so the band alone would pass it there.
if [[ ! "$freq_mhz" =~ ^[0-9]+\.[0-9]+$ ]] \
    || ! awk -v f="$freq_mhz" 'BEGIN { exit !(f + 0 >= 1 && f + 0 <= 100000) }'; then
    echo "FAIL: implausible calibrated frequency: ${freq_mhz} MHz (expected 1 .. 100000)" >&2
    exit 1
fi

# --- 3. The rate was not slept for -------------------------------------------
# At -v 4 proc_freq_resolve() prints "[D::proc_freq_resolve] tick rate <f> MHz
# from <source>". The tier banner from the same run says which source to
# expect: aarch64 (the neon tier) reads CNTFRQ_EL0, or measures over the run
# window where firmware left that register at 0; elsewhere main() opened the
# run window and the phiX run lasts well over its 1 ms floor. A sleep
# calibration fails either way.
"$bin" mem -v 4 "$ref" "$reads" > "$tmpdir/mem_v4.sam" 2> "$tmpdir/mem_v4.err" \
    || {
        echo "FAIL: bwa-mem3 mem -v 4 on phix.fa exited non-zero" >&2
        cat "$tmpdir/mem_v4.err" >&2
        exit 1
    }
if grep -q '^Executing in neon mode!!$' "$tmpdir/mem_v4.err"; then
    want_source="CNTFRQ_EL0 or the run window"
    want_re='^(CNTFRQ_EL0|the run window)$'
else
    want_source="the run window"
    want_re='^the run window$'
fi
source_line="$(grep -m1 '^\[D::proc_freq_resolve\] ' "$tmpdir/mem_v4.err" || true)"
got_source="${source_line##* MHz from }"
[[ -n "$source_line" && "$got_source" =~ $want_re ]] \
    || {
        echo "FAIL: expected the tick rate to come from ${want_source}; mem -v 4 reported:" >&2
        echo "      ${source_line:-<no [D::proc_freq_resolve] line>}" >&2
        exit 1
    }

# --- 4. `index` resolves the rate, and the time it prints is plausible -------
# `index` prints "Total time taken: <s>" from its own tick delta divided by
# proc_freq_resolve() (called in src/main.cpp, defined in src/profiling.cpp).
# proc_freq stays 0 until something resolves it, so a report that read the
# bare global would print `inf` here.
#
# The printed time is also checked against the wall time of the whole
# process, measured here with a sub-millisecond clock. The printed time covers
# part of the process, so it cannot exceed the wall time by more than read
# jitter: at the ~80% share it has on a laptop (the 1 Mb fixture's ~65 ms of
# indexing against process start-up), a rate more than ~1.5x too low fails.
# There is no lower bound: process start-up and a scheduler pause outside the
# indexing stage both add to the wall time and not to the printed time, so a
# small printed time is no evidence of a mis-scaled rate (a rate too high is
# pinned by check 2's band and by the unit test instead). The clock is
# CLOCK_MONOTONIC through perl's Time::HiRes, present on macOS and the Linux CI
# images, so a wall-clock step cannot skew it; without it the magnitude check
# is skipped (said so below), not passed. Index a copy so the fixture's own
# index files are left alone.
now_s() {
    perl -MTime::HiRes=clock_gettime,CLOCK_MONOTONIC \
        -e 'printf "%.6f\n", clock_gettime(CLOCK_MONOTONIC)' 2> /dev/null
}
idx_ref="$fixtures/synthetic_1mb.fa"
[[ -s "$idx_ref" ]] || {
    echo "FAIL: synthetic_1mb.fa missing at $idx_ref" >&2
    exit 1
}
cp "$idx_ref" "$tmpdir/idx.fa"
t_start="$(now_s || true)"
"$bin" index "$tmpdir/idx.fa" > /dev/null 2> "$tmpdir/index.err" \
    || {
        echo "FAIL: bwa-mem3 index on a copy of synthetic_1mb.fa exited non-zero" >&2
        cat "$tmpdir/index.err" >&2
        exit 1
    }
t_end="$(now_s || true)"
index_s="$(sed -nE 's/^Total time taken: ([^ ]+)$/\1/p' "$tmpdir/index.err" | head -n1)"
[[ -n "$index_s" ]] \
    || {
        echo "FAIL: 'Total time taken: ...' line missing from index stderr" >&2
        cat "$tmpdir/index.err" >&2
        exit 1
    }
[[ "$index_s" =~ ^[0-9]+\.[0-9]{4}$ ]] \
    || {
        echo "FAIL: malformed index 'Total time taken': ${index_s} (expected a finite, non-negative %0.4lf)" >&2
        exit 1
    }
if [[ -n "$t_start" && -n "$t_end" ]]; then
    wall_s="$(awk -v a="$t_start" -v b="$t_end" 'BEGIN { printf "%.4f", b - a }')"
    awk -v t="$index_s" -v w="$wall_s" 'BEGIN { exit !(t <= w * 1.1 + 0.01) }' \
        || {
            echo "FAIL: index 'Total time taken' ${index_s}s exceeds the ${wall_s}s wall time" >&2
            echo "      measured around it (expected at most wall*1.1+0.01): the tick rate is too low." >&2
            exit 1
        }
    index_note="index took ${index_s}s of ${wall_s}s wall"
else
    echo "SKIP: no sub-second clock (perl Time::HiRes); index time magnitude not checked" >&2
    index_note="index took ${index_s}s (magnitude not checked)"
fi

echo "PASS: startup does not sleep ($runs version calls in ${elapsed}s); calibrated @${freq_mhz} MHz from ${got_source}; ${index_note}"
