#!/usr/bin/env bash
# test/proc_freq_calibration_test.sh
#
# `proc_freq` (the __rdtsc() tick rate that every timing report divides by)
# is calibrated once at startup. It used to be measured with a literal
# `sleep(1)`, which charged a full second of wall time to EVERY invocation
# -- including `version`, `shm -l`, and usage errors, none of which ever
# read proc_freq. This test guards both halves of the replacement:
#
#   1. Startup is cheap: N back-to-back `version` calls finish in well
#      under N seconds. A reintroduced sleep(1) fails here immediately.
#   2. The calibration is still correct: a real `mem` run reports a
#      plausible processor frequency. A proc_freq of 0 prints @0.000000 MHz
#      and a broken window prints an absurd value; either fails here.
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
    echo "      Startup is sleeping again -- see calibrate_proc_freq() in src/main.cpp." >&2
    exit 1
fi

# --- 2. The calibration still produces a plausible frequency ----------------
# display_stats() prints "Processor is running @<freq> MHz" from proc_freq
# (src/profiling.cpp: `proc_freq*1.0/1e6` under %lf). proc_freq == 0 renders
# as @0.000000 MHz -- it is the *time* lines below it that divide by
# proc_freq and print `inf` -- and a mis-scaled window renders orders of
# magnitude off; require a value inside a deliberately generous band that no
# real tick source (a ~24 MHz arm64 CNTVCT_EL0 through a multi-GHz x86 TSC)
# falls outside of.
"$bin" mem "$ref" "$reads" > "$tmpdir/mem.sam" 2> "$tmpdir/mem.err" \
    || {
        echo "FAIL: bwa-mem3 mem on phix.fa exited non-zero" >&2
        cat "$tmpdir/mem.err" >&2
        exit 1
    }

# Capture the field verbatim ([^ ]+, not just digits) so whatever %lf
# actually printed -- 0.000000, an absurd magnitude, or a non-numeric
# nan/inf -- lands in the failure message instead of the match silently
# failing and reporting a missing line.
freq_mhz="$(sed -nE 's/^Processor is running @([^ ]+) MHz$/\1/p' "$tmpdir/mem.err" | head -n1)"
[[ -n "$freq_mhz" ]] \
    || {
        echo "FAIL: 'Processor is running @... MHz' line missing from mem stderr" >&2
        cat "$tmpdir/mem.err" >&2
        exit 1
    }

awk -v f="$freq_mhz" 'BEGIN { exit !(f + 0 >= 1 && f + 0 <= 100000) }' \
    || {
        echo "FAIL: implausible calibrated frequency: ${freq_mhz} MHz (expected 1 .. 100000)" >&2
        exit 1
    }

echo "PASS: startup does not sleep ($runs version calls in ${elapsed}s); calibrated @${freq_mhz} MHz"
