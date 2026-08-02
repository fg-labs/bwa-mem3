#!/usr/bin/env bash
# test/bam_compress_warn_test.sh
#
# `--bam=N` with N > 0 compresses BGZF on the single writer thread (bam_writer
# opens htslib with `wb<level>` and never calls hts_set_threads), so for large
# outputs that serial compression is usually what caps throughput. `mem` warns
# and points at `--bam=0` piped to a threaded compressor.
#
# The warning describes the *resolved* output setting, so it has to be emitted
# once, after option parsing, from the final opt->bam_mode / opt->bam_level --
# not from inside the parsing loop. The cases below pin that:
#
#   1. `--bam=6`             exactly one warning.
#   2. `--bam=6 --bam=0`     none: the last --bam wins and it is uncompressed,
#                            so there is nothing to warn about.
#   3. `--bam=0 --bam=6`     exactly one: the last --bam wins the other way.
#   4. `--bam=6 --bam=6`     exactly one, not one per occurrence.
#   5. `--bam=0`, bare       none.
#      `--bam`, no `--bam`
#   6. usage error           none: the run writes no output, so a warning about
#      (`--bam=6`, no reads) that output's compression level only buries the
#                            usage message.
#
# Every case also asserts the run's exit status, so a run that died before
# reaching the warning site cannot pass a "expected 0 warnings" case vacuously.
#
# Usage: test/bam_compress_warn_test.sh <bwa-mem3-binary> <fixtures-dir>

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

# Run `mem` with the given args and assert both the exit status and the number
# of compressed-BAM warnings on stderr.
#
#   expect <label> <want-rc-shape: zero|nonzero> <want-warnings> [mem args...]
#
# The exit status is asserted because it is what makes a zero-warning
# expectation meaningful: without it, a run that crashed during startup would
# also print zero warnings and satisfy cases 2, 5 and 6.
expect() {
    local label="$1" want_rc="$2" want_warns="$3"
    shift 3
    local err="$tmpdir/$label.err"
    local rc=0 warns grep_rc=0

    "$bin" mem "$@" > /dev/null 2> "$err" || rc=$?
    warns="$(grep -c 'WARNING: --bam=' "$err")" || grep_rc=$?

    # grep exits 1 when there are no matches -- an expected result here -- and
    # >1 only on a real error such as an unreadable capture file.
    if ((grep_rc > 1)); then
        echo "FAIL: [$label] could not scan $err for warnings (grep exit $grep_rc)" >&2
        exit 1
    fi

    case "$want_rc" in
        zero) ((rc == 0)) || {
            echo "FAIL: [$label] 'mem $*' exited $rc, expected 0" >&2
            cat "$err" >&2
            exit 1
        } ;;
        nonzero) ((rc != 0)) || {
            echo "FAIL: [$label] 'mem $*' exited 0, expected non-zero" >&2
            cat "$err" >&2
            exit 1
        } ;;
        *)
            echo "FAIL: [$label] bad want_rc '$want_rc' (use zero|nonzero)" >&2
            exit 1
            ;;
    esac

    if [[ "$warns" != "$want_warns" ]]; then
        echo "FAIL: [$label] 'mem $*' printed $warns compressed-BAM warning(s), expected $want_warns" >&2
        echo "---- stderr ----" >&2
        cat "$err" >&2
        exit 1
    fi
}

# --- 1. A single compressed --bam warns exactly once. -----------------------
expect single zero 1 --bam=6 "$ref" "$reads"

# --- 2/3. The last --bam wins, in both directions. --------------------------
expect overridden_off zero 0 --bam=6 --bam=0 "$ref" "$reads"
expect overridden_on zero 1 --bam=0 --bam=6 "$ref" "$reads"

# --- 4. Repeating the same compressed level warns once, not per occurrence. -
expect repeated zero 1 --bam=6 --bam=6 "$ref" "$reads"

# --- 5. Uncompressed BAM and SAM output are silent. -------------------------
expect level_zero zero 0 --bam=0 "$ref" "$reads"
expect bare_bam zero 0 --bam "$ref" "$reads"
expect sam zero 0 "$ref" "$reads"

# --- 6. A usage error produces no warning about output that never happens. --
# `--bam=6` with the reads file omitted fails the positional-argument check, so
# the run writes nothing and must exit non-zero without the warning.
expect usage_error nonzero 0 --bam=6 "$ref"

echo "PASS: compressed-BAM warning fires once for the resolved --bam level only"
