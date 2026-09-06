#!/usr/bin/env bash
# test/regression/bam_threads_byte_identity.sh
#
# Regression: the --bam-threads BGZF compression pool is byte-identical to the
# serial writer. htslib's ordered thread pool emits blocks in dispatch order, so
# --bam=6 --bam-threads=N must decode to exactly the same records as
# --bam=6 --bam-threads=0 (the same guarantee `samtools -@` relies on). This is
# the correctness contract that justifies the auto-parallelized default, so it
# is pinned here rather than left to manual checking.
#
# Also checks: the threaded output is deterministic (same compressed bytes across
# runs), and the --bam-threads parser rejects malformed / out-of-range values.
#
# Uses the chr22 holodeck fixture (multi-block BAM — a single BGZF block would
# not exercise the ordered-tpool block-ordering path).
#
# Inputs:
#   BWA_MEM3       — path to bwa-mem3 binary
#   CHR22_FA       — path to chr22.fa (pre-indexed with bwa-mem3 by caller)
#   CHR22_SIM_DIR  — directory containing holodeck reads.r[12].fastq.gz
set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${CHR22_FA:?CHR22_FA must be set}"
: "${CHR22_SIM_DIR:?CHR22_SIM_DIR must be set}"

cd "$CHR22_SIM_DIR"
[ -f reads.r1.fastq.gz ] && [ -f reads.r2.fastq.gz ] || {
    echo "FAIL: expected $CHR22_SIM_DIR/reads.r[12].fastq.gz" >&2
    exit 1
}

run() { # run <out.bam> <bam-threads>
    "$BWA_MEM3" mem --bam=6 --bam-threads="$2" "$CHR22_FA" \
        reads.r1.fastq.gz reads.r2.fastq.gz > "$1"
    samtools quickcheck "$1"
}

run serial.bam 0
run threaded_a.bam 4
run threaded_b.bam 4

# 1. Byte-identity of the decoded records (the contract). `samtools view` with no
#    -h emits records only, so the @PG command-line difference (--bam-threads=0
#    vs =4) does not enter the comparison.
serial_md5=$(samtools view serial.bam | md5sum | cut -d' ' -f1)
threaded_md5=$(samtools view threaded_a.bam | md5sum | cut -d' ' -f1)
serial_n=$(samtools view -c serial.bam)
threaded_n=$(samtools view -c threaded_a.bam)
if [ "$serial_md5" != "$threaded_md5" ] || [ "$serial_n" != "$threaded_n" ]; then
    echo "FAIL: --bam-threads=4 records differ from serial" \
        "(serial $serial_n/$serial_md5 vs threaded $threaded_n/$threaded_md5)" >&2
    exit 1
fi

# 2. Threaded output is deterministic: identical flags -> identical compressed
#    bytes (no block-ordering race in the tpool).
if ! cmp -s threaded_a.bam threaded_b.bam; then
    echo "FAIL: two --bam-threads=4 runs produced different compressed BAM (ordering race)" >&2
    exit 1
fi

# 3. The parser rejects malformed / out-of-range values (must exit non-zero and
#    write nothing usable). MAX_THREADS is 256.
for bad in abc -1 999 4x; do
    if "$BWA_MEM3" mem --bam=6 --bam-threads="$bad" "$CHR22_FA" \
        reads.r1.fastq.gz reads.r2.fastq.gz > /dev/null 2> /dev/null; then
        echo "FAIL: --bam-threads=$bad was accepted, expected a parse error" >&2
        exit 1
    fi
done

echo "PASS: --bam-threads byte-identical to serial ($serial_n records), deterministic, and range-validated"
