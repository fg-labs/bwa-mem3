#!/usr/bin/env bash
# test/regression/thread_determinism.sh
#
# Regression: bwa-mem3 -t 1 and -t 4 produce identical (sorted) output.
#
# Was: the "Thread-determinism smoke (phiX, -t 1 vs -t 4)" step inline in ci.yml.
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
"$BWA_MEM3" mem -t 1 "$CHR22_FA" \
    reads.r1.fastq.gz reads.r2.fastq.gz 2>/dev/null \
    | grep -v '^@PG' | sort > t1.sam
"$BWA_MEM3" mem -t 4 "$CHR22_FA" \
    reads.r1.fastq.gz reads.r2.fastq.gz 2>/dev/null \
    | grep -v '^@PG' | sort > t4.sam
if ! diff t1.sam t4.sam > /dev/null 2>&1; then
    echo "FAIL: -t 1 and -t 4 outputs differ after sort"
    diff t1.sam t4.sam | head -20
    exit 1
fi
echo "PASS: thread-determinism ($(grep -cv '^@' t1.sam || true) records)"
