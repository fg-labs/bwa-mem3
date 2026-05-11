#!/usr/bin/env bash
# test/regression/bam_roundtrip.sh
#
# Regression: bwa-mem3 mem --bam=6 produces a BAM that decodes cleanly
# and has the same record count as the SAM path.
#
# Was: the "--bam=6 roundtrip smoke (phiX)" step inline in ci.yml.
#
# Inputs:
#   BWA_MEM3       — path to bwa-mem3 binary
#   CHR22_FA       — path to chr22.fa (pre-indexed with bwa-mem3 by caller)
#   CHR22_SIM_DIR  — directory containing holodeck reads.r[12].fastq.gz
#                    and the bwamem3.sam written by chr22_parity.sh
set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${CHR22_FA:?CHR22_FA must be set}"
: "${CHR22_SIM_DIR:?CHR22_SIM_DIR must be set}"

cd "$CHR22_SIM_DIR"
# bwamem3.sam is produced by chr22_parity.sh and is required for the
# record-count comparison below. Without this guard the `|| true` on the
# grep would mask a missing-file error and report a misleading "0 vs N"
# mismatch.
[ -f bwamem3.sam ] || {
    echo "FAIL: expected $CHR22_SIM_DIR/bwamem3.sam; run chr22_parity.sh first" >&2
    exit 1
}
"$BWA_MEM3" mem --bam=6 "$CHR22_FA" \
    reads.r1.fastq.gz reads.r2.fastq.gz > bwamem3.bam
samtools quickcheck bwamem3.bam
# grep -c exits 1 on zero matches, which would abort the script under
# `set -euo pipefail` before we can report a real "0 vs 0" result —
# match the `|| true` pattern used elsewhere (thread_determinism.sh).
sam_records=$(grep -cv '^@' bwamem3.sam || true)
bam_records=$(samtools view -c bwamem3.bam)
if [ "$sam_records" != "$bam_records" ]; then
    echo "FAIL: SAM ($sam_records) vs --bam=6 BAM ($bam_records) record count mismatch"
    exit 1
fi
echo "PASS: --bam=6 roundtrip ($bam_records records)"
