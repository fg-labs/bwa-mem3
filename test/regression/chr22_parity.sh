#!/usr/bin/env bash
# test/regression/chr22_parity.sh
#
# Regression: bwa and bwa-mem3 produce identical SAM output on a ~50k PE
# holodeck simulation against hg38 chr22.
#
# Was: the "Index chr22 with bwa and align" + "Align chr22 with bwa-mem3"
# + "Compare chr22 bwa vs bwa-mem3 (parity)" steps inline in ci.yml.
#
# Inputs:
#   BWA_MEM3       — path to the bwa-mem3 binary under test
#   CHR22_FA       — path to chr22.fa (pre-indexed with bwa-mem3 by caller)
#   BWA_CHR22_FA   — path to a chr22.fa copy pre-indexed with bwa
#                    (sidecars: .amb .ann .bwt .pac .sa). Lives next to
#                    CHR22_FA so the bwa-mem3 sidecars don't collide.
#   CHR22_SIM_DIR  — directory containing holodeck-simulated
#                    reads.r1.fastq.gz / reads.r2.fastq.gz
set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${CHR22_FA:?CHR22_FA must be set}"
: "${BWA_CHR22_FA:?BWA_CHR22_FA must be set}"
: "${CHR22_SIM_DIR:?CHR22_SIM_DIR must be set}"

bwa mem -t 4 "$BWA_CHR22_FA" \
    "$CHR22_SIM_DIR/reads.r1.fastq.gz" \
    "$CHR22_SIM_DIR/reads.r2.fastq.gz" \
    > "$CHR22_SIM_DIR/bwa.sam" 2>"$CHR22_SIM_DIR/bwa.log"

"$BWA_MEM3" mem -t 4 "$CHR22_FA" \
    "$CHR22_SIM_DIR/reads.r1.fastq.gz" \
    "$CHR22_SIM_DIR/reads.r2.fastq.gz" \
    > "$CHR22_SIM_DIR/bwamem3.sam" 2>"$CHR22_SIM_DIR/bwamem3.log"

cd "$CHR22_SIM_DIR"
normalize() { grep -v '^@PG' "$1" | grep -v '^@HD' | sed 's/\tMQ:i:[0-9]*//' | sed 's/\tHN:i:[0-9]*//'; }
normalize bwa.sam     | sort > bwa.normalized.sam
normalize bwamem3.sam | sort > bwamem3.normalized.sam
echo "bwa records:      $(grep -cv '^@' bwa.normalized.sam)"
echo "bwa-mem3 records: $(grep -cv '^@' bwamem3.normalized.sam)"
if diff bwa.normalized.sam bwamem3.normalized.sam > /dev/null 2>&1; then
    echo "PASS: bwa == bwa-mem3 on chr22 (~50k PE holodeck reads)"
else
    echo "FAIL: chr22 parity diff; first 40 lines:"
    diff bwa.normalized.sam bwamem3.normalized.sam | head -40
    exit 1
fi
