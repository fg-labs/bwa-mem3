#!/usr/bin/env bash
# test/regression/bwtseed_lockstep_driver_identity.sh
#
# Regression: the third-pass bwtseed lockstep driver (the default) and the scalar
# driver (BWA3_BWTSEED_LOCKSTEP=0) produce identical alignment records.
#
# The lockstep driver is on by default at every thread count, so no other
# regression run exercises the scalar driver end to end. This pins both drivers
# against each other on real reads, oversubscribed (-t = 2x the runner's CPUs,
# the regime that used to select the scalar driver). The function-level parity
# harness (test/bwtseed_lockstep_parity_test.cpp) pins the per-read SMEM output;
# this is the whole-run check.
#
# Guarded against a vacuous pass: each run's startup line must report the driver
# it was asked for ("lockstep: on (default)" / "lockstep: off
# (BWA3_BWTSEED_LOCKSTEP)"), so a pin that stopped being honored fails here
# rather than comparing the default against itself.
#
# Inputs:
#   BWA_MEM3       — path to bwa-mem3 binary
#   CHR22_FA       — path to chr22.fa (pre-indexed with bwa-mem3 by caller)
#   CHR22_SIM_DIR  — directory containing holodeck reads.r[12].fastq.gz
set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${CHR22_FA:?CHR22_FA must be set}"
: "${CHR22_SIM_DIR:?CHR22_SIM_DIR must be set}"

ncpu=$(getconf _NPROCESSORS_ONLN 2> /dev/null || echo 2)
threads=$((ncpu * 2))

cd "$CHR22_SIM_DIR"
env -u BWA3_BWTSEED_LOCKSTEP "$BWA_MEM3" mem -t "$threads" "$CHR22_FA" \
    reads.r1.fastq.gz reads.r2.fastq.gz 2> lockstep_on.err \
    | grep -v '^@PG' > lockstep_on.sam
BWA3_BWTSEED_LOCKSTEP=0 "$BWA_MEM3" mem -t "$threads" "$CHR22_FA" \
    reads.r1.fastq.gz reads.r2.fastq.gz 2> lockstep_off.err \
    | grep -v '^@PG' > lockstep_off.sam

if ! grep -q 'third-pass bwtseed lockstep: on (default)' lockstep_on.err; then
    echo "FAIL: default run did not report the lockstep driver on by default"
    grep 'third-pass bwtseed lockstep' lockstep_on.err || true
    exit 1
fi
if ! grep -q 'third-pass bwtseed lockstep: off (BWA3_BWTSEED_LOCKSTEP)' lockstep_off.err; then
    echo "FAIL: BWA3_BWTSEED_LOCKSTEP=0 run did not report the scalar driver"
    grep 'third-pass bwtseed lockstep' lockstep_off.err || true
    exit 1
fi
# Header-only output on both sides would compare equal and pass having checked
# nothing, so require alignment records from each run before comparing them.
for sam in lockstep_on.sam lockstep_off.sam; do
    if [ "$(grep -cv '^@' "$sam" || true)" -eq 0 ]; then
        echo "FAIL: $sam has no alignment records -- nothing was compared"
        exit 1
    fi
done
if ! cmp -s lockstep_on.sam lockstep_off.sam; then
    echo "FAIL: lockstep and scalar bwtseed drivers differ at -t $threads"
    diff lockstep_on.sam lockstep_off.sam | head -20
    exit 1
fi
echo "PASS: bwtseed lockstep driver identity at -t $threads ($(grep -cv '^@' lockstep_on.sam || true) records)"
