#!/usr/bin/env bash
# test/regression/short_read_smoke.sh
#
# Regression: ASAN-instrumented bwa-mem3 mem on short, variable-length
# reads from a repeat-dense region of chr22 must complete cleanly,
# guarding the wsize_mem byte-vs-entry capacity-tracker class of bug
# fixed in PR #100.
#
# Bug recap. mmc->wsize_mem is set in *bytes* by mem_kernel1_core to
# size enc_qdb (= tot_len), then repurposed by mem_collect_smem as an
# *SMEM-entry count* to size matchArray. On subsequent batches the gate
# `if (tot_len >= wsize_mem)` compares bytes vs SMEM entries. Whenever
# the carried-over SMEM count exceeds the next batch's byte total,
# enc_qdb is left undersized and the next
# `memcpy(enc_qdb + offset, seq, l_seq)` in mem_collect_smem overflows
# its buffer.
#
# Two ingredients are required to surface the bug on CI-scale inputs:
#   1. High SMEM density per base — short reads with non-trivial error,
#      aligned to a repeat-rich slice (pericentromeric chr22:14M-20M).
#   2. Variable per-batch byte totals — variable read length (25-50 bp)
#      across the FASTQ so successive batches have very different
#      tot_len, letting one batch's stale wsize_mem (SMEM count) exceed
#      the next batch's tot_len (bytes).
#
# Uniform-length short reads from a low-error chr22 sim do NOT trigger
# the bug: every batch has the same tot_len, so the gate compares
# wsize_mem-as-bytes against tot_len-as-bytes correctly, and the
# misuse downstream in mem_collect_smem never lands on an undersized
# enc_qdb. The fixture was discovered the hard way; do not "simplify"
# it back to a single-length truncate of the canonical chr22 sim.
#
# Detection: build with `make ASAN=1` (forces USE_MIMALLOC=0 — mimalloc
# and asan can't coexist) so the heap-buffer-overflow in
# mem_collect_smem is caught directly. An un-instrumented build merely
# undersizes enc_qdb and the corruption-to-SIGSEGV translation depends
# on the allocator's adjacent-metadata layout, which is fragile across
# libcs and not portable enough for CI.
#
# Inputs (env vars):
#   BWA_MEM3       — path to the ASAN-instrumented bwa-mem3 binary
#   CHR22_FA       — path to chr22.fa (pre-indexed with bwa-mem3 by caller)
#   CHR22_SIM_DIR  — directory for fixture-private intermediates
#   PIXI_MANIFEST  — path to test/holodeck/pixi.toml (holodeck + samtools)

set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${CHR22_FA:?CHR22_FA must be set}"
: "${CHR22_SIM_DIR:?CHR22_SIM_DIR must be set}"
: "${PIXI_MANIFEST:?PIXI_MANIFEST must be set}"

mkdir -p "$CHR22_SIM_DIR"

# Surface the bwa-mem3 stderr log on any error. Without this, a non-zero
# exit from bwa-mem3 (the exact failure mode this fixture exists to
# catch) is killed by `set -e` before the ASAN report is echoed,
# leaving CI with only "exit 134" and no backtrace.
short_log="$CHR22_SIM_DIR/short_dense_var.log"
trap 'rc=$?; if [ $rc -ne 0 ] && [ -f "$short_log" ]; then
        echo "--- bwa-mem3 stderr ($short_log) ---" >&2
        tail -n 200 "$short_log" >&2
      fi; exit $rc' ERR

# Pericentromeric chr22:14M-20M is ~23% Ns interleaved with repeat-rich
# sequence; reads simulated from here hit many short SMEMs per base.
bed="$CHR22_SIM_DIR/chr22_dense.bed"
printf 'chr22\t14000000\t20000000\n' > "$bed"

# 1x coverage of 6 Mb at 45 bp = ~133k reads. The bug fires across
# successive 1024-read SMEM batches within a single mem_kernel1_core
# call, so coverage well below the "real workload" (45M reads in PR
# #100's reproducer) still produces plenty of inner batches per thread.
# 5-10% error rate keeps SMEMs short and numerous; --max-n-frac is
# loosened because the dense region is N-rich; --single-end emits R1
# only.
pixi run --frozen --manifest-path "$PIXI_MANIFEST" -- \
  holodeck simulate \
    -r "$CHR22_FA" \
    -b "$bed" \
    -o "$CHR22_SIM_DIR/short_dense" \
    --coverage 1 \
    --read-length 45 \
    --fragment-mean 45 --fragment-stddev 5 --min-fragment-length 30 \
    --seed 42 \
    --simple-names \
    --single-end \
    --min-error-rate 0.05 --max-error-rate 0.10 \
    --max-n-frac 0.5

# Truncate per-read to a deterministic length in [25, 50]. Variable
# lengths across the FASTQ are the second ingredient: a batch carrying
# wsize_mem from a long-batch tail can leave enc_qdb undersized for the
# next batch whose tot_len byte total is smaller but whose SMEM count
# has not yet shrunk below the carried-over value. Output is left
# uncompressed (bwa-mem3 is the sole consumer and decompresses anyway);
# mawk's END block writes the record count to a sidecar so we don't
# pay a second pass over the file just to count.
short_fq="$CHR22_SIM_DIR/short_dense_var.fq"
short_count="$CHR22_SIM_DIR/short_dense_var.count"
gunzip -c "$CHR22_SIM_DIR/short_dense.r1.fastq.gz" \
  | mawk -v count_out="$short_count" 'BEGIN{srand(42)}
      NR%4==1 { cur_len = 25 + int(26*rand()); print; next }
      NR%4==2 || NR%4==0 { print substr($0, 1, cur_len); next }
      { print }
      END { print NR/4 > count_out }' \
  > "$short_fq"

short_records=$(cat "$short_count")
echo "short-read fastq: $short_records records (25-50bp, pericentromeric)"

# Refuse to silently "pass" if holodeck produced essentially nothing —
# without this, an empty/near-empty FASTQ would let bwa-mem3 exit 0 on
# zero input and the >= record-count check below would trivially hold.
# Threshold is two orders of magnitude below the expected ~133k.
if [ "$short_records" -lt 1000 ]; then
    echo "FAIL: short-read fastq has $short_records records, expected ~133k (holodeck regression?)" >&2
    exit 1
fi

# Fail fast on the first ASAN/LSan report so the trap's tail of
# $short_log captures the actual overflow rather than downstream
# cascade noise. LSan (bundled with ASan on Linux) is left on its
# default so the fixture also catches per-thread leak regressions
# (e.g. the FMI_search lockstep SMEM caches fixed in #116).
short_sam="$CHR22_SIM_DIR/short_dense_var.sam"
ASAN_OPTIONS=abort_on_error=1:halt_on_error=1 \
"$BWA_MEM3" mem -t 4 "$CHR22_FA" "$short_fq" \
    > "$short_sam" 2>"$short_log"

# `|| true` because grep -c returns 1 on zero matches (a header-only SAM
# with set -e would die here instead of at the explicit FAIL below).
sam_records=$(grep -cv '^@' "$short_sam" || true)
echo "short-read SAM:   $sam_records records"

# Every input read must produce at least one primary record. The
# pre-fix corruption truncates the SAM mid-stream (un-instrumented
# builds) or aborts (ASAN builds), so the equality is a load-bearing
# check, not a tautology. Tolerate >= because secondary/supplementary
# alignments add extra lines on the chr22 pericentromere.
if [ "$sam_records" -lt "$short_records" ]; then
    echo "FAIL: short-read SE produced $sam_records SAM records, expected >= $short_records" >&2
    exit 1
fi

echo "PASS: short-read SE smoke ($sam_records records, 25-50bp pericentromeric)"
