#!/usr/bin/env bash
# test/regression/thread_determinism.sh
#
# Regression: bwa-mem3 output does not depend on -t.
#
#  1. -t 1 and -t 4 at the default batch size produce identical (sorted) output.
#     Was: the "Thread-determinism smoke (phiX, -t 1 vs -t 4)" step inline in
#     ci.yml.
#  2. At a FIXED -K, -t 1, -t 3 and -t 4 produce byte-identical UNSORTED output,
#     over several batches. -K pins the batches, but kt_for cuts the end of each
#     pass into smaller work items whose boundaries depend on -t, so which reads
#     share an item -- and a SIMD lane group in the batched extension and
#     mate-rescue kernels -- changes with -t. Identity here rests on every pair's
#     result being independent of what it is grouped with. Run on a copy of the
#     reads in which every fourth pair appears twice in a row, with
#     --dedup-reads on, so the duplicate-pair memo (whose representative and
#     duplicates can land in different items) is exercised too; on and off must
#     also agree.
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
    reads.r1.fastq.gz reads.r2.fastq.gz 2> /dev/null \
    | grep -v '^@PG' | sort > t1.sam
"$BWA_MEM3" mem -t 4 "$CHR22_FA" \
    reads.r1.fastq.gz reads.r2.fastq.gz 2> /dev/null \
    | grep -v '^@PG' | sort > t4.sam
if ! diff t1.sam t4.sam > /dev/null 2>&1; then
    echo "FAIL: -t 1 and -t 4 outputs differ after sort"
    diff t1.sam t4.sam | head -20
    exit 1
fi
echo "PASS: thread-determinism ($(grep -cv '^@' t1.sam || true) records)"

# (2) Fixed -K. The first 20000 pairs, every fourth one written twice in a row
# (25000 pairs, 5000 duplicated, ~7.5 Mbp); -K 1000000 cuts that into ~8
# batches of ~6700 reads, each a few full work items plus a -t dependent tail.
# `|| true` on the gzip side: head closes the pipe early (SIGPIPE).
fail() {
    echo "FAIL: $*"
    exit 1
}
dup_every4() { # <in.fastq.gz> <out.fq>
    { gzip -dc "$1" || true; } | head -n 80000 | paste - - - - \
        | awk '{ print } NR % 4 == 0 { print }' | tr '\t' '\n' > "$2"
}
dup_every4 reads.r1.fastq.gz fk.r1.fq
dup_every4 reads.r2.fastq.gz fk.r2.fq
[[ $(wc -l < fk.r1.fq) -eq 100000 && $(wc -l < fk.r2.fq) -eq 100000 ]] \
    || fail "fixed -K fixture is not 25000 pairs"
K=1000000
fk() { # <tag> <threads> <dedup mode>
    "$BWA_MEM3" mem -K "$K" -t "$2" --dedup-reads "$3" "$CHR22_FA" fk.r1.fq fk.r2.fq \
        2> "fk.$1.err" | grep -v '^@PG' > "fk.$1.sam" || fail "mem ($1)"
    [[ -s "fk.$1.sam" ]] || fail "empty output ($1)"
}
fk off.t1 1 off
for t in 1 3 4; do
    fk "on.t$t" "$t" on
    if ! cmp -s fk.off.t1.sam "fk.on.t$t.sam"; then
        diff fk.off.t1.sam "fk.on.t$t.sam" | head -20 || true
        fail "fixed -K $K: --dedup-reads on -t $t differs from off -t 1 (unsorted)"
    fi
done
nb=$(grep -c 'Calling mem_process_seqs' fk.on.t4.err || true)
[[ "$nb" -ge 4 ]] || fail "fixed -K run had $nb batches; expected several"
stats=$(BWAMEM3_DEDUP_READS_STATS=1 "$BWA_MEM3" mem -K "$K" -t 4 --dedup-reads on "$CHR22_FA" \
    fk.r1.fq fk.r2.fq 2>&1 > /dev/null | grep 'dedup-reads-stats' || true)
dup=$(awk -F'dup=' '{split($2, a, " "); print a[1]}' <<< "$stats")
[[ "${dup:-0}" -gt 0 ]] || fail "fixed -K fixture exercised no duplicate pairs"
echo "PASS: fixed -K $K, $nb batches: -t 1 == -t 3 == -t 4 unsorted, --dedup-reads on == off ($dup duplicate pairs, $(grep -cv '^@' fk.off.t1.sam || true) records)"

# (3) Output streaming. The process step writes each finished work item's SAM
# records while the worker_sam pass is still running, and the write step writes
# the rest (stream_emit_prefix in src/fastmap.cpp). On must equal off
# (BWA_MEM3_STREAM_OUTPUT=0) byte for byte, unsorted, over the same several
# batches: paired and single-end, at -t 2, 3 and 4, to a file and through a
# consumer slower than the aligner (about 10k lines/s), so the output stays
# backed up and a batch's early writes keep having to wait for the previous
# batch's write step -- dropping that wait interleaves two batches' records,
# which this catches. Also with every cohort sliced (BWA_MEM3_COHORT_SLICE_ALL),
# which emits through mem_pair_and_emit_cohort instead of mem_process_seqs.
# BWA_MEM3_STREAM_STATS=1 reports how many reads each batch wrote early: to a file, at least two
# batches must have (one would pass even if only the first batch ever got its
# turn), or the comparison proves nothing; at -t 1, which keeps its single
# pipeline thread, none may.
slow_reader() {
    perl -e 'while (<STDIN>) { print; select(undef, undef, undef, 0.005) unless $. % 50 }'
}
so() { # <tag> <sink: file|stall> <env...> -- <mem args...>
    local tag=$1 sink=$2
    shift 2
    local envs=()
    while [[ $1 != -- ]]; do
        envs+=("$1")
        shift
    done
    shift
    if [[ $sink == stall ]]; then
        env BWA_MEM3_STREAM_STATS=1 "${envs[@]}" "$BWA_MEM3" mem -K "$K" "$@" 2> "so.$tag.err" \
            | slow_reader | grep -v '^@PG' > "so.$tag.sam" || fail "mem ($tag)"
    else
        env BWA_MEM3_STREAM_STATS=1 "${envs[@]}" "$BWA_MEM3" mem -K "$K" -o "so.$tag.out" "$@" 2> "so.$tag.err" \
            || fail "mem ($tag)"
        grep -v '^@PG' "so.$tag.out" > "so.$tag.sam" || fail "empty output ($tag)"
    fi
    [[ -s "so.$tag.sam" ]] || fail "empty output ($tag)"
}
early() { awk '/stream_output:/ { s += $3 } END { print s + 0 }' "so.$1.err"; }
early_batches() { awk '/stream_output:/ && $3 > 0 { n++ } END { print n + 0 }' "so.$1.err"; }
for end in pe se; do
    if [[ $end == pe ]]; then in=(fk.r1.fq fk.r2.fq); else in=(fk.r1.fq); fi
    so "$end.off" file BWA_MEM3_STREAM_OUTPUT=0 -- -t 4 "$CHR22_FA" "${in[@]}"
    [[ $(grep -c 'stream_output:' "so.$end.off.err" || true) -eq 0 ]] \
        || fail "streaming $end: BWA_MEM3_STREAM_OUTPUT=0 still streamed"
    so "$end.t1" file BWA_MEM3_STREAM_OUTPUT=1 -- -t 1 "$CHR22_FA" "${in[@]}"
    cmp -s "so.$end.off.sam" "so.$end.t1.sam" || fail "streaming $end: -t 1 differs from off"
    [[ $(early "$end.t1") -eq 0 ]] || fail "streaming $end: -t 1 wrote early; it has one pipeline thread"
    for t in 2 3 4; do
        for sink in file stall; do
            tag="$end.t$t.$sink"
            so "$tag" "$sink" BWA_MEM3_STREAM_OUTPUT=1 -- -t "$t" "$CHR22_FA" "${in[@]}"
            if ! cmp -s "so.$end.off.sam" "so.$tag.sam"; then
                diff "so.$end.off.sam" "so.$tag.sam" | head -20 || true
                fail "streaming $tag differs from off"
            fi
            [[ $(early "$tag") -gt 0 ]] || fail "streaming $tag wrote nothing early"
            if [[ $sink == file ]]; then
                [[ $(early_batches "$tag") -ge 2 ]] \
                    || fail "streaming $tag: only $(early_batches "$tag") batch(es) wrote early"
            fi
        done
    done
done
for sink in file stall; do
    tag="pe.slice.$sink"
    so "$tag" "$sink" BWA_MEM3_COHORT_SLICE_ALL=1 -- -t 3 "$CHR22_FA" fk.r1.fq fk.r2.fq
    grep -q 'slicing every' "so.$tag.err" || fail "streaming $tag: the cohorts were not sliced"
    cmp -s so.pe.off.sam "so.$tag.sam" || fail "streaming $tag (every cohort sliced) differs from off"
    [[ $(early "$tag") -gt 0 ]] || fail "streaming $tag wrote nothing early"
done
# BAM output is written by the write step alone: its records are bam1_t, not
# the sam strings the process step would write.
BWA_MEM3_STREAM_STATS=1 "$BWA_MEM3" mem -K "$K" -t 4 --bam -o so.bam.out "$CHR22_FA" fk.r1.fq fk.r2.fq 2> so.bam.err \
    || fail "mem (--bam)"
[[ -s so.bam.out ]] || fail "empty output (--bam)"
[[ $(grep -c 'stream_output:' so.bam.err || true) -eq 0 ]] || fail "streaming: --bam output was streamed"
echo "PASS: output streaming on == off (paired and single-end, -t 1..4, file and slow pipe, sliced cohorts, not --bam; $(early pe.t4.stall) paired reads written early at -t 4)"
