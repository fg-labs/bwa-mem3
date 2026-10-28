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
# BWAMEM3_DEDUP_READS_STATS only counts (stderr, at exit), so the -t 4 run below
# also yields the duplicate-pair count instead of a run of its own.
for t in 1 3 4; do
    BWAMEM3_DEDUP_READS_STATS=1 fk "on.t$t" "$t" on
    if ! cmp -s fk.off.t1.sam "fk.on.t$t.sam"; then
        diff fk.off.t1.sam "fk.on.t$t.sam" | head -20 || true
        fail "fixed -K $K: --dedup-reads on -t $t differs from off -t 1 (unsorted)"
    fi
done
nb=$(grep -c 'Calling mem_process_seqs' fk.on.t4.err || true)
[[ "$nb" -ge 4 ]] || fail "fixed -K run had $nb batches; expected several"
stats=$(grep 'dedup-reads-stats' fk.on.t4.err || true)
dup=$(awk -F'dup=' '{split($2, a, " "); print a[1]}' <<< "$stats")
[[ "${dup:-0}" -gt 0 ]] || fail "fixed -K fixture exercised no duplicate pairs"
echo "PASS: fixed -K $K, $nb batches: -t 1 == -t 3 == -t 4 unsorted, --dedup-reads on == off ($dup duplicate pairs, $(grep -cv '^@' fk.off.t1.sam || true) records)"

# (3) Output streaming. The process step writes each finished work item's SAM
# records (and, at --bam=0, the last batch's BAM records) while the worker_sam
# pass is still running, and the write step writes the rest
# (stream_emit_prefix in src/fastmap.cpp). On must equal off
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
#
# Two inputs, each with its own off run. The file runs read all of fk (8 paired,
# 4 single-end batches at -K), which the at-least-two-batches check needs. The
# slow pipe runs at the consumer's pace, not the aligner's -- about 5 s per
# paired pass over fk -- so it, -t 1 and --bam read sk, the first 10000 pairs of
# fk: the same -K and so the same batches, 3 paired and 2 single-end, and every
# backed-up batch after the first still has to wait for the one before it.
head -n 40000 fk.r1.fq > sk.r1.fq
head -n 40000 fk.r2.fq > sk.r2.fq
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
set_in() { # <input: fk|sk> <end: pe|se>; sets in=(the fastqs)
    in=("$1.r1.fq")
    [[ $2 == se ]] || in+=("$1.r2.fq")
}
for end in pe se; do
    for fx in fk sk; do
        set_in "$fx" "$end"
        so "$end.$fx.off" file BWA_MEM3_STREAM_OUTPUT=0 -- -t 4 "$CHR22_FA" "${in[@]}"
        [[ $(grep -c 'stream_output:' "so.$end.$fx.off.err" || true) -eq 0 ]] \
            || fail "streaming $end: BWA_MEM3_STREAM_OUTPUT=0 still streamed"
    done
    set_in sk "$end"
    so "$end.t1" file BWA_MEM3_STREAM_OUTPUT=1 -- -t 1 "$CHR22_FA" "${in[@]}"
    cmp -s "so.$end.sk.off.sam" "so.$end.t1.sam" || fail "streaming $end: -t 1 differs from off"
    [[ $(early "$end.t1") -eq 0 ]] || fail "streaming $end: -t 1 wrote early; it has one pipeline thread"
    for t in 2 3 4; do
        for sink in file stall; do
            if [[ $sink == file ]]; then fx=fk; else fx=sk; fi
            set_in "$fx" "$end"
            tag="$end.t$t.$sink"
            so "$tag" "$sink" BWA_MEM3_STREAM_OUTPUT=1 -- -t "$t" "$CHR22_FA" "${in[@]}"
            if ! cmp -s "so.$end.$fx.off.sam" "so.$tag.sam"; then
                diff "so.$end.$fx.off.sam" "so.$tag.sam" | head -20 || true
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
    if [[ $sink == file ]]; then fx=fk; else fx=sk; fi
    tag="pe.slice.$sink"
    so "$tag" "$sink" BWA_MEM3_COHORT_SLICE_ALL=1 -- -t 3 "$CHR22_FA" "$fx.r1.fq" "$fx.r2.fq"
    grep -q 'slicing every' "so.$tag.err" || fail "streaming $tag: the cohorts were not sliced"
    cmp -s "so.pe.$fx.off.sam" "so.$tag.sam" || fail "streaming $tag (every cohort sliced) differs from off"
    [[ $(early "$tag") -gt 0 ]] || fail "streaming $tag wrote nothing early"
done
# --bam=0 streams the LAST batch only (the one whose read hit end of input);
# earlier batches, and every batch at --bam>=1, are written by the write step
# alone. The hard gate is the BAM file: streaming on must equal off byte for
# byte (same -o name, so @PG CL matches). BWA_MEM3_STREAM_STATS prints one line
# per batch armed to stream, so at --bam=0 there must be exactly one, for the
# last batch (its read count must be the last batch's). Whether that batch
# actually writes early is timing dependent (the previous batch's write must
# finish first), so it is reported, not required; a few attempts give it the
# chance, and every attempt is compared byte for byte.
bam_run() { # <lvl> <out> <K> <input prefix> [env...]
    local lvl=$1 out=$2 k=$3 fx=$4
    shift 4
    env BWA_MEM3_STREAM_STATS=1 "$@" "$BWA_MEM3" mem -K "$k" -t 4 --bam="$lvl" -o so.bam.cur \
        "$CHR22_FA" "$fx.r1.fq" "$fx.r2.fq" 2> "$out.err" || fail "mem (--bam=$lvl $out)"
    [[ -s so.bam.cur ]] || fail "empty output (--bam=$lvl $out)"
    mv so.bam.cur "$out"
}
stat_lines() { grep -c 'stream_output:' "$1.err" || true; }
# reads in the last batch (cohort): read_chunk lines tagged "(cohort slice)"
# accumulate into the cohort that the next untagged line completes
last_batch_n() {
    awk '/read_chunk:/ { n = $0; sub(/.*nseq: /, "", n); acc += n + 0
                         if ($0 !~ /cohort slice/) { if (acc > 0) last = acc; acc = 0 } }
         END { print last + 0 }' "$1.err"
}
bam0_leg() { # <tag> <K> <input prefix> [env...]: --bam=0 on == off, one stats line, for the last batch
    local tag=$1 k=$2 fx=$3
    shift 3
    bam_run 0 "so.$tag.off" "$k" "$fx" BWA_MEM3_STREAM_OUTPUT=0 "$@"
    [[ $(stat_lines "so.$tag.off") -eq 0 ]] || fail "$tag: streamed with BWA_MEM3_STREAM_OUTPUT=0"
    local attempt line early=0
    for attempt in 1 2 3 4 5; do
        bam_run 0 "so.$tag.on" "$k" "$fx" "$@"
        cmp -s "so.$tag.off" "so.$tag.on" || fail "$tag: --bam=0 streaming on differs from off (attempt $attempt)"
        [[ $(stat_lines "so.$tag.on") -eq 1 ]] \
            || fail "$tag: $(stat_lines "so.$tag.on") batches armed to stream at --bam=0, expected only the last"
        line=$(grep 'stream_output:' "so.$tag.on.err")
        [[ $(awk '{ print $5 }' <<< "$line") == "$(last_batch_n "so.$tag.on")" ]] \
            || fail "$tag: the batch armed to stream is not the last one ($line)"
        early=$(awk '{ print $3 }' <<< "$line")
        [[ $early -gt 0 ]] && break
    done
    echo "  $tag: last batch wrote $early reads early (attempt $attempt; informational)"
}
# Exact end of input: every batch is exactly full (fixed-length reads, -K a
# whole number of pairs that divides the input), so the last batch's read is
# not short and end of input only shows on the next, empty read. That batch is
# documented as not streamed; output must still match.
rl=$(awk 'NR % 4 == 2 { print length($0) }' sk.r1.fq sk.r2.fq | sort -u)
[[ $(wc -l <<< "$rl") -eq 1 ]] || fail "exact-boundary leg needs fixed-length reads (lengths: $(tr '\n' ' ' <<< "$rl"))"
npairs=$(($(wc -l < sk.r1.fq) / 4))
((npairs % 4 == 0)) || fail "exact-boundary leg needs a pair count divisible by 4 ($npairs)"
KX=$((npairs * 2 * rl / 4))
bam_run 0 so.bamx.off "$KX" sk BWA_MEM3_STREAM_OUTPUT=0
bam_run 0 so.bamx.on "$KX" sk
[[ $(grep 'read_chunk:' so.bamx.on.err | grep -v 'cohort slice' | grep -vc 'nseq: 0$' || true) -eq 4 ]] \
    || fail "exact-boundary leg: expected 4 full batches at -K $KX"
cmp -s so.bamx.off so.bamx.on || fail "exact-boundary: --bam=0 streaming on differs from off"
[[ $(stat_lines so.bamx.on) -eq 0 ]] \
    || fail "exact-boundary: a batch was armed to stream; an exactly full last batch is documented as not streamed"
echo "  exact end of input (-K $KX, 4 full batches): not streamed, BAM identical"
bam0_leg bam0 "$K" sk
# Every cohort sliced: the last cohort is paired and emitted through
# mem_pair_and_emit_cohort rather than mem_process_seqs.
bam0_leg bam0.slice "$K" fk BWA_MEM3_COHORT_SLICE_ALL=1
grep -q 'slicing every' so.bam0.slice.on.err || fail "bam0.slice: the cohorts were not sliced"
# --bam=6 (compressed) never streams: the level gates it.
bam_run 6 so.bam6.off "$K" sk BWA_MEM3_STREAM_OUTPUT=0
bam_run 6 so.bam6.on "$K" sk
cmp -s so.bam6.off so.bam6.on || fail "--bam=6 streaming on differs from off"
[[ $(stat_lines so.bam6.on) -eq 0 ]] || fail "--bam=6 was armed to stream (level 0 only)"
echo "PASS: output streaming on == off (paired and single-end, -t 1..4, file and slow pipe, sliced cohorts, --bam=0 last batch only; $(early pe.t4.stall) paired reads written early at -t 4)"
