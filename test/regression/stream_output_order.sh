#!/usr/bin/env bash
# test/regression/stream_output_order.sh
#
# Output streaming must never write a batch's records before the previous
# batch's write step has finished. The process step writes each finished work
# item's SAM records while the worker_sam pass runs (stream_emit_prefix in
# src/fastmap.cpp), but only once it is that batch's turn; otherwise it leaves
# the records to its own write step.
#
# Under ordinary timing the previous batch's write is usually done long before
# the next batch's records are ready, so dropping that wait rarely changes the
# output and the identity checks in thread_determinism.sh would not see it.
# This forces the case: the BWAMEM3_TESTING_WRITE_DELAY_MS hook holds every
# write step back, so each batch's pass ends before the previous batch's
# records are out. Then:
#   1. output must still equal BWA_MEM3_STREAM_OUTPUT=0 byte for byte
#      (paired and single-end, -t 2, 3 and 4) -- writing without the turn
#      would put a batch ahead of the one before it;
#   2. without the delay, every batch writes early (BWA_MEM3_STREAM_STATS=1
#      reports it per batch), so streaming is engaged and (1) is not vacuous.
#      BWAMEM3_TESTING_HOLD_LAST_ITEM keeps each pass's last item from leaving
#      until a prefix has been written, so this does not depend on the caller
#      running before every item has finished.
#   3. --bam=0 streams only its last batch: with the delay, the BAM file must
#      equal the BWA_MEM3_STREAM_OUTPUT=0 one byte for byte (same -o name, so
#      @PG CL matches) and exactly one batch is armed to stream; with the last
#      item held and no delay, that batch must write early.
# How many batches write early under the delay is not checked: a pass that
# outlasts the delay on a slow host may legitimately see the previous write
# finish and stream, so only the output order is asserted.
# The reference and reads are generated here (one random 300 kb contig, 8000
# read pairs), and -K cuts them into 8 batches of 2000 reads, each several
# work items, so a batch has finished items to offer while its pass runs.
#
# Inputs:
#   BWA_MEM3_TESTING — path to a binary built with `make TESTING_BUILD=1`

set -Eeuo pipefail
# A command that fails outside a `|| fail` guard (fixture generation, say)
# still ends the script with a FAIL: line.
trap 'echo "FAIL: line $LINENO: $BASH_COMMAND exited $?"' ERR

[[ -n "${BWA_MEM3_TESTING:-}" ]] || {
    echo "FAIL: BWA_MEM3_TESTING must be set (a TESTING_BUILD=1 binary)"
    exit 1
}
# The script works in a temporary directory, so a relative path (make
# test-injection passes ./bwa-mem3) must be resolved before the cd.
[[ -x "$BWA_MEM3_TESTING" ]] || {
    echo "FAIL: BWA_MEM3_TESTING=$BWA_MEM3_TESTING is not an executable"
    exit 1
}
BWA_MEM3_TESTING="$(cd "$(dirname "$BWA_MEM3_TESTING")" && pwd)/$(basename "$BWA_MEM3_TESTING")"

OUT_DIR="$(mktemp -d -t bwamem3-stream-order-XXXXXX)"
trap 'rm -rf "$OUT_DIR"' EXIT
cd "$OUT_DIR"

fail() { # <message> [stderr file of the failing command]
    echo "FAIL: $1"
    if [[ -n ${2:-} && -s $2 ]]; then
        echo "--- stderr ($2, last 40 lines) ---"
        tail -n 40 "$2"
    fi
    exit 1
}

awk 'BEGIN {
    srand(7); split("A C G T", b, " ")
    print ">chrT"
    for (i = 0; i < 300000; i++) {
        printf "%s", b[int(rand() * 4) + 1]
        if (i % 60 == 59) printf "\n"
    }
    printf "\n"
}' > ref.fa
awk 'function rc(s,    i, o, c) {
    o = ""
    for (i = length(s); i > 0; i--) {
        c = substr(s, i, 1)
        o = o (c == "A" ? "T" : c == "C" ? "G" : c == "G" ? "C" : "A")
    }
    return o
}
function mut(s,    p) { # one substitution in about a third of the reads
    if (rand() >= 0.33) return s
    p = int(rand() * length(s)) + 1
    return substr(s, 1, p - 1) (substr(s, p, 1) == "A" ? "C" : "A") substr(s, p + 1)
}
!/^>/ { ref = ref $0 }
END {
    srand(11); q = ""
    for (i = 0; i < 100; i++) q = q "I"
    for (k = 0; k < 8000; k++) {
        ins = 250 + int(rand() * 100)
        pos = int(rand() * (length(ref) - ins)) + 1
        r1 = mut(substr(ref, pos, 100))
        r2 = mut(rc(substr(ref, pos + ins - 100, 100)))
        printf "@p%d/1\n%s\n+\n%s\n", k, r1, q > "r1.fq"
        printf "@p%d/2\n%s\n+\n%s\n", k, r2, q > "r2.fq"
    }
}' ref.fa
"$BWA_MEM3_TESTING" index ref.fa 2> index.err || fail "index" index.err

K=200000
run() { # <tag> <env...> -- <mem args...>
    local tag=$1
    shift
    local envs=()
    while [[ $1 != -- ]]; do
        envs+=("$1")
        shift
    done
    shift
    # Defaults first, so an exported streaming switch or test hook cannot leak in;
    # each run's own settings (the =0 baseline included) come after and win.
    env BWA_MEM3_STREAM_STATS=1 BWA_MEM3_STREAM_OUTPUT=1 BWAMEM3_TESTING_WRITE_DELAY_MS=0 \
        BWAMEM3_TESTING_HOLD_LAST_ITEM=0 "${envs[@]}" "$BWA_MEM3_TESTING" mem -K "$K" -o "$tag.out" ref.fa "$@" 2> "$tag.err" \
        || fail "mem ($tag)" "$tag.err"
    grep -v '^@PG' "$tag.out" > "$tag.sam" || fail "empty output ($tag)"
}
batches() { grep -c 'stream_output:' "$1.err" || true; }
early_batches() { awk '/stream_output:/ && $3 > 0 { n++ } END { print n + 0 }' "$1.err"; }

for end in pe se; do
    if [[ $end == pe ]]; then in=(r1.fq r2.fq); else in=(r1.fq); fi
    run "$end.off" BWA_MEM3_STREAM_OUTPUT=0 -- -t 4 "${in[@]}"
    [[ $(wc -l < "$end.off.sam") -gt 8000 ]] || fail "$end: too few records"
    for t in 2 3 4; do
        tag="$end.delay.t$t"
        run "$tag" BWAMEM3_TESTING_WRITE_DELAY_MS=100 -- -t "$t" "${in[@]}"
        if ! cmp -s "$end.off.sam" "$tag.sam"; then
            diff "$end.off.sam" "$tag.sam" | head -20 || true
            fail "$tag: streamed output differs from unstreamed with the write step delayed"
        fi
        n=$(batches "$tag")
        [[ $n -ge 4 ]] || fail "$tag: only $n batches; expected 8 paired, 4 single-end"
    done
done
run pe.nodelay BWAMEM3_TESTING_WRITE_DELAY_MS=0 BWAMEM3_TESTING_HOLD_LAST_ITEM=1 -- -t 4 r1.fq r2.fq
cmp -s pe.off.sam pe.nodelay.sam || fail "pe.nodelay: streamed output differs from unstreamed"
n=$(batches pe.nodelay)
e=$(early_batches pe.nodelay)
[[ $n -ge 4 && $e -eq $n ]] || fail "pe.nodelay: $e of $n batches wrote early with the last item held; streaming is not engaging"
# (3) --bam=0, paired. One fixed -o name for every run so the headers match.
# K divides the 8000 pairs exactly, and an exactly full last batch is not
# streamed (see mem.md), so cut 7 batches with a short last one instead.
KB=230000
bam_run() { # <tag> <env...>: --bam=0 -t 4 paired into bam.cur, moved to <tag>.bam
    local tag=$1
    shift
    env BWA_MEM3_STREAM_STATS=1 BWA_MEM3_STREAM_OUTPUT=1 BWAMEM3_TESTING_WRITE_DELAY_MS=0 \
        BWAMEM3_TESTING_HOLD_LAST_ITEM=0 "$@" "$BWA_MEM3_TESTING" mem -K "$KB" -t 4 --bam=0 -o bam.cur ref.fa r1.fq r2.fq 2> "$tag.err" \
        || fail "mem ($tag)" "$tag.err"
    [[ -s bam.cur ]] || fail "empty output ($tag)"
    mv bam.cur "$tag.bam"
}
bam_run bam0.off BWA_MEM3_STREAM_OUTPUT=0
[[ $(batches bam0.off) -eq 0 ]] || fail "bam0.off: streamed with BWA_MEM3_STREAM_OUTPUT=0"
bam_run bam0.delay BWAMEM3_TESTING_WRITE_DELAY_MS=100
cmp -s bam0.off.bam bam0.delay.bam || fail "bam0.delay: streamed BAM differs from unstreamed with the write step delayed"
[[ $(batches bam0.delay) -eq 1 ]] || fail "bam0.delay: $(batches bam0.delay) batches armed to stream at --bam=0, expected only the last"
bam_run bam0.hold BWAMEM3_TESTING_HOLD_LAST_ITEM=1
cmp -s bam0.off.bam bam0.hold.bam || fail "bam0.hold: streamed BAM differs from unstreamed"
[[ $(batches bam0.hold) -eq 1 && $(early_batches bam0.hold) -eq 1 ]] \
    || fail "bam0.hold: the last batch did not write early with its last item held ($(early_batches bam0.hold) of $(batches bam0.hold))"
echo "PASS: stream_output_order ($(batches pe.delay.t4) batches; output is unchanged with the write step delayed; all $e batches write early without the delay; --bam=0 last batch unchanged and streamed)"
