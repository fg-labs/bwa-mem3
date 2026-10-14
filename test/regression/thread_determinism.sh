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
