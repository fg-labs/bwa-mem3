#!/usr/bin/env bash
# test/regression/chain_flat_replay_identity.sh
#
# Regression: seed chaining produces identical alignment records whether a read
# is chained on the flat sorted-key index (the default) or replayed through the
# kbtree (BWA3_CHAIN_FLAT_CAP), and the Pass-3 kept-set index reaches its
# rarely-taken branches on the same fixture.
#
# The flat index hands a read to the unchanged kbtree path when two chains would
# share a position (the kbtree's answer then depends on its node layout) or when
# the read has more chains than the cap, rewinding the seed buffer and freeing
# any seed arrays it heap-grew first. Neither happens often on ordinary reads, so
# this builds a fixture that forces each one and proves it did:
#
#   tie reads   A + F + A, where A is one 40-mer: A's two seeds hit the same
#               reference position 120 query bases apart -- beyond the band, so
#               they cannot join one chain -- and the second opens a chain at an
#               existing chain's position. BWA3_CHAIN_STATS must report
#               flat_declined_tie > 0 on the default run.
#   many-chain  six 30-mers from distant loci: six chains, so BWA3_CHAIN_FLAT_CAP=3
#   reads       declines each read mid-way, after the flat pass consumed seed-buffer
#               slots (flat_declined_cap > 0), exercising the rewind and replay.
#   long reads  36 kb with a substitution every 400 bases: one alignment spans the
#               read, so the Pass-3 index answers later seeds by scanning every
#               member rather than walking hash buckets (p3_linear_scans > 0).
#
# Each fixture run is repeated with BWA3_CHAIN_FLAT_CAP=0 (every read on the
# kbtree from its first seed) and =3, on the default seed order, on the opt-in
# reorder path (--seed-order local-longest, which replays from its own seed
# buffer), and under --meth; each must match its default run byte for byte. A
# malformed cap must be reported, which proves the knob is read at all.
#
# When CHR22_FA / CHR22_SIM_DIR are set, the same cap 0 / cap 3 identity is also
# checked on the chr22 holodeck reads. CI runs the fixture-only mode (unset
# CHR22_*) again under the ASAN build and the opt-in debug-macro build, where
# BWA_MEM3_DEBUG_P3_XCHECK / _UNGAPPED_XCHECK cross-check the fast paths on it.
#
# Inputs:
#   BWA_MEM3                  — the bwa-mem3 binary under test
#   CHR22_FA, CHR22_SIM_DIR   — optional: pre-indexed chr22.fa and the directory
#                               holding holodeck reads.r[12].fastq.gz
set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
case "$BWA_MEM3" in
    */*) BWA_MEM3="$(cd "$(dirname "$BWA_MEM3")" && pwd)/$(basename "$BWA_MEM3")" ;;
esac
for tool in mawk samtools; do
    command -v "$tool" > /dev/null 2>&1 || {
        echo "FAIL: $tool not on PATH (needed to build the fixture / read --meth BAM)"
        exit 1
    }
done

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
fail() {
    echo "FAIL: $*"
    exit 1
}

# align <out> <ref> <reads...> -- [env...] [-- flags...]: SAM minus @PG, stderr kept
align() {
    local out=$1 ref=$2
    shift 2
    local reads=()
    while [ $# -gt 0 ] && [ "$1" != -- ]; do
        reads+=("$1")
        shift
    done
    shift
    local envs=()
    while [ $# -gt 0 ] && [ "$1" != -- ]; do
        envs+=("$1")
        shift
    done
    [ $# -gt 0 ] && shift
    local rc=0
    env -u BWA3_CHAIN_FLAT_CAP -u BWA3_CHAIN_STATS ${envs[@]+"${envs[@]}"} \
        "$BWA_MEM3" mem -t 4 "$@" "$ref" "${reads[@]}" > "$out.raw" 2> "$out.err" || rc=$?
    if [ "$rc" -ne 0 ]; then
        tail -5 "$out.err"
        fail "bwa-mem3 mem exited $rc ($(basename "$out"))"
    fi
    samtools view -h "$out.raw" | grep -v '^@PG' > "$out.sam" \
        || fail "could not read the output of $(basename "$out")"
    [ -s "$out.sam" ] || fail "empty output from $(basename "$out")"
}
# Both helpers run under `set -o pipefail`, where a pipeline that finds a mismatch
# or no match would end the script before its FAIL: line; `|| true` keeps the
# pipeline's status from doing that, and the checks below decide instead.
same() { # same <label> <a.sam> <b.sam>
    cmp -s "$2" "$3" || {
        diff "$2" "$3" | head -20 || true
        fail "$1"
    }
}
stat_of() { # stat_of <key> <stderr file>: the value, or empty when absent
    { grep -o "$1=[0-9]*" "$2" || true; } | head -1 | cut -d= -f2
}

# --- the generated fixture ---------------------------------------------------
cd "$WORK"
mawk 'BEGIN {
    srand(7); b = "ACGT"; L = 300000; s = ""
    for (i = 0; i < L; i++) s = s substr(b, int(rand() * 4) + 1, 1)
    print ">g"
    for (i = 1; i <= L; i += 80) print substr(s, i, 80)
    q = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
    n = 0
    for (k = 0; k < 200; k++) {            # tie reads: A + F + A
        p = int(rand() * (L - 50000)) + 1; f = p + 20000
        r = substr(s, p, 40) substr(s, f, 80) substr(s, p, 40)
        printf "@tie%d\n%s\n+\n%s\n", k, r, substr(q q, 1, length(r)) > "reads.fq"
    }
    for (k = 0; k < 200; k++) {            # many-chain reads: six distant 30-mers
        r = ""
        for (j = 0; j < 6; j++) r = r substr(s, int(rand() * (L - 100)) + 1, 30)
        printf "@many%d\n%s\n+\n%s\n", k, r, substr(q q q, 1, length(r)) > "reads.fq"
    }
    for (k = 0; k < 4; k++) {              # long reads: 36 kb, a substitution every 400 bp
        p = int(rand() * (L - 40000)) + 1; r = substr(s, p, 36000); t = ""
        for (i = 1; i <= 36000; i += 400) {
            c = substr(r, i, 1); m = (c == "A") ? "C" : "A"
            t = t m substr(r, i + 1, 399)
        }
        qq = ""; while (length(qq) < 36000) qq = qq q
        printf "@long%d\n%s\n+\n%s\n", k, t, substr(qq, 1, 36000) > "reads.fq"
    }
}' > ref.fa 2> /dev/null
"$BWA_MEM3" index ref.fa > index.log 2>&1 || fail "bwa-mem3 index exited nonzero"
"$BWA_MEM3" index --meth ref.fa > index-meth.log 2>&1 || fail "bwa-mem3 index --meth exited nonzero"

# The default run carries the stats: each forced path must actually have run.
align def ref.fa reads.fq -- BWA3_CHAIN_STATS=1
ties=$(stat_of flat_declined_tie def.err)
linear=$(stat_of p3_linear_scans def.err)
flat=$(stat_of flat_reads def.err)
[ -n "$ties" ] || fail "no [chain-stats] line: is BWA3_CHAIN_STATS still read?"
[ "$ties" -gt 0 ] || fail "no read was handed to the kbtree on a tie (flat_declined_tie=0)"
[ "$linear" -gt 0 ] || fail "the Pass-3 index never took its full scan (p3_linear_scans=0)"
[ -n "$flat" ] || fail "no flat_reads count in the [chain-stats] line"
[ "$flat" -gt 0 ] || fail "no read completed on the flat index (flat_reads=0)"
align cap0 ref.fa reads.fq -- BWA3_CHAIN_FLAT_CAP=0 BWA3_CHAIN_STATS=1
[ "$(stat_of flat_reads cap0.err)" = 0 ] || fail "BWA3_CHAIN_FLAT_CAP=0 left reads on the flat index"
align cap3 ref.fa reads.fq -- BWA3_CHAIN_FLAT_CAP=3 BWA3_CHAIN_STATS=1
capped=$(stat_of flat_declined_cap cap3.err)
[ "$capped" -gt 0 ] || fail "BWA3_CHAIN_FLAT_CAP=3 declined no read mid-way (flat_declined_cap=0)"
align bad ref.fa reads.fq -- BWA3_CHAIN_FLAT_CAP=not-a-number
grep -q 'ERROR: BWA3_CHAIN_FLAT_CAP="not-a-number" is not a non-negative integer' bad.err \
    || fail "a malformed BWA3_CHAIN_FLAT_CAP was not reported; is the knob still read?"
same "flat index != kbtree (cap 0)" def.sam cap0.sam
same "flat index != mid-read replay (cap 3)" def.sam cap3.sam
same "malformed cap changed the output" def.sam bad.sam

align ord ref.fa reads.fq -- -- --seed-order local-longest
align ord3 ref.fa reads.fq -- BWA3_CHAIN_FLAT_CAP=3 -- --seed-order local-longest
align ord0 ref.fa reads.fq -- BWA3_CHAIN_FLAT_CAP=0 -- --seed-order local-longest
same "reorder path: flat index != mid-read replay (cap 3)" ord.sam ord3.sam
same "reorder path: flat index != kbtree (cap 0)" ord.sam ord0.sam

align meth ref.fa reads.fq -- -- --meth
align meth3 ref.fa reads.fq -- BWA3_CHAIN_FLAT_CAP=3 -- --meth
align meth0 ref.fa reads.fq -- BWA3_CHAIN_FLAT_CAP=0 -- --meth
same "--meth: flat index != mid-read replay (cap 3)" meth.sam meth3.sam
same "--meth: flat index != kbtree (cap 0)" meth.sam meth0.sam
fixture="fixture: $(grep -cv '^@' def.sam) records, $ties tie and $capped cap declines, $linear full Pass-3 scans"

# --- chr22 holodeck reads (optional) -----------------------------------------
if [ -n "${CHR22_FA:-}" ] && [ -n "${CHR22_SIM_DIR:-}" ]; then
    R1="$CHR22_SIM_DIR/reads.r1.fastq.gz"
    R2="$CHR22_SIM_DIR/reads.r2.fastq.gz"
    align c22 "$CHR22_FA" "$R1" "$R2" --
    align c22cap0 "$CHR22_FA" "$R1" "$R2" -- BWA3_CHAIN_FLAT_CAP=0
    align c22cap3 "$CHR22_FA" "$R1" "$R2" -- BWA3_CHAIN_FLAT_CAP=3
    same "chr22: flat index != kbtree (cap 0)" c22.sam c22cap0.sam
    same "chr22: flat index != mid-read replay (cap 3)" c22.sam c22cap3.sam
    echo "PASS: flat chaining index == kbtree replay ($fixture; chr22: $(grep -cv '^@' c22.sam) records)"
else
    echo "PASS: flat chaining index == kbtree replay ($fixture; chr22 not configured)"
fi
