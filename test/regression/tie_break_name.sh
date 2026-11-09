#!/usr/bin/env bash
# test/regression/tie_break_name.sh
#
# End-to-end behavior of --tie-break.
#
# When several alignments of a read (or pairings of a pair) score the same,
# bwa, bwa-mem2 and bwa-mem3 pick one with a hash of the read's position in the
# input. The position comes from a counter over the whole input, so -K and -t do
# not move it, but the same read can land on a different copy of a repeat when
# the input is reordered, subset or sharded. --tie-break name seeds the choice
# from a hash of the read name instead.
#
# This test pins the contract, on a generated fixture whose reads include exact
# repeats (MAPQ-0 ties):
#   * default == --tie-break index, byte for byte (the upstream-compatible path),
#     and the default is unchanged by -K/-t (its ids come from the global counter)
#   * --tie-break name gives the same records (after sorting) for the same reads
#     in a shuffled order, at a different -K and -t, and with every cohort taken
#     through the sliced align path, paired-end with -I fixed and single-end
#   * name mode keeps every read: same primary-record keys as the default
#   * the seed is the read NAME: renaming the reads (same order and sequences)
#     changes the tie-broken placements
#   * the default does NOT survive the shuffle on this fixture, and name mode
#     differs from the default, so none of the checks above is vacuous
#   * an unknown value is rejected, and name mode reports itself at startup
#
# Inputs:
#   BWA_MEM3 — path to the bwa-mem3 binary under test
set -euo pipefail

# Checked explicitly rather than with ${VAR:?}, which exits with a bare shell
# error (and before any ERR trap runs), so the failure would carry no FAIL: line.
if [[ -z ${BWA_MEM3:-} ]]; then
    echo "FAIL: BWA_MEM3 must be set" >&2
    exit 1
fi
# python3 generates the fixture. A missing python3 is a failure, not a skip: the
# CI step does not otherwise check for it, so a skip would turn it green untested.
command -v python3 > /dev/null 2>&1 || {
    echo "FAIL: python3 not on PATH (it generates the fixture)" >&2
    exit 1
}

case "$BWA_MEM3" in
    /*) BIN="$BWA_MEM3" ;;
    *) BIN="$PWD/$BWA_MEM3" ;;
esac
[[ -x "$BIN" ]] || {
    echo "FAIL: BWA_MEM3 ($BIN) is not executable" >&2
    exit 1
}

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
fail() {
    echo "FAIL: $*" >&2
    exit 1
}

# A 40 kb random reference with a 3 kb block copied exactly to two more places,
# so a pair drawn from inside the block aligns equally well at all three copies.
# Pairs (2x100, ~300 bp inserts) are drawn across the whole reference, two thirds
# of them from the block. The shuffled copy holds the same pairs in another order;
# the renamed copy holds the same pairs in the same order under other names.
python3 - << 'PY' || fail "fixture generation failed"
import random
rnd = random.Random(20261108)
L, BLOCK, LEN = 40000, 3000, 100
ref = [rnd.choice('ACGT') for _ in range(L)]
for at in (17000, 31000):
    ref[at:at + BLOCK] = ref[5000:5000 + BLOCK]
ref = ''.join(ref)
with open('ref.fa', 'w') as f:
    f.write('>ref\n')
    for i in range(0, L, 60):
        f.write(ref[i:i + 60] + '\n')
comp = str.maketrans('ACGT', 'TGCA')
pairs = []
for n in range(3000):
    insert = rnd.randint(260, 340)
    if n % 3:
        start = rnd.randint(5000, 5000 + BLOCK - insert)
    else:
        start = rnd.randint(0, L - insert)
    frag = ref[start:start + insert]
    r1, r2 = frag[:LEN], frag[-LEN:].translate(comp)[::-1]
    if rnd.random() < 0.5:
        r1, r2 = r2, r1
    name = 'p%05d' % n
    pairs.append((f'@{name}/1\n{r1}\n+\n{"I" * LEN}\n', f'@{name}/2\n{r2}\n+\n{"I" * LEN}\n'))
def write(tag, ps):
    with open(f'{tag}.R1.fq', 'w') as a, open(f'{tag}.R2.fq', 'w') as b:
        for x, y in ps:
            a.write(x)
            b.write(y)
write('orig', pairs)
write('ren', [(x.replace('@p', '@q', 1), y.replace('@p', '@q', 1)) for x, y in pairs])
rnd.shuffle(pairs)
write('shuf', pairs)
PY

"$BIN" index ref.fa > /dev/null 2>&1 || fail "index ref.fa"

# Pinned insert-size distribution: without it, paired-end output depends on which
# pairs share a batch (mem_pestat), whatever the tie-break.
INS=(-I "300,20,400,200")

# Run `bwa-mem3 mem` with the given args; stdout to $1.sam, stderr to $1.log.
run() {
    local out="$1"
    shift
    "$BIN" mem "$@" > "$out.sam" 2> "$out.log" || {
        cat "$out.log" >&2
        fail "mem $* exited nonzero"
    }
}
# Records only, sorted: the order-independent view of a run.
records() { grep -v '^@' "$1.sam" | LC_ALL=C sort; }
# Records without QNAME, sorted: compares placements across differently named runs.
placements() { grep -v '^@' "$1.sam" | cut -f2- | LC_ALL=C sort; }
# The primary-record keys (QNAME and the read-1/read-2 bits), sorted. Bit tests
# are arithmetic because not every awk implementation has and().
primary_keys() {
    awk -F'\t' '!/^@/ && int($2 / 256) % 2 == 0 && int($2 / 2048) % 2 == 0 { print $1 "\t" int($2 / 64) % 4 }' "$1.sam" \
        | LC_ALL=C sort
}
# Everything but @PG (its CL: differs between runs by construction).
no_pg() { grep -v '^@PG' "$1.sam"; }
# Assert that run $1 has $3 primary records and run $2 has the same primary keys.
same_reads() {
    local n
    n=$(primary_keys "$1" | wc -l | tr -d ' ')
    [[ $n -eq $3 ]] || fail "$1 has $n primary records, expected $3"
    cmp -s <(primary_keys "$1") <(primary_keys "$2") || fail "$2 does not keep the same reads as $1"
}

# --- default == --tie-break index, byte for byte -----------------------------
run def -t 2 "${INS[@]}" ref.fa orig.R1.fq orig.R2.fq
run idx -t 2 "${INS[@]}" --tie-break index ref.fa orig.R1.fq orig.R2.fq
cmp -s <(no_pg def) <(no_pg idx) || fail "--tie-break index changed the output"
# The default's ids come from the global read counter, so -K/-t must not move it;
# a call site that lost the global offset (a batch-local id) would fail here.
run def_k -t 1 -K 20000 "${INS[@]}" ref.fa orig.R1.fq orig.R2.fq
cmp -s <(records def) <(records def_k) || fail "paired-end default changed with -K/-t"
echo "PASS: --tie-break index is byte-identical to the default, which is unchanged by -K/-t"

# --- the fixture has ties: the default moves reads when the input is shuffled -
run def_shuf -t 2 "${INS[@]}" ref.fa shuf.R1.fq shuf.R2.fq
if cmp -s <(records def) <(records def_shuf); then
    fail "default output survived the shuffle; the fixture no longer exercises ties"
fi
n_mapq0=$(awk -F'\t' '!/^@/ && $5 == 0 { n++ } END { print n + 0 }' def.sam)
[[ $n_mapq0 -gt 0 ]] || fail "fixture produced no MAPQ-0 records"
echo "PASS: the default moves MAPQ-0 reads under a shuffle ($n_mapq0 MAPQ-0 records)"

# --- --tie-break name: order, batching and threads do not matter --------------
run name -t 2 "${INS[@]}" --tie-break name ref.fa orig.R1.fq orig.R2.fq
grep -q 'tie-break: name' name.log || fail "--tie-break name did not report itself at startup"
if grep -q 'tie-break: name' def.log; then
    fail "the default reported --tie-break name at startup"
fi
same_reads def name 6000
run name_shuf -t 2 "${INS[@]}" --tie-break name ref.fa shuf.R1.fq shuf.R2.fq
cmp -s <(records name) <(records name_shuf) || fail "paired-end --tie-break name changed under a shuffle"
# -K 20000 splits the 3000 pairs into many batches; -t 1 changes the work split.
run name_k -t 1 -K 20000 "${INS[@]}" --tie-break name ref.fa shuf.R1.fq shuf.R2.fq
cmp -s <(records name) <(records name_k) || fail "paired-end --tie-break name changed with -K/-t"
# A -K this small never reaches the sliced align path (slices have a 1 Mbase
# floor), so take every cohort through it explicitly.
BWA_MEM3_COHORT_SLICE_ALL=1 run name_slice -t 2 -K 20000 "${INS[@]}" --tie-break name ref.fa shuf.R1.fq shuf.R2.fq
grep -q 'BWA_MEM3_COHORT_SLICE_ALL set' name_slice.log || fail "BWA_MEM3_COHORT_SLICE_ALL was not honored"
cmp -s <(records name) <(records name_slice) || fail "paired-end --tie-break name changed on the sliced align path"
echo "PASS: paired-end --tie-break name is independent of input order, -K, -t and cohort slicing"

run se_def -t 2 ref.fa orig.R1.fq
run se_def_k -t 1 -K 20000 ref.fa orig.R1.fq
cmp -s <(records se_def) <(records se_def_k) || fail "single-end default changed with -K/-t"
run se_name -t 2 --tie-break name ref.fa orig.R1.fq
same_reads se_def se_name 3000
run se_name_shuf -t 1 -K 20000 --tie-break name ref.fa shuf.R1.fq
cmp -s <(records se_name) <(records se_name_shuf) || fail "single-end --tie-break name changed under a shuffle"
run se_def_shuf -t 2 ref.fa shuf.R1.fq
if cmp -s <(records se_def) <(records se_def_shuf); then
    fail "single-end default output survived the shuffle; the fixture no longer exercises ties"
fi
echo "PASS: single-end --tie-break name is independent of input order, -K and -t"

# --- name mode really changes which copy wins ---------------------------------
# Otherwise "independent of order" could be met by a no-op that happened to
# reproduce the default's placement for this input order.
if cmp -s <(records def) <(records name); then
    fail "--tie-break name produced the default's records exactly"
fi
echo "PASS: --tie-break name changes tie-broken placements relative to the default"

# --- the seed is the read name -----------------------------------------------
# Same reads, same order, other names: a seed taken from anything but the name
# (the sequence, say) would also be order-independent, so check the name matters.
run name_ren -t 2 "${INS[@]}" --tie-break name ref.fa ren.R1.fq ren.R2.fq
if cmp -s <(placements name) <(placements name_ren); then
    fail "paired-end --tie-break name placements did not change when the reads were renamed"
fi
run se_name_ren -t 2 --tie-break name ref.fa ren.R1.fq
if cmp -s <(placements se_name) <(placements se_name_ren); then
    fail "single-end --tie-break name placements did not change when the reads were renamed"
fi
echo "PASS: --tie-break name placements follow the read names"

# --- bad values are rejected ---------------------------------------------------
rc=0
"$BIN" mem --tie-break bogus ref.fa orig.R1.fq > bad.sam 2> bad.log || rc=$?
[[ $rc -ne 0 ]] || fail "--tie-break bogus was accepted"
# Rejected at option parsing: no SAM output and no reference load.
[[ ! -s bad.sam ]] || fail "--tie-break bogus produced SAM output"
if grep -qF '* Ref file:' bad.log; then
    fail "--tie-break bogus reached reference loading"
fi
grep -q "unknown --tie-break 'bogus'" bad.log || {
    cat bad.log >&2
    fail "--tie-break bogus did not name the bad value"
}
echo "PASS: an unknown --tie-break value is rejected"
