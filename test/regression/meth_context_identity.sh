#!/usr/bin/env bash
# test/regression/meth_context_identity.sh
#
# Regression: the identity and invariant gates for `--meth-context`. On inputs
# generated in this script (no committed data), all four hold:
#
#  (a) `--meth=5base --meth-context all` == `--meth=5base`, record for record.
#  (b) `--meth=5base --meth-context cpg` is identical across every
#      BWAMEM3_FORCE_TIER the host can run among avx2 and avx512bw (neon is the
#      one tier on arm64). Pre-AVX2 tiers are retired at runtime and are not
#      swept. With one runnable tier the cross-tier part reports SKIP:.
#  (c) the `cpg` records are identical with BWAMEM3_METH_BATCHED_RESCUE on
#      (batched kswv rescue) and off (scalar rescue). Same for plain 5base. On
#      pe5, turning mate rescue off (-S) must change the output, so the
#      comparison is known to cover rescued mates.
#  (d) without `--meth`, output equals the base commit's binary, passed as
#      BWA_MEM3_BASE. Reported SKIP: when BWA_MEM3_BASE is unset (CI has none).
#
# Every run of every leg also has to satisfy two invariants: no record has
# MAPQ > 60, and every mapped record carries AS:i. (A pair given a negative
# score once surfaced as MAPQ 191 on EM-like reads, which is why the EM-like set
# below is part of the fixture.) The legs must map most reads and, for the PE
# sets, produce proper pairs, so a run that silently maps nothing cannot pass.
#
# Fixture (python3, fixed seeds): a 120 kb two-contig reference, and three read
# sets cut from it:
#   pe5 -- 2x150 directional 5-base-like pairs, ~85% of CpG Cs (OT) or Gs (OB)
#          converted, 1.5% of the other Cs converted, 0.3% substitution errors,
#          one pair in 12 with a mate made unseedable by spaced mismatches so it
#          is placed by mate rescue.
#   em  -- 2x150 EM-like pairs, ~97% of every C (OT) or G (OB) converted.
#   se  -- 20-40 bp single reads, converted as pe5, either strand.
#
# Inputs:
#   BWA_MEM3      -- the bwa-mem3 binary under test
#   BWA_MEM3_BASE -- optional: the base commit's binary, for (d)
set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"
abspath() { printf '%s' "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
case "$BWA_MEM3" in
    */*) BIN="$(abspath "$BWA_MEM3")" ;;
    *) BIN="$(command -v "$BWA_MEM3")" ;;
esac
BASE_BIN=""
if [ -n "${BWA_MEM3_BASE:-}" ]; then
    case "$BWA_MEM3_BASE" in
        */*) BASE_BIN="$(abspath "$BWA_MEM3_BASE")" ;;
        *) BASE_BIN="$(command -v "$BWA_MEM3_BASE")" ;;
    esac
fi
for tool in samtools mawk python3; do
    command -v "$tool" > /dev/null 2>&1 || {
        echo "SKIP: $tool not on PATH"
        exit 0
    }
done

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
fail() {
    echo "FAIL: meth_context_identity: $*" >&2
    exit 1
}
[ -x "$BIN" ] || fail "BWA_MEM3 ($BIN) is not executable"
if [ -n "$BASE_BIN" ] && [ ! -x "$BASE_BIN" ]; then
    fail "BWA_MEM3_BASE ($BASE_BIN) is not executable"
fi

python3 -I - << 'PY'
import random

COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def rc(s):
    return ''.join(COMP[b] for b in reversed(s))

rng = random.Random(20261028)
def rand_seq(n):
    return ''.join(rng.choice('ACGT') for _ in range(n))

contigs = [rand_seq(70000), rand_seq(50000)]
with open('ref.fa', 'w') as f:
    for i, c in enumerate(contigs):
        f.write('>chr%d\n%s\n' % (i + 1, c))

def convert(frag, strand, p_cpg, p_other):
    """Convert top-strand fragment frag: strand 'T' converts C (OT), 'B' converts G (OB)."""
    out = list(frag)
    n = len(frag)
    for i, b in enumerate(frag):
        if strand == 'T' and b == 'C':
            cpg = i + 1 < n and frag[i + 1] == 'G'
            if rng.random() < (p_cpg if cpg else p_other):
                out[i] = 'T'
        elif strand == 'B' and b == 'G':
            cpg = i > 0 and frag[i - 1] == 'C'
            if rng.random() < (p_cpg if cpg else p_other):
                out[i] = 'A'
    return ''.join(out)

def err(seq, rate):
    out = list(seq)
    for i, b in enumerate(out):
        if rng.random() < rate:
            out[i] = rng.choice([x for x in 'ACGT' if x != b])
    return ''.join(out)

def unseedable(seq):
    # A mismatch every 14 bases leaves no 19-mer, so the mate is placed by rescue.
    out = list(seq)
    for i in range(7, len(out), 14):
        out[i] = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}[out[i]]
    return ''.join(out)

def pick(minlen, maxlen):
    c = rng.choice(contigs)
    n = rng.randint(minlen, maxlen)
    p = rng.randint(0, len(c) - n)
    return c[p:p + n]

def pairs(path, n, p_cpg, p_other, rescue_every):
    q = 'I' * 150
    r1, r2 = [], []
    for k in range(n):
        frag = pick(250, 400)
        strand = rng.choice('TB')
        cv = convert(frag, strand, p_cpg, p_other)
        a = err(cv[:150], 0.003)
        b = err(cv[-150:], 0.003)
        if rescue_every and k % rescue_every == 0:
            b = unseedable(b)
        # FR pair: read 1 is the fragment's left end for OT, its right end (rc) for OB.
        if strand == 'T':
            x, y = a, rc(b)
        else:
            x, y = rc(b), a
        r1.append('@p%d\n%s\n+\n%s\n' % (k, x, q))
        r2.append('@p%d\n%s\n+\n%s\n' % (k, y, q))
    with open(path + '_1.fq', 'w') as f: f.write(''.join(r1))
    with open(path + '_2.fq', 'w') as f: f.write(''.join(r2))

pairs('pe5', 2500, 0.85, 0.015, 12)
pairs('em', 1500, 0.97, 0.97, 0)

se = []
for k in range(3000):
    frag = pick(20, 40)
    strand = rng.choice('TB')
    s = err(convert(frag, strand, 0.85, 0.015), 0.005)
    if rng.random() < 0.5:
        s = rc(s)
    se.append('@s%d\n%s\n+\n%s\n' % (k, s, 'I' * len(s)))
with open('se.fq', 'w') as f: f.write(''.join(se))
PY
for f in ref.fa pe5_1.fq pe5_2.fq em_1.fq em_2.fq se.fq; do
    [ -s "$f" ] || fail "fixture generation failed ($f)"
done

"$BIN" index --meth ref.fa > /dev/null 2>&1 || fail "index --meth exit nonzero"

# Index for (d): plain, built by each binary, so neither depends on the other's format.
if [ -n "$BASE_BIN" ]; then
    mkdir plain_new plain_base
    cp ref.fa plain_new/ref.fa
    cp ref.fa plain_base/ref.fa
    "$BIN" index plain_new/ref.fa > /dev/null 2>&1 || fail "index (new) exit nonzero"
    "$BASE_BIN" index plain_base/ref.fa > /dev/null 2>&1 || fail "index (base) exit nonzero"
fi

# Tiers (b) can sweep: avx2 and avx512bw up to the host's, or the one native tier.
HOST_RAW="$(BWAMEM3_DEBUG_SIMD=1 "$BIN" 2>&1 || true)"
HOST_TIER="$(printf '%s\n' "$HOST_RAW" | sed -n 's/.*SIMD tier: \([a-z0-9]*\).*/\1/p' | head -1)"
case "$HOST_TIER" in
    neon) TIERS=(neon) ;;
    avx512bw) TIERS=(avx2 avx512bw) ;;
    avx2) TIERS=(avx2) ;;
    sse41 | sse42 | avx) TIERS=() ;;
    *) fail "could not detect the host SIMD tier (parsed '$HOST_TIER')" ;;
esac
echo "Host tier: $HOST_TIER; sweeping: ${TIERS[*]:-none}"

# run LABEL TIER|- ENV=VAL|- BINARY INDEX_REF SET MEM_ARGS... -> LABEL.sam (records only)
# SET: pe5 | em | se. The SE set takes shorter seeds and a lower minimum score.
run() {
    local out=$1 tier=$2 envkv=$3 bin=$4 ref=$5 set=$6
    shift 6
    local -a inputs=("${set}_1.fq" "${set}_2.fq") extra=()
    if [ "$set" = se ]; then
        inputs=(se.fq)
        extra=(-k 15 -T 15)
    fi
    local -a pre=(env)
    [ "$tier" != "-" ] && pre+=("BWAMEM3_FORCE_TIER=$tier")
    [ "$envkv" != "-" ] && pre+=("$envkv")
    "${pre[@]}" "$bin" mem -t 1 ${extra[@]+"${extra[@]}"} "$@" "$ref" "${inputs[@]}" \
        2> "$out.log" | samtools view - > "$out.sam" \
        || fail "$out: mem exited nonzero ($(tail -3 "$out.log" | tr '\n' ' '))"
    [ -s "$out.sam" ] || fail "$out: empty output"
    check_invariants "$out" "$set"
}

# Invariants every run must satisfy, whatever it is compared against.
check_invariants() {
    local out=$1 set=$2 stats
    stats="$(mawk -F'\t' '
        { n++ }
        int($2 / 4) % 2 == 0 { mapped++; if ($0 !~ /\tAS:i:/) noas++ }
        $5 > 60 { hi++ }
        int($2 / 2) % 2 == 1 { proper++ }
        END { printf "%d %d %d %d %d", n, mapped, noas, hi, proper }' "$out.sam")"
    local n mapped noas hi proper
    read -r n mapped noas hi proper <<< "$stats"
    [ "$hi" -eq 0 ] || fail "$out: $hi record(s) with MAPQ > 60"
    [ "$noas" -eq 0 ] || fail "$out: $noas mapped record(s) without AS:i"
    # Under cpg the EM-like set pays for every non-CpG conversion, so a good part
    # of it legitimately fails to map; it only has to map something.
    local floor=$((n * 7 / 10))
    [ "$set" = em ] && floor=$((n / 10))
    [ "$mapped" -ge "$floor" ] || fail "$out: only $mapped of $n records mapped"
    if [ "$set" != se ]; then
        local pfloor=$((n / 2))
        [ "$set" = em ] && pfloor=$((n / 10))
        [ "$proper" -ge "$pfloor" ] || fail "$out: only $proper of $n records properly paired"
    fi
}

compare() { # compare WHAT A B
    if ! cmp -s "$2.sam" "$3.sam"; then
        diff "$2.sam" "$3.sam" | head -10 >&2 || true
        fail "$1: $2 differs from $3"
    fi
}

NCMP=0
for set in pe5 em se; do
    # (a) context all == default
    run "a_${set}_default" - - "$BIN" ref.fa "$set" --meth=5base
    run "a_${set}_all" - - "$BIN" ref.fa "$set" --meth=5base --meth-context all
    compare "(a) $set" "a_${set}_default" "a_${set}_all"
    NCMP=$((NCMP + 1))

    # (b) cpg across tiers
    run "b_${set}_cpg_native" - - "$BIN" ref.fa "$set" --meth=5base --meth-context cpg
    # cpg must actually change something relative to all, or the set tests nothing
    if cmp -s "b_${set}_cpg_native.sam" "a_${set}_all.sam"; then
        fail "$set: cpg output equals all output; the fixture exercises no cpg difference"
    fi
    if [ "${#TIERS[@]}" -gt 1 ]; then
        for t in "${TIERS[@]}"; do
            run "b_${set}_cpg_$t" "$t" - "$BIN" ref.fa "$set" --meth=5base --meth-context cpg
            compare "(b) $set" "b_${set}_cpg_native" "b_${set}_cpg_$t"
            NCMP=$((NCMP + 1))
        done
    fi

    # (c) batched vs scalar rescue
    for ctx in "" cpg; do
        tag=${ctx:-5base}
        args=(--meth=5base)
        [ -n "$ctx" ] && args+=(--meth-context "$ctx")
        run "c_${set}_${tag}_scalar" - BWAMEM3_METH_BATCHED_RESCUE=0 "$BIN" ref.fa "$set" "${args[@]}"
        run "c_${set}_${tag}_batched" - BWAMEM3_METH_BATCHED_RESCUE=1 "$BIN" ref.fa "$set" "${args[@]}"
        compare "(c) $set $tag" "c_${set}_${tag}_scalar" "c_${set}_${tag}_batched"
        NCMP=$((NCMP + 1))
        # pe5 carries mates placed only by mate rescue: with rescue off (-S) the
        # output must change, or neither leg above reached the rescue kernels.
        if [ "$set" = pe5 ]; then
            run "c_${set}_${tag}_norescue" - - "$BIN" ref.fa "$set" -S "${args[@]}"
            if cmp -s "c_${set}_${tag}_norescue.sam" "c_${set}_${tag}_batched.sam"; then
                fail "(c) $set $tag: output is unchanged with -S; no mate reaches rescue"
            fi
        fi
    done
done

if [ "${#TIERS[@]}" -le 1 ]; then
    echo "SKIP: (b) cross-tier -- one runnable tier on this host (${TIERS[*]:-none}); nothing to compare"
fi

# (d) no --meth: identical to the base binary
if [ -n "$BASE_BIN" ]; then
    for set in pe5 em se; do
        run "d_${set}_new" - - "$BIN" plain_new/ref.fa "$set"
        run "d_${set}_base" - - "$BASE_BIN" plain_base/ref.fa "$set"
        compare "(d) $set" "d_${set}_new" "d_${set}_base"
        NCMP=$((NCMP + 1))
    done
else
    echo "SKIP: (d) no-meth identity to the base binary -- BWA_MEM3_BASE is not set"
fi

BASE_NOTE=skipped
[ -n "$BASE_BIN" ] && BASE_NOTE=checked
echo "PASS: meth_context_identity ($NCMP comparison(s); tiers: ${TIERS[*]:-none}; base: $BASE_NOTE)"
