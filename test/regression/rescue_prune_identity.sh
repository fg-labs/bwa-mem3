#!/usr/bin/env bash
# test/regression/rescue_prune_identity.sh
#
# Exact mate-rescue pruning (src/rescue_prune.h, wired in mem_matesw_batch_pre /
# _post), the banded rescue DP built on it (src/rescue_band.h) and the skipped
# or incremental post-rescue dedup (mem_matesw_batch_post) must leave the
# alignment output byte-identical to the full-window rescue with every dedup
# run. Every leg is the SAME binary; only the escape hatches differ:
#
#   REF_ENV (below): every BWA3_RESCUE_* shortcut off
#                        -> every rescue window computed in full by kswv with its
#                           original cell and its original NEON forms (USQADD,
#                           ROWPAIR, LAZYQE off), and every post-rescue dedup run in
#                           full (reference; every other leg runs the defaults)
#   BWA3_RESCUE_PRUNE=0 BWA3_RESCUE_BAND=0
#                        -> as the reference, but a dedup proven to be a no-op is
#                           skipped and a one-record insert done in O(n)
#   BWA3_RESCUE_BAND=0   -> proven failures dropped (B1), windows narrowed to a
#                           proven hull (B2), N windows kept whole; kswv on the hull
#   default              -> as above, plus the banded DP where the cost model picks it,
#                           and the banded start recovery (pass 1) of every eligible job
#   BWA3_RESCUE_BAND_COST=100000000
#                        -> the cost gate opened: every B2 hull that can be banded is,
#                           so the banded kernel runs whatever the cost model says
#
# The unit tests pin the filter's decisions against ksw_align2; this pins the
# plumbing around them, which no unit test reaches: the not-enqueued sentinel
# taking the ordinary failing-rescue path, the hull offset recorded in _pre and
# applied to the result in _post, and the mate staged from the filter's copy.
#
# The fixture is built programmatically (no committed test data): a generated
# reference, a bank of concordant FR pairs to seed the insert-size distribution
# (bwa refuses mate rescue until it can estimate one), and rescue pairs whose
# mate is made unseedable by spaced mismatches, so it only places via SW mate
# rescue -- some near the rescue threshold (narrowed windows), some unrelated to
# the reference (proven failures), and some with an N (full windows); anchors in a
# block present three times (one rescue per copy on the same mate, which is what
# exercises the dedup shortcuts); and 300 bp mates, long enough for the 16-bit
# kswv kernels (every 150 bp mate takes the 8-bit ones). It runs at -t 1 and -t 4.
#
# With RESCUE_TIERS set (a list of x86 tiers, e.g. "avx2 avx512bw"), the
# reference and default legs are also run under BWAMEM3_FORCE_TIER for each
# listed tier the host has, and must equal the unforced reference: the kswv
# rescue kernels differ per tier, so this is where each tier's 8- and 16-bit
# kernels are checked end to end. A listed tier above the host is reported SKIP:,
# the host's own tier is labelled as such (forcing it changes nothing), a force
# the dispatcher ignores fails, and a list that leaves nothing to force fails.
#
# The dedup shortcuts and the 16-bit kernels run on every architecture, so
# BWA3_RESCUE_PRUNE_STATS=1 must show skipped dedups, one-record inserts and
# 16-bit rescue jobs everywhere, and must show that each leg's switches took
# effect (nothing filtered or deduplicated early in the reference, nothing
# filtered in the dedup leg, nothing banded in the hull leg). Pruning and
# banding run wherever the build has a SIMD filter and band kernel: aarch64
# (NEON) and x86 builds whose SIMD floor is avx2 or avx512bw, at the fixture's
# default -k 19 (src/bwamem_pair.cpp rescue_prune_on, rescue_prune_runs,
# rescue_band_enabled). There, the stats must show proven failures, narrowed
# windows, fewer rows kept than examined and (with the gate opened) banded
# parents, or the identity would be vacuous; in a build without them (an x86
# floor below avx2) the pruning legs take the same path and their non-vacuity
# check is reported as skipped. A reference and a default leg at -k 25 must
# match too: the SIMD filters run at any threshold, but at the AVX-512BW kswv
# tier x86 prunes nothing there (its cost gate, rescue_prune_cost_ok) while it
# still bands pass 1, so the stats must show banded pass-1 jobs, and filtered
# jobs on every other tier and none at AVX-512BW. Two scorings other than the
# default must match the reference at the same scoring too: -B 6 and
# -O 8 -E 2, which the SIMD filters take, so their stats must show filtered
# jobs wherever they run; and -B 3, which the lemma refuses, so no filtered
# job. The band kernels take all three, so each must band pass 1.
# The same pairs as bisulfite reads,
# under --meth -B 4, --meth-scoring genomic and the default collapsed scoring,
# must match their references too (with samtools on PATH, since --meth emits
# BAM): filtered jobs on aarch64 at the first two, none at the collapsed one
# (refused) and none on x86.
#
# Inputs:
#   BWA_MEM3     — path to the bwa-mem3 binary under test
#   RESCUE_TIERS — optional: x86 tiers to force in addition (see above)
set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"
command -v python3 > /dev/null 2>&1 || {
    echo "SKIP: python3 not on PATH (it generates the fixture)"
    exit 0
}

case "$BWA_MEM3" in
    /*) BIN="$BWA_MEM3" ;;
    *) BIN="$PWD/$BWA_MEM3" ;;
esac
[ -x "$BIN" ] || {
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

python3 - << 'PY' || fail "fixture generation failed"
import random
rnd = random.Random(20260928)
L = 60000
ref = ''.join(rnd.choice('ACGT') for _ in range(L))
# A repeat: DUP_LEN bases at DUP_A copied to two more places, each copy with a
# substitution every 100 bases at its own phase, so an anchor inside it places at
# all three copies (within pen_unpaired) and mate rescue runs once per copy on the
# same mate, with a different score at each.
DUP_A, DUP_LEN = 10000, 600
for at, phase in ((30000, 50), (45000, 25)):
    dup = list(ref[DUP_A:DUP_A+DUP_LEN])
    for i in range(phase, DUP_LEN, 100):
        dup[i] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[dup[i]]
    ref = ref[:at] + ''.join(dup) + ref[at+DUP_LEN:]
# An exact repeat: DUP_LEN bases at DUP_B copied unchanged to two more places, so an
# anchor inside it rescues the same mate against three byte-identical windows.
DUP_B = 20000
for at in (38000, 52000):
    ref = ref[:at] + ref[DUP_B:DUP_B+DUP_LEN] + ref[at+DUP_LEN:]
with open('ref.fa', 'w') as f:
    f.write('>chrA\n')
    for i in range(0, L, 80):
        f.write(ref[i:i+80] + '\n')

def rc(s):
    c = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(c[b] for b in reversed(s))
def mism(b):
    return {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[b]

RL, INS = 150, 400
r1, r2 = [], []
def pair(name, anchor, mate):
    r1.append('@%s\n%s\n+\n%s\n' % (name, anchor, 'I' * len(anchor)))
    r2.append('@%s\n%s\n+\n%s\n' % (name, rc(mate), 'I' * len(mate)))

# Concordant pairs: only there to populate the insert-size distribution.
for k in range(400):
    p = rnd.randint(0, L - INS - 1)
    pair('c%d' % k, ref[p:p+RL], ref[p+INS-RL:p+INS])

for k in range(300):
    p = rnd.randint(0, L - INS - 1)
    anchor = ref[p:p+RL]
    mate = list(ref[p+INS-RL:p+INS])
    kind = k % 3
    if kind == 0:
        # Near the threshold: a mismatch every 14-17 bases leaves no exact 19-mer
        # (unseedable) while the SW score stays well above min_seed_len.
        i = rnd.randint(3, 10)
        while i < RL:
            mate[i] = mism(mate[i])
            i += rnd.randint(14, 17)
    elif kind == 1:
        # Unrelated to the window: the rescue SW is attempted and fails.
        mate = [rnd.choice('ACGT') for _ in range(RL)]
    else:
        # As kind 0 plus an N: the filter must keep the full window.
        i = rnd.randint(3, 10)
        while i < RL:
            mate[i] = mism(mate[i])
            i += rnd.randint(14, 17)
        mate[rnd.randint(0, RL - 1)] = 'N'
    pair('r%d_%d' % (kind, k), anchor, ''.join(mate))

# Anchors in the repeat, mates made unseedable as for kind 0. With the mate inside
# the repeat too, every copy rescues it (the later ones are one-record inserts
# into the deduped list, the third in O(n)); with the mate past the repeat's end,
# only the first copy does, and the later copies' dedups have nothing new to
# dedup (skips).
for k in range(60):
    p = DUP_A + (rnd.randint(0, DUP_LEN - INS) if k % 2 == 0 else rnd.randint(DUP_LEN - INS + 60, DUP_LEN - RL))
    anchor = ref[p:p+RL]
    mate = list(ref[p+INS-RL:p+INS])
    i = rnd.randint(3, 10)
    while i < RL:
        mate[i] = mism(mate[i])
        i += rnd.randint(14, 17)
    pair('d%d' % k, anchor, ''.join(mate))

# Anchors in the exact repeat, mates made unseedable as for kind 0 and inside the
# repeat: the second and third copies' rescues repeat the first byte for byte, so the
# filter answers them from its memo and _pre from the first copy's result.
for k in range(30):
    p = DUP_B + rnd.randint(0, DUP_LEN - INS)
    anchor = ref[p:p+RL]
    mate = list(ref[p+INS-RL:p+INS])
    i = rnd.randint(3, 10)
    while i < RL:
        mate[i] = mism(mate[i])
        i += rnd.randint(14, 17)
    pair('e%d' % k, anchor, ''.join(mate))

# 300 bp mates, made unseedable as for kind 0: a mate this long (times the match
# score, plus the kernel's shift) no longer fits a byte, so its rescue runs through
# the 16-bit kswv kernels.
ML = 300
for k in range(40):
    p = rnd.randint(0, L - INS - 1)
    anchor = ref[p:p+RL]
    mate = list(ref[p+INS-ML:p+INS])
    i = rnd.randint(3, 10)
    while i < ML:
        mate[i] = mism(mate[i])
        i += rnd.randint(14, 17)
    pair('l%d' % k, anchor, ''.join(mate))

with open('r1.fq', 'w') as f: f.write(''.join(r1))
with open('r2.fq', 'w') as f: f.write(''.join(r2))

# The same pairs as directional bisulfite reads for the --meth legs: every C of
# read 1 read as T (the OT strand, nothing methylated), every G of read 2 as A.
def conv(recs, f, t):
    out = []
    for r in recs:
        name, seq, plus, qual = r.split('\n')[:4]
        out.append('%s\n%s\n%s\n%s\n' % (name, seq.replace(f, t), plus, qual))
    return ''.join(out)
with open('m1.fq', 'w') as f: f.write(conv(r1, 'C', 'T'))
with open('m2.fq', 'w') as f: f.write(conv(r2, 'G', 'A'))
PY

"$BIN" index ref.fa > /dev/null 2>&1 || fail "index nonzero exit"

MEM_OPTS=() # extra `mem` options for the legs that follow; set per block
run_leg() { # $1 = threads, $2 = output stem, rest = env assignments
    local t="$1" stem="$2"
    shift 2
    env "$@" "$BIN" mem ${MEM_OPTS[@]+"${MEM_OPTS[@]}"} -t "$t" ref.fa r1.fq r2.fq 2> "$stem.err" \
        | grep -v '^@PG' > "$stem.sam" || fail "leg $stem nonzero exit"
    [ "$(grep -cv '^@' "$stem.sam" || true)" -gt 0 ] || fail "leg $stem produced no alignment records"
}

# The build's SIMD floor decides which shortcuts it carries: the SIMD filter and band kernel on
# neon, avx2 and avx512bw; neither below avx2 (an unrecognised floor is treated as carrying them,
# so its checks fail rather than skip).
floor="$("$BIN" version 2>&1 | sed -n 's/^SIMD floor: \([a-z0-9]*\).*/\1/p' | head -1 || true)"
[ -n "$floor" ] || fail "could not read the SIMD floor from \`$BIN version\` (no 'SIMD floor: <tier>' line)"
case "$floor" in
    sse41 | sse42 | avx | scalar) has_simd=0 ;;
    *) has_simd=1 ;;
esac
# The tier kswv runs at on this host (the dispatcher's choice, which the x86 cost gate keys on).
host_tier="$(BWAMEM3_DEBUG_SIMD=1 "$BIN" 2>&1 | sed -n 's/.*SIMD tier: \([a-z0-9]*\).*/\1/p' | head -1 || true)"

# The reference: every shortcut and alternative kernel form off. Defined once so
# the unforced and forced-tier references cannot drift apart.
REF_ENV=(BWA3_RESCUE_PRUNE=0 BWA3_RESCUE_BAND=0 BWA3_RESCUE_DEDUP_SKIP=0 BWA3_RESCUE_REPEAT=0
    BWA3_RESCUE_FSCAN=0 BWA3_RESCUE_USQADD=0 BWA3_RESCUE_ROWPAIR=0 BWA3_RESCUE_LAZYQE=0)

LEGS="dedup hull prune band"
for t in 1 4; do
    run_leg "$t" "full.t$t" BWA3_RESCUE_PRUNE_STATS=1 "${REF_ENV[@]}"
    run_leg "$t" "dedup.t$t" BWA3_RESCUE_PRUNE_STATS=1 BWA3_RESCUE_PRUNE=0 BWA3_RESCUE_BAND=0
    run_leg "$t" "hull.t$t" BWA3_RESCUE_PRUNE_STATS=1 BWA3_RESCUE_BAND=0
    run_leg "$t" "prune.t$t" BWA3_RESCUE_PRUNE_STATS=1
    run_leg "$t" "band.t$t" BWA3_RESCUE_PRUNE_STATS=1 BWA3_RESCUE_BAND_COST=100000000
    for leg in $LEGS; do
        if ! cmp -s "full.t$t.sam" "$leg.t$t.sam"; then
            echo "FAIL: rescue leg '$leg' differs from the reference rescue at -t $t:" >&2
            diff "full.t$t.sam" "$leg.t$t.sam" | head -20 >&2 || true
            exit 1
        fi
    done
done

# A seed length other than the default -k 19, so the filters run at another threshold: at the
# AVX-512BW kswv tier pruning's cost gate turns pruning off there (rescue_prune_cost_ok) while
# banded pass 1 still runs (rescue_band_runs), a combination no default-k leg reaches; elsewhere
# both run. Must equal the reference at the same -k.
MEM_OPTS=(-k 25)
run_leg 4 "k25full.t4" BWA3_RESCUE_PRUNE_STATS=1 "${REF_ENV[@]}"
run_leg 4 "k25.t4" BWA3_RESCUE_PRUNE_STATS=1
MEM_OPTS=()
if ! cmp -s k25full.t4.sam k25.t4.sam; then
    echo "FAIL: rescue at -k 25 differs from the reference at -k 25 (-t 4):" >&2
    diff k25full.t4.sam k25.t4.sam | head -20 >&2 || true
    exit 1
fi
kstats="$(grep '^\[RESCUE_PRUNE\]' k25.t4.err || true)"
[ -n "$kstats" ] || fail "BWA3_RESCUE_PRUNE_STATS=1 printed no [RESCUE_PRUNE] line (-k 25 leg)"
kjobs=$(printf '%s\n' "$kstats" | tr ' ' '\n' | sed -n 's/^jobs=//p')
kbstats="$(grep '^\[RESCUE_BAND\]' k25.t4.err || true)"
kp1=$(printf '%s\n' "$kbstats" | tr ' ' '\n' | sed -n 's/^pass1_banded=//p')
k25_note="; -k 25: jobs=${kjobs:-0} pass1_banded=${kp1:-0}"
if [ "$has_simd" = 1 ]; then
    [ "${kp1:-0}" -gt 0 ] || fail "no banded pass-1 job at -k 25 (SIMD floor '$floor'): $kstats $kbstats"
    [ -n "$host_tier" ] || fail "could not detect the kswv SIMD tier from BWAMEM3_DEBUG_SIMD output"
    case "$host_tier" in
        avx512bw) [ "${kjobs:-0}" -eq 0 ] || fail "x86 pruned at -k 25 at the AVX-512BW tier, past its cost gate: $kstats" ;;
        *) [ "${kjobs:-0}" -gt 0 ] || fail "pruning filtered no rescue job at -k 25 (kswv tier '$host_tier'): $kstats" ;;
    esac
fi

# Other scorings (rescue_prune_params), each against the reference at the same scoring: -B 6 (the
# default's bound weights) and -O 8 -E 2 (the per-diagonal charge c = 2, where a single-hit
# diagonal weighs a - c < 0), which the SIMD filters take, so they prune wherever those run; and
# -B 3, which the lemma refuses, so it prunes nowhere. The band kernels take all three, so each
# must band pass 1 (rescue_band_runs).
sc_note=""
check_scoring() { # $1 = label, $2 = where it prunes (simd | none), rest = mem options
    local label="$1" where="$2" stem st sj bst sp1
    shift 2
    stem="sc_$(printf '%s' "$label" | tr -c 'A-Za-z0-9' '_')"
    MEM_OPTS=("$@")
    run_leg 4 "$stem.full" BWA3_RESCUE_PRUNE_STATS=1 "${REF_ENV[@]}"
    run_leg 4 "$stem" BWA3_RESCUE_PRUNE_STATS=1
    MEM_OPTS=()
    if ! cmp -s "$stem.full.sam" "$stem.sam"; then
        echo "FAIL: rescue at $label differs from the reference at $label (-t 4):" >&2
        diff "$stem.full.sam" "$stem.sam" | head -20 >&2 || true
        exit 1
    fi
    st="$(grep '^\[RESCUE_PRUNE\]' "$stem.err" || true)"
    [ -n "$st" ] || fail "BWA3_RESCUE_PRUNE_STATS=1 printed no [RESCUE_PRUNE] line ($label leg)"
    sj=$(printf '%s\n' "$st" | tr ' ' '\n' | sed -n 's/^jobs=//p')
    bst="$(grep '^\[RESCUE_BAND\]' "$stem.err" || true)"
    sp1=$(printf '%s\n' "$bst" | tr ' ' '\n' | sed -n 's/^pass1_banded=//p')
    if [ "$has_simd" = 1 ]; then
        case "$where/$floor" in
            simd/*) [ "${sj:-0}" -gt 0 ] || fail "pruning filtered no rescue job at $label (SIMD floor '$floor'): $st" ;;
            *) [ "${sj:-0}" -eq 0 ] || fail "pruning filtered rescue jobs at $label, which the lemma refuses: $st" ;;
        esac
        [ "${sp1:-0}" -gt 0 ] || fail "no banded pass-1 job at $label: $bst"
    fi
    sc_note="$sc_note; $label: jobs=${sj:-0} pass1_banded=${sp1:-0}"
}
check_scoring "-B 6" simd -B 6
check_scoring "-O 8 -E 2" simd -O 8 -E 2
check_scoring "-B 3" none -B 3

# --meth (EM-seq chemistry, the default): the filter matches converted copies of the window and the
# mate (rescue_prune_params::set_meth). Each leg must equal the reference at the same options; the
# stats must show filtered jobs on aarch64 at --meth -B 4 and --meth-scoring genomic, and none at
# the default collapsed scoring (b = 2a, which the lemma refuses) or on x86 (its cost gate). --meth
# emits BAM, so these legs need samtools; without it they are reported SKIP:.
meth_note=""
if command -v samtools > /dev/null 2>&1; then
    "$BIN" index --meth ref.fa > /dev/null 2>&1 || fail "index --meth nonzero exit"
    meth_leg() { # $1 = output stem, $2 = options (one word per option), rest = env assignments
        local stem="$1" opts="$2"
        shift 2
        # shellcheck disable=SC2086 # opts is a list of option words
        env "$@" "$BIN" mem --meth $opts -t 4 ref.fa m1.fq m2.fq 2> "$stem.err" > "$stem.bam" \
            || fail "leg $stem nonzero exit"
        samtools view "$stem.bam" > "$stem.sam" || fail "leg $stem: samtools view failed"
        [ -s "$stem.sam" ] || fail "leg $stem produced no alignment records"
    }
    check_meth() { # $1 = label, $2 = prunes on aarch64 (1 | 0), $3 = options
        local label="$1" prunes="$2" opts="$3" stem st mj
        stem="meth_$(printf '%s' "$label" | tr -c 'A-Za-z0-9' '_')"
        meth_leg "$stem.full" "$opts" BWA3_RESCUE_PRUNE_STATS=1 "${REF_ENV[@]}"
        meth_leg "$stem" "$opts" BWA3_RESCUE_PRUNE_STATS=1
        if ! cmp -s "$stem.full.sam" "$stem.sam"; then
            echo "FAIL: rescue at --meth $label differs from the reference at the same options (-t 4):" >&2
            diff "$stem.full.sam" "$stem.sam" | head -20 >&2 || true
            exit 1
        fi
        st="$(grep '^\[RESCUE_PRUNE\]' "$stem.err" || true)"
        [ -n "$st" ] || fail "BWA3_RESCUE_PRUNE_STATS=1 printed no [RESCUE_PRUNE] line (--meth $label leg)"
        mj=$(printf '%s\n' "$st" | tr ' ' '\n' | sed -n 's/^jobs=//p')
        if [ "$has_simd" = 1 ]; then
            if [ "$prunes" = 1 ] && [ "$floor" = neon ]; then
                [ "${mj:-0}" -gt 0 ] || fail "pruning filtered no rescue job at --meth $label on aarch64: $st"
            else
                [ "${mj:-0}" -eq 0 ] || fail "pruning filtered rescue jobs at --meth $label (SIMD floor '$floor'), past its gate: $st"
            fi
        fi
        meth_note="$meth_note; --meth $label: jobs=${mj:-0}"
    }
    check_meth "-B 4" 1 "-B 4"
    check_meth "genomic" 1 "--meth-scoring genomic"
    check_meth "collapsed" 0 ""
else
    echo "SKIP: samtools not on PATH; the --meth legs did not run"
fi

# Forced tiers: the reference and default legs under each listed tier the host
# has must equal the unforced reference. Ranked against the dispatcher's own
# tier names (as in all_tiers_parity.sh): BWAMEM3_FORCE_TIER only downgrades.
tier_note=""
if [ -n "${RESCUE_TIERS:-}" ]; then
    [ -n "$host_tier" ] || fail "could not detect the host SIMD tier from BWAMEM3_DEBUG_SIMD output"
    rank() { # tier -> rank among the x86 tiers, or -1
        case "$1" in
            sse41) echo 0 ;; sse42) echo 1 ;; avx) echo 2 ;; avx2) echo 3 ;; avx512bw) echo 4 ;; *) echo -1 ;;
        esac
    }
    host_rank="$(rank "$host_tier")"
    [ "$host_rank" -ge 0 ] || fail "RESCUE_TIERS is set but the host tier '$host_tier' is not an x86 tier"
    ran=""
    for tier in $RESCUE_TIERS; do
        [ "$(rank "$tier")" -ge 0 ] || fail "RESCUE_TIERS names an unknown x86 tier '$tier'"
        if [ "$(rank "$tier")" -gt "$host_rank" ]; then
            echo "SKIP: tier $tier is above this host's ($host_tier); not forced"
            continue
        fi
        for t in 1 4; do
            run_leg "$t" "full.$tier.t$t" BWAMEM3_FORCE_TIER="$tier" "${REF_ENV[@]}"
            run_leg "$t" "prune.$tier.t$t" BWAMEM3_FORCE_TIER="$tier"
            for leg in full prune; do
                # An ignored force would run the host tier and compare equal.
                if grep -q 'ignoring BWAMEM3_FORCE_TIER' "$leg.$tier.t$t.err"; then
                    fail "the dispatcher ignored BWAMEM3_FORCE_TIER=$tier: $(grep 'ignoring BWAMEM3_FORCE_TIER' "$leg.$tier.t$t.err" | head -1)"
                fi
                if ! cmp -s "full.t$t.sam" "$leg.$tier.t$t.sam"; then
                    echo "FAIL: rescue leg '$leg' under BWAMEM3_FORCE_TIER=$tier differs from the reference at -t $t:" >&2
                    diff "full.t$t.sam" "$leg.$tier.t$t.sam" | head -20 >&2 || true
                    exit 1
                fi
            done
        done
        if [ "$tier" = "$host_tier" ]; then ran="$ran $tier(host tier)"; else ran="$ran $tier"; fi
    done
    [ -n "$ran" ] || fail "RESCUE_TIERS='$RESCUE_TIERS' left no tier to force on this host ($host_tier)"
    tier_note="; forced tiers:$ran"
fi

field() { printf '%s\n' "$stats" | tr ' ' '\n' | sed -n "s/^$1=//p"; }

# The reference leg's switches took effect: nothing filtered, no dedup skipped.
stats="$(grep '^\[RESCUE_PRUNE\]' full.t1.err || true)"
[ -n "$stats" ] || fail "BWA3_RESCUE_PRUNE_STATS=1 printed no [RESCUE_PRUNE] line (reference leg)"
[ "$(field jobs)" -eq 0 ] || fail "the reference leg filtered rescue jobs (BWA3_RESCUE_PRUNE=0 ignored): $stats"
[ "$(field dedup_skip)" -eq 0 ] && [ "$(field dedup_insert1)" -eq 0 ] \
    || fail "the reference leg took a dedup shortcut (BWA3_RESCUE_DEDUP_SKIP=0 ignored): $stats"

stats="$(grep '^\[RESCUE_PRUNE\]' dedup.t1.err || true)"
[ -n "$stats" ] || fail "BWA3_RESCUE_PRUNE_STATS=1 printed no [RESCUE_PRUNE] line (dedup leg)"
[ "$(field jobs)" -eq 0 ] || fail "the dedup leg filtered rescue jobs (BWA3_RESCUE_PRUNE=0 ignored): $stats"
jobs16=$(field jobs16)
[ "${jobs16:-0}" -gt 0 ] || fail "no 16-bit rescue job in the fixture (the 300 bp mates stopped reaching it): $stats"
dedup_skip=$(field dedup_skip) dedup_insert1=$(field dedup_insert1) dedup_fast=$(field dedup_insert1_fast)
[ "${dedup_skip:-0}" -gt 0 ] || fail "no post-rescue dedup was skipped as a no-op: $stats"
[ "${dedup_insert1:-0}" -gt 0 ] || fail "no post-rescue dedup took the one-record insert: $stats"
[ "${dedup_fast:-0}" -gt 0 ] || fail "no one-record insert was done in O(n): $stats"
dstats="dedup_skip=$dedup_skip dedup_insert1=$dedup_insert1 dedup_insert1_fast=$dedup_fast jobs16=$jobs16"

stats="$(grep '^\[RESCUE_PRUNE\]' prune.t1.err || true)"
[ -n "$stats" ] || fail "BWA3_RESCUE_PRUNE_STATS=1 printed no [RESCUE_PRUNE] line"
jobs=$(field jobs) full=$(field full) b1=$(field b1) b2=$(field b2) rows_in=$(field rows_in) rows_kept=$(field rows_kept)
if [ "$jobs" -eq 0 ]; then
    # Pruning is compiled in wherever the build's SIMD floor carries the SIMD filter, so no
    # filtered job there means it stopped engaging (a gate or fixture regression) and the identity
    # above held vacuously.
    [ "$has_simd" = 0 ] \
        || fail "pruning filtered no rescue job, though this build (SIMD floor '$floor') has the SIMD filter: $stats"
    echo "PASS: rescue_prune_identity (11-op cell, NEON-form defaults and dedup shortcuts == reference at -t 1 and -t 4; $dstats$tier_note$k25_note$sc_note$meth_note)"
    echo "SKIP: pruning does not run in this build (SIMD floor $floor, below the SIMD filter's avx2); its legs held trivially"
    exit 0
fi
[ "$full" -gt 0 ] || fail "no filtered job kept its full window (the N mates stopped reaching the filter): $stats"
[ "$b1" -gt 0 ] || fail "no proven-failure (B1) rescue in the fixture: $stats"
[ "$b2" -gt 0 ] || fail "no narrowed (B2) rescue in the fixture: $stats"
[ "$rows_kept" -lt "$rows_in" ] || fail "pruning kept every row: $stats"
# The anchors present three times rescue the same mate against identical windows: the repeats
# after the first must be answered from its result (BWA3_RESCUE_REPEAT), or that path went untested.
reused=$(field reused)
[ "${reused:-0}" -gt 0 ] || fail "no rescue job was answered from an identical earlier job's result: $stats"

bstats="$(grep '^\[RESCUE_BAND\]' band.t1.err || true)"
[ -n "$bstats" ] || fail "BWA3_RESCUE_PRUNE_STATS=1 printed no [RESCUE_BAND] line"
bfield() { printf '%s\n' "$bstats" | tr ' ' '\n' | sed -n "s/^$1=//p"; }
parents=$(bfield banded_parents)
[ "${parents:-0}" -gt 0 ] || fail "no banded rescue parent with the cost gate opened: $bstats"
# The hull leg's switch took effect: nothing banded in either pass.
hstats="$(grep '^\[RESCUE_BAND\]' hull.t1.err || true)"
for f in banded_parents pass1_banded; do
    v=$(printf '%s\n' "$hstats" | tr ' ' '\n' | sed -n "s/^$f=//p")
    [ "${v:-0}" -eq 0 ] || fail "the hull leg banded jobs (BWA3_RESCUE_BAND=0 ignored): $hstats"
done
p1stats="$(grep '^\[RESCUE_BAND\]' prune.t1.err || true)"
p1=$(printf '%s\n' "$p1stats" | tr ' ' '\n' | sed -n 's/^pass1_banded=//p')
[ "${p1:-0}" -gt 0 ] || fail "no banded pass-1 (start recovery) job at the defaults: $p1stats"

echo "PASS: rescue_prune_identity (11-op cell, NEON-form defaults, dedup shortcuts, hull, pruned and banded == reference at -t 1 and -t 4; $dstats$tier_note$k25_note$sc_note$meth_note; $stats; $bstats)"
