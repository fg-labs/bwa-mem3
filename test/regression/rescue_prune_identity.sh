#!/usr/bin/env bash
# test/regression/rescue_prune_identity.sh
#
# Exact mate-rescue pruning (src/rescue_prune.h, wired in mem_matesw_batch_pre /
# _post) and the banded rescue DP built on it (src/rescue_band.h) must leave the
# alignment output byte-identical to the full-window rescue. Every leg is the
# SAME binary; only the escape hatches differ:
#
#   BWA3_RESCUE_PRUNE=0 BWA3_RESCUE_BAND=0
#                        -> every rescue window computed in full by kswv (reference)
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
# the reference (proven failures), and some with an N (full windows). It runs at
# -t 1 and -t 4.
#
# Pruning and banding run only on aarch64 (src/bwamem_pair.cpp rescue_prune_on,
# rescue_band_enabled). There, BWA3_RESCUE_PRUNE_STATS=1 must show proven
# failures, narrowed windows, fewer rows kept than examined and (with the gate
# opened) banded parents, or the identity would be vacuous; elsewhere every leg
# takes the same path and the non-vacuity check is reported as skipped.
#
# Inputs:
#   BWA_MEM3 — path to the bwa-mem3 binary under test
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
q = 'I' * RL
r1, r2 = [], []
def pair(name, anchor, mate):
    r1.append('@%s\n%s\n+\n%s\n' % (name, anchor, q))
    r2.append('@%s\n%s\n+\n%s\n' % (name, rc(mate), q))

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

with open('r1.fq', 'w') as f: f.write(''.join(r1))
with open('r2.fq', 'w') as f: f.write(''.join(r2))
PY

"$BIN" index ref.fa > /dev/null 2>&1 || fail "index nonzero exit"

run_leg() { # $1 = threads, $2 = output stem, rest = env assignments
    local t="$1" stem="$2"
    shift 2
    env "$@" "$BIN" mem -t "$t" ref.fa r1.fq r2.fq 2> "$stem.err" | grep -v '^@PG' > "$stem.sam" \
        || fail "leg $stem nonzero exit"
    [ "$(grep -cv '^@' "$stem.sam" || true)" -gt 0 ] || fail "leg $stem produced no alignment records"
}

LEGS="hull prune band"
for t in 1 4; do
    run_leg "$t" "full.t$t" BWA3_RESCUE_PRUNE=0 BWA3_RESCUE_BAND=0
    run_leg "$t" "hull.t$t" BWA3_RESCUE_BAND=0
    run_leg "$t" "prune.t$t" BWA3_RESCUE_PRUNE_STATS=1
    run_leg "$t" "band.t$t" BWA3_RESCUE_PRUNE_STATS=1 BWA3_RESCUE_BAND_COST=100000000
    for leg in $LEGS; do
        if ! cmp -s "full.t$t.sam" "$leg.t$t.sam"; then
            echo "FAIL: rescue leg '$leg' differs from the full-window rescue at -t $t:" >&2
            diff "full.t$t.sam" "$leg.t$t.sam" | head -20 >&2 || true
            exit 1
        fi
    done
done

stats="$(grep '^\[RESCUE_PRUNE\]' prune.t1.err || true)"
[ -n "$stats" ] || fail "BWA3_RESCUE_PRUNE_STATS=1 printed no [RESCUE_PRUNE] line"
field() { printf '%s\n' "$stats" | tr ' ' '\n' | sed -n "s/^$1=//p"; }
jobs=$(field jobs) b1=$(field b1) b2=$(field b2) rows_in=$(field rows_in) rows_kept=$(field rows_kept)
if [ "$jobs" -eq 0 ]; then
    # Pruning is compiled in on aarch64, so no filtered job there means it stopped engaging (a
    # gate or fixture regression) and the identity above held vacuously.
    case "$(uname -m)" in
        aarch64 | arm64) fail "pruning filtered no rescue job on an aarch64 host: $stats" ;;
    esac
    echo "SKIP: pruning does not run on this host (aarch64 only); the identity above held trivially"
    exit 0
fi
[ "$b1" -gt 0 ] || fail "no proven-failure (B1) rescue in the fixture: $stats"
[ "$b2" -gt 0 ] || fail "no narrowed (B2) rescue in the fixture: $stats"
[ "$rows_kept" -lt "$rows_in" ] || fail "pruning kept every row: $stats"

bstats="$(grep '^\[RESCUE_BAND\]' band.t1.err || true)"
[ -n "$bstats" ] || fail "BWA3_RESCUE_PRUNE_STATS=1 printed no [RESCUE_BAND] line"
bfield() { printf '%s\n' "$bstats" | tr ' ' '\n' | sed -n "s/^$1=//p"; }
parents=$(bfield banded_parents)
[ "${parents:-0}" -gt 0 ] || fail "no banded rescue parent with the cost gate opened: $bstats"
p1stats="$(grep '^\[RESCUE_BAND\]' prune.t1.err || true)"
p1=$(printf '%s\n' "$p1stats" | tr ' ' '\n' | sed -n 's/^pass1_banded=//p')
[ "${p1:-0}" -gt 0 ] || fail "no banded pass-1 (start recovery) job at the defaults: $p1stats"

echo "PASS: rescue_prune_identity (hull, pruned and banded == full-window rescue at -t 1 and -t 4; $stats; $bstats)"
