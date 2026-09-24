#!/usr/bin/env bash
# test/regression/rescue_prune_identity.sh
#
# Exact mate-rescue pruning (src/rescue_prune.h, wired in mem_matesw_batch_pre /
# _post) must leave the alignment output byte-identical to the full-window
# rescue. Both legs are the SAME binary; only the BWA3_RESCUE_PRUNE escape hatch
# differs:
#
#   BWA3_RESCUE_PRUNE=0  -> every rescue window computed in full (reference)
#   default              -> proven failures dropped (B1), windows narrowed to a
#                           proven hull (B2), N windows kept whole (candidate)
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
# Pruning runs only on aarch64 (src/bwamem_pair.cpp rescue_prune_on). There,
# BWA3_RESCUE_PRUNE_STATS=1 must show proven failures, narrowed windows and fewer
# rows kept than examined, or the identity would be vacuous; elsewhere the two
# legs take the same path and the non-vacuity check is reported as skipped.
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

for t in 1 4; do
    run_leg "$t" "full.t$t" BWA3_RESCUE_PRUNE=0
    run_leg "$t" "prune.t$t" BWA3_RESCUE_PRUNE_STATS=1
    if ! cmp -s "full.t$t.sam" "prune.t$t.sam"; then
        echo "FAIL: pruned rescue differs from full-window rescue at -t $t:" >&2
        diff "full.t$t.sam" "prune.t$t.sam" | head -20 >&2
        exit 1
    fi
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

echo "PASS: rescue_prune_identity (pruned == full-window rescue at -t 1 and -t 4; $stats)"
