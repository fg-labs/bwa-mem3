#!/usr/bin/env bash
# test/regression/fixed_insert_slice_identity.sh
#
# Regression: with -I, reading a batch in slices must not change the output of
# the mate-concordant chain cap (--extend-mate-concordant auto, which --fast
# turns on).
#
# The auto window is the proper-pair insert high bound. Under -I that bound is
# known up front, and the align pass must use it from the first batch on. The
# whole-batch align path (mem_process_seqs) published it before aligning, but
# the sliced path (mem_align_cohort_slice) did not, so a first batch read in
# slices was capped with the 1000 bp fallback window instead. The same `--fast
# -I` run then gave different records with and without cohort slicing; on a
# real input only the first batch is sliced, so the drift hid in the first few
# million reads.
#
# The fixture is built so the window size matters: a 600 bp block repeated 26
# times, 800 bp apart, so a mate's other copies sit between -I's high bound
# (400) and the fallback (1000). The test asserts that (fixed windows of 400 and
# 1000 give different output), that the auto window with -I matches the fixed
# 400 bp window, and that the runs really slice, so no identity check can pass
# vacuously.
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
# python3 generates the fixture; without it the step would otherwise go green
# untested, so its absence is a failure.
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

# 40,000 2x100 pairs (8 Mbp) with 260-340 bp inserts over a 60 kb reference
# whose 600 bp block is copied 26 times, 800 bp apart, each copy with up to three
# substitutions.
python3 - << 'PY' || fail "fixture generation failed"
import random
rnd = random.Random(11)
L, BLOCK = 60000, 600
ref = [rnd.choice('ACGT') for _ in range(L)]
blk = ref[1000:1000 + BLOCK]
sub = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}
for k in range(26):
    at = 2000 + k * 800
    b = list(blk)
    for _ in range(rnd.randint(0, 3)):
        j = rnd.randrange(BLOCK)
        b[j] = sub[b[j]]
    ref[at:at + BLOCK] = b
ref = ''.join(ref)
with open('ref.fa', 'w') as f:
    f.write('>ref\n' + '\n'.join(ref[i:i + 60] for i in range(0, L, 60)) + '\n')
comp = str.maketrans('ACGT', 'TGCA')
with open('r1.fq', 'w') as a, open('r2.fq', 'w') as b:
    for n in range(40000):
        ins = rnd.randint(260, 340)
        s = rnd.randint(0, L - ins)
        frag = ref[s:s + ins]
        r1, r2 = frag[:100], frag[-100:].translate(comp)[::-1]
        a.write(f'@p{n}/1\n{r1}\n+\n{"I" * 100}\n')
        b.write(f'@p{n}/2\n{r2}\n+\n{"I" * 100}\n')
PY

"$BIN" index ref.fa > /dev/null 2>&1 || fail "index ref.fa"

# -K 4000000 makes two 4 Mbp batches, so slicing has a first batch to split.
INS=(-I "300,20,400,200")
COMMON=(-t 4 -K 4000000)

# Run `bwa-mem3 mem` with the given args; records (sorted) to $1.rec, stderr to $1.log.
run() {
    local out="$1"
    shift
    "$BIN" mem "$@" > "$out.sam" 2> "$out.log" || {
        cat "$out.log" >&2
        fail "mem $* exited nonzero"
    }
    # grep exits 1 when it selects nothing (a header-only SAM); let the empty
    # check below report that, and fail here only on a real read or sort error.
    { grep -v '^@' "$out.sam" || [[ $? -eq 1 ]]; } | LC_ALL=C sort > "$out.rec" \
        || fail "$out: extracting records from $out.sam failed"
    [[ -s "$out.rec" ]] || fail "$out produced no records"
}
slices() { grep -c '(cohort slice)$' "$1.log" || true; }

# --- non-vacuity: the mate window size changes this fixture's output ----------
run win400 "${COMMON[@]}" --max-extend-chains 20 --extend-mate-concordant=400 "${INS[@]}" --cohort-slices 0 ref.fa r1.fq r2.fq
run win1000 "${COMMON[@]}" --max-extend-chains 20 --extend-mate-concordant=1000 "${INS[@]}" --cohort-slices 0 ref.fa r1.fq r2.fq
if cmp -s win400.rec win1000.rec; then
    fail "fixed mate windows of 400 and 1000 gave the same output; the fixture cannot tell -I's bound from the fallback"
fi
echo "PASS: the mate-concordant window size changes this fixture's output"

# --- identity: unsliced == first batch sliced == every batch sliced -----------
# $1 = label, rest = the mode's options.
check_mode() {
    local label="$1"
    shift
    run "${label}_unsliced" "${COMMON[@]}" "$@" "${INS[@]}" --cohort-slices 0 ref.fa r1.fq r2.fq
    run "${label}_first" "${COMMON[@]}" "$@" "${INS[@]}" --cohort-slices 3 ref.fa r1.fq r2.fq
    BWA_MEM3_COHORT_SLICE_ALL=1 run "${label}_all" "${COMMON[@]}" "$@" "${INS[@]}" --cohort-slices 3 ref.fa r1.fq r2.fq
    [[ $(slices "${label}_unsliced") -eq 0 ]] || fail "$label: --cohort-slices 0 produced partial slices"
    [[ $(slices "${label}_first") -gt 0 ]] || fail "$label: --cohort-slices 3 produced no partial slices"
    [[ $(slices "${label}_all") -gt $(slices "${label}_first") ]] \
        || fail "$label: BWA_MEM3_COHORT_SLICE_ALL did not slice beyond the first batch"
    cmp -s "${label}_unsliced.rec" "${label}_first.rec" \
        || fail "$label: slicing the first batch changed $(diff "${label}_unsliced.rec" "${label}_first.rec" | grep -c '^<' || true) records"
    cmp -s "${label}_unsliced.rec" "${label}_all.rec" \
        || fail "$label: slicing every batch changed $(diff "${label}_unsliced.rec" "${label}_all.rec" | grep -c '^<' || true) records"
    echo "PASS: $label with -I is unchanged by cohort slicing"
}

check_mode fast --fast
check_mode auto_window --max-extend-chains 20 --extend-mate-concordant
# The auto window must be -I's high bound (400), not the 1000 bp fallback: if
# every auto run fell back, the slicing comparisons above would still agree.
cmp -s auto_window_unsliced.rec win400.rec \
    || fail "auto mate window with -I does not match the fixed 400 bp window ($(diff auto_window_unsliced.rec win400.rec | grep -c '^<' || true) records differ)"
echo "PASS: the auto mate window with -I matches the fixed 400 bp window"
# The default path has no chain cap; pin that slicing leaves it alone too.
check_mode default
