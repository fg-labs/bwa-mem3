#!/usr/bin/env bash
# test/regression/proper_pair_alt.sh
#
# End-to-end behavior of --proper-pair-from-emitted: which alignment the
# proper-pair FLAG bit (0x2) is derived from.
#
# bwa and bwa-mem2 both derive 0x2 from the top-scoring region a[0] even when
# the record they emit is a[which] (bwa bwamem_pair.c:411; bwa-mem2 inherited it
# verbatim). fg-labs/bwa-mem3#17 switched bwa-mem3 to a[which]; #362 moved that
# behind --proper-pair-from-emitted and restored the upstream default.
#
# The two derivations differ only when `which != 0`, which requires the read to
# have ALT hits -- so on any reference without a `.alt` sidecar the option is
# structurally inert, and no other fixture in this suite can reach the branch at
# all. That is the whole reason this script exists: without it, a consumer that
# ignored the option and hardcoded either index would pass every other
# regression here and every unit test, and would surface only as FLAG diffs on
# a real ALT-aware run -- 3,013 of 10,134,006 records on the one WGS slice
# measured for #362, which is that input's count and not a rate to expect
# elsewhere.
#
# Reaching the branch takes a read whose top PRIMARY region scores below `-T`
# while its top ALT region clears it. The fixture below builds exactly that: a
# read pair drawn from the ALT contig, with a 25 bp exact copy of each mate's 5'
# end planted far apart on the primary contig. Each mate then has a 120-point
# ALT hit and a 25-point primary hit (`-T` defaults to 30), so `which` selects
# the ALT region while a[0] stays the weak primary one. The two primary hits sit
# ~1,000 bp apart, outside the insert window the background pairs establish;
# the two ALT hits are 300 bp apart, inside it. So the bit is CLEAR by default
# (upstream's answer) and SET under the option, and the difference is the FLAG
# field alone.
#
# Fixture is generated here, not committed (repo convention). Sequence comes
# from python3's Mersenne Twister at a fixed seed, so it is identical on every
# run and platform.
#
# Inputs (env vars):
#   BWA_MEM3                  — path to the bwa-mem3 binary under test
#   PROPER_PAIR_ALT_WORK_DIR  — directory for fixture-private intermediates.
#                               This script DELETES the fixture files it owns
#                               inside it, so give it a directory of its own.

set -Eeuo pipefail

# Under `set -e` alone an unexpected failure -- python3, `index`, `mem` -- aborts
# with a bare nonzero status and no marker, which reads as "the harness broke"
# without naming what broke. This gives those the same `FAIL:` marker the
# deliberate assertions below use. `-E` is what makes the trap reach inside the
# helper functions; an explicit `exit` after a `FAIL:` echo does not re-trigger
# it, and a failure that is tested (`if ! cmd`, `cmd || true`) does not either,
# so the assertions stay single-reported.
# shellcheck disable=SC2154  # rc IS assigned by `rc=$?` below; shellcheck does
# not track assignments made inside a quoted trap body (false positive).
trap 'rc=$?; echo "FAIL: unexpected failure at line $LINENO (exit $rc): $BASH_COMMAND" >&2; exit "$rc"' ERR

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${PROPER_PAIR_ALT_WORK_DIR:?PROPER_PAIR_ALT_WORK_DIR must be set}"

mkdir -p "$PROPER_PAIR_ALT_WORK_DIR"
cd "$PROPER_PAIR_ALT_WORK_DIR"
# Remove only what a previous run of THIS script created. A blanket
# `rm -rf "$PROPER_PAIR_ALT_WORK_DIR"` would obey a mistyped or over-broad value
# of a caller-supplied variable; the `:?` guard above catches unset and empty,
# but not `/` or `$HOME`.
rm -f ref.fasta ref.fasta.* noalt.fasta noalt.fasta.* index.log \
    r1.fq r2.fq ./*.sam ./*.err

ALT_CONTIG=chr1_alt
PAIR=altpair # the one read pair that reaches the divergent branch

# --- Generate the ALT-aware fixture (deterministic; fixed PRNG seed). ---
python3 - << 'PY'
import random

random.seed(363)                      # fixed: fixture must be reproducible


def rnd(n):
    return "".join(random.choice("ACGT") for _ in range(n))


comp = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s[::-1].translate(comp)


CHR1_LEN, ALT_LEN = 30000, 6000
chr1 = list(rnd(CHR1_LEN))
alt = rnd(ALT_LEN)

# The divergent pair: both mates drawn from the ALT contig, FR, 300 bp insert.
R1 = alt[1000:1120]
R2 = rc(alt[1180:1300])

# Plant a 25 bp exact copy of each mate's 5' end on chr1, 1,000 bp apart. 25 is
# above the 19 bp minimum seed length (so the seed is found and an alnreg is
# built) and below the default -T of 30 (so `which` skips past it to the ALT
# region). R2 is already reverse-complemented, so plant rc(R2[:25]) to give it a
# forward-strand primary hit like the background pairs have.
chr1[25000:25025] = list(R1[:25])
chr1[26000:26025] = list(rc(R2[:25]))
chr1 = "".join(chr1)

with open("ref.fasta", "w") as f:
    for name, s in (("chr1", chr1), ("chr1_alt", alt)):
        f.write(">%s\n" % name)
        for i in range(0, len(s), 60):
            f.write(s[i:i + 60] + "\n")
# <prefix>.alt is a SAM file; bwa reads only column 1 (QNAME) to set is_alt.
with open("ref.fasta.alt", "w") as f:
    f.write("chr1_alt\t0\tchr1\t1\t60\t6000M\t*\t0\t0\t*\t*\n")

r1, r2 = [], []
# 30 background FR pairs off chr1 at a 300 bp insert. mem_pestat needs at least
# MIN_DIR_CNT (10) pairs in a direction before it will call that direction
# converged; without them pes[].failed is set and NOTHING is ever properly
# paired, which would make both modes agree for the wrong reason.
for i, off in enumerate(range(1000, 13000, 400)):
    r1.append(("bg%d" % i, chr1[off:off + 120]))
    r2.append(("bg%d" % i, rc(chr1[off + 180:off + 300])))
r1.append(("altpair", R1))
r2.append(("altpair", R2))
for path, recs in (("r1.fq", r1), ("r2.fq", r2)):
    with open(path, "w") as f:
        for n, s in recs:
            f.write("@%s\n%s\n+\n%s\n" % (n, s, "I" * len(s)))
PY

# The same reference WITHOUT the sidecar, to pin the inertness claim below.
cp ref.fasta noalt.fasta

"$BWA_MEM3" index ref.fasta > index.log 2>&1
"$BWA_MEM3" index noalt.fasta >> index.log 2>&1

# -t 1 keeps the batch partition (and so mem_pestat's sample) fixed; the insert
# statistics feed the window this test reads, so a thread-count-dependent
# partition would make the assertions below flaky rather than wrong.
run() { # <out-prefix> <ref> [extra flags...]
    local out=$1 ref=$2
    shift 2
    "$BWA_MEM3" mem -t 1 "$@" "$ref" r1.fq r2.fq > "$out.sam" 2> "$out.err"
}
run def ref.fasta
run emi ref.fasta --proper-pair-from-emitted

# --- 1. The fixture actually reaches the divergent branch. ----------------
# Assert the preconditions rather than trusting the construction: if a scoring
# or seeding change ever pushed the planted 25-mer above -T, or dropped it
# below the seed-length floor, `which` would be 0 and every assertion below
# would pass vacuously with both modes agreeing.
flags_of() { # <prefix> -- the pair's two FLAG fields, in file order
    grep "^$PAIR	" "$1.sam" | cut -f2
}
# Checked on BOTH streams, not just the default one. An unset array element is
# 0 in bash arithmetic rather than an error, so a run that emitted only one mate
# would not abort at the `mapfile` below -- it would reach the FLAG assertions
# with a phantom 0 and fail there, reporting a missing proper-pair bit when the
# real defect is a missing record. Assert the count up front so the diagnostic
# names the actual problem.
assert_pair_record_count() { # <prefix>
    local n
    n=$(grep -c "^$PAIR	" "$1.sam" || true)
    if [[ "$n" -ne 2 ]]; then
        echo "FAIL: expected exactly 2 records for $PAIR in $1.sam, found $n" >&2
        grep "^$PAIR	" "$1.sam" | cut -f1-6 >&2 || true
        exit 1
    fi
}
assert_pair_record_count def
assert_pair_record_count emi
# Both mates must EMIT on the ALT contig (that is a[which]) while carrying an
# XS from the weak primary hit (that is a[0], the one the default reads).
if [[ "$(grep -c "^$PAIR	.*	$ALT_CONTIG	" def.sam)" -ne 2 ]]; then
    echo "FAIL: $PAIR did not emit both mates on $ALT_CONTIG:" >&2
    grep "^$PAIR	" def.sam | cut -f1-6 >&2
    exit 1
fi
# Unanchored on the right on purpose. `XS:i:25` is a prefix of `XS:i:25X`, but
# the reads are 120 bp at a match score of 1, so no reachable XS exceeds 120 and
# 25 is the only value this can match. Anchoring to a trailing tab instead would
# pin the assertion to `HN:i` continuing to follow `XS:i`, which is a tag-order
# assumption this check has no reason to make.
if [[ "$(grep -c "^$PAIR	.*	XS:i:25" def.sam)" -ne 2 ]]; then
    echo "FAIL: $PAIR lost its sub-T (25 < 30) primary hit, so which == 0 and" >&2
    echo "      the a[0]-vs-a[which] branch is no longer exercised:" >&2
    grep "^$PAIR	" def.sam | cut -f1-6,12- >&2
    exit 1
fi
echo "PASS: fixture reaches the branch (ALT-emitted pair over a sub-T primary hit)"

# --- 2. The bit differs between the modes, by exactly 0x2. ----------------
mapfile -t def_flags < <(flags_of def)
mapfile -t emi_flags < <(flags_of emi)
for i in 0 1; do
    if ((def_flags[i] & 2)); then
        echo "FAIL: mate $i has FLAG 0x2 set by DEFAULT (${def_flags[i]}); the" >&2
        echo "      default must derive it from a[0], whose hits are outside" >&2
        echo "      the insert window" >&2
        exit 1
    fi
    if ! ((emi_flags[i] & 2)); then
        echo "FAIL: mate $i lacks FLAG 0x2 under --proper-pair-from-emitted" >&2
        echo "      (${emi_flags[i]}); the emitted ALT hits are 300 bp apart," >&2
        echo "      inside the window" >&2
        exit 1
    fi
    # Nothing else in the FLAG may move: this option selects a source for one
    # bit, and a change that also flipped, say, 0x100 would be a different bug.
    if (((def_flags[i] ^ emi_flags[i]) != 2)); then
        echo "FAIL: mate $i FLAGs differ by more than 0x2:" >&2
        echo "      default ${def_flags[i]}, emitted ${emi_flags[i]}" >&2
        exit 1
    fi
done
echo "PASS: --proper-pair-from-emitted flips FLAG 0x2, and only 0x2, on the ALT pair"

# --- 3. Nothing else in the output moves. --------------------------------
# The pair above is the ENTIRE difference. @PG is excluded because its CL:
# records the argv, which necessarily differs between the two runs.
if ! diff <(grep -v '^@PG' def.sam | grep -v "^$PAIR	") \
    <(grep -v '^@PG' emi.sam | grep -v "^$PAIR	") > /dev/null; then
    echo "FAIL: the option moved records other than $PAIR:" >&2
    diff <(grep -v '^@PG' def.sam | grep -v "^$PAIR	") \
        <(grep -v '^@PG' emi.sam | grep -v "^$PAIR	") | head -20 >&2
    exit 1
fi
# ...and the two mates differ in the FLAG field alone, not in position, MAPQ,
# CIGAR, or any tag.
if ! diff <(grep "^$PAIR	" def.sam | cut -f1,3-) \
    <(grep "^$PAIR	" emi.sam | cut -f1,3-) > /dev/null; then
    echo "FAIL: $PAIR differs outside the FLAG field:" >&2
    diff <(grep "^$PAIR	" def.sam | cut -f1,3-) \
        <(grep "^$PAIR	" emi.sam | cut -f1,3-) >&2
    exit 1
fi
echo "PASS: the ALT pair's FLAG is the only difference between the two modes"

# --- 4. Without a `.alt` sidecar the option is inert. ---------------------
# The structural claim the default rests on -- `which != 0` requires ALT hits,
# and is_alt is never set without a sidecar. Same reads, same sequence, no
# sidecar: the two modes must be byte-identical, ALT pair included. `@PG` is
# excluded for the same reason as in section 3 -- its CL: field records the
# argv, so it differs whenever the flag is passed, on any reference. That is
# the one record the flag always moves; the claim here is over everything else.
run noalt_def noalt.fasta
run noalt_emi noalt.fasta --proper-pair-from-emitted
if ! diff <(grep -v '^@PG' noalt_def.sam) <(grep -v '^@PG' noalt_emi.sam) > /dev/null; then
    echo "FAIL: --proper-pair-from-emitted changed output on a reference with no" >&2
    echo "      .alt sidecar; the option is documented as inert there" >&2
    diff <(grep -v '^@PG' noalt_def.sam) <(grep -v '^@PG' noalt_emi.sam) | head -20 >&2
    exit 1
fi
echo "PASS: with no .alt sidecar the option is inert (byte-identical output apart from @PG)"

echo "PASS: proper-pair ALT derivation regression"
