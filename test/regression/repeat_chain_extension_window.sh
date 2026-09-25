#!/usr/bin/env bash
set -euo pipefail
# test/regression/repeat_chain_extension_window.sh
#
# Regression: extension staging must not copy a chain's whole reference window
# once PER SEED.
#
# A short read from a long homopolymer has an SMEM with tens of thousands of
# hits. Seeding keeps a max_occ (500) sample of them, spaced about
# (run length / 500) apart. When that spacing is within the band (-w 100), the
# samples all chain, giving one chain of ~500 seeds whose reference window spans
# the whole run. Each seed's left extension targets the window prefix up to the
# seed, and its right extension the suffix after it. Staging a private copy of
# those per seed costs ~500 x (run length) bytes per read, so one batch of a few
# hundred such 150 bp reads grew the per-thread extension buffers past the
# int32-addressable limit and the run aborted with
#   ERROR: seqBufRef in stage_seed_extension cannot grow past the int32 offset range
# The chain's left targets are all suffixes of one reversed prefix, and its
# right targets all suffixes of one forward suffix, so they are now staged once
# per chain.
#
# The fixture is a 60 kb reference: 20 kb random, a 20 kb poly-A run, 20 kb
# random. The reads are 600 poly-A 150-mers, each with up to three G/T
# substitutions (positions are drawn independently and may coincide) so that no
# single SMEM spans the read and every seed needs extension. Both the
# default path and --meth are checked (--meth doubles the hits across its two
# converted strands and disables the ungapped fast path, so it stages every
# seed). Each run must exit 0 and emit one primary record per read, mapped
# on the poly-A run. Each run also counts the banded-extension jobs
# (BWAMEM3_DEDUP_STATS) and requires at least 100 per read: a read on a short
# chain stages a handful, so this is what shows the fixture still reaches the
# ~500-seed chains the bug needs, rather than passing because the chains got
# short.
#
# Not covered here: the read-side (query) window sharing. 150 bp reads cannot
# fill the query buffers even with per-seed copies; that takes kilobase reads,
# which align too slowly on this fixture for a CI test.
#
# Inputs:
#   BWA_MEM3 - the bwa-mem3 binary under test: a bare command name resolved via
#              PATH, or a path to the binary (a relative path is made absolute
#              before the script changes directory).
if [ -z "${BWA_MEM3:-}" ]; then
    echo "FAIL: BWA_MEM3 must be set" >&2
    exit 1
fi
case "$BWA_MEM3" in
    */*) BWA_MEM3="$(cd "$(dirname "$BWA_MEM3")" && pwd)/$(basename "$BWA_MEM3")" ;;
esac
command -v mawk > /dev/null 2>&1 || {
    echo "SKIP: mawk not on PATH (required to generate the fixture and inspect records)"
    exit 0
}

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
fail() {
    echo "FAIL: $*" >&2
    exit 1
}

FLANK=20000
RUN=20000
N_READS=600
READ_LEN=150

# Park-Miller LCG: exact in mawk's doubles (16807 * 2^31 < 2^53), so the
# fixture is identical on every host.
mawk -v flank="$FLANK" -v run="$RUN" 'BEGIN {
    x = 4242; split("A C G T", b, " ")
    s = ""
    for (i = 0; i < flank; i++) { x = (x * 16807) % 2147483647; s = s b[x % 4 + 1] }
    for (i = 0; i < run; i++) s = s "A"
    for (i = 0; i < flank; i++) { x = (x * 16807) % 2147483647; s = s b[x % 4 + 1] }
    print ">poly"
    for (i = 1; i <= length(s); i += 80) print substr(s, i, 80)
}' > ref.fa

mawk -v n="$N_READS" -v len="$READ_LEN" 'BEGIN {
    x = 1111; q = ""
    for (i = 0; i < len; i++) q = q "I"
    for (r = 0; r < n; r++) {
        for (i = 1; i <= len; i++) s[i] = "A"
        for (k = 0; k < 3; k++) {
            x = (x * 16807) % 2147483647; p = x % len + 1
            x = (x * 16807) % 2147483647; s[p] = (x % 2) ? "G" : "T"
        }
        seq = ""
        for (i = 1; i <= len; i++) seq = seq s[i]
        printf "@r%d\n%s\n+\n%s\n", r, seq, q
    }
}' > reads.fq

"$BWA_MEM3" index ref.fa > index.log 2>&1 || {
    tail -5 index.log >&2
    fail "index exited nonzero"
}
"$BWA_MEM3" index --meth ref.fa > index_meth.log 2>&1 || {
    tail -5 index_meth.log >&2
    fail "index --meth exited nonzero"
}

# Run one mode and check its records. $1 is a label; the rest are mem options.
check_mode() {
    local label="$1"
    shift
    if ! BWAMEM3_DEDUP_STATS=1 "$BWA_MEM3" mem "$@" -t 1 ref.fa reads.fq > "$label.sam" 2> "$label.err"; then
        grep -E '^ERROR' "$label.err" >&2 || tail -5 "$label.err" >&2
        fail "$label: mem exited nonzero on $N_READS homopolymer reads"
    fi
    # Primary records: count, how many are mapped, and how many of those
    # overlap the poly-A run (1-based [FLANK+1, FLANK+RUN]). A read may start a
    # few bases into the left flank when the flank happens to end in A.
    local counts
    counts=$(mawk -v lo=$((FLANK - READ_LEN + 2)) -v hi=$((FLANK + RUN)) '
        /^@/ { next }
        int($2 / 256) % 2 == 1 || int($2 / 2048) % 2 == 1 { next }
        { n++; if (int($2 / 4) % 2 == 0) { m++; if ($4 >= lo && $4 <= hi) inrun++ } }
        END { printf "%d %d %d\n", n, m, inrun }' "$label.sam")
    local n m inrun
    read -r n m inrun <<< "$counts"
    [ "$n" -eq "$N_READS" ] || fail "$label: expected $N_READS primary records, got $n"
    [ "$m" -eq "$N_READS" ] || fail "$label: expected all $N_READS reads mapped, got $m"
    [ "$inrun" -eq "$N_READS" ] || fail "$label: expected all reads placed on the poly-A run, got $inrun"
    # Extension jobs per read: the long-chain guard (see the header).
    local jobs
    jobs=$(sed -n 's/^\[dedup-stats\] total_jobs=\([0-9]*\).*/\1/p' "$label.err")
    [ -n "$jobs" ] || fail "$label: no [dedup-stats] line in stderr (BWAMEM3_DEDUP_STATS)"
    [ "$jobs" -ge $((N_READS * 100)) ] \
        || fail "$label: only $jobs extension jobs for $N_READS reads; the fixture no longer builds long chains"
    echo "  ok: $label ($n records, all mapped on the poly-A run, $jobs extension jobs)"
}

check_mode default
check_mode meth --meth

echo "PASS: repeat_chain_extension_window (a long homopolymer chain stages its reference window once, default and --meth)"
