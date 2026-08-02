#!/usr/bin/env bash
# test/regression/cohort_slice_identity.sh
#
# Regression: reading a batch in several slices must not change the output.
#
# A batch is also the mem_pestat cohort -- the read set the insert-size
# distribution is inferred from -- and those percentile bounds feed pairing,
# mate rescue and MAPQ. Slicing exists purely so compute can start before the
# whole first batch has been read; it must never move the batch BOUNDARY.
#
# The subtle failure this guards against: both readers stop at the first
# whole-record, even-n boundary at or PAST the requested size, so slices
# overshoot independently. Without clamping each slice's target to the bases the
# cohort still needs, T/8 + T/4 + T/2 + T/8 lands on a different final record
# than a single read of T -- a different cohort, a different mem_pestat, and
# quietly different alignments. The drift is a record or two per slice, so it
# shows up as neither a crash nor a record-count change; only a byte comparison
# catches it.
#
# Two things make this test actually powerful rather than merely present:
#
#   1. BWA_MEM3_COHORT_SLICE_ALL slices EVERY cohort, not just the first (which
#      is all production does). Measured against a build with the clamp removed:
#      slicing only the first cohort differs by 4 lines, slicing every cohort
#      differs by 24912. The same bug, ~6000x the signal.
#
#   2. The input must have a genuinely variable insert size. With a constant
#      insert, mem_pestat returns the same bounds for every possible subset of
#      reads, so cohort composition cannot affect output and the comparison is
#      vacuous no matter what the code does. The test asserts non-vacuity below
#      by requiring the inferred percentiles to differ between cohorts.
#
# Inputs (env vars):
#   BWA_MEM3       -- path to the bwa-mem3 binary under test
#   CHR22_FA       -- path to chr22.fa, pre-indexed with bwa-mem3 by the caller
#   CHR22_SIM_DIR  -- directory with reads.r1.fastq.gz / reads.r2.fastq.gz
set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${CHR22_FA:?CHR22_FA must be set}"
: "${CHR22_SIM_DIR:?CHR22_SIM_DIR must be set}"

R1="$CHR22_SIM_DIR/reads.r1.fastq.gz"
R2="$CHR22_SIM_DIR/reads.r2.fastq.gz"
[ -f "$R1" ] && [ -f "$R2" ] || {
    echo "FAIL: missing $R1 / $R2" >&2
    exit 1
}

WORK=$(mktemp -d "${TMPDIR:-/tmp}/cohort_slice.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

# -K is deliberately small so the input spans MANY cohorts: that exercises the
# accumulator, the cohort-id advance, the realloc growth path and the boundary
# clamp repeatedly rather than once.
K=2000000

# Partial slices are the reads that did NOT end their cohort; step 0 marks
# exactly those lines with a trailing " (cohort slice)".
slice_count() { # $1 = a run's stderr log
    grep -c '(cohort slice)$' "$1" || true
}

align() { # $1 = tag, $2 = --cohort-slices value, $3 = slice-every-cohort (0/1)
    local tag="$1" slices="$2" all="$3"
    # The aligner and the header strip run as separate steps so their exit
    # statuses stay distinguishable. Piped together, a run that emitted headers
    # but no records makes `grep -v` exit 1 (it selected nothing), and under
    # `set -o pipefail` that ends the script before the explicit "produced no
    # records" check below -- the failure this function most needs to report.
    BWA_MEM3_COHORT_SLICE_ALL="$all" \
        "$BWA_MEM3" mem -t 4 -K "$K" --cohort-slices "$slices" \
        "$CHR22_FA" "$R1" "$R2" 2> "$WORK/$tag.err" > "$WORK/$tag.raw"
    grep -v '^@' "$WORK/$tag.raw" > "$WORK/$tag.sam" || true
    local n
    n=$(wc -l < "$WORK/$tag.sam")
    [ "$n" -gt 0 ] || {
        echo "FAIL: $tag produced no records" >&2
        tail -15 "$WORK/$tag.err" >&2
        exit 1
    }
    echo "  $tag: $n records, $(grep -c 'Inferring insert size' "$WORK/$tag.err" || true) cohorts, $(slice_count "$WORK/$tag.err") partial slice(s)"
}

align unsliced 0 0
align first3 3 0
align all3 3 1
align all8 8 1

# The all-cohort runs have to actually slice more than the first cohort, or the
# comparisons below prove nothing about the path they exist to cover. If
# BWA_MEM3_COHORT_SLICE_ALL were ignored -- misspelled, dropped from the option
# block, or parsed as off -- all3/all8 would degenerate into first3/first8:
# still slicing cohort 0, still byte-identical, still passing every cmp.
unsliced_slices=$(slice_count "$WORK/unsliced.err")
first3_slices=$(slice_count "$WORK/first3.err")
all3_slices=$(slice_count "$WORK/all3.err")
all8_slices=$(slice_count "$WORK/all8.err")
if [ "$unsliced_slices" -ne 0 ]; then
    echo "FAIL: --cohort-slices 0 produced $unsliced_slices partial slice(s); it must produce none." >&2
    exit 1
fi
if [ "$first3_slices" -le 0 ]; then
    echo "FAIL: --cohort-slices 3 produced no partial slices, so nothing below is" >&2
    echo "      exercising the accumulator at all." >&2
    exit 1
fi
if [ "$all3_slices" -le "$first3_slices" ] || [ "$all8_slices" -le "$first3_slices" ]; then
    echo "FAIL: all-cohort slicing did not produce slices beyond the first cohort" >&2
    echo "      (first3=$first3_slices all3=$all3_slices all8=$all8_slices)." >&2
    echo "      BWA_MEM3_COHORT_SLICE_ALL is not taking effect, so the every-cohort" >&2
    echo "      path is untested. See the slice_all gate in kt_pipeline step 0." >&2
    exit 1
fi
echo "  slice-all effective: first3=$first3_slices all3=$all3_slices all8=$all8_slices partial slice(s)"

# Non-vacuity. If every cohort infers the same insert-size bounds then cohort
# composition cannot influence the output and an IDENTICAL result below would
# prove nothing about the clamp.
# `|| true` on both greps: zero "percentile" lines is precisely the vacuous
# fixture the check below exists to report, but grep exits 1 when it matches
# nothing (as does `grep -c` on a zero count), so pipefail would end the run
# before the diagnostic could print.
distinct=$({ grep -h "percentile" "$WORK/unsliced.err" || true; } | sort -u | wc -l | tr -d ' ')
total=$({ grep -hc "percentile" "$WORK/unsliced.err" || true; } | tr -d ' ')
echo "  non-vacuity: $distinct distinct percentile tuples across $total cohorts"
if [ "${distinct:-0}" -lt 2 ]; then
    echo "FAIL: mem_pestat inferred identical bounds for every cohort, so this" >&2
    echo "      fixture cannot detect cohort divergence. Use reads with a" >&2
    echo "      variable insert size." >&2
    exit 1
fi

fail=0
for t in first3 all3 all8; do
    if cmp -s "$WORK/unsliced.sam" "$WORK/$t.sam"; then
        echo "  unsliced vs $t: IDENTICAL"
    else
        echo "  unsliced vs $t: DIFFERS" >&2
        diff "$WORK/unsliced.sam" "$WORK/$t.sam" > "$WORK/$t.diff" 2>&1 || true
        echo "    differing lines: $(grep -c '^[<>]' "$WORK/$t.diff" || true)" >&2
        head -4 "$WORK/$t.diff" >&2
        fail=1
    fi
done

if [ "$fail" -ne 0 ]; then
    echo "FAIL: cohort slicing changed the output; the batch boundary moved." >&2
    echo "      See the slice-target clamp in kt_pipeline step 0 (src/fastmap.cpp)." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# EOF landing exactly ON a slice boundary.
#
# A partial slice does not end its cohort, so when the input runs out right after
# one, step 0 is asked for the next slice and gets zero reads while a cohort is
# still accumulating. Retiring the pipeline at that point discards every read
# already accumulated AND leaks the cohort's string arena -- silently, with exit
# status 0, so only a record count catches it.
#
# The full-length runs above cannot detect this: whichever slice consumes the
# last record is short, which ends its cohort normally. The input has to stop
# exactly where a slice stopped.
#
# The truncation point is taken from the binary's own slice log rather than
# computed, so it stays correct regardless of read length, -K or the ramp shape.
# ---------------------------------------------------------------------------
first_slice_reads=$(awk '/\(cohort slice\)$/ { sub(/.*nseq: /, ""); sub(/ .*/, ""); print; exit }' \
    "$WORK/all3.err")
if [ -z "${first_slice_reads:-}" ] || [ "$first_slice_reads" -lt 2 ]; then
    echo "FAIL: could not read the first partial slice's record count from the" >&2
    echo "      all3 log, so the EOF-on-a-slice-boundary case cannot be built." >&2
    exit 1
fi

pairs=$((first_slice_reads / 2))
expected=$((pairs * 2))
# `sed -n '1,Np'` rather than `head -n N`, so pipefail can stay on.
#
# head closes the pipe as soon as it has its lines; the SIGPIPE that kills the
# upstream gzip is what forced `set +o pipefail` here. But that window also hid
# a genuine decompression failure: a corrupt R1/R2 whose first N lines happen to
# inflate cleanly would still produce both fixtures, and this identity check
# would pass on truncated input -- exactly the outcome it exists to catch. sed
# drains the decompressor instead, so a real gzip failure surfaces and fails the
# fixture build loudly.
eof_lines=$((pairs * 4))
if ! gzip -cd "$R1" | sed -n "1,${eof_lines}p" | gzip > "$WORK/eof.r1.fastq.gz"; then
    echo "FAIL: could not create EOF R1 fixture" >&2
    exit 1
fi
if ! gzip -cd "$R2" | sed -n "1,${eof_lines}p" | gzip > "$WORK/eof.r2.fastq.gz"; then
    echo "FAIL: could not create EOF R2 fixture" >&2
    exit 1
fi

echo "  EOF-on-boundary: truncating to $pairs pairs ($expected records), the first"
echo "                   partial slice's exact size"

align_eof() { # $1 = tag, $2 = --cohort-slices value, $3 = slice-every-cohort (0/1)
    local tag="$1" slices="$2" all="$3" n rc
    set +e
    BWA_MEM3_COHORT_SLICE_ALL="$all" \
        "$BWA_MEM3" mem -t 4 -K "$K" --cohort-slices "$slices" \
        "$CHR22_FA" "$WORK/eof.r1.fastq.gz" "$WORK/eof.r2.fastq.gz" \
        > "$WORK/$tag.full.sam" 2> "$WORK/$tag.err"
    rc=$?
    set -e
    [ "$rc" -eq 0 ] || {
        echo "FAIL: $tag exited $rc" >&2
        tail -15 "$WORK/$tag.err" >&2
        exit 1
    }
    grep -v '^@' "$WORK/$tag.full.sam" > "$WORK/$tag.sam" || true
    n=$(wc -l < "$WORK/$tag.sam" | tr -d ' ')
    if [ "$n" -ne "$expected" ]; then
        echo "FAIL: $tag emitted $n records, expected $expected." >&2
        echo "      Reads accumulated in an in-flight cohort were dropped when the" >&2
        echo "      input ended on a slice boundary. See the zero-read/EOF early" >&2
        echo "      return in kt_pipeline step 0 (src/fastmap.cpp): it must still" >&2
        echo "      flush a cohort that step 1 is holding." >&2
        exit 1
    fi
    echo "  $tag: $n records"
}

align_eof eof_unsliced 0 0
align_eof eof_first3 3 0
align_eof eof_all3 3 1
align_eof eof_all8 8 1

for t in eof_first3 eof_all3 eof_all8; do
    if cmp -s "$WORK/eof_unsliced.sam" "$WORK/$t.sam"; then
        echo "  eof_unsliced vs $t: IDENTICAL"
    else
        echo "  eof_unsliced vs $t: DIFFERS" >&2
        diff "$WORK/eof_unsliced.sam" "$WORK/$t.sam" > "$WORK/$t.diff" 2>&1 || true
        echo "    differing lines: $(grep -c '^[<>]' "$WORK/$t.diff" || true)" >&2
        head -4 "$WORK/$t.diff" >&2
        fail=1
    fi
done

if [ "$fail" -ne 0 ]; then
    echo "FAIL: cohort slicing changed the output when the input ended exactly on" >&2
    echo "      a slice boundary." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# The same EOF-on-a-boundary input, now checking what the cohort RESERVED.
#
# Capacity never changes the output, so the identity checks above pass either
# way and cannot see this. The cohort accumulator projects its final read count
# from a slice's mean read length, which is only sound while more input is
# coming -- and on the first slice that is not yet known. A file ending exactly
# on the boundary returns a full-size slice with seqs != NULL, so it fails every
# cohort_complete test and reads as mid-cohort; EOF only surfaces on the next
# read. Projecting there reserves a whole task_size for a cohort that turns out
# to be one slice (measured: 40002 slots held 9766 reads at -K 4000000, and the
# ratio grows with task_size -- ~330 MB of never-used bseq1_t at a default t=64).
#
# Doubling alone never exceeds 2x what it holds, so that is the bound: it passes
# for any growth by doubling and fails for a task_size projection. -v 4 because
# the reservation line is opt-in; default verbosity is unchanged.
# ---------------------------------------------------------------------------
BWA_MEM3_COHORT_SLICE_ALL=1 \
    "$BWA_MEM3" mem -t 4 -K "$K" --cohort-slices 3 -v 4 \
    "$CHR22_FA" "$WORK/eof.r1.fastq.gz" "$WORK/eof.r2.fastq.gz" \
    > /dev/null 2> "$WORK/eof_reserve.err"

reserve_line=$(grep -m1 'cohort_reserve:' "$WORK/eof_reserve.err" || true)
if [ -z "$reserve_line" ]; then
    echo "FAIL: no 'cohort_reserve:' line under -v 4, so the reservation could" >&2
    echo "      not be checked. The accumulator must still log it." >&2
    exit 1
fi
res_cap=$(echo "$reserve_line" | sed -E 's/.*cap: ([0-9]+).*/\1/')
res_held=$(echo "$reserve_line" | sed -E 's/.*held: ([0-9]+).*/\1/')

if [ "$res_held" -ge 1024 ] && [ "$res_cap" -ge $(( res_held * 2 )) ]; then
    echo "FAIL: the cohort reserved $res_cap slots while holding $res_held reads." >&2
    echo "      Growth by doubling never exceeds 2x, so this is a task_size" >&2
    echo "      projection fired on a cohort whose input had already ended." >&2
    echo "      See the cohort_n > 0 gate in kt_pipeline step 1 (src/fastmap.cpp):" >&2
    echo "      a second slice is what proves more input exists." >&2
    exit 1
fi
echo "  EOF-on-boundary reservation: cap $res_cap for $res_held reads (under 2x)"

echo "PASS: output is identical with and without cohort slicing, including with"
echo "      every cohort sliced 8 ways, and when the input ends exactly on a"
echo "      slice boundary -- which also reserves no more than doubling would."
