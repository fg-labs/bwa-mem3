#!/usr/bin/env bash
# test/regression/profile_slice_cpu.sh
#
# Regression: --profile must account for the compute CPU of a PARTIAL cohort.
#
# Cohort slicing gave step 1 a second exit. A cohort that is still accumulating
# slices returns early -- it has seeded and extended its slice, but pairing and
# SAM must wait for the whole cohort -- and that path originally saved only
# proc_wall. Because g_ktfor is reset at every step-1 entry, the slice's
# proc_cpu / compute / encode / thr_busy_* were not merely omitted from the row,
# they were destroyed: nothing else ever reads them.
#
# The consequence is a profiler that quietly lies in proportion to how much
# slicing is happening. At -t 64 the ramp covers ~40% of the run's bases, so
# ~40% of the compute CPU vanished from the report -- and since the aggregate
# row derives from the per-chunk rows, the in-PROCESS() budget could not be made
# to balance. Nothing crashes and no alignment changes, so only an explicit
# check on the profile output catches it.
#
# What makes this test sharp rather than merely present:
#
#   1. BWA_MEM3_COHORT_SLICE_ALL=1 slices EVERY cohort rather than only the
#      first (which is all production does), so a short run produces dozens of
#      partial slices instead of a handful.
#
#   2. The assertion is a boolean, not a threshold: every chunk that read bases
#      must report proc_cpu > 0. A slice that ran kt_for cannot legitimately
#      have consumed zero CPU, so there is nothing to tune and nothing to flake.
#      Before the fix every ramp slice reports exactly 0.0000.
#
# Inputs (env vars):
#   BWA_MEM3 -- path to a bwa-mem3 binary built with STAGE_PROF=1 (so --profile
#               exists). A default build compiles profiling out entirely.
#   FIXTURES -- directory containing synthetic_1mb.fa (default: test/fixtures)
set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
FIXTURES="${FIXTURES:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../fixtures" && pwd)}"

[[ -x "$BWA_MEM3" ]] || { echo "FAIL: binary not executable: $BWA_MEM3" >&2; exit 1; }

src_ref="$FIXTURES/synthetic_1mb.fa"
[[ -s "$src_ref" ]] || { echo "FAIL: synthetic_1mb.fa missing: $src_ref" >&2; exit 1; }

WORK=$(mktemp -d "${TMPDIR:-/tmp}/profile_slice_cpu.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

# --profile only exists in a STAGE_PROF=1 build; without it there is nothing to
# assert on, and silently passing would make this test worthless in CI.
# `mem` with no arguments prints usage and exits nonzero, so the output has to be
# captured before grepping -- piping it straight into grep trips pipefail.
"$BWA_MEM3" mem > "$WORK/usage.txt" 2>&1 || true
if ! grep -q -- '--profile' "$WORK/usage.txt"; then
    echo "FAIL: this binary has no --profile; build with STAGE_PROF=1" >&2
    exit 1
fi

ref="$WORK/ref.fa"
cp "$src_ref" "$ref"
"$BWA_MEM3" index "$ref" > "$WORK/index.log" 2>&1 \
    || { echo "FAIL: could not index the reference" >&2; tail -20 "$WORK/index.log" >&2; exit 1; }

# Reads are generated here rather than committed. They must come FROM the
# reference so they align and the aligner does real per-slice work -- a slice
# whose reads all fail to seed would report near-zero CPU for honest reasons.
READS="$WORK/reads.fa"
READ_LEN=100
N_READS=20000
awk -v len="$READ_LEN" -v n="$N_READS" '
    /^>/  { next }
          { seq = seq $0 }
    END {
        if (length(seq) < len * 2) { print "reference too short" > "/dev/stderr"; exit 1 }
        # Walk the reference with a stride, wrapping around, so the reads are
        # spread over the whole sequence instead of piling onto one locus.
        span   = length(seq) - len
        stride = 7
        pos    = 1
        for (i = 0; i < n; i++) {
            printf(">r%d\n%s\n", i, substr(seq, pos, len))
            pos += stride
            if (pos > span) pos = (pos % span) + 1
        }
    }
' "$ref" > "$READS"
[[ -s "$READS" ]] || { echo "FAIL: read generation produced nothing" >&2; exit 1; }

# -K is deliberately small so the input spans many cohorts, each of which
# SLICE_ALL then splits -- exercising the partial-cohort exit repeatedly.
K=200000

align() {   # $1 = tag, $2 = --cohort-slices, $3 = slice-every-cohort (0/1)
    local tag="$1" slices="$2" all="$3"
    BWA_MEM3_COHORT_SLICE_ALL="$all" \
        "$BWA_MEM3" mem -t 4 -K "$K" --cohort-slices "$slices" \
        --profile "$WORK/$tag.tsv" "$ref" "$READS" \
        > "$WORK/$tag.sam" 2> "$WORK/$tag.err" \
        || { echo "FAIL: $tag exited nonzero" >&2; tail -20 "$WORK/$tag.err" >&2; exit 1; }
    [[ -s "$WORK/$tag.tsv" ]] || { echo "FAIL: $tag wrote no profile TSV" >&2; exit 1; }
}

align unsliced 0 0
align sliced   4 1

# Resolve columns by header name: the TSV schema is append-friendly and a
# hard-coded index would rot the next time a field is added.
#
# A missing header is fatal rather than an empty index, because an empty index
# does NOT reliably fail: awk reads `$""` as `$0`, so `$bp+0 > 0` silently
# degrades into "the whole line coerced to a number", and every assertion below
# would pass on a broken build. It is also implementation-dependent -- mawk and
# gawk (what CI runs) take the `$0` reading, while the BWK awk on macOS errors
# out -- so checking here is what makes the behaviour the same everywhere.
col() {   # $1 = tsv, $2 = header name -> 1-based column index; fails if absent
    local idx
    idx=$(awk -F'\t' -v want="$2" 'NR==1 { for (i=1;i<=NF;i++) if ($i==want) { print i; exit } }' "$1")
    if [[ -z "$idx" ]]; then
        echo "FAIL: $1 has no '$2' column; header is: $(head -1 "$1")" >&2
        return 1
    fi
    printf '%s\n' "$idx"
}

# Rows that read bases (n_bp > 0), excluding the ALL aggregate.
data_rows() {   # $1 = tsv
    local f="$1" c_chunk c_bp
    c_chunk=$(col "$f" chunk) || return 1
    c_bp=$(col "$f" n_bp)     || return 1
    awk -F'\t' -v ch="$c_chunk" -v bp="$c_bp" \
        'NR>1 && $ch!="ALL" && $bp+0 > 0' "$f"
}

n_unsliced=$(data_rows "$WORK/unsliced.tsv" | wc -l | tr -d ' ')
n_sliced=$(data_rows "$WORK/sliced.tsv" | wc -l | tr -d ' ')
partials=$(grep -c 'cohort slice' "$WORK/sliced.err" || true)

echo "  unsliced: $n_unsliced chunk rows"
echo "  sliced:   $n_sliced chunk rows, $partials partial slice(s) reported by the reader"

# Non-vacuity. If slicing did not actually split any cohort then every row takes
# the cohort-complete exit, the partial-cohort path is never entered, and the
# assertion below would pass on a broken build.
if [ "$n_sliced" -le "$n_unsliced" ] || [ "${partials:-0}" -lt 3 ]; then
    echo "FAIL: slicing did not produce partial cohorts, so this run cannot" >&2
    echo "      detect the bug. Expected more rows than the $n_unsliced" >&2
    echo "      unsliced rows and at least 3 partial slices, got $n_sliced" >&2
    echo "      rows and ${partials:-0} partial slices." >&2
    exit 1
fi

# The assertion: a chunk that read bases ran kt_for, so it cannot have used
# zero compute CPU. `compute` is derived from proc_cpu, so it must be populated
# too (an empty cell is NaN -- see sp_chunk_init).
check() {   # $1 = tsv, $2 = tag
    local f="$1" tag="$2" c_chunk c_bp c_cpu c_comp
    c_chunk=$(col "$f" chunk)    || return 1
    c_bp=$(col "$f" n_bp)        || return 1
    c_cpu=$(col "$f" proc_cpu)   || return 1
    c_comp=$(col "$f" compute)   || return 1
    awk -F'\t' -v ch="$c_chunk" -v bp="$c_bp" -v cpu="$c_cpu" -v comp="$c_comp" -v tag="$tag" '
        NR>1 && $ch!="ALL" && $bp+0 > 0 {
            rows++
            if ($cpu+0 <= 0) { printf("  %s chunk %s: %s bases but proc_cpu=%s\n", tag, $ch, $bp, $cpu); bad++ }
            else if ($comp == "")   { printf("  %s chunk %s: proc_cpu=%s but compute is empty\n", tag, $ch, $cpu); bad++ }
            else total += $cpu
        }
        END {
            printf("  %s: %d/%d rows report compute CPU, %.4f CPU sec total\n", tag, rows-bad, rows, total)
            exit (bad > 0)
        }
    ' "$f"
}

fail=0
check "$WORK/unsliced.tsv" unsliced || fail=1
check "$WORK/sliced.tsv"   sliced   || fail=1

if [ "$fail" -ne 0 ]; then
    # Two causes reach here, and the diagnostics above say which: an unresolved
    # TSV column (a schema change -- col() names the missing header), or the
    # regression itself. Don't attribute the former to the latter.
    echo "FAIL: --profile did not account for every chunk that read bases." >&2
    echo "      If a column could not be resolved, the TSV schema changed and" >&2
    echo "      this script needs updating. Otherwise a partial cohort failed to" >&2
    echo "      harvest g_ktfor before returning; see sp_harvest_proc() and" >&2
    echo "      step 1 of kt_pipeline (src/fastmap.cpp)." >&2
    exit 1
fi

echo "PASS: every chunk that read bases reports its compute CPU, including the"
echo "      partial cohort slices that return before pairing and SAM."
