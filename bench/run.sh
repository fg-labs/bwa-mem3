#!/usr/bin/env bash
# Run bwa-mem3 mem on a fixed harness input and append measurements to results.csv.
# Usage: bench/run.sh <tag>
#   <tag>  a label for this run (e.g. "baseline", "pr-15", "main-after-lto").
#
# Reads bench/config.env. Runs:
#   - 1 single-thread trial (for golden md5)
#   - BENCH_TRIALS multi-thread trials (for wall-clock + RSS)
#
# Each trial appends one CSV row to BENCH_RESULTS.

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 <tag>" >&2
    exit 2
fi

TAG="$1"
HERE="$(cd "$(dirname "$0")" && pwd)"

# Load config.
CONFIG="$HERE/config.env"
if [[ ! -f "$CONFIG" ]]; then
    echo "error: $CONFIG not found. copy from config.env.example and edit." >&2
    exit 2
fi
# shellcheck disable=SC1090
. "$CONFIG"

# Resolve binary by tag using an explicit allow-list. Unknown tags fail
# fast (rather than silently picking the candidate path) so a typo like
# "baseliine" can't quietly run with the wrong binary.
case "$TAG" in
    baseline) BIN="${BWA_MEM3_BASELINE:?BWA_MEM3_BASELINE is required (set in bench/config.env)}" ;;
    candidate | pr-* | perf-* | main-*) BIN="${BWA_MEM3_CANDIDATE:?BWA_MEM3_CANDIDATE is required (set in bench/config.env)}" ;;
    *)
        echo "error: unknown TAG '$TAG' — expected 'baseline', 'candidate', 'pr-*', 'perf-*', or 'main-*'." >&2
        exit 2
        ;;
esac
echo "TAG=$TAG -> BIN=$BIN" >&2

# BENCH_INDEX is a bwa-mem3 index *prefix* (the loader probes for
# <prefix>.0123, <prefix>.bwt.2bit.64, etc.), so the bare prefix may not
# itself exist as a file. Probe the .0123 sidecar — it's the only one
# bwa-mem3 mem unconditionally requires — but report the prefix name so
# the user sees the variable they configured.
[[ -e "$BIN" ]] || {
    echo "error: missing path: $BIN" >&2
    exit 2
}
[[ -e "${BENCH_INDEX}.0123" ]] || {
    echo "error: missing path: $BENCH_INDEX (.0123 sidecar not found)" >&2
    exit 2
}
[[ -e "$BENCH_R1" ]] || {
    echo "error: missing path: $BENCH_R1" >&2
    exit 2
}
[[ -e "$BENCH_R2" ]] || {
    echo "error: missing path: $BENCH_R2" >&2
    exit 2
}

mkdir -p "$BENCH_OUTDIR"

HOST="$(hostname -s)"
ARCH="$(uname -m)"
OS="$(uname -s)"

# Pick md5 tool per OS (macOS ships `md5`, Linux ships `md5sum`).
if command -v md5sum > /dev/null 2>&1; then
    md5_hash() { md5sum "$1" | awk '{print $1}'; }
    md5_stream() { md5sum | awk '{print $1}'; }
else
    md5_hash() { md5 -q "$1"; }
    md5_stream() { md5 -q; }
fi

# Pick /usr/bin/time invocation per OS. Both emit max RSS; units differ.
# shellcheck disable=SC2016  # $1/$2/$NF below are awk fields, not shell vars
if [[ "$OS" == "Darwin" ]]; then
    TIME_CMD=(/usr/bin/time -l)
    # macOS time -l:
    #   "        1.01 real         0.00 user         0.00 sys"
    #   "             5767168  maximum resident set size"
    # RSS is in bytes; divide by 1024 → kB.
    WALL_AWK='$2=="real" && $4=="user" && $6=="sys" {print $1; exit}'
    RSS_AWK='/maximum resident set size/ {print $1; exit}'
    RSS_SCALE=1024
else
    TIME_CMD=(/usr/bin/time -v)
    # GNU time -v:
    #   "        Elapsed (wall clock) time (h:mm:ss or m:ss): 1:23.45"
    #   "        Maximum resident set size (kbytes): 4567890"
    WALL_AWK='/Elapsed \(wall clock\) time/ {print $NF; exit}'
    RSS_AWK='/Maximum resident set size/ {print $NF; exit}'
    RSS_SCALE=1
fi

# Ensure results header exists.
RESULTS="$BENCH_RESULTS"
if [[ ! -f "$RESULTS" ]]; then
    mkdir -p "$(dirname "$RESULTS")"
    echo "tag,host,arch,binary,threads,trial,wall_s,max_rss_kb,md5" > "$RESULTS"
fi

BIN_FP="$(md5_hash "$BIN" | cut -c1-12)"

run_trial() {
    local trial="$1" threads="$2"
    local time_file="$BENCH_OUTDIR/time.$TAG.$trial.$$"
    # Only the golden trial consumes the SAM (md5 below). Perf trials
    # ("perf$i") would otherwise pay full SAM-write I/O — for the smoke1M
    # PE150 input that's hundreds of MB per trial, folding filesystem-write
    # latency (and page-cache state across runs) into wall_s and masking
    # exactly the CPU-side wins this harness is designed to measure.
    # Redirect non-golden trials to /dev/null while keeping the SAM-build
    # code path live.
    local sam_file
    if [[ "$trial" == "golden" ]]; then
        sam_file="$BENCH_OUTDIR/sam.$TAG.$trial.$$.sam"
    else
        sam_file="/dev/null"
    fi
    local md5=""

    # Drop the time_file (and the golden SAM, if created) if the trial aborts
    # before reaching the explicit cleanup below.
    #
    # RETURN on its own is not enough. It fires for normal returns and for
    # errors under `set -e`, but a SIGINT or SIGTERM arriving during the timed
    # `mem` run kills the shell without ever running it, stranding a
    # hundreds-of-MB golden SAM in BENCH_OUTDIR. EXIT covers the remaining
    # abort paths; INT and TERM clean up and then re-raise on a cleared trap so
    # the caller still observes a signal death rather than a plain exit.
    #
    # The trap body is a string eval'd by the shell when the trap fires, so
    # paths are escaped via printf %q (which produces a shell-safe quoting of
    # arbitrary bytes) rather than wrapped in single quotes — single quotes
    # don't survive an embedded apostrophe in BENCH_OUTDIR.
    local time_file_q sam_file_q rm_cmd
    printf -v time_file_q '%q' "$time_file"
    printf -v sam_file_q '%q' "$sam_file"
    if [[ "$sam_file" != "/dev/null" ]]; then
        rm_cmd="rm -f -- $time_file_q $sam_file_q"
    else
        rm_cmd="rm -f -- $time_file_q"
    fi
    # shellcheck disable=SC2064  # bake the quoted paths in now, not at trap time
    trap "$rm_cmd" RETURN EXIT
    # shellcheck disable=SC2064
    trap "$rm_cmd; trap - INT; kill -INT \$\$" INT
    # shellcheck disable=SC2064
    trap "$rm_cmd; trap - TERM; kill -TERM \$\$" TERM

    # Run. stdout → SAM (or /dev/null); stderr → time_file (includes bwa-mem3 log + /usr/bin/time output).
    # shellcheck disable=SC2086
    "${TIME_CMD[@]}" "$BIN" mem -t "$threads" "$BENCH_INDEX" "$BENCH_R1" "$BENCH_R2" \
        > "$sam_file" 2> "$time_file"

    # Wall clock.
    local wall
    wall="$(awk "$WALL_AWK" "$time_file")"
    # GNU time may emit h:mm:ss or m:ss — normalize.
    case "$wall" in
        *:*:*) wall=$(awk -F: '{print $1*3600+$2*60+$3}' <<< "$wall") ;;
        *:*) wall=$(awk -F: '{print $1*60+$2}' <<< "$wall") ;;
    esac
    [[ -n "$wall" ]] || wall=NA

    # Max RSS.
    local rss_raw rss_kb
    rss_raw="$(awk "$RSS_AWK" "$time_file")"
    if [[ -n "$rss_raw" ]]; then
        rss_kb=$((rss_raw / RSS_SCALE))
    else
        rss_kb=NA
    fi

    # Refuse to record a trial that exited 0 having aligned nothing.
    #
    # Perf trials send their SAM to /dev/null, so the stream itself cannot be
    # inspected without putting a consumer back into the timed path — which is
    # the very I/O the redirect above exists to keep out of wall_s. The
    # aligner's own progress log is already captured in $time_file for the
    # wall/RSS parse, so count the reads it reports processing instead: that
    # costs nothing, perturbs no measurement, and covers both trial kinds.
    local processed
    processed="$(awk '/Processed [0-9]+ reads/ {
        for (i = 1; i <= NF; i++) if ($i == "Processed") { n += $(i + 1); break }
    } END {print n + 0}' "$time_file")"
    if [[ "$processed" -le 0 ]]; then
        echo "error: trial '$trial' exited 0 but processed no reads; refusing to record it" >&2
        tail -20 "$time_file" >&2
        return 1
    fi

    # Golden md5 only on single-thread trial.
    if [[ "$trial" == "golden" ]]; then
        md5="$(grep -v '^@PG' "$sam_file" | md5_stream)"
    fi

    echo "$TAG,$HOST,$ARCH,$BIN_FP,$threads,$trial,$wall,$rss_kb,$md5" | tee -a "$RESULTS"

    if [[ "${BENCH_KEEP_LOGS:-0}" == "1" ]]; then
        mv "$time_file" "$BENCH_OUTDIR/time.$TAG.$trial.kept" || true
    else
        rm -f "$time_file"
    fi
    if [[ "$sam_file" != "/dev/null" ]]; then
        rm -f "$sam_file"
    fi

    # Bash traps are global (no functrace / `local -` here), so leaving these
    # set would clobber any caller's handlers and persist into the next
    # run_trial invocation. Clear them now that cleanup is done; the next trial
    # registers its own against its own paths.
    trap - RETURN EXIT INT TERM
}

echo "# run.sh tag=$TAG bin=$BIN host=$HOST arch=$ARCH" >&2

# 1 golden, single-thread.
run_trial golden 1

# BENCH_TRIALS perf trials, multi-thread.
for i in $(seq 1 "$BENCH_TRIALS"); do
    run_trial "perf$i" "$BENCH_THREADS"
done
