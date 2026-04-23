#!/usr/bin/env bash
# Run bwa-mem2 mem on a fixed harness input and append measurements to results.csv.
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

# Resolve binary by tag. "baseline" -> BWA_MEM2_BASELINE; anything else -> BWA_MEM2_CANDIDATE.
case "$TAG" in
  baseline) BIN="$BWA_MEM2_BASELINE" ;;
  *)        BIN="$BWA_MEM2_CANDIDATE" ;;
esac

for path in "$BIN" "$BENCH_INDEX" "$BENCH_R1" "$BENCH_R2"; do
  [[ -e "$path" ]] || { echo "error: missing path: $path" >&2; exit 2; }
done

mkdir -p "$BENCH_OUTDIR"

HOST="$(hostname -s)"
ARCH="$(uname -m)"
OS="$(uname -s)"

# Pick /usr/bin/time invocation per OS. Both emit max RSS; units differ.
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

BIN_FP="$(md5sum "$BIN" | awk '{print $1}' | cut -c1-12)"

run_trial() {
  local trial="$1" threads="$2"
  local time_file="$BENCH_OUTDIR/time.$TAG.$trial.$$"
  local sam_file="$BENCH_OUTDIR/sam.$TAG.$trial.$$.sam"
  local md5=""

  # Warm filesystem cache on first trial only (optional — comment out for cold-cache runs).
  # shellcheck disable=SC2034
  local _warm_done=0

  # Run. stdout → SAM; stderr → time_file (includes bwa-mem2 log + /usr/bin/time output).
  # shellcheck disable=SC2086
  "${TIME_CMD[@]}" "$BIN" mem -t "$threads" "$BENCH_INDEX" "$BENCH_R1" "$BENCH_R2" \
      > "$sam_file" 2> "$time_file"

  # Wall clock.
  local wall
  wall="$(awk "$WALL_AWK" "$time_file")"
  # GNU time may emit h:mm:ss or m:ss — normalize.
  case "$wall" in
    *:*:*) wall=$(awk -F: '{print $1*3600+$2*60+$3}' <<<"$wall") ;;
    *:*)   wall=$(awk -F: '{print $1*60+$2}' <<<"$wall") ;;
  esac
  [[ -n "$wall" ]] || wall=NA

  # Max RSS.
  local rss_raw rss_kb
  rss_raw="$(awk "$RSS_AWK" "$time_file")"
  if [[ -n "$rss_raw" ]]; then
    rss_kb=$(( rss_raw / RSS_SCALE ))
  else
    rss_kb=NA
  fi

  # Golden md5 only on single-thread trial.
  if [[ "$trial" == "golden" ]]; then
    md5="$(grep -v '^@PG' "$sam_file" | md5sum | awk '{print $1}')"
  fi

  echo "$TAG,$HOST,$ARCH,$BIN_FP,$threads,$trial,$wall,$rss_kb,$md5" | tee -a "$RESULTS"

  if [[ "${BENCH_KEEP_LOGS:-0}" == "1" ]]; then
    mv "$time_file" "$BENCH_OUTDIR/time.$TAG.$trial.kept" || true
  else
    rm -f "$time_file"
  fi
  rm -f "$sam_file"
}

echo "# run.sh tag=$TAG bin=$BIN host=$HOST arch=$ARCH" >&2

# 1 golden, single-thread.
run_trial golden 1

# BENCH_TRIALS perf trials, multi-thread.
for i in $(seq 1 "$BENCH_TRIALS"); do
  run_trial "perf$i" "$BENCH_THREADS"
done
