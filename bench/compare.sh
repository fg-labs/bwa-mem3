#!/usr/bin/env bash
# Compare two tags in bench/results.csv.
# Usage: bench/compare.sh <tag_a> <tag_b>
# Emits a short report: golden md5 equality, median wall delta, median RSS delta.

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: $0 <tag_a> <tag_b>" >&2
  exit 2
fi

TAG_A="$1"; TAG_B="$2"
HERE="$(cd "$(dirname "$0")" && pwd)"

# Load config for BENCH_RESULTS path.
if [[ ! -f "$HERE/config.env" ]]; then
  echo "no config at $HERE/config.env (copy from config.env.example)" >&2
  exit 2
fi
# shellcheck disable=SC1091
. "$HERE/config.env"

RESULTS="$BENCH_RESULTS"
[[ -f "$RESULTS" ]] || { echo "no results at $RESULTS" >&2; exit 2; }

# Pass tag values through the environment (ENVIRON["TAG"]) instead of
# `awk -v t="$1"` — the latter interprets backslash escapes in the value
# (\n, \t, etc.), which would silently corrupt comparisons if a tag ever
# contains a backslash. ENVIRON treats values as literal bytes.
md5_for() {
  TAG="$1" awk -F, '$1==ENVIRON["TAG"] && $6=="golden" {print $9}' "$RESULTS" | tail -1
}

# Extract median wall-clock over perf trials (simple sort-middle).
# NR>1 skips the CSV header row; the prior `$6!="trial"` predicate did the
# same thing implicitly (the header's 6th column is the literal string
# "trial") but read like a trial-name filter at the call site.
# bench/run.sh writes `NA` to wall_s / max_rss_kb when /usr/bin/time output
# can't be parsed. `sort -n` collates `NA` as 0, so a single failed trial
# would silently drag the median toward zero — exactly the kind of
# attribution bug this harness exists to avoid. Drop NA rows before sort.
median_wall() {
  TAG="$1" awk -F, 'NR>1 && $1==ENVIRON["TAG"] && $6!="golden" && $7!="NA" {print $7}' "$RESULTS" \
    | sort -n \
    | awk '{a[NR]=$1} END{if(NR==0){print "NA"; exit} print (NR%2==1?a[(NR+1)/2]:(a[NR/2]+a[NR/2+1])/2)}'
}

median_rss() {
  TAG="$1" awk -F, 'NR>1 && $1==ENVIRON["TAG"] && $6!="golden" && $8!="NA" {print $8}' "$RESULTS" \
    | sort -n \
    | awk '{a[NR]=$1} END{if(NR==0){print "NA"; exit} print (NR%2==1?a[(NR+1)/2]:(a[NR/2]+a[NR/2+1])/2)}'
}

MD5_A="$(md5_for "$TAG_A")"
MD5_B="$(md5_for "$TAG_B")"
WA="$(median_wall "$TAG_A")"
WB="$(median_wall "$TAG_B")"
RA="$(median_rss "$TAG_A")"
RB="$(median_rss "$TAG_B")"

echo "tag A: $TAG_A"
echo "tag B: $TAG_B"
echo ""
echo "golden md5 A: $MD5_A"
echo "golden md5 B: $MD5_B"
if [[ -z "$MD5_A" || -z "$MD5_B" ]]; then
  missing=()
  [[ -z "$MD5_A" ]] && missing+=("$TAG_A")
  [[ -z "$MD5_B" ]] && missing+=("$TAG_B")
  echo "golden: NA — no golden md5 recorded for: ${missing[*]}"
elif [[ "$MD5_A" == "$MD5_B" ]]; then
  echo "golden: MATCH (byte-identical SAM)"
else
  echo "golden: MISMATCH — output changed"
fi
echo ""
echo "median wall A: ${WA}s"
echo "median wall B: ${WB}s"
# Numeric zero check inside awk so values like "0.0" or "0e0" don't slip
# past a literal-string `[[ "$WA" != "0" ]]` and trip a divide-by-zero.
if [[ "$WA" != "NA" && "$WB" != "NA" ]]; then
  awk -v a="$WA" -v b="$WB" 'BEGIN{
    if (a+0 == 0) { print "wall delta: n/a (A baseline is 0)"; exit }
    printf "wall delta: %+.2f%% (B vs A)\n", (b-a)/a*100
  }'
fi
echo ""
echo "median RSS A: ${RA} kb"
echo "median RSS B: ${RB} kb"
if [[ "$RA" != "NA" && "$RB" != "NA" ]]; then
  awk -v a="$RA" -v b="$RB" 'BEGIN{
    if (a+0 == 0) { print "RSS delta:  n/a (A baseline is 0)"; exit }
    printf "RSS delta:  %+.2f%% (B vs A)\n", (b-a)/a*100
  }'
fi
