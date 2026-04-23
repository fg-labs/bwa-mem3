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
# shellcheck disable=SC1091
. "$HERE/config.env"

RESULTS="$BENCH_RESULTS"
[[ -f "$RESULTS" ]] || { echo "no results at $RESULTS" >&2; exit 2; }

# Extract golden md5 for each tag (last occurrence wins, in case of reruns).
md5_for() {
  awk -F, -v t="$1" '$1==t && $6=="golden" {print $9}' "$RESULTS" | tail -1
}

# Extract median wall-clock over perf trials (simple sort-middle).
median_wall() {
  awk -F, -v t="$1" '$1==t && $6!="golden" && $6!="trial" {print $7}' "$RESULTS" \
    | sort -n \
    | awk '{a[NR]=$1} END{if(NR==0){print "NA"; exit} print (NR%2==1?a[(NR+1)/2]:(a[NR/2]+a[NR/2+1])/2)}'
}

median_rss() {
  awk -F, -v t="$1" '$1==t && $6!="golden" && $6!="trial" {print $8}' "$RESULTS" \
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
if [[ -n "$MD5_A" && "$MD5_A" == "$MD5_B" ]]; then
  echo "golden: MATCH (byte-identical SAM)"
else
  echo "golden: MISMATCH — output changed"
fi
echo ""
echo "median wall A: ${WA}s"
echo "median wall B: ${WB}s"
if [[ "$WA" != "NA" && "$WB" != "NA" ]]; then
  awk -v a="$WA" -v b="$WB" 'BEGIN{d=(b-a)/a*100; printf "wall delta: %+.2f%% (B vs A)\n", d}'
fi
echo ""
echo "median RSS A: ${RA} kb"
echo "median RSS B: ${RB} kb"
if [[ "$RA" != "NA" && "$RB" != "NA" && "$RA" != "0" ]]; then
  awk -v a="$RA" -v b="$RB" 'BEGIN{d=(b-a)/a*100; printf "RSS delta:  %+.2f%% (B vs A)\n", d}'
fi
