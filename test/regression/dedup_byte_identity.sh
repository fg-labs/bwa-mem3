#!/usr/bin/env bash
# Regression: extension-DP job dedup must be byte-identical to no-dedup in
# EVERY mode. The dedup wrapper reuses a scored result only after a full
# byte-compare of (target, query, h0), so off/on/auto must produce identical
# SAM streams (modulo @PG CL, which embeds argv). Checked on a dup-rich and a
# low-dup fixture, single- and multi-threaded.
set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${DEDUP_PHIX_FA:?DEDUP_PHIX_FA must be set}"
: "${DEDUP_WORK_DIR:?DEDUP_WORK_DIR must be set}"
HERE="$(cd "$(dirname "$0")" && pwd)"
fail() { echo "FAIL: $*" >&2; exit 1; }
ok()   { echo "OK:   $*"; }
W="$DEDUP_WORK_DIR"; mkdir -p "$W"
cp "$DEDUP_PHIX_FA" "$W/phix.fa"
"$BWA_MEM3" index "$W/phix.fa" > /dev/null 2>&1 || fail "index"
"$HERE/../fixtures/make_dedup_reads.sh" "$W/phix.fa" "$W/dup"  600 8   # dup-rich (4800 pairs)
"$HERE/../fixtures/make_dedup_reads.sh" "$W/phix.fa" "$W/uniq" 4800 1  # low-dup  (4800 pairs)
run() { # run <tag> <mode-env> <threads> <prefix>
    local tag=$1 mode=$2 t=$3 p=$4
    BWAMEM3_DEDUP=$mode "$BWA_MEM3" mem -t "$t" "$W/phix.fa" \
        "${p}_1.fq" "${p}_2.fq" 2> "$W/$tag.err" | grep -v '^@PG' > "$W/$tag.sam" \
        || fail "$tag: mem run failed"
    [[ -s "$W/$tag.sam" ]] || fail "$tag: empty output"
}
for t in 1 4; do
  for p in dup uniq; do
    run "off_${p}_t${t}"  off  "$t" "$W/$p"
    run "on_${p}_t${t}"   on   "$t" "$W/$p"
    run "auto_${p}_t${t}" auto "$t" "$W/$p"
    cmp "$W/off_${p}_t${t}.sam" "$W/on_${p}_t${t}.sam"   || fail "on!=off  ($p, -t $t)"
    cmp "$W/off_${p}_t${t}.sam" "$W/auto_${p}_t${t}.sam" || fail "auto!=off ($p, -t $t)"
    ok "byte-identity off==on==auto ($p, -t $t)"
  done
done
# line-count sanity: every input pair produced >=2 alignment records
n=$(grep -vc '^@' "$W/off_dup_t1.sam"); [[ "$n" -ge 9600 ]] || fail "record count $n < 9600"
ok "record-count sanity ($n records)"
