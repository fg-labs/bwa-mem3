#!/usr/bin/env bash
# Regression: whole-read-pair memoization (--dedup-reads / BWAMEM3_DEDUP_READS)
# must be byte-identical to no-memoization in EVERY mode. A duplicate PAIR skips
# seed->chain->extend and copies its representative's regs, then replays the
# per-read SAM stage (pairing / MAPQ / tie-break / QNAME / QUAL), so off/on/auto
# must produce identical SAM (modulo @PG CL, which embeds argv).
#
# Why this needs its own gate: --dedup-reads is on-by-default (auto), but the
# net-cycles controller stays OFF on cheap alignments (phix included), so the
# other byte-identity gates run the flag at its default and never engage the
# memoize path -- a PASS there says nothing about it. This test forces the path
# ON and asserts non-vacuity via BWAMEM3_DEDUP_READS_STATS. Checked on a dup-rich
# and a low-dup fixture, on an N+lowercase fixture (the DUP-read 2-bit conversion
# path, which kernel1 would otherwise skip), single- and multi-threaded.
set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${DEDUP_READS_PHIX_FA:?DEDUP_READS_PHIX_FA must be set}"
: "${DEDUP_READS_WORK_DIR:?DEDUP_READS_WORK_DIR must be set}"
HERE="$(cd "$(dirname "$0")" && pwd)"
fail() {
    echo "FAIL: $*" >&2
    exit 1
}
ok() { echo "PASS: $*"; }
W="$DEDUP_READS_WORK_DIR"
mkdir -p "$W" || fail "create work directory"
cp "$DEDUP_READS_PHIX_FA" "$W/phix.fa" || fail "copy reference"
"$BWA_MEM3" index "$W/phix.fa" > /dev/null 2>&1 || fail "index"

# dup-rich: each 97 bp concordant pair repeated 8x under distinct read names, so
# every window is one representative + 7 byte-identical duplicate PAIRS.
"$HERE/../fixtures/make_dedup_reads.sh" "$W/phix.fa" "$W/dup" 600 8 \
    || fail "generate duplicate fixture" # 4800 pairs
# low-dup control: every pair unique (no memoization opportunity).
"$HERE/../fixtures/make_dedup_reads.sh" "$W/phix.fa" "$W/uniq" 4800 1 \
    || fail "generate unique fixture" # 4800 pairs

# D1 fixture: N + lowercase bases in duplicate pairs. DUP reads skip kernel1's
# in-place ASCII->2-bit conversion, so worker_bwt_memo must convert them up
# front or the SAM stage sees garbage bases. Derive from the dup fixture with a
# uniform, deterministic per-base mangle so a window's repeats stay identical
# (and therefore still detected as duplicate pairs). Length is preserved (97).
mangle() { # <in.fq> <out.fq> : lowercase bases 10-15, set bases 50 & 52 to N
    awk 'NR % 4 == 2 {
        s = substr($0, 1, 9) tolower(substr($0, 10, 6)) substr($0, 16)
        s = substr(s, 1, 49) "N" substr(s, 51, 1) "N" substr(s, 53)
        print s
        next
    } { print }' "$1" > "$2"
}
mangle "$W/dup_1.fq" "$W/nlow_1.fq" || fail "generate N/lowercase R1"
mangle "$W/dup_2.fq" "$W/nlow_2.fq" || fail "generate N/lowercase R2"

run() { # run <tag> <mode> <threads> <prefix>
    local tag=$1 mode=$2 t=$3 p=$4
    "$BWA_MEM3" mem --dedup-reads "$mode" -t "$t" "$W/phix.fa" \
        "${p}_1.fq" "${p}_2.fq" 2> "$W/$tag.err" | grep -v '^@PG' > "$W/$tag.sam" \
        || fail "$tag: mem run failed"
    [[ -s "$W/$tag.sam" ]] || fail "$tag: empty output"
}

for t in 1 4; do
    for p in dup uniq nlow; do
        run "off_${p}_t${t}" off "$t" "$W/$p"
        run "on_${p}_t${t}" on "$t" "$W/$p"
        run "auto_${p}_t${t}" auto "$t" "$W/$p"
        cmp "$W/off_${p}_t${t}.sam" "$W/on_${p}_t${t}.sam" || fail "on!=off  ($p, -t $t)"
        cmp "$W/off_${p}_t${t}.sam" "$W/auto_${p}_t${t}.sam" || fail "auto!=off ($p, -t $t)"
        ok "byte-identity off==on==auto ($p, -t $t)"
    done
done

# Thread determinism of the ON path specifically: the compaction and copy pass
# run under kt_for, so -t1 must equal -t4 with the memoize path live.
for p in dup nlow; do
    cmp "$W/on_${p}_t1.sam" "$W/on_${p}_t4.sam" || fail "on -t1 != -t4 ($p)"
    ok "thread determinism --dedup-reads on ($p, -t1 == -t4)"
done

# Non-vacuity: the memoize path must ACTUALLY be exercised, or on==off is a
# vacuous pass. In on-mode, BWAMEM3_DEDUP_READS_STATS must report dup>0 and that
# the memo engaged (final=ON). Capture the pipeline status: under `set -e` +
# `pipefail` a bare `stats=$(...)` would abort on a nonzero `mem` exit before
# fail() could report.
stats=""
if ! stats=$(BWAMEM3_DEDUP_READS_STATS=1 "$BWA_MEM3" mem --dedup-reads on -t 1 \
    "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" 2>&1 > /dev/null \
    | grep 'dedup-reads-stats'); then
    fail "dedup-reads-stats run failed"
fi
dup=$(awk -F'dup=' '{split($2, a, " "); print a[1]}' <<< "$stats")
fin=$(awk -F'final=' '{split($2, a, " "); print a[1]}' <<< "$stats")
[[ "${dup:-0}" -gt 0 ]] || fail "memoize path not exercised (dup=${dup:-0}); fixture is vacuous"
[[ "$fin" == "ON" ]] || fail "on-mode did not engage the memoize path (final=$fin)"
ok "memoize path exercised (${dup} duplicate pairs, final=ON)"

# CLI validation: a bad value and an explicit-but-empty value are both fatal
# (the empty form must NOT silently inherit the env), and usage advertises the
# flag. Mirrors dedup_byte_identity.sh's --dedup checks.
if "$BWA_MEM3" mem --dedup-reads bogus "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" \
    > /dev/null 2> "$W/bogus.err"; then
    fail "--dedup-reads bogus should be fatal"
fi
grep -q "expected off|on|auto" "$W/bogus.err" || fail "--dedup-reads bogus: wrong error"
if BWAMEM3_DEDUP_READS=on "$BWA_MEM3" mem --dedup-reads= "$W/phix.fa" \
    "$W/dup_1.fq" "$W/dup_2.fq" > /dev/null 2> "$W/empty.err"; then
    fail "--dedup-reads= (empty) should be fatal"
fi
grep -q "expected off|on|auto" "$W/empty.err" || fail "--dedup-reads= empty: wrong error"
# `mem` with no args exits non-zero after printing usage() to stderr; capture
# first so `set -o pipefail` does not fail the pipeline on that exit code (the
# all_tiers_parity.sh / dedup_byte_identity.sh house pattern).
USAGE_OUT="$("$BWA_MEM3" mem 2>&1 || true)"
grep -q -- '--dedup-reads STR' <<< "$USAGE_OUT" || fail "usage() missing --dedup-reads"
ok "CLI flags (--dedup-reads / validation / usage)"
