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
fail() {
    echo "FAIL: $*" >&2
    exit 1
}
ok() { echo "PASS: $*"; }
W="$DEDUP_WORK_DIR"
mkdir -p "$W" || fail "create work directory"
cp "$DEDUP_PHIX_FA" "$W/phix.fa" || fail "copy reference"
"$BWA_MEM3" index "$W/phix.fa" > /dev/null 2>&1 || fail "index"
"$HERE/../fixtures/make_dedup_reads.sh" "$W/phix.fa" "$W/dup" 600 8 \
    || fail "generate duplicate fixture" # dup-rich (4800 pairs)
"$HERE/../fixtures/make_dedup_reads.sh" "$W/phix.fa" "$W/uniq" 4800 1 \
    || fail "generate unique fixture" # low-dup  (4800 pairs)
run() {                               # run <tag> <mode-env> <threads> <prefix>
    local tag=$1 mode=$2 t=$3 p=$4
    BWAMEM3_DEDUP=$mode "$BWA_MEM3" mem -t "$t" "$W/phix.fa" \
        "${p}_1.fq" "${p}_2.fq" 2> "$W/$tag.err" | grep -v '^@PG' > "$W/$tag.sam" \
        || fail "$tag: mem run failed"
    [[ -s "$W/$tag.sam" ]] || fail "$tag: empty output"
}
for t in 1 4; do
    for p in dup uniq; do
        run "off_${p}_t${t}" off "$t" "$W/$p"
        run "on_${p}_t${t}" on "$t" "$W/$p"
        run "auto_${p}_t${t}" auto "$t" "$W/$p"
        cmp "$W/off_${p}_t${t}.sam" "$W/on_${p}_t${t}.sam" || fail "on!=off  ($p, -t $t)"
        cmp "$W/off_${p}_t${t}.sam" "$W/auto_${p}_t${t}.sam" || fail "auto!=off ($p, -t $t)"
        ok "byte-identity off==on==auto ($p, -t $t)"
    done
done
# line-count sanity: every input pair produced >=2 alignment records.
# grep -vc returns 1 for zero non-header lines (header-only SAM) and >1 on a real
# read error; under `set -e` the bare assignment would abort before fail() could
# report. Capture the status so header-only yields n=0 (which then trips the
# <9600 guard with a diagnostic) while a genuine grep read error still fails loud.
rc=0
n=$(grep -vc '^@' "$W/off_dup_t1.sam") || rc=$?
[[ "$rc" -gt 1 ]] && fail "off_dup_t1.sam: grep read error (status $rc)"
[[ "$n" -ge 9600 ]] || fail "record count $n < 9600"
ok "record-count sanity ($n records)"

# Guard: the dedup path must ACTUALLY be exercised. The fixture's 3 bp deletion
# forces gapped extension so reads reach banded SW; if that ever regresses to an
# exact-substring / ungapped-fast-path fixture, total_jobs would be 0 and the
# byte-identity checks above would be vacuously true. Assert it is non-zero.
# Capture the pipeline status: under `set -e` + `pipefail` a bare `tj=$(...)`
# assignment would abort the whole script on a nonzero `mem` exit, before the
# guard below can call fail().
tj=""
if ! tj=$(BWAMEM3_DEDUP=on BWAMEM3_DEDUP_STATS=1 "$BWA_MEM3" mem -t 1 "$W/phix.fa" \
    "$W/dup_1.fq" "$W/dup_2.fq" 2>&1 > /dev/null \
    | awk -F'total_jobs=' '/dedup-stats/{split($2,a," "); print a[1]}'); then
    fail "dedup-stats run failed"
fi
[[ "${tj:-0}" -gt 0 ]] || fail "dedup path not exercised (total_jobs=${tj:-0}); fixture is vacuous"
ok "dedup path exercised (${tj} banded-SW jobs deduped)"

# CLI flags: --dedup must mirror the env, flag wins over env, bad values fatal
"$BWA_MEM3" mem -t 1 --dedup on "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" 2> /dev/null \
    | grep -v '^@PG' > "$W/cli_on.sam" || fail "--dedup on run failed"
cmp "$W/off_dup_t1.sam" "$W/cli_on.sam" || fail "--dedup on != off (byte-identity)"
# Flag wins over env, observably: with the env asking for 'off' and the flag for
# 'on', the dedup path must still run (stats are emitted only from that path, so
# total_jobs>0 discriminates the winner -- a plain cmp is vacuous here because
# output is byte-identical in every mode).
prec_jobs=""
if ! prec_jobs=$(BWAMEM3_DEDUP=off BWAMEM3_DEDUP_STATS=1 "$BWA_MEM3" mem -t 1 --dedup on \
    "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" 2>&1 > /dev/null \
    | awk -F'total_jobs=' '/dedup-stats/{split($2,a," "); print a[1]}'); then
    fail "flag/env precedence run failed"
fi
[[ "${prec_jobs:-0}" -gt 0 ]] || fail "flag/env precedence: --dedup on did not win over BWAMEM3_DEDUP=off"
# The flag also short-circuits env MODE parsing entirely: an invalid env value must not be fatal.
BWAMEM3_DEDUP=bogus "$BWA_MEM3" mem -t 1 --dedup on "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" \
    2> /dev/null | grep -v '^@PG' > "$W/cli_prec.sam" || fail "flag-over-env run failed"
cmp "$W/off_dup_t1.sam" "$W/cli_prec.sam" || fail "flag/env precedence output mismatch"
if "$BWA_MEM3" mem --dedup bogus "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" > /dev/null 2> "$W/bogus.err"; then
    fail "--dedup bogus should be fatal"
fi
grep -q "expected off|on|auto" "$W/bogus.err" || fail "--dedup bogus: wrong error"
# An explicit-but-empty CLI value must be fatal, NOT silently inherit the env.
# With BWAMEM3_DEDUP=on set, `--dedup=` would otherwise fall through to the env
# in mem_dedup_configure() (which treats an empty mode_arg as absent); assert it
# still exits non-zero.
if BWAMEM3_DEDUP=on "$BWA_MEM3" mem --dedup= "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" \
    > /dev/null 2> "$W/empty.err"; then
    fail "--dedup= (empty) should be fatal"
fi
grep -q "expected off|on|auto" "$W/empty.err" || fail "--dedup= empty: wrong error"
# The env-only expert knobs must full-string parse, not take a leading prefix:
# atof("2x")==2.0 and atoll("12M")==12 both slip past a bare sign check and
# silently mis-configure the controller. A trailing-junk value must be fatal.
if BWAMEM3_DEDUP_Z=2x "$BWA_MEM3" mem "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" > /dev/null 2> "$W/z.err"; then
    fail "BWAMEM3_DEDUP_Z=2x should be fatal"
fi
grep -q "BWAMEM3_DEDUP_Z" "$W/z.err" || fail "BWAMEM3_DEDUP_Z=2x: wrong error"
if BWAMEM3_DEDUP_REPROBE=12M "$BWA_MEM3" mem "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" > /dev/null 2> "$W/rp.err"; then
    fail "BWAMEM3_DEDUP_REPROBE=12M should be fatal"
fi
grep -q "BWAMEM3_DEDUP_REPROBE" "$W/rp.err" || fail "BWAMEM3_DEDUP_REPROBE=12M: wrong error"
# Explicit-but-empty env values must be fatal, NOT silently inherit the default.
# `BWAMEM3_DEDUP=`/`BWAMEM3_DEDUP_Z=`/`BWAMEM3_DEDUP_REPROBE=` are set-but-empty:
# a bare `if (getenv(x) && *x)` guard treats them as absent and uses the default,
# hiding a cleared/misconfigured variable. Each empty form must exit non-zero.
if BWAMEM3_DEDUP='' "$BWA_MEM3" mem "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" > /dev/null 2> "$W/me.err"; then
    fail "BWAMEM3_DEDUP= (empty) should be fatal"
fi
grep -q "expected off|on|auto" "$W/me.err" || fail "BWAMEM3_DEDUP= empty: wrong error"
if BWAMEM3_DEDUP_Z='' "$BWA_MEM3" mem "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" > /dev/null 2> "$W/ze.err"; then
    fail "BWAMEM3_DEDUP_Z= (empty) should be fatal"
fi
grep -q "BWAMEM3_DEDUP_Z" "$W/ze.err" || fail "BWAMEM3_DEDUP_Z= empty: wrong error"
if BWAMEM3_DEDUP_REPROBE='' "$BWA_MEM3" mem "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" > /dev/null 2> "$W/re.err"; then
    fail "BWAMEM3_DEDUP_REPROBE= (empty) should be fatal"
fi
grep -q "BWAMEM3_DEDUP_REPROBE" "$W/re.err" || fail "BWAMEM3_DEDUP_REPROBE= empty: wrong error"
# The override exception: a non-empty --dedup must still win over an empty
# BWAMEM3_DEDUP= (the empty env is ignored, not fatal, when the flag is set).
BWAMEM3_DEDUP='' "$BWA_MEM3" mem -t 1 --dedup on "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" \
    2> /dev/null | grep -v '^@PG' > "$W/cli_over_empty.sam" || fail "--dedup on over empty env failed"
cmp "$W/off_dup_t1.sam" "$W/cli_over_empty.sam" || fail "--dedup on over empty env: output mismatch"
# `mem` with no args exits non-zero after printing usage() to stderr; under
# `set -o pipefail` a direct `... | grep -q` would fail the pipeline on that
# exit code even when grep matches, so capture first (matches the
# all_tiers_parity.sh / cohort_slice_identity.sh `|| true` house pattern).
USAGE_OUT="$("$BWA_MEM3" mem 2>&1 || true)"
grep -q -- '--dedup STR' <<< "$USAGE_OUT" || fail "usage() missing --dedup"
ok "CLI flags (--dedup / precedence / validation / usage)"
