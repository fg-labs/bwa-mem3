#!/usr/bin/env bash
# Regression: cross-read (k,s) SA-interval dedup must be byte-identical to
# no-dedup in EVERY mode. A duplicate (k,s) interval resolves to exactly the
# same coordinates in the same order as its representative, so the memcpy that
# replicates the rep's slot-run is a by-construction identity: off/on/auto must
# produce identical SAM streams (modulo @PG CL, which embeds argv). Checked on a
# dup-rich and a low-dup fixture, single- and multi-threaded.
set -euo pipefail
# The tool under test is a runtime dependency, not harness config: an unset,
# missing, or non-executable BWA_MEM3 means "cannot run here", not a failure.
# Emit SKIP: and exit 0 so it is not misclassified (the raw `:?` expansion below
# would print a shell error, and a nonexistent path would later report FAIL:
# index). The fixture/work-dir vars ARE harness config, so they keep `:?`.
if [[ -z "${BWA_MEM3:-}" || ! -x "$BWA_MEM3" ]]; then
    echo "SKIP: BWA_MEM3 is unavailable"
    exit 0
fi
: "${KS_DEDUP_PHIX_FA:?KS_DEDUP_PHIX_FA must be set}"
: "${KS_DEDUP_WORK_DIR:?KS_DEDUP_WORK_DIR must be set}"
HERE="$(cd "$(dirname "$0")" && pwd)"
fail() {
    echo "FAIL: $*" >&2
    exit 1
}
ok() { echo "PASS: $*"; }
W="$KS_DEDUP_WORK_DIR"
mkdir -p "$W" || fail "create work directory"
cp "$KS_DEDUP_PHIX_FA" "$W/phix.fa" || fail "copy reference"
"$BWA_MEM3" index "$W/phix.fa" > /dev/null 2>&1 || fail "index"
# Reuse the shared dedup fixture generator: N distinct templates, each repeated
# COPIES times. Identical read-pairs carry identical SMEMs and therefore
# identical (k,s) intervals within a chunk -- exactly the cross-read duplicates
# the SA-interval dedup collapses.
"$HERE/../fixtures/make_dedup_reads.sh" "$W/phix.fa" "$W/dup" 600 8 \
    || fail "generate duplicate fixture" # dup-rich (4800 pairs)
"$HERE/../fixtures/make_dedup_reads.sh" "$W/phix.fa" "$W/uniq" 4800 1 \
    || fail "generate unique fixture" # low-dup  (4800 pairs)
run() {                               # run <tag> <mode-env> <threads> <prefix>
    local tag=$1 mode=$2 t=$3 p=$4
    BWA3_KS_DEDUP=$mode "$BWA_MEM3" mem -t "$t" "$W/phix.fa" \
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

# Guard: the dedup path must ACTUALLY collapse cross-read duplicates, else the
# byte-identity checks above are vacuously true. BWA3_KS_DEDUP_STATS emits, at
# exit, `[ks-dedup-stats] total_pos=T distinct_pos=D saved=...` from the dedup
# path only. Assert T>0 (path exercised) AND D<T (duplicates were found and
# skipped). Capture the pipeline status: under `set -e`+`pipefail` a bare
# assignment would abort on the nonzero `mem` exit before the guard can fail().
stats=""
if ! stats=$(BWA3_KS_DEDUP=on BWA3_KS_DEDUP_STATS=1 "$BWA_MEM3" mem -t 1 "$W/phix.fa" \
    "$W/dup_1.fq" "$W/dup_2.fq" 2>&1 > /dev/null \
    | grep 'ks-dedup-stats'); then
    fail "ks-dedup-stats run failed or emitted no stats line"
fi
tot=$(sed -n 's/.*total_pos=\([0-9]*\).*/\1/p' <<< "$stats")
dis=$(sed -n 's/.*distinct_pos=\([0-9]*\).*/\1/p' <<< "$stats")
[[ "${tot:-0}" -gt 0 ]] || fail "dedup path not exercised (total_pos=${tot:-0}); fixture is vacuous"
[[ "${dis:-0}" -lt "${tot:-0}" ]] \
    || fail "dedup collapsed nothing (distinct_pos=${dis:-0} == total_pos=${tot:-0}); fixture is vacuous"
ok "dedup path exercised and non-vacuous (distinct_pos=$dis < total_pos=$tot)"

# CLI flags: --ks-dedup must mirror the env, flag wins over env, bad values fatal.
"$BWA_MEM3" mem -t 1 --ks-dedup on "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" 2> /dev/null \
    | grep -v '^@PG' > "$W/cli_on.sam" || fail "--ks-dedup on run failed"
cmp "$W/off_dup_t1.sam" "$W/cli_on.sam" || fail "--ks-dedup on != off (byte-identity)"
# Flag wins over env, observably: with the env asking for 'off' and the flag for
# 'on', the dedup path must still run (stats emit only from that path, so
# total_pos>0 discriminates the winner -- a plain cmp is vacuous here because
# output is byte-identical in every mode).
prec_tot=""
if ! prec_tot=$(BWA3_KS_DEDUP=off BWA3_KS_DEDUP_STATS=1 "$BWA_MEM3" mem -t 1 --ks-dedup on \
    "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" 2>&1 > /dev/null \
    | sed -n 's/.*total_pos=\([0-9]*\).*/\1/p'); then
    fail "flag/env precedence run failed"
fi
[[ "${prec_tot:-0}" -gt 0 ]] || fail "flag/env precedence: --ks-dedup on did not win over BWA3_KS_DEDUP=off"
# The flag also short-circuits env MODE parsing entirely: an invalid env value must not be fatal.
BWA3_KS_DEDUP=bogus "$BWA_MEM3" mem -t 1 --ks-dedup on "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" \
    2> /dev/null | grep -v '^@PG' > "$W/cli_prec.sam" || fail "flag-over-env run failed"
cmp "$W/off_dup_t1.sam" "$W/cli_prec.sam" || fail "flag/env precedence output mismatch"
if "$BWA_MEM3" mem --ks-dedup bogus "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" > /dev/null 2> "$W/bogus.err"; then
    fail "--ks-dedup bogus should be fatal"
fi
grep -q "expected off|on|auto" "$W/bogus.err" || fail "--ks-dedup bogus: wrong error"
# An explicit-but-empty CLI value must be fatal, NOT silently inherit the env.
if BWA3_KS_DEDUP=on "$BWA_MEM3" mem --ks-dedup= "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" \
    > /dev/null 2> "$W/empty.err"; then
    fail "--ks-dedup= (empty) should be fatal"
fi
grep -q "expected off|on|auto" "$W/empty.err" || fail "--ks-dedup= empty: wrong error"
# The env-only expert knobs must full-string parse, not take a leading prefix:
# atof("2x")==2.0 and atoll("12M")==12 both slip past a bare sign check and
# silently mis-configure the controller. A trailing-junk value must be fatal.
if BWA3_KS_DEDUP_Z=2x "$BWA_MEM3" mem "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" > /dev/null 2> "$W/z.err"; then
    fail "BWA3_KS_DEDUP_Z=2x should be fatal"
fi
grep -q "BWA3_KS_DEDUP_Z" "$W/z.err" || fail "BWA3_KS_DEDUP_Z=2x: wrong error"
if BWA3_KS_DEDUP_REPROBE=12M "$BWA_MEM3" mem "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" > /dev/null 2> "$W/rp.err"; then
    fail "BWA3_KS_DEDUP_REPROBE=12M should be fatal"
fi
grep -q "BWA3_KS_DEDUP_REPROBE" "$W/rp.err" || fail "BWA3_KS_DEDUP_REPROBE=12M: wrong error"
# Explicit-but-empty env values must be fatal, NOT silently inherit the default.
if BWA3_KS_DEDUP='' "$BWA_MEM3" mem "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" > /dev/null 2> "$W/me.err"; then
    fail "BWA3_KS_DEDUP= (empty) should be fatal"
fi
grep -q "expected off|on|auto" "$W/me.err" || fail "BWA3_KS_DEDUP= empty: wrong error"
if BWA3_KS_DEDUP_Z='' "$BWA_MEM3" mem "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" > /dev/null 2> "$W/ze.err"; then
    fail "BWA3_KS_DEDUP_Z= (empty) should be fatal"
fi
grep -q "BWA3_KS_DEDUP_Z" "$W/ze.err" || fail "BWA3_KS_DEDUP_Z= empty: wrong error"
if BWA3_KS_DEDUP_REPROBE='' "$BWA_MEM3" mem "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" > /dev/null 2> "$W/re.err"; then
    fail "BWA3_KS_DEDUP_REPROBE= (empty) should be fatal"
fi
grep -q "BWA3_KS_DEDUP_REPROBE" "$W/re.err" || fail "BWA3_KS_DEDUP_REPROBE= empty: wrong error"
# The override exception: a non-empty --ks-dedup must still win over an empty
# BWA3_KS_DEDUP= (the empty env is ignored, not fatal, when the flag is set).
BWA3_KS_DEDUP='' "$BWA_MEM3" mem -t 1 --ks-dedup on "$W/phix.fa" "$W/dup_1.fq" "$W/dup_2.fq" \
    2> /dev/null | grep -v '^@PG' > "$W/cli_over_empty.sam" || fail "--ks-dedup on over empty env failed"
cmp "$W/off_dup_t1.sam" "$W/cli_over_empty.sam" || fail "--ks-dedup on over empty env: output mismatch"
# `mem` with no args exits non-zero after printing usage() to stderr; under
# `set -o pipefail` a direct `... | grep -q` would fail the pipeline on that
# exit code even when grep matches, so capture first (matches the
# all_tiers_parity.sh / cohort_slice_identity.sh `|| true` house pattern).
USAGE_OUT="$("$BWA_MEM3" mem 2>&1 || true)"
grep -q -- '--ks-dedup STR' <<< "$USAGE_OUT" || fail "usage() missing --ks-dedup"
ok "CLI flags (--ks-dedup / precedence / validation / usage)"
