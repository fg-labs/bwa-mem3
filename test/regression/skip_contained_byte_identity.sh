#!/usr/bin/env bash
# test/regression/skip_contained_byte_identity.sh
#
# Regression: contained-seed extension skipping is ON by default, and it must be
# byte-identical to the reference extension path (`--keep-contained-ext`).
#
# The default two-wave skip defers a seed that is contained (same diagonal) in a
# longer in-chain seed past the main extension batch, then skips its banded-SW
# only when the post-extension containment purge predicate (PE18) confirms it
# against the container's REAL extended alnreg; a seed it does not confirm is
# extended in a second batch. The skip set is therefore a proven subset of the
# post-extension purge and the predicate is scoring-independent, so the default
# output is bit-for-bit the reference path's. This pins that contract in CI —
# the whole justification for making it the default rests on it.
#
# Also asserts the two contract properties the CLI change introduced:
#   - `--keep-contained-ext` opts out to the reference path (== default output).
#   - `--skip-contained-ext` is a deprecated accepted no-op (== default output),
#     and prints a deprecation notice to stderr.
#
# Ordering is contractually deterministic across thread counts, so every
# comparison runs at -t 1 and -t 4. Checked on a dup-rich and a low-dup fixture.
# (--meth parity is validated separately on production bisulfite data; a meth CI
# A/B is tracked as a follow-up.)
#
# Fixtures are generated here, not committed (repo convention).
#
# Inputs (env vars):
#   BWA_MEM3                    — path to the bwa-mem3 binary under test
#   SKIP_CONTAINED_PHIX_FA      — phix reference FASTA
#   SKIP_CONTAINED_WORK_DIR     — private scratch directory for this script

set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${SKIP_CONTAINED_PHIX_FA:?SKIP_CONTAINED_PHIX_FA must be set}"
: "${SKIP_CONTAINED_WORK_DIR:?SKIP_CONTAINED_WORK_DIR must be set}"

HERE="$(cd "$(dirname "$0")" && pwd)"
fail() {
    echo "FAIL: $*" >&2
    exit 1
}
ok() { echo "PASS: $*"; }

W="$SKIP_CONTAINED_WORK_DIR"
mkdir -p "$W" || fail "create work directory"
cp "$SKIP_CONTAINED_PHIX_FA" "$W/phix.fa" || fail "copy reference"
"$BWA_MEM3" index "$W/phix.fa" > /dev/null 2>&1 || fail "index"

"$HERE/../fixtures/make_dedup_reads.sh" "$W/phix.fa" "$W/dup" 600 8 \
    || fail "generate duplicate fixture" # dup-rich (4800 pairs)
"$HERE/../fixtures/make_dedup_reads.sh" "$W/phix.fa" "$W/uniq" 4800 1 \
    || fail "generate unique fixture" # low-dup  (4800 pairs)

run() { # run <tag> <threads> <prefix> [extra flags...]
    local tag=$1 t=$2 p=$3
    shift 3
    "$BWA_MEM3" mem -t "$t" "$@" "$W/phix.fa" "${p}_1.fq" "${p}_2.fq" \
        2> "$W/$tag.err" | grep -v '^@PG' > "$W/$tag.sam" \
        || fail "$tag: mem run failed"
    [[ -s "$W/$tag.sam" ]] || fail "$tag: empty output"
}

for t in 1 4; do
    for p in dup uniq; do
        run "def_${p}_t${t}" "$t" "$W/$p"                       # default: skip ON
        run "keep_${p}_t${t}" "$t" "$W/$p" --keep-contained-ext # reference path
        run "dep_${p}_t${t}" "$t" "$W/$p" --skip-contained-ext  # deprecated no-op

        cmp "$W/def_${p}_t${t}.sam" "$W/keep_${p}_t${t}.sam" \
            || fail "default(skip) != --keep-contained-ext ($p, -t $t)"
        cmp "$W/def_${p}_t${t}.sam" "$W/dep_${p}_t${t}.sam" \
            || fail "--skip-contained-ext (deprecated) != default ($p, -t $t)"
        ok "byte-identity default == --keep-contained-ext == --skip-contained-ext ($p, -t $t)"
    done
done

# The deprecated flag must announce itself (so users learn to migrate).
grep -q 'skip-contained-ext is deprecated' "$W/dep_dup_t1.err" \
    || fail "--skip-contained-ext did not print a deprecation notice"
ok "--skip-contained-ext prints a deprecation notice"

# Record-count sanity: every input pair produced >=2 alignment records.
# grep -vc returns 1 for a header-only SAM and >1 on a real read error; capture
# the status so header-only trips the count guard with a diagnostic while a
# genuine grep error still fails loud.
rc=0
n=$(grep -vc '^@' "$W/def_dup_t1.sam") || rc=$?
[[ "$rc" -gt 1 ]] && fail "def_dup_t1.sam: grep read error (status $rc)"
[[ "$n" -ge 9600 ]] || fail "record count $n < 9600"
ok "record count sanity ($n records)"
