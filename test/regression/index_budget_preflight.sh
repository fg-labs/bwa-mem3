#!/usr/bin/env bash
# test/regression/index_budget_preflight.sh
#
# A plain (non --meth) `index` must check its memory budget before it writes
# any output file. Rebuilding over a prefix that already holds an index is the
# case that matters: if the budget check runs only once the FM-index build
# starts, a rejected rebuild has already replaced <prefix>.pac/.ann/.amb with the
# new reference's, leaving them beside the old <prefix>.bwt.2bit.64. When the
# reference changed, that prefix is a mix of two indexes and aligns silently
# wrong.
#
# Cases:
#   1. a budget-exceeding rebuild over an existing index is rejected with exit
#      code 3 and the standard shortfall message, every existing index file is
#      byte-identical afterwards, and no scratch file is left beside them;
#   2. an over-budget build into an empty directory is rejected the same way
#      and writes nothing there;
#   3. an empty reference over an existing index is rejected without touching
#      it (the FM-index build used to reject it only after the replacement);
#   4. a rebuild within budget still replaces the index, byte-identical to a
#      fresh build of the same reference at a new prefix.
#
# The references are generated here (fixed-seed awk), not committed.
#
# Environment: BWA_MEM3 is required (the binary under test: a path, or a
# command on PATH). No staged inputs.

set -euo pipefail

[[ -n ${BWA_MEM3:-} ]] || {
    echo "FAIL: BWA_MEM3 must be set to the bwa-mem3 binary under test" >&2
    exit 1
}
bin="$BWA_MEM3"
[[ $bin == */* ]] || bin="$(command -v "$bin")" || {
    echo "FAIL: BWA_MEM3=$BWA_MEM3 is not on PATH" >&2
    exit 1
}
[[ -f "$bin" && -x "$bin" ]] || {
    echo "FAIL: bwa-mem3 binary not executable at $bin" >&2
    exit 1
}

fail() { # $1 = stderr file to dump, $2... = message
    local err="$1"
    shift
    echo "FAIL [index budget preflight]: $*" >&2
    if [[ -s "$err" ]]; then
        echo "--- index stderr:" >&2
        cat "$err" >&2
    fi
    exit 1
}

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

# Emit a FASTA of `n_contigs` random contigs of `len` bases each, one sequence
# line per contig, from a fixed seed so every run sees the same reference.
gen_ref() { # $1 = seed, $2 = n_contigs, $3 = len, $4 = output path
    awk -v seed="$1" -v n="$2" -v len="$3" 'BEGIN {
        srand(seed)
        split("A C G T", b, " ")
        for (c = 1; c <= n; c++) {
            printf(">ctg%d\n", c)
            s = ""
            for (i = 0; i < len; i++) s = s b[int(rand() * 4) + 1]
            print s
        }
    }' > "$4"
}

# Two different references, so a mixed prefix is detectable by content and by
# length: B is longer than A.
gen_ref 11 4 25000 "$tmp/a.fa"
gen_ref 29 5 40000 "$tmp/b.fa"

files=(pac ann amb bwt.2bit.64)
prefix="$tmp/idx/ref"
mkdir -p "$tmp/idx" "$tmp/snap"

# A ~200 kbp reference needs ~3 MiB; 64M clears it with room to spare, and 1M
# cannot.
"$bin" index --max-memory 64M -p "$prefix" "$tmp/a.fa" > /dev/null 2> "$tmp/build_a.err" \
    || fail "$tmp/build_a.err" "baseline build of reference A failed"
for f in "${files[@]}"; do
    cp "$prefix.$f" "$tmp/snap/ref.$f"
done
before_listing="$(cd "$tmp/idx" && ls -A)"

# --- 1: rejected rebuild leaves the existing index untouched ---------------
rc=0
"$bin" index --max-memory 1M -p "$prefix" --force "$tmp/b.fa" > /dev/null 2> "$tmp/reject.err" || rc=$?
[[ $rc -eq 3 ]] || fail "$tmp/reject.err" "over-budget rebuild exited $rc, expected 3"
grep -q 'ERROR: libsais needs ~.* to index a .* reference; --max-memory is 1.0 MiB' "$tmp/reject.err" \
    || fail "$tmp/reject.err" "over-budget rebuild did not print the budget shortfall"
# --max-memory was given, so the remedy is raising it, not a bigger host.
grep -q 'Raise --max-memory to at least' "$tmp/reject.err" \
    || fail "$tmp/reject.err" "user-budget rejection did not suggest raising --max-memory"
if grep -q 'Retry on a host' "$tmp/reject.err"; then
    fail "$tmp/reject.err" "user-budget rejection printed the host-size hint"
fi
for f in "${files[@]}"; do
    cmp -s "$prefix.$f" "$tmp/snap/ref.$f" \
        || fail "$tmp/reject.err" "rejected rebuild modified $f: the prefix now mixes two indexes"
done
after_listing="$(cd "$tmp/idx" && ls -A)"
[[ "$after_listing" == "$before_listing" ]] \
    || fail "$tmp/reject.err" "rejected rebuild changed the output directory:
before: $before_listing
after:  $after_listing"
echo "PASS: over-budget rebuild rejected (exit 3), existing index byte-identical"

# --- 2: rejected build into an empty directory writes nothing --------------
mkdir "$tmp/empty"
rc=0
"$bin" index --max-memory 1M -p "$tmp/empty/ref" "$tmp/b.fa" > /dev/null 2> "$tmp/reject_empty.err" || rc=$?
[[ $rc -eq 3 ]] || fail "$tmp/reject_empty.err" "over-budget fresh build exited $rc, expected 3"
leftover="$(ls -A "$tmp/empty")"
[[ -z "$leftover" ]] || fail "$tmp/reject_empty.err" "over-budget fresh build left files behind: $leftover"
echo "PASS: over-budget fresh build rejected (exit 3), nothing written"

# --- 3: empty reference leaves the existing index untouched -----------------
: > "$tmp/none.fa"
rc=0
"$bin" index --max-memory 64M -p "$prefix" --force "$tmp/none.fa" > /dev/null 2> "$tmp/none.err" || rc=$?
[[ $rc -eq 1 ]] || fail "$tmp/none.err" "empty-reference rebuild exited $rc, expected 1"
grep -q "malformed reference" "$tmp/none.err" \
    || fail "$tmp/none.err" "empty-reference rebuild did not report the empty reference"
for f in "${files[@]}"; do
    cmp -s "$prefix.$f" "$tmp/snap/ref.$f" \
        || fail "$tmp/none.err" "empty-reference rebuild modified $f"
done
[[ "$(cd "$tmp/idx" && ls -A)" == "$before_listing" ]] \
    || fail "$tmp/none.err" "empty-reference rebuild changed the output directory"
echo "PASS: empty-reference rebuild rejected (exit 1), existing index byte-identical"

# --- 4: within-budget rebuild still replaces the index ----------------------
"$bin" index --max-memory 64M -p "$prefix" --force "$tmp/b.fa" > /dev/null 2> "$tmp/rebuild.err" \
    || fail "$tmp/rebuild.err" "within-budget rebuild of reference B failed"
"$bin" index --max-memory 64M -p "$tmp/fresh" "$tmp/b.fa" > /dev/null 2> "$tmp/fresh.err" \
    || fail "$tmp/fresh.err" "fresh build of reference B failed"
for f in "${files[@]}"; do
    cmp -s "$prefix.$f" "$tmp/fresh.$f" \
        || fail "$tmp/rebuild.err" "rebuilt $f differs from a fresh build of the same reference"
done
cmp -s "$prefix.pac" "$tmp/snap/ref.pac" \
    && fail "$tmp/rebuild.err" "within-budget rebuild left the old .pac in place"
after_rebuild="$(cd "$tmp/idx" && ls -A)"
[[ "$after_rebuild" == "$before_listing" ]] \
    || fail "$tmp/rebuild.err" "within-budget rebuild left extra files:
expected: $before_listing
got:      $after_rebuild"
echo "PASS: within-budget rebuild matches a fresh build"

echo "PASS: index_budget_preflight"
