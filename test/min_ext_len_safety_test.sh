#!/usr/bin/env bash
# test/min_ext_len_safety_test.sh
#
# Guards the --min-ext-len recall-safety invariant (the "non-emptying" filter,
# mem_chain_drop_short_seeds): a chain with no seed >= min_ext_len is left
# untouched, so the filter can never empty a chain and never lose a read. A
# min_ext_len larger than any possible seed therefore makes EVERY chain all-short
# -> every chain untouched -> output identical to default.
#
# Regression for the smoke-1M collapse: the previous (emptying) filter dropped
# every seed < min_ext_len unconditionally, so reads whose only seeds were short
# (low-mappability / repetitive) went unmapped -- 63% of a real 151bp WGS sample.
# Under that behavior a huge --min-ext-len would unmap (nearly) everything; under
# the fixed filter it is a no-op.
#
# Usage: test/min_ext_len_safety_test.sh <bwa-mem3-binary> <fixtures-dir>
set -euo pipefail

[[ $# -eq 2 ]] || {
    echo "usage: $0 <bwa-mem3-binary> <fixtures-dir>" >&2
    exit 2
}
bin="$1"
fixtures="$2"
ref="$fixtures/phix.fa"
reads="$fixtures/reads.fa"

[[ -x "$bin" ]] || {
    echo "FAIL: binary not executable: $bin" >&2
    exit 1
}
[[ -s "$ref" ]] || {
    echo "FAIL: phix.fa missing: $ref" >&2
    exit 1
}
[[ -s "$reads" ]] || {
    echo "FAIL: reads.fa missing: $reads" >&2
    exit 1
}

if [[ ! -s "$ref.bwt.2bit.64" || ! -s "$ref.amb" || ! -s "$ref.ann" || ! -s "$ref.pac" ]]; then
    "$bin" index "$ref" > /dev/null 2>&1 || {
        echo "FAIL: index phix.fa" >&2
        exit 1
    }
fi

d="$(mktemp -d)"
trap 'rm -rf "$d"' EXIT

# Alignment records only (strip @ header lines, whose @PG CL: differs by flags).
# -t 1 for deterministic output. A min-ext-len far above any seed length makes
# every chain all-short, which the fixed filter must leave untouched.
#
# The aligner and the header strip run as separate steps so their exit statuses
# stay distinguishable. Piped together, a run that emitted headers but no
# alignments makes `grep -v` exit 1 (it selected nothing), and under
# `set -o pipefail` that ends the script before the explicit "produced no
# alignment records" check below -- the one failure this test most needs to
# report. Split apart, a genuine aligner failure still aborts via `set -e`
# while an empty record set reaches its diagnostic.
"$bin" mem -t 1 "$ref" "$reads" 2> /dev/null > "$d/default.raw"
"$bin" mem -t 1 --min-ext-len 1000000 "$ref" "$reads" 2> /dev/null > "$d/huge.raw"
grep -v '^@' "$d/default.raw" > "$d/default.sam" || true
grep -v '^@' "$d/huge.raw" > "$d/huge.sam" || true

[[ -s "$d/default.sam" ]] || {
    echo "FAIL: default run produced no alignment records" >&2
    exit 1
}

# 1. Core invariant: no read is lost. Count mapped records (FLAG bit 0x4 clear);
#    int($2/4)%2==0 is pure arithmetic so it works in any awk (no and()).
dmap=$(awk 'int($2/4)%2==0' "$d/default.sam" | wc -l | tr -d ' ')
hmap=$(awk 'int($2/4)%2==0' "$d/huge.sam" | wc -l | tr -d ' ')
[[ "$dmap" -gt 0 ]] || {
    echo "FAIL: default mapped 0 records -- fixture cannot exercise the test" >&2
    exit 1
}
[[ "$dmap" == "$hmap" ]] \
    || {
        echo "FAIL: huge --min-ext-len dropped reads ($dmap mapped at default vs $hmap)" >&2
        exit 1
    }
echo "OK:   huge --min-ext-len keeps all $dmap mapped records (filter never empties a chain)"

# 2. Stronger: the whole alignment block is byte-identical to default, since every
#    chain is all-short and therefore untouched.
if diff -q "$d/default.sam" "$d/huge.sam" > /dev/null; then
    echo "OK:   huge --min-ext-len output byte-identical to default (all-short chains untouched)"
else
    echo "FAIL: huge --min-ext-len diverged from default output" >&2
    # `|| true`: this runs in the "files differ" branch, so diff exits 1 by
    # construction and pipefail would abort before the explicit exit below.
    diff "$d/default.sam" "$d/huge.sam" | head -10 >&2 || true
    exit 1
fi

echo "PASS: min_ext_len_safety_test"
