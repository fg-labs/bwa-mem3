#!/usr/bin/env bash
# test/regression/aln_index_parity.sh
#
# Regression: bwa-mem3's index agrees row for row with a stock `bwa index`
# (0.7.19), so the `aln` port's FM-index adapter (src/aln/aln_bwt) returns the
# SA intervals and positions stock bwa would.
#
# For each reference -- the committed phix.fa and a generated FASTA (three
# contigs covering an N run, a reverse-complement copy and a tandem repeat, plus
# an .alt file naming one of them) -- the script indexes it with the oracle and
# with bwa-mem3 into separate directories, then runs aln_index_parity_test,
# which compares occ, SA, primary and L2 on every row.
#
# Inputs (env vars):
#   BWA_MEM3          -- path to the bwa-mem3 binary under test
#   BWA_ORACLE        -- path to a stock bwa 0.7.19 binary
#                        (test/tools/build_bwa_oracle.sh); unset => SKIP, unless
#   BWA_ORACLE_REQUIRED=1  -- (CI) an unset or non-executable BWA_ORACLE is a FAIL
#   FIXTURES          -- directory holding phix.fa (default: test/fixtures)
#   ALN_PARITY_BIN    -- path to aln_index_parity_test (default: repo root)

set -Eeuo pipefail
# shellcheck disable=SC2154  # rc is assigned by `rc=$?` at the start of the trap body
trap 'rc=$?; echo "FAIL: unexpected command failure at line $LINENO (exit $rc)" >&2; exit "$rc"' ERR

if [ "${BWA_ORACLE_REQUIRED:-0}" = 1 ] && [ ! -x "${BWA_ORACLE:-}" ]; then
    echo "FAIL: BWA_ORACLE_REQUIRED=1 but BWA_ORACLE ('${BWA_ORACLE:-}') is not an executable" >&2
    exit 1
fi
if [ -z "${BWA_ORACLE:-}" ]; then
    echo "SKIP: BWA_ORACLE unset"
    exit 0
fi
if [ -z "${BWA_MEM3:-}" ]; then
    echo "FAIL: BWA_MEM3 must be set" >&2
    exit 1
fi

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fixtures="${FIXTURES:-$repo_root/test/fixtures}"
parity_bin="${ALN_PARITY_BIN:-$repo_root/aln_index_parity_test}"

work="$(mktemp -d "${TMPDIR:-/tmp}/aln_index_parity.XXXXXX")"
cleanup() { rm -rf "$work"; }
trap cleanup EXIT

# Generated reference: the contigs the unit tests use.
"$parity_bin" --write-fasta "$work/synthetic.fa"
printf 'c2\n' > "$work/synthetic.fa.alt"

check() {
    local name="$1" fa="$2"
    mkdir -p "$work/$name/stock" "$work/$name/mem3"
    local d
    for d in stock mem3; do
        cp "$fa" "$work/$name/$d/ref.fa"
        if [ -f "$fa.alt" ]; then cp "$fa.alt" "$work/$name/$d/ref.fa.alt"; fi
    done
    "$BWA_ORACLE" index "$work/$name/stock/ref.fa" > "$work/$name/stock.log" 2>&1
    "$BWA_MEM3" index "$work/$name/mem3/ref.fa" > "$work/$name/mem3.log" 2>&1
    echo "aln_index_parity: $name"
    "$parity_bin" "$work/$name/stock/ref.fa" "$work/$name/mem3/ref.fa"
}

check phix "$fixtures/phix.fa"
check synthetic "$work/synthetic.fa"

echo "PASS: aln_index_parity.sh (occ, SA, primary, L2 identical to stock bwa 0.7.19 on phix and a synthetic reference)"
