#!/usr/bin/env bash
# Local unit-test harness for the five C++ binaries under test/.
# Builds, indexes the phiX fixture, runs each binary, asserts exit 0
# and non-empty output. No hash pinning — that lives in the CI workflow.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
BWAMEM2="$ROOT/bwa-mem2"
FIXTURES="$HERE/fixtures"

fail() { echo "FAIL: $*" >&2; exit 1; }
ok()   { echo "OK:   $*"; }

[[ -x "$BWAMEM2" ]] || fail "bwa-mem2 not built at $BWAMEM2 (run 'make' or 'make arm64' first)"

# Build the five unit binaries.
(cd "$HERE" && make) || fail "test/ make failed"

# Build the phiX FMI index if not already present. Check all five artifacts
# that `bwa-mem2 index` produces so a corrupt/partial prior run is re-indexed.
if [[ ! -s "$FIXTURES/phix.fa.bwt.2bit.64" || \
      ! -s "$FIXTURES/phix.fa.0123"        || \
      ! -s "$FIXTURES/phix.fa.amb"         || \
      ! -s "$FIXTURES/phix.fa.ann"         || \
      ! -s "$FIXTURES/phix.fa.pac" ]]; then
    "$BWAMEM2" index "$FIXTURES/phix.fa" >/dev/null 2>&1 || fail "bwa-mem2 index on phix.fa failed"
fi

# --- fmi_test --------------------------------------------------------------
OUT="$(cd "$HERE" && ./fmi_test "$FIXTURES/phix.fa" "$FIXTURES/reads.fa" 10 19 1 2>&1)"
echo "$OUT" | grep -q 'numReads = 10'       || fail "fmi_test: numReads line missing"
echo "$OUT" | grep -q 'totalSmems ='        || fail "fmi_test: totalSmems line missing"
echo "$OUT" | grep -qE '\[[0-9]+,[0-9]+\]'  || fail "fmi_test: no SMEM output (PRINT_OUTPUT patch not applied?)"
ok "fmi_test"

# --- smem2_test ------------------------------------------------------------
OUT="$(cd "$HERE" && ./smem2_test "$FIXTURES/phix.fa" "$FIXTURES/smem2_input.txt" 5 50 19 1 2>&1)"
echo "$OUT" | grep -q 'numReads = 5'        || fail "smem2_test: numReads line missing"
echo "$OUT" | grep -q 'totalSmems ='        || fail "smem2_test: totalSmems line missing"
ok "smem2_test"

# --- bwt_seed_strategy_test ------------------------------------------------
# Note: unlike fmi_test and smem2_test, this binary doesn't emit a
# "numReads = N" line. Only assert the post-SW totalSmems marker.
OUT="$(cd "$HERE" && ./bwt_seed_strategy_test "$FIXTURES/phix.fa" "$FIXTURES/bwt_seed_input.fa" 5 50 19 20 2>&1)"
echo "$OUT" | grep -q 'minSeedLen ='        || fail "bwt_seed_strategy_test: minSeedLen line missing (read-parsing broken?)"
echo "$OUT" | grep -q 'totalSmems ='        || fail "bwt_seed_strategy_test: totalSmems line missing"
ok "bwt_seed_strategy_test"

# --- sa2ref_test -----------------------------------------------------------
OUT_FILE="$(mktemp)"
FKSW="$HERE/fksw.txt"
trap 'rm -f "$OUT_FILE" "$FKSW"' EXIT
(cd "$HERE" && ./sa2ref_test "$FIXTURES/phix.fa" "$FIXTURES/sa2ref_input.txt" "$OUT_FILE" 20 >/dev/null 2>&1) || fail "sa2ref_test crashed"
[[ -s "$OUT_FILE" ]] || fail "sa2ref_test: output file empty"
ok "sa2ref_test (wrote $(wc -l < "$OUT_FILE" | tr -d ' ') coords)"

# --- xeonbsw (main_banded) -------------------------------------------------
(cd "$HERE" && ./xeonbsw -pairs "$FIXTURES/pairs.txt" >/dev/null 2>&1) || fail "xeonbsw crashed"
[[ -s "$FKSW" ]] || fail "xeonbsw: fksw.txt empty (main_banded fksw-write patch not applied?)"
LINES="$(wc -l < "$FKSW" | tr -d ' ')"
[[ "$LINES" -eq 3 ]] || fail "xeonbsw: expected 3 lines in fksw.txt, got $LINES"
ok "xeonbsw ($LINES pair scores emitted)"

# --- pg_cl_escape_test ----------------------------------------------------
# Regression for issue #45 / upstream #293: tabs inside `-R` must not
# bleed into the @PG CL: value. Uses the same phiX fixture as the rest
# of the harness (indexed above if needed).
"$HERE/pg_cl_escape_test.sh" "$BWAMEM2" "$FIXTURES" || fail "pg_cl_escape_test failed"
ok "pg_cl_escape_test"

echo "ALL UNIT TESTS PASSED"
