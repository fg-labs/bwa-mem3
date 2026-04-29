#!/usr/bin/env bash
# test/help_prescan_test.sh
#
# `bwa-mem3 mem` pre-scans argv for `--help` so the AVX/SA banner and the
# post-run profiling trailer are skipped on a real help request. The scan
# must NOT match `--help` when it appears as the value of an option that
# takes an argument (e.g. `-R --help`, `-o --help`,
# `--set-as-failed --help`). Without the fix, those legitimate runs
# wrongly short-circuit to `main_mem` and suppress the banner.
#
# Usage: test/help_prescan_test.sh <bwa-mem3-binary> <fixtures-dir>

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 <bwa-mem3-binary> <fixtures-dir>" >&2
    exit 2
fi

# Resolve to absolute paths so the `cd "$tmpdir"` block in test 2 still
# finds the binary, reference, and reads regardless of how the script
# was invoked.
abspath() { (cd "$(dirname "$1")" && printf '%s/%s\n' "$(pwd)" "$(basename "$1")"); }
bin="$(abspath "$1")"
fixtures="$(abspath "$2")"
ref="$fixtures/phix.fa"
reads="$fixtures/reads.fa"

[[ -x "$bin" ]]   || { echo "FAIL: bwa-mem3 binary not executable at $bin" >&2; exit 1; }
[[ -s "$ref" ]]   || { echo "FAIL: phix.fa missing at $ref" >&2; exit 1; }
[[ -s "$reads" ]] || { echo "FAIL: reads.fa missing at $reads" >&2; exit 1; }

# Build the phiX FMI index if not already present.
if [[ ! -s "$ref.bwt.2bit.64" || ! -s "$ref.0123" || ! -s "$ref.amb" \
      || ! -s "$ref.ann"       || ! -s "$ref.pac" ]]; then
    "$bin" index "$ref" >/dev/null 2>&1 \
        || { echo "FAIL: bwa-mem3 index on phix.fa failed" >&2; exit 1; }
fi

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

# --- 1. Real `mem --help` short-circuits: option list, no banner, exit 0. ---
set +e
"$bin" mem --help > "$tmpdir/help.out" 2> "$tmpdir/help.err"
rc=$?
set -e
[[ $rc -eq 0 ]] \
    || { echo "FAIL: mem --help exited with $rc, expected 0" >&2
         echo "---- stderr ----" >&2; cat "$tmpdir/help.err" >&2; exit 1; }
grep -q 'Algorithm options' "$tmpdir/help.err" \
    || { echo "FAIL: mem --help did not print the usage option list" >&2
         echo "---- stderr ----" >&2; cat "$tmpdir/help.err" >&2; exit 1; }
if grep -q 'Executing in' "$tmpdir/help.err"; then
    echo "FAIL: mem --help printed the AVX/SA banner; pre-scan should suppress it" >&2
    exit 1
fi

# --- 2. `mem -o --help <ref> <reads>` is a legitimate run whose -o value -----
# is the literal string `--help`. The banner must be printed, the run must
# exit 0, and SAM must be written to a file literally named `--help`. Run
# from the tmp dir so the output file lands there even though `--help`
# looks like an option.
set +e
( cd "$tmpdir" && "$bin" mem -o "--help" "$ref" "$reads" > run.out 2> run.err )
rc=$?
set -e
[[ $rc -eq 0 ]] \
    || { echo "FAIL: 'mem -o --help' exited $rc, expected 0" >&2
         echo "---- stderr ----" >&2; cat "$tmpdir/run.err" >&2; exit 1; }
[[ -s "$tmpdir/--help" ]] \
    || { echo "FAIL: mem -o --help did not write SAM to file '--help'" >&2
         echo "---- stderr ----" >&2; cat "$tmpdir/run.err" >&2; exit 1; }
grep -q 'Executing in' "$tmpdir/run.err" \
    || { echo "FAIL: banner missing from 'mem -o --help' invocation; pre-scan wrongly matched the option's value" >&2
         echo "---- stderr ----" >&2; cat "$tmpdir/run.err" >&2; exit 1; }

# --- 3. `mem -R --help <ref> <reads>` -- here -R rejects `--help` as an -----
# invalid @RG, but the banner must still be printed first and the run must
# exit non-zero (so the test can't pass if `--help` were ever wrongly
# consumed and the run completed successfully).
set +e
"$bin" mem -R "--help" "$ref" "$reads" > "$tmpdir/rg.out" 2> "$tmpdir/rg.err"
rc=$?
set -e
[[ $rc -ne 0 ]] \
    || { echo "FAIL: 'mem -R --help' exited 0; expected non-zero (invalid @RG)" >&2
         echo "---- stderr ----" >&2; cat "$tmpdir/rg.err" >&2; exit 1; }
grep -q 'Executing in' "$tmpdir/rg.err" \
    || { echo "FAIL: banner missing from 'mem -R --help' invocation; pre-scan wrongly matched the option's value" >&2
         echo "---- stderr ----" >&2; cat "$tmpdir/rg.err" >&2; exit 1; }

# --- 4. `mem --set-as-failed --help <ref> <reads>` -- long option with a ----
# required argument; the banner must still be printed and the run must
# exit non-zero.
set +e
"$bin" mem --set-as-failed "--help" "$ref" "$reads" \
    > "$tmpdir/saf.out" 2> "$tmpdir/saf.err"
rc=$?
set -e
[[ $rc -ne 0 ]] \
    || { echo "FAIL: 'mem --set-as-failed --help' exited 0; expected non-zero (--help is not a valid 'f'|'r' value)" >&2
         echo "---- stderr ----" >&2; cat "$tmpdir/saf.err" >&2; exit 1; }
grep -q 'Executing in' "$tmpdir/saf.err" \
    || { echo "FAIL: banner missing from 'mem --set-as-failed --help' invocation; pre-scan wrongly matched the option's value" >&2
         echo "---- stderr ----" >&2; cat "$tmpdir/saf.err" >&2; exit 1; }

echo "PASS: mem --help short-circuits; mem <opt-with-arg> --help does not"
