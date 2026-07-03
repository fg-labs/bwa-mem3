#!/usr/bin/env bash
# test/fast_preset_test.sh
#
# Asserts that `bwa-mem3 mem --fast` resolves the characterized speed levers
# (-m 10 -y 0 --min-ext-len 30 --smem-dedup --skip-contained-ext
# --max-extend-chains 5, plus -s 2, --max-extend-chains 10, and
# --extend-mate-concordant under --meth), that explicit user flags override the
# preset, and that the default path is untouched when --fast is absent.
# --skip-contained-ext no-ops under --meth (internal gate), so it is omitted from
# the meth audit line; --max-extend-chains applies under --meth too but at a
# higher cap of 10 (non-meth keeps 5). --extend-mate-concordant is meth-only
# (recovers the bisulfite pairing regression the chain cap otherwise causes) and
# must be absent from the non-meth audit line.
#
# The assertion surface is the audit line main_mem prints to stderr when
# --fast is active: it reports the *resolved* mem_opt_t values, so an explicit
# override is visible in the line itself.
#
# Usage: test/fast_preset_test.sh <bwa-mem3-binary> <fixtures-dir>
set -euo pipefail

[[ $# -eq 2 ]] || { echo "usage: $0 <bwa-mem3-binary> <fixtures-dir>" >&2; exit 2; }
bin="$1"; fixtures="$2"
src_ref="$fixtures/phix.fa"; reads="$fixtures/reads.fa"

[[ -x "$bin" ]]     || { echo "FAIL: binary not executable: $bin" >&2; exit 1; }
[[ -s "$src_ref" ]] || { echo "FAIL: phix.fa missing: $src_ref" >&2; exit 1; }
[[ -s "$reads" ]]   || { echo "FAIL: reads.fa missing: $reads" >&2; exit 1; }

err="$(mktemp)"; mdir="$(mktemp -d)"
trap 'rm -f "$err"; rm -rf "$mdir"' EXIT

# Index a private copy of phiX in temp space so the test never writes derived
# index files into the (possibly read-only or shared) fixtures tree. The meth
# path below does the same; keep both hermetic. The copy is always fresh, so
# index unconditionally.
ref="$mdir/phix.fa"
cp "$src_ref" "$ref"
"$bin" index "$ref" >/dev/null 2>&1 || { echo "FAIL: index phix.fa" >&2; exit 1; }

# Run mem with extra args; echo the resolved --fast audit line (empty if none).
fast_line() {
    "$bin" mem "$@" "$ref" "$reads" >/dev/null 2>"$err" \
        || { echo "FAIL: mem exited nonzero (args: $*)" >&2; cat "$err" >&2; exit 1; }
    grep -E '^\[M::main_mem\] --fast:' "$err" || true
}

# 1. Bundle correctness (non-meth).
line="$(fast_line --fast)"
[[ "$line" == *"-m 10"* && "$line" == *"-y 0"* \
   && "$line" == *"--min-ext-len 30"* && "$line" == *"--smem-dedup"* \
   && "$line" == *"--skip-contained-ext"* \
   && "$line" == *"--max-extend-chains 5"* ]] \
    || { echo "FAIL: --fast bundle wrong: '$line'" >&2; exit 1; }
[[ "$line" != *"-s "* ]] \
    || { echo "FAIL: non-meth --fast must not set -s: '$line'" >&2; exit 1; }
[[ "$line" != *"--extend-mate-concordant"* ]] \
    || { echo "FAIL: non-meth --fast must not enable --extend-mate-concordant: '$line'" >&2; exit 1; }
echo "OK:   --fast bundle resolves -m 10 -y 0 --min-ext-len 30 --smem-dedup --skip-contained-ext --max-extend-chains 5 (no -s, no mate-concordant)"

# 2. Override precedence: an explicit flag wins *only* for its own field; the
#    rest of the preset (including the unconditional --smem-dedup) must survive.
line="$(fast_line --fast -m 30)"
[[ "$line" == *"-m 30"* && "$line" == *"-y 0"* \
   && "$line" == *"--min-ext-len 30"* && "$line" == *"--smem-dedup"* \
   && "$line" == *"--skip-contained-ext"* \
   && "$line" == *"--max-extend-chains 5"* ]] \
    || { echo "FAIL: explicit -m 30 should only override -m: '$line'" >&2; exit 1; }
line="$(fast_line --fast -y 5)"
[[ "$line" == *"-m 10"* && "$line" == *"-y 5"* \
   && "$line" == *"--min-ext-len 30"* && "$line" == *"--smem-dedup"* \
   && "$line" == *"--skip-contained-ext"* \
   && "$line" == *"--max-extend-chains 5"* ]] \
    || { echo "FAIL: explicit -y 5 should only override -y: '$line'" >&2; exit 1; }
line="$(fast_line --fast --min-ext-len 45)"
[[ "$line" == *"-m 10"* && "$line" == *"-y 0"* \
   && "$line" == *"--min-ext-len 45"* && "$line" == *"--smem-dedup"* \
   && "$line" == *"--skip-contained-ext"* \
   && "$line" == *"--max-extend-chains 5"* ]] \
    || { echo "FAIL: explicit --min-ext-len 45 should only override min-ext-len: '$line'" >&2; exit 1; }
# --max-extend-chains is overridable under --fast (respects an explicit value).
line="$(fast_line --fast --max-extend-chains 8)"
[[ "$line" == *"-m 10"* && "$line" == *"-y 0"* \
   && "$line" == *"--min-ext-len 30"* && "$line" == *"--smem-dedup"* \
   && "$line" == *"--max-extend-chains 8"* ]] \
    || { echo "FAIL: explicit --max-extend-chains 8 should only override that lever: '$line'" >&2; exit 1; }
echo "OK:   explicit -m/-y/--min-ext-len/--max-extend-chains override only their field; rest of preset survives"

# 3. Default contract: no --fast => no audit line at all.
"$bin" mem "$ref" "$reads" >/dev/null 2>"$err" || { echo "FAIL: plain mem nonzero" >&2; exit 1; }
! grep -qE '^\[M::main_mem\] --fast:' "$err" \
    || { echo "FAIL: audit line present without --fast" >&2; exit 1; }
echo "OK:   no --fast => no audit line (default path untouched)"

# 4. Meth path: --fast --meth additionally sets -s 2. Build a meth index on a
#    copy of phiX; if meth indexing is unavailable in this build, SKIP.
cp "$ref" "$mdir/ref.fa"
if "$bin" index --meth "$mdir/ref.fa" >/dev/null 2>&1; then
    "$bin" mem --meth --fast -t 1 "$mdir/ref.fa" "$reads" >/dev/null 2>"$err" \
        || { echo "FAIL: mem --meth --fast nonzero" >&2; cat "$err" >&2; exit 1; }
    line="$(grep -E '^\[M::main_mem\] --fast:' "$err" || true)"
    [[ "$line" == *"-s 2"* ]] \
        || { echo "FAIL: --fast --meth should resolve -s 2: '$line'" >&2; exit 1; }
    [[ "$line" != *"--skip-contained-ext"* ]] \
        || { echo "FAIL: --skip-contained-ext no-ops under --meth; must be absent from audit line: '$line'" >&2; exit 1; }
    [[ "$line" == *"--max-extend-chains 10"* ]] \
        || { echo "FAIL: --max-extend-chains applies under --meth (cap 10); must be in audit line: '$line'" >&2; exit 1; }
    [[ "$line" == *"--extend-mate-concordant"* ]] \
        || { echo "FAIL: --fast --meth must enable --extend-mate-concordant: '$line'" >&2; exit 1; }
    echo "OK:   --fast --meth additionally sets -s 2 and --extend-mate-concordant (skip-contained-ext omitted meth-gated, --max-extend-chains raised to 10)"
    # Explicit -s wins even under --meth (src/fastmap.cpp: -s 2 is gated on !opt0.split_width).
    "$bin" mem --meth --fast -s 7 -t 1 "$mdir/ref.fa" "$reads" >/dev/null 2>"$err" \
        || { echo "FAIL: mem --meth --fast -s 7 nonzero" >&2; cat "$err" >&2; exit 1; }
    line="$(grep -E '^\[M::main_mem\] --fast:' "$err" || true)"
    [[ "$line" == *"-s 7"* ]] \
        || { echo "FAIL: explicit -s 7 should win under --meth: '$line'" >&2; exit 1; }
    echo "OK:   explicit -s 7 overrides --fast --meth's -s 2"
else
    echo "SKIP: index --meth unavailable; meth -s 2 case not checked" >&2
fi

echo "PASS: fast_preset_test"
