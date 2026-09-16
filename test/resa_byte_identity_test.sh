#!/usr/bin/env bash
# test/resa_byte_identity_test.sh
#
# Asserts the correctness contract of `bwa-mem3 re-sa`: resampling an existing
# index's on-disk SA sample table produces a `.bwt.2bit.64` byte-identical to one
# built directly at the target rate with `index -u`, in both directions, and
# leaves alignments unchanged.
#
# Why byte-identity is the gate: re-sa recovers each added SA sample by the same
# LF-walk the resolver uses (densify) or keeps a verbatim subset of the stored
# samples (coarsen). Either way the result must equal a freshly built index at
# that rate down to the byte -- the cp_occ block, sentinel, count[] and trailing
# sa_compx tag included. A byte diff is the earliest, strictest signal that the
# on-disk layout or a synthesised sample is wrong.
#
# Cases:
#   A  densify (stride 8 -> stride 4): re-sa -u 2  ==  index -u 2
#   B  coarsen (stride 8 -> stride 16): re-sa -u 4  ==  index -u 4
#   C  thread determinism: re-sa -u 1 -t 4  ==  re-sa -u 1 -t 1  ==  index -u 1
#   D  no-op: re-sa to the current rate exits 0 and does not touch the file
#   E  alignment equivalence: mem against the original vs the re-sa'd index
#      emits identical records
#   F  inspection (no -u) reports the current rate without touching the file
#   G  argument errors exit non-zero
#
# Usage: test/resa_byte_identity_test.sh <bwa-mem3-binary> <fixtures-dir>

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 <bwa-mem3-binary> <fixtures-dir>" >&2
    exit 2
fi

bin="$1"
fixtures="$2"
[[ -x "$bin" ]] || { echo "FAIL: bwa-mem3 binary not executable at $bin" >&2; exit 1; }
ref_src="$fixtures/synthetic_1mb.fa"
reads="$fixtures/reads.fa"
[[ -s "$ref_src" ]] || { echo "FAIL: missing fixture $ref_src" >&2; exit 1; }
[[ -s "$reads" ]]   || { echo "FAIL: missing fixture $reads" >&2; exit 1; }

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

fail() { echo "FAIL [re-sa byte-identity]: $*" >&2; exit 1; }

# Build a fresh index of the fixture in a private dir and echo the resulting
# .bwt.2bit.64 path. Extra args (e.g. -u 2) are passed through to `index`.
build_index() { # $1 = dir, $2... = extra index args
    local dir="$1"; shift
    mkdir -p "$dir"
    cp "$ref_src" "$dir/ref.fa"
    "$bin" index "$@" "$dir/ref.fa" > "$dir/index.log" 2>&1 \
        || { cat "$dir/index.log" >&2; fail "index $* failed"; }
    echo "$dir/ref.fa.bwt.2bit.64"
}

# --- Case A: densify (stride 8 -> stride 4) --------------------------------
a_re="$(build_index "$tmp/a_re")"                 # default rate (u=3, stride 8)
"$bin" re-sa -u 2 "$tmp/a_re/ref.fa" > "$tmp/a_re/resa.log" 2>&1 \
    || { cat "$tmp/a_re/resa.log" >&2; fail "case A: re-sa -u 2 failed"; }
a_direct="$(build_index "$tmp/a_direct" -u 2)"
cmp -s "$a_re" "$a_direct" \
    || fail "case A: re-sa -u 2 differs from index -u 2 (densify not byte-identical)"

# --- Case B: coarsen (stride 8 -> stride 16) -------------------------------
b_re="$(build_index "$tmp/b_re")"
"$bin" re-sa -u 4 "$tmp/b_re/ref.fa" > "$tmp/b_re/resa.log" 2>&1 \
    || { cat "$tmp/b_re/resa.log" >&2; fail "case B: re-sa -u 4 failed"; }
b_direct="$(build_index "$tmp/b_direct" -u 4)"
cmp -s "$b_re" "$b_direct" \
    || fail "case B: re-sa -u 4 differs from index -u 4 (coarsen not byte-identical)"

# --- Case B2: resample FROM a non-default source stride -------------------
# Every other case starts at the default stride-8 (u=3); this exercises
# shift-dependent logic at disk_sa_compx != 3 in both directions.
b2d="$(build_index "$tmp/b2d" -u 2)"            # source = stride-4 (u=2)
"$bin" re-sa -u 1 "$tmp/b2d/ref.fa" > "$tmp/b2d/resa.log" 2>&1 || fail "case B2: densify from stride-4 failed"
cmp -s "$b2d" "$(build_index "$tmp/b2d_oracle" -u 1)" \
    || fail "case B2: re-sa -u 1 from a stride-4 source differs from index -u 1"
b2c="$(build_index "$tmp/b2c" -u 2)"            # source = stride-4 (u=2)
"$bin" re-sa -u 5 "$tmp/b2c/ref.fa" > "$tmp/b2c/resa.log" 2>&1 || fail "case B2: coarsen from stride-4 failed"
cmp -s "$b2c" "$(build_index "$tmp/b2c_oracle" -u 5)" \
    || fail "case B2: re-sa -u 5 from a stride-4 source differs from index -u 5"

# --- Case C: thread determinism --------------------------------------------
c_t4="$(build_index "$tmp/c_t4")"
c_t1="$(build_index "$tmp/c_t1")"
"$bin" re-sa -u 1 -t 4 "$tmp/c_t4/ref.fa" > "$tmp/c_t4/resa.log" 2>&1 || fail "case C: re-sa -t4 failed"
"$bin" re-sa -u 1 -t 1 "$tmp/c_t1/ref.fa" > "$tmp/c_t1/resa.log" 2>&1 || fail "case C: re-sa -t1 failed"
cmp -s "$c_t4" "$c_t1" || fail "case C: densify output depends on thread count"
c_direct="$(build_index "$tmp/c_direct" -u 1)"
cmp -s "$c_t4" "$c_direct" || fail "case C: threaded densify differs from index -u 1"

# --- Case D: no-op ---------------------------------------------------------
d_idx="$(build_index "$tmp/d")"                   # default rate = u3
before="$(cksum "$d_idx")"
"$bin" re-sa -u 3 "$tmp/d/ref.fa" > "$tmp/d/resa.log" 2>&1 \
    || fail "case D: re-sa to the current rate should exit 0"
grep -q "nothing to do" "$tmp/d/resa.log" \
    || fail "case D: re-sa to the current rate should report a no-op"
after="$(cksum "$d_idx")"
[[ "$before" == "$after" ]] || fail "case D: no-op re-sa mutated the index file"

# --- Case E: alignment equivalence -----------------------------------------
e_orig="$(build_index "$tmp/e_orig")"             # stride 8
e_dense="$(build_index "$tmp/e_dense")"
"$bin" re-sa -u 2 "$tmp/e_dense/ref.fa" > "$tmp/e_dense/resa.log" 2>&1 || fail "case E: re-sa failed"
# Compare alignment RECORDS only: the @PG header's CL: field embeds the index
# path, which necessarily differs between the two index dirs.
"$bin" mem "$tmp/e_orig/ref.fa"  "$reads" 2>/dev/null | grep -v '^@' > "$tmp/e_orig.sam"
"$bin" mem "$tmp/e_dense/ref.fa" "$reads" 2>/dev/null | grep -v '^@' > "$tmp/e_dense.sam"
[[ -s "$tmp/e_orig.sam" ]] || fail "case E: no alignment records produced"
cmp -s "$tmp/e_orig.sam" "$tmp/e_dense.sam" \
    || fail "case E: alignments differ between the original and the re-sa'd index"

# --- Case F: inspection (no -u) --------------------------------------------
f_idx="$(build_index "$tmp/f")"                   # default rate = stride 8
before_f="$(cksum "$f_idx")"
"$bin" re-sa "$tmp/f/ref.fa" > "$tmp/f/inspect.log" 2>&1 \
    || fail "case F: inspection (no -u) should exit 0"
grep -qE "SA sample rate 1/8 \(shift 3\)" "$tmp/f/inspect.log" \
    || { cat "$tmp/f/inspect.log" >&2; fail "case F: inspection must report the current rate (1/8)"; }
[[ "$before_f" == "$(cksum "$f_idx")" ]] || fail "case F: inspection mutated the index file"

# --- Case G: argument errors -----------------------------------------------
"$bin" re-sa -u 9 "$tmp/f/ref.fa"   > /dev/null 2>&1 && fail "case G: out-of-range -u must fail"
"$bin" re-sa -u 2 /nonexistent/idx  > /dev/null 2>&1 && fail "case G: missing index must fail"
"$bin" re-sa                        > /dev/null 2>&1 && fail "case G: missing idxbase must fail"

echo "PASS: re-sa produces a byte-identical index to \`index -u\` (densify + coarsen), is thread-deterministic, no-ops on the current rate, preserves alignments, reports the rate in inspection mode, and rejects bad arguments"
