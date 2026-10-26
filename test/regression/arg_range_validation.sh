#!/usr/bin/env bash
# test/regression/arg_range_validation.sh
#
# Asserts that `mem` rejects, at parse time and with a clear ERROR: message, a
# non-positive gap-extension penalty (-E), a non-positive insert-size standard
# deviation (-I mean,std), a negative gap-open penalty (-O) and a band width
# (-w) outside 0..536870911:
#   - -E 0 divides by zero in cal_max_gap; the (int) cast of the inf/NaN result
#     is UB that resolves differently per architecture (arch-divergent bands).
#   - -I mean,0 divides by zero in mem_pair; the NaN collapses every pair score
#     to 0 and silently disables pair-aware placement.
#   - -O -1 makes a gap open cheaper than an extension, which the SIMD
#     Smith-Waterman kernels do not compute exactly (see the -O checks below).
#   - -w -1 inverts the reference window and aborts; -w -50 lets the banded
#     kernels store past their rows (see the -w checks below).
#   - a scoring the aligner cannot carry (see the scoring checks below): an -A
#     outside 1..127, an effective -B outside 0..128, a negative -d or -U, or
#     an effective gap open plus extension sum out of range.
#   - a -T that is not a complete int (atoi read garbage as 0 and wrapped
#     out-of-range values).
#
# Validation fires during option parsing, before the reference is opened, so
# /dev/null placeholders suffice for the rejections. The message (not just a
# non-zero exit) is asserted, so a non-zero exit from a later missing-input
# error cannot masquerade as the rejection under test. The last check aligns
# one read against the phiX fixture (FIXTURES, default test/fixtures) to show
# an accepted -T reaches the output filter.

set -euo pipefail

BWA_MEM3="${BWA_MEM3:-$(dirname "$0")/../../bwa-mem3}"
if [[ ! -x "$BWA_MEM3" ]]; then
    echo "FAIL: BWA_MEM3 ('$BWA_MEM3') is not an executable binary" >&2
    exit 1
fi

OUT_DIR="$(mktemp -d -t bwamem3-argval-XXXXXX)"
trap 'rm -rf "$OUT_DIR"' EXIT

# check_reject <description> <expected-stderr-substring> <mem-args...>
check_reject() {
    local desc=$1 pattern=$2
    shift 2
    local rc=0
    "$BWA_MEM3" mem "$@" /dev/null /dev/null /dev/null \
        > "$OUT_DIR/out" 2> "$OUT_DIR/err" || rc=$?
    if [[ $rc -eq 0 ]]; then
        echo "FAIL: $desc: expected non-zero exit, got 0" >&2
        cat "$OUT_DIR/err" >&2
        exit 1
    fi
    if ! grep -qF "$pattern" "$OUT_DIR/err"; then
        echo "FAIL: $desc: stderr does not contain '$pattern'" >&2
        cat "$OUT_DIR/err" >&2
        exit 1
    fi
}

check_reject "-E 0" "ERROR: -E gap-extension penalty must be a positive integer" -E 0
check_reject "-E 5,0" "ERROR: -E gap-extension penalty must be a positive integer" -E 5,0
# An out-of-range -E value must be rejected before it is narrowed to int: 2^32+1
# narrows to 1, which would pass the old >= 1 check and silently apply a
# gap-extension penalty different from the one requested (and feed cal_max_gap's
# divide with a wrong value). Reject on the parsed long, not the wrapped int.
check_reject "-E 4294967297" "ERROR: -E gap-extension penalty must be a positive integer" -E 4294967297
check_reject "-E 5,4294967297" "ERROR: -E gap-extension penalty must be a positive integer" -E 5,4294967297
check_reject "-I 300,0" "ERROR: -I standard deviation must be a positive number" -I 300,0

# A negative gap-open penalty makes opening a gap cheaper than extending one
# (o + e < e), which the SIMD Smith-Waterman kernels do not support: the
# mate-rescue kernels (kswv) return scores, ends and suboptimal scores that
# differ from the exact affine-gap DP, and the striped ksw_u8/ksw_i16 return
# differing scores on some of the same inputs. Upstream
# bwa/bwa-mem2 accept the value and return those inexact results; reject it at
# parse instead. -O 0 stays valid. The parse is the -E one: a complete token,
# no sign, no ERANGE, at most INT_MAX, so 8,-1 (whose second value the old
# parser silently dropped, leaving the insertion open at 8) is rejected too.
check_reject "-O -1" "ERROR: -O gap-open penalty must be a non-negative integer" -O -1
check_reject "-O 8,-1" "ERROR: -O gap-open penalty must be a non-negative integer" -O 8,-1
check_reject "-O -1,6" "ERROR: -O gap-open penalty must be a non-negative integer" -O -1,6
check_reject "-O 6x" "ERROR: -O gap-open penalty must be a non-negative integer" -O 6x
check_reject "-O 4294967297" "ERROR: -O gap-open penalty must be a non-negative integer" -O 4294967297
# strtol() returns 0 for a token with no digits, and 0 is a valid -O, so an
# empty or missing first value must be rejected on the parse itself: "" would
# otherwise run with both opens at 0, and ",1" with a deletion open of 0.
check_reject "-O ''" "ERROR: -O gap-open penalty must be a non-negative integer" -O ''
check_reject "-O ,1" "ERROR: -O gap-open penalty must be a non-negative integer" -O ,1

# A negative band width turns cal_max_gap's `w << 1` cap negative: -w -1 inverts
# the reference window and aborts, -w -50 lets the banded kernels store past
# their rows and still exits 0, and -w -2147483648 segfaults. The widest
# w-scaled int is `w << 2`, so the top of the range is INT_MAX >> 2 = 536870911.
check_reject "-w -1" "ERROR: -w band width must be an integer in 0..536870911" -w -1
check_reject "-w -50" "ERROR: -w band width must be an integer in 0..536870911" -w -50
check_reject "-w -2147483648" "ERROR: -w band width must be an integer in 0..536870911" -w -2147483648
check_reject "-w 536870912" "ERROR: -w band width must be an integer in 0..536870911" -w 536870912
check_reject "-w 100x" "ERROR: -w band width must be an integer in 0..536870911" -w 100x
check_reject "-w ''" "ERROR: -w band width must be an integer in 0..536870911" -w ''
# The ends of the range still parse: the run then fails on the /dev/null
# reference, but not on -w. Require that index-load error too, so a run that
# failed on something else before reaching the index proves nothing.
for w in 0 536870911; do
    rc=0
    "$BWA_MEM3" mem -w "$w" /dev/null /dev/null /dev/null \
        > "$OUT_DIR/out" 2> "$OUT_DIR/err" || rc=$?
    if grep -qF "ERROR: -w band width" "$OUT_DIR/err"; then
        echo "FAIL: -w $w: rejected a band width inside 0..536870911 (exit $rc)" >&2
        cat "$OUT_DIR/err" >&2
        exit 1
    fi
    if ! grep -qF "Unable to open the file: /dev/null" "$OUT_DIR/err"; then
        echo "FAIL: -w $w: the run never reached loading the /dev/null index, so the parse is unproven (exit $rc)" >&2
        cat "$OUT_DIR/err" >&2
        exit 1
    fi
done

# A non-finite -I mean or std (nan/inf, or an overflowing token that strtod maps
# to inf) reaches an (int) cast whose result is undefined and arch-divergent;
# reject it at parse rather than casting inf/NaN downstream.
check_reject "-I nan,1" "ERROR: -I mean insert size must be a finite number" -I nan,1
check_reject "-I inf" "ERROR: -I mean insert size must be a finite number" -I inf
check_reject "-I 300,1e400" "ERROR: -I standard deviation must be a positive number" -I 300,1e400
check_reject "-I 300,1,1e400" "ERROR: -I insert-size max must be a finite number" -I 300,1,1e400

# strtol()/strtod() stop at the first non-numeric byte and the old parsers
# applied the valid prefix, silently dropping the remainder. Require complete
# consumption of the token so a typo fails loudly instead of running with a
# value the user did not write.
check_reject "-E 5abc" "ERROR: -E gap-extension penalty must be a positive integer" -E 5abc
check_reject "-E 5,3xyz" "ERROR: -E gap-extension penalty must be a positive integer" -E 5,3xyz
check_reject "-I 300,50xy" "ERROR: -I expects mean" -I 300,50xy
check_reject "-I 300,50,junk" "ERROR: -I expects mean" -I 300,50,junk
# strtod() returns 0 for a token with no digits, so a missing mean (",50") would
# run with a mean insert size of 0 and the given std; require the mean itself.
check_reject "-I ,50" "ERROR: -I expects mean" -I ,50

# A mean (or explicit max/min) large enough that the rounded insert-size bound
# overflows int makes the (int) cast UB; reject on the double before casting.
check_reject "-I 3000000000" "ERROR: -I mean/std imply an insert-size bound outside" -I 3000000000
check_reject "-I 300,50,4000000000" "ERROR: -I insert-size max must be in" -I 300,50,4000000000

# -c is the seed max-occurrence cap and a divisor in the SA-resolve step math
# (`p->s / opt->max_occ`), so max_occ == 0 is a hard divide-by-zero and a
# negative value silently yields no seeds; atoi() would have accepted both, plus
# trailing garbage and an out-of-range value that narrows to a valid-looking
# int. Reject all at parse, before SA resolution.
check_reject "-c 0" "ERROR: -c max occurrences must be a positive integer" -c 0
check_reject "-c -1" "ERROR: -c max occurrences must be a positive integer" -c -1
check_reject "-c 5abc" "ERROR: -c max occurrences must be a positive integer" -c 5abc
check_reject "-c 4294967297" "ERROR: -c max occurrences must be a positive integer" -c 4294967297

# check_accept <description> <mem-args...>: the options pass validation. The
# run still fails (the /dev/null index does not exist), so require evidence it
# got as far as loading the index -- an error naming one of the index files
# (or, under --meth, the seed-index probe) -- and that no scoring rejection
# fired on the way.
check_accept() {
    local desc=$1
    shift
    local rc=0
    "$BWA_MEM3" mem "$@" /dev/null /dev/null /dev/null \
        > "$OUT_DIR/out" 2> "$OUT_DIR/err" || rc=$?
    if grep -qE "ERROR: (-[ABdTU] |the effective )" "$OUT_DIR/err"; then
        echo "FAIL: $desc: a valid scoring was rejected" >&2
        cat "$OUT_DIR/err" >&2
        exit 1
    fi
    if ! grep -qE "/dev/null\.[a-z0-9]|seed index '/dev/null" "$OUT_DIR/err"; then
        echo "FAIL: $desc: the run never reached index loading (exit $rc)" >&2
        cat "$OUT_DIR/err" >&2
        exit 1
    fi
}

# Scoring. The kernels store the match score and the negated mismatch penalty
# in int8_t (bwa_fill_scmat, kswv, the --meth matrices), so a wider value
# wraps: -A 128 becomes -128, and -A 256 zeroes kswv's maximum substitution
# score, which its score2 window divides by (SIGFPE at the AVX-512BW tier). A
# match score of 0 or less is refused too: nothing can score above zero, so
# the same divide traps on x86, and the pairing and mapping-quality formulas
# overflow on negative scores, differently per SIMD tier. -A is held to 1..127
# at parse, which also keeps update_a()'s scaled defaults inside int and their
# signs unchanged.
A_RANGE="ERROR: -A match score must be an integer in 1..127"
check_reject "-A 128" "$A_RANGE" -A 128
check_reject "-A 256" "$A_RANGE" -A 256
check_reject "-A 20000000" "$A_RANGE" -A 20000000
check_reject "-A 2x" "$A_RANGE" -A 2x
check_reject "-A 0" "$A_RANGE" -A 0
check_reject "-A -1" "$A_RANGE" -A -1
check_reject "-A 0 (all pinned)" "$A_RANGE" -A 0 -B 4 -O 6 -E 1 -d 100 -L 5 -U 17
check_reject "-A -1 (all pinned)" "$A_RANGE" -A -1 -B 4 -O 6 -E 1 -d 100 -L 5 -U 17
check_reject "-A -128" "$A_RANGE" -A -128
# The mismatch penalty is checked after -A scales it, so the message reports
# the effective value: -A 33 scales the default -B 4 to 132. A negative
# penalty scores a mismatch above zero, which the 8-bit kernels' admission
# tests and the mapping-quality model do not allow.
B_RANGE="outside 0..128"
check_reject "-B 129" "ERROR: the effective mismatch penalty is 129, $B_RANGE" -B 129
check_reject "-B -1" "ERROR: the effective mismatch penalty is -1, $B_RANGE" -B -1
check_reject "-B -4" "ERROR: the effective mismatch penalty is -4, $B_RANGE" -B -4
check_reject "-B -127" "ERROR: the effective mismatch penalty is -127, $B_RANGE" -B -127
check_reject "-A 4 -B -4" "ERROR: the effective mismatch penalty is -4, $B_RANGE" -A 4 -B -4
check_reject "-A 33" "ERROR: the effective mismatch penalty is 132, $B_RANGE" -A 33
check_reject "-A 127" "ERROR: the effective mismatch penalty is 508, $B_RANGE" -A 127
check_reject "-B 4x" "ERROR: -B mismatch penalty must be an integer" -B 4x
# atoi() mapped garbage to 0 (output every alignment) and wrapped out-of-range
# values; -T now takes a complete int.
check_reject "-T 30x" "ERROR: -T minimum score to output must be an integer" -T 30x
check_reject "-T abc" "ERROR: -T minimum score to output must be an integer" -T abc
check_reject "-T 2147483648" "ERROR: -T minimum score to output must be an integer" -T 2147483648
check_reject "-T -2147483649" "ERROR: -T minimum score to output must be an integer" -T -2147483649
# The 8-bit banded extension asserts zdrop >= 0 and aborts.
check_reject "-d -1" "ERROR: -d z-drop must be a non-negative integer" -d -1
check_reject "-d 5x" "ERROR: -d z-drop must be a non-negative integer" -d 5x
# Paired-end scoring subtracts -U from alignment scores; near INT_MIN that
# overflows, and any negative value rewards leaving a pair unpaired.
check_reject "-U -1" "ERROR: -U unpaired penalty must be a non-negative integer" -U -1
check_reject "-U -2147483647" "ERROR: -U unpaired penalty must be a non-negative integer" -U -2147483647
check_reject "-U 5x" "ERROR: -U unpaired penalty must be a non-negative integer" -U 5x
# Cells outside ksw_global2's band hold -2^30, and a gap state derived from one
# subtracts o + e twice more, so each open + extension sum must be <= 2^29 - 1.
check_reject "-O 536870911" "ERROR: the effective gap open plus extension penalty is 536870912 for deletions" -O 536870911
check_reject "-O 6,536870911" "and 536870912 for insertions" -O 6,536870911
check_reject "-E 536870906" "ERROR: the effective gap open plus extension penalty is 536870912 for deletions" -E 536870906

# The last accepted values, the shipped presets and --meth's scaled defaults
# pass validation.
check_accept "defaults"
check_accept "-A 1" -A 1
check_accept "-A 2" -A 2
check_accept "-A 32" -A 32
check_accept "-A 127 -B 4" -A 127 -B 4
check_accept "-B 0" -B 0
check_accept "-B 128" -B 128
check_accept "-d 0 -U 0" -d 0 -U 0
check_accept "-T 0" -T 0
check_accept "-T -5" -T -5
check_accept "-T 2147483647" -T 2147483647
check_accept "-T -2147483648" -T -2147483648
check_accept "-O 536870910" -O 536870910
check_accept "-x pacbio" -x pacbio
check_accept "-x ont2d" -x ont2d
check_accept "-x intractg" -x intractg
check_accept "--meth" --meth
check_accept "--meth -A 2" --meth -A 2
check_accept "--meth-scoring genomic -A 32" --meth --meth-scoring genomic -A 32

# An accepted -T must reach the output filter, not just get past parsing: a
# 100 bp read copied from phiX scores exactly 100 at the default scoring, so it
# is reported at -T 100 and unmapped (flag 0x4) at -T 101.
FIXTURES="${FIXTURES:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../fixtures" && pwd)}"
if [[ ! -s "$FIXTURES/phix.fa" ]]; then
    echo "FAIL: phix.fa missing: $FIXTURES/phix.fa" >&2
    exit 1
fi
cp "$FIXTURES/phix.fa" "$OUT_DIR/phix.fa"
if ! "$BWA_MEM3" index "$OUT_DIR/phix.fa" > /dev/null 2>&1; then
    echo "FAIL: could not index $OUT_DIR/phix.fa" >&2
    exit 1
fi
read_seq="$(grep -v '^>' "$OUT_DIR/phix.fa" | tr -d '\n' | cut -c1001-1100)"
if [[ ${#read_seq} -ne 100 ]]; then
    echo "FAIL: could not cut a 100 bp read from phix.fa (got ${#read_seq} bp)" >&2
    exit 1
fi
printf '>r1\n%s\n' "$read_seq" > "$OUT_DIR/read.fa"
# Print the FLAG and AS:i of the read's primary record at threshold $1.
primary_flag_and_score() {
    "$BWA_MEM3" mem -T "$1" "$OUT_DIR/phix.fa" "$OUT_DIR/read.fa" 2> "$OUT_DIR/err" > "$OUT_DIR/out.sam"
    awk -F'\t' '!/^@/ && $2 < 256 { as = "none"; for (i = 12; i <= NF; i++) if ($i ~ /^AS:i:/) as = $i; print $2, as; exit }' "$OUT_DIR/out.sam"
}
got="$(primary_flag_and_score 100)"
if [[ $got != "0 AS:i:100" ]]; then
    echo "FAIL: -T 100: expected the exact 100 bp read reported forward with AS:i:100, got '$got'" >&2
    cat "$OUT_DIR/err" >&2
    exit 1
fi
got="$(primary_flag_and_score 101)"
if [[ ${got%% *} != 4 ]]; then
    echo "FAIL: -T 101: expected the read reported unmapped (flag 4), got '$got'" >&2
    cat "$OUT_DIR/err" >&2
    exit 1
fi

echo "PASS: arg_range_validation (-E, -I, and -c reject non-positive, -O rejects negative, -w rejects values outside 0..536870911, -A, -B, -d, -U and the effective gap penalties are range-checked, -T must be an int and is applied to output, -O and -I reject a first value with no digits, and all reject non-finite, out-of-range, and trailing-garbage values)"
