#!/usr/bin/env bash
# test/regression/rescue_kmer_options.sh
#
# Regression: --rescue-kmer / --rescue-band must reject malformed and
# out-of-range values instead of silently substituting a working default.
#
# Both options were originally parsed with atoi(), which maps every unparseable
# string to 0. For --rescue-kmer, 0 means "off" AND counts as an explicitly-set
# value, so `--fast --rescue-kmer=x` would disable the single largest lever in
# --fast while looking accepted. For --rescue-band, both unparseable text and a
# negative value produced a number the kernel's `band > 0 ? band : 50` guard
# then turned back into the default -- so a typo read as a deliberate band width
# but quietly ran the default one. Neither failure is visible in the output:
# --rescue-kmer only moves MAPQ on rescued reads.
#
# This test pins the contract:
#   * --rescue-kmer  -> integer in 0..16 (0 = off); anything else is a hard error
#   * --rescue-band  -> integer in 1..1000000 bp; anything else is a hard error
#   * the K bound matches the uint32 k-mer encoder, so K=17.. is REJECTED rather
#     than silently clamped to 16 and reported as a distinct setting
#   * the value must be attached with `=`: --rescue-kmer takes an *optional*
#     argument, so `--rescue-kmer 6` selects the bare default and leaves 6 as a
#     positional argument. This is the whole reason the docs and --help spell the
#     --fast expansion as --rescue-kmer=6.
#
# Inputs:
#   BWA_MEM3 — path to bwa-mem3 binary
#   FIXTURES — directory containing phix.fa and reads.fa (default: test/fixtures)
set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
FIXTURES="${FIXTURES:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../fixtures" && pwd)}"

src_ref="$FIXTURES/phix.fa"
reads="$FIXTURES/reads.fa"
[[ -x "$BWA_MEM3" ]] || {
    echo "FAIL: binary not executable: $BWA_MEM3" >&2
    exit 1
}
[[ -s "$src_ref" ]] || {
    echo "FAIL: phix.fa missing: $src_ref" >&2
    exit 1
}
[[ -s "$reads" ]] || {
    echo "FAIL: reads.fa missing: $reads" >&2
    exit 1
}

# Index a private copy so the test never writes into the fixtures tree.
mdir="$(mktemp -d)"
err="$mdir/err.log"
trap 'rm -rf "$mdir"' EXIT
ref="$mdir/phix.fa"
cp "$src_ref" "$ref"
"$BWA_MEM3" index "$ref" > /dev/null 2>&1 || {
    echo "FAIL: index phix.fa" >&2
    exit 1
}

fails=0

# $1 = full option word (e.g. --rescue-kmer=abc), $2 = expected diagnostic substring.
# Always uses the `=` form so a leading '-' in the value is not taken for another
# option.
reject() {
    local optword="$1" want="$2" rc=0
    "$BWA_MEM3" mem "$optword" "$ref" "$reads" > /dev/null 2> "$err" || rc=$?
    if [[ "$rc" -eq 0 ]]; then
        echo "FAIL: '$optword' was accepted; expected a non-zero exit"
        fails=$((fails + 1))
    elif ! grep -q "$want" "$err"; then
        echo "FAIL: '$optword' exited $rc but without the expected diagnostic ('$want')"
        fails=$((fails + 1))
    else
        echo "  ok: '$optword' rejected"
    fi
}

# $1 = full option word that must be ACCEPTED (the run itself must succeed).
accept() {
    local optword="$1"
    if "$BWA_MEM3" mem "$optword" "$ref" "$reads" > /dev/null 2> "$err"; then
        echo "  ok: '$optword' accepted"
    else
        echo "FAIL: '$optword' was rejected; expected it to run"
        cat "$err" >&2
        fails=$((fails + 1))
    fi
}

kmer_err='ERROR: --rescue-kmer requires an integer in 0..16'
band_err='ERROR: --rescue-band requires an integer in 1..1000000 bp'

# --- --rescue-kmer -----------------------------------------------------------
reject "--rescue-kmer=abc" "$kmer_err" # not a number at all -> would have meant "off"
reject "--rescue-kmer=6x" "$kmer_err"  # trailing junk: partial parse must not pass
reject "--rescue-kmer=-1" "$kmer_err"  # negative
# Above the uint32 k-mer code width. Rejected rather than clamped: a clamped K=20
# behaves exactly as K=16, so accepting it would report a setting that does not
# exist. This is also what the documented 1..16 range promises.
reject "--rescue-kmer=17" "$kmer_err"
reject "--rescue-kmer=20" "$kmer_err"

accept "--rescue-kmer"    # bare: selects the default K
accept "--rescue-kmer=0"  # explicit off
accept "--rescue-kmer=1"  # low edge
accept "--rescue-kmer=16" # high edge -- must be INSIDE the accepted range

# --- --rescue-band -----------------------------------------------------------
reject "--rescue-band=abc" "$band_err"
reject "--rescue-band=50bp" "$band_err" # trailing suffix: the option takes plain bp
reject "--rescue-band=-5" "$band_err"   # negative
# 0 is rejected rather than accepted-and-ignored: it reads as "no band" but the
# kernel's `band > 0 ? band : 50` guard would have silently run the 50 bp default.
reject "--rescue-band=0" "$band_err"
reject "--rescue-band=1000001" "$band_err" # above the documented cap

accept "--rescue-band=1"
accept "--rescue-band=1000000"

# --- the value must be attached with `=` -------------------------------------
# --rescue-kmer takes an optional argument, so getopt_long never consumes a
# following word: `mem --rescue-kmer 6 ref reads` leaves THREE positionals, and
# `6` is taken as the index base. Pinning this keeps the docs and the --help
# --fast expansion honest about the `=` form.
rc=0
"$BWA_MEM3" mem --rescue-kmer 6 "$ref" "$reads" > /dev/null 2> "$err" || rc=$?
if [[ "$rc" -eq 0 ]]; then
    echo "FAIL: '--rescue-kmer 6' (space-separated) unexpectedly succeeded;"
    echo "      it should leave 6 as a positional argument"
    fails=$((fails + 1))
elif ! grep -q '6\.bwt\.2bit\.64' "$err"; then
    echo "FAIL: '--rescue-kmer 6' should have tried to open '6' as the index base"
    cat "$err" >&2
    fails=$((fails + 1))
else
    echo "  ok: '--rescue-kmer 6' leaves 6 as a positional (the '=' form is required)"
fi

if [[ "$fails" -ne 0 ]]; then
    echo "FAIL: rescue-kmer option validation ($fails failure(s))"
    exit 1
fi
echo "PASS: --rescue-kmer/--rescue-band reject malformed and out-of-range values"
