#!/usr/bin/env bash
# test/regression/meth_5base_alias.sh
#
# Regression: --meth=5base (and its alias 5-base) selects the TAPS chemistry;
# see docs/src/methylation/5base.md for why. Every record must be identical to
# --meth=taps.
#
# Byte-identity alone would also pass if both values were silently mapped to
# em-seq, so the fixture pins the polarity too: one pair from each original
# strand, each mate carrying exactly one converted CpG, which must be that
# mate's only methylated CpG call ('Z') under 5base and its only unmethylated
# one ('z') under bare --meth. That check needs mawk and reports SKIP for itself
# alone when mawk is missing; the identity and rejection checks always run, and
# so does the check that every run (5base, its spellings, the em-seq control
# and both last-wins orders) gives the same four placed records (POS and CIGAR), so two
# header-only outputs cannot compare equal: taps and the em-seq control are
# pinned to them, and every other run is compared record for record with one
# of those two.
#
# The retained CpGs pin it from the other side: each is unmethylated ('z')
# under 5base and methylated ('Z') under bare --meth, mate for mate. And
# `--meth=taps --meth-scoring neutral` must equal plain --meth=taps, so a change
# to TAPS's default scoring (which 5base inherits) fails here.
#
# The last --meth on the command line wins, as for any option: a bare --meth
# after --meth=5base is em-seq (it used to keep 5base), and --meth=5base after a
# bare --meth is 5base. Both orders are checked record for record.
#
# Also asserted: an unknown chemistry is rejected with a message naming 5base,
# and so is `--meth=illumina` (Illumina also sells bisulfite kits, whose calls
# 5base would invert). The advice for `--meth 5base` (no '=') is checked with
# the other orphaned values in meth_flag_footguns.sh.
#
# Inputs:
#   BWA_MEM3 — the bwa-mem3 binary under test: a bare command name resolved via
#              PATH, or a path to the binary (a relative one is made absolute
#              before the `cd "$WORK"` below)
set -euo pipefail
# shellcheck source=test/framework/dna.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../framework" && pwd)/dna.sh"
# Checked explicitly rather than with ${VAR:?}, which exits with a bare shell
# error (and before any ERR trap runs), so the failure would carry no FAIL: line.
if [[ -z ${BWA_MEM3:-} ]]; then
    echo "FAIL: BWA_MEM3 must be set" >&2
    exit 1
fi
case "$BWA_MEM3" in
    */*) BWA_MEM3="$(cd "$(dirname "$BWA_MEM3")" && pwd)/$(basename "$BWA_MEM3")" ;;
esac

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
fail() {
    echo "FAIL: $*" >&2
    exit 1
}

# Deterministic 1500 bp reference (same generator as meth_output_integrity).
REF=TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTAGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGAATCGTCCCTTCTAACGTAGTTCGCATAGTTCCCGTCCGTAGCCGGACTATTCGAACACCCAGTATTCGATTAACTCGGGCTTGACGTATTAGAGGCGTTAGTGTGCCAGGTAAGATACGCCAACGGAATTAACCTCTGTGACACTCCGCGGAGCCTTCGGACATATAAGTGATCGGGTCTACGTTTGTTAGACTTGAGACGTCTGTTAAGAGTTGGGTCTAATAAATCGCCTACACGTGGAGTCTAACGGGGAAGCGTCGAATCCTGATACATCATATAATGGAGCGTGTTATGAAAAAAGAGCATTCCATTGTACGAGCCGTGCCAGAAACGGCTTGACTACGTGAGCGTAGTGTTAGATAAACAGGAAACACTGACGCGGTTAGAAGGCGGATTGCCGGTAGGTTTTGGAAACATAAATACACACGGTATCATGTTGGGTCACGATTCCTATCACCGCACAGGGCCAACCATAGAAGAACTGAAAGAACTAATCTGGCGGCGGGCTCGGTGCTTATATTTTCCACCCAACATCGTGCACATTAGGCTCACCGCGCCCTACGGGCGAAGGGTGCGTACGGTGTTTATAAGGCGTGACGGCCCCAAGTAGAGGGTAATTCTGTGAAAGAATCTCAGGACGGTGGCATGAATTCAATTCCTTTTAAACCTATCGTTCCGACCTTATGCAATCCTTCAATGAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTACTCTTCACTGCCCAAGCGAGGCTAACCTTTTGAAATGTCACAGTCGAAGCATATCTCCCGTACATCTTTTTCGGAGATCGCAGCTCGCGGAGCTATAAGCGACTTAAGCCCTTGTGTCGGTGATCCCAAGGGTCTGACTCCTGTACCAGGGTTACTGTTTCGCTTTACGGAGTAGCCTGTGAGGTGAACTGAAAGGAGCATATTTGAGATCTAAGATAGGGTCCTCCTCTGCGTCTACGTTCTCTCCGTTACGTACGGCTTCGCACCGGAGTGCATCTTGGCCCCGAAACGCACTGTGTGTGCTGATACAGCGTCCCTGGCCGGCCATGGGTTCAGAACTCCCGGGAACGCTTTTCAACTTAGAGGAACCCCGTCATGGAAGTAGATCGCGTCGAATGAGGGAGTTAGTCCTCGTTCCAGCTGGTAATTGTTTTACCGCTTGGGACCACTATAGGCCGCGGGTAGAAGTTGCTGGGTGTTGATTCCAACCCTCGAACCACGATACGACCTGCCATTTATGGCACAGTAAGGTTCAAACAGCATAATGAATACAGTTATAGTAACTTCCTCACGTACGATTAGGACGCAGCCTTG
printf '>chrA\n%s\n' "$REF" > ref.fa
"$BWA_MEM3" index --meth ref.fa > /dev/null 2>&1 || fail "index --meth nonzero exit"

# Two 300 bp fragments, methylated at the first CpG of each 100 bp end. The
# library is directional: R1 reads the converted ORIGINAL strand and R2 the
# reverse complement of its far end, so both mates carry that strand's C->T.
# On the top strand a converted CpG reads TG; on the bottom strand its C is
# the top strand's G, so in top-strand coordinates the CpG reads CA.
qual=$(printf 'I%.0s' $(seq 1 100))
add_pair() { # $1 = name, $2 = R1, $3 = R2
    printf '@%s\n%s\n+\n%s\n' "$1" "$2" "$qual" >> r1.fq
    printf '@%s\n%s\n+\n%s\n' "$1" "$3" "$qual" >> r2.fq
}
ot_left=${REF:200:100}
ot_right=${REF:400:100}
add_pair ot "${ot_left/CG/TG}" "$(dna_revcomp "${ot_right/CG/TG}")"
ob_left=${REF:600:100}
ob_right=${REF:800:100}
add_pair ob "$(dna_revcomp "${ob_right/CG/CA}")" "${ob_left/CG/CA}"

align() { # $1 = output, rest = --meth spelling
    local out="$1"
    shift
    "$BWA_MEM3" mem "$@" -t 1 ref.fa r1.fq r2.fq > "$out" 2> "$out.err" \
        || fail "mem $* nonzero exit: $(tail -3 "$out.err" | tr '\n' ' ')"
}
records() { grep -v '^@' "$1"; }
# QNAME FLAG RNAME POS CIGAR of every record, in output order, '|'-joined (no
# awk: this must hold where mawk is missing).
placements() { records "$1" | cut -f 1-4,6 | tr '\t' ' ' | tr '\n' '|'; }
# Both pairs placed mate for mate: ot at 201/401, ob at 801/601, each end
# aligned end to end (100M). Two pairs give no insert-size statistics, so
# neither is flagged proper (no 0x2 in FLAG).
PLACED="ot 97 chrA 201 100M|ot 145 chrA 401 100M|ob 81 chrA 801 100M|ob 161 chrA 601 100M|"
check_placed() { # $1 = SAM, $2 = label
    [ "$(placements "$1")" = "$PLACED" ] \
        || fail "$2: expected records '$PLACED', got '$(placements "$1")'"
}
# Count of character $2 in each record's XM:Z, as "<name>/<mate>=<n>" sorted.
cpg_calls() {
    mawk -v want="$2" '!/^@/ {
        for (i = 12; i <= NF; i++) if ($i ~ /^XM:Z:/) {
            xm = substr($i, 6)
            print $1 "/" ((int($2 / 64) % 2) ? 1 : 2) "=" gsub(want, "", xm)
        }
    }' "$1" | sort | tr '\n' ' ' | sed 's/ $//'
}
ONE_EACH="ob/1=1 ob/2=1 ot/1=1 ot/2=1"

align taps.sam --meth=taps
check_placed taps.sam --meth=taps

# --- 1. every 5-base spelling is --meth=taps, record for record ------------
for chem in 5base 5-base; do
    align "$chem.sam" --meth="$chem"
    diff <(records taps.sam) <(records "$chem.sam") > "$chem.diff" \
        || fail "--meth=$chem records differ from --meth=taps: $(head -4 "$chem.diff" | cut -c1-200 | tr '\n' ' ')"
done
# TAPS (and so 5base) defaults to neutral scoring.
align taps-neutral.sam --meth=taps --meth-scoring neutral
diff <(records taps.sam) <(records taps-neutral.sam) > taps-neutral.diff \
    || fail "--meth=taps does not default to --meth-scoring neutral: $(head -4 taps-neutral.diff | cut -c1-200 | tr '\n' ' ')"

align emseq.sam --meth
check_placed emseq.sam "--meth (em-seq control)"

# --- 1b. the last --meth wins -----------------------------------------------
align last-bare.sam --meth=5base --meth
diff <(records emseq.sam) <(records last-bare.sam) > last-bare.diff \
    || fail "--meth=5base --meth: a later bare --meth must select em-seq, but records differ from --meth: $(head -4 last-bare.diff | cut -c1-200 | tr '\n' ' ')"
align last-chem.sam --meth --meth=5base
diff <(records 5base.sam) <(records last-chem.sam) > last-chem.diff \
    || fail "--meth --meth=5base: records differ from --meth=5base: $(head -4 last-chem.diff | cut -c1-200 | tr '\n' ' ')"

# --- 2. polarity: the converted CpG is the methylated call -----------------
polarity=checked
if command -v mawk > /dev/null 2>&1; then
    [ "$(cpg_calls 5base.sam Z)" = "$ONE_EACH" ] \
        || fail "--meth=5base: want exactly one methylated CpG call (Z) per mate, got '$(cpg_calls 5base.sam Z)'"
    [ "$(cpg_calls emseq.sam z)" = "$ONE_EACH" ] \
        || fail "--meth (em-seq control): want exactly one unmethylated CpG call (z) per mate, got '$(cpg_calls emseq.sam z)'"
    # The retained CpGs: unmethylated under 5base, methylated under em-seq.
    [ -n "$(cpg_calls 5base.sam z)" ] && [ "$(cpg_calls 5base.sam z)" = "$(cpg_calls emseq.sam Z)" ] \
        || fail "retained CpGs: 5base 'z' calls '$(cpg_calls 5base.sam z)' differ from em-seq 'Z' calls '$(cpg_calls emseq.sam Z)'"
else
    polarity=SKIPPED
    echo "SKIP: polarity check only (mawk not on PATH, required to inspect SAM tags)"
fi

# --- 3. bad spellings are diagnosed -----------------------------------------
if "$BWA_MEM3" mem --meth=6base ref.fa r1.fq > /dev/null 2> err.txt; then
    fail "--meth=6base exited 0"
fi
grep -q "'5base'" err.txt || fail "unknown-chemistry error does not list 5base; got: $(head -1 err.txt)"
if "$BWA_MEM3" mem --meth=illumina ref.fa r1.fq > /dev/null 2> err.txt; then
    fail "--meth=illumina exited 0; it must not alias 5base"
fi
grep -q "got 'illumina'" err.txt || fail "--meth=illumina: not rejected as an unknown chemistry; got: $(head -1 err.txt)"
echo "PASS: meth_5base_alias (every run gives the same four placed records; 5base/5-base == taps record for record; taps defaults to neutral; the last --meth wins; 5mC->T read as methylated and retained C as unmethylated: $polarity; illumina and unknown chemistries rejected)"
