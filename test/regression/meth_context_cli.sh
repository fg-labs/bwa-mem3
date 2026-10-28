#!/usr/bin/env bash
# test/regression/meth_context_cli.sh
#
# Regression: the --meth-context command line. Pins the contract that the
# option parser, the cross-option validation and the stray-value diagnosis keep,
# and that each accepted spelling and scoring reaches the alignment (by AS).
#
# Asserted here:
#  1. Each of the four rejections prints its exact message and exits non-zero:
#     an unknown value, cpg under em-seq/bisulfite, cpg with collapsed scoring,
#     and --meth-context without --meth.
#  2. `--meth=5base --meth-context all` is record-for-record identical to
#     `--meth=5base` alone (all is the default), and the `--meth-context=cpg`
#     spelling with '=' gives the separated spelling's records, on a read that
#     cpg and all score differently.
#  3. cpg is accepted with taps, 5base and 5-base, and the read still maps
#     forward at its locus; an explicit --meth-scoring genomic is honoured
#     under cpg (a freed CpG conversion scores +A rather than neutral's 0).
#  4. A stray `cpg` (`--meth=5base cpg ref r1`) is diagnosed by the stray-value
#     check, naming --meth-context, rather than reported as a missing index.
#
# Inputs:
#   BWA_MEM3 — the bwa-mem3 binary under test: a bare command name resolved via
#              PATH, or a path to the binary (a relative one is made absolute
#              before the `cd "$WORK"` below)
set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"
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

# One 100 bp read with a converted CpG (CG -> TG), as a 5-base library gives.
read=${REF:200:100}
printf '@r\n%s\n+\n%s\n' "${read/CG/TG}" "$(printf 'I%.0s' $(seq 1 100))" > r.fq
# One 100 bp read with a non-CpG conversion mid-read (CA -> TA), which cpg scores as a
# mismatch and all frees, so the two contexts give it different records.
read2=${REF:600:100}
tail2=${read2:30}
case "$tail2" in *CA*) ;; *) fail "no CA after offset 30 in read2" ;; esac
head2=${tail2%%CA*}
nc_pos=$((30 + ${#head2})) # the first C of a CA at or after offset 30
read2="${read2:0:$nc_pos}T${read2:$((nc_pos + 1))}"
printf '@r2\n%s\n+\n%s\n' "$read2" "$(printf 'I%.0s' $(seq 1 100))" > r2.fq

records() { grep -v '^@' "$1"; }

# --- 1. the four rejections, exact copy ------------------------------------
expect_error() { # $1 = label, $2 = exact message, rest = mem argv
    local label="$1" want="$2"
    shift 2
    if "$BWA_MEM3" mem "$@" -t 1 ref.fa r.fq > /dev/null 2> err.txt; then
        fail "$label: exited 0"
    fi
    grep -Fxq -- "$want" err.txt \
        || fail "$label: message not found verbatim; want '$want', got '$(head -3 err.txt | tr '\n' ' ')'"
}
expect_error "unknown value" \
    "ERROR: --meth-context accepts 'all' (default) or 'cpg', got 'chg'" \
    --meth=5base --meth-context chg
expect_error "cpg + emseq" \
    "ERROR: --meth-context cpg requires --meth=taps or --meth=5base: under EM-seq/bisulfite every unmethylated C converts, so conversions are not confined to CpG" \
    --meth --meth-context cpg
expect_error "cpg + bisulfite" \
    "ERROR: --meth-context cpg requires --meth=taps or --meth=5base: under EM-seq/bisulfite every unmethylated C converts, so conversions are not confined to CpG" \
    --meth=bisulfite --meth-context cpg
expect_error "cpg + collapsed" \
    "ERROR: --meth-context cpg is incompatible with --meth-scoring collapsed: collapsed also frees reference T against read C, which has no CpG meaning" \
    --meth=5base --meth-context cpg --meth-scoring collapsed
expect_error "no --meth" \
    "ERROR: --meth-context requires --meth" \
    --meth-context all

# --- 2. all is the default --------------------------------------------------
"$BWA_MEM3" mem --meth=5base -t 1 ref.fa r.fq > plain.sam 2> plain.err || fail "--meth=5base nonzero exit"
"$BWA_MEM3" mem --meth=5base --meth-context all -t 1 ref.fa r.fq > all.sam 2> all.err \
    || fail "--meth-context all nonzero exit"
[ "$(records plain.sam | wc -l | tr -d ' ')" = "1" ] || fail "--meth=5base: expected 1 record"
diff <(records plain.sam) <(records all.sam) > all.diff \
    || fail "--meth-context all differs from the default: $(head -4 all.diff | cut -c1-200 | tr '\n' ' ')"

# --- 3. cpg is accepted wherever it is valid --------------------------------
for chem in taps 5base 5-base; do
    "$BWA_MEM3" mem --meth="$chem" --meth-context cpg -t 1 ref.fa r.fq > "cpg-$chem.sam" 2> "cpg-$chem.err" \
        || fail "--meth=$chem --meth-context cpg rejected: $(head -2 "cpg-$chem.err" | tr '\n' ' ')"
    [ "$(records "cpg-$chem.sam" | wc -l | tr -d ' ')" = "1" ] || fail "--meth=$chem --meth-context cpg: expected 1 record"
    # The read is chrA:201-300 with one converted CpG: it must still map, there.
    loc="$(records "cpg-$chem.sam" | cut -f 2-4 | tr '\t' ' ')"
    [ "$loc" = "0 chrA 201" ] || fail "--meth=$chem --meth-context cpg: want FLAG 0 at chrA:201, got '$loc'"
done
as_of() { records "$1" | grep -o 'AS:i:-\{0,1\}[0-9]*' | cut -d: -f3; }
"$BWA_MEM3" mem --meth=taps --meth-context=cpg -t 1 ref.fa r2.fq > eq.sam 2> eq.err \
    || fail "--meth-context=cpg (with '=') rejected: $(head -2 eq.err | tr '\n' ' ')"
"$BWA_MEM3" mem --meth=taps --meth-context cpg -t 1 ref.fa r2.fq > sep.sam 2> sep.err \
    || fail "--meth-context cpg (separated) on r2 rejected: $(head -2 sep.err | tr '\n' ' ')"
"$BWA_MEM3" mem --meth=taps --meth-context all -t 1 ref.fa r2.fq > all2.sam 2> all2.err \
    || fail "--meth-context all on r2 rejected: $(head -2 all2.err | tr '\n' ' ')"
diff <(records sep.sam) <(records eq.sam) > eq.diff \
    || fail "--meth-context=cpg differs from --meth-context cpg: $(head -4 eq.diff | cut -c1-200 | tr '\n' ' ')"
# The comparison only means something if r2 tells the contexts apart: cpg charges
# its non-CpG conversion as a mismatch, so its AS is below all's.
as_cpg=$(as_of sep.sam)
as_all=$(as_of all2.sam)
[ -n "$as_cpg" ] && [ -n "$as_all" ] && [ "$as_cpg" -lt "$as_all" ] \
    || fail "r2 does not separate cpg from all: AS cpg='$as_cpg' all='$as_all'"
"$BWA_MEM3" mem --meth=taps --meth-context cpg --meth-scoring genomic -t 1 ref.fa r.fq > g.sam 2> g.err \
    || fail "cpg with explicit genomic scoring rejected: $(head -2 g.err | tr '\n' ' ')"
# r's one conversion is at a CpG, freed under cpg: genomic scores it +A (1), neutral
# (the taps default, cpg-taps.sam above) scores it 0, so AS is exactly one higher.
as_g=$(as_of g.sam)
as_n=$(as_of cpg-taps.sam)
[ -n "$as_g" ] && [ -n "$as_n" ] && [ "$as_g" -eq $((as_n + 1)) ] \
    || fail "cpg + genomic: want AS one above neutral's, got genomic='$as_g' neutral='$as_n'"

# --- 4. a stray cpg is diagnosed ---------------------------------------------
if "$BWA_MEM3" mem --meth=5base cpg ref.fa r.fq > /dev/null 2> stray.txt; then
    fail "a stray 'cpg' exited 0"
fi
grep -q "'cpg' was taken as a positional argument" stray.txt \
    || fail "stray cpg: error does not name the token; got: $(head -2 stray.txt | tr '\n' ' ')"
grep -Fq -- "the value for --meth-context and no such file exists" stray.txt \
    || fail "stray cpg: error does not attribute it to --meth-context; got: $(head -2 stray.txt | tr '\n' ' ')"
grep -Fq -- "--meth-context cpg." stray.txt \
    || fail "stray cpg: advice does not name '--meth-context cpg'; got: $(tail -2 stray.txt | tr '\n' ' ')"

echo "PASS: meth_context_cli (unknown value, cpg+emseq, cpg+collapsed and no --meth rejected with exact copy; all == default; --meth-context=cpg == separated; cpg accepted with taps/5base/5-base; genomic honoured under cpg; stray cpg diagnosed)"
