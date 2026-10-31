#!/usr/bin/env bash
# test/regression/aln_parity.sh
#
# Regression: `bwa-mem3 aln` writes the same SAM as stock bwa 0.7.19's `aln` followed by
# `samse` (single-end) or `sampe` (paired-end) -- every record, byte for byte; only the @PG line differs (it names
# the program that ran) and is dropped from both sides before the diff.
#
# The reference is phiX plus the synthetic contigs the unit tests use (an N run, a
# reverse-complement copy, a tandem repeat), indexed once by stock bwa and once by
# bwa-mem3, and a second copy of both indices with an .alt file naming one contig. Reads
# come from test/regression/aln_reads.awk (exact copies, substitutions, 1-bp indels,
# ambiguous bases, reads from the repeat and the reverse-complement copy, both strands),
# plus hand-made reads: lower case, IUPAC codes R and Y, a '-', an all-N read.
#
# (A read with a '-' makes stock write a NUL byte into the SAM -- its base table maps '-' to
# 5 and "ACGTN"[5] is the string's terminator -- so every grep and diff here is binary-safe.)
#
# The paired-end section does the same with `aln` on each file and `sampe` (pairs come from
# test/regression/aln_pairs.awk), against `bwa-mem3 aln <opts> [sampe options] ref r1.fq r2.fq`.
#
# Each case runs with -t 1 and -t 4. Stock is `aln <opts> ref r.fq > r.sai` then
# `samse [-n K] [-r RG] ref r.sai r.fq`; ours is `aln <opts> [--max-xa-hits K] [-r RG] ref r.fq`.
#
# Inputs (env vars):
#   BWA_MEM3          -- path to the bwa-mem3 binary under test
#   BWA_ORACLE        -- path to a stock bwa 0.7.19 binary
#                        (test/tools/build_bwa_oracle.sh); unset => SKIP
#   FIXTURES          -- directory holding phix.fa (default: test/fixtures)
#   ALN_SAI_DUMP      -- path to aln_sai_dump, used only to write the synthetic contigs
#                        (default: repo root)

set -Eeuo pipefail
# shellcheck disable=SC2154  # rc is assigned by `rc=$?` at the start of the trap body
trap 'rc=$?; echo "FAIL: unexpected command failure at line $LINENO (exit $rc)" >&2; exit "$rc"' ERR

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
dump="${ALN_SAI_DUMP:-$repo_root/aln_sai_dump}"

work="$(mktemp -d "${TMPDIR:-/tmp}/aln_parity.XXXXXX")"
cleanup() { rm -rf "$work"; }
trap cleanup EXIT

# Reference: phiX followed by the synthetic contigs; "plain" has no .alt, "alt" marks c2.
"$dump" --write-fasta "$work/synthetic.fa"
cat "$fixtures/phix.fa" "$work/synthetic.fa" > "$work/ref.fa"
for d in plain alt; do
    for tool in stock mem3; do
        mkdir -p "$work/$d/$tool"
        cp "$work/ref.fa" "$work/$d/$tool/ref.fa"
        if [ "$d" = alt ]; then printf 'c2\n' > "$work/$d/$tool/ref.fa.alt"; fi
    done
    "$BWA_ORACLE" index "$work/$d/stock/ref.fa" > "$work/$d.stock.index.log" 2>&1
    "$BWA_MEM3" index "$work/$d/mem3/ref.fa" > "$work/$d.mem3.index.log" 2>&1
done

gen_reads() {
    awk -v seed="$1" -v n="$2" -v lens="$3" -v tail="$4" -f "$repo_root/test/regression/aln_reads.awk" "$work/ref.fa"
}

# Rewrites every eighth read: lower case, an R, a Y, a '-', an N (the rest are untouched).
mutate_reads() {
    awk 'NR % 4 == 2 {
        i = (NR - 2) / 4
        k = i % 8
        if (k == 1) $0 = tolower($0)
        else if (k == 2) $0 = substr($0, 1, 9) "R" substr($0, 11)
        else if (k == 3) $0 = substr($0, 1, 11) "Y" substr($0, 13)
        else if (k == 4) $0 = substr($0, 1, 4) "-" substr($0, 6)
        else if (k == 5) $0 = substr($0, 1, 6) "N" substr($0, 8)
    } { print }'
}

n_cfg=0
mapped_total=0
xa_total=0
ah_total=0

# compare <tag> <refdir> <fq> <xa> <rg> [aln options...]
# xa and rg are the samse -n / --max-xa-hits value and the -r line; "" leaves them at the default.
compare() {
    local tag="$1" d="$2" fq="$3" xa="$4" rg="$5"
    shift 5
    local t
    for t in 1 4; do
        n_cfg=$((n_cfg + 1))
        local samse_args=() ours_args=()
        if [ -n "$xa" ]; then
            samse_args+=(-n "$xa")
            ours_args+=(--max-xa-hits "$xa")
        fi
        if [ -n "$rg" ]; then
            samse_args+=(-r "$rg")
            ours_args+=(-r "$rg")
        fi
        "$BWA_ORACLE" aln "$@" -t "$t" "$work/$d/stock/ref.fa" "$fq" > "$work/stock.sai" 2> "$work/stock.aln.log"
        "$BWA_ORACLE" samse ${samse_args[@]+"${samse_args[@]}"} "$work/$d/stock/ref.fa" "$work/stock.sai" "$fq" \
            > "$work/stock.sam" 2> "$work/stock.samse.log"
        "$BWA_MEM3" aln "$@" -t "$t" ${ours_args[@]+"${ours_args[@]}"} "$work/$d/mem3/ref.fa" "$fq" \
            > "$work/mem3.sam" 2> "$work/mem3.aln.log"
        grep -av '^@PG' "$work/stock.sam" > "$work/stock.nopg" || true
        grep -av '^@PG' "$work/mem3.sam" > "$work/mem3.nopg" || true
        if ! diff -a "$work/stock.nopg" "$work/mem3.nopg" > "$work/sam.diff"; then
            echo "FAIL: aln_parity: $tag -t $t $* differs from stock samse; first lines of diff (< stock, > ours):" >&2
            head -n 10 "$work/sam.diff" >&2
            exit 1
        fi
        if ! grep -q '^@PG' "$work/mem3.sam"; then
            echo "FAIL: aln_parity: $tag -t $t has no @PG line" >&2
            exit 1
        fi
        mapped_total=$((mapped_total + $({ grep -av '^@' "$work/mem3.nopg" || true; } | mawk_or_awk '$3 != "*"' | wc -l)))
        xa_total=$((xa_total + $(grep -ac 'XA:Z:' "$work/mem3.nopg" || true)))
        ah_total=$((ah_total + $(grep -ac 'AH:\*' "$work/mem3.nopg" || true)))
        echo "aln_parity: ok $tag -t $t $* ($(grep -avc '^@' "$work/mem3.nopg" || true) records)"
    done
}

mawk_or_awk() { if command -v mawk > /dev/null 2>&1; then mawk "$@"; else awk "$@"; fi; }

# need_tag <tag> <pattern>: the case's last output must contain the pattern, so a case that
# silently stopped exercising its feature fails instead of passing on empty output.
need_tag() {
    if ! grep -aq -- "$2" "$work/mem3.nopg"; then
        echo "FAIL: aln_parity: $1 produced no '$2' in its output" >&2
        exit 1
    fi
}

# run_case <name> <seed> <n> <lens> <tail> <xa> <rg> [aln options...]
run_case() {
    local name="$1" seed="$2" n="$3" lens="$4" tail="$5" xa="$6" rg="$7"
    shift 7
    gen_reads "$seed" "$n" "$lens" "$tail" | mutate_reads_if "$name" > "$work/$name.fq"
    compare "$name" plain "$work/$name.fq" "$xa" "$rg" "$@"
}
# Only the ambiguity cases rewrite bases; the others keep the generator's reads as they are.
mutate_reads_if() { if [ "$1" = ambig ]; then mutate_reads; else cat; fi; }

run_case default36 201 2000 36 0 "" ""
need_tag default36 'XA:Z:'
run_case adna 202 2000 20-40 0 "" "" -l 1024 -n 0.04 -o 2
run_case indel100 203 1000 100 0 "" "" -n 0.04 -o 2
need_tag indel100 'MD:Z:[0-9]*\^'
run_case ambig 204 1500 36 0 "" ""
# Reads that overlap the 50-base N run of c2 by 1-4 bases at either end, on both strands: the
# alignments cover reference N bases, so stock prints an XN tag.
awk '
function rc(s,   i, r, c) {
    r = ""
    for (i = length(s); i >= 1; i--) {
        c = substr(s, i, 1)
        r = r (c == "A" ? "T" : c == "C" ? "G" : c == "G" ? "C" : c == "T" ? "A" : c)
    }
    return r
}
/^>/ { in_c2 = ($1 == ">c2"); next }
in_c2 { c2 = c2 $0 }
END {
    n = 0
    for (k = 1; k <= 4; k++) {
        for (side = 0; side < 2; side++) {
            t = side ? substr(c2, 351 - k, 36) : substr(c2, 300 - 35 + k, 36)
            if (n % 2) t = rc(t)
            q = ""
            for (i = 1; i <= length(t); i++) q = q "I"
            printf "@nrun%d\n%s\n+\n%s\n", n++, t, q
        }
    }
}' "$work/synthetic.fa" > "$work/nrun.fq"
compare nrun plain "$work/nrun.fq" "" "" -n 4
need_tag nrun 'XN:i:'
run_case qtrim 205 500 100 1 "" "" -q 20
need_tag qtrim 'XC:i:'
run_case xa0 206 1500 36 0 0 ""
run_case xa10 207 1500 24-36 0 10 "" -n 2
run_case readgroup 208 800 36 0 "" '@RG\tID:x\tSM:y'
need_tag readgroup 'RG:Z:x'
run_case nonstop 209 600 36 0 "" "" -N

# One 1,000 bp read and twenty 250 bp reads in a single batch (the batch-wide gap-open cap
# follows the longest read), plus an all-N read and a read with a '-'.
{
    gen_reads 301 1 1000 0
    gen_reads 302 20 250 0
    printf '@alln\n%s\n+\n%s\n' "$(printf 'N%.0s' $(seq 1 36))" "$(printf 'I%.0s' $(seq 1 36))"
} > "$work/long.fq"
compare long plain "$work/long.fq" "" "" -n 0.04 -o 2

# An empty FASTQ gives a header and nothing else.
: > "$work/empty.fq"
compare empty plain "$work/empty.fq" "" ""
if grep -aqv '^@' "$work/mem3.nopg"; then
    echo "FAIL: aln_parity: the empty FASTQ produced records" >&2
    exit 1
fi

# A reference with an .alt file: @SQ carries AH:* for the alt contig.
gen_reads 401 1000 36 0 > "$work/alt.fq"
compare alt alt "$work/alt.fq" "" ""
need_tag alt 'AH:\*'

if [ "$mapped_total" -le 0 ] || [ "$xa_total" -le 0 ] || [ "$ah_total" -le 0 ]; then
    echo "FAIL: aln_parity: the cases found no mapped reads, XA tags or AH:* markers" >&2
    exit 1
fi
echo "PASS: aln_parity SE ($((n_cfg / 2)) cases x 2 thread counts)"

# ---- paired end: aln R1, aln R2, sampe ------------------------------------------------
# Reads come from aln_pairs.awk (fragments of 250-349 bp, so most pairs are concordant; a few
# have a far second mate or a random one; "rescue" gives some second mates too many scattered
# mismatches for the search, so only the mate Smith-Waterman can place them).
pe_cfg=0
pe_xt_m=0
pe_ref=plain

gen_pairs() { # seed n lens rescue names out1 out2 [lens2 heavy1]
    awk -v seed="$1" -v n="$2" -v lens="$3" -v rescue="$4" -v names="$5" -v out1="$6" -v out2="$7" \
        -v lens2="${8:-}" -v heavy1="${9:-}" \
        -f "$repo_root/test/regression/aln_pairs.awk" "$work/ref.fa"
}

# compare_pe <tag> <fq1> <fq2> <sampe args> <ours args> [aln options...]
# The two argument strings are split on spaces (-r values with spaces are not used here).
compare_pe() {
    local tag="$1" fq1="$2" fq2="$3" sampe_s="$4" ours_s="$5"
    shift 5
    local t sampe_args=() ours_args=()
    read -r -a sampe_args <<< "$sampe_s" || true
    read -r -a ours_args <<< "$ours_s" || true
    for t in 1 4; do
        pe_cfg=$((pe_cfg + 1))
        "$BWA_ORACLE" aln "$@" -t "$t" "$work/$pe_ref/stock/ref.fa" "$fq1" > "$work/stock.1.sai" 2> "$work/stock.aln1.log"
        "$BWA_ORACLE" aln "$@" -t "$t" "$work/$pe_ref/stock/ref.fa" "$fq2" > "$work/stock.2.sai" 2> "$work/stock.aln2.log"
        "$BWA_ORACLE" sampe ${sampe_args[@]+"${sampe_args[@]}"} "$work/$pe_ref/stock/ref.fa" "$work/stock.1.sai" "$work/stock.2.sai" "$fq1" "$fq2" \
            > "$work/stock.sam" 2> "$work/stock.sampe.log"
        "$BWA_MEM3" aln "$@" -t "$t" ${ours_args[@]+"${ours_args[@]}"} "$work/$pe_ref/mem3/ref.fa" "$fq1" "$fq2" \
            > "$work/mem3.sam" 2> "$work/mem3.aln.log"
        grep -av '^@PG' "$work/stock.sam" > "$work/stock.nopg" || true
        grep -av '^@PG' "$work/mem3.sam" > "$work/mem3.nopg" || true
        if ! diff -a "$work/stock.nopg" "$work/mem3.nopg" > "$work/sam.diff"; then
            echo "FAIL: aln_parity: PE $tag -t $t $* differs from stock sampe; first lines of diff (< stock, > ours):" >&2
            head -n 10 "$work/sam.diff" >&2
            exit 1
        fi
        if ! grep -q '^@PG' "$work/mem3.sam"; then
            echo "FAIL: aln_parity: PE $tag -t $t has no @PG line" >&2
            exit 1
        fi
        pe_xt_m=$((pe_xt_m + $(grep -ac 'XT:A:M' "$work/mem3.nopg" || true)))
        echo "aln_parity: ok PE $tag -t $t $* ($(grep -avc '^@' "$work/mem3.nopg" || true) records)"
    done
}

# run_pe <name> <seed> <n> <lens> <rescue> <names> <sampe args> <ours args> [aln options...]
run_pe() {
    local name="$1" seed="$2" n="$3" lens="$4" rescue="$5" names="$6" sampe_s="$7" ours_s="$8"
    shift 8
    gen_pairs "$seed" "$n" "$lens" "$rescue" "$names" "$work/$name.1.fq" "$work/$name.2.fq"
    compare_pe "$name" "$work/$name.1.fq" "$work/$name.2.fq" "$sampe_s" "$ours_s" "$@"
}

# need_pe <tag> <pattern>: as need_tag, on the last PE case.
need_pe() { need_tag "PE $1" "$2"; }

run_pe pe36 501 2000 36 0 plain "" ""
need_pe pe36 'XT:A:U'
run_pe pe150 502 800 150 0 plain "" ""
run_pe pe_adna 503 1500 20-40 0 plain "" "" -l 1024 -n 0.04 -o 2
# Files of different read lengths (R1 36 bp, R2 150 bp) under fractional -n with a large -o and a
# small -m: each file is searched as its own batch, as stock's separate aln runs do. (The
# batch-wide gap-open cap turned out not to change any output -- see the aln README -- so this
# case guards the mixed-length path, it does not tell a per-file search from a combined one.)
gen_pairs 521 1500 36 0 plain "$work/pe_mixlen.1.fq" "$work/pe_mixlen.2.fq" 150 30
compare_pe pe_mixlen "$work/pe_mixlen.1.fq" "$work/pe_mixlen.2.fq" "" "" -n 0.04 -o 5 -m 100
run_pe pe_rescue 504 1500 100 30 plain "" ""
need_pe pe_rescue 'XT:A:M'
run_pe pe_rescue150 505 600 150 30 plain "" ""
run_pe pe_maxocc 506 1500 36 0 plain "-o 3" "--max-occ 3"
run_pe pe_xa 507 1500 36 0 plain "-n 5 -N 20" "--max-xa-hits 5 --max-xa-hits-discordant 20"
run_pe pe_a300 508 1500 36 0 plain "-a 300" "-a 300"
run_pe pe_a300_A 520 1500 36 0 plain "-A -a 300" "-A -a 300"
run_pe pe_s 509 1500 100 30 plain "-s" "-s"
run_pe pe_A 510 1500 100 30 plain "-A" "-A"
run_pe pe_c 511 1500 36 0 plain "-c 0.1" "-c 0.1"
run_pe pe_slash 512 800 36 0 slash "" ""
run_pe pe_comment 513 800 36 0 comment "" ""
run_pe pe_rg 514 500 36 0 plain "-r @RG\\tID:x\\tSM:y" "-r @RG\\tID:x\\tSM:y"
need_pe pe_rg 'RG:Z:x'

# Repeats: contigs of five and of twelve exact copies of a 200 bp unit, 40 bp apart, so each end
# has five (twelve) equally good places and XA lists depend on --max-xa-hits and
# --max-xa-hits-discordant. Reads come from those contigs only; the reference is phiX plus them.
awk 'BEGIN {
    srand(7)
    split("5 12", ncopy, " ")
    for (m = 1; m <= 2; m++) {
        unit = ""
        for (i = 0; i < 200; i++) unit = unit substr("ACGT", int(rand() * 4) + 1, 1)
        s = ""
        for (c = 0; c < ncopy[m]; c++) {
            s = s unit
            for (i = 0; i < 40; i++) s = s substr("ACGT", int(rand() * 4) + 1, 1)
        }
        printf ">m%d\n", ncopy[m]
        for (i = 1; i <= length(s); i += 60) print substr(s, i, 60)
    }
}' > "$work/multi.fa"
for tool in stock mem3; do
    mkdir -p "$work/multi/$tool"
    cat "$fixtures/phix.fa" "$work/multi.fa" > "$work/multi/$tool/ref.fa"
done
"$BWA_ORACLE" index "$work/multi/stock/ref.fa" > "$work/multi.stock.index.log" 2>&1
"$BWA_MEM3" index "$work/multi/mem3/ref.fa" > "$work/multi.mem3.index.log" 2>&1
pe_ref=multi
for spec in "multi_default||" "multi_xa5|-n 5|--max-xa-hits 5" "multi_xa20|-N 20|--max-xa-hits-discordant 20" "multi_occ|-o 3|--max-occ 3"; do
    IFS='|' read -r name sampe_s ours_s <<< "$spec"
    awk -v seed=517 -v n=800 -v lens=36 -v rescue=0 -v uniform=1 -v names=plain -v out1="$work/$name.1.fq" -v out2="$work/$name.2.fq" \
        -f "$repo_root/test/regression/aln_pairs.awk" "$work/multi.fa"
    compare_pe "$name" "$work/$name.1.fq" "$work/$name.2.fq" "$sampe_s" "$ours_s"
    if [ "$name" = multi_xa5 ]; then need_pe multi_xa5 'XA:Z:'; fi
done
pe_ref=plain

# R1 holds three reads more than R2: sampe pairs the first R2-count reads and stops; we do
# the same and say so on stderr.
gen_pairs 515 600 36 0 plain "$work/pe_extra.1.fq" "$work/pe_extra.2.fq"
gen_pairs 516 3 36 0 plain "$work/pe_extra_tail.1.fq" "$work/pe_extra_tail.2.fq"
sed 's/^@r/@x/' "$work/pe_extra_tail.1.fq" >> "$work/pe_extra.1.fq"
compare_pe pe_extra "$work/pe_extra.1.fq" "$work/pe_extra.2.fq" "" ""
if ! grep -q 'read counts differ' "$work/mem3.aln.log"; then
    echo "FAIL: aln_parity: PE pe_extra printed no warning about the unequal read counts" >&2
    exit 1
fi

if [ "$pe_xt_m" -le 0 ]; then
    echo "FAIL: aln_parity: no PE case produced a mate-rescued record (XT:A:M)" >&2
    exit 1
fi
echo "PASS: aln_parity PE ($((pe_cfg / 2)) cases x 2 thread counts, $pe_xt_m rescued mates)"
