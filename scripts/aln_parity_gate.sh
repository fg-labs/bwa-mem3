#!/usr/bin/env bash
# scripts/aln_parity_gate.sh
#
# The large parity gate for `bwa-mem3 aln`: it must write the same SAM as stock bwa 0.7.19's
# `aln` followed by `samse` (single end) or `sampe` (paired end), every record byte for byte;
# only the @PG line differs (it names the program that ran) and is dropped from both sides.
#
# This is too large for CI (hundreds of thousands of reads per case, several minutes), so it
# lives here and NOT under test/regression/: the CI coverage lint does not apply to it. The
# small, CI-sized version of the same check is test/regression/aln_parity.sh. Run this by hand
# on each host (Linux glibc and macOS) before a release that touches `aln`.
#
# Reads are simulated with holodeck (a seeded run, so the same inputs give the same reads);
# nothing is committed. Stock is run once per case at THREADS threads (its output does not
# depend on the thread count); bwa-mem3 runs at every thread count in the case's list, by
# default 1 and THREADS.
#
# Cases (select with CASES, default all that have their inputs):
#   1   SE and PE, 700,000 reads / pairs, spanning three 262,144-read batches
#   2   PE where the middle batch is random sequence absent from the reference, so stock's
#       infer_isize fails there and the third batch reuses the previous batch's estimate
#   3   mixed 20-150 bp reads, -n 0.04 -o 2, one batch (SE and PE)
#   4   20-40 bp reads, -l 1024 -n 0.04 -o 2 (SE and PE)
#   5   2 x 150 bp PE, default options
#   6   read names ending /1, /2, /3 and names with comments (SE and PE)
#   7   the option sweep: each aln option alone, then the samse and sampe options
#   8   a reference with an .alt file (ALT_REF), SE and PE
#   9   real ancient-DNA reads (ADNA_DIR) against a whole-genome index (REF_STOCK and REF_MEM3),
#       at default options and at -l 1024 -n 0.04 -o 2; SE, plus PE on the merged-library pairs
#   10  PE where R2 ends before R1 at, just before, just after and well across a batch boundary
#
# Inputs (environment variables):
#   BWA_ORACLE   stock bwa 0.7.19 binary (test/tools/build_bwa_oracle.sh)            required
#   BWA_MEM3     the bwa-mem3 binary under test                                      required
#   REF          reference FASTA the simulated reads come from (cases 1-7, 10)       required
#   WORK         work directory (created; per-case files are deleted after a PASS)  required
#   THREADS      thread count for the multi-threaded runs (default: all CPUs)
#   REF_STOCK    stock index prefix for REF (default: index REF into WORK with stock bwa)
#   REF_MEM3     bwa-mem3 index prefix for REF (default: index REF into WORK with bwa-mem3)
#   ALT_REF      FASTA of a reference with alt contigs, with ALT_REF.alt next to it (case 8)
#   ADNA_DIR     directory holding the case 9 reads (file names are in the case 9 section)
#   ADNA_THREADS thread counts for case 9 (default "4 $THREADS": -t 1 is too slow there)
#   CASES        space-separated case numbers (default: 1 2 3 4 5 6 7 8 9 10)
#   SEED         base seed for the simulated reads (default 20261030)
#
# Requires holodeck, samtools, gzip and mawk or awk on PATH. Exits 0 only if no case FAILed
# (a case with missing inputs is reported SKIP and does not fail the run, so check the
# summary). On a FAIL the first differing record from each side is printed.

set -Eeuo pipefail

: "${BWA_ORACLE:?BWA_ORACLE must be set}"
: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${WORK:?WORK must be set}"
SEED="${SEED:-20261030}"
CASES="${CASES:-1 2 3 4 5 6 7 8 9 10}"
if [ -z "${THREADS:-}" ]; then
    THREADS="$(nproc 2> /dev/null || sysctl -n hw.ncpu)"
fi
ADNA_THREADS="${ADNA_THREADS:-4 $THREADS}"
for tool in holodeck samtools gzip; do
    command -v "$tool" > /dev/null 2>&1 || {
        echo "FAIL: $tool not found on PATH" >&2
        exit 2
    }
done

mkdir -p "$WORK/run" "$WORK/reads"
RUN="$WORK/run"
RESULTS="$WORK/results.tsv"
: > "$RESULTS"
T_LIST="1 $THREADS"
if [ "$THREADS" = 1 ]; then T_LIST=1; fi
n_fail=0

mawk_or_awk() { if command -v mawk > /dev/null 2>&1; then mawk "$@"; else awk "$@"; fi; }
say() { printf '%s\n' "$*"; }

# result <PASS|FAIL|SKIP> <case tag> <detail>
result() {
    printf '%s\t%s\t%s\n' "$1" "$2" "$3" >> "$RESULTS"
    say "$1 $2: $3"
    if [ "$1" = FAIL ]; then n_fail=$((n_fail + 1)); fi
}

# ---- references -------------------------------------------------------------------------
# prepare_ref <name> <fasta> [stock prefix] [mem3 prefix]
# Sets, for later use, REF_<name>_FA (the FASTA with its .fai), _STOCK and _MEM3 (index prefixes).
prepare_ref() {
    local name="$1" fa="$2" stock="${3:-}" mem3="${4:-}" d="$WORK/ref_$1"
    mkdir -p "$d"
    if [ ! -e "$d/ref.fa" ]; then ln -s "$fa" "$d/ref.fa"; fi
    if [ ! -s "$d/ref.fa.fai" ]; then samtools faidx "$d/ref.fa"; fi
    if [ -z "$stock" ]; then
        mkdir -p "$d/stock" "$d/mem3"
        for tool in stock mem3; do
            cp "$fa" "$d/$tool/ref.fa"
            if [ -f "$fa.alt" ]; then cp "$fa.alt" "$d/$tool/ref.fa.alt"; fi
        done
        if [ ! -s "$d/stock/ref.fa.bwt" ]; then "$BWA_ORACLE" index "$d/stock/ref.fa" > "$d/stock.index.log" 2>&1; fi
        if [ ! -s "$d/mem3/ref.fa.bwt.2bit.64" ]; then "$BWA_MEM3" index "$d/mem3/ref.fa" > "$d/mem3.index.log" 2>&1; fi
        stock="$d/stock/ref.fa"
        mem3="$d/mem3/ref.fa"
    fi
    printf -v "REF_${name}_FA" '%s' "$d/ref.fa"
    printf -v "REF_${name}_STOCK" '%s' "$stock"
    printf -v "REF_${name}_MEM3" '%s' "$mem3"
}

# use_ref <name>: select the reference the next comparisons run against.
use_ref() {
    local s="REF_$1_STOCK" m="REF_$1_MEM3" f="REF_$1_FA"
    RS="${!s}"
    RM="${!m}"
    RFA="${!f}"
}

# ---- read generation ----------------------------------------------------------------------
# gen <out prefix> <n reads or pairs> <read length> <pe 0|1> <seed offset> [holodeck options]
# Writes <prefix>.fq (SE) or <prefix>.1.fq and <prefix>.2.fq (PE), exactly n records each.
gen() {
    local prefix="$1" n="$2" len="$3" pe="$4" seed="$5"
    shift 5
    local glen cov d="$RUN/holo.$$" i got
    glen="$(mawk_or_awk '{s += $2} END {print s}' "$RFA.fai")"
    # generous: a part of the reference is N, and a few pairs are rejected and resampled
    cov="$(mawk_or_awk -v n="$n" -v len="$len" -v g="$glen" -v pe="$pe" 'BEGIN { printf "%.4f", 1.6 * n * len * (pe ? 2 : 1) / g }')"
    rm -rf "$d"
    mkdir -p "$d"
    local pe_args=(--single-end)
    if [ "$pe" = 1 ]; then pe_args=(); fi
    holodeck simulate -r "$RFA" -o "$d/h" -c "$cov" -l "$len" --seed "$((SEED + seed))" \
        ${pe_args[@]+"${pe_args[@]}"} "$@" > "$d/holo.log" 2>&1
    if [ "$pe" = 1 ]; then
        for i in 1 2; do
            gzip -dc "$d/h.r$i.fastq.gz" | mawk_or_awk -v lim="$((4 * n))" 'NR <= lim' > "$prefix.$i.fq"
        done
    else
        local f="$d/h.fastq.gz"
        if [ ! -e "$f" ]; then f="$d/h.r1.fastq.gz"; fi
        gzip -dc "$f" | mawk_or_awk -v lim="$((4 * n))" 'NR <= lim' > "$prefix.fq"
    fi
    rm -rf "$d"
    for got in "$prefix.fq" "$prefix.1.fq" "$prefix.2.fq"; do
        if [ -e "$got" ] && [ "$(($(wc -l < "$got") / 4))" -ne "$n" ]; then
            say "FAIL: holodeck gave too few reads for $got (wanted $n)" >&2
            exit 2
        fi
    done
}

# Filters (stdin to stdout) over four-line FASTQ records.
# trim_lengths <lo> <hi> <seed>: cut every read to a random length in [lo, hi] (if longer).
trim_lengths() {
    mawk_or_awk -v lo="$1" -v hi="$2" -v seed="$3" 'BEGIN { srand(seed) }
        NR % 4 == 2 { l = lo + int(rand() * (hi - lo + 1)); if (l < length($0)) $0 = substr($0, 1, l); n = length($0) }
        NR % 4 == 0 { $0 = substr($0, 1, n) }
        { print }'
}
# randomize_range <first record> <last record> <seed>: replace the bases of those records (1-based)
# with random sequence.
randomize_range() {
    mawk_or_awk -v lo="$1" -v hi="$2" -v seed="$3" 'BEGIN { srand(seed) }
        NR % 4 == 2 { r = int((NR - 2) / 4) + 1
            if (r >= lo && r <= hi) { s = ""; for (i = 0; i < length($0); i++) s = s substr("ACGT", int(rand() * 4) + 1, 1); $0 = s } }
        { print }'
}
# rename_se: names ending /1, /2, /3, or carrying a comment, in turn.
rename_se() {
    mawk_or_awk 'NR % 4 == 1 { i = int((NR - 1) / 4) % 4
        if (i == 0) $0 = $0 "/1"; else if (i == 1) $0 = $0 "/2"; else if (i == 2) $0 = $0 "/3"
        else $0 = $0 " 1:N:0:ACGT extra comment" } { print }'
}
# rename_pe <mate 1|2>: pairs named base/1 + base/2, base/3 + base/3, base + comment, base/1 + comment.
rename_pe() {
    mawk_or_awk -v m="$1" 'NR % 4 == 1 { i = int((NR - 1) / 4) % 4
        if (i == 0) $0 = $0 "/" m; else if (i == 1) $0 = $0 "/3"
        else if (i == 2) $0 = $0 " " m ":N:0:ACGT"; else $0 = $0 "/" m " mate comment" } { print }'
}
# phred64: shift the qualities from phred+33 to phred+64.
phred64() {
    mawk_or_awk 'BEGIN { for (i = 33; i < 127; i++) ord[sprintf("%c", i)] = i }
        NR % 4 == 0 { s = ""; for (i = 1; i <= length($0); i++) s = s sprintf("%c", ord[substr($0, i, 1)] + 31); $0 = s }
        { print }'
}
# casava: add a Casava 1.8 comment to every read, filtered (Y) on every other one.
casava() { mawk_or_awk 'NR % 4 == 1 { $0 = $0 " 1:" ((int((NR - 1) / 4) % 2) ? "Y" : "N") ":0:ACGT" } { print }'; }

# ---- comparison -----------------------------------------------------------------------------
n_records() { grep -avc '^@' "$1" || true; }

# same_sam <stock sam> <ours sam>: 0 if equal apart from @PG; else print the first differing
# record from each side and return 1.
same_sam() {
    if cmp -s <(grep -av '^@PG' "$1" || true) <(grep -av '^@PG' "$2" || true); then return 0; fi
    local line
    line="$(cmp <(grep -av '^@PG' "$1" || true) <(grep -av '^@PG' "$2" || true) 2>&1 | mawk_or_awk '{ for (i = 1; i <= NF; i++) if ($i == "line") { gsub(",", "", $(i + 1)); print $(i + 1); exit } }' || true)"
    say "  first difference at line ${line:-?} of the SAM without @PG lines"
    say "  stock: $(grep -av '^@PG' "$1" | sed -n "${line:-1}p" | cut -c1-600 || true)"
    say "  ours:  $(grep -av '^@PG' "$2" | sed -n "${line:-1}p" | cut -c1-600 || true)"
    return 1
}

# run_tool <log> <command...>: run quietly; on failure print the log tail and return 1.
run_tool() {
    local log="$1"
    shift
    if ! "$@" 2> "$log"; then
        say "  command failed: $*"
        tail -n 5 "$log" | sed 's/^/  | /'
        return 1
    fi
}

# split_words <string>: split on spaces into the global array WORDS (backslashes kept).
split_words() {
    WORDS=()
    read -r -a WORDS <<< "$1" || true
}

# ours_vs_stock <SE|PE> <tag> <stock seconds> <ours args string> <n inputs> <inputs...> [aln options...]
# Runs bwa-mem3 at every thread count in T_LIST on the inputs, compares each run with
# $RUN/stock.sam (already written by the caller) and reports the case. The part compare_se and
# compare_pe share; they differ only in the stock pipeline and the number of inputs.
ours_vs_stock() {
    local kind="$1" tag="$2" w_stock="$3" ours_s="$4" n_in="$5" t t0 w_ours ok=1 detail
    shift 5
    local inputs=("${@:1:$n_in}")
    shift "$n_in"
    local ours_args=()
    split_words "$ours_s"
    ours_args=(${WORDS[@]+"${WORDS[@]}"})
    detail="stock ${w_stock}s"
    for t in $T_LIST; do
        t0=$SECONDS
        if ! run_tool "$RUN/ours.err" "$BWA_MEM3" aln "$@" -t "$t" ${ours_args[@]+"${ours_args[@]}"} "$RM" "${inputs[@]}" > "$RUN/ours.sam"; then
            ok=0
            detail="$detail; -t $t: bwa-mem3 failed"
            continue
        fi
        w_ours=$((SECONDS - t0))
        if same_sam "$RUN/stock.sam" "$RUN/ours.sam" > "$RUN/diff.txt"; then
            detail="$detail; -t $t ${w_ours}s"
        else
            ok=0
            detail="$detail; -t $t DIFFERS"
            cat "$RUN/diff.txt"
        fi
    done
    if [ "$ok" = 1 ]; then
        result PASS "$kind $tag $*" "$(n_records "$RUN/ours.sam") records; $detail"
    else
        result FAIL "$kind $tag $*" "$detail"
    fi
}

# compare_se <tag> <fq> <samse args> <ours args> [aln options...]
compare_se() {
    local tag="$1" fq="$2" samse_s="$3" ours_s="$4" t0
    shift 4
    local samse_args=()
    split_words "$samse_s"
    samse_args=(${WORDS[@]+"${WORDS[@]}"})
    t0=$SECONDS
    if ! { run_tool "$RUN/stock.aln.err" "$BWA_ORACLE" aln "$@" -t "$THREADS" "$RS" "$fq" > "$RUN/stock.sai" \
        && run_tool "$RUN/stock.samse.err" "$BWA_ORACLE" samse ${samse_args[@]+"${samse_args[@]}"} "$RS" "$RUN/stock.sai" "$fq" > "$RUN/stock.sam"; }; then
        result FAIL "SE $tag" "stock failed"
        return
    fi
    ours_vs_stock SE "$tag" $((SECONDS - t0)) "$ours_s" 1 "$fq" "$@"
    rm -f "$RUN"/stock.sai "$RUN"/stock.sam "$RUN"/ours.sam
}

# compare_pe <tag> <fq1> <fq2> <sampe args> <ours args> [aln options...]
# Prints, into $RUN/stock.sampe.err, stock's messages (cases look for infer_isize failures).
compare_pe() {
    local tag="$1" fq1="$2" fq2="$3" sampe_s="$4" ours_s="$5" t0
    shift 5
    local sampe_args=()
    split_words "$sampe_s"
    sampe_args=(${WORDS[@]+"${WORDS[@]}"})
    t0=$SECONDS
    if ! { run_tool "$RUN/stock.aln1.err" "$BWA_ORACLE" aln "$@" -t "$THREADS" "$RS" "$fq1" > "$RUN/stock.1.sai" \
        && run_tool "$RUN/stock.aln2.err" "$BWA_ORACLE" aln "$@" -t "$THREADS" "$RS" "$fq2" > "$RUN/stock.2.sai" \
        && run_tool "$RUN/stock.sampe.err" "$BWA_ORACLE" sampe ${sampe_args[@]+"${sampe_args[@]}"} "$RS" "$RUN/stock.1.sai" "$RUN/stock.2.sai" "$fq1" "$fq2" > "$RUN/stock.sam"; }; then
        result FAIL "PE $tag" "stock failed"
        return
    fi
    ours_vs_stock PE "$tag" $((SECONDS - t0)) "$ours_s" 2 "$fq1" "$fq2" "$@"
    rm -f "$RUN"/stock.1.sai "$RUN"/stock.2.sai "$RUN"/stock.sam "$RUN"/ours.sam
}

wanted() { case " $CASES " in *" $1 "*) return 0 ;; *) return 1 ;; esac }

# ---- the cases ------------------------------------------------------------------------------
if [ -n "${REF:-}" ]; then
    prepare_ref main "$REF" "${REF_STOCK:-}" "${REF_MEM3:-}"
    use_ref main
fi
need_ref() {
    if [ -z "${REF:-}" ]; then
        result SKIP "case $1" "REF is not set"
        return 1
    fi
}

BATCH=262144
# begin_case <n>: a fresh read directory for the case; end_case removes it.
begin_case() {
    R="$WORK/reads/c$1"
    rm -rf "$R"
    mkdir -p "$R"
}
end_case() { rm -rf "$R"; }

if wanted 1 && need_ref 1; then
    begin_case 1
    say "== case 1: three batches (700,000 reads / pairs)"
    gen "$R/c1" 700000 100 0 1 --max-error-rate 0.03
    compare_se c1 "$R/c1.fq" "" ""
    gen "$R/c1pe" 700000 100 1 2 --max-error-rate 0.03
    compare_pe c1 "$R/c1pe.1.fq" "$R/c1pe.2.fq" "" ""
    end_case
fi

if wanted 2 && need_ref 2; then
    begin_case 2
    say "== case 2: PE, middle batch random (infer_isize fails there; batch 3 reuses the estimate)"
    gen "$R/c2" 700000 100 1 3 --max-error-rate 0.03
    randomize_range "$((BATCH + 1))" "$((2 * BATCH))" 11 < "$R/c2.1.fq" > "$R/c2r.1.fq"
    randomize_range "$((BATCH + 1))" "$((2 * BATCH))" 12 < "$R/c2.2.fq" > "$R/c2r.2.fq"
    compare_pe c2 "$R/c2r.1.fq" "$R/c2r.2.fq" "" ""
    n_isize="$(grep -c 'fail to infer insert size' "$RUN/stock.sampe.err" || true)"
    if [ "$n_isize" -lt 1 ]; then
        result FAIL "PE c2 setup" "stock did not hit an infer_isize failure, so the carry-over is not exercised"
    else
        result PASS "PE c2 setup" "stock hit $n_isize infer_isize failure(s)"
    fi
    end_case
fi

if wanted 3 && need_ref 3; then
    begin_case 3
    say "== case 3: mixed 20-150 bp, -n 0.04 -o 2, one batch"
    gen "$R/c3" 200000 150 0 4 --max-error-rate 0.03
    trim_lengths 20 150 31 < "$R/c3.fq" > "$R/c3m.fq"
    compare_se c3 "$R/c3m.fq" "" "" -n 0.04 -o 2
    gen "$R/c3pe" 100000 150 1 5 --max-error-rate 0.03
    trim_lengths 20 150 32 < "$R/c3pe.1.fq" > "$R/c3m.1.fq"
    trim_lengths 20 150 33 < "$R/c3pe.2.fq" > "$R/c3m.2.fq"
    compare_pe c3 "$R/c3m.1.fq" "$R/c3m.2.fq" "" "" -n 0.04 -o 2
    end_case
fi

if wanted 4 && need_ref 4; then
    begin_case 4
    say "== case 4: 20-40 bp, -l 1024 -n 0.04 -o 2"
    gen "$R/c4" 300000 40 0 6 --max-error-rate 0.03
    trim_lengths 20 40 41 < "$R/c4.fq" > "$R/c4m.fq"
    compare_se c4 "$R/c4m.fq" "" "" -l 1024 -n 0.04 -o 2
    gen "$R/c4pe" 200000 40 1 7 --max-error-rate 0.03
    trim_lengths 20 40 42 < "$R/c4pe.1.fq" > "$R/c4m.1.fq"
    trim_lengths 20 40 43 < "$R/c4pe.2.fq" > "$R/c4m.2.fq"
    compare_pe c4 "$R/c4m.1.fq" "$R/c4m.2.fq" "" "" -l 1024 -n 0.04 -o 2
    end_case
fi

if wanted 5 && need_ref 5; then
    begin_case 5
    say "== case 5: 2 x 150 bp PE, default options"
    gen "$R/c5" 300000 150 1 8 --max-error-rate 0.03
    compare_pe c5 "$R/c5.1.fq" "$R/c5.2.fq" "" ""
    end_case
fi

if wanted 6 && need_ref 6; then
    begin_case 6
    say "== case 6: names ending /1 /2 /3, names with comments"
    gen "$R/c6" 100000 100 0 9
    rename_se < "$R/c6.fq" > "$R/c6n.fq"
    compare_se c6 "$R/c6n.fq" "" ""
    gen "$R/c6pe" 100000 100 1 10
    rename_pe 1 < "$R/c6pe.1.fq" > "$R/c6n.1.fq"
    rename_pe 2 < "$R/c6pe.2.fq" > "$R/c6n.2.fq"
    compare_pe c6 "$R/c6n.1.fq" "$R/c6n.2.fq" "" ""
    end_case
fi

if wanted 7 && need_ref 7; then
    begin_case 7
    say "== case 7: option sweep"
    gen "$R/c7" 150000 100 0 12 --max-error-rate 0.04
    # name|samse args|ours args|aln options
    sweep_se=(
        'n3|||-n 3' 'n001|||-n 0.01' 'o0|||-o 0' 'e2|||-e 2' 'i10|||-i 10' 'd5|||-d 5'
        'l20k1|||-l 20 -k 1' 'q20|||-q 20' 'B4|||-B 4' 'L|||-L -n 0.04 -o 2' 'R10|||-R 10'
        'MOE|||-M 2 -O 9 -E 3 -n 0.04 -o 2' 'N|||-N' 'xa0|-n 0|--max-xa-hits 0|'
        'xa10|-n 10|--max-xa-hits 10|-n 2' 'rg|-r @RG\tID:x\tSM:y|-r @RG\tID:x\tSM:y|'
    )
    for spec in "${sweep_se[@]}"; do
        IFS='|' read -r name samse_s ours_s opts <<< "$spec"
        split_words "$opts"
        compare_se "sweep_$name" "$R/c7.fq" "$samse_s" "$ours_s" ${WORDS[@]+"${WORDS[@]}"}
    done
    trim_lengths 20 100 51 < "$R/c7.fq" > "$R/c7m.fq"
    compare_se sweep_m50 "$R/c7m.fq" "" "" -m 50
    phred64 < "$R/c7.fq" > "$R/c7i.fq"
    compare_se sweep_I "$R/c7i.fq" "" "" -I -q 20
    casava < "$R/c7.fq" > "$R/c7y.fq"
    compare_se sweep_Y "$R/c7y.fq" "" "" -Y
    gen "$R/c7pe" 80000 100 1 13 --max-error-rate 0.04
    # name|sampe args|ours args
    sweep_pe=(
        'default||' 'maxocc|-o 3|--max-occ 3' 'xa|-n 5 -N 20|--max-xa-hits 5 --max-xa-hits-discordant 20'
        'a300|-a 300|-a 300' 'a300A|-A -a 300|-A -a 300' 's|-s|-s' 'A|-A|-A' 'c|-c 0.1|-c 0.1'
        'rg|-r @RG\tID:x\tSM:y|-r @RG\tID:x\tSM:y'
    )
    for spec in "${sweep_pe[@]}"; do
        IFS='|' read -r name sampe_s ours_s <<< "$spec"
        compare_pe "sweep_$name" "$R/c7pe.1.fq" "$R/c7pe.2.fq" "$sampe_s" "$ours_s"
    done
    compare_pe sweep_adna "$R/c7pe.1.fq" "$R/c7pe.2.fq" "" "" -l 1024 -n 0.04 -o 2
    end_case
fi

if wanted 8; then
    if [ -z "${ALT_REF:-}" ] || [ ! -f "${ALT_REF}.alt" ]; then
        result SKIP "case 8" "ALT_REF (with ALT_REF.alt) is not set"
    else
        begin_case 8
        say "== case 8: reference with alt contigs"
        prepare_ref alt "$ALT_REF"
        use_ref alt
        gen "$R/c8" 200000 100 0 14 --max-error-rate 0.03
        compare_se c8 "$R/c8.fq" "" ""
        head -n 400 "$R/c8.fq" > "$R/c8h.fq"
        "$BWA_MEM3" aln -t 1 "$RM" "$R/c8h.fq" > "$R/c8h.sam" 2> /dev/null
        if ! grep -aq 'AH:\*' "$R/c8h.sam"; then
            result FAIL "SE c8 setup" "no AH:* marker in the header (the .alt file was not read)"
        fi
        compare_se c8x "$R/c8.fq" "-n 100" "--max-xa-hits 100" -n 2
        gen "$R/c8pe" 100000 100 1 15 --max-error-rate 0.03
        compare_pe c8 "$R/c8pe.1.fq" "$R/c8pe.2.fq" "" ""
        end_case
        if [ -n "${REF:-}" ]; then use_ref main; fi
    fi
fi

if wanted 9; then
    if [ -z "${ADNA_DIR:-}" ] || [ -z "${REF_STOCK:-}" ] || [ -z "${REF_MEM3:-}" ]; then
        result SKIP "case 9" "ADNA_DIR, REF_STOCK and REF_MEM3 are not all set"
    else
        begin_case 9
        say "== case 9: real ancient-DNA reads (threads: $ADNA_THREADS)"
        RS="$REF_STOCK"
        RM="$REF_MEM3"
        T_LIST="$ADNA_THREADS"
        adna=("-" "-l 1024 -n 0.04 -o 2")
        for lib in ERR5343086.1M ERR5343087.1M ERR13475328.1M SRR957738.1M; do
            for o in "${adna[@]}"; do
                if [ "$o" = "-" ]; then o=""; fi
                split_words "$o"
                compare_se "adna_$lib" "$ADNA_DIR/$lib.fq" "" "" ${WORDS[@]+"${WORDS[@]}"}
            done
        done
        for o in "${adna[@]}"; do
            if [ "$o" = "-" ]; then o=""; fi
            split_words "$o"
            compare_pe adna_ERR13475328.pe "$ADNA_DIR/ERR13475328.pe.1M_1.fq" "$ADNA_DIR/ERR13475328.pe.1M_2.fq" "" "" ${WORDS[@]+"${WORDS[@]}"}
        done
        compare_pe adna_ERR668531 "$ADNA_DIR/ERR668531.1M_1.fq" "$ADNA_DIR/ERR668531.1M_2.fq" "" ""
        T_LIST="1 $THREADS"
        if [ "$THREADS" = 1 ]; then T_LIST=1; fi
        if [ -n "${REF:-}" ]; then use_ref main; fi
    fi
fi

if wanted 10 && need_ref 10; then
    begin_case 10
    say "== case 10: R2 ends before R1 at a batch boundary"
    gen "$R/c10" "$((BATCH + 8000))" 100 1 16 --max-error-rate 0.03
    for n2 in "$((BATCH - 1))" "$BATCH" "$((BATCH + 1))" "$((BATCH + 5000))"; do
        mawk_or_awk -v lim="$((4 * n2))" 'NR <= lim' "$R/c10.2.fq" > "$R/c10s.2.fq"
        compare_pe "c10_r2_${n2}" "$R/c10.1.fq" "$R/c10s.2.fq" "" ""
        if ! grep -Eq 'read 2 file ended before read 1|read counts differ' "$RUN/ours.err" 2> /dev/null; then
            result FAIL "PE c10_r2_${n2} warning" "no warning about the shorter R2 on stderr"
        fi
    done
    end_case
fi

# ---- summary --------------------------------------------------------------------------------
say ""
say "== summary ($(grep -c '^PASS' "$RESULTS" || true) PASS, $n_fail FAIL, $(grep -c '^SKIP' "$RESULTS" || true) SKIP)"
cat "$RESULTS"
if [ "$n_fail" -gt 0 ]; then
    say "FAIL: aln_parity_gate: $n_fail case(s) differ from stock"
    exit 1
fi
say "PASS: aln_parity_gate"
