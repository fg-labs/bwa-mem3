#!/usr/bin/env bash
# test/regression/sam_bam_header_parity.sh
#
# Regression: `--bam` is a container choice, not a content change, so a --bam
# run carries the header and records the SAM text of the same run does. The two
# writers used to build their headers separately, and disagreed:
#
#   - a -H @SQ block (which, as in bwa, replaces the index's @SQ) worked in SAM
#     text, but --bam added the generated @SQ too and htslib refused the
#     duplicate, so the run failed;
#   - an @SQ block listing the contigs in another order than the index (from the
#     index's .hdr / .dict sidecar, or -H) made --bam label every record with
#     the wrong contig: its records carried the index's contig numbers, which
#     index the header's list;
#   - an @HD that was not the first line of the sidecar or -H stayed where it
#     was in SAM text, not first as the SAM spec requires (and as --bam had it).
#
# Each case runs paired-end reads through both writers -- most pairs on one
# contig, a few with their mates on different contigs so RNEXT names a contig
# too -- and requires the same header (bar the @PG command line) and the
# same records:
#
#   plain            the index's own @SQ
#   sidecar order    a .hdr sidecar listing the contigs in reverse
#   -H @PG           a supplied @PG, including a tag following its CL field
#   -H order         a -H @SQ block listing them in reverse
#   sidecar @HD      a sidecar whose @HD follows an @CO; @HD is the first line
#   -H @HD           -H text whose @HD follows an @RG; @HD is the first line
#   --compat -H @HD  the same -H under --compat=bwa-mem warns (SAM text only)
#   -H partial       a -H @SQ block missing a contig: SAM text still runs, and
#                    --bam fails naming the contig rather than mislabel records
#   -H wrong LN      a -H @SQ block with a contig's LN short of its length: SAM
#                    text still runs, and --bam fails naming the contig
#
# The header comparison ignores only @PG CL fields (the command lines). Every
# compared case also requires one bwa-mem3 @PG in each header with a nonempty CL
# and its VN equal to the binary's version.
#
# Inputs (env vars):
#   BWA_MEM3              path to the bwa-mem3 binary under test
#   HDR_PARITY_WORK_DIR   directory for the fixture and outputs
# Needs samtools (to read the --bam output) and python3.

set -euo pipefail

trap 'echo "FAIL: sam_bam_header_parity: line $LINENO: $BASH_COMMAND" >&2' ERR

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${HDR_PARITY_WORK_DIR:?HDR_PARITY_WORK_DIR must be set}"

command -v samtools > /dev/null 2>&1 || {
    echo "FAIL: samtools not on PATH; the --bam output cannot be read" >&2
    exit 1
}
command -v python3 > /dev/null 2>&1 || {
    echo "FAIL: python3 not on PATH; the fixture cannot be generated" >&2
    exit 1
}

# The script cds into the work directory, so make BWA_MEM3 absolute first: a
# bare name is looked up on PATH (whose entries may be relative), and a
# relative path is taken against the current directory.
if [[ $BWA_MEM3 != */* ]]; then
    BWA_MEM3="$(command -v "$BWA_MEM3")" || {
        echo "FAIL: sam_bam_header_parity: BWA_MEM3 is not on PATH" >&2
        exit 1
    }
fi
case "$BWA_MEM3" in
    /*) ;;
    *) BWA_MEM3="$(CDPATH='' cd -- "$(dirname -- "$BWA_MEM3")" && pwd)/$(basename -- "$BWA_MEM3")" ;;
esac

mkdir -p "$HDR_PARITY_WORK_DIR"
cd "$HDR_PARITY_WORK_DIR"
rm -f ref.fa ref.fa.* r1.fq r2.fq ./*.sam ./*.bam ./*.txt ./*.err ./*.log

# Three contigs of random sequence; 120 pairs, 40 from each contig, then 6
# cross-contig pairs (x0..x5) whose read 1 is on chrA and read 2 on chrC.
python3 - << 'PY'
import random

random.seed(653)
B = "ACGT"
ctg = {"chrA": 9000, "chrB": 7000, "chrC": 5000}
seq = {n: "".join(random.choice(B) for _ in range(l)) for n, l in ctg.items()}
comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
with open("ref.fa", "w") as fh:
    for n, s in seq.items():
        fh.write(">%s\n%s\n" % (n, "\n".join(s[i:i + 60] for i in range(0, len(s), 60))))
with open("r1.fq", "w") as f1, open("r2.fq", "w") as f2:
    k = 0
    for n, s in seq.items():
        for _ in range(40):
            st = random.randrange(0, len(s) - 400)
            r1 = s[st:st + 100]
            r2 = "".join(comp[b] for b in reversed(s[st + 300:st + 400]))
            f1.write("@p%d\n%s\n+\n%s\n" % (k, r1, "I" * 100))
            f2.write("@p%d\n%s\n+\n%s\n" % (k, r2, "I" * 100))
            k += 1
    for k in range(6):
        a = random.randrange(0, len(seq["chrA"]) - 100)
        c = random.randrange(0, len(seq["chrC"]) - 100)
        r1 = seq["chrA"][a:a + 100]
        r2 = "".join(comp[b] for b in reversed(seq["chrC"][c:c + 100]))
        f1.write("@x%d\n%s\n+\n%s\n" % (k, r1, "I" * 100))
        f2.write("@x%d\n%s\n+\n%s\n" % (k, r2, "I" * 100))
PY

"$BWA_MEM3" index ref.fa > index.log 2>&1

# The version every @PG VN must carry: the first line of `bwa-mem3 version`.
"$BWA_MEM3" version > version.txt 2>&1
bin_version=$(sed -n 1p version.txt)
[[ -n $bin_version ]] || {
    echo "FAIL: sam_bam_header_parity: could not read the binary's version" >&2
    exit 1
}

fails=0
checks=0
fail() {
    echo "FAIL: $*" >&2
    fails=$((fails + 1))
}

# run_both <case> <mem args...>: SAM text to <case>.sam, --bam to <case>.bam.
run_both() {
    local c="$1"
    shift
    "$BWA_MEM3" mem -t 1 "$@" ref.fa r1.fq r2.fq > "$c.sam" 2> "$c.sam.err" \
        || {
            fail "$c: mem (SAM) exited non-zero: $(tail -1 "$c.sam.err")"
            return 1
        }
    "$BWA_MEM3" mem -t 1 --bam=0 "$@" ref.fa r1.fq r2.fq > "$c.bam" 2> "$c.bam.err" \
        || {
            fail "$c: mem --bam exited non-zero: $(tail -1 "$c.bam.err")"
            return 1
        }
}

# same <case>: compare the full headers (ignoring only @PG CL fields) and
# alignment records, and validate the generated @PG identity.
same() {
    local c="$1"
    checks=$((checks + 1))
    grep '^@' "$c.sam" | sed $'/^@PG\t/s/\tCL:[^\t]*//' > "$c.sam.hdr.txt"
    samtools view -H --no-PG "$c.bam" | sed $'/^@PG\t/s/\tCL:[^\t]*//' > "$c.bam.hdr.txt"
    grep -v '^@' "$c.sam" > "$c.sam.rec.txt"
    samtools view "$c.bam" > "$c.bam.rec.txt"
    # Compare every header record, ignoring only @PG CL fields. Also require
    # bwa-mem3's own @PG in each container, with a command line and the binary's
    # version, so a shared omission cannot pass the parity check.
    grep $'^@PG\tID:bwa-mem3\t' "$c.sam" > "$c.sam.pg.txt" || true
    samtools view -H --no-PG "$c.bam" | grep $'^@PG\tID:bwa-mem3\t' > "$c.bam.pg.txt" || true
    local sam_pg bam_pg sam_pg_nocl bam_pg_nocl
    sam_pg=$(wc -l < "$c.sam.pg.txt" | tr -d ' ')
    bam_pg=$(wc -l < "$c.bam.pg.txt" | tr -d ' ')
    sam_pg_nocl=$(sed $'s/\tCL:[^\t]*//' "$c.sam.pg.txt")
    bam_pg_nocl=$(sed $'s/\tCL:[^\t]*//' "$c.bam.pg.txt")
    if [[ $sam_pg != 1 || $bam_pg != 1 ]]; then
        fail "$c: expected one bwa-mem3 @PG line in each header, got $sam_pg (SAM) and $bam_pg (--bam)"
    elif ! grep -Eq $'\tCL:[^\t]+' "$c.sam.pg.txt" || ! grep -Eq $'\tCL:[^\t]+' "$c.bam.pg.txt"; then
        fail "$c: bwa-mem3 @PG line has no nonempty CL (command line) field"
    elif [[ $sam_pg_nocl != "$bam_pg_nocl" ]]; then
        fail "$c: the @PG lines differ outside CL: '$sam_pg_nocl' (SAM) vs '$bam_pg_nocl' (--bam)"
    elif [[ $(tr '\t' '\n' <<< "$sam_pg_nocl" | sed -n 's/^VN://p') != "$bin_version" ]]; then
        fail "$c: @PG VN is not the binary's version '$bin_version': '$sam_pg_nocl'"
    elif ! cmp -s "$c.sam.hdr.txt" "$c.bam.hdr.txt"; then
        fail "$c: --bam header differs from the SAM text's"
        diff "$c.sam.hdr.txt" "$c.bam.hdr.txt" | head -10 >&2 || true
    elif ! cmp -s "$c.sam.rec.txt" "$c.bam.rec.txt"; then
        fail "$c: --bam records differ from the SAM text's"
        diff "$c.sam.rec.txt" "$c.bam.rec.txt" | head -6 >&2 || true
    else
        echo "ok: $c: --bam header and $(wc -l < "$c.sam.rec.txt" | tr -d ' ') records equal the SAM text's"
    fi
}

# expect <label> <want> <got>
expect() {
    checks=$((checks + 1))
    if [[ $2 == "$3" ]]; then
        echo "ok: $1"
    else
        fail "$1: expected '$2', got '$3'"
    fi
}

# Every same-contig pair's RNAME (and its mate's, RNEXT "=") must be its own
# contig, and each cross-contig pair must name the other mate's contig in RNEXT.
contigs_right() {
    local c="$1"
    expect "$c: every pair on its own contig" "chrA 80 chrB 80 chrC 80" \
        "$(awk -F '\t' '!/^@/ && $7 == "=" { n[$3]++ } END { printf "chrA %d chrB %d chrC %d", n["chrA"], n["chrB"], n["chrC"] }' "$c.sam")"
    expect "$c: cross-contig mates name each other's contig" "chrA>chrC 6 chrC>chrA 6" \
        "$(awk -F '\t' '!/^@/ && $1 ~ /^x/ && $2 < 256 { n[$3 ">" $7]++ } END { printf "chrA>chrC %d chrC>chrA %d", n["chrA>chrC"], n["chrC>chrA"] }' "$c.sam")"
}

# sq_order <case>: one check that the SAM header lists the contigs in reverse.
sq_order() {
    expect "$1: @SQ order" "SN:chrC SN:chrB SN:chrA" \
        "$(grep '^@SQ' "$1.sam" | cut -f2 | paste -sd' ' -)"
}

rev_sq=$'@SQ\tSN:chrC\tLN:5000\n@SQ\tSN:chrB\tLN:7000\n@SQ\tSN:chrA\tLN:9000'

if run_both plain; then
    same plain
    contigs_right plain
    expect "plain: exactly one @HD" 1 "$(grep -c $'^@HD\t' plain.sam || true)"
    expect "plain: default @HD first" $'@HD\tVN:1.5\tSO:unsorted\tGO:query' "$(head -1 plain.sam)"
fi

printf '%s\n' "$rev_sq" > ref.fa.hdr
if run_both sidecar_order; then
    same sidecar_order
    contigs_right sidecar_order
    sq_order sidecar_order
fi
rm -f ref.fa.hdr

# CL deliberately precedes VN so normalization must retain following tags.
user_pg=$'@PG\tID:prep\tPN:prep\tCL:prep reads\tVN:1.0'
if run_both h_pg -H "$user_pg"; then
    same h_pg
    contigs_right h_pg
    expect "-H @PG: supplied record preserved" "$user_pg" "$(grep $'^@PG\tID:prep\t' h_pg.sam)"
fi

if run_both h_order -H "$rev_sq"; then
    same h_order
    contigs_right h_order
    sq_order h_order
fi

printf '@CO\tfrom the sidecar\n@HD\tVN:1.6\tSO:unsorted\n%s\n' "$rev_sq" > ref.fa.hdr
if run_both sidecar_hd; then
    same sidecar_hd
    contigs_right sidecar_hd
    expect "sidecar @HD: the first line" $'@HD\tVN:1.6\tSO:unsorted' "$(head -1 sidecar_hd.sam)"
fi
rm -f ref.fa.hdr

if run_both h_hd -H $'@RG\tID:rg1\tSM:s1\n@HD\tVN:1.6\tSO:unsorted'; then
    same h_hd
    contigs_right h_hd
    expect "-H @HD: the first line" $'@HD\tVN:1.6\tSO:unsorted' "$(head -1 h_hd.sam)"
fi

# Under --compat the moved @HD makes the header differ from the target's, so
# an @HD that was not the first line of -H warns too.
checks=$((checks + 1))
if ! "$BWA_MEM3" mem -t 1 --compat=bwa-mem -H $'@RG\tID:rg1\tSM:s1\n@HD\tVN:1.6\tSO:unsorted' \
    ref.fa r1.fq r2.fq > compat_hd.sam 2> compat_hd.sam.err; then
    fail "--compat -H @HD: mem exited non-zero: $(tail -1 compat_hd.sam.err)"
elif ! grep -q 'with an @HD from -H' compat_hd.sam.err; then
    fail "--compat -H @HD: a non-leading -H @HD did not warn"
else
    expect "--compat -H @HD: the first line" $'@HD\tVN:1.6\tSO:unsorted' "$(head -1 compat_hd.sam)"
    echo "ok: --compat -H @HD: a non-leading -H @HD warns"
fi

checks=$((checks + 1))
partial=$'@SQ\tSN:chrA\tLN:9000\n@SQ\tSN:chrB\tLN:7000'
if ! "$BWA_MEM3" mem -t 1 -H "$partial" ref.fa r1.fq r2.fq > h_partial.sam 2> h_partial.sam.err; then
    fail "-H partial: SAM text should still run"
elif [[ $(grep $'^@SQ\t' h_partial.sam) != "$partial" ]]; then
    fail "-H partial: SAM @SQ block did not preserve the supplied records"
elif [[ $(grep -cv '^@' h_partial.sam || true) != 252 ]]; then
    fail "-H partial: SAM text did not preserve all 252 alignment records"
elif "$BWA_MEM3" mem -t 1 --bam=0 -H "$partial" ref.fa r1.fq r2.fq > h_partial.bam 2> h_partial.bam.err; then
    fail "-H partial: --bam should fail when an index contig has no @SQ record"
elif ! grep -q "contig 'chrC' of the index has no @SQ record" h_partial.bam.err; then
    fail "-H partial: --bam failed without naming the missing contig: $(head -1 h_partial.bam.err)"
else
    echo "ok: -H partial: SAM text runs; --bam fails naming chrC"
fi

checks=$((checks + 1))
short_ln=$'@SQ\tSN:chrA\tLN:9000\n@SQ\tSN:chrB\tLN:7000\n@SQ\tSN:chrC\tLN:4000'
if ! "$BWA_MEM3" mem -t 1 -H "$short_ln" ref.fa r1.fq r2.fq > h_short_ln.sam 2> h_short_ln.sam.err; then
    fail "-H wrong LN: SAM text should still run"
elif [[ $(grep $'^@SQ\t' h_short_ln.sam) != "$short_ln" ]]; then
    fail "-H wrong LN: SAM @SQ block did not preserve the supplied records"
elif [[ $(grep -cv '^@' h_short_ln.sam || true) != 252 ]]; then
    fail "-H wrong LN: SAM text did not preserve all 252 alignment records"
elif "$BWA_MEM3" mem -t 1 --bam=0 -H "$short_ln" ref.fa r1.fq r2.fq > h_short_ln.bam 2> h_short_ln.bam.err; then
    fail "-H wrong LN: --bam should fail when an @SQ LN is not the index contig's length"
elif ! grep -q "contig 'chrC' is 5000 bp in the index but its @SQ record in the output header has LN:4000" h_short_ln.bam.err; then
    fail "-H wrong LN: --bam failed without naming the contig and both lengths: $(head -1 h_short_ln.bam.err)"
else
    echo "ok: -H wrong LN: SAM text runs; --bam fails naming chrC and both lengths"
fi

if [[ $fails -gt 0 ]]; then
    echo "FAIL: sam_bam_header_parity: $fails of $checks checks failed"
    exit 1
fi
echo "PASS: sam_bam_header_parity: $checks checks"
