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
# Each case runs paired-end reads (so mate contigs are covered too) through
# both writers and requires the same header (bar the @PG command line) and the
# same records:
#
#   plain            the index's own @SQ
#   sidecar order    a .hdr sidecar listing the contigs in reverse
#   -H order         a -H @SQ block listing them in reverse
#   sidecar @HD      a sidecar whose @HD follows an @CO; @HD is the first line
#   -H @HD           -H text whose @HD follows an @RG; @HD is the first line
#   -H partial       a -H @SQ block missing a contig: SAM text still runs, and
#                    --bam fails naming the contig rather than mislabel records
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

mkdir -p "$HDR_PARITY_WORK_DIR"
cd "$HDR_PARITY_WORK_DIR"
rm -f ref.fa ref.fa.* r1.fq r2.fq ./*.sam ./*.bam ./*.txt ./*.err ./*.log

# Three contigs of random sequence; 120 pairs, 40 from each contig.
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
PY

"$BWA_MEM3" index ref.fa > index.log 2>&1

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

# same <case>: one check that the SAM and --bam outputs carry the same header
# (the @PG line, whose CL records the differing command line, aside) and the
# same records.
same() {
    local c="$1"
    checks=$((checks + 1))
    grep '^@' "$c.sam" | grep -v '^@PG' > "$c.sam.hdr.txt"
    samtools view -H --no-PG "$c.bam" | grep -v '^@PG' > "$c.bam.hdr.txt"
    grep -v '^@' "$c.sam" > "$c.sam.rec.txt"
    samtools view "$c.bam" > "$c.bam.rec.txt"
    if ! cmp -s "$c.sam.hdr.txt" "$c.bam.hdr.txt"; then
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

# Every read's RNAME (and its mate's, RNEXT "=") must be its own contig.
contigs_right() {
    local c="$1"
    expect "$c: every pair on its own contig" "chrA 80 chrB 80 chrC 80" \
        "$(awk -F '\t' '!/^@/ && $7 == "=" { n[$3]++ } END { printf "chrA %d chrB %d chrC %d", n["chrA"], n["chrB"], n["chrC"] }' "$c.sam")"
}

rev_sq=$'@SQ\tSN:chrC\tLN:5000\n@SQ\tSN:chrB\tLN:7000\n@SQ\tSN:chrA\tLN:9000'

if run_both plain; then same plain; fi

printf '%s\n' "$rev_sq" > ref.fa.hdr
if run_both sidecar_order; then
    same sidecar_order
    contigs_right sidecar_order
fi
rm -f ref.fa.hdr

if run_both h_order -H "$rev_sq"; then
    same h_order
    contigs_right h_order
fi

printf '@CO\tfrom the sidecar\n@HD\tVN:1.6\tSO:unsorted\n%s\n' "$rev_sq" > ref.fa.hdr
if run_both sidecar_hd; then
    same sidecar_hd
    expect "sidecar @HD: the first line" $'@HD\tVN:1.6\tSO:unsorted' "$(head -1 sidecar_hd.sam)"
fi
rm -f ref.fa.hdr

if run_both h_hd -H $'@RG\tID:rg1\tSM:s1\n@HD\tVN:1.6\tSO:unsorted'; then
    same h_hd
    expect "-H @HD: the first line" $'@HD\tVN:1.6\tSO:unsorted' "$(head -1 h_hd.sam)"
fi

checks=$((checks + 1))
partial=$'@SQ\tSN:chrA\tLN:9000\n@SQ\tSN:chrB\tLN:7000'
if ! "$BWA_MEM3" mem -t 1 -H "$partial" ref.fa r1.fq r2.fq > h_partial.sam 2> h_partial.sam.err; then
    fail "-H partial: SAM text should still run"
elif "$BWA_MEM3" mem -t 1 --bam=0 -H "$partial" ref.fa r1.fq r2.fq > h_partial.bam 2> h_partial.bam.err; then
    fail "-H partial: --bam should fail when an index contig has no @SQ record"
elif ! grep -q "contig 'chrC' of the index has no @SQ record" h_partial.bam.err; then
    fail "-H partial: --bam failed without naming the missing contig: $(head -1 h_partial.bam.err)"
else
    echo "ok: -H partial: SAM text runs; --bam fails naming chrC"
fi

if [[ $fails -gt 0 ]]; then
    echo "FAIL: sam_bam_header_parity: $fails of $checks checks failed"
    exit 1
fi
echo "PASS: sam_bam_header_parity: $checks checks"
