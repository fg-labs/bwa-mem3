#!/usr/bin/env bash
# test/regression/meth_rescue_strand.sh
#
# Regression: a --meth mate brought in by SW mate rescue must carry the original
# strand of its own fragment (XG:Z) whichever strand its anchor mapped to.
#
# The rescue scores the mate against a window on one strand half of the doubled
# reference, so the matrix it picks is relative to that half. That value is right
# for the SW but is not the genome-strand hypothesis the rest of the aligner
# stores (mem_reg2aln flips the hypothesis itself for reverse regions).
# Recording it unchanged mislabelled every mate rescued off a reverse-strand
# anchor: an OT R1 forward came out XG:Z:GA and an OB R2 forward XG:Z:CT, and
# the final CIGAR/NM/XM were regenerated under the wrong conversion, so every
# C->T (or G->A) became a mismatch. Near the end of the last contig the window
# can also be clamped onto the half opposite the anchor's, so the strand must
# come from where the mate landed, not from the anchor.
#
# Four fragments, one per (fragment strand, anchor strand) pair, plus one pair at
# the end of the contig whose window lands on the other half. In each, the
# anchor maps cleanly and the mate carries spaced non-bisulfite errors that
# leave no seed of minimum length, so it can only be placed by mate rescue (a
# leg run with -S proves that: every such mate is unmapped there).
#
#   OT fragment, anchor R2 reverse, rescued R1 forward   -> XG:Z:CT
#   OB fragment, anchor R1 reverse, rescued R2 forward   -> XG:Z:GA
#   OT fragment, anchor R1 forward, rescued R2 reverse   -> XG:Z:CT
#   OB fragment, anchor R2 forward, rescued R1 reverse   -> XG:Z:GA
#   contig end,  anchor R1 forward, rescued R2 forward   -> XG:Z:GA
#
# Under the default (collapsed) and genomic scoring each rescued mate must land
# at its true position on its true strand with a full-length CIGAR, NM equal to
# the planted errors only (its conversions are scored free), and XM calling
# exactly its converted cytosines, unmethylated. TAPS (--meth=taps) shares the
# rescue path; its neutral scoring may clip a planted error at a read end, so
# that leg checks strand and XG, and that NM counts no conversion: MD is literal,
# so in every leg NM must equal MD's mismatches at the planted errors' reference
# A/T plus any indel bases, and the full legs also require MD to list each
# converted cytosine.
#
# Inputs:
#   BWA_MEM3 — path to the bwa-mem3 binary under test
set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"
command -v samtools > /dev/null 2>&1 || {
    echo "SKIP: samtools not on PATH (--meth emits BAM)"
    exit 0
}
for tool in python3 mawk; do
    command -v "$tool" > /dev/null 2>&1 || {
        echo "SKIP: $tool not on PATH (needed to build and read the fixture)"
        exit 0
    }
done

# Resolve BWA_MEM3 to an absolute path before we cd into the temp workdir.
case "$BWA_MEM3" in
    /*) BIN="$BWA_MEM3" ;;
    *) BIN="$PWD/$BWA_MEM3" ;;
esac
[ -x "$BIN" ] || {
    echo "FAIL: BWA_MEM3 ($BIN) is not executable" >&2
    exit 1
}

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
fail() {
    echo "FAIL: $*" >&2
    exit 1
}

# Deterministic 1500 bp reference (same PRNG seed 4242 as the sibling meth tests).
REF=TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTAGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGAATCGTCCCTTCTAACGTAGTTCGCATAGTTCCCGTCCGTAGCCGGACTATTCGAACACCCAGTATTCGATTAACTCGGGCTTGACGTATTAGAGGCGTTAGTGTGCCAGGTAAGATACGCCAACGGAATTAACCTCTGTGACACTCCGCGGAGCCTTCGGACATATAAGTGATCGGGTCTACGTTTGTTAGACTTGAGACGTCTGTTAAGAGTTGGGTCTAATAAATCGCCTACACGTGGAGTCTAACGGGGAAGCGTCGAATCCTGATACATCATATAATGGAGCGTGTTATGAAAAAAGAGCATTCCATTGTACGAGCCGTGCCAGAAACGGCTTGACTACGTGAGCGTAGTGTTAGATAAACAGGAAACACTGACGCGGTTAGAAGGCGGATTGCCGGTAGGTTTTGGAAACATAAATACACACGGTATCATGTTGGGTCACGATTCCTATCACCGCACAGGGCCAACCATAGAAGAACTGAAAGAACTAATCTGGCGGCGGGCTCGGTGCTTATATTTTCCACCCAACATCGTGCACATTAGGCTCACCGCGCCCTACGGGCGAAGGGTGCGTACGGTGTTTATAAGGCGTGACGGCCCCAAGTAGAGGGTAATTCTGTGAAAGAATCTCAGGACGGTGGCATGAATTCAATTCCTTTTAAACCTATCGTTCCGACCTTATGCAATCCTTCAATGAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTACTCTTCACTGCCCAAGCGAGGCTAACCTTTTGAAATGTCACAGTCGAAGCATATCTCCCGTACATCTTTTTCGGAGATCGCAGCTCGCGGAGCTATAAGCGACTTAAGCCCTTGTGTCGGTGATCCCAAGGGTCTGACTCCTGTACCAGGGTTACTGTTTCGCTTTACGGAGTAGCCTGTGAGGTGAACTGAAAGGAGCATATTTGAGATCTAAGATAGGGTCCTCCTCTGCGTCTACGTTCTCTCCGTTACGTACGGCTTCGCACCGGAGTGCATCTTGGCCCCGAAACGCACTGTGTGTGCTGATACAGCGTCCCTGGCCGGCCATGGGTTCAGAACTCCCGGGAACGCTTTTCAACTTAGAGGAACCCCGTCATGGAAGTAGATCGCGTCGAATGAGGGAGTTAGTCCTCGTTCCAGCTGGTAATTGTTTTACCGCTTGGGACCACTATAGGCCGCGGGTAGAAGTTGCTGGGTGTTGATTCCAACCCTCGAACCACGATACGACCTGCCATTTATGGCACAGTAAGGTTCAAACAGCATAATGAATACAGTTATAGTAACTTCCTCACGTACGATTAGGACGCAGCCTTG

printf '>chrA\n%s\n' "$REF" > ref.fa
"$BIN" index --meth ref.fa > /dev/null 2>&1 || fail "index --meth nonzero exit"

# Writes r1.fq/r2.fq and expect.tsv, one line per rescued mate:
#   name  mate-flag(64|128)  POS(1-based)  rev(0|1)  XG  converted-offsets  errors
# converted-offsets are the 0-based offsets in SEQ of the mate's converted
# cytosines, comma-separated.
python3 - "$REF" << 'PY' || fail "fixture generation failed"
import random
import sys

ref = sys.argv[1]
L = len(ref)
READ_LEN = 60
INS = 250
QUAL = "I" * READ_LEN
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def rc(s):
    return "".join(COMP[b] for b in reversed(s))


r1, r2, expect = [], [], []

# Plain concordant FR pairs: bwa estimates the insert-size distribution from
# these and declines mate rescue until it can (a lone pair leaves pes failed).
# Their inserts spread around INS so the fixtures' INS sits inside the rescue
# window rather than on its edge, where the window would cut off the mate's
# last base.
rnd = random.Random(7)
for k in range(80):
    ins = rnd.randint(INS - 30, INS + 30)
    p = rnd.randint(40, L - ins - 40)
    r1.append(ref[p : p + READ_LEN])
    r2.append(rc(ref[p + ins - READ_LEN : p + ins]))
names = [f"c{k}" for k in range(80)]


def convert(span, src, dst):
    """The span's bases in genomic letters with src->dst conversion and A<->T
    errors at reference A/T, each at least 12 bases past the previous one, so
    that no error-free stretch reaches a 19-base seed: the mate has to come in
    through mate rescue. Returns (bases, number of errors)."""
    bases = list(ref[span[0] : span[1]].replace(src, dst))
    planted = []
    target = span[0] + 6
    while target < span[1] - 6:
        g = target
        while g < span[1] and ref[g] not in "AT":
            g += 1
        if g >= span[1]:
            break
        bases[g - span[0]] = "T" if bases[g - span[0]] == "A" else "A"
        planted.append(g)
        target = g + 12
    edges = [span[0] - 1] + planted + [span[1]]
    assert max(b - a - 1 for a, b in zip(edges, edges[1:])) < 19, "mate is still seedable"
    return "".join(bases), len(planted)


def add(name, r1_seq, r2_seq, rescued_mate, span, rev, strand, errors):
    src = "C" if strand == "OT" else "G"
    offsets = [i for i in range(READ_LEN) if ref[span[0] + i] == src]
    assert offsets, "mate has no conversions to call"
    names.append(name)
    r1.append(r1_seq)
    r2.append(r2_seq)
    flag = 64 if rescued_mate == 1 else 128
    xg = "CT" if strand == "OT" else "GA"
    expect.append(
        f"{name}\t{flag}\t{span[0] + 1}\t{rev}\t{xg}\t{','.join(map(str, offsets))}\t{errors}"
    )


def fragment(start, strand, rescued_mate):
    """One fragment at [start, start + INS) from original strand OT or OB.

    The fragment's converted strand is written in genomic letters (OT: C->T,
    OB: G->A); R1 reads that strand from its 5' end and R2 reads its copy, so
    for OT R1 is forward and R2 reverse, and for OB the reverse.
    """
    end = start + INS
    src, dst = ("C", "T") if strand == "OT" else ("G", "A")
    fwd_span = (start, start + READ_LEN)
    rev_span = (end - READ_LEN, end)
    r1_span, r1_rev = (fwd_span, 0) if strand == "OT" else (rev_span, 1)
    r2_span, r2_rev = (rev_span, 1) if strand == "OT" else (fwd_span, 0)
    span, rev = (r1_span, r1_rev) if rescued_mate == 1 else (r2_span, r2_rev)
    errors = 0
    reads = {}
    for mate, (s, is_rev) in ((1, (r1_span, r1_rev)), (2, (r2_span, r2_rev))):
        if mate == rescued_mate:
            bases, errors = convert(s, src, dst)
        else:
            bases = ref[s[0] : s[1]].replace(src, dst)
        reads[mate] = rc(bases) if is_rev else bases
    add(f"{strand}_r{rescued_mate}", reads[1], reads[2], rescued_mate, span, rev, strand, errors)


fragment(100, "OT", rescued_mate=1)  # anchor R2 reverse
fragment(500, "OB", rescued_mate=2)  # anchor R1 reverse
fragment(900, "OT", rescued_mate=2)  # anchor R1 forward
fragment(1150, "OB", rescued_mate=1)  # anchor R2 forward

# Contig end: a forward R1 anchor whose FR rescue window runs past the end of
# the last contig, so the clamp puts the window on the reverse half, where it
# finds R2 forward (G->A, its read-number conversion) at the very end.
anchor = (L - 150, L - 90)
mate = (L - READ_LEN, L)
mate_bases, mate_errors = convert(mate, "G", "A")
add(
    "edge_r2",
    ref[anchor[0] : anchor[1]].replace("C", "T"),
    mate_bases,
    2,
    mate,
    0,
    "OB",
    mate_errors,
)

for path, reads in (("r1.fq", r1), ("r2.fq", r2)):
    with open(path, "w") as out:
        for name, seq in zip(names, reads):
            out.write(f"@{name}\n{seq}\n+\n{QUAL}\n")
with open("expect.tsv", "w") as out:
    out.write("\n".join(expect) + "\n")
PY

# Print "FLAG POS CIGAR NM XG CALLS OTHER MISMATCHES_AT_AT MISMATCHES_AT_CG INDELS" for
# the primary record of mate $3 of read $2 in $1, where CALLS are the comma-separated SEQ
# offsets of the XM calls in $4 (a bracket expression), OTHER counts the calls of either
# case not in $4, MISMATCHES_AT_AT / _AT_CG count the mismatches MD lists at a reference
# A/T and C/G, and INDELS the inserted (CIGAR) and deleted (MD) bases.
get() {
    mawk -v name="$2" -v bit="$3" -v want="$4" '
    $1 == name && (int($2/bit) % 2) == 1 && (int($2/256) % 2) == 0 && (int($2/2048) % 2) == 0 {
        nm = ""; xg = ""; calls = ""; other = 0; mm_at = 0; mm_cg = 0; indels = 0
        cig = $6
        while (match(cig, /[0-9]+[MIDNSHP=X]/)) {
            op = substr(cig, RSTART + RLENGTH - 1, 1)
            if (op == "I") indels += substr(cig, RSTART, RLENGTH - 1)
            cig = substr(cig, RSTART + RLENGTH)
        }
        for (i = 12; i <= NF; i++) {
            if ($i ~ /^MD:Z:/) {
                md = substr($i, 6); deleted = 0
                for (j = 1; j <= length(md); j++) {
                    c = substr(md, j, 1)
                    if (c ~ /[0-9]/) deleted = 0
                    else if (c == "^") deleted = 1
                    else if (deleted) indels++
                    else if (c ~ /[CGcg]/) mm_cg++
                    else mm_at++
                }
            }
            if ($i ~ /^NM:i:/) nm = substr($i, 6)
            if ($i ~ /^XG:Z:/) xg = substr($i, 6)
            if ($i ~ /^XM:Z:/) {
                xm = substr($i, 6)
                for (j = 1; j <= length(xm); j++) {
                    c = substr(xm, j, 1)
                    if (c ~ want) calls = calls (calls == "" ? "" : ",") (j - 1)
                    else if (c ~ /[zxhuZXHU]/) other++
                }
            }
        }
        out = $2 " " $4 " " $6 " " nm " " xg " " (calls == "" ? "-" : calls) " " other " " mm_at " " mm_cg " " indels
    }
    END { if (out != "") print out }' "$1"
}

# One leg per scoring mode: $1 label, $2 "full" or "taps", then the mem options.
run_leg() {
    local leg="$1" checks="$2"
    shift 2
    "$BIN" mem "$@" -t 1 ref.fa r1.fq r2.fq 2> /dev/null | samtools view - > "$leg.sam" \
        || fail "$leg: mem $* nonzero exit"
    "$BIN" mem "$@" -S -t 1 ref.fa r1.fq r2.fq 2> /dev/null | samtools view - > "$leg.norescue.sam" \
        || fail "$leg: mem $* -S nonzero exit"
    [ -s "$leg.sam" ] || fail "$leg: mem $* produced no records"
    [ -s "$leg.norescue.sam" ] || fail "$leg: mem $* -S produced no records"

    # EM-seq calls a converted cytosine unmethylated (lower case); TAPS methylated.
    local want='[zxhu]'
    [ "$checks" = taps ] && want='[ZXHU]'
    local n_checked=0 name bit pos rev xg offsets errors label fields
    local flag got_pos cigar nm got_xg calls other mm_at mm_cg indels
    local want_checked=5
    [ "$checks" = taps ] && want_checked=4
    while IFS=$'\t' read -r name bit pos rev xg offsets errors; do
        label="$leg: $name (mate flag $bit)"
        # The contig-end mate is EM-seq-dense in conversions, which TAPS's
        # neutral scoring scores 0, leaving too little for rescue to place it.
        [ "$checks" = taps ] && [ "$name" = edge_r2 ] && continue

        fields="$(get "$leg.norescue.sam" "$name" "$bit" "$want")" || fail "$label: reading the -S run failed"
        [ -n "$fields" ] || fail "$label: no record in the -S run"
        read -r flag _ <<< "$fields"
        (((flag & 4) != 0)) || fail "$label: mapped without mate rescue (flag $flag); the fixture no longer exercises rescue"

        fields="$(get "$leg.sam" "$name" "$bit" "$want")" || fail "$label: reading $leg.sam failed"
        [ -n "$fields" ] || fail "$label: no record"
        read -r flag got_pos cigar nm got_xg calls other mm_at mm_cg indels <<< "$fields"
        (((flag & 4) == 0)) || fail "$label: unmapped (flag $flag); mate rescue did not place it"
        [ "$(((flag >> 4) & 1))" = "$rev" ] || fail "$label: strand rev=$(((flag >> 4) & 1)), want $rev"
        [ "$got_xg" = "$xg" ] || fail "$label: XG:Z:$got_xg, want XG:Z:$xg (rescued mate's original strand)"
        # MD is literal under --meth, so it lists both the planted errors (at reference A/T)
        # and the conversions (at reference C/G). NM must count exactly the former plus the
        # indel bases: a conversion counted in NM fails this wherever the alignment is clipped.
        [ "$nm" -eq $((mm_at + indels)) ] || fail "$label: NM $nm, want $((mm_at + indels)) (MD's $mm_at mismatches at reference A/T + $indels indel bases; conversions must not count)"
        if [ "$checks" = taps ]; then
            [ "$nm" -le "$errors" ] || fail "$label: NM $nm, want at most $errors (only the planted errors; conversions are free)"
        else
            [ "$got_pos" = "$pos" ] || fail "$label: POS $got_pos, want $pos"
            [ "$cigar" = "60M" ] || fail "$label: CIGAR $cigar, want 60M"
            [ "$nm" = "$errors" ] || fail "$label: NM $nm, want $errors (only the planted errors; conversions are free)"
            [ "$calls" = "$offsets" ] || fail "$label: XM calls converted cytosines at $calls, want $offsets"
            [ "$mm_cg" = "$(($(tr -cd , <<< "$offsets" | wc -c) + 1))" ] \
                || fail "$label: MD lists $mm_cg mismatch(es) at a reference C/G, want one per converted cytosine ($offsets)"
            [ "$other" = 0 ] || fail "$label: XM has $other call(s) that are not converted cytosines"
        fi
        n_checked=$((n_checked + 1))
    done < expect.tsv
    [ "$n_checked" -eq "$want_checked" ] || fail "$leg: checked $n_checked rescued mates, want $want_checked"
}

run_leg collapsed full --meth
run_leg genomic full --meth --meth-scoring genomic
run_leg taps taps --meth=taps

echo "PASS: meth_rescue_strand (rescued mates keep their own strand off forward and reverse anchors and at a contig end)"
