#!/usr/bin/env bash
# test/regression/aln_meth_mapq.sh
#
# `bwa-mem3 aln --meth --meth-mapq`: the conversion-blind MAPQ cap, on a seeded synthetic reference
# whose contigs each hold two copies of a 40-base segment that differ only where one copy has the
# converted base's source (C, or G for the G->A strand) and the other already has the converted
# base (T, or A):
#
#   ab  copies A (pos 401) and B (pos 541) differ at two C/T positions (offsets 10 and 25)
#   ef  copies E (pos 401) and F (pos 841) differ at one C/T position (offset 20)
#   gh  copies G (pos 401) and H (pos 841) differ at one G/A position (offset 20)
#
# A read 1 that kept such a C (a methylated C) and converted the rest is unique under genomic
# scoring, because at the other copy the C sits over a T (the mirror cell) and costs a mismatch.
# With the mirror cell free the copies cost the same, so the conversion-blind MAPQ is lower and
# YQ:i keeps the genomic MAPQ. Under aln's collapsed scoring (the mirror cell free, same -M) the
# single-end MAPQ is the independent oracle for the cap.
#
# Checks: the default equals an explicit conversion-blind; single-end reads (one and two retained
# Cs on the C->T strand, the second one past the genomic cut; one retained C on the G->A strand)
# are capped to collapsed's MAPQ with YQ:i, and controls (a fully converted ambiguous read, a unique
# read) are not; only MAPQ and YQ:i differ between the modes; under -n 1 (or a seed's -k 0) a copy
# over the limit only by its mirror mismatches still caps, as collapsed keeps it;
# a locus with several rows enters the cap, single and paired end, at its cheapest row with the
# mirror cell free; the default for collapsed, neutral, TAPS and 5-base is scoring, and collapsed accepts conversion-blind as a no-op; paired end: each
# record's conversion-blind MAPQ is collapsed's, capped at the genomic one, where collapsed places
# it, else 0 (including where collapsed starts an end at a tied mirror-free locus, and where the
# genomic pairing places read 1 at a genomic tie past collapsed's mirror-free cut, and where
# collapsed leaves a pair improper and rescues a mate), a read 1 that collapsed's pairing places at
# the other copy gets 0, a mate the genomic run rescues and collapsed leaves unmapped gets 0, a mate
# that moves read 1 to the other copy anchors it (no cap) unless read 1's mirror-free loci exceed
# --max-occ (then collapsed does not pair it, and its cap is 0), and no MAPQ is ever raised; and
# the option's refusals.
#
# Env:
#   BWA_MEM3  the binary (default: the repository's ./bwa-mem3)
set -euo pipefail
trap 'echo "FAIL: line $LINENO: $BASH_COMMAND" >&2' ERR
here=$(cd "$(dirname "$0")" && pwd)
bin="${BWA_MEM3:-$here/../../bwa-mem3}"
[[ $bin = /* ]] || bin=$PWD/$bin
[[ -x $bin ]] || {
    echo "FAIL: BWA_MEM3 ('$bin') is not an executable binary" >&2
    exit 1
}
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# The reference and reads, from the shared seeded generator helpers (aln_rng.awk).
cat > "$tmp/gen.awk" << 'EOF'
function dna(n,   s, i) { s = ""; for (i = 0; i < n; i++) s = s substr("ACGT", ri(4) + 1, 1); return s }
function put(s, i, c) { return substr(s, 1, i) c substr(s, i + 2) }
# a base other than c, never T (so it cannot match a converted read C)
function alt(c) { return c == "G" ? "A" : "G" }
# from read base b to t everywhere except at the 0-based offsets in keep (" 10 25 ")
function conv(s, b, t, keep,   i, c, out) {
    out = ""
    for (i = 1; i <= length(s); i++) {
        c = substr(s, i, 1)
        out = out ((c == b && index(keep, " " (i - 1) " ") == 0) ? t : c)
    }
    return out
}
function fq(file, name, s,   q) { q = s; gsub(/./, "I", q); printf "@%s\n%s\n+\n%s\n", name, s, q > file }
# a read 2 (G->A) from the reverse strand of contig s at 0-based start st, length n
function mate(s, st, n) { return conv(rc(substr(s, st + 1, n)), "G", "A", "") }
BEGIN {
    seed_rng(20261121)
    a = put(put(dna(40), 10, "C"), 25, "C"); b = put(put(a, 10, "T"), 25, "T")
    ab = dna(400) a dna(100) b dna(400)
    e = put(dna(40), 20, "C"); f = put(e, 20, "T")
    ef = dna(400) e dna(400) f dna(400)
    g = put(dna(40), 20, "G"); h = put(g, 20, "A")
    gh = dna(400) g dna(400) h dna(400)
    # copy Z is copy A with one mismatch; copy Y has a T at A's C and one mismatch, so with the
    # mirror cell free it ties Z (a locus only the mirror-free strata keep)
    ka = put(dna(40), 10, "C")
    z = put(ka, 5, alt(substr(ka, 6, 1)))
    y = put(put(ka, 10, "T"), 33, alt(substr(ka, 34, 1)))
    ayz = dna(400) ka dna(400) y dna(400) z dna(400)
    printf ">ab\n%s\n>ef\n%s\n>gh\n%s\n>ayz\n%s\n", ab, ef, gh, ayz > (dir "/ref.fa")
    # single end
    fq(dir "/se.fq", "ab2", conv(a, "C", "T", " 10 25 "))
    fq(dir "/se.fq", "ab1", conv(a, "C", "T", " 10 "))
    fq(dir "/se.fq", "ef1", conv(e, "C", "T", " 20 "))
    fq(dir "/se.fq", "gh1", conv(rc(g), "C", "T", " 19 "))
    fq(dir "/se.fq", "ef0", conv(e, "C", "T", ""))
    fq(dir "/se.fq", "uniq", conv(substr(ef, 101, 40), "C", "T", ""))
    # paired end, one batch on ab (enough pairs for an insert-size estimate): read 1 keeps both Cs
    # of copy A, read 2 starts at st
    r1 = conv(a, "C", "T", " 10 25 ")
    for (st = 440; st <= 900; st += 20) { fq(dir "/p1.fq", "m" st, r1); fq(dir "/p2.fq", "m" st, mate(ab, st, 60)) }
    # a read 2 with too many mismatches for the search, so mate rescue places it next to copy A
    m = mate(ab, 600, 60)
    for (k = 5; k <= 55; k += 10) m = put(m, k, other(substr(m, k + 1, 1)))
    fq(dir "/p1.fq", "resc", r1); fq(dir "/p2.fq", "resc", m)
    # paired end, a single pair on ef: read 2 from the flank after copy F moves read 1 there
    fq(dir "/q1.fq", "anchor", conv(e, "C", "T", " 20 ")); fq(dir "/q2.fq", "anchor", mate(ef, 1180, 60))
    # paired end, a single pair on ayz: read 2 from the flank after copy Z moves read 1 there
    fq(dir "/o1.fq", "maxocc", conv(ka, "C", "T", " 10 ")); fq(dir "/o2.fq", "maxocc", mate(ayz, 1580, 60))
}
EOF
awk -v dir="$tmp" -f "$here/aln_rng.awk" -f "$tmp/gen.awk" < /dev/null
"$bin" index --meth "$tmp/ref.fa" > "$tmp/index.out" 2> "$tmp/index.err" || {
    echo "FAIL: index --meth exited nonzero" >&2
    cat "$tmp/index.err" >&2
    exit 1
}

# grep -v exits 1 when there are no records: that is an empty result, not an error
records() { grep -v '^@' "$1" || [ $? -eq 1 ]; }
# name FLAG RNAME POS MAPQ YQ (or -) of each record, one per line
summary() { awk '{ yq = "-"; for (i = 12; i <= NF; ++i) if (index($i, "YQ:i:") == 1) yq = substr($i, 6); print $1, $2, $3, $4, $5, yq }' "$1"; }
# the records with MAPQ and YQ:i removed
strip() { awk 'BEGIN { OFS = "\t" } { $5 = "."; out = $1; for (i = 2; i <= NF; ++i) if (index($i, "YQ:i:") != 1) out = out "\t" $i; print out }' "$1"; }
run() {
    local out=$1
    shift
    # the EXIT trap removes $tmp: show the aligner's stderr before it does
    "$bin" aln "$@" > "$tmp/$out.sam" 2> "$tmp/$out.err" || {
        echo "FAIL: $out: aln exited nonzero" >&2
        cat "$tmp/$out.err" >&2
        exit 1
    }
    records "$tmp/$out.sam" > "$tmp/$out.rec"
    summary "$tmp/$out.rec" > "$tmp/$out.sum"
}
pick() { awk -v n="$2" '$1 == n' "$1"; }
expect() {
    local name=$1 file=$2 want=$3
    [ "$(cat "$file")" = "$want" ] || {
        echo "FAIL: $name: expected" >&2
        printf '%s\n' "$want" >&2
        echo "got" >&2
        cat "$file" >&2
        exit 1
    }
}

# The conversion-blind MAPQ observes collapsed's: each record (paired runs in the same order) has
# the lower of the scoring and collapsed MAPQs where collapsed places it at the same position and
# strand, else 0 (collapsed's pairing and mate rescue are rerun for the cap). Reads
# $tmp/<name>.{blind,scoring,collapsed}.sum.
observer() {
    local label=$1 name=$2 rc=0
    paste -d ' ' "$tmp/$name.blind.sum" "$tmp/$name.scoring.sum" "$tmp/$name.collapsed.sum" \
        | awk '{
            # a line pairs the same record of each run: its QNAME and read 1/read 2 flags agree
            if ($1 != $7 || $1 != $13 || int($2 / 64) % 4 != int($8 / 64) % 4 || int($2 / 64) % 4 != int($14 / 64) % 4) {
                ident = 1; print "  " $0 > "/dev/stderr"; next
            }
            same = $3 == $15 && $4 == $16 && int($2 / 16) % 2 == int($14 / 16) % 2
            want = !same ? 0 : $11 < $17 ? $11 : $17
            if (int($2 / 4) % 2 == 0 && $5 != want) { bad = 1; print "  " $0 > "/dev/stderr" }
        } END { exit ident ? 2 : bad }' || rc=$?
    [ "$rc" -ne 2 ] || {
        echo "FAIL: $label: the runs do not list the same records (QNAME, read 1/read 2) in the same order" >&2
        exit 1
    }
    [ "$rc" -eq 0 ] || {
        echo "FAIL: $label: a conversion-blind MAPQ is not collapsed's, capped at the genomic one" >&2
        exit 1
    }
}

# Single end.
run se.blind --meth "$tmp/ref.fa" "$tmp/se.fq"
run se.explicit --meth --meth-mapq conversion-blind "$tmp/ref.fa" "$tmp/se.fq"
run se.scoring --meth --meth-mapq scoring "$tmp/ref.fa" "$tmp/se.fq"
run se.collapsed --meth --meth-scoring collapsed "$tmp/ref.fa" "$tmp/se.fq"
cmp "$tmp/se.blind.rec" "$tmp/se.explicit.rec" || {
    echo "FAIL: default: differs from an explicit --meth-mapq conversion-blind" >&2
    exit 1
}
echo 'PASS: default-is-conversion-blind'

# ab2's second copy costs two mismatches under genomic, past the genomic cut (X1:i:0): the cap
# must count it from the candidates before the cut.
grep -E '^ab2	' "$tmp/se.scoring.rec" | grep -qE '	X0:i:1	X1:i:0	' || {
    echo "FAIL: fixture: ab2's other copy is not past the genomic cut" >&2
    exit 1
}
expect se-scoring "$tmp/se.scoring.sum" "ab2 0 ab 401 37 -
ab1 0 ab 401 23 -
ef1 0 ef 401 23 -
gh1 16 gh 401 23 -
ef0 0 ef 841 0 -
uniq 0 ef 101 37 -"
expect se-conversion-blind "$tmp/se.blind.sum" "ab2 0 ab 401 0 37
ab1 0 ab 401 0 23
ef1 0 ef 401 0 23
gh1 16 gh 401 0 23
ef0 0 ef 841 0 -
uniq 0 ef 101 37 -"
# the oracle: collapsed's MAPQ (placement can differ, at an equal-cost copy)
paste -d ' ' "$tmp/se.blind.sum" "$tmp/se.collapsed.sum" | awk '$5 != $11 { bad = 1; print "  " $0 > "/dev/stderr" } END { exit bad }' || {
    echo "FAIL: single end: the conversion-blind MAPQ differs from collapsed's" >&2
    exit 1
}
echo 'PASS: single-end-capped-to-collapsed'

cmp <(strip "$tmp/se.blind.rec") <(strip "$tmp/se.scoring.rec") || {
    echo "FAIL: modes differ outside MAPQ and YQ:i" >&2
    exit 1
}
echo 'PASS: only-mapq-and-yq-differ'

# Under -n 1, ab2's other copy (two mirror mismatches) is over the limit genomically, so the search
# drops it, but within it once the mirror cell is free: the cap must still count it, as collapsed
# keeps it as a locus.
run se.n1.blind --meth -n 1 "$tmp/ref.fa" "$tmp/se.fq"
run se.n1.collapsed --meth --meth-scoring collapsed -n 1 "$tmp/ref.fa" "$tmp/se.fq"
expect se-n1-capped <(pick "$tmp/se.n1.blind.sum" ab2) "ab2 0 ab 401 0 37"
paste -d ' ' "$tmp/se.n1.blind.sum" "$tmp/se.n1.collapsed.sum" | awk '$5 != $11 { bad = 1; print "  " $0 > "/dev/stderr" } END { exit bad }' || {
    echo "FAIL: single end, -n 1: the conversion-blind MAPQ differs from collapsed's" >&2
    exit 1
}
echo 'PASS: single-end-mirror-free-n-limit'

# The same for the seed's -k: with a 20-base seed and -k 0, ab2's other copy has one mismatch in the
# seed (its mirror mismatch at offset 10), so the search drops it, but none once the mirror cell is
# free: the cap must still count it, as collapsed keeps it as a locus.
run se.k0.blind --meth -l 20 -k 0 "$tmp/ref.fa" "$tmp/se.fq"
run se.k0.collapsed --meth --meth-scoring collapsed -l 20 -k 0 "$tmp/ref.fa" "$tmp/se.fq"
expect se-k0-capped <(pick "$tmp/se.k0.blind.sum" ab2 | awk '{ print $1, $2, $3, $4, $5 }') "ab2 0 ab 401 0"
paste -d ' ' "$tmp/se.k0.blind.sum" "$tmp/se.k0.collapsed.sum" | awk '$5 != $11 { bad = 1; print "  " $0 > "/dev/stderr" } END { exit bad }' || {
    echo "FAIL: single end, -k 0: the conversion-blind MAPQ differs from collapsed's" >&2
    exit 1
}
echo 'PASS: single-end-mirror-free-k-limit'

# A locus holds several rows (search paths with a gap in different places); the genomic search keeps
# the genomically cheapest, collapsed the cheapest with the mirror cell free, and the cap takes each
# locus at that row. On a reference of its own (row.fa), read 1 'tie' (50 bp) has ACA at offsets
# 16-18; contig ta holds it with ATAA there (one base longer) and a mismatch at offset 40, so the
# deletion gives two rows of genomic cost 17: two mismatches (the lower CIGAR, which the genomic
# search keeps) or one and a read C over a T (mirror-free 14). Contig tb holds it with a two-base
# deletion (15). The genomic MAPQ is 23 (ta within tb + 3); with ta at 14 the cap is 0, and
# collapsed places the read at ta. Read 'rep' has CAC there; contig ra holds TCTA, whose rows cost
# 17 (two mismatches) and 20 (two mirror mismatches, mirror-free 14): the first replaces the second.
# Contig rb has four mismatches (12): the genomic MAPQ is 37 (ra past rb + 3); with ra at 14 the cap
# is collapsed's 23.
cat > "$tmp/row.awk" << 'EOF'
function dna(n,   s, i) { s = ""; for (i = 0; i < n; i++) s = s substr("ACGT", ri(4) + 1, 1); return s }
function put(s, i, c) { return substr(s, 1, i) c substr(s, i + 2) }
function alt(c) { return c == "G" ? "A" : "G" }
function fq(file, name, s,   q) { q = s; gsub(/./, "I", q); printf "@%s\n%s\n+\n%s\n", name, s, q > file }
BEGIN {
    seed_rng(20261123)
    t = dna(50); t = substr(t, 1, 16) "ACA" substr(t, 20)
    ta = substr(t, 1, 16) "ATAA" substr(t, 20); ta = put(ta, 41, alt(substr(t, 41, 1)))
    x = alt(substr(t, 31, 1)); tb = substr(t, 1, 30) x alt(x) substr(t, 31)
    r = dna(50); r = substr(r, 1, 16) "CAC" substr(r, 20)
    ra = substr(r, 1, 16) "TCTA" substr(r, 20)
    rb = r
    split("5 25 33 44", at, " ")
    for (k = 1; k <= 4; ++k) rb = put(rb, at[k], alt(substr(r, at[k] + 1, 1)))
    printf ">ta\n%s\n>tb\n%s\n>ra\n%s\n>rb\n%s\n", dna(100) ta dna(100), dna(100) tb dna(100), dna(100) ra dna(100), dna(100) rb dna(100) > (dir "/row.fa")
    fq(dir "/row.fq", "tie", t)
    fq(dir "/row.fq", "rep", r)
}
EOF
awk -v dir="$tmp" -f "$here/aln_rng.awk" -f "$tmp/row.awk" < /dev/null
"$bin" index --meth "$tmp/row.fa" > "$tmp/row-index.out" 2> "$tmp/row-index.err" || {
    echo "FAIL: index --meth (row.fa) exited nonzero" >&2
    cat "$tmp/row-index.err" >&2
    exit 1
}
run row.scoring --meth --meth-mapq scoring -n 5 -l 1024 "$tmp/row.fa" "$tmp/row.fq"
run row.blind --meth -n 5 -l 1024 "$tmp/row.fa" "$tmp/row.fq"
run row.collapsed --meth --meth-scoring collapsed -n 5 -l 1024 "$tmp/row.fa" "$tmp/row.fq"
expect row-scoring "$tmp/row.scoring.sum" "tie 0 tb 101 23 -
rep 0 rb 101 37 -"
# the genomic search keeps ta's two-mismatch row (the lower CIGAR), which XA shows
grep -E '^tie	' "$tmp/row.scoring.rec" | grep -qF '	XA:Z:ta,+101,17M1D33M,3;' || {
    echo "FAIL: fixture: the genomic search does not keep ta's two-mismatch row" >&2
    exit 1
}
expect row-conversion-blind "$tmp/row.blind.sum" "tie 0 tb 101 0 23
rep 0 rb 101 23 37"
expect row-collapsed "$tmp/row.collapsed.sum" "tie 0 ta 101 23 -
rep 0 rb 101 23 -"
cmp <(strip "$tmp/row.blind.rec") <(strip "$tmp/row.scoring.rec") || {
    echo "FAIL: row: modes differ outside MAPQ and YQ:i" >&2
    exit 1
}
echo 'PASS: single-end-lowest-mirror-free-row'

# Paired end, on a reference of its own (rowpe.fa): read 1 is 'tie' again, at contigs ta (one
# deletion row of genomic 17 and mirror-free 14) and tb (a two-base deletion, 15), and read 2 comes
# from after tb, so the genomic pairing (scoring and conversion-blind alike) places the pair at tb.
# The shadow pairing must see ta at 14, where collapsed starts read 1 (its only cheapest locus
# mirror-free). Collapsed's insert size comes from those starts, all on another contig than read 2,
# so it has none and pairs by -a 500: the 11 pairs within it move read 1 to tb (MAPQ 29 and 37), the
# other 13 leave it at ta (MAPQ 23; read 2 37). So the conversion-blind MAPQs are 29 and 37, and 0
# for read 1 of the 13. At ta's genomic row (17) the shadow is the real pairing and nothing is
# capped.
cat > "$tmp/rowpe.awk" << 'EOF'
function dna(n,   s, i) { s = ""; for (i = 0; i < n; i++) s = s substr("ACGT", ri(4) + 1, 1); return s }
function put(s, i, c) { return substr(s, 1, i) c substr(s, i + 2) }
function alt(c) { return c == "G" ? "A" : "G" }
function fq(file, name, s,   q) { q = s; gsub(/./, "I", q); printf "@%s\n%s\n+\n%s\n", name, s, q > file }
function ga(s,   i, c, out) { out = ""; for (i = 1; i <= length(s); i++) { c = substr(s, i, 1); out = out (c == "G" ? "A" : c) } return out }
BEGIN {
    seed_rng(20261123)
    t = dna(50); t = substr(t, 1, 16) "ACA" substr(t, 20)
    ta = substr(t, 1, 16) "ATAA" substr(t, 20); ta = put(ta, 41, alt(substr(t, 41, 1)))
    x = alt(substr(t, 31, 1)); tb = substr(t, 1, 30) x alt(x) substr(t, 31)
    ca = dna(100) ta dna(100)
    cb = dna(100) tb dna(900)
    printf ">ta\n%s\n>tb\n%s\n", ca, cb > (dir "/rowpe.fa")
    for (st = 340; st <= 800; st += 20) {
        fq(dir "/rowpe1.fq", "m" st, t)
        fq(dir "/rowpe2.fq", "m" st, ga(rc(substr(cb, st + 1, 60))))
    }
}
EOF
awk -v dir="$tmp" -f "$here/aln_rng.awk" -f "$tmp/rowpe.awk" < /dev/null
"$bin" index --meth "$tmp/rowpe.fa" > "$tmp/rowpe-index.out" 2> "$tmp/rowpe-index.err" || {
    echo "FAIL: index --meth (rowpe.fa) exited nonzero" >&2
    cat "$tmp/rowpe-index.err" >&2
    exit 1
}
run rowpe.scoring --meth --meth-mapq scoring -n 5 -l 1024 "$tmp/rowpe.fa" "$tmp/rowpe1.fq" "$tmp/rowpe2.fq"
run rowpe.blind --meth -n 5 -l 1024 "$tmp/rowpe.fa" "$tmp/rowpe1.fq" "$tmp/rowpe2.fq"
run rowpe.collapsed --meth --meth-scoring collapsed -n 5 -l 1024 "$tmp/rowpe.fa" "$tmp/rowpe1.fq" "$tmp/rowpe2.fq"
# FLAG RNAME MAPQ YQ and how many of the 48 records have them
tally() { awk '{ n[$2 " " $3 " " $5 " " $6]++ } END { for (k in n) print n[k], k }' "$1" | LC_ALL=C sort; }
expect rowpe-scoring <(tally "$tmp/rowpe.scoring.sum") "24 147 tb 60 -
24 99 tb 60 -"
expect rowpe-collapsed <(tally "$tmp/rowpe.collapsed.sum") "11 147 tb 37 -
11 99 tb 29 -
13 145 tb 37 -
13 97 ta 23 -"
expect rowpe-conversion-blind <(tally "$tmp/rowpe.blind.sum") "11 99 tb 29 60
13 99 tb 0 60
24 147 tb 37 60"
observer rowpe rowpe
cmp <(strip "$tmp/rowpe.blind.rec") <(strip "$tmp/rowpe.scoring.rec") || {
    echo "FAIL: rowpe: modes differ outside MAPQ and YQ:i" >&2
    exit 1
}
echo 'PASS: paired-lowest-mirror-free-row'

# Paired end, collapsed's starts (on a reference of its own, start.fa): pairing's MAPQ depends on
# where the ends start, and collapsed starts each end at its own draw among the loci tied with the
# mirror cell free. Each pair has a home contig (read 1's 40-base segment, read 2's 300 bases on)
# and copies of a segment on contigs of their own: exact, or with a T under the C the read keeps
# (one mismatch genomically, none mirror-free). Sixteen copies of each pair draw with sixteen seeds;
# -A leaves out the insert-size estimate (and mate rescue).
#   tie   read 2 has an exact and a mirror copy: where collapsed starts read 2 off home, pairing
#         moves it there (29) and the genomic pair has 36
#   low   read 2 has a mirror copy: likewise 29 where the genomic pair has 60
#   both  each read has a mirror copy: where collapsed starts both off home both move (9), and the
#         genomic pair has 46
cat > "$tmp/start.awk" << 'EOF'
function dna(n,   s, i) { s = ""; for (i = 0; i < n; i++) s = s substr("ACGT", ri(4) + 1, 1); return s }
function put(s, i, c) { return substr(s, 1, i) c substr(s, i + 2) }
function fq(file, name, s,   q) { q = s; gsub(/./, "I", q); printf "@%s\n%s\n+\n%s\n", name, s, q > file }
function ga(s,   i, c, out) { out = ""; for (i = 1; i <= length(s); i++) { c = substr(s, i, 1); out = out (c == "G" ? "A" : c) } return out }
function ct(s,   i, c, out) { out = ""; for (i = 1; i <= length(s); i++) { c = substr(s, i, 1); out = out (c == "C" ? "T" : c) } return out }
# a contig of its own holding seg between two random flanks
function copy(name, seg) { printf ">%s\n%s\n", name, dna(100) seg dna(100) > (dir "/start.fa") }
# pair name: its home contig, and an exact copy of read 2's segment (dup2), a mirror copy of read
# 1's (mir1) and one of read 2's (mir2); read 1 keeps the C at offset 12, read 2 the one at 20
function pair(name, dup2, mir1, mir2,   s1, s2, a, b, k) {
    s1 = put(dna(40), 12, "C")
    s2 = put(dna(40), 20, "C")
    printf ">%s\n%s\n", name, dna(100) s1 dna(260) s2 dna(100) > (dir "/start.fa")
    if (dup2) copy(name "_dup2", s2)
    if (mir1) copy(name "_mir1", put(s1, 12, "T"))
    if (mir2) copy(name "_mir2", put(s2, 20, "T"))
    a = put(ct(s1), 12, "C")
    b = put(ga(rc(s2)), 19, "G")
    for (k = 0; k < 16; ++k) {
        fq(dir "/start1.fq", name k, a)
        fq(dir "/start2.fq", name k, b)
    }
}
BEGIN {
    seed_rng(20261124)
    pair("tie", 1, 0, 1)
    pair("low", 0, 0, 1)
    pair("both", 0, 1, 1)
}
EOF
awk -v dir="$tmp" -f "$here/aln_rng.awk" -f "$tmp/start.awk" < /dev/null
"$bin" index --meth "$tmp/start.fa" > "$tmp/start-index.out" 2> "$tmp/start-index.err" || {
    echo "FAIL: index --meth (start.fa) exited nonzero" >&2
    cat "$tmp/start-index.err" >&2
    exit 1
}
run start.scoring --meth --meth-mapq scoring -A "$tmp/start.fa" "$tmp/start1.fq" "$tmp/start2.fq"
run start.blind --meth -A "$tmp/start.fa" "$tmp/start1.fq" "$tmp/start2.fq"
run start.collapsed --meth --meth-scoring collapsed -A "$tmp/start.fa" "$tmp/start1.fq" "$tmp/start2.fq"
cmp <(strip "$tmp/start.blind.rec") <(strip "$tmp/start.scoring.rec") || {
    echo "FAIL: start: modes differ outside MAPQ and YQ:i" >&2
    exit 1
}
observer start start
# each pair (its name less the copy number) with the genomic and conversion-blind MAPQs of a record
shapes() { paste -d ' ' "$tmp/start.scoring.sum" "$tmp/start.blind.sum" | awk '{ sub(/[0-9]+$/, "", $1); print $1, $5, $11 }' | LC_ALL=C sort -u; }
for want in 'tie 36 29' 'low 60 29' 'both 46 9'; do
    shapes | grep -qx "$want" || {
        echo "FAIL: start: no record with (pair, genomic, conversion-blind) MAPQs $want" >&2
        shapes >&2
        exit 1
    }
done
echo 'PASS: paired-collapsed-starts'

# Paired end, a genomic hit past the mirror-free cut (cut.fa): read 1 keeps two Cs; contig cx has a
# T under both (two mismatches genomically, none mirror-free) and contig ch two real mismatches, so
# the two tie genomically, but with the mirror cell free ch is past cx + 3 and is not one of
# collapsed's loci. Read 2 comes from after ch: the genomic pairing places read 1 there (MAPQ 29),
# collapsed cannot (it leaves read 1 at cx), so read 1's cap is 0.
cat > "$tmp/cut.awk" << 'EOF'
function dna(n,   s, i) { s = ""; for (i = 0; i < n; i++) s = s substr("ACGT", ri(4) + 1, 1); return s }
function put(s, i, c) { return substr(s, 1, i) c substr(s, i + 2) }
function alt(c) { return c == "G" ? "A" : "G" }
function fq(file, name, s,   q) { q = s; gsub(/./, "I", q); printf "@%s\n%s\n+\n%s\n", name, s, q > file }
function ga(s,   i, c, out) { out = ""; for (i = 1; i <= length(s); i++) { c = substr(s, i, 1); out = out (c == "G" ? "A" : c) } return out }
function ct(s,   i, c, out) { out = ""; for (i = 1; i <= length(s); i++) { c = substr(s, i, 1); out = out (c == "C" ? "T" : c) } return out }
BEGIN {
    seed_rng(7)
    s = put(put(dna(40), 10, "C"), 28, "C")
    x = put(put(s, 10, "T"), 28, "T")
    h = put(put(s, 5, alt(substr(s, 6, 1))), 20, alt(substr(s, 21, 1)))
    cx = dna(100) x dna(900)
    ch = dna(100) h dna(900)
    printf ">cx\n%s\n>ch\n%s\n", cx, ch > (dir "/cut.fa")
    a = put(put(ct(s), 10, "C"), 28, "C")
    for (st = 340; st <= 800; st += 20) {
        fq(dir "/cut1.fq", "m" st, a)
        fq(dir "/cut2.fq", "m" st, ga(rc(substr(ch, st + 1, 60))))
    }
}
EOF
awk -v dir="$tmp" -f "$here/aln_rng.awk" -f "$tmp/cut.awk" < /dev/null
"$bin" index --meth "$tmp/cut.fa" > "$tmp/cut-index.out" 2> "$tmp/cut-index.err" || {
    echo "FAIL: index --meth (cut.fa) exited nonzero" >&2
    cat "$tmp/cut-index.err" >&2
    exit 1
}
run cut.scoring --meth --meth-mapq scoring "$tmp/cut.fa" "$tmp/cut1.fq" "$tmp/cut2.fq"
run cut.blind --meth "$tmp/cut.fa" "$tmp/cut1.fq" "$tmp/cut2.fq"
run cut.collapsed --meth --meth-scoring collapsed "$tmp/cut.fa" "$tmp/cut1.fq" "$tmp/cut2.fq"
cmp <(strip "$tmp/cut.blind.rec") <(strip "$tmp/cut.scoring.rec") || {
    echo "FAIL: cut: modes differ outside MAPQ and YQ:i" >&2
    exit 1
}
observer cut cut
# read 1 records the genomic pairing places at ch with a nonzero MAPQ, and their conversion-blind MAPQs
past=$(awk '$2 % 256 >= 64 && $2 % 256 < 128 && $3 == "ch" && $5 > 0 { print $1 }' "$tmp/cut.scoring.sum")
[ -n "$past" ] || {
    echo "FAIL: fixture: the genomic pairing places no read 1 at ch with a nonzero MAPQ" >&2
    exit 1
}
awk '$2 % 256 >= 64 && $2 % 256 < 128 && $3 == "ch" && $5 != 0 { bad = 1; print "  " $0 > "/dev/stderr" } END { exit bad }' "$tmp/cut.blind.sum" || {
    echo "FAIL: cut: a read 1 at ch, past collapsed's cut, is not capped at 0" >&2
    exit 1
}
echo 'PASS: paired-past-mirror-free-cut'

# Paired end, collapsed's mate rescue (mresc.fa): read 1 is unique on contig home; read 2's segment
# lies 280 bases on with two real mismatches there, and on contig copy with a T under each of the
# three Cs read 2 keeps (three mismatches in the seed genomically, so -k drops it; none mirror-free).
# The genomic pair is proper at home (MAPQ 60). Collapsed places read 2 at copy (home is past its
# cut), pairs nothing, and rescues read 2 home from read 1, which lowers both ends to 11. So both
# conversion-blind MAPQs are 11. Pairs on contigs of their own give both runs an insert size.
cat > "$tmp/mresc.awk" << 'EOF'
function dna(n,   s, i) { s = ""; for (i = 0; i < n; i++) s = s substr("ACGT", ri(4) + 1, 1); return s }
function put(s, i, c) { return substr(s, 1, i) c substr(s, i + 2) }
function alt(c) { return c == "G" ? "A" : "G" }
function fq(file, name, s,   q) { q = s; gsub(/./, "I", q); printf "@%s\n%s\n+\n%s\n", name, s, q > file }
function ga(s,   i, c, out) { out = ""; for (i = 1; i <= length(s); i++) { c = substr(s, i, 1); out = out (c == "G" ? "A" : c) } return out }
function ct(s,   i, c, out) { out = ""; for (i = 1; i <= length(s); i++) { c = substr(s, i, 1); out = out (c == "C" ? "T" : c) } return out }
BEGIN {
    seed_rng(20261124)
    for (k = 0; k < 24; ++k) {
        s1 = dna(40); s2 = dna(40)
        printf ">bg%d\n%s\n", k, dna(100) s1 dna(200 + 3 * k) s2 dna(100) > (dir "/mresc.fa")
        fq(dir "/mresc1.fq", "bg" k, ct(s1)); fq(dir "/mresc2.fq", "bg" k, ga(rc(s2)))
    }
    s1 = dna(40)
    s2 = put(put(put(dna(40), 10, "C"), 20, "C"), 30, "C")
    h2 = s2
    for (p = 5; p <= 35; p += 30) h2 = put(h2, p, alt(substr(ct(s2), p + 1, 1)))
    printf ">home\n%s\n", dna(100) s1 dna(240) h2 dna(100) > (dir "/mresc.fa")
    printf ">copy\n%s\n", dna(100) put(put(put(s2, 10, "T"), 20, "T"), 30, "T") dna(100) > (dir "/mresc.fa")
    fq(dir "/mresc1.fq", "target", ct(s1))
    fq(dir "/mresc2.fq", "target", put(put(put(ga(rc(s2)), 29, "G"), 19, "G"), 9, "G"))
}
EOF
awk -v dir="$tmp" -f "$here/aln_rng.awk" -f "$tmp/mresc.awk" < /dev/null
"$bin" index --meth "$tmp/mresc.fa" > "$tmp/mresc-index.out" 2> "$tmp/mresc-index.err" || {
    echo "FAIL: index --meth (mresc.fa) exited nonzero" >&2
    cat "$tmp/mresc-index.err" >&2
    exit 1
}
run mresc.scoring --meth --meth-mapq scoring "$tmp/mresc.fa" "$tmp/mresc1.fq" "$tmp/mresc2.fq"
run mresc.blind --meth "$tmp/mresc.fa" "$tmp/mresc1.fq" "$tmp/mresc2.fq"
run mresc.collapsed --meth --meth-scoring collapsed "$tmp/mresc.fa" "$tmp/mresc1.fq" "$tmp/mresc2.fq"
cmp <(strip "$tmp/mresc.blind.rec") <(strip "$tmp/mresc.scoring.rec") || {
    echo "FAIL: mresc: modes differ outside MAPQ and YQ:i" >&2
    exit 1
}
observer mresc mresc
expect collapsed-rescue-genomic <(pick "$tmp/mresc.scoring.sum" target) "target 99 home 101 60 -
target 147 home 381 60 -"
grep -E '^target	147	' "$tmp/mresc.collapsed.rec" | grep -q '	XT:A:M	' || {
    echo "FAIL: fixture: collapsed does not rescue target's read 2 (XT:A:M)" >&2
    exit 1
}
expect collapsed-rescue <(pick "$tmp/mresc.blind.sum" target) "target 99 home 101 11 60
target 147 home 381 11 60"
echo 'PASS: paired-collapsed-mate-rescue'

# The default is scoring outside em-seq genomic; collapsed accepts conversion-blind (a no-op).
for args in '--meth --meth-scoring collapsed' '--meth=taps' '--meth=5base' '--meth --meth-scoring neutral'; do
    read -r -a options <<< "$args"
    run other "${options[@]}" "$tmp/ref.fa" "$tmp/se.fq"
    run other.scoring "${options[@]}" --meth-mapq scoring "$tmp/ref.fa" "$tmp/se.fq"
    cmp "$tmp/other.rec" "$tmp/other.scoring.rec" || {
        echo "FAIL: $args: the default is not --meth-mapq scoring" >&2
        exit 1
    }
done
run other.blind --meth --meth-scoring collapsed --meth-mapq conversion-blind "$tmp/ref.fa" "$tmp/se.fq"
cmp "$tmp/se.collapsed.rec" "$tmp/other.blind.rec" || {
    echo "FAIL: collapsed: conversion-blind is not a no-op" >&2
    exit 1
}
echo 'PASS: scoring-outside-emseq-genomic'

# Paired end.
run pe.blind --meth "$tmp/ref.fa" "$tmp/p1.fq" "$tmp/p2.fq"
run pe.scoring --meth --meth-mapq scoring "$tmp/ref.fa" "$tmp/p1.fq" "$tmp/p2.fq"
run pe.collapsed --meth --meth-scoring collapsed "$tmp/ref.fa" "$tmp/p1.fq" "$tmp/p2.fq"
grep -q 'infer_isize] (25, 50, 75) percentile' "$tmp/pe.blind.err" || {
    echo "FAIL: fixture: the paired batch gave no insert-size estimate" >&2
    exit 1
}
cmp <(strip "$tmp/pe.blind.rec") <(strip "$tmp/pe.scoring.rec") || {
    echo "FAIL: paired: modes differ outside MAPQ and YQ:i" >&2
    exit 1
}
paste -d ' ' "$tmp/pe.blind.sum" "$tmp/pe.scoring.sum" | awk '$5 > $11 { bad = 1; print "  " $0 > "/dev/stderr" } END { exit bad }' || {
    echo "FAIL: paired: conversion-blind raised a MAPQ" >&2
    exit 1
}
echo 'PASS: paired-only-lowers'

observer paired pe
# m660: both copies pair with read 2; collapsed draws read 1 at copy A and pairs it there, with the
# MAPQ of a pair with a close suboptimal pair, and read 2 keeps its own.
expect paired-near <(pick "$tmp/pe.blind.sum" m660) "m660 99 ab 401 7 60
m660 147 ab 661 37 60"
[ "$(pick "$tmp/pe.collapsed.sum" m660 | awk '$2 == 99 { print $4 }')" = 401 ] || {
    echo "FAIL: fixture: collapsed does not place m660's read 1 at copy A" >&2
    exit 1
}
# m880: collapsed's pairing places read 1 at copy B, so its cap is 0.
expect paired-elsewhere <(pick "$tmp/pe.blind.sum" m880) "m880 99 ab 401 0 60
m880 147 ab 881 37 60"
[ "$(pick "$tmp/pe.collapsed.sum" m880 | awk '$2 == 99 { print $4 }')" = 541 ] || {
    echo "FAIL: fixture: collapsed does not place m880's read 1 at copy B" >&2
    exit 1
}
# resc: read 2 is found only by mate rescue (XT:A:M) next to read 1, so there is no pairing. Read 1
# gets collapsed's MAPQ, 0 (copy B ties once the mirror cell is free), too low to anchor collapsed's
# rescue, so collapsed leaves read 2 unmapped and its cap is 0.
grep -E '^resc	' "$tmp/pe.blind.rec" | grep -q '	XT:A:M	' || {
    echo "FAIL: fixture: resc's read 2 is not rescued (XT:A:M)" >&2
    exit 1
}
expect paired-rescued <(pick "$tmp/pe.blind.sum" resc) "resc 99 ab 401 0 29
resc 147 ab 601 0 29"
echo 'PASS: paired-capped'

# A single pair: read 2 lies after copy F, so pairing moves read 1 there in both modes; the mate
# anchors it and nothing is capped.
run q.blind --meth "$tmp/ref.fa" "$tmp/q1.fq" "$tmp/q2.fq"
run q.scoring --meth --meth-mapq scoring "$tmp/ref.fa" "$tmp/q1.fq" "$tmp/q2.fq"
cmp "$tmp/q.blind.rec" "$tmp/q.scoring.rec" || {
    echo "FAIL: anchored pair: conversion-blind changed it" >&2
    exit 1
}
expect paired-anchor "$tmp/q.blind.sum" "anchor 99 ef 841 29 -
anchor 147 ef 1181 37 -"
echo 'PASS: paired-mate-anchors'

# The same move with read 1's three mirror-free loci over --max-occ 2 (its two stock loci are not):
# the pairing still moves read 1 to copy Z, but collapsed would not pair it and would leave it at its
# single-end draw, so its cap is 0. Under the default --max-occ collapsed's pairing anchors it too.
run o.blind --meth "$tmp/ref.fa" "$tmp/o1.fq" "$tmp/o2.fq"
expect maxocc-default "$tmp/o.blind.sum" "maxocc 99 ayz 1281 29 -
maxocc 147 ayz 1581 37 -"
run o2.blind --meth --max-occ 2 "$tmp/ref.fa" "$tmp/o1.fq" "$tmp/o2.fq"
expect maxocc-moved "$tmp/o2.blind.sum" "maxocc 99 ayz 1281 0 29
maxocc 147 ayz 1581 37 -"
echo 'PASS: paired-max-occ-moved'

check_error() {
    local name=$1 pattern=$2 rc=0
    shift 2
    "$bin" "$@" > "$tmp/$name.out" 2> "$tmp/$name.err" || rc=$?
    [ "$rc" = 1 ] && grep -F -- "$pattern" "$tmp/$name.err" > /dev/null && [ ! -s "$tmp/$name.out" ] || {
        echo "FAIL: $name: exit $rc, expected a refusal with '$pattern'" >&2
        sed -n '1,10p' "$tmp/$name.err" >&2
        exit 1
    }
    echo "PASS: $name"
}
check_error bad-value "--meth-mapq accepts 'conversion-blind' or 'scoring', got 'nope'" aln --meth --meth-mapq nope "$tmp/ref.fa" "$tmp/se.fq"
check_error without-meth '--meth-mapq requires --meth' aln --meth-mapq scoring "$tmp/ref.fa" "$tmp/se.fq"
check_error taps '--meth-mapq conversion-blind applies to em-seq' aln --meth=taps --meth-mapq conversion-blind "$tmp/ref.fa" "$tmp/se.fq"
check_error 5base '--meth-mapq conversion-blind applies to em-seq' aln --meth=5base --meth-mapq conversion-blind "$tmp/ref.fa" "$tmp/se.fq"
check_error neutral '--meth-mapq conversion-blind applies to em-seq' aln --meth --meth-scoring neutral --meth-mapq conversion-blind "$tmp/ref.fa" "$tmp/se.fq"
