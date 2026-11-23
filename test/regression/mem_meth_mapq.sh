#!/usr/bin/env bash
# test/regression/mem_meth_mapq.sh
#
# `bwa-mem3 mem --meth --meth-mapq`: the conversion-blind MAPQ cap, on a seeded synthetic reference
# whose contigs each hold two copies of a 40-base segment that differ only where one copy has the
# converted base's source (C, or G for the G->A strand) and the other already has the converted
# base (T, or A):
#
#   ab  copies A (pos 401) and B (pos 541) differ at two C/T positions (offsets 10 and 25)
#   ef  copies E (pos 401) and F (pos 841) differ at one C/T position (offset 20)
#   gh  copies G (pos 401) and H (pos 841) differ at one G/A position (offset 20)
#   ij  copies I and J differ at one C/T position (offset 12) and by one base J has after offset 35,
#       so the competing alignment is gapped; kl holds such a pair (C/T at offset 30, the extra
#       base after offset 10) reverse-complemented
#   xy, pq, mn, ed/ex/ey  a mirror mismatch past where the competitor's score peaks (no gain), a
#       competitor with an insertion and a deletion over equal spans (none), one clipped before its
#       retained Cs (gains only by extending), and one split across a contig boundary
#   rs  adjacent copies U1 and U2 differ at three C/T positions, for a mate found only by rescue
#   zt26, zt40  as rs, but U2 overlaps U1 and ends inside (26) or outside (40) the rescue zone
#
# A read 1 that kept such a C (a methylated C) and converted the rest is unique under genomic
# scoring, because at the other copy the C sits over a T (the mirror cell) and costs a mismatch.
# With the mirror cell free the copies score the same, so the conversion-blind MAPQ is 0 and YQ:i
# keeps the genomic MAPQ. Collapsed scoring at genomic's -B 4 (the mirror cell free, same
# penalties) is the independent oracle for the single-end cap.
#
# Checks: the default equals an explicit conversion-blind; single-end reads (one and two retained
# Cs on the C->T strand, one retained G on the G->A strand, and a read whose competitor needs a
# gap, on either strand, capped part-way) are capped to collapsed -B 4's MAPQ with YQ:i, and
# controls (a fully converted ambiguous read, a unique read) are not; only MAPQ and
# YQ:i differ between the modes, and YQ:i is the scoring MAPQ of exactly the lowered records; the default for collapsed, neutral, TAPS and 5-base is scoring,
# and collapsed accepts conversion-blind as a no-op; paired end, in a batch with an insert-size
# estimate: only MAPQ and YQ:i (and the MAPQ MQ:i and SA:Z repeat) differ, no MAPQ is raised, a
# pair whose MAPQ read 1's retained C raised is capped, unpaired ends take their single-end cap, a
# rescued mate is capped by its rescue window's second-best alignment (outside the best's zone), and
# every MAPQ equals collapsed -B 4's; --meth-mapq-max-len caps only reads within it (0: every
# length), at the boundary inclusive, and a whole pair when one end is; and the options' refusals.
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

# The reference and reads: a Park-Miller generator (16807 mod 2^31-1, exact in awk's doubles).
cat > "$tmp/gen.awk" << 'EOF'
function rnd() { st = (st * 16807) % 2147483647; return st }
function ri(n) { return rnd() % n }
function rc(s,   i, r, c) {
    r = ""
    for (i = length(s); i >= 1; i--) {
        c = substr(s, i, 1)
        r = r (c == "A" ? "T" : c == "C" ? "G" : c == "G" ? "C" : c == "T" ? "A" : c)
    }
    return r
}
function dna(n,   s, i) { s = ""; for (i = 0; i < n; i++) s = s substr("ACGT", ri(4) + 1, 1); return s }
function put(s, i, c) { return substr(s, 1, i) c substr(s, i + 2) }
function other(c,   k) { do { k = substr("ACGT", ri(4) + 1, 1) } while (k == c); return k }
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
# a read 2 (G->A) from the reverse strand of contig s at 0-based start p, length n
function mate(s, p, n) { return conv(rc(substr(s, p + 1, n)), "G", "A", "") }
BEGIN {
    st = 20261121
    a = put(put(dna(40), 10, "C"), 25, "C"); b = put(put(a, 10, "T"), 25, "T")
    ab = dna(400) a dna(100) b dna(400)
    e = put(dna(40), 20, "C"); f = put(e, 20, "T")
    ef = dna(400) e dna(400) f dna(400)
    g = put(dna(40), 20, "G"); h = put(g, 20, "A")
    gh = dna(400) g dna(400) h dna(400)
    # gapped competitors: copy J is copy I with a T at I's C (offset 12) and one more base after
    # offset 35, so a read from I aligns to J only with a deletion; kl holds a second such pair
    # (the T at offset 30, the extra base after offset 10, before it) reverse-complemented, so its
    # read aligns to the reverse strand and its mirror mismatch lies past the gap
    i = put(dna(50), 12, "C"); j = put(i, 12, "T"); j = substr(j, 1, 35) "A" substr(j, 36)
    ij = dna(300) i dna(300) j dna(300)
    k = put(dna(50), 30, "C"); l = put(k, 30, "T"); l = substr(l, 1, 10) "A" substr(l, 11)
    kl = rc(dna(300) k dna(300) l dna(300))
    # tail: copy Y is X with real mismatches at offsets 10 and 46 and a T at X's C 47, past where
    # Y's alignment score peaks; freeing that cell does not raise Y's score
    x = put(dna(50), 47, "C"); y = put(put(put(x, 10, other(substr(x, 11, 1))), 46, other(substr(x, 47, 1))), 47, "T")
    xy = dna(300) x dna(300) y dna(300)
    # indel: copy B is A without the base at offset 15 and with an extra base after offset 40, so its
    # alignment holds an insertion and a deletion over equal spans; A keeps Cs between the two
    p2 = put(put(put(dna(60), 20, "C"), 25, "C"), 30, "C"); q2 = substr(p2, 1, 15) substr(p2, 17, 25) "G" substr(p2, 42)
    pq = dna(300) p2 dna(300) q2 dna(300)
    # clipped: copy N is M with a real mismatch at offset 30 and Ts at M's Cs 33..45 after it, so N's
    # genomic alignment stops at 30, and only freed does it run on
    mm = dna(50)
    for (t = 33; t <= 45; t += 3) mm = put(mm, t, "C")
    nn = put(mm, 30, other(substr(mm, 31, 1)))
    for (t = 33; t <= 45; t += 3) nn = put(nn, t, "T")
    mn = dna(300) mm dna(300) nn dna(300)
    # edge: copy D (D0 with a T at offset 20) is split across a contig boundary: its first 6 bases
    # end contig ex, the rest starts contig ey; an alignment cannot run across it
    d0 = put(dna(40), 20, "C"); d1 = put(d0, 20, "T")
    printf ">ab\n%s\n>ef\n%s\n>gh\n%s\n>ij\n%s\n>kl\n%s\n", ab, ef, gh, ij, kl > (dir "/ref.fa")
    printf ">xy\n%s\n>pq\n%s\n>mn\n%s\n", xy, pq, mn > (dir "/ref.fa")
    printf ">ed\n%s\n>ex\n%s\n>ey\n%s\n", dna(300) d0 dna(300), dna(300) substr(d1, 1, 6), substr(d1, 7) dna(300) > (dir "/ref.fa")
    # gclip: as mn, but copy N2 also has one more base after offset 3, so its genomic alignment is
    # gapped and stops at the real mismatch at 30; only freed does it run on past its span
    m2 = dna(50)
    for (t = 33; t <= 45; t += 3) m2 = put(m2, t, "C")
    n2 = put(m2, 30, other(substr(m2, 31, 1)))
    for (t = 33; t <= 45; t += 3) n2 = put(n2, t, "T")
    n2 = substr(n2, 1, 3) "A" substr(n2, 4)
    printf ">gc\n%s\n", dna(300) m2 dna(300) n2 dna(300) > (dir "/ref.fa")
    # lo: copies L1 and L2 of a 70-base segment differ at one C/T position (offset 30); its reads
    # are longer than --meth-mapq-max-len 60
    l1 = put(dna(70), 30, "C"); l2 = put(l1, 30, "T")
    lo = dna(300) l1 dna(300) l2 dna(300)
    printf ">lo\n%s\n", lo > (dir "/ref.fa")
    fq(dir "/long.fq", "long1", conv(l1, "C", "T", " 30 "))
    # a mixed-length pair: a unique 40-base read 1 and a 70-base read 2 from L1's reverse strand
    # that kept L1's C (a G at offset 39 of read 2), so only read 1 is within a limit of 60
    fq(dir "/lp1.fq", "mix", conv(substr(lo, 101, 40), "C", "T", ""))
    fq(dir "/lp2.fq", "mix", conv(rc(l1), "G", "A", " 39 "))
    # chim: ef1's read followed by 80 converted bases found nowhere in the reference (no C, so
    # freeing the mirror cell does not extend into them), so its longest aligned run is under
    # --chimera-qc's 44% of the read
    fq(dir "/chim.fq", "chim", conv(e, "C", "T", " 20 ") conv(dna(80), "C", "T", ""))
    fq(dir "/se.fq", "ab2", conv(a, "C", "T", " 10 25 "))
    fq(dir "/se.fq", "ab1", conv(a, "C", "T", " 10 "))
    fq(dir "/se.fq", "ef1", conv(e, "C", "T", " 20 "))
    fq(dir "/se.fq", "gh1", conv(rc(g), "C", "T", " 19 "))
    fq(dir "/se.fq", "ef0", conv(e, "C", "T", ""))
    fq(dir "/se.fq", "uniq", conv(substr(ef, 101, 40), "C", "T", ""))
    fq(dir "/se.fq", "ij1", conv(i, "C", "T", " 12 "))
    fq(dir "/se.fq", "kl1", conv(k, "C", "T", " 30 "))
    fq(dir "/se.fq", "tail", conv(x, "C", "T", " 47 "))
    fq(dir "/se.fq", "indel", conv(p2, "C", "T", " 20 25 30 "))
    fq(dir "/se.fq", "clipped", conv(mm, "C", "T", " 33 36 39 42 45 "))
    fq(dir "/se.fq", "edge", conv(d0, "C", "T", " 20 "))
    fq(dir "/se.fq", "gclip", conv(m2, "C", "T", " 33 36 39 42 45 "))
    # paired end, one batch on ef (enough pairs for an insert-size estimate): read 1 keeps copy E's
    # C, read 2 starts at p in the flank after copy E
    r1 = conv(e, "C", "T", " 20 ")
    for (p = 460; p <= 780; p += 10) { fq(dir "/p1.fq", "m" p, r1); fq(dir "/p2.fq", "m" p, mate(ef, p, 60)) }
    # rescue: copies U1 and U2 of a 60-base unit sit side by side, both within the insert-size window
    # of a read 1 from the unique flank before them, and differ only at three C/T positions (U1's Cs
    # at offsets 10, 30 and 50). Read 2 comes from U1's reverse strand with a real mismatch every 13
    # bases, so it has no seed and is found only by mate rescue, with U2 as the second-best
    # alignment in the rescue window (its csub). r1kept keeps U1's three Cs (Gs in read 2): U2 is
    # then three mirror mismatches below U1, and tied with it once they are free; r1conv converts
    # them, a tie under either scoring.
    u1 = put(put(put(dna(60), 10, "C"), 30, "C"), 50, "C")
    u2 = put(put(put(u1, 10, "T"), 30, "T"), 50, "T")
    rs = dna(300) u1 u2 dna(300)
    printf ">rs\n%s\n", rs > (dir "/ref.fa")
    for (t = 0; t < 2; ++t) {
        m = conv(rc(u1), "G", "A", t ? "" : " 9 29 49 ")
        for (o = 4; o < 60; o += 13) { c = substr(m, o + 1, 1); m = put(m, o, c == "A" ? "C" : "A") }
        fq(dir "/p1.fq", t ? "r1conv" : "r1kept", conv(substr(rs, 111, 40), "C", "T", ""))
        fq(dir "/p2.fq", t ? "r1conv" : "r1kept", m)
    }
    # unique pairs from gh's first flank, so the batch has pairs to estimate the insert size from
    # (an end with a close second hit, such as every read 1 above, is left out of the estimate)
    for (p = 0; p <= 290; p += 10) {
        fq(dir "/p1.fq", "u" p, conv(substr(gh, p + 1, 40), "C", "T", ""))
        fq(dir "/p2.fq", "u" p, mate(gh, p + 180 + (p / 10 % 7) * 10, 60))
    }
    # rescue zone: as rs, but U2 overlaps U1. W repeats a P-base unit to 60 + P bases; U1 = W's
    # first 60 bases, with Cs at offsets 5, 15 and 25 where the unit has Ts, and U2 = W from P on,
    # the same bases with those Ts. Mate rescue's second-best skips alignments ending within
    # score / A columns of the best (the zone, 36 here), and U2 ends P columns after U1. At P = 26
    # U2 is inside the zone and the second-best is a shorter alignment to the repeat; at P = 40 it
    # is U2. Generated last, so the fixtures above keep their bases.
    for (P = 26; P <= 40; P += 14) {
        v = put(put(put(dna(P), 5, "T"), 15, "T"), 25, "T")
        w = ""
        while (length(w) < 60 + P) w = w v
        w = put(put(put(substr(w, 1, 60 + P), 5, "C"), 15, "C"), 25, "C")
        zt = dna(300) w dna(300)
        printf ">zt%d\n%s\n", P, zt > (dir "/ref.fa")
        for (t = 0; t < 2; ++t) {
            m = conv(rc(substr(w, 1, 60)), "G", "A", t ? "" : " 34 44 54 ")
            for (o = 4; o < 60; o += 13) { c = substr(m, o + 1, 1); m = put(m, o, c == "A" ? "C" : "A") }
            fq(dir "/p1.fq", "z" P (t ? "conv" : "kept"), conv(substr(zt, 111, 40), "C", "T", ""))
            fq(dir "/p2.fq", "z" P (t ? "conv" : "kept"), m)
        }
    }
}
EOF
awk -v dir="$tmp" -f "$tmp/gen.awk" < /dev/null
"$bin" index --meth "$tmp/ref.fa" > "$tmp/index.out" 2> "$tmp/index.err" || {
    echo "FAIL: index --meth exited nonzero" >&2
    cat "$tmp/index.err" >&2
    exit 1
}

records() { grep -v '^@' "$1"; }
# name FLAG RNAME POS MAPQ YQ (or -) of each record, one per line
summary() { awk '{ yq = "-"; for (i = 12; i <= NF; ++i) if (index($i, "YQ:i:") == 1) yq = substr($i, 6); print $1, $2, $3, $4, $5, yq }' "$1"; }
# the records without MAPQ, YQ:i and the MAPQ that MQ:i and SA:Z repeat (MQ:i itself is kept)
strip() {
    awk 'BEGIN { OFS = "\t" } {
        $5 = "."; out = $1
        for (i = 2; i <= NF; ++i) {
            if (index($i, "YQ:i:") == 1) continue
            if (index($i, "MQ:i:") == 1) $i = "MQ:i:."
            if (index($i, "SA:Z:") == 1) { # rname,pos,strand,CIGAR,MAPQ,NM; per entry
                n = split(substr($i, 6), e, ";"); v = "SA:Z:"
                for (k = 1; k < n; ++k) { split(e[k], f, ","); f[5] = "."; v = v f[1] "," f[2] "," f[3] "," f[4] "," f[5] "," f[6] ";" }
                $i = v
            }
            out = out "\t" $i
        }
        print out
    }' "$1"
}
run() {
    local out=$1
    shift
    # the EXIT trap removes $tmp: show the aligner's stderr before it does
    "$bin" mem -k 14 -T 18 "$@" > "$tmp/$out.sam" 2> "$tmp/$out.err" || {
        echo "FAIL: $out: mem exited nonzero" >&2
        cat "$tmp/$out.err" >&2
        exit 1
    }
    records "$tmp/$out.sam" > "$tmp/$out.rec"
    summary "$tmp/$out.rec" > "$tmp/$out.sum"
}
# each record of a conversion-blind run ($1) against the same record of the scoring run ($2), both
# summaries: YQ:i is the scoring MAPQ when the cap lowered the MAPQ, and absent otherwise
yq_is_scoring() {
    paste -d ' ' "$1" "$2" | awk '$1 != $7 || $2 != $8 || ($5 < $11 && $6 != $11) || ($5 >= $11 && $6 != "-") { bad = 1; print "  " $0 > "/dev/stderr" } END { exit bad }' || {
        echo "FAIL: $3: YQ:i is not the scoring MAPQ of exactly the lowered records" >&2
        exit 1
    }
}
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

# Single end.
run se.blind --meth "$tmp/ref.fa" "$tmp/se.fq"
run se.explicit --meth --meth-mapq conversion-blind "$tmp/ref.fa" "$tmp/se.fq"
run se.scoring --meth --meth-mapq scoring "$tmp/ref.fa" "$tmp/se.fq"
run se.oracle --meth --meth-scoring collapsed -B 4 "$tmp/ref.fa" "$tmp/se.fq"
cmp "$tmp/se.blind.rec" "$tmp/se.explicit.rec" || {
    echo "FAIL: default: differs from an explicit --meth-mapq conversion-blind" >&2
    exit 1
}
echo 'PASS: default-is-conversion-blind'

expect se-scoring "$tmp/se.scoring.sum" "ab2 0 ab 401 60 -
ab1 0 ab 401 27 -
ef1 0 ef 401 27 -
gh1 16 gh 401 27 -
ef0 0 ef 401 0 -
uniq 0 ef 101 60 -
ij1 0 ij 301 42 -
kl1 16 kl 652 42 -
tail 0 xy 301 32 -
indel 0 pq 301 48 -
clipped 0 mn 301 60 -
edge 0 ed 301 60 -
gclip 0 gc 301 60 -"
expect se-conversion-blind "$tmp/se.blind.sum" "ab2 0 ab 401 0 60
ab1 0 ab 401 0 27
ef1 0 ef 401 0 27
gh1 16 gh 401 0 27
ef0 0 ef 401 0 -
uniq 0 ef 101 60 -
ij1 0 ij 301 22 42
kl1 16 kl 652 22 42
tail 0 xy 301 32 -
indel 0 pq 301 48 -
clipped 0 mn 301 15 60
edge 0 ed 301 33 60
gclip 0 gc 301 28 60"
# the oracle: collapsed -B 4's MAPQ (placement can differ, at an equal-scoring copy, so each row is
# checked to be the same read, and for a pair the same read end, but not the same placement)
paste -d ' ' "$tmp/se.blind.sum" "$tmp/se.oracle.sum" | awk 'NF != 12 || $1 != $7 || int($2 / 64) % 4 != int($8 / 64) % 4 || $5 != $11 { bad = 1; print "  " $0 > "/dev/stderr" } END { exit bad }' || {
    echo "FAIL: single end: a row is not the same read in both runs, or the conversion-blind MAPQ differs from collapsed -B 4's" >&2
    exit 1
}
echo 'PASS: single-end-capped-to-collapsed-B4'

cmp <(strip "$tmp/se.blind.rec") <(strip "$tmp/se.scoring.rec") || {
    echo "FAIL: modes differ outside MAPQ and YQ:i" >&2
    exit 1
}
yq_is_scoring "$tmp/se.blind.sum" "$tmp/se.scoring.sum" 'single end'
echo 'PASS: only-mapq-and-yq-differ'

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
run other --meth --meth-scoring collapsed "$tmp/ref.fa" "$tmp/se.fq"
run other.blind --meth --meth-scoring collapsed --meth-mapq conversion-blind "$tmp/ref.fa" "$tmp/se.fq"
cmp "$tmp/other.rec" "$tmp/other.blind.rec" || {
    echo "FAIL: collapsed: conversion-blind is not a no-op" >&2
    exit 1
}
echo 'PASS: scoring-outside-emseq-genomic'

# Paired end.
run pe.blind --meth "$tmp/ref.fa" "$tmp/p1.fq" "$tmp/p2.fq"
run pe.scoring --meth --meth-mapq scoring "$tmp/ref.fa" "$tmp/p1.fq" "$tmp/p2.fq"
grep -q 'mean and std.dev' "$tmp/pe.blind.err" || {
    echo "FAIL: fixture: the paired batch gave no insert-size estimate" >&2
    exit 1
}
cmp <(strip "$tmp/pe.blind.rec") <(strip "$tmp/pe.scoring.rec") || {
    echo "FAIL: paired: modes differ outside MAPQ and YQ:i" >&2
    exit 1
}
paste -d ' ' "$tmp/pe.blind.sum" "$tmp/pe.scoring.sum" | awk 'NF != 12 || $1 != $7 || int($2 / 64) % 4 != int($8 / 64) % 4 || $5 > $11 { bad = 1; print "  " $0 > "/dev/stderr" } END { exit bad }' || {
    echo "FAIL: paired: a row is not the same read end in both runs, or conversion-blind raised a MAPQ" >&2
    exit 1
}
yq_is_scoring "$tmp/pe.blind.sum" "$tmp/pe.scoring.sum" paired
echo 'PASS: paired-only-lowers'
pick() { awk -v n="$2" '$1 == n' "$1"; }
# m600: a proper pair. Read 1's single-end MAPQ is capped to 0, and mem's pairing then gives it at
# most 40 above that (min(q_pe, q_se + 40)), as collapsed -B 4 does; read 2 keeps its own 60.
expect paired-capped <(pick "$tmp/pe.blind.sum" m600) "m600 99 ef 401 40 60
m600 147 ef 601 60 -"
# m460: read 2 is too close for the insert-size estimate, so the ends are not paired and read 1
# takes its single-end cap.
expect paired-unpaired <(pick "$tmp/pe.blind.sum" m460) "m460 97 ef 401 0 27
m460 145 ef 461 60 -"
# r1kept: read 2, found by mate rescue, is capped to 0 by its rescue window's second-best alignment
# (U2), which ties U1 once the mirror cell is free, as r1conv's does under either scoring
{
    pick "$tmp/pe.blind.sum" r1kept
    pick "$tmp/pe.blind.sum" r1conv
} > "$tmp/rescue.sum"
expect paired-rescue-csub "$tmp/rescue.sum" "r1kept 99 rs 111 60 -
r1kept 147 rs 305 0 60
r1conv 99 rs 111 60 -
r1conv 147 rs 305 0 -"
# the rescue zone: at P = 26 U2 is inside it, so z26kept's read 2 is capped only to the MAPQ the
# shorter alignment leaves once its mirror cell is free, as z26conv's genomic MAPQ; at P = 40 U2 is
# outside it, and ties U1 once the mirror cell is free
{
    pick "$tmp/pe.blind.sum" z26kept
    pick "$tmp/pe.blind.sum" z26conv
    pick "$tmp/pe.blind.sum" z40kept
    pick "$tmp/pe.blind.sum" z40conv
} > "$tmp/zone.sum"
expect paired-rescue-zone "$tmp/zone.sum" "z26kept 99 zt26 111 60 -
z26kept 147 zt26 305 48 60
z26conv 99 zt26 111 60 -
z26conv 147 zt26 305 48 -
z40kept 99 zt40 111 60 -
z40kept 147 zt40 305 0 60
z40conv 99 zt40 111 60 -
z40conv 147 zt40 305 0 -"
# in both modes, the records come as mates in order (read 1 and read 2 of one name; an odd count
# fails), and each carries MQ:i, which repeats its mate's emitted (capped) MAPQ
for mode in blind scoring; do
    awk '{ mq[NR] = "none"; for (i = 12; i <= NF; ++i) if (index($i, "MQ:i:") == 1) mq[NR] = substr($i, 6); q[NR] = $5
          name[NR] = $1; end[NR] = int($2 / 64) % 4 } # read end: 1 (read 1) or 2 (read 2)
        END { if (NR % 2) { print "  an unpaired trailing record" > "/dev/stderr"; bad = 1 }
              for (r = 1; r < NR; r += 2)
                  if (name[r] != name[r + 1] || end[r] * end[r + 1] != 2 || mq[r] != q[r + 1] || mq[r + 1] != q[r]) { print "  pair " (r + 1) / 2 > "/dev/stderr"; bad = 1 }
              exit bad }' "$tmp/pe.$mode.rec" || {
        echo "FAIL: paired ($mode): the records are not mates in order, or MQ:i is missing or does not repeat the mate's emitted MAPQ" >&2
        exit 1
    }
done
run pe.oracle --meth --meth-scoring collapsed -B 4 "$tmp/ref.fa" "$tmp/p1.fq" "$tmp/p2.fq"
paste -d ' ' "$tmp/pe.blind.sum" "$tmp/pe.oracle.sum" | awk 'NF != 12 || $1 != $7 || int($2 / 64) % 4 != int($8 / 64) % 4 || $5 != $11 { bad = 1; print "  " $0 > "/dev/stderr" } END { exit bad }' || {
    echo "FAIL: paired: a row is not the same read end in both runs, or the conversion-blind MAPQ differs from collapsed -B 4's" >&2
    exit 1
}
echo 'PASS: paired-capped'

# --meth-mapq-max-len: only reads of at most that length are capped (a pair: both ends, when either
# end is); 0, the default, caps every length. Every se.fq and p1/p2 read is at most 60 bp. Paired
# read 1 is 40 bp and read 2 60 bp, so at 50 only read 1 is within it and the whole pair is still
# capped.
run pe.len50 --meth --meth-mapq-max-len 50 "$tmp/ref.fa" "$tmp/p1.fq" "$tmp/p2.fq"
cmp "$tmp/pe.len50.rec" "$tmp/pe.blind.rec" || {
    echo "FAIL: paired: --meth-mapq-max-len between the ends' lengths does not cap the whole pair" >&2
    exit 1
}
run se.len0 --meth --meth-mapq-max-len 0 "$tmp/ref.fa" "$tmp/se.fq"
run pe.len0 --meth --meth-mapq-max-len 0 "$tmp/ref.fa" "$tmp/p1.fq" "$tmp/p2.fq"
run se.len10 --meth --meth-mapq-max-len 10 "$tmp/ref.fa" "$tmp/se.fq"
run pe.len10 --meth --meth-mapq-max-len 10 "$tmp/ref.fa" "$tmp/p1.fq" "$tmp/p2.fq"
for end in se pe; do
    cmp "$tmp/$end.len0.rec" "$tmp/$end.blind.rec" || {
        echo "FAIL: $end: --meth-mapq-max-len 0 differs from the default" >&2
        exit 1
    }
    cmp "$tmp/$end.len10.rec" "$tmp/$end.scoring.rec" || {
        echo "FAIL: $end: --meth-mapq-max-len below every read's length still caps" >&2
        exit 1
    }
done
# at the boundary: ab2 (capped 60 -> 0) is capped at its own length, and not one base below it
ab2_len=$(awk 'NR % 4 == 1 && $1 == "@ab2" { getline; print length($0) }' "$tmp/se.fq")
run se.at --meth --meth-mapq-max-len "$ab2_len" "$tmp/ref.fa" "$tmp/se.fq"
run se.below --meth --meth-mapq-max-len "$((ab2_len - 1))" "$tmp/ref.fa" "$tmp/se.fq"
expect len-at <(pick "$tmp/se.at.sum" ab2) "ab2 0 ab 401 0 60"
expect len-below <(pick "$tmp/se.below.sum" ab2) "ab2 0 ab 401 60 -"
# a 70-base read is capped by default, and not past a limit of 60
run long.blind --meth "$tmp/ref.fa" "$tmp/long.fq"
run long.len60 --meth --meth-mapq-max-len 60 "$tmp/ref.fa" "$tmp/long.fq"
expect len-long-default "$tmp/long.blind.sum" "long1 0 lo 301 0 12"
expect len-long-60 "$tmp/long.len60.sum" "long1 0 lo 301 12 -"
# a pair is capped as a whole when its shorter end is within the limit: read 2 (70 bases) is
# lowered at 60 because read 1 (40 bases) is within it, and not below read 1's length
run mix.len60 --meth --meth-mapq-max-len 60 "$tmp/ref.fa" "$tmp/lp1.fq" "$tmp/lp2.fq"
run mix.short --meth --meth-mapq-max-len 39 "$tmp/ref.fa" "$tmp/lp1.fq" "$tmp/lp2.fq"
expect len-mixed-60 <(pick "$tmp/mix.len60.sum" mix | sed -n 2p) "mix 145 lo 301 0 12"
expect len-mixed-below <(pick "$tmp/mix.short.sum" mix | sed -n 2p) "mix 145 lo 301 12 -"
echo 'PASS: max-len'

# --chimera-qc lowers a flagged record's MAPQ to at most 1 in both modes; YQ:i is the uncapped MAPQ
# under the same ceiling, so it appears only when the cap lowered the MAPQ below that.
run chim.scoring --meth --chimera-qc --meth-mapq scoring "$tmp/ref.fa" "$tmp/chim.fq"
run chim.blind --meth --chimera-qc --meth-mapq-max-len 0 "$tmp/ref.fa" "$tmp/chim.fq"
expect chimera-qc-scoring "$tmp/chim.scoring.sum" "chim 512 ef 401 1 -"
expect chimera-qc-blind "$tmp/chim.blind.sum" "chim 512 ef 401 1 -"
yq_is_scoring "$tmp/chim.blind.sum" "$tmp/chim.scoring.sum" chimera-qc
# and in every run above, a YQ:i is above its record's MAPQ
for f in "$tmp"/*.sum; do
    awk '$6 != "-" && $6 + 0 <= $5 + 0 { print "  " FILENAME ": " $0 > "/dev/stderr"; bad = 1 } END { exit bad }' "$f" || {
        echo "FAIL: a YQ:i is not above its record's MAPQ" >&2
        exit 1
    }
done
echo 'PASS: chimera-qc-yq'

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
check_error bad-value "--meth-mapq accepts 'conversion-blind' or 'scoring', got 'nope'" mem --meth --meth-mapq nope "$tmp/ref.fa" "$tmp/se.fq"
check_error without-meth '--meth-mapq requires --meth' mem --meth-mapq scoring "$tmp/ref.fa" "$tmp/se.fq"
check_error bad-max-len "--meth-mapq-max-len accepts a read length >= 0 (0: every length), got '-1'" mem --meth --meth-mapq-max-len -1 "$tmp/ref.fa" "$tmp/se.fq"
check_error max-len-without-meth '--meth-mapq-max-len requires --meth' mem --meth-mapq-max-len 50 "$tmp/ref.fa" "$tmp/se.fq"
check_error taps '--meth-mapq conversion-blind applies to em-seq' mem --meth=taps --meth-mapq conversion-blind "$tmp/ref.fa" "$tmp/se.fq"
check_error 5base '--meth-mapq conversion-blind applies to em-seq' mem --meth=5base --meth-mapq conversion-blind "$tmp/ref.fa" "$tmp/se.fq"
check_error neutral '--meth-mapq conversion-blind applies to em-seq' mem --meth --meth-scoring neutral --meth-mapq conversion-blind "$tmp/ref.fa" "$tmp/se.fq"
