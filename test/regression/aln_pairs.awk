# Seeded read-pair generator for aln_parity.sh (paired-end section).
# usage: awk -v seed=S -v n=N -v lens=L -v rescue=R -v names=M -v out1=F1 -v out2=F2 -f aln_rng.awk -f aln_pairs.awk ref.fa
#   lens   -- one read length, or MIN-MAX for a uniform range
#   rescue -- percent of pairs whose second mate carries len/12 spread mismatches (too many for
#             the search, within reach of the mate Smith-Waterman)
#   lens2  -- optional read length of the reverse-strand mate of each pair (default: as lens)
#   heavy1 -- optional percent of pairs whose forward mate gets a mismatch every fifth base, so its
#             search runs through its whole entry budget
#   uniform -- 1 draws every pair from a uniformly chosen contig (default: mostly the first)
#   names  -- "plain" (r<i>), "slash" (r<i>/1, r<i>/2) or "comment" (r<i> c1 / r<i> c2)
# Fragments are 250-349 bp from the contigs; 3% of pairs get a second mate from elsewhere and
# 3% a random second mate. rnd, ri, rc and other come from aln_rng.awk (shared with aln_reads.awk).
function edit(t,   e, p, len) {
    len = length(t)
    e = ri(100)
    p = ri(len) + 1
    if (e < 25) t = substr(t, 1, p - 1) other(substr(t, p, 1)) substr(t, p + 1)
    else if (e < 35) {
        t = substr(t, 1, p - 1) other(substr(t, p, 1)) substr(t, p + 1)
        p = ri(len) + 1
        t = substr(t, 1, p - 1) other(substr(t, p, 1)) substr(t, p + 1)
    }
    else if (e < 42) t = substr(t, 1, p - 1) substr(t, p + 1)
    else if (e < 49) t = substr(t, 1, p - 1) substr("ACGT", ri(4) + 1, 1) substr(t, p)
    return t
}
function spread(t,   k, p, step, len) {
    len = length(t)
    step = 12
    for (p = 6; p <= len; p += step) t = substr(t, 1, p - 1) other(substr(t, p, 1)) substr(t, p + 1)
    return t
}
function heavy(t,   p, len) {
    len = length(t)
    for (p = 4; p <= len; p += 5) t = substr(t, 1, p - 1) other(substr(t, p, 1)) substr(t, p + 1)
    return t
}
function randseq(len,   s, i) {
    s = ""
    for (i = 0; i < len; i++) s = s substr("ACGT", ri(4) + 1, 1)
    return s
}
function qual(m,   q, i) { q = ""; for (i = 1; i <= m; i++) q = q "I"; return q }
function pick(   x) {
    if (uniform) return ri(nc) + 1
    x = ri(10)
    if (x == 6) return nc - 2
    if (x == 7) return nc - 1
    if (x == 8) return nc
    return 1
}
/^>/ { nc++; next }
{ seq[nc] = seq[nc] toupper($0) }
END {
    st = seed
    for (i = 0; i < 20; i++) rnd()
    for (r = 0; r < n; r++) {
        if (index(lens, "-")) { split(lens, mm, "-"); len = mm[1] + ri(mm[2] - mm[1] + 1) } else len = lens + 0
        len2 = (lens2 == "") ? len : lens2 + 0
        c = pick(); s = seq[c]; L = length(s)
        f = 250 + ri(100)
        if (f > L) f = L
        if (len > f) len = f
        pos = 1 + ri(L - f + 1)
        frag = substr(s, pos, f)
        fw = substr(frag, 1, len)
        if (len2 > f) len2 = f
        rv = rc(substr(frag, f - len2 + 1, len2))
        a = edit(fw); b = edit(rv)
        if (ri(100) < rescue) b = spread(rv)
        if (heavy1 != "" && ri(100) < heavy1) a = heavy(fw) # too many mismatches for any search budget
        x = ri(100)
        if (x < 3) {
            c2 = pick(); s2 = seq[c2]; L2 = length(s2)
            l2 = (len2 > L2) ? L2 : len2
            b = rc(substr(s2, 1 + ri(L2 - l2 + 1), l2))
        } else if (x < 6) b = randseq(len2)
        if (lens2 == "" && ri(2)) { t = a; a = b; b = t }
        if (names == "slash") { n1 = "r" r "/1"; n2 = "r" r "/2" }
        else if (names == "comment") { n1 = "r" r " c1"; n2 = "r" r " c2" }
        else { n1 = "r" r; n2 = "r" r }
        printf "@%s\n%s\n+\n%s\n", n1, a, qual(length(a)) > out1
        printf "@%s\n%s\n+\n%s\n", n2, b, qual(length(b)) > out2
    }
}
