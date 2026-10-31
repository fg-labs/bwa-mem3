# Seeded read generator shared by aln_sai_parity.sh and aln_parity.sh.
# usage: awk -v seed=S -v n=N -v lens=L -v tail=T -f aln_reads.awk ref.fa
#   lens  -- comma-separated read lengths, or MIN-MAX for a uniform range
#   tail  -- 1 gives a fifth of the reads a low-quality 30-base tail
# Park-Miller generator (16807 mod 2^31-1): products stay below 2^46, so the
# arithmetic is exact in any awk's doubles.
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
function other(c,   b, k) {
    b = "ACGT"
    do { k = substr(b, ri(4) + 1, 1) } while (k == c)
    return k
}
/^>/ { nc++; name[nc] = substr($1, 2); next }
{ seq[nc] = seq[nc] toupper($0) }
END {
    st = seed
    for (i = 0; i < 20; i++) rnd()
    nl = split(lens, lv, ",")
    for (r = 0; r < n; r++) {
        if (nl == 1 && index(lv[1], "-")) {
            split(lv[1], mm, "-"); len = mm[1] + ri(mm[2] - mm[1] + 1)
        } else len = lv[ri(nl) + 1]
        # contig: a tenth of the reads come from the tandem repeat, another tenth
        # from the reverse-complement copy, the rest uniformly over contigs
        x = ri(10)
        if (x == 0) c = nc; else if (x == 1) c = nc - 1; else c = ri(nc) + 1
        s = seq[c]
        L = length(s)
        if (c == nc - 1 && x == 1) { lo = 1; span = 300 } else { lo = 1; span = L }
        if (len > span) len = span
        pos = lo + ri(span - len + 1)
        t = substr(s, pos, len)
        if (ri(2)) t = rc(t)
        # edits
        e = ri(100)
        p = ri(len > 0 ? len : 1) + 1
        if (e < 25) t = substr(t, 1, p - 1) other(substr(t, p, 1)) substr(t, p + 1)
        else if (e < 35) {
            t = substr(t, 1, p - 1) other(substr(t, p, 1)) substr(t, p + 1)
            p = ri(len) + 1
            t = substr(t, 1, p - 1) other(substr(t, p, 1)) substr(t, p + 1)
        }
        else if (e < 50) t = substr(t, 1, p - 1) substr(t, p + 1)
        else if (e < 60) t = substr(t, 1, p - 1) substr("ACGT", ri(4) + 1, 1) substr(t, p)
        else if (e < 65) t = substr(t, 1, p - 1) "N" substr(t, p + 1)
        n2 = length(t)
        q = ""
        for (i = 1; i <= n2; i++) q = q ((tail && r % 5 == 0 && i > n2 - 30) ? "#" : "I")
        printf "@r%d\n%s\n+\n%s\n", r, t, q
    }
}
