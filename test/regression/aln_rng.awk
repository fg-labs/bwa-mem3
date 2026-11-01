# Helpers shared by the seeded generators aln_reads.awk and aln_pairs.awk; load it first:
#   awk -f aln_rng.awk -f aln_reads.awk ...
# rnd/ri: Park-Miller generator (16807 mod 2^31-1) over the global `st`, which the caller seeds;
# products stay below 2^46, so the arithmetic is exact in any awk's doubles. rc: reverse
# complement. other(c): a random base other than c.
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
