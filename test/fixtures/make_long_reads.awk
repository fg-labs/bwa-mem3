# test/fixtures/make_long_reads.awk
#
# Deterministic long single-end reads with scattered errors, sliced from one
# contig of a FASTA reference. Used by the contained-seed skip regressions to
# drive the skip's second extension batch: long reads are rescored by
# mem_flt_chained_seeds, so container/contained visit order can invert and the
# deferral guard sends a seed to the second batch.
#
# Each read is a window of LMIN..LMAX bases with ~1% substitutions and ~0.1%
# each of 1 bp deletions and 1 bp insertions, drawn from a Park-Miller
# generator. Its arithmetic is exact in awk's doubles, so the output is
# byte-identical across awk implementations.
#
# Usage:
#   awk -v NAME=<contig> -v SEED=<int> -v NREADS=<n> -v LMIN=<bp> -v LMAX=<bp> \
#       -f make_long_reads.awk ref.fa > reads.fq
function rnd() { X = (X * 16807) % 2147483647; return X / 2147483647 }
function rot(b) { return b == "A" ? "C" : b == "C" ? "G" : b == "G" ? "T" : "A" }
/^>/ { cur = substr($1, 2); next }
cur == NAME { seq = seq toupper($0) }
END {
    X = SEED + 0
    for (k = 0; k < NREADS; k++) {
        L = LMIN + int(rnd() * (LMAX - LMIN + 1))
        if (L > length(seq)) L = length(seq)
        off = int(rnd() * (length(seq) - L + 1))
        s = substr(seq, off + 1, L)
        o = ""
        for (i = 1; i <= L; i++) {
            c = substr(s, i, 1)
            u = rnd()
            if (u < 0.010) o = o rot(c)
            else if (u < 0.011) continue
            else if (u < 0.012) o = o c rot(c)
            else o = o c
        }
        q = ""
        for (i = 1; i <= length(o); i++) q = q "I"
        printf "@%s_long_%d\n%s\n+\n%s\n", NAME, k, o, q
    }
}
