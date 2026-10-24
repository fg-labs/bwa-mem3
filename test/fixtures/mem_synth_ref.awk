# test/fixtures/mem_synth_ref.awk
#
# Generator (no committed data) for the synthetic reference the `mem` CLI
# regression scripts (test/regression/mem_*.sh) and the matching C++ fixture
# (test/framework/mem_fixture.h) align against. Run with no input:
#
#     awk -f test/fixtures/mem_synth_ref.awk > mem_synth.fa
#
# Bases come from a Park-Miller generator (x = 48271 * x mod 2^31-1, seed
# 20261024) whose products stay below 2^47, so the arithmetic is exact in awk's
# doubles and the output is byte-identical under gawk, mawk and BSD awk. Base
# i of the stream is "ACGT"[int(x / 2^29)] after the i-th step. The contigs are
# drawn from one stream in order, then the features below are written over
# them. Coordinates are 1-based and inclusive, as in SAM.
#
#   ctg1  5000 bp  unique random sequence, except:
#                  2001-2300  block R (300 bp; the original copy)
#                  3001-3500  block S (500 bp; the original copy)
#   ctg2  4000 bp  1001-1300  exact copy of R, same strand
#                  2001-2500  copy of S with a substitution at block offsets
#                             10, 30, 50, ..., 490 (25 sites, 5%): each base
#                             becomes the next one in ACGT order (A>C>G>T>A)
#   ctg3  3000 bp   501-620   (AC)x60 dinucleotide repeat
#                  1501-1540  40 N
#                  2201-2250  poly-A (50 A)
#   ctg4   600 bp   101-400   reverse complement of R
#
# So a read inside R has three equally good hits (ctg1 +, ctg2 +, ctg4 -), a
# read inside S has a best hit on ctg1 and a 5%-diverged second hit on ctg2,
# and every 19-mer outside the planted features is unique in the reference, on
# both strands (checked when this generator was written; a change here must
# recheck it).
#
# The POSIX `cksum` of the output is pinned in test/framework/mem_cli.sh and
# test/framework/mem_fixture.h, so a change to this file, or an awk that
# computes it differently, fails loudly instead of moving every expected value.

function draw(n, s, i) {
    s = ""
    for (i = 0; i < n; i++) {
        x = (x * 48271) % 2147483647
        s = s substr("ACGT", int(x / 536870912) + 1, 1)
    }
    return s
}

function put(s, pos, blk) {
    return substr(s, 1, pos - 1) blk substr(s, pos + length(blk))
}

function revcomp(s, r, i, c) {
    r = ""
    for (i = length(s); i > 0; i--) {
        c = substr(s, i, 1)
        r = r (c == "A" ? "T" : c == "C" ? "G" : c == "G" ? "C" : c == "T" ? "A" : "N")
    }
    return r
}

function diverge(s, r, i, c) {
    r = ""
    for (i = 1; i <= length(s); i++) {
        c = substr(s, i, 1)
        if (i % 20 == 10)
            c = c == "A" ? "C" : c == "C" ? "G" : c == "G" ? "T" : "A"
        r = r c
    }
    return r
}

function run(unit, n, s, i) {
    s = ""
    for (i = 0; i < n; i++)
        s = s unit
    return s
}

function emit(name, s, i) {
    print ">" name
    for (i = 1; i <= length(s); i += 60)
        print substr(s, i, 60)
}

BEGIN {
    x = 20261024
    c1 = draw(5000)
    c2 = draw(4000)
    c3 = draw(3000)
    c4 = draw(600)
    R = substr(c1, 2001, 300)
    S = substr(c1, 3001, 500)
    c2 = put(c2, 1001, R)
    c2 = put(c2, 2001, diverge(S))
    c3 = put(c3, 501, run("AC", 60))
    c3 = put(c3, 1501, run("N", 40))
    c3 = put(c3, 2201, run("A", 50))
    c4 = put(c4, 101, revcomp(R))
    emit("ctg1", c1)
    emit("ctg2", c2)
    emit("ctg3", c3)
    emit("ctg4", c4)
}
