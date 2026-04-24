// test/framework/seqpair_gen.cpp
#include "seqpair_gen.h"

#include <cassert>

namespace bwa_tests {

TestPair gen_random_pair(std::mt19937 &rng, int qlen, int rlen) {
    std::uniform_int_distribution<int> base(0, 3);
    TestPair p;
    p.ref.resize(rlen);
    p.qry.resize(qlen);
    p.tag = "random";
    for (int i = 0; i < rlen; i++) p.ref[i] = static_cast<uint8_t>(base(rng));
    for (int i = 0; i < qlen; i++) p.qry[i] = static_cast<uint8_t>(base(rng));
    return p;
}

TestPair gen_exact_match_pair(int len) {
    TestPair p;
    p.ref.assign(len * 2, 0);
    p.qry.assign(len, 0);
    p.tag = "exact_match";
    for (int i = 0; i < len * 2; i++) {
        p.ref[i] = static_cast<uint8_t>((i * 7 + 3) & 3);
    }
    for (int i = 0; i < len; i++) {
        p.qry[i] = p.ref[len / 2 + i];
    }
    return p;
}

TestPair gen_all_mismatch_pair(int len) {
    TestPair p;
    p.ref.assign(len, 0);
    p.qry.assign(len, 1);
    p.tag = "all_mismatch";
    return p;
}

TestPair gen_homopolymer_pair(int len, uint8_t base) {
    TestPair p;
    p.ref.assign(len, base);
    p.qry.assign(len, base);
    p.tag = "homopolymer";
    return p;
}

TestPair gen_sub_cluster_pair(std::mt19937 &rng, int qlen, int rlen,
                              int run_pos, int run_len) {
    assert(run_pos >= 0 && run_len >= 0 && run_pos + run_len <= qlen &&
           "gen_sub_cluster_pair: run_pos + run_len must fit within qlen");
    TestPair p = gen_random_pair(rng, qlen, rlen);
    p.tag = "sub_cluster";
    std::uniform_int_distribution<int> base(0, 3);
    for (int i = 0; i < run_len; i++) {
        p.qry[run_pos + i] = static_cast<uint8_t>(base(rng));
    }
    return p;
}

TestPair gen_with_n_bases_pair(std::mt19937 &rng, int qlen, int rlen,
                               int n_frac_pct) {
    TestPair p = gen_random_pair(rng, qlen, rlen);
    p.tag = "with_n_bases";
    std::uniform_int_distribution<int> pct(0, 99);
    for (int i = 0; i < qlen; i++) {
        if (pct(rng) < n_frac_pct) p.qry[i] = 4;
    }
    for (int i = 0; i < rlen; i++) {
        if (pct(rng) < n_frac_pct) p.ref[i] = 4;
    }
    return p;
}

} // namespace bwa_tests
