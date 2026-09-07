// Regression test: mem_chain2aln_across_reads_V2's srtgg seed-order buffer
// must grow to fit the request in ONE resize, no matter how many doublings
// that takes.
//
// Bug: the growth branch doubled `fac` exactly once:
//
//     if ((spos + c->n) > seeds_per_read_eff * fac * nseq) {
//         fac <<= 1;
//         srtgg = realloc(srtgg, nseq * seeds_per_read_eff * fac * sizeof(uint32_t));
//         ...
//     }
//
// `c->n` (the seed count of a single chain) is not bounded by
// seeds_per_read_eff, so a single doubling can still leave
// (spos + c->n) past the new capacity -- the very next loop,
// `srtgg[spos++] = srt[i]`, then writes past the just-grown allocation.
//
// This test builds ONE chain whose seed count is large enough that a single
// doubling of the initial capacity (nseq * SEEDS_PER_READ * FAC) is not
// enough, so the fixed code must loop the doubling. Every seed spans the
// full read (qbeg == 0 and qbeg + len == l_query), which routes both the
// left- and right-extension branches to their skip paths (no SW/ungapped
// work, no seqPairArray/seqBuf setup needed) -- so the harness only has to
// get the function to the srtgg growth site, not drive a full extension.
//
// Meaningful under AddressSanitizer (`make ASAN=1 srtgg_grow_multi_test`):
// pre-fix, this reports a heap-buffer-overflow WRITE at the
// `srtgg[spos++] = srt[i]` line; post-fix it exits 0.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>

extern "C" {
#include "bntseq.h"
}
#include "bwamem.h"

// _set_pac is already defined by bntseq.h (transitively included above).

int main() {
    // Single-contig reference, long enough to hold the read window with
    // room to spare. Base content is irrelevant -- every seed skips
    // extension, so rseq is fetched but never scored.
    const int64_t l_pac = 200000;
    std::vector<uint8_t> pac((size_t)(l_pac + 3) / 4 + 1, 0);
    for (int64_t i = 0; i < l_pac; ++i) _set_pac(pac.data(), i, (uint8_t)(i & 3));

    bntann1_t anns[1];
    std::memset(anns, 0, sizeof(anns));
    anns[0].offset = 0;
    anns[0].len    = (int32_t)l_pac;
    anns[0].name   = const_cast<char *>("chr_test");
    anns[0].anno   = const_cast<char *>("");

    bntseq_t bns;
    std::memset(&bns, 0, sizeof(bns));
    bns.l_pac  = l_pac;
    bns.n_seqs = 1;
    bns.anns   = anns;
    bns.n_holes = 0;
    bns.ambs    = nullptr;

    mem_opt_t *opt = mem_opt_init();
    if (opt == nullptr) { fprintf(stderr, "mem_opt_init failed\n"); return 2; }

    // One read, one chain, with a single-diagonal chain of seeds so
    // heavily populated that the srtgg buffer's initial capacity
    // (nseq * SEEDS_PER_READ * FAC) needs more than one doubling to fit
    // it. SEEDS_PER_READ=500, FAC=8 -> initial capacity 4000; doubling once
    // gives 8000, so a seed count comfortably past that (10000) forces a
    // second doubling under the fixed code, and overflows the first-doubling
    // buffer under the buggy code.
    const int l_query  = 60;
    const int rbeg      = 100000;
    const int n_seeds   = 10000;

    const int nseq = 1;
    bseq1_t seq;
    std::memset(&seq, 0, sizeof(seq));
    seq.l_seq = l_query;
    std::vector<char> seqbuf((size_t)l_query, 'A');
    seq.seq = seqbuf.data();

    std::vector<mem_seed_t> seeds(n_seeds);
    for (int i = 0; i < n_seeds; ++i) {
        mem_seed_t s;  // default ctor: rbeg=qbeg=len=score=aln=0, n_hits=1
        s.rbeg  = rbeg;
        s.qbeg  = 0;
        s.len   = l_query;             // spans the whole read: skips both
                                        // left- and right-extension branches
        s.score = l_query;
        seeds[i] = s;
    }

    mem_chain_t chain;
    std::memset(&chain, 0, sizeof(chain));
    chain.seqid = 0;
    chain.rid   = 0;
    chain.n     = n_seeds;
    chain.m     = n_seeds;
    chain.first = 0;
    chain.meth_hypothesis = -1;
    chain.frac_rep = 0.0f;
    chain.seeds = seeds.data();

    mem_chain_v chain_v;
    std::memset(&chain_v, 0, sizeof(chain_v));
    chain_v.n = 1;
    chain_v.m = 1;
    chain_v.a = &chain;

    mem_alnreg_v av;
    std::memset(&av, 0, sizeof(av));  // av.a is calloc'd inside the function

    mem_cache mmc;
    std::memset(&mmc, 0, sizeof(mmc));
    // Every mmc pointer/counter for tid=0 starts NULL/zero -- every buffer
    // this call touches (seqPairArray*, seqBuf*) is lazily allocated on
    // first use, and this test's full-length seeds never trigger that.
    //
    // wsize[tid] (*wsize_pair) is the exception: the function's trailing
    // sanity check requires numPairsLeft/Right < *wsize_pair, and with every
    // seed skipping extension numPairsLeft/Right stay 0. A wsize_pair of 0
    // would make that check "0 >= 0" and abort, so give it a nonzero value;
    // seqPairArrayLeft128/Right128/Aux stay NULL, which is fine because
    // nothing dereferences them when numPairsLeft/Right are 0.
    mmc.wsize[0] = 1024;

    int32_t lim_g[BATCH_SIZE + 2];
    std::memset(lim_g, 0, sizeof(lim_g));
    mmc.lim[0] = lim_g;

    // ref_string == NULL routes bns_fetch_seq_v2 through the pac-fetch
    // scratch path (on-demand unpack from `pac`), so no materialized
    // .0123 buffer is needed.
    mem_chain2aln_across_reads_V2(opt, &bns, pac.data(), &seq, nseq,
                                  &chain_v, &av, &mmc, /*ref_string=*/nullptr,
                                  /*tid=*/0);

    if (av.a) free(av.a);
    free(opt);

    fprintf(stderr, "srtgg_grow_multi_test: OK (n_seeds=%d survived without "
                    "overrunning srtgg)\n", n_seeds);
    return 0;
}
