#include "seed_order.h"
#include <cassert>
#include <string.h>
#include <vector>

// Cap for the O(n^2) genomic modes (absorb-count, most-absorb). Above it, both
// fall back to the O(n) global-longest ordering (which performs ~equivalently),
// so neither matrix memory nor O(n^2) CPU blows up on pathological repeat reads.
#define GENOMIC_ORDER_MAX_N 1024

static const char *kNames[] = {
    "off", "global-longest", "local-longest", "absorb-count", "most-absorb"
};

seed_order_t seed_order_from_str(const char *s) {
    if (s)
        for (int i = 0; i <= (int)SEED_ORDER_MOST_ABSORB; ++i)
            if (strcmp(s, kNames[i]) == 0) return (seed_order_t)i;
    return (seed_order_t)-1;
}

const char *seed_order_to_str(seed_order_t m) {
    if ((int)m < 0 || (int)m > (int)SEED_ORDER_MOST_ABSORB) return "?";
    return kNames[(int)m];
}

// Stable counting sort of recs[0..n) by key[i] in [0,maxk]. desc descends
// (via maxk-key, preserving stability). scratch is reused across passes.
static void counting_sort_by(seed_rec_t *recs, int n, std::vector<int> &key,
                             int maxk, bool desc, std::vector<seed_rec_t> &scratch) {
    if (desc) for (int i = 0; i < n; ++i) key[i] = maxk - key[i];
    std::vector<int> cnt(maxk + 2, 0);
    for (int i = 0; i < n; ++i) ++cnt[key[i] + 1];
    for (int k = 1; k <= maxk + 1; ++k) cnt[k] += cnt[k - 1];   // start offsets
    scratch.resize(n);
    for (int i = 0; i < n; ++i) scratch[cnt[key[i]]++] = recs[i];  // stable, ascending
    for (int i = 0; i < n; ++i) recs[i] = scratch[i];
}

static inline bool absorbs(const seed_rec_t &A, const seed_rec_t &B) {
    if (A.rid != B.rid) return false;  // seeds on different chromosomes/contigs can't absorb each other
    const mem_seed_t &a = A.seed, &b = B.seed;
    bool qin = (b.qbeg >= a.qbeg) && (b.qbeg + b.len <= a.qbeg + a.len);
    bool rin = (b.rbeg >= a.rbeg) && (b.rbeg + b.len <= a.rbeg + a.len);
    if (!qin || !rin) return false;
    if (a.len != b.len) return a.len > b.len;            // strictly larger
    return A.orig_ix < B.orig_ix;                        // equal-size: lower ix wins
}

static int max_of(seed_rec_t *recs, int n, bool use_qbeg) {
    int m = 0;
    for (int i = 0; i < n; ++i) { int v = use_qbeg ? recs[i].seed.qbeg : recs[i].seed.len; if (v > m) m = v; }
    return m;
}

void order_seeds(seed_rec_t *recs, int64_t n, seed_order_t mode) {
    if (mode == SEED_ORDER_OFF || n < 2) return;
    assert(n <= INT_MAX);
    int ni = (int)n;
    std::vector<seed_rec_t> scratch;
    std::vector<int> key(ni);
    switch (mode) {
    case SEED_ORDER_GLOBAL_LONGEST: {
        int mk = max_of(recs, ni, false);
        for (int i = 0; i < ni; ++i) key[i] = recs[i].seed.len;
        counting_sort_by(recs, ni, key, mk, /*desc=*/true, scratch);
        break;
    }
    case SEED_ORDER_LOCAL_LONGEST: {
        int mklen = max_of(recs, ni, false), mkq = max_of(recs, ni, true);
        for (int i = 0; i < ni; ++i) key[i] = recs[i].seed.len;          // LSD: secondary first
        counting_sort_by(recs, ni, key, mklen, /*desc=*/true, scratch);
        for (int i = 0; i < ni; ++i) key[i] = recs[i].seed.qbeg;         // then primary qbeg asc
        counting_sort_by(recs, ni, key, mkq, false, scratch);
        break;
    }
    case SEED_ORDER_ABSORB_COUNT: {
        if (ni > GENOMIC_ORDER_MAX_N) {   // O(n^2) count too costly; longest-first is ~equivalent and O(n)
            order_seeds(recs, n, SEED_ORDER_GLOBAL_LONGEST);
            break;
        }
        for (int i = 0; i < ni; ++i) {
            int c = 0;
            for (int j = 0; j < ni; ++j)
                if (j != i && absorbs(recs[i], recs[j])) ++c;
            key[i] = c;
        }
        counting_sort_by(recs, ni, key, ni - 1, /*desc=*/true, scratch);
        break;
    }
    case SEED_ORDER_MOST_ABSORB: {
        if (ni > GENOMIC_ORDER_MAX_N) {        // guard: n x n matrix + O(n^2) too costly; fall back O(n)
            order_seeds(recs, n, SEED_ORDER_GLOBAL_LONGEST);
            break;
        }
        std::vector<unsigned char> M((size_t)ni * ni, 0);
        std::vector<int> rowsum(ni, 0);
        for (int i = 0; i < ni; ++i)
            for (int j = 0; j < ni; ++j)
                if (j != i && absorbs(recs[i], recs[j])) { M[(size_t)i*ni + j] = 1; ++rowsum[i]; }
        std::vector<char> consumed(ni, 0);
        std::vector<seed_rec_t> out; out.reserve(ni);
        for (int picked = 0; picked < ni; ++picked) {
            int best = -1;
            // Phase A assigns orig_ix in append order, so the buffer is in orig_ix order on entry;
            // iterating i=0..ni-1 therefore breaks ties by orig_ix (lower orig_ix wins), as required.
            for (int i = 0; i < ni; ++i) {       // max rowsum, tiebreak orig_ix (i is in ix order)
                if (consumed[i]) continue;
                if (best < 0 || rowsum[i] > rowsum[best]) best = i;
            }
            if (best < 0) break;
            out.push_back(recs[best]);
            consumed[best] = 1;
            for (int j = 0; j < ni; ++j) {       // consume everything best absorbs
                if (consumed[j] || !M[(size_t)best*ni + j]) continue;
                consumed[j] = 1;
                out.push_back(recs[j]);          // absorbed seeds still chained (dropped by chainer)
                for (int i = 0; i < ni; ++i)     // decrement live absorbers of j (column j)
                    if (!consumed[i] && M[(size_t)i*ni + j]) --rowsum[i];
            }
        }
        assert((int)out.size() == ni);
        for (int i = 0; i < ni; ++i) recs[i] = out[i];
        break;
    }
    default: assert(0 && "unhandled seed_order_t in order_seeds"); break;
    }
}
