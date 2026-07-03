#include "seed_order.h"
#include <cassert>
#include <climits>
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

namespace {

/* Reusable per-thread scratch. order_seeds runs once per read on each worker
 * thread, so the buffers below would otherwise be malloc'd/freed millions of
 * times over a run. thread_local keeps one set per worker, grown to the
 * high-water mark and reused; nothing here is shared across threads. */
struct OrderScratch {
    std::vector<int>           perm;      // index permutation being sorted (4 B/elem)
    std::vector<int>           tmp;       // stable-scatter target for perm
    std::vector<int>           key;       // per-seed sort key, indexed by orig position
    std::vector<int>           cnt;       // counting-sort histogram
    std::vector<seed_rec_t>    applied;   // final materialization of the permutation
    std::vector<unsigned char> mat;       // genomic modes: n x n absorb matrix
    std::vector<int>           rowsum;    // genomic modes: per-seed absorb count
    std::vector<char>          consumed;  // most-absorb: greedy consumed flags
};
thread_local OrderScratch g;

/* Stable counting sort of the index permutation g.perm[0..n) by g.key[perm[i]]
 * in [0,maxk]. desc descends via (maxk-key) while preserving stability. Only
 * 4-byte indices move here; order_seeds applies the final permutation to recs
 * once, so a full seed_rec_t is copied a single time regardless of pass count. */
void sort_perm(int n, int maxk, bool desc) {
    int *perm = g.perm.data();
    int *key  = g.key.data();
    if (desc) for (int i = 0; i < n; ++i) key[i] = maxk - key[i];
    g.cnt.assign(maxk + 2, 0);
    int *cnt = g.cnt.data();
    for (int i = 0; i < n; ++i) ++cnt[key[perm[i]] + 1];
    for (int k = 1; k <= maxk + 1; ++k) cnt[k] += cnt[k - 1];   // start offsets
    g.tmp.resize(n);
    int *tmp = g.tmp.data();
    for (int i = 0; i < n; ++i) tmp[cnt[key[perm[i]]]++] = perm[i];  // stable, ascending
    for (int i = 0; i < n; ++i) perm[i] = tmp[i];
}

inline bool absorbs(const seed_rec_t &A, const seed_rec_t &B) {
    if (A.rid != B.rid) return false;  // seeds on different chromosomes/contigs can't absorb each other
    const mem_seed_t &a = A.seed, &b = B.seed;
    bool qin = (b.qbeg >= a.qbeg) && (b.qbeg + b.len <= a.qbeg + a.len);
    bool rin = (b.rbeg >= a.rbeg) && (b.rbeg + b.len <= a.rbeg + a.len);
    if (!qin || !rin) return false;
    if (a.len != b.len) return a.len > b.len;            // strictly larger
    return A.orig_ix < B.orig_ix;                        // equal-size: lower ix wins
}

int max_of(const seed_rec_t *recs, int n, bool use_qbeg) {
    int m = 0;
    for (int i = 0; i < n; ++i) { int v = use_qbeg ? recs[i].seed.qbeg : recs[i].seed.len; if (v > m) m = v; }
    return m;
}

} // namespace

void order_seeds(seed_rec_t *recs, int64_t n, seed_order_t mode) {
    if (mode == SEED_ORDER_OFF || n < 2) return;
    assert(n <= INT_MAX);
    int ni = (int)n;

    // perm starts as the identity => recs' current (orig_ix) order, which every
    // stable pass below preserves as the final tiebreak.
    g.perm.resize(ni);
    for (int i = 0; i < ni; ++i) g.perm[i] = i;
    g.key.resize(ni);

    switch (mode) {
    case SEED_ORDER_GLOBAL_LONGEST: {
        int mk = max_of(recs, ni, false);
        for (int i = 0; i < ni; ++i) g.key[i] = recs[i].seed.len;
        sort_perm(ni, mk, /*desc=*/true);
        break;
    }
    case SEED_ORDER_LOCAL_LONGEST: {
        int mklen = max_of(recs, ni, false), mkq = max_of(recs, ni, true);
        for (int i = 0; i < ni; ++i) g.key[i] = recs[i].seed.len;          // LSD: secondary first
        sort_perm(ni, mklen, /*desc=*/true);
        for (int i = 0; i < ni; ++i) g.key[i] = recs[i].seed.qbeg;         // then primary qbeg asc
        sort_perm(ni, mkq, /*desc=*/false);
        break;
    }
    case SEED_ORDER_ABSORB_COUNT: {
        if (ni > GENOMIC_ORDER_MAX_N) {   // O(n^2) count too costly; longest-first is ~equivalent and O(n)
            order_seeds(recs, n, SEED_ORDER_GLOBAL_LONGEST);
            return;
        }
        for (int i = 0; i < ni; ++i) {
            int c = 0;
            for (int j = 0; j < ni; ++j)
                if (j != i && absorbs(recs[i], recs[j])) ++c;
            g.key[i] = c;
        }
        sort_perm(ni, ni - 1, /*desc=*/true);
        break;
    }
    case SEED_ORDER_MOST_ABSORB: {
        if (ni > GENOMIC_ORDER_MAX_N) {        // guard: n x n matrix + O(n^2) too costly; fall back O(n)
            order_seeds(recs, n, SEED_ORDER_GLOBAL_LONGEST);
            return;
        }
        g.mat.assign((size_t)ni * ni, 0);
        g.rowsum.assign(ni, 0);
        unsigned char *M = g.mat.data();
        int *rowsum = g.rowsum.data();
        for (int i = 0; i < ni; ++i)
            for (int j = 0; j < ni; ++j)
                if (j != i && absorbs(recs[i], recs[j])) { M[(size_t)i*ni + j] = 1; ++rowsum[i]; }
        g.consumed.assign(ni, 0);
        char *consumed = g.consumed.data();
        int np = 0;  // greedy pick order is written straight into perm
        for (int picked = 0; picked < ni; ++picked) {
            int best = -1;
            // recs is in orig_ix order, so iterating i=0..ni-1 breaks rowsum ties
            // by orig_ix (lower wins), as required.
            for (int i = 0; i < ni; ++i) {
                if (consumed[i]) continue;
                if (best < 0 || rowsum[i] > rowsum[best]) best = i;
            }
            if (best < 0) break;
            g.perm[np++] = best;
            consumed[best] = 1;
            // No need to decrement column `best` (live absorbers of best): at
            // selection time no live seed can absorb best. Such a seed would
            // absorb best plus everything best absorbs (containment + length
            // order are transitive), giving rowsum >= rowsum[best]+1 > rowsum[best],
            // so it would have been picked as best instead. rowsum thus stays an
            // exact live-count for every still-live seed without this decrement.
            for (int j = 0; j < ni; ++j) {       // consume everything best absorbs
                if (consumed[j] || !M[(size_t)best*ni + j]) continue;
                consumed[j] = 1;
                g.perm[np++] = j;                // absorbed seeds still chained (dropped by chainer)
                for (int i = 0; i < ni; ++i)     // decrement live absorbers of j (column j)
                    if (!consumed[i] && M[(size_t)i*ni + j]) --rowsum[i];
            }
        }
        assert(np == ni);
        break;
    }
    default: assert(0 && "unhandled seed_order_t in order_seeds"); return;
    }

    // Apply the permutation to recs exactly once — the only full seed_rec_t move.
    g.applied.resize(ni);
    for (int i = 0; i < ni; ++i) g.applied[i] = recs[g.perm[i]];
    for (int i = 0; i < ni; ++i) recs[i] = g.applied[i];
}
