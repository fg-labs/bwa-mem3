// test/unit/test_mem_seed_oracle.cpp
//
// The FM index and the three seeding rounds of `bwa-mem3 mem`, checked against
// a brute-force suffix array of the text the index searches. The oracle is
// written from the model in the bwa papers and the comments of
// src/FMI_search.h, not from the code under test:
//
//   - The indexed text is the forward strand of every contig as the index
//     stores it (the 2-bit .pac, an N replaced by some base), followed by its
//     reverse complement, with a sentinel `$` that sorts before every base. Row
//     r of the BWT is the r-th suffix of that text in sorted order, so row 0 is
//     the suffix `$`. count[b] (after the loader's +1) is the row where the
//     suffixes starting with base b begin, and the SA interval (k, s) of a
//     string is the run of rows whose suffixes start with it.
//   - Round 1 (getSMEMsAllPosOneThread, min_intv 1) emits every maximal exact
//     match of length >= k: an interval [m, n] of the read that occurs in the
//     text and can be extended neither left nor right without losing every
//     occurrence. With a threshold t in place of 1 (round 2's min_intv = s + 1,
//     getSMEMsOnePosOneThread_lockstep), a call at position x emits the
//     intervals through x that are maximal under "occurs at least t times".
//     The brute force finds, for every start i, the last end R(i) that still
//     occurs t times; [i, R(i)] is maximal when R(i - 1) < R(i). Every case
//     uses a minimum length of at least 2: the forward phase records the
//     interval it extends before testing the extension against t, starting
//     with the lone base at x, so with a minimum length of 1 a base occurring
//     fewer than t times can be emitted, which the maximal-interval model does
//     not describe.
//   - Round 3 (bwtSeedStrategyAllPosOneThread and its lockstep twin) walks x
//     from 0: it extends [x, j] to the right until fewer than max_intv
//     occurrences remain and the length is at least the floor, emits the
//     interval when it still occurs, and restarts after it (also when it no
//     longer occurs: the walk then goes on to the length floor, emits nothing
//     and restarts there); an N restarts the walk past it, and running off the
//     read ends it.
//   - mem_collect_smem runs the three rounds (round 2 inside each round-1 SMEM
//     at least split_len long and occurring at most -s times, at its middle
//     position), sorts the result by read and then by (m, n), and with
//     --smem-dedup drops adjacent duplicates. split_len is k * r rounded to
//     the nearest integer with a half rounded down (bwa's -r); each case below
//     states its value, worked out by hand.
//
// Every SMEM is compared on (rid, m, n, k, s). The reverse-complement interval
// start `l` is excluded: the backward phase leaves it stale by design
// (backwardExt_konly in FMI_search.h) and nothing downstream reads it. The
// order in which a seeding round writes its SMEMs is not asserted, only the
// set: mem_collect_smem sorts them by (rid, m, n) before anyone reads them,
// and two SMEMs equal on (rid, m, n) are equal on (k, s). mem_collect_smem's
// own output is compared in order.
//
// The lockstep drivers interleave reads in slots, so every round is run at
// every lockstep width from 1 to the compile-time maximum; the SA-sample
// resolver is checked on indexes built at every SA stride `index -u` accepts
// (1 to 64 rows), under each --ks-dedup mode and both sentinel arms. The
// reference is a small random text with planted exact, diverged, reverse-
// complement and tandem repeats, an N run and three contigs; the reads are
// slices of it with substitutions, indels, N bases, chimeras and random
// sequence, of lengths 0 to 150, and two identical reads side by side (so
// --smem-dedup meets equal SMEMs of different reads).

#include "doctest/doctest.h"
#include "mem_fixture.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "FMI_search.h"
#include "bntseq.h"
#include "bwamem.h"
#include "lockstep_width.h"

// Defined in src/bwamem.cpp, where mem_kernel1_core is its only caller; there
// is no header declaration to include.
SMEM *mem_collect_smem(FMI_search *fmi, const mem_opt_t *opt, const bseq1_t *seq_, int nseq,
                       SMEM *matchArray, int32_t *min_intv_ar, int32_t *query_pos_ar,
                       uint8_t *enc_qdb, int32_t *rid, mem_cache *mmc, int64_t &tot_smem, int tid);

namespace {

using bwa_tests::Contig;

// ---------------------------------------------------------------------------
// The reference and the reads.

char next_base(char c) { return c == 'A' ? 'C' : c == 'C' ? 'G' : c == 'G' ? 'T' : 'A'; }

std::string random_bases(std::mt19937 &rng, size_t n) {
    std::string s(n, 'A');
    for (size_t i = 0; i < n; ++i) s[i] = "ACGT"[rng() % 4];
    return s;
}

// 0-based, half-open coordinates:
//   ctgA 900 bp: X = [100, 220) copied exactly to [600, 720); (AC)x15 at [300, 330).
//   ctgB 700 bp: reverse complement of X at [50, 170); Y = ctgA[400, 500) at
//                [250, 350) with every 12th base changed (11-bp exact windows);
//                12 N at [400, 412); 25 A at [500, 525); Z = ctgA[800, 830) at
//                [560, 590), [600, 630), [640, 670), and its reverse
//                complement at [680, 710).
//   ctgC  60 bp: random, so a read can span the ctgB/ctgC junction.
// So in the indexed text (both strands together) X occurs three times, Z five
// times, and Y once exactly plus once diverged.
std::vector<Contig> oracle_reference() {
    std::mt19937 rng(20261024u);
    std::string a = random_bases(rng, 900), b = random_bases(rng, 700), c = random_bases(rng, 60);
    const std::string X = a.substr(100, 120);
    a.replace(600, X.size(), X);
    for (int i = 0; i < 15; ++i) a.replace((size_t)(300 + 2 * i), 2, "AC");
    b.replace(50, X.size(), bwa_tests::revcomp(X));
    std::string Y = a.substr(400, 100);
    for (size_t i = 11; i < Y.size(); i += 12) Y[i] = next_base(Y[i]);
    b.replace(250, Y.size(), Y);
    b.replace(400, 12, std::string(12, 'N'));
    b.replace(500, 25, std::string(25, 'A'));
    const std::string Z = a.substr(800, 30);
    b.replace(560, 30, Z);
    b.replace(600, 30, Z);
    b.replace(640, 30, Z);
    b.replace(680, 30, bwa_tests::revcomp(Z));
    return {{"ctgA", a}, {"ctgB", b}, {"ctgC", c}};
}

// Reads as text (A C G T N). F is the contigs concatenated, so a slice of it
// can span a contig junction, as a seed in the indexed text can.
std::vector<std::string> oracle_reads(const std::vector<Contig> &ref) {
    std::mt19937 rng(7u);
    std::string F;
    for (const Contig &c : ref) F += c.seq;
    const size_t nF = F.size();
    auto slice = [&](size_t len) { return F.substr(rng() % (nF - len + 1), len); };
    auto mutate = [&](std::string s, int n_sub) -> std::string {
        for (int i = 0; i < n_sub; ++i) {
            const size_t p = rng() % s.size();
            s[p] = next_base(s[p] == 'N' ? 'A' : s[p]);
        }
        return s;
    };
    std::vector<std::string> reads;
    reads.push_back(F.substr(0, 40));               // the first bases of the text (SA value 0)
    reads.push_back(F.substr(nF - 30));             // the last forward bases
    reads.push_back("");                            // an empty read
    for (int i = 0; i < 20; ++i) reads.push_back(slice(30 + rng() % 121));
    for (int i = 0; i < 15; ++i) reads.push_back(bwa_tests::revcomp(slice(30 + rng() % 121)));
    for (int i = 0; i < 20; ++i) {
        // One draw per statement: argument evaluation order is unspecified.
        const std::string s = slice(40 + rng() % 111);
        reads.push_back(mutate(s, 1 + (int)(rng() % 4)));
    }
    for (int i = 0; i < 10; ++i) {                  // an insertion or a deletion of 1-3 bases
        std::string s = slice(60 + rng() % 91);
        const size_t p = 10 + rng() % (s.size() - 20), k = 1 + rng() % 3;
        if (i % 2 == 0) s.erase(p, k);
        else s.insert(p, random_bases(rng, k));
        reads.push_back(i % 4 == 1 ? bwa_tests::revcomp(s) : s);
    }
    for (int i = 0; i < 8; ++i) {                   // chimeras of ctgA and ctgB
        const std::string left = ref[0].seq.substr(rng() % 800, 50);
        reads.push_back(left + ref[1].seq.substr(rng() % 600, 50));
    }
    for (int i = 0; i < 10; ++i) {                  // N bases inside, first or last
        std::string s = slice(30 + rng() % 91);
        const int n = 1 + (int)(rng() % 3);
        for (int j = 0; j < n; ++j) s[rng() % s.size()] = 'N';
        if (i == 0) s[0] = 'N';
        if (i == 1) s[s.size() - 1] = 'N';
        reads.push_back(s);
    }
    reads.push_back(std::string(20, 'N'));
    for (int i = 0; i < 6; ++i) reads.push_back(random_bases(rng, 30 + rng() % 91));
    const int short_lens[] = {1, 2, 3, 5, 8, 13, 18, 19, 20, 24};
    for (int len : short_lens) reads.push_back(slice((size_t)len));
    const std::string &A = ref[0].seq, &B = ref[1].seq;
    reads.push_back(A.substr(120, 80));             // inside X
    reads.push_back(A.substr(80, 60));              // across X's left edge
    reads.push_back(A.substr(690, 60));             // across X's second copy's right edge
    reads.push_back(B.substr(40, 140));             // the reverse-complement copy of X, with flanks
    reads.push_back(A.substr(400, 100));            // Y, which has a diverged copy
    reads.push_back(A.substr(430, 50));
    reads.push_back(A.substr(790, 50));             // Z with flanks
    reads.push_back(B.substr(555, 120));            // three copies of Z
    reads.push_back(A.substr(295, 40));             // the (AC)x15 repeat
    reads.push_back(B.substr(495, 35));             // the poly-A run
    reads.push_back(B.substr(390, 40));             // across the N run
    reads.push_back(B.substr(650, 50) + ref[2].seq.substr(0, 40));  // across the ctgB/ctgC junction
    reads.push_back(mutate(A.substr(120, 90), 2));
    reads.push_back(mutate(B.substr(560, 110), 3));
    // Two identical unique 20-mers: each has one SMEM (offsets 0-19, found by
    // rounds 1 and 3), so after sorting, the first read's last SMEM and the
    // second read's first SMEM are equal in everything but the read.
    reads.push_back(A.substr(30, 20));
    reads.push_back(A.substr(30, 20));
    return reads;
}

// ---------------------------------------------------------------------------
// The brute-force suffix array.

struct SuffixLess {
    const std::vector<uint8_t> *t;
    int64_t n;
    bool operator()(int64_t a, int64_t b) const {
        if (a == b) return false;
        for (;;) {
            const int ca = a < n ? (*t)[(size_t)a] : -1, cb = b < n ? (*t)[(size_t)b] : -1;
            if (ca != cb) return ca < cb;
            ++a, ++b;
        }
    }
};

struct Brute {
    std::vector<uint8_t> fwd;   // forward strand, 2-bit codes, as the .pac holds it
    std::vector<uint8_t> t;     // fwd followed by its reverse complement; `$` follows
    int64_t n = 0;              // t.size()
    std::vector<int64_t> sa;    // n + 1 rows, row 0 = the suffix `$` (position n)
    std::vector<int64_t> rank;  // rank[p] = row of the suffix at position p, p in [0, n]

    // Character d of row's suffix: 0-3, or -1 for `$`.
    int at(int64_t row, int64_t d) const {
        const int64_t p = sa[(size_t)row] + d;
        return p < n ? t[(size_t)p] : -1;
    }

    // [lo, hi) holds the rows whose suffixes start with some string of length
    // d; narrow it to those whose next character is c.
    void narrow(int64_t &lo, int64_t &hi, int64_t d, int c) const {
        int64_t a = lo, b = hi;
        while (a < b) {
            const int64_t m = a + (b - a) / 2;
            if (at(m, d) < c) a = m + 1; else b = m;
        }
        int64_t a2 = a, b2 = hi;
        while (a2 < b2) {
            const int64_t m = a2 + (b2 - a2) / 2;
            if (at(m, d) <= c) a2 = m + 1; else b2 = m;
        }
        lo = a;
        hi = a2;
    }

    // SA interval (k = first row, s = row count) of q[0, len).
    void interval(const uint8_t *q, int len, int64_t &k, int64_t &s) const {
        int64_t lo = 0, hi = n + 1;
        for (int d = 0; d < len && lo < hi; ++d) narrow(lo, hi, d, q[d]);
        k = lo;
        s = hi - lo;
    }
};

Brute brute_of(const uint8_t *pac, int64_t l_pac) {
    Brute b;
    b.fwd.resize((size_t)l_pac);
    for (int64_t i = 0; i < l_pac; ++i) b.fwd[(size_t)i] = (uint8_t)_get_pac(pac, i);
    b.t = b.fwd;
    for (int64_t i = l_pac - 1; i >= 0; --i) b.t.push_back((uint8_t)(3 - b.fwd[(size_t)i]));
    b.n = (int64_t)b.t.size();
    b.sa.resize((size_t)b.n + 1);
    for (int64_t i = 0; i <= b.n; ++i) b.sa[(size_t)i] = i;
    SuffixLess less = {&b.t, b.n};
    std::sort(b.sa.begin(), b.sa.end(), less);
    b.rank.resize((size_t)b.n + 1);
    for (int64_t r = 0; r <= b.n; ++r) b.rank[(size_t)b.sa[(size_t)r]] = r;
    return b;
}

// ---------------------------------------------------------------------------
// Fixtures shared by the cases: the reference, one index per SA stride, the
// brute force, and the reads (built once; the index build is the costly part).

struct Fixture {
    std::vector<Contig> ref;
    std::unique_ptr<bwa_tests::TempIndex> ix[7];   // sa_compx 0..6 (SA stride 1..64)
    Brute bf;
    std::vector<std::vector<uint8_t>> reads;       // 0-4 codes
    // maximal_intervals(reads[rid], t), computed once per (rid, t).
    std::map<std::pair<size_t, int64_t>, std::vector<std::pair<int, int>>> maximal;
};

Fixture &fixture() {
    static std::unique_ptr<Fixture> f;
    if (!f) {
        f.reset(new Fixture);
        f->ref = oracle_reference();
        f->ix[3].reset(new bwa_tests::TempIndex(f->ref, 3));
        f->bf = brute_of(f->ix[3]->pac(), f->ix[3]->bns()->l_pac);
        for (const std::string &r : oracle_reads(f->ref)) f->reads.push_back(bwa_tests::nt4(r));
    }
    return *f;
}

FMI_search *index_at(int sa_compx) {
    Fixture &f = fixture();
    if (!f.ix[sa_compx]) f.ix[sa_compx].reset(new bwa_tests::TempIndex(f.ref, sa_compx));
    return f.ix[sa_compx]->fmi();
}

// ---------------------------------------------------------------------------
// SMEMs and the three oracles.

struct Smem {
    uint32_t rid, m, n;
    int64_t k, s;
    bool operator<(const Smem &o) const {
        if (rid != o.rid) return rid < o.rid;
        if (m != o.m) return m < o.m;
        if (n != o.n) return n < o.n;
        if (k != o.k) return k < o.k;
        return s < o.s;
    }
    bool operator==(const Smem &o) const {
        return rid == o.rid && m == o.m && n == o.n && k == o.k && s == o.s;
    }
};

std::string show(const std::vector<Smem> &v) {
    std::ostringstream os;
    for (const Smem &x : v) os << "(" << x.rid << "," << x.m << "," << x.n << "," << x.k << "," << x.s << ") ";
    return os.str();
}

// "" when equal, else the first differing position (or the length mismatch).
std::string first_diff(const std::vector<Smem> &got, const std::vector<Smem> &want) {
    if (got == want) return "";
    size_t i = 0;
    while (i < got.size() && i < want.size() && got[i] == want[i]) ++i;
    std::ostringstream os;
    os << "sizes " << got.size() << " vs " << want.size() << "; first difference at " << i << ": got "
       << (i < got.size() ? show(std::vector<Smem>(1, got[i])) : "nothing") << " want "
       << (i < want.size() ? show(std::vector<Smem>(1, want[i])) : "nothing");
    return os.str();
}

std::vector<Smem> observable(const SMEM *a, int64_t n, bool sort_them) {
    std::vector<Smem> v;
    for (int64_t i = 0; i < n; ++i) v.push_back({a[i].rid, a[i].m, a[i].n, a[i].k, a[i].s});
    if (sort_them) std::sort(v.begin(), v.end());
    return v;
}

Smem make_smem(const Brute &bf, uint32_t rid, const std::vector<uint8_t> &q, int m, int e) {
    Smem x = {rid, (uint32_t)m, (uint32_t)e, 0, 0};
    bf.interval(q.data() + m, e - m + 1, x.k, x.s);
    return x;
}

// For each start i, the last end R[i] such that q[i..R[i]] occurs at least t
// times (i - 1 when even q[i] does not); then the maximal intervals.
std::vector<std::pair<int, int>> maximal_intervals(const Brute &bf, const std::vector<uint8_t> &q, int64_t t) {
    const int len = (int)q.size();
    std::vector<int> R((size_t)len);
    for (int i = 0; i < len; ++i) {
        int e = i - 1;
        int64_t lo = 0, hi = bf.n + 1;
        for (int j = i; j < len && q[(size_t)j] < 4; ++j) {
            bf.narrow(lo, hi, j - i, q[(size_t)j]);
            if (hi - lo < t) break;
            e = j;
        }
        R[(size_t)i] = e;
    }
    std::vector<std::pair<int, int>> out;
    for (int i = 0; i < len; ++i)
        if (R[(size_t)i] >= i && (i == 0 || R[(size_t)i - 1] < R[(size_t)i])) out.push_back({i, R[(size_t)i]});
    return out;
}

const std::vector<std::pair<int, int>> &maximal_of(Fixture &f, size_t rid, int64_t t) {
    const std::pair<size_t, int64_t> key(rid, t);
    std::map<std::pair<size_t, int64_t>, std::vector<std::pair<int, int>>>::iterator it = f.maximal.find(key);
    if (it == f.maximal.end()) it = f.maximal.insert(std::make_pair(key, maximal_intervals(f.bf, f.reads[rid], t))).first;
    return it->second;
}

// Round 1 for read rid: every maximal exact match of length >= min_len.
std::vector<Smem> oracle_round1(Fixture &f, size_t rid, int min_len) {
    std::vector<Smem> out;
    for (const std::pair<int, int> &iv : maximal_of(f, rid, 1))
        if (iv.second - iv.first + 1 >= min_len)
            out.push_back(make_smem(f.bf, (uint32_t)rid, f.reads[rid], iv.first, iv.second));
    return out;
}

// Round 2 on read rid at position x with threshold t: the maximal intervals through x.
std::vector<Smem> oracle_round2(Fixture &f, size_t rid, int x, int64_t t, int min_len) {
    std::vector<Smem> out;
    for (const std::pair<int, int> &iv : maximal_of(f, rid, t))
        if (iv.first <= x && x <= iv.second && iv.second - iv.first + 1 >= min_len)
            out.push_back(make_smem(f.bf, (uint32_t)rid, f.reads[rid], iv.first, iv.second));
    return out;
}

// Round 3 for one read: the forward walk described in the header.
std::vector<Smem> oracle_round3(const Brute &bf, uint32_t rid, const std::vector<uint8_t> &q, int64_t max_intv,
                                int min_len) {
    std::vector<Smem> out;
    const int len = (int)q.size();
    int x = 0;
    while (x < len) {
        if (q[(size_t)x] > 3) { ++x; continue; }
        int64_t lo = 0, hi = bf.n + 1;
        bf.narrow(lo, hi, 0, q[(size_t)x]);
        int next = len;
        for (int j = x + 1; j < len; ++j) {
            if (q[(size_t)j] > 3) { next = j + 1; break; }
            if (lo < hi) bf.narrow(lo, hi, j - x, q[(size_t)j]);
            if (hi - lo < max_intv && j - x + 1 >= min_len) {
                if (hi > lo) out.push_back({rid, (uint32_t)x, (uint32_t)j, lo, hi - lo});
                next = j + 1;
                break;
            }
        }
        x = next;
    }
    return out;
}

// ---------------------------------------------------------------------------
// Running the code under test.

// A batch of reads in the layout the seeding functions take: the 0-4 codes
// concatenated (enc_qdb), each read's offset into it, and bseq1_t records whose
// `seq` holds the same codes (mem_collect_smem copies them from there).
struct Batch {
    std::vector<std::vector<char>> codes;
    std::vector<bseq1_t> seqs;
    std::vector<uint8_t> enc;
    std::vector<int32_t> cum;
    int32_t max_len = 0;

    explicit Batch(const std::vector<std::vector<uint8_t>> &reads) {
        codes.resize(reads.size());
        seqs.resize(reads.size());
        for (size_t i = 0; i < reads.size(); ++i) {
            codes[i].assign(reads[i].begin(), reads[i].end());
            codes[i].push_back(0);
            std::memset(&seqs[i], 0, sizeof(bseq1_t));
            seqs[i].l_seq = (int)reads[i].size();
            seqs[i].seq = codes[i].data();
            cum.push_back((int32_t)enc.size());
            enc.insert(enc.end(), reads[i].begin(), reads[i].end());
            max_len = std::max(max_len, (int32_t)reads[i].size());
        }
        enc.push_back(4);   // never read; keeps data() valid for an all-empty batch
    }
    int32_t size() const { return (int32_t)seqs.size(); }
};

// Restores the lockstep widths and the third-round switch on scope exit.
struct LockstepGuard {
    int32_t smem_n = g_smem_lockstep_n, bwt_n = g_bwtseed_lockstep_n, bwt_on = g_bwtseed_lockstep;
    ~LockstepGuard() {
        g_smem_lockstep_n = smem_n;
        g_bwtseed_lockstep_n = bwt_n;
        g_bwtseed_lockstep = bwt_on;
    }
};

// Lockstep widths to run: every width 1..max when `every` (a lane-mixing fault
// can show at one width only), else a spread of small, odd, power-of-two and
// maximal ones.
std::vector<int32_t> widths(int32_t max_width, bool every) {
    std::vector<int32_t> w;
    if (every) {
        for (int32_t i = 1; i <= max_width; ++i) w.push_back(i);
    } else {
        const int32_t some[] = {1, 2, 3, 5, 8, 13, 16, 24, 31, 32, 33, 47, 63, 64};
        for (int32_t x : some)
            if (x <= max_width) w.push_back(x);
    }
    return w;
}

std::vector<Smem> run_round1(FMI_search *fmi, Batch &b, int min_len) {
    const int32_t nr = b.size();
    std::vector<int32_t> min_intv((size_t)nr, 1), rid((size_t)nr);
    for (int32_t i = 0; i < nr; ++i) rid[(size_t)i] = i;
    std::vector<SMEM> out((size_t)nr * (size_t)b.max_len + 1);
    int64_t n_out = 0;
    fmi->getSMEMsAllPosOneThread(b.enc.data(), min_intv.data(), rid.data(), nr, nr, b.seqs.data(), b.cum.data(),
                                 b.max_len, min_len, out.data(), &n_out);
    REQUIRE(n_out <= (int64_t)out.size());
    return observable(out.data(), n_out, true);
}

struct Probe {
    int32_t rid, x, t;
};

std::vector<Smem> run_round2(FMI_search *fmi, Batch &b, const std::vector<Probe> &probes, int min_len) {
    const int32_t np = (int32_t)probes.size();
    std::vector<int32_t> qpos, min_intv, rid;
    for (const Probe &p : probes) {
        qpos.push_back(p.x);
        min_intv.push_back(p.t);
        rid.push_back(p.rid);
    }
    std::vector<SMEM> out((size_t)np * (size_t)b.max_len + 1);
    int64_t n_out = 0;
    fmi->getSMEMsOnePosOneThread_lockstep(b.enc.data(), qpos.data(), min_intv.data(), rid.data(), np, np,
                                          b.seqs.data(), b.cum.data(), b.max_len, min_len, out.data(), &n_out);
    REQUIRE(n_out <= (int64_t)out.size());
    return observable(out.data(), n_out, true);
}

std::vector<Smem> run_round3(FMI_search *fmi, Batch &b, int32_t max_intv, int min_len, bool lockstep) {
    const int32_t nr = b.size();
    std::vector<int32_t> mi((size_t)nr, max_intv);
    std::vector<SMEM> out((size_t)nr * (size_t)b.max_len + 1);
    const int64_t n_out =
        lockstep ? fmi->bwtSeedStrategyAllPosOneThread_lockstep(b.enc.data(), mi.data(), nr, b.seqs.data(),
                                                                 b.cum.data(), min_len, out.data(), b.max_len)
                 : fmi->bwtSeedStrategyAllPosOneThread(b.enc.data(), mi.data(), nr, b.seqs.data(), b.cum.data(),
                                                       min_len, out.data());
    REQUIRE(n_out <= (int64_t)out.size());
    return observable(out.data(), n_out, true);
}

// mem_collect_smem with the buffers sized the way mem_kernel1_core sizes them
// (tot_len entries), so an overrun is an overrun here too (and an AddressSanitizer
// build reports it). The result is the returned array; freeing the per-thread
// cache follows the current layout of mem_cache.
std::vector<Smem> run_collect(FMI_search *fmi, Batch &b, const mem_opt_t *opt) {
    int64_t tot_len = 0;
    for (const bseq1_t &s : b.seqs) tot_len += s.l_seq;
    REQUIRE(tot_len > 0);
    std::unique_ptr<mem_cache, void (*)(void *)> mmc((mem_cache *)calloc(1, sizeof(mem_cache)), free);
    REQUIRE(mmc != nullptr);
    mmc->wsize_mem[0] = tot_len;
    mmc->matchArray[0] = (SMEM *)_mm_malloc((size_t)tot_len * sizeof(SMEM), 64);
    REQUIRE(mmc->matchArray[0] != nullptr);
    std::vector<int32_t> min_intv((size_t)tot_len), qpos((size_t)tot_len), rid((size_t)tot_len);
    std::vector<uint8_t> enc((size_t)tot_len);
    int64_t tot = 0;
    SMEM *res = mem_collect_smem(fmi, opt, b.seqs.data(), b.size(), mmc->matchArray[0], min_intv.data(),
                                 qpos.data(), enc.data(), rid.data(), mmc.get(), tot, 0);
    REQUIRE(res != nullptr);
    std::vector<Smem> v = observable(res, tot, false);
    _mm_free(mmc->matchArray[0]);
    _mm_free(mmc->smem_sort_scratch[0].cnt);
    _mm_free(mmc->smem_sort_scratch[0].tmp);
    _mm_free(mmc->smem_sort_scratch[0].pairs);
    return v;
}

// The oracle for mem_collect_smem: rounds 1-3 as the header describes, sorted
// by (rid, m, n) (ties are identical), adjacent duplicates removed under
// --smem-dedup.
// n_round[] counts what each round contributed, so a case can show that the
// round it varies ran.
std::vector<Smem> oracle_collect(Fixture &f, const mem_opt_t *opt, int split_len, size_t n_round[3]) {
    std::vector<Smem> all;
    n_round[0] = n_round[1] = n_round[2] = 0;
    for (size_t r = 0; r < f.reads.size(); ++r) {
        const std::vector<Smem> r1 = oracle_round1(f, r, opt->min_seed_len);
        all.insert(all.end(), r1.begin(), r1.end());
        n_round[0] += r1.size();
        for (const Smem &p : r1) {
            if ((int)(p.n + 1 - p.m) < split_len || p.s > opt->split_width) continue;
            const std::vector<Smem> r2 = oracle_round2(f, r, (int)((p.m + p.n + 1) >> 1), p.s + 1, opt->min_seed_len);
            all.insert(all.end(), r2.begin(), r2.end());
            n_round[1] += r2.size();
        }
        if (opt->max_mem_intv > 0) {
            const std::vector<Smem> r3 =
                oracle_round3(f.bf, (uint32_t)r, f.reads[r], (int64_t)opt->max_mem_intv, opt->min_seed_len + 1);
            all.insert(all.end(), r3.begin(), r3.end());
            n_round[2] += r3.size();
        }
    }
    std::sort(all.begin(), all.end());
    if (opt->smem_dedup) all.erase(std::unique(all.begin(), all.end()), all.end());
    return all;
}

int64_t count_base(const Brute &bf, int c) {
    int64_t n = 0;
    for (uint8_t x : bf.t) n += (x == c);
    return n;
}

}  // namespace

// ---------------------------------------------------------------------------

TEST_CASE("mem seeding oracle: the index's text is the FASTA on both strands"
          * doctest::test_suite("unit/mem_seed")) {
    Fixture &f = fixture();
    std::string F;
    for (const Contig &c : f.ref) F += c.seq;
    REQUIRE(f.bf.fwd.size() == F.size());
    int64_t n_bases = 0;
    for (size_t i = 0; i < F.size(); ++i) {
        if (F[i] == 'N') continue;   // an N is stored as some base
        CAPTURE(i);
        REQUIRE(f.bf.fwd[i] == bwa_tests::nt4(F.substr(i, 1))[0]);
        ++n_bases;
    }
    CHECK(n_bases == (int64_t)F.size() - 12);
    // The planted repeats are where the reference comment says (guards the oracle's inputs).
    const std::string &A = f.ref[0].seq, &B = f.ref[1].seq;
    CHECK(A.substr(100, 120) == A.substr(600, 120));
    CHECK(B.substr(50, 120) == bwa_tests::revcomp(A.substr(100, 120)));
    CHECK(B.substr(560, 30) == A.substr(800, 30));
    CHECK(B.substr(400, 12) == std::string(12, 'N'));
    // And occur as often as the reference comment says, in the indexed text.
    int64_t k, n_occ;
    const std::vector<uint8_t> x = bwa_tests::nt4(A.substr(130, 25)), z = bwa_tests::nt4(A.substr(805, 20)),
                               y = bwa_tests::nt4(A.substr(400, 100)), w = bwa_tests::nt4(A.substr(412, 11));
    f.bf.interval(x.data(), (int)x.size(), k, n_occ);
    CHECK(n_occ == 3);
    f.bf.interval(z.data(), (int)z.size(), k, n_occ);
    CHECK(n_occ == 5);
    f.bf.interval(y.data(), (int)y.size(), k, n_occ);
    CHECK(n_occ == 1);
    f.bf.interval(w.data(), (int)w.size(), k, n_occ);   // an 11-bp window of Y, also in its diverged copy
    CHECK(n_occ == 2);
}

TEST_CASE("FMI_search load: count[] and the sentinel row and every SA row match a brute-force suffix array at SA strides 1-64"
          * doctest::test_suite("unit/mem_seed")) {
    Fixture &f = fixture();
    const Brute &bf = f.bf;
    for (int compx = 0; compx <= 6; ++compx) {
        CAPTURE(compx);
        FMI_search *fmi = index_at(compx);
        // The loader detects the stride from the file's trailing tag.
        CHECK(fmi->sa_compx == compx);
        CHECK(fmi->reference_seq_len == bf.n + 1);
        // count[b] = 1 (the `$` row) + the number of bases of the text below b.
        const int64_t *count = fmi->count_data();
        int64_t below = 0;
        for (int b = 0; b < 4; ++b) {
            CAPTURE(b);
            CHECK(count[b] == 1 + below);
            below += count_base(bf, b);
        }
        CHECK(count[4] == 1 + bf.n);
        // The row whose BWT symbol is `$` is the row of the suffix at position 0.
        CHECK(fmi->sentinel_index == bf.rank[0]);
        // Every row 1..n resolved through the SA-sample pipeline (one interval
        // covering them all, max_occ = n so the stride is 1).
        SMEM all;
        std::memset(&all, 0, sizeof(all));
        all.k = 1;
        all.s = bf.n;
        std::vector<int64_t> coord((size_t)bf.n);
        int64_t cnt = 0, id = 0;
        fmi->get_sa_entries_prefetch(&all, coord.data(), &cnt, 1, (int32_t)bf.n, 0, id, 0);
        REQUIRE(cnt == bf.n);
        CHECK(id == bf.n);
        int64_t bad = 0;
        for (int64_t r = 1; r <= bf.n; ++r)
            if (coord[(size_t)r - 1] != bf.sa[(size_t)r]) {
                if (bad++ < 5) {
                    CAPTURE(r);
                    CHECK(coord[(size_t)r - 1] == bf.sa[(size_t)r]);
                }
            }
        CHECK(bad == 0);
        // The pac and the contig table are the same at every stride.
        CHECK(fixture().ix[compx]->bns()->l_pac == (int64_t)bf.fwd.size());
        for (size_t i = 0; i < bf.fwd.size(); ++i)
            if ((uint8_t)_get_pac(fixture().ix[compx]->pac(), (int64_t)i) != bf.fwd[i]) {
                CAPTURE(i);
                FAIL("pac differs from the stride-8 index");
            }
    }
}

TEST_CASE("get_sa_entries_prefetch: coordinates follow the max_occ stride rule under every --ks-dedup mode and sentinel arm"
          * doctest::test_suite("unit/mem_seed")) {
    Fixture &f = fixture();
    const Brute &bf = f.bf;

    // Intervals of strings that occur once (including near text position 0,
    // whose LF walks can reach the `$` row), a few times (the planted repeats)
    // and many times (single bases, dinucleotides, the AC repeat). Some are
    // repeated, adjacent and apart, so --ks-dedup has duplicates to collapse.
    std::vector<std::vector<uint8_t>> pats;
    for (int p = 0; p < 16; ++p) pats.push_back(std::vector<uint8_t>(bf.t.begin() + p, bf.t.begin() + p + 12));
    for (int b = 0; b < 4; ++b) pats.push_back(std::vector<uint8_t>(1, (uint8_t)b));
    pats.push_back({0, 1});
    pats.push_back({2, 3, 3});
    pats.push_back(bwa_tests::nt4(f.ref[0].seq.substr(130, 25)));   // X: 3 occurrences in the text
    pats.push_back(bwa_tests::nt4(f.ref[0].seq.substr(805, 20)));   // Z: 5 occurrences
    pats.push_back(bwa_tests::nt4(f.ref[0].seq.substr(300, 12)));   // inside (AC)x15
    pats.push_back(bwa_tests::nt4(f.ref[1].seq.substr(500, 10)));   // inside the poly-A run
    for (int p = 0; p < 30; ++p) {
        const size_t pos = (size_t)(37 * p + 5) % (bf.t.size() - 30);
        pats.push_back(std::vector<uint8_t>(bf.t.begin() + (long)pos, bf.t.begin() + (long)pos + 9 + p % 13));
    }
    // Intervals that start at the same row with different sizes: a string and
    // its extension by the next base of its first row's suffix. (k, s), not k
    // alone, identifies an interval.
    const size_t n_base_pats = pats.size();
    int n_same_k = 0;
    for (size_t i = 0; i < n_base_pats; ++i) {
        std::vector<uint8_t> q = pats[i];
        int64_t k, s;
        bf.interval(q.data(), (int)q.size(), k, s);
        for (int ext = 0; ext < 6 && s > 1; ++ext) {
            const int c = bf.at(k, (int64_t)q.size());
            if (c < 0) break;
            q.push_back((uint8_t)c);
            int64_t k2, s2;
            bf.interval(q.data(), (int)q.size(), k2, s2);
            REQUIRE(k2 == k);
            if (s2 < s) {
                pats.push_back(q);
                ++n_same_k;
            }
            s = s2;
        }
    }
    REQUIRE(n_same_k > 0);
    std::vector<SMEM> smems;
    for (const std::vector<uint8_t> &q : pats) {
        SMEM s;
        std::memset(&s, 0, sizeof(s));
        bf.interval(q.data(), (int)q.size(), s.k, s.s);
        REQUIRE(s.s >= 1);
        smems.push_back(s);
    }
    const size_t n_distinct = smems.size();
    for (size_t i = 0; i < n_distinct; i += 3) smems.push_back(smems[i]);
    smems.push_back(smems[n_distinct - 1]);
    smems.push_back(smems[n_distinct - 1]);

    const int32_t max_occs[] = {1, 2, 3, 7, 40, 500};
    // --ks-dedup auto reports its latch decision at the default verbosity.
    const int saved_verbose = bwa_verbose;
    bwa_verbose = 1;
    const char *modes[] = {"off", "on", "auto"};
    for (int compx = 0; compx <= 6; compx += 3) {
        FMI_search *fmi = index_at(compx);
        const int64_t mask = (1 << compx) - 1;
        // The resolved coordinate of row r. With drop_sentinel_offset, a walk
        // that reaches the row of position 0 before any sampled row reports 0
        // (bwa-mem2's behaviour, kept for --compat=bwa-mem2); otherwise every
        // row resolves to its suffix-array value.
        int64_t sentinel_rows = 0;
        auto expect_row = [&](int64_t row, int drop) -> int64_t {
            int64_t r = row, p = bf.sa[(size_t)row];
            for (;;) {
                if ((r & mask) == 0) return bf.sa[(size_t)row];
                if (p == 0) return drop ? 0 : bf.sa[(size_t)row];
                --p;
                r = bf.rank[(size_t)p];
            }
        };
        for (int32_t max_occ : max_occs) {
            for (const char *mode : modes) {
                ks_dedup_configure(mode);
                for (int drop = 0; drop <= 1; ++drop) {
                    CAPTURE(compx);
                    CAPTURE(max_occ);
                    CAPTURE(mode);
                    CAPTURE(drop);
                    std::vector<int64_t> want;
                    for (const SMEM &s : smems) {
                        const int64_t step = s.s > max_occ ? s.s / max_occ : 1;
                        int32_t c = 0;
                        for (int64_t j = s.k; j < s.k + s.s && c < max_occ; j += step, ++c) {
                            want.push_back(expect_row(j, drop));
                            if (drop && want.back() != expect_row(j, 0)) ++sentinel_rows;
                        }
                    }
                    std::vector<int64_t> got(want.size() + 8, -7);
                    int64_t cnt = 0, id = 0;
                    fmi->get_sa_entries_prefetch(smems.data(), got.data(), &cnt, (int64_t)smems.size(), max_occ, 0,
                                                 id, drop);
                    CHECK(cnt == (int64_t)want.size());
                    CHECK(id == (int64_t)want.size());
                    got.resize(want.size() + 8);
                    for (size_t i = want.size(); i < got.size(); ++i) CHECK(got[i] == -7);   // no overrun
                    got.resize(want.size());
                    CHECK(got == want);
                }
            }
        }
        // At stride 64 some of the rows near position 0 do reach the `$` row
        // unsampled, so the drop_sentinel_offset arm was exercised.
        if (compx == 6) CHECK(sentinel_rows > 0);
    }
    ks_dedup_configure(NULL);   // back to the environment / default resolution
    bwa_verbose = saved_verbose;
}

TEST_CASE("getSMEMsAllPosOneThread: round-1 SMEMs are the maximal exact matches at every lockstep width"
          * doctest::test_suite("unit/mem_seed")) {
    Fixture &f = fixture();
    FMI_search *fmi = index_at(3);
    Batch b(f.reads);
    LockstepGuard guard;
    const int min_lens[] = {2, 19};
    for (int min_len : min_lens) {
        std::vector<Smem> want;
        for (size_t r = 0; r < f.reads.size(); ++r) {
            const std::vector<Smem> v = oracle_round1(f, r, min_len);
            want.insert(want.end(), v.begin(), v.end());
        }
        std::sort(want.begin(), want.end());
        REQUIRE(want.size() > f.reads.size());
        for (int32_t w : widths(SMEM_LOCKSTEP_N_MAX, min_len == 19)) {
            CAPTURE(min_len);
            CAPTURE(w);
            g_smem_lockstep_n = w;
            const std::vector<Smem> got = run_round1(fmi, b, min_len);
            if (got != want) {
                CHECK(first_diff(got, want) == "");
                break;
            }
        }
    }
}

TEST_CASE("getSMEMsOnePosOneThread_lockstep: round-2 SMEMs are the maximal intervals through x with min_intv occurrences at every width"
          * doctest::test_suite("unit/mem_seed")) {
    Fixture &f = fixture();
    FMI_search *fmi = index_at(3);
    Batch b(f.reads);
    LockstepGuard guard;
    std::mt19937 rng(11u);
    const int32_t thresholds[] = {1, 2, 3, 4, 5, 6, 7, 12, 40, 200, 5000};
    std::vector<Probe> probes;
    for (int32_t r = 0; r < b.size(); ++r) {
        const int len = (int)f.reads[(size_t)r].size();
        for (int i = 0; len > 0 && i < 2; ++i)
            probes.push_back({r, (int32_t)(rng() % (uint32_t)len), thresholds[rng() % 11]});
    }
    const int min_lens[] = {2, 12, 19};
    for (int min_len : min_lens) {
        std::vector<Smem> want;
        for (const Probe &p : probes) {
            const std::vector<Smem> v = oracle_round2(f, (size_t)p.rid, p.x, p.t, min_len);
            want.insert(want.end(), v.begin(), v.end());
        }
        std::sort(want.begin(), want.end());
        REQUIRE(want.size() > 20);
        for (int32_t w : widths(SMEM_LOCKSTEP_N_MAX, min_len == 19)) {
            CAPTURE(min_len);
            CAPTURE(w);
            g_smem_lockstep_n = w;
            const std::vector<Smem> got = run_round2(fmi, b, probes, min_len);
            if (got != want) {
                CHECK(first_diff(got, want) == "");
                break;
            }
        }
    }
}

TEST_CASE("bwtSeedStrategyAllPosOneThread: round-3 seeds are the forward walk's in the scalar and lockstep drivers at every width"
          * doctest::test_suite("unit/mem_seed")) {
    Fixture &f = fixture();
    FMI_search *fmi = index_at(3);
    Batch b(f.reads);
    LockstepGuard guard;
    struct Case { int32_t max_intv; int min_len; };
    const Case cases[] = {{20, 20}, {5, 13}, {100, 3}, {2, 31}, {1, 20}};
    for (const Case &c : cases) {
        std::vector<Smem> want;
        for (size_t r = 0; r < f.reads.size(); ++r) {
            const std::vector<Smem> v = oracle_round3(f.bf, (uint32_t)r, f.reads[r], c.max_intv, c.min_len);
            want.insert(want.end(), v.begin(), v.end());
        }
        std::sort(want.begin(), want.end());
        if (c.max_intv > 1) REQUIRE(want.size() > 20);
        else REQUIRE(want.empty());   // s < 1 means no occurrence, and an empty interval is never emitted
        CAPTURE(c.max_intv);
        CAPTURE(c.min_len);
        CHECK(first_diff(run_round3(fmi, b, c.max_intv, c.min_len, false), want) == "");
        for (int32_t w : widths(BWTSEED_LOCKSTEP_N_MAX, c.max_intv == 20)) {
            CAPTURE(w);
            g_bwtseed_lockstep_n = w;
            const std::vector<Smem> got = run_round3(fmi, b, c.max_intv, c.min_len, true);
            if (got != want) {
                CHECK(first_diff(got, want) == "");
                break;
            }
        }
    }
}

TEST_CASE("mem_collect_smem: the three rounds sorted equal the oracle under -k -r -y -s and --smem-dedup"
          * doctest::test_suite("unit/mem_seed")) {
    Fixture &f = fixture();
    FMI_search *fmi = index_at(3);
    Batch b(f.reads);
    LockstepGuard guard;
    std::unique_ptr<mem_opt_t, void (*)(void *)> opt(mem_opt_init(), free);
    // The defaults (-k 19 -r 1.5 -y 20 -s 10), then each knob moved. split:
    // the re-seed length k * r, a half rounded down (19 * 1.5 = 28.5 -> 28,
    // 25 * 1.5 = 37.5 -> 37, 12 * 1.0 = 12, 19 * 10 = 190, 5 * 1.2 = 6).
    // want_r2 / want_r3: whether round 2 / round 3 must contribute (1), must
    // not (0), or either (-1), so a case that varies a round shows the round ran.
    struct Case { int k; float r; int split; int y; int s; int dedup; int want_r2; int want_r3; };
    const Case cases[] = {
        {19, 1.5f, 28, 20, 10, 0, 1, 1},    {19, 1.5f, 28, 20, 10, 1, 1, 1},  {12, 1.0f, 12, 20, 10, 0, 1, 1},
        {25, 1.5f, 37, 20, 10, 0, -1, 1},   {19, 10.0f, 190, 20, 10, 0, 0, 1}, {19, 1.5f, 28, 0, 10, 0, 1, 0},
        {19, 1.5f, 28, 3, 10, 0, 1, 1},     {19, 1.5f, 28, 20, 1, 0, 1, 1},   {12, 1.0f, 12, 5, 1000, 1, 1, 1},
        {5, 1.2f, 6, 50, 1000, 0, 1, 1},
    };
    size_t n_dups_removed = 0, n_cross_read_equal = 0;
    for (const Case &c : cases) {
        opt->min_seed_len = c.k;
        opt->split_factor = c.r;
        opt->max_mem_intv = (uint64_t)c.y;
        opt->split_width = c.s;
        opt->smem_dedup = c.dedup;
        CAPTURE(c.k);
        CAPTURE(c.r);
        CAPTURE(c.y);
        CAPTURE(c.s);
        CAPTURE(c.dedup);
        size_t n_round[3];
        const std::vector<Smem> want = oracle_collect(f, opt.get(), c.split, n_round);
        CHECK(n_round[0] > 0);
        if (c.want_r2 >= 0) CHECK((n_round[1] > 0) == (c.want_r2 == 1));
        if (c.want_r3 >= 0) CHECK((n_round[2] > 0) == (c.want_r3 == 1));
        if (c.dedup) {
            n_dups_removed += n_round[0] + n_round[1] + n_round[2] - want.size();
            // Neighbours equal in all but the read, which --smem-dedup must keep.
            for (size_t i = 1; i < want.size(); ++i)
                n_cross_read_equal += want[i].rid != want[i - 1].rid && want[i].m == want[i - 1].m
                                      && want[i].n == want[i - 1].n && want[i].k == want[i - 1].k
                                      && want[i].s == want[i - 1].s;
        }
        g_smem_lockstep_n = guard.smem_n;
        g_bwtseed_lockstep_n = guard.bwt_n;
        g_bwtseed_lockstep = 1;
        CHECK(first_diff(run_collect(fmi, b, opt.get()), want) == "");
        g_bwtseed_lockstep = 0;
        CHECK(first_diff(run_collect(fmi, b, opt.get()), want) == "");
        // Narrow and odd widths for both lockstep drivers.
        g_bwtseed_lockstep = 1;
        const int32_t widths[] = {1, 3, 64};
        for (int32_t w : widths) {
            CAPTURE(w);
            g_smem_lockstep_n = w;
            g_bwtseed_lockstep_n = w;
            CHECK(first_diff(run_collect(fmi, b, opt.get()), want) == "");
        }
    }
    // --smem-dedup had duplicates to remove (round 2 or 3 re-finding an SMEM),
    // and equal SMEMs of different reads to keep.
    CHECK(n_dups_removed > 0);
    CHECK(n_cross_read_equal > 0);
}
