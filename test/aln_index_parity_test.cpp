// Row-for-row parity of bwa-mem3's index against a stock `bwa index` (0.7.19).
//
// usage: aln_index_parity_test <stock-prefix> <mem3-prefix>
//        aln_index_parity_test --write-fasta <path>
//
// The `aln` port answers bwa's FM-index queries (bwt_occ, bwt_sa, L2, primary)
// from bwa-mem3's own index through src/aln/aln_bwt. That is only correct if the
// two indexes agree on every row. This test loads the stock .bwt/.sa with
// test-local copies of upstream's reader and query functions, then compares, for
// every row k in [0, seq_len] and every base c: occ(k, c), SA[k], and the
// header scalars primary and L2[0..4]. At the first mismatch it prints
// `FAIL: <field> differs at k=<k> c=<c>: stock=<v> mem3=<v>` and exits 1.
//
// --write-fasta writes the shared synthetic reference (test/unit/aln_test_ref.h)
// so the regression script indexes the same contigs the unit tests use.
//
// The stock functions below are copied from bwa's Apache2 branch (see
// src/aln/README.md), bwt.c, under the upstream notice below. They are MODIFIED
// from the original: placed in namespace stock, die()/fread_exact() wrappers
// replace the err_* I/O helpers, std:: prefixes, stock_-prefixed layout macros,
// the bwt_size computation uses std::ftell,
// the bwt_t struct is local, and the bwt_cal_sa/bwt_gen_cnt_table paths are omitted.
//
// The MIT License
//
// Copyright (c) 2008 Genome Research Ltd (GRL).
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "FMI_search.h"
#include "aln/aln_bwt.h"
#include "unit/aln_test_ref.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

namespace stock {

typedef uint64_t bwtint_t;
typedef unsigned char ubyte_t;

#define OCC_INTV_SHIFT 7
#define OCC_INTERVAL (1LL << OCC_INTV_SHIFT)
#define OCC_INTV_MASK (OCC_INTERVAL - 1)

typedef struct {
    bwtint_t primary;
    bwtint_t L2[5];
    bwtint_t seq_len;
    bwtint_t bwt_size;
    uint32_t *bwt;
    uint32_t cnt_table[256];
    int sa_intv;
    bwtint_t n_sa;
    bwtint_t *sa;
} bwt_t;

// The upstream layout is only correct for OCC_INTERVAL == 0x80.
#define stock_bwt_bwt(b, k) ((b)->bwt[((k) >> 7 << 4) + sizeof(bwtint_t) + (((k) & 0x7f) >> 4)])
#define stock_bwt_occ_intv(b, k) ((b)->bwt + ((k) >> 7 << 4))
#define stock_bwt_B0(b, k) (stock_bwt_bwt(b, k) >> ((~(k) & 0xf) << 1) & 3)

static void die(const char *msg, const char *fn)
{
    std::fprintf(stderr, "FAIL: %s: %s\n", msg, fn);
    std::exit(2);
}

static FILE *open_or_die(const char *fn)
{
    FILE *fp = std::fopen(fn, "rb");
    if (!fp) die("cannot open", fn);
    return fp;
}

static void fread_exact(void *p, size_t sz, size_t n, FILE *fp, const char *fn)
{
    if (std::fread(p, sz, n, fp) != n) die("short read", fn);
}

static inline bwtint_t bwt_invPsi(const bwt_t *bwt, bwtint_t k);

static inline int __occ_aux(uint64_t y, int c)
{
    // reduce nucleotide counting to bits counting
    y = ((c & 2) ? y : ~y) >> 1 & ((c & 1) ? y : ~y) & 0x5555555555555555ull;
    // count the number of 1s in y
    y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
    return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, ubyte_t c)
{
    bwtint_t n;
    uint32_t *p, *end;

    if (k == bwt->seq_len) return bwt->L2[c + 1] - bwt->L2[c];
    if (k == (bwtint_t)(-1)) return 0;
    k -= (k >= bwt->primary); // because $ is not in bwt

    // retrieve Occ at k/OCC_INTERVAL
    n = ((bwtint_t *)(p = stock_bwt_occ_intv(bwt, k)))[c];
    p += sizeof(bwtint_t); // jump to the start of the first BWT cell

    // calculate Occ up to the last k/32
    end = p + (((k >> 5) - ((k & ~OCC_INTV_MASK) >> 5)) << 1);
    for (; p < end; p += 2) n += __occ_aux((uint64_t)p[0] << 32 | p[1], c);

    // calculate Occ
    n += __occ_aux(((uint64_t)p[0] << 32 | p[1]) & ~((1ull << ((~k & 31) << 1)) - 1), c);
    if (c == 0) n -= ~k & 31; // corrected for the masked bits

    return n;
}

static inline bwtint_t bwt_invPsi(const bwt_t *bwt, bwtint_t k) // compute inverse CSA
{
    bwtint_t x = k - (k > bwt->primary);
    x = stock_bwt_B0(bwt, x);
    x = bwt->L2[x] + bwt_occ(bwt, k, x);
    return k == bwt->primary ? 0 : x;
}

bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t k)
{
    bwtint_t sa = 0, mask = bwt->sa_intv - 1;
    while (k & mask) {
        ++sa;
        k = bwt_invPsi(bwt, k);
    }
    /* without setting bwt->sa[0] = -1, the following line should be
       changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
    return sa + bwt->sa[k / bwt->sa_intv];
}

void bwt_restore_sa(const char *fn, bwt_t *bwt)
{
    char skipped[256];
    bwtint_t primary;

    FILE *fp = open_or_die(fn);
    fread_exact(&primary, sizeof(bwtint_t), 1, fp, fn);
    if (primary != bwt->primary) die("SA-BWT inconsistency: primary is not the same", fn);
    fread_exact(skipped, sizeof(bwtint_t), 4, fp, fn); // skip
    fread_exact(&bwt->sa_intv, sizeof(bwtint_t), 1, fp, fn);
    fread_exact(&primary, sizeof(bwtint_t), 1, fp, fn);
    if (primary != bwt->seq_len) die("SA-BWT inconsistency: seq_len is not the same", fn);

    bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
    bwt->sa = (bwtint_t *)calloc(bwt->n_sa, sizeof(bwtint_t));
    bwt->sa[0] = -1;

    fread_exact(bwt->sa + 1, sizeof(bwtint_t), bwt->n_sa - 1, fp, fn);
    std::fclose(fp);
}

bwt_t *bwt_restore_bwt(const char *fn)
{
    bwt_t *bwt = (bwt_t *)calloc(1, sizeof(bwt_t));
    FILE *fp = open_or_die(fn);
    std::fseek(fp, 0, SEEK_END);
    bwt->bwt_size = (std::ftell(fp) - sizeof(bwtint_t) * 5) >> 2;
    bwt->bwt = (uint32_t *)calloc(bwt->bwt_size, 4);
    std::fseek(fp, 0, SEEK_SET);
    fread_exact(&bwt->primary, sizeof(bwtint_t), 1, fp, fn);
    fread_exact(bwt->L2 + 1, sizeof(bwtint_t), 4, fp, fn);
    fread_exact(bwt->bwt, 4, bwt->bwt_size, fp, fn);
    bwt->seq_len = bwt->L2[4];
    std::fclose(fp);
    return bwt;
}

} // namespace stock

static std::string fmt(uint64_t v)
{
    return v == (uint64_t)(-1) ? std::string("-1") : std::to_string((unsigned long long)v);
}

[[noreturn]] static void mismatch(const char *field, uint64_t k, int c, uint64_t stock_v, uint64_t mem3_v)
{
    std::printf("FAIL: %s differs at k=%s c=%d: stock=%s mem3=%s\n", field,
                fmt(k).c_str(), c, fmt(stock_v).c_str(), fmt(mem3_v).c_str());
    std::exit(1);
}

int main(int argc, char *argv[])
{
    if (argc == 3 && std::strcmp(argv[1], "--write-fasta") == 0)
        return aln_test::write_test_fasta(argv[2]) ? 0 : 1;
    if (argc != 3) {
        std::fprintf(stderr, "usage: %s <stock-prefix> <mem3-prefix>\n       %s --write-fasta <path>\n",
                     argv[0], argv[0]);
        return 2;
    }

    const std::string stock_prefix = argv[1];
    stock::bwt_t *s = stock::bwt_restore_bwt((stock_prefix + ".bwt").c_str());
    stock::bwt_restore_sa((stock_prefix + ".sa").c_str(), s);

    FMI_search fmi(argv[2]);
    fmi.load_index(false, 1);
    bwa_aln::bwt_t *m = bwa_aln::bwt_from_fmi(&fmi);

    if (s->primary != m->primary) mismatch("primary", 0, 0, s->primary, m->primary);
    for (int i = 0; i < 5; ++i)
        if (s->L2[i] != m->L2[i]) mismatch("L2", 0, i, s->L2[i], m->L2[i]);
    if (s->seq_len != m->seq_len) mismatch("seq_len", 0, 0, s->seq_len, m->seq_len);

    for (uint64_t k = 0; k <= s->seq_len; ++k) {
        for (int c = 0; c < 4; ++c) {
            uint64_t a = stock::bwt_occ(s, k, (stock::ubyte_t)c);
            uint64_t b = bwa_aln::bwt_occ(m, k, (bwa_aln::ubyte_t)c);
            if (a != b) mismatch("occ", k, c, a, b);
        }
        uint64_t a = stock::bwt_sa(s, k);
        uint64_t b = bwa_aln::bwt_sa(m, k);
        if (a != b) mismatch("sa", k, 0, a, b);
    }

    std::printf("PASS: aln_index_parity (%llu rows)\n", (unsigned long long)(s->seq_len + 1));
    bwa_aln::bwt_destroy(m);
    return 0;
}
