/* alt_lift.cpp -- see alt_lift.h. The table and its interval arithmetic, plus
 * mem_alt_project, which applies them to a read's regions; it lives here rather
 * than beside dedup in bwamem.cpp, sharing dedup's redundancy threshold through
 * mem_overlap_redundant (bwamem.h). */

/* htslib headers before any bwa-mem3 header that pulls in kstring.h (they share
 * the KSTRING_H include guard). Only the CIGAR op constants and tables are used. */
#include "htslib/sam.h"

#include "alt_lift.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <vector>

#include "bwamem.h"
#include "khash.h"
#include "utils.h"

KHASH_MAP_INIT_STR(alt_name, int)

/* Op classes of a .alt CIGAR. H counts as query: the .alt query is the whole ALT
 * contig, hard-clipped bases included. */
static inline int op_is_aligned(int op) { return bam_cigar_type(op) == 3; }             /* M, =, X */
static inline int op_is_clip(int op) { return op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP; }
static inline int op_eats_query(int op) { return (bam_cigar_type(op) & 1) || op == BAM_CHARD_CLIP; }
static inline int op_eats_ref(int op) { return (bam_cigar_type(op) & 2) != 0; }

/* The ops of a CIGAR the .alt sidecar carries, BAM-encoded (len << 4 | BAM_C*).
 * Not sam_parse_cigar: it accepts zero-length ops and logs every malformed line
 * through hts_log, while a malformed .alt line here must be skipped and counted
 * under the loader's own warning cap. P is dropped (it consumes nothing). */
static int parse_cigar(const char *s, uint32_t **out, int32_t *n_out, int32_t *qlen_out)
{
    int n = 0, cap = 8;
    int64_t qlen = 0;
    uint32_t *c = (uint32_t *) malloc(cap * sizeof(uint32_t));
    xassert(c != NULL, "out of memory: .alt cigar");
    while (*s) {
        /* digits only before each op: strtol alone would take a sign or
         * leading whitespace ("10M+5I", "10M 5I") */
        if (*s < '0' || *s > '9') { free(c); return -1; }
        char *end;
        long len = strtol(s, &end, 10);
        if (end == s || len <= 0 || len >= (1L << 28) || !*end) { free(c); return -1; }
        const int op = bam_cigar_table[(unsigned char) *end];
        if (op < 0 || op == BAM_CBACK) { free(c); return -1; }
        s = end + 1;
        if (op == BAM_CPAD) continue;
        if (n == cap) {
            cap <<= 1;
            uint32_t *t = (uint32_t *) realloc(c, cap * sizeof(uint32_t));
            xassert(t != NULL, "out of memory: .alt cigar");
            c = t;
        }
        c[n++] = bam_cigar_gen(len, op);
        if (op_eats_query(op)) qlen += len;
    }
    if (n == 0 || qlen > INT32_MAX) { free(c); return -1; }
    *out = c; *n_out = n; *qlen_out = (int32_t) qlen;
    return 0;
}

/* Parse a whole SAM integer field (FLAG or POS): decimal digits only, no sign,
 * nothing after them. Returns 0 and sets *out, or -1. */
static int parse_uint_field(const char *s, long *out)
{
    if (*s < '0' || *s > '9') return -1;
    char *end;
    errno = 0;
    const long v = strtol(s, &end, 10);
    if (errno != 0 || *end != 0) return -1;
    *out = v;
    return 0;
}

/* Lines named one by one in a warning before the summary count takes over, so a
 * .alt with a systematic problem does not flood stderr. */
#define ALT_LIFT_MAX_LINE_WARNINGS 10

static void build_offsets(alt_lift1_t *a)
{
    const int n = a->n_cigar;
    a->qoff = (int64_t *) malloc((size_t) (n + 1) * 2 * sizeof(int64_t));
    xassert(a->qoff != NULL, "out of memory: .alt offsets");
    a->roff = a->qoff + (n + 1);
    a->qoff[0] = 0; a->roff[0] = a->tpos;
    for (int i = 0; i < n; ++i) {
        const int op = (int) bam_cigar_op(a->cigar[i]);
        const int64_t len = (int64_t) bam_cigar_oplen(a->cigar[i]);
        a->qoff[i + 1] = a->qoff[i] + (op_eats_query(op) ? len : 0);
        a->roff[i + 1] = a->roff[i] + (op_eats_ref(op) ? len : 0);
    }
}

/* bns_restore reads the same file for the ALT contig names alone (to set is_alt)
 * and keeps nothing else; this pass keeps the alignments, so it reads it again. */
alt_lift_t *alt_lift_load(const char *path, const bntseq_t *bns)
{
    FILE *fp = fopen(path, "r");
    if (fp == NULL) return NULL;
    khash_t(alt_name) *h = kh_init(alt_name);
    xassert(h != NULL, "out of memory: .alt name hash");
    for (int i = 0; i < bns->n_seqs; ++i) {
        int absent;
        khint_t k = kh_put(alt_name, h, bns->anns[i].name, &absent);
        kh_val(h, k) = i;
    }

    alt_lift_t *l = (alt_lift_t *) calloc(1, sizeof(alt_lift_t));
    xassert(l != NULL, "out of memory: alt_lift_t");
    l->n_seqs = bns->n_seqs;
    l->a = (alt_lift1_t *) calloc((size_t) bns->n_seqs, sizeof(alt_lift1_t));
    xassert(l->a != NULL, "out of memory: alt_lift_t");
    for (int i = 0; i < bns->n_seqs; ++i) l->a[i].trid = -1;

    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    long lineno = 0;
    alt_lift_skips_t *sk = &l->skip;
    /* A malformed line: named with its reason while under the warning cap, then
     * counted. */
#define ALT_LIFT_BAD(...) do { \
        if (sk->bad < ALT_LIFT_MAX_LINE_WARNINGS) { \
            fprintf(stderr, "[W::%s] %s line %ld: skipping ", __func__, path, lineno); \
            fprintf(stderr, __VA_ARGS__); \
            fputc('\n', stderr); \
        } \
        ++sk->bad; \
    } while (0)
    while ((len = getline(&line, &cap, fp)) > 0) {
        ++lineno;
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) line[--len] = 0;
        if (len == 0 || line[0] == '@') continue;
        char *f[6];
        int nf = 0;
        char *p = line;
        while (nf < 6) {
            f[nf++] = p;
            char *t = strchr(p, '\t');
            if (t == NULL) break;
            *t = 0; p = t + 1;
        }
        if (nf < 6) { ALT_LIFT_BAD("a line with fewer than 6 fields"); continue; }
        khint_t q = kh_get(alt_name, h, f[0]);
        if (q == kh_end(h)) { ++sk->no_ctg; continue; }
        long flag, pos;
        if (parse_uint_field(f[1], &flag) != 0 || parse_uint_field(f[3], &pos) != 0) {
            ALT_LIFT_BAD("%s (FLAG or POS is not a non-negative integer)", f[0]);
            continue;
        }
        if (flag & 0x4 || pos < 1 || strcmp(f[2], "*") == 0 || strcmp(f[5], "*") == 0) { ++sk->unmapped; continue; }
        if (flag & 0x900) { ++sk->not_primary; continue; }    /* secondary / supplementary */
        khint_t t = kh_get(alt_name, h, f[2]);
        if (t == kh_end(h)) { ++sk->no_target; continue; }
        const int qrid = kh_val(h, q), trid = kh_val(h, t);
        const int64_t tlen = bns->anns[trid].len;
        if (bns->anns[trid].is_alt) { ++sk->alt_target; continue; }
        alt_lift1_t *a = &l->a[qrid];
        if (a->trid >= 0) { ++sk->dup; continue; }              /* first usable line wins */
        /* POS can be as large as LONG_MAX; reject it before any arithmetic on it
         * (the reasons above do none, so they are counted first) */
        if (pos > tlen) {
            ALT_LIFT_BAD("%s (POS %ld is past the end of %s, length %lld)", f[0], pos, f[2], (long long) tlen);
            continue;
        }
        if (parse_cigar(f[5], &a->cigar, &a->n_cigar, &a->qlen) != 0) {
            ALT_LIFT_BAD("%s (bad CIGAR)", f[0]);
            continue;
        }
        int64_t rspan = 0;
        for (int i = 0; i < a->n_cigar; ++i)
            if (op_eats_ref((int) bam_cigar_op(a->cigar[i]))) rspan += (int64_t) bam_cigar_oplen(a->cigar[i]);
        /* 1 <= pos <= tlen, so tlen - (pos - 1) is at least 1 and cannot overflow */
        int ok = 1;
        if (a->qlen != bns->anns[qrid].len) {
            ALT_LIFT_BAD("%s (its CIGAR has %d query bases, the contig has %d)", f[0], (int) a->qlen,
                         (int) bns->anns[qrid].len);
            ok = 0;
        } else if (rspan > tlen - (pos - 1)) {
            ALT_LIFT_BAD("%s (its alignment ends at %lld, past the end of %s, length %lld)", f[0],
                         (long long) (pos - 1 + rspan), f[2], (long long) tlen);
            ok = 0;
        }
        if (!ok) { free(a->cigar); a->cigar = NULL; a->n_cigar = 0; a->qlen = 0; continue; }
        a->trid = trid;
        a->rev = (flag & 0x10) ? 1 : 0;
        a->tpos = pos - 1;
        a->tlen = tlen;
        build_offsets(a);
        ++l->n_lift;
    }
#undef ALT_LIFT_BAD
    free(line);
    fclose(fp);
    kh_destroy(alt_name, h);
    if (sk->bad > ALT_LIFT_MAX_LINE_WARNINGS)
        fprintf(stderr, "[W::%s] %s: %ld more malformed line(s) not named above\n",
                __func__, path, sk->bad - ALT_LIFT_MAX_LINE_WARNINGS);
    const long n_skip = sk->no_ctg + sk->no_target + sk->alt_target + sk->unmapped + sk->not_primary
        + sk->dup + sk->bad;
    if (n_skip)
        fprintf(stderr, "[W::%s] %s: skipped %ld line(s): %ld ALT contig not in the index, "
                "%ld target not in the index, %ld target is an ALT contig, %ld unmapped, "
                "%ld secondary/supplementary, %ld duplicate, %ld malformed\n",
                __func__, path, n_skip, sk->no_ctg, sk->no_target, sk->alt_target, sk->unmapped,
                sk->not_primary, sk->dup, sk->bad);
    return l;
}

void alt_lift_destroy(alt_lift_t *l)
{
    if (l == NULL) return;
    for (int i = 0; i < l->n_seqs; ++i) { free(l->a[i].cigar); free(l->a[i].qoff); }
    free(l->a);
    free(l);
}

/* The op containing query coordinate x (strict: the first op i with
 * qoff[i + 1] > x, for an interval's start) or the op that ends at or after it
 * (not strict: the first i with qoff[i + 1] >= x, for its end). For
 * 0 <= x < qlen (strict) or 0 < x <= qlen it consumes query, since qoff[i] <= x
 * (strict) or qoff[i] < x stays below qoff[i + 1]. */
static int find_op(const alt_lift1_t *a, int64_t x, int strict)
{
    int lo = 0, hi = a->n_cigar;
    while (lo < hi) {
        const int mid = (lo + hi) >> 1;
        if (strict ? a->qoff[mid + 1] > x : a->qoff[mid + 1] >= x) hi = mid;
        else lo = mid + 1;
    }
    return lo;
}

/* Reference position of query coordinate x within op i. */
static inline int64_t ref_at(const alt_lift1_t *a, int i, int64_t x)
{
    return a->roff[i] + (op_is_aligned((int) bam_cigar_op(a->cigar[i])) ? x - a->qoff[i] : 0);
}

int alt_lift_interval(const alt_lift_t *l, int rid, int64_t s, int64_t e, alt_lift_iv_t *out)
{
    if (l == NULL || rid < 0 || rid >= l->n_seqs) return -1;
    const alt_lift1_t *a = &l->a[rid];
    if (a->trid < 0 || a->tlen < 1 || e <= s || s < 0 || e > a->qlen) return -1;
    /* query coordinates in the CIGAR's orientation (the primary's forward strand) */
    const int64_t qs = a->rev ? a->qlen - e : s, qe = a->rev ? a->qlen - s : e;
    const int i0 = find_op(a, qs, 1), i1 = find_op(a, qe, 0);
    xassert(i0 <= i1 && i1 < a->n_cigar, "alt_lift_interval: interval outside the .alt CIGAR");
    int64_t b = ref_at(a, i0, qs), en = ref_at(a, i1, qe);
    /* alt_lift_load keeps only alignments that end within the target, so
     * tpos <= b <= en <= tlen; the start reaches tlen (and the interval is
     * empty on the target) only when it lies wholly in a trailing clip or
     * insertion of an alignment that ends at the target's end. */
    xassert(b >= 0 && en <= a->tlen, "alt_lift_interval: .alt alignment past its target");
    if (b >= a->tlen) b = a->tlen - 1;
    if (en <= b) en = b + 1;
    /* One walk over the ops from i0 to i1: the interval's query bases in aligned
     * (M/=/X) and clipped (S/H) ops, and the diagonal (inserted minus deleted
     * bases since the interval's start) at each op end inside it, the band the
     * read's alignment to the primary needs to follow the ALT. Query-only ops (I,
     * and S/H, which the callers reject anyway) move the diagonal up,
     * reference-only ops (D, N) down; D/N ops lie strictly between i0 and i1,
     * since both of those consume query. */
    int64_t d = 0, max_diag = 0, n_aln = 0, n_clip = 0;
    for (int k = i0; k <= i1; ++k) {
        const int op = (int) bam_cigar_op(a->cigar[k]);
        if (op_eats_query(op)) {
            const int64_t lo = qs > a->qoff[k] ? qs : a->qoff[k];
            const int64_t hi = qe < a->qoff[k + 1] ? qe : a->qoff[k + 1];
            if (op_is_aligned(op)) { n_aln += hi - lo; continue; }
            if (op_is_clip(op)) n_clip += hi - lo;
            d += hi - lo;
        } else if (op_eats_ref(op)) {
            d -= (int64_t) bam_cigar_oplen(a->cigar[k]);
        }
        const int64_t ad = d < 0 ? -d : d;
        if (ad > max_diag) max_diag = ad;
    }
    out->trid = a->trid;
    out->rev = a->rev;
    out->ps = b;
    out->pe = en;
    out->n_aln = n_aln;
    out->n_clip = n_clip;
    out->max_diag = max_diag;
    out->edge_ins = bam_cigar_op(a->cigar[i0]) == BAM_CINS || bam_cigar_op(a->cigar[i1]) == BAM_CINS;
    return 0;
}

/* --alt-project: are lifted region x and region y (either may be the lifted one)
 * the same locus? The reference test is dedup's redundancy threshold
 * (mem_overlap_redundant): they overlap on the reference by more than
 * mask_level_redun (0.95) of the shorter reference span, so a placement one
 * repeat unit away in a tandem repeat stays a separate hit, as ordinary dedup
 * keeps it, while every pair dedup would call redundant is one locus. The read
 * intervals must also share more than half of the shorter one. The rid test is
 * only a fast reject: an overlap in bwa's doubled coordinates already means the
 * same contig and strand. */
static inline int alt_project_same_locus(const mem_opt_t *opt, const mem_alnreg_t *x, const mem_alnreg_t *y)
{
    if (x->rid != y->rid) return 0;
    const int64_t or_ = (x->re < y->re ? x->re : y->re) - (x->rb > y->rb ? x->rb : y->rb);
    const int64_t mr = (x->re - x->rb < y->re - y->rb) ? x->re - x->rb : y->re - y->rb;
    if (or_ <= 0 || !mem_overlap_redundant(opt, or_, mr)) return 0;
    const int oq = (x->qe < y->qe ? x->qe : y->qe) - (x->qb > y->qb ? x->qb : y->qb);
    const int mq = (x->qe - x->qb < y->qe - y->qb) ? x->qe - x->qb : y->qe - y->qb;
    return oq > 0 && 2 * oq > mq;
}

/* --alt-project: may lifted region l replace unlifted region p at its locus?
 * Only if it outscores p by more than p's read bases outside l's read interval
 * are worth at the match score (-A): l's record covers l's read interval only,
 * so those bases would lose the alignment p gave them. An exact tie of the
 * adjusted scores keeps p, so a tie goes to the primary assembly. */
static inline int alt_project_replaces(const mem_opt_t *opt, const mem_alnreg_t *l, const mem_alnreg_t *p)
{
    const int64_t n_out = (l->qb > p->qb ? l->qb - p->qb : 0) + (p->qe > l->qe ? p->qe - l->qe : 0);
    return (int64_t) l->score > (int64_t) p->score + (int64_t) opt->a * n_out;
}

namespace {
/* One region in mem_alt_project's grouping: its index in a[], the contig it had
 * before the lift (the lift overwrites rid, and the old one breaks ties between
 * lifted regions), and where the grouping put it. */
struct alt_project_ent_t {
    int i, rid0;
    int kind;   /* AP_UNLIFTED, or lifted and AP_CAN_REPLACE or AP_CANNOT_REPLACE */
    int keep;
};
enum { AP_UNLIFTED = 0, AP_CAN_REPLACE = 1, AP_CANNOT_REPLACE = 2 };
} // namespace

/* --alt-project (see alt_lift.h). A region on ALT contig `rid` covers forward
 * contig interval [fb, fe) (a reverse-strand region, rb >= l_pac, is the mirror
 * of [2*l_pac - re, 2*l_pac - rb)); alt_lift_interval maps it to [ps, pe) on the
 * primary contig, and the strand flips when the ALT contig is reverse
 * complemented. qb/qe, score, truesc and w are kept: the ALT score ranks the read
 * (and is the record's AS), and mem_reg2aln aligns the read to the primary
 * reference at the lifted span when the record is written, passing
 * bwa_gen_cigar4 a band of at least the ALT's diagonal excursion in the interval
 * (carried as alt_bw). MAPQ's identity term (mem_approx_mapq_se) uses the longer
 * of the query span and the lifted reference span against the ALT score, so an
 * ALT deletion inside the read lengthens it and lowers MAPQ, more the longer the
 * deletion, while an ALT insertion shortens the reference span only and leaves
 * the term unchanged.
 *
 * A region is lifted only when the lift describes it: none of it lies in a clip
 * of the .alt CIGAR (ALT sequence with no primary counterpart), it neither
 * starts nor ends inside an ALT insertion (the read's end would have no primary
 * base to anchor on, so its alignment there would begin or end with an
 * insertion), at least half of it lies in aligned (M/=/X) ops, and the
 * alignment to the primary can follow the ALT's indels inside it: their widest
 * diagonal excursion is at most the band width opt->w, at most MEM_ALT_BW_MAX
 * (the most alt_bw can carry; this caps the excursion, not the band), and at
 * most the band bwa_gen_cigar4 will actually use for this region
 * (bwa_gen_cigar_band: it caps the band by the query length, to about 35 for a
 * 150-base region at the default -A/-O/-E, then floors it at the length
 * difference plus 3). Any other region stays on the ALT contig, with ALT-aware
 * handling, as without the flag.
 *
 * Then regions a lift made redundant are dropped, only ever in pairs where one
 * of the two is lifted and the two are the same locus (alt_project_same_locus);
 * pairs of unlifted regions are left as the default leaves them. Regions are
 * visited in one priority order: dedup's by-score order (mem_alnreg_score_lt:
 * score descending, rb, qb), then, for regions tied on all three, the lower
 * contig index before the lift; regions tied on all four keep their input
 * (dedup) order. (The rb before the lift would never break a tie of all four:
 * a lifted and an unlifted region with one contig index are on different
 * contigs, and two regions lifted from one ALT contig and strand with the same
 * rb after the lift had the same rb before it, as a lift maps distinct aligned
 * ALT bases to distinct primary bases and a region starting or ending in an ALT
 * insertion or clip is not lifted.) Which steps a lifted region takes part in
 * is decided first, by whether it can replace (step 0), and the order only
 * ranks regions within a step. So of two lifted regions at one locus that tie
 * on score, rb and qb and that both can, or both cannot, replace the unlifted
 * regions there, the one from the ALT contig that comes first in the index
 * wins; where only one of them can, that one is kept, whatever its index.
 *  0. A lifted region can replace the unlifted regions at its locus only if it
 *     replaces each of them (alt_project_replaces); tested against every
 *     unlifted region, so it does not depend on order.
 *  1. Which unlifted regions survive: visited in priority order, an unlifted
 *     region is dropped when a kept lifted region that can replace is at its
 *     locus (it comes first, as it outscores the unlifted one), and a lifted
 *     region that can replace is dropped when a kept one is at its locus.
 *  2. A lifted region that cannot replace is tested only against what step 1
 *     kept: visited in priority order, it is dropped when any kept region (an
 *     unlifted survivor, a lifted region kept in step 1, or one kept before it
 *     in this step) is at its locus, and kept otherwise. It never drops a
 *     region, so a primary-assembly hit is never replaced by a lifted hit that
 *     cannot replace every primary-assembly hit at its locus.
 * The score credit of alt_project_replaces applies only when a lifted region
 * would replace an unlifted one. Between two lifted regions at one locus the
 * higher score wins (steps 1 and 2), as in ordinary dedup, with no credit for
 * the read bases one aligns and the other does not. So whether an unlifted
 * region is at the locus can change which lifted region is kept: if the
 * higher-scoring lifted region cannot replace it and a lower-scoring one that
 * aligns more of the read can, the lower-scoring one is kept; with no unlifted
 * region there, the higher-scoring one is.
 * A region is dropped only for a kept region at its locus, so a region dropped
 * in favour of X never also takes down the regions X alone overlapped, whatever
 * the order of their pairs; an unlifted region is dropped only for a kept
 * lifted region that replaces it; and no kept lifted region shares a locus
 * with another kept region.
 * Neither existing suppression pass fits: dedup's window pass ran before the
 * lift (when the two were on different contigs) and has no notion of the
 * primary-wins-a-tie rule, and mem_mark_primary_se drops nothing (it marks
 * secondaries by query overlap alone, across loci, which is what keeps a
 * different locus counting for MAPQ).
 *
 * The survivors are left in that priority order, a refinement of dedup's
 * by-score order (lifting moves rb). Patching (mem_patch_reg) is not re-run: a
 * lifted region and a colinear primary one with no overlapping read bases stay
 * two regions, as the ALT piece and the primary piece of such a read do without
 * the flag. */
int mem_alt_project(const mem_opt_t *opt, const bntseq_t *bns, int n, mem_alnreg_t *a)
{
    const alt_lift_t *L = opt->alt_lift;
    if (L == NULL || n <= 0) return n;
    const int64_t l_pac = bns->l_pac;
    /* filled at the first lift; one buffer per thread, reused across reads */
    static thread_local std::vector<alt_project_ent_t> e;
    e.clear();
    for (int i = 0; i < n; ++i) {
        mem_alnreg_t *p = &a[i];
        if (!p->is_alt || p->rid < 0) continue;
        const int is_rev = p->rb >= l_pac;
        const int64_t fb = is_rev ? (l_pac << 1) - p->re : p->rb;
        const int64_t fe = is_rev ? (l_pac << 1) - p->rb : p->re;
        const int64_t off = bns->anns[p->rid].offset;
        alt_lift_iv_t iv;
        if (alt_lift_interval(L, p->rid, fb - off, fe - off, &iv) != 0) continue;
        if (iv.n_clip > 0 || iv.edge_ins || 2 * iv.n_aln < fe - fb
            || iv.max_diag > opt->w || iv.max_diag > MEM_ALT_BW_MAX) continue;
        /* mem_reg2aln passes a band of at least max_diag, and a wider band never
         * narrows bwa_gen_cigar4's, so this is the narrowest it can use */
        if (bwa_gen_cigar_band(opt->mat[0], opt->o_del, opt->e_del, opt->o_ins, opt->e_ins,
                               (int) iv.max_diag, p->qe - p->qb, iv.pe - iv.ps) < iv.max_diag) continue;
        if (e.empty()) {   /* nothing is lifted yet: every rid and rb is its own */
            e.resize(n);
            for (int k = 0; k < n; ++k) e[k] = {k, a[k].rid, AP_UNLIFTED, 0};
        }
        const int64_t toff = bns->anns[iv.trid].offset;
        if (is_rev ^ iv.rev) { p->rb = (l_pac << 1) - (toff + iv.pe); p->re = (l_pac << 1) - (toff + iv.ps); }
        else                 { p->rb = toff + iv.ps;                 p->re = toff + iv.pe; }
        p->rid = iv.trid;
        p->is_alt = 0;
        p->alt_bw = (unsigned) iv.max_diag + 1;
        p->dedup_re_rank = 0;  /* its `re` moved: unranked for the next dedup-only call */
        e[i].kind = AP_CAN_REPLACE;
    }
    if (e.empty()) return n;
    /* Step 0 */
    for (int i = 0; i < n; ++i) {
        if (e[i].kind == AP_UNLIFTED) continue;
        for (int j = 0; j < n; ++j)
            if (e[j].kind == AP_UNLIFTED && alt_project_same_locus(opt, &a[i], &a[j])
                && !alt_project_replaces(opt, &a[i], &a[j])) {
                e[i].kind = AP_CANNOT_REPLACE;
                break;
            }
    }
    std::stable_sort(e.begin(), e.end(), [a](const alt_project_ent_t &x, const alt_project_ent_t &y) {
        const mem_alnreg_t *p = &a[x.i], *q = &a[y.i];
        if (mem_alnreg_score_lt(p, q)) return true;
        if (mem_alnreg_score_lt(q, p)) return false;
        return x.rid0 < y.rid0;
    });
    /* Step 1 */
    for (int r = 0; r < n; ++r) {
        if (e[r].kind == AP_CANNOT_REPLACE) continue;
        e[r].keep = 1;
        for (int k = 0; k < r; ++k)
            if (e[k].keep && e[k].kind == AP_CAN_REPLACE && alt_project_same_locus(opt, &a[e[k].i], &a[e[r].i])) {
                e[r].keep = 0;
                break;
            }
    }
    /* Step 2: e[k].keep is final for every region but the lifted ones that
     * cannot replace, which are still 0 when visited after r */
    for (int r = 0; r < n; ++r) {
        if (e[r].kind != AP_CANNOT_REPLACE) continue;
        e[r].keep = 1;
        for (int k = 0; k < n; ++k)
            if (k != r && e[k].keep && alt_project_same_locus(opt, &a[e[k].i], &a[e[r].i])) {
                e[r].keep = 0;
                break;
            }
    }
    /* survivors in priority order; one buffer per thread, reused across reads */
    static thread_local std::vector<mem_alnreg_t> out;
    out.clear();
    for (int r = 0; r < n; ++r)
        if (e[r].keep) out.push_back(a[e[r].i]);
    std::copy(out.begin(), out.end(), a);
    return (int) out.size();
}
