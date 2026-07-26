#include "fast_reader_bseq.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fr_fastq.h"
#include "utils.h"   /* err_fatal */
#include "stage_prof.h"
#include "read_arena.h"

/* This adapter turns fr_fastq record slices into the bseq1_t array the pipeline
 * consumes. fr_fastq replaced the kseq layer (see fr_fastq.c): records are now
 * sliced once out of the reader's own buffer rather than copied through a
 * kstring_t intermediate. The public handle names below are kept for the
 * pipeline's call sites; the opaque handle is now an fr_fastq_t*. */

/* Strip a trailing /<digit> from the name, matching bwa's trim_readno. It used
 * to mutate the kseq kstring_t in place before the copy; here it just shortens
 * the length we copy out of the (const) record slice. */
static inline size_t fr_trim_readno_len(const char *name, size_t l)
{
    if (l > 2 && name[l - 2] == '/' && isdigit((unsigned char)name[l - 1])) return l - 2;
    return l;
}

/* Length-aware field copy: malloc(len+1)+memcpy using the parser-known length,
 * avoiding strdup's internal strlen scan on the hot path. The source slice is
 * not assumed NUL-terminated; we set dst[len] explicitly. Aborts on OOM to
 * match the house style of failing loudly on allocation failure (see the
 * realloc in bseq_read_fast below). */
static inline char *fr_dup_field(const char *src, size_t len)
{
    char *dst = (char *)malloc(len + 1);
    if (dst == NULL)
        err_fatal(__func__, "failed to allocate %zu bytes for a read field", len + 1);
    memcpy(dst, src, len);
    dst[len] = '\0';
    return dst;
}

/* Copy one parsed record into a bseq1_t, producing output identical to the old
 * kseq path: per-field copy with comment/qual==NULL when empty, l_seq from the
 * sequence length, and the readno suffix trimmed from the name.
 *
 * The leading memset is retained on purpose: bseq_read_fast grows `seqs` with
 * realloc, which leaves new entries uninitialized, and the output loop in
 * fastmap.cpp free()s sam/bams unconditionally — so every field must start
 * well-defined. Zeroing the whole struct (rather than hand-listing
 * sam/bams/n_bams/cap_bams) stays correct if a field is ever added to bseq1_t.
 * id is overwritten by the caller. */
static inline void fr_rec_to_bseq1(const fr_fastq_rec_t *r, bseq1_t *s, read_arena_t *arena)
{
    memset(s, 0, sizeof(*s));
    /* Honor the bseq1_t.meth_base_ot -1 sentinel ("non-meth") from bwa.h: the
     * memset above would otherwise leave it 0, which the seed-chemistry filter
     * reads as OB. --meth ingest overwrites it with the read-number 0/1. */
    s->meth_base_ot = -1;
    size_t name_l = fr_trim_readno_len(r->name, r->name_l);
    /* PIPE-F6: name/seq/qual are carved from the per-chunk bump arena instead
     * of individually malloc'd. read_arena_dup does the same length-known
     * malloc+memcpy+NUL as fr_dup_field, so the bytes, lengths, and NUL
     * termination are identical — only the backing allocation changes. comment
     * stays heap-owned (fr_dup_field): --meth (fastmap step 0) frees and
     * reassigns it, and the non-copy_comment path frees it early, so its
     * ownership is not uniform enough to live in the arena. */
    s->name    = read_arena_dup(arena, r->name, name_l);
    s->comment = r->comment_l ? fr_dup_field(r->comment, r->comment_l) : 0;
    s->seq     = read_arena_dup(arena, r->seq, r->seq_l);
    s->qual    = r->qual_l ? read_arena_dup(arena, r->qual, r->qual_l) : 0;
    s->l_seq   = (int)r->seq_l;
}

void *fast_kseq_init(fast_reader_t *fr) { return fr_fastq_init(fr); }

void fast_kseq_destroy(void *p) { fr_fastq_destroy((fr_fastq_t *)p); }

bseq1_t *bseq_read_fast(int64_t chunk_size, int *n_, void *ks1_, void *ks2_, int64_t *s,
                        read_arena_t **arena_out)
{
    fr_fastq_t *p1 = (fr_fastq_t *)ks1_, *p2 = (fr_fastq_t *)ks2_;
    int64_t size = 0, m, n;
    bseq1_t *seqs;
    m = n = 0; seqs = 0;
    /* `arena_out` is in/out: NULL in means "create a fresh arena", non-NULL means
     * "keep carving from this one". Cohort slicing reads one batch in several
     * calls but consumes all of its fields together at output, so the arena's
     * lifetime is the COHORT, not the individual read -- carrying it across
     * slices is what keeps that one arena, destroyed once, correct. */
    const int arena_created_here = (*arena_out == NULL);
    read_arena_t *arena = arena_created_here ? read_arena_create() : *arena_out;
    for (;;) {
        /* Tokenize the next record(s). fr_fastq_next scans record boundaries and
         * pulls bytes through the codec layer, which charges its read()/inflate
         * time to the diskwait/decompress timers from inside this call. Time the
         * whole region and subtract that IO delta so read_parse accounts for the
         * tokenization (boundary scan + any multi-line compaction). The IO
         * interval is strictly nested inside the wall bracket, so the remainder
         * is non-negative. */
        double _tk0 = sp_enabled() ? sp_wall() : 0.0, _io0 = 0.0;
        if (sp_enabled()) { double d, c; sp_read_get(&d, &c, NULL); _io0 = d + c; }
        fr_fastq_rec_t r1, r2;
        int g1 = fr_fastq_next(p1, &r1);
        int g2 = (g1 == 1 && p2) ? fr_fastq_next(p2, &r2) : 0;
        if (sp_enabled()) {
            double d, c; sp_read_get(&d, &c, NULL);
            sp_read_add(2, (sp_wall() - _tk0) - ((d + c) - _io0));
        }
        if (g1 != 1) break;                              /* clean EOF or malformed 1st file */
        if (p2 && g2 != 1) {                             /* 2nd file has fewer */
            fprintf(stderr, "[W::%s] the 2nd file has fewer sequences.\n", __func__);
            break;
        }
        if (n >= m) {
            m = m ? m << 1 : 256;
            /* Grow via a temp so a failed realloc doesn't leak the old buffer
             * (cppcheck memleakOnRealloc). Abort loudly on OOM rather than
             * returning a short batch: the pipeline reads n_seqs==0 as clean
             * EOF and would silently truncate the alignment. Matches bwa-mem3's
             * house style of aborting on allocation failure. */
            bseq1_t *tmp = (bseq1_t *)realloc(seqs, m * sizeof(bseq1_t));
            if (tmp == NULL)
                err_fatal(__func__, "failed to grow read buffer to %ld records", (long)m);
            seqs = tmp;
        }
        double _tp = sp_enabled() ? sp_wall() : 0.0;
        fr_rec_to_bseq1(&r1, &seqs[n], arena);
        seqs[n].id = n;
        size += seqs[n++].l_seq;
        if (p2) {
            fr_rec_to_bseq1(&r2, &seqs[n], arena);
            seqs[n].id = n;
            size += seqs[n++].l_seq;
        }
        if (sp_enabled()) sp_read_add(2, sp_wall() - _tp);
        if (size >= chunk_size && (n & 1) == 0) break;   /* even-parity cut, all modes */
    }
    if (size == 0) {                                      /* 1st file has fewer */
        fr_fastq_rec_t rr;
        if (p2 && fr_fastq_next(p2, &rr) == 1)
            fprintf(stderr, "[W::%s] the 1st file has fewer sequences.\n", __func__);
    }
    /* PIPE-F6: hand the arena (backing every read's name/seq/qual) to the caller,
     * which destroys it in the write stage once everything it backs is consumed.
     * An empty batch (n == 0 → seqs == NULL, read as clean EOF by the pipeline)
     * carved nothing, so an arena created HERE is freed here and reported as
     * none. A CARRIED arena is never destroyed on this path: it still backs the
     * earlier slices of an in-flight cohort, and it is not ours to free. */
    if (n == 0 && arena_created_here) { read_arena_destroy(arena); arena = NULL; }
    *arena_out = arena;
    *n_ = n;
    *s = size;
    return seqs;
}
