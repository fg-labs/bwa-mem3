#include "fast_reader_bseq.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kstring.h"
#include "kseq.h"
#include "utils.h"   /* err_fatal */
#include "stage_prof.h"

/* kseq instantiated over the fast_reader byte source. This is the only TU that
 * uses this handle type, so there is no kseq_t name collision with the gzFile
 * instantiation used elsewhere. */
KSEQ_INIT(fast_reader_t *, fast_reader_read)

/* fr_trim_readno mirrors the static helper in bwa.cpp exactly: strip a trailing
 * /<digit> for any digit. It runs BEFORE the copy and shortens ks->name.l, so
 * fr_kseq2bseq1 must read the live kseq length rather than re-scanning. */
static inline void fr_trim_readno(kstring_t *s)
{
    if (s->l > 2 && s->s[s->l - 2] == '/' && isdigit((unsigned char)s->s[s->l - 1]))
        s->l -= 2, s->s[s->l] = 0;
}

/* Length-aware field copy: malloc(len+1)+memcpy using the kseq-known length,
 * avoiding strdup's internal strlen scan on the hot path. The kseq source is
 * already NUL-terminated at .l; we set dst[len] explicitly so the copy is
 * NUL-terminated regardless. Aborts on OOM to match the house style of failing
 * loudly on allocation failure (see the realloc in bseq_read_fast below). */
static inline char *fr_dup_field(const char *src, size_t len)
{
    char *dst = (char *)malloc(len + 1);
    if (dst == NULL)
        err_fatal(__func__, "failed to allocate %zu bytes for a read field", len + 1);
    memcpy(dst, src, len);
    dst[len] = '\0';
    return dst;
}

/* fr_kseq2bseq1 produces output identical to bwa.cpp's kseq2bseq1, but copies
 * each field with its kseq-known length (ks->*.l) so the copy never scans for a
 * terminator, and takes l_seq from ks->seq.l directly. qual==NULL when empty.
 *
 * The leading memset is retained on purpose and is NOT the cost the strlen
 * removal targets: bseq_read_fast grows `seqs` with realloc, which leaves new
 * entries uninitialized, and the output loop in fastmap.cpp free()s sam/bams
 * unconditionally — so every field must start well-defined. Zeroing the whole
 * struct (rather than hand-listing sam/bams/n_bams/cap_bams) stays correct if a
 * field is ever added to bseq1_t, matching the documented rationale on the
 * bwa.cpp sibling. id is overwritten by the caller. */
static inline void fr_kseq2bseq1(const kseq_t *ks, bseq1_t *s)
{
    memset(s, 0, sizeof(*s));
    s->name    = fr_dup_field(ks->name.s, ks->name.l);
    s->comment = ks->comment.l ? fr_dup_field(ks->comment.s, ks->comment.l) : 0;
    s->seq     = fr_dup_field(ks->seq.s, ks->seq.l);
    s->qual    = ks->qual.l ? fr_dup_field(ks->qual.s, ks->qual.l) : 0;
    s->l_seq   = (int)ks->seq.l;
}

void *fast_kseq_init(fast_reader_t *fr) { return kseq_init(fr); }

void fast_kseq_destroy(void *ks) { if (ks) kseq_destroy((kseq_t *)ks); }

bseq1_t *bseq_read_fast(int64_t chunk_size, int *n_, void *ks1_, void *ks2_, int64_t *s)
{
    kseq_t *ks = (kseq_t *)ks1_, *ks2 = (kseq_t *)ks2_;
    int64_t size = 0, m, n;
    bseq1_t *seqs;
    m = n = 0; seqs = 0;
    for (;;) {
        /* Tokenize the next record(s). kseq_read does the FASTQ tokenization
         * (newline scan + copy of each field into its kstring_t) and pulls
         * bytes through the codec layer, which charges its own read()/inflate
         * time to the diskwait/decompress timers from inside this call. To make
         * read_parse account for the tokenization itself (the dominant, formerly
         * uninstrumented cost), time the whole kseq_read region and subtract the
         * IO delta the codec already counted. The IO interval is strictly nested
         * inside this wall bracket, so the remainder is non-negative. */
        double _tk0 = sp_enabled() ? sp_wall() : 0.0, _io0 = 0.0;
        if (sp_enabled()) { double d, c; sp_read_get(&d, &c, NULL); _io0 = d + c; }
        int r1 = kseq_read(ks);
        int r2 = (r1 >= 0 && ks2) ? kseq_read(ks2) : 0;
        if (sp_enabled()) {
            double d, c; sp_read_get(&d, &c, NULL);
            sp_read_add(2, (sp_wall() - _tk0) - ((d + c) - _io0));
        }
        if (r1 < 0) break;                               /* clean EOF on 1st file */
        if (ks2 && r2 < 0) {                             /* 2nd file has fewer */
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
        fr_trim_readno(&ks->name);
        fr_kseq2bseq1(ks, &seqs[n]);
        seqs[n].id = n;
        size += seqs[n++].l_seq;
        if (ks2) {
            fr_trim_readno(&ks2->name);
            fr_kseq2bseq1(ks2, &seqs[n]);
            seqs[n].id = n;
            size += seqs[n++].l_seq;
        }
        if (sp_enabled()) sp_read_add(2, sp_wall() - _tp);
        if (size >= chunk_size && (n & 1) == 0) break;   /* even-parity cut, all modes */
    }
    if (size == 0) {                                      /* 1st file has fewer */
        if (ks2 && kseq_read(ks2) >= 0)
            fprintf(stderr, "[W::%s] the 1st file has fewer sequences.\n", __func__);
    }
    *n_ = n;
    *s = size;
    return seqs;
}
