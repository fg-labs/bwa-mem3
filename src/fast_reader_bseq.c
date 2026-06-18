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

/* trim_readno / kseq2bseq1 mirror the static helpers in bwa.cpp exactly:
 * strip a trailing /<digit> for any digit; per-field strdup; qual==NULL when
 * empty; zero the struct first (the output loop free()s every field). */
static inline void fr_trim_readno(kstring_t *s)
{
    if (s->l > 2 && s->s[s->l - 2] == '/' && isdigit((unsigned char)s->s[s->l - 1]))
        s->l -= 2, s->s[s->l] = 0;
}

static inline void fr_kseq2bseq1(const kseq_t *ks, bseq1_t *s)
{
    memset(s, 0, sizeof(*s));
    s->name    = strdup(ks->name.s);
    s->comment = ks->comment.l ? strdup(ks->comment.s) : 0;
    s->seq     = strdup(ks->seq.s);
    s->qual    = ks->qual.l ? strdup(ks->qual.s) : 0;
    s->l_seq   = strlen(s->seq);
}

void *fast_kseq_init(fast_reader_t *fr) { return kseq_init(fr); }

void fast_kseq_destroy(void *ks) { if (ks) kseq_destroy((kseq_t *)ks); }

bseq1_t *bseq_read_fast(int64_t chunk_size, int *n_, void *ks1_, void *ks2_, int64_t *s)
{
    kseq_t *ks = (kseq_t *)ks1_, *ks2 = (kseq_t *)ks2_;
    int64_t size = 0, m, n;
    bseq1_t *seqs;
    m = n = 0; seqs = 0;
    while (kseq_read(ks) >= 0) {
        if (ks2 && kseq_read(ks2) < 0) {                 /* 2nd file has fewer */
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
