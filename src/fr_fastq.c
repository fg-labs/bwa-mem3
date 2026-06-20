/* fr_fastq: single-copy FASTQ/FASTA parser over a fast_reader byte source.
 *
 * Strategy: keep a large refillable buffer and guarantee the *whole* current
 * record is resident before handing back field slices. For the dominant case
 * (single-line FASTQ/FASTA) every field is then a contiguous run in that buffer,
 * so the caller copies each field exactly once (no kstring_t intermediate).
 * Multi-line seq/qual are the only case that needs concatenation; those are
 * compacted into a small scratch buffer, matching what kseq did anyway.
 *
 * The grammar mirrors kseq exactly for realistic inputs (see the differential
 * test test/fr_fastq_diff_test.c). Pathological shapes that never occur in real
 * FASTQ (e.g. a bare "\r" line inside a sequence) are out of scope; the
 * --legacy-reader path retains exact kseq behavior for anyone who needs it.
 */
#include "fr_fastq.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#define FR_FASTQ_BUF_INIT (1u << 20)   /* 1 MiB; a record almost always fits */
#define FR_FASTQ_READ_MAX (1u << 20)   /* cap a single fast_reader_read */

/* scan_record return sentinel: the record is not yet fully buffered. */
#define FR_UNDERRUN (-3)

struct fr_fastq {
    fast_reader_t *fr;             /* borrowed byte source */
    unsigned char *buf;            /* refillable read buffer */
    size_t cap, beg, end;          /* capacity, parse cursor, valid byte count */
    int    eof;                    /* fast_reader_read returned 0 */
    int    err;                    /* fast_reader_read returned -1 (decode error) */
    /* scratch for assembling multi-line seq/qual (allocation persists across
     * records; the used length is recomputed per record). */
    unsigned char *seq_scr;  size_t seq_scr_cap;
    unsigned char *qual_scr; size_t qual_scr_cap;
};

static void *xrealloc(void *p, size_t n)
{
    void *q = realloc(p, n);
    if (q == NULL) { abort(); }   /* OOM on the read path is unrecoverable */
    return q;
}

fr_fastq_t *fr_fastq_init(fast_reader_t *fr)
{
    fr_fastq_t *p = (fr_fastq_t *)calloc(1, sizeof(*p));
    if (p == NULL) abort();
    p->fr  = fr;
    p->cap = FR_FASTQ_BUF_INIT;
    /* Test-only: shrink the initial buffer so the differential test can drive
     * mid-record underruns and buffer growth at every byte boundary. Never set
     * in production; the buffer still grows to fit any record either way. */
    const char *tc = getenv("FR_FASTQ_TEST_BUFSZ");
    if (tc && *tc) { long v = atol(tc); if (v > 0) p->cap = (size_t)v; }
    p->buf = (unsigned char *)xrealloc(NULL, p->cap);
    return p;
}

void fr_fastq_destroy(fr_fastq_t *p)
{
    if (!p) return;
    free(p->buf); free(p->seq_scr); free(p->qual_scr);
    free(p);
}

/* Pull more bytes in: drop the consumed prefix [0, beg), then read into the
 * tail, growing the buffer if a single record fills it. Returns bytes read
 * (0 at EOF, -1 on decode error). */
static int fr_refill(fr_fastq_t *p)
{
    if (p->beg > 0) {
        memmove(p->buf, p->buf + p->beg, p->end - p->beg);
        p->end -= p->beg;
        p->beg  = 0;
    }
    if (p->end == p->cap) {        /* record larger than the buffer: grow */
        p->cap *= 2;
        p->buf  = (unsigned char *)xrealloc(p->buf, p->cap);
    }
    size_t want = p->cap - p->end;
    if (want > FR_FASTQ_READ_MAX) want = FR_FASTQ_READ_MAX;
    int got = fast_reader_read(p->fr, p->buf + p->end, (int)want);
    if (got < 0) { p->err = 1; return -1; }
    if (got == 0) p->eof = 1;
    else          p->end += (size_t)got;
    return got;
}

/* Append n bytes to a scratch buffer, growing as needed. */
static void scr_put(unsigned char **bufp, size_t *capp, size_t *lenp,
                    const unsigned char *src, size_t n)
{
    if (*lenp + n > *capp) {
        size_t ncap = *capp ? *capp : 256;
        while (ncap < *lenp + n) ncap *= 2;
        *bufp = (unsigned char *)xrealloc(*bufp, ncap);
        *capp = ncap;
    }
    memcpy(*bufp + *lenp, src, n);
    *lenp += n;
}

/* Length of a line's content with a single trailing '\r' dropped, matching
 * kseq's KS_SEP_LINE rule: strip only when the copied run (which includes the
 * '\r') is longer than one byte. */
static size_t line_content_len(const unsigned char *s, size_t raw_len)
{
    if (raw_len > 1 && s[raw_len - 1] == '\r') return raw_len - 1;
    return raw_len;
}

/* Index of the next '\n' in b[from..end), or `end` if none. memchr is the
 * libc-vectorised newline scan (NEON/SSE/AVX) — the sequence and quality lines
 * are the dominant scan cost, so this is where the speedup lives. */
static inline size_t find_nl(const unsigned char *b, size_t from, size_t end)
{
    const unsigned char *nl = (const unsigned char *)memchr(b + from, '\n', end - from);
    return nl ? (size_t)(nl - b) : end;
}

/* Attempt to parse one record entirely out of the resident buffer.
 *   1           -> *rec filled, *new_beg = byte after the record
 *   0           -> clean EOF
 *   -2          -> malformed (matches kseq_read's -2)
 *   FR_UNDERRUN -> the record is not fully buffered yet (caller refills) */
static int scan_record(fr_fastq_t *p, fr_fastq_rec_t *rec, size_t *new_beg)
{
    const unsigned char *b = p->buf;
    const size_t end = p->end;

    /* Phase 1: skip to a header line. Junk/blank lines before it are consumed
     * for good, so advance beg past them. */
    while (p->beg < end && b[p->beg] != '>' && b[p->beg] != '@') p->beg++;
    if (p->beg >= end) return p->eof ? 0 : FR_UNDERRUN;

    const int is_fastq_marker = (b[p->beg] == '@');
    size_t i = p->beg + 1;        /* first byte after the header char */

    /* name: up to the first whitespace (kseq KS_SEP_SPACE). */
    size_t name_s = i;
    while (i < end && !isspace(b[i])) i++;
    if (i >= end && !p->eof) return FR_UNDERRUN;
    size_t name_e = i;

    /* comment: present iff the name stopped at a non-newline whitespace. */
    size_t com_s = name_e, com_e = name_e;
    size_t after_header;          /* first byte of the line below the header */
    if (i < end && b[i] != '\n') {
        size_t cs = i + 1, j = find_nl(b, cs, end);
        if (j >= end && !p->eof) return FR_UNDERRUN;
        com_s = cs;
        com_e = cs + line_content_len(b + cs, j - cs);
        after_header = (j < end) ? j + 1 : end;
    } else {
        after_header = (i < end) ? i + 1 : end;
    }

    /* seq: sequence lines (newlines stripped) until the next header or '+'. */
    size_t k = after_header;
    size_t seq_s = after_header, seq_e = after_header;  /* single-line fast path */
    size_t seq_len = 0; int seq_lines = 0;
    unsigned char term = 0;       /* 0 = EOF, else '>' '@' '+' */
    for (;;) {
        if (k >= end) { if (p->eof) break; else return FR_UNDERRUN; }
        unsigned char c = b[k];
        if (c == '>' || c == '@' || c == '+') { term = c; break; }
        if (c == '\n') { k++; continue; }        /* skip empty lines */
        size_t ls = k;
        k = find_nl(b, k, end);
        if (k >= end && !p->eof) return FR_UNDERRUN;
        size_t clen = line_content_len(b + ls, k - ls);
        if (seq_lines == 0) { seq_s = ls; seq_e = ls + clen; }
        else {
            if (seq_lines == 1) {                /* promote first line to scratch */
                size_t l = 0;
                scr_put(&p->seq_scr, &p->seq_scr_cap, &l, b + seq_s, seq_e - seq_s);
                seq_len = l;
            }
            scr_put(&p->seq_scr, &p->seq_scr_cap, &seq_len, b + ls, clen);
        }
        if (seq_lines == 0) seq_len = clen;
        seq_lines++;
        if (k < end) k++;                        /* step past '\n' */
    }

    rec->name    = (const char *)(b + name_s);
    rec->name_l  = name_e - name_s;
    rec->comment = (const char *)(b + com_s);
    rec->comment_l = com_e - com_s;
    if (seq_lines <= 1) { rec->seq = (const char *)(b + seq_s); }
    else                { rec->seq = (const char *)p->seq_scr; }
    rec->seq_l = seq_len;

    if (term != '+') {            /* FASTA (or EOF): no quality */
        (void)is_fastq_marker;
        rec->qual = NULL; rec->qual_l = 0;
        *new_beg = k;             /* leave the next header unconsumed */
        return 1;
    }

    /* skip the '+' separator line. */
    size_t pl = find_nl(b, k, end);
    if (pl >= end && !p->eof) return FR_UNDERRUN;
    if (pl >= end) return -2;     /* '+' line with no newline and no quality */
    size_t m = pl + 1;

    /* qual: quality lines accumulated until they reach the sequence length
     * (kseq reads at least one line, then stops once qual.l >= seq.l). */
    size_t qual_s = m, qual_e = m, qual_len = 0; int qual_lines = 0;
    for (;;) {
        if (m >= end) { if (p->eof) break; else return FR_UNDERRUN; }
        size_t ls = m;
        m = find_nl(b, m, end);
        if (m >= end && !p->eof) return FR_UNDERRUN;
        size_t clen = line_content_len(b + ls, m - ls);
        if (qual_lines == 0) { qual_s = ls; qual_e = ls + clen; }
        else {
            if (qual_lines == 1) {
                size_t l = 0;
                scr_put(&p->qual_scr, &p->qual_scr_cap, &l, b + qual_s, qual_e - qual_s);
                qual_len = l;
            }
            scr_put(&p->qual_scr, &p->qual_scr_cap, &qual_len, b + ls, clen);
        }
        if (qual_lines == 0) qual_len = clen;
        qual_lines++;
        if (m < end) m++;                        /* step past '\n' */
        if (qual_len >= seq_len) break;
    }
    if (qual_len != seq_len) return -2;          /* kseq: seq.l != qual.l */

    if (qual_lines <= 1) { rec->qual = (const char *)(b + qual_s); }
    else                 { rec->qual = (const char *)p->qual_scr; }
    rec->qual_l = qual_len;
    *new_beg = m;
    return 1;
}

int fr_fastq_next(fr_fastq_t *p, fr_fastq_rec_t *rec)
{
    for (;;) {
        size_t new_beg = p->beg;
        int r = scan_record(p, rec, &new_beg);
        if (r != FR_UNDERRUN) {
            if (r == 1) p->beg = new_beg;
            return r;
        }
        if (fr_refill(p) < 0) return -2;         /* decode error: stop the batch */
        /* on EOF, scan_record will now resolve (no more underruns possible). */
    }
}
