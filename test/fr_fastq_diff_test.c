/* Differential test: fr_fastq must parse byte-identically to kseq.
 *
 * kseq is the oracle. For each payload in a deliberately adversarial corpus
 * (single/multi-line FASTQ and FASTA, comments, CRLF, blank lines, missing
 * trailing newline, leading junk, truncated quality), we parse the same bytes
 * through kseq and through fr_fastq and assert that every record's
 * name/comment/seq/qual (length AND bytes) matches, and that both stop after
 * the same number of records (so the -2 malformed cut lines up too).
 *
 * Build (standalone — only fr_fastq + fast_reader, no bwa-mem3 link):
 *   cc -I src test/fr_fastq_diff_test.c src/fr_fastq.c src/fast_reader.c \
 *      -lz -ldeflate -o fr_fastq_diff_test
 */
#include "fast_reader.h"
#include "fr_fastq.h"
#include "kseq.h"

#include <assert.h>
#include <fcntl.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

KSEQ_INIT(fast_reader_t *, fast_reader_read)

/* fast_reader.c is instrumented with stage_prof hooks; stub them out so this
 * standalone parser test links without the bwa-mem3 object graph. sp_enabled()
 * returning 0 disables every timing path inside fast_reader.c. */
int    sp_enabled(void) { return 0; }
double sp_wall(void) { return 0.0; }
void   sp_read_add(int which, double seconds) { (void)which; (void)seconds; }
void   sp_read_bytes(long fd_bytes, long bgzf_blocks) { (void)fd_bytes; (void)bgzf_blocks; }

static int g_fail = 0;

/* A field copied out of a parser (both kseq and fr_fastq reuse their buffers
 * across records, so each record's bytes must be snapshotted immediately). */
typedef struct { char *s; size_t l; int present; } field_t;
typedef struct { field_t name, comment, seq, qual; } rec_t;

static field_t dup_field(const char *s, size_t l, int present)
{
    field_t f;
    f.l = l; f.present = present;
    f.s = (char *)malloc(l + 1);
    assert(f.s != NULL);
    if (l) memcpy(f.s, s, l);
    f.s[l] = '\0';
    return f;
}

static void free_rec(rec_t *r)
{
    free(r->name.s); free(r->comment.s); free(r->seq.s); free(r->qual.s);
}

/* Double the record array with overflow- and realloc-failure checks. Both
 * parsers grow `recs` identically, so this keeps the growth path safe in one
 * place: a doubling that wraps int, or a byte count that exceeds SIZE_MAX, or a
 * realloc that fails, all abort rather than under-allocate or leak the original. */
static rec_t *grow_recs(rec_t *recs, int *cap)
{
    int next = (*cap == 0) ? 16 : (*cap * 2);
    if (next <= *cap) abort();                            /* int overflow */
    if ((size_t)next > SIZE_MAX / sizeof(*recs)) abort(); /* size overflow */
    rec_t *tmp = (rec_t *)realloc(recs, (size_t)next * sizeof(*recs));
    if (!tmp) abort();                                    /* realloc failure */
    *cap = next;
    return tmp;
}

/* Parse via kseq. Mirrors bseq_read_fast's `kseq_read(...) >= 0` stop. */
static rec_t *parse_kseq(const char *path, int *n_out)
{
    int fd = open(path, O_RDONLY);
    assert(fd >= 0);
    const char *err = NULL;
    fast_reader_t *fr = fast_reader_dopen(fd, &err);
    assert(fr != NULL);
    kseq_t *ks = kseq_init(fr);

    int cap = 16, n = 0;
    rec_t *recs = (rec_t *)malloc((size_t)cap * sizeof(*recs));
    assert(recs != NULL);
    while (kseq_read(ks) >= 0) {
        if (n == cap) recs = grow_recs(recs, &cap);
        recs[n].name    = dup_field(ks->name.s, ks->name.l, 1);
        recs[n].comment = dup_field(ks->comment.s ? ks->comment.s : "", ks->comment.l, ks->comment.l > 0);
        recs[n].seq     = dup_field(ks->seq.s, ks->seq.l, 1);
        recs[n].qual    = dup_field(ks->qual.s ? ks->qual.s : "", ks->qual.l, ks->qual.l > 0);
        n++;
    }
    kseq_destroy(ks);
    fast_reader_close(fr);
    *n_out = n;
    return recs;
}

/* Parse via fr_fastq. Mirrors the same stop (next() != 1). */
static rec_t *parse_fr(const char *path, int *n_out)
{
    int fd = open(path, O_RDONLY);
    assert(fd >= 0);
    const char *err = NULL;
    fast_reader_t *fr = fast_reader_dopen(fd, &err);
    assert(fr != NULL);
    fr_fastq_t *p = fr_fastq_init(fr);

    int cap = 16, n = 0;
    rec_t *recs = (rec_t *)malloc((size_t)cap * sizeof(*recs));
    assert(recs != NULL);
    fr_fastq_rec_t rec;
    while (fr_fastq_next(p, &rec) == 1) {
        if (n == cap) recs = grow_recs(recs, &cap);
        recs[n].name    = dup_field(rec.name, rec.name_l, 1);
        recs[n].comment = dup_field(rec.comment ? rec.comment : "", rec.comment_l, rec.comment_l > 0);
        recs[n].seq     = dup_field(rec.seq, rec.seq_l, 1);
        recs[n].qual    = dup_field(rec.qual ? rec.qual : "", rec.qual_l, rec.qual_l > 0);
        n++;
    }
    fr_fastq_destroy(p);
    fast_reader_close(fr);
    *n_out = n;
    return recs;
}

static int field_eq(const field_t *a, const field_t *b)
{
    /* Only length + bytes matter for byte-identical bseq1_t: a 0-length field
     * and an absent field both produce a NULL bseq1 pointer downstream. */
    if (a->l != b->l) return 0;
    return memcmp(a->s, b->s, a->l) == 0;
}

static void report_field(const char *which, int rec, const field_t *k, const field_t *f)
{
    if (field_eq(k, f)) return;
    fprintf(stderr, "    rec %d %-7s MISMATCH: kseq(l=%zu)='%s'  fr(l=%zu)='%s'\n",
            rec, which, k->l, k->s, f->l, f->s);
}

static void run_case(const char *name, const char *payload, size_t len)
{
    char path[] = "/tmp/fr_fastq_diff_XXXXXX";
    int fd = mkstemp(path);
    assert(fd >= 0);
    { size_t off = 0; while (off < len) { ssize_t w = write(fd, payload + off, len - off); assert(w > 0); off += (size_t)w; } }
    close(fd);

    int nk = 0, nf = 0;
    rec_t *k = parse_kseq(path, &nk);
    rec_t *f = parse_fr(path, &nf);

    int ok = 1;
    if (nk != nf) { ok = 0; fprintf(stderr, "    record count: kseq=%d fr=%d\n", nk, nf); }
    int n = nk < nf ? nk : nf;
    for (int i = 0; i < n; i++) {
        int rec_ok = field_eq(&k[i].name, &f[i].name) && field_eq(&k[i].comment, &f[i].comment)
                  && field_eq(&k[i].seq, &f[i].seq) && field_eq(&k[i].qual, &f[i].qual);
        if (!rec_ok) {
            ok = 0;
            report_field("name", i, &k[i].name, &f[i].name);
            report_field("comment", i, &k[i].comment, &f[i].comment);
            report_field("seq", i, &k[i].seq, &f[i].seq);
            report_field("qual", i, &k[i].qual, &f[i].qual);
        }
    }
    fprintf(stderr, "%s: %s (%d records)\n", ok ? "ok  " : "FAIL", name, nk);
    if (!ok) g_fail = 1;

    for (int i = 0; i < nk; i++) free_rec(&k[i]);
    for (int i = 0; i < nf; i++) free_rec(&f[i]);
    free(k); free(f);
    unlink(path);
}

typedef struct { const char *name; const char *payload; size_t len; } case_t;

#define LIT(name, lit) { (name), (lit), sizeof(lit) - 1 }

static const case_t CASES[] = {
    LIT("basic single-line FASTQ",
        "@r1\nACGTACGT\n+\nIIIIIIII\n"
        "@r2\nTTTTGGGG\n+\nFFFFFFFF\n"
        "@r3\nCCCCAAAA\n+\n########\n"),
    LIT("readno suffixes /1 /2",
        "@read1/1\nACGT\n+\nIIII\n@read1/2\nTTGG\n+\nFFFF\n"),
    LIT("with comment",
        "@r1 here is a comment\nACGT\n+\nIIII\n@r2 1:N:0:ACGT\nTTGG\n+\nFFFF\n"),
    LIT("multiple spaces before comment", "@r1   spaced comment\nACGT\n+\nIIII\n"),
    LIT("tab before comment",             "@r1\tfield2\tfield3\nACGT\n+\nIIII\n"),
    LIT("empty comment (trailing space)", "@r1 \nACGT\n+\nIIII\n"),
    LIT("FASTA single-line",              ">s1\nACGTACGT\n>s2\nTTTTGGGG\n"),
    LIT("FASTA multi-line seq",           ">s1 desc\nACGT\nTTGG\nCCAA\n>s2\nGGGG\nTTTT\n"),
    LIT("multi-line FASTQ seq and qual",  "@r1\nACGT\nTTGG\n+\nIIII\nFFFF\n@r2\nCCCC\n+\nHHHH\n"),
    LIT("CRLF line endings",              "@r1\r\nACGT\r\n+\r\nIIII\r\n@r2\r\nTTGG\r\n+\r\nFFFF\r\n"),
    LIT("blank lines between records",    "@r1\nACGT\n+\nIIII\n\n@r2\nTTGG\n+\nFFFF\n"),
    LIT("no trailing newline",            "@r1\nACGT\n+\nIIII\n@r2\nTTGG\n+\nFFFF"),
    LIT("leading junk before first header", "; comment line\n# another\n@r1\nACGT\n+\nIIII\n"),
    LIT("lowercase and N in seq",         "@r1\nacgtNNNNACGT\n+\nIIIIIIIIIIII\n"),
    LIT("qual chars include @ and +",     "@r1\nACGTACGT\n+\n@@++!!~~\n@r2\nTTTTGGGG\n+\n++++@@@@\n"),
    LIT("plus line with trailing text",   "@r1\nACGT\n+r1 same name\nIIII\n"),
    LIT("truncated final qual (kseq -2)", "@r1\nACGT\n+\nIIII\n@r2\nACGTACGT\n+\nIIII\n"),
    LIT("single record, no comment, exact", "@only\nA\n+\nI\n"),
};

/* Deterministic LCG so the fuzz corpus is reproducible across runs. */
static unsigned long g_rng = 0x9e3779b97f4a7c15UL;
static unsigned rnd(unsigned mod) { g_rng = g_rng * 6364136223846793005UL + 1; return (unsigned)((g_rng >> 33) % mod); }

/* Build a payload of random-but-valid FASTQ/FASTA records: varied name, seq and
 * qual lengths (qual always matches seq), optional comments, occasional FASTA
 * records. Exercises every field-length boundary against kseq. */
static size_t make_fuzz_payload(char *out, size_t cap)
{
    static const char BASES[] = "ACGTNacgtn";
    static const char QUALS[] = "!#$%&()*+,-./0123456789:;<=>?@ABCDEFGHIIIII";
    size_t n = 0;
    while (n + 1024 < cap) {
        int fasta = (rnd(5) == 0);
        out[n++] = fasta ? '>' : '@';
        int namelen = 1 + (int)rnd(20);
        for (int j = 0; j < namelen; j++) out[n++] = (char)('a' + rnd(26));
        if (rnd(2)) { out[n++] = '/'; out[n++] = (char)('1' + rnd(2)); } /* sometimes /1 or /2 */
        if (rnd(3) == 0) { out[n++] = ' '; int cl = (int)rnd(15); for (int j = 0; j < cl; j++) out[n++] = (char)('A' + rnd(40)); }
        out[n++] = '\n';
        int seqlen = (int)rnd(40);
        for (int j = 0; j < seqlen; j++) out[n++] = BASES[rnd(sizeof BASES - 1)];
        out[n++] = '\n';
        if (!fasta) {
            out[n++] = '+'; out[n++] = '\n';
            for (int j = 0; j < seqlen; j++) out[n++] = QUALS[rnd(sizeof QUALS - 1)];
            out[n++] = '\n';
        }
    }
    return n;
}

int main(void)
{
    const int NCASES = (int)(sizeof(CASES) / sizeof(CASES[0]));

    /* A moderate payload spanning many records, used to stress buffer refills
     * (small enough to stay fast even at a 1-byte buffer). */
    size_t big_cap = 200u << 10;
    char *big = (char *)malloc(big_cap);
    assert(big != NULL);
    size_t big_n = 0;
    for (int i = 0; big_n + 256 < big_cap; i++)
        big_n += (size_t)snprintf(big + big_n, big_cap - big_n,
            "@read%d/1 comment %d here\nACGTACGTNNNNacgtACGT\n+\nIIIIIIIIIIIIIIIIIIII\n", i, i);

    /* Run the whole corpus at the default buffer, then at a ladder of tiny
     * buffer sizes (via the FR_FASTQ_TEST_BUFSZ hook) so every record straddles
     * the buffer end and the underrun/compact/grow path is exercised at every
     * byte alignment. */
    const long bufszs[] = { 0, 1, 2, 3, 5, 7, 16, 64, 4096 };
    for (size_t bi = 0; bi < sizeof(bufszs) / sizeof(bufszs[0]); bi++) {
        char val[32];
        if (bufszs[bi] == 0) unsetenv("FR_FASTQ_TEST_BUFSZ");
        else { snprintf(val, sizeof val, "%ld", bufszs[bi]); setenv("FR_FASTQ_TEST_BUFSZ", val, 1); }
        fprintf(stderr, "--- buffer size: %s ---\n", bufszs[bi] ? val : "default(1MiB)");
        for (int i = 0; i < NCASES; i++)
            run_case(CASES[i].name, CASES[i].payload, CASES[i].len);
        run_case("refill-stress payload", big, big_n);
    }

    free(big);

    /* Deterministic fuzz: many randomly-sized records, diffed at the default
     * buffer and a couple of tiny ones so field-length boundaries land at every
     * offset relative to a refill. */
    {
        size_t fcap = 256u << 10;
        char *fz = (char *)malloc(fcap);
        assert(fz != NULL);
        size_t fn = make_fuzz_payload(fz, fcap);
        const long fbufs[] = { 0, 1, 3, 7, 64 };
        for (size_t bi = 0; bi < sizeof(fbufs) / sizeof(fbufs[0]); bi++) {
            char val[32];
            if (fbufs[bi] == 0) unsetenv("FR_FASTQ_TEST_BUFSZ");
            else { snprintf(val, sizeof val, "%ld", fbufs[bi]); setenv("FR_FASTQ_TEST_BUFSZ", val, 1); }
            fprintf(stderr, "--- fuzz @ buffer %s ---\n", fbufs[bi] ? val : "default");
            run_case("fuzz corpus", fz, fn);
        }
        free(fz);
    }

    fprintf(stderr, "\n%s\n", g_fail ? "DIFFERENTIAL TEST FAILED" : "ALL CASES IDENTICAL");
    return g_fail ? 1 : 0;
}
