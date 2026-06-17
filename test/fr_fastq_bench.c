/* Read+parse microbenchmark: kseq path vs fr_fastq path, on a real FASTQ.
 *
 * Isolates exactly what the single-copy parser changed — decompress (identical
 * in both modes) + tokenize + field copy into bseq1-style allocations — with no
 * index load or alignment, so it runs in ~1s/rep and the delta between modes is
 * the parse change alone.
 *
 *   mode kseq    : kseq_read + (trim_readno + per-field copy), then free
 *   mode frfastq : fr_fastq_next + (trim_readno + per-field copy), then free
 *
 * Build:
 *   cc -O2 -Wall -Wextra -Isrc test/fr_fastq_bench.c src/fr_fastq.c src/fast_reader.c \
 *      -lz -ldeflate -o fr_fastq_bench
 * Run:  ./fr_fastq_bench {kseq|frfastq} reads.fq.gz [reps]
 */
#include "fast_reader.h"
#include "fr_fastq.h"
#include "kseq.h"

#include <ctype.h>
#include <fcntl.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

KSEQ_INIT(fast_reader_t *, fast_reader_read)

int    sp_enabled(void) { return 0; }
double sp_wall(void) { return 0.0; }
void   sp_read_add(int w, double s) { (void)w; (void)s; }
void   sp_read_bytes(long a, long b) { (void)a; (void)b; }

static double now(void)
{
    struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + (double)t.tv_nsec * 1e-9;
}

/* Mirror the real per-field copy: malloc(l+1)+memcpy, no strlen. */
static char *dup_field(const char *s, size_t l)
{
    char *d = (char *)malloc(l + 1);
    if (!d) abort();
    memcpy(d, s, l); d[l] = '\0';
    return d;
}
static size_t trim_readno_len(const char *name, size_t l)
{
    if (l > 2 && name[l - 2] == '/' && isdigit((unsigned char)name[l - 1])) return l - 2;
    return l;
}

/* Consume one batch's worth of allocations exactly as the pipeline would: name,
 * comment, seq, qual malloc'd then freed. Returns reads + bases via *nb. */
static long run_kseq(const char *path, long *bases)
{
    int fd = open(path, O_RDONLY);
    if (fd < 0) { perror("open"); exit(2); }
    const char *err = NULL;
    fast_reader_t *fr = fast_reader_dopen(fd, &err);
    if (!fr) { fprintf(stderr, "open fail: %s\n", err ? err : "?"); exit(2); }
    kseq_t *ks = kseq_init(fr);
    long n = 0, bp = 0;
    while (kseq_read(ks) >= 0) {
        size_t nl = trim_readno_len(ks->name.s, ks->name.l);
        char *name = dup_field(ks->name.s, nl);
        char *comment = ks->comment.l ? dup_field(ks->comment.s, ks->comment.l) : NULL;
        char *seq = dup_field(ks->seq.s, ks->seq.l);
        char *qual = ks->qual.l ? dup_field(ks->qual.s, ks->qual.l) : NULL;
        bp += (long)ks->seq.l; n++;
        free(name); free(comment); free(seq); free(qual);
    }
    kseq_destroy(ks); fast_reader_close(fr);
    *bases = bp; return n;
}

/* malloc(l+1) + free, no memcpy: isolates allocation cost. A single touch
 * keeps the optimizer from eliding the allocation. */
static volatile char g_sink;
static void alloc_field(size_t l)
{
    char *d = (char *)malloc(l + 1);
    if (!d) abort();
    d[0] = (char)l; g_sink = d[0];
    free(d);
}

/* level: 0 = scan only, 1 = scan + alloc (no copy), 2 = scan + alloc + copy. */
static long run_frfastq(const char *path, long *bases, int level)
{
    int fd = open(path, O_RDONLY);
    if (fd < 0) { perror("open"); exit(2); }
    const char *err = NULL;
    fast_reader_t *fr = fast_reader_dopen(fd, &err);
    if (!fr) { fprintf(stderr, "open fail: %s\n", err ? err : "?"); exit(2); }
    fr_fastq_t *p = fr_fastq_init(fr);
    long n = 0, bp = 0;
    fr_fastq_rec_t r;
    while (fr_fastq_next(p, &r) == 1) {
        size_t nl = trim_readno_len(r.name, r.name_l);
        if (level == 2) {
            char *name = dup_field(r.name, nl);
            char *comment = r.comment_l ? dup_field(r.comment, r.comment_l) : NULL;
            char *seq = dup_field(r.seq, r.seq_l);
            char *qual = r.qual_l ? dup_field(r.qual, r.qual_l) : NULL;
            free(name); free(comment); free(seq); free(qual);
        } else if (level == 1) {
            alloc_field(nl);
            if (r.comment_l) alloc_field(r.comment_l);
            alloc_field(r.seq_l);
            if (r.qual_l) alloc_field(r.qual_l);
        } else {
            g_sink = (char)(nl ^ r.seq_l ^ r.qual_l); /* force the slices to be read */
        }
        bp += (long)r.seq_l; n++;
    }
    fr_fastq_destroy(p); fast_reader_close(fr);
    *bases = bp; return n;
}

int main(int argc, char **argv)
{
    if (argc < 3) { fprintf(stderr, "usage: %s {kseq|frfastq} reads.fq[.gz] [reps]\n", argv[0]); return 2; }
    const char *mode = argv[1], *path = argv[2];
    int reps = 5;
    if (argc > 3) {
        char *end = NULL;
        long v = strtol(argv[3], &end, 10);
        if (!end || *end != '\0' || v < 1 || v > INT_MAX) {
            fprintf(stderr, "invalid reps: %s\n", argv[3]);
            return 2;
        }
        reps = (int)v;
    }
    /* modes: kseq | frfastq(=copy) | fr_scan | fr_alloc | fr_copy */
    int use_kseq = strcmp(mode, "kseq") == 0;
    int level = 2;
    if (strcmp(mode, "fr_scan") == 0)  level = 0;
    else if (strcmp(mode, "fr_alloc") == 0) level = 1;
    else if (strcmp(mode, "fr_copy") == 0 || strcmp(mode, "frfastq") == 0) level = 2;
    else if (!use_kseq) { fprintf(stderr, "bad mode\n"); return 2; }

    double best = 1e30; long n = 0, bp = 0;
    for (int i = 0; i < reps; i++) {
        double t0 = now();
        n = use_kseq ? run_kseq(path, &bp) : run_frfastq(path, &bp, level);
        double dt = now() - t0;
        if (dt < best) best = dt;
        fprintf(stderr, "  rep %d: %.4fs  (%.1f Mreads/s)\n", i, dt, n / dt / 1e6);
    }
    fprintf(stdout, "%-8s best=%.4fs  reads=%ld  bp=%ld  %.2f Mreads/s  %.1f MB/s\n",
            mode, best, n, bp, n / best / 1e6, bp / best / 1e6);
    return 0;
}
