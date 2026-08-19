/* Regression: bseq_read_fast must not write past the read buffer in paired-end
 * mode when the initial capacity estimate is odd.
 *
 * The batch reader sizes its bseq1_t buffer once from the first read's length:
 *   m = chunk_size / r1.seq_l + 256   (floored at 256)
 * and, in paired-end mode, writes TWO records per loop iteration (seqs[n] and
 * seqs[n+1]) while checking capacity only once at the top of the loop. `n`
 * therefore only ever takes even values in PE mode, but `m` is ODD whenever
 * chunk_size / r1.seq_l is odd. When n reaches m-1 (even, since m is odd) the
 * guard `n >= m` is false, so the iteration writes seqs[m-1] AND seqs[m] — and
 * seqs[m] is one past the end of an m-element buffer.
 *
 * This test crafts exactly that shape: a long first read fixes a small, odd
 * quotient (chunk_size / 1000 = 3 -> m = 259), and all subsequent reads are
 * tiny so the byte budget (chunk_size) is not reached until well after n passes
 * m-1. Built with ASan, the overflowing write to seqs[259] aborts on the
 * pre-fix code and passes once the capacity guard accounts for the paired write.
 *
 * Build (standalone — only the reader TUs, no bwa-mem3 link):
 *   cc -I src -fsanitize=address test/bseq_read_pe_oob_test.c \
 *      src/fast_reader_bseq.c src/fr_fastq.c src/fast_reader.c \
 *      -lz -ldeflate -o bseq_read_pe_oob_test
 */
#include "bwa.h"
#include "fast_reader.h"
#include "fast_reader_bseq.h"

#include <assert.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* fast_reader.c / fast_reader_bseq.c are instrumented with stage_prof hooks;
 * stub them so this standalone test links without the bwa-mem3 object graph.
 * sp_enabled() == 0 disables every timing path. err_fatal is bwa's abort-loudly
 * helper (utils.c); stub it to a plain abort so the test needs no libbwa link. */
int    sp_enabled(void) { return 0; }
double sp_wall(void) { return 0.0; }
void   sp_read_add(int which, double seconds) { (void)which; (void)seconds; }
void   sp_read_bytes(long fd_bytes, long bgzf_blocks) { (void)fd_bytes; (void)bgzf_blocks; }
void   sp_read_get(double *disk, double *decomp, double *parse)
{ if (disk) *disk = 0; if (decomp) *decomp = 0; if (parse) *parse = 0; }
void   err_fatal(const char *header, const char *fmt, ...)
{ (void)header; (void)fmt; fprintf(stderr, "err_fatal\n"); abort(); }

/* Write a FASTQ file whose first record has `first_len` bases and the remaining
 * `n_records - 1` records have `rest_len` bases each. Returns a heap path the
 * caller must unlink()+free(). */
static char *write_fastq(int n_records, int first_len, int rest_len)
{
    char *path = strdup("/tmp/bseq_oob_XXXXXX");
    int fd = mkstemp(path);
    assert(fd >= 0);
    FILE *f = fdopen(fd, "w");
    assert(f != NULL);
    char *seq = (char *)malloc((size_t)first_len + 1);
    char *qual = (char *)malloc((size_t)first_len + 1);
    assert(seq && qual);
    for (int i = 0; i < n_records; i++) {
        int len = (i == 0) ? first_len : rest_len;
        memset(seq, 'A', (size_t)len);
        memset(qual, 'I', (size_t)len);
        seq[len] = qual[len] = '\0';
        fprintf(f, "@r%d\n%s\n+\n%s\n", i, seq, qual);
    }
    free(seq);
    free(qual);
    fclose(f);
    return path;
}

int main(void)
{
    /* quotient = chunk_size / first_len = 3000 / 1000 = 3 (odd) -> initial
     * m = 3 + 256 = 259 (odd). With 1-base tail reads the byte budget (3000)
     * is not reached until ~3000 records, so n climbs past m-1 = 258 first and
     * the paired write at n=258 touches seqs[259]. 200 records per mate is more
     * than the 130 pairs needed to reach n=258. */
    const int64_t chunk_size = 3000;
    const int n_records = 200, first_len = 1000, rest_len = 1;

    char *p1 = write_fastq(n_records, first_len, rest_len);
    char *p2 = write_fastq(n_records, first_len, rest_len);

    int fd1 = open(p1, O_RDONLY), fd2 = open(p2, O_RDONLY);
    assert(fd1 >= 0 && fd2 >= 0);
    const char *err = NULL;
    fast_reader_t *fr1 = fast_reader_dopen(fd1, &err);
    fast_reader_t *fr2 = fast_reader_dopen(fd2, &err);
    assert(fr1 && fr2);
    void *ks1 = fast_kseq_init(fr1), *ks2 = fast_kseq_init(fr2);

    int n = 0;
    int64_t size = 0;
    read_arena_t *arena = NULL;
    bseq1_t *seqs = bseq_read_fast(chunk_size, &n, ks1, ks2, &size, &arena, 1);

    /* Reaching here without an ASan abort means the paired write stayed in
     * bounds. Sanity-check the batch is non-trivial (n crossed the odd initial
     * capacity) so the test can't pass by reading nothing. */
    assert(n > 258 && "batch too small to exercise the odd-capacity boundary");

    for (int i = 0; i < n; i++) free(seqs[i].comment);
    free(seqs);
    if (arena) read_arena_destroy(arena);
    fast_kseq_destroy(ks1);
    fast_kseq_destroy(ks2);
    fast_reader_close(fr1);
    fast_reader_close(fr2);
    unlink(p1); unlink(p2);
    free(p1); free(p2);

    fprintf(stderr, "ok: bseq_read_fast PE odd-capacity boundary stayed in bounds (n=%d)\n", n);
    return 0;
}
