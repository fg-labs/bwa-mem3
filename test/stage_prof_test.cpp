/* Unit test for the pure helpers of stage_prof (stats + clocks + NaN init).
 * Builds standalone: no htslib, no index, no pipeline. */
#include "../src/stage_prof.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Copy the 0-based tab-delimited field `col` of `line` into `out`. Unlike
 * strtok, this preserves empty fields (consecutive tabs), which is essential
 * for distinguishing a blank TSV cell from a "0". Trailing newline stripped.
 * Returns the number of fields seen (so callers can detect a short row). */
static int tsv_field(const char *line, int col, char *out, size_t out_sz) {
    out[0] = '\0';
    int idx = 0;
    const char *p = line;
    while (1) {
        const char *tab = strchr(p, '\t');
        const char *end = tab ? tab : p + strlen(p);
        if (idx == col) {
            size_t n = (size_t)(end - p);
            while (n > 0 && (p[n - 1] == '\n' || p[n - 1] == '\r')) n--;   /* strip EOL on last field */
            if (n >= out_sz) n = out_sz - 1;
            memcpy(out, p, n);
            out[n] = '\0';
        }
        idx++;
        if (!tab) break;
        p = tab + 1;
    }
    return idx;
}

/* Read the aggregate ("ALL") row of a stage_prof TSV and return its tab field
 * at the given 0-based column index into `out` (empty string if the cell is
 * blank). Returns 1 on success, 0 if no ALL row was found. The aggregate row
 * carries "ALL" in the chunk column (index 8). */
static int agg_field(const char *path, int col, char *out, size_t out_sz) {
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    /* generous: the input column alone can be up to 4096 chars in real output */
    char line[16384], chunk[64];
    int found = 0;
    while (fgets(line, sizeof line, f)) {
        if (tsv_field(line, 8, chunk, sizeof chunk) > 8 && strcmp(chunk, "ALL") == 0) {
            tsv_field(line, col, out, out_sz);
            found = 1;
            break;
        }
    }
    fclose(f);
    return found;
}

int main(void) {
    /* sp_chunk_init sets the maybe-N/A doubles to NaN */
    prof_chunk_t c;
    sp_chunk_init(&c);
    assert(isnan(c.read_parse));
    assert(isnan(c.thr_busy_mean));
    assert(isnan(c.write_compress));
    assert(c.n_reads == 0 && c.write_bytes == 0);

    /* sp_thread_stats over {1,3,3,1}: min 1, max 3, mean 2, stdev 1 (population) */
    double busy[4] = {1.0, 3.0, 3.0, 1.0};
    sp_thread_stats(&c, busy, 4);
    assert(c.thr_busy_min == 1.0);
    assert(c.thr_busy_max == 3.0);
    assert(fabs(c.thr_busy_mean - 2.0) < 1e-12);
    assert(fabs(c.thr_busy_stdev - 1.0) < 1e-9);

    /* single-element stats are degenerate but well-defined */
    prof_chunk_t c1; sp_chunk_init(&c1);
    double one[1] = {5.0};
    sp_thread_stats(&c1, one, 1);
    assert(c1.thr_busy_min == 5.0 && c1.thr_busy_max == 5.0);
    assert(c1.thr_busy_mean == 5.0 && c1.thr_busy_stdev == 0.0);

    /* clocks: monotonic and nonnegative */
    double w0 = sp_wall(), w1 = sp_wall();
    assert(w1 >= w0);
    assert(sp_thread_cpu() >= 0.0);

    /* profiling is off until sp_init with a real path */
    assert(sp_enabled() == 0);
    sp_init("", "t", "v", "x86_64", 4, "sam", -1, "in");   /* empty path -> stays off */
    assert(sp_enabled() == 0);

    /* Aggregate N/A semantics: when every chunk reports read_bytes_in/bgzf_blocks
     * as N/A (-1), the aggregate row must stay blank, not collapse to a false 0. */
    const char *tsv = "/tmp/stage_prof_test.tsv";
    remove(tsv);
    sp_init(tsv, "t", "v", "x86_64", 4, "sam", -1, "in");
    assert(sp_enabled() == 1);
    prof_chunk_t a; sp_chunk_init(&a); a.chunk = 0; a.n_reads = 10; a.n_bp = 1500;
    prof_chunk_t b; sp_chunk_init(&b); b.chunk = 1; b.n_reads = 20; b.n_bp = 3000;
    sp_add_chunk(&a);   /* both leave read_bytes_in/bgzf_blocks at the -1 N/A default */
    sp_add_chunk(&b);
    sp_finish(1.0, 0.0, 0.0);

    char field[64];
    assert(agg_field(tsv, 10, field, sizeof field) && strcmp(field, "30") == 0);  /* n_reads sums */
    assert(agg_field(tsv, 16, field, sizeof field) && field[0] == '\0');           /* read_bytes_in blank */
    assert(agg_field(tsv, 17, field, sizeof field) && field[0] == '\0');           /* bgzf_blocks blank */
    remove(tsv);

    printf("stage_prof_test OK\n");
    return 0;
}
