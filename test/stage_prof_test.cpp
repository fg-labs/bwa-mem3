/* Unit test for the pure helpers of stage_prof (stats + clocks + NaN init).
 * Builds standalone: no htslib, no index, no pipeline. */
#include "../src/stage_prof.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>   /* dup/dup2/fileno for the stderr-capture in the OOM-drop check */

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

/* Verify the retained (non-aggregate) chunk rows of a stage_prof TSV are exactly
 * the generated set after a forced drop: chunk ids 0..n-1, each present once,
 * with n_reads (col 10) and n_bp (col 11) intact at 1. A bare aggregate-count
 * check (n_reads sums to 64) still passes if the drop path corrupts, replaces,
 * or duplicates a retained record while leaving the sum unchanged; asserting
 * per-row identity is what catches that. Skips the header (col 8 == "chunk")
 * and the aggregate row (col 8 == "ALL"). */
static void assert_retained_chunks(const char *path, int n) {
    FILE *f = fopen(path, "r");
    assert(f != NULL);
    char line[16384], chunk[64], reads[64], bp[64];
    int seen[64] = {0};
    assert(n <= (int)(sizeof seen / sizeof seen[0]));
    int rows = 0;
    while (fgets(line, sizeof line, f)) {
        if (tsv_field(line, 8, chunk, sizeof chunk) <= 8) continue;  /* short/blank line */
        if (strcmp(chunk, "chunk") == 0 || strcmp(chunk, "ALL") == 0) continue;
        char *end = NULL;
        long id = strtol(chunk, &end, 10);
        assert(end != chunk && *end == '\0');   /* chunk id parses cleanly */
        assert(id >= 0 && id < n);              /* within the generated range */
        assert(!seen[id]);                      /* no duplicate id */
        seen[id] = 1;
        assert(tsv_field(line, 10, reads, sizeof reads) > 10 && strcmp(reads, "1") == 0);  /* n_reads intact */
        assert(tsv_field(line, 11, bp, sizeof bp) > 11 && strcmp(bp, "1") == 0);           /* n_bp intact */
        rows++;
    }
    fclose(f);
    assert(rows == n);                          /* exactly n retained rows, no extras */
    for (int i = 0; i < n; i++) assert(seen[i]); /* every generated id present */
}

/* Run sp_finish(w,i,c) with stderr redirected to `errpath`, then restore the
 * real stderr so any later assert stays visible. Lets a caller inspect the
 * incomplete-report diagnostic sp_finish writes to stderr. */
static void finish_capturing(const char *errpath, double w, double i, double c) {
    fflush(stderr);
    int saved_err = dup(fileno(stderr));
    assert(saved_err >= 0);
    FILE *ef = freopen(errpath, "w", stderr);
    assert(ef != NULL);
    sp_finish(w, i, c);
    fflush(stderr);
    assert(dup2(saved_err, fileno(stderr)) >= 0);
    close(saved_err);
    clearerr(stderr);
}

/* Read the whole file at `path` into `buf` (NUL-terminated, truncated to fit). */
static void slurp(const char *path, char *buf, size_t bufsz) {
    FILE *f = fopen(path, "r");
    assert(f != NULL);
    size_t n = fread(buf, 1, bufsz - 1, f);
    buf[n] = '\0';
    fclose(f);
}

/* Fill the chunk buffer to its current capacity, then arm one forced realloc
 * failure so the next add crosses the grow boundary and is dropped. Leaves
 * g_n_dropped == 1 for the current run without relying on host OOM. */
static void force_one_drop(void) {
    /* g_cap seeds at 64 and doubles; once seeded (by any prior run in this
     * process) it stays >= 64, so filling to it then adding one more always
     * crosses the boundary. */
    prof_chunk_t g;
    for (int i = 0; i < 64; i++) {
        sp_chunk_init(&g); g.chunk = i; g.n_reads = 1; g.n_bp = 1;
        sp_add_chunk(&g);   /* all succeed: capacity is not yet exhausted */
    }
    sp_test_arm_realloc_fail(1);   /* the next growth realloc behaves as if it returned NULL */
    sp_chunk_init(&g); g.chunk = 64; g.n_reads = 1; g.n_bp = 1;
    sp_add_chunk(&g);   /* crosses the boundary -> forced failure -> chunk dropped */
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

    /* sp_init must reset run-scoped accumulators so a second run started in
     * the same process (a host embedding bwa-mem3 as a library, or -- as here
     * -- a second sp_init/sp_finish cycle in one test binary) reports only
     * its own totals, not a mix carried over from a prior run. This exercises
     * the chunk count and idle time directly; the OOM-drop counter resets via
     * the same statement (see sp_init) and is forced nonzero deterministically
     * in the dedicated OOM-drop block below. */
    const char *tsv2 = "/tmp/stage_prof_test2.tsv";
    const char *err1 = "/tmp/stage_prof_test1.err";
    const char *err2 = "/tmp/stage_prof_test2.err";
    char errbuf[1024];
    remove(tsv);
    remove(tsv2);
    remove(err1);
    remove(err2);

    /* Run 1: accrue idle time AND force a dropped chunk, so this run ends with
     * g_n_dropped > 0 and its sp_finish flags the report incomplete. Capturing
     * run 1's stderr is the independent oracle that the drop really happened. */
    sp_init(tsv, "t", "v", "x86_64", 4, "sam", -1, "in");
    force_one_drop();
    sp_add_idle(0, 5.0);   /* run 1 accrues idle time */
    finish_capturing(err1, 1.0, 0.0, 0.0);
    slurp(err1, errbuf, sizeof errbuf);
    assert(strstr(errbuf, "incomplete") != NULL);   /* run 1 dropped a chunk */

    /* Run 2: a fresh run in the same process adds one small chunk and neither
     * drops nor accrues idle. sp_init must reset the run-scoped accumulators
     * (chunk count, idle, AND the dropped-chunk counter) so run 2 reports only
     * its own totals -- not state carried over from run 1. */
    sp_init(tsv2, "t", "v", "x86_64", 4, "sam", -1, "in");   /* fresh run, same process */
    prof_chunk_t r2; sp_chunk_init(&r2); r2.chunk = 0; r2.n_reads = 7; r2.n_bp = 7;
    sp_add_chunk(&r2);
    finish_capturing(err2, 1.0, 0.0, 0.0);
    slurp(err2, errbuf, sizeof errbuf);

    /* g_n_dropped is reset by sp_init, so run 2 must NOT inherit run 1's drop:
     * its report is complete and its diagnostic mentions no dropped chunk. */
    assert(strstr(errbuf, "incomplete") == NULL);
    assert(strstr(errbuf, "dropped") == NULL);
    /* run 2's aggregate must reflect only its own chunk (7), not the 64 rows
     * that landed in run 1, inherited via an un-reset chunk count. */
    assert(agg_field(tsv2, 10, field, sizeof field) && strcmp(field, "7") == 0);
    /* run 2's idle-time column must be 0, not run 1's leftover 5.0. */
    assert(agg_field(tsv2, 31, field, sizeof field) && strcmp(field, "0.0000") == 0);
    remove(tsv);
    remove(tsv2);
    remove(err1);
    remove(err2);

    /* Deterministic OOM-drop coverage: force the buffer past capacity with one
     * armed realloc failure so the next add takes the drop path. The dropped
     * chunk must NOT appear in the report, the already-stored rows must remain
     * valid, and sp_finish must emit the "incomplete" diagnostic naming the
     * exact dropped count -- all without relying on host OOM. */
    const char *tsv3 = "/tmp/stage_prof_test3.tsv";
    const char *err3 = "/tmp/stage_prof_test3.err";
    remove(tsv3);
    remove(err3);
    sp_init(tsv3, "t", "v", "x86_64", 4, "sam", -1, "in");
    force_one_drop();   /* fills to capacity, then drops exactly one chunk */

    finish_capturing(err3, 1.0, 0.0, 0.0);

    /* The stored rows are intact: aggregate n_reads sums only the 64 that landed,
     * not the dropped 65th (each chunk carried n_reads == 1). */
    assert(agg_field(tsv3, 10, field, sizeof field) && strcmp(field, "64") == 0);

    /* Aggregate identity is not enough: assert the 64 retained rows are exactly
     * chunk ids 0..63, each present once, with n_reads/n_bp intact -- so a drop
     * path that corrupts, replaces, or duplicates a record (leaving the sum at
     * 64) is still caught. */
    assert_retained_chunks(tsv3, 64);

    /* sp_finish reported the exact drop count and flagged the report incomplete.
     * Asserting the exact "1 chunk(s) dropped" -- not merely the "dropped"
     * keyword -- catches a miscount a substring check would let pass. */
    slurp(err3, errbuf, sizeof errbuf);
    assert(strstr(errbuf, "1 chunk(s) dropped") != NULL);
    assert(strstr(errbuf, "incomplete") != NULL);
    remove(tsv3);
    remove(err3);

    printf("stage_prof_test OK\n");
    return 0;
}
