#include "bwa_thp.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__linux__)
#include <errno.h>
#include <fcntl.h>
#include <time.h>
#include <unistd.h>
#endif

extern int bwa_verbose;   /* declared in bwa.h; not included here to keep this file dependency-free */

#define THP_SYSFS "/sys/kernel/mm/transparent_hugepage/"

/* The anonymous mTHP sizes below the 2 MB PMD size on a 4 KB-page kernel
 * (orders 2-8; Linux >= 6.8). Absent on older kernels. */
static const char *const thp_other_knobs[] = {
    THP_SYSFS "hugepages-16kB/enabled",  THP_SYSFS "hugepages-32kB/enabled",
    THP_SYSFS "hugepages-64kB/enabled",  THP_SYSFS "hugepages-128kB/enabled",
    THP_SYSFS "hugepages-256kB/enabled", THP_SYSFS "hugepages-512kB/enabled",
    THP_SYSFS "hugepages-1024kB/enabled",
};
#define THP_N_OTHER ((int)(sizeof thp_other_knobs / sizeof thp_other_knobs[0]))

/* Value of the `key` line ("Key:   123 kB") at the start of any line of `text`;
 * -1 if absent or malformed. */
static long long thp_field(const char *text, const char *key)
{
    const size_t klen = strlen(key);
    for (const char *p = text; p != NULL && *p != '\0'; ) {
        if (strncmp(p, key, klen) == 0) {
            char *end = NULL;
            long long v = strtoll(p + klen, &end, 10);
            if (end == p + klen || v < 0) return -1;
            return v;
        }
        p = strchr(p, '\n');
        if (p != NULL) ++p;
    }
    return -1;
}

bwa_thp_mode_t bwa_thp_parse_enabled(const char *text)
{
    if (text == NULL) return BWA_THP_MODE_UNKNOWN;
    const char *open = strchr(text, '[');
    if (open == NULL) return BWA_THP_MODE_UNKNOWN;
    const char *close = strchr(open, ']');
    if (close == NULL) return BWA_THP_MODE_UNKNOWN;
    const size_t len = (size_t)(close - open - 1);
    if (len == 6 && strncmp(open + 1, "always", 6) == 0)  return BWA_THP_MODE_ALWAYS;
    if (len == 7 && strncmp(open + 1, "madvise", 7) == 0) return BWA_THP_MODE_MADVISE;
    if (len == 5 && strncmp(open + 1, "never", 5) == 0)   return BWA_THP_MODE_NEVER;
    if (len == 7 && strncmp(open + 1, "inherit", 7) == 0) return BWA_THP_MODE_INHERIT;
    return BWA_THP_MODE_UNKNOWN;
}

bwa_thp_config_t bwa_thp_resolve_config(bwa_thp_mode_t global, bwa_thp_mode_t pmd,
                                        const bwa_thp_mode_t *other, int n_other)
{
    bwa_thp_config_t c = { BWA_THP_MODE_UNKNOWN, NULL };
    if (global != BWA_THP_MODE_ALWAYS && global != BWA_THP_MODE_MADVISE &&
        global != BWA_THP_MODE_NEVER)
        return c;
    if (pmd == BWA_THP_MODE_UNKNOWN) return c;
    for (int i = 0; i < n_other && other != NULL; ++i) {
        const bwa_thp_mode_t o = other[i] == BWA_THP_MODE_INHERIT ? global : other[i];
        if (o == BWA_THP_MODE_UNKNOWN || o == BWA_THP_MODE_ALWAYS || o == BWA_THP_MODE_MADVISE)
            return c;
    }
    if (pmd == BWA_THP_MODE_INHERIT) {
        c.pmd_mode = global;
        c.knob = BWA_THP_KNOB_GLOBAL;
    } else {
        c.pmd_mode = pmd;
        c.knob = BWA_THP_KNOB_PMD;
    }
    return c;
}

void bwa_thp_smaps_init(bwa_thp_smaps_parser_t *p)
{
    memset(p, 0, sizeof *p);
    p->cur_anon = p->cur_huge = -1;
    p->cur_elig = -1;
}

/* Close the current VMA record into the totals. */
static void thp_smaps_commit(bwa_thp_smaps_parser_t *p)
{
    if (!p->in_vma) return;
    const long long a = p->cur_anon > 0 ? p->cur_anon : 0;
    const long long h = p->cur_huge > 0 ? p->cur_huge : 0;
    p->u.vmas++;
    p->u.anon_kb += a;
    p->u.anon_huge_kb += h;
    if (p->cur_elig >= 0) p->u.elig_fields++;
    if (p->cur_elig == 1) {
        p->u.elig_kb += a;
        p->u.elig_huge_kb += h;
    }
    p->cur_anon = p->cur_huge = -1;
    p->cur_elig = -1;
}

static void thp_smaps_line(bwa_thp_smaps_parser_t *p)
{
    const char *s = p->line;
    const char c = s[0];
    /* VMA header: "<start>-<end> perms ..." in lowercase hex; field names
     * start with an upper-case letter. */
    if ((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f')) {
        thp_smaps_commit(p);
        p->in_vma = 1;
        return;
    }
    if (!p->in_vma) return;
    if (strncmp(s, "Anonymous:", 10) == 0)           p->cur_anon = thp_field(s, "Anonymous:");
    else if (strncmp(s, "AnonHugePages:", 14) == 0)  p->cur_huge = thp_field(s, "AnonHugePages:");
    else if (strncmp(s, "THPeligible:", 12) == 0) {
        const long long v = thp_field(s, "THPeligible:");
        p->cur_elig = v == 0 || v == 1 ? (int)v : -1;
    }
}

void bwa_thp_smaps_feed(bwa_thp_smaps_parser_t *p, const char *data, size_t n)
{
    while (n > 0) {
        const char *nl = (const char *)memchr(data, '\n', n);
        const size_t take = nl != NULL ? (size_t)(nl - data) : n;
        const size_t room = BWA_THP_LINE_MAX - 1 - p->len;
        const size_t keep = take < room ? take : room;
        memcpy(p->line + p->len, data, keep);
        p->len += keep;
        if (nl == NULL) return;
        p->line[p->len] = '\0';
        thp_smaps_line(p);
        p->len = 0;
        data = nl + 1;
        n -= take + 1;
    }
}

void bwa_thp_smaps_finish(bwa_thp_smaps_parser_t *p, bwa_thp_usage_t *out)
{
    if (p->len > 0) {
        p->line[p->len] = '\0';
        thp_smaps_line(p);
        p->len = 0;
    }
    thp_smaps_commit(p);
    if (out != NULL) *out = p->u;
}

bwa_thp_status_t bwa_thp_parse_status(const char *text)
{
    bwa_thp_status_t s = { -1, -1 };
    if (text == NULL) return s;
    s.rss_anon_kb = thp_field(text, "RssAnon:");
    const long long e = thp_field(text, "THP_enabled:");
    s.thp_enabled = e == 0 || e == 1 ? (int)e : -1;
    return s;
}

int bwa_thp_sample_complete(long long smaps_anon_kb, long long rss_anon_kb)
{
    if (rss_anon_kb < 0) return 1;
    return smaps_anon_kb * 10 >= rss_anon_kb * 9;
}

long long bwa_thp_short_kb(long long anon_kb, long long huge_kb)
{
    return huge_kb < anon_kb ? anon_kb - huge_kb : 0;
}

/* total_anon: the process size gate; base and short_kb: the population judged. */
static int thp_material(long long total_anon, long long base, long long short_kb)
{
    return total_anon >= BWA_THP_MIN_ANON_KB && short_kb >= BWA_THP_MIN_SHORT_KB &&
           short_kb * 100 >= (long long)BWA_THP_MIN_SHORT_PCT * base;
}

bwa_thp_verdict_t bwa_thp_assess(bwa_thp_mode_t pmd_mode, int thp_enabled,
                                 const bwa_thp_usage_t *u)
{
    if (u == NULL) return BWA_THP_OK;
    if (pmd_mode != BWA_THP_MODE_ALWAYS && pmd_mode != BWA_THP_MODE_MADVISE &&
        pmd_mode != BWA_THP_MODE_NEVER)
        return BWA_THP_OK;
    const long long all_short = bwa_thp_short_kb(u->anon_kb, u->anon_huge_kb);
    if (pmd_mode == BWA_THP_MODE_NEVER)
        return thp_material(u->anon_kb, u->anon_kb, all_short) ? BWA_THP_WARN_DISABLED : BWA_THP_OK;
    if (thp_enabled == 0)
        return thp_material(u->anon_kb, u->anon_kb, all_short) ? BWA_THP_WARN_PROCESS_DISABLED
                                                               : BWA_THP_OK;
    if (u->elig_fields == 0) return BWA_THP_OK;   /* no THPeligible (Linux < 5.0) */
    const long long elig_short = bwa_thp_short_kb(u->elig_kb, u->elig_huge_kb);
    return thp_material(u->anon_kb, u->elig_kb, elig_short) ? BWA_THP_WARN_FALLBACK : BWA_THP_OK;
}

int bwa_thp_should_sample(long batch)
{
    return batch > 0 && (batch & (batch - 1)) == 0;
}

int bwa_thp_format_warning(bwa_thp_verdict_t v, const bwa_thp_usage_t *u,
                           const char *knob, const char *func, char *buf, size_t n)
{
    if (buf == NULL || n == 0) return 0;
    buf[0] = '\0';
    if (v == BWA_THP_OK || u == NULL) return 0;
    const double kb_per_gib = 1024.0 * 1024.0;
    const char *f = func != NULL ? func : "mem";
    static const char docs[] =
        "see docs/src/best-practices/optimization-checklist.md, \"Keep transparent huge pages working\".";
    static const char cost[] = "on human WGS this costs on the order of 1% wall time per GiB";
    if (v == BWA_THP_WARN_FALLBACK) {
        const double base_g  = (double)u->elig_kb / kb_per_gib;
        const double short_g = (double)bwa_thp_short_kb(u->elig_kb, u->elig_huge_kb) / kb_per_gib;
        return snprintf(buf, n,
            "[W::%s] %.1f GiB of %.1f GiB of THP-eligible anonymous memory is on 4 KB pages, "
            "not transparent huge pages (typically because host memory was full of page "
            "cache or fragmented); %s.\n"
            "[W::%s] start with free memory (e.g. drop the page cache) or run in a "
            "memory-limited cgroup; %s\n",
            f, short_g, base_g, cost, f, docs);
    }
    const double anon_g  = (double)u->anon_kb / kb_per_gib;
    const double short_g = (double)bwa_thp_short_kb(u->anon_kb, u->anon_huge_kb) / kb_per_gib;
    if (v == BWA_THP_WARN_DISABLED) {
        return snprintf(buf, n,
            "[W::%s] transparent huge pages are disabled (%s is [never]): %.1f GiB of "
            "%.1f GiB of anonymous memory is on 4 KB pages; %s (about 28%% on 5M-pair human WGS).\n"
            "[W::%s] set it to madvise or always; %s\n",
            f, knob != NULL ? knob : BWA_THP_KNOB_GLOBAL, short_g, anon_g, cost, f, docs);
    }
    return snprintf(buf, n,
        "[W::%s] transparent huge pages are disabled for this process (/proc/self/status "
        "THP_enabled: 0, set by prctl(PR_SET_THP_DISABLE), e.g. MIMALLOC_ALLOW_THP=0 or "
        "inherited from the launcher): %.1f GiB of %.1f GiB of anonymous memory is on "
        "4 KB pages; %s.\n"
        "[W::%s] unset MIMALLOC_ALLOW_THP=0 or launch without PR_SET_THP_DISABLE; %s\n",
        f, short_g, anon_g, cost, f, docs);
}

const char *bwa_thp_stop_reason(bwa_thp_stop_t s)
{
    switch (s) {
    case BWA_THP_STOP_PAGE_SIZE: return "base page is not 4 KB";
    case BWA_THP_STOP_CONFIG:    return "THP mode unreadable, unrecognised, or a multi-size THP size other than 2 MB enabled";
    case BWA_THP_STOP_SMAPS:     return "/proc/self/smaps unreadable";
    default:                     return "";
    }
}

bwa_thp_verdict_t bwa_thp_step(bwa_thp_monitor_t *m, int verbose,
                               const bwa_thp_io_t *io, bwa_thp_report_t *rep)
{
    bwa_thp_report_t r;
    memset(&r, 0, sizeof r);
    r.status.rss_anon_kb = -1;
    r.status.thp_enabled = -1;
    bwa_thp_verdict_t v = BWA_THP_OK;
    if (m == NULL || io == NULL || m->done) goto out;
    if (verbose < 2) { m->done = 1; goto out; }
    ++m->batches;
    {
        const int scheduled = bwa_thp_should_sample(m->batches);
        if (!scheduled && !m->retry_pending) goto out;
        m->retry_pending = 0;
        r.sampled = 1;
        r.retry = !scheduled;

        if (m->cfg.pmd_mode == BWA_THP_MODE_UNKNOWN) {
            /* The messages speak of 4 KB pages and 2 MB THP; other base page
             * sizes (e.g. a 64 KB arm64 kernel, 512 MB PMD THP) are out of scope. */
            if (io->page_size(io->ctx) != 4096) {
                m->done = 1;
                r.stop = BWA_THP_STOP_PAGE_SIZE;
                goto out;
            }
            const bwa_thp_mode_t g = bwa_thp_parse_enabled(io->read_text(io->ctx, BWA_THP_KNOB_GLOBAL));
            const char *t = io->read_text(io->ctx, BWA_THP_KNOB_PMD);
            const bwa_thp_mode_t pm = t != NULL ? bwa_thp_parse_enabled(t) : BWA_THP_MODE_INHERIT;
            bwa_thp_mode_t om[THP_N_OTHER];
            for (int i = 0; i < THP_N_OTHER; ++i) {
                t = io->read_text(io->ctx, thp_other_knobs[i]);
                om[i] = t != NULL ? bwa_thp_parse_enabled(t) : BWA_THP_MODE_NEVER;
            }
            m->cfg = bwa_thp_resolve_config(g, pm, om, THP_N_OTHER);
            if (m->cfg.pmd_mode == BWA_THP_MODE_UNKNOWN) {
                m->done = 1;
                r.stop = BWA_THP_STOP_CONFIG;
                goto out;
            }
        }

        while (r.walks < BWA_THP_MAX_WALKS && !r.complete) {
            ++r.walks;
            if (io->walk_smaps(io->ctx, &r.usage) != 0) {
                m->done = 1;
                r.stop = BWA_THP_STOP_SMAPS;
                goto out;
            }
            r.status = bwa_thp_parse_status(io->read_text(io->ctx, "/proc/self/status"));
            r.complete = bwa_thp_sample_complete(r.usage.anon_kb, r.status.rss_anon_kb);
        }
        if (!r.complete) {
            /* At most one deferred retry per scheduled sample. */
            ++m->incomplete;
            m->retry_pending = scheduled;
            goto out;
        }
        v = bwa_thp_assess(m->cfg.pmd_mode, r.status.thp_enabled, &r.usage);
        if (v != BWA_THP_OK) m->done = 1;
    }
out:
    if (rep != NULL) *rep = r;
    return v;
}

#if defined(__linux__)
/* A separate reader from system.cpp's read_small_file(): that one is an
 * anonymous-namespace C++ helper returning std::string and refuses files over
 * 4 KiB, while /proc/self/status can exceed that (a long Groups: line) and
 * smaps is megabytes, so this module reads to EOF with a growing buffer
 * (status, sysfs) or streams in chunks (smaps). */
typedef struct {
    char  *buf;
    size_t cap;
} thp_real_ctx_t;

#define THP_READ_CAP (1u << 20)   /* status and sysfs knobs are far below this */

static long thp_real_page_size(void *ctx)
{
    (void)ctx;
    return sysconf(_SC_PAGESIZE);
}

static const char *thp_real_read_text(void *vctx, const char *path)
{
    thp_real_ctx_t *c = (thp_real_ctx_t *)vctx;
    const int fd = open(path, O_RDONLY | O_CLOEXEC);
    if (fd < 0) return NULL;
    size_t len = 0;
    for (;;) {
        if (len + 1 >= c->cap) {
            const size_t ncap = c->cap ? c->cap * 2 : 4096;
            char *nb = ncap <= THP_READ_CAP ? (char *)realloc(c->buf, ncap) : NULL;
            if (nb == NULL) { close(fd); return NULL; }
            c->buf = nb;
            c->cap = ncap;
        }
        const ssize_t r = read(fd, c->buf + len, c->cap - 1 - len);
        if (r < 0) {
            if (errno == EINTR) continue;
            close(fd);
            return NULL;
        }
        if (r == 0) break;
        len += (size_t)r;
    }
    close(fd);
    c->buf[len] = '\0';
    return c->buf;
}

static int thp_real_walk_smaps(void *vctx, bwa_thp_usage_t *out)
{
    (void)vctx;
    enum { CHUNK = 64 * 1024 };
    char *chunk = (char *)malloc(CHUNK);
    if (chunk == NULL) return -1;
    const int fd = open("/proc/self/smaps", O_RDONLY | O_CLOEXEC);
    if (fd < 0) { free(chunk); return -1; }
    bwa_thp_smaps_parser_t p;
    bwa_thp_smaps_init(&p);
    int rc = 0;
    for (;;) {
        const ssize_t r = read(fd, chunk, CHUNK);
        if (r < 0) {
            if (errno == EINTR) continue;
            rc = -1;
            break;
        }
        if (r == 0) break;
        bwa_thp_smaps_feed(&p, chunk, (size_t)r);
    }
    close(fd);
    free(chunk);
    if (rc == 0) bwa_thp_smaps_finish(&p, out);
    return rc;
}
#endif

void bwa_thp_check_batch(bwa_thp_monitor_t *m)
{
#if defined(__linux__)
    if (m == NULL || m->done) return;
    thp_real_ctx_t ctx = { NULL, 0 };
    const bwa_thp_io_t io = { &ctx, thp_real_page_size, thp_real_read_text, thp_real_walk_smaps };
    const int dbg = bwa_verbose >= 4;
    struct timespec t0 = { 0, 0 }, t1 = { 0, 0 };
    if (dbg) clock_gettime(CLOCK_MONOTONIC, &t0);
    bwa_thp_report_t rep;
    const bwa_thp_verdict_t v = bwa_thp_step(m, bwa_verbose, &io, &rep);
    if (dbg) clock_gettime(CLOCK_MONOTONIC, &t1);
    free(ctx.buf);
    /* Named like the other `mem` warnings ([W::main_mem]); this runs from its pipeline. */
    static const char func[] = "main_mem";
    if (dbg && rep.sampled && rep.stop != BWA_THP_STOP_NONE)
        fprintf(stderr, "[D::%s] THP check batch %ld: check disabled: %s\n", func, m->batches,
                bwa_thp_stop_reason(rep.stop));
    else if (dbg && rep.sampled)
        fprintf(stderr, "[D::%s] THP check batch %ld%s: Anonymous %lld kB (AnonHugePages %lld kB), "
                        "THP-eligible %lld kB (AnonHugePages %lld kB), %ld VMAs, RssAnon %lld kB, "
                        "THP_enabled %d, %d walk(s)%s, %ld incomplete sample(s) so far, %.3f ms\n",
                func, m->batches, rep.retry ? " (retry)" : "",
                rep.usage.anon_kb, rep.usage.anon_huge_kb, rep.usage.elig_kb,
                rep.usage.elig_huge_kb, rep.usage.vmas, rep.status.rss_anon_kb,
                rep.status.thp_enabled, rep.walks, rep.complete ? "" : " (incomplete; skipped)",
                m->incomplete,
                (double)(t1.tv_sec - t0.tv_sec) * 1e3 + (double)(t1.tv_nsec - t0.tv_nsec) / 1e6);
    if (v != BWA_THP_OK) {
        char msg[1024];
        bwa_thp_format_warning(v, &rep.usage, m->cfg.knob, func, msg, sizeof msg);
        fputs(msg, stderr);
    }
#else
    (void)m;
#endif
}
