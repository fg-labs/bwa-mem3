#ifndef BWA_THP_H
#define BWA_THP_H

#include <stddef.h>

/* Transparent-huge-page (THP) fallback check for `mem` (Linux).
 *
 * The index arrays are madvised for THP (2 MB pages), and the allocator's
 * arenas are too where it does so, but the kernel only backs them with huge
 * pages when it can find free, unfragmented 2 MB blocks at fault time. On a
 * host whose RAM is full of page cache it often cannot, and part of the hot
 * working set silently lands on 4 KB pages, at a measurable wall-time cost.
 * The fix is environmental, so `mem` only says so on stderr. The costs, walk
 * timings and the user-facing description are in
 * docs/src/best-practices/optimization-checklist.md, "Keep transparent huge
 * pages working".
 *
 * What counts as a fallback: only memory the kernel considered eligible for
 * PMD (2 MB) THP. The check walks /proc/self/smaps and sums `Anonymous:` and
 * `AnonHugePages:` over the VMAs whose `THPeligible:` is 1. On the 6.1 kernel
 * (AL2023) that bit is hugepage_vma_check(): the VMA is anonymous (or THP-able
 * file), not MADV_NOHUGEPAGE, its process has no PR_SET_THP_DISABLE, it can
 * hold at least one aligned 2 MB range, and the sysfs mode allows it (`always`,
 * or `madvise` and the VMA was madvised). So eligible-but-not-huge memory is
 * memory the kernel would have backed with 2 MB pages and did not: a fallback.
 * Unadvised memory in `madvise` mode (e.g. a USE_MIMALLOC=0 build's arenas) is
 * not eligible and never counts; under `always` every anonymous VMA that can
 * hold an aligned 2 MB range is eligible, advised or not. The unaligned head/tail of an eligible VMA
 * (< 2 MB each) also counts as not huge; the 1 GiB floor absorbs that.
 * Kernel caveats: THPeligible exists since Linux 5.0 (absent: the fallback
 * verdict is skipped); from 6.8 it also reports eligibility for smaller
 * multi-size THP (mTHP) orders, so any mTHP size other than 2 MB being enabled
 * puts the run out of scope (see bwa_thp_resolve_config).
 *
 * Three warnings, at most one per run:
 *   DISABLED          effective 2 MB THP mode is `never` (sysfs).
 *   PROCESS_DISABLED  /proc/self/status `THP_enabled: 0` (Linux >= 5.0):
 *                     prctl(PR_SET_THP_DISABLE), e.g. MIMALLOC_ALLOW_THP=0, or
 *                     inherited from the launcher across fork/exec.
 *   FALLBACK          eligible memory materially short of huge pages.
 *
 * Everything but the /proc and /sys reads is a pure function (the reads come in
 * through bwa_thp_io_t), so the decisions and the per-batch state machine are
 * unit-tested with synthetic text. Nothing here touches alignment output. */

#ifdef __cplusplus
extern "C" {
#endif

/* A bracketed mode in /sys/kernel/mm/transparent_hugepage/enabled or in a
 * per-size .../hugepages-<N>kB/enabled knob (Linux >= 6.8). */
typedef enum {
    BWA_THP_MODE_UNKNOWN = 0, /* unreadable, no bracketed token, or out of scope */
    BWA_THP_MODE_ALWAYS,
    BWA_THP_MODE_MADVISE,
    BWA_THP_MODE_NEVER,
    BWA_THP_MODE_INHERIT      /* per-size knobs only: follow the global mode */
} bwa_thp_mode_t;

/* Parse the bracketed mode, e.g. "always [madvise] never\n" -> MADVISE. */
bwa_thp_mode_t bwa_thp_parse_enabled(const char *text);

/* The sysfs knobs that set the 2 MB mode. */
#define BWA_THP_KNOB_GLOBAL "/sys/kernel/mm/transparent_hugepage/enabled"
#define BWA_THP_KNOB_PMD    "/sys/kernel/mm/transparent_hugepage/hugepages-2048kB/enabled"

/* The effective 2 MB (PMD) THP configuration. */
typedef struct {
    bwa_thp_mode_t pmd_mode; /* ALWAYS, MADVISE or NEVER; UNKNOWN = do not check */
    const char *knob;        /* the knob that set pmd_mode (BWA_THP_KNOB_*); NULL if UNKNOWN */
} bwa_thp_config_t;

/* Resolve the effective PMD mode from the parsed knobs: the global `enabled`
 * mode, the hugepages-2048kB/enabled mode (pass INHERIT when the knob is absent,
 * i.e. a kernel without per-size knobs; INHERIT follows the global mode), and
 * the other per-size knobs' modes (pass NEVER for an absent knob). If any other
 * size resolves to always/madvise, memory in those folios is neither on 4 KB
 * pages nor in AnonHugePages and THPeligible covers those orders too, so the
 * run is out of scope: UNKNOWN. Any UNKNOWN input is UNKNOWN as well. */
bwa_thp_config_t bwa_thp_resolve_config(bwa_thp_mode_t global, bwa_thp_mode_t pmd,
                                        const bwa_thp_mode_t *other, int n_other);

/* Per-walk totals over /proc/self/smaps, in kB. Shmem (a `bwa-mem3 shm` index)
 * and HugeTLB (a `--huge-pages` reservation) are not `Anonymous`. */
typedef struct {
    long long anon_kb;      /* Anonymous: over all VMAs */
    long long anon_huge_kb; /* AnonHugePages: (PMD-mapped THP) over all VMAs */
    long long elig_kb;      /* Anonymous: over VMAs with THPeligible: 1 */
    long long elig_huge_kb; /* AnonHugePages: over VMAs with THPeligible: 1 */
    long vmas;              /* VMAs seen */
    long elig_fields;       /* VMAs that carried a THPeligible: line (0 on < 5.0) */
} bwa_thp_usage_t;

/* Streaming /proc/self/smaps parser: smaps is one record per VMA and can be
 * megabytes, so it is fed in chunks of any size (lines may span chunks). Only
 * a line's first BWA_THP_LINE_MAX-1 bytes are kept, which holds every field
 * the parser reads; a VMA header line is recognised by its leading hex digit. */
#define BWA_THP_LINE_MAX 64
typedef struct {
    bwa_thp_usage_t u;
    long long cur_anon, cur_huge;
    int  cur_elig, in_vma;
    size_t len;
    char line[BWA_THP_LINE_MAX];
} bwa_thp_smaps_parser_t;

void bwa_thp_smaps_init(bwa_thp_smaps_parser_t *p);
void bwa_thp_smaps_feed(bwa_thp_smaps_parser_t *p, const char *data, size_t n);
/* Flush the last line and VMA; copy the totals to *out. */
void bwa_thp_smaps_finish(bwa_thp_smaps_parser_t *p, bwa_thp_usage_t *out);

/* The two /proc/self/status fields the check uses; -1 when absent. */
typedef struct {
    long long rss_anon_kb; /* RssAnon: */
    int thp_enabled;       /* THP_enabled: (Linux >= 5.0), 0 or 1 */
} bwa_thp_status_t;
bwa_thp_status_t bwa_thp_parse_status(const char *text);

/* An smaps walk is not a snapshot: the kernel releases mmap_lock between
 * read() calls, and allocator mmap/munmap from the worker threads can make a
 * walk undercount (one smaps_rollup walk reported 10 of 24 GiB). RssAnon is an
 * O(1) counter, so a walk is accepted as complete when its total Anonymous is
 * >= 90% of RssAnon read right after it. An unknown RssAnon (< 0) accepts. */
int bwa_thp_sample_complete(long long smaps_anon_kb, long long rss_anon_kb);

typedef enum {
    BWA_THP_OK = 0,              /* nothing to report */
    BWA_THP_WARN_DISABLED,       /* 2 MB THP mode is `never` and memory is large */
    BWA_THP_WARN_PROCESS_DISABLED, /* THP_enabled: 0 for this process, memory large */
    BWA_THP_WARN_FALLBACK        /* eligible memory materially short of huge pages */
} bwa_thp_verdict_t;

/* Thresholds (kB). A verdict is material only when all three hold:
 *   total Anonymous >= BWA_THP_MIN_ANON_KB              (a large-genome run)
 *   shortfall       >= BWA_THP_MIN_SHORT_KB             (>= 1 GiB not huge)
 *   shortfall * 100 >= BWA_THP_MIN_SHORT_PCT * base
 * where shortfall/base are over the eligible VMAs for FALLBACK and over all
 * anonymous memory for the two DISABLED verdicts. A fully-THP human WGS run
 * measured ~0.02 GiB short of ~28 GiB; runs on page-cache-full hosts were
 * 1.2-10 GiB short. */
#define BWA_THP_MIN_ANON_KB   (4LL * 1024 * 1024)
#define BWA_THP_MIN_SHORT_KB  (1LL * 1024 * 1024)
#define BWA_THP_MIN_SHORT_PCT 4

/* anon_kb - huge_kb, clamped at 0 (AnonHugePages can momentarily exceed
 * Anonymous across lines of one read). Used by assess and format alike. */
long long bwa_thp_short_kb(long long anon_kb, long long huge_kb);

/* The decision. UNKNOWN mode never warns. NEVER takes precedence over a
 * process-level disable (THP_enabled 0; -1 = unknown, ignored), which takes
 * precedence over FALLBACK. FALLBACK needs THPeligible lines (Linux >= 5.0). */
bwa_thp_verdict_t bwa_thp_assess(bwa_thp_mode_t pmd_mode, int thp_enabled,
                                 const bwa_thp_usage_t *u);

/* Sampling schedule over 1-based completed-batch numbers: 1, 2, 4, 8, ...
 * Bounded (log2 of the batch count) and front-loaded: the index and most of
 * the per-batch buffers are touched by the end of the first batch. A walk
 * visits every PTE, so on 4 KB pages it costs tens of ms or more (timings in
 * the docs section named at the top); it must not run per batch. */
int bwa_thp_should_sample(long batch);

/* Format the two warning lines (each ending in '\n', each prefixed
 * "[W::<func>] ", func NULL -> "mem") for a non-OK verdict into buf. `knob` is
 * the sysfs path named by DISABLED (NULL -> BWA_THP_KNOB_GLOBAL). Returns the
 * snprintf length (output truncated if >= n); returns 0 and writes "" for
 * BWA_THP_OK or a NULL u. */
int bwa_thp_format_warning(bwa_thp_verdict_t v, const bwa_thp_usage_t *u,
                           const char *knob, const char *func, char *buf, size_t n);

/* The reads the per-batch step needs, injectable for tests. */
typedef struct {
    void *ctx;
    long (*page_size)(void *ctx);
    /* Whole file as NUL-terminated text, valid until the next call; NULL if
     * absent or unreadable. Used for the sysfs knobs and /proc/self/status. */
    const char *(*read_text)(void *ctx, const char *path);
    /* One full walk of /proc/self/smaps into *out; 0 on success, -1 unreadable. */
    int (*walk_smaps)(void *ctx, bwa_thp_usage_t *out);
} bwa_thp_io_t;

/* Why the check stopped for good before reaching a verdict (for -v 4). */
typedef enum {
    BWA_THP_STOP_NONE = 0,
    BWA_THP_STOP_PAGE_SIZE,  /* base page is not 4 KB */
    BWA_THP_STOP_CONFIG,     /* THP mode unreadable, unrecognised, or mTHP out of scope */
    BWA_THP_STOP_SMAPS       /* /proc/self/smaps unreadable */
} bwa_thp_stop_t;

/* Short text for a stop reason, e.g. "base page is not 4 KB"; "" for NONE. */
const char *bwa_thp_stop_reason(bwa_thp_stop_t s);

/* Per-run state, zero-initialised. Touched only from the write step of the
 * `mem` pipeline, which runs one batch at a time in order. */
typedef struct {
    long batches;          /* completed batches seen */
    int  done;             /* warned, or checking impossible: stop */
    int  retry_pending;    /* the last scheduled sample was incomplete: one retry next batch */
    long incomplete;       /* samples skipped as incomplete */
    bwa_thp_config_t cfg;  /* resolved on the first sample; pmd_mode UNKNOWN until then */
} bwa_thp_monitor_t;

/* What one step did, for the -v 4 debug line. */
typedef struct {
    int sampled;              /* a sample was taken this batch */
    int walks;                /* smaps walks made (<= 3) */
    int complete;             /* the last walk passed bwa_thp_sample_complete */
    int retry;                /* this sample was the deferred retry */
    bwa_thp_stop_t stop;      /* set when this step stopped the check for good */
    bwa_thp_usage_t usage;    /* the last walk */
    bwa_thp_status_t status;  /* status read after the last walk */
} bwa_thp_report_t;

#define BWA_THP_MAX_WALKS 3

/* Two retry mechanisms, because an smaps walk races two kinds of churn: the
 * back-to-back re-walks (BWA_THP_MAX_WALKS) ride out a brief mmap/munmap
 * burst that overlaps one walk, while the deferred retry covers a sample that
 * lands inside a longer allocation phase, which would otherwise wait for the
 * next power-of-two batch (e.g. 64 -> 128) and so may never be judged.
 *
 * One completed batch. Pure apart from `io`: with verbose < 2 it stops; it
 * samples on the bwa_thp_should_sample() schedule, plus one deferred retry on
 * the next batch when a scheduled sample had no complete walk in
 * BWA_THP_MAX_WALKS; the first sample resolves the config and stops when the
 * base page is not 4 KB or the config is UNKNOWN; an unreadable smaps stops.
 * Each of those stops sets rep->stop. Returns the verdict to print (at most
 * once: a non-OK verdict sets done). */
bwa_thp_verdict_t bwa_thp_step(bwa_thp_monitor_t *m, int verbose,
                               const bwa_thp_io_t *io, bwa_thp_report_t *rep);

/* Call once per completed output batch: bwa_thp_step() over the real /proc and
 * /sys, printing the warning and (at -v 4) a [D::] line per sample. No-op
 * off Linux. */
void bwa_thp_check_batch(bwa_thp_monitor_t *m);

#ifdef __cplusplus
}
#endif

#endif /* BWA_THP_H */
