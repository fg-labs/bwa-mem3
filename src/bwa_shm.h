/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
         Heng Li <hli@jimmy.harvard.edu>
*****************************************************************************************/

#ifndef BWA_SHM_H_
#define BWA_SHM_H_

#include <stdint.h>
#include <stddef.h>
#include <limits.h>   /* PATH_MAX */

/* "BWAMEM2\0" interpreted as a little-endian uint64_t:
 *   'B'=0x42 'W'=0x57 'A'=0x41 'M'=0x4D 'E'=0x45 'M'=0x4D '2'=0x32 '\0'=0x00
 * little-endian => bytes 0x42 0x57 0x41 0x4D 0x45 0x4D 0x32 0x00 read as u64. */
#define BWA_SHM_MAGIC          0x00324D454D415742ull
/* Bumped 1 -> 2 when FMI_SCALARS grew from 56 to 64 bytes (sa_compx appended).
 * A stale v1 segment would otherwise pass bwa_shm_attach() and then hard-fail
 * the 64-byte FMI_SCALARS size check in load_index_from_shm() with no disk
 * fallback; rejecting it at attach lets load_index() fall back to disk.
 * Bumped 2 -> 3 when FMI_SA_MS changed from one byte per SA sample to the packed
 * high parts (FMI_search::sa_hi_size_bytes()).
 * Bumped 3 -> 4 when the header gained cp_layout (the compact checkpoint
 * layout, `shm --compact-occ`): a version-3 reader would take a compact
 * segment's header for a default one, and a version-4 reader cannot tell a
 * version-3 segment's layout, so each rejects the other at attach and loads
 * from disk.
 * Bump it whenever the header or any packed section layout (e.g.
 * BWA_SHM_FMI_SCALARS_BYTES) changes. */
#define BWA_SHM_VERSION        4u
/* The occurrence-checkpoint layout a segment holds (bwa_shm_header_t.cp_layout).
 * DEFAULT: the file's 64-byte records in FMI_CP_OCC. COMPACT (`shm
 * --compact-occ`): the 32-byte records in FMI_CP_PLANES instead, plus
 * FMI_CP_BASE and FMI_CP_LAYOUT. Any other value is refused at attach. */
#define BWA_SHM_CP_LAYOUT_DEFAULT 0u
#define BWA_SHM_CP_LAYOUT_COMPACT 1u
#define BWA_SHM_CTL_NAME       "/bwactl"
/* Named POSIX semaphore guarding read-modify-write sequences against
 * /bwactl. A stager holds it only for the few milliseconds of registry and
 * segment bookkeeping, never while packing the index. Created on first use and
 * unlinked by bwa_shm_destroy alongside the registry, so `bwa-mem3 shm -d`
 * recovers from a stager killed inside that window (POSIX semaphores have no
 * SEM_UNDO). Public so a test that cleans up its own staged index takes the
 * same lock. */
#define BWA_SHM_LOCK_NAME      "/bwactl_lock"
#define BWA_SHM_IDX_PREFIX     "/bwaidx-"
#define BWA_SHM_CTL_SIZE       0x10000   /* 64 KiB registry segment */

/* FMI_SCALARS section: int64_t reference_seq_len, count[5], sentinel_index,
 * sa_compx. sa_compx (the SA sample-rate shift) is carried so a non-default
 * `-u`-built index attached via shm is sized correctly rather than assuming
 * the compile-time default. */
#define BWA_SHM_FMI_SCALARS_BYTES (sizeof(int64_t) * 8)

/* Section kinds — see implementation plan for what each holds. */
#define BWA_SHM_SEC_FMI_SCALARS  1u
#define BWA_SHM_SEC_FMI_CP_OCC   2u
#define BWA_SHM_SEC_FMI_SA_MS    3u
#define BWA_SHM_SEC_FMI_SA_LS    4u
#define BWA_SHM_SEC_BNS_STRUCT   5u
#define BWA_SHM_SEC_BNS_AMBS     6u
#define BWA_SHM_SEC_BNS_ANNS     7u
#define BWA_SHM_SEC_BNS_NAMES    8u
#define BWA_SHM_SEC_PAC          9u
#define BWA_SHM_SEC_REF_STRING   10u
/* The compact layout's sections (BWA_SHM_CP_LAYOUT_COMPACT). FMI_CP_PLANES:
 * CP_OCC_PLANES[] (32 bytes per 64 rows). FMI_CP_BASE: int64_t[4] per
 * 2^span_shift rows, the full counts the span-relative ones are taken from;
 * empty unless the counts are wide. FMI_CP_LAYOUT: BWA_SHM_FMI_CP_LAYOUT_BYTES,
 * int64_t span_shift and wide. */
#define BWA_SHM_SEC_FMI_CP_PLANES 11u
#define BWA_SHM_SEC_FMI_CP_BASE   12u
#define BWA_SHM_SEC_FMI_CP_LAYOUT 13u
#define BWA_SHM_FMI_CP_LAYOUT_BYTES (sizeof(int64_t) * 2)

/* 1 GB staging (`bwa-mem3 shm --huge-pages`). /dev/shm backs a segment with
 * 2 MB pages at most, so the index goes into a System V segment created with
 * SHM_HUGETLB | SHM_HUGE_1GB, and the POSIX segment `/bwaidx-<name>` holds only
 * a stub: one bwa_shm_header_t whose reserved[2] is the System V segment's
 * size in bytes (non-zero marks a stub) and whose `magic`, published last, is
 * BWA_SHM_HUGE_MAGIC ("BWAHUGE\0", little-endian). The registry entry, the
 * staging PID and token, and so the waiting, reclaim and ownership rules, are
 * those of any other stage. The System V key is derived from the token.
 *
 * The stub is registered before the System V segment is created and removed
 * after it, so every segment a stage creates is named by a registered stub,
 * even when the stage is killed, and `shm -d` or a restage frees its pages.
 * Only the segment that the stub's own stager created is ever removed: the
 * one whose id the stub records or, for a stager killed before it recorded the
 * id, the one whose kernel-recorded creator is the stub's PID and whose size
 * is the stub's. System V keys
 * are per IPC namespace, so the stub also records the stager's: seen from
 * another namespace (a container that shares /dev/shm but not IPC) the stub is
 * unusable and is left alone, never taken for one whose segment is gone. The
 * kernel reports a creator PID as seen from the asking process's PID
 * namespace, so the stub records the stager's PID namespace too, and the
 * creator of an unrecorded segment is compared only from that namespace. A
 * stub's total_size is that of the index in its segment (informational). */
#define BWA_SHM_HUGE_MAGIC     0x0045475548415742ull

/* A segment is complete once `magic` is BWA_SHM_MAGIC (a 1 GB stub: once it is
 * BWA_SHM_HUGE_MAGIC, published after its segment's). The stager writes the
 * header last and publishes `magic` with a release store, so a reader that
 * loads it with acquire ordering and sees BWA_SHM_MAGIC sees the whole index.
 * While staging, `magic` is 0, reserved[0] holds the stager's PID and
 * reserved[1] a random per-stage token; the stager keeps both in the final
 * header. A segment with magic 0 whose stager is no longer alive was left by
 * an interrupted stage. */
typedef struct {
	uint64_t magic;          /* BWA_SHM_MAGIC (a 1 GB stub: BWA_SHM_HUGE_MAGIC),
	                          * or 0 while staging */
	uint32_t version;        /* BWA_SHM_VERSION */
	uint16_t n_sections;
	uint16_t cp_layout;      /* BWA_SHM_CP_LAYOUT_*; 0 in a 1 GB stub */
	uint64_t total_size;     /* bytes including this header */
	uint64_t reserved[5];    /* [0] stager PID, [1] stage token; 1 GB stub only:
	                          * [2] System V segment bytes, [3] its id + 1,
	                          * [4] the stager's IPC namespace in the low
	                          *     32 bits, its PID namespace in the high
	                          *     32 (0: unknown);
	                          * rest zero */
} bwa_shm_header_t;           /* 64 bytes; first object in every segment */
#ifdef __cplusplus
static_assert(sizeof(bwa_shm_header_t) == 64, "bwa_shm_header_t is 64 bytes");
#endif

typedef struct {
	uint32_t kind;           /* one of BWA_SHM_SEC_* */
	uint32_t flags;          /* zero-fill for now */
	uint64_t offset;         /* bytes from segment start to payload */
	uint64_t size;           /* payload bytes (does not count alignment padding) */
	uint64_t reserved;       /* zero-fill for forward-compat */
} bwa_shm_section_t;          /* 32 bytes */

#ifdef __cplusplus
extern "C" {
#endif

	/* 1 if `prefix`'s basename is registered (which includes a stage still in
	 * progress or one that was interrupted), 0 if not, -1 on error. */
	int       bwa_shm_test(const char *prefix);
	/* 1 if `prefix`'s basename is registered and its segment is complete (the
	 * check bwa_shm_stage makes to find an index already staged), 0 if not, -1
	 * on error. A registered stage still in progress or interrupted, or one
	 * with no segment or an unusable one, gives 0. Complete does not mean a
	 * given load can attach it: bwa_shm_attach also refuses, e.g., a seed-only
	 * stage when the load needs the packed reference. */
	int       bwa_shm_test_complete(const char *prefix);
	/* bns_only=true stages a SEED-only segment (D3 --meth): the FM-index + BNS
	 * are packed, but PAC and REF_STRING (.0123) are omitted (zero-length
	 * sections) because `mem --meth` never reads the seed pac/.0123. The seed
	 * `.0123` need not even exist on disk. Saves ~14.5 GB of shm on hg38. */
	/* target_sa_compx: if in [0, disk_sa_compx), the staged SA sample table is
	 * densified to this shift on the fly (see bwa-mem3 shm -u); -1 (default)
	 * stages the disk sample rate unchanged. Densifying trades a one-time
	 * staging cost + extra shm for faster per-query SA resolution, without
	 * rebuilding the on-disk index.
	 * n_threads: worker count for the densify pass (the per-sample LF-walks are
	 * independent) and the compact conversion; <=0 or 1 runs them serially.
	 * compact: stage the 32-byte occurrence-checkpoint layout
	 * (BWA_SHM_CP_LAYOUT_COMPACT). Not combinable with a target_sa_compx (any
	 * value >= 0), nor with huge.
	 * huge: stage on 1 GB huge pages (`shm --huge-pages`, Linux only; see
	 * BWA_SHM_HUGE_MAGIC). Fails, staging nothing, when they cannot be had; on
	 * an AMD x86 host, where 1 GB staging measured no gain (unless
	 * BWAMEM3_SHM_HUGE_PAGES_AMD=1 is set); where kernel.shm_rmid_forced
	 * would remove the segment when the stager exits, or cannot be read; and
	 * when the index is
	 * already staged on other pages (it is left as it is). */
	int       bwa_shm_stage(const char *prefix, bool bns_only = false,
	                        int target_sa_compx = -1,
	                        int n_threads = 1,
	                        bool huge = false,
	                        bool compact = false); /* loads from disk, packs, stages */
	/* Test-only hook (see bwa_shm.cpp): bwa_shm_stage calls it with its
	 * writable index mapping after it has registered the segment and released
	 * the stager lock, before the pack fills it, so a unit test can inspect
	 * that mapping or act while the stage is in progress. NULL on every
	 * production path. */
	extern void (*bwa_shm_stage_map_hook)(const void *base, size_t len);
	/* Drops all (matches v1 -d). -1 if a segment could not be removed (a 1 GB
	 * segment, a header that could not be read, or, when the stager lock could
	 * not be taken, a 1 GB stage in progress); that index's segment, its entry
	 * and the registry are then kept. */
	int       bwa_shm_destroy(void);
	int       bwa_shm_list(void);                                    /* prints staged indices to stdout */
	int       main_shm(int argc, char *argv[]);                      /* CLI entry — see src/main.cpp dispatch */

	/* Attach helpers used by FMI_search / indexEle / fastmap.
	 * *len_out is the packed index size (the header's total_size), the bound
	 * for every section lookup. *map_len_out (optional) is the length of the
	 * mapping, which can be longer and is what bwa_shm_detach takes. An
	 * attached index is released with bwa_shm_detach, never with a bare
	 * munmap: a 1 GB index is a System V attachment. */
	/* need_pac: the caller loads the packed reference (PAC), so a stage
	 * without a usable one (a seed-only `shm --meth` stage) is declined (NULL)
	 * rather than attached for a load that would exit on it.
	 * cp_layout: a BWA_SHM_CP_LAYOUT_* the caller can use, so a segment of the
	 * other layout is declined (NULL, with a note); -1 takes either, and the
	 * caller reads which with bwa_shm_segment_compact. */
	uint8_t  *bwa_shm_attach(const char *prefix, size_t *len_out,
	                         size_t *map_len_out = NULL,
	                         bool need_pac = false,
	                         int cp_layout = -1);                    /* NULL if not staged */
	void      bwa_shm_detach(uint8_t *base, size_t map_len);
	int       bwa_shm_section_find(const uint8_t *base, uint32_t kind,
	                               uint64_t *offset_out, uint64_t *size_out);
	/* 1 if the attached segment at `base` holds the compact checkpoint layout
	 * (header cp_layout BWA_SHM_CP_LAYOUT_COMPACT), 0 if the default one. */
	int       bwa_shm_segment_compact(const uint8_t *base);

	/* The huge-page warning that `shm` (after the pack) and bwa_shm_attach (for a
	 * segment it returns) print, at most once per process, when the segment is not
	 * on huge pages. Exposed for its unit tests (test/unit/test_shm_hugepage_warn.cpp);
	 * nothing else calls them.
	 * bwa_shm_hugepage_warn_wanted: whether to warn for a segment of seg_len bytes.
	 * on_huge is bwa_smaps_vma_finish's verdict for its mapping (bwa_thp.h);
	 * pmd_bytes the PMD (transparent huge page) size; shmem_enabled the text of
	 * /sys/kernel/mm/transparent_hugepage/shmem_enabled, or NULL when it cannot be
	 * read. 1 only when the remount the warning advises could put the segment on
	 * huge pages: the mapping is resident and none of it is huge (on_huge == 0),
	 * the segment holds at least one PMD page (a smaller one is never PMD-mapped,
	 * whatever the mount), and the kernel has shmem THP not set to [deny] (which
	 * turns it off for every tmpfs mount, whatever the mount option).
	 * bwa_shm_format_hugepage_warning: the one-line "[W::<func>]" message for a
	 * segment on page_bytes pages, newline-terminated, into buf (n bytes; NULL
	 * and 0 to size it). Returns the snprintf length. */
	int       bwa_shm_hugepage_warn_wanted(int on_huge, uint64_t seg_len, uint64_t pmd_bytes,
	                                       const char *shmem_enabled);
	int       bwa_shm_format_hugepage_warning(char *buf, size_t n, const char *func,
	                                          const char *name, uint64_t page_bytes,
	                                          uint64_t pmd_bytes);

	/* Layout of a packed bwa-mem3 index segment. Computed by bwa_shm_compute
	 * from a prefix on disk; consumed by bwa_shm_pack_into and bwa_shm_stage.
	 * Owns a heap-loaded bntseq_t held across compute->pack_into so we don't
	 * read .ann/.amb twice. Free with bwa_shm_layout_free. */
	typedef struct {
		char     prefix[PATH_MAX];       /* matches the disk loader's PATH_MAX buffers */
		void    *bns;                    /* bntseq_t *; void* to keep bntseq.h optional here */
		int64_t  reference_seq_len;
		int64_t  count[5];               /* +1-adjusted, ready to write */
		int64_t  sentinel_index;
		int64_t  sa_compx;                /* SA sample-rate shift STAGED into shm. Equals
		                                   * disk_sa_compx unless `shm -u` densified it. */
		int64_t  disk_sa_compx;           /* SA shift on disk (tail-detected). Source stride
		                                   * for a densify pass; == sa_compx when not densifying. */
		int       densify_threads;        /* worker count for the SA densify pass and the
		                                   * compact conversion (>=1); 1 == serial. */
		bool      compact;                /* stage the compact checkpoint layout */
		int64_t   cp_span_shift;          /* compact only: rows per count span, log2 */
		int64_t   cp_wide;                /* compact only: counts are span-relative */
		int64_t  ref_string_len;         /* file size of <prefix>.0123 */
		uint64_t total_size;
		uint32_t n_sections;
		bwa_shm_section_t sections[16];  /* 10 (default layout) or 12 (compact) used */
		uint64_t stager_pid;             /* header reserved[0]; 0 outside bwa_shm_stage */
		uint64_t stage_token;            /* header reserved[1]; 0 outside bwa_shm_stage */
	} bwa_shm_layout_t;

	/* Compute layout and load BNS. Returns 0 on success, -1 on error. When
	 * bns_only=true, the PAC and REF_STRING sections are sized to zero and the
	 * `.0123` is not stat'd (the seed `.0123` need not exist) — see bwa_shm_stage. */
	int  bwa_shm_compute(const char *prefix, bwa_shm_layout_t *layout,
	                     bool bns_only = false, int target_sa_compx = -1,
	                     int n_threads = 1, bool compact = false);

	/* Pack the index described by `layout` into `dest`, which must be at least
	 * layout->total_size bytes. Compact layouts require 32-byte alignment of
	 * dest; a misaligned compact destination is rejected before any writes.
	 * Default layouts retain their natural alignment requirements (8 bytes).
	 * Streams cp_occ / sa_* / pac / ref_string from disk straight into dest
	 * (no intermediate heap copy). Returns 0/-1. */
	int  bwa_shm_pack_into(const bwa_shm_layout_t *layout, uint8_t *dest);

	/* Release the layout's loaded BNS. Safe to call after compute fails. */
	void bwa_shm_layout_free(bwa_shm_layout_t *layout);

	/* Internal: load index components from disk and pack them into a single
	 * self-describing buffer, as bwa_shm_stage packs a segment. Used by
	 * tests. On success, *buf_out is heap-allocated, with the padding
	 * between sections zero as in a staged segment, and *len_out is its
	 * size; caller frees with free(). Returns 0 on success, -1 on error. */
	int bwa_shm_pack_from_disk(const char *prefix, uint8_t **buf_out, uint64_t *len_out);

#ifdef __cplusplus
}
#endif

#endif /* BWA_SHM_H_ */
