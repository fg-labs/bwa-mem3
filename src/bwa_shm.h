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
 * Increment whenever any packed section layout (e.g. BWA_SHM_FMI_SCALARS_BYTES)
 * changes. */
#define BWA_SHM_VERSION        2u
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

/* A segment is complete once `magic` is BWA_SHM_MAGIC. The stager writes the
 * header last and publishes `magic` with a release store, so a reader that
 * loads it with acquire ordering and sees BWA_SHM_MAGIC sees the whole index.
 * While staging, `magic` is 0, reserved[0] holds the stager's PID and
 * reserved[1] a random per-stage token; the stager keeps both in the final
 * header. A segment with magic 0 whose stager is no longer alive was left by
 * an interrupted stage. */
typedef struct {
	uint64_t magic;          /* BWA_SHM_MAGIC, or 0 while staging */
	uint32_t version;        /* BWA_SHM_VERSION */
	uint32_t n_sections;
	uint64_t total_size;     /* bytes including this header */
	uint64_t reserved[5];    /* [0] stager PID, [1] stage token; rest zero */
} bwa_shm_header_t;           /* 64 bytes; first object in every segment */

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
	 * independent); <=0 or 1 runs it serially. Ignored when not densifying. */
	/* hugetlb_dir: NULL (default) stages into POSIX shared memory. Non-NULL
	 * asks for the index on 1 GB hugetlb pages (Linux only; see "1 GB hugetlb
	 * staging" below): "" finds a writable hugetlbfs mount with pagesize=1G in
	 * /proc/mounts, any other value names such a mount. When that is not
	 * possible the index is staged in POSIX shared memory as if hugetlb_dir
	 * were NULL, with an [M::] line saying why. */
	int       bwa_shm_stage(const char *prefix, bool bns_only = false,
	                        int target_sa_compx = -1,
	                        int n_threads = 1,
	                        const char *hugetlb_dir = NULL); /* loads from disk, packs, stages */
	/* Test-only hook (see bwa_shm.cpp): bwa_shm_stage calls it with its
	 * writable index mapping after it has registered the segment and released
	 * the stager lock, before the pack fills it, so a unit test can inspect
	 * that mapping or act while the stage is in progress. NULL on every
	 * production path. */
	extern void (*bwa_shm_stage_map_hook)(const void *base, size_t len);
	int       bwa_shm_destroy(void);                                 /* drops all (matches v1 -d) */
	int       bwa_shm_list(void);                                    /* prints staged indices to stdout */
	int       main_shm(int argc, char *argv[]);                      /* CLI entry — see src/main.cpp dispatch */

	/* Attach helpers used by FMI_search / indexEle / fastmap.
	 * *len_out is the packed index size (the header's total_size), the bound
	 * for every section lookup. *map_len_out (optional) is the length of the
	 * mapping, which is what the caller passes to munmap. */
	uint8_t  *bwa_shm_attach(const char *prefix, size_t *len_out,
	                         size_t *map_len_out = NULL);            /* NULL if not staged */
	int       bwa_shm_section_find(const uint8_t *base, uint32_t kind,
	                               uint64_t *offset_out, uint64_t *size_out);

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
		int       densify_threads;        /* worker count for the SA densify pass (>=1);
		                                   * 1 == serial. Ignored when not densifying. */
		int64_t  ref_string_len;         /* file size of <prefix>.0123 */
		uint64_t total_size;
		uint32_t n_sections;
		bwa_shm_section_t sections[16];  /* 10 used today; spare for forward-compat */
		uint64_t stager_pid;             /* header reserved[0]; 0 outside bwa_shm_stage */
		uint64_t stage_token;            /* header reserved[1]; 0 outside bwa_shm_stage */
	} bwa_shm_layout_t;

	/* Compute layout and load BNS. Returns 0 on success, -1 on error. When
	 * bns_only=true, the PAC and REF_STRING sections are sized to zero and the
	 * `.0123` is not stat'd (the seed `.0123` need not exist) — see bwa_shm_stage. */
	int  bwa_shm_compute(const char *prefix, bwa_shm_layout_t *layout,
	                     bool bns_only = false, int target_sa_compx = -1,
	                     int n_threads = 1);

	/* Pack the index described by `layout` into `dest`, which must be at least
	 * layout->total_size bytes. Streams cp_occ / sa_* / pac / ref_string from
	 * disk straight into dest (no intermediate heap copy). Returns 0/-1. */
	int  bwa_shm_pack_into(const bwa_shm_layout_t *layout, uint8_t *dest);

	/* Release the layout's loaded BNS. Safe to call after compute fails. */
	void bwa_shm_layout_free(bwa_shm_layout_t *layout);

	/* Internal: load index components from disk and pack them into a single
	 * self-describing buffer. Used by bwa_shm_stage (Phase 2) and exercised
	 * by tests. On success, *buf_out is malloc'd and *len_out is its size;
	 * caller frees with free(). Returns 0 on success, -1 on error. */
	int bwa_shm_pack_from_disk(const char *prefix, uint8_t **buf_out, uint64_t *len_out);

	/* ---- 1 GB hugetlb staging (`bwa-mem3 shm --huge-pages`) -------------
	 *
	 * /dev/shm is a tmpfs, which backs a segment with 2 MB pages at most. To
	 * put a staged index on 1 GB pages, the stager writes it to the file
	 * `<dir>/bwaidx-<name>` on a hugetlbfs mount with pagesize=1G, and the
	 * POSIX segment `/bwaidx-<name>` holds a redirect stub
	 * (bwa_shm_hugetlb_stub_t) naming that file. The registry entry is the
	 * usual one, so `mem`, `shm -l` and `shm -d` find the index by the same
	 * name, and bwa_shm_attach follows the stub.
	 *
	 * The stub's first 64 bytes read as a bwa_shm_header_t: magic 0 while
	 * staging, the stager PID and stage token where the header keeps them, and
	 * its own magic published last. A redirect stage is therefore waited on,
	 * reclaimed and checked for ownership exactly like any other. The stub is
	 * written before the file is created and removed after the file, so every
	 * file a stage creates is named by a registered stub, even one left by a
	 * killed stage, and a restage or `shm -d` frees its pages.
	 *
	 * A stub is followed only to `bwaidx-<name>` for its own registry name, as
	 * a regular file (not a symlink) owned by the stub's owner in a directory
	 * on a hugetlbfs mount with 1 GB pages; the attach also requires the
	 * file's header to carry the stub's stage token and size. */
#define BWA_SHM_HUGETLB_STUB_MAGIC   0x0045475548415742ull /* "BWAHUGE\0" little-endian */
#define BWA_SHM_HUGETLB_STUB_VERSION 1u
#define BWA_SHM_HUGETLB_KIND_1G      1u   /* non-zero marks a stub while staging */
#define BWA_SHM_HUGETLB_PATH_MAX     4096 /* fixed, so the stub size does not vary by OS */
#define BWA_SHM_HUGETLB_FILE_PREFIX  "bwaidx-"

	typedef struct {
		uint64_t magic;          /* 0 while staging; BWA_SHM_HUGETLB_STUB_MAGIC once complete */
		uint32_t version;        /* BWA_SHM_HUGETLB_STUB_VERSION */
		uint32_t kind;           /* BWA_SHM_HUGETLB_KIND_1G (header n_sections) */
		uint64_t total_size;     /* packed index bytes (header total_size) */
		uint64_t stager_pid;     /* header reserved[0] */
		uint64_t stage_token;    /* header reserved[1] */
		uint64_t file_size;      /* header reserved[2]: the file's size, whole 1 GiB
		                          * pages; 0 until the stager has created the file */
		uint64_t reserved[2];    /* zero */
		char     path[BWA_SHM_HUGETLB_PATH_MAX]; /* absolute path of the file */
	} bwa_shm_hugetlb_stub_t;     /* 64 + 4096 bytes */

	/* A complete stub (magic set, PID and token zero) for `path` holding an
	 * index of `total_size` bytes in a file of `file_size` bytes. Returns 0,
	 * or -1 when the fields would not decode. Pure. */
	int bwa_shm_hugetlb_stub_encode(const char *path, uint64_t total_size,
	                                uint64_t file_size, bwa_shm_hugetlb_stub_t *out);

	/* Validate `len` bytes at `buf` as a complete stub and copy it to *out.
	 * Returns 0 when the length, magic, version and kind match, the path is
	 * absolute and NUL-terminated, and 0 < total_size <= file_size with
	 * file_size a whole number of 1 GiB pages; -1 otherwise. The binding of
	 * the path to a registry name is the caller's check. Pure. */
	int bwa_shm_hugetlb_stub_decode(const void *buf, size_t len,
	                                bwa_shm_hugetlb_stub_t *out);

	/* Parse one line of /proc/mounts ("<src> <dir> <fstype> <opts> ..."). If it
	 * is a hugetlbfs mount whose page size is 1 GiB (pagesize=1024M, 1G or
	 * 1048576K), copy its directory, with the octal escapes of /proc/mounts
	 * (\040 for a space, ...) decoded, into dir_out and return 1. Returns 0 for
	 * any other line, a malformed one, or a directory that does not fit in
	 * dir_len bytes. Pure. */
	int bwa_shm_parse_hugetlbfs_1g_mount(const char *line, char *dir_out, size_t dir_len);

#ifdef __cplusplus
}

static_assert(sizeof(bwa_shm_header_t) == 64, "index header is 64 bytes");
static_assert(sizeof(bwa_shm_hugetlb_stub_t) == 64 + BWA_SHM_HUGETLB_PATH_MAX,
              "stub is a 64-byte header and a fixed-size path");
static_assert(offsetof(bwa_shm_hugetlb_stub_t, kind) ==
              offsetof(bwa_shm_header_t, n_sections), "stub kind overlays n_sections");
static_assert(offsetof(bwa_shm_hugetlb_stub_t, total_size) ==
              offsetof(bwa_shm_header_t, total_size), "stub total_size overlays the header's");
static_assert(offsetof(bwa_shm_hugetlb_stub_t, stager_pid) ==
              offsetof(bwa_shm_header_t, reserved), "stub PID overlays header reserved[0]");
static_assert(offsetof(bwa_shm_hugetlb_stub_t, stage_token) ==
              offsetof(bwa_shm_header_t, reserved) + sizeof(uint64_t),
              "stub token overlays header reserved[1]");
static_assert(offsetof(bwa_shm_hugetlb_stub_t, file_size) ==
              offsetof(bwa_shm_header_t, reserved) + 2 * sizeof(uint64_t),
              "stub file_size overlays header reserved[2]");
static_assert(offsetof(bwa_shm_hugetlb_stub_t, path) == sizeof(bwa_shm_header_t),
              "stub path follows the 64-byte header");
static_assert(PATH_MAX <= BWA_SHM_HUGETLB_PATH_MAX, "a stub holds any path");
static_assert(BWA_SHM_VERSION != BWA_SHM_HUGETLB_STUB_VERSION,
              "a staging stub is told from a staging index by its version");
#endif

#endif /* BWA_SHM_H_ */
