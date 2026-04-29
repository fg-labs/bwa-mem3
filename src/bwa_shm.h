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
#define BWA_SHM_VERSION        1u
#define BWA_SHM_CTL_NAME       "/bwactl"
#define BWA_SHM_IDX_PREFIX     "/bwaidx-"
#define BWA_SHM_CTL_SIZE       0x10000   /* 64 KiB. Must equal bwa.h::BWA_CTL_SIZE
                                          * (the legacy v1 stub); a static_assert in
                                          * bwa_shm.cpp will enforce equality once
                                          * that translation unit is added. */

/* FMI_SCALARS section: int64_t reference_seq_len, count[5], sentinel_index. */
#define BWA_SHM_FMI_SCALARS_BYTES (sizeof(int64_t) * 7)

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

typedef struct {
	uint64_t magic;          /* BWA_SHM_MAGIC */
	uint32_t version;        /* BWA_SHM_VERSION */
	uint32_t n_sections;
	uint64_t total_size;     /* bytes including this header */
	uint64_t reserved[5];    /* zero-fill for forward-compat */
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

	int       bwa_shm_test(const char *prefix);                      /* 1 if staged, 0 if not, -1 on error */
	int       bwa_shm_stage(const char *prefix);                     /* loads from disk, packs, stages */
	int       bwa_shm_destroy(void);                                 /* drops all (matches v1 -d) */
	int       bwa_shm_list(void);                                    /* prints staged indices to stdout */
	int       main_shm(int argc, char *argv[]);                      /* CLI entry — see src/main.cpp dispatch */

	/* Attach helpers used by FMI_search / indexEle / fastmap. */
	uint8_t  *bwa_shm_attach(const char *prefix, size_t *len_out);   /* NULL if not staged */
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
		int64_t  ref_string_len;         /* file size of <prefix>.0123 */
		uint64_t total_size;
		uint32_t n_sections;
		bwa_shm_section_t sections[16];  /* 10 used today; spare for forward-compat */
	} bwa_shm_layout_t;

	/* Compute layout and load BNS. Returns 0 on success, -1 on error. */
	int  bwa_shm_compute(const char *prefix, bwa_shm_layout_t *layout);

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

#ifdef __cplusplus
}
#endif

#endif /* BWA_SHM_H_ */
