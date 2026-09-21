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

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>.
*****************************************************************************************/

#ifndef __INDEX_ELE_HPP
#define __INDEX_ELE_HPP

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "utils.h"
#include "bntseq.h"
#include "macro.h"

#define BWA_IDX_BWT 0x1
#define BWA_IDX_BNS 0x2
#define BWA_IDX_PAC 0x4
#define BWA_IDX_ALL 0x7

typedef struct {
	bntseq_t *bns; // information on the reference sequences
	uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base
	int    is_shm;
	int64_t l_mem;
	uint8_t  *mem;
} bwaidx_fm_t;


class indexEle {
	
public:
	bwaidx_fm_t *idx;
	
	indexEle();
	~indexEle();
	/* Load BNS (and, when `which` includes BWA_IDX_PAC, the 2-bit packed
	 * reference) from disk. `pread_workers` is the already-resolved worker count
	 * for the parallel .pac slurp (as returned by index_load_threads); 1 reads
	 * serially. The loaded bytes are identical regardless of worker count. */
	void bwa_idx_load_ele(const char *hint, int which, int pread_workers = 1);
	char *bwa_idx_infer_prefix(const char *hint);

	/* Attach BNS + PAC from a packed bwa-mem3 index segment produced by
	 * bwa_shm_pack_from_disk. `base`/`len` describe the segment (typically
	 * from bwa_shm_attach). The bntseq_t struct and anns[] array are
	 * heap-copied (mapper mutates them); ambs[], name/anno strings, and
	 * pac are aliased into the segment. Sets idx->is_shm=1 and idx->mem,
	 * idx->l_mem so the existing destructor gate skips the wrong frees.
	 * The segment's lifetime belongs to the loader process; shm pages are
	 * never freed by this object.
	 *
	 * load_pac=false skips the PAC section entirely (idx->pac=NULL): the D3
	 * --meth seed segment is staged bns_only (no PAC section), and `mem --meth`
	 * extends against the ORIGINAL pac, so the seed pac is never needed.
	 */
	void bwa_idx_load_ele_from_shm(uint8_t *base, size_t len, bool load_pac = true);
};

/* Slurp the whole of an open .pac stream (`*fp_pac`, positioned at its start as
 * bns_restore leaves it) into the caller-allocated `dst` (`pac_bytes` bytes),
 * reading in parallel across `pread_workers` pread() workers (a resolved count,
 * as from index_load_threads; 1 reads serially), then close the stream and NULL
 * `*fp_pac`. Advises MADV_HUGEPAGE on `dst` first, matching the FM-index arrays.
 * A read error or short file is fatal (aborts the process). Shared by the
 * seed-index loader (bwa_idx_load_ele) and the --meth original-reference loader
 * so the sizing/read/close logic lives in one place. */
void pac_slurp_and_close(FILE **fp_pac, uint8_t *dst, int64_t pac_bytes, int pread_workers);
#endif
