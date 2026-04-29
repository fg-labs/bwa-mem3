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

#ifndef FASTMAP_HPP
#define FASTMAP_HPP

#include <ctype.h>
#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include <fstream>
#include "bwa.h"
#include "bwamem.h"
#include "kthread.h"
#include "kvec.h"
#include "utils.h"
#include "bntseq.h"
#include "kseq.h"
#include "profiling.h"

KSEQ_DECLARE(gzFile)

/* Forward-declared — keeps htslib out of this header. */
struct bam_writer_s;

typedef struct {
	kseq_t *ks, *ks2;
	mem_opt_t *opt;
	mem_pestat_t *pes0;
	int64_t n_processed;
	int copy_comment;
	int64_t my_ntasks;
	int64_t ntasks;
	int64_t task_size;
	int64_t actual_chunk_size;
	FILE *fp;
	uint8_t *ref_string;
	int      ref_string_is_shm;   /* 1 if ref_string aliases shm pages; do not _mm_free. */
	FMI_search *fmi;
	uint8_t *shm_base;            /* if non-NULL, the active /bwaidx-<base> mapping */
	struct bam_writer_s *bam_writer;  /* non-NULL when opt->bam_mode is set */
} ktp_aux_t;

typedef struct {
	ktp_aux_t *aux;
	int n_seqs;
	bseq1_t *seqs;
} ktp_data_t;

    
void *kopen(const char *fn, int *_fd);
int kclose(void *a);
int main_mem(int argc, char *argv[]);

// Allocate all per-worker scratch buffers (chaining arrays, BSW buffers, BWT
// scratch) on `w`, sized for `nreads` and `nthreads`. Exposed so that
// consumers of libbwa.a that build their own worker_t (e.g. language
// bindings) can reuse the exact same allocation layout as the bwa-mem3
// pipeline and stay in sync across future changes. Records `nthreads` on
// `w.nthreads` so the matching worker_free can validate the pairing. Prints
// a small summary to stderr. Asserts on allocation failure (matching the
// rest of bwa-mem3).
void worker_alloc(const mem_opt_t *opt, worker_t &w, int32_t nreads, int32_t nthreads);

// Release all per-worker scratch buffers previously allocated by
// worker_alloc. `nthreads` must match the value passed to the paired
// worker_alloc call (asserted against `w.nthreads`) so the per-thread loops
// iterate over exactly the slots that were populated. Exposed as the public
// counterpart to worker_alloc so external consumers don't need to duplicate
// the teardown logic.
void worker_free(worker_t &w, int32_t nthreads);

#endif
