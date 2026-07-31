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
#include "stage_prof.h"
#include "kvec.h"
#include "utils.h"
#include "bntseq.h"
#include "kseq.h"
#include "profiling.h"

KSEQ_DECLARE(gzFile)

/* Forward-declared — keeps htslib out of this header. */
struct bam_writer_s;

#include "fast_reader.h"

typedef struct {
	kseq_t *ks, *ks2;
	/* Fast-path reader (default). When legacy_reader is set, the gzFile/kseq
	 * path (ks, ks2, fp) is used; otherwise the fast_reader path below. */
	int legacy_reader;
	fast_reader_t *fr1, *fr2;
	void *frks, *frks2;   /* kseq_t* over fast_reader; opaque in this header */
	mem_opt_t *opt;
	mem_pestat_t *pes0;
	int64_t n_processed;
	int copy_comment;
	int64_t my_ntasks;
	int64_t ntasks;
	int64_t task_size;
	int64_t actual_chunk_size;
	/* --- pestat cohort slicing -------------------------------------------
	 * A "cohort" is one task_size batch: exactly the read set bwa-mem2 would
	 * hand mem_pestat. It may be read as several smaller slices so that the
	 * first read of a run does not stall the compute pipeline, but the cohort
	 * boundary itself is never moved -- that is what keeps output identical.
	 *
	 * Owned by the read step (step 0) and the process step (step 1), both of
	 * which kt_pipeline serialises, so none of this needs a lock. */
	int64_t cohort_slices;      /* slices for the FIRST cohort; 0 = off (one slice) */
	int     cohort_slice_all;   /* stress knob: slice EVERY cohort (identity testing) */
	int64_t cohort_bases;       /* bases of the current cohort read so far */
	int64_t cohort_index;       /* how many cohorts have been started */
	int     cohort_slice;       /* slice number within the current cohort */
	/* The in-flight cohort's string arena, carried across its slices. Every slice
	 * carves name/seq/qual from this one arena; the accumulator shares those
	 * bytes by pointer, so they must live until the whole cohort is written. The
	 * completing slice hands it to the write stage (which destroys it) and this
	 * goes back to NULL so the next cohort creates a fresh one. See read_arena.h
	 * for why the cohort, not the read call, is the arena's unit. */
	read_arena_t *cohort_arena;
	/* How fast the ramp grows per slice, and the previous slice's requested size
	 * to grow from. The ratio is what decides whether step 1 ever waits: while
	 * step 1 computes slice k the reader must deliver slice k+1, so the ramp is
	 * only free while ratio <= (read seconds per base) / (compute seconds per
	 * base). Both are step-0-only state, which kt_pipeline serialises. */
	double  cohort_ramp_ratio;  /* growth per slice; see --cohort-ramp-ratio */
	int64_t cohort_ramp_first;  /* first slice, in bases; 0 = fraction of task_size */
	int64_t ramp_prev_target;   /* previous slice's REQUESTED size, in bases */
	/* Accumulator for a multi-slice cohort, owned by step 1. */
	bseq1_t *cohort_seqs;       /* contiguous reads for the cohort being built */
	int      cohort_n;          /* reads accumulated so far */
	int      cohort_cap;        /* capacity of cohort_seqs, in reads */
	int64_t  cohort_first_id;   /* global id of the cohort's first read */
	FILE *fp;
	uint8_t *ref_string;
	int      ref_string_is_shm;   /* 1 if ref_string aliases shm pages; do not _mm_free. */
	FMI_search *fmi;
	uint8_t *shm_base;            /* if non-NULL, the active /bwaidx-<base> mapping */
	/* Generic (non-meth) BAM writer: non-NULL only when opt->bam_mode is set
	 * AND opt->meth_mode is not -- under --meth the meth writer owns output
	 * (g_meth_bam_writer) in both containers, and this stays NULL. */
	struct bam_writer_s *bam_writer;
	/* D3 (--meth) ORIGINAL-reference handles, resident alongside (and distinct
	 * from) the seed FM-index in `fmi`. The seed index (`fmi->idx->bns/pac`) is
	 * the f/r-doubled converted reference used for candidate generation; these
	 * are the un-converted original `.bns`/`.pac` (real chrom names, N contigs)
	 * loaded for the future extension/scoring phase (D3 spec §5.2/§6). Both must
	 * stay resident. NULL outside --meth (and when no original prefix resolved).
	 * Load-only in this step — no consumer yet (extension is behind the
	 * meth-mode checkpoint). */
	bntseq_t *meth_orig_bns;
	uint8_t  *meth_orig_pac;
	/* D3 (--meth, PR-3) ORIGINAL unpacked reference (`<orig>.0123`), the
	 * extension/dedup ref bases for the (original-coord) remapped seeds. Distinct
	 * from `ref_string` above, which is the SEED `.0123`. Always heap-allocated
	 * via _mm_malloc here (never shm) and _mm_free'd on teardown. NULL outside
	 * --meth (or when no original prefix resolved). */
	uint8_t  *meth_orig_ref_string;
} ktp_aux_t;

typedef struct {
	ktp_aux_t *aux;
	int n_seqs;
	bseq1_t *seqs;
	/* PIPE-F6: per-chunk bump arena backing every read's name/seq/qual fields.
	 * Set by the reader in step 0; destroyed once in step 2 after the chunk's
	 * output is written. NULL when the chunk carried no reads (clean EOF). */
	read_arena_t *read_arena;
	prof_chunk_t prof;   /* stage_prof: per-chunk read/process/write timing (--profile) */
	/* Set by the read step when this slice completes a pestat cohort (either it
	 * reached task_size or the input ran out). Until then the process step
	 * aligns the slice and returns an empty item, which the pipeline passes
	 * through harmlessly -- returning NULL would retire the worker. */
	int cohort_complete;
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
//
// Buffer lifetimes, which callers driving the kernels themselves MUST respect:
//   * `w.regs` is sized by `nreads` and indexed by the read's position in the
//     chunk (`seq_id`). It must stay live from the align pass through pairing
//     and SAM emission.
//   * `w.chain_scratch` and `w.seed_scratch` are PER-THREAD scratch, sized by
//     `nthreads * BATCH_SIZE` and indexed by `tid`, NOT by `seq_id`. Thread
//     `tid` owns `chain_scratch[tid * BATCH_SIZE ..]` and
//     `seed_scratch[tid * BATCH_SIZE * AVG_SEEDS_PER_READ ..]` and may reuse
//     that window for every work item it runs, because mem_kernel2_core frees
//     every chain in the item before returning. Indexing either by `seq_id`
//     (as bwa-mem3 did before these buffers were shrunk) now runs off the end
//     of the allocation and corrupts other threads' windows.
void worker_alloc(const mem_opt_t *opt, worker_t &w, int32_t nreads, int32_t nthreads);

// Release all per-worker scratch buffers previously allocated by
// worker_alloc. `nthreads` must match the value passed to the paired
// worker_alloc call (asserted against `w.nthreads`) so the per-thread loops
// iterate over exactly the slots that were populated. Exposed as the public
// counterpart to worker_alloc so external consumers don't need to duplicate
// the teardown logic.
void worker_free(worker_t &w, int32_t nthreads);

#endif
