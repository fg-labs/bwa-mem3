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

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bwa_madvise.h"
#if NUMA_ENABLED
#include <numa.h>
#endif
#include <sstream>
#include <getopt.h>
#include "fastmap.h"
#include "FMI_search.h"
#include "bam_writer.h"
#include "meth_bam.h"
#include "stage_prof.h"
#include "seed_order.h"
#include "version.h"
#include <sys/resource.h>
#include "bwa_shm.h"
#include "fast_reader_bseq.h"

#if AFF && (__linux__)
#include <sys/sysinfo.h>
int affy[256];
#endif

// --------------
extern uint64_t tprof[LIM_R][LIM_C];
// ---------------

#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
/* ARM/Apple Silicon - no CPUID instruction, use system calls */

#ifdef __APPLE__
#include <sys/sysctl.h>
#endif

void __cpuid(unsigned int i, unsigned int cpuid[4]) {
    /* ARM doesn't have CPUID - return zeros */
    cpuid[0] = cpuid[1] = cpuid[2] = cpuid[3] = 0;
}

/* Get L2 cache size in bytes (for dynamic tuning) */
static int64_t get_l2_cache_size() {
#ifdef __APPLE__
    int64_t l2_size = 0;
    size_t size = sizeof(l2_size);
    if (sysctlbyname("hw.l2cachesize", &l2_size, &size, NULL, 0) == 0) {
        return l2_size;
    }
#endif
    return 4 * 1024 * 1024;  /* Default 4MB */
}

/* Get optimal batch size based on L2 cache */
int get_dynamic_batch_size() {
    int64_t l2_size = get_l2_cache_size();
    /* Heuristic: ~1KB working set per read pair, use 1/4 of L2 for batching
     * to leave room for other data structures */
    int batch_size = (int)(l2_size / 4 / 1024);
    /* Clamp to reasonable range */
    if (batch_size < 256) batch_size = 256;
    if (batch_size > 4096) batch_size = 4096;
    /* Round down to power of 2 for alignment */
    int pow2 = 256;
    while (pow2 * 2 <= batch_size) pow2 *= 2;
    return pow2;
}

int HTStatus()
{
    /* ARM/Apple Silicon doesn't have hyperthreading in the x86 sense.
     * Apple Silicon has P-cores and E-cores, which we handle differently.
     * Return 0 to indicate no HT (we handle core types via QoS instead).
     */
#ifdef __APPLE__
    int pcore_count = 0, ecore_count = 0;
    size_t size = sizeof(int);

    /* Try to get P-core and E-core counts on Apple Silicon */
    if (sysctlbyname("hw.perflevel0.physicalcpu", &pcore_count, &size, NULL, 0) == 0) {
        sysctlbyname("hw.perflevel1.physicalcpu", &ecore_count, &size, NULL, 0);
        fprintf(stderr, "Platform vendor: Apple Silicon.\n");
        fprintf(stderr, "P-cores: %d, E-cores: %d\n", pcore_count, ecore_count);
    } else {
        /* Fallback for older macOS or non-Apple Silicon */
        int total_cores = 0;
        sysctlbyname("hw.physicalcpu", &total_cores, &size, NULL, 0);
        fprintf(stderr, "Platform: ARM64, Physical CPUs: %d\n", total_cores);
    }

    /* Report L2 cache and dynamic batch size */
    int64_t l2_size = get_l2_cache_size();
    int dynamic_batch = get_dynamic_batch_size();
    fprintf(stderr, "L2 Cache: %lld MB, Dynamic batch size: %d (compile-time: %d)\n",
            l2_size / (1024 * 1024), dynamic_batch, BATCH_SIZE);
#else
    fprintf(stderr, "Platform vendor: ARM64.\n");
#endif
    return 0;  /* No hyperthreading on ARM */
}

#else
/* x86 CPUID implementation */

void __cpuid(unsigned int i, unsigned int cpuid[4]) {
#ifdef _WIN32
    __cpuid((int *) cpuid, (int)i);

#else
    asm volatile
        ("cpuid" : "=a" (cpuid[0]), "=b" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
            : "0" (i), "2" (0));
#endif
}


int HTStatus()
{
    unsigned int cpuid[4];
    char platform_vendor[12];
    __cpuid(0, cpuid);
    ((unsigned int *)platform_vendor)[0] = cpuid[1]; // B
    ((unsigned int *)platform_vendor)[1] = cpuid[3]; // D
    ((unsigned int *)platform_vendor)[2] = cpuid[2]; // C
    std::string platform = std::string(platform_vendor, 12);

    __cpuid(1, cpuid);
    unsigned int platform_features = cpuid[3]; //D

    // __cpuid(1, cpuid);
    unsigned int num_logical_cpus = (cpuid[1] >> 16) & 0xFF; // B[23:16]
    // fprintf(stderr, "#logical cpus: ", num_logical_cpus);

    unsigned int num_cores = -1;
    if (platform == "GenuineIntel") {
        __cpuid(4, cpuid);
        num_cores = ((cpuid[0] >> 26) & 0x3f) + 1; //A[31:26] + 1
        fprintf(stderr, "Platform vendor: Intel.\n");
    } else  {
        fprintf(stderr, "Platform vendor unknown.\n");
    }

    // fprintf(stderr, "#physical cpus: ", num_cores);

    int ht = platform_features & (1 << 28) && num_cores < num_logical_cpus;
    if (ht)
        fprintf(stderr, "CPUs support hyperThreading !!\n");

    return ht;
}

#endif /* ARM vs x86 */


/*** Memory pre-allocations ***/
// Core allocation routine, parameterised only on mem_opt_t. Exposed so that
// library consumers which build up a worker_t themselves (e.g. language
// bindings that don't construct a ktp_aux_t) can reuse the exact same
// allocation sequence instead of re-implementing it and drifting out of sync.
void worker_alloc(const mem_opt_t *opt, worker_t &w, int32_t nreads, int32_t nthreads)
{
    assert(opt != NULL);
    assert(nreads >= 0);
    assert(nthreads > 0);

    // Record the thread count on the worker so worker_free can validate the
    // paired call and the per-thread loops can never walk past the slots
    // populated here.
    w.nthreads = nthreads;

    int32_t memSize = nreads;

    /* Mem allocation section for core kernels */
    w.regs = NULL; w.chain_ar = NULL; w.seedBuf = NULL;

    w.regs = (mem_alnreg_v *) calloc(memSize, sizeof(mem_alnreg_v));
    w.chain_ar = (mem_chain_v*) malloc (memSize * sizeof(mem_chain_v));
    w.seedBuf = (mem_seed_t *) calloc(sizeof(mem_seed_t),  memSize * AVG_SEEDS_PER_READ);

    assert(w.seedBuf  != NULL);
    assert(w.regs     != NULL);
    assert(w.chain_ar != NULL);

    w.seedBufSize = BATCH_SIZE * AVG_SEEDS_PER_READ;

    /*** printing ***/
    int64_t allocMem = memSize * sizeof(mem_alnreg_v) +
        memSize * sizeof(mem_chain_v) +
        sizeof(mem_seed_t) * memSize * AVG_SEEDS_PER_READ;
    fprintf(stderr, "------------------------------------------\n");
    fprintf(stderr, "1. Memory pre-allocation for Chaining: %0.4lf MB\n", allocMem/1e6);


    /* SWA mem allocation */
    int64_t wsize = BATCH_SIZE * SEEDS_PER_READ;
    for(int l=0; l<nthreads; l++)
    {
        w.mmc.seqBufLeftRef[l*CACHE_LINE]  = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN, 64);
        w.mmc.seqBufLeftQer[l*CACHE_LINE]  = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN, 64);
        w.mmc.seqBufRightRef[l*CACHE_LINE] = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN, 64);
        w.mmc.seqBufRightQer[l*CACHE_LINE] = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN, 64);

        w.mmc.wsize_buf_ref[l*CACHE_LINE] = wsize * MAX_SEQ_LEN_REF;
        w.mmc.wsize_buf_qer[l*CACHE_LINE] = wsize * MAX_SEQ_LEN_QER;

        assert(w.mmc.seqBufLeftRef[l*CACHE_LINE]  != NULL);
        assert(w.mmc.seqBufLeftQer[l*CACHE_LINE]  != NULL);
        assert(w.mmc.seqBufRightRef[l*CACHE_LINE] != NULL);
        assert(w.mmc.seqBufRightQer[l*CACHE_LINE] != NULL);
    }

    for(int l=0; l<nthreads; l++) {
        w.mmc.seqPairArrayAux[l]      = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
        w.mmc.seqPairArrayLeft128[l]  = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
        w.mmc.seqPairArrayRight128[l] = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
        w.mmc.wsize[l] = wsize;

        assert(w.mmc.seqPairArrayAux[l] != NULL);
        assert(w.mmc.seqPairArrayLeft128[l] != NULL);
        assert(w.mmc.seqPairArrayRight128[l] != NULL);
    }


    allocMem = (wsize * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN) * opt->n_threads * 2+
        (wsize * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN) * opt->n_threads  * 2 +
        wsize * sizeof(SeqPair) * opt->n_threads * 3;
    fprintf(stderr, "2. Memory pre-allocation for BSW: %0.4lf MB\n", allocMem/1e6);

    // SMEM buffers (matchArray / min_intv_ar / query_pos_ar / enc_qdb / rid)
    // and the lockstep-batch slot buffers are sized from the observed max
    // read length on each batch in mem_collect_smem; they're NULL here and
    // grow on first use. `lim` is still a fixed BATCH_SIZE+32 allocation
    // because its size does not depend on read length.
    for (int l=0; l<nthreads; l++)
    {
        w.mmc.wsize_mem[l]     = 0;
        w.mmc.wsize_mem_s[l]   = 0;
        w.mmc.wsize_mem_r[l]   = 0;
        w.mmc.wsize_qdb[l]     = 0;
        w.mmc.matchArray[l]    = NULL;
        w.mmc.min_intv_ar[l]   = NULL;
        w.mmc.query_pos_ar[l]  = NULL;
        w.mmc.enc_qdb[l]       = NULL;
        w.mmc.rid[l]           = NULL;
        w.mmc.lim[l]           = (int32_t *) _mm_malloc((BATCH_SIZE + 32) * sizeof(int32_t), 64);

        w.mmc.lockstep_prev[l]      = NULL;
        w.mmc.lockstep_match_buf[l] = NULL;
        w.mmc.lockstep_buf_cap[l]   = 0;
    }

    allocMem = nthreads * (BATCH_SIZE + 32) * sizeof(int32_t);
    fprintf(stderr, "3. Memory pre-allocation for BWT (lazy SMEM buffers, %ld B fixed): %0.4lf MB\n",
            (long)(nthreads * (BATCH_SIZE + 32) * sizeof(int32_t)), allocMem/1e6);
    fprintf(stderr, "------------------------------------------\n");
}

// Release every per-worker scratch buffer allocated by worker_alloc. The
// nthreads argument must match the value passed to the paired worker_alloc
// call so the per-thread loops iterate over exactly the slots that were
// populated.
void worker_free(worker_t &w, int32_t nthreads)
{
    assert(nthreads > 0);
    // Catch mismatched alloc/free pairs before they drive out-of-bounds frees.
    assert(w.nthreads == nthreads);

    free(w.chain_ar);
    free(w.regs);
    free(w.seedBuf);

    for(int l=0; l<nthreads; l++) {
        _mm_free(w.mmc.seqBufLeftRef[l*CACHE_LINE]);
        _mm_free(w.mmc.seqBufRightRef[l*CACHE_LINE]);
        _mm_free(w.mmc.seqBufLeftQer[l*CACHE_LINE]);
        _mm_free(w.mmc.seqBufRightQer[l*CACHE_LINE]);
    }

    for(int l=0; l<nthreads; l++) {
        free(w.mmc.seqPairArrayAux[l]);
        free(w.mmc.seqPairArrayLeft128[l]);
        free(w.mmc.seqPairArrayRight128[l]);
    }

    // NULL-safe: SMEM buffers are now allocated lazily on first batch;
    // workers that never ran a batch leave them as NULL. _mm_free / free
    // are both well-defined on NULL.
    for(int l=0; l<nthreads; l++) {
        _mm_free(w.mmc.matchArray[l]);
        free(w.mmc.min_intv_ar[l]);
        free(w.mmc.query_pos_ar[l]);
        free(w.mmc.enc_qdb[l]);
        free(w.mmc.rid[l]);
        _mm_free(w.mmc.lim[l]);

        _mm_free(w.mmc.lockstep_prev[l]);
        _mm_free(w.mmc.lockstep_match_buf[l]);
    }
}

// Back-compat wrapper used by the bwa-mem3 pipeline.
void memoryAlloc(ktp_aux_t *aux, worker_t &w, int32_t nreads, int32_t nthreads)
{
    worker_alloc(aux->opt, w, nreads, nthreads);
}

ktp_data_t *kt_pipeline(void *shared, int step, void *data, mem_opt_t *opt, worker_t &w)
{
    ktp_aux_t *aux = (ktp_aux_t*) shared;
    ktp_data_t *ret = (ktp_data_t*) data;

    if (step == 0)
    {
        ktp_data_t *ret = (ktp_data_t *) calloc(1, sizeof(ktp_data_t));
        assert(ret != NULL);
        uint64_t tim = __rdtsc();
        double sp_r0 = 0.0;
        if (sp_enabled()) {
            sp_chunk_init(&ret->prof); sp_read_reset();
            ret->prof.chunk_start = sp_run_elapsed();   /* timeline anchor */
            sp_r0 = sp_wall();
        }

        /* Read "reads" from input file (fread) */
        int64_t sz = 0;
        ret->seqs = aux->legacy_reader
            ? bseq_read_orig(aux->task_size, &ret->n_seqs, aux->ks, aux->ks2, &sz)
            : bseq_read_fast(aux->task_size, &ret->n_seqs, aux->frks, aux->frks2, &sz);

        tprof[READ_IO][0] += __rdtsc() - tim;

        if (sp_enabled()) {
            ret->prof.read_wall = sp_wall() - sp_r0;
            /* Only the fast reader is instrumented; the legacy reader leaves the
             * sub-splits at their NaN init (blank) rather than reporting a fake 0. */
            if (!aux->legacy_reader) {
                sp_read_get(&ret->prof.read_diskwait, &ret->prof.read_decompress, &ret->prof.read_parse);
                sp_read_get_bytes(&ret->prof.read_bytes_in, &ret->prof.bgzf_blocks);
            }
            ret->prof.n_reads = ret->n_seqs;
            ret->prof.n_bp = sz;
        }

        fprintf(stderr, "[0000] read_chunk: %lld, work_chunk_size: %lld, nseq: %d\n",
                (long long)aux->task_size, (long long)sz, ret->n_seqs);

        if (ret->seqs == 0) {
            free(ret);
            return 0;
        }
        if (!aux->copy_comment){
            for (int i = 0; i < ret->n_seqs; ++i){
                free(ret->seqs[i].comment);
                ret->seqs[i].comment = 0;
            }
        }

        /* Inline bwameth-style c2t conversion on read ingest. Matches
         * `bwameth.py c2t R1.fq R2.fq`: R1 reads get C→T, R2 reads get G→A.
         * The original sequence is stashed as a YS:Z comment tag and the
         * conversion type as YC:Z — these pass through to SAM via
         * copy_comment (which --meth sets).
         *
         * R1/R2 classification:
         *   - Non-smart-pair PE (two FASTQs): records are strictly
         *     interleaved R1/R2/R1/R2… in ret->seqs, so record parity
         *     (i & 1) identifies R2.
         *   - Smart-pair (-p, single stream): the stream can contain
         *     orphans, so parity mis-tags them. Classify by adjacent-name
         *     pairing (same rule as bseq_classify) — if this read's name
         *     matches the previous unmatched read's name, it is R2;
         *     otherwise it starts a new pair as R1.
         *   - SE: always R1. */
        if (aux->opt->meth_mode) {
            int is_pe    = (aux->opt->flag & MEM_F_PE) != 0;
            int is_smart = (aux->opt->flag & MEM_F_SMARTPE) != 0;
            int prev_is_r1 = 0;
            for (int i = 0; i < ret->n_seqs; ++i) {
                bseq1_t *s = &ret->seqs[i];
                int is_r2;
                if (!is_pe) {
                    is_r2 = 0;
                } else if (!is_smart) {
                    is_r2 = (i & 1);
                } else if (prev_is_r1 && i > 0
                           && strcmp(s->name, ret->seqs[i-1].name) == 0) {
                    is_r2 = 1;
                    prev_is_r1 = 0;
                } else {
                    is_r2 = 0;
                    prev_is_r1 = 1;
                }
                const char *yc = is_r2 ? "GA" : "CT";
                char from = is_r2 ? 'G' : 'C';
                char from_lo = is_r2 ? 'g' : 'c';
                char to   = is_r2 ? 'A' : 'T';
                int l = s->l_seq;
                /* Build the YS:Z/YC:Z comment. Preserve any prior FASTQ
                 * comment (e.g. -C carries barcode/UMI SAM tags) by appending
                 * it after YC; otherwise --meth silently strips -C metadata. */
                const char *prior = s->comment;
                size_t prior_len = prior ? strlen(prior) : 0;
                size_t yslen = (size_t)l + 32 + (prior_len ? prior_len + 1 : 0);
                char *comment = (char *)malloc(yslen);
                assert(comment != NULL);
                int off = snprintf(comment, yslen, "YS:Z:");
                memcpy(comment + off, s->seq, (size_t)l);
                off += l;
                off += snprintf(comment + off, yslen - off, "\tYC:Z:%s", yc);
                if (prior_len)
                    snprintf(comment + off, yslen - off, "\t%s", prior);
                free(s->comment);
                s->comment = comment;
                /* Retain the ORIGINAL (unconverted) read bases as a first-class
                 * field BEFORE the in-place projection below overwrites s->seq.
                 * Same orientation/order as s->seq (original read order, ASCII);
                 * downstream consumers must RC it wherever they RC s->seq — see
                 * the bseq1_t.meth_orig_seq orientation contract in bwa.h.
                 * strdup is fine for the draft; freed in the per-batch free
                 * loop below alongside s->seq. */
                s->meth_orig_seq = strdup(s->seq);
                assert(s->meth_orig_seq != NULL);
                /* Project in place. */
                for (int j = 0; j < l; ++j) {
                    char c = s->seq[j];
                    if (c == from || c == from_lo) s->seq[j] = to;
                }
            }
        }
        {
            int64_t size = 0;
            for (int i = 0; i < ret->n_seqs; ++i) size += ret->seqs[i].l_seq;

            fprintf(stderr, "\t[0000][ M::%s] read %d sequences (%ld bp)...\n",
                    __func__, ret->n_seqs, (long)size);
        }

        return ret;
    } // Step 0
    else if (step == 1)  /* Step 2: Main processing-engine */
    {
        static int task = 0;
        if (w.nreads < ret->n_seqs)
        {
            fprintf(stderr, "[0000] Reallocating initial memory allocations!!\n");
            free(w.regs); free(w.chain_ar); free(w.seedBuf);
            w.nreads = ret->n_seqs;
            w.regs = (mem_alnreg_v *) calloc(w.nreads, sizeof(mem_alnreg_v));
            w.chain_ar = (mem_chain_v*) malloc (w.nreads * sizeof(mem_chain_v));
            w.seedBuf = (mem_seed_t *) calloc(sizeof(mem_seed_t), w.nreads * AVG_SEEDS_PER_READ);
            assert(w.regs != NULL); assert(w.chain_ar != NULL); assert(w.seedBuf != NULL);
        }

        fprintf(stderr, "[0000] Calling mem_process_seqs.., task: %d\n", task++);

        double sp_p0 = 0.0;
        if (sp_enabled()) {
            sp_chunk_init(&g_ktfor);          /* reset balance + encode accumulator */
            g_ktfor.encode = 0.0;             /* accumulate (sp_chunk_init left it NaN) */
            sp_p0 = sp_wall();
        }

        uint64_t tim = __rdtsc();
        if (opt->flag & MEM_F_SMARTPE)
        {
            bseq1_t *sep[2];
            int n_sep[2];
            mem_opt_t tmp_opt = *opt;

            bseq_classify(ret->n_seqs, ret->seqs, n_sep, sep);

            fprintf(stderr, "[M::%s] %d single-end sequences; %d paired-end sequences.....\n",
                    __func__, n_sep[0], n_sep[1]);

            if (n_sep[0]) {
                tmp_opt.flag &= ~MEM_F_PE;
                /* single-end sequences, in the mixture */
                mem_process_seqs(&tmp_opt,
                                 aux->n_processed,
                                 n_sep[0],
                                 sep[0],
                                 0,
                                 w);

                for (int i = 0; i < n_sep[0]; ++i) {
                    bseq1_t *dst = &ret->seqs[sep[0][i].id];
                    bseq1_t *src = &sep[0][i];
                    dst->sam      = src->sam;      src->sam      = NULL;
                    dst->bams     = src->bams;     src->bams     = NULL;
                    dst->n_bams   = src->n_bams;   src->n_bams   = 0;
                    dst->cap_bams = src->cap_bams; src->cap_bams = 0;
                }
            }
            if (n_sep[1]) {
                tmp_opt.flag |= MEM_F_PE;
                /* paired-end sequences, in the mixture */
                mem_process_seqs(&tmp_opt,
                                 aux->n_processed + n_sep[0],
                                 n_sep[1],
                                 sep[1],
                                 aux->pes0,
                                 w);

                for (int i = 0; i < n_sep[1]; ++i) {
                    bseq1_t *dst = &ret->seqs[sep[1][i].id];
                    bseq1_t *src = &sep[1][i];
                    dst->sam      = src->sam;      src->sam      = NULL;
                    dst->bams     = src->bams;     src->bams     = NULL;
                    dst->n_bams   = src->n_bams;   src->n_bams   = 0;
                    dst->cap_bams = src->cap_bams; src->cap_bams = 0;
                }
            }
            free(sep[0]); free(sep[1]);
        }
        else {
            /* pure (single/paired-end), reads processing */
            mem_process_seqs(opt,
                             aux->n_processed,
                             ret->n_seqs,
                             ret->seqs,
                             aux->pes0,
                             w);
        }
        tprof[MEM_PROCESS2][0] += __rdtsc() - tim;

        if (sp_enabled()) {
            ret->prof.proc_wall      = sp_wall() - sp_p0;
            ret->prof.proc_cpu       = g_ktfor.proc_cpu;
            ret->prof.thr_busy_min   = g_ktfor.thr_busy_min;
            ret->prof.thr_busy_max   = g_ktfor.thr_busy_max;
            ret->prof.thr_busy_mean  = g_ktfor.thr_busy_mean;
            ret->prof.thr_busy_stdev = g_ktfor.thr_busy_stdev;
            /* encode = SAM/BAM-build CPU (accurate, summed over compute threads);
             * compute = the rest of the alignment CPU. Same clock, so subtractable. */
            ret->prof.encode  = g_ktfor.encode;
            ret->prof.compute = (g_ktfor.proc_cpu > g_ktfor.encode)
                                ? g_ktfor.proc_cpu - g_ktfor.encode : NAN;
        }

        aux->n_processed += ret->n_seqs;
        return ret;
    }
    /* Step 3: Write output */
    else if (step == 2)
    {
        uint64_t tim = __rdtsc();
        double sp_w0 = sp_enabled() ? sp_wall() : 0.0;
        long sp_wbytes = 0;

        for (int i = 0; i < ret->n_seqs; )
        {
            int group_size = 1;
            if (aux->opt->meth_mode && (aux->opt->flag & MEM_F_PE)
                && i + 1 < ret->n_seqs
                && strcmp(ret->seqs[i].name, ret->seqs[i+1].name) == 0) {
                group_size = 2;
            }

            if (aux->opt->meth_mode && g_meth_bam_writer != NULL) {
                /* Gather all bam1_t* in the QNAME group, propagate QC fail, emit. */
                int total = 0;
                for (int k = 0; k < group_size; ++k) total += ret->seqs[i+k].n_bams;
                if (total > 0) {
                    struct bam1_t **group = (struct bam1_t **)malloc(total * sizeof(struct bam1_t *));
                    if (group == NULL)
                        err_fatal(__func__, "out of memory gathering meth BAM group of %d", total);
                    int idx = 0;
                    for (int k = 0; k < group_size; ++k) {
                        for (int j = 0; j < ret->seqs[i+k].n_bams; ++j) {
                            group[idx++] = (struct bam1_t *)ret->seqs[i+k].bams[j];
                        }
                    }
                    meth_bam_group_propagate_qcfail(group, total);
                    for (int j = 0; j < total; ++j) {
#ifndef DISABLE_OUTPUT
                        if (meth_bam_writer_write(g_meth_bam_writer, group[j]) < 0)
                            err_fatal(__func__, "failed to write meth BAM record");
#endif
                        bam_writer_free(group[j]);
                    }
                    free(group);
                }
            } else if (aux->bam_writer != NULL) {
                for (int k = 0; k < group_size; ++k) {
                    for (int j = 0; j < ret->seqs[i+k].n_bams; ++j) {
#ifndef DISABLE_OUTPUT
                        if (bam_writer_write(aux->bam_writer, (struct bam1_t *)ret->seqs[i+k].bams[j]) < 0)
                            err_fatal(__func__, "failed to write BAM record");
#endif
                        bam_writer_free((struct bam1_t *)ret->seqs[i+k].bams[j]);
                    }
                }
            } else {
                for (int k = 0; k < group_size; ++k) {
                    if (ret->seqs[i+k].sam) {
#ifndef DISABLE_OUTPUT
                        fputs(ret->seqs[i+k].sam, aux->fp);
#endif
                    }
                    /* meth_mode populates bams[] regardless of writer state;
                     * under DISABLE_OUTPUT both writers are forced NULL and
                     * we land here, so the per-record bam1_t allocations
                     * would otherwise leak (only the pointer array below is
                     * freed). Profile-build only, but it skews the very
                     * Maximum RSS that bench/run.sh records. */
                    for (int j = 0; j < ret->seqs[i+k].n_bams; ++j) {
                        bam_writer_free((struct bam1_t *)ret->seqs[i+k].bams[j]);
                    }
                }
            }

            for (int k = 0; k < group_size; ++k) {
                if (sp_enabled() && ret->seqs[i+k].sam) sp_wbytes += (long)strlen(ret->seqs[i+k].sam);
                free(ret->seqs[i+k].name);
                free(ret->seqs[i+k].comment);
                free(ret->seqs[i+k].seq);
                free(ret->seqs[i+k].qual);
                free(ret->seqs[i+k].sam);
                free(ret->seqs[i+k].bams);
                free(ret->seqs[i+k].meth_orig_seq); /* NULL outside --meth; free() is NULL-safe */
            }
            i += group_size;
        }
        if (sp_enabled()) {
            ret->prof.write_wall = sp_wall() - sp_w0;
            ret->prof.write_bytes = sp_wbytes;
            if (aux->opt->bam_mode) {            /* htslib fuses compress+diskwrite */
                ret->prof.write_compress = ret->prof.write_wall;   /* diskwrite stays NaN */
            } else {
                ret->prof.write_diskwrite = ret->prof.write_wall;
                ret->prof.write_compress = 0.0;
            }
            static long g_sp_chunk = 0;
            ret->prof.chunk = __sync_fetch_and_add(&g_sp_chunk, 1);
            sp_add_chunk(&ret->prof);
        }
        free(ret->seqs);
        free(ret);
        tprof[SAM_IO][0] += __rdtsc() - tim;

        return 0;
    } // step 2

    return 0;
}

static void *ktp_worker(void *data)
{
    ktp_worker_t *w = (ktp_worker_t*) data;
    ktp_t *p = w->pl;

    while (w->step < p->n_steps) {
        // test whether we can kick off the job with this worker
        double sp_i0 = sp_enabled() ? sp_wall() : 0.0;   // idle = time waiting for our turn
        int pthread_ret = pthread_mutex_lock(&p->mutex);
        assert(pthread_ret == 0);
        for (;;) {
            int i;
            // test whether another worker is doing the same step
            for (i = 0; i < p->n_workers; ++i) {
                if (w == &p->workers[i]) continue; // ignore itself
                if (p->workers[i].step <= w->step && p->workers[i].index < w->index)
                    break;
            }
            if (i == p->n_workers) break; // no workers with smaller indices are doing w->step or the previous steps
            pthread_ret = pthread_cond_wait(&p->cv, &p->mutex);
            assert(pthread_ret == 0);
        }
        pthread_ret = pthread_mutex_unlock(&p->mutex);
        assert(pthread_ret == 0);
        if (sp_enabled()) sp_add_idle(w->step, sp_wall() - sp_i0);   /* idle attributed to next step */

        // working on w->step
        w->data = kt_pipeline(p->shared, w->step, w->step? w->data : 0, w->opt, *(w->w)); // for the first step, input is NULL

        // update step and let other workers know
        pthread_ret = pthread_mutex_lock(&p->mutex);
        assert(pthread_ret == 0);
        w->step = w->step == p->n_steps - 1 || w->data? (w->step + 1) % p->n_steps : p->n_steps;

        if (w->step == 0) w->index = p->index++;
        pthread_ret = pthread_cond_broadcast(&p->cv);
        assert(pthread_ret == 0);
        pthread_ret = pthread_mutex_unlock(&p->mutex);
        assert(pthread_ret == 0);
    }
    pthread_exit(0);
}

static int process(void *shared, gzFile gfp, gzFile gfp2, int pipe_threads)
{
    ktp_aux_t   *aux = (ktp_aux_t*) shared;
    worker_t     w;
    mem_opt_t   *opt = aux->opt;
    int32_t nthreads = opt->n_threads; // global variable for profiling!
    w.nthreads = opt->n_threads;

#if NUMA_ENABLED
    int  deno = 1;
    int tc = numa_num_task_cpus();
    int tn = numa_num_task_nodes();
    int tcc = numa_num_configured_cpus();
    fprintf(stderr, "num_cpus: %d, num_numas: %d, configured cpus: %d\n", tc, tn, tcc);
    int ht = HTStatus();
    if (ht) deno = 2;

    if (nthreads < tcc/tn/deno) {
        fprintf(stderr, "Enabling single numa domain...\n\n");
        // numa_set_preferred(0);
        // bitmask mask(0);
        struct bitmask *mask = numa_bitmask_alloc(numa_num_possible_nodes());
        numa_bitmask_clearall(mask);
        numa_bitmask_setbit(mask, 0);
        numa_bind(mask);
        numa_bitmask_free(mask);
    }
#else
    /* Report platform info on non-NUMA systems (e.g., macOS/Apple Silicon) */
#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
    HTStatus();
#endif
#endif
#if AFF && (__linux__)
    { // Affinity/HT stuff
        unsigned int cpuid[4];
        asm volatile
            ("cpuid" : "=a" (cpuid[0]), "=b" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
             : "0" (0xB), "2" (1));
        int num_logical_cpus = cpuid[1] & 0xFFFF;

        asm volatile
            ("cpuid" : "=a" (cpuid[0]), "=b" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
             : "0" (0xB), "2" (0));
        int num_ht = cpuid[1] & 0xFFFF;
        int num_total_logical_cpus = get_nprocs_conf();
        int num_sockets = num_total_logical_cpus / num_logical_cpus;
        fprintf(stderr, "#sockets: %d, #cores/socket: %d, #logical_cpus: %d, #ht/core: %d\n",
                num_sockets, num_logical_cpus/num_ht, num_total_logical_cpus, num_ht);

        for (int i=0; i<num_total_logical_cpus; i++) affy[i] = i;
        int slookup[256] = {-1};

        if (num_ht == 2 && num_sockets == 2)  // generalize it for n sockets
        {
            for (int i=0; i<num_total_logical_cpus; i++) {
                std::ostringstream ss;
                ss << i;
                std::string str = "/sys/devices/system/cpu/cpu"+ ss.str();
                str = str +"/topology/thread_siblings_list";
                // std::cout << str << std::endl;
                // std::string str = "cpu.txt";
                FILE *fp = fopen(str.c_str(), "r");
                if (fp == NULL) {
                    fprintf("Error: Cant open the file..\n");
                    break;
                }
                else {
                    int a, b, v;
                    char ch[10] = {'\0'};
                    fgets(ch, 10, fp);
                    v = sscanf(ch, "%u,%u",&a,&b);
                    if (v == 1) v = sscanf(ch, "%u-%u",&a,&b);
                    if (v == 1) {
                        fprintf(stderr, "Mis-match between HT and threads_sibling_list...%s\n", ch);
                        fprintf(stderr, "Continuing with default affinity settings..\n");
                        break;
                    }
                    slookup[a] = 1;
                    slookup[b] = 2;
                    fclose(fp);
                }
            }
            int a = 0, b = num_total_logical_cpus / num_ht;
            for (int i=0; i<num_total_logical_cpus; i++) {
                if (slookup[i] == -1) {
                    fprintf(stderr, "Unseen cpu topology..\n");
                    break;
                }
                if (slookup[i] == 1) affy[a++] = i;
                else affy[b++] = i;
            }
        }
    }
#endif

    int32_t nreads = aux->actual_chunk_size / NREADS_ESTIMATE_AVG_BASES + 10;

    /* All memory allocation */
    memoryAlloc(aux, w, nreads, nthreads);
    fprintf(stderr, "* Threads used (compute): %d\n", nthreads);

    /* pipeline using pthreads */
    ktp_t aux_;
    int p_nt = pipe_threads; // 2;
    int n_steps = 3;

    w.ref_string = aux->ref_string;
    // Mirror into mem_cache so helpers that take only `mmc` (e.g.
    // mem_matesw_batch_pre/post) can call bns_fetch_seq_v2 without
    // threading ref_string through every signature on the way down.
    w.mmc.ref_string = aux->ref_string;
    w.fmi = aux->fmi;
    /* D3 (--meth, PR-3): hand the ORIGINAL bns/pac/.0123 to the worker so the
     * (remapped, original-coord) seeds chain/extend/pair/output against the
     * original reference. NULL outside --meth → mem_aln_* fall back to the seed
     * index and the non-meth path is unchanged. The batched mate-rescue helpers
     * read the ref via w.mmc.ref_string, so point THAT at the original .0123 too
     * in --meth (else mate rescue would fetch SEED bases at original coords). */
    w.meth_orig_bns        = aux->meth_orig_bns;
    w.meth_orig_pac        = aux->meth_orig_pac;
    w.meth_orig_ref_string = aux->meth_orig_ref_string;
    if (aux->opt->meth_mode && aux->meth_orig_ref_string != NULL)
        w.mmc.ref_string = aux->meth_orig_ref_string;
    w.nreads  = nreads;
    // w.memSize = nreads;

    aux_.n_workers = p_nt;
    aux_.n_steps = n_steps;
    aux_.shared = aux;
    aux_.index = 0;
    int pthread_ret = pthread_mutex_init(&aux_.mutex, 0);
    assert(pthread_ret == 0);
    pthread_ret = pthread_cond_init(&aux_.cv, 0);
    assert(pthread_ret == 0);

    fprintf(stderr, "* No. of pipeline threads: %d\n\n", p_nt);
    aux_.workers = (ktp_worker_t*) malloc(p_nt * sizeof(ktp_worker_t));
    assert(aux_.workers != NULL);

    for (int i = 0; i < p_nt; ++i) {
        ktp_worker_t *wr = &aux_.workers[i];
        wr->step = 0; wr->pl = &aux_; wr->data = 0;
        wr->index = aux_.index++;
        wr->i = i;
        wr->opt = opt;
        wr->w = &w;
    }

    pthread_t *ptid = (pthread_t *) calloc(p_nt, sizeof(pthread_t));
    assert(ptid != NULL);

    for (int i = 0; i < p_nt; ++i)
        pthread_create(&ptid[i], 0, ktp_worker, (void*) &aux_.workers[i]);

    for (int i = 0; i < p_nt; ++i)
        pthread_join(ptid[i], 0);

    pthread_ret = pthread_mutex_destroy(&aux_.mutex);
    assert(pthread_ret == 0);
    pthread_ret = pthread_cond_destroy(&aux_.cv);
    assert(pthread_ret == 0);

    free(ptid);
    free(aux_.workers);
    /***** pipeline ends ******/

    /* Retire the kt_for() worker pool (created lazily on the first chunk).
     * No kt_for() calls remain past this point. */
    kt_pool_destroy();

    fprintf(stderr, "[0000] Computation ends..\n");

    /* Dealloc per-worker scratch buffers allocated in the header section */
    worker_free(w, nthreads);

    return 0;
}

static void update_a(mem_opt_t *opt, const mem_opt_t *opt0)
{
    if (opt0->a) { // matching score is changed
        if (!opt0->b) opt->b *= opt->a;
        if (!opt0->T) opt->T *= opt->a;
        if (!opt0->o_del) opt->o_del *= opt->a;
        if (!opt0->e_del) opt->e_del *= opt->a;
        if (!opt0->o_ins) opt->o_ins *= opt->a;
        if (!opt0->e_ins) opt->e_ins *= opt->a;
        if (!opt0->zdrop) opt->zdrop *= opt->a;
        if (!opt0->pen_clip5) opt->pen_clip5 *= opt->a;
        if (!opt0->pen_clip3) opt->pen_clip3 *= opt->a;
        if (!opt0->pen_unpaired) opt->pen_unpaired *= opt->a;
    }
}

static void usage(const mem_opt_t *opt)
{
    fprintf(stderr, "Usage: bwa-mem3 mem [options] <idxbase> <in1.fq> [in2.fq]\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  Algorithm options:\n");
    fprintf(stderr, "    -o STR        Output SAM file name\n");
    fprintf(stderr, "    --bam[=N]     Emit BAM instead of SAM text. N=0 (default) = uncompressed;\n");
    fprintf(stderr, "                  1..9 = BGZF deflate levels. Writes to stdout; redirect with `>`.\n");
    fprintf(stderr, "    -t INT        number of threads [%d]\n", opt->n_threads);
    fprintf(stderr, "    -k INT        minimum seed length [%d]\n", opt->min_seed_len);
    fprintf(stderr, "    -w INT        band width for banded alignment [%d]\n", opt->w);
    fprintf(stderr, "    -d INT        off-diagonal X-dropoff [%d]\n", opt->zdrop);
    fprintf(stderr, "    -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [%g]\n", opt->split_factor);
    fprintf(stderr, "    -y INT        seed occurrence for the 3rd round seeding [%ld]\n", (long)opt->max_mem_intv);
    fprintf(stderr, "    -c INT        skip seeds with more than INT occurrences [%d]\n", opt->max_occ);
    fprintf(stderr, "    --smem-dedup  dedup identical SMEMs before chaining: fewer SA lookups, ~10%% fewer; opt-in, NOT byte-identical (changes XS/secondary on a small fraction of reads) [off]\n");
    fprintf(stderr, "    --skip-contained-ext  skip banded-SW extension of seeds contained (same diagonal) in a longer in-chain seed; byte-identical (non-meth), ~10%% less alignment CPU; no effect under --meth [off]\n");
    fprintf(stderr, "    --max-extend-chains INT  cap chains extended per read to the top-INT by weight; ~23%% less alignment CPU, high-confidence placement unaffected; ignored for reads with >4096 chains; opt-in, NOT byte-identical (0 = off) [%d]\n", opt->max_extend_chains);
    fprintf(stderr, "    --adaptive-band  adaptive banded-SW: start tight and expand each pair to its chain-geometry band on long-extension reads; long-read speedup (~1.3x on SBX), no-op on short reads; opt-in, NOT byte-identical [%s]\n", opt->band_start? "on":"off");
    fprintf(stderr, "    --extend-mate-concordant[=INT]  when --max-extend-chains caps a PE read, also keep any chain concordant (same contig, FR, within INT bp) with a mate chain; recovers the true pair's low-weight chain the cap would drop (mainly --meth). Bare = auto (window = estimated proper-pair insert high bound); =INT = fixed bp; =0 = off. Opt-in, NOT byte-identical [%s]\n", opt->mate_concordant_window? (opt->mate_concordant_window<0? "auto":"fixed") : "off");
    fprintf(stderr, "    -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [%.2f]\n", opt->drop_ratio);
    fprintf(stderr, "    -W INT        discard a chain if seeded bases shorter than INT [0]\n");
    fprintf(stderr, "    -m INT        perform at most INT rounds of mate rescues for each read [%d]\n", opt->max_matesw);
    fprintf(stderr, "    -S            skip mate rescue\n");
    fprintf(stderr, "    -P            skip pairing; mate rescue performed unless -S also in use\n");
    fprintf(stderr, "    --fast        speed preset: -m 10 -y 0 --min-ext-len 30 --smem-dedup\n");
    fprintf(stderr, "                  --skip-contained-ext --max-extend-chains 5 --adaptive-band (and\n");
    fprintf(stderr, "                  -s 2 --extend-mate-concordant under --meth). Opt-in; explicit\n");
    fprintf(stderr, "                  flags override where applicable; --smem-dedup,\n");
    fprintf(stderr, "                  --skip-contained-ext and --adaptive-band are always enabled.\n");
    fprintf(stderr, "                  NOT byte-identical to the default (divergence confined to the\n");
    fprintf(stderr, "                  low-confidence tail).\n");
    fprintf(stderr, "Scoring options:\n");
    fprintf(stderr, "   -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [%d]\n", opt->a);
    fprintf(stderr, "   -B INT        penalty for a mismatch [%d]\n", opt->b);
    fprintf(stderr, "   -O INT[,INT]  gap open penalties for deletions and insertions [%d,%d]\n", opt->o_del, opt->o_ins);
    fprintf(stderr, "   -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [%d,%d]\n", opt->e_del, opt->e_ins);
    fprintf(stderr, "   -L INT[,INT]  penalty for 5'- and 3'-end clipping [%d,%d]\n", opt->pen_clip5, opt->pen_clip3);
    fprintf(stderr, "   -U INT        penalty for an unpaired read pair [%d]\n", opt->pen_unpaired);
//  fprintf(stderr, "   -x STR        read type. Setting -x changes multiple parameters unless overriden [null]\n");
//  fprintf(stderr, "                 pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)\n");
//  fprintf(stderr, "                 ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)\n");
//  fprintf(stderr, "                 intractg: -B9 -O16 -L5  (intra-species contigs to ref)\n");
    fprintf(stderr, "Input/output options:\n");
    fprintf(stderr, "   -p            smart pairing (ignoring in2.fq)\n");
    fprintf(stderr, "   -R STR        read group header line such as '@RG\\tID:foo\\tSM:bar' [null]\n");
    fprintf(stderr, "   -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]\n");
    fprintf(stderr, "   -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)\n");
    fprintf(stderr, "   -5            for split alignment, take the alignment with the smallest coordinate as primary\n");
    fprintf(stderr, "   -q            don't modify mapQ of supplementary alignments\n");
    fprintf(stderr, "   -K INT        process INT input bases in each batch regardless of nThreads (for reproducibility) []\n");
    fprintf(stderr, "   -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
    fprintf(stderr, "   -T INT        minimum score to output [%d]\n", opt->T);
    fprintf(stderr, "   -h INT[,INT]  if there are <INT hits with score >%.2f%% of the max score, output all in XA [%d,%d]\n",
            opt->XA_drop_ratio * 100.0, opt->max_XA_hits, opt->max_XA_hits_alt);
    fprintf(stderr, "   -z FLOAT      the fraction of the max score to use with -h [%.2f]\n", opt->XA_drop_ratio);
    fprintf(stderr, "   -u            output XB instead of XA; XB is XA with the alignment score and mapping quality added\n");
    fprintf(stderr, "   -a            output all alignments for SE or unpaired PE\n");
    fprintf(stderr, "   -C            append FASTA/FASTQ comment to SAM output\n");
    fprintf(stderr, "   -V            output the reference FASTA header in the XR tag\n");
    fprintf(stderr, "   -Y            use soft clipping for supplementary alignments\n");
    fprintf(stderr, "   -M            mark shorter split hits as secondary\n");
    fprintf(stderr, "   -I FLOAT[,FLOAT[,INT[,INT]]]\n");
    fprintf(stderr, "                 specify the mean, standard deviation (10%% of the mean if absent), max\n");
    fprintf(stderr, "                 (4 sigma from the mean if absent) and min of the insert size distribution.\n");
    fprintf(stderr, "                 FR orientation only. [inferred]\n");
    fprintf(stderr, "Bisulfite (--meth) options:\n");
    fprintf(stderr, "   --meth        enable inline bwameth-style C→T/G→A read conversion + meth-aware BAM\n");
    fprintf(stderr, "                 emission. Implies --bam. Requires the reference to have been built\n");
    fprintf(stderr, "                 with `bwa-mem3 index --meth` (emits the original index plus a\n");
    fprintf(stderr, "                 ref.fa.meth.* converted seed index).\n");
    fprintf(stderr, "   --meth-scoring collapsed|genomic\n");
    fprintf(stderr, "                 bisulfite scoring mode [collapsed]. collapsed: C/T (and G/A)\n");
    fprintf(stderr, "                 interchangeable, bwameth-compatible placement (sets -B 2).\n");
    fprintf(stderr, "                 genomic: free only the conversion direction, keep variants as\n");
    fprintf(stderr, "                 mismatches (variant-aware, truthful NM/MD; -B 4).\n");
    fprintf(stderr, "   --set-as-failed f|r\n");
    fprintf(stderr, "                 flag alignments to the matching strand ('f' or 'r') as QC-fail (0x200)\n");
    fprintf(stderr, "   --chimera-qc\n");
    fprintf(stderr, "                 enable the bwameth.py-style longest-match <44%% chimera heuristic\n");
    fprintf(stderr, "                 (sets 0x200, clears 0x2, caps MAPQ at 1). Off by default; not in Bismark.\n");
    fprintf(stderr, "Supplementary MAPQ rescoring (fg-labs extension):\n");
    fprintf(stderr, "   --supp-rep-hard-cap INT\n");
    fprintf(stderr, "                 force MAPQ=0 for supplementary alignments whose chain contains any seed\n");
    fprintf(stderr, "                 with >=INT genome occurrences (i.e. the supp region is repetitive on its\n");
    fprintf(stderr, "                 own). 0 disables (default). Typical values 5-20; lower = more aggressive.\n");
    fprintf(stderr, "                 Primary MAPQ is unaffected.\n");
    fprintf(stderr, "Seed ordering (fg-labs extension):\n");
    fprintf(stderr, "   --seed-order STR\n");
    fprintf(stderr, "                 seed emission order before chaining: off|local-longest [off]\n");
    fprintf(stderr, "                 (advanced modes: global-longest, absorb-count, most-absorb; see docs)\n");
    fprintf(stderr, "Input reader:\n");
    fprintf(stderr, "   --legacy-reader\n");
    fprintf(stderr, "                 use the legacy gzFile/kseq input reader instead of the default\n");
    fprintf(stderr, "                 content-detecting fast reader (escape hatch / A-B baseline).\n");
#ifdef STAGE_PROF
    fprintf(stderr, "Profiling:\n");
    fprintf(stderr, "   --profile FILE\n");
    fprintf(stderr, "                 write per-chunk stage profiling TSV to FILE (off by default).\n");
#endif
    fprintf(stderr, "Help:\n");
    fprintf(stderr, "   --help        print this help message and exit\n");
    fprintf(stderr, "Note: Please read the man page for detailed description of the command line and options.\n");
}

/* D3 (--meth) only: load the ORIGINAL reference's bns + pac (un-converted, real
 * chrom names) from `prefix` as resident handles for the future extension/scoring
 * phase, distinct from the seed FM-index. Mirrors indexEle::bwa_idx_load_ele's
 * disk path (bns_restore then slurp the full .pac into memory and close fp_pac).
 * On success writes *bns_out / *pac_out and returns 0; on failure frees any
 * partial allocation, leaves the out-params NULL, and returns -1. */
static int meth_orig_ref_load_handles(const char *prefix,
                                      bntseq_t **bns_out, uint8_t **pac_out)
{
    *bns_out = NULL;
    *pac_out = NULL;
    bntseq_t *bns = bns_restore(prefix);
    if (bns == NULL) {
        fprintf(stderr,
                "ERROR: --meth could not load the original reference bns from "
                "'%s.{amb,ann,pac}'\n", prefix);
        return -1;
    }
    int64_t pac_bytes = bns->l_pac / 4 + 1;
    uint8_t *pac = (uint8_t*) calloc(pac_bytes, 1);
    if (pac == NULL) {
        fprintf(stderr, "ERROR: --meth failed to allocate %lld bytes for the "
                "original reference pac\n", (long long)pac_bytes);
        bns_destroy(bns);
        return -1;
    }
    bwamem_madv_hugepage(pac, pac_bytes);
    /* bns_restore left .pac open in bns->fp_pac; slurp it whole, then close. */
    err_fread_noeof(pac, 1, pac_bytes, bns->fp_pac);
    err_fclose(bns->fp_pac);
    bns->fp_pac = NULL;

    *bns_out = bns;
    *pac_out = pac;
    return 0;
}

/* Free the original-reference handles loaded by meth_orig_ref_load_handles and
 * NULL them. Idempotent (NULL-safe) so it can sit on every exit path the seed
 * index is freed on without double-free risk. */
static void meth_orig_ref_free_handles(ktp_aux_t *aux)
{
    if (aux->meth_orig_pac != NULL) { free(aux->meth_orig_pac); aux->meth_orig_pac = NULL; }
    if (aux->meth_orig_bns != NULL) { bns_destroy(aux->meth_orig_bns); aux->meth_orig_bns = NULL; }
    /* D3 (--meth): meth_orig_ref_string is NULL under pac-fetch (the only path
     * now — the original reference is unpacked from `.pac` on demand, never
     * materialized). This NULL-safe free is retained defensively; if it is ever
     * non-NULL it would be an _mm_malloc'd buffer, freed with _mm_free. */
    if (aux->meth_orig_ref_string != NULL) {
        _mm_free(aux->meth_orig_ref_string);
        aux->meth_orig_ref_string = NULL;
    }
}

int main_mem(int argc, char *argv[])
{
    int          i, c, ignore_alt = 0, no_mt_io = 0;
    int          fixed_chunk_size          = -1;
    char        *p, *rg_line               = 0, *hdr_line = 0;
    const char  *mode                      = 0;
    int          fast                      = 0;

    mem_opt_t    *opt, opt0;
    gzFile        fp = 0, fp2 = 0;
    void         *ko = 0, *ko2 = 0;
    int           fd, fd2;
    mem_pestat_t  pes[4];
    ktp_aux_t     aux;
    bool          is_o       = 0;
    const char   *out_path   = NULL;   /* -o/-f path; opened lazily so --bam can claim it */
    bool          out_opened = false;  /* true iff aux.fp is a real fopen()'d FILE* we own */
    uint8_t      *ref_string;

    memset(&aux, 0, sizeof(ktp_aux_t));
    memset(pes, 0, 4 * sizeof(mem_pestat_t));
    for (i = 0; i < 4; ++i) pes[i].failed = 1;

    // opterr = 0;
    aux.fp = stdout;
    aux.opt = opt = mem_opt_init();
    memset(&opt0, 0, sizeof(mem_opt_t));

    /* Parse input arguments */
    // comment: added option '5' in the list
    //
    // Long-only options for bisulfite mode (bwa-mem3 meth fork):
    //   --meth              Enable inline bwameth-style c2t + post-processing + BAM output.
    //                       Expects a reference built with `bwa-mem3 index --meth`.
    //   --set-as-failed f|r Flag alignments to this strand as QC-fail (0x200)
    //   --chimera-qc        Enable the bwameth.py-style longest-M <44% chimera heuristic
    //                       (off by default; not part of Bismark)
    enum {
        OPT_BAM = 1000,
        OPT_METH,
        OPT_METH_SCORING,
        OPT_METH_SET_AS_FAILED,
        OPT_METH_CHIMERA_QC,
        OPT_SUPP_REP_HARD_CAP,
        OPT_LEGACY_READER,
        OPT_MIN_EXT_LEN,
        OPT_MAX_EXTEND_CHAINS,
        OPT_SEED_ORDER,
        OPT_SMEM_DEDUP,
        OPT_FAST,
        OPT_SKIP_CONTAINED_EXT,
        OPT_ADAPTIVE_BAND,
        OPT_EXTEND_MATE_CONCORDANT,
#ifdef STAGE_PROF
        OPT_PROFILE,
#endif
        OPT_HELP,
    };
    static struct option long_opts[] = {
        {"bam",                      optional_argument, 0, OPT_BAM},
        {"min-ext-len",              required_argument, 0, OPT_MIN_EXT_LEN},
        {"max-extend-chains",        required_argument, 0, OPT_MAX_EXTEND_CHAINS},
        {"smem-dedup",               no_argument,       0, OPT_SMEM_DEDUP},
        {"fast",                     no_argument,       0, OPT_FAST},
        {"skip-contained-ext",       no_argument,       0, OPT_SKIP_CONTAINED_EXT},
        {"adaptive-band",            no_argument,       0, OPT_ADAPTIVE_BAND},
        {"extend-mate-concordant",   optional_argument, 0, OPT_EXTEND_MATE_CONCORDANT},
        {"meth",                     no_argument,       0, OPT_METH},
        {"meth-scoring",             required_argument, 0, OPT_METH_SCORING},
        {"set-as-failed",            required_argument, 0, OPT_METH_SET_AS_FAILED},
        {"chimera-qc",               no_argument,       0, OPT_METH_CHIMERA_QC},
        {"supp-rep-hard-cap",        required_argument, 0, OPT_SUPP_REP_HARD_CAP},
        {"seed-order",               required_argument, 0, OPT_SEED_ORDER},
        {"legacy-reader",            no_argument,       0, OPT_LEGACY_READER},
#ifdef STAGE_PROF
        {"profile",                  required_argument, 0, OPT_PROFILE},
#endif
        {"help",                     no_argument,       0, OPT_HELP},
        {0, 0, 0, 0}
    };
#ifdef STAGE_PROF
    const char *profile_path = NULL;   /* --profile <path>: stage_prof TSV output */
#endif
    while ((c = getopt_long(argc, argv, "51qpaMCSPVYjuk:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:W:x:G:h:y:K:X:H:o:f:z:",
                            long_opts, NULL)) >= 0)
    {
        if (c == 'k') opt->min_seed_len = atoi(optarg), opt0.min_seed_len = 1;
        else if (c == OPT_MIN_EXT_LEN) opt->min_ext_len = atoi(optarg), opt0.min_ext_len = 1;
        else if (c == OPT_MAX_EXTEND_CHAINS) opt->max_extend_chains = atoi(optarg), opt0.max_extend_chains = 1;
        else if (c == '1') no_mt_io = 1;
        else if (c == 'x') mode = optarg;
        else if (c == 'w') opt->w = atoi(optarg), opt0.w = 1;
        else if (c == 'A') opt->a = atoi(optarg), opt0.a = 1, assert(opt->a >= INT_MIN && opt->a <= INT_MAX);
        else if (c == 'B') opt->b = atoi(optarg), opt0.b = 1, assert(opt->b >= INT_MIN && opt->b <= INT_MAX);
        else if (c == 'T') opt->T = atoi(optarg), opt0.T = 1, assert(opt->T >= INT_MIN && opt->T <= INT_MAX);
        else if (c == 'U')
            opt->pen_unpaired = atoi(optarg), opt0.pen_unpaired = 1, assert(opt->pen_unpaired >= INT_MIN && opt->pen_unpaired <= INT_MAX);
        else if (c == 't')
            opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1, assert(opt->n_threads >= INT_MIN && opt->n_threads <= INT_MAX);
        else if (c == 'o' || c == 'f')
        {
            /* Capture the path; defer opening until after --bam is parsed so
             * BAM mode can hand the path to htslib instead of truncating it
             * here as a SAM-text FILE*. */
            is_o = 1;
            out_path = optarg;
        }
        else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;
        else if (c == 'a') opt->flag |= MEM_F_ALL;
        else if (c == 'p') opt->flag |= MEM_F_PE | MEM_F_SMARTPE;
        else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;
        else if (c == 'S') opt->flag |= MEM_F_NO_RESCUE;
        else if (c == 'Y') opt->flag |= MEM_F_SOFTCLIP;
        else if (c == 'V') opt->flag |= MEM_F_REF_HDR;
        else if (c == '5') opt->flag |= MEM_F_PRIMARY5 | MEM_F_KEEP_SUPP_MAPQ; // always apply MEM_F_KEEP_SUPP_MAPQ with -5
        else if (c == 'q') opt->flag |= MEM_F_KEEP_SUPP_MAPQ;
        else if (c == 'u') opt->flag |= MEM_F_XB;
        else if (c == 'c') opt->max_occ = atoi(optarg), opt0.max_occ = 1;
        else if (c == 'd') opt->zdrop = atoi(optarg), opt0.zdrop = 1;
        else if (c == 'v') bwa_verbose = atoi(optarg);
        else if (c == 'j') ignore_alt = 1;
        else if (c == 'r')
            opt->split_factor = atof(optarg), opt0.split_factor = 1.;
        else if (c == 'D') opt->drop_ratio = atof(optarg), opt0.drop_ratio = 1.;
        else if (c == 'm') opt->max_matesw = atoi(optarg), opt0.max_matesw = 1;
        else if (c == 's') opt->split_width = atoi(optarg), opt0.split_width = 1;
        else if (c == 'G')
            opt->max_chain_gap = atoi(optarg), opt0.max_chain_gap = 1;
        else if (c == 'N')
            opt->max_chain_extend = atoi(optarg), opt0.max_chain_extend = 1;
        else if (c == 'W')
            opt->min_chain_weight = atoi(optarg), opt0.min_chain_weight = 1;
        else if (c == 'y')
            opt->max_mem_intv = atol(optarg), opt0.max_mem_intv = 1;
        else if (c == 'C') aux.copy_comment = 1;
        else if (c == 'K') fixed_chunk_size = atoi(optarg);
        else if (c == 'X') opt->mask_level = atof(optarg);
        else if (c == 'h')
        {
            opt0.max_XA_hits = opt0.max_XA_hits_alt = 1;
            opt->max_XA_hits = opt->max_XA_hits_alt = strtol(optarg, &p, 10);
            if (*p != 0 && ispunct(*p) && isdigit(p[1]))
                opt->max_XA_hits_alt = strtol(p+1, &p, 10);
        }
        else if (c == 'z') opt->XA_drop_ratio = atof(optarg);
        else if (c == 'Q')
        {
            opt0.mapQ_coef_len = 1;
            opt->mapQ_coef_len = atoi(optarg);
            opt->mapQ_coef_fac = opt->mapQ_coef_len > 0? log(opt->mapQ_coef_len) : 0;
        }
        else if (c == 'O')
        {
            opt0.o_del = opt0.o_ins = 1;
            opt->o_del = opt->o_ins = strtol(optarg, &p, 10);
            if (*p != 0 && ispunct(*p) && isdigit(p[1]))
                opt->o_ins = strtol(p+1, &p, 10);
        }
        else if (c == 'E')
        {
            opt0.e_del = opt0.e_ins = 1;
            opt->e_del = opt->e_ins = strtol(optarg, &p, 10);
            if (*p != 0 && ispunct(*p) && isdigit(p[1]))
                opt->e_ins = strtol(p+1, &p, 10);
        }
        else if (c == 'L')
        {
            opt0.pen_clip5 = opt0.pen_clip3 = 1;
            opt->pen_clip5 = opt->pen_clip3 = strtol(optarg, &p, 10);
            if (*p != 0 && ispunct(*p) && isdigit(p[1]))
                opt->pen_clip3 = strtol(p+1, &p, 10);
        }
        else if (c == 'R')
        {
            if ((rg_line = bwa_set_rg(optarg)) == 0) {
                free(opt);
                if (out_opened)
                    fclose(aux.fp);
                return 1;
            }
        }
        else if (c == 'H')
        {
            if (optarg[0] != '@')
            {
                FILE *fp;
                if ((fp = fopen(optarg, "r")) != 0)
                {
                    hdr_line = bwa_insert_header_file(fp, hdr_line);
                    fclose(fp);
                }
            } else hdr_line = bwa_insert_header(optarg, hdr_line);
        }
        else if (c == OPT_BAM) {
            opt->bam_mode = 1;
            opt->bam_level = 0;
            if (optarg != NULL && optarg[0] != '\0') {
                int lvl = atoi(optarg);
                if (lvl < 0 || lvl > 9) {
                    fprintf(stderr, "ERROR: --bam level must be 0..9 (got '%s')\n", optarg);
                    free(opt);
                    if (out_opened) fclose(aux.fp);
                    return 1;
                }
                opt->bam_level = lvl;
            }
        }
#ifdef STAGE_PROF
        else if (c == OPT_PROFILE) {
            profile_path = optarg;
        }
#endif
        else if (c == OPT_METH) {
            opt->meth_mode = 1;
            opt->bam_mode = 1;  /* meth implies BAM output */
        }
        else if (c == OPT_METH_SCORING) {
            if (optarg != NULL && strcmp(optarg, "collapsed") == 0) {
                opt->meth_scoring = MEM_METH_SCORING_COLLAPSED;
            } else if (optarg != NULL && strcmp(optarg, "genomic") == 0) {
                opt->meth_scoring = MEM_METH_SCORING_GENOMIC;
            } else {
                fprintf(stderr, "ERROR: --meth-scoring requires 'collapsed' or 'genomic'\n");
                free(opt);
                if (out_opened) fclose(aux.fp);
                return 1;
            }
        }
        else if (c == OPT_METH_SET_AS_FAILED) {
            if (optarg == NULL || !(optarg[0] == 'f' || optarg[0] == 'r') || optarg[1] != '\0') {
                fprintf(stderr, "ERROR: --set-as-failed requires 'f' or 'r'\n");
                free(opt);
                if (out_opened) fclose(aux.fp);
                return 1;
            }
            opt->meth_set_as_failed = optarg[0];
        }
        else if (c == OPT_METH_CHIMERA_QC) {
            opt->meth_chimera_qc = 1;
        }
        else if (c == OPT_SUPP_REP_HARD_CAP) {
            char *end = NULL;
            errno = 0;
            long v = strtol(optarg, &end, 10);
            if (end == optarg || end == NULL || *end != '\0' ||
                errno == ERANGE || v < 0 || v > INT_MAX) {
                fprintf(stderr, "ERROR: --supp-rep-hard-cap requires a non-negative integer\n");
                free(opt);
                if (out_opened) fclose(aux.fp);
                return 1;
            }
            opt->supp_rep_hard_cap = (int)v;
        }
        else if (c == OPT_LEGACY_READER) aux.legacy_reader = 1;
        else if (c == OPT_SEED_ORDER) {
            opt->seed_emit_order = seed_order_from_str(optarg);
            if ((int)opt->seed_emit_order < 0) {
                fprintf(stderr, "[E::%s] unknown --seed-order '%s' (off|local-longest; "
                        "see docs for advanced modes)\n", __func__, optarg);
                free(opt);
                if (out_opened) fclose(aux.fp);
                return 1;
            }
        }
        else if (c == OPT_SMEM_DEDUP) opt->smem_dedup = 1;
        else if (c == OPT_FAST) fast = 1;
        else if (c == OPT_SKIP_CONTAINED_EXT) opt->skip_contained_ext = 1;
        else if (c == OPT_ADAPTIVE_BAND) opt->band_start = ADAPTIVE_BAND_START;
        else if (c == OPT_EXTEND_MATE_CONCORDANT) {
            /* bare flag = auto (-1, use the estimated insert-size high bound);
             * =INT = fixed window in bp; =0 = off. */
            opt->mate_concordant_window = optarg ? atoi(optarg) : -1;
            opt0.mate_concordant_window = 1;
        }
        else if (c == OPT_HELP) {
            usage(opt);
            free(opt);
            free(hdr_line);
            free(rg_line);
            if (out_opened) fclose(aux.fp);
            return 0;
        }
        else if (c == 'I')
        {
            aux.pes0 = pes;
            pes[1].failed = 0;
            pes[1].avg = strtod(optarg, &p);
            pes[1].std = pes[1].avg * .1;
            if (*p != 0 && ispunct(*p) && isdigit(p[1]))
                pes[1].std = strtod(p+1, &p);
            pes[1].high = (int)(pes[1].avg + 4. * pes[1].std + .499);
            pes[1].low  = (int)(pes[1].avg - 4. * pes[1].std + .499);
            if (pes[1].low < 1) pes[1].low = 1;
            if (*p != 0 && ispunct(*p) && isdigit(p[1]))
                pes[1].high = (int)(strtod(p+1, &p) + .499);
            if (*p != 0 && ispunct(*p) && isdigit(p[1]))
                pes[1].low  = (int)(strtod(p+1, &p) + .499);
        }
        else {
            free(opt);
            if (out_opened)
                fclose(aux.fp);
            return 1;
        }
    }

    /* Check output file name */
    if (rg_line)
    {
        hdr_line = bwa_insert_header(rg_line, hdr_line);
        free(rg_line);
    }

    if (opt->n_threads < 1) opt->n_threads = 1;
    if (optind + 2 != argc && optind + 3 != argc) {
        usage(opt);
        free(opt);
        if (out_opened)
            fclose(aux.fp);
        return 1;
    }

    /* Further input parsing */
    if (mode)
    {
        fprintf(stderr, "WARNING: bwa-mem3 doesn't work well with long reads or contigs; please use minimap2 instead.\n");
        if (strcmp(mode, "intractg") == 0)
        {
            if (!opt0.o_del) opt->o_del = 16;
            if (!opt0.o_ins) opt->o_ins = 16;
            if (!opt0.b) opt->b = 9;
            if (!opt0.pen_clip5) opt->pen_clip5 = 5;
            if (!opt0.pen_clip3) opt->pen_clip3 = 5;
        }
        else if (strcmp(mode, "pacbio") == 0 || strcmp(mode, "pbref") == 0 || strcmp(mode, "ont2d") == 0)
        {
            if (!opt0.o_del) opt->o_del = 1;
            if (!opt0.e_del) opt->e_del = 1;
            if (!opt0.o_ins) opt->o_ins = 1;
            if (!opt0.e_ins) opt->e_ins = 1;
            if (!opt0.b) opt->b = 1;
            if (opt0.split_factor == 0.) opt->split_factor = 10.;
            if (strcmp(mode, "ont2d") == 0)
            {
                if (!opt0.min_chain_weight) opt->min_chain_weight = 20;
                if (!opt0.min_seed_len) opt->min_seed_len = 14;
                if (!opt0.pen_clip5) opt->pen_clip5 = 0;
                if (!opt0.pen_clip3) opt->pen_clip3 = 0;
            }
            else
            {
                if (!opt0.min_chain_weight) opt->min_chain_weight = 40;
                if (!opt0.min_seed_len) opt->min_seed_len = 17;
                if (!opt0.pen_clip5) opt->pen_clip5 = 0;
                if (!opt0.pen_clip3) opt->pen_clip3 = 0;
            }
        }
        else
        {
            fprintf(stderr, "[E::%s] unknown read type '%s'\n", __func__, mode);
            free(opt);
            if (out_opened)
                fclose(aux.fp);
            return 1;
        }
    } else update_a(opt, &opt0);

    /* --fast: one-flag shorthand for the characterized speed levers
     *   -m 10  -y 0  --min-ext-len 30  --smem-dedup  --skip-contained-ext
     *   --max-extend-chains 5  --adaptive-band
     *   (under --meth: also adds -s 2 and --extend-mate-concordant).
     * Mirrors the -x preset: each lever is applied only when the user did not
     * set it explicitly (opt0), so explicit flags win where applicable. The
     * exceptions are --smem-dedup and --skip-contained-ext, which are plain
     * on/off booleans forced on unconditionally (no opt-out flag exists).
     * --skip-contained-ext is byte-identical on non-meth SE/PE and no-ops under
     * --meth via its own internal gate (see bwamem.cpp), so forcing it on here is
     * safe for --fast --meth too.
     * Output is NOT byte-identical to the default; divergence is confined to the
     * low-confidence tail (see docs/best-practices/settings-profiles.md).
     * meth_mode is already resolved here (parsed in the getopt loop above). */
    if (fast) {
        if (!opt0.max_matesw)   opt->max_matesw   = 10;  /* -m 10 */
        if (!opt0.max_mem_intv) opt->max_mem_intv = 0;   /* -y 0  */
        if (!opt0.min_ext_len)  opt->min_ext_len  = 30;  /* --min-ext-len 30 */
        /* --max-extend-chains: 5 for non-meth; 10 under --meth. A 7-point ablation
         * ({0,5,10,20,50,100,1000}) on 1M sim-meth PE pairs (with mate-concordant
         * rescue on, below) shows chr-accuracy flat (0.9908) at every cap but the
         * confident wrong-chromosome rate is U-shaped, minimized at 10 (cap 5: 592
         * MAPQ>=30 mismaps; cap 10: 382; uncapped: 1056), for +0.7s wall (20.2->20.9s,
         * still -6% vs uncapped). Non-meth keeps 5 (its placement is cap-insensitive
         * and 5 is the pure-speed pick). */
        if (!opt0.max_extend_chains) opt->max_extend_chains = opt->meth_mode ? 10 : 5;
        opt->smem_dedup = 1;                             /* --smem-dedup (plain on/off) */
        opt->skip_contained_ext = 1;                     /* --skip-contained-ext (plain on/off;
                                                          * meth-gated internally) */
        opt->band_start = ADAPTIVE_BAND_START;           /* --adaptive-band: no-op on short reads
                                                          * (8-bit tier untouched), ~25% faster on
                                                          * long-read (SBX/HiFi/ONT) runs. */
        /* --extend-mate-concordant (meth only): the top-5 chain cap regresses
         * bisulfite PE placement. Mechanism (instrumented on 50k sim-meth-place
         * pairs vs truth): NOT chain-dropping -- in 89% of regressions the read's
         * true alignment is still a candidate, but the collapsed 3-letter alphabet
         * flattens chain weights so the read carries many chains, and capping to 5
         * starves PE pairing/mate-rescue of the secondary anchors that let the true
         * concordant pair win; both mates then flip together to a wrong concordant
         * locus (99% proper-pair). Keeping any capped chain that is concordant with
         * a mate chain retains exactly the true pair's low-weight chain while still
         * dropping the far/redundant ones, recovering placement to default parity
         * (97.64% -> 98.08%, == cap-off 98.09%). Non-meth --fast keeps the plain
         * cap (WGS placement is already unaffected and the exemption would erode
         * the speedup). Auto (-1) sizes the concordance window to the estimated
         * proper-pair insert bound so only genuine pair anchors are retained. */
        if (opt->meth_mode && !opt0.mate_concordant_window)
            opt->mate_concordant_window = -1;
        if (opt->meth_mode && !opt0.split_width)
            opt->split_width = 2;                        /* -s 2 (meth only): light Pass-2 reseed.
                                                          * -s 0 (no reseed) inflates MAPQ on bisulfite
                                                          * reads (interior-repeat competitors go unfound);
                                                          * -s 2 reseeds the occurrence-1 SMEMs that inflate,
                                                          * recovering MAPQ+placement at ~the same speed. */
    }

    /* Meth-mode default tuning. bwameth.py runs bwa as
     * `bwa mem -T 40 -B 2 -L 10 -CM`, adding `-U 100 -p` for paired-end. We adopt
     * the soft-clip (-L 10), unpaired (-U 100), output-threshold (-T 40), -M and
     * -C defaults for BOTH --meth-scoring modes; the ONLY mode-dependent knob is
     * the mismatch penalty (the leniency gate):
     *   COLLAPSED (bwameth drop-in): -B 2 — bwameth's lenient mismatch. Combined
     *     with the two-cell matrix (C/T and G/A free both ways) this reproduces
     *     bwameth's collapsed-space placement.
     *   GENOMIC (variant-aware): keep bwa's default -B 4 — the full-hg38 variant
     *     A/B with the asymmetric matrix showed b=4 places better and is better
     *     MAPQ-calibrated than b=2 (placement 92.6 vs 92.5, discordant MAPQ
     *     1.8 vs 2.1).
     * pen_unpaired is only consulted for paired-end rescue, so setting it
     * unconditionally is a no-op for single-end. -A/-B always override and reach
     * the matrices (mem_opt_fill_meth_mat below). */
    if (opt->meth_mode) {
        if (!opt0.pen_clip5)    opt->pen_clip5   = 10;
        if (!opt0.pen_clip3)    opt->pen_clip3   = 10;
        if (!opt0.pen_unpaired) opt->pen_unpaired = 100;  /* bwameth -U 100 (paired) */
        if (!opt0.T)            opt->T           = 40;
        opt->flag |= MEM_F_NO_MULTI;   /* -M */
        aux.copy_comment = 1;          /* -C, needed for YS:Z/YC:Z passthrough */
        if (opt->meth_scoring == MEM_METH_SCORING_COLLAPSED) {
            if (!opt0.b) opt->b = 2;   /* bwameth's lenient mismatch */
        }
        /* GENOMIC keeps bwa's default b=4 (variant-aware). */
    }

    /* Matrix for SWA */
    bwa_fill_scmat(opt->a, opt->b, opt->mat);
    /* D3 (--meth): re-derive the per-hypothesis asymmetric matrices from the matrix
     * we just rebuilt, so -A/-B and the -x presets reach meth scoring (they set
     * opt->a/opt->b above; without this the meth matrices keep init-time defaults). */
    mem_opt_fill_meth_mat(opt);

    /* In --meth (D3) the canonical UX is "bwa-mem3 mem --meth ref.fa": we
     * auto-append ".meth" to find the converted SEED FM-index built by
     * "bwa-mem3 index --meth" (the original-alphabet index lives at the bare
     * prefix and supplies BNS+PAC via FMI_search::set_meth_ref_prefix below).
     * If the user already passed the ".meth" path directly, use it as-is. */
    char c2t_ref[PATH_MAX];
    char orig_ref_buf[PATH_MAX];
    /* In --meth, the ORIGINAL (un-projected) reference prefix: supplies BNS+PAC
     * for extension/coords (dual-index, set_meth_ref_prefix) and the
     * .hdr/.dict sidecar for @SQ M5/UR enrichment. NULL outside --meth. */
    const char *meth_orig_ref_prefix = NULL;
    const char *ref_prefix = argv[optind];
    if (opt->meth_mode) {
        const char *suffix = ".meth";
        size_t slen = strlen(suffix);
        size_t alen = strlen(argv[optind]);
        int already_c2t = (alen >= slen) &&
                          (strcmp(argv[optind] + alen - slen, suffix) == 0);
        if (!already_c2t) {
            int n = snprintf(c2t_ref, sizeof(c2t_ref), "%s%s", argv[optind], suffix);
            if (n <= 0 || (size_t)n >= sizeof(c2t_ref)) {
                fprintf(stderr, "ERROR: ref path too long for --meth\n");
                exit(EXIT_FAILURE);
            }
            ref_prefix = c2t_ref;
            meth_orig_ref_prefix = argv[optind];
        } else {
            /* User passed the ".meth" seed-index path directly; recover the
             * original reference prefix by stripping the suffix so its
             * sidecar (not the seed index's) supplies @SQ identity tags. */
            size_t base = alen - slen;
            if (base < sizeof(orig_ref_buf)) {
                memcpy(orig_ref_buf, argv[optind], base);
                orig_ref_buf[base] = '\0';
                meth_orig_ref_prefix = orig_ref_buf;
            } else {
                /* Degenerate: the prefix doesn't fit. The original sidecar is
                 * optional (the c2t index alone aligns), so warn and continue
                 * without @SQ M5/UR enrichment rather than aborting. */
                fprintf(stderr,
                        "[bwa-mem3:--meth] WARNING: reference path too long to "
                        "derive the original prefix; @SQ M5/UR enrichment "
                        "skipped.\n");
            }
        }
    }

    /* D3: fail fast with an actionable message if the `.meth` seed index is
     * absent — distinguishing a never-built index from a stale D1 `.bwameth.c2t`
     * index (the format changed; the user must re-run `index --meth`). */
    if (opt->meth_mode) {
        char probe[PATH_MAX];
        snprintf(probe, sizeof(probe), "%s%s", ref_prefix, ".bwt.2bit.64");
        if (access(probe, F_OK) != 0) {
            const char *orig = (meth_orig_ref_prefix != NULL) ? meth_orig_ref_prefix
                                                              : argv[optind];
            char old_probe[PATH_MAX];
            snprintf(old_probe, sizeof(old_probe), "%s.bwameth.c2t.bwt.2bit.64", orig);
            if (access(old_probe, F_OK) == 0)
                fprintf(stderr,
                        "ERROR: --meth found a legacy '.bwameth.c2t' index, but the "
                        "format changed. Re-run: bwa-mem3 index --meth %s\n", orig);
            else
                fprintf(stderr,
                        "ERROR: --meth seed index '%s.*' not found. Run: "
                        "bwa-mem3 index --meth %s\n", ref_prefix, orig);
            free(opt);
            if (out_opened) fclose(aux.fp);
            return 1;
        }
    }

    if (opt->seed_emit_order != SEED_ORDER_OFF)
        fprintf(stderr, "[M::%s] seed order: %s\n", __func__, seed_order_to_str(opt->seed_emit_order));
    if (fast) {
        if (opt->meth_mode)
            /* --skip-contained-ext is set but no-ops under --meth (internal gate), so it is
             * intentionally omitted from the meth audit line to reflect the effective levers.
             * --adaptive-band is set unconditionally and applies under --meth, so it stays. */
            fprintf(stderr, "[M::%s] --fast: -m %d -y %ld --min-ext-len %d --smem-dedup --max-extend-chains %d --adaptive-band -s %d --extend-mate-concordant\n",
                    __func__, opt->max_matesw, (long)opt->max_mem_intv, opt->min_ext_len, opt->max_extend_chains, opt->split_width);
        else
            fprintf(stderr, "[M::%s] --fast: -m %d -y %ld --min-ext-len %d --smem-dedup --skip-contained-ext --max-extend-chains %d --adaptive-band\n",
                    __func__, opt->max_matesw, (long)opt->max_mem_intv, opt->min_ext_len, opt->max_extend_chains);
    }

    /* Load bwt2/FMI index */
    uint64_t tim = __rdtsc();

    fprintf(stderr, "* Ref file: %s\n", ref_prefix);
    aux.fmi = new FMI_search(ref_prefix);
    /* D3 dual-index: the FM-index AND its BNS/PAC come from the `.meth` SEED prefix.
     * The seed BNS is required to decode seed positions into (seed contig, local pos,
     * strand): the seed reference has the f/r-doubled contig layout (r0,f0,r1,f1,...).
     * Seed contigs are remapped to ORIGINAL coordinates arithmetically
     * (orig_tid = seed_rid/2; hypothesis = seed_rid & 1; pos preserved). The ORIGINAL
     * reference's BNS/PAC are loaded separately as the remap/extension target in the
     * extension phase (NOT a replacement of the seed BNS). */
    /* D3 --meth: load the SEED index's FM + bns but NOT its pac. The seed pac is
     * never read in --meth (extension/scoring/mate-rescue use meth_orig_pac);
     * skipping it saves ~1.6 GB on hg38. Outside --meth, load the pac as before. */
    aux.fmi->load_index(/*load_pac=*/!opt->meth_mode);
    aux.shm_base = aux.fmi->shm_attached_base();
    tprof[FMI][0] += __rdtsc() - tim;

    /* D3: load the ORIGINAL reference's bns/pac as resident handles for the
     * (future) extension/scoring phase — distinct from the seed FM-index above.
     * The seed BNS (aux.fmi->idx->bns) is the f/r-doubled converted reference
     * used to decode seed positions; these handles are the un-converted original
     * (real chrom names, N contigs) that extension will score against. No
     * consumer yet (extension is behind the meth-mode checkpoint below); this is
     * a load-only building block. Freed on every exit path the seed index is. */
    if (opt->meth_mode && meth_orig_ref_prefix != NULL) {
        if (meth_orig_ref_load_handles(meth_orig_ref_prefix,
                                       &aux.meth_orig_bns, &aux.meth_orig_pac) != 0) {
            delete aux.fmi;
            free(opt);
            if (out_opened) fclose(aux.fp);
            return 1;
        }
        fprintf(stderr,
                "[bwa-mem3:--meth] original reference bns/pac loaded for "
                "extension (%d contig(s)).\n", aux.meth_orig_bns->n_seqs);
    }

    /* D3 (--meth, PR-3): the early-return seeding checkpoint is REMOVED — the
     * pipeline now runs end to end. Seeds generated against the `.meth` seed
     * FM-index are remapped to ORIGINAL coordinates + an OT/OB hypothesis inside
     * mem_chain_seeds (see meth_seed_to_orig), and chaining/extension/pairing/
     * mate-rescue/output all operate against the ORIGINAL bns/pac carried on
     * worker_t::meth_orig_* (the original reference bases are pac-fetched from
     * `.pac` on demand; no unpacked `.0123` is loaded for either the seed or the
     * original — see below).
     * Extension and mate-rescue score the ORIGINAL read against the original
     * reference with the per-hypothesis asymmetric OT/OB matrix (PR-4 + A1: the
     * batched extension partitions a mixed PE batch by hypothesis; mate rescue
     * routes meth pairs through the scalar ksw_align2 path). The projected read is
     * used only for seeding against the `.meth` FM-index. */
    /* pac-fetch: reconstruct the ORIGINAL reference from its `.pac` on demand
     * (bns_get_seq_v2's ref_string==NULL path) instead of loading the unpacked
     * `.0123` (~6.4 GB on hg38). This is the only reference path: `index` no
     * longer builds `.0123` and `mem` never reads it (byte-identical to the
     * historical `.0123` load, verified on plain + --meth incl. batched rescue). */
    if (opt->meth_mode && meth_orig_ref_prefix != NULL) {
        fprintf(stderr,
                "[bwa-mem3:--meth] pac-fetch: original reference unpacked from "
                ".pac on demand (.0123 not loaded).\n");
    }

    // reference bases are pac-fetched from .pac on demand; no .0123 is loaded
    tim = __rdtsc();
    uint64_t timer;
    if (opt->meth_mode) {
        /* D3 --meth: the SEED `.0123` is dead weight (~13 GB on hg38). Every
         * downstream consumer reads the ORIGINAL unpacked reference via
         * mem_aln_ref_string()/mmc.ref_string (= meth_orig_ref_string, loaded
         * above); seeding uses the FM-index, not the unpacked seed `.0123`. So
         * skip loading it entirely and poison aux.ref_string to NULL — any
         * consumer that bypasses the mem_aln_* helpers will then crash loudly
         * rather than silently read seed bases at original coordinates. */
        ref_string             = NULL;
        aux.ref_string         = NULL;
        aux.ref_string_is_shm  = 0;
        timer = __rdtsc();
        fprintf(stderr, "* [--meth] seed reference `.0123` not loaded "
                "(extension uses the original reference)\n");
    } else {
        /* plain pac-fetch: unpack the original reference from `.pac` on demand
         * (bns_get_seq_v2's ref_string==NULL path). idx->pac is resident on both
         * the disk path (load_pac=true) and the shm-attach path (aliases the
         * staged PAC section). aux.ref_string==NULL routes every consumer —
         * extension + batched mate rescue — to pac-fetch (the rescue copies each
         * window into seqBufRef in-iteration, so the single-live-window contract
         * holds). −6.4 GB on hg38, byte-identical. */
        ref_string             = NULL;
        aux.ref_string         = NULL;
        aux.ref_string_is_shm  = 0;
        timer = __rdtsc();
        fprintf(stderr, "* pac-fetch: reference `.0123` not loaded; "
                "unpacking original bases from .pac on demand\n");
    }
    tprof[REF_IO][0] += timer - tim;

    if (ignore_alt)
        for (i = 0; i < aux.fmi->idx->bns->n_seqs; ++i)
            aux.fmi->idx->bns->anns[i].is_alt = 0;
    /* D3 (--meth, PR-3): --ignore-alt clears is_alt on the SEED bns above, but
     * --meth chaining/extension/output read is_alt from the ORIGINAL bns
     * (aux.meth_orig_bns) — mirror the clear there so --ignore-alt takes effect
     * in --meth. (Original genome indexes rarely carry ALT contigs, but keep the
     * two views consistent.) */
    if (ignore_alt && opt->meth_mode && aux.meth_orig_bns != NULL)
        for (i = 0; i < aux.meth_orig_bns->n_seqs; ++i)
            aux.meth_orig_bns->anns[i].is_alt = 0;

    /* READS file operations */
    ko = kopen(argv[optind + 1], &fd);
	if (ko == 0) {
		fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 1]);
        free(opt);
        if (out_opened)
            fclose(aux.fp);
        delete aux.fmi;
        // kclose(ko);
        return 1;
    }
    // fp = gzopen(argv[optind + 1], "r");
    if (aux.legacy_reader) {
        fp = gzdopen(fd, "r");
        aux.ks = kseq_init(fp);
    } else {
        const char *fr_err = NULL;
        aux.fr1 = fast_reader_dopen(fd, &fr_err);
        if (aux.fr1 == NULL) {
            fprintf(stderr, "[E::%s] %s\n", __func__, fr_err ? fr_err : "failed to open input");
            free(opt);
            if (out_opened) fclose(aux.fp);
            delete aux.fmi;
            kclose(ko);
            return 1;
        }
        aux.frks = fast_kseq_init(aux.fr1);
    }

    // PAIRED_END
    /* Handling Paired-end reads */
    aux.ks2 = 0;
    if (optind + 2 < argc) {
        if (opt->flag & MEM_F_PE) {
            fprintf(stderr, "[W::%s] when '-p' is in use, the second query file is ignored.\n",
                    __func__);
        }
        else
        {
            ko2 = kopen(argv[optind + 2], &fd2);
            if (ko2 == 0) {
                fprintf(stderr, "[E::%s] failed to open file `%s'.\n", __func__, argv[optind + 2]);
                free(opt);
                free(ko);
                if (aux.legacy_reader) { err_gzclose(fp); kseq_destroy(aux.ks); }
                else { fast_kseq_destroy(aux.frks); fast_reader_close(aux.fr1); }
                if (out_opened)
                    fclose(aux.fp);
                delete aux.fmi;
                kclose(ko);
                // kclose(ko2);
                return 1;
            }
            if (aux.legacy_reader) {
                fp2 = gzdopen(fd2, "r");
                aux.ks2 = kseq_init(fp2);
                assert(aux.ks2 != 0);
            } else {
                const char *fr_err2 = NULL;
                aux.fr2 = fast_reader_dopen(fd2, &fr_err2);
                if (aux.fr2 == NULL) {
                    fprintf(stderr, "[E::%s] %s\n", __func__, fr_err2 ? fr_err2 : "failed to open input");
                    free(opt);
                    fast_kseq_destroy(aux.frks); fast_reader_close(aux.fr1);
                    if (out_opened) fclose(aux.fp);
                    delete aux.fmi; kclose(ko); kclose(ko2);
                    return 1;
                }
                aux.frks2 = fast_kseq_init(aux.fr2);
            }
            opt->flag |= MEM_F_PE;
        }
    }

    /* Look up optional per-index header records (<prefix>.hdr or
     * <baseprefix>.dict) once and route them into both output paths. The
     * SAM text path merges these with user -H per lh3/bwa#348 precedence;
     * the --bam path forwards them to htslib's sam_hdr_add_lines so the
     * rich @SQ (AS/M5/SP/AH/…) also makes it into the BAM header. */
    char *idx_hdr_lines = bwa_load_hdr_from_index(ref_prefix);
    /* --meth only: the original (pre-c2t) reference's .hdr/.dict sidecar, for
     * @SQ M5/UR enrichment and @CO/@PG/@RG pass-through in meth_bam_writer_open.
     * NULL outside --meth or when the original has no sidecar. */
    char *meth_orig_hdr_lines = (meth_orig_ref_prefix != NULL)
                                ? bwa_load_hdr_from_index(meth_orig_ref_prefix)
                                : NULL;

    /* Output path:
     *  - --meth: open meth_bam_writer with strand-consolidated SQ headers.
     *    Honors -o/-f (target path) or stdout ("-").
     *  - --bam (no --meth): open generic bam_writer; htslib writes its own
     *    @HD + @SQ + @PG header. Honors -o/-f or stdout.
     *  - SAM text: open -o/-f path (if any) as a FILE*; bwa_print_sam_hdr2. */
    bam_writer_t *bam_writer = NULL;
#ifdef DISABLE_OUTPUT
    /* profile-build (-DDISABLE_OUTPUT) skips ALL filesystem-touching output
     * code so wall-clock measurements aren't gated on disk speed. The
     * per-record write sites below are also #ifndef-guarded; this block
     * additionally skips writer open + header emit so a non-writable -o
     * (or read-only fs) doesn't fail before compute begins. */
    aux.bam_writer = NULL;
    g_meth_bam_writer = NULL;
    if (opt->meth_mode) {
        /* D3 (PR-5): output is native original-alphabet. RNAME/POS/XM/XG all
         * derive from the ORIGINAL bns/pac (loaded into aux.meth_orig_bns/pac
         * by PR-1) and the per-alignment hypothesis — no f/r chrom map, no
         * un-converted ref view. The per-record path needs only the original
         * pac global (the original bns reaches it via mem_aln_bns()). With
         * output disabled there is no writer to open, but the global must still
         * be set so meth_mem_aln_to_bam builds correct XM/coords. Mirror the
         * non-DISABLE_OUTPUT guard below: a NULL original bns/pac (e.g. the
         * degenerate "prefix too long" path never loaded the handles) would
         * otherwise let chaining run against a NULL original reference and
         * silently corrupt coordinates. */
        if (aux.meth_orig_bns == NULL || aux.meth_orig_pac == NULL) {
            fprintf(stderr, "ERROR: meth: original reference (bns/pac) not loaded\n");
            free(opt);
            delete aux.fmi;
            return 1;
        }
        g_meth_orig_pac = aux.meth_orig_pac;
    }
    (void)is_o;
    (void)hdr_line;
    (void)idx_hdr_lines;
    (void)meth_orig_hdr_lines;
    (void)out_path;
#else
    if (opt->meth_mode) {
        /* D3 (PR-5): native original-alphabet output. The @SQ header is built
         * straight from the ORIGINAL (un-converted) bns (aux.meth_orig_bns),
         * and the per-record path consults the original pac (g_meth_orig_pac)
         * for XM:Z; alignments already carry original rids/coords (PR-3) and the
         * hypothesis (XG strand), so there is no f/r chrom map and no
         * un-converted ref fold. */
        if (aux.meth_orig_bns == NULL || aux.meth_orig_pac == NULL) {
            fprintf(stderr, "ERROR: meth: original reference (bns/pac) not loaded\n");
            free(opt);
            delete aux.fmi;
            return 1;
        }
        g_meth_orig_pac = aux.meth_orig_pac;
        const char *meth_out_path = is_o ? out_path : "-";
        extern char *bwa_pg;
        /* meth_orig_hdr_lines is the *original* reference's .hdr/.dict sidecar;
         * the writer enriches each @SQ with its M5/UR tags and forwards
         * @CO/@PG/@RG provenance. */
        g_meth_bam_writer = meth_bam_writer_open(meth_out_path, aux.meth_orig_bns,
                                                 bwa_pg, NULL,
                                                 hdr_line, meth_orig_hdr_lines,
                                                 opt->bam_level);
        if (g_meth_bam_writer == NULL) {
            fprintf(stderr, "ERROR: meth: failed to open BAM writer for '%s'\n", meth_out_path);
            g_meth_orig_pac = NULL;
            free(opt);
            delete aux.fmi;
            return 1;
        }
    } else if (opt->bam_mode) {
        const char *bam_path = is_o ? out_path : "-";
        extern char *bwa_pg;
        /* Suppress idx .hdr/.dict records entirely when the user's -H
         * supplies any @SQ, matching bwa_print_sam_hdr2's SAM precedence. */
        const char *bam_idx_hdr = idx_hdr_lines;
        if (hdr_line != NULL) {
            if (strncmp(hdr_line, "@SQ\t", 4) == 0 ||
                strstr(hdr_line, "\n@SQ\t") != NULL)
                bam_idx_hdr = NULL;
        }
        bam_writer = bam_writer_open(bam_path, aux.fmi->idx->bns,
                                     bam_idx_hdr, hdr_line,
                                     bwa_pg, opt->bam_level);
        if (bam_writer == NULL) {
            fprintf(stderr, "ERROR: failed to open BAM writer at '%s'\n", bam_path);
            free(opt);
            delete aux.fmi;
            return 1;
        }
        aux.bam_writer = bam_writer;
    } else {
        aux.bam_writer = NULL;
        if (is_o) {
            aux.fp = fopen(out_path, "w");
            if (aux.fp == NULL) {
                fprintf(stderr, "Error: can't open %s output file\n", out_path);
                free(opt);
                delete aux.fmi;
                return 1;
            }
            out_opened = true;
        }
        bwa_print_sam_hdr2(aux.fmi->idx->bns, idx_hdr_lines, hdr_line, aux.fp);
    }
#endif

    if (fixed_chunk_size > 0) {
        /* -K: the user pinned the batch size (for reproducibility) — honor it
         * exactly, never cap. */
        aux.task_size = fixed_chunk_size;
    } else {
        /* Default batch size is chunk_size (~10M bases) per thread. That keeps
         * each thread well-fed, but at very high -t it makes a single chunk
         * enormous (10M * 192 ~= 1.9G bases), so the input is only ~3-4 chunks
         * and the pipeline starves: the first chunk's read and the last chunk's
         * write don't overlap anything, leaving cores idle (fill/drain). Cap the
         * default so high -t still produces enough chunks to keep read/compute/
         * write overlapped, while keeping each chunk far above the ~33k-pairs
         * floor below which per-chunk overhead (pestat/barriers) starts to bite.
         * Output stays identical for -t small enough that the cap doesn't engage
         * (scaled <= cap); above that, batch composition changes exactly as -K
         * would (validated to leave proper-pair rate unchanged). The cap is
         * overridable via BWA_MEM3_CHUNK_CAP (bases; <=0 disables, for sweeps). */
        int64_t cap = 256000000;
        const char *cap_env = getenv("BWA_MEM3_CHUNK_CAP");
        if (cap_env && *cap_env) cap = (int64_t)atoll(cap_env);
        int64_t scaled = (int64_t)opt->chunk_size * (int64_t)opt->n_threads;
        aux.task_size = (cap > 0 && scaled > cap) ? cap : scaled;
    }
    tprof[MISC][1] = opt->chunk_size = aux.actual_chunk_size = aux.task_size;

    /* Pipeline depth. The 3-step pipeline (read / process / write) is gated by
     * how many of those steps can run at once. With 2 workers only 2 of the 3
     * run concurrently, so the single I/O-side worker serialises read(N+1) and
     * write(N-1) around the compute worker; at high thread counts that
     * read+write chain, not compute, binds the wall. A 3rd worker lets
     * read || process || write triple-overlap. It is not oversubscription: the
     * step mutex still admits only one worker into the compute step at a time,
     * so the extra worker only ever does single-threaded I/O. Cost is one more
     * chunk in flight. -1 (no_mt_io) forces single-threaded I/O as before. */
    const int pipe_workers = no_mt_io ? 1 : (opt->n_threads > 2 ? 3 : opt->n_threads);

    double sp_t0 = 0.0;
#ifdef STAGE_PROF
    /* stage_prof: arm per-chunk read/process/write profiling if --profile given */
    if (profile_path && profile_path[0]) {
#if defined(__x86_64__) || defined(_M_X64)
        const char *sp_arch = "x86_64";
#elif defined(__aarch64__) || defined(__arm64__)
        const char *sp_arch = "arm64";
#else
        const char *sp_arch = "unknown";
#endif
        const char *sp_in = (optind + 1 < argc) ? argv[optind + 1]
                          : (optind < argc ? argv[optind] : "");
        sp_init(profile_path, "bwa-mem3", PACKAGE_VERSION, sp_arch, opt->n_threads,
                opt->bam_mode ? "bam" : "sam",
                opt->bam_mode ? opt->bam_level : -1, sp_in);
        sp_set_workers(pipe_workers);   /* pipeline depth (process() arg) */
        sp_t0 = sp_wall();
    }
#endif

    tim = __rdtsc();

    /* Relay process function */
    process(&aux, fp, fp2, pipe_workers);

    tprof[PROCESS][0] += __rdtsc() - tim;

    if (sp_enabled()) {
        struct rusage ru; getrusage(RUSAGE_SELF, &ru);
#ifdef __linux__
        double rss_mb = (double)ru.ru_maxrss / 1024.0;          /* Linux: ru_maxrss is KiB */
#else
        double rss_mb = (double)ru.ru_maxrss / (1024.0*1024.0); /* macOS/BSD: ru_maxrss is bytes */
#endif
        /* mean_cores_busy is computed inside sp_finish from per-chunk proc_cpu */
        sp_finish(sp_wall() - sp_t0, 0.0, rss_mb);
    }

    /* Close meth BAM writer BEFORE free(opt) — opt->meth_mode is checked here. */
    int meth_mode_local = opt->meth_mode;

    // free memory
    int32_t nt = aux.opt->n_threads;
    if (!aux.ref_string_is_shm) {
        _mm_free(ref_string);
    }
    /* When ref_string aliases shm, the segment is owned by the loader
     * process; the kernel reclaims the mapping at process exit. */
    free(hdr_line);
    free(idx_hdr_lines);
    free(meth_orig_hdr_lines);
    free(opt);
    if (aux.legacy_reader) { kseq_destroy(aux.ks); err_gzclose(fp); }
    else { fast_kseq_destroy(aux.frks); fast_reader_close(aux.fr1); }
    kclose(ko);

    // PAIRED_END
    if (aux.ks2 || aux.fr2) {
        if (aux.legacy_reader) { kseq_destroy(aux.ks2); err_gzclose(fp2); }
        else { fast_kseq_destroy(aux.frks2); fast_reader_close(aux.fr2); }
        kclose(ko2);
    }

    /* BGZF flush + EOF marker errors surface only on close. Propagate to the
     * exit code so a truncated BAM doesn't masquerade as a successful run. */
    int exit_code = 0;
    if (meth_mode_local && g_meth_bam_writer != NULL) {
        int rc = meth_bam_writer_close(g_meth_bam_writer);
        if (rc != 0) {
            fprintf(stderr, "[meth] ERROR: BAM writer close rc=%d\n", rc);
            exit_code = 1;
        }
        g_meth_bam_writer = NULL;
    }
    /* D3 (PR-5): the original pac global is a borrowed pointer into
     * aux.meth_orig_pac (freed by meth_orig_ref_free_handles below); just
     * clear it. The retired f/r chrom map and un-converted ref fold no longer
     * exist. */
    g_meth_orig_pac = NULL;
    if (out_opened) {
        fclose(aux.fp);
    }

    if (bam_writer != NULL) {
        if (bam_writer_close(bam_writer) != 0) {
            fprintf(stderr, "[bam_writer] ERROR: close returned non-zero\n");
            exit_code = 1;
        }
        aux.bam_writer = NULL;
    }

    // new bwt/FMI
    /* D3 (--meth, PR-3): free the original-reference handles (bns/pac/.0123)
     * alongside the seed FM-index. NULL/no-op outside --meth; idempotent and
     * NULL-safe. Now reached on every run since the seeding checkpoint is gone. */
    meth_orig_ref_free_handles(&aux);
    delete(aux.fmi);

    /* Display runtime profiling stats */
    tprof[MEM][0] = __rdtsc() - tprof[MEM][0];
    display_stats(nt);

    return exit_code;
}
