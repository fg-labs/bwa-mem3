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
#include "bwa_shm.h"

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

        /* Read "reads" from input file (fread) */
        int64_t sz = 0;
        ret->seqs = bseq_read_orig(aux->task_size,
                                   &ret->n_seqs,
                                   aux->ks, aux->ks2,
                                   &sz);

        tprof[READ_IO][0] += __rdtsc() - tim;

        fprintf(stderr, "[0000] read_chunk: %ld, work_chunk_size: %ld, nseq: %d\n",
                aux->task_size, sz, ret->n_seqs);

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
                
        aux->n_processed += ret->n_seqs;
        return ret;
    }
    /* Step 3: Write output */
    else if (step == 2)
    {
        uint64_t tim = __rdtsc();

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
                free(ret->seqs[i+k].name);
                free(ret->seqs[i+k].comment);
                free(ret->seqs[i+k].seq);
                free(ret->seqs[i+k].qual);
                free(ret->seqs[i+k].sam);
                free(ret->seqs[i+k].bams);
            }
            i += group_size;
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
    w.fmi = aux->fmi;
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
    fprintf(stderr, "    -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [%.2f]\n", opt->drop_ratio);
    fprintf(stderr, "    -W INT        discard a chain if seeded bases shorter than INT [0]\n");
    fprintf(stderr, "    -m INT        perform at most INT rounds of mate rescues for each read [%d]\n", opt->max_matesw);
    fprintf(stderr, "    -S            skip mate rescue\n");
    fprintf(stderr, "    -P            skip pairing; mate rescue performed unless -S also in use\n");
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
    fprintf(stderr, "                 with `bwa-mem3 index --meth` (emits ref.fa.bwameth.c2t).\n");
    fprintf(stderr, "   --set-as-failed f|r\n");
    fprintf(stderr, "                 flag alignments to the matching strand ('f' or 'r') as QC-fail (0x200)\n");
    fprintf(stderr, "   --do-not-penalize-chimeras\n");
    fprintf(stderr, "                 disable the longest-match <44%% chimera heuristic (no 0x200 / MAPQ cap)\n");
    fprintf(stderr, "Supplementary MAPQ rescoring (fg-labs extension):\n");
    fprintf(stderr, "   --supp-rep-hard-cap INT\n");
    fprintf(stderr, "                 force MAPQ=0 for supplementary alignments whose chain contains any seed\n");
    fprintf(stderr, "                 with >=INT genome occurrences (i.e. the supp region is repetitive on its\n");
    fprintf(stderr, "                 own). 0 disables (default). Typical values 5-20; lower = more aggressive.\n");
    fprintf(stderr, "                 Primary MAPQ is unaffected.\n");
    fprintf(stderr, "Help:\n");
    fprintf(stderr, "   --help        print this help message and exit\n");
    fprintf(stderr, "Note: Please read the man page for detailed description of the command line and options.\n");
}

/* Resolve `<prefix>.0123` to a `uint8_t *` view of the ref string. If
 * `shm_base` is non-NULL, the FMI loader already attached this segment;
 * resolve REF_STRING via that single mapping. Otherwise fall back to the
 * disk slurp. The two paths must agree: if FMI came from disk, ref string
 * must come from disk too, so we never call `bwa_shm_attach` here. */
static uint8_t *load_ref_string(const char *prefix, uint8_t *shm_base,
                                int64_t *rlen_out, int *is_shm_out)
{
    *is_shm_out = 0;
    *rlen_out   = 0;

    if (shm_base != NULL) {
        uint64_t off = 0, sz = 0;
        if (bwa_shm_section_find(shm_base, BWA_SHM_SEC_REF_STRING, &off, &sz) != 0) {
            fprintf(stderr,
                "ERROR: shm segment for '%s' is missing REF_STRING; aborting.\n"
                "       The segment was staged by an older bwa-mem3; drop and re-stage.\n",
                prefix);
            exit(EXIT_FAILURE);
        }
        *is_shm_out = 1;
        *rlen_out   = (int64_t)sz;
        fprintf(stderr, "* Reference genome attached from shm: %lld bytes\n",
                (long long)sz);
        return shm_base + off;
    }

    /* Disk path. */
    char binary_seq_file[PATH_MAX];
    strcpy_s(binary_seq_file, PATH_MAX, prefix);
    strcat_s(binary_seq_file, PATH_MAX, ".0123");

    fprintf(stderr, "* Binary seq file = %s\n", binary_seq_file);
    FILE *fr = fopen(binary_seq_file, "r");
    if (fr == NULL) {
        fprintf(stderr, "Error: can't open %s input file\n", binary_seq_file);
        return NULL;
    }
    int64_t rlen = 0;
    if (fseek(fr, 0, SEEK_END) != 0) {
        fprintf(stderr, "Error: fseek failed on %s\n", binary_seq_file);
        fclose(fr);
        return NULL;
    }
    rlen = ftell(fr);
    if (rlen <= 0) {
        fprintf(stderr, "Error: %s is empty or unseekable (ftell=%lld)\n",
                binary_seq_file, (long long)rlen);
        fclose(fr);
        return NULL;
    }
    uint8_t *buf = (uint8_t*) _mm_malloc(rlen, 64);
    assert_not_null(buf, rlen, rlen);
    bwamem_madv_hugepage(buf, rlen);
    rewind(fr);
    err_fread_noeof(buf, 1, rlen, fr);
    fclose(fr);

    *rlen_out = rlen;
    return buf;
}

int main_mem(int argc, char *argv[])
{
    int          i, c, ignore_alt = 0, no_mt_io = 0;
    int          fixed_chunk_size          = -1;
    char        *p, *rg_line               = 0, *hdr_line = 0;
    const char  *mode                      = 0;

    mem_opt_t    *opt, opt0;
    gzFile        fp, fp2 = 0;
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
    //   --meth                      Enable inline bwameth-style c2t + post-processing + BAM output.
    //                               Expects a reference built with `bwa-mem3 index --meth`.
    //   --set-as-failed f|r         Flag alignments to this strand as QC-fail (0x200)
    //   --do-not-penalize-chimeras  Skip the longest-match <44% chimera heuristic
    enum {
        OPT_BAM = 1000,
        OPT_METH,
        OPT_METH_SET_AS_FAILED,
        OPT_METH_NO_CHIMERA,
        OPT_SUPP_REP_HARD_CAP,
        OPT_HELP,
    };
    static struct option long_opts[] = {
        {"bam",                      optional_argument, 0, OPT_BAM},
        {"meth",                     no_argument,       0, OPT_METH},
        {"set-as-failed",            required_argument, 0, OPT_METH_SET_AS_FAILED},
        {"do-not-penalize-chimeras", no_argument,       0, OPT_METH_NO_CHIMERA},
        {"supp-rep-hard-cap",        required_argument, 0, OPT_SUPP_REP_HARD_CAP},
        {"help",                     no_argument,       0, OPT_HELP},
        {0, 0, 0, 0}
    };
    while ((c = getopt_long(argc, argv, "51qpaMCSPVYjuk:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:W:x:G:h:y:K:X:H:o:f:z:",
                            long_opts, NULL)) >= 0)
    {
        if (c == 'k') opt->min_seed_len = atoi(optarg), opt0.min_seed_len = 1;
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
        else if (c == OPT_METH) {
            opt->meth_mode = 1;
            opt->bam_mode = 1;  /* meth implies BAM output */
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
        else if (c == OPT_METH_NO_CHIMERA) {
            opt->meth_no_chim = 1;
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

    /* Meth-mode default tuning. bwameth.py runs bwa-mem3 with
     * -B 2 -L 10 -U 100 -T 40 -CM — these reduce mismatch and soft-clip
     * penalties so BS reads get long un-clipped alignments, raise the
     * output score threshold, mark shorter hits as secondary, and
     * pass the YS/YC comment tags through to SAM. We apply the same
     * defaults when --meth is set, modulo explicit CLI overrides. */
    if (opt->meth_mode) {
        if (!opt0.b)            opt->b           = 2;
        if (!opt0.pen_clip5)    opt->pen_clip5   = 10;
        if (!opt0.pen_clip3)    opt->pen_clip3   = 10;
        if (!opt0.pen_unpaired) opt->pen_unpaired= 100;
        if (!opt0.T)            opt->T           = 40;
        opt->flag |= MEM_F_NO_MULTI;   /* -M */
        aux.copy_comment = 1;          /* -C, needed for YS:Z/YC:Z passthrough */
    }

    /* Matrix for SWA */
    bwa_fill_scmat(opt->a, opt->b, opt->mat);

    /* In --meth the canonical UX is "bwa-mem3 mem --meth ref.fa" and
     * we auto-append ".bwameth.c2t" to find the index built by
     * "bwa-mem3 index --meth". If the user (or bwameth.py's internal
     * invocation) already passed the ".bwameth.c2t" path directly, use
     * it as-is rather than double-appending. */
    char c2t_ref[PATH_MAX];
    const char *ref_prefix = argv[optind];
    if (opt->meth_mode) {
        const char *suffix = ".bwameth.c2t";
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
        }
    }

    /* Load bwt2/FMI index */
    uint64_t tim = __rdtsc();

    fprintf(stderr, "* Ref file: %s\n", ref_prefix);
    aux.fmi = new FMI_search(ref_prefix);
    aux.fmi->load_index();
    aux.shm_base = aux.fmi->shm_attached_base();
    tprof[FMI][0] += __rdtsc() - tim;

    // reading ref string (from shm if FMI attached, else from .0123 file)
    tim = __rdtsc();
    fprintf(stderr, "* Reading reference genome..\n");
    int64_t rlen = 0;
    int     ref_is_shm = 0;
    ref_string = load_ref_string(ref_prefix, aux.shm_base, &rlen, &ref_is_shm);
    if (ref_string == NULL) {
        exit(EXIT_FAILURE);
    }
    aux.ref_string         = ref_string;
    aux.ref_string_is_shm  = ref_is_shm;
    uint64_t timer = __rdtsc();
    tprof[REF_IO][0] += timer - tim;
    fprintf(stderr, "* Reference genome size: %ld bp\n", (long)rlen);
    fprintf(stderr, "* Done reading reference genome !!\n\n");

    if (ignore_alt)
        for (i = 0; i < aux.fmi->idx->bns->n_seqs; ++i)
            aux.fmi->idx->bns->anns[i].is_alt = 0;

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
    fp = gzdopen(fd, "r");
    aux.ks = kseq_init(fp);

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
                err_gzclose(fp);
                kseq_destroy(aux.ks);
                if (out_opened)
                    fclose(aux.fp);
                delete aux.fmi;
                kclose(ko);
                // kclose(ko2);
                return 1;
            }
            // fp2 = gzopen(argv[optind + 2], "r");
            fp2 = gzdopen(fd2, "r");
            aux.ks2 = kseq_init(fp2);
            opt->flag |= MEM_F_PE;
            assert(aux.ks2 != 0);
        }
    }

    /* Look up optional per-index header records (<prefix>.hdr or
     * <baseprefix>.dict) once and route them into both output paths. The
     * SAM text path merges these with user -H per lh3/bwa#348 precedence;
     * the --bam path forwards them to htslib's sam_hdr_add_lines so the
     * rich @SQ (AS/M5/SP/AH/…) also makes it into the BAM header. */
    char *idx_hdr_lines = bwa_load_hdr_from_index(ref_prefix);

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
        g_meth_cmap = meth_chrom_map_build_from_bns(aux.fmi->idx->bns);
        /* g_meth_cmap is consulted by per-record paths even with output
         * disabled; build it so meth_mode tagging stays consistent. */
        if (g_meth_cmap == NULL) {
            fprintf(stderr, "ERROR: meth: failed to build chrom map\n");
            free(opt);
            delete aux.fmi;
            return 1;
        }
    }
    (void)is_o;
    (void)hdr_line;
    (void)idx_hdr_lines;
    (void)out_path;
#else
    if (opt->meth_mode) {
        g_meth_cmap = meth_chrom_map_build_from_bns(aux.fmi->idx->bns);
        if (g_meth_cmap == NULL) {
            fprintf(stderr, "ERROR: meth: failed to build chrom map\n");
            free(opt);
            delete aux.fmi;
            return 1;
        }
        const char *meth_out_path = is_o ? out_path : "-";
        extern char *bwa_pg;
        g_meth_bam_writer = meth_bam_writer_open(meth_out_path, g_meth_cmap, bwa_pg, NULL,
                                                 opt->bam_level);
        if (g_meth_bam_writer == NULL) {
            fprintf(stderr, "ERROR: meth: failed to open BAM writer for '%s'\n", meth_out_path);
            meth_chrom_map_free(g_meth_cmap); g_meth_cmap = NULL;
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

    if (fixed_chunk_size > 0)
        aux.task_size = fixed_chunk_size;
    else {
        //aux.task_size = 10000000 * opt->n_threads; //aux.actual_chunk_size;
        aux.task_size = opt->chunk_size * opt->n_threads; //aux.actual_chunk_size;
    }
    tprof[MISC][1] = opt->chunk_size = aux.actual_chunk_size = aux.task_size;

    tim = __rdtsc();

    /* Relay process function */
    process(&aux, fp, fp2, no_mt_io? 1:2);

    tprof[PROCESS][0] += __rdtsc() - tim;

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
    free(opt);
    kseq_destroy(aux.ks);
    err_gzclose(fp); kclose(ko);

    // PAIRED_END
    if (aux.ks2) {
        kseq_destroy(aux.ks2);
        err_gzclose(fp2); kclose(ko2);
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
    /* Free g_meth_cmap independently of g_meth_bam_writer: under
     * -DDISABLE_OUTPUT the writer is never opened, but the chrom map is
     * still built so per-record paths see consistent tagging. The
     * branch above only fires when the writer exists, so freeing the
     * map there would leak on the DISABLE_OUTPUT path. */
    if (meth_mode_local && g_meth_cmap != NULL) {
        meth_chrom_map_free(g_meth_cmap);
        g_meth_cmap = NULL;
    }
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
    delete(aux.fmi);

    /* Display runtime profiling stats */
    tprof[MEM][0] = __rdtsc() - tprof[MEM][0];
    display_stats(nt);

    return exit_code;
}
