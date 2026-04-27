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

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>
*****************************************************************************************/

#ifndef _MACRO_HPP
#define _MACRO_HPP

#include <stdio.h>

#define VER 0
#define printf_(x,y...)								\
	{												\
		if(x)										\
			fprintf(stderr, y);						\
	}

/* Branch-prediction hints for hot paths. Available to every TU that
 * includes macro.h. ksw.cpp and kswv.h predate this and have their own
 * local copies — guarded so there's no redefinition warning. */
#ifndef LIKELY
#  if defined(__GNUC__)
#    define LIKELY(x)   __builtin_expect(!!(x), 1)
#    define UNLIKELY(x) __builtin_expect(!!(x), 0)
#  else
#    define LIKELY(x)   (x)
#    define UNLIKELY(x) (x)
#  endif
#endif

/* Note: BSW-specific macros are in src/bandedSWA.h file */

#define H0_ -99
#define SEEDS_PER_READ 500           /* Avg seeds per read */
#define MAX_SEEDS_PER_READ 500       /* Max seeds per read */
#define AVG_SEEDS_PER_READ 64        /* Used for storing seeds in chains*/

// Average bases per read, used as a coarse upper-bound estimate (#reads
// per chunk) for per-run regs/chain_ar/seedBuf allocation. Not a
// correctness bound — the SMEM and related per-thread buffers are sized
// at batch time from the observed maximum read length and grown on demand.
#define NREADS_ESTIMATE_AVG_BASES 100

/* BWAMEM_BATCHED_MATESW:
 *   1 -> worker_sam takes the batched mate-rescue SW path
 *        (mem_sam_pe_batch_pre / mem_sam_pe_batch / mem_sam_pe_batch_post,
 *        feeding kswv::getScores8 / getScores16).
 *   0 -> worker_sam takes the legacy scalar mem_sam_pe + ksw_align2 path.
 *
 * Historically gated on __AVX512BW__ only, which routed non-AVX-512 builds
 * (ARM, AVX2-only x86) to the scalar path even though batched kernels can
 * be implemented for those architectures. As of the NEON + AVX2 ports this
 * gate accepts any arch with a batched kswv kernel.
 * DISABLE_BATCHED_MATESW is an escape hatch for the A/B test in CI. */
#ifndef BWAMEM_BATCHED_MATESW
  #if DISABLE_BATCHED_MATESW
    #define BWAMEM_BATCHED_MATESW 0
  #elif __AVX512BW__ || __AVX2__ \
        || defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
    #define BWAMEM_BATCHED_MATESW 1
  #else
    #define BWAMEM_BATCHED_MATESW 0
  #endif
#endif

/* Apple Silicon has larger L2 caches (4-16MB per cluster) and benefits from
 * larger batch sizes to better utilize cache locality. M1/M2/M3/M4 all have
 * significantly more L2 cache per core than typical x86 consumer CPUs. */
#if defined(APPLE_SILICON) || defined(__aarch64__) || defined(__ARM_NEON)
#define BATCH_SIZE 1024              /* Larger batch for Apple Silicon's big L2 cache */
#define BATCH_MUL 40                 /* Also scale BATCH_MUL proportionally */
#else
#define BATCH_SIZE 512               /* Block of reads allocated to a thread for processing */
#define BATCH_MUL 20
#endif /* APPLE_SILICON batch size */

#define SEEDS_PER_CHAIN 1
#define N_SMEM_KERNEL 3

#define SEQ_LEN8 128   // redundant??

#define MAX_LINE_LEN 256
#define CACHE_LINE 16        // 16 INT32
#define ALIGN_OFF 1

#define MAX_THREADS 256
/* LIM_R must stay >= UGP_L_CAT_FIN_END (the high-watermark UGP counter
 * index — currently 279, defined further down in this file). Bumped
 * 128 → 256 → 384 as UGP_* instrumentation grew; kept on a 64-multiple
 * so memset-init stays cheap. The static_assert near UGP_L_CAT_FIN_END
 * enforces the invariant at compile time — bump this number if the
 * assert ever fires. */
#define LIM_R 384
/* LIM_C is the per-thread (column) dimension of tprof[][]. Must stay
 * >= MAX_THREADS or any user invocation with -t > LIM_C corrupts
 * adjacent globals via tprof[][tid] writes from the UGP_* and
 * per-stage profiling sites. Sized to MAX_THREADS so the two move
 * together. */
#define LIM_C 256

#define SA_COMPRESSION 1
#define SA_COMPX 03 // (= power of 2)
#define SA_COMPX_MASK 0x7    // 0x7 or 0x3 or 0x1

/*** Runtime profiling macros ***/
#define INDEX 0
#define MEM 1
#define MEM2 2
#define MEM3 4
#define SAM1 5
#define SAM2 3
#define SAM3 7
#define MPI_TIME 8
#define MEM_PROCESS10 9
#define MEM_PROCESS2 10
#define READ_IO 11
#define PROCESS 12
#define REF_IO 13
#define PREPROCESS 14
#define CONVERT 15
#define MPI_TIME_SUM 16
#define OUTPUT 17
#define MPI_TIME_MIN 18
#define POST_SWA 19
#define MPI_TIME_MAX 20
#define SAM_IO 21
#define ALIGN1 22

#define KT_FOR 24
#define KTF_WORKER 26
#define WORKER20 28
#define WORKER21 30
#define WORKER10 32
#define WORKER11 34
#define MEM_ALN 36
#define MEM_CHAIN 38
#define MEM_COLLECT 40
#define BWT_REVERSE 41
#define BWA_BUILD 42
#define PACKED 43
#define SA 44
#define BWT_REVERSE_A 45
#define BWT_REVERSE_B 46
#define MEM_SA 47
#define MEM_ALN2 48
#define MEM_ALN2_A 49
#define MEM_ALN2_B 50
#define MEM_ALN2_C 51
#define EXTEND 52
#define FORWARD 53
#define MEM_CHAIN1 54
#define MEM_CHAIN2 55
#define SMEM1 56
#define SMEM2 57
#define SMEM3 58
#define BWT_FORWARD_A 59
#define STR 60
#define MISC 61
#define MEM_ALN2_UP 62
#define BWT_FORWARD_B 63
#define CLEFT 64
#define CRIGHT 65
#define MEM_ALN_M1 66
#define MEM_ALN_M2 67
#define MEM_SA_BLOCK 68
#define SEQ_FETCH 69
#define MEM_ALN2_PRE 70
#define QLEN 71
#define TLEN 72
#define CNT 73
#define WAVG 74
#define WCNT 75
#define WMAX 76
#define WMIN 77
#define KSW 78
#define PE 79
#define PESW 80
#define PESORT 81
#define INTROSORT 82
#define PE1 83
#define PE3 84
#define PE2 85
#define PE4 86
#define PE5 87
#define PE6 88
#define PE7 89
#define PE8 90
#define PE11 91
#define PE12 92
#define PE13 93
#define PE14 94
#define PE15 95
#define PE16 96
#define PE17 97
#define PE18 98
#define PE19 99
#define PE20 100
#define PE21 101
#define PE22 102
#define PE23 103
#define MEM_ALN2_DOWN 104
#define MEM_ALN2_DOWN1 105
#define SORT 106
#define FMI 107
#define MEM_ALN2_D 108
#define MEM_ALN2_E 109
#define PE24 110
#define PE25 111
#define PE26 112

// PR 26a/c.1/e diagnostics. Per-thread, aggregated in display_stats().
#define UGP_L_ATTEMPT     113
#define UGP_L_HIT         114
#define UGP_L_TIGHT       115
#define UGP_L_TB_1_8      116
#define UGP_L_TB_9_32     117
#define UGP_L_TB_33_MAX   118
#define UGP_L_DISP_NARROW 119
#define UGP_R_ATTEMPT     120
#define UGP_R_HIT         121
#define UGP_R_TIGHT       122
#define UGP_R_TB_1_8      123
#define UGP_R_TB_9_32     124
#define UGP_R_TB_33_MAX   125
#define UGP_R_DISP_NARROW 126

/* PR 26c.2 design instrumentation. Bucketed counters accessed via base+offset
 * to avoid macro sprawl. Each group covers LEFT (side=0) then RIGHT (side=1).
 *
 *   Group A: Fine tight_band distribution (informs adaptive banding).
 *     10 bins per side: {0, 1-2, 3-4, 5-8, 9-16, 17-32, 33-48, 49-64, 65-80, 81-100}.
 *     Index: UGP_FINE_BASE + side*UGP_FINE_NBINS + bin_idx.
 *
 *   Group B: Per-tier band breakdown (informs 16-bit tier bucketing decision).
 *     4 bins per (side, tier ∈ {8-bit, 16-bit}): {tb=0, tb 1-8, tb 9-32, tb 33+}.
 *     Index: UGP_TIER_TB_BASE + side*8 + tier*4 + bin_idx.
 *
 *   Group C: Per-worker-batch narrow bucket size (informs MIN_BUCKET threshold).
 *     8 bins per side: {0, 1-15, 16-31, 32-63, 64-127, 128-255, 256-511, 512+}.
 *     Index: UGP_NARROW_SZ_BASE + side*UGP_NARROW_SZ_NBINS + bin_idx.
 */
#define UGP_FINE_NBINS         10
#define UGP_FINE_BASE          127
#define UGP_FINE_END           (UGP_FINE_BASE + 2 * UGP_FINE_NBINS)            /* 147 */

#define UGP_TIER_TB_BASE       UGP_FINE_END                                    /* 147 */
#define UGP_TIER_TB_END        (UGP_TIER_TB_BASE + 2 * 2 * 4)                  /* 163 */

#define UGP_NARROW_SZ_NBINS    8
#define UGP_NARROW_SZ_BASE     UGP_TIER_TB_END                                 /* 163 */
#define UGP_NARROW_SZ_END      (UGP_NARROW_SZ_BASE + 2 * UGP_NARROW_SZ_NBINS)  /* 179 */

/* Q1+Q2 instrumentation. Counts the FINAL (committed) outcome of each
 * extension — incremented at the per-seed HIT site and at every post-SW
 * retry-collect commit branch. Pairs that loop back for a wider band are
 * NOT counted at the intermediate iteration; only when their decision is
 * final.
 *
 *   Q1 (UGP_OUTCOME_BASE): ungapped vs gapped. Proxy for SW result is
 *     sp->qle == sp->tle (no net query/target offset). Not 100% rigorous
 *     — an I+D of equal length cancels — but accurate for short reads.
 *     HIT path is always ungapped by construction.
 *     Index: UGP_OUTCOME_BASE + side * 2 + (gapped ? 1 : 0)
 *
 *   Q2 (UGP_SCORE_HIST_BASE): per-side alignment-score histogram.
 *     8 bins: {0-10, 11-25, 26-50, 51-75, 76-100, 101-125, 126-150, 151+}.
 *     Score is a->score (HIT) or sp->score (SW commit) — both include h0.
 *     Index: UGP_SCORE_HIST_BASE + side * UGP_SCORE_HIST_NBINS + bin_idx
 */
#define UGP_OUTCOME_BASE       UGP_NARROW_SZ_END                                 /* 179 */
#define UGP_OUTCOME_END        (UGP_OUTCOME_BASE + 2 * 2)                        /* 183 */
#define UGP_L_UNGAPPED         (UGP_OUTCOME_BASE + 0)
#define UGP_L_GAPPED           (UGP_OUTCOME_BASE + 1)
#define UGP_R_UNGAPPED         (UGP_OUTCOME_BASE + 2)
#define UGP_R_GAPPED           (UGP_OUTCOME_BASE + 3)

#define UGP_SCORE_HIST_NBINS   8
#define UGP_SCORE_HIST_BASE    UGP_OUTCOME_END                                   /* 183 */
#define UGP_SCORE_HIST_END     (UGP_SCORE_HIST_BASE + 2 * UGP_SCORE_HIST_NBINS)  /* 199 */

/* Q3 instrumentation: per-category delta-from-perfect histograms for LEFT
 * extensions split into 5 categories.
 *
 *   cat 0  ALL                            every LEFT extension
 *   cat 1  UNGAP_FINAL                    HIT, or sp->qle == sp->tle on SW commit
 *   cat 2  GAPPED_FINAL                   sp->qle != sp->tle on SW commit
 *   cat 3  HIT                            ungapped fast-path returned HIT
 *   cat 4  UNGAP_FINAL_NOT_HIT            cat1 ∩ ~cat3
 *
 * Each pair contributes its delta = (h0 + a*len2) − aln_score, where
 * (h0 + a*len2) is the score of a hypothetical perfect ungapped extension.
 * delta == 0 means HIT with all matches; delta grows with mismatches/gaps.
 *
 * Two histograms per category, 8 bins each (ugp_delta_bin: 0, 1-5, 6-10,
 * 11-25, 26-50, 51-75, 76-100, 101+):
 *
 *   UGP_L_CAT_UNG  delta on would-be ungapped extension score
 *                  HIT contributes (perfect − fp_score); non-HIT contributes
 *                  (perfect − sp->ugp_walk_score) using the walk score
 *                  computed at LEFT queue point.
 *
 *   UGP_L_CAT_FIN  delta on final committed alignment score
 *                  HIT contributes (perfect − fp_score); non-HIT contributes
 *                  (perfect − sp->score).
 */
#define UGP_CAT_NBINS          8
#define UGP_CAT_NCAT           5
#define UGP_L_CAT_UNG_BASE     UGP_SCORE_HIST_END                                  /* 199 */
#define UGP_L_CAT_UNG_END      (UGP_L_CAT_UNG_BASE + UGP_CAT_NCAT * UGP_CAT_NBINS) /* 239 */
#define UGP_L_CAT_FIN_BASE     UGP_L_CAT_UNG_END                                   /* 239 */
#define UGP_L_CAT_FIN_END      (UGP_L_CAT_FIN_BASE + UGP_CAT_NCAT * UGP_CAT_NBINS) /* 279 */

/* Compile-time tripwire: any new UGP_* counter group that pushes the high
 * watermark past LIM_R would silently overflow tprof[]/prof[] at runtime.
 * UGP_L_CAT_FIN_END is the largest UGP index in use today (=279).
 *
 * Note: this header is included only from C++ TUs (every consumer ends in
 * .cpp), so the C++11 `static_assert` keyword resolves cleanly. If a C TU
 * ever pulls macro.h in, swap to `_Static_assert` (C11) or include
 * <assert.h> first. */
static_assert(UGP_L_CAT_FIN_END <= LIM_R,
              "LIM_R too small for UGP counters; bump LIM_R in macro.h");

#endif

