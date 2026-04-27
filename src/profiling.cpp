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

#include <stdio.h>
#include "macro.h"
#include <stdint.h>
#include <assert.h>
#include "bwa.h"        /* bwa_verbose: gates the per-category UGP histograms */
#include "profiling.h"

// Profiling globals. Defined here rather than in main.cpp so they are part of
// libbwa.a and visible to consumers that link the library without pulling in
// main.o (e.g. language bindings).
uint64_t proc_freq = 0;
uint64_t tprof[LIM_R][LIM_C] = {{0}};
uint64_t prof[LIM_R] = {0};

int find_opt(uint64_t *a, int len, uint64_t *max, uint64_t *min, double *avg)
{
    *max = 0;
    *min = 1e15;
    *avg = 0;

    int i=0;
    for (i=0; i<len; i++)
    {
        if (a[i] > *max) *max = a[i];
        if (a[i] < *min) *min = a[i];
        *avg += a[i];
    }
    *avg /= len;

    return 1;
}

int display_stats(int nthreads)
{
    uint64_t max, min;
    double avg;
    fprintf(stderr, "No. of OMP threads: %d\n", nthreads);
    fprintf(stderr, "Processor is running @%lf MHz\n", proc_freq*1.0/1e6);
    fprintf(stderr, "Runtime profile:\n");

    fprintf(stderr, "\n\tTime taken for main_mem function: %0.2lf sec\n\n",
            tprof[MEM][0]*1.0/proc_freq);

    fprintf(stderr, "\tIO times (sec) :\n");
    find_opt(tprof[READ_IO], 1, &max, &min, &avg);
    fprintf(stderr, "\tReading IO time (reads) avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

    find_opt(tprof[SAM_IO], 1, &max, &min, &avg);
    fprintf(stderr, "\tWriting IO time (SAM) avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

    find_opt(tprof[REF_IO], 1, &max, &min, &avg);
    fprintf(stderr, "\tReading IO time (Reference Genome) avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

    find_opt(tprof[FMI], 1, &max, &min, &avg);
    fprintf(stderr, "\tIndex read time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

    fprintf(stderr, "\n\tOverall time (sec) (Excluding Index reading time):\n");
    // find_opt(tprof[PROCESS], 1, &max, &min, &avg);
    fprintf(stderr, "\tPROCESS() (Total compute time + (read + SAM) IO time) : %0.2lf\n",
            tprof[PROCESS][0]*1.0/proc_freq);

    find_opt(tprof[MEM_PROCESS2], 1, &max, &min, &avg);
    fprintf(stderr, "\tMEM_PROCESS_SEQ() (Total compute time (Kernel + SAM)), avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

    fprintf(stderr, "\n\t SAM Processing time (sec):\n");
    find_opt(tprof[WORKER20], 1, &max, &min, &avg);
    fprintf(stderr, "\t--WORKER_SAM avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    
#if HIDE
    find_opt(tprof[SAM1], 1, &max, &min, &avg);
    fprintf(stderr, "\t\tWORKER_SAM1 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    find_opt(tprof[SAM2], 1, &max, &min, &avg);
    fprintf(stderr, "\t\tWORKER_SAM2 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    
    find_opt(tprof[ALIGN1], 1, &max, &min, &avg);
    fprintf(stderr, "\t\t\tWORKER_ALIGN1 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    
    find_opt(tprof[SAM3], 1, &max, &min, &avg);
    fprintf(stderr, "\t\tWORKER_SAM3 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
#endif
    
    fprintf(stderr, "\n\tKernels' compute time (sec):\n");
    find_opt(tprof[WORKER10], 1, &max, &min, &avg);
    fprintf(stderr, "\tTotal kernel (smem+sal+bsw) time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    
#if HIDE
    find_opt(tprof[MEM_ALN_M1], nthreads, &max, &min, &avg);
    fprintf(stderr, "\t\tMEM_ALN_CHAIN_FLT avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    
    find_opt(tprof[MEM_ALN_M2], nthreads, &max, &min, &avg);
    fprintf(stderr, "\t\tMEM_ALN_CHAIN_SEED avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
#endif
    
    find_opt(tprof[MEM_COLLECT], nthreads, &max, &min, &avg);
    fprintf(stderr, "\t\tSMEM compute avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

    find_opt(tprof[MEM_CHAIN], nthreads, &max, &min, &avg);
    fprintf(stderr, "\t\tMEM_CHAIN avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    
    find_opt(tprof[MEM_SA_BLOCK], nthreads, &max, &min, &avg);
    fprintf(stderr, "\t\tSAL compute avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    
    #if 1 //HIDE
    find_opt(tprof[MEM_SA], nthreads, &max, &min, &avg);
    fprintf(stderr, "\t\t\t\tMEM_SA avg: %0.2lf, (%0.2lf, %0.2lf)\n\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    #endif
    
    // printf("\n\t BSW compute time (sec):\n");
    find_opt(tprof[MEM_ALN2], nthreads, &max, &min, &avg);
    fprintf(stderr, "\t\tBSW time, avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

    // PR 26a/c.1/e/26c.2/Q1/Q2/Q3 diagnostics. Gated behind bwa_verbose >= 4
    // — these are histogram-heavy profile-build sections, useful for tuning
    // but noisy in production stderr. Match the convention used elsewhere
    // for high-verbosity diagnostics in bwa.cpp.
    if (bwa_verbose >= 4) {
    // PR 26a/c.1/e diagnostics.
    {
        uint64_t la=0, lh=0, lt=0, ltb1=0, ltb2=0, ltb3=0, ln=0;
        uint64_t ra=0, rh=0, rt=0, rtb1=0, rtb2=0, rtb3=0, rn=0;
        for (int t = 0; t < nthreads; t++) {
            la   += tprof[UGP_L_ATTEMPT][t];
            lh   += tprof[UGP_L_HIT][t];
            lt   += tprof[UGP_L_TIGHT][t];
            ltb1 += tprof[UGP_L_TB_1_8][t];
            ltb2 += tprof[UGP_L_TB_9_32][t];
            ltb3 += tprof[UGP_L_TB_33_MAX][t];
            ln   += tprof[UGP_L_DISP_NARROW][t];
            ra   += tprof[UGP_R_ATTEMPT][t];
            rh   += tprof[UGP_R_HIT][t];
            rt   += tprof[UGP_R_TIGHT][t];
            rtb1 += tprof[UGP_R_TB_1_8][t];
            rtb2 += tprof[UGP_R_TB_9_32][t];
            rtb3 += tprof[UGP_R_TB_33_MAX][t];
            rn   += tprof[UGP_R_DISP_NARROW][t];
        }
        if (la + ra > 0) {
            fprintf(stderr, "\n\tUngapped fast-path / tight_band diagnostics:\n");
            fprintf(stderr, "\t  LEFT   att=%lu  hit=%lu (%.2f%%)  tight=%lu  tb[1..8]=%lu  tb[9..32]=%lu  tb[33..]=%lu  disp_narrow=%lu\n",
                    (unsigned long)la, (unsigned long)lh,
                    la ? 100.0 * (double)lh / (double)la : 0.0,
                    (unsigned long)lt, (unsigned long)ltb1,
                    (unsigned long)ltb2, (unsigned long)ltb3,
                    (unsigned long)ln);
            fprintf(stderr, "\t  RIGHT  att=%lu  hit=%lu (%.2f%%)  tight=%lu  tb[1..8]=%lu  tb[9..32]=%lu  tb[33..]=%lu  disp_narrow=%lu\n",
                    (unsigned long)ra, (unsigned long)rh,
                    ra ? 100.0 * (double)rh / (double)ra : 0.0,
                    (unsigned long)rt, (unsigned long)rtb1,
                    (unsigned long)rtb2, (unsigned long)rtb3,
                    (unsigned long)rn);
        }
    }

    /* PR 26c.2 design instrumentation. */
    {
        const char *fine_lbl[UGP_FINE_NBINS] = {
            "0", "1-2", "3-4", "5-8", "9-16", "17-32", "33-48", "49-64", "65-80", "81-100"
        };
        const char *tier_lbl[2]   = { "8bit ", "16bit" };
        const char *side_lbl[2]   = { "LEFT ", "RIGHT" };
        const char *band_lbl[4]   = { "tb=0", "tb 1-8", "tb 9-32", "tb 33+" };
        const char *narrow_lbl[UGP_NARROW_SZ_NBINS] = {
            "0", "1-15", "16-31", "32-63", "64-127", "128-255", "256-511", "512+"
        };

        /* Group A: fine tight_band histogram per side. */
        for (int side = 0; side < 2; side++) {
            uint64_t total = 0;
            uint64_t bins[UGP_FINE_NBINS] = {0};
            for (int b = 0; b < UGP_FINE_NBINS; b++) {
                for (int t = 0; t < nthreads; t++)
                    bins[b] += tprof[UGP_FINE_BASE + side * UGP_FINE_NBINS + b][t];
                total += bins[b];
            }
            if (total == 0) continue;
            fprintf(stderr, "\t  %s tight_band fine: total=%lu  ",
                    side_lbl[side], (unsigned long)total);
            for (int b = 0; b < UGP_FINE_NBINS; b++)
                fprintf(stderr, "%s=%lu ", fine_lbl[b], (unsigned long)bins[b]);
            fprintf(stderr, "\n");
        }

        /* Group B: tier × band per side. The 2 tiers and 4 bands here mirror
         * the producer-side layout written into UGP_TIER_TB_BASE; keep these
         * named so any future bin/tier-count tweak stays in sync. */
        const int tier_n = 2;
        const int band_n = 4;
        for (int side = 0; side < 2; side++) {
            for (int tier = 0; tier < tier_n; tier++) {
                uint64_t bins[band_n] = {0};
                uint64_t total = 0;
                for (int b = 0; b < band_n; b++) {
                    for (int t = 0; t < nthreads; t++)
                        bins[b] += tprof[UGP_TIER_TB_BASE
                                         + side * (tier_n * band_n)
                                         + tier * band_n + b][t];
                    total += bins[b];
                }
                if (total == 0) continue;
                fprintf(stderr, "\t  %s %s tier:       total=%lu  ",
                        side_lbl[side], tier_lbl[tier], (unsigned long)total);
                for (int b = 0; b < band_n; b++)
                    fprintf(stderr, "%s=%lu ", band_lbl[b], (unsigned long)bins[b]);
                fprintf(stderr, "\n");
            }
        }

        /* Group C: per-batch narrow bucket size histogram per side. */
        for (int side = 0; side < 2; side++) {
            uint64_t bins[UGP_NARROW_SZ_NBINS] = {0};
            uint64_t total = 0;
            for (int b = 0; b < UGP_NARROW_SZ_NBINS; b++) {
                for (int t = 0; t < nthreads; t++)
                    bins[b] += tprof[UGP_NARROW_SZ_BASE + side * UGP_NARROW_SZ_NBINS + b][t];
                total += bins[b];
            }
            if (total == 0) continue;
            fprintf(stderr, "\t  %s narrow bucket sz/batch: total=%lu  ",
                    side_lbl[side], (unsigned long)total);
            for (int b = 0; b < UGP_NARROW_SZ_NBINS; b++)
                fprintf(stderr, "%s=%lu ", narrow_lbl[b], (unsigned long)bins[b]);
            fprintf(stderr, "\n");
        }
    }

    /* Q1+Q2 instrumentation: post-SW ungapped/gapped final outcome and
     * per-side alignment-score histogram. */
    {
        const char *side_lbl[2]  = { "LEFT ", "RIGHT" };
        const char *score_lbl[UGP_SCORE_HIST_NBINS] = {
            "0-10", "11-25", "26-50", "51-75", "76-100", "101-125", "126-150", "151+"
        };
        for (int side = 0; side < 2; side++) {
            uint64_t ung = 0, gap = 0;
            for (int t = 0; t < nthreads; t++) {
                ung += tprof[UGP_OUTCOME_BASE + side * 2 + 0][t];
                gap += tprof[UGP_OUTCOME_BASE + side * 2 + 1][t];
            }
            uint64_t total = ung + gap;
            if (total == 0) continue;
            fprintf(stderr, "\t  %s outcome: total=%lu  ungapped=%lu (%.2f%%)  gapped=%lu (%.2f%%)\n",
                    side_lbl[side], (unsigned long)total,
                    (unsigned long)ung, 100.0 * (double)ung / (double)total,
                    (unsigned long)gap, 100.0 * (double)gap / (double)total);
        }
        for (int side = 0; side < 2; side++) {
            uint64_t bins[UGP_SCORE_HIST_NBINS] = {0};
            uint64_t total = 0;
            for (int b = 0; b < UGP_SCORE_HIST_NBINS; b++) {
                for (int t = 0; t < nthreads; t++)
                    bins[b] += tprof[UGP_SCORE_HIST_BASE + side * UGP_SCORE_HIST_NBINS + b][t];
                total += bins[b];
            }
            if (total == 0) continue;
            fprintf(stderr, "\t  %s score hist:  total=%lu  ",
                    side_lbl[side], (unsigned long)total);
            for (int b = 0; b < UGP_SCORE_HIST_NBINS; b++)
                fprintf(stderr, "%s=%lu ", score_lbl[b], (unsigned long)bins[b]);
            fprintf(stderr, "\n");
        }
    }

    /* Q3 instrumentation: per-category delta-from-perfect histograms (LEFT).
     * delta = (h0 + a*len2) - aln_score; 0 = perfect ungapped extension. */
    {
        const char *cat_lbl[UGP_CAT_NCAT] = {
            "ALL                 ",
            "UNGAP_FINAL         ",
            "GAPPED_FINAL        ",
            "HIT                 ",
            "UNGAP_FINAL_NOT_HIT "
        };
        const char *cat_delta_lbl[UGP_CAT_NBINS] = {
            "0", "1-5", "6-10", "11-25", "26-50", "51-75", "76-100", "101+"
        };
        uint64_t any_ung = 0, any_fin = 0;
        for (int c = 0; c < UGP_CAT_NCAT; c++)
            for (int b = 0; b < UGP_CAT_NBINS; b++)
                for (int t = 0; t < nthreads; t++) {
                    any_ung += tprof[UGP_L_CAT_UNG_BASE + c * UGP_CAT_NBINS + b][t];
                    any_fin += tprof[UGP_L_CAT_FIN_BASE + c * UGP_CAT_NBINS + b][t];
                }
        if (any_ung > 0 || any_fin > 0) {
            fprintf(stderr, "\n\tQ3 LEFT per-category delta-from-perfect histograms:\n");
        }
        if (any_ung > 0) {
            fprintf(stderr, "\t  -- delta on would-be ungapped extension score --\n");
            for (int c = 0; c < UGP_CAT_NCAT; c++) {
                uint64_t bins[UGP_CAT_NBINS] = {0};
                uint64_t total = 0;
                for (int b = 0; b < UGP_CAT_NBINS; b++) {
                    for (int t = 0; t < nthreads; t++)
                        bins[b] += tprof[UGP_L_CAT_UNG_BASE + c * UGP_CAT_NBINS + b][t];
                    total += bins[b];
                }
                fprintf(stderr, "\t  cat%d %s ung total=%lu  ",
                        c, cat_lbl[c], (unsigned long)total);
                for (int b = 0; b < UGP_CAT_NBINS; b++)
                    fprintf(stderr, "%s=%lu ", cat_delta_lbl[b], (unsigned long)bins[b]);
                fprintf(stderr, "\n");
            }
        }
        if (any_fin > 0) {
            fprintf(stderr, "\t  -- delta on final committed alignment score --\n");
            for (int c = 0; c < UGP_CAT_NCAT; c++) {
                uint64_t bins[UGP_CAT_NBINS] = {0};
                uint64_t total = 0;
                for (int b = 0; b < UGP_CAT_NBINS; b++) {
                    for (int t = 0; t < nthreads; t++)
                        bins[b] += tprof[UGP_L_CAT_FIN_BASE + c * UGP_CAT_NBINS + b][t];
                    total += bins[b];
                }
                fprintf(stderr, "\t  cat%d %s fin total=%lu  ",
                        c, cat_lbl[c], (unsigned long)total);
                for (int b = 0; b < UGP_CAT_NBINS; b++)
                    fprintf(stderr, "%s=%lu ", cat_delta_lbl[b], (unsigned long)bins[b]);
                fprintf(stderr, "\n");
            }
        }
    }
    } /* end bwa_verbose >= 4 gate for diagnostic histograms */

    #if HIDE
    int agg1 = 0, agg2 = 0, agg3 = 0;
    for (int i=0; i<nthreads; i++) {
        agg1 += tprof[PE11][i];
        agg2 += tprof[PE12][i];
        agg3 += tprof[PE13][i];
    }
    if (agg1 != agg3) 
        fprintf(stderr, "There is a discrepancy re-allocs, plz rectify!!\n");

    if(agg2 > 0)
    {
        fprintf(stderr, "\n\tTotal re-allocs: %d out of total requests: %d, Rate: %0.2f\n",
                agg1, agg2, agg1*1.0/agg2);
    }

    double res = 0.0, max_ = 0, min_=1e10;
    for (int i=0; i<nthreads; i++) {
        double val = (tprof[ALIGN1][i]*1.0) / tprof[MEM_CHAIN][i];
        res += val;
        if (max_ < val) max_ = val;
        if (min_ > val) min_ = val;
    }
    fprintf(stderr, "\tAvg. FM-index traversal per get_sa_entry(): avg: %lf, max: %lf, min: %lf\n",
            res/nthreads, max_, min_);

    int64_t tot_inst1 = 0, tot_inst2 = 0;
    for (int i=0; i<nthreads; i++) {
        tot_inst1 += tprof[SAM1][i];
        tot_inst2 += tprof[SAM2][i];
    }
    
    fprintf(stderr, "\ttot_inst1: %ld, tot_inst2: %ld, over %d threads\n",
            tot_inst1, tot_inst2, nthreads);
    #endif
    
#if HIDE
    fprintf(stderr, "\n BSW Perf.:\n");
    find_opt(tprof[MEM_ALN2_B], 1, &max, &min, &avg);
    fprintf(stderr, "\tLeft 16-bit time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    find_opt(tprof[MEM_ALN2_D], 1, &max, &min, &avg);
    fprintf(stderr, "\tLeft 8-bit time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    find_opt(tprof[MEM_ALN2_C], 1, &max, &min, &avg);
    fprintf(stderr, "\tRight 16-bit time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    find_opt(tprof[MEM_ALN2_E], 1, &max, &min, &avg);
    fprintf(stderr, "\tRight 8-bit time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
#endif

#if HIDE
    fprintf(stderr, "\nSTATSV\n");
    fprintf(stderr, "%0.2lf\n", tprof[PROCESS][0]*1.0/proc_freq);
    find_opt(tprof[MEM_PROCESS2], 1, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    find_opt(tprof[WORKER10], 1, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    find_opt(tprof[MEM_COLLECT], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    find_opt(tprof[MEM_SA_BLOCK], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    find_opt(tprof[MEM_SA], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    double val = 0;
    find_opt(tprof[MEM_ALN2_UP], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    val += avg;
    find_opt(tprof[CLEFT], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    val += avg;
    find_opt(tprof[CRIGHT], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    val += avg;
    find_opt(tprof[MEM_ALN2_DOWN], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    val += avg;
    fprintf(stderr, "%0.2lf\n", val*1.0/proc_freq);
    fprintf(stderr, "%0.2lf\n", (tprof[REF_IO][0] + tprof[FMI][0])*1.0/proc_freq);
    fprintf(stderr, "%0.2lf\n", tprof[READ_IO][0]*1.0/proc_freq);
#endif
    // printf("\tMemory usage (GB):\n");
    // printf("\tAvg: %0.2lf, Peak: %0.2lf\n", tprof[PE21][0]*1.0/1e9, tprof[PE22][0]*1.0/1e9);

    return 1;
}

