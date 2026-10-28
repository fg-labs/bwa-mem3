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

#ifndef _PROFILE_HPP
#define _PROFILE_HPP

#include <stdint.h>
#include "macro.h"

/* One tprof row: a per-thread counter for each tid, indexed exactly like the
 * old uint64_t[LIM_C] row (tprof[row][tid]). Each tid's counter sits alone on a
 * 64-byte cache line, so worker threads bumping the same row on a hot path
 * never write to a shared line (eight adjacent tids used to share one). */
#define TPROF_SLOT_U64 (64 / sizeof(uint64_t))
struct tprof_row_t {
    uint64_t slot[LIM_C][TPROF_SLOT_U64];
    uint64_t &operator[](int tid) { return slot[tid][0]; }
    const uint64_t &operator[](int tid) const { return slot[tid][0]; }
};

int display_stats(int );
extern uint64_t proc_freq;

/* proc_freq is the __rdtsc() tick rate (ticks per second) every timing report
 * divides by. It is measured over the process's own run rather than by sleeping
 * at start-up: proc_freq_mark_start() records a (CLOCK_MONOTONIC, __rdtsc())
 * pair, and proc_freq_resolve() sets proc_freq from a second pair the first
 * time a report needs it (display_stats(), the `index` "Total time taken"
 * line), falling back to a short sleep when that window is too short. A
 * nonzero proc_freq set by the caller is kept. Returns proc_freq. */
void proc_freq_mark_start(void);
uint64_t proc_freq_resolve(void);
extern tprof_row_t tprof[LIM_R];
extern uint64_t prof[LIM_R];
#endif
