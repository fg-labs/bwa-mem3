/* SPDX-License-Identifier: MIT */
#ifndef BWAMEM3_HTS_PERIODIC_SYNC_H
#define BWAMEM3_HTS_PERIODIC_SYNC_H

/* Spread a BGZF output's fdatasync across the run instead of paying it at close.
 *
 * htslib's bgzf_close() ends with hflush(), which fdatasync()s an output that
 * is a regular file. htslib's multithreaded BGZF writer (--bam-threads > 0)
 * calls hflush() every 512 blocks while it writes, so the syncs are spread
 * through the run on its writer thread. The single-threaded writer, the
 * default, never does, so the whole output is synced at close: ~3.5 GB of
 * uncompressed BAM took ~6 s after the last record was written on a 500 MB/s
 * EBS volume, with every compute thread idle.
 *
 * hts_periodic_sync() gives the single-threaded writer the same 512-block
 * cadence. It runs on the thread that writes records -- the pipeline's write
 * step, which overlaps the next batch's compute -- so close is left with at
 * most one interval to sync. It only calls hflush(), which writes out what is
 * already buffered: the output bytes are unchanged. On a pipe or other
 * non-file output the sync is a no-op (htslib ignores EINVAL from fdatasync).
 *
 * Include after the htslib headers (htslib/sam.h pulls in hts.h and bgzf.h);
 * see bam_writer.h for why they cannot be included from a bwa-mem3 header
 * that callers also see. */

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#include "htslib/hts.h"

/* Compressed bytes between syncs: 512 maximum-size BGZF blocks (32 MiB), the
 * interval htslib's multithreaded writer uses (bgzf.c, "Frequency of 1/512"). */
#define HTS_PERIODIC_SYNC_BYTES ((int64_t)512 * BGZF_MAX_BLOCK_SIZE)

typedef struct {
    int64_t synced_at; /* block_address at the last sync (0 = none yet) */
    int64_t n_syncs;   /* syncs issued, for tests */
    int     off;       /* BWA_MEM3_BAM_PERIODIC_SYNC=0 */
} hts_periodic_sync_t;

/* Read BWA_MEM3_BAM_PERIODIC_SYNC into `s` (zero-initialized = on). `0` turns
 * the periodic sync off, for A/B timing; `1` or an empty value keeps the
 * default. Any other value would read as "on" and silently not do what was
 * asked, so it is reported and the default kept, as BWA_MEM3_STREAM_OUTPUT is.
 * Call once per writer, when it is opened. */
static inline void hts_periodic_sync_init(hts_periodic_sync_t *s)
{
    const char *env = getenv("BWA_MEM3_BAM_PERIODIC_SYNC");
    s->synced_at = 0;
    s->n_syncs = 0;
    s->off = 0;
    if (env == NULL || *env == '\0') return;
    if (strcmp(env, "0") == 0) s->off = 1;
    else if (strcmp(env, "1") != 0)
        fprintf(stderr, "WARNING: BWA_MEM3_BAM_PERIODIC_SYNC='%s' is not 0 or 1; "
                        "ignoring it and syncing periodically\n", env);
}

/* Sync `fp` if HTS_PERIODIC_SYNC_BYTES of BGZF output have been written since
 * the last sync. Call after each record write. A no-op for non-BGZF output
 * (SAM text) and when htslib's multithreaded writer is attached, which syncs on
 * its own. Returns 0, or -1 if the flush or sync failed (e.g. ENOSPC, EIO), which
 * the caller reports as a write error. */
static inline int hts_periodic_sync(hts_periodic_sync_t *s, htsFile *fp)
{
    if (s == NULL || s->off || fp == NULL || !fp->is_bgzf) return 0;
    BGZF *bgzf = fp->fp.bgzf;
    if (bgzf == NULL || bgzf->mt != NULL) return 0;
    if (bgzf->block_address - s->synced_at < HTS_PERIODIC_SYNC_BYTES) return 0;
    s->synced_at = bgzf->block_address;
    ++s->n_syncs;
    return hflush(bgzf->fp) == 0 ? 0 : -1;
}

#endif /* BWAMEM3_HTS_PERIODIC_SYNC_H */
