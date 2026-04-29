/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2026  Fulcrum Genomics, contributors of bwa-mem3.

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
*****************************************************************************************/

#include "bwa_shm.h"
#include "bwa.h"           /* BWA_CTL_SIZE — see static_assert below */
#include "FMI_search.h"    /* CP_FILENAME_SUFFIX, CP_OCC, CP_SHIFT */
#include "macro.h"         /* SA_COMPX */
#include "bntseq.h"        /* bntseq_t, bns_restore, bns_destroy */
#include "utils.h"         /* err_fread_noeof, err_fclose, xopen */
#include "safestringlib.h" /* strcpy_s, strcat_s */

#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <sys/mman.h>      /* shm_open, mmap, munmap */
#include <sys/stat.h>      /* mode constants, stat */
#include <fcntl.h>         /* O_* */
#include <unistd.h>        /* close, ftruncate */
#include <getopt.h>        /* getopt_long for `bwa-mem3 shm --meth ...` */
#include <errno.h>
#include <limits.h>        /* PATH_MAX */

/* /bwactl layout: 16-bit n_entries, 16-bit next_write_offset, then packed
 * (int64_t l_mem, char[] basename) entries. The header is 4 bytes total
 * (two uint16_t fields). next_write_offset starts at 4 on a fresh segment. */
#define BWA_SHM_CTL_HEADER_BYTES 4

/* .bwt.2bit.64 prefix: int64_t reference_seq_len, int64_t count[5]. */
#define BWA_BWT_2BIT_HEADER_BYTES (sizeof(int64_t) * 6)

/* Sanity bound on reference_seq_len read from disk. The bwa-mem2 FM-index
 * has long enforced reference_seq_len <= 0x7fffffffff (about 549 Gbp); see
 * FMI_search::load_index. Anything beyond this is a malformed file. */
#define BWA_SHM_MAX_REF_SEQ_LEN 0x7fffffffffLL

/* Sanity bound on the total packed segment size we'll create. 1 TiB is far
 * larger than any realistic reference genome's index and well below what
 * 64-bit arithmetic can express, so overflow checks against this max are
 * effectively impossible to trip in practice but cheap to enforce. */
#define BWA_SHM_MAX_TOTAL_SIZE (1ull << 40)

/* The two constants must agree. bwa.h::BWA_CTL_SIZE is a v1-era stub that
 * we'll remove in a separate cleanup commit; until then, this assert
 * prevents silent drift if either definition is edited. */
static_assert(BWA_SHM_CTL_SIZE == BWA_CTL_SIZE,
              "BWA_SHM_CTL_SIZE must equal bwa.h::BWA_CTL_SIZE");

/* Return a pointer to the filename component of `hint` — everything after
 * the last '/'. If `hint` has no '/', returns `hint`. */
static const char *prefix_basename(const char *hint)
{
    if (hint == NULL || hint[0] == '\0') return hint;
    const char *p = hint + std::strlen(hint) - 1;
    while (p >= hint && *p != '/') --p;
    return p + 1;
}

/* Build "<a><b>" into `out` (PATH_MAX-bounded) using safestringlib helpers,
 * matching the path-building convention used elsewhere in the codebase
 * (see bntseq.cpp, FMI_search.cpp, fastmap.cpp). */
static void path_concat2(char out[PATH_MAX], const char *a, const char *b)
{
    strcpy_s(out, PATH_MAX, a);
    strcat_s(out, PATH_MAX, b);
}

/* Open `/bwactl` read-write, creating it if necessary and zero-initializing
 * a fresh segment with next-write offset = BWA_SHM_CTL_HEADER_BYTES.
 *
 * On success returns the mmap'd address (a `BWA_SHM_CTL_SIZE`-byte writable
 * mapping) and `*fd_out` holds the still-open file descriptor — caller
 * closes it after the mapping is no longer needed (typical pattern: close
 * immediately, since the mapping survives close).
 *
 * NOTE: there is no advisory lock around the registry. POSIX shm fds on
 * macOS reject flock(2) with EOPNOTSUPP, so a portable lock would require
 * a side-channel lockfile. We accept the v1 race window here: two
 * concurrent `bwa-mem3 shm <prefix>` invocations may both pass the
 * pre-stage `bwa_shm_test` and append duplicate entries (or have one
 * `O_EXCL`-create the segment and the other discover it via EEXIST). The
 * EEXIST branch in `bwa_shm_stage` treats that as already-staged. */
static void *ctl_open_rw(int *fd_out)
{
    int created = 0;
    int fd = shm_open(BWA_SHM_CTL_NAME, O_RDWR, 0);
    if (fd < 0) {
        if (errno != ENOENT) return NULL;
        fd = shm_open(BWA_SHM_CTL_NAME, O_CREAT | O_RDWR | O_EXCL, 0644);
        if (fd < 0) {
            fd = shm_open(BWA_SHM_CTL_NAME, O_RDWR, 0);
            if (fd < 0) return NULL;
        } else {
            created = 1;
        }
    }
    if (ftruncate(fd, BWA_SHM_CTL_SIZE) < 0) { close(fd); return NULL; }
    void *m = mmap(NULL, BWA_SHM_CTL_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (m == MAP_FAILED) { close(fd); return NULL; }
    if (created) {
        std::memset(m, 0, BWA_SHM_CTL_SIZE);
        uint16_t next_write = BWA_SHM_CTL_HEADER_BYTES;
        std::memcpy((uint8_t *)m + 2, &next_write, sizeof(uint16_t));
    }
    *fd_out = fd;
    return m;
}

static void ctl_close(void *m, int fd)
{
    if (m != NULL) munmap(m, BWA_SHM_CTL_SIZE);
    if (fd >= 0)   close(fd);
}

/* Walk the registry, calling `cb(l_mem, name, ctx)` for each entry. cb may
 * return 0 to continue or non-zero to stop. Returns the cb's last return
 * value (or 0 if the registry is absent or empty), -1 on open/mmap error. */
typedef int (*bwa_shm_ctl_cb_t)(int64_t l_mem, const char *name, void *ctx);

static int ctl_walk(bwa_shm_ctl_cb_t cb, void *ctx)
{
    int fd = shm_open(BWA_SHM_CTL_NAME, O_RDONLY, 0);
    if (fd < 0) return errno == ENOENT ? 0 : -1;
    /* A previous failed stage can leave a zero-sized registry; on macOS an
     * unlinked-but-still-referenced segment also lingers at size 0. mmap
     * rejects 0-byte mappings with EINVAL. ctl_open_rw always ftruncates
     * the segment to the full BWA_SHM_CTL_SIZE on creation, so anything
     * smaller is junk left behind by an aborted stage — treat as empty. */
    struct stat st;
    if (fstat(fd, &st) < 0 || (size_t)st.st_size < BWA_SHM_CTL_SIZE) {
        close(fd);
        return 0;
    }
    void *m = mmap(NULL, BWA_SHM_CTL_SIZE, PROT_READ, MAP_SHARED, fd, 0);
    close(fd);
    if (m == MAP_FAILED) return -1;

    const uint8_t *base = (const uint8_t *)m;
    uint16_t n_entries = 0;
    std::memcpy(&n_entries, base, sizeof(uint16_t));
    const uint8_t *p   = base + BWA_SHM_CTL_HEADER_BYTES;
    const uint8_t *end = base + BWA_SHM_CTL_SIZE;
    int rc = 0;
    for (uint16_t i = 0; i < n_entries; ++i) {
        if (p + sizeof(int64_t) >= end) break;
        int64_t l_mem = 0;
        std::memcpy(&l_mem, p, sizeof(int64_t));
        const char *name = (const char *)(p + sizeof(int64_t));
        size_t name_max  = (size_t)(end - (const uint8_t *)name);
        size_t namelen   = strnlen(name, name_max);
        if (namelen == name_max) break;
        rc = cb(l_mem, name, ctx);
        if (rc != 0) break;
        p = (const uint8_t *)name + namelen + 1;
    }

    munmap(m, BWA_SHM_CTL_SIZE);
    return rc;
}

extern "C" {

int bwa_shm_compute(const char *prefix, bwa_shm_layout_t *layout)
{
    if (prefix == NULL || layout == NULL) return -1;

    memset(layout, 0, sizeof(*layout));
    size_t plen = strlen(prefix);
    if (plen + 1 > sizeof(layout->prefix)) {
        fprintf(stderr, "[E::%s] prefix too long\n", __func__);
        return -1;
    }
    memcpy(layout->prefix, prefix, plen + 1);

    /* 1. Load BNS (.ann + .amb + name/anno strings; opens .pac fp which we
     *    don't need here). */
    bntseq_t *bns = bns_restore(prefix);
    if (bns == NULL) {
        fprintf(stderr, "[E::%s] bns_restore(%s) failed\n", __func__, prefix);
        return -1;
    }
    /* Close the .pac stream that bns_restore left open; pack_into re-opens
     * it via xopen on the .pac path. Zero the field so bns_destroy doesn't
     * double-close. */
    if (bns->fp_pac != NULL) {
        err_fclose(bns->fp_pac);
        bns->fp_pac = NULL;
    }
    layout->bns = bns;

    if (bns->n_seqs < 0 || bns->n_holes < 0 || bns->l_pac <= 0) {
        fprintf(stderr,
            "[E::%s] malformed BNS in %s (n_seqs=%d n_holes=%d l_pac=%lld)\n",
            __func__, prefix, bns->n_seqs, bns->n_holes, (long long)bns->l_pac);
        bwa_shm_layout_free(layout);
        return -1;
    }

    /* 2. Peek FMI scalars + sentinel_index out of <prefix>.bwt.2bit.64. */
    char cp_path[PATH_MAX];
    path_concat2(cp_path, prefix, CP_FILENAME_SUFFIX);
    FILE *cp = xopen(cp_path, "rb");
    err_fread_noeof(&layout->reference_seq_len, sizeof(int64_t), 1, cp);
    err_fread_noeof(layout->count, sizeof(int64_t), 5, cp);

    if (layout->reference_seq_len <= 0 ||
        layout->reference_seq_len > BWA_SHM_MAX_REF_SEQ_LEN) {
        fprintf(stderr,
            "[E::%s] %s: reference_seq_len=%lld out of bounds (1..%lld)\n",
            __func__, cp_path, (long long)layout->reference_seq_len,
            (long long)BWA_SHM_MAX_REF_SEQ_LEN);
        err_fclose(cp);
        bwa_shm_layout_free(layout);
        return -1;
    }
    /* count[i] read from disk is the cumulative occurrence; +1 mirrors the
     * adjustment in FMI_search::load_index. Each value must be in
     * [0, reference_seq_len] after adjustment. */
    for (int i = 0; i < 5; ++i) {
        if (layout->count[i] < 0 || layout->count[i] > layout->reference_seq_len) {
            fprintf(stderr,
                "[E::%s] %s: count[%d]=%lld out of bounds (0..%lld)\n",
                __func__, cp_path, i, (long long)layout->count[i],
                (long long)layout->reference_seq_len);
            err_fclose(cp);
            bwa_shm_layout_free(layout);
            return -1;
        }
        layout->count[i] += 1;
    }

    int64_t cp_occ_count = (layout->reference_seq_len >> CP_SHIFT) + 1;
    int64_t sa_count     = (layout->reference_seq_len >> SA_COMPX) + 1;
    int64_t cp_occ_bytes = cp_occ_count * (int64_t)sizeof(CP_OCC);
    int64_t sa_ms_bytes  = sa_count     * (int64_t)sizeof(int8_t);
    int64_t sa_ls_bytes  = sa_count     * (int64_t)sizeof(uint32_t);

    int64_t sentinel_off = (int64_t)BWA_BWT_2BIT_HEADER_BYTES
                         + cp_occ_bytes + sa_ms_bytes + sa_ls_bytes;
    if (fseek(cp, (long)sentinel_off, SEEK_SET) != 0) {
        fprintf(stderr, "[E::%s] fseek(%lld) failed in %s\n",
                __func__, (long long)sentinel_off, cp_path);
        err_fclose(cp);
        bwa_shm_layout_free(layout);
        return -1;
    }
    err_fread_noeof(&layout->sentinel_index, sizeof(int64_t), 1, cp);
    err_fclose(cp);

    if (layout->sentinel_index < 0 ||
        layout->sentinel_index >= layout->reference_seq_len) {
        fprintf(stderr,
            "[E::%s] %s: sentinel_index=%lld out of bounds (0..%lld)\n",
            __func__, cp_path, (long long)layout->sentinel_index,
            (long long)layout->reference_seq_len);
        bwa_shm_layout_free(layout);
        return -1;
    }

    /* 3. Stat <prefix>.0123 for its size. */
    char zer_path[PATH_MAX];
    path_concat2(zer_path, prefix, ".0123");
    struct stat zst;
    if (stat(zer_path, &zst) != 0) {
        fprintf(stderr, "[E::%s] stat(%s) failed: %s\n",
                __func__, zer_path, strerror(errno));
        bwa_shm_layout_free(layout);
        return -1;
    }
    if (zst.st_size <= 0) {
        fprintf(stderr, "[E::%s] %s is empty\n", __func__, zer_path);
        bwa_shm_layout_free(layout);
        return -1;
    }
    layout->ref_string_len = (int64_t)zst.st_size;

    /* 4. BNS_NAMES section size: walk anns[]. */
    int64_t bns_names_bytes = 0;
    for (int i = 0; i < bns->n_seqs; ++i) {
        bns_names_bytes += (int64_t)strlen(bns->anns[i].name) + 1
                         + (int64_t)strlen(bns->anns[i].anno) + 1;
    }

    /* 5. Compute section sizes (matching the historical layout). */
    int64_t fmi_scalars_bytes = (int64_t)BWA_SHM_FMI_SCALARS_BYTES;
    int64_t bns_struct_bytes  = (int64_t)sizeof(bntseq_t);
    int64_t bns_ambs_bytes    = (int64_t)bns->n_holes * (int64_t)sizeof(bntamb1_t);
    int64_t bns_anns_bytes    = (int64_t)bns->n_seqs  * (int64_t)sizeof(bntann1_t);
    int64_t pac_bytes         = bns->l_pac / 4 + 1;

    /* 6. Lay out sections at 64-byte aligned offsets. */
    static const uint32_t N_SECTIONS = 10;
    uint64_t cursor = sizeof(bwa_shm_header_t)
                    + (uint64_t)N_SECTIONS * sizeof(bwa_shm_section_t);
    cursor = (cursor + 63ull) & ~63ull;

    bwa_shm_section_t *sec = layout->sections;
    auto place = [&](uint32_t kind, int64_t size, uint32_t i) {
        sec[i].kind   = kind;
        sec[i].offset = cursor;
        sec[i].size   = (uint64_t)size;
        cursor = (cursor + (uint64_t)size + 63ull) & ~63ull;
    };

    place(BWA_SHM_SEC_FMI_SCALARS, fmi_scalars_bytes, 0);
    place(BWA_SHM_SEC_FMI_CP_OCC,  cp_occ_bytes,      1);
    place(BWA_SHM_SEC_FMI_SA_MS,   sa_ms_bytes,       2);
    place(BWA_SHM_SEC_FMI_SA_LS,   sa_ls_bytes,       3);
    place(BWA_SHM_SEC_BNS_STRUCT,  bns_struct_bytes,  4);
    place(BWA_SHM_SEC_BNS_AMBS,    bns_ambs_bytes,    5);
    place(BWA_SHM_SEC_BNS_ANNS,    bns_anns_bytes,    6);
    place(BWA_SHM_SEC_BNS_NAMES,   bns_names_bytes,   7);
    place(BWA_SHM_SEC_PAC,         pac_bytes,         8);
    place(BWA_SHM_SEC_REF_STRING,  layout->ref_string_len, 9);

    if (cursor > BWA_SHM_MAX_TOTAL_SIZE) {
        fprintf(stderr,
            "[E::%s] computed segment size %llu exceeds max %llu\n",
            __func__, (unsigned long long)cursor,
            (unsigned long long)BWA_SHM_MAX_TOTAL_SIZE);
        bwa_shm_layout_free(layout);
        return -1;
    }

    layout->n_sections = N_SECTIONS;
    layout->total_size = cursor;
    return 0;
}

int bwa_shm_pack_into(const bwa_shm_layout_t *layout, uint8_t *dest)
{
    if (layout == NULL || dest == NULL || layout->bns == NULL) return -1;
    const bntseq_t *bns = (const bntseq_t *)layout->bns;
    const bwa_shm_section_t *sec = layout->sections;

    /* 1. Header. */
    bwa_shm_header_t hdr;
    memset(&hdr, 0, sizeof(hdr));
    hdr.magic       = BWA_SHM_MAGIC;
    hdr.version     = BWA_SHM_VERSION;
    hdr.n_sections  = layout->n_sections;
    hdr.total_size  = layout->total_size;
    memcpy(dest, &hdr, sizeof(hdr));

    /* 2. Section table. */
    memcpy(dest + sizeof(bwa_shm_header_t),
           layout->sections,
           (size_t)layout->n_sections * sizeof(bwa_shm_section_t));

    /* 3. FMI_SCALARS (assembled on stack). */
    {
        uint8_t scratch[BWA_SHM_FMI_SCALARS_BYTES];
        memcpy(scratch,        &layout->reference_seq_len, sizeof(int64_t));
        memcpy(scratch + 8,    layout->count,              sizeof(int64_t) * 5);
        memcpy(scratch + 48,   &layout->sentinel_index,    sizeof(int64_t));
        memcpy(dest + sec[0].offset, scratch, sizeof(scratch));
    }

    /* 4. CP_OCC + SA_MS + SA_LS: one fp, three contiguous freads. */
    char cp_path[PATH_MAX];
    path_concat2(cp_path, layout->prefix, CP_FILENAME_SUFFIX);
    FILE *cp = xopen(cp_path, "rb");
    if (fseek(cp, (long)BWA_BWT_2BIT_HEADER_BYTES, SEEK_SET) != 0) {
        fprintf(stderr, "[E::%s] fseek on %s failed\n", __func__, cp_path);
        err_fclose(cp);
        return -1;
    }
    err_fread_noeof(dest + sec[1].offset, 1, (size_t)sec[1].size, cp);
    err_fread_noeof(dest + sec[2].offset, 1, (size_t)sec[2].size, cp);
    err_fread_noeof(dest + sec[3].offset, 1, (size_t)sec[3].size, cp);
    err_fclose(cp);

    /* 5. BNS_STRUCT, BNS_AMBS, BNS_ANNS — all small. */
    memcpy(dest + sec[4].offset, bns,       (size_t)sec[4].size);
    memcpy(dest + sec[5].offset, bns->ambs, (size_t)sec[5].size);
    memcpy(dest + sec[6].offset, bns->anns, (size_t)sec[6].size);

    /* 6. BNS_NAMES: concatenate name\0 anno\0 for each seq. The bns_anns
     *    payload's name/anno pointers are not patched here — the unpacker
     *    points them into BNS_NAMES at attach time once the segment base
     *    address is known. */
    {
        uint8_t *dst = dest + sec[7].offset;
        for (int i = 0; i < bns->n_seqs; ++i) {
            size_t n = strlen(bns->anns[i].name) + 1;
            memcpy(dst, bns->anns[i].name, n); dst += n;
            size_t a = strlen(bns->anns[i].anno) + 1;
            memcpy(dst, bns->anns[i].anno, a); dst += a;
        }
    }

    /* 7. PAC: stream <prefix>.pac. */
    {
        char pac_path[PATH_MAX];
        path_concat2(pac_path, layout->prefix, ".pac");
        FILE *fp = xopen(pac_path, "rb");
        err_fread_noeof(dest + sec[8].offset, 1, (size_t)sec[8].size, fp);
        err_fclose(fp);
    }

    /* 8. REF_STRING: stream <prefix>.0123. */
    {
        char zer_path[PATH_MAX];
        path_concat2(zer_path, layout->prefix, ".0123");
        FILE *fp = xopen(zer_path, "rb");
        err_fread_noeof(dest + sec[9].offset, 1, (size_t)sec[9].size, fp);
        err_fclose(fp);
    }

    return 0;
}

void bwa_shm_layout_free(bwa_shm_layout_t *layout)
{
    if (layout == NULL) return;
    if (layout->bns != NULL) {
        bns_destroy((bntseq_t *)layout->bns);
        layout->bns = NULL;
    }
}

int bwa_shm_pack_from_disk(const char *prefix, uint8_t **buf_out, uint64_t *len_out)
{
    *buf_out = NULL;
    *len_out = 0;

    bwa_shm_layout_t layout;
    if (bwa_shm_compute(prefix, &layout) != 0) return -1;

    uint8_t *buf = (uint8_t *)malloc(layout.total_size);
    if (buf == NULL) {
        bwa_shm_layout_free(&layout);
        return -1;
    }
    if (bwa_shm_pack_into(&layout, buf) != 0) {
        free(buf);
        bwa_shm_layout_free(&layout);
        return -1;
    }

    *buf_out = buf;
    *len_out = layout.total_size;
    bwa_shm_layout_free(&layout);
    return 0;
}

static int test_cb(int64_t /*l_mem*/, const char *name, void *ctx)
{
    return std::strcmp(name, (const char *)ctx) == 0 ? 1 : 0;
}

int bwa_shm_test(const char *prefix)
{
    if (prefix == NULL || prefix[0] == '\0') return -1;
    return ctl_walk(test_cb, (void *)prefix_basename(prefix));
}

static int list_cb(int64_t l_mem, const char *name, void * /*ctx*/)
{
    std::printf("%s\t%lld\n", name, (long long)l_mem);
    return 0;
}

int bwa_shm_list(void)
{
    int rc = ctl_walk(list_cb, NULL);
    return rc < 0 ? rc : 0;
}

static int destroy_cb(int64_t /*l_mem*/, const char *name, void * /*ctx*/)
{
    char idx_path[PATH_MAX];
    path_concat2(idx_path, BWA_SHM_IDX_PREFIX, name);
    shm_unlink(idx_path);   /* ignore errors (entry may already be gone) */
    return 0;
}

int bwa_shm_destroy(void)
{
    int rc = ctl_walk(destroy_cb, NULL);
    if (rc < 0) return rc;
    shm_unlink(BWA_SHM_CTL_NAME);
    return 0;
}

int bwa_shm_stage(const char *prefix)
{
    if (prefix == NULL || prefix[0] == '\0') {
        std::fprintf(stderr, "[E::%s] empty prefix\n", __func__);
        return -1;
    }
    int already = bwa_shm_test(prefix);
    if (already < 0) {
        std::fprintf(stderr, "[E::%s] failed to probe registry\n", __func__);
        return -1;
    }
    if (already == 1) {
        std::fprintf(stderr, "[M::%s] index '%s' is already in shm\n",
                     __func__, prefix);
        return 0;
    }

    const char *name = prefix_basename(prefix);
    if (name[0] == '\0') {
        std::fprintf(stderr, "[E::%s] prefix has no basename component\n", __func__);
        return -1;
    }

    /* 1. Compute the layout (loads BNS; peeks scalars from the FMI file). */
    bwa_shm_layout_t layout;
    if (bwa_shm_compute(prefix, &layout) != 0) {
        std::fprintf(stderr, "[E::%s] failed to compute layout for '%s'\n",
                     __func__, prefix);
        return -1;
    }

    /* 2. Open the control segment R/W. There is no advisory lock — see the
     *    note on ctl_open_rw for why. We re-check the registry below to
     *    catch the case where another stager appended this prefix between
     *    our pre-stage bwa_shm_test and now. */
    int   ctl_fd  = -1;
    void *ctl_map = ctl_open_rw(&ctl_fd);
    if (ctl_map == NULL) {
        std::fprintf(stderr, "[E::%s] failed to open control segment %s: %s\n",
                     __func__, BWA_SHM_CTL_NAME, std::strerror(errno));
        bwa_shm_layout_free(&layout);
        return -1;
    }

    /* 3. Re-scan the live registry for this basename. */
    uint16_t n_entries = 0, next_write = 0;
    std::memcpy(&n_entries,  (uint8_t *)ctl_map,     sizeof(uint16_t));
    std::memcpy(&next_write, (uint8_t *)ctl_map + 2, sizeof(uint16_t));
    {
        const uint8_t *p   = (const uint8_t *)ctl_map + BWA_SHM_CTL_HEADER_BYTES;
        const uint8_t *end = (const uint8_t *)ctl_map + BWA_SHM_CTL_SIZE;
        for (uint16_t i = 0; i < n_entries; ++i) {
            if (p + sizeof(int64_t) >= end) break;
            const char *entry_name = (const char *)(p + sizeof(int64_t));
            size_t name_max = (size_t)(end - (const uint8_t *)entry_name);
            size_t namelen  = strnlen(entry_name, name_max);
            if (namelen == name_max) break;
            if (std::strcmp(entry_name, name) == 0) {
                std::fprintf(stderr,
                    "[M::%s] index '%s' was staged by another process\n",
                    __func__, prefix);
                ctl_close(ctl_map, ctl_fd);
                bwa_shm_layout_free(&layout);
                return 0;
            }
            p = (const uint8_t *)entry_name + namelen + 1;
        }
    }

    /* 4. Bound-check the new entry. next_write is a uint16_t and
     *    BWA_SHM_CTL_SIZE is 0x10000 (65536), which uint16_t cannot represent.
     *    Reject the exact-fill case so the post-write cast back to uint16_t
     *    cannot wrap to 0 and have the next stage overwrite the header. */
    size_t name_bytes  = std::strlen(name) + 1;
    size_t entry_bytes = sizeof(int64_t) + name_bytes;
    if ((size_t)next_write + entry_bytes >= BWA_SHM_CTL_SIZE) {
        std::fprintf(stderr,
                     "[E::%s] control segment full (cannot stage '%s')\n",
                     __func__, name);
        ctl_close(ctl_map, ctl_fd);
        bwa_shm_layout_free(&layout);
        return -1;
    }

    /* 5. Create the per-index segment. EEXIST means another stager raced us
     *    to a segment with the same basename — treat that as already-staged
     *    rather than unlinking and re-creating, which would invalidate any
     *    process currently attached to it. */
    char idx_path[PATH_MAX];
    path_concat2(idx_path, BWA_SHM_IDX_PREFIX, name);
    int idx_fd = shm_open(idx_path, O_CREAT | O_RDWR | O_EXCL, 0644);
    if (idx_fd < 0 && errno == EEXIST) {
        std::fprintf(stderr,
            "[M::%s] segment %s already exists; treating as already staged\n",
            __func__, idx_path);
        ctl_close(ctl_map, ctl_fd);
        bwa_shm_layout_free(&layout);
        return 0;
    }
    if (idx_fd < 0) {
        std::fprintf(stderr, "[E::%s] shm_open(%s) failed: %s\n",
                     __func__, idx_path, std::strerror(errno));
        ctl_close(ctl_map, ctl_fd);
        bwa_shm_layout_free(&layout);
        return -1;
    }
    if (ftruncate(idx_fd, (off_t)layout.total_size) < 0) {
        std::fprintf(stderr, "[E::%s] ftruncate(%s, %llu) failed: %s\n",
                     __func__, idx_path,
                     (unsigned long long)layout.total_size, std::strerror(errno));
        close(idx_fd);
        shm_unlink(idx_path);
        ctl_close(ctl_map, ctl_fd);
        bwa_shm_layout_free(&layout);
        return -1;
    }
    void *idx_map = mmap(NULL, (size_t)layout.total_size,
                         PROT_READ | PROT_WRITE, MAP_SHARED, idx_fd, 0);
    if (idx_map == MAP_FAILED) {
        std::fprintf(stderr, "[E::%s] mmap(%s) failed: %s\n",
                     __func__, idx_path, std::strerror(errno));
        close(idx_fd);
        shm_unlink(idx_path);
        ctl_close(ctl_map, ctl_fd);
        bwa_shm_layout_free(&layout);
        return -1;
    }

    /* 6. Stream-pack directly into the mmap'd region. No heap intermediate. */
    if (bwa_shm_pack_into(&layout, (uint8_t *)idx_map) != 0) {
        std::fprintf(stderr, "[E::%s] pack_into failed for %s\n", __func__, prefix);
        munmap(idx_map, (size_t)layout.total_size);
        close(idx_fd);
        shm_unlink(idx_path);
        ctl_close(ctl_map, ctl_fd);
        bwa_shm_layout_free(&layout);
        return -1;
    }
    munmap(idx_map, (size_t)layout.total_size);
    close(idx_fd);

    /* 7. Append the entry to the control segment. */
    uint8_t *entry = (uint8_t *)ctl_map + next_write;
    int64_t l_mem = (int64_t)layout.total_size;
    std::memcpy(entry,                       &l_mem, sizeof(int64_t));
    std::memcpy(entry + sizeof(int64_t),     name,   name_bytes);
    n_entries  = (uint16_t)(n_entries + 1);
    next_write = (uint16_t)(next_write + entry_bytes);
    std::memcpy((uint8_t *)ctl_map,     &n_entries,  sizeof(uint16_t));
    std::memcpy((uint8_t *)ctl_map + 2, &next_write, sizeof(uint16_t));

    ctl_close(ctl_map, ctl_fd);

    std::fprintf(stderr, "[M::%s] staged '%s' (%llu bytes) as %s\n",
                 __func__, name,
                 (unsigned long long)layout.total_size, idx_path);
    bwa_shm_layout_free(&layout);
    return 0;
}

uint8_t *bwa_shm_attach(const char *prefix, size_t *len_out)
{
    if (len_out != NULL) *len_out = 0;
    if (prefix == NULL || prefix[0] == '\0') return NULL;

    /* Cheap pre-check: skip the open if there's no registry entry. */
    if (bwa_shm_test(prefix) != 1) return NULL;

    const char *name = prefix_basename(prefix);
    char idx_path[PATH_MAX];
    path_concat2(idx_path, BWA_SHM_IDX_PREFIX, name);

    int fd = shm_open(idx_path, O_RDONLY, 0);
    if (fd < 0) {
        std::fprintf(stderr, "[E::%s] shm_open(%s) failed: %s\n",
                     __func__, idx_path, std::strerror(errno));
        return NULL;
    }

    struct stat st;
    if (fstat(fd, &st) < 0) {
        std::fprintf(stderr, "[E::%s] fstat(%s) failed: %s\n",
                     __func__, idx_path, std::strerror(errno));
        close(fd);
        return NULL;
    }
    if (st.st_size < (off_t)sizeof(bwa_shm_header_t)) {
        std::fprintf(stderr, "[E::%s] segment %s too small (%lld bytes)\n",
                     __func__, idx_path, (long long)st.st_size);
        close(fd);
        return NULL;
    }
    size_t len = (size_t)st.st_size;

    void *m = mmap(NULL, len, PROT_READ, MAP_SHARED, fd, 0);
    close(fd);
    if (m == MAP_FAILED) {
        std::fprintf(stderr, "[E::%s] mmap(%s, %llu) failed: %s\n",
                     __func__, idx_path, (unsigned long long)len, std::strerror(errno));
        return NULL;
    }

    /* Validate the on-segment header. */
    bwa_shm_header_t hdr;
    std::memcpy(&hdr, m, sizeof(hdr));
    if (hdr.magic != BWA_SHM_MAGIC) {
        std::fprintf(stderr, "[E::%s] bad magic in %s\n", __func__, idx_path);
        munmap(m, len);
        return NULL;
    }
    if (hdr.version != BWA_SHM_VERSION) {
        std::fprintf(stderr, "[E::%s] unsupported segment version %u in %s (expected %u)\n",
                     __func__, hdr.version, idx_path, BWA_SHM_VERSION);
        munmap(m, len);
        return NULL;
    }
    /* total_size is the exact data written; st_size may be page-rounded up.
     * Require total_size <= segment size to tolerate ftruncate's rounding. */
    if (hdr.total_size > (uint64_t)len) {
        std::fprintf(stderr, "[E::%s] header total_size=%llu exceeds segment size %llu in %s\n",
                     __func__, (unsigned long long)hdr.total_size,
                     (unsigned long long)len, idx_path);
        munmap(m, len);
        return NULL;
    }

    if (len_out != NULL) *len_out = (size_t)hdr.total_size;
    return (uint8_t *)m;
}

int bwa_shm_section_find(const uint8_t *base, uint32_t kind,
                         uint64_t *offset_out, uint64_t *size_out)
{
    if (base == NULL) return -1;
    const bwa_shm_header_t *hdr = (const bwa_shm_header_t *)base;
    if (hdr->magic   != BWA_SHM_MAGIC)   return -1;
    if (hdr->version != BWA_SHM_VERSION) return -1;

    /* Bound-check the section table itself before walking it. A corrupt
     * n_sections could otherwise drive us past the mapped region. */
    if (hdr->total_size < sizeof(bwa_shm_header_t)) return -1;
    const uint64_t table_bytes =
        (uint64_t)hdr->n_sections * (uint64_t)sizeof(bwa_shm_section_t);
    if (table_bytes > hdr->total_size - sizeof(bwa_shm_header_t)) return -1;

    const bwa_shm_section_t *table =
        (const bwa_shm_section_t *)(base + sizeof(bwa_shm_header_t));
    for (uint32_t i = 0; i < hdr->n_sections; ++i) {
        if (table[i].kind == kind) {
            /* Reject sections whose extent escapes the segment. */
            if (table[i].size > hdr->total_size ||
                table[i].offset > hdr->total_size - table[i].size) {
                return -1;
            }
            if (offset_out != NULL) *offset_out = table[i].offset;
            if (size_out   != NULL) *size_out   = table[i].size;
            return 0;
        }
    }
    return -1;
}

/* Resolve the on-disk prefix for `bwa-mem3 shm --meth <idxbase>` so that
 * the registry basename matches what `bwa-mem3 mem --meth <idxbase>` will
 * later look up. Mirrors fastmap.cpp's c2t-suffix resolution: append
 * `.bwameth.c2t` unless `prefix` already ends with it. */
static int resolve_meth_prefix(const char *prefix, char out[PATH_MAX])
{
    static const char SUFFIX[] = ".bwameth.c2t";
    static const size_t SUFLEN = sizeof(SUFFIX) - 1;
    size_t plen = std::strlen(prefix);
    int already = (plen >= SUFLEN) &&
                  (std::strcmp(prefix + plen - SUFLEN, SUFFIX) == 0);
    if (already) {
        if (plen + 1 > (size_t)PATH_MAX) {
            std::fprintf(stderr, "[E::%s] prefix too long\n", __func__);
            return -1;
        }
        strcpy_s(out, PATH_MAX, prefix);
    } else {
        if (plen + SUFLEN + 1 > (size_t)PATH_MAX) {
            std::fprintf(stderr, "[E::%s] prefix too long for --meth\n", __func__);
            return -1;
        }
        path_concat2(out, prefix, SUFFIX);
    }
    return 0;
}

/* If `<prefix>.0123` is absent but `<prefix>.bwameth.c2t.0123` is present,
 * the user likely meant `--meth`. Print a hint and return 1. Returns 0
 * when no hint applies (let the regular stage path produce its own error). */
static int hint_missing_meth_flag(const char *prefix)
{
    char zer_path[PATH_MAX];
    path_concat2(zer_path, prefix, ".0123");
    struct stat st;
    if (stat(zer_path, &st) == 0) return 0;        /* literal index present */
    if (errno != ENOENT) return 0;                 /* I/O error: defer */

    char c2t_zer[PATH_MAX];
    path_concat2(c2t_zer, prefix, ".bwameth.c2t.0123");
    if (stat(c2t_zer, &st) != 0) return 0;         /* no meth index either */

    std::fprintf(stderr,
        "[E::main_shm] no FMI at '%s.0123' but a meth index is present at "
        "'%s.bwameth.c2t.*'\n"
        "             did you mean `bwa-mem3 shm --meth %s`?\n",
        prefix, prefix, prefix);
    return 1;
}

static void print_shm_usage(void)
{
    std::fprintf(stderr,
        "\nUsage: bwa-mem3 shm [-d|-l|--help] [--meth] [idxbase]\n\n"
        "Options:\n"
        "  -d        destroy all indices in shared memory (matches bwa v1 behavior)\n"
        "  -l        list names of indices in shared memory\n"
        "  --meth    stage a `bwa-mem3 index --meth` index — auto-appends\n"
        "            `.bwameth.c2t` to <idxbase>, mirroring `mem --meth`\n"
        "  -h --help print this help and exit\n\n"
        "Stage with no flags: `bwa-mem3 shm <idxbase>` loads the index into\n"
        "POSIX shared memory; subsequent `bwa-mem3 mem <idxbase> ...` runs\n"
        "auto-attach instead of re-reading from disk. For meth indices, pass\n"
        "the same plain `<idxbase>` to all three commands plus `--meth` on\n"
        "`index`, `shm`, and `mem` (the c2t suffix is auto-appended).\n\n"
        "Footgun: if you re-build the index, run `bwa-mem3 shm -d` first.\n"
        "There is no staleness check -- a stale segment will silently mis-align.\n\n"
        "macOS: POSIX shm has implementation-defined per-segment caps; large\n"
        "       indices may simply fail to stage. Prefer Linux for production.\n"
        "Linux: /dev/shm defaults to ~50%% of RAM on bare metal; in containers\n"
        "       it is often much smaller and may need raising via --shm-size\n"
        "       (Docker) or an emptyDir tmpfs (Kubernetes).\n\n");
}

int main_shm(int argc, char *argv[])
{
    int c, to_list = 0, to_drop = 0, meth = 0, ret = 0;

    /* Pre-scan for --help / -h. getopt_long would treat -h as unknown
     * (no short opt declared) and print a generic error, but the top-level
     * usage advertises `bwa-mem3 <command> --help`, so we honor both. */
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--help") == 0 || std::strcmp(argv[i], "-h") == 0) {
            print_shm_usage();
            return 0;
        }
    }

    /* getopt's global state may carry over from prior subcommand parsing.
     * Reset for a clean reparse. POSIX requires optind=1; macOS/BSD also
     * needs optreset=1 to clear the static `place` pointer inside getopt. */
    optind = 1;
    opterr = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    extern int optreset;
    optreset = 1;
#endif
    static struct option long_opts[] = {
        {"meth", no_argument, 0, 1000},
        {0,      0,           0, 0   }
    };
    while ((c = getopt_long(argc, argv, "ld", long_opts, NULL)) >= 0) {
        if      (c == 'l')   to_list = 1;
        else if (c == 'd')   to_drop = 1;
        else if (c == 1000)  meth    = 1;
        else                 return 1;   /* getopt printed the error */
    }
    if (optind == argc && !to_list && !to_drop) {
        print_shm_usage();
        return 1;
    }
    if (optind < argc && (to_list || to_drop)) {
        std::fprintf(stderr,
            "[E::%s] -l or -d cannot be combined with idxbase\n", __func__);
        return 1;
    }
    if (meth && (to_list || to_drop)) {
        std::fprintf(stderr,
            "[E::%s] --meth cannot be combined with -l or -d\n", __func__);
        return 1;
    }
    if (optind + 1 < argc) {
        std::fprintf(stderr,
            "[E::%s] too many positional arguments; expected at most one idxbase\n",
            __func__);
        return 1;
    }

    /* Stage. */
    if (optind < argc) {
        const char *user_prefix = argv[optind];
        char       resolved[PATH_MAX];
        if (meth) {
            if (resolve_meth_prefix(user_prefix, resolved) != 0) return 1;
            user_prefix = resolved;
        } else if (hint_missing_meth_flag(user_prefix)) {
            return 1;
        }
        if (bwa_shm_stage(user_prefix) < 0) {
            std::fprintf(stderr,
                "[E::%s] failed to stage '%s' in shared memory\n",
                __func__, user_prefix);
            ret = 1;
        }
    }
    if (to_list)  ret |= (bwa_shm_list()    < 0) ? 1 : 0;
    if (to_drop)  ret |= (bwa_shm_destroy() < 0) ? 1 : 0;
    return ret;
}

} /* extern "C" */
