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

#include "read_index_ele.h"
#include "safestringlib.h"

#include "bwa_madvise.h"
#include "bwa_shm.h"

#include <cstring>     /* memcpy */

indexEle::indexEle()
{
    idx = (bwaidx_fm_t*) calloc(1, sizeof(bwaidx_fm_t));
    assert(idx != NULL);
}

indexEle::~indexEle()
{
    if (idx == 0) return;
    if (idx->mem == 0)
    {
        if (idx->bns) bns_destroy(idx->bns);
        if (idx->pac) free(idx->pac);
    } else {
        free(idx->bns->anns); free(idx->bns);
        if (!idx->is_shm) free(idx->mem);
    }
    free(idx);
}

void indexEle::bwa_idx_load_ele(const char *hint, int which)
{
    char *prefix;
    int l_hint = strlen(hint);
    prefix = (char *) malloc(l_hint + 3 + 4 + 1);
    assert(prefix != NULL);
    strcpy_s(prefix, l_hint + 3 + 4 + 1, hint);

    fprintf(stderr, "* Index prefix: %s\n", prefix);
    
    // idx = (bwaidx_fm_t*) calloc(1, sizeof(bwaidx_fm_t));
    if (which & BWA_IDX_BNS) {
        int i, c;
        idx->bns = bns_restore(prefix);
        if (idx->bns == 0) {
            printf("Error!! : [%s] bns is NULL!!\n", __func__);
            exit(EXIT_FAILURE);
        }
        for (i = c = 0; i < idx->bns->n_seqs; ++i)
            if (idx->bns->anns[i].is_alt) ++c;
        
        fprintf(stderr, "* Read %d ALT contigs\n", c);
        
        if (which & BWA_IDX_PAC)
        {
            int64_t pac_bytes = idx->bns->l_pac/4+1;
            idx->pac = (uint8_t*) calloc(pac_bytes, 1);
            assert(idx->pac != NULL);
            bwamem_madv_hugepage(idx->pac, pac_bytes);
            err_fread_noeof(idx->pac, 1, pac_bytes, idx->bns->fp_pac); // concatenated 2-bit encoded sequence
            err_fclose(idx->bns->fp_pac);
            idx->bns->fp_pac = 0;
        }
    }
    free(prefix);
}

void indexEle::bwa_idx_load_ele_from_shm(uint8_t *base, size_t len)
{
    if (base == NULL) {
        fprintf(stderr, "ERROR! bwa_idx_load_ele_from_shm called with NULL base\n");
        exit(EXIT_FAILURE);
    }

    uint64_t off_bns = 0, sz_bns = 0;
    uint64_t off_ambs = 0, sz_ambs = 0;
    uint64_t off_anns = 0, sz_anns = 0;
    uint64_t off_names = 0, sz_names = 0;
    uint64_t off_pac = 0, sz_pac = 0;

    if (bwa_shm_section_find(base, BWA_SHM_SEC_BNS_STRUCT, &off_bns,  &sz_bns)  != 0
        || bwa_shm_section_find(base, BWA_SHM_SEC_BNS_AMBS,  &off_ambs, &sz_ambs) != 0
        || bwa_shm_section_find(base, BWA_SHM_SEC_BNS_ANNS,  &off_anns, &sz_anns) != 0
        || bwa_shm_section_find(base, BWA_SHM_SEC_BNS_NAMES, &off_names, &sz_names) != 0
        || bwa_shm_section_find(base, BWA_SHM_SEC_PAC,       &off_pac,  &sz_pac)  != 0)
    {
        fprintf(stderr, "ERROR! shm segment missing one or more BNS/PAC sections\n");
        exit(EXIT_FAILURE);
    }

    if (sz_bns != sizeof(bntseq_t)) {
        fprintf(stderr, "ERROR! BNS_STRUCT size %llu != sizeof(bntseq_t) %llu\n",
                (unsigned long long)sz_bns, (unsigned long long)sizeof(bntseq_t));
        exit(EXIT_FAILURE);
    }

    /* Heap-copy bns so the mapper can mutate fp_pac and friends without
     * writing into the shared segment. */
    idx->bns = (bntseq_t *) malloc(sizeof(bntseq_t));
    if (idx->bns == NULL) {
        fprintf(stderr, "ERROR! malloc(bntseq_t) failed in %s\n", __func__);
        exit(EXIT_FAILURE);
    }
    memcpy(idx->bns, base + off_bns, sizeof(bntseq_t));
    /* fp_pac is closed by the disk loader after slurping into idx->pac;
     * the segment was packed after that close, but zeroing here is cheap
     * insurance against a future writer that forgets. */
    idx->bns->fp_pac = NULL;

    /* Validate the trusted-by-default scalar fields immediately after the
     * memcpy from the segment. A corrupt or maliciously crafted segment
     * could otherwise drive negative/overflowing size computations below.
     * Bounds match the writer-side checks in bwa_shm_compute. */
    if (idx->bns->n_seqs < 0 || idx->bns->n_holes < 0 || idx->bns->l_pac <= 0) {
        fprintf(stderr,
            "ERROR! shm BNS_STRUCT scalars out of bounds "
            "(n_seqs=%d n_holes=%d l_pac=%lld)\n",
            idx->bns->n_seqs, idx->bns->n_holes, (long long)idx->bns->l_pac);
        exit(EXIT_FAILURE);
    }

    if ((int64_t)sz_ambs != (int64_t)idx->bns->n_holes * (int64_t)sizeof(bntamb1_t)) {
        fprintf(stderr, "ERROR! BNS_AMBS size %llu != n_holes %d * sizeof(bntamb1_t) %llu\n",
                (unsigned long long)sz_ambs, idx->bns->n_holes,
                (unsigned long long)sizeof(bntamb1_t));
        exit(EXIT_FAILURE);
    }
    if ((int64_t)sz_anns != (int64_t)idx->bns->n_seqs * (int64_t)sizeof(bntann1_t)) {
        fprintf(stderr, "ERROR! BNS_ANNS size %llu != n_seqs %d * sizeof(bntann1_t) %llu\n",
                (unsigned long long)sz_anns, idx->bns->n_seqs,
                (unsigned long long)sizeof(bntann1_t));
        exit(EXIT_FAILURE);
    }

    /* ambs is aliased into shm — read-only is fine. */
    idx->bns->ambs = (bntamb1_t *)(base + off_ambs);

    /* anns is heap-copied so the mapper can mutate is_alt etc. Each anns[i]
     * has a name and anno pointer that we patch to point into the on-segment
     * BNS_NAMES section (the original pointers came from the loader process
     * heap and are stale). */
    idx->bns->anns = (bntann1_t *) malloc((size_t)sz_anns);
    if (idx->bns->anns == NULL) {
        fprintf(stderr, "ERROR! malloc(%llu bytes for anns[]) failed in %s\n",
                (unsigned long long)sz_anns, __func__);
        exit(EXIT_FAILURE);
    }
    memcpy(idx->bns->anns, base + off_anns, (size_t)sz_anns);

    /* Walk BNS_NAMES once, patching anns[i].name/anno pointers. memchr
     * bounds each strlen so a malformed/truncated segment can't read past
     * the section. Walk must consume the section exactly. */
    {
        const char *cursor = (const char *)(base + off_names);
        const char *end    = cursor + sz_names;
        for (int i = 0; i < idx->bns->n_seqs; ++i) {
            const char *name_end =
                (const char *)memchr(cursor, '\0', (size_t)(end - cursor));
            if (name_end == NULL) {
                fprintf(stderr,
                    "ERROR! shm segment BNS_NAMES: unterminated name for seq %d\n", i);
                exit(EXIT_FAILURE);
            }
            idx->bns->anns[i].name = (char *)cursor;
            cursor = name_end + 1;

            const char *anno_end =
                (const char *)memchr(cursor, '\0', (size_t)(end - cursor));
            if (anno_end == NULL) {
                fprintf(stderr,
                    "ERROR! shm segment BNS_NAMES: unterminated anno for seq %d\n", i);
                exit(EXIT_FAILURE);
            }
            idx->bns->anns[i].anno = (char *)cursor;
            cursor = anno_end + 1;
        }
        if (cursor != end) {
            fprintf(stderr,
                "ERROR! shm segment BNS_NAMES walk consumed %lld of %llu bytes\n",
                (long long)(cursor - (const char *)(base + off_names)),
                (unsigned long long)sz_names);
            exit(EXIT_FAILURE);
        }
    }

    /* pac aliased into shm. The on-disk loader allocates a dedicated buffer;
     * here we just point at the segment. The size sanity check ensures the
     * segment matches what bns->l_pac demands. */
    int64_t expected_pac = idx->bns->l_pac / 4 + 1;
    if ((int64_t)sz_pac != expected_pac) {
        fprintf(stderr, "ERROR! BNS_PAC size %llu != l_pac/4+1 %lld\n",
                (unsigned long long)sz_pac, (long long)expected_pac);
        exit(EXIT_FAILURE);
    }
    idx->pac = (uint8_t *)(base + off_pac);

    /* The destructor's existing is_shm gate skips free(idx->mem). We set
     * idx->mem = base / idx->l_mem = len so the gate works as designed.
     * Note: idx->bns->ambs and idx->pac point inside [base, base+len); they
     * must NOT be freed individually. Likewise the patched anns[i].name /
     * anns[i].anno pointers. The destructor only frees idx->bns and
     * idx->bns->anns (both heap-copies above), gated by is_shm correctly. */
    idx->is_shm = 1;
    idx->mem    = base;
    idx->l_mem  = (int64_t)len;

    fprintf(stderr, "* BNS+PAC attached from shm: n_seqs=%d n_holes=%d l_pac=%ld\n",
            idx->bns->n_seqs, idx->bns->n_holes, (long)idx->bns->l_pac);
}

#include <sys/file.h>
char* indexEle::bwa_idx_infer_prefix(const char *hint)
{
    char *prefix;
    int l_hint;
    FILE *fp;
    l_hint = strlen(hint);
    prefix = (char *) malloc(l_hint + 3 + 4 + 1);
    assert(prefix != NULL);
    strcpy_s(prefix, l_hint + 3 + 4 + 1, hint);
    strcpy_s(prefix + l_hint, 8, ".64.bwt");
    if ((fp = fopen(prefix, "rb")) != 0)
    {
        fclose(fp);
        prefix[l_hint + 3] = 0;
        return prefix;
    } else {
        strcpy_s(prefix + l_hint, 8, ".bwt");
        if ((fp = fopen(prefix, "rb")) == 0)
        {
            free(prefix);
            return 0;
        } else {
            //flock(fileno(fp), 1);
            //flock(fileno(fp), 1);  // Unlock the file
            fclose(fp);
            prefix[l_hint] = 0;
            return prefix;
        }
    }
}

#if TEST
//int main(int argc, char* argv[])
//{
//  printf("Testing read_index_ele...\n");
//  indexEle *bwaEle = new indexEle();
//  
//  bwaEle->bwa_idx_load_ele("/projects/PCL-GBB/wasim/read_and_ref_data_1/hgaa.fa",
//                          BWA_IDX_ALL);
//
//  delete bwaEle;
//}
#endif
