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

Authors: Sanchit Misra <sanchit.misra@intel.com>; Vasimuddin Md <vasimuddin.md@intel.com>;
*****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "sais.h"
#include "FMI_search.h"
#include "memcpy_bwamem.h"
#include "profiling.h"

// Structured-exit on SMEM-buffer overflow. Called from every matchArray
// and match_buf write site when a bounds check fails. Prints a clear
// message and exits so memory corruption never propagates.
static void smem_overflow_die(const char *where, int64_t cap,
                              int64_t needed, int32_t read_len)
{
    fprintf(stderr,
        "FATAL: SMEM buffer overflow in %s: capacity=%" PRId64
        ", attempted write at index=%" PRId64 ", read_len=%d. This is a "
        "bwa-mem2 bug — please file an issue at "
        "https://github.com/fg-labs/bwa-mem2/issues with the input that "
        "triggered it.\n",
        where, cap, needed, read_len);
    exit(1);
}

#include "safestringlib.h"

FMI_search::FMI_search(const char *fname)
{
    fprintf(stderr, "* Entering FMI_search\n");
    //strcpy(file_name, fname);
    strcpy_s(file_name, PATH_MAX, fname);
    reference_seq_len = 0;
    sentinel_index = 0;
    index_alloc = 0;
    sa_ls_word = NULL;
    sa_ms_byte = NULL;
    cp_occ = NULL;
    one_hot_mask_array = NULL;
}

FMI_search::~FMI_search()
{
    if(sa_ms_byte)
        _mm_free(sa_ms_byte);
    if(sa_ls_word)
        _mm_free(sa_ls_word);
    if(cp_occ)
        _mm_free(cp_occ);
    if(one_hot_mask_array)
        _mm_free(one_hot_mask_array);
}

int64_t FMI_search::pac_seq_len(const char *fn_pac)
{
	FILE *fp;
	int64_t pac_len;
	uint8_t c;
	fp = xopen(fn_pac, "rb");
	err_fseek(fp, -1, SEEK_END);
	pac_len = err_ftell(fp);
	err_fread_noeof(&c, 1, 1, fp);
	err_fclose(fp);
	return (pac_len - 1) * 4 + (int)c;
}

void FMI_search::pac2nt(const char *fn_pac, std::string &reference_seq)
{
	uint8_t *buf2;
	int64_t i, pac_size, seq_len;
	FILE *fp;

	// initialization
	seq_len = pac_seq_len(fn_pac);
    assert(seq_len > 0);
    assert(seq_len <= 0x7fffffffffL);
	fp = xopen(fn_pac, "rb");

	// prepare sequence
	pac_size = (seq_len>>2) + ((seq_len&3) == 0? 0 : 1);
	buf2 = (uint8_t*)calloc(pac_size, 1);
    assert(buf2 != NULL);
	err_fread_noeof(buf2, 1, pac_size, fp);
	err_fclose(fp);
	for (i = 0; i < seq_len; ++i) {
		int nt = buf2[i>>2] >> ((3 - (i&3)) << 1) & 3;
        switch(nt)
        {
            case 0:
                reference_seq += "A";
            break;
            case 1:
                reference_seq += "C";
            break;
            case 2:
                reference_seq += "G";
            break;
            case 3:
                reference_seq += "T";
            break;
            default:
                fprintf(stderr, "ERROR! Value of nt is not in 0,1,2,3!");
                exit(EXIT_FAILURE);
        }
	}
    for(i = seq_len - 1; i >= 0; i--)
    {
        char c = reference_seq[i];
        switch(c)
        {
            case 'A':
                reference_seq += "T";
            break;
            case 'C':
                reference_seq += "G";
            break;
            case 'G':
                reference_seq += "C";
            break;
            case 'T':
                reference_seq += "A";
            break;
        }
    }
	free(buf2);
}

int FMI_search::build_fm_index(const char *ref_file_name, char *binary_seq, int64_t ref_seq_len, int64_t *sa_bwt, int64_t *count) {
    printf("ref_seq_len = %ld\n", ref_seq_len);
    fflush(stdout);

    char outname[PATH_MAX];

    strcpy_s(outname, PATH_MAX, ref_file_name);
    strcat_s(outname, PATH_MAX, CP_FILENAME_SUFFIX);
    //sprintf(outname, "%s.bwt.2bit.%d", ref_file_name, CP_BLOCK_SIZE);

    std::fstream outstream (outname, std::ios::out | std::ios::binary);
    outstream.seekg(0);	

    printf("count = %ld, %ld, %ld, %ld, %ld\n", count[0], count[1], count[2], count[3], count[4]);
    fflush(stdout);

    uint8_t *bwt;

    ref_seq_len++;
    outstream.write((char *)(&ref_seq_len), 1 * sizeof(int64_t));
    outstream.write((char*)count, 5 * sizeof(int64_t));

    int64_t i;
    int64_t ref_seq_len_aligned = ((ref_seq_len + CP_BLOCK_SIZE - 1) / CP_BLOCK_SIZE) * CP_BLOCK_SIZE;
    int64_t size = ref_seq_len_aligned * sizeof(uint8_t);
    bwt = (uint8_t *)_mm_malloc(size, 64);
    assert_not_null(bwt, size, index_alloc);

    int64_t sentinel_index = -1;
    for(i=0; i< ref_seq_len; i++)
    {
        if(sa_bwt[i] == 0)
        {
            bwt[i] = 4;
            printf("BWT[%ld] = 4\n", i);
            sentinel_index = i;
        }
        else
        {
            char c = binary_seq[sa_bwt[i]-1];
            switch(c)
            {
                case 0: bwt[i] = 0;
                          break;
                case 1: bwt[i] = 1;
                          break;
                case 2: bwt[i] = 2;
                          break;
                case 3: bwt[i] = 3;
                          break;
                default:
                        fprintf(stderr, "ERROR! i = %ld, c = %c\n", i, c);
                        exit(EXIT_FAILURE);
            }
        }
    }
    for(i = ref_seq_len; i < ref_seq_len_aligned; i++)
        bwt[i] = DUMMY_CHAR;


    printf("CP_SHIFT = %d, CP_MASK = %d\n", CP_SHIFT, CP_MASK);
    printf("sizeof CP_OCC = %ld\n", sizeof(CP_OCC));
    fflush(stdout);
    // create checkpointed occ
    int64_t cp_occ_size = (ref_seq_len >> CP_SHIFT) + 1;
    CP_OCC *cp_occ = NULL;

    size = cp_occ_size * sizeof(CP_OCC);
    cp_occ = (CP_OCC *)_mm_malloc(size, 64);
    assert_not_null(cp_occ, size, index_alloc);
    memset(cp_occ, 0, cp_occ_size * sizeof(CP_OCC));
    int64_t cp_count[16];

    memset(cp_count, 0, 16 * sizeof(int64_t));
    for(i = 0; i < ref_seq_len; i++)
    {
        if((i & CP_MASK) == 0)
        {
            CP_OCC cpo;
            cpo.cp_count[0] = cp_count[0];
            cpo.cp_count[1] = cp_count[1];
            cpo.cp_count[2] = cp_count[2];
            cpo.cp_count[3] = cp_count[3];

			int32_t j;
            cpo.one_hot_bwt_str[0] = 0;
            cpo.one_hot_bwt_str[1] = 0;
            cpo.one_hot_bwt_str[2] = 0;
            cpo.one_hot_bwt_str[3] = 0;

			for(j = 0; j < CP_BLOCK_SIZE; j++)
			{
                cpo.one_hot_bwt_str[0] = cpo.one_hot_bwt_str[0] << 1;
                cpo.one_hot_bwt_str[1] = cpo.one_hot_bwt_str[1] << 1;
                cpo.one_hot_bwt_str[2] = cpo.one_hot_bwt_str[2] << 1;
                cpo.one_hot_bwt_str[3] = cpo.one_hot_bwt_str[3] << 1;
				uint8_t c = bwt[i + j];
                //printf("c = %d\n", c);
                if(c < 4)
                {
                    cpo.one_hot_bwt_str[c] += 1;
                }
			}

            cp_occ[i >> CP_SHIFT] = cpo;
        }
        cp_count[bwt[i]]++;
    }
    outstream.write((char*)cp_occ, cp_occ_size * sizeof(CP_OCC));
    _mm_free(cp_occ);
    _mm_free(bwt);

    #if SA_COMPRESSION  

    size = ((ref_seq_len >> SA_COMPX)+ 1)  * sizeof(uint32_t);
    uint32_t *sa_ls_word = (uint32_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ls_word, size, index_alloc);
    size = ((ref_seq_len >> SA_COMPX) + 1) * sizeof(int8_t);
    int8_t *sa_ms_byte = (int8_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ms_byte, size, index_alloc);
    int64_t pos = 0;
    for(i = 0; i < ref_seq_len; i++)
    {
        if ((i & SA_COMPX_MASK) == 0)
        {
            sa_ls_word[pos] = sa_bwt[i] & 0xffffffff;
            sa_ms_byte[pos] = (sa_bwt[i] >> 32) & 0xff;
            pos++;
        }
    }
    fprintf(stderr, "pos: %d, ref_seq_len__: %ld\n", pos, ref_seq_len >> SA_COMPX);
    outstream.write((char*)sa_ms_byte, ((ref_seq_len >> SA_COMPX) + 1) * sizeof(int8_t));
    outstream.write((char*)sa_ls_word, ((ref_seq_len >> SA_COMPX) + 1) * sizeof(uint32_t));
    
    #else
    
    size = ref_seq_len * sizeof(uint32_t);
    uint32_t *sa_ls_word = (uint32_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ls_word, size, index_alloc);
    size = ref_seq_len * sizeof(int8_t);
    int8_t *sa_ms_byte = (int8_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ms_byte, size, index_alloc);
    for(i = 0; i < ref_seq_len; i++)
    {
        sa_ls_word[i] = sa_bwt[i] & 0xffffffff;
        sa_ms_byte[i] = (sa_bwt[i] >> 32) & 0xff;
    }
    outstream.write((char*)sa_ms_byte, ref_seq_len * sizeof(int8_t));
    outstream.write((char*)sa_ls_word, ref_seq_len * sizeof(uint32_t));
    
    #endif

    outstream.write((char *)(&sentinel_index), 1 * sizeof(int64_t));
    outstream.close();
    printf("max_occ_ind = %ld\n", i >> CP_SHIFT);    
    fflush(stdout);

    _mm_free(sa_ms_byte);
    _mm_free(sa_ls_word);
    return 0;
}

int FMI_search::build_index() {

    char *prefix = file_name;
    uint64_t startTick;
    startTick = __rdtsc();
    index_alloc = 0;

    std::string reference_seq;
    char pac_file_name[PATH_MAX];
    strcpy_s(pac_file_name, PATH_MAX, prefix);
    strcat_s(pac_file_name, PATH_MAX, ".pac");
    //sprintf(pac_file_name, "%s.pac", prefix);
    pac2nt(pac_file_name, reference_seq);
	int64_t pac_len = reference_seq.length();
    int status;
    int64_t size = pac_len * sizeof(char);
    char *binary_ref_seq = (char *)_mm_malloc(size, 64);
    index_alloc += size;
    assert_not_null(binary_ref_seq, size, index_alloc);
    char binary_ref_name[PATH_MAX];
    strcpy_s(binary_ref_name, PATH_MAX, prefix);
    strcat_s(binary_ref_name, PATH_MAX, ".0123");
    //sprintf(binary_ref_name, "%s.0123", prefix);
    std::fstream binary_ref_stream (binary_ref_name, std::ios::out | std::ios::binary);
    binary_ref_stream.seekg(0);
    fprintf(stderr, "init ticks = %llu\n", __rdtsc() - startTick);
    startTick = __rdtsc();
    int64_t i, count[16];
	memset(count, 0, sizeof(int64_t) * 16);
    for(i = 0; i < pac_len; i++)
    {
        switch(reference_seq[i])
        {
            case 'A':
            binary_ref_seq[i] = 0, ++count[0];
            break;
            case 'C':
            binary_ref_seq[i] = 1, ++count[1];
            break;
            case 'G':
            binary_ref_seq[i] = 2, ++count[2];
            break;
            case 'T':
            binary_ref_seq[i] = 3, ++count[3];
            break;
            default:
            binary_ref_seq[i] = 4;

        }
    }
    count[4]=count[0]+count[1]+count[2]+count[3];
    count[3]=count[0]+count[1]+count[2];
    count[2]=count[0]+count[1];
    count[1]=count[0];
    count[0]=0;
    fprintf(stderr, "ref seq len = %ld\n", pac_len);
    binary_ref_stream.write(binary_ref_seq, pac_len * sizeof(char));
    fprintf(stderr, "binary seq ticks = %llu\n", __rdtsc() - startTick);
    startTick = __rdtsc();

    size = (pac_len + 2) * sizeof(int64_t);
    int64_t *suffix_array=(int64_t *)_mm_malloc(size, 64);
    index_alloc += size;
    assert_not_null(suffix_array, size, index_alloc);
    startTick = __rdtsc();
	//status = saisxx<const char *, int64_t *, int64_t>(reference_seq.c_str(), suffix_array + 1, pac_len, 4);
	status = saisxx(reference_seq.c_str(), suffix_array + 1, pac_len);
	suffix_array[0] = pac_len;
    fprintf(stderr, "build suffix-array ticks = %llu\n", __rdtsc() - startTick);
    startTick = __rdtsc();

	build_fm_index(prefix, binary_ref_seq, pac_len, suffix_array, count);
    fprintf(stderr, "build fm-index ticks = %llu\n", __rdtsc() - startTick);
    _mm_free(binary_ref_seq);
    _mm_free(suffix_array);
    return 0;
}

void FMI_search::load_index()
{
    one_hot_mask_array = (uint64_t *)_mm_malloc(64 * sizeof(uint64_t), 64);
    one_hot_mask_array[0] = 0;
    uint64_t base = 0x8000000000000000L;
    one_hot_mask_array[1] = base;
    int64_t i = 0;
    for(i = 2; i < 64; i++)
    {
        one_hot_mask_array[i] = (one_hot_mask_array[i - 1] >> 1) | base;
    }

    char *ref_file_name = file_name;
    //beCalls = 0;
    char cp_file_name[PATH_MAX];
    strcpy_s(cp_file_name, PATH_MAX, ref_file_name);
    strcat_s(cp_file_name, PATH_MAX, CP_FILENAME_SUFFIX);

    // Read the BWT and FM index of the reference sequence
    FILE *cpstream = NULL;
    cpstream = fopen(cp_file_name,"rb");
    if (cpstream == NULL)
    {
        fprintf(stderr, "ERROR! Unable to open the file: %s\n", cp_file_name);
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stderr, "* Index file found. Loading index from %s\n", cp_file_name);
    }

    err_fread_noeof(&reference_seq_len, sizeof(int64_t), 1, cpstream);
    assert(reference_seq_len > 0);
    assert(reference_seq_len <= 0x7fffffffffL);

    fprintf(stderr, "* Reference seq len for bi-index = %ld\n", reference_seq_len);

    // create checkpointed occ
    int64_t cp_occ_size = (reference_seq_len >> CP_SHIFT) + 1;
    cp_occ = NULL;

    err_fread_noeof(&count[0], sizeof(int64_t), 5, cpstream);
    if ((cp_occ = (CP_OCC *)_mm_malloc(cp_occ_size * sizeof(CP_OCC), 64)) == NULL) {
        fprintf(stderr, "ERROR! unable to allocated cp_occ memory\n");
        exit(EXIT_FAILURE);
    }

    err_fread_noeof(cp_occ, sizeof(CP_OCC), cp_occ_size, cpstream);
    int64_t ii = 0;
    for(ii = 0; ii < 5; ii++)// update read count structure
    {
        count[ii] = count[ii] + 1;
    }

    #if SA_COMPRESSION

    int64_t reference_seq_len_ = (reference_seq_len >> SA_COMPX) + 1;
    sa_ms_byte = (int8_t *)_mm_malloc(reference_seq_len_ * sizeof(int8_t), 64);
    sa_ls_word = (uint32_t *)_mm_malloc(reference_seq_len_ * sizeof(uint32_t), 64);
    err_fread_noeof(sa_ms_byte, sizeof(int8_t), reference_seq_len_, cpstream);
    err_fread_noeof(sa_ls_word, sizeof(uint32_t), reference_seq_len_, cpstream);
    
    #else
    
    sa_ms_byte = (int8_t *)_mm_malloc(reference_seq_len * sizeof(int8_t), 64);
    sa_ls_word = (uint32_t *)_mm_malloc(reference_seq_len * sizeof(uint32_t), 64);
    err_fread_noeof(sa_ms_byte, sizeof(int8_t), reference_seq_len, cpstream);
    err_fread_noeof(sa_ls_word, sizeof(uint32_t), reference_seq_len, cpstream);

    #endif

    sentinel_index = -1;
    #if SA_COMPRESSION
    err_fread_noeof(&sentinel_index, sizeof(int64_t), 1, cpstream);
    fprintf(stderr, "* sentinel-index: %ld\n", sentinel_index);
    #endif
    fclose(cpstream);

    int64_t x;
    #if !SA_COMPRESSION
    for(x = 0; x < reference_seq_len; x++)
    {
        // fprintf(stderr, "x: %ld\n", x);
        #if SA_COMPRESSION
        if(get_sa_entry_compressed(x) == 0) {
            sentinel_index = x;
            break;
        }
        #else
        if(get_sa_entry(x) == 0) {
            sentinel_index = x;
            break;
        }
        #endif
    }
    fprintf(stderr, "\nsentinel_index: %ld\n", x);    
    #endif

    fprintf(stderr, "* Count:\n");
    for(x = 0; x < 5; x++)
    {
        fprintf(stderr, "%ld,\t%lu\n", x, (unsigned long)count[x]);
    }
    fprintf(stderr, "\n");  

    fprintf(stderr, "* Reading other elements of the index from files %s\n",
            ref_file_name);
    bwa_idx_load_ele(ref_file_name, BWA_IDX_ALL);

    fprintf(stderr, "* Done reading Index!!\n");
}

void FMI_search::getSMEMsOnePosOneThread(uint8_t *enc_qdb,
                                         int16_t *query_pos_array,
                                         int32_t *min_intv_array,
                                         int32_t *rid_array,
                                         int32_t numReads,
                                         int32_t batch_size,
                                         const bseq1_t *seq_,
                                         int32_t *query_cum_len_ar,
                                         int32_t max_readlength,
                                         int32_t minSeedLen,
                                         SMEM *matchArray,
                                         int64_t max_smem,
                                         int64_t *__numTotalSmem)
{
    int64_t numTotalSmem = *__numTotalSmem;
    // Heap-allocate prevArray instead of using a stack VLA. max_readlength
    // is now driven by batch content (no compile-time cap), so on long-read
    // workloads a stack array would risk exhaustion. Sized identically to
    // the lockstep path's per-slot buffer (max_readlength SMEMs).
    SMEM *prevArray = (SMEM *)_mm_malloc((size_t)max_readlength * sizeof(SMEM), 64);
    assert(prevArray != NULL);

    uint32_t i;
    // Perform SMEM for original reads
    for(i = 0; i < numReads; i++)
    {
        int x = query_pos_array[i];
        int32_t rid = rid_array[i];
        int next_x = x + 1;

        int readlength = seq_[rid].l_seq;
        int offset = query_cum_len_ar[rid];
        // uint8_t a = enc_qdb[rid * readlength + x];
        uint8_t a = enc_qdb[offset + x];

        if(a < 4)
        {
            SMEM smem;
            smem.rid = rid;
            smem.m = x;
            smem.n = x;
            smem.k = count[a];
            smem.l = count[3 - a];
            smem.s = count[a+1] - count[a];
            int numPrev = 0;
            
            int j;
            for(j = x + 1; j < readlength; j++)
            {
                // a = enc_qdb[rid * readlength + j];
                a = enc_qdb[offset + j];
                next_x = j + 1;
                if(a < 4)
                {
                    SMEM smem_ = smem;

                    // Forward extension is backward extension with the BWT of reverse complement
                    smem_.k = smem.l;
                    smem_.l = smem.k;
                    SMEM newSmem_ = backwardExt(smem_, 3 - a);
                    //SMEM newSmem_ = forwardExt(smem_, 3 - a);
                    SMEM newSmem = newSmem_;
                    newSmem.k = newSmem_.l;
                    newSmem.l = newSmem_.k;
                    newSmem.n = j;

                    int32_t s_neq_mask = newSmem.s != smem.s;

                    prevArray[numPrev] = smem;
                    numPrev += s_neq_mask;
                    if(newSmem.s < min_intv_array[i])
                    {
                        next_x = j;
                        break;
                    }
                    smem = newSmem;
#ifdef ENABLE_PREFETCH
                    _mm_prefetch((const char *)(&cp_occ[(smem.k) >> CP_SHIFT]), _MM_HINT_T0);
                    _mm_prefetch((const char *)(&cp_occ[(smem.l) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                }
                else
                {
                    break;
                }
            }
            if(smem.s >= min_intv_array[i])
            {

                prevArray[numPrev] = smem;
                numPrev++;
            }

            SMEM *prev;
            prev = prevArray;

            int p;
            for(p = 0; p < (numPrev/2); p++)
            {
                SMEM temp = prev[p];
                prev[p] = prev[numPrev - p - 1];
                prev[numPrev - p - 1] = temp;
            }

            // Backward search
            int cur_j = readlength;
            for(j = x - 1; j >= 0; j--)
            {
                int numCurr = 0;
                int curr_s = -1;
                // a = enc_qdb[rid * readlength + j];
                a = enc_qdb[offset + j];

                if(a > 3)
                {
                    break;
                }
                for(p = 0; p < numPrev; p++)
                {
                    SMEM smem = prev[p];
                    SMEM newSmem = backwardExt(smem, a);
                    newSmem.m = j;

                    if((newSmem.s < min_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen))
                    {
                        cur_j = j;

                        if (numTotalSmem >= max_smem)
                            smem_overflow_die("getSMEMsOnePosOneThread",
                                              max_smem, numTotalSmem,
                                              readlength);
                        matchArray[numTotalSmem++] = smem;
                        break;
                    }
                    if((newSmem.s >= min_intv_array[i]) && (newSmem.s != curr_s))
                    {
                        curr_s = newSmem.s;
                        prev[numCurr++] = newSmem;
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                        break;
                    }
                }
                p++;
                for(; p < numPrev; p++)
                {
                    SMEM smem = prev[p];

                    SMEM newSmem = backwardExt(smem, a);
                    newSmem.m = j;


                    if((newSmem.s >= min_intv_array[i]) && (newSmem.s != curr_s))
                    {
                        curr_s = newSmem.s;
                        prev[numCurr++] = newSmem;
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                    }
                }
                numPrev = numCurr;
                if(numCurr == 0)
                {
                    break;
                }
            }
            if(numPrev != 0)
            {
                SMEM smem = prev[0];
                if(((smem.n - smem.m + 1) >= minSeedLen))
                {

                    if (numTotalSmem >= max_smem)
                        smem_overflow_die("getSMEMsOnePosOneThread",
                                          max_smem, numTotalSmem, readlength);
                    matchArray[numTotalSmem++] = smem;
                }
                numPrev = 0;
            }
        }
        query_pos_array[i] = next_x;
    }
    (*__numTotalSmem) = numTotalSmem;
    _mm_free(prevArray);
}

// ===== Lockstep SMEM batching =====
// Per-slot state for the lockstep SMEM walk. One instance per in-flight read
// in the batch. Every field mirrors a per-read local in the scalar
// getSMEMsOnePosOneThread body, so parity with the scalar path follows from
// composing the existing primitives (backwardExt, count[], cp_occ) in the
// same sequence the scalar would.
//
// The prev[] and match_buf[] buffers are heap-owned by the caller (mmc in
// mem_collect_smem), sized from the batch's max_readlength, and handed to
// each slot as stride-offset pointers at driver entry. No compile-time cap.

enum LockstepPhase : uint8_t {
    PH_FWD      = 0,  // forward extension inner loop active
    PH_BWD_INIT = 1,  // between-phases housekeeping pending
    PH_BWD      = 2,  // backward search outer loop active
    PH_DONE     = 3   // slot finished; match_buf ready for flush
};

struct FMI_search::BatchSlot {
    // Input identity — copied at init, never mutated thereafter.
    int32_t input_idx;           // index into the caller's input arrays
    int32_t rid;                 // rid_array[input_idx]
    int16_t start_pos;           // query_pos_array[input_idx] (saved for bwd init)
    int32_t min_intv;            // min_intv_array[input_idx]
    int32_t readlength;          // seq_[rid].l_seq
    int32_t offset;              // query_cum_len_ar[rid]

    // Output for query_pos_array[input_idx] write-back at flush time.
    int16_t next_x;

    // Walk state (mirrors scalar locals).
    SMEM smem;                   // current SA interval
    int32_t j;                   // current query position in the active phase's inner loop
    int32_t cur_j;               // backward phase bookkeeping (scalar's cur_j)
    LockstepPhase phase;

    // Backward-search state. numCurr/curr_s are loop-local in
    // ls_advance_backward_step — they don't live across outer-j steps,
    // so they don't appear here.
    int32_t numPrev;
    SMEM   *prev;                // borrowed; caller-owned heap buffer, buf_cap SMEMs
    int32_t buf_cap;

    // Per-slot match buffer + reorder-buffer ready flag for in-order flush.
    SMEM   *match_buf;           // borrowed; caller-owned, buf_cap SMEMs
    int32_t match_count;
    bool    ready;
};

void FMI_search::ls_prefetch_cp_occ(const BatchSlot *s)
{
#ifdef ENABLE_PREFETCH
    _mm_prefetch((const char *)(&cp_occ[(s->smem.k) >> CP_SHIFT]), _MM_HINT_T0);
    _mm_prefetch((const char *)(&cp_occ[(s->smem.l) >> CP_SHIFT]), _MM_HINT_T0);
    _mm_prefetch((const char *)(&cp_occ[(s->smem.k + s->smem.s) >> CP_SHIFT]), _MM_HINT_T0);
    _mm_prefetch((const char *)(&cp_occ[(s->smem.l + s->smem.s) >> CP_SHIFT]), _MM_HINT_T0);
#else
    (void)s;
#endif
}

// Populate a slot from the caller's input arrays at the given input_idx.
// After this call either:
//   phase == PH_FWD  — slot is ready to step the forward-extension inner loop
//   phase == PH_DONE — the first base is non-ACGT; scalar skips this read, so
//                      we match that (zero matches emitted, ready to flush).
// prev_buf/match_buf_buf are borrowed stride views into a caller-owned
// per-thread heap allocation (buf_cap SMEMs each).
void FMI_search::ls_init_slot(BatchSlot *s,
                              int32_t input_idx,
                              const int16_t *query_pos_array,
                              const int32_t *min_intv_array,
                              const int32_t *rid_array,
                              const bseq1_t *seq_,
                              const int32_t *query_cum_len_ar,
                              const uint8_t *enc_qdb,
                              SMEM *prev_buf,
                              SMEM *match_buf_buf,
                              int32_t buf_cap)
{
    s->input_idx   = input_idx;
    s->rid         = rid_array[input_idx];
    s->start_pos   = query_pos_array[input_idx];
    s->min_intv    = min_intv_array[input_idx];
    s->readlength  = seq_[s->rid].l_seq;
    s->offset      = query_cum_len_ar[s->rid];
    s->next_x      = s->start_pos + 1;
    s->numPrev     = 0;
    s->prev        = prev_buf;
    s->match_buf   = match_buf_buf;
    s->buf_cap     = buf_cap;
    s->match_count = 0;
    s->ready       = false;

    int32_t x = s->start_pos;
    uint8_t a = enc_qdb[s->offset + x];
    if (a < 4) {
        s->smem.rid = s->rid;
        s->smem.m   = x;
        s->smem.n   = x;
        s->smem.k   = count[a];
        s->smem.l   = count[3 - a];
        s->smem.s   = count[a + 1] - count[a];
        s->j        = x + 1;
        s->phase    = PH_FWD;
    } else {
        // Scalar path skips the whole read when the first base is N.
        // Match that by going straight to DONE with zero matches.
        s->phase = PH_DONE;
        s->ready = true;
    }
}
// Advance one slot through one step of forward extension.
// Mirrors the inner j-loop body of the scalar getSMEMsOnePosOneThread
// (src/FMI_search.cpp — the for(j = x+1; j < readlength; j++) loop).
// On the step that exits forward-ext (end-of-read, non-ACGT at j, or
// s < min_intv), transitions phase to PH_BWD_INIT.
void FMI_search::ls_advance_forward_step(BatchSlot *s, const uint8_t *enc_qdb)
{
    if (s->j >= s->readlength) {
        // Ran off the end of the read; keep the still-valid smem if any.
        if (s->smem.s >= s->min_intv) {
            s->prev[s->numPrev++] = s->smem;
        }
        s->phase = PH_BWD_INIT;
        return;
    }

    uint8_t a = enc_qdb[s->offset + s->j];
    s->next_x = s->j + 1;
    if (a >= 4) {
        // Non-ACGT base terminates forward ext.
        if (s->smem.s >= s->min_intv) {
            s->prev[s->numPrev++] = s->smem;
        }
        s->phase = PH_BWD_INIT;
        return;
    }

    SMEM smem_ = s->smem;
    smem_.k = s->smem.l;
    smem_.l = s->smem.k;
    SMEM newSmem_ = backwardExt(smem_, 3 - a);
    SMEM newSmem  = newSmem_;
    newSmem.k = newSmem_.l;
    newSmem.l = newSmem_.k;
    newSmem.n = s->j;

    int32_t s_neq_mask = (newSmem.s != s->smem.s);
    s->prev[s->numPrev] = s->smem;
    s->numPrev += s_neq_mask;

    if (newSmem.s < s->min_intv) {
        s->next_x = s->j;
        s->phase = PH_BWD_INIT;
        return;
    }

    s->smem = newSmem;
    s->j++;
#ifdef ENABLE_PREFETCH
    _mm_prefetch((const char *)(&cp_occ[(s->smem.k) >> CP_SHIFT]), _MM_HINT_T0);
    _mm_prefetch((const char *)(&cp_occ[(s->smem.l) >> CP_SHIFT]), _MM_HINT_T0);
#endif
}

// Between-phases housekeeping: reverse prev[] and set up backward-phase
// cursors. After this call the slot is either PH_BWD (ready for
// ls_advance_backward_step) or PH_DONE (nothing to do — backward outer
// loop would terminate on the first step, so we short-circuit to the
// final prev[0] emit and done).
void FMI_search::ls_prepare_backward(BatchSlot *s)
{
    // Reverse prev[] in place (matches scalar behavior before backward loop).
    for (int p = 0; p < (s->numPrev / 2); p++) {
        SMEM tmp = s->prev[p];
        s->prev[p] = s->prev[s->numPrev - p - 1];
        s->prev[s->numPrev - p - 1] = tmp;
    }
    s->cur_j = s->readlength;
    s->j = s->start_pos - 1;  // first position for the backward outer loop

    if (s->numPrev == 0) {
        // Nothing to emit; scalar's final `if (numPrev != 0)` block is a no-op.
        s->phase = PH_DONE;
        s->ready = true;
        return;
    }
    if (s->j < 0) {
        // Backward outer loop cannot execute (start_pos == 0). Transition to
        // PH_BWD so the first ls_advance_backward_step call sees j < 0, jumps
        // to DONE, and runs the final prev[0] emit with the correct minSeedLen.
        s->phase = PH_BWD;
        return;
    }
    s->phase = PH_BWD;
}

// Advance one slot through ONE outer-j iteration of backward search.
// Mirrors the body of the scalar `for (j = x-1; j >= 0; j--)` outer loop
// in getSMEMsOnePosOneThread. On terminating conditions (numCurr == 0,
// j < 0 after decrement, or non-ACGT at j), runs the scalar's final
// prev[0] emit and transitions to PH_DONE.
void FMI_search::ls_advance_backward_step(BatchSlot *s,
                                          const uint8_t *enc_qdb,
                                          int32_t minSeedLen)
{
    if (s->j < 0) goto DONE;

    {
        int32_t numCurr = 0;
        int32_t curr_s = -1;
        uint8_t a = enc_qdb[s->offset + s->j];
        if (a > 3) goto DONE;

        int p;
        for (p = 0; p < s->numPrev; p++) {
            SMEM smem = s->prev[p];
            SMEM newSmem = backwardExt(smem, a);
            newSmem.m = s->j;
            if ((newSmem.s < s->min_intv) && ((smem.n - smem.m + 1) >= minSeedLen)) {
                s->cur_j = s->j;
                if (s->match_count >= s->buf_cap)
                    smem_overflow_die("ls_advance_backward_step/match_buf",
                                      s->buf_cap, s->match_count,
                                      s->readlength);
                s->match_buf[s->match_count++] = smem;
                break;
            }
            if ((newSmem.s >= s->min_intv) && (newSmem.s != curr_s)) {
                curr_s = newSmem.s;
                s->prev[numCurr++] = newSmem;
#ifdef ENABLE_PREFETCH
                _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                _mm_prefetch((const char *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                break;
            }
        }
        p++;
        for (; p < s->numPrev; p++) {
            SMEM smem = s->prev[p];
            SMEM newSmem = backwardExt(smem, a);
            newSmem.m = s->j;
            if ((newSmem.s >= s->min_intv) && (newSmem.s != curr_s)) {
                curr_s = newSmem.s;
                s->prev[numCurr++] = newSmem;
#ifdef ENABLE_PREFETCH
                _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                _mm_prefetch((const char *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
            }
        }
        s->numPrev = numCurr;
        s->j--;                // advance for the next outer step
        if (numCurr == 0) goto DONE;
        return;                // stay in PH_BWD; next call continues
    }

DONE:
    // Scalar end-of-function final prev[0] emit (lines ~650-659 of scalar).
    if (s->numPrev != 0) {
        SMEM smem = s->prev[0];
        if ((smem.n - smem.m + 1) >= minSeedLen) {
            if (s->match_count >= s->buf_cap)
                smem_overflow_die("ls_advance_backward_step/match_buf(final)",
                                  s->buf_cap, s->match_count, s->readlength);
            s->match_buf[s->match_count++] = smem;
        }
        s->numPrev = 0;
    }
    s->phase = PH_DONE;
    s->ready = true;
}
// ===== End lockstep scaffolding =====

void FMI_search::getSMEMsAllPosOneThread(uint8_t *enc_qdb,
                                         int32_t *min_intv_array,
                                         int32_t *rid_array,
                                         int32_t numReads,
                                         int32_t batch_size,
                                         const bseq1_t *seq_,
                                         int32_t *query_cum_len_ar,
                                         int32_t max_readlength,
                                         int32_t minSeedLen,
                                         SMEM *matchArray,
                                         int64_t max_smem,
                                         SMEM *lockstep_prev_base,
                                         SMEM *lockstep_match_buf_base,
                                         int64_t *__numTotalSmem)
{
    int16_t *query_pos_array = (int16_t *)_mm_malloc(numReads * sizeof(int16_t), 64);

    int32_t i;
    for(i = 0; i < numReads; i++)
        query_pos_array[i] = 0;

    int32_t numActive = numReads;
    (*__numTotalSmem) = 0;

    do
    {
        int32_t head = 0;
        int32_t tail = 0;
        for(head = 0; head < numActive; head++)
        {
            int readlength = seq_[rid_array[head]].l_seq;
            if(query_pos_array[head] < readlength)
            {
                rid_array[tail] = rid_array[head];
                query_pos_array[tail] = query_pos_array[head];
                min_intv_array[tail] = min_intv_array[head];
                tail++;
            }
        }
#if SMEM_LOCKSTEP_N > 1
        getSMEMsOnePosOneThread_lockstep(enc_qdb,
                                         query_pos_array,
                                         min_intv_array,
                                         rid_array,
                                         tail,
                                         batch_size,
                                         seq_,
                                         query_cum_len_ar,
                                         max_readlength,
                                         minSeedLen,
                                         matchArray,
                                         max_smem,
                                         lockstep_prev_base,
                                         lockstep_match_buf_base,
                                         __numTotalSmem);
#else
        (void)lockstep_prev_base;
        (void)lockstep_match_buf_base;
        getSMEMsOnePosOneThread(enc_qdb,
                                query_pos_array,
                                min_intv_array,
                                rid_array,
                                tail,
                                batch_size,
                                seq_,
                                query_cum_len_ar,
                                max_readlength,
                                minSeedLen,
                                matchArray,
                                max_smem,
                                __numTotalSmem);
#endif
        numActive = tail;
    } while(numActive > 0);

    _mm_free(query_pos_array);
}

void FMI_search::getSMEMsOnePosOneThread_lockstep(uint8_t *enc_qdb,
                                                   int16_t *query_pos_array,
                                                   int32_t *min_intv_array,
                                                   int32_t *rid_array,
                                                   int32_t numReads,
                                                   int32_t batch_size,
                                                   const bseq1_t *seq_,
                                                   int32_t *query_cum_len_ar,
                                                   int32_t max_readlength,
                                                   int32_t minSeedLen,
                                                   SMEM *matchArray,
                                                   int64_t max_smem,
                                                   SMEM *lockstep_prev_base,
                                                   SMEM *lockstep_match_buf_base,
                                                   int64_t *__numTotalSmem)
{
    (void)batch_size;

    if (numReads == 0) return;

    const int32_t N = SMEM_LOCKSTEP_N;
    // Per-slot prev[] and match_buf[] are borrowed stride views into a
    // caller-owned per-thread heap allocation, sized `max_readlength` SMEMs
    // each. BatchSlot itself (a handful of ints + pointers) lives on the
    // stack — small and bounded, no compile-time read-length cap.
    BatchSlot slots[SMEM_LOCKSTEP_N];
    const int32_t slot_cap = max_readlength;

    auto slot_prev      = [&](int32_t s) -> SMEM * {
        return lockstep_prev_base      + (int64_t)s * slot_cap;
    };
    auto slot_match_buf = [&](int32_t s) -> SMEM * {
        return lockstep_match_buf_base + (int64_t)s * slot_cap;
    };

    // Seed the first min(N, numReads) slots.
    int32_t initial = (numReads < N) ? numReads : N;
    for (int32_t s = 0; s < initial; s++) {
        ls_init_slot(&slots[s], s, query_pos_array, min_intv_array,
                     rid_array, seq_, query_cum_len_ar, enc_qdb,
                     slot_prev(s), slot_match_buf(s), slot_cap);
        if (slots[s].phase == PH_FWD) ls_prefetch_cp_occ(&slots[s]);
    }
    // Any unused slots (numReads < N) are marked DONE with invalid input_idx
    // so the stepping pass skips them and the flush pass ignores them.
    for (int32_t s = initial; s < N; s++) {
        slots[s].phase = PH_DONE;
        slots[s].ready = false;
        slots[s].input_idx = -1;
    }

    int32_t next_input   = initial;
    int32_t flush_cursor = 0;
    int64_t numTotalSmem = *__numTotalSmem;

    while (flush_cursor < numReads) {
        // --- Stepping pass: advance each non-DONE slot by one phase-step. ---
        for (int32_t s = 0; s < N; s++) {
            switch (slots[s].phase) {
                case PH_FWD:
                    ls_advance_forward_step(&slots[s], enc_qdb);
                    // If forward just transitioned to PH_BWD_INIT, run the
                    // between-phases housekeeping immediately. prepare_backward
                    // may transition straight to PH_DONE (e.g. empty prev[]).
                    if (slots[s].phase == PH_BWD_INIT) {
                        ls_prepare_backward(&slots[s]);
                    }
                    break;
                case PH_BWD_INIT:
                    ls_prepare_backward(&slots[s]);
                    break;
                case PH_BWD:
                    ls_advance_backward_step(&slots[s], enc_qdb, minSeedLen);
                    break;
                case PH_DONE:
                    break;
            }
        }

        // --- Flush pass: in-order emit + slot recycle. ---
        bool progress = true;
        while (progress) {
            progress = false;
            for (int32_t s = 0; s < N; s++) {
                if (slots[s].phase == PH_DONE &&
                    slots[s].ready &&
                    slots[s].input_idx == flush_cursor) {
                    for (int32_t m = 0; m < slots[s].match_count; m++) {
                        if (numTotalSmem >= max_smem)
                            smem_overflow_die("getSMEMsOnePosOneThread_lockstep",
                                              max_smem, numTotalSmem,
                                              slots[s].readlength);
                        matchArray[numTotalSmem++] = slots[s].match_buf[m];
                    }
                    query_pos_array[slots[s].input_idx] = slots[s].next_x;
                    flush_cursor++;

                    if (next_input < numReads) {
                        ls_init_slot(&slots[s], next_input,
                                     query_pos_array, min_intv_array,
                                     rid_array, seq_, query_cum_len_ar, enc_qdb,
                                     slot_prev(s), slot_match_buf(s), slot_cap);
                        if (slots[s].phase == PH_FWD) ls_prefetch_cp_occ(&slots[s]);
                        next_input++;
                    } else {
                        slots[s].input_idx = -1;  // retired
                        slots[s].ready = false;
                    }
                    progress = true;
                }
            }
        }
    }

    *__numTotalSmem = numTotalSmem;
}

int64_t FMI_search::bwtSeedStrategyAllPosOneThread(uint8_t *enc_qdb,
                                                   int32_t *max_intv_array,
                                                   int32_t numReads,
                                                   const bseq1_t *seq_,
                                                   int32_t *query_cum_len_ar,
                                                   int32_t minSeedLen,
                                                   SMEM *matchArray,
                                                   int64_t max_smem)
{
    int32_t i;

    int64_t numTotalSeed = 0;

    for(i = 0; i < numReads; i++)
    {
        int readlength = seq_[i].l_seq;
        int16_t x = 0;
        while(x < readlength)
        {
            int next_x = x + 1;

            // Forward search
            SMEM smem;
            smem.rid = i;
            smem.m = x;
            smem.n = x;
            
            int offset = query_cum_len_ar[i];
            uint8_t a = enc_qdb[offset + x];
            // uint8_t a = enc_qdb[i * readlength + x];

            if(a < 4)
            {
                smem.k = count[a];
                smem.l = count[3 - a];
                smem.s = count[a+1] - count[a];


                int j;
                for(j = x + 1; j < readlength; j++)
                {
                    next_x = j + 1;
                    // a = enc_qdb[i * readlength + j];
                    a = enc_qdb[offset + j];
                    if(a < 4)
                    {
                        SMEM smem_ = smem;

                        // Forward extension is backward extension with the BWT of reverse complement
                        smem_.k = smem.l;
                        smem_.l = smem.k;
                        SMEM newSmem_ = backwardExt(smem_, 3 - a);
                        //SMEM smem = backwardExt(smem, 3 - a);
                        //smem.n = j;
                        SMEM newSmem = newSmem_;
                        newSmem.k = newSmem_.l;
                        newSmem.l = newSmem_.k;
                        newSmem.n = j;
                        smem = newSmem;
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(smem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(smem.l) >> CP_SHIFT]), _MM_HINT_T0);
#endif


                        if((smem.s < max_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen))
                        {

                            if(smem.s > 0)
                            {
                                if (numTotalSeed >= max_smem)
                                    smem_overflow_die("bwtSeedStrategyAllPosOneThread",
                                                      max_smem, numTotalSeed,
                                                      readlength);
                                matchArray[numTotalSeed++] = smem;
                            }
                            break;
                        }
                    }
                    else
                    {

                        break;
                    }
                }

            }
            x = next_x;
        }
    }
    return numTotalSeed;
}


void FMI_search::getSMEMs(uint8_t *enc_qdb,
        int32_t numReads,
        int32_t batch_size,
        int32_t readlength,
        int32_t minSeedLen,
        int32_t nthreads,
        SMEM *matchArray,
        int64_t *numTotalSmem)
{
    SMEM *prevArray = (SMEM *)_mm_malloc(nthreads * readlength * sizeof(SMEM), 64);
    SMEM *currArray = (SMEM *)_mm_malloc(nthreads * readlength * sizeof(SMEM), 64);


// #pragma omp parallel num_threads(nthreads)
    {
        int tid = 0; //omp_get_thread_num();   // removed omp
        numTotalSmem[tid] = 0;
        SMEM *myPrevArray = prevArray + tid * readlength;
        SMEM *myCurrArray = prevArray + tid * readlength;

        int32_t perThreadQuota = (numReads + (nthreads - 1)) / nthreads;
        int32_t first = tid * perThreadQuota;
        int32_t last  = (tid + 1) * perThreadQuota;
        if(last > numReads) last = numReads;
        SMEM *myMatchArray = matchArray + first * readlength;

        uint32_t i;
        // Perform SMEM for original reads
        for(i = first; i < last; i++)
        {
            int x = readlength - 1;
            int numPrev = 0;
            int numSmem = 0;

            while (x >= 0)
            {
                // Forward search
                SMEM smem;
                smem.rid = i;
                smem.m = x;
                smem.n = x;
                uint8_t a = enc_qdb[i * readlength + x];

                if(a > 3)
                {
                    x--;
                    continue;
                }
                smem.k = count[a];
                smem.l = count[3 - a];
                smem.s = count[a+1] - count[a];

                int j;
                for(j = x + 1; j < readlength; j++)
                {
                    a = enc_qdb[i * readlength + j];
                    if(a < 4)
                    {
                        SMEM smem_ = smem;

                        // Forward extension is backward extension with the BWT of reverse complement
                        smem_.k = smem.l;
                        smem_.l = smem.k;
                        SMEM newSmem_ = backwardExt(smem_, 3 - a);
                        SMEM newSmem = newSmem_;
                        newSmem.k = newSmem_.l;
                        newSmem.l = newSmem_.k;
                        newSmem.n = j;

                        if(newSmem.s != smem.s)
                        {
                            myPrevArray[numPrev] = smem;
                            numPrev++;
                        }
                        smem = newSmem;
                        if(newSmem.s == 0)
                        {
                            break;
                        }
                    }
                    else
                    {
                        myPrevArray[numPrev] = smem;
                        numPrev++;
                        break;
                    }
                }
                if(smem.s != 0)
                {
                    myPrevArray[numPrev++] = smem;
                }

                SMEM *curr, *prev;
                prev = myPrevArray;
                curr = myCurrArray;

                int p;
                for(p = 0; p < (numPrev/2); p++)
                {
                    SMEM temp = prev[p];
                    prev[p] = prev[numPrev - p - 1];
                    prev[numPrev - p - 1] = temp;
                }

                int next_x = x - 1;

                // Backward search
                int cur_j = readlength;
                for(j = x - 1; j >= 0; j--)
                {
                    int numCurr = 0;
                    int curr_s = -1;
                    a = enc_qdb[i * readlength + j];
                    //printf("a = %d\n", a);
                    if(a > 3)
                    {
                        next_x = j - 1;
                        break;
                    }
                    for(p = 0; p < numPrev; p++)
                    {
                        SMEM smem = prev[p];
                        SMEM newSmem = backwardExt(smem, a);
                        newSmem.m = j;

                        if(newSmem.s == 0)
                        {
                            if((numCurr == 0) && (j < cur_j))
                            {
                                cur_j = j;
                                if((smem.n - smem.m + 1) >= minSeedLen)
                                    myMatchArray[numTotalSmem[tid] + numSmem++] = smem;
                            }
                        }
                        if((newSmem.s != 0) && (newSmem.s != curr_s))
                        {
                            curr_s = newSmem.s;
                            curr[numCurr++] = newSmem;
                        }
                    }
                    SMEM *temp = prev;
                    prev = curr;
                    curr = temp;
                    numPrev = numCurr;
                    if(numCurr == 0)
                    {
                        next_x = j;
                        break;
                    }
                    else
                    {
                        next_x = j - 1;
                    }
                }
                if(numPrev != 0)
                {
                    SMEM smem = prev[0];
                    if((smem.n - smem.m + 1) >= minSeedLen)
                        myMatchArray[numTotalSmem[tid] + numSmem++] = smem;
                    numPrev = 0;
                }
                x = next_x;
            }
            numTotalSmem[tid] += numSmem;
        }
    }

    _mm_free(prevArray);
    _mm_free(currArray);
}


int compare_smem(const void *a, const void *b)
{
    SMEM *pa = (SMEM *)a;
    SMEM *pb = (SMEM *)b;

    if(pa->rid < pb->rid)
        return -1;
    if(pa->rid > pb->rid)
        return 1;

    if(pa->m < pb->m)
        return -1;
    if(pa->m > pb->m)
        return 1;
    if(pa->n > pb->n)
        return -1;
    if(pa->n < pb->n)
        return 1;
    return 0;
}

void FMI_search::sortSMEMs(SMEM *matchArray,
        int64_t numTotalSmem[],
        int32_t numReads,
        int32_t readlength,
        int nthreads)
{
    int tid;
    int32_t perThreadQuota = (numReads + (nthreads - 1)) / nthreads;
    for(tid = 0; tid < nthreads; tid++)
    {
        int32_t first = tid * perThreadQuota;
        SMEM *myMatchArray = matchArray + first * readlength;
        qsort(myMatchArray, numTotalSmem[tid], sizeof(SMEM), compare_smem);
    }
}


SMEM FMI_search::backwardExt(SMEM smem, uint8_t a)
{
    //beCalls++;
    uint8_t b;

    int64_t k[4], l[4], s[4];
    for(b = 0; b < 4; b++)
    {
        int64_t sp = (int64_t)(smem.k);
        int64_t ep = (int64_t)(smem.k) + (int64_t)(smem.s);
        GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);
        GET_OCC(ep, b, occ_id_ep, y_ep, occ_ep, one_hot_bwt_str_c_ep, match_mask_ep);
        k[b] = count[b] + occ_sp;
        s[b] = occ_ep - occ_sp;
    }

    int64_t sentinel_offset = 0;
    if((smem.k <= sentinel_index) && ((smem.k + smem.s) > sentinel_index)) sentinel_offset = 1;
    l[3] = smem.l + sentinel_offset;
    l[2] = l[3] + s[3];
    l[1] = l[2] + s[2];
    l[0] = l[1] + s[1];

    smem.k = k[a];
    smem.l = l[a];
    smem.s = s[a];
    return smem;
}

int64_t FMI_search::get_sa_entry(int64_t pos)
{
    int64_t sa_entry = sa_ms_byte[pos];
    sa_entry = sa_entry << 32;
    sa_entry = sa_entry + sa_ls_word[pos];
    return sa_entry;
}

void FMI_search::get_sa_entries(int64_t *posArray, int64_t *coordArray, uint32_t count, int32_t nthreads)
{
    uint32_t i;
// #pragma omp parallel for num_threads(nthreads)
    for(i = 0; i < count; i++)
    {
        int64_t pos = posArray[i];
        int64_t sa_entry = sa_ms_byte[pos];
        sa_entry = sa_entry << 32;
        sa_entry = sa_entry + sa_ls_word[pos];
        //_mm_prefetch((const char *)(sa_ms_byte + pos + SAL_PFD), _MM_HINT_T0);
        coordArray[i] = sa_entry;
    }
}

void FMI_search::get_sa_entries(SMEM *smemArray, int64_t *coordArray, int32_t *coordCountArray, uint32_t count, int32_t max_occ)
{
    uint32_t i;
    int32_t totalCoordCount = 0;
    for(i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for(j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
        {
            int64_t pos = j;
            int64_t sa_entry = sa_ms_byte[pos];
            sa_entry = sa_entry << 32;
            sa_entry = sa_entry + sa_ls_word[pos];
            //_mm_prefetch((const char *)(sa_ms_byte + pos + SAL_PFD * step), _MM_HINT_T0);
            coordArray[totalCoordCount + c] = sa_entry;
        }
        coordCountArray[i] = c;
        totalCoordCount += c;
    }
}

// sa_compression
int64_t FMI_search::get_sa_entry_compressed(int64_t pos, int tid)
{
    if ((pos & SA_COMPX_MASK) == 0) {
        
        #if  SA_COMPRESSION
        int64_t sa_entry = sa_ms_byte[pos >> SA_COMPX];
        #else
        int64_t sa_entry = sa_ms_byte[pos];     // simulation
        #endif
        
        sa_entry = sa_entry << 32;
        
        #if  SA_COMPRESSION
        sa_entry = sa_entry + sa_ls_word[pos >> SA_COMPX];
        #else
        sa_entry = sa_entry + sa_ls_word[pos];   // simulation
        #endif
        
        return sa_entry;        
    }
    else {
        // tprof[MEM_CHAIN][tid] ++;
        int64_t offset = 0; 
        int64_t sp = pos;
        while(true)
        {
            int64_t occ_id_pp_ = sp >> CP_SHIFT;
            int64_t y_pp_ = CP_BLOCK_SIZE - (sp & CP_MASK) - 1; 
            uint64_t *one_hot_bwt_str = cp_occ[occ_id_pp_].one_hot_bwt_str;
            uint8_t b;

            if((one_hot_bwt_str[0] >> y_pp_) & 1)
                b = 0;
            else if((one_hot_bwt_str[1] >> y_pp_) & 1)
                b = 1;
            else if((one_hot_bwt_str[2] >> y_pp_) & 1)
                b = 2;
            else if((one_hot_bwt_str[3] >> y_pp_) & 1)
                b = 3;
            else
                b = 4;

            if (b == 4) {
                return offset;
            }

            GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);

            sp = count[b] + occ_sp;
            
            offset ++;
            // tprof[ALIGN1][tid] ++;
            if ((sp & SA_COMPX_MASK) == 0) break;
        }
        // assert((reference_seq_len >> SA_COMPX) - 1 >= (sp >> SA_COMPX));
        #if  SA_COMPRESSION
        int64_t sa_entry = sa_ms_byte[sp >> SA_COMPX];
        #else
        int64_t sa_entry = sa_ms_byte[sp];      // simultion
        #endif
        
        sa_entry = sa_entry << 32;

        #if  SA_COMPRESSION
        sa_entry = sa_entry + sa_ls_word[sp >> SA_COMPX];
        #else
        sa_entry = sa_entry + sa_ls_word[sp];      // simulation
        #endif
        
        sa_entry += offset;
        return sa_entry;
    }
}

void FMI_search::get_sa_entries(SMEM *smemArray, int64_t *coordArray, int32_t *coordCountArray, uint32_t count, int32_t max_occ, int tid)
{
    
    uint32_t i;
    int32_t totalCoordCount = 0;
    for(i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for(j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
        {
            int64_t pos = j;
            int64_t sa_entry = get_sa_entry_compressed(pos, tid);
            coordArray[totalCoordCount + c] = sa_entry;
        }
        // coordCountArray[i] = c;
        *coordCountArray += c;
        totalCoordCount += c;
    }
}

// SA_COPMRESSION w/ PREFETCH
int64_t FMI_search::call_one_step(int64_t pos, int64_t &sa_entry, int64_t &offset)
{
    if ((pos & SA_COMPX_MASK) == 0) {        
        sa_entry = sa_ms_byte[pos >> SA_COMPX];        
        sa_entry = sa_entry << 32;        
        sa_entry = sa_entry + sa_ls_word[pos >> SA_COMPX];        
        // return sa_entry;
        return 1;
    }
    else {
        // int64_t offset = 0; 
        int64_t sp = pos;

        int64_t occ_id_pp_ = sp >> CP_SHIFT;
        int64_t y_pp_ = CP_BLOCK_SIZE - (sp & CP_MASK) - 1; 
        uint64_t *one_hot_bwt_str = cp_occ[occ_id_pp_].one_hot_bwt_str;
        uint8_t b;

        if((one_hot_bwt_str[0] >> y_pp_) & 1)
            b = 0;
        else if((one_hot_bwt_str[1] >> y_pp_) & 1)
            b = 1;
        else if((one_hot_bwt_str[2] >> y_pp_) & 1)
            b = 2;
        else if((one_hot_bwt_str[3] >> y_pp_) & 1)
            b = 3;
        else
            b = 4;
        if (b == 4) {
            sa_entry = 0;
            return 1;
        }
        
        GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);
        
        sp = count[b] + occ_sp;
        
        offset ++;
        if ((sp & SA_COMPX_MASK) == 0) {
    
            sa_entry = sa_ms_byte[sp >> SA_COMPX];        
            sa_entry = sa_entry << 32;
            sa_entry = sa_entry + sa_ls_word[sp >> SA_COMPX];
            
            sa_entry += offset;
            // return sa_entry;
            return 1;
        }
        else {
            sa_entry = sp;
            return 0;
        }
    } // else
}

void FMI_search::get_sa_entries_prefetch(SMEM *smemArray, int64_t *coordArray,
                                         int64_t *coordCountArray, int64_t count,
                                         const int32_t max_occ, int tid, int64_t &id_)
{
    
    // uint32_t i;
    int32_t totalCoordCount = 0;
    int32_t mem_lim = 0, id = 0;
    
    for(int i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        mem_lim += smem.s;
    }

    int64_t *pos_ar = (int64_t *) _mm_malloc( mem_lim * sizeof(int64_t), 64);
    int64_t *map_ar = (int64_t *) _mm_malloc( mem_lim * sizeof(int64_t), 64);

    for(int i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for(j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
        {
            int64_t pos = j;
             pos_ar[id]  = pos;
             map_ar[id++] = totalCoordCount + c;
            // int64_t sa_entry = get_sa_entry_compressed(pos, tid);
            // coordArray[totalCoordCount + c] = sa_entry;
        }
        //coordCountArray[i] = c;
        *coordCountArray += c;
        totalCoordCount += c;
    }
    
    id_ += id;
    
    const int32_t sa_batch_size = 20;
    int64_t working_set[sa_batch_size], map_pos[sa_batch_size];;
    int64_t offset[sa_batch_size] = {-1};
    
    int i = 0, j = 0;    
    while(i<id && j<sa_batch_size)
    {
        int64_t pos =  pos_ar[i];
        working_set[j] = pos;
        map_pos[j] = map_ar[i];
        offset[j] = 0;
        
        if (pos & SA_COMPX_MASK == 0) {
            _mm_prefetch(&sa_ms_byte[pos >> SA_COMPX], _MM_HINT_T0);
            _mm_prefetch(&sa_ls_word[pos >> SA_COMPX], _MM_HINT_T0);
        }
        else {
            int64_t occ_id_pp_ = pos >> CP_SHIFT;
            _mm_prefetch(&cp_occ[occ_id_pp_], _MM_HINT_T0);
        }
        i++;
        j++;
    }
        
    int lim = j, all_quit = 0;
    while (all_quit < id)
    {
        
        for (int k=0; k<lim; k++)
        {
            int64_t sp = 0, pos = 0;
            bool quit;
            if (offset[k] >= 0) {
                quit = call_one_step(working_set[k], sp, offset[k]);
            }
            else
                continue;
            
            if (quit) {
                coordArray[map_pos[k]] = sp;
                all_quit ++;
                
                if (i < id)
                {
                    pos = pos_ar[i];
                    working_set[k] = pos;
                    map_pos[k] = map_ar[i++];
                    offset[k] = 0;
                    
                    if (pos & SA_COMPX_MASK == 0) {
                        _mm_prefetch(&sa_ms_byte[pos >> SA_COMPX], _MM_HINT_T0);
                        _mm_prefetch(&sa_ls_word[pos >> SA_COMPX], _MM_HINT_T0);
                    }
                    else {
                        int64_t occ_id_pp_ = pos >> CP_SHIFT;
                        _mm_prefetch(&cp_occ[occ_id_pp_], _MM_HINT_T0);
                    }
                }
                else
                    offset[k] = -1;
            }
            else {
                working_set[k] = sp;
                if (sp & SA_COMPX_MASK == 0) {
                    _mm_prefetch(&sa_ms_byte[sp >> SA_COMPX], _MM_HINT_T0);
                    _mm_prefetch(&sa_ls_word[sp >> SA_COMPX], _MM_HINT_T0);
                }
                else {
                    int64_t occ_id_pp_ = sp >> CP_SHIFT;
                    _mm_prefetch(&cp_occ[occ_id_pp_], _MM_HINT_T0);
                }                
            }
        }
    }
    
    _mm_free(pos_ar);
    _mm_free(map_ar);
}
