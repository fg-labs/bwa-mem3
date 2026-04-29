/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019 Intel Corporation, Heng Li.

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

Contacts: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
                                Heng Li <hli@jimmy.harvard.edu> 
*****************************************************************************************/

// ----------------------------------
#include "main.h"
#include "version.h"
#include "bwa_shm.h"

#ifdef USE_MIMALLOC
#include <mimalloc.h>
#endif


// ----------------------------------
// Profiling globals are now defined in profiling.cpp so they live in
// libbwa.a and are visible to library consumers that don't link main.o.
// ----------------------------------

int usage()
{
    fprintf(stderr, "Usage: bwa-mem3 <command> <arguments>\n");
    fprintf(stderr, "Commands:\n");
    fprintf(stderr, "  index         create index (add --meth to build a bwameth-style doubled c2t reference)\n");
    fprintf(stderr, "  mem           alignment (add --meth for bisulfite-seq: inline c2t + BAM output)\n");
    fprintf(stderr, "  shm           load/list/drop the index in POSIX shared memory\n");
    fprintf(stderr, "  version       print version number\n");
    fprintf(stderr, "Run `bwa-mem3 <command> --help` for command-specific options.\n");
    return 1;
}

// Append a single argv token to the @PG CL: value. Tabs inside the token
// (common when the caller passes e.g. `-R $'@RG\tID:x\tSM:y'`) would
// otherwise bleed into the @PG line as extra tag-separators, producing
// SAM that strict validators reject (issue #45 / upstream bwa-mem2#293).
// Newlines and carriage returns would be worse still: a literal \n
// terminates the @PG record mid-line and corrupts the whole header.
// Replace any of these with a single space; do not mutate argv itself.
static void append_pg_cl_arg(kstring_t *pg, const char *arg)
{
    kputc(' ', pg);
    for (const char *c = arg; *c != '\0'; ++c) {
        kputc((*c == '\t' || *c == '\n' || *c == '\r') ? ' ' : *c, pg);
    }
}

int main(int argc, char* argv[])
{
        
    // ---------------------------------    
    uint64_t tim = __rdtsc();
    sleep(1);
    proc_freq = __rdtsc() - tim;

    int ret = -1;
    if (argc < 2) return usage();

    if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        usage();
        return 0;
    }

    if (strcmp(argv[1], "index") == 0)
    {
         uint64_t tim = __rdtsc();
         ret = bwa_index(argc-1, argv+1);
         fprintf(stderr, "Total time taken: %0.4lf\n", (__rdtsc() - tim)*1.0/proc_freq);
         return ret;
    }
    else if (strcmp(argv[1], "mem") == 0)
    {
        // Short-circuit `mem --help` so we skip the AVX/SA banner and the
        // post-run profiling trailer (printed below when ret==0). Skip
        // tokens that are the value of an option that takes an argument
        // -- otherwise `mem -R --help ...`, `mem -o --help ...`, or
        // `mem --set-as-failed --help ...` would treat the option's
        // value as a help request and suppress the banner. Keep the
        // short-option string and long-option list in sync with the
        // optstring and long_opts table in fastmap.cpp::main_mem.
        static const char *const MEM_SHORT_OPTS_WITH_ARG =
            "kcvsrtRABOEUwLdTQDmINWxGhyKXHofz";
        static const char *const MEM_LONG_OPTS_WITH_ARG[] = {
            "--set-as-failed", "--supp-rep-hard-cap", NULL,
        };
        for (int i = 2; i < argc; ++i) {
            const char *t = argv[i];
            if (strcmp(t, "--") == 0) break;
            if (strcmp(t, "--help") == 0) return main_mem(argc-1, argv+1);
            // Long option (no `=`): next argv token is its value.
            int matched_long = 0;
            for (int j = 0; MEM_LONG_OPTS_WITH_ARG[j]; ++j) {
                if (strcmp(t, MEM_LONG_OPTS_WITH_ARG[j]) == 0) {
                    matched_long = 1; break;
                }
            }
            if (matched_long) { ++i; continue; }
            // Short-option bundle `-abc[value]`: walk chars; the first
            // arg-taking option consumes the rest of the bundle (if any)
            // or the next argv token (if the bundle ends at it).
            if (t[0] != '-' || t[1] == '\0' || t[1] == '-') continue;
            for (const char *p = t + 1; *p; ++p) {
                if (strchr(MEM_SHORT_OPTS_WITH_ARG, *p)) {
                    if (*(p + 1) == '\0') ++i;
                    break;
                }
            }
        }

        tprof[MEM][0] = __rdtsc();
        kstring_t pg = {0,0,0};
        extern char *bwa_pg;

        fprintf(stderr, "-----------------------------\n");
#if __AVX512BW__
        fprintf(stderr, "Executing in AVX512 mode!!\n");
#elif __AVX2__
        fprintf(stderr, "Executing in AVX2 mode!!\n");
#elif __AVX__
        fprintf(stderr, "Executing in AVX mode!!\n");        
#elif __SSE4_2__
        fprintf(stderr, "Executing in SSE4.2 mode!!\n");
#elif __SSE4_1__
        fprintf(stderr, "Executing in SSE4.1 mode!!\n");        
#endif
        fprintf(stderr, "-----------------------------\n");

        #if SA_COMPRESSION
        fprintf(stderr, "* SA compression enabled with xfactor: %d\n", 0x1 << SA_COMPX);
        #endif
        
        ksprintf(&pg, "@PG\tID:bwa-mem3\tPN:bwa-mem3\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);

        for (int i = 1; i < argc; ++i) append_pg_cl_arg(&pg, argv[i]);
        ksprintf(&pg, "\n");
        bwa_pg = pg.s;
        ret = main_mem(argc-1, argv+1);
        free(bwa_pg);
        
        /** Enable this return to avoid printing of the runtime profiling **/
        //return ret;
    }
    else if (strcmp(argv[1], "version") == 0)
    {
        puts(PACKAGE_VERSION);
#ifdef USE_MIMALLOC
        {
            int mv = mi_version();
            fprintf(stderr, "mimalloc %d.%d.%d\n", mv / 10000, (mv / 100) % 100, mv % 100);
        }
#endif
        return 0;
    }
    else if (strcmp(argv[1], "shm") == 0)
    {
        return main_shm(argc - 1, argv + 1);
    }
    else {
        fprintf(stderr, "ERROR: unknown command '%s'\n", argv[1]);
        return 1;
    }

    if (ret == 0) {
        fprintf(stderr, "\nImportant parameter settings: \n");
        fprintf(stderr, "\tBATCH_SIZE: %d\n", BATCH_SIZE);
        fprintf(stderr, "\tMAX_SEQ_LEN_REF: %d\n", MAX_SEQ_LEN_REF);
        fprintf(stderr, "\tMAX_SEQ_LEN_QER: %d\n", MAX_SEQ_LEN_QER);
        fprintf(stderr, "\tMAX_SEQ_LEN8: %d\n", MAX_SEQ_LEN8);
        fprintf(stderr, "\tSEEDS_PER_READ: %d\n", SEEDS_PER_READ);
        fprintf(stderr, "\tSIMD_WIDTH8 X: %d\n", SIMD_WIDTH8);
        fprintf(stderr, "\tSIMD_WIDTH16 X: %d\n", SIMD_WIDTH16);
        fprintf(stderr, "\tAVG_SEEDS_PER_READ: %d\n", AVG_SEEDS_PER_READ);
    }
    
    return ret;
}
