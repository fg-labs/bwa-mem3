/* The `bwa-mem3 aln` subcommand: bwa's aln + samse/sampe, fused. New code. */
#ifndef BWA_MEM3_ALN_MAIN_H
#define BWA_MEM3_ALN_MAIN_H

#include "aln_out.h"
#include "bwtaln.h"

namespace bwa_aln {

struct aln_opt_t {
	gap_opt_t gap;              // stock aln's options (n_threads included)
	int max_xa_hits;            // --max-xa-hits: stock samse -n (default 3)
	int batch_size;             // reads per batch (default 0x40000); a test hook, not a CLI option
};

/* Stock defaults: gap_init_opt(), max_xa_hits 3, batch_size 0x40000. */
void aln_opt_init(aln_opt_t *o);

/* Single-end run: reads fn_fq in batches and appends each batch's SAM records to `out`
 * (the header is the caller's). Seeds drand48 once from bns->seed. `pac` is the packed
 * forward reference. Returns 0, or 1 on an input or output error. */
int aln_run_se(const aln_opt_t *o, const bwt_t *bwt, const bntseq_t *bns, ubyte_t *pac,
			   const char *fn_fq, aln_out_t *out);

} // namespace bwa_aln

/* argv[0] is "aln". */
int main_aln(int argc, char *argv[]);

#endif
