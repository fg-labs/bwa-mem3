/* SPDX-License-Identifier: MIT */
/*
 * The --meth option family: long-option table fragments, value parsing,
 * cross-option validation, chemistry defaults and the help text.
 *
 * Lives apart from fastmap.cpp so every subcommand that takes these options
 * parses, validates and documents them identically.
 */
#ifndef BWA_MEM3_METH_OPTS_H
#define BWA_MEM3_METH_OPTS_H

#include <getopt.h>
#include <stdio.h>

#include "bwamem.h"

/* getopt_long ids for the family. A range of its own, clear of every
 * subcommand's local enum (those start at 1000). */
enum {
    OPT_METH = 1900,
    OPT_METH_CONTEXT,
    OPT_METH_SCORING,
    OPT_METH_SEED_PRUNE,
    OPT_METH_TAGS,
    OPT_METH_SET_AS_FAILED,
    OPT_METH_CHIMERA_QC
};

/* Where a methylation conversion is free (--meth-context). Stored in
 * mem_opt_t::meth_context as an int. */
enum meth_ctx_t { METH_CTX_ALL = 0, METH_CTX_CPG = 1 };

/* Which family options the user gave explicitly: a default must not override
 * them. Zero-initialise before parsing. */
struct meth_cli_t {
    int set_meth, set_scoring, set_context, set_seed_prune, set_tags;
};

/* If opt_id belongs to the family, apply optarg to *opt, record it in *cli and
 * return 1; return 0 if it does not. A bad value prints its error and exits 1. */
int meth_opts_parse(int opt_id, const char *optarg, mem_opt_t *opt, meth_cli_t *cli);

/* The family's cross-option checks, judged as meth_opts_finalize judges them:
 * the first error's message, without an "ERROR: " prefix, or NULL when the
 * options are consistent. It does not print, exit or change *opt, so a program
 * embedding bwa-mem3 can report the error itself and call meth_opts_finalize
 * only on NULL. The returned string is static.
 *
 * Call it BEFORE meth_opts_finalize: it judges the scoring mode finalize would
 * resolve (TAPS without cli->set_scoring counts as NEUTRAL, whatever
 * opt->meth_scoring holds), so after finalize it is accurate only if cli still
 * records every scoring change made since. With meth_mode set, it also rejects
 * a meth_context, meth_chem or meth_scoring value outside its enum (the parser
 * never produces one). That range check reads the stored values, before the
 * scoring is resolved, so an out-of-range meth_scoring is rejected even under
 * TAPS without cli->set_scoring, whose resolved scoring would be NEUTRAL. With
 * meth_mode unset, it checks only "--meth-context requires --meth" (a
 * meth_context other than METH_CTX_ALL, or cli->set_context) and does not judge
 * the stored meth_chem or meth_scoring. The other family fields (seed prune,
 * tags, --set-as-failed) are validated only by meth_opts_parse. Not an error,
 * so not reported here: the em-seq + neutral scoring warning finalize prints
 * (meth_opts_emseq_neutral tests for it). */
const char *meth_opts_check(const mem_opt_t *opt, const meth_cli_t *cli);

/* Cross-option checks (meth_opts_check; exit 1 with its error on a bad
 * combination) and the chemistry defaults (TAPS scores NEUTRAL,
 * --meth-seed-prune defaults to spec30). Prints a warning to stderr, and
 * continues, when meth_opts_emseq_neutral holds. Call once, after every option
 * is parsed and before mem_opt_apply_meth_defaults. */
void meth_opts_finalize(mem_opt_t *opt, const meth_cli_t *cli);

/* 1 when the options score em-seq's dense conversions as 0 (--meth=emseq with
 * --meth-scoring neutral), which places short reads poorly; meth_opts_finalize
 * prints a warning for it. */
int meth_opts_emseq_neutral(const mem_opt_t *opt);

/* Help text: the options every subcommand shares, then the ones only `mem` has. */
void meth_opts_usage_shared(FILE *fp);
void meth_opts_usage_mem_only(FILE *fp);

/* The closed-vocabulary values of the family's options, NULL-terminated. A
 * positional that equals one of them (and names no file) is a value orphaned
 * from its flag. */
const char *const *meth_opts_stray_values(void);

/* The family flag that owns `s` (case-insensitive), or NULL. When `canon` is
 * non-NULL it receives the vocabulary's own spelling, or NULL for the
 * `^`-prefixed --meth-tags exclusion, which has none. */
const char *meth_opts_stray_flag(const char *s, const char **canon);

/* struct option initializer fragments (each entry ends in a comma). */
#define METH_LONG_OPTIONS_SHARED \
    {"meth",         optional_argument, 0, OPT_METH}, \
    {"meth-context", required_argument, 0, OPT_METH_CONTEXT}, \
    {"meth-scoring", required_argument, 0, OPT_METH_SCORING}, \
    {"meth-tags",    required_argument, 0, OPT_METH_TAGS},

#define METH_LONG_OPTIONS_MEM_ONLY \
    {"meth-seed-prune", optional_argument, 0, OPT_METH_SEED_PRUNE}, \
    {"set-as-failed",   required_argument, 0, OPT_METH_SET_AS_FAILED}, \
    {"chimera-qc",      no_argument,       0, OPT_METH_CHIMERA_QC},

#endif /* BWA_MEM3_METH_OPTS_H */
