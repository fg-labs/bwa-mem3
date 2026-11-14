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

/* Which family options the user gave explicitly: a default must not override
 * them. Zero-initialise before parsing. */
struct meth_cli_t {
    int set_meth, set_scoring, set_context, set_seed_prune, set_tags, set_failed;
};

/* If opt_id belongs to the family, apply optarg to *opt, record it in *cli and
 * return 1; return 0 if it does not. A bad value prints its error and exits 1.
 * mem's parser: the shared options (meth_opts_parse_shared) and mem's own. */
int meth_opts_parse(int opt_id, const char *optarg, mem_opt_t *opt, meth_cli_t *cli);

/* The shared options only (--meth, --meth-context, --meth-scoring, --meth-tags,
 * --set-as-failed), into *meth: what every subcommand that takes --meth parses. */
int meth_opts_parse_shared(int opt_id, const char *optarg, meth_opt_t *meth, meth_cli_t *cli);

/* --set-as-failed's value ('f' or 'r') into *out; exits 1 on anything else. */
int meth_parse_set_as_failed(const char *optarg, char *out);

/* Cross-option checks (exit 1 with the error on a bad combination) and the
 * chemistry defaults (TAPS scores NEUTRAL, --meth-seed-prune defaults to
 * spec30). Call once, after every option is parsed and before
 * mem_opt_apply_meth_defaults. */
void meth_opts_finalize(mem_opt_t *opt, const meth_cli_t *cli);

/* 1 when the options score em-seq's dense conversions as 0 (--meth=emseq with
 * --meth-scoring neutral), which places short reads poorly; meth_opts_finalize_shared
 * prints a warning for it. */
int meth_opts_emseq_neutral(const meth_opt_t *meth);
/* The shared part of meth_opts_finalize: the cross-option checks and the
 * chemistry defaults (TAPS scores NEUTRAL) on the shared options alone. */
void meth_opts_finalize_shared(meth_opt_t *meth, const meth_cli_t *cli);

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
    {"meth-tags",    required_argument, 0, OPT_METH_TAGS}, \
    {"set-as-failed", required_argument, 0, OPT_METH_SET_AS_FAILED},

#define METH_LONG_OPTIONS_MEM_ONLY \
    {"meth-seed-prune", optional_argument, 0, OPT_METH_SEED_PRUNE}, \
    {"chimera-qc",      no_argument,       0, OPT_METH_CHIMERA_QC},

#endif /* BWA_MEM3_METH_OPTS_H */
