/* SPDX-License-Identifier: MIT */

#include "meth_opts.h"

#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "meth_xm.h"

/* Every --meth=CHEM spelling and the chemistry it selects. The one list both
 * the --meth parser and the stray-value check read, so a new alias cannot be
 * accepted by one and missed by the other. 5base/5-base (Illumina 5-base)
 * convert only 5mC to T, directionally: the same aligner-visible change, call
 * polarity and sparsity as TAPS, so they are TAPS here. The one chemistry
 * difference, 5hmC (TAPS converts it, 5-base does not), changes what a call
 * measures, not how a read aligns. */
static const struct { const char *name; meth_chem_t chem; } meth_chem_names[] = {
    { "emseq",  METH_CHEM_EMSEQ }, { "em-seq", METH_CHEM_EMSEQ }, { "bisulfite", METH_CHEM_EMSEQ },
    { "taps",   METH_CHEM_TAPS  }, { "5base",  METH_CHEM_TAPS  }, { "5-base",    METH_CHEM_TAPS  },
};

/* Closed-vocabulary values of the family's other options and the flag each
 * belongs to. `all` is also a --meth-context value; it is listed once, under
 * --meth-tags, which is where an orphaned `all` has always been attributed. */
static const struct { const char *value, *flag; } meth_stray_known[] = {
    { "collapsed", "--meth-scoring" }, { "genomic",   "--meth-scoring" },
    { "neutral",   "--meth-scoring" },
    { "cpg",       "--meth-context" },
    { "XR",        "--meth-tags"    }, { "XG",        "--meth-tags" },
    { "XM",        "--meth-tags"    }, { "all",       "--meth-tags" },
    { "none",      "--meth-tags"    },
    { "spec30",    "--meth-seed-prune" }, { "baseline",  "--meth-seed-prune" },
    { "off",       "--meth-seed-prune" },
};

#define N_CHEM  (sizeof(meth_chem_names) / sizeof(meth_chem_names[0]))
#define N_KNOWN (sizeof(meth_stray_known) / sizeof(meth_stray_known[0]))

const char *const *meth_opts_stray_values(void)
{
    static const char *values[N_CHEM + N_KNOWN + 1];
    if (values[0] == NULL) {
        size_t n = 0;
        for (size_t i = 0; i < N_CHEM; ++i) values[n++] = meth_chem_names[i].name;
        for (size_t i = 0; i < N_KNOWN; ++i) values[n++] = meth_stray_known[i].value;
        values[n] = NULL;
    }
    return values;
}

const char *meth_opts_stray_flag(const char *s, const char **canon)
{
    if (canon != NULL) *canon = NULL;
    if (s == NULL) return NULL;
    /* A `^`-prefixed exclusion can only have come from --meth-tags. Its `-XM`
     * synonym needs no case here: `X` takes an argument in the short optstring,
     * so getopt binds `-XM` as `-X M` and it never reaches a positional slot. */
    if (s[0] == '^') return "--meth-tags";
    for (size_t i = 0; i < N_CHEM; ++i)
        if (strcasecmp(s, meth_chem_names[i].name) == 0) {
            if (canon != NULL) *canon = meth_chem_names[i].name;
            return "--meth";
        }
    for (size_t i = 0; i < N_KNOWN; ++i)
        if (strcasecmp(s, meth_stray_known[i].value) == 0) {
            if (canon != NULL) *canon = meth_stray_known[i].value;
            return meth_stray_known[i].flag;
        }
    return NULL;
}

int meth_opts_parse(int opt_id, const char *optarg, mem_opt_t *opt, meth_cli_t *cli)
{
    switch (opt_id) {
    case OPT_METH:
        opt->meth_mode = 1;
        cli->set_meth = 1;
        /* --meth selects bisulfite ALIGNMENT semantics only; the output
         * container is --bam's job here exactly as it is without --meth.
         * (Through 0.7.x this set bam_mode=1, which made text SAM
         * unreachable under --meth — there was no flag that undid it.) */
        /* Optional chemistry argument. getopt_long's optional_argument only
         * accepts the `--meth=taps` form (a separate word is treated as a
         * positional), so a bare `--meth` leaves optarg NULL => em-seq.
         * Set em-seq first so the last --meth wins: a bare --meth after
         * --meth=taps is em-seq, not a leftover taps. */
        opt->meth_chem = METH_CHEM_EMSEQ;
        if (optarg != NULL) {
            size_t k = 0;
            while (k < N_CHEM && strcmp(optarg, meth_chem_names[k].name) != 0) ++k;
            if (k == N_CHEM) {
                fprintf(stderr, "ERROR: --meth accepts 'emseq' (default; aliases em-seq, bisulfite), "
                                "'taps', or '5base' (alias 5-base), got '%s'\n", optarg);
                exit(1);
            }
            opt->meth_chem = meth_chem_names[k].chem;
        }
        return 1;
    case OPT_METH_CONTEXT:
        if (optarg != NULL && strcmp(optarg, "all") == 0) {
            opt->meth_context = METH_CTX_ALL;
        } else if (optarg != NULL && strcmp(optarg, "cpg") == 0) {
            opt->meth_context = METH_CTX_CPG;
        } else {
            fprintf(stderr, "ERROR: --meth-context accepts 'all' (default) or 'cpg', got '%s'\n",
                    optarg != NULL ? optarg : "");
            exit(1);
        }
        cli->set_context = 1;
        return 1;
    case OPT_METH_SCORING:
        if (optarg != NULL && strcmp(optarg, "collapsed") == 0) {
            opt->meth_scoring = MEM_METH_SCORING_COLLAPSED;
        } else if (optarg != NULL && strcmp(optarg, "genomic") == 0) {
            opt->meth_scoring = MEM_METH_SCORING_GENOMIC;
        } else if (optarg != NULL && strcmp(optarg, "neutral") == 0) {
            opt->meth_scoring = MEM_METH_SCORING_NEUTRAL;
        } else {
            fprintf(stderr, "ERROR: --meth-scoring requires 'collapsed', 'genomic', or 'neutral'\n");
            exit(1);
        }
        cli->set_scoring = 1;
        return 1;
    case OPT_METH_SEED_PRUNE:
        /* --meth-seed-prune[=spec30|baseline|off]: prune the 3-letter
         * over-seeding before SA resolution (~30% faster --meth at ~0 truth
         * accuracy cost; NOT byte-identical). ON by default under --meth
         * (SPEC30, set in meth_opts_finalize); pass =off to restore
         * byte-identical seeding. getopt_long's optional_argument only binds
         * with '=' (a bare word is a positional), so a plain --meth-seed-prune
         * leaves optarg NULL => spec30 (the recommended rule). */
        if (optarg == NULL || strcmp(optarg, "spec30") == 0) {
            opt->meth_seed_prune = MEM_METH_PRUNE_SPEC30;
        } else if (strcmp(optarg, "baseline") == 0) {
            opt->meth_seed_prune = MEM_METH_PRUNE_BASELINE;
        } else if (strcmp(optarg, "off") == 0) {
            opt->meth_seed_prune = MEM_METH_PRUNE_OFF;
        } else {
            /* optarg is non-NULL here (NULL took the spec30 branch), so the
             * value was bound with '=' -- the orphaned-space case cannot
             * reach this branch, and a "use =, not a space" note would
             * misdirect. Just name the accepted vocabulary. */
            fprintf(stderr, "ERROR: --meth-seed-prune accepts 'spec30' (default), 'baseline', or 'off'\n");
            exit(1);
        }
        cli->set_seed_prune = 1;   /* explicit: wins over the --meth spec30 default */
        return 1;
    case OPT_METH_TAGS: {
        const char *tag_err = NULL;
        if (mem_opt_parse_meth_tags(optarg, &opt->meth_tags, &tag_err) != 0) {
            fprintf(stderr, "ERROR: --meth-tags '%s': %s\n"
                            "       expected 'all', 'none', a comma-separated list "
                            "(e.g. XR,XG), or ^-prefixed exclusions (e.g. ^XM)\n",
                    optarg != NULL ? optarg : "", tag_err);
            exit(1);
        }
        cli->set_tags = 1;
        return 1;
    }
    case OPT_METH_SET_AS_FAILED:
        if (optarg == NULL || !(optarg[0] == 'f' || optarg[0] == 'r') || optarg[1] != '\0') {
            fprintf(stderr, "ERROR: --set-as-failed requires 'f' or 'r'\n");
            exit(1);
        }
        opt->meth_set_as_failed = optarg[0];
        return 1;
    case OPT_METH_CHIMERA_QC:
        opt->meth_chimera_qc = 1;
        return 1;
    default:
        return 0;
    }
}

void meth_opts_finalize(mem_opt_t *opt, const meth_cli_t *cli)
{
    if (cli->set_context && !opt->meth_mode) {
        fprintf(stderr, "ERROR: --meth-context requires --meth\n");
        exit(1);
    }
    if (!opt->meth_mode) return;
    /* TAPS (and 5base, which is TAPS here) defaults to NEUTRAL scoring.
     * TAPS converts only the METHYLATED cytosines, so conversions are ~20-30x rarer than under em-seq (measured:
     * 3.2% of C vs 94.6% on chr22 @ 12x). The collapsed 3-letter alphabet buys
     * little at that density and costs specificity. NEUTRAL goes one step
     * further: it scores the conversion cell as 0 rather than a full match,
     * so a sparse TAPS conversion is tolerated without over-crediting spurious
     * C->T alignments. Measured on 4.07M simulated TAPS reads vs full hg38,
     * across three methylation loads: NEUTRAL placed 95.95-96.01% (essentially
     * the 96.02% unconverted ceiling) vs GENOMIC 95.68-95.73% and COLLAPSED
     * 95.37-95.43% -- a robust +0.24-0.28 pp over genomic -- with IDENTICAL
     * MAPQ calibration (60+ bin 99.99% both) and NM (1.966 vs 1.969), so the
     * gain is not bought with over-confidence or inflated edit distance.
     * NEUTRAL's freed cell is not rank-1, but the generalized kswv freed-cell
     * blend scores it to its matrix value, so mate rescue stays on the batched
     * kernel on the freed-capable tiers (NEON/AVX2/AVX512BW) and falls back to
     * scalar ksw_align2 only on the freed-less x86 tiers (sse41/sse42/avx),
     * exactly as GENOMIC and COLLAPSED do. An explicit --meth-scoring still wins.
     * Set before mem_opt_apply_meth_defaults so its COLLAPSED -B 2 branch keys
     * off the resolved scoring mode. */
    if (opt->meth_chem == METH_CHEM_TAPS && !cli->set_scoring)
        opt->meth_scoring = MEM_METH_SCORING_NEUTRAL;
    /* --meth-seed-prune defaults to spec30 under --meth. Unlike non-meth, --meth
     * is a bwa-mem3-native feature -- bwa/bwa-mem2 never had a bisulfite mode, so
     * there is no upstream reference output to stay byte-identical to. spec30 is
     * ~30% faster on em-seq at ~0 truth-based accuracy cost (mapq>=20 accuracy
     * identical; 8x fewer newly-unmapped than the `baseline` rule). An explicit
     * --meth-seed-prune=off|baseline (or the BWAMEM3_METH_SEED_PRUNE env) still
     * wins; no effect outside --meth (the prune is gated on the meth remap). */
    if (!cli->set_seed_prune)
        opt->meth_seed_prune = MEM_METH_PRUNE_SPEC30;

    if (opt->meth_context == METH_CTX_CPG) {
        if (opt->meth_chem == METH_CHEM_EMSEQ) {
            fprintf(stderr, "ERROR: --meth-context cpg requires --meth=taps or --meth=5base: under "
                            "EM-seq/bisulfite every unmethylated C converts, so conversions are not "
                            "confined to CpG\n");
            exit(1);
        }
        if (opt->meth_scoring == MEM_METH_SCORING_COLLAPSED) {
            fprintf(stderr, "ERROR: --meth-context cpg is incompatible with --meth-scoring collapsed: "
                            "collapsed also frees reference T against read C, which has no CpG meaning\n");
            exit(1);
        }
    }
}

void meth_opts_usage_shared(FILE *fp)
{
    fprintf(fp, "   --meth[=CHEM] enable inline bwameth-style C→T/G→A read conversion + meth-aware\n");
    fprintf(fp, "                 record emission (XM:Z/XG:Z/XR:Z). Output is SAM text by default,\n");
    fprintf(fp, "                 as without --meth; add --bam for BAM. Requires the reference to\n");
    fprintf(fp, "                 have been built with `bwa-mem3 index --meth` (emits the original\n");
    fprintf(fp, "                 index plus a ref.fa.meth.* converted seed index).\n");
    fprintf(fp, "                 CHEM selects the chemistry and thus the XM:Z call polarity:\n");
    fprintf(fp, "                   emseq  bisulfite/EM-seq, UNmethylated C converts [default]\n");
    fprintf(fp, "                          (aliases: em-seq, bisulfite)\n");
    fprintf(fp, "                   taps   TET-assisted pyridine borane, METHYLATED C converts;\n");
    fprintf(fp, "                          also defaults --meth-scoring to neutral\n");
    fprintf(fp, "                   5base  Illumina 5-base, METHYLATED C (5mC) converts; same\n");
    fprintf(fp, "                          aligner semantics as taps (alias: 5-base)\n");
    fprintf(fp, "                 NOTE: use --meth=taps (with '='), not --meth taps. Running TAPS\n");
    fprintf(fp, "                 or 5-base data without =taps or =5base inverts every\n");
    fprintf(fp, "                 methylation call.\n");
    fprintf(fp, "   --meth-context all|cpg\n");
    fprintf(fp, "                 where a conversion is free. all (default): at every C, as the\n");
    fprintf(fp, "                 chemistry converts it. cpg: only at a CpG, so a C->T at a\n");
    fprintf(fp, "                 non-CpG C scores as a mismatch. Needs --meth=taps or --meth=5base\n");
    fprintf(fp, "                 (EM-seq/bisulfite convert every unmethylated C) and a --meth-scoring\n");
    fprintf(fp, "                 other than collapsed.\n");
    fprintf(fp, "   --meth-scoring collapsed|genomic|neutral\n");
    fprintf(fp, "                 scoring mode. Default depends on chemistry: collapsed for\n");
    fprintf(fp, "                 --meth/--meth=emseq, neutral for --meth=taps|5base (their\n");
    fprintf(fp, "                 conversions are sparse, so collapsing costs specificity it\n");
    fprintf(fp, "                 can't repay).\n");
    fprintf(fp, "                 collapsed: C/T (and G/A) interchangeable, bwameth-compatible\n");
    fprintf(fp, "                 placement (sets -B 2).\n");
    fprintf(fp, "                 genomic: free only the conversion direction, scored as a full\n");
    fprintf(fp, "                 match, keep variants as mismatches (variant-aware: variants\n");
    fprintf(fp, "                 outside the conversion direction counted in NM; -B 4).\n");
    fprintf(fp, "                 neutral: free only the conversion direction but score it 0\n");
    fprintf(fp, "                 (tolerated, not rewarded); best for TAPS (variant-aware:\n");
    fprintf(fp, "                 variants counted in NM as in genomic; -B 4).\n");
    fprintf(fp, "                 In genomic/neutral a real variant in the conversion direction\n");
    fprintf(fp, "                 itself (C->T at a reference C) is indistinguishable from a\n");
    fprintf(fp, "                 conversion and is not counted in NM; under --meth-context cpg\n");
    fprintf(fp, "                 only at a CpG (a C->T at any other C counts).\n");
    fprintf(fp, "   --meth-tags SPEC\n");
    fprintf(fp, "                 which Bismark tags to emit: 'all' (default), 'none', a\n");
    fprintf(fp, "                 comma-separated list (XR,XG), or ^-prefixed exclusions (^XM).\n");
    fprintf(fp, "                 Comma-separated, NOT space-separated (a second word becomes\n");
    fprintf(fp, "                 the reference). An exclusion may be written ^XM or -XM; prefer\n");
    fprintf(fp, "                 -XM in scripts (bare ^ is a negated glob in zsh EXTENDED_GLOB).\n");
    fprintf(fp, "                 Unselected tags are not computed. XM:Z is a read-length string\n");
    fprintf(fp, "                 and dominates the BAM's aux payload; '^XM' drops it for callers\n");
    fprintf(fp, "                 that recompute calls from the reference and the read.\n");
    fprintf(fp, "                 Keep XM for Bismark-family tools (bismark_methylation_extractor)\n");
    fprintf(fp, "                 and other callers that read the per-base calls from XM:Z.\n");
}

void meth_opts_usage_mem_only(FILE *fp)
{
    fprintf(fp, "   --meth-seed-prune[=spec30|baseline|off]\n");
    fprintf(fp, "                 prune the 3-letter alphabet's short, repetitive spurious SMEMs\n");
    fprintf(fp, "                 before SA resolution: ~30%% faster --meth at ~0 accuracy cost\n");
    fprintf(fp, "                 vs truth (measured WGS/WES/em-seq, mapq>=20 accuracy identical).\n");
    fprintf(fp, "                 ON by default under --meth (spec30) -- --meth is bwa-mem3-native,\n");
    fprintf(fp, "                 so there is no upstream reference output to preserve. spec30:\n");
    fprintf(fp, "                 per-read two-regime rule; baseline: drop len<25 & SA-count>1;\n");
    fprintf(fp, "                 off: restore the un-pruned seeding (byte-identical). NOT\n");
    fprintf(fp, "                 byte-identical vs off. Env BWAMEM3_METH_SEED_PRUNE (off|spec30|\n");
    fprintf(fp, "                 baseline) overrides this flag for A/B testing.\n");
    fprintf(fp, "   --set-as-failed f|r\n");
    fprintf(fp, "                 flag alignments to the matching strand ('f' or 'r') as QC-fail (0x200)\n");
    fprintf(fp, "   --chimera-qc\n");
    fprintf(fp, "                 enable the bwameth.py-style longest-match <44%% chimera heuristic\n");
    fprintf(fp, "                 (sets 0x200, clears 0x2, caps MAPQ at 1). Off by default; not in Bismark.\n");
}
