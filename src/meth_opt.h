/* SPDX-License-Identifier: MIT */
/*
 * The options every --meth subcommand shares (`mem` and `aln`), their value
 * types, and their defaults (meth_opt_init): set by meth_opts_parse_shared /
 * meth_opts_finalize_shared (meth_opts.h) and read by the matrix builder
 * meth_fill_mats (meth_model.h). mem embeds the struct in mem_opt_t as `meth`;
 * mem-only --meth options (--meth-seed-prune, --chimera-qc) stay in mem_opt_t.
 * Kept free of mem's headers so a subcommand can take the options without them.
 */
#ifndef BWA_MEM3_METH_OPT_H
#define BWA_MEM3_METH_OPT_H

/* The --meth SW scoring model, selected by --meth-scoring. All
 * three free the conversion cell (OT: ref-C x read-T; OB: ref-G x read-A); they
 * differ in to what value, and whether the mirror cell is freed too:
 *   COLLAPSED: conversion + mirror both freed to +a (C/T interchangeable,
 *              bwameth-compatible placement, -B 2).
 *   GENOMIC:   conversion freed to +a, mirror kept as a mismatch (variant-aware,
 *              -B 4, rank-1 batched-expressible).
 *   NEUTRAL:   conversion freed to 0 (tolerated but not rewarded), mirror kept a
 *              mismatch (-B 4). For TAPS, whose conversions are sparse (~3% of C),
 *              a full-match reward over-credits spurious C->T alignments; scoring
 *              the conversion neutrally measured ~+0.25 pp placement over GENOMIC
 *              at every methylation load. Not rank-1 (the freed value is neither
 *              match nor mismatch), but the generalized kswv freed-cell blend
 *              scores the cell to its matrix value, so mate rescue is batched on
 *              the freed-capable tiers (NEON/AVX2/AVX512BW) exactly like
 *              GENOMIC/COLLAPSED, and falls back to scalar ksw_align2 on the
 *              freed-less x86 tiers (sse41/sse42/avx) as all three modes do. */
enum mem_meth_scoring { MEM_METH_SCORING_COLLAPSED = 0, MEM_METH_SCORING_GENOMIC = 1,
                        MEM_METH_SCORING_NEUTRAL = 2 };

/* Bismark tag selection under --meth (--meth-tags). XR:Z and XG:Z are two-byte
 * strand labels; XM:Z is a read-length methylation-call string and dominates the
 * BAM's aux payload, so it is the one worth being able to turn off. */
#define MEM_METH_TAG_XR   0x1
#define MEM_METH_TAG_XG   0x2
#define MEM_METH_TAG_XM   0x4
#define MEM_METH_TAGS_ALL (MEM_METH_TAG_XR | MEM_METH_TAG_XG | MEM_METH_TAG_XM)

/* Where a methylation conversion is free (--meth-context). */
enum meth_ctx_t { METH_CTX_ALL = 0, METH_CTX_CPG = 1 };

/* The methylation chemistry (--meth=CHEM): what a converted base means (XM:Z
 * call polarity, meth_xm.h). 5base is TAPS. */
typedef enum {
    METH_CHEM_EMSEQ = 0,  /* default: bisulfite / em-seq */
    METH_CHEM_TAPS  = 1
} meth_chem_t;

struct meth_opt_t {
    int meth_mode;    // 1 = bisulfite mode (--meth)
    int meth_scoring; // bisulfite matrix mode (--meth-scoring):
                      // MEM_METH_SCORING_{COLLAPSED,GENOMIC,NEUTRAL}
    int meth_chem;    // methylation chemistry (--meth=emseq|taps|5base): meth_chem_t.
                      // Selects the XM:Z call polarity and the defaults that
                      // follow from conversion density (--meth-scoring, the
                      // rescue-prune mode); seeding, the converted index and
                      // each mode's scoring matrix are shared, because every
                      // chemistry produces the same C->T / G->A base change
                      // as far as the aligner can see. 5base is TAPS.
    int meth_tags;    // bitmask of Bismark tags to emit (--meth-tags):
                      // MEM_METH_TAG_{XR,XG,XM}. Tags not selected are
                      // neither computed nor emitted -- clearing XM skips
                      // the per-read meth_build_xm() pass entirely.
    int meth_context; // --meth-context: where a conversion is free, a meth_ctx_t
                      // value: METH_CTX_ALL (default; every C) or
                      // METH_CTX_CPG (only at a CpG).
    char meth_set_as_failed; // --set-as-failed: 'f' (OT), 'r' (OB) or 0; records
                             // on that strand get the QC-fail flag 0x200
};

/* The defaults, in one place for every subcommand. */
static inline void meth_opt_init(struct meth_opt_t *m)
{
    m->meth_mode = 0;
    m->meth_scoring = MEM_METH_SCORING_COLLAPSED; /* bwameth-compatible */
    m->meth_chem = METH_CHEM_EMSEQ;               /* bisulfite/em-seq */
    m->meth_tags = MEM_METH_TAGS_ALL;             /* XR, XG and XM (the Bismark set) */
    m->meth_context = METH_CTX_ALL;               /* every C converts */
    m->meth_set_as_failed = 0;
}

#endif
