/* Test-side conveniences over src/rescue_prune.h, kept out of the production header because only
 * the unit tests and the band harness use them: the default-scoring overloads of the filter entry
 * points (production always passes a rescue_prune_params) and a default-scoring predicate. */
#ifndef BWA_MEM3_TEST_RESCUE_PRUNE_TEST_H
#define BWA_MEM3_TEST_RESCUE_PRUNE_TEST_H

#include "rescue_prune.h"

/* Whether p is the default scoring (-A 1 -B 4 -O 6 -E 1) without a --meth conversion or relation. */
static inline bool rescue_prune_is_default_scoring(const rescue_prune_params &p)
{
    return p.conv_from < 0 && p.relx < 0 && p.a == 1 && p.b == 4 && p.o_del == 6 && p.e_del == 1 && p.o_ins == 6
           && p.e_ins == 1;
}

/* rescue_prune_window / rescue_prune_window_scalar at the default scoring and threshold minsc. */
static inline int rescue_prune_window(const uint8_t *ref, int len1, const uint8_t *q, int len2, int minsc,
                                      int max_hits, int *hb, int *he, rescue_prune_view *view = nullptr)
{
    return rescue_prune_window(ref, len1, q, len2, rescue_prune_params::defaults(minsc), max_hits, hb, he,
                               view);
}
static inline int rescue_prune_window_scalar(const uint8_t *ref, int len1, const uint8_t *q, int len2,
                                             int minsc, int max_hits, rescue_prune_scratch &s, int *hb,
                                             int *he)
{
    return rescue_prune_window_scalar(ref, len1, q, len2, rescue_prune_params::defaults(minsc), max_hits,
                                      s, hb, he);
}

#endif
