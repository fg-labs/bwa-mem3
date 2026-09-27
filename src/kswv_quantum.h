/* Query-padding quantum of the batched mate-rescue kernels (kswv).
 *
 * A lane's query occupies columns [0, len2) and the batch wrapper pads it out to
 * the quantum below (see the query-padding contract in kswv.cpp). The quantum is
 * the 8-bit SSE lane width (16) for the 8-bit kernels and the 16-bit one (8) for
 * the 16-bit kernels, matching bwa-mem2.
 *
 * It is its own header because more than the kernels depend on it: the exact
 * rescue-pruning filter (rescue_prune.h) bounds how far past its last hit an
 * alignment can extend by the pad columns, which score like matches, so its hull
 * must use at least the quantum the kernel pads to. It uses the 8-bit quantum on
 * jobs of both widths, which is never less than the 16-bit one. Both read it from
 * here. */
#ifndef BWA_MEM3_KSWV_QUANTUM_H
#define BWA_MEM3_KSWV_QUANTUM_H

static inline int kswv_query_quantum8(int len2)  { return ((len2 + 16 - 1) / 16) * 16; }
static inline int kswv_query_quantum16(int len2) { return ((len2 + 8 - 1) / 8) * 8; }

#endif
