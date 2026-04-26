#ifndef BWA_INDEX_PRELUDE_H
#define BWA_INDEX_PRELUDE_H

#include <cstdint>

class PackedText;

// Parallel scan of a PackedText into `<prefix>.0123` (one byte per base,
// values 0..3). Matches bwa-mem2's existing .0123 format.
int emit_0123(const PackedText& pac, const char* prefix, int num_threads);

// Five-entry cumulative count array used by the FM index: count[0] = 0,
// count[i] = #{bases in pac with value < i} for i = 1..4. Bases with value
// > 3 are skipped. Parallelises the scan across `num_threads` stripes and
// reduces per-stripe histograms.
void compute_counts(const PackedText& pac, int64_t count[5], int num_threads = 1);

#endif
