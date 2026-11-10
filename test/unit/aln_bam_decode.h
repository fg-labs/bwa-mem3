#ifndef BWA_TEST_ALN_BAM_DECODE_H
#define BWA_TEST_ALN_BAM_DECODE_H

#include <string>

namespace bwa_tests {
// Decode all records and require clean EOF; isolated from the aligner's kstring type.
std::string decoded_aln_records(const std::string &path);
}

#endif
