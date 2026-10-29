/* Shared by the lockstep parity tests and seed_jump_table_test: reads around
 * the seed jump table's depths. */
#ifndef SEED_JUMP_TEST_UTIL_H
#define SEED_JUMP_TEST_UTIL_H

#include <algorithm>
#include <string>
#include <vector>

/* Reads from a phiX window: every length 1..14, a 14-base read with an N at
 * each position, and a reverse-complement window. A walk's table run is cut
 * short by the read end or the N at every length from 1 base up past the
 * default depth. */
static inline std::vector<std::string> jump_depth_reads()
{
    const std::string src = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAA";
    std::vector<std::string> v;
    for (size_t L = 1; L <= 14; L++) v.push_back(src.substr(5, L));
    for (size_t p = 0; p < 14; p++) {
        std::string r = src.substr(30, 14);
        r[p] = 'N';
        v.push_back(r);
    }
    std::string rc = src.substr(10, 14);
    std::reverse(rc.begin(), rc.end());
    for (char &c : rc) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
    v.push_back(rc);
    return v;
}

#endif
