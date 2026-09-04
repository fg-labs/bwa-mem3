// test/nst_nt4_decode_test.cpp -- pins the contract of nst_nt4_decode(), the
// helper that unified the nt4 re-encoding idiom. The load-bearing case is the
// signed-char one: a byte >= 0x80 arriving as a signed `char` must map to N
// (nt4 code 4), never be treated as a value below the guard and stored verbatim
// (the OOB-into-the-scoring-matrix bug the helper's unsigned-char parameter
// fixes at the call boundary). A regression that narrows the parameter back to
// `char` fails this on signed-char platforms.

#include "bntseq.h"

#include <cstdio>
#include <cstdlib>

#define CHECK(c) do { if (!(c)) { \
    std::fprintf(stderr, "CHECK failed at %s:%d: %s\n", __FILE__, __LINE__, #c); \
    std::abort(); } } while (0)

int main() {
    // Values below the limit pass through unchanged.
    CHECK(nst_nt4_decode(0, 4) == 0);
    CHECK(nst_nt4_decode(3, 4) == 3);
    CHECK(nst_nt4_decode(4, 5) == 4);   // N passes through when the caller uses limit 5
    CHECK(nst_nt4_decode(4, 4) == nst_nt4_table[4]);   // 4 is not < 4, so it maps via the table

    // A high-bit byte arriving as a signed char: the unsigned-char parameter
    // converts it value-mod-256 (0xC8 -> 200), so it is NOT < limit and maps to
    // N via the table -- rather than being stored as a negative value.
    char high_bit = (char)0xC8;                 // -56 on a signed-char platform
    CHECK(nst_nt4_decode(high_bit, 4) == 4);    // must be N, not a stored negative
    CHECK(nst_nt4_decode(high_bit, 5) == 4);
    CHECK(nst_nt4_decode((unsigned char)0x80, 4) == nst_nt4_table[0x80]);

    // ASCII bases map through the table (>= limit).
    CHECK(nst_nt4_decode((unsigned char)'A', 4) == nst_nt4_table[(unsigned char)'A']);
    CHECK(nst_nt4_decode((unsigned char)'N', 4) == nst_nt4_table[(unsigned char)'N']);

    std::fprintf(stderr, "nst_nt4_decode_test OK\n");
    return 0;
}
