// Shared synthetic reference for the aln port's tests.
//
// Three contigs chosen to cover what an FM-index adapter can trip over:
//   c1  600 bp of deterministic pseudo-random ACGT
//   c2  c1[100:400] reverse-complemented, 50 N, then 200 bp of random: a
//       reverse-complement copy and a run of ambiguous bases
//   c3  64 x "AC": a tandem repeat (long runs of equal suffixes)
// The sequence is a pure function of the LCG below, so it is identical on
// every host.

#ifndef BWA_MEM3_TEST_ALN_TEST_REF_H
#define BWA_MEM3_TEST_ALN_TEST_REF_H

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

namespace aln_test {

inline std::string lcg_dna(uint32_t &state, int n)
{
    static const char kBases[] = "ACGT";
    std::string s;
    for (int i = 0; i < n; ++i) {
        state = state * 1103515245u + 12345u;
        s += kBases[(state >> 16) & 3];
    }
    return s;
}

inline std::string revcomp(const std::string &s)
{
    std::string r(s.rbegin(), s.rend());
    for (char &ch : r) {
        switch (ch) {
        case 'A': ch = 'T'; break;
        case 'C': ch = 'G'; break;
        case 'G': ch = 'C'; break;
        case 'T': ch = 'A'; break;
        default: break;
        }
    }
    return r;
}

// A mkdtemp/mkstemp template under $TMPDIR (or /tmp), as the other unit tests do: CI and
// sandboxed hosts point TMPDIR somewhere writable. `stem` is e.g. "aln_se_unit".
inline std::string tmp_template(const char *stem)
{
    const char *d = std::getenv("TMPDIR");
    return std::string(d && *d ? d : "/tmp") + "/" + stem + "_XXXXXX";
}

// (name, sequence) for c1, c2, c3.
inline std::vector<std::pair<std::string, std::string>> test_contigs()
{
    uint32_t state = 11;
    std::string c1 = lcg_dna(state, 600);
    std::string c2 = revcomp(c1.substr(100, 300)) + std::string(50, 'N') + lcg_dna(state, 200);
    std::string c3;
    for (int i = 0; i < 64; ++i) c3 += "AC";
    return {{"c1", c1}, {"c2", c2}, {"c3", c3}};
}

// Writes `contigs` as a FASTA (60 columns); returns false on I/O failure.
inline bool write_test_fasta(const std::string &path, const std::vector<std::pair<std::string, std::string>> &contigs)
{
    FILE *fp = std::fopen(path.c_str(), "w");
    if (!fp) return false;
    for (const auto &c : contigs) {
        std::fprintf(fp, ">%s\n", c.first.c_str());
        for (size_t i = 0; i < c.second.size(); i += 60)
            std::fprintf(fp, "%s\n", c.second.substr(i, 60).c_str());
    }
    return std::fclose(fp) == 0;
}

// The test contigs (test_contigs()) as a FASTA.
inline bool write_test_fasta(const std::string &path) { return write_test_fasta(path, test_contigs()); }

} // namespace aln_test

#endif
