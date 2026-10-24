// test/framework/mem_fixture.h
//
// Header-only fixtures for doctest cases on the `mem` path:
//
//   - mem_synth_reference(): the synthetic reference of
//     test/fixtures/mem_synth_ref.awk, generated in memory with the same
//     Park-Miller stream, byte for byte (mem_synth_fasta() is the awk script's
//     exact output, and its POSIX cksum is pinned in kMemSynthCksum, the value
//     test/framework/mem_cli.sh pins for the awk output). Coordinates of the
//     planted repeats are in the awk script's header.
//   - nt4() / revcomp(): the 0-4 encoding `mem` aligns in (A C G T N).
//   - TempDir / TempIndex: a scratch directory under $TMPDIR (removed with
//     everything in it on destruction), and a `bwa-mem3 index` of a set of
//     contigs built there and loaded back through FMI_search, for tests that
//     need a real FM index, bns and pac. Building and loading the synthetic
//     reference takes a few milliseconds.
//
// Unit tests use these only with generated sequence; nothing is read from
// test/fixtures/.

#ifndef BWA_TESTS_MEM_FIXTURE_H
#define BWA_TESTS_MEM_FIXTURE_H

#include <dirent.h>
#include <fcntl.h>
#include <unistd.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "FMI_search.h"
#include "bwa.h"

namespace bwa_tests {

struct Contig {
    std::string name;
    std::string seq;
};

// POSIX `cksum` of mem_synth_fasta(): {CRC, byte count}. Same value as
// MEM_CLI_REF_CKSUM in test/framework/mem_cli.sh.
constexpr uint32_t kMemSynthCksum = 1290712259u;
constexpr uint64_t kMemSynthBytes = 12835u;

namespace mem_fixture_detail {

inline std::string draw(uint64_t &x, int n) {
    std::string s((size_t)n, 'A');
    for (int i = 0; i < n; ++i) {
        x = (x * 48271u) % 2147483647u;
        s[(size_t)i] = "ACGT"[x / 536870912u];
    }
    return s;
}

inline void put(std::string &s, int pos1, const std::string &blk) {
    s.replace((size_t)(pos1 - 1), blk.size(), blk);
}

inline char next_base(char c) {
    return c == 'A' ? 'C' : c == 'C' ? 'G' : c == 'G' ? 'T' : 'A';
}

}  // namespace mem_fixture_detail

inline std::string revcomp(const std::string &s) {
    std::string r(s.rbegin(), s.rend());
    for (char &c : r)
        c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : 'N';
    return r;
}

// A C G T -> 0 1 2 3; anything else -> 4 (N), as bwa's nst_nt4_table does for
// the bases `mem` sees.
inline std::vector<uint8_t> nt4(const std::string &s) {
    std::vector<uint8_t> v(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        switch (s[i]) {
            case 'A': case 'a': v[i] = 0; break;
            case 'C': case 'c': v[i] = 1; break;
            case 'G': case 'g': v[i] = 2; break;
            case 'T': case 't': v[i] = 3; break;
            default: v[i] = 4; break;
        }
    }
    return v;
}

// The synthetic reference of test/fixtures/mem_synth_ref.awk (see its header):
// ctg1 5000 bp, ctg2 4000 bp, ctg3 3000 bp, ctg4 600 bp.
inline std::vector<Contig> mem_synth_reference() {
    using namespace mem_fixture_detail;
    uint64_t x = 20261024u;
    std::string c1 = draw(x, 5000), c2 = draw(x, 4000), c3 = draw(x, 3000), c4 = draw(x, 600);
    const std::string R = c1.substr(2000, 300), S = c1.substr(3000, 500);
    std::string Sd = S;
    for (size_t i = 1; i <= Sd.size(); ++i)
        if (i % 20 == 10) Sd[i - 1] = next_base(Sd[i - 1]);
    std::string ac;
    for (int i = 0; i < 60; ++i) ac += "AC";
    put(c2, 1001, R);
    put(c2, 2001, Sd);
    put(c3, 501, ac);
    put(c3, 1501, std::string(40, 'N'));
    put(c3, 2201, std::string(50, 'A'));
    put(c4, 101, revcomp(R));
    return {{"ctg1", c1}, {"ctg2", c2}, {"ctg3", c3}, {"ctg4", c4}};
}

// FASTA text, `width` bases per line (the awk script writes 60).
inline std::string fasta_text(const std::vector<Contig> &contigs, size_t width = 60) {
    std::string out;
    for (const Contig &c : contigs) {
        out += ">" + c.name + "\n";
        for (size_t i = 0; i < c.seq.size(); i += width) out += c.seq.substr(i, width) + "\n";
    }
    return out;
}

inline std::string mem_synth_fasta() { return fasta_text(mem_synth_reference()); }

// POSIX cksum: CRC-32 (polynomial 0x04C11DB7, MSB first) over the bytes and then
// the length's bytes, least significant first, complemented. {crc, length}.
inline std::pair<uint32_t, uint64_t> posix_cksum(const std::string &data) {
    uint32_t crc = 0;
    auto feed = [&crc](uint8_t byte) {
        crc ^= (uint32_t)byte << 24;
        for (int k = 0; k < 8; ++k) crc = (crc & 0x80000000u) ? (crc << 1) ^ 0x04C11DB7u : crc << 1;
    };
    for (unsigned char c : data) feed(c);
    for (uint64_t n = data.size(); n != 0; n >>= 8) feed((uint8_t)(n & 0xff));
    return {~crc, (uint64_t)data.size()};
}

// A fresh directory under $TMPDIR (default /tmp), removed with every regular
// file in it when the object goes out of scope. Flat: tests that need
// subdirectories remove them themselves.
class TempDir {
  public:
    TempDir() {
        const char *t = getenv("TMPDIR");
        std::string tmpl = (t != NULL && t[0] != '\0') ? t : "/tmp";
        if (tmpl.back() != '/') tmpl += '/';
        tmpl += "bwa_mem3_test.XXXXXX";
        std::vector<char> buf(tmpl.begin(), tmpl.end());
        buf.push_back('\0');
        if (mkdtemp(buf.data()) == NULL) throw std::runtime_error("mkdtemp failed for " + tmpl);
        path_ = buf.data();
    }
    ~TempDir() {
        if (DIR *d = opendir(path_.c_str())) {
            while (struct dirent *e = readdir(d)) {
                if (strcmp(e->d_name, ".") == 0 || strcmp(e->d_name, "..") == 0) continue;
                unlink((path_ + "/" + e->d_name).c_str());
            }
            closedir(d);
        }
        rmdir(path_.c_str());
    }
    TempDir(const TempDir &) = delete;
    TempDir &operator=(const TempDir &) = delete;

    const std::string &path() const { return path_; }
    std::string file(const std::string &name) const { return path_ + "/" + name; }

    void write(const std::string &name, const std::string &text) const {
        FILE *fp = fopen(file(name).c_str(), "wb");
        if (fp == NULL) throw std::runtime_error("cannot write " + file(name));
        const bool ok = fwrite(text.data(), 1, text.size(), fp) == text.size();
        if (fclose(fp) != 0 || !ok) throw std::runtime_error("short write to " + file(name));
    }

  private:
    std::string path_;
};

// Redirects stderr to /dev/null for its lifetime: the index builder and loader
// report progress there unconditionally.
class QuietStderr {
  public:
    QuietStderr() {
        fflush(stderr);
        saved_ = dup(fileno(stderr));
        const int devnull = open("/dev/null", O_WRONLY);
        if (devnull >= 0) {
            dup2(devnull, fileno(stderr));
            close(devnull);
        }
    }
    ~QuietStderr() {
        fflush(stderr);
        if (saved_ >= 0) {
            dup2(saved_, fileno(stderr));
            close(saved_);
        }
    }
    QuietStderr(const QuietStderr &) = delete;
    QuietStderr &operator=(const QuietStderr &) = delete;

  private:
    int saved_ = -1;
};

// `bwa-mem3 index` of `contigs` (written as ref.fa in a TempDir), loaded back
// with FMI_search::load_index. sa_compx is the SA sampling shift (`index -u`,
// default 3). fmi() is the FM index; bns() and pac() are its reference
// handles (pac is the 2-bit forward strand, N as a random base, as `mem` sees
// it).
class TempIndex {
  public:
    explicit TempIndex(const std::vector<Contig> &contigs, int sa_compx = 3) {
        dir_.write("ref.fa", fasta_text(contigs));
        prefix_ = dir_.file("ref.fa");
        QuietStderr quiet;
        if (bwa_idx_build(prefix_.c_str(), prefix_.c_str(), 0, sa_compx) != 0)
            throw std::runtime_error("bwa_idx_build failed for " + prefix_);
        fmi_.reset(new FMI_search(prefix_.c_str()));
        fmi_->load_index(/*load_pac=*/true, /*n_threads=*/1);
    }

    FMI_search *fmi() const { return fmi_.get(); }
    const bntseq_t *bns() const { return fmi_->idx->bns; }
    const uint8_t *pac() const { return fmi_->idx->pac; }
    const std::string &prefix() const { return prefix_; }

  private:
    TempDir dir_;
    std::string prefix_;
    std::unique_ptr<FMI_search> fmi_;
};

}  // namespace bwa_tests

#endif
