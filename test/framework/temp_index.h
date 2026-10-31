// test/framework/temp_index.h
//
// A real `bwa-mem3 index` of a few test contigs, built in a TempDir and loaded
// back through FMI_search, for unit tests that need a genuine FM index, bns
// and pac (seeding, SA lookup, reference fetches, chaining). Building and
// loading a reference of a few kilobases takes a few milliseconds. Also the
// sequence helpers such tests use to build their reads: nt4() and revcomp().

#ifndef BWA_TESTS_TEMP_INDEX_H
#define BWA_TESTS_TEMP_INDEX_H

#include <errno.h>
#include <fcntl.h>
#include <unistd.h>

#include <cstdint>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "FMI_search.h"
#include "bntseq.h"
#include "bwa.h"
#include "stderr_capture.h"
#include "temp_dir.h"

namespace bwa_tests {

struct Contig {
    std::string name;
    std::string seq;
};

// Reverse complement over A C G T in either case (case is kept); anything else
// becomes N.
inline std::string revcomp(const std::string &s) {
    std::string r(s.rbegin(), s.rend());
    for (char &c : r) {
        switch (c) {
            case 'A': c = 'T'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            case 'T': c = 'A'; break;
            case 'a': c = 't'; break;
            case 'c': c = 'g'; break;
            case 'g': c = 'c'; break;
            case 't': c = 'a'; break;
            default: c = 'N'; break;
        }
    }
    return r;
}

// The 0-4 encoding `mem` aligns in (A C G T, either case, -> 0 1 2 3; anything
// else -> 4), through bwa's own nst_nt4_table.
inline std::vector<uint8_t> nt4(const std::string &s) {
    std::vector<uint8_t> v(s.size());
    for (size_t i = 0; i < s.size(); ++i) v[i] = nst_nt4_table[(unsigned char)s[i]];
    return v;
}

// FASTA text, `width` (> 0) bases per line.
inline std::string fasta_text(const std::vector<Contig> &contigs, size_t width = 60) {
    if (width == 0) throw std::invalid_argument("fasta_text: width must be positive");
    std::string out;
    for (const Contig &c : contigs) {
        out += ">" + c.name + "\n";
        for (size_t i = 0; i < c.seq.size(); i += width) out += c.seq.substr(i, width) + "\n";
    }
    return out;
}

namespace temp_index_detail {

// While an index is being built or loaded, its log and a copy of the real
// stderr. bwa_idx_build and load_index end the process with exit() on most
// failures, which would otherwise discard the log unseen; the atexit handler
// replays it to the real stderr.
struct ActiveLog {
    std::string path;
    int real_stderr = -1;
};

inline ActiveLog &active_log() {
    static ActiveLog log;
    return log;
}

// Writes all n bytes, retrying an interrupted write and continuing after a
// partial one; gives up on any other error.
inline void write_all(int fd, const char *p, size_t n) {
    while (n > 0) {
        const ssize_t w = write(fd, p, n);
        if (w < 0 && errno == EINTR) continue;
        if (w <= 0) return;
        p += w;
        n -= (size_t)w;
    }
}

inline void replay_active_log() {
    ActiveLog &log = active_log();
    if (log.real_stderr < 0) return;
    fflush(stderr);
    const int in = open(log.path.c_str(), O_RDONLY);
    if (in >= 0) {
        static const char kHead[] = "TempIndex: index build or load exited the process; its log follows:\n";
        write_all(log.real_stderr, kHead, sizeof(kHead) - 1);
        char buf[4096];
        ssize_t got;
        while ((got = read(in, buf, sizeof(buf))) > 0) write_all(log.real_stderr, buf, (size_t)got);
        close(in);
    }
}

// Redirects stderr into `path` for its lifetime, arming replay_active_log().
class LoggedStderr {
public:
    explicit LoggedStderr(const std::string &path)
        : fd_(open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0600)), real_(dup(fileno(stderr))),
          redirect_(checked(fd_.get(), real_.get())) {
        // Construct the static log first: atexit handlers and static
        // destructors run in reverse order of registration, so the replay
        // must be registered after the log it reads, or it runs on a
        // destroyed one.
        (void)active_log();
        static const bool registered = (atexit(replay_active_log) == 0);
        if (!registered) throw std::runtime_error("cannot register the TempIndex log replay");
        active_log().path = path;
        active_log().real_stderr = real_.get();
    }
    ~LoggedStderr() { disarm(); }
    LoggedStderr(const LoggedStderr &) = delete;
    LoggedStderr &operator=(const LoggedStderr &) = delete;

    void disarm() {
        active_log().real_stderr = -1;
        redirect_.release();
    }
private:
    static int checked(int log_fd, int real_fd) {
        if (log_fd < 0 || real_fd < 0) throw std::runtime_error("cannot set up the TempIndex log");
        return log_fd;
    }
    OwnedFd fd_;
    OwnedFd real_;              // the real stderr, taken before the redirect
    StderrRedirect redirect_;   // declared last: restored before the fds close
};

} // namespace temp_index_detail

// `contigs` written as a FASTA in a TempDir, indexed with bwa_idx_build and
// loaded with FMI_search::load_index. sa_compx is the SA sampling shift
// (`index -u`, default 3). fmi() is the FM index; bns() and pac() are its
// reference handles (pac is the 2-bit forward strand, N as a random base, as
// `mem` sees it).
//
// The FASTA's basename is the TempDir's own unique leaf: load_index first
// tries to attach an index staged in shared memory (`bwa-mem3 shm`), matched
// by basename alone, so a fixed name such as ref.fa would silently load any
// staged index of that name instead of the one just built. The build's
// progress output goes to index.log in the TempDir; on a failure that returns,
// the log is part of the exception, and on one that exits, it is replayed to
// stderr.
class TempIndex {
public:
    explicit TempIndex(const std::vector<Contig> &contigs, int sa_compx = 3) {
        const std::string &dir = dir_.path();
        const std::string fasta = dir.substr(dir.rfind('/') + 1) + ".fa";
        dir_.write(fasta, fasta_text(contigs));
        prefix_ = dir_.file(fasta);
        const std::string log = dir_.file("index.log");
        temp_index_detail::LoggedStderr logged(log);
        if (bwa_idx_build(prefix_.c_str(), prefix_.c_str(), 0, sa_compx) != 0) {
            logged.disarm();
            throw std::runtime_error("bwa_idx_build failed for " + prefix_ + "; log:\n" + read_file(log));
        }
        fmi_.reset(new FMI_search(prefix_.c_str()));
        fmi_->load_index(/*load_pac=*/true, /*n_threads=*/1);
    }

    FMI_search *fmi() const { return fmi_.get(); }
    const bntseq_t *bns() const { return fmi_->idx->bns; }
    const uint8_t *pac() const { return fmi_->idx->pac; }
    const std::string &prefix() const { return prefix_; }
private:
    static std::string read_file(const std::string &path) {
        std::string out;
        if (FILE *fp = fopen(path.c_str(), "rb")) {
            char buf[4096];
            size_t got;
            while ((got = fread(buf, 1, sizeof(buf), fp)) > 0) out.append(buf, got);
            fclose(fp);
        }
        return out;
    }
    TempDir dir_;
    std::string prefix_;
    std::unique_ptr<FMI_search> fmi_;
};

} // namespace bwa_tests

#endif // BWA_TESTS_TEMP_INDEX_H
