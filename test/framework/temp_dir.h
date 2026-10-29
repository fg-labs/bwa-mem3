// test/framework/temp_dir.h
//
// RAII scratch directory for unit tests that need files on disk (an index
// built for the test, a corrupt header, captured output).

#ifndef BWA_TESTS_TEMP_DIR_H
#define BWA_TESTS_TEMP_DIR_H

#include <errno.h>
#include <unistd.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include "temp_path.h"

namespace bwa_tests {

// The scratch root (bwa_test_temp_root: $TMPDIR, else /tmp), with a trailing
// slash.
inline std::string temp_root() {
    std::string root = bwa_test_temp_root();
    if (root.back() != '/') root += '/';
    return root;
}

// An open, already-unlinked scratch file under temp_root() (so nothing is left
// behind however the test ends), or -1 if mkstemp or the unlink fails. A
// failed unlink closes the descriptor and names the leftover file on stderr.
inline int anon_temp_fd(const std::string &prefix = "bwa_mem3_test") {
    const std::string tmpl = temp_root() + prefix + ".XXXXXX";
    std::vector<char> buf(tmpl.begin(), tmpl.end());
    buf.push_back('\0');
    const int fd = mkstemp(buf.data());
    if (fd < 0) return -1;
    if (unlink(buf.data()) != 0) {
        const int err = errno;
        close(fd);
        std::fprintf(stderr, "anon_temp_fd: cannot unlink %s (%s); it is left behind\n",
                     buf.data(), std::strerror(err));
        return -1;
    }
    return fd;
}

// A fresh directory under temp_root(), unique per process and per object (so
// concurrent unit-test shards never collide), removed with every regular file
// in it on destruction. Flat: a test that makes subdirectories removes them;
// a directory left non-empty is reported on stderr rather than leaked silently.
class TempDir {
public:
    explicit TempDir(const std::string &prefix = "bwa_mem3_test") {
        const std::string tmpl = temp_root() + prefix + ".XXXXXX";
        std::vector<char> buf(tmpl.begin(), tmpl.end());
        buf.push_back('\0');
        if (mkdtemp(buf.data()) == NULL) throw std::runtime_error("mkdtemp failed for " + tmpl);
        path_ = buf.data();
    }
    ~TempDir() {
        if (bwa_test_remove_flat_dir(path_.c_str()) != 0)
            fprintf(stderr, "warning: TempDir could not remove %s: %s\n", path_.c_str(), strerror(errno));
    }
    TempDir(const TempDir &) = delete;
    TempDir &operator=(const TempDir &) = delete;

    const std::string &path() const { return path_; }
    std::string file(const std::string &name) const { return path_ + "/" + name; }

    // Writes `text` to file(name), replacing it; throws on any short write.
    void write(const std::string &name, const std::string &text) const {
        FILE *fp = fopen(file(name).c_str(), "wb");
        if (fp == NULL) throw std::runtime_error("cannot write " + file(name));
        const bool ok = fwrite(text.data(), 1, text.size(), fp) == text.size();
        if (fclose(fp) != 0 || !ok) throw std::runtime_error("short write to " + file(name));
    }
private:
    std::string path_;
};

} // namespace bwa_tests

#endif // BWA_TESTS_TEMP_DIR_H
