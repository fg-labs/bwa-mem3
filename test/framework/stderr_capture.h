// test/framework/stderr_capture.h
//
// RAII redirection of the process's stderr for unit tests: StderrCapture keeps
// what the code under test writes there so a case can assert on a diagnostic,
// and QuietStderr discards it.

#ifndef BWA_TESTS_STDERR_CAPTURE_H
#define BWA_TESTS_STDERR_CAPTURE_H

#include <fcntl.h>
#include <unistd.h>

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>

#include "temp_dir.h"

namespace bwa_tests {

// Points stderr at `fd` for its lifetime and restores it on destruction (or on
// release()). Does not take ownership of `fd`.
class StderrRedirect {
public:
    explicit StderrRedirect(int fd) {
        fflush(stderr);
        saved_ = dup(fileno(stderr));
        if (saved_ < 0 || dup2(fd, fileno(stderr)) < 0) {
            if (saved_ >= 0) close(saved_);
            throw std::runtime_error("cannot redirect stderr");
        }
    }
    ~StderrRedirect() { release(); }
    StderrRedirect(const StderrRedirect &) = delete;
    StderrRedirect &operator=(const StderrRedirect &) = delete;

    void release() {
        if (saved_ < 0) return;
        fflush(stderr);
        const int ok = dup2(saved_, fileno(stderr));
        close(saved_);
        saved_ = -1;
        if (ok < 0) abort();   // stderr would stay redirected for the rest of the binary
    }
private:
    int saved_ = -1;
};

// An fd closed on destruction; the first member of the classes below, so a
// redirect that throws in their constructors does not leak it.
class OwnedFd {
public:
    explicit OwnedFd(int fd) : fd_(fd) {}
    ~OwnedFd() {
        if (fd_ >= 0) close(fd_);
    }
    OwnedFd(const OwnedFd &) = delete;
    OwnedFd &operator=(const OwnedFd &) = delete;
    int get() const { return fd_; }
private:
    int fd_;
};

// stderr to /dev/null for the object's lifetime.
class QuietStderr {
public:
    QuietStderr() : devnull_(checked(open("/dev/null", O_WRONLY), "cannot open /dev/null")),
                    redirect_(devnull_.get()) {}
    QuietStderr(const QuietStderr &) = delete;
    QuietStderr &operator=(const QuietStderr &) = delete;
private:
    static int checked(int fd, const char *what) {
        if (fd < 0) throw std::runtime_error(what);
        return fd;
    }
    OwnedFd devnull_;
    StderrRedirect redirect_;   // declared after devnull_: restored before it closes
    friend class StderrCapture;
};

// stderr to an unlinked temporary file under temp_root() until text() is
// called, which restores stderr and returns everything written meanwhile.
class StderrCapture {
public:
    StderrCapture()
        : fd_(QuietStderr::checked(anon_temp_fd("bwa_mem3_stderr"), "cannot create a stderr capture file")),
          redirect_(fd_.get()) {}
    StderrCapture(const StderrCapture &) = delete;
    StderrCapture &operator=(const StderrCapture &) = delete;

    std::string text() {
        redirect_.release();
        std::string out;
        if (lseek(fd_.get(), 0, SEEK_SET) != 0) throw std::runtime_error("cannot rewind captured stderr");
        char buf[4096];
        ssize_t got;
        while ((got = read(fd_.get(), buf, sizeof(buf))) > 0) out.append(buf, (size_t)got);
        if (got < 0) throw std::runtime_error("cannot read captured stderr");
        return out;
    }
private:
    OwnedFd fd_;
    StderrRedirect redirect_;   // declared after fd_: restored before it closes
};

} // namespace bwa_tests

#endif // BWA_TESTS_STDERR_CAPTURE_H
