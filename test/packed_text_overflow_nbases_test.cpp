// test/packed_text_overflow_nbases_test.cpp -- the PackedText ctor must reject
// an oversized n_bases (one for which n_bases + 3 overflows int64_t) with a
// dedicated, accurate diagnostic *before* computing expected_bytes, instead of
// relying on the incidental (and technically undefined-behavior) wraparound
// of the (n_bases + 3) >> 2 expression to trip the "too small" size check.
//
// The correct behavior is a fatal abort naming the oversized n_bases; this
// forks: the child constructs a PackedText with n_bases = INT64_MAX and must
// die with a message that says the length is too large. A regression (guard
// removed) still aborts on the same input -- because the overflowed
// expected_bytes happens to be enormous on this compiler/optimization level
// -- but with the generic "too small" message instead, which the parent
// detects as the failure (the guard is not actually preventing the signed
// overflow, it is masked by an accident of the current codegen).

#include "packed_text.h"

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/wait.h>
#include <unistd.h>

int main() {
    char path[] = "/tmp/pt_overflow_nbases_XXXXXX";
    int fd = mkstemp(path);
    if (fd < 0) { perror("mkstemp"); return 1; }
    const char pac[8] = {0};
    if (write(fd, pac, sizeof(pac)) != (ssize_t)sizeof(pac)) { perror("write"); return 1; }
    close(fd);

    int pipefd[2];
    if (pipe(pipefd) != 0) { perror("pipe"); return 1; }
    pid_t pid = fork();
    if (pid < 0) { perror("fork"); return 1; }
    if (pid == 0) {
        dup2(pipefd[1], STDERR_FILENO);
        close(pipefd[0]); close(pipefd[1]);
        PackedText pt(path, INT64_MAX);   // n_bases + 3 overflows int64_t
        _exit(0);
    }
    close(pipefd[1]);
    char buf[4096]; size_t n = 0; ssize_t r;
    while (n < sizeof(buf) - 1 && (r = read(pipefd[0], buf + n, sizeof(buf) - 1 - n)) > 0)
        n += (size_t)r;
    buf[n] = '\0'; close(pipefd[0]);
    int status = 0; waitpid(pid, &status, 0);
    unlink(path);
    bool aborted = !WIFEXITED(status) || WEXITSTATUS(status) != 0;
    bool has_msg = strstr(buf, "too large") != nullptr;
    if (aborted && has_msg) { fprintf(stderr, "packed_text_overflow_nbases_test OK\n"); return 0; }
    fprintf(stderr, "FAIL: aborted=%d has_msg=%d (overflow guard did not fire as expected)\n"
                    "  child stderr: %s\n", aborted, has_msg, buf);
    return 1;
}
