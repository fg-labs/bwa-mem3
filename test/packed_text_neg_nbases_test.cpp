// test/packed_text_neg_nbases_test.cpp -- the PackedText ctor must fail closed
// (err_fatal "negative n_bases") on a negative n_bases before the
// (n_bases + 3) >> 2 expected-size computation, instead of proceeding with a
// nonsensical size.
//
// The correct behavior is a fatal abort, so this forks: the child constructs a
// PackedText with n_bases = -1 and must die with that message; the parent
// verifies. A regression (guard removed) computes expected_bytes = 0, passes
// the size check, and constructs successfully -- a clean exit(0) with no
// message -- which the parent detects as the failure.

#include "packed_text.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/wait.h>
#include <unistd.h>

int main() {
    char path[] = "/tmp/pt_neg_nbases_XXXXXX";
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
        PackedText pt(path, -1);   // must err_fatal and abort
        _exit(0);                  // reached only if the guard is gone
    }
    close(pipefd[1]);
    char buf[4096]; size_t n = 0; ssize_t r;
    while (n < sizeof(buf) - 1 && (r = read(pipefd[0], buf + n, sizeof(buf) - 1 - n)) > 0)
        n += (size_t)r;
    buf[n] = '\0'; close(pipefd[0]);
    int status = 0; waitpid(pid, &status, 0);
    unlink(path);
    bool aborted = !WIFEXITED(status) || WEXITSTATUS(status) != 0;
    bool has_msg = strstr(buf, "negative n_bases") != nullptr;
    if (aborted && has_msg) { fprintf(stderr, "packed_text_neg_nbases_test OK\n"); return 0; }
    fprintf(stderr, "FAIL: aborted=%d has_msg=%d (guard did not fire as expected)\n"
                    "  child stderr: %s\n", aborted, has_msg, buf);
    return 1;
}
