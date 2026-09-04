// test/read_arena_overflow_test.cpp -- read_arena_dup(len == SIZE_MAX) must fail
// closed via err_fatal ("read field length overflows size_t") rather than wrap
// len+1 to 0 and memcpy SIZE_MAX bytes into a zero-length allocation.
//
// The correct behavior is a fatal abort, so this forks: the child triggers the
// guard and must die with that message; the parent verifies both. A regression
// (guard removed) instead segfaults on the memcpy -- also a nonzero exit, but
// WITHOUT the message -- so the message check is what distinguishes the two.

#include "read_arena.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/wait.h>
#include <unistd.h>

int main() {
    int pipefd[2];
    if (pipe(pipefd) != 0) { perror("pipe"); return 1; }
    pid_t pid = fork();
    if (pid < 0) { perror("fork"); return 1; }
    if (pid == 0) {
        dup2(pipefd[1], STDERR_FILENO);
        close(pipefd[0]); close(pipefd[1]);
        read_arena_t *a = read_arena_create();
        char c = 'x';
        (void) read_arena_dup(a, &c, SIZE_MAX);   // must err_fatal and abort
        _exit(0);                                 // reached only if the guard is gone
    }
    close(pipefd[1]);
    char buf[4096]; size_t n = 0; ssize_t r;
    while (n < sizeof(buf) - 1 && (r = read(pipefd[0], buf + n, sizeof(buf) - 1 - n)) > 0)
        n += (size_t)r;
    buf[n] = '\0'; close(pipefd[0]);
    int status = 0; waitpid(pid, &status, 0);
    bool aborted = !WIFEXITED(status) || WEXITSTATUS(status) != 0;
    bool has_msg = strstr(buf, "read field length overflows size_t") != nullptr;
    if (aborted && has_msg) { fprintf(stderr, "read_arena_overflow_test OK\n"); return 0; }
    fprintf(stderr, "FAIL: aborted=%d has_msg=%d (guard did not fire as expected)\n"
                    "  child stderr: %s\n", aborted, has_msg, buf);
    return 1;
}
