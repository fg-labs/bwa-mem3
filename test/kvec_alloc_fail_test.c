/* test/kvec_alloc_fail_test.c -- kvec.h growth macros must abort loudly, not
 * leak-and-NULL-deref, when the backing realloc fails.
 *
 * Historically kv_push/kv_pushp/kv_resize/kv_a assigned realloc()'s result
 * straight onto (v).a. On failure realloc returns NULL: the old buffer leaks
 * and the very next element write dereferences NULL. This test forces that
 * failure and asserts the fixed contract -- a loud, controlled abort() carrying
 * an out-of-memory diagnostic; never a segfault, never silent continuation.
 *
 * Standalone (header-only; links no bwa-mem3 objects) and forked: the abort has
 * to be contained so the parent can inspect *how* the child died. A buggy
 * kvec.h kills the child with SIGSEGV; the fix kills it with SIGABRT. That
 * signal is the whole discriminator, so the parent asserts SIGABRT specifically
 * and checks the child's stderr for the diagnostic to prove the child aborted
 * through kvec's own guard rather than some unrelated fault.
 *
 * The realloc seam: test_realloc() is defined *before* `#define realloc`, so its
 * own realloc() call binds to libc; kvec.h, included after, routes every resize
 * through it. No allocator interposition and no linker --wrap, so it is portable
 * to the macOS and Linux CI rows alike.
 */
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

/* Armed by the child immediately before the growth that must fail; every other
 * allocation (the initial grow, all startup allocations) takes the real path. */
static int g_fail_next_realloc = 0;

static void *test_realloc(void *p, size_t n)
{
    if (g_fail_next_realloc) return NULL;
    return realloc(p, n);
}

#define realloc(p, n) test_realloc((p), (n))

#include "kvec.h"

/* Which growth macro the child exercises. kv_push (a do/while store) and
 * kv_pushp (a comma expression yielding a slot pointer) are the two named in
 * the issue and expand differently; kv_resize and kv_a route through the same
 * kv_realloc_or_die() call, so guarding these two guards all four. */
enum grow_macro { GROW_PUSH, GROW_PUSHP };

/* Grow a kvec until a realloc is forced to fail. On a fixed kvec.h this aborts
 * inside the third growth and never returns; on a buggy one it NULL-derefs. */
static void child_force_oom(int wfd, enum grow_macro which)
{
    /* Route this child's stderr into the parent's pipe so the parent can read
     * the diagnostic the abort path prints. */
    dup2(wfd, STDERR_FILENO);

    kvec_t(int) v;
    kv_init(v);
    if (which == GROW_PUSH) {
        kv_push(int, v, 1);        /* n=0,m=0 -> grow to m=2 (real realloc)   */
        kv_push(int, v, 2);        /* n=1<m  -> no grow, n=2                  */
        g_fail_next_realloc = 1;   /* arm the failure for the next grow       */
        kv_push(int, v, 3);        /* n=2==m -> grow -> forced NULL -> abort  */
    } else {
        *kv_pushp(int, v) = 1;     /* n=0,m=0 -> grow to m=2 (real realloc)   */
        *kv_pushp(int, v) = 2;     /* n=1<m  -> no grow, n=2                  */
        g_fail_next_realloc = 1;   /* arm the failure for the next grow       */
        *kv_pushp(int, v) = 3;     /* n=2==m -> grow -> forced NULL -> abort  */
    }

    /* Reached only if the growth path did NOT abort on the failed realloc. The
     * exit code here is immaterial: the parent's discriminator is a SIGABRT
     * carrying the OOM diagnostic, which only the guarded fix produces and which
     * no optimization level can elide (abort() is a side-effecting noreturn).
     * Unguarded code instead reaches this line (a compiler may forward the dead
     * store past the NULL it produced) or faults -- either way, not SIGABRT. */
    _exit(0);
}

/* Fork a child that grows a kvec through `which` until the armed realloc fails,
 * and assert it terminated the fixed way: SIGABRT carrying the OOM diagnostic.
 * Returns 0 on pass, 1 on failure (with a reason on stderr). */
static int run_case(const char *name, enum grow_macro which)
{
    int pipefd[2];
    if (pipe(pipefd) != 0) { perror("pipe"); return 1; }

    /* Flush before forking: a prior case's PASS line may still sit in this
     * process's stdout buffer, and the child would otherwise inherit and (on
     * abort) re-emit it, double-printing. */
    fflush(NULL);

    pid_t pid = fork();
    if (pid < 0) { perror("fork"); return 1; }
    if (pid == 0) {
        close(pipefd[0]);
        child_force_oom(pipefd[1], which);
        _exit(0);              /* unreachable on both fixed and buggy code */
    }
    close(pipefd[1]);

    /* Drain the child's stderr so a full pipe can never deadlock it. */
    char buf[512];
    size_t got = 0;
    ssize_t r;
    while (got < sizeof(buf) - 1 &&
           (r = read(pipefd[0], buf + got, sizeof(buf) - 1 - got)) > 0) {
        got += (size_t)r;
    }
    buf[got] = '\0';
    close(pipefd[0]);

    int status = 0;
    if (waitpid(pid, &status, 0) < 0) { perror("waitpid"); return 1; }

    if (!WIFSIGNALED(status)) {
        fprintf(stderr, "FAIL[%s]: child exited normally (status %d); the OOM "
                        "growth path did not abort -- leak-and-NULL-deref not "
                        "fixed\n", name, WEXITSTATUS(status));
        return 1;
    }
    if (WTERMSIG(status) != SIGABRT) {
        fprintf(stderr, "FAIL[%s]: child died from signal %d, expected SIGABRT "
                        "(%d). SIGSEGV here means the growth macro "
                        "NULL-dereferenced instead of guarding the failed "
                        "realloc.\n", name, WTERMSIG(status), SIGABRT);
        return 1;
    }
    if (strstr(buf, "out of memory") == NULL) {
        fprintf(stderr, "FAIL[%s]: child aborted (SIGABRT) but without the "
                        "expected out-of-memory diagnostic; its stderr was:\n%s\n",
                        name, buf);
        return 1;
    }

    printf("PASS[%s]: kvec growth aborts with an OOM diagnostic when realloc "
           "fails\n", name);
    return 0;
}

int main(void)
{
    /* Prove the realloc seam is actually in effect before relying on it: with
     * the failure armed, the interposed realloc must hand back NULL. If some
     * build quirk bypassed the `#define`, catch it here instead of letting the
     * test pass for the wrong reason. */
    g_fail_next_realloc = 1;
    void *probe = realloc(NULL, 8);
    if (probe != NULL) {
        fprintf(stderr, "FAIL: realloc seam not in effect (armed realloc "
                        "returned non-NULL); cannot exercise the OOM path\n");
        free(probe);
        return 1;
    }
    g_fail_next_realloc = 0;

    int rc = 0;
    rc |= run_case("kv_push", GROW_PUSH);
    rc |= run_case("kv_pushp", GROW_PUSHP);
    return rc;
}
