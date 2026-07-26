// test/bns_zero_calloc_test.cpp — a zero-sequence index must load on an
// allocator where calloc(0, ...) returns NULL.
//
// bns_restore_core sizes bns->anns from the sequence count in the .ann header
// and then guards the allocation. That guard has to tell two different NULLs
// apart:
//
//   * calloc(n, ...) returning NULL for n > 0 -- genuinely out of memory.
//   * calloc(0,  ...) returning NULL          -- explicitly permitted by C99
//     7.20.3 for a zero-size request, and not a failure at all.
//
// An unconditional non-NULL guard conflates them, and because the guard is an
// xassert rather than an assert it fires in release builds too -- which is the
// whole point of the change it belongs to.
//
// Why this is a standalone binary rather than a case in the doctest suite:
// glibc and macOS libmalloc both return a unique non-NULL pointer for
// calloc(0, ...), so on those allocators the buggy and fixed code behave
// identically and no test that uses the real allocator can tell them apart.
// The only way to make the regression detectable is to supply an allocator that
// takes the other legal option -- and interposing calloc inside the shared
// doctest binary would redirect every allocation of all ~180 other test cases
// to it. Here the blast radius is one executable that does nothing else.
//
// The interposition works because bntseq.o is linked statically into this
// executable, so its call to calloc binds to the definition below at link time.
// A regression aborts the process on the xassert, so `make test` sees a nonzero
// exit; that is the whole assertion.
//
// The .ann/.amb/.pac fixture is generated at run time -- no test data files.

#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "bntseq.h"

/* Armed only around the call under test, so fixture setup and every startup
 * allocation take the ordinary path. Not thread-safe and does not need to be:
 * this binary is single-threaded by construction. */
static int g_zero_calloc_returns_null = 0;

/* Interposed calloc.
 *
 * When armed, a zero-size request takes the legal-but-less-common branch and
 * returns NULL. Otherwise this is an ordinary calloc: malloc + zero-fill, with
 * calloc's own overflow check preserved, since delegating to the real calloc
 * would mean resolving it through dlsym and dlsym allocates. */
extern "C" void *calloc(size_t nmemb, size_t size)
{
    if (nmemb == 0 || size == 0) {
        if (g_zero_calloc_returns_null) return NULL;
        /* Unarmed zero-size request: hand back a unique freeable pointer, which
         * is the other behaviour the standard allows and what the host
         * allocators do. */
        return malloc(1);
    }
    if (nmemb > SIZE_MAX / size) {   /* what calloc must reject rather than wrap */
        errno = ENOMEM;
        return NULL;
    }
    const size_t total = nmemb * size;
    void *p = malloc(total);
    if (p != NULL) memset(p, 0, total);
    return p;
}

namespace {

/* Write `contents` to `path`, or die. */
void write_file(const char *path, const char *contents)
{
    FILE *fp = fopen(path, "w");
    if (fp == NULL) {
        fprintf(stderr, "FAIL: cannot create %s: %s\n", path, strerror(errno));
        exit(1);
    }
    if (fputs(contents, fp) < 0 || fclose(fp) != 0) {
        fprintf(stderr, "FAIL: cannot write %s: %s\n", path, strerror(errno));
        exit(1);
    }
}

}  // namespace

int main(void)
{
    /* Honour TMPDIR, as test/unit/test_fmi_pread_from_stream.cpp does: CI
     * runners and sandboxes often point it somewhere with room, and /tmp is
     * RAM-backed tmpfs on several modern distributions. */
    const char *tmpdir = getenv("TMPDIR");
    if (tmpdir == NULL || tmpdir[0] == '\0') tmpdir = "/tmp";
    char dir[PATH_MAX];
    const int dir_len = snprintf(dir, sizeof(dir), "%s%sbns_zero_calloc_XXXXXX",
                                 tmpdir, tmpdir[strlen(tmpdir) - 1] == '/' ? "" : "/");
    if (dir_len < 0 || (size_t)dir_len >= sizeof(dir)) {
        fprintf(stderr, "FAIL: TMPDIR path too long: %s\n", tmpdir);
        return 1;
    }
    if (mkdtemp(dir) == NULL) {
        fprintf(stderr, "FAIL: mkdtemp(%s): %s\n", dir, strerror(errno));
        return 1;
    }

    char ann[PATH_MAX], amb[PATH_MAX], pac[PATH_MAX];
    snprintf(ann, sizeof(ann), "%s/ref.ann", dir);
    snprintf(amb, sizeof(amb), "%s/ref.amb", dir);
    snprintf(pac, sizeof(pac), "%s/ref.pac", dir);

    /* .ann header is "l_pac n_seqs seed"; .amb is "l_pac n_seqs n_holes". The
     * two must agree on l_pac and n_seqs or bns_restore_core rejects the pair.
     * With n_seqs and n_holes both 0, neither file has a body. */
    write_file(ann, "0 0 11\n");
    write_file(amb, "0 0 0\n");
    write_file(pac, "");   /* opened, but a zero-length genome reads nothing */

    /* Prove the interposition is live before relying on it: if calloc were NOT
     * bound to the definition above, this returns non-NULL and the real test
     * below would silently pass for the wrong reason. */
    g_zero_calloc_returns_null = 1;
    void *probe = calloc(0, sizeof(int));
    if (probe != NULL) {
        fprintf(stderr, "FAIL: calloc interposition is not in effect "
                        "(calloc(0, 4) returned non-NULL); this test cannot "
                        "prove anything -- check the link line\n");
        free(probe);
        return 1;
    }

    /* The call under test. With the guard wrong this aborts inside xassert and
     * never returns, which is exactly the regression being detected. */
    bntseq_t *bns = bns_restore_core(ann, amb, pac);
    g_zero_calloc_returns_null = 0;

    int rc = 0;
    if (bns == NULL) {
        fprintf(stderr, "FAIL: bns_restore_core returned NULL for a "
                        "zero-sequence index\n");
        rc = 1;
    } else {
        if (bns->n_seqs != 0) {
            fprintf(stderr, "FAIL: n_seqs is %d, expected 0\n", bns->n_seqs);
            rc = 1;
        }
        if (bns->l_pac != 0) {
            fprintf(stderr, "FAIL: l_pac is %lld, expected 0\n",
                    (long long)bns->l_pac);
            rc = 1;
        }
        bns_destroy(bns);   /* must also survive anns == NULL */
    }

    unlink(ann);
    unlink(amb);
    unlink(pac);
    rmdir(dir);

    if (rc == 0)
        printf("PASS: a zero-sequence index loads when calloc(0, ...) is NULL\n");
    return rc;
}
