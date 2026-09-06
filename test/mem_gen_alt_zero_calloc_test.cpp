// test/mem_gen_alt_zero_calloc_test.cpp -- mem_gen_alt must not treat
// calloc(0) == NULL as OOM. calloc(0) may legally return NULL; when a->n == 0
// the cnt buffer is unused, so the guard is xassert(cnt != NULL || a->n == 0).
//
// Like test/bns_zero_calloc_test.cpp, this interposes calloc so a zero-size
// request returns NULL (the legal-but-less-common branch that glibc/macOS
// libmalloc do not take), which is the only way to make the regression
// observable. The fixed code returns normally (exit 0); the pre-fix code
// aborts on the unconditional xassert(cnt != NULL), which `make test` sees as a
// nonzero exit. mem_gen_alt dereferences none of opt/bns/pac/query when
// a->n == 0 (the region loop and everything past the tot==0 early-out are
// skipped), so zero-initialized structs suffice -- no index fixture needed.

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "bwamem.h"

// Armed only around the call under test so startup allocations take the normal path.
static int g_zero_calloc_returns_null = 0;

extern "C" void *calloc(size_t nmemb, size_t size) {
    if (nmemb == 0 || size == 0) {
        if (g_zero_calloc_returns_null) return NULL;
        return malloc(1);   // unarmed: a unique freeable pointer, as the host allocators do
    }
    if (nmemb > SIZE_MAX / size) { errno = ENOMEM; return NULL; }
    const size_t total = nmemb * size;
    void *p = malloc(total);
    if (p != NULL) memset(p, 0, total);
    return p;
}

int main(void) {
    mem_opt_t opt;      memset(&opt, 0, sizeof opt);
    bntseq_t  bns;      memset(&bns, 0, sizeof bns);
    uint8_t   pac = 0;
    mem_alnreg_v a;     memset(&a, 0, sizeof a);   // n = 0, m = 0, a = NULL
    int *hn = NULL;

    g_zero_calloc_returns_null = 1;
    char **xa = mem_gen_alt(&opt, &bns, &pac, &a, 0, NULL, &hn, NULL);
    g_zero_calloc_returns_null = 0;

    // Reaching here means no abort: the calloc(0)==NULL case was handled. With
    // a->n == 0 there is nothing to emit, so XA and the hit-count array are empty.
    if (xa != NULL) { fprintf(stderr, "FAIL: expected NULL XA for a->n==0\n"); return 1; }
    free(hn);   // NULL here; free(NULL) is a no-op
    fprintf(stderr, "mem_gen_alt_zero_calloc_test OK\n");
    return 0;
}
