/* Verifies that bwa_shm_destroy unlinks the named POSIX semaphore that
 * guards the /bwactl registry. Together with the lock acquired in
 * bwa_shm_stage, this is the recovery path for a stager that segfaults or
 * is kill -9'd while holding the lock — without sem_unlink, the stuck
 * semaphore would block every subsequent stage forever (POSIX named
 * semaphores have no SEM_UNDO equivalent).
 *
 * Tests:
 *   1. After bwa_shm_destroy, sem_open with O_CREAT|O_EXCL succeeds — i.e.
 *      the prior name no longer exists.
 *   2. The sequence is repeatable (force the semaphore in, destroy, again).
 */

#include "bwa_shm.h"

#include <semaphore.h>
#include <fcntl.h>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define BWA_SHM_LOCK_NAME "/bwactl_lock"

#define CHECK(cond) do { \
    if (!(cond)) { \
        std::fprintf(stderr, "CHECK failed at %s:%d: %s (errno=%d %s)\n", \
                     __FILE__, __LINE__, #cond, errno, std::strerror(errno)); \
        std::fflush(stderr); \
        std::abort(); \
    } \
} while (0)

static void clear_lock_name(void)
{
    /* Best-effort. ENOENT is fine; we just want a clean baseline. */
    sem_unlink(BWA_SHM_LOCK_NAME);
}

static void force_lock_exists(void)
{
    sem_t *s = sem_open(BWA_SHM_LOCK_NAME, O_CREAT, 0644, 1);
    CHECK(s != SEM_FAILED);
    CHECK(sem_close(s) == 0);
}

static void check_lock_absent(void)
{
    /* O_CREAT|O_EXCL fails with EEXIST iff the name is currently linked. */
    sem_t *t = sem_open(BWA_SHM_LOCK_NAME, O_CREAT | O_EXCL, 0644, 1);
    CHECK(t != SEM_FAILED);
    sem_close(t);
    sem_unlink(BWA_SHM_LOCK_NAME);
}

int main(void)
{
    /* Start from a known-clean state. The registry segment may or may not
     * exist; either way bwa_shm_destroy is idempotent. */
    clear_lock_name();
    CHECK(bwa_shm_destroy() == 0);

    /* 1. Force the semaphore name into existence (simulates a prior stager
     *    leaving the lock object behind), then destroy must remove it. */
    force_lock_exists();
    CHECK(bwa_shm_destroy() == 0);
    check_lock_absent();

    /* 2. Repeatable. */
    force_lock_exists();
    CHECK(bwa_shm_destroy() == 0);
    check_lock_absent();

    /* 3. Destroy on an already-clean state is still 0 (no semaphore to
     *    unlink). */
    CHECK(bwa_shm_destroy() == 0);

    std::printf("shm_lock_destroy_test: OK\n");
    return 0;
}
