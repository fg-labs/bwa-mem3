// test/framework/shm_registry_guard.h — refuse to destroy a registry a test
// does not own.
//
// bwa_shm_destroy drops every index registered on the host, so a test that
// calls it must first check that the registry holds nothing but its own
// fixture. Header-only so standalone tests can include it without a link
// dependency.

#ifndef BWA_TESTS_SHM_REGISTRY_GUARD_H
#define BWA_TESTS_SHM_REGISTRY_GUARD_H

#include "bwa_shm.h"   /* BWA_SHM_CTL_SIZE */

#include <errno.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cstdint>
#include <cstring>

namespace bwa_tests {

/* True when the shm registry `registry` (normally BWA_SHM_CTL_NAME) is absent
 * or lists no index other than `fixture`. Parses it as ctl_walk does: a uint16
 * entry count, a uint16 next-write offset, then packed (int64 size,
 * NUL-terminated basename) entries. A registry smaller than BWA_SHM_CTL_SIZE is
 * empty, as ctl_walk treats it. Any registry this cannot read counts as holding
 * another index, so the caller never runs bwa_shm_destroy over indexes it does
 * not own. */
inline bool shm_registry_holds_only(const char *registry, const char *fixture)
{
    const int fd = shm_open(registry, O_RDONLY, 0);
    if (fd < 0) return errno == ENOENT;
    struct stat st;
    if (fstat(fd, &st) < 0) { close(fd); return false; }
    if ((size_t)st.st_size < BWA_SHM_CTL_SIZE) { close(fd); return true; }
    void *m = mmap(NULL, BWA_SHM_CTL_SIZE, PROT_READ, MAP_SHARED, fd, 0);
    close(fd);
    if (m == MAP_FAILED) return false;
    const uint8_t *p = (const uint8_t *)m + 4;   /* past the two uint16 fields */
    const uint8_t *end = (const uint8_t *)m + BWA_SHM_CTL_SIZE;
    uint16_t n_entries = 0;
    std::memcpy(&n_entries, m, sizeof(n_entries));
    bool only = true;
    for (uint16_t i = 0; i < n_entries && only; i++) {
        if ((size_t)(end - p) <= sizeof(int64_t)) { only = false; break; }
        const char *name = (const char *)(p + sizeof(int64_t));
        const size_t name_max = (size_t)(end - (const uint8_t *)name);
        const size_t len = strnlen(name, name_max);
        only = len < name_max && std::strcmp(name, fixture) == 0;
        p = (const uint8_t *)name + len + 1;
    }
    munmap(m, BWA_SHM_CTL_SIZE);
    return only;
}

}  // namespace bwa_tests

#endif  // BWA_TESTS_SHM_REGISTRY_GUARD_H
