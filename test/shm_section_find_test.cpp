#include "bwa_shm.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define CHECK(cond) do { \
    if (!(cond)) { \
        std::fprintf(stderr, "CHECK failed at %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        std::fflush(stderr); \
        std::abort(); \
    } \
} while (0)

int main() {
    uint8_t buf[1024];
    std::memset(buf, 0, sizeof(buf));

    bwa_shm_header_t hdr = {};
    hdr.magic = BWA_SHM_MAGIC;
    hdr.version = BWA_SHM_VERSION;
    hdr.n_sections = 2;
    hdr.total_size = sizeof(buf);
    std::memcpy(buf, &hdr, sizeof(hdr));

    bwa_shm_section_t s1 = { BWA_SHM_SEC_FMI_SCALARS, 0, 128, 56, 0 };
    bwa_shm_section_t s2 = { BWA_SHM_SEC_PAC,         0, 192, 789, 0 };
    std::memcpy(buf + sizeof(hdr),                    &s1, sizeof(s1));
    std::memcpy(buf + sizeof(hdr) + sizeof(s1),       &s2, sizeof(s2));

    uint64_t off = 0, sz = 0;

    CHECK(bwa_shm_section_find(buf, BWA_SHM_SEC_FMI_SCALARS, &off, &sz) == 0);
    CHECK(off == 128 && sz == 56);

    off = 0; sz = 0;
    CHECK(bwa_shm_section_find(buf, BWA_SHM_SEC_PAC, &off, &sz) == 0);
    CHECK(off == 192 && sz == 789);

    CHECK(bwa_shm_section_find(buf, BWA_SHM_SEC_BNS_AMBS, &off, &sz) == -1);
    CHECK(bwa_shm_section_find(NULL, BWA_SHM_SEC_PAC, &off, &sz) == -1);

    /* Bad magic. */
    buf[0] ^= 0xFF;
    CHECK(bwa_shm_section_find(buf, BWA_SHM_SEC_PAC, &off, &sz) == -1);
    buf[0] ^= 0xFF;

    /* Bad version. */
    bwa_shm_header_t bad_ver = hdr;
    bad_ver.version = 99u;
    std::memcpy(buf, &bad_ver, sizeof(bad_ver));
    CHECK(bwa_shm_section_find(buf, BWA_SHM_SEC_PAC, &off, &sz) == -1);
    std::memcpy(buf, &hdr, sizeof(hdr));

    /* Sanity: lookup works again after restore. */
    CHECK(bwa_shm_section_find(buf, BWA_SHM_SEC_PAC, &off, &sz) == 0);

    std::printf("shm_section_find_test: OK\n");
    return 0;
}
