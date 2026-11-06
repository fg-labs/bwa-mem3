/* Exercise the real destroy implementation using only per-run owned names.
 * A missing semaphore unlink must fail without consulting the host registry.
 */

// Include declarations before replacing the fixed names: bwa_shm.h is guarded,
// so its inclusion inside bwa_shm.cpp cannot reset them to the production names.
#define bwa_shm_destroy test_private_shm_destroy
#include "bwa_shm.h"
#include "framework/temp_dir.h"

#include <semaphore.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

static char private_control[32], private_lock[32], private_index[32];
#undef BWA_SHM_CTL_NAME
#undef BWA_SHM_LOCK_NAME
#undef BWA_SHM_IDX_PREFIX
#define BWA_SHM_CTL_NAME private_control
#define BWA_SHM_LOCK_NAME private_lock
#define BWA_SHM_IDX_PREFIX private_index

// Defines a distinct destroy symbol, so an accidental archive fallback cannot
// satisfy the call. The other implementation symbols keep libbwa's global-name
// object out of this executable, matching the private seed-table fixture.
#include "../src/bwa_shm.cpp"
#undef BWA_SHM_CTL_NAME
#undef BWA_SHM_LOCK_NAME
#undef BWA_SHM_IDX_PREFIX
#undef bwa_shm_destroy

#define CHECK(cond) do { \
    if (!(cond)) throw std::runtime_error(std::string("CHECK failed: ") + #cond); \
} while (0)

namespace {

// Unlink only names for which this object successfully obtained O_EXCL ownership.
// This also removes our private leftovers when a destroy assertion throws.
class OwnedShm {
public:
    explicit OwnedShm(const std::string &name) : name_(name) {
        const int fd = shm_open(name_.c_str(), O_CREAT | O_EXCL | O_RDWR, 0600);
        CHECK(fd >= 0);
        close(fd);
    }
    ~OwnedShm() { shm_unlink(name_.c_str()); }
    OwnedShm(const OwnedShm &) = delete;
    OwnedShm &operator=(const OwnedShm &) = delete;
    const char *name() const { return name_.c_str(); }
private:
    std::string name_;
};

class OwnedSemaphore {
public:
    explicit OwnedSemaphore(const std::string &name) : name_(name) {
        sem_t *s = sem_open(name_.c_str(), O_CREAT | O_EXCL, 0600, 1);
        CHECK(s != SEM_FAILED);
        CHECK(sem_close(s) == 0);
    }
    ~OwnedSemaphore() { sem_unlink(name_.c_str()); }
    OwnedSemaphore(const OwnedSemaphore &) = delete;
    OwnedSemaphore &operator=(const OwnedSemaphore &) = delete;
    const char *name() const { return name_.c_str(); }
private:
    std::string name_;
};

bool shm_exists(const char *name) {
    const int fd = shm_open(name, O_RDONLY, 0);
    if (fd < 0) { CHECK(errno == ENOENT); return false; }
    CHECK(close(fd) == 0);
    return true;
}

bool semaphore_exists(const char *name) {
    sem_t *s = sem_open(name, 0);
    if (s == SEM_FAILED) { CHECK(errno == ENOENT); return false; }
    CHECK(sem_close(s) == 0);
    return true;
}

// Registry format is a count/write-offset header followed by size/name entries.
void write_registry(const OwnedShm &registry) {
    const int fd = shm_open(registry.name(), O_RDWR, 0);
    CHECK(fd >= 0);
    const bool sized = ftruncate(fd, BWA_SHM_CTL_SIZE) == 0;
    void *map = sized ? mmap(NULL, BWA_SHM_CTL_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0) : MAP_FAILED;
    close(fd);
    CHECK(map != MAP_FAILED);
    memset(map, 0, BWA_SHM_CTL_SIZE);
    const uint16_t count = 2;
    uint16_t offset = 4;
    memcpy(map, &count, sizeof count);
    for (const char *name : {"first", "second"}) {
        const int64_t bytes = 1;
        memcpy(static_cast<char *>(map) + offset, &bytes, sizeof bytes);
        offset += sizeof bytes;
        memcpy(static_cast<char *>(map) + offset, name, strlen(name) + 1);
        offset += strlen(name) + 1;
    }
    memcpy(static_cast<char *>(map) + 2, &offset, sizeof offset);
    CHECK(munmap(map, BWA_SHM_CTL_SIZE) == 0);
}

void check_destroy() {
    bwa_tests::TempDir names("shm_destroy");
    const std::string token = names.path().substr(names.path().rfind('.') + 1);
    CHECK(token.size() == 6);
    snprintf(private_control, sizeof private_control, "/sd%sc", token.c_str());
    snprintf(private_lock, sizeof private_lock, "/sd%sl", token.c_str());
    snprintf(private_index, sizeof private_index, "/sd%s-", token.c_str());

    // Different private registry/lock, plus an unregistered same-prefix segment.
    // They must remain intact across every destroy, including repeated no-ops.
    OwnedShm sentinel_registry("/sd" + token + "u");
    OwnedSemaphore sentinel_lock("/sd" + token + "v");
    OwnedShm sentinel_index(std::string(private_index) + "unrelated");
    auto check_sentinels = [&] {
        CHECK(shm_exists(sentinel_registry.name()));
        CHECK(semaphore_exists(sentinel_lock.name()));
        CHECK(shm_exists(sentinel_index.name()));
    };

    // Reserve each name exclusively before passing it to the isolated destroy.
    // The first pass covers a missing registry with an existing stager lock.
    {
        OwnedShm registry(private_control);
        OwnedSemaphore lock(private_lock);
        CHECK(shm_unlink(registry.name()) == 0);
        CHECK(test_private_shm_destroy() == 0);
        CHECK(!semaphore_exists(private_lock));
        CHECK(!shm_exists(private_control));
        check_sentinels();
        CHECK(test_private_shm_destroy() == 0);
    }
    for (int repeat = 0; repeat < 2; ++repeat) {
        OwnedShm registry(private_control);
        OwnedSemaphore lock(private_lock);
        OwnedShm first(std::string(private_index) + "first");
        OwnedShm second(std::string(private_index) + "second");
        write_registry(registry);
        CHECK(bwa_shm_test("first") == 1);
        CHECK(bwa_shm_test("second") == 1);
        CHECK(test_private_shm_destroy() == 0);
        CHECK(!shm_exists(first.name()));
        CHECK(!shm_exists(second.name()));
        CHECK(!shm_exists(registry.name()));
        CHECK(!semaphore_exists(lock.name()));
        CHECK(bwa_shm_test("first") == 0);
        CHECK(bwa_shm_test("second") == 0);
        // Successful exclusive creation proves the semaphore name was unlinked.
        OwnedSemaphore replacement(private_lock);
        CHECK(test_private_shm_destroy() == 0);
        CHECK(!semaphore_exists(private_lock));
        CHECK(test_private_shm_destroy() == 0);
        check_sentinels();
    }
    check_sentinels();
}

} // namespace

int main() {
    try {
        check_destroy();
        std::printf("shm_lock_destroy_test: private destroy, repeatability and unrelated preservation OK\n");
        return 0;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "%s (errno=%d %s)\n", e.what(), errno, std::strerror(errno));
        return 1;
    }
}
