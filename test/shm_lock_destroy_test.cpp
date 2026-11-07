/* Exercise the real destroy implementation using only per-run owned names.
 * A missing semaphore unlink must fail without consulting the host registry.
 * Also checks the registry guard the legacy global-registry tests use, on a
 * private registry name.
 */

// Include declarations before replacing the fixed names: bwa_shm.h is guarded,
// so its inclusion inside bwa_shm.cpp cannot reset them to the production names.
#define bwa_shm_destroy test_private_shm_destroy
#include "bwa_shm.h"
#include "framework/shm_registry_guard.h"
#include "framework/temp_dir.h"

#include <semaphore.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <functional>
#include <stdexcept>
#include <initializer_list>
#include <string>
#include <thread>
#include <vector>

static char private_control[32], private_lock[32], private_index[32];
// How long destroy waits for a held stager lock; a check shortens it.
static int g_destroy_wait_ms = 5000;
#define BWA_SHM_DESTROY_LOCK_WAIT_MS g_destroy_wait_ms
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
void write_registry(const OwnedShm &registry,
                    std::initializer_list<const char *> names = {"first", "second"}) {
    const int fd = shm_open(registry.name(), O_RDWR, 0);
    CHECK(fd >= 0);
    const bool sized = ftruncate(fd, BWA_SHM_CTL_SIZE) == 0;
    void *map = sized ? mmap(NULL, BWA_SHM_CTL_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0) : MAP_FAILED;
    close(fd);
    CHECK(map != MAP_FAILED);
    memset(map, 0, BWA_SHM_CTL_SIZE);
    const uint16_t count = static_cast<uint16_t>(names.size());
    uint16_t offset = 4;
    memcpy(map, &count, sizeof count);
    for (const char *name : names) {
        const int64_t bytes = 1;
        memcpy(static_cast<char *>(map) + offset, &bytes, sizeof bytes);
        offset += sizeof bytes;
        memcpy(static_cast<char *>(map) + offset, name, strlen(name) + 1);
        offset += strlen(name) + 1;
    }
    memcpy(static_cast<char *>(map) + 2, &offset, sizeof offset);
    CHECK(munmap(map, BWA_SHM_CTL_SIZE) == 0);
}

// A raw registry: `fill` writes every byte of the zeroed BWA_SHM_CTL_SIZE map.
void write_raw_registry(const OwnedShm &registry, const std::function<void(uint8_t *)> &fill) {
    const int fd = shm_open(registry.name(), O_RDWR, 0);
    CHECK(fd >= 0);
    const bool sized = ftruncate(fd, BWA_SHM_CTL_SIZE) == 0;
    void *map = sized ? mmap(NULL, BWA_SHM_CTL_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0) : MAP_FAILED;
    close(fd);
    CHECK(map != MAP_FAILED);
    memset(map, 0, BWA_SHM_CTL_SIZE);
    fill(static_cast<uint8_t *>(map));
    CHECK(munmap(map, BWA_SHM_CTL_SIZE) == 0);
}

// Append `name` to the existing registry, as a stager does under the lock.
// False if the registry is gone.
bool append_entry(const char *name) {
    int fd = -1;
    void *m = ctl_open_rw(&fd, /*create=*/false);
    if (m == NULL) return false;
    ctl_header h = ctl_header_load(m);
    const int64_t bytes = 1;
    memcpy(static_cast<uint8_t *>(m) + h.next_write, &bytes, sizeof bytes);
    memcpy(static_cast<uint8_t *>(m) + h.next_write + sizeof bytes, name, strlen(name) + 1);
    h.n_entries++;
    h.next_write = static_cast<uint16_t>(h.next_write + sizeof bytes + strlen(name) + 1);
    ctl_header_store(m, h);
    ctl_close(m, fd);
    return true;
}

// Point the private names at a fresh per-run token; `names` owns the token.
std::string use_private_names(const bwa_tests::TempDir &names) {
    const std::string token = names.path().substr(names.path().rfind('.') + 1);
    CHECK(token.size() == 6);
    snprintf(private_control, sizeof private_control, "/sd%sc", token.c_str());
    snprintf(private_lock, sizeof private_lock, "/sd%sl", token.c_str());
    snprintf(private_index, sizeof private_index, "/sd%s-", token.c_str());
    return token;
}

void check_destroy() {
    bwa_tests::TempDir names("shm_destroy");
    const std::string token = use_private_names(names);

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
    // Destroy takes the stager lock before it walks the registry: a stager
    // holding it registers "second" before releasing, and destroy still
    // removes it. A lockless walk would run first and leave "second" behind.
    {
        OwnedShm registry(private_control);
        OwnedSemaphore lock(private_lock);
        OwnedShm first(std::string(private_index) + "first");
        OwnedShm second(std::string(private_index) + "second");
        write_registry(registry, {"first"});
        sem_t *s = sem_open(private_lock, 0);
        CHECK(s != SEM_FAILED);
        CHECK(sem_wait(s) == 0);
        std::atomic<bool> added{false};
        std::thread stager([&] {
            std::this_thread::sleep_for(std::chrono::milliseconds(300));
            added = append_entry("second");   // fails if destroy ran first
            sem_post(s);
        });
        const auto t0 = std::chrono::steady_clock::now();
        const int rc = test_private_shm_destroy();
        const auto waited = std::chrono::steady_clock::now() - t0;
        stager.join();
        sem_close(s);
        CHECK(rc == 0);
        CHECK(added);
        CHECK(waited >= std::chrono::milliseconds(250));
        CHECK(!shm_exists(first.name()));
        CHECK(!shm_exists(second.name()));
        CHECK(!shm_exists(registry.name()));
        CHECK(!semaphore_exists(lock.name()));
        check_sentinels();
    }
    // A lock that stays held (a stager killed inside its registry update) does
    // not block destroy: past the wait it proceeds and unlinks the lock.
    {
        OwnedShm registry(private_control);
        OwnedSemaphore lock(private_lock);
        OwnedShm first(std::string(private_index) + "first");
        write_registry(registry, {"first"});
        sem_t *s = sem_open(private_lock, 0);
        CHECK(s != SEM_FAILED);
        CHECK(sem_wait(s) == 0);
        sem_close(s);   // never posted
        g_destroy_wait_ms = 100;
        const int rc = test_private_shm_destroy();
        g_destroy_wait_ms = 5000;
        CHECK(rc == 0);
        CHECK(!shm_exists(first.name()));
        CHECK(!shm_exists(registry.name()));
        CHECK(!semaphore_exists(lock.name()));
        check_sentinels();
    }
    check_sentinels();
}

// The refusal predicate the legacy global-registry tests use, run on a private
// registry name so the host registry is never read or written.
void check_registry_guard() {
    bwa_tests::TempDir names("shm_guard");
    use_private_names(names);
    CHECK(bwa_tests::shm_registry_holds_only(private_control, "phix.fa"));   // absent
    {
        OwnedShm registry(private_control);                                  // 0 bytes
        CHECK(bwa_tests::shm_registry_holds_only(private_control, "phix.fa"));
    }
    // macOS sizes a shm object only once, so each registry is created afresh.
    auto holds_only = [](std::initializer_list<const char *> entries, const char *fixture) {
        OwnedShm registry(private_control);
        write_registry(registry, entries);
        return bwa_tests::shm_registry_holds_only(private_control, fixture);
    };
    CHECK(holds_only({}, "phix.fa"));
    CHECK(holds_only({"phix.fa"}, "phix.fa"));
    CHECK(!holds_only({"phix.fa"}, "other.fa"));
    CHECK(!holds_only({"phix.fa", "hg38.fa"}, "phix.fa"));
    CHECK(!holds_only({"hg38.fa"}, "phix.fa"));
    // Malformed registries are refused. One name fills the registry so that
    // only 4 bytes follow it: as the sole entry it is accepted, but a second
    // entry would start past the end and must be refused, not read.
    const size_t long_len = BWA_SHM_CTL_SIZE - BWA_SHM_CTL_HEADER_BYTES - sizeof(int64_t) - 1 - 4;
    const std::string long_name(long_len, 'a');
    auto bounded = [&](uint16_t count) {
        OwnedShm registry(private_control);
        write_raw_registry(registry, [&](uint8_t *m) {
            memcpy(m, &count, sizeof count);
            const uint16_t end = static_cast<uint16_t>(BWA_SHM_CTL_SIZE - 4);
            memcpy(m + 2, &end, sizeof end);
            memcpy(m + BWA_SHM_CTL_HEADER_BYTES + sizeof(int64_t), long_name.c_str(), long_len + 1);
        });
        return bwa_tests::shm_registry_holds_only(private_control, long_name.c_str());
    };
    CHECK(bounded(1));
    CHECK(!bounded(2));
    // A name with no terminating NUL before the end of the registry is refused.
    {
        OwnedShm registry(private_control);
        write_raw_registry(registry, [](uint8_t *m) {
            const uint16_t count = 1;
            memcpy(m, &count, sizeof count);
            uint8_t *name = m + BWA_SHM_CTL_HEADER_BYTES + sizeof(int64_t);
            memset(name, 'a', (size_t)(m + BWA_SHM_CTL_SIZE - name));
        });
        CHECK(!bwa_tests::shm_registry_holds_only(private_control, "phix.fa"));
    }
    CHECK(!shm_exists(private_control));
}

} // namespace

// With no argument every check runs. One argument (destroy or registry-guard)
// runs only that check, so a failure in one cannot mask another.
int main(int argc, char *argv[]) {
    const std::string only = argc > 1 ? argv[1] : "";
    if (argc > 2 || (!only.empty() && only != "destroy" && only != "registry-guard")) {
        std::fprintf(stderr, "usage: %s [destroy|registry-guard]\n", argv[0]);
        return 2;
    }
    try {
        if (only.empty() || only == "destroy") check_destroy();
        if (only.empty() || only == "registry-guard") check_registry_guard();
        std::printf("shm_lock_destroy_test: %s OK\n", only.empty() ? "all checks" : only.c_str());
        return 0;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "%s (errno=%d %s)\n", e.what(), errno, std::strerror(errno));
        return 1;
    }
}
