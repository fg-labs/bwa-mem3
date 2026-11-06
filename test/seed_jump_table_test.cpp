/* Exhaustive gate for the seed jump table (FMI_search::seed_jump).
 *
 * The lockstep parity tests only read the entries of strings that occur in
 * their reads. This test checks every entry of every level against the
 * interval the seeding walks reach by stepping: the 1-base interval
 * (count[a], count[3 - a], count[a + 1] - count[a]), then one forward step per
 * base (swap k <-> l, backwardExt with 3 - c, swap back), s == 0 strings and
 * strings whose interval holds the sentinel row included. It does so at the
 * depths on both sides of the builder's serial prefix (2, 3, 4) and at the
 * default depth, each at several builder thread counts, so a wrong level
 * offset, child index or per-prefix partition fails here whatever the reads.
 * At depths 2-4 it also checks every entry's occurrence count s, and the
 * interval starts k and l, against a brute-force count and a naive suffix
 * array of the index text (the .pac forward strand plus its reverse
 * complement), so a defect that backwardExt or the 1-base interval shares
 * with the stepped walk above fails too.
 * It checks that the lockstep SMEM and bwtSeed drivers take their first
 * steps from the table at the default depth: with the table's entries
 * replaced by other valid intervals, their output changes.
 * It also checks the depth BWA3_BWTSEED_JUMP selects (unset, empty, 0,
 * in-range and malformed values), the build_seed_jump=false opt-out on
 * load_index and on fmi_seed_open, the table an index attached from a packed
 * shm segment builds, and load_index's shm-attach path with and without the
 * opt-out (when a phiX segment can be staged on this host).
 *
 * Usage: seed_jump_table_test <bwa-mem3 index prefix>
 */

#define private public
#include "FMI_search.h"
#undef private
#include "bwa.h"
#include "bwa_shm.h"
#include "bntseq.h"
#include "fmi_seed_api.h"
#include "seed_jump_test_util.h"

#include <algorithm>
#include <exception>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <semaphore.h>
#include <limits.h>
#include <poll.h>
#include <signal.h>
#include <sys/wait.h>
#include "framework/temp_dir.h"
#include "framework/stderr_capture.h"

// Own every IPC name before the real implementation can open it. Reserved
// index descriptors are handed to the stager's O_EXCL open without releasing
// their names, closing the reservation-to-stage race. Packing/mapping stay real.
// Exactly one process owns a name, from its exclusive creation to its unlink;
// ownership never crosses fork. A forked child inherits only a read-only view.
#define NON_OWNER_REMOVE_EXIT 93
class PrivateObjects;
// The owner of the main namespace, if any. Its exit handler cleans up a plain
// exit() in that process only; inherited copies in children do nothing.
static PrivateObjects *exit_owner = NULL;
class PrivateObjects {
    struct Object { std::string name; bool semaphore; int fd; };
    std::vector<Object> objects;
    const pid_t owner = getpid();
public:
    PrivateObjects() = default;
    PrivateObjects(const PrivateObjects &) = delete;
    PrivateObjects &operator=(const PrivateObjects &) = delete;
    // Unregister on destruction so the exit handler never sees a dead owner.
    ~PrivateObjects() { clear(); if (exit_owner == this) exit_owner = NULL; }
    pid_t owner_pid() const { return owner; }
    bool owns(const std::string &name) const {
        for (const auto &o : objects) if (o.name == name) return true;
        return false;
    }
    int descriptor(const std::string &name) const {
        for (const auto &o : objects) if (o.name == name && !o.semaphore) return o.fd;
        return -1;
    }
    bool acquire(const std::string &name, bool semaphore = false) {
        if (owns(name)) return false;
        int fd = -1;
        if (semaphore) {
            sem_t *lock = sem_open(name.c_str(), O_CREAT | O_EXCL, 0600, 1);
            if (lock == SEM_FAILED) return false;
            sem_close(lock);
        } else {
            fd = shm_open(name.c_str(), O_CREAT | O_EXCL | O_RDWR, 0600);
            if (fd < 0) return false;
        }
        objects.push_back({name, semaphore, fd});
        return true;
    }
    void remove(const std::string &name) {
        if (getpid() != owner) {
            // A non-owner can never unlink, so this is a test defect, not a no-op.
            fprintf(stderr, "[FATAL] PrivateObjects::remove(%s) from pid %d; owner is pid %d\n",
                    name.c_str(), (int)getpid(), (int)owner);
            fflush(stderr);
            _exit(NON_OWNER_REMOVE_EXIT);
        }
        for (auto it = objects.begin(); it != objects.end(); ++it) {
            if (it->name != name) continue;
            if (it->semaphore) sem_unlink(name.c_str());
            else { shm_unlink(name.c_str()); close(it->fd); }
            objects.erase(it);
            return;
        }
    }
    void clear() {
        if (getpid() != owner) return;
        while (!objects.empty()) remove(objects.back().name);
    }
    bool reserve(const std::string &control, const std::string &lock,
                 const std::vector<std::string> &segments) {
        if (!acquire(control) || !acquire(lock, true)) { clear(); return false; }
        for (const auto &name : segments)
            if (!acquire(name)) { clear(); return false; }
        return true;
    }
};

static void owner_exit_cleanup()
{
    if (exit_owner != NULL && getpid() == exit_owner->owner_pid()) exit_owner->clear();
}

static char private_control[32], private_lock[32], private_index[32];
static PrivateObjects *stage_objects = NULL;
static int reserved_stage_open(const char *name, int flags, mode_t mode) {
    if (!stage_objects || !stage_objects->owns(name)) { errno = EPERM; return -1; }
    if ((flags & O_EXCL) && stage_objects && stage_objects->owns(name)) {
        const int fd = stage_objects->descriptor(name);
        if (fd >= 0) return dup(fd);
    }
    return shm_open(name, flags, mode);
}
// When set, each suppressed unlink writes its name and a newline here, so a test
// can show the guard was exercised rather than the call never being made.
static int suppressed_unlink_fd = -1;
// Production failure paths unlink the reserved index name from the stage child.
// Only the owner may unlink an owned name, so a non-owner's unlink of one is a
// no-op (the owner unlinks it later); any other name passes through.
static bool stage_owned_elsewhere(const char *name)
{
    const bool suppressed = stage_objects && getpid() != stage_objects->owner_pid()
                            && stage_objects->owns(name);
    if (suppressed && suppressed_unlink_fd >= 0) {
        const std::string line = std::string(name) + "\n";
        if (write(suppressed_unlink_fd, line.data(), line.size()) != (ssize_t)line.size()) _exit(92);
    }
    return suppressed;
}
static int reserved_stage_unlink(const char *name)
{
    return stage_owned_elsewhere(name) ? 0 : shm_unlink(name);
}
static int reserved_sem_unlink(const char *name)
{
    return stage_owned_elsewhere(name) ? 0 : sem_unlink(name);
}
#undef BWA_SHM_CTL_NAME
#undef BWA_SHM_LOCK_NAME
#undef BWA_SHM_IDX_PREFIX
#define BWA_SHM_CTL_NAME private_control
#define BWA_SHM_LOCK_NAME private_lock
#define BWA_SHM_IDX_PREFIX private_index
#define shm_open reserved_stage_open
#define shm_unlink reserved_stage_unlink
#define sem_unlink reserved_sem_unlink
#include "../src/bwa_shm.cpp"
#undef sem_unlink
#undef shm_unlink
#undef shm_open
#undef BWA_SHM_CTL_NAME
#undef BWA_SHM_LOCK_NAME
#undef BWA_SHM_IDX_PREFIX

static int32_t total_cases = 0;
static int32_t passed_cases = 0;

static void report(bool ok, const char *name)
{
    total_cases++;
    if (ok) passed_cases++;
    fprintf(stderr, "[%s] %s\n", ok ? "PASS" : "FAIL", name);
}

/* Set BWA3_BWTSEED_JUMP (NULL = unset). */
static void set_jump_env(const char *v)
{
    if (v == NULL) unsetenv("BWA3_BWTSEED_JUMP");
    else setenv("BWA3_BWTSEED_JUMP", v, 1);
}

/* Registry names are injectable so cleanup checks never use the host registry. */
struct ShmRegistry {
    const char *control;
    const char *lock;
    const char *index_prefix;
    PrivateObjects *owner;
};
static ShmRegistry fixture_registry = {
    private_control, private_lock, private_index, NULL
};
typedef void *(*RegistryMap)(void *, size_t, int, int, int, off_t);
typedef int (*RegistryTryWait)(sem_t *);

enum RegistryState { REGISTRY_MISSING, REGISTRY_READABLE, REGISTRY_UNREADABLE };

/* Parse ctl_walk's uint16 count/offset header and packed int64 size + basename
 * entries. Only a readable registry permits cleanup; a short registry is empty,
 * matching ctl_walk, but an incomplete entry makes the registry unreadable. */
static RegistryState read_shm_registry(std::vector<std::string> *names,
                                       const ShmRegistry &registry = fixture_registry,
                                       RegistryMap map_registry = mmap)
{
    names->clear();
    const int fd = shm_open(registry.control, O_RDONLY, 0);
    if (fd < 0) return errno == ENOENT ? REGISTRY_MISSING : REGISTRY_UNREADABLE;
    struct stat st;
    if (fstat(fd, &st) < 0) {
        close(fd);
        return REGISTRY_UNREADABLE;
    }
    if ((size_t)st.st_size < BWA_SHM_CTL_SIZE) {
        close(fd);
        return REGISTRY_READABLE;
    }
    void *m = map_registry(NULL, BWA_SHM_CTL_SIZE, PROT_READ, MAP_SHARED, fd, 0);
    close(fd);
    if (m == MAP_FAILED) return REGISTRY_UNREADABLE;
    const uint8_t *base = (const uint8_t *)m;
    const uint8_t *end = base + BWA_SHM_CTL_SIZE;
    uint16_t n_entries = 0;
    memcpy(&n_entries, base, sizeof(uint16_t));
    bool valid = true;
    const uint8_t *p = base + 4;   /* past n_entries and next_write */
    for (uint16_t i = 0; i < n_entries; i++) {
        if ((size_t)(end - p) <= sizeof(int64_t)) { valid = false; break; }
        const char *name = (const char *)(p + sizeof(int64_t));
        const size_t name_max = (size_t)(end - (const uint8_t *)name);
        const size_t len = strnlen(name, name_max);
        if (len == name_max) { valid = false; break; }
        names->push_back(std::string(name, len));
        p = (const uint8_t *)name + len + 1;
    }
    munmap(m, BWA_SHM_CTL_SIZE);
    return valid ? REGISTRY_READABLE : REGISTRY_UNREADABLE;
}

/* Undo this test's staging of `prefix` without touching any other index.
 * Unlike bwa_shm_destroy, which unlinks every staged index, the registry and
 * its lock, this holds the lock the stager takes (private_lock), so no
 * stage can run while it reads the registry, and drops anything only when the
 * registry holds only the fixture basenames owned by this run, or no entries
 * (an aborted stage's leftover). Then it unlinks the fixture segments and registry,
 * retaining the lock so waiting stagers cannot acquire a different semaphore.
 * The optional second fixture is owned by this test too. If another index is
 * registered, or the registry cannot be read, it drops nothing and warns.
 * With no registry (a stage that failed before creating one) it does nothing. */
static void drop_fixture_shm(const char *prefix,
                             const ShmRegistry &registry = fixture_registry,
                             RegistryMap map_registry = mmap,
                             const char *second_fixture = NULL,
                             RegistryTryWait try_wait = sem_trywait)
{
    const char *slash = strrchr(prefix, '/');
    const char *fixture = slash != NULL ? slash + 1 : prefix;
    const std::string first_segment = std::string(registry.index_prefix) + fixture;
    const std::string second_segment = second_fixture ? std::string(registry.index_prefix) + second_fixture : "";
    if (!registry.owner || !registry.owner->owns(registry.control)
        || !registry.owner->owns(registry.lock) || !registry.owner->owns(first_segment)
        || (second_fixture && !registry.owner->owns(second_segment))) {
        fprintf(stderr, "[WARN] private cleanup lacks ownership: registry %s lock %s segment %s %s; dropped nothing\n",
                registry.control, registry.lock, first_segment.c_str(), second_segment.c_str());
        return;
    }
    std::vector<std::string> names;

    sem_t *lock = sem_open(registry.lock, 0);
    if (lock == SEM_FAILED) {
        /* Every stager creates the lock before the registry. */
        if (read_shm_registry(&names, registry, map_registry) != REGISTRY_MISSING)
            fprintf(stderr, "[WARN] shm registry %s exists but its lock %s does not; "
                    "left private objects for their owning fixture cleanup\n",
                    registry.control, registry.lock);
        return;
    }
    /* Bounded wait: a crashed private stager may leave its semaphore held.
     * The parent owner can unlink its own objects after the child is reaped. */
    bool locked = false;
    int err = 0;
    for (int tries = 0; tries < 1000 && !locked; tries++) {
        locked = try_wait(lock) == 0;
        err = locked ? 0 : errno;
        if (!locked && err != EAGAIN && err != EINTR) break;
        if (!locked) usleep(10000);
    }
    if (!locked) {
        fprintf(stderr, "[WARN] could not take the shm lock %s (%s); left %s staged "
                "(private registry %s; owning fixture retains cleanup responsibility)\n", registry.lock,
                err == EAGAIN || err == EINTR ? "still held after 10 s" : strerror(err), fixture, registry.control);
        sem_close(lock);
        return;
    }
    const RegistryState state = read_shm_registry(&names, registry, map_registry);
    if (state != REGISTRY_READABLE) {
        if (state == REGISTRY_UNREADABLE)
            fprintf(stderr, "[WARN] could not read shm registry %s; dropped nothing\n",
                    registry.control);
        sem_post(lock);
        sem_close(lock);
        return;
    }
    const bool only_fixtures = std::all_of(names.begin(), names.end(), [&](const std::string &n) {
        return n == fixture || (second_fixture != NULL && n == second_fixture);
    });
    if (only_fixtures) {
        const std::string idx = std::string(registry.index_prefix) + fixture;
        registry.owner->remove(idx);   /* absent after an aborted stage */
        if (second_fixture != NULL)
            registry.owner->remove(std::string(registry.index_prefix) + second_fixture);
        registry.owner->remove(registry.control);
        sem_post(lock);
        sem_close(lock);
        return;
    }
    sem_post(lock);
    sem_close(lock);
    std::string staged;
    for (const std::string &n : names) staged += (staged.empty() ? "" : ", ") + n;
    fprintf(stderr, "[WARN] another index was staged in shm while this test ran; dropped "
            "nothing (private segment %s%s, registry %s, lock %s, entries: %s); "
            "only the owning fixture may clean these names\n", registry.index_prefix, fixture, registry.control, registry.lock, staged.c_str());
}

/* These private names exercise cleanup even when staging a packed index is not
 * supported, without changing the host's real staging registry. */
static bool shm_name_exists(const char *name)
{
    int fd = shm_open(name, O_RDONLY, 0);
    if (fd < 0) return false;
    close(fd);
    return true;
}

static void *fail_registry_map(void *, size_t, int, int, int, off_t)
{
    errno = ENOMEM;
    return MAP_FAILED;
}

// The child announces a real EAGAIN at the lock while the writer holds it ('r'),
// then its acquisition after that wait ('l'), then each registry map ('m'). The
// notification wait and the cleaner's lock acquisition are both bounded.
static int lock_reached_fd = -1;
static bool lock_waited = false;
static void announce(char event)
{
    if (lock_reached_fd >= 0 && write(lock_reached_fd, &event, 1) != 1) _exit(91);
}
static bool pipe_ready(int fd)
{
    struct pollfd p = {fd, POLLIN, 0};
    int rc;
    do { rc = poll(&p, 1, 5000); } while (rc < 0 && errno == EINTR);
    return rc == 1 && (p.revents & POLLIN) != 0;
}

static int synchronized_trywait(sem_t *lock)
{
    const int rc = sem_trywait(lock);
    const int err = errno;
    if (rc < 0 && err == EAGAIN && !lock_waited) {
        lock_waited = true;
        announce('r');
    } else if (rc == 0 && lock_waited) {
        announce('l');
    }
    errno = err;
    return rc;
}

static void *announced_registry_map(void *addr, size_t len, int prot, int flags, int fd, off_t off)
{
    announce('m');
    return mmap(addr, len, prot, flags, fd, off);
}

static void check_read_after_lock(void *map, sem_t *lock, const ShmRegistry &registry,
                                  const char *fixture, const char *other)
{
    // Start fixture-only. While the cleaner demonstrably waits on the lock,
    // publish the unrelated entry. A pre-lock snapshot would delete both names.
    const uint16_t one = 1, two = 2;
    memcpy(map, &one, sizeof one);
    int reached[2];
    if (pipe(reached) != 0) { report(false, "lock ordering: notification pipe"); return; }
    if (sem_trywait(lock) != 0) {
        close(reached[0]); close(reached[1]);
        report(false, "lock ordering: writer takes private lock"); return;
    }
    fflush(NULL);
    const pid_t child = fork();
    if (child == 0) {
        close(reached[0]);
        lock_reached_fd = reached[1];
        drop_fixture_shm(fixture, registry, announced_registry_map, NULL, synchronized_trywait);
        _exit(0);
    }
    close(reached[1]);
    const bool waiting = child > 0 && pipe_ready(reached[0]);
    char ready;
    const bool notified = waiting && read(reached[0], &ready, 1) == 1;
    memcpy(map, &two, sizeof two);
    sem_post(lock);
    if (!notified && child > 0) kill(child, SIGKILL);
    int status = 0;
    const bool ended = child > 0 && waitpid(child, &status, 0) == child;
    // After the child exits, its write end is closed: read every later event.
    std::string events;
    char event;
    while (ended && read(reached[0], &event, 1) == 1) events += event;
    close(reached[0]);
    std::vector<std::string> names;
    report(notified && ended && events == "lm",
           "fixture cleanup acquires the stager lock after waiting, then maps the registry once");
    report(notified && ended && WIFEXITED(status) && WEXITSTATUS(status) == 0
               && read_shm_registry(&names, registry) == REGISTRY_READABLE
               && names == std::vector<std::string>{fixture, other}
               && shm_name_exists((std::string(registry.index_prefix) + fixture).c_str())
               && shm_name_exists((std::string(registry.index_prefix) + other).c_str()),
           "fixture cleanup reads the private registry after taking the stager lock");
}

static void check_fixture_cleanup()
{
    PrivateObjects owned;
    std::string ctl, lock_name, idx;
    bool reserved = false;
    for (int attempt = 0; attempt < 20 && !reserved; ++attempt) {
        bwa_tests::TempDir private_names("sj");
        const std::string token = "/sh" + private_names.path().substr(private_names.path().rfind('.') + 1);
        ctl = token + "c"; lock_name = token + "l"; idx = token + "-";
        reserved = owned.reserve(ctl, lock_name, {idx + "fixture", idx + "other"});
    }
    report(reserved, "helper reserves complete private namespace without foreign reuse");
    if (!reserved) return;
    const ShmRegistry registry = {ctl.c_str(), lock_name.c_str(), idx.c_str(), &owned};
    const char *fixture = "fixture", *other = "other";
    const std::string fixture_idx = idx + fixture, other_idx = idx + other;
    std::vector<std::string> names;
    report(read_shm_registry(&names, registry) == REGISTRY_READABLE && names.empty(),
           "fixture registry starts with an exclusively owned empty registry");
    sem_t *lock = sem_open(registry.lock, 0);
    int fd = dup(owned.descriptor(ctl));
    if (fd >= 0)
        report(read_shm_registry(&names, registry) == REGISTRY_READABLE && names.empty(),
               "fixture registry treats a short registry as empty");
    int fixture_fd = dup(owned.descriptor(fixture_idx));
    int other_fd = dup(owned.descriptor(other_idx));
    bool ready = lock != SEM_FAILED && fd >= 0 && fixture_fd >= 0 && other_fd >= 0;
    void *map = MAP_FAILED;
    if (ready && ftruncate(fd, BWA_SHM_CTL_SIZE) == 0)
        map = mmap(NULL, BWA_SHM_CTL_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    ready = ready && map != MAP_FAILED;
    if (ready) {
        memset(map, 0, BWA_SHM_CTL_SIZE);
        uint16_t count = 2, offset = 4;
        memcpy(map, &count, sizeof(count));
        for (const char *name : {fixture, other}) {
            int64_t bytes = 1;
            memcpy((char *)map + offset, &bytes, sizeof(bytes));
            offset += sizeof(bytes);
            memcpy((char *)map + offset, name, strlen(name) + 1);
            offset += strlen(name) + 1;
        }
        memcpy((char *)map + 2, &offset, sizeof(offset));
        check_read_after_lock(map, lock, registry, fixture, other);
        std::vector<uint8_t> before((uint8_t *)map, (uint8_t *)map + BWA_SHM_CTL_SIZE);
        drop_fixture_shm(fixture, registry);
        report(read_shm_registry(&names, registry) == REGISTRY_READABLE && names.size() == 2
                   && names[0] == fixture && names[1] == other
                   && shm_name_exists(fixture_idx.c_str()) && shm_name_exists(other_idx.c_str()),
               "fixture cleanup preserves another registered index and both segments");
        drop_fixture_shm(fixture, registry, fail_registry_map);
        report(read_shm_registry(&names, registry) == REGISTRY_READABLE && names.size() == 2
                   && names[0] == fixture && names[1] == other
                   && memcmp(map, before.data(), before.size()) == 0
                   && shm_name_exists(fixture_idx.c_str()) && shm_name_exists(other_idx.c_str()),
               "fixture cleanup preserves an unreadable registry and its segments");
        /* A valid fixture entry followed by an unterminated entry must not be
         * mistaken for a fixture-only registry and removed. */
        const size_t malformed_offset = 4 + sizeof(int64_t) + strlen(fixture) + 1;
        memset((char *)map + malformed_offset, 'x', BWA_SHM_CTL_SIZE - malformed_offset);
        count = 2;
        memcpy(map, &count, sizeof(count));
        drop_fixture_shm(fixture, registry);
        report(read_shm_registry(&names, registry) == REGISTRY_UNREADABLE
                   && shm_name_exists(registry.control) && shm_name_exists(fixture_idx.c_str())
                   && shm_name_exists(other_idx.c_str()),
               "fixture cleanup preserves a malformed registry and its segments");
        // Restore the valid two-entry registry, then remove both owned names.
        memcpy(map, before.data(), before.size());
        const std::string unrelated_idx = idx + "unrelated";
        PrivateObjects unrelated;
        report(unrelated.acquire(unrelated_idx), "two-owned cleanup: unrelated segment created");
        drop_fixture_shm(fixture, registry, mmap, other);
        report(read_shm_registry(&names, registry) == REGISTRY_MISSING && names.empty()
                   && !shm_name_exists(fixture_idx.c_str()) && !shm_name_exists(other_idx.c_str())
                   && shm_name_exists(unrelated_idx.c_str()),
               "two-owned cleanup removes both segments and registry entries, preserves unrelated segment");
        unrelated.clear();
        // Recreate a one-owned registry to check lock identity independently.
        const int replacement = owned.acquire(ctl) ? dup(owned.descriptor(ctl)) : -1;
        bool replaced = replacement >= 0 && ftruncate(replacement, BWA_SHM_CTL_SIZE) == 0;
        void *replacement_map = replaced ? mmap(NULL, BWA_SHM_CTL_SIZE, PROT_READ | PROT_WRITE,
                                                MAP_SHARED, replacement, 0) : MAP_FAILED;
        replaced = replaced && replacement_map != MAP_FAILED;
        report(replaced, "one-owned cleanup: replacement registry created");
        if (replaced) {
            memcpy(replacement_map, before.data(), before.size());
            const uint16_t one = 1;
            memcpy(replacement_map, &one, sizeof one);
            munmap(replacement_map, BWA_SHM_CTL_SIZE);
        }
        if (replacement >= 0) close(replacement);
        const int recreated_fixture = owned.acquire(fixture_idx) ? dup(owned.descriptor(fixture_idx)) : -1;
        const bool recreated = recreated_fixture >= 0 && shm_name_exists(fixture_idx.c_str());
        report(recreated, "one-owned cleanup requires an exclusively recreated existing segment");
        if (recreated_fixture >= 0) close(recreated_fixture);
        /* Remove only the recreated fixture; the other name is already absent. */
        if (recreated && replaced) drop_fixture_shm(fixture, registry);
        sem_t *reopened = sem_open(registry.lock, 0);
        bool same_lock = false;
        if (reopened != SEM_FAILED) {
            if (sem_trywait(lock) == 0) {
                const int rc = sem_trywait(reopened);
                same_lock = rc < 0 && errno == EAGAIN;
                if (rc == 0) sem_post(reopened);
                sem_post(lock);
            }
            sem_close(reopened);
        }
        report(same_lock, "fixture cleanup keeps the same stager lock linked and released");
        report(!shm_name_exists(registry.control) && !shm_name_exists(fixture_idx.c_str())
                   && !shm_name_exists(other_idx.c_str()),
               "fixture-only cleanup removes its registry and segment");
    } else {
        fprintf(stderr, "[SKIP] fixture cleanup checks: private POSIX shm unavailable (%s)\n", strerror(errno));
    }
    if (map != MAP_FAILED) munmap(map, BWA_SHM_CTL_SIZE);
    if (fd >= 0) close(fd);
    if (fixture_fd >= 0) close(fixture_fd);
    if (other_fd >= 0) close(other_fd);
    if (lock != SEM_FAILED) sem_close(lock);
}

/* Alias the generated disk fixture under short per-run names: POSIX staging
 * keys on basename, not content, and macOS limits shm names to 31 characters. */
struct ShmDiskFixtures {
    bwa_tests::TempDir directory{"sj"};
    std::string first, second;
    bool ready = false;

    explicit ShmDiskFixtures(const char *source)
    {
        const std::string source_path(source);
        const size_t slash = source_path.rfind('/');
        const std::string parent = slash == std::string::npos ? "." :
                                   slash == 0 ? "/" : source_path.substr(0, slash);
        const std::string basename = source_path.substr(slash == std::string::npos ? 0 : slash + 1);
        char resolved_parent[PATH_MAX];
        if (basename.empty() || realpath(parent.c_str(), resolved_parent) == NULL) return;
        const std::string resolved = std::string(resolved_parent) + "/" + basename;
        first = directory.file("first");
        second = directory.file("second");
        for (const std::string &prefix : {first, second}) {
            for (const char *suffix : {".ann", ".amb", ".pac", ".bwt.2bit.64"}) {
                const std::string dest = prefix + suffix;
                const std::string src = std::string(resolved) + suffix;
                if (symlink(src.c_str(), dest.c_str()) != 0) return;
            }
        }
        ready = true;
    }

};

/* Check the host facility before staging. After this succeeds, every stage
 * error is unexpected and fails with the implementation's captured diagnostic. */
static bool staging_prerequisite()
{
#if defined(__linux__)
    PrivateObjects probe;
    const std::string probe_name = std::string(private_control) + "p";
    const int fd = probe.acquire(probe_name) ? dup(probe.descriptor(probe_name)) : -1;
    if (fd < 0) {
        fprintf(stderr, "[SKIP] real staging prerequisite: private control shm_open: %s\n", strerror(errno));
        return false;
    }
    bool ok = ftruncate(fd, BWA_SHM_CTL_SIZE) == 0;
    void *map = ok ? mmap(NULL, BWA_SHM_CTL_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0) : MAP_FAILED;
    ok = ok && map != MAP_FAILED;
    const int err = errno;
    if (map != MAP_FAILED) munmap(map, BWA_SHM_CTL_SIZE);
    close(fd);
    probe.clear();
    if (!ok) fprintf(stderr, "[SKIP] real staging prerequisite: private control mapping: %s\n", strerror(err));
    return ok;
#else
    fprintf(stderr, "[SKIP] real staging prerequisite: supported Linux POSIX shm required; "
                    "portable private-registry checks still run\n");
    return false;
#endif
}

#if defined(__GLIBC__) && (defined(__SANITIZE_ADDRESS__) || defined(__has_feature))
#if defined(__SANITIZE_ADDRESS__)
#define PRIVATE_STAGE_ASAN 1
#elif __has_feature(address_sanitizer)
#define PRIVATE_STAGE_ASAN 1
#endif
#endif
#if defined(PRIVATE_STAGE_ASAN)
static void stage_child_exit(int status, void *) { fflush(NULL); _exit(status); }
#endif

static bool stage_required(const char *prefix, const char *name)
{
    bwa_tests::StderrCapture capture;
    fflush(NULL);
    const pid_t child = fork();
    if (child == 0) {
#if defined(PRIVATE_STAGE_ASAN)
        on_exit(stage_child_exit, NULL);
#endif
        const int rc = bwa_shm_stage(prefix);
        fflush(NULL);
        _exit(rc == 0 ? 0 : 1);
    }
    int status = 0;
    pid_t ended = -1;
    if (child > 0) do { ended = waitpid(child, &status, 0); } while (ended < 0 && errno == EINTR);
    const bool ok = child > 0 && ended == child && WIFEXITED(status) && WEXITSTATUS(status) == 0;
    const std::string diagnostic = capture.text();
    if (!ok)
        fprintf(stderr, "[ERROR] %s: unexpected staging child status %d:\n%s", name, status, diagnostic.c_str());
    report(ok, name);
    return ok;
}

// Deliberate collisions are exclusively ours, never borrowed host objects.
static void check_namespace_collisions()
{
    for (const char *prefix : {"sm", "sh", "sj"}) {
        bwa_tests::TempDir temp("collision");
        const std::string token = temp.path().substr(temp.path().rfind('.') + 1);
        const std::string base = std::string("/") + prefix + token;
        const std::string ctl = base + "c", lock = base + "l";
        const std::vector<std::string> segments = {base + "-first", base + "-second"};
        for (int collision = 0; collision < 4; ++collision) {
            PrivateObjects foreign;
            const std::string target = collision == 0 ? ctl : collision == 1 ? lock : segments[collision - 2];
            const bool created = foreign.acquire(target, collision == 1);
            report(created, "collision sentinel exclusively created");
            if (!created) continue;
            {
                PrivateObjects contender;
                const bool reserved = contender.reserve(ctl, lock, segments);
                bool kept = shm_name_exists(target.c_str());
                // Close the probe without unlinking the foreign semaphore.
                if (collision == 1) { sem_t *probe = sem_open(lock.c_str(), 0); kept = probe != SEM_FAILED; if (kept) sem_close(probe); }
                report(!reserved && !contender.owns(target) && kept,
                       "partial namespace collision rolls back only acquired objects");
            }
            bool preserved = shm_name_exists(target.c_str());
            if (collision == 1) {
                sem_t *probe = sem_open(lock.c_str(), 0);
                preserved = probe != SEM_FAILED;
                if (preserved) sem_close(probe);
            }
            report(preserved, "normal cleanup preserves foreign collision object");
            for (const auto &name : {ctl, segments[0], segments[1]})
                if (name != target) report(!shm_name_exists(name.c_str()), "partial acquisition leaves no owned segment or registry");
        }
    }
    // Cleanup must reject an unowned second segment even with an owned first.
    bwa_tests::TempDir temp("unowned_second");
    const std::string base = "/su" + temp.path().substr(temp.path().rfind('.') + 1);
    const std::string ctl = base + "c", lock = base + "l", idx = base + "-";
    PrivateObjects owner, foreign;
    const bool ready = owner.reserve(ctl, lock, {idx + "first"}) && foreign.acquire(idx + "second");
    report(ready, "unowned second-segment sentinel prepared");
    if (ready) {
        const ShmRegistry registry = {ctl.c_str(), lock.c_str(), idx.c_str(), &owner};
        drop_fixture_shm("first", registry, mmap, "second");
        report(shm_name_exists((idx + "second").c_str()) && shm_name_exists(ctl.c_str()),
               "cleanup refuses an unowned second segment and preserves registry");
    }
}

static void check_fatal_stage()
{
    const std::string saved_control(private_control), saved_lock(private_lock), saved_index(private_index);
    PrivateObjects *saved_owner = stage_objects;
    for (bool truncated : {false, true}) {
        bwa_tests::TempDir input("fatal_stage");
        const std::string token = input.path().substr(input.path().rfind('.') + 1);
        const std::string ctl = "/sf" + token + "c", lock = "/sf" + token + "l", idx = "/sf" + token + "-";
        PrivateObjects owned, sentinel;
        const bool ready = owned.reserve(ctl, lock, {idx + "first", idx + "second"})
                           && sentinel.acquire(idx + "unrelated") && sentinel.acquire(ctl + "u")
                           && sentinel.acquire(lock + "u", true);
        report(ready, "fatal stage parent owns complete namespace and unrelated sentinel");
        if (!ready) continue;
        snprintf(private_control, sizeof private_control, "%s", ctl.c_str());
        snprintf(private_lock, sizeof private_lock, "%s", lock.c_str());
        snprintf(private_index, sizeof private_index, "%s", idx.c_str());
        stage_objects = &owned;
        if (truncated) input.write("broken.ann", "");
        const int before_total = total_cases, before_passed = passed_cases;
        bwa_tests::StderrCapture capture;
        const bool staged = stage_required(input.file(truncated ? "broken" : "missing").c_str(), "fatal private stage");
        const std::string text = capture.text();
        const bool failed = !staged && total_cases == before_total + 1 && passed_cases == before_passed;
        total_cases = before_total; passed_cases = before_passed;
        report(failed && text.find("[ERROR]") != std::string::npos && text.find("[FAIL]") != std::string::npos
                   && text.find("status 256") != std::string::npos
                   && text.find(truncated ? "Unexpected end of file" : "fail to open file") != std::string::npos,
               truncated ? "truncated index fatal exit reports diagnostic and FAIL in parent" :
                           "missing index fatal exit reports diagnostic and FAIL in parent");
        owned.clear();
        sem_t *probe = sem_open(lock.c_str(), 0);
        const bool missing_lock = probe == SEM_FAILED && errno == ENOENT;
        if (probe != SEM_FAILED) sem_close(probe);
        sem_t *sentinel_probe = sem_open((lock + "u").c_str(), 0);
        const bool sentinel_lock = sentinel_probe != SEM_FAILED;
        if (sentinel_lock) sem_close(sentinel_probe);
        report(sentinel_lock && shm_name_exists((ctl + "u").c_str())
                   && !shm_name_exists(ctl.c_str()) && missing_lock && !shm_name_exists((idx + "first").c_str())
                   && !shm_name_exists((idx + "second").c_str()) && shm_name_exists((idx + "unrelated").c_str()),
               "parent cleans owned objects after fatal child exit and preserves unrelated sentinel");
    }
    snprintf(private_control, sizeof private_control, "%s", saved_control.c_str());
    snprintf(private_lock, sizeof private_lock, "%s", saved_lock.c_str());
    snprintf(private_index, sizeof private_index, "%s", saved_index.c_str());
    stage_objects = saved_owner;
}

struct Walk { int64_t k, l, s; };

/* One forward step of the seeding walks (bsd_advance_step,
 * ls_advance_forward_step) from w with base c. */
static Walk forward_step(FMI_search *fmi, const Walk &w, int c)
{
    SMEM t;
    memset(&t, 0, sizeof(t));
    t.k = w.l;
    t.l = w.k;
    t.s = w.s;
    const SMEM r = fmi->backwardExt(t, (uint8_t)(3 - c));
    Walk o = { r.l, r.k, r.s };
    return o;
}

/* Depth-first over every string of length L..D extending `code`, comparing
 * each table entry with the stepped interval w. */
static void check_subtree(FMI_search *fmi, int L, uint32_t code, const Walk &w, int D,
                          int64_t *checked, int64_t *bad)
{
    const FMI_search::SeedJumpEntry &e = fmi->seed_jump[FMI_search::seed_jump_level_off(L) + code];
    (*checked)++;
    if (e.k != w.k || e.l != w.l || e.s != w.s) {
        if (*bad < 5)
            fprintf(stderr, "  level %d code %u: table (%lld, %lld, %lld) != stepped (%lld, %lld, %lld)\n",
                    L, code, (long long)e.k, (long long)e.l, (long long)e.s,
                    (long long)w.k, (long long)w.l, (long long)w.s);
        (*bad)++;
    }
    if (L == D) return;
    for (int c = 0; c < 4; c++)
        check_subtree(fmi, L + 1, code * 4 + (uint32_t)c, forward_step(fmi, w, c), D, checked, bad);
}

/* Every entry of fmi's table at depth D. */
static bool check_table(FMI_search *fmi, int D)
{
    if (fmi->seed_jump == NULL || fmi->seed_jump_depth != D) {
        fprintf(stderr, "  table %p, depth %d; expected depth %d\n",
                (void *)fmi->seed_jump, fmi->seed_jump_depth, D);
        return false;
    }
    int64_t checked = 0, bad = 0;
    for (int a = 0; a < 4; a++) {
        Walk w = { fmi->count[a], fmi->count[3 - a], fmi->count[a + 1] - fmi->count[a] };
        check_subtree(fmi, 1, (uint32_t)a, w, D, &checked, &bad);
    }
    const int64_t want = FMI_search::seed_jump_level_off(D + 1);
    if (checked != want) {
        fprintf(stderr, "  checked %lld entries, expected %lld\n", (long long)checked, (long long)want);
        return false;
    }
    return bad == 0;
}

/* ---- Brute-force oracle over the index text -------------------------------
 * The index is built over the .pac forward strand followed by its reverse
 * complement, with a sentinel ($) after it that sorts below every base. */
struct TextOracle {
    std::vector<uint8_t> t;      // 2-bit bases, forward then reverse complement
    std::vector<int32_t> sa;     // suffix array of t + $ (the $ suffix first)
};

static bool load_text(const char *prefix, TextOracle *o)
{
    bntseq_t *bns = bns_restore(prefix);
    if (bns == NULL) return false;
    const int64_t l_pac = bns->l_pac;
    bns_destroy(bns);
    std::string pac_name = std::string(prefix) + ".pac";
    FILE *fp = fopen(pac_name.c_str(), "rb");
    if (fp == NULL) return false;
    std::vector<uint8_t> pac((size_t)(l_pac / 4 + 1));
    const size_t got = fread(pac.data(), 1, pac.size(), fp);
    fclose(fp);
    if ((int64_t)got * 4 < l_pac) return false;
    o->t.resize((size_t)(2 * l_pac));
    for (int64_t i = 0; i < l_pac; i++) {
        const uint8_t b = (pac[(size_t)(i >> 2)] >> ((~i & 3) << 1)) & 3;
        o->t[(size_t)i] = b;
        o->t[(size_t)(2 * l_pac - 1 - i)] = (uint8_t)(3 - b);
    }
    const int32_t n = (int32_t)o->t.size();
    o->sa.resize((size_t)n + 1);
    for (int32_t i = 0; i <= n; i++) o->sa[(size_t)i] = i;
    const uint8_t *t = o->t.data();
    std::sort(o->sa.begin(), o->sa.end(), [t, n](int32_t a, int32_t b) {
        while (a < n && b < n) {
            if (t[a] != t[b]) return t[a] < t[b];
            a++; b++;
        }
        return a == n && b != n;   // the shorter suffix (reaching $) sorts first
    });
    return true;
}

/* Occurrences of p in t, counted position by position. */
static int64_t brute_count(const TextOracle &o, const std::vector<uint8_t> &p)
{
    int64_t c = 0;
    const size_t n = o.t.size(), L = p.size();
    for (size_t i = 0; i + L <= n; i++)
        if (std::equal(p.begin(), p.end(), o.t.begin() + (long)i)) c++;
    return c;
}

/* Rank of the first suffix of t + $ that is not below p (p's SA interval start). */
static int64_t sa_lower(const TextOracle &o, const std::vector<uint8_t> &p)
{
    const int32_t n = (int32_t)o.t.size();
    const uint8_t *t = o.t.data();
    auto it = std::lower_bound(o.sa.begin(), o.sa.end(), 0, [&](int32_t suf, int) {
        for (size_t i = 0; i < p.size(); i++) {
            if (suf + (int32_t)i >= n) return true;   // $ sorts below every base
            const uint8_t b = t[suf + i];
            if (b != p[i]) return b < p[i];
        }
        return false;   // p is a prefix of the suffix: not below p
    });
    return (int64_t)(it - o.sa.begin());
}

/* Every entry of fmi's table (depth D) against the text oracle: s is p's
 * occurrence count; when s > 0, k is p's SA interval start and l is the start
 * of p's reverse complement (both forward-strand and reverse-complement
 * occurrences are in the text, so the two counts agree). */
static bool check_table_text(FMI_search *fmi, const TextOracle &o, int D)
{
    if (fmi->seed_jump == NULL || fmi->seed_jump_depth != D) return false;
    int64_t checked = 0, bad = 0;
    for (int L = 1; L <= D; L++) {
        for (uint32_t code = 0; code < (1u << (2 * L)); code++) {
            std::vector<uint8_t> p((size_t)L), rc((size_t)L);
            for (int i = 0; i < L; i++) p[(size_t)i] = (code >> (2 * (L - 1 - i))) & 3;
            for (int i = 0; i < L; i++) rc[(size_t)i] = (uint8_t)(3 - p[(size_t)(L - 1 - i)]);
            const FMI_search::SeedJumpEntry &e = fmi->seed_jump[FMI_search::seed_jump_level_off(L) + code];
            const int64_t s = brute_count(o, p);
            bool ok = e.s == s;
            if (ok && s > 0) ok = e.k == sa_lower(o, p) && e.l == sa_lower(o, rc);
            checked++;
            if (!ok) {
                if (bad < 5)
                    fprintf(stderr, "  level %d code %u: table (%lld, %lld, %lld), text s %lld k %lld l %lld\n",
                            L, code, (long long)e.k, (long long)e.l, (long long)e.s, (long long)s,
                            (long long)sa_lower(o, p), (long long)sa_lower(o, rc));
                bad++;
            }
        }
    }
    return checked == FMI_search::seed_jump_level_off(D + 1) && bad == 0;
}

/* ---- Does a lockstep driver read the table? --------------------------------
 * Reads: phiX windows (forward and reverse-complement text) and the jump-depth
 * reads, 2-bit encoded, N = 4. */
struct Reads {
    std::vector<uint8_t> enc;
    std::vector<bseq1_t> seq;
    std::vector<int32_t> cum;
    int32_t max_len = 0;
};

static Reads make_reads(const TextOracle &o)
{
    std::vector<std::vector<uint8_t>> rs;
    for (size_t st = 0; st + 60 <= o.t.size() && rs.size() < 40; st += o.t.size() / 40)
        rs.push_back(std::vector<uint8_t>(o.t.begin() + (long)st, o.t.begin() + (long)st + 60));
    for (const std::string &r : jump_depth_reads()) {
        std::vector<uint8_t> e;
        for (char c : r) e.push_back(c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : 4);
        rs.push_back(e);
    }
    Reads R;
    R.seq.resize(rs.size());
    memset(R.seq.data(), 0, rs.size() * sizeof(bseq1_t));
    for (size_t i = 0; i < rs.size(); i++) {
        R.cum.push_back((int32_t)R.enc.size());
        R.enc.insert(R.enc.end(), rs[i].begin(), rs[i].end());
        R.seq[i].l_seq = (int)rs[i].size();
        R.max_len = std::max(R.max_len, (int32_t)rs[i].size());
    }
    R.cum.push_back((int32_t)R.enc.size());
    return R;
}

static bool same_smems(const std::vector<SMEM> &a, const std::vector<SMEM> &b)
{
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++)
        if (a[i].rid != b[i].rid || a[i].m != b[i].m || a[i].n != b[i].n
            || a[i].k != b[i].k || a[i].l != b[i].l || a[i].s != b[i].s
            || a[i].info != b[i].info)
            return false;
    return true;
}

/* The SMEM lockstep driver's output: one forward walk from x = 0 per read. */
static std::vector<SMEM> run_smem_lockstep(FMI_search *fmi, Reads &R, int32_t min_intv)
{
    const int32_t n = (int32_t)R.seq.size();
    std::vector<int32_t> qpa((size_t)n, 0), mia((size_t)n, min_intv), rid((size_t)n);
    for (int32_t i = 0; i < n; i++) rid[(size_t)i] = i;
    std::vector<SMEM> out((size_t)n * (size_t)R.max_len);
    int64_t n_out = 0;
    fmi->getSMEMsOnePosOneThread_lockstep(R.enc.data(), qpa.data(), mia.data(), rid.data(),
                                          n, n, R.seq.data(), R.cum.data(), R.max_len,
                                          /*minSeedLen=*/19, out.data(), &n_out);
    out.resize((size_t)n_out);
    return out;
}

/* The bwtSeed lockstep driver's output at the default minSeedLen. */
static std::vector<SMEM> run_bwtseed_lockstep(FMI_search *fmi, Reads &R, int32_t max_intv)
{
    const int32_t n = (int32_t)R.seq.size();
    std::vector<int32_t> mia((size_t)n, max_intv);
    std::vector<SMEM> out((size_t)n * (size_t)R.max_len);
    const int64_t n_out = fmi->bwtSeedStrategyAllPosOneThread_lockstep(
        R.enc.data(), mia.data(), n, R.seq.data(), R.cum.data(), /*minSeedLen=*/19,
        out.data(), R.max_len);
    out.resize((size_t)n_out);
    return out;
}

/* Replace every entry of levels 2..depth with its level's code-0 entry (the
 * interval of a run of As): a valid interval, so a walk that reads it steps on
 * safely, but the wrong one for every other string. */
static void poison_table(FMI_search *fmi)
{
    for (int L = 2; L <= fmi->seed_jump_depth; L++) {
        FMI_search::SeedJumpEntry *lev = fmi->seed_jump + FMI_search::seed_jump_level_off(L);
        for (int64_t c = 1; c < ((int64_t)1 << (2 * L)); c++) lev[c] = lev[0];
    }
}

/* FmiSeed is opaque in fmi_seed_api.h; this is its definition in
 * fmi_seed_api.cpp (identical, as the one-definition rule requires), so the
 * test can reach the handle's FMI_search. */
struct FmiSeed { FMI_search *fmi; };

// Reserve and initialize the main private namespace in this process. Called only
// after every fatal-capable disk load, so an invalid index exits owning nothing.
static bool reserve_main_namespace(PrivateObjects &owned)
{
    bool reserved = false;
    for (int attempt = 0; attempt < 20 && !reserved; ++attempt) {
        bwa_tests::TempDir private_names("seed_jump");
        const std::string token = private_names.path().substr(private_names.path().rfind('.') + 1);
        snprintf(private_control, sizeof private_control, "/sm%sc", token.c_str());
        snprintf(private_lock, sizeof private_lock, "/sm%sl", token.c_str());
        snprintf(private_index, sizeof private_index, "/sm%s-", token.c_str());
        reserved = owned.reserve(private_control, private_lock,
                                 {std::string(private_index) + "first", std::string(private_index) + "second"});
    }
    report(reserved, "main reserves complete private namespace after disk loads, before staging");
    if (!reserved) return false;
    exit_owner = &owned;
    // A signal, _exit or abort after this point leaks these names; name them.
    fprintf(stderr, "[INFO] private namespace reserved: %s %s %sfirst %ssecond (pid %d)\n",
            private_control, private_lock, private_index, private_index, (int)getpid());
    const int ctl = owned.descriptor(private_control);
    const uint16_t header[2] = {0, 4};
    void *control_map = ftruncate(ctl, BWA_SHM_CTL_SIZE) == 0 ?
        mmap(NULL, BWA_SHM_CTL_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, ctl, 0) : MAP_FAILED;
    const bool initialized = control_map != MAP_FAILED;
    if (initialized) { memcpy(control_map, header, sizeof header); munmap(control_map, BWA_SHM_CTL_SIZE); }
    report(initialized, "reserved registry initialized for real staging");
    if (!initialized) return false;
    fixture_registry.owner = &owned;
    stage_objects = &owned;
    return true;
}

// The owner unlinks its whole namespace itself, then checks every name is gone.
static bool release_main_namespace(PrivateObjects &owned)
{
    owned.clear();
    exit_owner = NULL;
    stage_objects = NULL;
    fixture_registry.owner = NULL;
    sem_t *lock_probe = sem_open(private_lock, 0);
    const bool no_lock = lock_probe == SEM_FAILED && errno == ENOENT;
    if (lock_probe != SEM_FAILED) sem_close(lock_probe);
    const bool cleaned = no_lock && !shm_name_exists(private_control)
                         && !shm_name_exists((std::string(private_index) + "first").c_str())
                         && !shm_name_exists((std::string(private_index) + "second").c_str());
    report(cleaned, "owner removes complete main namespace after full workload");
    return cleaned;
}

static int fixture_workload(const char *prefix, bool portable_checks)
{
    if (portable_checks) {
    check_namespace_collisions();
    check_fatal_stage();
    check_fixture_cleanup();
    {
        // An actual rejected stage must count as FAIL, not a host-unavailable SKIP.
        const int before_total = total_cases, before_passed = passed_cases;
        bwa_tests::StderrCapture capture;
        const bool staged = stage_required(NULL, "second private fixture staging");
        const std::string diagnostic = capture.text();
        const bool failed = !staged && total_cases == before_total + 1 && passed_cases == before_passed;
        total_cases = before_total;
        passed_cases = before_passed;
        report(failed && diagnostic.find("[ERROR]") != std::string::npos
                   && diagnostic.find("[FAIL]") != std::string::npos
                   && diagnostic.find("empty prefix") != std::string::npos
                   && diagnostic.find("[SKIP]") == std::string::npos,
               "unexpected second-stage errors fail with the real diagnostic, never SKIP");
    }
    {
        bwa_tests::TempDir suffix_only("suffix_only");
        for (const char *suffix : {".ann", ".amb", ".pac", ".bwt.2bit.64"})
            suffix_only.write(std::string("index") + suffix, "fixture");
        ShmDiskFixtures aliases(suffix_only.file("index").c_str());
        report(aliases.ready && access(aliases.first.c_str(), F_OK) != 0
                   && access((aliases.first + ".ann").c_str(), R_OK) == 0,
               "suffix-only prefix resolves through directory and basename");
    }
    }
    if (strcmp(prefix, "--fixtures-only") == 0)
        return passed_cases == total_cases ? 0 : 1;
    char name[160];

    /* Default load: depth 10, built by load_index. */
    set_jump_env(NULL);
    FMI_search *fmi = new FMI_search(prefix);
    fmi->load_index();
    report(check_table(fmi, 10), "load_index builds the default depth-10 table");

    /* Rebuilds at the depths around the serial prefix (3) and the default,
     * at several thread counts, including more threads than the 64 prefixes
     * the parallel levels are split by (so the builder clamps the count). */
    const struct { int depth; int threads; } builds[] = {
        {2, 1}, {2, 4}, {3, 1}, {3, 4}, {4, 1}, {4, 3}, {5, 64}, {5, 100}, {10, 1}, {10, 7},
    };
    for (const auto &b : builds) {
        char v[8];
        snprintf(v, sizeof v, "%d", b.depth);
        set_jump_env(v);
        fmi->build_seed_jump_table(b.threads);
        snprintf(name, sizeof name, "depth %d, %d builder thread(s): every entry equals the stepped walk",
                 b.depth, b.threads);
        report(check_table(fmi, b.depth), name);
    }

    /* Depths 2-4 against the brute-force text oracle (no backwardExt). */
    TextOracle text;
    const bool have_text = load_text(prefix, &text);
    report(have_text, "index text and naive suffix array loaded from the .pac");
    for (int D = 2; D <= 4; D++) {
        char v[8];
        snprintf(v, sizeof v, "%d", D);
        set_jump_env(v);
        fmi->build_seed_jump_table(2);
        snprintf(name, sizeof name, "depth %d: every entry's s, k and l equal the index text's", D);
        report(have_text && check_table_text(fmi, text, D), name);
    }

    /* The lockstep drivers read the table at the default depth: replacing its
     * entries with other valid intervals must change their output. Without the
     * table both drivers step every walk, and the output is the same either
     * way. */
    if (have_text) {
        set_jump_env(NULL);
        Reads R = make_reads(text);
        fmi->build_seed_jump_table(2);
        const std::vector<SMEM> smem_ok = run_smem_lockstep(fmi, R, 1);
        const std::vector<SMEM> seed_ok = run_bwtseed_lockstep(fmi, R, 1000);
        poison_table(fmi);
        const std::vector<SMEM> smem_bad = run_smem_lockstep(fmi, R, 1);
        const std::vector<SMEM> seed_bad = run_bwtseed_lockstep(fmi, R, 1000);
        report(fmi->seed_jump_depth == 10 && !smem_ok.empty() && !same_smems(smem_ok, smem_bad),
               "SMEM lockstep driver takes its walks' first steps from the depth-10 table");
        report(fmi->seed_jump_depth == 10 && !seed_ok.empty() && !same_smems(seed_ok, seed_bad),
               "bwtSeed lockstep driver takes its walks' first steps from the depth-10 table");

        set_jump_env("0");
        fmi->build_seed_jump_table(2);
        report(fmi->seed_jump == NULL && same_smems(run_smem_lockstep(fmi, R, 1), smem_ok)
                   && same_smems(run_bwtseed_lockstep(fmi, R, 1000), seed_ok),
               "BWA3_BWTSEED_JUMP=0: both lockstep drivers' output equals the table's");
    } else {
        report(false, "lockstep drivers read the table: no index text");
    }

    /* BWA3_BWTSEED_JUMP=0: no table. */
    set_jump_env("0");
    fmi->build_seed_jump_table(4);
    report(fmi->seed_jump == NULL && fmi->seed_jump_depth == 0, "BWA3_BWTSEED_JUMP=0: no table");

    /* Empty, malformed and out-of-range values: the default depth. */
    const char *defaults[] = {"", "1", "13", "-1", "-0", "+5", " 5", "5 ", "10x", "abc",
                              "99999999999999999999"};
    for (const char *v : defaults) {
        set_jump_env(v);
        fmi->build_seed_jump_table(2);
        snprintf(name, sizeof name, "BWA3_BWTSEED_JUMP=\"%s\": default depth 10", v);
        report(fmi->seed_jump != NULL && fmi->seed_jump_depth == 10, name);
    }
    set_jump_env(NULL);

    /* The build_seed_jump=false opt-out, on load_index and on the facade that
     * passes it (fmi_seed_open). */
    {
        FMI_search *np = new FMI_search(prefix);
        np->load_index(/*load_pac=*/false, /*n_threads=*/1, /*build_seed_jump=*/false);
        report(np->seed_jump == NULL && np->seed_jump_depth == 0,
               "load_index(build_seed_jump=false): no table");
        delete np;
        FmiSeed *h = fmi_seed_open(prefix);
        report(h != NULL && h->fmi->seed_jump == NULL && h->fmi->seed_jump_depth == 0,
               "fmi_seed_open (disk or already-staged attach): no table");
        fmi_seed_close(h);
    }

    /* An index attached from a packed segment: load_index_from_shm builds no
     * table, and build_seed_jump_table over the attached arrays builds one
     * equal to the disk load's. The packed segment is copied into an anonymous
     * mapping because the attached FMI_search munmaps its base on destruction. */
    {
        uint8_t *buf = NULL;
        uint64_t len = 0;
        bool ok = bwa_shm_pack_from_disk(prefix, &buf, &len) == 0 && buf != NULL;
        void *map = MAP_FAILED;
        if (ok) {
            map = mmap(NULL, (size_t)len, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0);
            ok = map != MAP_FAILED;
        }
        bool none = false;
        if (ok) {
            memcpy(map, buf, (size_t)len);
            fmi->build_seed_jump_table(1);   // the disk load's default table again
            FMI_search *at = new FMI_search(prefix);
            at->load_index_from_shm((uint8_t *)map, (size_t)len);
            none = at->seed_jump == NULL && at->seed_jump_depth == 0;
            at->build_seed_jump_table(3);
            const size_t bytes = (size_t)FMI_search::seed_jump_level_off(11) * sizeof(FMI_search::SeedJumpEntry);
            ok = check_table(at, 10) && fmi->seed_jump != NULL
                 && memcmp(at->seed_jump, fmi->seed_jump, bytes) == 0;
            delete at;   // munmaps map
        }
        report(none, "load_index_from_shm alone: no table");
        report(ok, "attached index: table equals the disk load's");
        free(buf);
    }

    /* Real staging uses the exclusively reserved private namespace.
     * Alias prefixes have only the required suffix files; .0123 is not staged.
     * The scratch fixture is built first: if it throws, no IPC name is owned yet. */
    ShmDiskFixtures fixture(prefix);
    PrivateObjects owned;
    if (!reserve_main_namespace(owned)) {
        release_main_namespace(owned);
        delete fmi;
        return 1;
    }
    report(fixture.ready, "real staging: suffix-only disk aliases prepared");
    if (fixture.ready && staging_prerequisite()) {
        const char *staged_prefix = fixture.first.c_str();
        const char *second_prefix = fixture.second.c_str();
        const std::string first_name = strrchr(staged_prefix, '/') + 1;
        const std::string second_name = strrchr(second_prefix, '/') + 1;
        const bool staged = stage_required(staged_prefix, "first private fixture staging");
        if (staged) {
            FMI_search *at = new FMI_search(staged_prefix);
            at->load_index();
            report(at->shm_base != NULL && check_table(at, 10),
                   "load_index (suffix-only shm attach): builds the default table");
            delete at;
            at = new FMI_search(staged_prefix);
            at->load_index(false, 1, false);
            report(at->shm_base != NULL && at->seed_jump == NULL && at->seed_jump_depth == 0,
                   "load_index(build_seed_jump=false) (shm attach): no table");
            delete at;
            FmiSeed *h = fmi_seed_open(staged_prefix);
            report(h != NULL && h->fmi->shm_base != NULL && h->fmi->seed_jump == NULL
                       && h->fmi->seed_jump_depth == 0,
                   "fmi_seed_open (already-staged shm attach): no table");
            fmi_seed_close(h);
            if (stage_required(second_prefix, "second private fixture staging")) {
                drop_fixture_shm(staged_prefix);
                std::vector<std::string> names;
                report(read_shm_registry(&names) == REGISTRY_READABLE
                           && names == std::vector<std::string>{first_name, second_name}
                           && bwa_shm_test(staged_prefix) == 1 && bwa_shm_test(second_prefix) == 1
                           && shm_name_exists((std::string(private_index) + first_name).c_str())
                           && shm_name_exists((std::string(private_index) + second_name).c_str()),
                       "fixture cleanup preserves a second staged index and its segment");
            }
        }
        drop_fixture_shm(staged_prefix, fixture_registry, mmap, second_name.c_str());
        std::vector<std::string> names;
        report(read_shm_registry(&names) == REGISTRY_MISSING && names.empty()
                   && !shm_name_exists((std::string(private_index) + first_name).c_str())
                   && !shm_name_exists((std::string(private_index) + second_name).c_str()),
               "real staging cleanup removes both owned segments and their registry entries");
    }

    release_main_namespace(owned);
    delete fmi;
    fprintf(stderr, "%d / %d cases passed\n", passed_cases, total_cases);
    return passed_cases == total_cases ? 0 : 1;
}

static bool lock_exists(const std::string &name)
{
    sem_t *probe = sem_open(name.c_str(), 0);
    if (probe == SEM_FAILED) return false;
    sem_close(probe);
    return true;
}

static int wait_status(pid_t child)
{
    int status = 0;
    pid_t ended;
    do { ended = waitpid(child, &status, 0); } while (ended < 0 && errno == EINTR);
    return ended == child && WIFEXITED(status) ? WEXITSTATUS(status) : -1;
}

// Unique private names for one regression; the test exclusively creates each.
static std::string private_base(const char *tag, const char *prefix)
{
    bwa_tests::TempDir names(tag);
    return std::string("/") + prefix + names.path().substr(names.path().rfind('.') + 1);
}

// A missing or truncated index exits from the disk loader before the main
// namespace is reserved. The workload runs in a child here only as a harness for
// that fatal exit: the child owns nothing when it dies, which is what is proven.
static void check_invalid_full_workload()
{
    for (bool corrupt : {false, true}) {
        bwa_tests::TempDir input("invalid_full");
        // A valid first scalar followed by a truncated file reaches a different
        // fatal read than the missing-file exit, without constructing an index.
        if (corrupt) {
            const int64_t length = 101;
            input.write("broken.bwt.2bit.64", std::string(reinterpret_cast<const char *>(&length), sizeof length));
        }
        PrivateObjects unrelated;
        const std::string base = private_base("invalid_full", "si");
        const std::string ctl = base + "c", lock = base + "l", segment = base + "-sentinel";
        const bool ready = unrelated.reserve(ctl, lock, {segment});
        report(ready, "invalid full-input unrelated sentinels reserved");
        if (!ready) continue;
        const std::string prefix = input.file(corrupt ? "broken" : "missing");
        bwa_tests::StderrCapture capture;
        fflush(NULL);
        const pid_t child = fork();
        if (child == 0) {
#if defined(PRIVATE_STAGE_ASAN)
            on_exit(stage_child_exit, NULL);
#endif
            const int rc = fixture_workload(prefix.c_str(), false);
            fflush(NULL);
            _exit(rc);
        }
        const int rc = child > 0 ? wait_status(child) : -1;
        const std::string diagnostic = capture.text();
        report(rc > 0 && diagnostic.find(corrupt ? "Unexpected end of file" : "Unable to open the file") != std::string::npos
                   && diagnostic.find("[INFO] private namespace reserved") == std::string::npos
                   && shm_name_exists(ctl.c_str()) && lock_exists(lock) && shm_name_exists(segment.c_str()),
               corrupt ? "corrupt full input exits with its diagnostic before reserving the main namespace" :
                         "missing full input exits with its diagnostic before reserving the main namespace");
    }
}

// A forked child runs the real bwa_shm_destroy over an owned private namespace:
// its production shm_unlink/sem_unlink calls are suppressed, the owner keeps
// every name and later removes them itself. Unrelated names are untouched.
static void check_child_unlink_suppressed()
{
    const std::string base = private_base("child_unlink", "sc");
    const std::string ctl = base + "c", lock = base + "l", idx = base + "-";
    const std::string segment = idx + "owned", foreign_name = base + "-foreign";
    PrivateObjects owner, foreign;
    bool ready = owner.reserve(ctl, lock, {segment}) && foreign.acquire(foreign_name);
    // Register the owned segment so destroy's registry walk unlinks it too.
    void *map = ready && ftruncate(owner.descriptor(ctl), BWA_SHM_CTL_SIZE) == 0 ?
        mmap(NULL, BWA_SHM_CTL_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, owner.descriptor(ctl), 0) : MAP_FAILED;
    ready = ready && map != MAP_FAILED;
    if (ready) {
        const uint16_t count = 1, next = (uint16_t)(4 + sizeof(int64_t) + strlen("owned") + 1);
        const int64_t bytes = 1;
        memset(map, 0, BWA_SHM_CTL_SIZE);
        memcpy(map, &count, sizeof count);
        memcpy((char *)map + 2, &next, sizeof next);
        memcpy((char *)map + 4, &bytes, sizeof bytes);
        memcpy((char *)map + 4 + sizeof bytes, "owned", strlen("owned") + 1);
        munmap(map, BWA_SHM_CTL_SIZE);
    }
    report(ready, "child unlink: owned private namespace and foreign sentinel exclusively created");
    if (!ready) return;
    const std::string saved_control(private_control), saved_lock(private_lock), saved_index(private_index);
    PrivateObjects *saved = stage_objects;
    snprintf(private_control, sizeof private_control, "%s", ctl.c_str());
    snprintf(private_lock, sizeof private_lock, "%s", lock.c_str());
    snprintf(private_index, sizeof private_index, "%s", idx.c_str());
    stage_objects = &owner;
    int attempts[2] = {-1, -1};
    const bool piped = pipe(attempts) == 0;
    fflush(NULL);
    const pid_t child = piped ? fork() : -1;
    if (child == 0) {
        close(attempts[0]);
        suppressed_unlink_fd = attempts[1];
        _exit(bwa_shm_destroy() == 0 ? 0 : 1);
    }
    if (piped) close(attempts[1]);
    const int rc = child > 0 ? wait_status(child) : -1;
    // The child has exited, so its write end is closed: read every attempt.
    std::string attempted;
    char buf[256];
    ssize_t got;
    while (piped && (got = read(attempts[0], buf, sizeof buf)) > 0) attempted.append(buf, (size_t)got);
    if (piped) close(attempts[0]);
    report(attempted == segment + "\n" + ctl + "\n" + lock + "\n",
           "child destroy attempts to unlink the segment, registry and lock, in that order");
    report(rc == 0 && shm_name_exists(segment.c_str()) && shm_name_exists(ctl.c_str()) && lock_exists(lock)
               && owner.owns(segment) && owner.owns(ctl) && owner.owns(lock),
           "child production unlinks of owned names are suppressed; names still exist and are owned");
    owner.clear();
    stage_objects = saved;
    snprintf(private_control, sizeof private_control, "%s", saved_control.c_str());
    snprintf(private_lock, sizeof private_lock, "%s", saved_lock.c_str());
    snprintf(private_index, sizeof private_index, "%s", saved_index.c_str());
    report(!shm_name_exists(segment.c_str()) && !shm_name_exists(ctl.c_str()) && !lock_exists(lock)
               && shm_name_exists(foreign_name.c_str()),
           "owner cleanup removes its names and preserves a foreign sentinel");
}

// remove() from a non-owner is a hard failure that unlinks nothing.
static void check_non_owner_remove_fatal()
{
    const std::string segment = private_base("non_owner", "sn") + "-owned";
    PrivateObjects owner;
    const bool ready = owner.acquire(segment);
    report(ready, "non-owner remove: owner segment exclusively created");
    if (!ready) return;
    bwa_tests::StderrCapture capture;
    fflush(NULL);
    const pid_t child = fork();
    if (child == 0) { owner.remove(segment); _exit(0); }
    const int rc = child > 0 ? wait_status(child) : -1;
    const std::string diagnostic = capture.text();
    report(rc == NON_OWNER_REMOVE_EXIT && diagnostic.find("[FATAL] PrivateObjects::remove") != std::string::npos
               && shm_name_exists(segment.c_str()) && owner.owns(segment),
           "non-owner remove exits with its distinct code and leaves the name");
    owner.clear();
    report(!shm_name_exists(segment.c_str()), "owner removes the name after a rejected non-owner remove");
}

// exit() in an owner runs its exit handler, which removes its names. A child of
// that owner inherits the handler, but its plain exit must not unlink anything.
static void check_exit_handler()
{
    const std::string base = private_base("exit_handler", "se");
    const std::string ctl = base + "c", lock = base + "l", segment = base + "-first";
    for (bool grandchild : {false, true}) {
        int ready[2];
        if (pipe(ready) != 0) { report(false, "exit handler: notification pipe"); return; }
        fflush(NULL);
        const pid_t child = fork();
        if (child == 0) {
            close(ready[0]);
            // This child becomes the only owner of its own names. No destructor
            // runs on exit(), so only main's registered handler can clean up; the
            // global exit_owner keeps the object reachable.
            exit_owner = new PrivateObjects();
            if (!exit_owner->reserve(ctl, lock, {segment})) _exit(95);
            if (!grandchild) exit(41);
            fflush(NULL);
            const pid_t inner = fork();
            if (inner == 0) exit(0);   // inherited handler: PID guard must skip
            const int inner_rc = inner > 0 ? wait_status(inner) : -1;
            const bool kept = shm_name_exists(ctl.c_str()) && lock_exists(lock) && shm_name_exists(segment.c_str());
            const char verdict = inner_rc == 0 && kept ? 'k' : 'x';
            if (write(ready[1], &verdict, 1) != 1) _exit(97);
            exit(42);
        }
        close(ready[1]);
        char verdict = 0;
        const bool notified = child > 0 && read(ready[0], &verdict, 1) == 1;
        close(ready[0]);
        const int rc = child > 0 ? wait_status(child) : -1;
        const bool absent = !shm_name_exists(ctl.c_str()) && !lock_exists(lock) && !shm_name_exists(segment.c_str());
        if (grandchild)
            report(notified && verdict == 'k' && rc == 42 && absent,
                   "owner's child exit(0) keeps the owner's names; owner exit then removes them");
        else
            report(rc == 41 && absent, "owner exit(41) runs the exit handler and leaves no names");
    }
}

// Destroying the registered owner unregisters it, so a later exit() cannot
// dereference it after an exception unwinds the owning frame.
static void check_owner_destructor_unregisters()
{
    PrivateObjects *const saved = exit_owner;
    {
        PrivateObjects scoped;
        exit_owner = &scoped;
    }
    report(exit_owner == NULL, "destroying the exit owner unregisters it from the exit handler");
    exit_owner = saved;
}

static int run_checks(const char *prefix)
{
    check_child_unlink_suppressed();
    check_non_owner_remove_fatal();
    check_exit_handler();
    check_owner_destructor_unregisters();
    check_invalid_full_workload();
    return fixture_workload(prefix, true);
}

int main(int argc, char *argv[])
{
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <bwa-mem3 index prefix>\n", argv[0]);
        return 2;
    }
    // Covers a plain exit() in the owner after reservation; PID-guarded.
    if (atexit(owner_exit_cleanup) != 0) return 1;
    // Catching here unwinds every stack owner, so its destructor unlinks its
    // names; an uncaught exception may terminate without unwinding.
    try {
        return run_checks(argv[1]);
    } catch (const std::exception &e) {
        fprintf(stderr, "[FAIL] unexpected exception: %s\n", e.what());
    } catch (...) {
        fprintf(stderr, "[FAIL] unexpected non-standard exception\n");
    }
    return 1;
}
