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
};
static const ShmRegistry fixture_registry = {
    BWA_SHM_CTL_NAME, BWA_SHM_LOCK_NAME, BWA_SHM_IDX_PREFIX
};
typedef void *(*RegistryMap)(void *, size_t, int, int, int, off_t);

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
 * its lock, this holds the lock the stager takes (BWA_SHM_LOCK_NAME), so no
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
                             const char *second_fixture = NULL)
{
    const char *slash = strrchr(prefix, '/');
    const char *fixture = slash != NULL ? slash + 1 : prefix;
    std::vector<std::string> names;

    sem_t *lock = sem_open(registry.lock, 0);
    if (lock == SEM_FAILED) {
        /* Every stager creates the lock before the registry. */
        if (read_shm_registry(&names, registry, map_registry) != REGISTRY_MISSING)
            fprintf(stderr, "[WARN] shm registry %s exists but its lock %s does not; "
                    "left it staged (`bwa-mem3 shm -d` removes it)\n",
                    registry.control, registry.lock);
        return;
    }
    /* Up to 10 s, polled, rather than sem_wait: a stager that crashed holding
     * the lock leaves it taken until `bwa-mem3 shm -d`, and the test must not
     * hang on it. */
    bool locked = false;
    int err = 0;
    for (int tries = 0; tries < 1000 && !locked; tries++) {
        locked = sem_trywait(lock) == 0;
        err = locked ? 0 : errno;
        if (!locked && err != EAGAIN && err != EINTR) break;
        if (!locked) usleep(10000);
    }
    if (!locked) {
        fprintf(stderr, "[WARN] could not take the shm lock %s (%s); left %s staged "
                "(`bwa-mem3 shm -d` removes it)\n", registry.lock,
                err == EAGAIN || err == EINTR ? "still held after 10 s" : strerror(err), fixture);
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
        shm_unlink(idx.c_str());   /* absent after an aborted stage */
        if (second_fixture != NULL)
            shm_unlink((std::string(registry.index_prefix) + second_fixture).c_str());
        shm_unlink(registry.control);
        sem_post(lock);
        sem_close(lock);
        return;
    }
    sem_post(lock);
    sem_close(lock);
    std::string staged;
    for (const std::string &n : names) staged += (staged.empty() ? "" : ", ") + n;
    fprintf(stderr, "[WARN] another index was staged in shm while this test ran; dropped "
            "nothing, so %s may still be staged (registry: %s). `bwa-mem3 shm -d` removes "
            "every staged index\n", fixture, staged.c_str());
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

static void check_fixture_cleanup()
{
    char unique[] = "/tmp/bwa-sj-XXXXXX";
    char *dir = mkdtemp(unique);
    if (dir == NULL) { report(false, "cleanup fixture: unique names allocated"); return; }
    const std::string token = std::string("/sj") + (strrchr(dir, '-') + 1);
    const std::string ctl = token + "c", lock_name = token + "l", idx = token + "-";
    const ShmRegistry registry = {ctl.c_str(), lock_name.c_str(), idx.c_str()};
    const char *fixture = "fixture", *other = "other";
    const std::string fixture_idx = idx + fixture, other_idx = idx + other;
    std::vector<std::string> names;
    report(read_shm_registry(&names, registry) == REGISTRY_MISSING && names.empty(),
           "fixture registry distinguishes a missing registry");
    sem_t *lock = sem_open(registry.lock, O_CREAT | O_EXCL, 0600, 1);
    int fd = shm_open(registry.control, O_CREAT | O_EXCL | O_RDWR, 0600);
    if (fd >= 0)
        report(read_shm_registry(&names, registry) == REGISTRY_READABLE && names.empty(),
               "fixture registry treats a short registry as empty");
    int fixture_fd = shm_open(fixture_idx.c_str(), O_CREAT | O_EXCL | O_RDWR, 0600);
    int other_fd = shm_open(other_idx.c_str(), O_CREAT | O_EXCL | O_RDWR, 0600);
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
        /* Restore a fixture-only registry even when the old cleanup removed its
         * name, so the lock-lifetime regression remains independently testable. */
        int restored = shm_open(registry.control, O_CREAT | O_RDWR, 0600);
        if (restored >= 0) close(restored);
        count = 1;
        memcpy(map, &count, sizeof(count));
        drop_fixture_shm(fixture, registry);
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
                   && shm_name_exists(other_idx.c_str()),
               "fixture-only cleanup removes its registry and segment, not other segments");
    } else {
        fprintf(stderr, "[SKIP] fixture cleanup checks: private POSIX shm unavailable (%s)\n", strerror(errno));
    }
    if (map != MAP_FAILED) munmap(map, BWA_SHM_CTL_SIZE);
    if (fd >= 0) close(fd);
    if (fixture_fd >= 0) close(fixture_fd);
    if (other_fd >= 0) close(other_fd);
    if (lock != SEM_FAILED) sem_close(lock);
    if (fd >= 0) shm_unlink(registry.control);
    if (fixture_fd >= 0) shm_unlink(fixture_idx.c_str());
    if (other_fd >= 0) shm_unlink(other_idx.c_str());
    if (lock != SEM_FAILED) sem_unlink(registry.lock);
    rmdir(dir);
}

/* Alias the generated disk fixture under short per-run names: POSIX staging
 * keys on basename, not content, and macOS limits shm names to 31 characters. */
struct ShmDiskFixtures {
    char directory[64] = "/tmp/bwa-sj-XXXXXX";
    std::string first, second;
    std::vector<std::string> links;
    bool ready = false;

    explicit ShmDiskFixtures(const char *source)
    {
        char resolved[PATH_MAX];
        if (realpath(source, resolved) == NULL || mkdtemp(directory) == NULL) return;
        const std::string token = strrchr(directory, '-') + 1;
        first = std::string(directory) + "/sj" + token + "a";
        second = std::string(directory) + "/sj" + token + "b";
        for (const std::string &prefix : {first, second}) {
            for (const char *suffix : {".ann", ".amb", ".pac", ".bwt.2bit.64"}) {
                const std::string dest = prefix + suffix;
                const std::string src = std::string(resolved) + suffix;
                if (symlink(src.c_str(), dest.c_str()) != 0) return;
                links.push_back(dest);
            }
        }
        ready = true;
    }

    ~ShmDiskFixtures()
    {
        for (const std::string &link : links) unlink(link.c_str());
        if (!first.empty()) rmdir(directory);
    }
};

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

int main(int argc, char *argv[])
{
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <bwa-mem3 index prefix>\n", argv[0]);
        return 2;
    }
    check_fixture_cleanup();
    const char *prefix = argv[1];
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
        report(h != NULL && h->fmi->shm_base == NULL && h->fmi->seed_jump == NULL
                   && h->fmi->seed_jump_depth == 0,
               "fmi_seed_open (disk load): no table");
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

    /* Use unique basenames even if another process stages an index while the
     * test runs. Cleanup leaves the existing stager lock linked and preserves
     * any registry containing an index not owned by this run. */
    ShmDiskFixtures fixture(prefix);
    int ctl_fd = shm_open(BWA_SHM_CTL_NAME, O_RDONLY, 0);
    if (ctl_fd >= 0) {
        close(ctl_fd);
        fprintf(stderr, "[SKIP] load_index's shm-attach path: a registry already exists; not changing it\n");
    } else if (errno != ENOENT || !fixture.ready) {
        fprintf(stderr, "[SKIP] load_index's shm-attach path: cannot prepare an isolated fixture\n");
    } else if (bwa_shm_stage(fixture.first.c_str()) == 0) {
        const char *staged_prefix = fixture.first.c_str();
        FMI_search *at = new FMI_search(staged_prefix);
        at->load_index();
        report(at->shm_base != NULL && check_table(at, 10),
               "load_index (shm attach): builds the default table");
        delete at;
        at = new FMI_search(staged_prefix);
        at->load_index(/*load_pac=*/false, /*n_threads=*/1, /*build_seed_jump=*/false);
        report(at->shm_base != NULL && at->seed_jump == NULL && at->seed_jump_depth == 0,
               "load_index(build_seed_jump=false) (shm attach): no table");
        delete at;
        FmiSeed *h = fmi_seed_open(staged_prefix);
        report(h != NULL && h->fmi->shm_base != NULL && h->fmi->seed_jump == NULL
                   && h->fmi->seed_jump_depth == 0,
               "fmi_seed_open (shm attach): no table");
        fmi_seed_close(h);
        const char *second_prefix = fixture.second.c_str();
        if (bwa_shm_stage(second_prefix) == 0) {
            drop_fixture_shm(staged_prefix);
            std::vector<std::string> names;
            const std::string first_name = strrchr(staged_prefix, '/') + 1;
            const std::string second_name = strrchr(second_prefix, '/') + 1;
            report(read_shm_registry(&names) == REGISTRY_READABLE
                       && std::find(names.begin(), names.end(), first_name) != names.end()
                       && std::find(names.begin(), names.end(), second_name) != names.end()
                       && bwa_shm_test(staged_prefix) == 1 && bwa_shm_test(second_prefix) == 1
                       && shm_name_exists((std::string(BWA_SHM_IDX_PREFIX) + first_name).c_str())
                       && shm_name_exists((std::string(BWA_SHM_IDX_PREFIX) + second_name).c_str()),
                   "fixture cleanup preserves a second staged index and its segment");
        } else {
            fprintf(stderr, "[SKIP] second-index preservation: staging unavailable on this host\n");
        }
        /* Remove both test-owned entries only if no user entry has appeared. */
        drop_fixture_shm(staged_prefix, fixture_registry, mmap, strrchr(second_prefix, '/') + 1);
    } else {
        drop_fixture_shm(fixture.first.c_str());
        fprintf(stderr, "[SKIP] load_index's shm-attach path: could not stage the unique fixture on this host\n");
    }

    delete fmi;
    fprintf(stderr, "%d / %d cases passed\n", passed_cases, total_cases);
    return passed_cases == total_cases ? 0 : 1;
}
