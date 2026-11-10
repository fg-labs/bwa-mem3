/* bwa_shm_stage holds the stager lock only for registry bookkeeping, marks a
 * segment complete only after the pack, and recovers from an interrupted or
 * failed stage. Every check stages a small index under private per-run
 * registry, lock and segment names, so the host registry is never touched.
 *
 * The `shm --huge-pages` checks (h1..h13) run the redirect paths on any OS:
 * an ordinary directory stands in for a hugetlbfs mount with 1 GB pages (the
 * hugetlbfs check accepts only it), the page prefault is a no-op, and the
 * minimum index size is 0, so the "1 GB page file" is a sparse 1 GiB file.
 * On Linux with BWA_TEST_HUGETLB_DIR set to a hugetlbfs mount with
 * pagesize=1G and at least one free 1 GB page, a further check stages on real
 * 1 GB pages; otherwise it prints SKIP.
 *
 * Usage: shm_stage_test <indexed-prefix>   (e.g. test/fixtures/phix.fa)
 */

// Distinct stage/destroy symbols, so an accidental archive fallback cannot
// satisfy those calls with the global-name implementation. The other symbols
// defined here keep libbwa's bwa_shm.o out of the link, as in
// shm_lock_destroy_test.
#define bwa_shm_stage   test_private_shm_stage
#define bwa_shm_destroy test_private_shm_destroy
#include "bwa_shm.h"
#include "bwa_hugepages.h"
#include "framework/temp_dir.h"

#include <algorithm>
#include <dirent.h>
#include <semaphore.h>
#include <signal.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/statvfs.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <unistd.h>
#include <cerrno>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <stdexcept>
#include <string>
#include <vector>

static char private_control[32], private_lock[32], private_index[32];
// The wait deadline and the /dev/shm probe are test-controlled: a child sets
// g_wait_s before staging, and g_fake_free_bytes >= 0 makes the probe report
// that many free bytes (otherwise it fails, so the preflight is skipped).
static int g_wait_s = 30;
static long long g_fake_free_bytes = -1;
static int fake_statvfs(const char *, struct statvfs *v) {
    if (g_fake_free_bytes < 0) return -1;
    memset(v, 0, sizeof *v);
    v->f_frsize = 1;
    v->f_bavail = (fsblkcnt_t)g_fake_free_bytes;
    return 0;
}
// The 1 GB page seams. The hugetlbfs check passes only for the directory
// g_huge_dev/g_huge_ino; with g_fs_ok_budget >= 0, only that many more checks
// can pass. g_prefault_fail makes the prefault fail as the real one does
// (unmapping the buffer). g_real_hugetlb restores both real implementations.
static bool g_huge_set = false;
static dev_t g_huge_dev;
static ino_t g_huge_ino;
static int g_fs_ok_budget = -1;
static bool g_prefault_fail = false;
static bool g_real_hugetlb = false;
static bool hugetlbfs_1g_dir(int fd);   // the real check, in bwa_shm.cpp
static bool fake_fs_ok(int fd) {
    if (g_real_hugetlb) return hugetlbfs_1g_dir(fd);
    if (g_fs_ok_budget == 0) return false;
    if (g_fs_ok_budget > 0) --g_fs_ok_budget;
    struct stat st;
    return g_huge_set && fstat(fd, &st) == 0 && S_ISDIR(st.st_mode) &&
           st.st_dev == g_huge_dev && st.st_ino == g_huge_ino;
}
static int fake_prefault(bwamem_hugetlb_buf_t *buf, int nthreads) {
    if (g_real_hugetlb) return bwamem_hugetlb_prefault(buf, nthreads);
    if (!g_prefault_fail) return 0;
    munmap(buf->base, buf->map_bytes);
    buf->base = NULL;
    buf->map_bytes = 0;
    return 1;
}
#define BWA_SHM_STAGE_WAIT_S g_wait_s
#define BWA_SHM_STATVFS fake_statvfs
#define BWA_SHM_HUGETLB_FS_OK(dirfd) fake_fs_ok(dirfd)
#define BWA_SHM_HUGETLB_PREFAULT(buf, nthreads) fake_prefault(buf, nthreads)
#define BWA_SHM_HUGETLB_MIN_BYTES 0
#undef BWA_SHM_CTL_NAME
#undef BWA_SHM_LOCK_NAME
#undef BWA_SHM_IDX_PREFIX
#define BWA_SHM_CTL_NAME private_control
#define BWA_SHM_LOCK_NAME private_lock
#define BWA_SHM_IDX_PREFIX private_index
#include "../src/bwa_shm.cpp"
#undef BWA_SHM_CTL_NAME
#undef BWA_SHM_LOCK_NAME
#undef BWA_SHM_IDX_PREFIX
#undef bwa_shm_stage
#undef bwa_shm_destroy

#define CHECK(cond) do { \
    if (!(cond)) throw std::runtime_error(std::string("CHECK failed: ") + #cond + \
                                          " (line " + std::to_string(__LINE__) + ")"); \
} while (0)

namespace {

const char *g_source = NULL;   // indexed prefix from argv[1]

// One case's private namespace: a fresh token, with the registry and lock
// reserved by O_EXCL. Destruction unlinks the token's registry, lock and the
// segments for `bases` (all names only this run can have created).
// Remove the files in directory `path` (not recursive), then the directory.
void remove_dir(const std::string &path) {
    DIR *d = opendir(path.c_str());
    if (d == NULL) return;
    while (struct dirent *e = readdir(d)) {
        if (std::strcmp(e->d_name, ".") != 0 && std::strcmp(e->d_name, ".."))
            unlink((path + "/" + e->d_name).c_str());
    }
    closedir(d);
    rmdir(path.c_str());
}

struct Names {
    bwa_tests::TempDir dir{"ss"};
    std::string segment;                 // the first base's segment
    std::vector<std::string> segments;
    std::string huge;                    // the stand-in hugetlbfs mount (resolved)
    std::vector<std::string> subdirs;    // removed with their files
    explicit Names(std::initializer_list<const char *> bases) {
        const std::string token = dir.path().substr(dir.path().rfind('.') + 1);
        CHECK(token.size() == 6);
        snprintf(private_control, sizeof private_control, "/ss%sc", token.c_str());
        snprintf(private_lock, sizeof private_lock, "/ss%sl", token.c_str());
        snprintf(private_index, sizeof private_index, "/ss%s-", token.c_str());
        for (const char *base : bases) segments.push_back(std::string(private_index) + base);
        segment = segments.at(0);
        const int fd = shm_open(private_control, O_CREAT | O_EXCL | O_RDWR, 0600);
        CHECK(fd >= 0);
        close(fd);
        sem_t *s = sem_open(private_lock, O_CREAT | O_EXCL, 0600, 1);
        CHECK(s != SEM_FAILED);
        sem_close(s);
        for (const std::string &seg : segments)
            CHECK(shm_open(seg.c_str(), O_RDONLY, 0) < 0 && errno == ENOENT);
        huge = subdir("huge");
        struct stat st;
        CHECK(stat(huge.c_str(), &st) == 0);
        g_huge_dev = st.st_dev;
        g_huge_ino = st.st_ino;
        g_huge_set = true;
    }
    ~Names() {
        for (const std::string &sub : subdirs) remove_dir(sub);
        g_huge_set = false;
        for (const std::string &seg : segments) shm_unlink(seg.c_str());
        shm_unlink(private_control);
        sem_unlink(private_lock);
    }
    // A new subdirectory of this case's directory, by its resolved path (the
    // path a stage records, since it resolves --huge-pages=DIR).
    std::string subdir(const char *base) {
        CHECK(mkdir(dir.file(base).c_str(), 0700) == 0);
        char resolved[PATH_MAX];
        CHECK(realpath(dir.file(base).c_str(), resolved) != NULL);
        subdirs.push_back(resolved);
        return resolved;
    }
    // The 1 GB page file a --huge-pages stage of `base` creates.
    std::string huge_file(const char *base) const {
        return huge + "/" BWA_SHM_HUGETLB_FILE_PREFIX + base;
    }
    // An index prefix named `base` in this case's directory: symlinks to the
    // source index files, or copies when the case alters them.
    std::string prefix(const char *base, bool copy = false) const {
        const std::string p = dir.file(base);
        for (const char *suffix : {".ann", ".amb", ".pac", ".bwt.2bit.64"}) {
            const std::string src = std::string(g_source) + suffix, dst = p + suffix;
            if (!copy) { CHECK(symlink(src.c_str(), dst.c_str()) == 0); continue; }
            FILE *in = fopen(src.c_str(), "rb"), *out = fopen(dst.c_str(), "wb");
            CHECK(in != NULL && out != NULL);
            char buf[1 << 16];
            size_t n;
            while ((n = fread(buf, 1, sizeof buf, in)) > 0) CHECK(fwrite(buf, 1, n, out) == n);
            fclose(in);
            CHECK(fclose(out) == 0);
        }
        return p;
    }
};

bool segment_exists(const std::string &name) {
    const int fd = shm_open(name.c_str(), O_RDONLY, 0);
    if (fd < 0) { CHECK(errno == ENOENT); return false; }
    close(fd);
    return true;
}

// Registry entry names, via the implementation's own parser.
std::vector<std::string> entries() {
    std::vector<std::string> out;
    CHECK(ctl_walk([](int64_t, const char *name, void *ctx) {
        static_cast<std::vector<std::string> *>(ctx)->push_back(name);
        return 0;
    }, &out) == 0);
    return out;
}

// The stager lock is free: a non-blocking take succeeds (and is undone).
bool lock_free() {
    sem_t *s = sem_open(private_lock, 0);
    if (s == SEM_FAILED) return errno == ENOENT;
    const bool got = sem_trywait(s) == 0;
    if (got) sem_post(s);
    sem_close(s);
    return got;
}

// bwa_shm_list's output, captured through a pipe.
std::string list_output() {
    fflush(stdout);
    int p[2];
    CHECK(pipe(p) == 0);
    const int saved = dup(STDOUT_FILENO);
    dup2(p[1], STDOUT_FILENO);
    close(p[1]);
    const int rc = bwa_shm_list();
    fflush(stdout);
    dup2(saved, STDOUT_FILENO);
    close(saved);
    std::string out;
    char buf[512];
    ssize_t n;
    while ((n = read(p[0], buf, sizeof buf)) > 0) out.append(buf, (size_t)n);
    close(p[0]);
    CHECK(rc == 0);
    return out;
}

const uint64_t GIB = 1ull << 30;

// Attach succeeds and its bytes equal bwa_shm_pack_from_disk's. With `huge`,
// the attach must be of a --huge-pages stage: whole 1 GiB pages are mapped.
void check_attach_matches(const std::string &prefix, bool huge = false) {
    size_t len = 0, map_len = 0;
    uint8_t *m = bwa_shm_attach(prefix.c_str(), &len, &map_len);
    CHECK(m != NULL);
    CHECK(map_len >= len);
    CHECK(huge == (map_len % GIB == 0));
    uint8_t *disk = NULL;
    uint64_t disk_len = 0;
    CHECK(bwa_shm_pack_from_disk(prefix.c_str(), &disk, &disk_len) == 0);
    CHECK(len == disk_len);
    // The staged header also carries the stager PID and token in reserved[0..1];
    // the rest of reserved stays zero.
    const size_t r0 = offsetof(bwa_shm_header_t, reserved);
    CHECK(std::memcmp(m, disk, r0) == 0);
    const bwa_shm_header_t *h = (const bwa_shm_header_t *)m;
    CHECK(h->reserved[0] != 0 && (h->reserved[1] & 1) == 1);
    CHECK(h->reserved[2] == 0 && h->reserved[3] == 0 && h->reserved[4] == 0);
    // BNS_STRUCT and BNS_ANNS are memcpys of the packer's structs, so their
    // pointer fields are that process's heap addresses. Blank both sections in
    // copies, compare every other byte, then the scalars of the two structs.
    std::vector<uint8_t> a(m, m + len), b(disk, disk + len);
    uint64_t soff = 0, ssize = 0, aoff = 0, asize = 0;
    CHECK(bwa_shm_section_find(disk, BWA_SHM_SEC_BNS_STRUCT, &soff, &ssize) == 0);
    CHECK(bwa_shm_section_find(disk, BWA_SHM_SEC_BNS_ANNS, &aoff, &asize) == 0);
    CHECK(ssize == sizeof(bntseq_t) && asize % sizeof(bntann1_t) == 0);
    // bwa_shm_pack_from_disk packs into malloc'd memory, so the alignment
    // padding between sections is whatever the heap held; compare only the
    // header, the section table and the section payloads.
    std::vector<bool> covered(len, false);
    const bwa_shm_header_t *dh = (const bwa_shm_header_t *)disk;
    const size_t table_end = sizeof(bwa_shm_header_t) + dh->n_sections * sizeof(bwa_shm_section_t);
    for (size_t i = 0; i < table_end && i < len; ++i) covered[i] = true;
    for (uint32_t k = 0; k < dh->n_sections; ++k) {
        const bwa_shm_section_t *t = (const bwa_shm_section_t *)(disk + sizeof(bwa_shm_header_t)) + k;
        for (uint64_t i = t->offset; i < t->offset + t->size && i < len; ++i) covered[i] = true;
    }
    for (uint8_t *v : {a.data(), b.data()}) {
        for (size_t i = 0; i < len; ++i) if (!covered[i]) v[i] = 0;
        std::memset(v + r0, 0, 16);
        std::memset(v + soff, 0, ssize);
        std::memset(v + aoff, 0, asize);
    }
    CHECK(a == b);
    bntseq_t sb, db;
    std::memcpy(&sb, m + soff, sizeof sb);
    std::memcpy(&db, disk + soff, sizeof db);
    CHECK(sb.l_pac == db.l_pac && sb.n_seqs == db.n_seqs && sb.n_holes == db.n_holes
          && sb.seed == db.seed);
    for (uint64_t k = 0; k < asize / sizeof(bntann1_t); ++k) {
        bntann1_t sa, da;
        std::memcpy(&sa, m + aoff + k * sizeof sa, sizeof sa);
        std::memcpy(&da, disk + aoff + k * sizeof da, sizeof da);
        CHECK(sa.offset == da.offset && sa.len == da.len && sa.n_ambs == da.n_ambs
              && sa.gi == da.gi && sa.is_alt == da.is_alt);
    }
    free(disk);
    munmap(m, map_len);
}

// Create segment `name` by hand, `size` bytes, with a header of `magic`, PID
// `pid`, token `token` and total_size `declared` (0: `size`) when size allows.
// With `fill`, every byte is written, so tmpfs allocates every page.
void make_segment(const std::string &name, size_t size, uint64_t magic, uint64_t pid,
                  uint64_t token, uint64_t declared = 0, bool fill = false) {
    const int fd = shm_open(name.c_str(), O_CREAT | O_EXCL | O_RDWR, 0600);
    CHECK(fd >= 0);
    if (size > 0) {
        CHECK(ftruncate(fd, (off_t)size) == 0);
        void *m = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        CHECK(m != MAP_FAILED);
        if (fill) memset(m, 0xa5, size);
        bwa_shm_header_t hdr;
        memset(&hdr, 0, sizeof hdr);
        hdr.magic = magic;
        hdr.version = BWA_SHM_VERSION;
        hdr.total_size = declared != 0 ? declared : size;
        hdr.reserved[0] = pid;
        hdr.reserved[1] = token;
        memcpy(m, &hdr, sizeof hdr);
        munmap(m, size);
    }
    close(fd);
}

// The header of segment `name` (CHECKs it exists and is header-sized).
bwa_shm_header_t read_header(const std::string &name) {
    const int fd = shm_open(name.c_str(), O_RDONLY, 0);
    CHECK(fd >= 0);
    void *m = mmap(NULL, sizeof(bwa_shm_header_t), PROT_READ, MAP_SHARED, fd, 0);
    close(fd);
    CHECK(m != MAP_FAILED);
    bwa_shm_header_t hdr;
    memcpy(&hdr, m, sizeof hdr);
    munmap(m, sizeof hdr);
    return hdr;
}

// Append a registry entry by hand, through the implementation's own writer
// path (under the private lock).
void register_entry(const char *name, int64_t l_mem) {
    sem_t *lock = ctl_lock_acquire();
    CHECK(lock != NULL);
    int fd = -1;
    void *m = ctl_open_rw(&fd);
    CHECK(m != NULL);
    ctl_header h = ctl_header_load(m);
    memcpy((uint8_t *)m + h.next_write, &l_mem, sizeof l_mem);
    memcpy((uint8_t *)m + h.next_write + sizeof l_mem, name, strlen(name) + 1);
    h.n_entries++;
    h.next_write = (uint16_t)(h.next_write + sizeof l_mem + strlen(name) + 1);
    ctl_header_store(m, h);
    ctl_close(m, fd);
    ctl_lock_release(lock);
}

// Per-child hook behaviour; the child's map hook dispatches on it.
enum HookMode { HOOK_NONE, HOOK_PAUSE, HOOK_SIGKILL, HOOK_FIFO_INPUT, HOOK_FATAL_INPUT,
                HOOK_DESTROY };
HookMode g_mode = HOOK_NONE;
int g_ready_fd = -1, g_go_fd = -1;
// Parent-side pipe ends of live paused stagers; every stage child closes them,
// so only the parent can release (EOF) a paused stager.
std::vector<int> g_child_close;
std::string g_break_path;

void stage_hook(const void *, size_t) {
    switch (g_mode) {
    case HOOK_PAUSE: {
        char c = 'r';
        if (write(g_ready_fd, &c, 1) != 1) _exit(90);
        if (read(g_go_fd, &c, 1) < 0) _exit(91);   // EOF releases it
        break;
    }
    case HOOK_SIGKILL: kill(getpid(), SIGKILL); break;
    case HOOK_FIFO_INPUT:
        // A FIFO in place of the FM-index file: the pack opens it (this child
        // holds it open for writing, so the open does not block) and its
        // fseek fails with ESPIPE, the pack's returning error path.
        if (unlink(g_break_path.c_str()) != 0 || mkfifo(g_break_path.c_str(), 0600) != 0
            || open(g_break_path.c_str(), O_RDWR) < 0) _exit(92);
        break;
    case HOOK_FATAL_INPUT:
        // A missing .pac: the pack's open of it exits the process midway.
        if (unlink(g_break_path.c_str()) != 0) _exit(92);
        break;
    case HOOK_DESTROY: if (test_private_shm_destroy() != 0) _exit(93); break;
    case HOOK_NONE: break;
    }
}

// Fork a child that stages `prefix` with hook `mode` under alarm(20). Its
// exit status is 0 for a stage returning 0, 1 for -1.
pid_t spawn_stage(const std::string &prefix, HookMode mode, const char *huge = NULL) {
    fflush(NULL);
    const pid_t pid = fork();
    CHECK(pid >= 0);
    if (pid == 0) {
        for (int fd : g_child_close) close(fd);
        alarm(20);
        g_mode = mode;
        bwa_shm_stage_map_hook = stage_hook;
        const int rc = test_private_shm_stage(prefix.c_str(), false, -1, 1, huge);
        fflush(NULL);
        _exit(rc == 0 ? 0 : 1);
    }
    return pid;
}

int wait_status(pid_t pid) {
    int status = 0;
    pid_t r;
    do { r = waitpid(pid, &status, 0); } while (r < 0 && errno == EINTR);
    CHECK(r == pid);
    return status;
}

bool exited(int status, int code) { return WIFEXITED(status) && WEXITSTATUS(status) == code; }

// Seconds on the monotonic clock.
double mono_now() {
    struct timespec ts;
    CHECK(clock_gettime(CLOCK_MONOTONIC, &ts) == 0);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

// A child paused inside its stage's map hook. Destruction releases and reaps
// it, so a failed CHECK cannot leave it blocked.
struct PausedStager {
    pid_t pid = -1;
    int ready[2] = {-1, -1}, go[2] = {-1, -1};
    explicit PausedStager(const std::string &prefix, const char *huge = NULL) {
        CHECK(pipe(ready) == 0 && pipe(go) == 0);
        g_ready_fd = ready[1];
        g_go_fd = go[0];
        g_child_close.push_back(ready[0]);
        g_child_close.push_back(go[1]);
        pid = spawn_stage(prefix, HOOK_PAUSE, huge);
        close(ready[1]); ready[1] = -1;
        close(go[0]); go[0] = -1;
        char c = 0;
        CHECK(read(ready[0], &c, 1) == 1 && c == 'r');   // child is mid-stage
    }
    void forget(int &fd) {
        g_child_close.erase(std::remove(g_child_close.begin(), g_child_close.end(), fd),
                            g_child_close.end());
        close(fd);
        fd = -1;
    }
    int finish() {
        forget(go[1]);   // EOF releases the hook
        const int status = wait_status(pid);
        pid = -1;
        return status;
    }
    ~PausedStager() {
        if (go[1] >= 0) forget(go[1]);
        if (ready[0] >= 0) forget(ready[0]);
        if (pid > 0) { kill(pid, SIGKILL); waitpid(pid, NULL, 0); }
    }
};

// The child `pid` has not exited yet.
bool running(pid_t pid) { return waitpid(pid, NULL, WNOHANG) == 0; }

// 1. The lock is free while the pack runs; a second stage of the same index
//    waits for it (holding no lock) and then succeeds.
void check_lock_not_held_during_pack() {
    Names n({"live"});
    const std::string prefix = n.prefix("live");
    PausedStager a(prefix);
    CHECK(lock_free());
    CHECK(bwa_shm_attach(prefix.c_str(), NULL) == NULL);
    const std::string listed = list_output();
    CHECK(listed.rfind("live\t", 0) == 0 && listed.find("\tincomplete\n") != std::string::npos);
    const pid_t b = spawn_stage(prefix, HOOK_NONE);
    usleep(1500000);                 // B is in its wait loop by now
    CHECK(running(b));
    CHECK(lock_free());              // and holds no lock while waiting
    CHECK(exited(a.finish(), 0));
    CHECK(exited(wait_status(b), 0));
    CHECK(list_output().find("incomplete") == std::string::npos);
    check_attach_matches(prefix);
    CHECK(entries() == std::vector<std::string>{"live"});
    std::printf("ok 1 lock not held during pack; a concurrent stage waits, then succeeds\n");
}

// 1b. B waits behind A; A is killed; B reclaims and stages.
void check_wait_then_reclaim() {
    Names n({"waitkill"});
    const std::string prefix = n.prefix("waitkill");
    PausedStager a(prefix);
    const pid_t b = spawn_stage(prefix, HOOK_NONE);
    usleep(1500000);
    CHECK(running(b));
    CHECK(kill(a.pid, SIGKILL) == 0);
    const int a_status = wait_status(a.pid);
    a.pid = -1;
    CHECK(WIFSIGNALED(a_status) && WTERMSIG(a_status) == SIGKILL);
    CHECK(exited(wait_status(b), 0));
    check_attach_matches(prefix);
    CHECK(entries() == std::vector<std::string>{"waitkill"});
    std::printf("ok 1b a waiting stage reclaims after the stager is killed\n");
}

// 2. A stager killed mid-pack leaves an incomplete segment the next stage reclaims.
void check_interrupted_stage() {
    Names n({"killed"});
    const std::string prefix = n.prefix("killed");
    const int status = wait_status(spawn_stage(prefix, HOOK_SIGKILL));
    CHECK(WIFSIGNALED(status) && WTERMSIG(status) == SIGKILL);
    CHECK(lock_free());
    CHECK(bwa_shm_attach(prefix.c_str(), NULL) == NULL);
    CHECK(list_output().find("killed\t") != std::string::npos);
    CHECK(list_output().find("\tincomplete\n") != std::string::npos);
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE)), 0));   // reclaimed
    check_attach_matches(prefix);
    CHECK(entries() == std::vector<std::string>{"killed"});
    std::printf("ok 2 interrupted stage reclaimed\n");
}

// 3. A pack that returns failure removes the entry and the segment. A pack that
//    exits the process instead (err_fatal on a missing file) cannot clean up,
//    so it leaves an incomplete segment that the next stage reclaims.
void check_pack_failure() {
    {
        Names n({"fifo"});
        const std::string prefix = n.prefix("fifo", /*copy=*/true);
        g_break_path = prefix + ".bwt.2bit.64";
        CHECK(exited(wait_status(spawn_stage(prefix, HOOK_FIFO_INPUT)), 1));
        CHECK(entries().empty());
        CHECK(!segment_exists(n.segment));
        CHECK(lock_free());
    }
    {
        Names n({"fatal"});
        const std::string prefix = n.prefix("fatal", /*copy=*/true);
        g_break_path = prefix + ".pac";
        CHECK(exited(wait_status(spawn_stage(prefix, HOOK_FATAL_INPUT)), 1));
        CHECK(lock_free());
        CHECK(entries() == std::vector<std::string>{"fatal"});
        CHECK(list_output().find("\tincomplete\n") != std::string::npos);
        CHECK(bwa_shm_attach(prefix.c_str(), NULL) == NULL);
        CHECK(symlink((std::string(g_source) + ".pac").c_str(), g_break_path.c_str()) == 0);
        CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE)), 0));   // reclaimed
        check_attach_matches(prefix);
        CHECK(entries() == std::vector<std::string>{"fatal"});
    }
    std::printf("ok 3 pack failure cleans up (returned) or is reclaimable (fatal exit)\n");
}

// 4. `shm -d` during a stage: the stage fails and leaves no segment or entry.
void check_destroy_during_stage() {
    Names n({"dropped"});
    const std::string prefix = n.prefix("dropped");
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_DESTROY)), 1));
    CHECK(!segment_exists(n.segment));
    CHECK(entries().empty());
    // Step 4 reopens without O_CREAT, so the removed registry and lock stay gone.
    CHECK(shm_open(private_control, O_RDONLY, 0) < 0 && errno == ENOENT);
    CHECK(sem_open(private_lock, 0) == SEM_FAILED && errno == ENOENT);
    std::printf("ok 4 destroy during stage\n");
}

// 5. Restaging a complete index is a no-op.
void check_restage_noop() {
    Names n({"again"});
    const std::string prefix = n.prefix("again");
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE)), 0));
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE)), 0));
    CHECK(entries() == std::vector<std::string>{"again"});
    check_attach_matches(prefix);
    std::printf("ok 5 restage of a complete index is a no-op\n");
}

// 6. A live stager that never finishes: the waiting stage gives up at the
//    deadline (2 s here) and leaves the segment alone.
void check_wait_deadline() {
    Names n({"stuck"});
    const std::string prefix = n.prefix("stuck");
    make_segment(n.segment, sizeof(bwa_shm_header_t), 0, (uint64_t)getpid(), 7);
    register_entry("stuck", 64);
    const double t0 = mono_now();
    g_wait_s = 2;
    const int status = wait_status(spawn_stage(prefix, HOOK_NONE));
    g_wait_s = 30;
    const double elapsed = mono_now() - t0;
    CHECK(exited(status, 1));
    // The deadline is whole seconds of CLOCK_MONOTONIC from the stage's start,
    // so a 2 s wait lasts more than 1 s; an immediate failure would not.
    CHECK(elapsed > 1.0);
    CHECK(elapsed <= 4.0);
    const bwa_shm_header_t h = read_header(n.segment);
    CHECK(h.magic == 0 && h.reserved[0] == (uint64_t)getpid() && h.reserved[1] == 7);
    CHECK(entries() == std::vector<std::string>{"stuck"});
    std::printf("ok 6 waiting stage gives up at the deadline\n");
}

// 7. The name was taken over mid-stage (another stager reclaimed it): the
//    stage fails and leaves the replacement and its entry alone.
void check_not_ours() {
    Names n({"taken"});
    const std::string prefix = n.prefix("taken");
    PausedStager a(prefix);
    CHECK(shm_unlink(n.segment.c_str()) == 0);
    make_segment(n.segment, sizeof(bwa_shm_header_t), 0, (uint64_t)getpid(), 2);
    CHECK(exited(a.finish(), 1));
    const bwa_shm_header_t h = read_header(n.segment);
    CHECK(h.magic == 0 && h.reserved[1] == 2);
    CHECK(entries() == std::vector<std::string>{"taken"});
    std::printf("ok 7 a taken-over name is left to its new owner\n");
}

// 8. A segment with a foreign magic, or a complete header declaring more bytes
//    than the segment holds, is listed `unusable`, refused and untouched; a
//    registered 0-byte segment (never written) is reclaimed.
void check_bad_and_undersized() {
    Names n({"foreign", "empty", "short"});
    make_segment(n.segments[0], sizeof(bwa_shm_header_t), 0x1234, 0, 0);
    CHECK(exited(wait_status(spawn_stage(n.prefix("foreign"), HOOK_NONE)), 1));
    CHECK(read_header(n.segments[0]).magic == 0x1234);
    const uint64_t declared = 1u << 20;
    make_segment(n.segments[2], sizeof(bwa_shm_header_t), BWA_SHM_MAGIC, 0, 0, declared);
    register_entry("short", (int64_t)declared);
    CHECK(list_output().find("short\t1048576\tunusable\n") != std::string::npos);
    CHECK(bwa_shm_attach(n.prefix("short").c_str(), NULL) == NULL);
    CHECK(exited(wait_status(spawn_stage(n.dir.file("short"), HOOK_NONE)), 1));
    CHECK(read_header(n.segments[2]).total_size == declared);
    CHECK(entries() == std::vector<std::string>{"short"});
    make_segment(n.segments[1], 0, 0, 0, 0);
    register_entry("empty", 64);
    const std::string prefix = n.prefix("empty");
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE)), 0));
    check_attach_matches(prefix);
    CHECK(entries() == (std::vector<std::string>{"short", "empty"}));
    std::printf("ok 8 foreign and truncated segments refused; undersized segment reclaimed\n");
}

// 9. Reclaiming one of two entries leaves the other intact; removing the first
//    of three entries keeps the rest in order with the right write cursor.
void check_two_entries() {
    {
        Names n({"keep", "redo"});
        const std::string keep = n.prefix("keep"), redo = n.prefix("redo");
        CHECK(exited(wait_status(spawn_stage(keep, HOOK_NONE)), 0));
        const int status = wait_status(spawn_stage(redo, HOOK_SIGKILL));
        CHECK(WIFSIGNALED(status));
        CHECK(exited(wait_status(spawn_stage(redo, HOOK_NONE)), 0));
        CHECK(entries() == (std::vector<std::string>{"keep", "redo"}));
        check_attach_matches(keep);
        check_attach_matches(redo);
    }
    {
        Names m({"x"});   // a fresh registry for the removal check
        register_entry("aa", 1);
        register_entry("bbb", 2);
        register_entry("c", 3);
        int fd = -1;
        void *map = ctl_open_rw(&fd);
        CHECK(map != NULL);
        ctl_remove(map, "aa");
        const ctl_header h = ctl_header_load(map);
        ctl_close(map, fd);
        CHECK(h.n_entries == 2);
        CHECK(h.next_write == BWA_SHM_CTL_HEADER_BYTES + 2 * sizeof(int64_t) + 4 + 2);
        CHECK(entries() == (std::vector<std::string>{"bbb", "c"}));
    }
    std::printf("ok 9 multi-entry registry survives reclaim and removal\n");
}

// 10. Free-space preflight counts other in-progress stages as used.
void check_capacity() {
    Names n({"cap", "other"});
    const std::string prefix = n.prefix("cap");
    uint8_t *disk = NULL;
    uint64_t need = 0;
    CHECK(bwa_shm_pack_from_disk(prefix.c_str(), &disk, &need) == 0);
    free(disk);
    // Another live stager has reserved `need` bytes; free space covers one index.
    make_segment(n.segments[1], sizeof(bwa_shm_header_t), 0, (uint64_t)getpid(), 9);
    register_entry("other", (int64_t)need);
    g_fake_free_bytes = (long long)need + (long long)need / 2;
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE)), 1));
    CHECK(!segment_exists(n.segments[0]));
    g_fake_free_bytes = 2 * (long long)need;
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE)), 0));
    g_fake_free_bytes = -1;
    std::printf("ok 10 preflight subtracts in-progress stages\n");
}

// 10b. Pages an in-progress stage has already written are in use, so the
//      preflight reserves only its unwritten remainder. Where st_blocks is not
//      reported (macOS) the whole stage stays reserved.
void check_capacity_allocated() {
    Names n({"cap2", "full"});
    const std::string prefix = n.prefix("cap2");
    uint8_t *disk = NULL;
    uint64_t need = 0;
    CHECK(bwa_shm_pack_from_disk(prefix.c_str(), &disk, &need) == 0);
    free(disk);
    make_segment(n.segments[1], (size_t)need, 0, (uint64_t)getpid(), 9, 0, /*fill=*/true);
    register_entry("full", (int64_t)need);
    const int fd = shm_open(n.segments[1].c_str(), O_RDONLY, 0);
    CHECK(fd >= 0);
    struct stat st;
    CHECK(fstat(fd, &st) == 0);
    close(fd);
    const bool reported = st.st_blocks > 0;
    if (reported) CHECK((uint64_t)st.st_blocks * 512u >= need);
    const long long reserved = reported ? 0 : (long long)need;
    g_fake_free_bytes = (long long)need + reserved - 1;
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE)), 1));
    CHECK(!segment_exists(n.segments[0]));
    g_fake_free_bytes = (long long)need + reserved;
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE)), 0));
    g_fake_free_bytes = -1;
    std::printf("ok 10b preflight reserves only the unwritten part of a stage%s\n",
                reported ? "" : " (st_blocks not reported: whole stage reserved)");
}

// 11. Restaging a staged index reads no disk files: with them deleted it still
//     returns 0 (the unlocked early exit).
void check_early_exit() {
    Names n({"gone"});
    const std::string prefix = n.prefix("gone");
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE)), 0));
    for (const char *suffix : {".ann", ".amb", ".pac", ".bwt.2bit.64"})
        CHECK(unlink((prefix + suffix).c_str()) == 0);
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE)), 0));
    std::printf("ok 11 restage of a staged index reads no disk files\n");
}

// ---- shm --huge-pages ---------------------------------------------------

bool file_exists(const std::string &path) {
    struct stat st;
    return lstat(path.c_str(), &st) == 0;
}

// The bytes of file `path` ("" if it cannot be opened).
std::string file_text(const std::string &path) {
    std::string text;
    FILE *fp = fopen(path.c_str(), "rb");
    if (fp == NULL) return text;
    char buf[256];
    size_t n;
    while ((n = fread(buf, 1, sizeof buf, fp)) > 0) text.append(buf, n);
    fclose(fp);
    return text;
}

void write_text(const std::string &path, const std::string &text) {
    FILE *fp = fopen(path.c_str(), "wb");
    CHECK(fp != NULL);
    CHECK(fwrite(text.data(), 1, text.size(), fp) == text.size());
    CHECK(fclose(fp) == 0);
}

off_t segment_size(const std::string &name) {
    const int fd = shm_open(name.c_str(), O_RDONLY, 0);
    CHECK(fd >= 0);
    struct stat st;
    CHECK(fstat(fd, &st) == 0);
    close(fd);
    return st.st_size;
}

bwa_shm_hugetlb_stub_t read_stub(const std::string &name) {
    // macOS rounds a shm object up to whole pages.
    CHECK(segment_size(name) >= (off_t)sizeof(bwa_shm_hugetlb_stub_t));
    const int fd = shm_open(name.c_str(), O_RDONLY, 0);
    CHECK(fd >= 0);
    void *m = mmap(NULL, sizeof(bwa_shm_hugetlb_stub_t), PROT_READ, MAP_SHARED, fd, 0);
    close(fd);
    CHECK(m != MAP_FAILED);
    bwa_shm_hugetlb_stub_t stub;
    memcpy(&stub, m, sizeof stub);
    munmap(m, sizeof stub);
    return stub;
}

// Create segment `name` holding a stub for `path` by hand.
void make_stub(const std::string &name, const std::string &path, uint64_t magic, uint64_t pid,
               uint64_t token, uint64_t total_size, uint64_t file_size) {
    bwa_shm_hugetlb_stub_t stub;
    memset(&stub, 0, sizeof stub);
    stub.magic = magic;
    stub.version = BWA_SHM_HUGETLB_STUB_VERSION;
    stub.kind = BWA_SHM_HUGETLB_KIND_1G;
    stub.total_size = total_size;
    stub.stager_pid = pid;
    stub.stage_token = token;
    stub.file_size = file_size;
    CHECK(path.size() < sizeof stub.path);
    memcpy(stub.path, path.c_str(), path.size() + 1);
    const int fd = shm_open(name.c_str(), O_CREAT | O_EXCL | O_RDWR, 0600);
    CHECK(fd >= 0);
    CHECK(ftruncate(fd, sizeof stub) == 0);
    void *m = mmap(NULL, sizeof stub, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    close(fd);
    CHECK(m != MAP_FAILED);
    memcpy(m, &stub, sizeof stub);
    munmap(m, sizeof stub);
}

// Write a 1 GiB sparse file at `path` holding `prefix`'s packed index with
// stage token `token` (and magic `magic`) in its header; returns total_size.
uint64_t make_index_file(const std::string &path, const std::string &prefix, uint64_t token,
                         uint64_t magic = BWA_SHM_MAGIC) {
    uint8_t *packed = NULL;
    uint64_t len = 0;
    CHECK(bwa_shm_pack_from_disk(prefix.c_str(), &packed, &len) == 0);
    bwa_shm_header_t *h = (bwa_shm_header_t *)packed;
    h->reserved[1] = token;
    h->magic = magic;
    const int fd = open(path.c_str(), O_CREAT | O_EXCL | O_RDWR, 0644);
    CHECK(fd >= 0);
    CHECK(ftruncate(fd, (off_t)GIB) == 0);
    CHECK(pwrite(fd, packed, len, 0) == (ssize_t)len);
    close(fd);
    free(packed);
    return len;
}

// A PID that is no longer running: a child that has exited and been reaped.
uint64_t dead_pid() {
    const pid_t pid = fork();
    CHECK(pid >= 0);
    if (pid == 0) _exit(0);
    CHECK(exited(wait_status(pid), 0));
    return (uint64_t)pid;
}

// h1. A --huge-pages stage puts the index in a whole-page file under the mount
//     and a complete stub in the segment; attach maps the file, and `shm -d`
//     removes the file and then the stub.
void check_huge_round_trip() {
    Names n({"h1"});
    const std::string prefix = n.prefix("h1"), file = n.huge_file("h1");
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE, n.huge.c_str())), 0));
    struct stat st;
    CHECK(lstat(file.c_str(), &st) == 0 && S_ISREG(st.st_mode));
    CHECK(st.st_size > 0 && (uint64_t)st.st_size % GIB == 0);
    const bwa_shm_hugetlb_stub_t stub = read_stub(n.segment);
    CHECK(stub.magic == BWA_SHM_HUGETLB_STUB_MAGIC && std::string(stub.path) == file);
    CHECK(stub.file_size == (uint64_t)st.st_size && stub.stager_pid != 0);
    check_attach_matches(prefix, /*huge=*/true);
    const std::string listed = list_output();
    CHECK(listed.rfind("h1\t", 0) == 0 && listed.find("\t", 3) == std::string::npos);
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE, n.huge.c_str())), 0));   // no-op
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE)), 0));   // stays where it is
    CHECK(test_private_shm_destroy() == 0);
    CHECK(!file_exists(file));
    CHECK(!segment_exists(n.segment));
    std::printf("ok h1 --huge-pages stage, attach and destroy\n");
}

// h2. Mid-pack the lock is free, the stub is registered with magic 0 and names
//     its file, attach falls back, and a second stager waits, then succeeds.
void check_huge_in_progress() {
    Names n({"h2"});
    const std::string prefix = n.prefix("h2");
    PausedStager a(prefix, n.huge.c_str());
    CHECK(lock_free());
    const bwa_shm_hugetlb_stub_t stub = read_stub(n.segment);
    CHECK(stub.magic == 0 && stub.stager_pid == (uint64_t)a.pid && stub.file_size == GIB);
    CHECK(file_exists(n.huge_file("h2")));
    CHECK(bwa_shm_attach(prefix.c_str(), NULL) == NULL);
    CHECK(entries() == std::vector<std::string>{"h2"});
    CHECK(list_output().find("\tincomplete\n") != std::string::npos);
    const pid_t b = spawn_stage(prefix, HOOK_NONE, n.huge.c_str());
    usleep(1500000);
    CHECK(running(b));
    CHECK(lock_free());
    CHECK(exited(a.finish(), 0));
    CHECK(exited(wait_status(b), 0));
    check_attach_matches(prefix, true);
    std::printf("ok h2 a --huge-pages stage in progress is waited on\n");
}

// h3. A --huge-pages stage killed mid-pack leaves an incomplete stub naming
//     its file; a restage with or without --huge-pages, or `shm -d`, removes
//     the file.
void check_huge_interrupted() {
    for (int how = 0; how < 3; ++how) {
        Names n({"h3"});
        const std::string prefix = n.prefix("h3"), file = n.huge_file("h3");
        const int status = wait_status(spawn_stage(prefix, HOOK_SIGKILL, n.huge.c_str()));
        CHECK(WIFSIGNALED(status) && WTERMSIG(status) == SIGKILL);
        CHECK(file_exists(file));
        CHECK(list_output().find("\tincomplete\n") != std::string::npos);
        if (how == 0) {
            CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE)), 0));
            CHECK(!file_exists(file));
            CHECK(read_header(n.segment).magic == BWA_SHM_MAGIC);   // an index, not a stub
            check_attach_matches(prefix);
        } else if (how == 1) {
            const bwa_shm_hugetlb_stub_t old = read_stub(n.segment);
            CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE, n.huge.c_str())), 0));
            CHECK(read_stub(n.segment).stage_token != old.stage_token);
            check_attach_matches(prefix, true);
        } else {
            CHECK(test_private_shm_destroy() == 0);
            CHECK(!file_exists(file));
            CHECK(!segment_exists(n.segment));
        }
    }
    std::printf("ok h3 an interrupted --huge-pages stage is reclaimed with its file\n");
}

// h4. A dead stub that never recorded a file (killed between creating the
//     stub and the file) is reclaimed without touching a file at its path.
void check_huge_dead_stub_without_file() {
    Names n({"h4"});
    const std::string prefix = n.prefix("h4"), file = n.huge_file("h4");
    write_text(file, "keep");
    make_stub(n.segment, file, 0, dead_pid(), 5, 64, 0);
    register_entry("h4", 64);
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE)), 0));
    CHECK(file_text(file) == "keep");
    check_attach_matches(prefix);
    std::printf("ok h4 a dead stub with no file recorded leaves the path alone\n");
}

// h5. A file already at the path, named by no stub, is never adopted or
//     removed: the stage fails and leaves no stub or entry behind.
void check_huge_orphan_file() {
    Names n({"h5"});
    const std::string prefix = n.prefix("h5"), file = n.huge_file("h5");
    write_text(file, "keep");
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE, n.huge.c_str())), 1));
    CHECK(file_text(file) == "keep");
    CHECK(!segment_exists(n.segment));
    CHECK(entries().empty());
    CHECK(lock_free());
    std::printf("ok h5 an unnamed file at the path is an error, left untouched\n");
}

// h6. No mount (the directory is missing), or the mount check failing once
//     the file exists: the stage falls back to POSIX shared memory, and no
//     file is left.
void check_huge_reserve_fallback() {
    {
        Names n({"h6"});
        const std::string prefix = n.prefix("h6"), missing = n.huge + "/missing";
        CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE, missing.c_str())), 0));
        check_attach_matches(prefix);
    }
    {
        Names n({"h6"});
        const std::string prefix = n.prefix("h6");
        g_fs_ok_budget = 1;   // the check before the lock passes, the one after create fails
        const int status = wait_status(spawn_stage(prefix, HOOK_NONE, n.huge.c_str()));
        g_fs_ok_budget = -1;
        CHECK(exited(status, 0));
        CHECK(!file_exists(n.huge_file("h6")));
        CHECK(read_header(n.segment).magic == BWA_SHM_MAGIC);   // an index, not a stub
        check_attach_matches(prefix);
        CHECK(entries() == std::vector<std::string>{"h6"});
    }
    std::printf("ok h6 a mount that cannot be used falls back to /dev/shm\n");
}

// h7. A prefault failure (an old kernel) rolls the registered stage back and
//     stages again in POSIX shared memory.
void check_huge_prefault_fallback() {
    Names n({"h7"});
    const std::string prefix = n.prefix("h7");
    g_prefault_fail = true;
    const int status = wait_status(spawn_stage(prefix, HOOK_NONE, n.huge.c_str()));
    g_prefault_fail = false;
    CHECK(exited(status, 0));
    CHECK(!file_exists(n.huge_file("h7")));
    check_attach_matches(prefix);
    CHECK(entries() == std::vector<std::string>{"h7"});
    std::printf("ok h7 a prefault failure falls back to /dev/shm\n");
}

// h8. A --huge-pages pack that returns failure removes entry, file and stub.
void check_huge_pack_failure() {
    Names n({"h8"});
    const std::string prefix = n.prefix("h8", /*copy=*/true);
    g_break_path = prefix + ".bwt.2bit.64";
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_FIFO_INPUT, n.huge.c_str())), 1));
    CHECK(entries().empty());
    CHECK(!segment_exists(n.segment));
    CHECK(!file_exists(n.huge_file("h8")));
    CHECK(lock_free());
    std::printf("ok h8 a failed --huge-pages pack removes entry, file and stub\n");
}

// h9. `shm -d` during a --huge-pages stage removes the file and the stub; the
//     stage fails and recreates nothing.
void check_huge_destroy_during_stage() {
    Names n({"h9"});
    const std::string prefix = n.prefix("h9");
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_DESTROY, n.huge.c_str())), 1));
    CHECK(!segment_exists(n.segment));
    CHECK(!file_exists(n.huge_file("h9")));
    CHECK(shm_open(private_control, O_RDONLY, 0) < 0 && errno == ENOENT);
    CHECK(sem_open(private_lock, 0) == SEM_FAILED && errno == ENOENT);
    std::printf("ok h9 destroy during a --huge-pages stage\n");
}

// h10. The name was taken over mid-stage (a reclaimer replaced stub and file at
//      the same path): the first stager fails and leaves both alone.
void check_huge_not_ours() {
    Names n({"h10"});
    const std::string prefix = n.prefix("h10"), file = n.huge_file("h10");
    PausedStager a(prefix, n.huge.c_str());
    CHECK(unlink(file.c_str()) == 0);
    CHECK(shm_unlink(n.segment.c_str()) == 0);
    write_text(file, "new!");
    make_stub(n.segment, file, 0, (uint64_t)getpid(), 2, 64, GIB);
    CHECK(exited(a.finish(), 1));
    CHECK(file_text(file) == "new!");
    CHECK(read_stub(n.segment).stage_token == 2);
    CHECK(entries() == std::vector<std::string>{"h10"});
    std::printf("ok h10 a taken-over --huge-pages stage is left to its new owner\n");
}

// h11. A complete stub is followed only to its own index: never to a file of
//      another name, a symlink, a directory that is not the mount, a file of
//      another size, or one whose header is not the stub's. Such a stub is
//      `unusable`, attach refuses it, a stage of its name refuses to replace
//      it, and `shm -d` removes the stub but not the file. A stub that checks
//      out (s0) is followed by all three.
void check_huge_stub_safety() {
    Names n({"s0", "s1", "s2", "s3", "s4", "s5", "s6"});
    const std::string plain = n.subdir("plain");
    const std::string p0 = n.prefix("s0");
    CHECK(exited(wait_status(spawn_stage(p0, HOOK_NONE, n.huge.c_str())), 0));
    std::vector<std::string> keep;   // files that must survive `shm -d`
    // s1: the file of another name.
    std::string f = n.huge_file("other");
    uint64_t total = make_index_file(f, n.prefix("s1"), 11);
    make_stub(n.segments[1], f, BWA_SHM_HUGETLB_STUB_MAGIC, 1, 11, total, GIB);
    keep.push_back(f);
    // s2: a symlink to a good file outside the mount.
    const std::string real = n.dir.file("real-s2");
    total = make_index_file(real, n.prefix("s2"), 12);
    CHECK(symlink(real.c_str(), n.huge_file("s2").c_str()) == 0);
    make_stub(n.segments[2], n.huge_file("s2"), BWA_SHM_HUGETLB_STUB_MAGIC, 1, 12, total, GIB);
    keep.push_back(real);
    keep.push_back(n.huge_file("s2"));
    // s3: a directory the hugetlbfs check rejects.
    f = plain + "/" BWA_SHM_HUGETLB_FILE_PREFIX "s3";
    total = make_index_file(f, n.prefix("s3"), 13);
    make_stub(n.segments[3], f, BWA_SHM_HUGETLB_STUB_MAGIC, 1, 13, total, GIB);
    keep.push_back(f);
    // s4: a file of another size than the stub records.
    f = n.huge_file("s4");
    total = make_index_file(f, n.prefix("s4"), 14);
    make_stub(n.segments[4], f, BWA_SHM_HUGETLB_STUB_MAGIC, 1, 14, total, 2 * GIB);
    keep.push_back(f);
    // s5: a file with another stage's token.
    f = n.huge_file("s5");
    total = make_index_file(f, n.prefix("s5"), 7);
    make_stub(n.segments[5], f, BWA_SHM_HUGETLB_STUB_MAGIC, 1, 8, total, GIB);
    keep.push_back(f);
    // s6: a file still being filled (header magic 0).
    f = n.huge_file("s6");
    total = make_index_file(f, n.prefix("s6"), 16, /*magic=*/0);
    make_stub(n.segments[6], f, BWA_SHM_HUGETLB_STUB_MAGIC, 1, 16, total, GIB);
    keep.push_back(f);
    for (const char *s : {"s1", "s2", "s3", "s4", "s5", "s6"}) register_entry(s, (int64_t)total);

    const std::string listed = list_output();
    CHECK(listed.find("s0\t" + std::to_string(total) + "\n") != std::string::npos);
    for (const char *s : {"s1", "s2", "s3", "s4", "s5", "s6"}) {
        CHECK(listed.find(std::string(s) + "\t" + std::to_string(total) + "\tunusable\n")
              != std::string::npos);
        CHECK(bwa_shm_attach(n.dir.file(s).c_str(), NULL) == NULL);
    }
    check_attach_matches(p0, true);
    CHECK(exited(wait_status(spawn_stage(n.dir.file("s5"), HOOK_NONE)), 1));
    CHECK(exited(wait_status(spawn_stage(n.dir.file("s5"), HOOK_NONE, n.huge.c_str())), 1));
    CHECK(test_private_shm_destroy() == 0);
    for (const std::string &k : keep) CHECK(file_exists(k));
    CHECK(!file_exists(n.huge_file("s0")));
    for (const std::string &seg : n.segments) CHECK(!segment_exists(seg));
    std::printf("ok h11 a stub is followed only to its own index\n");
}

// h11b. A stub is never followed to a file owned by another user (needs root
//       to make one).
void check_huge_stub_owner() {
    if (geteuid() != 0) {
        std::printf("SKIP h11b stub owner check (needs root to chown a file)\n");
        return;
    }
    Names n({"o1"});
    const std::string f = n.huge_file("o1");
    const uint64_t total = make_index_file(f, n.prefix("o1"), 21);
    CHECK(chown(f.c_str(), 1, (gid_t)-1) == 0);
    make_stub(n.segment, f, BWA_SHM_HUGETLB_STUB_MAGIC, 1, 21, total, GIB);
    register_entry("o1", (int64_t)total);
    CHECK(list_output().find("\tunusable\n") != std::string::npos);
    CHECK(bwa_shm_attach(n.dir.file("o1").c_str(), NULL) == NULL);
    CHECK(test_private_shm_destroy() == 0);
    CHECK(file_exists(f));
    std::printf("ok h11b a stub is not followed to another user's file\n");
}

// h12. A complete stub whose file is gone is `missing`; a restage reclaims it.
void check_huge_dangling() {
    Names n({"h12"});
    const std::string prefix = n.prefix("h12");
    make_stub(n.segment, n.huge_file("h12"), BWA_SHM_HUGETLB_STUB_MAGIC, 1, 3, 64, GIB);
    register_entry("h12", 64);
    CHECK(list_output() == "h12\t64\tmissing\n");
    CHECK(bwa_shm_attach(prefix.c_str(), NULL) == NULL);
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE, n.huge.c_str())), 0));
    check_attach_matches(prefix, true);
    CHECK(entries() == std::vector<std::string>{"h12"});
    std::printf("ok h12 a stub whose file is gone is missing, then reclaimed\n");
}

// h13. A --huge-pages stage in progress takes nothing from /dev/shm, so a
//      POSIX stage that exactly fits the free space still succeeds.
void check_huge_capacity() {
    Names n({"cap3", "big"});
    const std::string prefix = n.prefix("cap3");
    uint8_t *disk = NULL;
    uint64_t need = 0;
    CHECK(bwa_shm_pack_from_disk(prefix.c_str(), &disk, &need) == 0);
    free(disk);
    PausedStager a(n.prefix("big"), n.huge.c_str());
    g_fake_free_bytes = (long long)need;
    const int status = wait_status(spawn_stage(prefix, HOOK_NONE));
    g_fake_free_bytes = -1;
    CHECK(exited(status, 0));
    CHECK(exited(a.finish(), 0));
    std::printf("ok h13 a --huge-pages stage in progress is not charged to /dev/shm\n");
}

// The hugetlb pool's free and reserved 1 GB page counts.
std::pair<long, long> hugetlb_counts() {
    auto read_long = [](const char *path) {
        FILE *fp = fopen(path, "r");
        long v = -1;
        if (fp != NULL) { if (fscanf(fp, "%ld", &v) != 1) v = -1; fclose(fp); }
        return v;
    };
    const char *dir = "/sys/kernel/mm/hugepages/hugepages-1048576kB/";
    return {read_long((std::string(dir) + "free_hugepages").c_str()),
            read_long((std::string(dir) + "resv_hugepages").c_str())};
}

// KernelPageSize in kB of the /proc/self/smaps entry holding `addr`, or -1.
long kernel_page_kb(const void *addr) {
    FILE *fp = fopen("/proc/self/smaps", "r");
    if (fp == NULL) return -1;
    char line[512];
    bool in = false;
    long kb = -1;
    const uintptr_t a = (uintptr_t)addr;
    while (fgets(line, sizeof(line), fp) != NULL) {
        unsigned long lo, hi;
        if (sscanf(line, "%lx-%lx ", &lo, &hi) == 2) { in = lo <= a && a < hi; continue; }
        if (in && sscanf(line, "KernelPageSize: %ld kB", &kb) == 1) break;
    }
    fclose(fp);
    return in ? kb : -1;
}

// hp. Linux, BWA_TEST_HUGETLB_DIR: a stage on real 1 GB pages is mapped with
//     them, and its pages return to the pool after `shm -d`, and after a
//     killed stage is reclaimed or destroyed.
void check_huge_real_pool() {
    const char *hdir = getenv("BWA_TEST_HUGETLB_DIR");
#if defined(__linux__)
    if (hdir == NULL || hdir[0] == '\0') {
        std::printf("SKIP hp real 1 GB pages (BWA_TEST_HUGETLB_DIR not set)\n");
        return;
    }
    g_real_hugetlb = true;
    const std::pair<long, long> start = hugetlb_counts();
    {
        Names n({"hpreal"});
        const std::string prefix = n.prefix("hpreal");
        CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE, hdir)), 0));
        size_t len = 0, map_len = 0;
        uint8_t *m = bwa_shm_attach(prefix.c_str(), &len, &map_len);
        CHECK(m != NULL && map_len == GIB);
        CHECK(kernel_page_kb(m) == 1048576);
        munmap(m, map_len);
        check_attach_matches(prefix, true);
        CHECK(test_private_shm_destroy() == 0);
        CHECK(hugetlb_counts() == start);
    }
    for (int how = 0; how < 2; ++how) {
        Names n({"hpreal"});
        const std::string prefix = n.prefix("hpreal");
        const int status = wait_status(spawn_stage(prefix, HOOK_SIGKILL, hdir));
        CHECK(WIFSIGNALED(status));
        if (how == 0) {
            CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE)), 0));
            CHECK(hugetlb_counts() == start);
        }
        CHECK(test_private_shm_destroy() == 0);
        CHECK(hugetlb_counts() == start);
    }
    g_real_hugetlb = false;
    std::printf("ok hp real 1 GB pages: mapped on them, returned after destroy and reclaim\n");
#else
    (void)hdir;
    std::printf("SKIP hp real 1 GB pages (not Linux)\n");
#endif
}

} // namespace

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: %s <indexed-prefix>\n", argv[0]);
        return 2;
    }
    char resolved[PATH_MAX];
    if (realpath(argv[1], resolved) == NULL) {
        std::fprintf(stderr, "cannot resolve %s: %s\n", argv[1], std::strerror(errno));
        return 2;
    }
    g_source = resolved;
    try {
        check_lock_not_held_during_pack();
        check_wait_then_reclaim();
        check_interrupted_stage();
        check_pack_failure();
        check_destroy_during_stage();
        check_restage_noop();
        check_wait_deadline();
        check_not_ours();
        check_bad_and_undersized();
        check_two_entries();
        check_capacity();
        check_capacity_allocated();
        check_early_exit();
        check_huge_round_trip();
        check_huge_in_progress();
        check_huge_interrupted();
        check_huge_dead_stub_without_file();
        check_huge_orphan_file();
        check_huge_reserve_fallback();
        check_huge_prefault_fallback();
        check_huge_pack_failure();
        check_huge_destroy_during_stage();
        check_huge_not_ours();
        check_huge_stub_safety();
        check_huge_stub_owner();
        check_huge_dangling();
        check_huge_capacity();
        check_huge_real_pool();
        std::printf("shm_stage_test: OK\n");
        return 0;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "%s (errno=%d %s)\n", e.what(), errno, std::strerror(errno));
        return 1;
    }
}
